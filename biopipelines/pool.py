# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Pool tool for gathering N StandardizedOutputs from parallel runs of the
same upstream tool into one combined StandardizedOutput.

Typical use case::

    runs = []
    with Parallel():
        for i in range(10):
            Resources(gpu="A100")
            runs.append(RFdiffusion(contigs="A1-100", num_designs=10))
    Resources(...)
    combined = Pool(runs=runs)   # 100 dense ids; pool.path column tracks origin

Pool is type-agnostic: it iterates ``runs[0].streams.items()`` and treats
every DataStream uniformly, regardless of whether it carries structures,
sequences, compounds, or anything else. Same for tables. All inputs must
expose the same set of stream names and the same set of table names; a
mismatch raises with a descriptive message.

ID renumbering: by default, outputs are ``<orig_id>_<pool_idx>`` (1-based)
where ``pool_idx`` is the source run's position in ``Pool(runs=runs)``.
Identical original ids across runs (the typical case — same tool, same
parameters) no longer collide. The ``pool.path`` column in every emitted
map_table records the same pool index for downstream queries.

Pass ``recount_prefix="name"`` to switch to a flat 1-based renumber across
the combined rows: ids become ``name_1``, ``name_2``, ..., ``name_N``,
the original ids are kept in an ``original.id`` column, and ``pool.path``
still records the source run. Counted at config time when every input run
has fully-resolved ids; otherwise the framework emits a lazy
``name_[<N>]`` pattern that pipe_pool resolves at runtime.
"""

import os
import json
from typing import Dict, List, Any, Optional

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream, create_map_table
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream, create_map_table


class Pool(BaseConfig):
    """
    Pool tool for gathering N parallel-run StandardizedOutputs into one.

    Concatenates every shared stream and table across the inputs with
    composite ids ``<orig_id>_<pool_idx>`` (1-based) and a ``pool.path``
    provenance column tracking each id's source run.
    """

    TOOL_NAME = "Pool"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== Pool ==="
echo "Uses biopipelines environment (no additional installation needed)."
touch "$INSTALL_SUCCESS"
echo "=== Pool ready ==="
"""

    pool_config_json = Path(lambda self: self.configuration_path("pool_config.json"))
    pool_py = Path(lambda self: self.pipe_script_path("pipe_pool.py"))

    def __init__(self, runs: List[StandardizedOutput],
                 recount_prefix: Optional[str] = None,
                 **kwargs):
        """
        Args:
            runs: List of two or more StandardizedOutput objects, all from
                runs of the same upstream tool. Must expose identical
                stream-name sets and identical table-name sets, with
                matching formats per shared stream.
            recount_prefix: If set, replace the default per-run id suffix
                (`<orig_id>_<pool_idx>`) with a flat 1-based renumber across
                all rows of the pool, producing ids `<recount_prefix>_1`,
                `<recount_prefix>_2`, ..., `<recount_prefix>_N`. Counted at
                config time when every input run has fully-resolved ids;
                otherwise the framework emits a lazy `<recount_prefix>_[<N>]`
                pattern that pipe_pool resolves at runtime.
        """
        if not isinstance(runs, (list, tuple)):
            raise ValueError(
                f"Pool 'runs' must be a list of StandardizedOutput, got "
                f"{type(runs).__name__}. Pass the runs as a list, e.g. "
                "Pool(runs=[run_a, run_b]) or Pool(runs=runs_list)."
            )
        if len(runs) < 2:
            raise ValueError(
                f"Pool requires at least 2 inputs, got {len(runs)}. "
                "Pooling a single run is a no-op; pass the run through "
                "directly or use Panda for filtering."
            )
        for i, r in enumerate(runs):
            if not isinstance(r, StandardizedOutput):
                raise ValueError(
                    f"Pool input #{i + 1} is {type(r).__name__}, expected "
                    "StandardizedOutput. Pool only accepts the output of "
                    "tool calls (e.g., the value returned by RFdiffusion(...))."
                )

        self.recount_prefix = recount_prefix
        self.runs = list(runs)

        # Validate identical stream / table names + matching formats.
        first = self.runs[0]
        first_streams = self._stream_names(first)
        first_tables = self._table_names(first)
        first_formats = self._stream_formats(first)
        for i, r in enumerate(self.runs[1:], start=2):
            other_streams = self._stream_names(r)
            if other_streams != first_streams:
                raise ValueError(
                    f"Pool input #{i} has stream names {sorted(other_streams)} "
                    f"but input #1 has {sorted(first_streams)}. All inputs "
                    "must come from the same upstream tool."
                )
            other_tables = self._table_names(r)
            if other_tables != first_tables:
                raise ValueError(
                    f"Pool input #{i} has table names {sorted(other_tables)} "
                    f"but input #1 has {sorted(first_tables)}. All inputs "
                    "must come from the same upstream tool."
                )
            other_formats = self._stream_formats(r)
            for name, fmt in first_formats.items():
                if other_formats.get(name) != fmt:
                    raise ValueError(
                        f"Pool input #{i} stream '{name}' has format "
                        f"{other_formats.get(name)!r} but input #1 has {fmt!r}."
                    )

        # Validate that every run uses the same storage shape (shared-file
        # vs per-id) for a given stream. Mixed shapes can't be combined
        # coherently — slice+merge requires shared sources on every side,
        # and per-id copy requires per-id sources. Catch the mismatch here
        # rather than silently dropping the odd-one-out at pipe_pool time.
        for name in first_streams:
            shapes = [
                bool(getattr(r.streams.get(name), "is_shared_file", False))
                for r in self.runs
                if r.streams.get(name) is not None
            ]
            if len(set(shapes)) > 1:
                shared_idxs = [i for i, s in enumerate(shapes, start=1) if s]
                perid_idxs = [i for i, s in enumerate(shapes, start=1) if not s]
                raise ValueError(
                    f"Pool stream '{name}' has mixed storage across runs: "
                    f"shared-file in runs {shared_idxs}, per-id in runs "
                    f"{perid_idxs}. All runs must agree."
                )

        # Cache per-stream / per-table info for use in get_output_files.
        self._shared_streams = sorted(first_streams)
        self._shared_tables = sorted(first_tables)
        self._stream_format_map = first_formats

        super().__init__(**kwargs)

    # ── helpers: enumerate streams / tables of an upstream output ────────────

    @staticmethod
    def _stream_names(run: StandardizedOutput) -> set:
        """Names of non-empty DataStreams on a StandardizedOutput."""
        names = set()
        for k, v in run.streams.items():
            if isinstance(v, DataStream) and len(v) > 0:
                names.add(k)
        return names

    @staticmethod
    def _stream_formats(run: StandardizedOutput) -> Dict[str, str]:
        """Map stream name -> declared format for non-empty streams."""
        out: Dict[str, str] = {}
        for k, v in run.streams.items():
            if isinstance(v, DataStream) and len(v) > 0:
                out[k] = v.format
        return out

    @staticmethod
    def _table_names(run: StandardizedOutput) -> set:
        """Names of TableInfo entries on a StandardizedOutput."""
        if hasattr(run, 'tables') and hasattr(run.tables, '_tables'):
            return set(run.tables._tables.keys())
        return set()

    # ── BaseConfig hooks ─────────────────────────────────────────────────────

    def validate_params(self):
        """All validation is in __init__; nothing additional at config time."""
        pass

    def configure_inputs(self, pipeline_folders):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"INPUTS: {len(self.runs)} runs")
        lines.append(f"STREAMS: {', '.join(self._shared_streams)}")
        if self._shared_tables:
            lines.append(f"TABLES: {', '.join(self._shared_tables)}")
        return lines

    def generate_script(self, script_path: str) -> str:
        """Emit a thin bash script that hands off to pipe_pool.py."""
        # Per-input run config: serialise enough metadata that pipe_pool.py
        # can read each run's streams + tables without re-importing the
        # source tool's Python class.
        runs_cfg: List[Dict[str, Any]] = []
        for r in self.runs:
            run_streams: List[Dict[str, Any]] = []
            for name in self._shared_streams:
                ds = r.streams.get(name)
                # Preserve shared-file form (str) — wrapping a path string in
                # list() would shatter it into a list of characters.
                files_payload = ds.files if ds.is_shared_file else list(ds.files)
                run_streams.append({
                    "name": ds.name,
                    "stream_key": name,
                    "format": ds.format,
                    "ids": list(ds.ids),
                    "files": files_payload,
                    "map_table": ds.map_table,
                    "is_shared_file": bool(ds.is_shared_file),
                })
            run_tables: List[Dict[str, Any]] = []
            for name in self._shared_tables:
                t = r.tables._tables[name]
                run_tables.append({
                    "name": name,
                    "path": t.info.path,
                    "columns": list(t.info.columns) if t.info.columns else [],
                })
            runs_cfg.append({"streams": run_streams, "tables": run_tables})

        # Pre-compute the per-stream output stream dirs and per-table output
        # paths so pipe_pool.py doesn't need to know about FolderManager.
        out_streams: Dict[str, Dict[str, str]] = {}
        for name in self._shared_streams:
            out_streams[name] = {
                "stream_dir": self.stream_folder(name),
                "map_table": self.stream_map_path(name),
                "format": self._stream_format_map[name],
            }
        out_tables: Dict[str, str] = {
            name: self.table_path(name) for name in self._shared_tables
        }

        config = {
            "runs": runs_cfg,
            "shared_streams": self._shared_streams,
            "shared_tables": self._shared_tables,
            "out_streams": out_streams,
            "out_tables": out_tables,
            "output_folder": self.output_folder,
            "recount_prefix": self.recount_prefix,
        }

        with open(self.pool_config_json, 'w') as f:
            json.dump(config, f, indent=2)

        script = "#!/bin/bash\n"
        script += "# Pool: gather N upstream runs into one combined output\n\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += 'echo "=== Pool: gathering runs ==="\n'
        script += f'python "{self.pool_py}" "{self.pool_config_json}"\n'
        script += 'echo "=== Pool complete ==="\n\n'
        script += self.generate_completion_check_footer()
        return script

    def get_output_files(self) -> Dict[str, Any]:
        """Build the combined StandardizedOutput shape.

        Each shared stream gathers ids of the form ``<orig_id>_<pool_idx>``
        across runs (1-based pool index). Each shared table is emitted
        with the matching renumbered ids. The pool.path column is
        materialised at runtime by pipe_pool.py.
        """
        output: Dict[str, Any] = {}

        # Decide once whether any input stream carries lazy ids. If yes,
        # config-time counting for recount mode is impossible — we have to
        # emit a lazy `<prefix>_[<N>]` pattern and let pipe_pool resolve.
        any_lazy = any(
            getattr(r.streams.get(name), "is_lazy", False)
            for name in self._shared_streams
            for r in self.runs
            if r.streams.get(name) is not None
        )
        recount = self.recount_prefix is not None

        for name in self._shared_streams:
            fmt = self._stream_format_map[name]
            stream_dir = self.stream_folder(name)
            map_path = self.stream_map_path(name)

            # Shared-file streams: every run owns one artifact. We slice each
            # run (with the per-run rename to composite ids) and merge them
            # at execution time; here we just predict the combined path.
            any_shared = any(
                getattr(r.streams.get(name), "is_shared_file", False)
                for r in self.runs
                if r.streams.get(name) is not None
            )

            new_ids: List[str] = []
            new_files: List[str] = []
            pool_paths: List[int] = []
            # Original ids preserved per row for the map_table when
            # recounting, so the renumber stays auditable.
            orig_id_col: List[str] = []
            # Only populated in the lazy-recount path; placeholder otherwise
            # so the DataStream construction below can reference it freely.
            lazy_pattern_ids: List[str] = []

            for pool_idx, r in enumerate(self.runs, start=1):
                ds = r.streams.get(name)
                orig_ids = list(ds.ids_expanded) if hasattr(ds, "ids_expanded") and ds.ids_expanded else list(ds.ids)
                # Shared-file streams: files_expanded would return a single-path
                # list but len != len(orig_ids), and we don't want it threaded
                # through the per-id append branch below anyway. Skip files here.
                if getattr(ds, "is_shared_file", False):
                    orig_files = []
                else:
                    orig_files = list(ds.files_expanded) if hasattr(ds, "files_expanded") and ds.files_expanded else list(ds.files)
                for j, oid in enumerate(orig_ids):
                    if recount and not any_lazy:
                        composite = f"{self.recount_prefix}_{len(new_ids) + 1}"
                    else:
                        composite = f"{oid}_{pool_idx}"
                    new_ids.append(composite)
                    pool_paths.append(pool_idx)
                    orig_id_col.append(oid)
                    if orig_files and len(orig_files) == len(orig_ids):
                        src = orig_files[j] or ""
                        # Skip wildcard extensions ('.*') — the upstream
                        # tool didn't resolve a concrete extension at
                        # config time. pipe_pool.py picks the real one
                        # at runtime by reading the upstream map_table.
                        if src.endswith(".*") or "*" in os.path.basename(src):
                            ext = ""
                        else:
                            ext = os.path.splitext(src)[1]
                        if ext:
                            new_files.append(os.path.join(stream_dir, f"{composite}{ext}"))

            # In recount mode with lazy upstream ids, drop the per-row id
            # list and emit a single lazy pattern. pipe_pool fills the
            # actual count and rewrites the map at runtime.
            if recount and any_lazy:
                lazy_pattern_ids = [f"{self.recount_prefix}_[<N>]"]

            # Always materialise the map_table at config time with the
            # composite ids and pool.path column. pipe_pool.py overwrites
            # it at runtime once the actual files have been copied (so
            # post-run inspection sees the same shape as config-time
            # prediction). If any composite id has no resolved extension,
            # drop the files list entirely — runtime fills it in from the
            # upstream map_tables. With lazy recount we skip the config-time
            # map (no concrete ids to write).
            if recount and any_lazy:
                pass
            else:
                files_for_map = new_files if (new_files and len(new_files) == len(new_ids)) else None
                additional = {"pool.path": pool_paths}
                if recount:
                    additional["original.id"] = orig_id_col
                create_map_table(
                    map_path,
                    ids=new_ids,
                    files=files_for_map,
                    additional_columns=additional,
                )

            # The DataStream's `files` declaration drives the framework's
            # completion-check glob. For file-bearing streams we use the
            # wildcard '<id>.*' form because the runtime extension can
            # disagree with whatever the upstream tool declared at config
            # time (Ligand declares `.csv` for its compound table but
            # renders `.pdb`; PDB declares `.*` outright). The wildcard
            # makes the completion check tolerant of the actual extension
            # pipe_pool.py picked. For value-only streams (no per-id
            # files upstream), keep `files=[]` so the completion check
            # doesn't expect anything that won't be on disk.
            any_input_has_files = any(
                bool(r.streams.get(name).files) for r in self.runs
                if r.streams.get(name) is not None
            )
            ids_for_stream = lazy_pattern_ids if (recount and any_lazy) else new_ids
            # Shared-file streams: emit a single combined artifact under the
            # stream folder. Pick the first run's basename as canonical; pipe
            # pool slices each run with rename + merges into this path.
            if any_shared:
                first_shared = next(
                    (r.streams.get(name).files for r in self.runs
                     if r.streams.get(name) is not None
                     and getattr(r.streams.get(name), "is_shared_file", False)),
                    None,
                )
                if first_shared:
                    stream_files = os.path.join(stream_dir, os.path.basename(first_shared))
                else:
                    stream_files = []
            elif ids_for_stream and any_input_has_files:
                stream_files = [os.path.join(stream_dir, "<id>.*")]
            else:
                stream_files = []

            output[name] = DataStream(
                name=name,
                ids=ids_for_stream,
                files=stream_files,
                map_table=map_path,
                format=fmt,
            )

        tables: Dict[str, TableInfo] = {}
        for name in self._shared_tables:
            # Carry forward the column schema from input #1; pipe_pool.py
            # appends "pool.path" to the actual rows at runtime.
            first_t = self.runs[0].tables._tables[name]
            cols = list(first_t.info.columns) if first_t.info.columns else []
            if "pool.path" not in cols:
                cols.append("pool.path")
            tables[name] = TableInfo(
                name=name,
                path=self.table_path(name),
                columns=cols,
                description=f"Pooled {name} from {len(self.runs)} runs"
            )
        if tables:
            output["tables"] = tables

        output["output_folder"] = self.output_folder
        return output

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                "n_runs": len(self.runs),
                "shared_streams": self._shared_streams,
                "shared_tables": self._shared_tables,
            }
        })
        return base_dict
