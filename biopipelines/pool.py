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
    from .datastream import DataStream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream


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
                 streams: Optional[List[str]] = None,
                 **kwargs):
        """
        Args:
            runs: List of two or more StandardizedOutput objects. With the
                default ``streams=None`` all inputs must expose identical
                stream-name sets (typically because they come from runs of the
                same tool). Pass ``streams=[...]`` to pool only those named
                streams when inputs differ (e.g. one structure carries a
                ``compounds`` stream and another doesn't).
            recount_prefix: If set, replace the default per-run id suffix
                (`<orig_id>_<pool_idx>`) with a flat 1-based renumber across
                all rows of the pool, producing ids `<recount_prefix>_1`,
                `<recount_prefix>_2`, ..., `<recount_prefix>_N`. Counted at
                config time when every input run has fully-resolved ids;
                otherwise the framework emits a lazy `<recount_prefix>_[<N>]`
                pattern that pipe_pool resolves at runtime.
            streams: Stream names to pool. ``None`` (default) is strict: every
                input must expose the same stream set. A list keeps exactly
                those streams (each must exist in every input, with matching
                format/shape) and drops the rest. Tables are always pooled as
                the intersection across inputs, regardless of this setting.
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

        run_streams = [self._stream_names(r) for r in self.runs]
        run_tables = [self._table_names(r) for r in self.runs]

        if streams is None:
            # Strict: every input must expose the same stream set.
            first_streams = run_streams[0]
            for i, other_streams in enumerate(run_streams[1:], start=2):
                if other_streams != first_streams:
                    raise ValueError(
                        f"Pool input #{i} has stream names {sorted(other_streams)} "
                        f"but input #1 has {sorted(first_streams)}. All inputs "
                        "must come from the same upstream tool, or pass "
                        "streams=[...] to pool only the streams they share."
                    )
        else:
            # Explicit: pool exactly the named streams; each must exist everywhere.
            requested = set(streams)
            for i, other_streams in enumerate(run_streams, start=1):
                missing = requested - other_streams
                if missing:
                    raise ValueError(
                        f"Pool input #{i} is missing requested stream(s) "
                        f"{sorted(missing)} (has {sorted(other_streams)})."
                    )
            first_streams = requested

        # Tables are always pooled as the intersection across inputs.
        first_tables = set.intersection(*run_tables) if run_tables else set()

        # Formats must be COMPATIBLE across runs for every pooled stream. A
        # multi-format stream (e.g. "pdb|cif", the form RCSB uses before the
        # extension is known) is compatible with any run whose tokens overlap;
        # an exact-match check would wrongly reject pooling a "pdb" stream with
        # a "pdb|cif" one. The pooled stream advertises the union of tokens, so
        # per-file renderers/consumers resolve the actual extension at runtime.
        def _tokens(fmt):
            return {t.strip().lower() for t in (fmt or "").split("|") if t.strip()}

        first_formats = {}
        for name in first_streams:
            token_union = set()
            for i, r in enumerate(self.runs, start=1):
                toks = _tokens(self._stream_formats(r).get(name))
                if token_union and toks and not (token_union & toks):
                    raise ValueError(
                        f"Pool input #{i} stream '{name}' has format "
                        f"{self._stream_formats(r).get(name)!r}, incompatible "
                        f"with {'|'.join(sorted(token_union))!r} from earlier inputs."
                    )
                token_union |= toks
            first_formats[name] = "|".join(sorted(token_union))

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

        # Content-bearing streams: a stream OWNS a content table when its
        # map_table is the same file as the TableInfo of the same name. A
        # stream that merely reuses a sibling stream's content CSV as a
        # metadata lookup (e.g. Sequence's `fasta` reusing `sequences.csv`)
        # is NOT content-bearing — matching by table name (not bare path)
        # keeps that many-to-one reuse a lineage-only stream.
        self._content_streams: set = set()
        for name in first_streams:
            if self._is_content_stream(self.runs[0], name):
                self._content_streams.add(name)

        # Classification must agree across all runs.
        for i, r in enumerate(self.runs[1:], start=2):
            for name in first_streams:
                is_content = self._is_content_stream(r, name)
                if is_content != (name in self._content_streams):
                    raise ValueError(
                        f"Pool stream '{name}' is content-bearing in some "
                        f"runs but not others (input #{i} disagrees with "
                        "input #1). All inputs must come from the same "
                        "upstream tool."
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

    @staticmethod
    def _is_content_stream(run: StandardizedOutput, name: str) -> bool:
        """A stream is content-bearing when it OWNS the content table of the
        same name: a TableInfo named ``name`` whose path is the stream's
        own ``map_table``. Streams that reuse a sibling's content CSV only as
        a metadata lookup (different table name) are lineage-only."""
        ds = run.streams.get(name)
        mt = getattr(ds, "map_table", None) if ds is not None else None
        if not mt:
            return False
        table = run.tables._tables.get(name)
        return table is not None and table.info.path == mt

    def _content_file_col(self, stream_name: str) -> Optional[str]:
        """Name of the file column in a content stream's table, or None when
        the content stream carries no per-id files (value-based). Read from
        the upstream TableInfo so the pooled table keeps the same schema
        (PDB's ``structures`` uses ``file_path``, not the generic ``file``)."""
        cols = self.runs[0].tables._tables.get(stream_name).info.columns or []
        for cand in ("file_path", "file"):
            if cand in cols:
                return cand
        return None

    def _content_path(self, stream_name: str) -> str:
        """Combined-file path for a content-bearing stream, following the
        upstream ``<stream>/<basename>.csv`` convention (reusing the upstream
        map_table's basename)."""
        upstream_mt = self.runs[0].streams.get(stream_name).map_table
        return self.stream_path(stream_name, os.path.basename(upstream_mt))

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
        content_file_cols: Dict[str, str] = {}
        for name in self._shared_streams:
            is_content = name in self._content_streams
            out_streams[name] = {
                "stream_dir": self.stream_folder(name),
                "map_table": self._content_path(name) if is_content else self.stream_map_path(name),
                "format": self._stream_format_map[name],
            }
            if is_content:
                fc = self._content_file_col(name)
                if fc:
                    content_file_cols[name] = fc
        # A content table shares its file with the stream of the same name;
        # other tables get the canonical tables/<name>.csv path.
        out_tables: Dict[str, str] = {}
        for name in self._shared_tables:
            if name in self._content_streams:
                out_tables[name] = self._content_path(name)
            else:
                out_tables[name] = self.table_path(name)

        config = {
            "runs": runs_cfg,
            "shared_streams": self._shared_streams,
            "shared_tables": self._shared_tables,
            "out_streams": out_streams,
            "out_tables": out_tables,
            "content_bearing_streams": sorted(self._content_streams),
            "content_file_cols": content_file_cols,
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
            # Content streams share one file with their TableInfo; others use
            # the lineage-only <stream>_map.csv.
            is_content = name in self._content_streams
            map_path = self._content_path(name) if is_content else self.stream_map_path(name)

            # Shared-file streams: every run owns one artifact. We slice each
            # run (with the per-run rename to composite ids) and merge them
            # at execution time; here we just predict the combined path.
            any_shared = any(
                getattr(r.streams.get(name), "is_shared_file", False)
                for r in self.runs
                if r.streams.get(name) is not None
            )

            new_ids: List[str] = []
            lazy_pattern_ids: List[str] = []

            for pool_idx, r in enumerate(self.runs, start=1):
                ds = r.streams.get(name)
                orig_ids = list(ds.ids_expanded) if hasattr(ds, "ids_expanded") and ds.ids_expanded else list(ds.ids)
                for oid in orig_ids:
                    if recount and not any_lazy:
                        composite = f"{self.recount_prefix}_{len(new_ids) + 1}"
                    else:
                        composite = f"{oid}_{pool_idx}"
                    new_ids.append(composite)

            # In recount mode with lazy upstream ids, emit a single lazy
            # pattern; pipe_pool fills the count at runtime.
            if recount and any_lazy:
                lazy_pattern_ids = [f"{self.recount_prefix}_[<N>]"]

            # No config-time map: pipe_pool writes every map_table at runtime
            # (declarative get_output_files, per the Map Table Contract).

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
            # Content tables (name matches a content-bearing stream) share the
            # one file with that stream; others get tables/<name>.csv.
            if name in self._content_streams:
                table_path = self._content_path(name)
            else:
                table_path = self.table_path(name)
            tables[name] = TableInfo(
                name=name,
                path=table_path,
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
