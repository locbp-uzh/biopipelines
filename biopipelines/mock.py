# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Mock tool only meant to be used for testing the functionalities of BioPipeliens.
It emits stub outputs (empty files, value streams, CSV tables) for
pipeline-wiring tests and local dry-runs. Follows the real tool contract:
compact IDs + file templates at config time, map_table materialized at
runtime (by default) or at config time when requested.

Design:
    - `source`: one or more DataStreams (wrap in Bundle(...) / Each(...) to
      control combinatorics, same as real tools like Boltz2). `source` is
      mutually exclusive with `ids`.
    - `ids`: explicit output IDs when there is no upstream source.
    - `children`: optional suffix pattern appended to each parent ID.
        - Deterministic (e.g. "<1..3>", "<A B>"): output IDs stay expandable;
          config-time count is exact.
        - Lazy (e.g. "[_<N><A V>]"): output IDs remain lazy; `produce` is
          REQUIRED so the runtime knows what files to actually emit per
          parent.
    - `produce`: list of expanded suffix strings applied per parent at
      runtime (required whenever `children` is lazy).
    - `streams`: dict[name -> {format, file, values}].
    - `tables`:  dict[name -> {columns, fill, rows}].
    - `map_table_strategy`: "runtime" (default) | "config" | "both".
    - `missing`: IDs pretended to have failed (applied post-expansion).

Examples:

    # Simple stub with explicit IDs and a per-ID file stream.
    Mock(
        ids=["a", "b", "c"],
        streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
        tables={"scores": {"columns": ["score"], "fill": {"score": 0.0}}},
    )

    # Fan-out — each parent ID produces 3 deterministic children.
    Mock(
        source=parent_tool.streams.structures,
        children="<1..3>",
        streams={"designs": {"format": "pdb", "file": "<id>.pdb"}},
    )

    # Lazy fan-out — parents spawn an unknown-at-config-time set of children,
    # but the test declares exactly what the runtime should emit.
    Mock(
        source=parent_tool.streams.structures,
        children="[_<N><A V>]",
        produce=["_1A", "_1V", "_2A"],
        streams={"mutants": {"format": "pdb", "file": "<id>.pdb"}},
    )

    # Combinatorics via Bundle/Each (same wrappers as real tools).
    Mock(
        source=[Each(lib_a), Each(lib_b)],
        streams={"predictions": {"format": "pdb", "file": "<id>.pdb"}},
    )
"""

import json
import os
from typing import Any, Dict, List, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream, create_map_table
    from .combinatorics import (
        Bundle,
        Each,
        predict_output_ids_with_provenance,
    )
    from . import id_patterns
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream, create_map_table
    from combinatorics import Bundle, Each, predict_output_ids_with_provenance
    import id_patterns


MapTableStrategy = str  # "runtime" | "config" | "both"


def _collect_datastreams(source: Any) -> List[DataStream]:
    """Flatten Bundle/Each/list into a flat DataStream list."""
    if source is None:
        return []
    if isinstance(source, DataStream):
        return [source]
    if isinstance(source, (Bundle, Each)):
        out = []
        for s in source.sources:
            out.extend(_collect_datastreams(s))
        return out
    if isinstance(source, list):
        out = []
        for s in source:
            out.extend(_collect_datastreams(s))
        return out
    raise ValueError(
        f"Mock source must be DataStream, Bundle, Each, or list thereof; got {type(source)}"
    )


class Mock(BaseConfig):
    """Stub tool for tests. See module docstring for usage."""

    TOOL_NAME = "Mock"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== Mock ==="
echo "No installation required."
touch "$INSTALL_SUCCESS"
echo "=== Mock ready ==="
"""

    config_file = Path(lambda self: self.configuration_path("mock_config.json"))
    mock_py = Path(lambda self: self.pipe_script_path("pipe_mock.py"))

    def __init__(self,
                 source: Optional[Union[DataStream, Bundle, Each, List]] = None,
                 ids: Optional[Union[str, List[str]]] = None,
                 children: Optional[str] = None,
                 produce: Optional[List[str]] = None,
                 streams: Optional[Dict[str, Dict[str, Any]]] = None,
                 tables: Optional[Dict[str, Dict[str, Any]]] = None,
                 map_table_strategy: MapTableStrategy = "runtime",
                 missing: Optional[List[str]] = None,
                 **kwargs):
        # ── mutex + basic validation ──────────────────────────────────────────
        if source is None and ids is None:
            raise ValueError("Mock requires either 'source' or 'ids'")
        if source is not None and ids is not None:
            raise ValueError("Mock: 'source' and 'ids' are mutually exclusive")

        if map_table_strategy not in ("runtime", "config", "both"):
            raise ValueError(
                f"Mock: map_table_strategy must be 'runtime'|'config'|'both', got {map_table_strategy!r}"
            )

        if children is not None and id_patterns.is_lazy(children) and not produce:
            raise ValueError(
                "Mock: 'children' is a lazy pattern; you must also pass 'produce' "
                "with the suffixes the runtime should emit per parent"
            )
        if produce is not None and children is None:
            raise ValueError("Mock: 'produce' only makes sense together with 'children'")

        # ── resolve parent ids + provenance ───────────────────────────────────
        self.source_raw = source
        self.source_streams = _collect_datastreams(source)

        if ids is not None:
            self.parent_ids: List[str] = [ids] if isinstance(ids, str) else list(ids)
            self.provenance: Dict[str, List[str]] = {}
            self.axis_names: List[str] = []
        else:
            parent_ids, provenance, axis_names = self._resolve_from_source(source)
            self.parent_ids = parent_ids
            self.provenance = provenance
            self.axis_names = axis_names

        if not self.parent_ids:
            raise ValueError("Mock: resolved parent ID list is empty")

        # ── compact output IDs (keep patterns; NEVER expand) ──────────────────
        self.children = children
        self.produce = list(produce) if produce else None
        if children is None:
            self.output_ids = list(self.parent_ids)
        else:
            self.output_ids = id_patterns.append_suffix(self.parent_ids, children)

        # ── stream/table specs ────────────────────────────────────────────────
        self.stream_specs: Dict[str, Dict[str, Any]] = streams or {}
        self.table_specs: Dict[str, Dict[str, Any]] = tables or {}
        self.map_table_strategy = map_table_strategy
        self.missing: List[str] = list(missing) if missing else []

        self._validate_specs()

        super().__init__(**kwargs)

    # ── helpers ───────────────────────────────────────────────────────────────

    @staticmethod
    def _resolve_from_source(source):
        """Compute parent IDs + provenance from a DataStream / Bundle / Each / list.

        Uses `predict_output_ids_with_provenance` so Mock's axis handling
        matches the real tool contract. DataStreams are wrapped in
        StandardizedOutput(s) so the combinatorics machinery can reach them
        via `.streams.<stream_name>`.
        """
        from .base_config import StandardizedOutput

        def _wrap(value, stream_name):
            """Replace bare DataStreams with StandardizedOutputs keyed by stream_name."""
            if isinstance(value, DataStream):
                return StandardizedOutput({stream_name: value})
            if isinstance(value, Bundle):
                return Bundle(*[_wrap(s, stream_name) for s in value.sources])
            if isinstance(value, Each):
                return Each(*[_wrap(s, stream_name) for s in value.sources])
            return value

        named = {}

        def _axis(value):
            inner = _collect_datastreams(value)
            if not inner:
                raise ValueError("Mock: source axis has no DataStreams")
            axis_name = inner[0].name or f"axis{len(named)}"
            wrapped = _wrap(value, axis_name)
            named[axis_name] = (wrapped, axis_name)

        if isinstance(source, (DataStream, Bundle, Each)):
            _axis(source)
        elif isinstance(source, list):
            for s in source:
                _axis(s)
        else:
            raise ValueError(f"Mock: unsupported source type {type(source)}")

        parent_ids, provenance = predict_output_ids_with_provenance(**named)
        return parent_ids, provenance, list(named.keys())

    def _validate_specs(self):
        n = len(self.parent_ids)
        for name, spec in self.stream_specs.items():
            vals = spec.get("values")
            if vals is not None and len(vals) != n:
                raise ValueError(
                    f"Mock stream '{name}': values has {len(vals)} entries "
                    f"but {n} parent IDs were declared"
                )
        for name, spec in self.table_specs.items():
            rows = spec.get("rows")
            if rows is not None and len(rows) != n:
                raise ValueError(
                    f"Mock table '{name}': rows has {len(rows)} entries but "
                    f"{n} parent IDs were declared"
                )

    def validate_params(self):
        pass

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_id_provenance(self) -> Dict[str, Dict[str, str]]:
        """Synthesize provenance from in-memory config data.

        The default implementation reads `{axis}.id` columns from the
        stream map_tables — but with `map_table_strategy="runtime"` those
        files don't exist until the pipeline actually runs. The lineage
        writer runs at `pipeline.save()` time (config time), so we build
        the same {axis: {child_id: parent_id}} shape from `self.provenance`
        (upstream axes) + the children fan-out (when present).

        After execution, runtime map_tables carry the ground-truth IDs
        (expanded patterns, dropped `missing` rows, etc.), so we defer to
        the superclass reader whenever any map_table is present on disk."""
        for name, spec in self.stream_specs.items():
            if not (spec.get("file") or spec.get("values")):
                continue
            if os.path.exists(self.stream_map_path(name)):
                runtime = super().get_id_provenance()
                if runtime:
                    return runtime
                break

        result: Dict[str, Dict[str, str]] = {}

        # Upstream axes: provenance is indexed per parent_id; project onto
        # the compact output_ids so downstream tracers can match.
        for axis, prov_list in (self.provenance or {}).items():
            if not prov_list:
                continue
            mapping: Dict[str, str] = {}
            for out_id, parent_id in zip(self.output_ids, prov_list):
                if parent_id:
                    mapping[out_id] = parent_id
            if mapping:
                result[axis] = mapping

        # Fan-out axis: when `children` is set, add a self-named link from
        # each output_id back to the parent it grew from. Use the first
        # stream name as the axis label (matches the runtime column).
        if self.children and self.stream_specs:
            axis_name = next(iter(self.stream_specs.keys()))
            mapping = {}
            for out_id, parent_id in zip(self.output_ids, self.parent_ids):
                mapping[out_id] = parent_id
            if mapping and axis_name not in result:
                result[axis_name] = mapping

        # Source passthrough (no children, but has a source): identity
        # provenance so a single-hop chain connects visually.
        if not self.children and self.source_streams and self.axis_names:
            for axis in self.axis_names:
                if axis in result:
                    continue
                # parent_ids == output_ids in passthrough; still useful for
                # the lineage tracer since axis names link tools.
                result[axis] = {oid: oid for oid in self.output_ids}

        return result

    # ── path helpers ──────────────────────────────────────────────────────────

    def _stream_file_paths(self, name: str, spec: Dict[str, Any]) -> List[str]:
        """Mock-specific: resolve a ``streams={...}`` spec's ``file`` template
        into an absolute path under the stream's folder. Other tools build
        their own file lists in ``get_output_files()`` — this helper only
        exists because Mock's test-facing spec dict carries a template
        string rather than a list of paths."""
        template = spec.get("file")
        if not template:
            return []
        return [self.stream_path(name, template)]

    def _table_columns(self, spec: Dict[str, Any]) -> List[str]:
        cols = list(spec.get("columns") or [])
        if "id" not in cols:
            cols = ["id"] + cols
        return cols

    # ── output wiring ─────────────────────────────────────────────────────────

    def get_output_files(self) -> Dict[str, Any]:
        output: Dict[str, Any] = {}

        # Pre-write map_tables at config time when requested. For lazy outputs
        # we only write what can be expanded (deterministic prefix); the pipe
        # script overwrites with the full table at runtime.
        for name, spec in self.stream_specs.items():
            fmt = spec.get("format", "txt")
            files = self._stream_file_paths(name, spec)
            map_path = self.stream_map_path(name) if (files or spec.get("values")) else ""

            if map_path and self.map_table_strategy in ("config", "both"):
                # create_map_table() handles mkdir for the target path.
                prov = self.provenance if self.provenance else None
                kwargs = {"files": files or None}
                if spec.get("values"):
                    kwargs["values"] = list(spec["values"])
                if prov:
                    kwargs["provenance"] = prov
                create_map_table(map_path, ids=list(self.output_ids), **kwargs)

            output[name] = DataStream(
                name=name,
                ids=list(self.output_ids),
                files=files,
                map_table=map_path,
                format=fmt,
            )

        tables = {}
        for name, spec in self.table_specs.items():
            tables[name] = TableInfo(
                name=name,
                path=self.table_path(name),
                columns=self._table_columns(spec),
                description=f"Mock table '{name}'"
            )
        if tables:
            output["tables"] = tables

        output["output_folder"] = self.output_folder
        return output

    # ── script emission ───────────────────────────────────────────────────────

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"PARENT IDS: {len(self.parent_ids)}")
        if self.children:
            kind = "lazy" if id_patterns.is_lazy(self.children) else "deterministic"
            produce = f", produce={self.produce}" if self.produce else ""
            lines.append(f"CHILDREN: {self.children} ({kind}{produce})")
        if self.stream_specs:
            lines.append(f"STREAMS: {', '.join(self.stream_specs.keys())}")
        if self.table_specs:
            lines.append(f"TABLES: {', '.join(self.table_specs.keys())}")
        lines.append(f"MAP_TABLE STRATEGY: {self.map_table_strategy}")
        if self.missing:
            lines.append(f"MISSING: {self.missing}")
        return lines

    def generate_script(self, script_path: str) -> str:
        streams_cfg = {}
        for name, spec in self.stream_specs.items():
            streams_cfg[name] = {
                "format": spec.get("format", "txt"),
                "file": spec.get("file"),
                "values": list(spec["values"]) if spec.get("values") else None,
                "map_table": self.stream_map_path(name)
                if (spec.get("file") or spec.get("values"))
                else None,
                # Per-stream output directory. The pipe script writes
                # <stream_folder>/<file-template-with-id-substituted> here
                # instead of joining templates to output_folder directly.
                "stream_folder": self.stream_folder(name),
            }

        tables_cfg = {}
        for name, spec in self.table_specs.items():
            tables_cfg[name] = {
                "path": self.table_path(name),
                "columns": self._table_columns(spec),
                "fill": spec.get("fill"),
                "rows": spec.get("rows"),
            }

        source_stream_jsons = []
        if self.source_streams:
            for idx, ds in enumerate(self.source_streams):
                json_path = self.configuration_path(f"source_{idx}_{ds.name or 'stream'}.json")
                ds.save_json(json_path)
                source_stream_jsons.append({"name": ds.name, "path": json_path})

        config_data = {
            "output_folder": self.output_folder,
            "parent_ids": self.parent_ids,
            "output_ids": self.output_ids,
            "provenance": self.provenance,
            "axis_names": self.axis_names,
            "children": self.children,
            "produce": self.produce,
            "streams": streams_cfg,
            "tables": tables_cfg,
            "map_table_strategy": self.map_table_strategy,
            "missing": self.missing,
            "source_streams": source_stream_jsons,
        }
        with open(self.config_file, "w") as f:
            json.dump(config_data, f, indent=2)

        script = "#!/bin/bash\n# Mock stub-output generator\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += f'echo "=== Mock: generating stub outputs ==="\n'
        script += f'python "{self.mock_py}" "{self.config_file}"\n'
        script += f'echo "=== Mock complete ==="\n'
        script += self.generate_completion_check_footer()
        return script

    def to_dict(self) -> Dict[str, Any]:
        base = super().to_dict()
        base.update({
            "tool_params": {
                "parent_ids": self.parent_ids,
                "output_ids": self.output_ids,
                "children": self.children,
                "produce": self.produce,
                "streams": list(self.stream_specs.keys()),
                "tables": list(self.table_specs.keys()),
                "map_table_strategy": self.map_table_strategy,
                "missing": self.missing,
            }
        })
        return base
