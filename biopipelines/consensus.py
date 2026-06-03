# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Consensus - per-group aggregation of a resi-csv stream.

Collapses an N-file (per-id) ``resi-csv`` stream into per-(chain, resi) aggregate
values within each group, and emits a new ``resi-csv`` stream carrying the
aggregate columns. It is value-agnostic: it does not know about selections,
distances, or any specific column meaning - it groups input ids by ``groups=``
(using the framework's id matching) and applies aggregation operations across
each group, the same way Panda applies table operations.

Output is one file per group, keyed by the group id, so a downstream stream-consuming tool
(e.g. Selection thresholding ``include="frequency>=0.5"``) stays id-keyed and
matches per id exactly. The output stream keeps the input stream's name
(Panda/Pool convention), so ``Consensus(dsel.streams.distances, ...)`` is read
back as ``consensus.streams.distances``.

Example::

    consensus = Consensus(
        dsel.streams.distances,            # resi-csv: id, chain, resi, distance
        groups=parents,                    # stream whose ids define the partition
        operations=[Consensus.fraction("distance<=6", name="frequency"),
                    Consensus.min("distance")],
    )
    # consensus.streams.distances -> resi-csv: id, chain, resi, frequency, distance_min, n_group
    surface = Selection(Selection.add(consensus.streams.distances, include="frequency>=0.5"),
                        Selection.invert(), structures=poses)
"""

import os
import json
from dataclasses import dataclass, field
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput
    from .file_paths import Path
    from .datastream import DataStream
except ImportError:
    import sys
    sys.path.insert(0, os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput
    from file_paths import Path
    from datastream import DataStream


@dataclass
class ConsensusOp:
    """A single aggregation applied per (chain, resi) across a group's ids."""

    type: str
    params: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {"type": self.type, "params": self.params}


class Consensus(BaseConfig):
    """
    Per-group aggregation of a ``resi-csv`` stream.

    Takes a ``resi-csv`` stream (one CSV per id, each row a ``(chain, resi)``
    with value columns), partitions the input ids into groups via ``groups=``,
    and for every ``(chain, resi)`` computes aggregate columns across the ids of
    each group. The result is written back as a ``resi-csv`` stream, one file per
    group keyed by the group id. The output stream keeps the input stream's name
    (Panda/Pool convention).

    Operations (each adds one column to the output rows):
        - ``Consensus.fraction(predicate, name=...)`` - fraction of the group's
          ids whose row at that ``(chain, resi)`` satisfies ``predicate``
          (``"column op value"``, e.g. ``"distance<=6"``). Ids lacking that
          residue count toward the denominator.
        - ``Consensus.mean(column, name=...)`` - mean of ``column`` across the
          group's rows present at that ``(chain, resi)``.
        - ``Consensus.min(column, name=...)`` / ``Consensus.max(column, name=...)``
        - ``Consensus.count(name=...)`` - number of ids in the group with a row
          at that ``(chain, resi)``.

    Every output row also carries ``n_group`` (the group size).
    """

    TOOL_NAME = "Consensus"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== Consensus ==="
echo "Uses biopipelines environment (no additional installation needed)."
touch "$INSTALL_SUCCESS"
echo "=== Consensus ready ==="
"""

    config_json = Path(lambda self: self.configuration_path("config.json"))
    stream_json = Path(lambda self: self.configuration_path(".input_stream.json"))
    groups_json = Path(lambda self: self.configuration_path(".groups.json"))
    consensus_map = Path(lambda self: self.stream_map_path(self.output_stream_name))
    missing_csv = Path(lambda self: self.table_path("missing"))
    consensus_py = Path(lambda self: self.pipe_script_path("pipe_consensus.py"))

    # ── Operation factories ──

    @staticmethod
    def fraction(predicate: str, name: str = "frequency") -> ConsensusOp:
        """Fraction of a group's ids whose row satisfies ``predicate``
        (``"column op value"``, e.g. ``"distance<=6"``)."""
        return ConsensusOp("fraction", {"predicate": predicate, "name": name})

    @staticmethod
    def mean(column: str, name: Optional[str] = None) -> ConsensusOp:
        return ConsensusOp("mean", {"column": column, "name": name or f"{column}_mean"})

    @staticmethod
    def min(column: str, name: Optional[str] = None) -> ConsensusOp:
        return ConsensusOp("min", {"column": column, "name": name or f"{column}_min"})

    @staticmethod
    def max(column: str, name: Optional[str] = None) -> ConsensusOp:
        return ConsensusOp("max", {"column": column, "name": name or f"{column}_max"})

    @staticmethod
    def count(name: str = "count") -> ConsensusOp:
        return ConsensusOp("count", {"name": name})

    # ── Constructor ──

    def __init__(self,
                 stream: Union[DataStream, StandardizedOutput],
                 *,
                 groups: Union[DataStream, StandardizedOutput],
                 operations: List[ConsensusOp],
                 **kwargs):
        # Keep original inputs for upstream missing-table detection.
        self.stream_input = stream
        self.groups_input = groups
        self.input_stream = self._resolve_resi_csv(stream, "stream")
        self.groups_stream = self._resolve_group_stream(groups)
        self.operations: List[ConsensusOp] = list(operations) if operations else []
        # Preserve the input stream's name on the output (Panda/Pool convention).
        self.output_stream_name = self.input_stream.name
        super().__init__(**kwargs)

    @staticmethod
    def _resolve_resi_csv(obj, argname: str) -> DataStream:
        if isinstance(obj, DataStream):
            return obj
        if isinstance(obj, StandardizedOutput):
            for _, ds in obj.streams.items():
                if ds is not None and len(ds) > 0 and ds.format == "resi-csv":
                    return ds
            raise ValueError(f"{argname}: no resi-csv stream found in output")
        raise ValueError(f"{argname} must be a resi-csv DataStream or StandardizedOutput, got {type(obj)}")

    @staticmethod
    def _resolve_group_stream(obj) -> DataStream:
        if isinstance(obj, DataStream):
            return obj
        if isinstance(obj, StandardizedOutput):
            for name in ("sequences", "structures", "compounds"):
                ds = obj.streams.get(name)
                if ds is not None and len(ds) > 0:
                    return ds
            raise ValueError("groups: no non-empty sequences/structures/compounds stream found")
        raise ValueError(f"groups must be a DataStream or StandardizedOutput, got {type(obj)}")

    def validate_params(self):
        if not self.input_stream or len(self.input_stream) == 0:
            raise ValueError("stream is required and must not be empty")
        if self.input_stream.format != "resi-csv":
            raise ValueError(f"stream must be a resi-csv stream, got format '{self.input_stream.format}'")
        if not self.groups_stream or len(self.groups_stream) == 0:
            raise ValueError("groups is required and must not be empty")
        if not self.operations:
            raise ValueError("at least one operation is required")

        seen_names = set()
        for op in self.operations:
            if op.type not in ("fraction", "mean", "min", "max", "count"):
                raise ValueError(f"unknown operation type '{op.type}'")
            name = op.params.get("name")
            if not name or not isinstance(name, str):
                raise ValueError(f"operation '{op.type}' is missing an output column name")
            if name in ("id", "chain", "resi", "n_group"):
                raise ValueError(f"operation name '{name}' collides with a reserved column")
            if name in seen_names:
                raise ValueError(f"duplicate operation output column '{name}'")
            seen_names.add(name)
            if op.type == "fraction" and not op.params.get("predicate"):
                raise ValueError("fraction operation requires a predicate")
            elif op.type in ("mean", "min", "max") and not op.params.get("column"):
                raise ValueError(f"operation '{op.type}' requires a column")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"INPUT STREAM: {len(self.input_stream)} files (resi-csv)")
        lines.append(f"GROUPS: {len(self.groups_stream)} ids")
        op_summaries = []
        for op in self.operations:
            if op.type == "fraction":
                op_summaries.append(f"fraction({op.params['predicate']} -> {op.params['name']})")
            elif op.type == "count":
                op_summaries.append(f"count -> {op.params['name']}")
            else:
                op_summaries.append(f"{op.type}({op.params['column']} -> {op.params['name']})")
        lines.append(f"OPERATIONS: {', '.join(op_summaries)}")
        return lines

    def generate_script(self, script_path: str) -> str:
        self.input_stream.save_json(self.stream_json)
        self.groups_stream.save_json(self.groups_json)

        config_data = {
            "stream_json": self.stream_json,
            "groups_json": self.groups_json,
            "operations": [op.to_dict() for op in self.operations],
            "consensus_dir": self.stream_folder(self.output_stream_name),
            "consensus_map_csv": self.consensus_map,
        }
        with open(self.config_json, "w") as f:
            json.dump(config_data, f, indent=2)

        script = "#!/bin/bash\n"
        script += "# Consensus execution script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += f"""echo "Aggregating {len(self.input_stream)} resi-csv files over {len(self.groups_stream)} groups"

python "{self.consensus_py}" "{self.config_json}"

echo "Consensus written to: {self.stream_folder(self.output_stream_name)}"
"""
        # A group whose members were all dropped upstream produces no <id>.csv;
        # propagate the upstream manifest so the completion check excuses it.
        script += self.generate_missing_propagation(
            self.stream_input, self.groups_input, missing_csv=self.missing_csv
        )
        script += self.generate_completion_check_footer()
        return script

    def get_output_files(self) -> Dict[str, Any]:
        name = self.output_stream_name
        # One file per group, keyed by the group id (matches the groups stream).
        consensus_stream = DataStream(
            name=name,
            ids=self.groups_stream.ids,
            files=[self.stream_path(name, "<id>.csv")],
            map_table=self.consensus_map,
            format="resi-csv",
        )
        tables = {}
        # Excuse groups whose members were all dropped upstream: declare a
        # `missing` table when any input axis carries one, so the completion
        # check doesn't fail the never-written <group>.csv.
        if self._collect_upstream_missing_paths(self.stream_input, self.groups_input):
            tables["missing"] = self.missing_table_info(self.missing_csv)

        return {
            name: consensus_stream,
            "tables": tables,
            "output_folder": self.output_folder,
        }
