# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Selection tool with composable, operations-based API.

Combines and modifies PyMOL-formatted selection strings using a sequence of
operations (add, subtract, expand, shrink, shift, invert) applied left-to-right
to a running selection set.

Example::

    # Union two columns
    sel = Selection(Selection.add(fuse.tables.sequences.L1,
                                  fuse.tables.sequences.L2))

    # Sequential: union, subtract, then expand
    sel = Selection(
        Selection.add(fuse.tables.sequences.L1,
                      fuse.tables.sequences.L2),
        Selection.subtract(other.tables.structures.col8),
        Selection.expand(1),
        structures=rfd
    )

    # Single column with expand + invert
    sel = Selection(
        Selection.add(rfd.tables.structures.designed),
        Selection.expand(5),
        Selection.invert(),
        structures=rfd
    )
"""

import os
import json
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
    from .biopipelines_io import TableReference
except ImportError:
    import sys
    sys.path.insert(0, os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream
    from biopipelines_io import TableReference


class SelectionOp:
    """Lightweight descriptor for a selection operation."""

    def __init__(self, op_type: str, refs=None, value=None,
                 stream: Optional[DataStream] = None,
                 filter_expr: Optional[str] = None,
                 direction: Optional[str] = None):
        self.op_type = op_type      # "add" | "subtract" | "expand" | "shrink" | "shift" | "invert"
        self.refs = refs or []      # list of TableReference (for table-based add/subtract)
        self.value = value          # int (for expand/shrink/shift) or None
        self.stream = stream        # DataStream (for stream-based add/subtract)
        self.filter_expr = filter_expr    # e.g. "propensity>0.5", "rmsf<=0.3"
        self.direction = direction  # "n" | "c" | "nc" | None (for expand/shrink)

    def to_dict(self) -> Dict[str, Any]:
        d = {"op": self.op_type}
        if self.refs:
            d["refs"] = [{"table": ref.path, "column": ref.column} for ref in self.refs]
        if self.value is not None:
            d["value"] = self.value
        if self.stream is not None:
            json_path = getattr(self.stream, "_json_path", None)
            if json_path is None:
                raise RuntimeError(
                    "SelectionOp.to_dict() called before generate_script() — "
                    "stream JSON path has not been set yet."
                )
            d["stream_json"] = json_path
        if self.filter_expr is not None:
            d["filter_expr"] = self.filter_expr
        if self.direction is not None:
            d["direction"] = self.direction
        return d

    PDB_AWARE_OPS = {"expand", "shrink", "shift", "invert",
                     "n_terminus", "c_terminus", "termini", "all_residues", "gaps"}


class Selection(BaseConfig):
    """
    Selection tool with composable, operations-based API.

    Takes a sequence of ``SelectionOp`` objects as positional arguments.
    Operations are applied left-to-right to a running selection (set of
    residues per ID).

    Operations:
        - ``Selection.add(*refs)``                        — union table columns into running selection
        - ``Selection.add(ds, include="propensity>0.5")``  — union residues from per-residue stream
        - ``Selection.subtract(*refs)``                   — remove table columns from running selection
        - ``Selection.subtract(ds, include="propensity>0.5")`` — remove residues from per-residue stream
        - ``Selection.expand(n, direction="nc")``  — PDB-aware: expand by *n* residues (n/c/nc)
        - ``Selection.shrink(n, direction="nc")``  — PDB-aware: shrink by *n* residues (n/c/nc)
        - ``Selection.dilate`` / ``Selection.erode``       — aliases for expand / shrink
        - ``Selection.close(n, direction="nc")``   — expand then shrink (fills gaps <= *n*)
        - ``Selection.open(n, direction="nc")``    — shrink then expand (removes protrusions <= *n*)
        - ``Selection.shift(n)``        — PDB-aware: shift intervals by *n* residues
        - ``Selection.invert()``        — PDB-aware: replace with complement
        - ``Selection.n_terminus(n=1)`` — PDB-aware: select first *n* residues of each chain
        - ``Selection.c_terminus(n=1)`` — PDB-aware: select last *n* residues of each chain
        - ``Selection.termini(n=1)``    — PDB-aware: select first and last *n* residues of each chain
        - ``Selection.all_residues()``  — PDB-aware: select all residues in the structure
        - ``Selection.gaps()``          — PDB-aware: select residues flanking structural gaps

    Direction parameter (for expand/shrink/close/open):
        - ``"n"``:  toward N-terminus (lower residue numbers)
        - ``"c"``:  toward C-terminus (higher residue numbers)
        - ``"nc"``: both directions (default)

    Stream-based filtering:
        When a DataStream is passed, it must have format ``per-residue-values-csv``
        with a ``resi`` column. The ``include=`` or ``exclude=`` expression
        specifies both the column and threshold: ``"column op value"``, e.g.
        ``include="propensity>0.5"`` keeps rows where propensity > 0.5.
        Supported operators: ``>``, ``<``, ``>=``, ``<=``, ``==``, ``!=``.
    """

    TOOL_NAME = "Selection"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== Selection ==="
echo "Requires ProteinEnv (installed with PyMOL.install())"
echo "No additional installation needed."
touch "$INSTALL_SUCCESS"
echo "=== Selection ready ==="
"""

    # Lazy path descriptors
    selections_csv = Path(lambda self: self.table_path("selections"))
    config_json = Path(lambda self: self.configuration_path("config.json"))
    structures_ds_json = Path(lambda self: self.configuration_path("structures.json"))
    selection_editor_py = Path(lambda self: self.pipe_script_path("pipe_selection.py"))

    # ── Static factory methods (return SelectionOp, NOT StandardizedOutput) ──

    @staticmethod
    def add(*args, include: Optional[str] = None, exclude: Optional[str] = None) -> SelectionOp:
        """Union table refs or per-residue stream rows into the running selection.

        Table-based (pass TableReference objects)::

            Selection.add(fuse.tables.sequences.L1, fuse.tables.sequences.L2)

        Stream-based (pass a DataStream of format per-residue-values-csv)::

            Selection.add(cpred.streams.propensities, include="propensity>0.5")

        Args:
            *args:    One or more TableReference objects, or a single DataStream.
            include:  Filter expression ``"column op value"`` for rows to include,
                      e.g. ``"propensity>0.5"``, ``"rmsf<=1.2"``.
            exclude:  Filter expression ``"column op value"`` for rows to exclude.
        """
        if len(args) == 1 and isinstance(args[0], DataStream):
            if include is not None and exclude is not None:
                raise ValueError("Selection.add: use either include= or exclude=, not both")
            filter_expr = include if include is not None else exclude
            return SelectionOp("add", stream=args[0], filter_expr=filter_expr)
        if include is not None or exclude is not None:
            raise ValueError("Selection.add: include=/exclude= are only valid with a DataStream argument")
        return SelectionOp("add", refs=list(args))

    @staticmethod
    def subtract(*args, include: Optional[str] = None, exclude: Optional[str] = None) -> SelectionOp:
        """Remove table refs or per-residue stream rows from the running selection.

        Table-based (pass TableReference objects)::

            Selection.subtract(other.tables.structures.col8)

        Stream-based (pass a DataStream of format per-residue-values-csv)::

            Selection.subtract(cpred.streams.propensities, include="propensity>0.8")

        Args:
            *args:    One or more TableReference objects, or a single DataStream.
            include:  Filter expression ``"column op value"`` for rows to subtract,
                      e.g. ``"propensity>0.8"``.
            exclude:  Filter expression ``"column op value"`` for rows to subtract (inverted).
        """
        if len(args) == 1 and isinstance(args[0], DataStream):
            if include is not None and exclude is not None:
                raise ValueError("Selection.subtract: use either include= or exclude=, not both")
            filter_expr = include if include is not None else exclude
            return SelectionOp("subtract", stream=args[0], filter_expr=filter_expr)
        if include is not None or exclude is not None:
            raise ValueError("Selection.subtract: include=/exclude= are only valid with a DataStream argument")
        return SelectionOp("subtract", refs=list(args))

    @staticmethod
    def expand(n: int, direction: str = "nc") -> SelectionOp:
        """Expand selection by *n* residues. Direction: ``"n"`` (toward N-term),
        ``"c"`` (toward C-term), or ``"nc"`` (both, default)."""
        if direction not in ("n", "c", "nc"):
            raise ValueError(f"direction must be 'n', 'c', or 'nc', got '{direction}'")
        return SelectionOp("expand", value=n, direction=direction)

    @staticmethod
    def shrink(n: int, direction: str = "nc") -> SelectionOp:
        """Shrink selection by *n* residues. Direction: ``"n"`` (from N-term side),
        ``"c"`` (from C-term side), or ``"nc"`` (both, default)."""
        if direction not in ("n", "c", "nc"):
            raise ValueError(f"direction must be 'n', 'c', or 'nc', got '{direction}'")
        return SelectionOp("shrink", value=n, direction=direction)

    # Morphological aliases
    dilate = expand
    erode = shrink

    @staticmethod
    def close(n: int, direction: str = "nc") -> List[SelectionOp]:
        """Morphological close: expand then shrink. Fills gaps of size <= *n*."""
        if direction not in ("n", "c", "nc"):
            raise ValueError(f"direction must be 'n', 'c', or 'nc', got '{direction}'")
        return [SelectionOp("expand", value=n, direction=direction),
                SelectionOp("shrink", value=n, direction=direction)]

    @staticmethod
    def open(n: int, direction: str = "nc") -> List[SelectionOp]:
        """Morphological open: shrink then expand. Removes protrusions of size <= *n*."""
        if direction not in ("n", "c", "nc"):
            raise ValueError(f"direction must be 'n', 'c', or 'nc', got '{direction}'")
        return [SelectionOp("shrink", value=n, direction=direction),
                SelectionOp("expand", value=n, direction=direction)]

    @staticmethod
    def shift(n: int) -> SelectionOp:
        return SelectionOp("shift", value=n)

    @staticmethod
    def invert() -> SelectionOp:
        return SelectionOp("invert")

    @staticmethod
    def n_terminus(n: int = 1) -> SelectionOp:
        """Select the first *n* residues of each chain (N-terminal region).
        Requires ``structures=``."""
        return SelectionOp("n_terminus", value=n)

    @staticmethod
    def c_terminus(n: int = 1) -> SelectionOp:
        """Select the last *n* residues of each chain (C-terminal region).
        Requires ``structures=``."""
        return SelectionOp("c_terminus", value=n)

    @staticmethod
    def termini(n: int = 1) -> SelectionOp:
        """Select the first and last *n* residues of each chain.
        Requires ``structures=``."""
        return SelectionOp("termini", value=n)

    @staticmethod
    def all_residues() -> SelectionOp:
        """Select all residues present in the structure.
        Requires ``structures=``."""
        return SelectionOp("all_residues")

    @staticmethod
    def gaps() -> SelectionOp:
        """Select the residues immediately flanking each structural gap
        (one residue on each side of each unresolved region).
        Requires ``structures=``."""
        return SelectionOp("gaps")

    # ── Constructor ──

    def __init__(self,
                 *operations: SelectionOp,
                 structures=None,
                 **kwargs):
        """
        Initialize Selection tool.

        Args:
            *operations: Sequence of SelectionOp objects describing the pipeline.
            structures: DataStream or StandardizedOutput providing PDB files.
                Required when any PDB-aware operation (expand/shrink/shift/invert)
                is used.
            **kwargs: Additional parameters forwarded to BaseConfig.

        Output:
            Tables:
                selections: id | selection
        """
        # Build ordered operations list (flatten lists from close/open composites)
        self.operations: List[SelectionOp] = []

        for op in operations:
            if isinstance(op, list):
                for sub_op in op:
                    if not isinstance(sub_op, SelectionOp):
                        raise TypeError(
                            f"Composite operations must contain SelectionOp objects, "
                            f"got {type(sub_op).__name__}"
                        )
                    self.operations.append(sub_op)
            elif isinstance(op, SelectionOp):
                self.operations.append(op)
            else:
                # Implicit add(): bare TableReferences, DataStreams, or selection
                # strings are wrapped in Selection.add(x). Lets
                # `Selection(rfd.tables.structures.designed)` stand in for
                # `Selection(Selection.add(rfd.tables.structures.designed))`.
                self.operations.append(Selection.add(op))

        # Resolve structures (needed for PDB-aware ops)
        self.structures_stream = None
        if isinstance(structures, StandardizedOutput):
            self.structures_stream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        elif structures is not None:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate Selection parameters."""
        # First operation must initialise the running selection
        VALID_FIRST_OPS = {"add", "n_terminus", "c_terminus", "termini", "all_residues", "gaps"}
        if not self.operations or self.operations[0].op_type not in VALID_FIRST_OPS:
            raise ValueError(
                "The first operation must be add() or a predefined pattern "
                "(n_terminus, c_terminus, termini, all_residues, gaps)."
            )

        # Each add/subtract must have either refs or a stream (not neither)
        for op in self.operations:
            if op.op_type in ("add", "subtract"):
                if not op.refs and op.stream is None:
                    raise ValueError(
                        f"Selection.{op.op_type}() requires either TableReference args or a DataStream"
                    )

        # PDB-aware ops require structures
        has_pdb_op = any(op.op_type in SelectionOp.PDB_AWARE_OPS for op in self.operations)
        if has_pdb_op and self.structures_stream is None:
            raise ValueError(
                "structures= is required when using PDB-aware operations "
                "(expand, shrink, shift, invert)."
            )

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input sources."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        if self.structures_stream:
            config_lines.append(f"INPUT STRUCTURES: {len(self.structures_stream)} files")

        # Summarise operations
        op_summaries = []
        for op in self.operations:
            if op.op_type in ("add", "subtract"):
                if op.stream is not None:
                    op_summaries.append(f"{op.op_type}(stream, {op.filter_expr})")
                else:
                    cols = ", ".join(ref.column for ref in op.refs)
                    op_summaries.append(f"{op.op_type}({cols})")
            elif op.value is not None:
                dir_str = f", {op.direction}" if op.direction and op.direction != "nc" else ""
                op_summaries.append(f"{op.op_type}({op.value}{dir_str})")
            else:
                op_summaries.append(f"{op.op_type}()")

        config_lines.append(f"OPERATIONS: {' -> '.join(op_summaries)}")
        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate execution script."""
        # Serialize structures DataStream if available
        structures_json_path = None
        if self.structures_stream:
            self.structures_stream.save_json(self.structures_ds_json)
            structures_json_path = self.structures_ds_json

        # Serialize any per-residue streams used in operations
        for i, op in enumerate(self.operations):
            if op.stream is not None:
                stream_json = os.path.join(self.output_folder, f".stream_op{i}.json")
                op.stream.save_json(stream_json)
                op.stream._json_path = stream_json

        # Build config for pipe_script
        config_data = {
            "operations": [op.to_dict() for op in self.operations],
            "output_csv": self.selections_csv
        }
        if structures_json_path:
            config_data["structures_json"] = structures_json_path

        with open(self.config_json, 'w') as f:
            json.dump(config_data, f, indent=2)

        # Build display string for echo
        op_display = " -> ".join(
            op.op_type + (f"({op.value})" if op.value is not None else
                          f"({', '.join(r.column for r in op.refs)})" if op.refs else "()")
            for op in self.operations
        )
        num_structures = len(self.structures_stream) if self.structures_stream else "N/A"

        script_content = "#!/bin/bash\n"
        script_content += "# SelectionEditor execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()

        script_content += f"""echo "Modifying selections"
echo "Structures: {num_structures}"
echo "Operations: {op_display}"

# Run selection pipeline
python {self.selection_editor_py} \\
    "{self.config_json}"

if [ $? -eq 0 ]; then
    echo "Selection modification completed successfully"
else
    echo "Selection modification failed"
fi

echo "Modified selections saved to: {self.selections_csv}"

"""
        script_content += self.generate_completion_check_footer()
        return script_content

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files."""
        tables = {
            "selections": TableInfo(
                name="selections",
                path=self.selections_csv,
                columns=["id", "selection"],
                description="Modified selections"
            )
        }

        return {
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "selection_params": {
                "operations": [op.to_dict() for op in self.operations],
                "has_structures": self.structures_stream is not None
            }
        })
        return base_dict
