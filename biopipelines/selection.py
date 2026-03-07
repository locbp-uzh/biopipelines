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
from typing import Dict, List, Any, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
    from .biopipelines_io import TableReference
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream
    from biopipelines_io import TableReference


class SelectionOp:
    """Lightweight descriptor for a selection operation."""

    def __init__(self, op_type: str, refs=None, value=None):
        self.op_type = op_type      # "add" | "subtract" | "expand" | "shrink" | "shift" | "invert"
        self.refs = refs or []      # list of TableReference (for add/subtract)
        self.value = value          # int (for expand/shrink/shift) or None

    def to_dict(self) -> Dict[str, Any]:
        d = {"op": self.op_type}
        if self.refs:
            d["refs"] = [{"table": ref.path, "column": ref.column} for ref in self.refs]
        if self.value is not None:
            d["value"] = self.value
        return d

    PDB_AWARE_OPS = {"expand", "shrink", "shift", "invert"}


class Selection(BaseConfig):
    """
    Selection tool with composable, operations-based API.

    Takes a sequence of ``SelectionOp`` objects as positional arguments.
    Operations are applied left-to-right to a running selection (set of
    residues per ID).

    Operations:
        - ``Selection.add(*refs)``      — union referenced columns into running selection
        - ``Selection.subtract(*refs)`` — remove referenced columns from running selection
        - ``Selection.expand(n)``       — PDB-aware: expand intervals by *n* residues each side
        - ``Selection.shrink(n)``       — PDB-aware: shrink intervals by *n* residues each side
        - ``Selection.shift(n)``        — PDB-aware: shift intervals by *n* residues
        - ``Selection.invert()``        — PDB-aware: replace with complement
    """

    TOOL_NAME = "SelectionEditor"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== SelectionEditor ==="
echo "Requires ProteinEnv (installed with PyMOL.install())"
echo "No additional installation needed."
echo "=== SelectionEditor ready ==="
"""

    # Lazy path descriptors
    selections_csv = Path(lambda self: os.path.join(self.output_folder, "selections.csv"))
    config_json = Path(lambda self: os.path.join(self.output_folder, "config.json"))
    structures_ds_json = Path(lambda self: os.path.join(self.output_folder, "structures.json"))
    selection_editor_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_selection.py"))

    # ── Static factory methods (return SelectionOp, NOT StandardizedOutput) ──

    @staticmethod
    def add(*refs) -> SelectionOp:
        return SelectionOp("add", refs=list(refs))

    @staticmethod
    def subtract(*refs) -> SelectionOp:
        return SelectionOp("subtract", refs=list(refs))

    @staticmethod
    def expand(n: int) -> SelectionOp:
        return SelectionOp("expand", value=n)

    @staticmethod
    def shrink(n: int) -> SelectionOp:
        return SelectionOp("shrink", value=n)

    @staticmethod
    def shift(n: int) -> SelectionOp:
        return SelectionOp("shift", value=n)

    @staticmethod
    def invert() -> SelectionOp:
        return SelectionOp("invert")

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
        # Build ordered operations list
        self.operations: List[SelectionOp] = []

        for op in operations:
            if not isinstance(op, SelectionOp):
                raise TypeError(
                    f"Positional arguments must be SelectionOp objects (use Selection.add(), "
                    f"Selection.expand(), etc.), got {type(op).__name__}"
                )
            self.operations.append(op)

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
        # First operation must be add to initialise the running selection
        if not self.operations or self.operations[0].op_type != "add":
            raise ValueError(
                "The first operation must be add() to provide an initial selection value."
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
                cols = ", ".join(ref.column for ref in op.refs)
                op_summaries.append(f"{op.op_type}({cols})")
            elif op.value is not None:
                op_summaries.append(f"{op.op_type}({op.value})")
            else:
                op_summaries.append(f"{op.op_type}()")

        config_lines.append(f"OPERATIONS: {' -> '.join(op_summaries)}")
        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        # Serialize structures DataStream if available
        structures_json_path = None
        if self.structures_stream:
            self.structures_stream.save_json(self.structures_ds_json)
            structures_json_path = self.structures_ds_json

        # Build config for HelpScript
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
    exit 1
fi

echo "Modified selections saved to: {self.selections_csv}"

"""
        script_content += self.generate_completion_check_footer()
        return script_content

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files."""
        count = len(self.structures_stream) if self.structures_stream else "variable"

        tables = {
            "selections": TableInfo(
                name="selections",
                path=self.selections_csv,
                columns=["id", "selection"],
                description="Modified selections",
                count=count
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
