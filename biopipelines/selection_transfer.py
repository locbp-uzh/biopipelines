# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
SelectionTransfer: carry a residue selection from the structures it is numbered on onto structures from a later stage of the pipeline.

Residue correspondence comes from a per-chain sequence alignment against the doubled reference sequence, so trims, extensions, shifts, renumbering, chain renames, point mutations and circular permutations all map. Residues with no counterpart are reported, never guessed.

Example::

    st = SelectionTransfer(selection=rfd.tables.structures.designed,
                           reference=rfd,
                           structures=permuted_af)
    ProteinMPNN(structures=permuted_af, redesigned=st.tables.selections.selection)
"""

import os
import re
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


_SPAN = re.compile(r"^[A-Z]?\d+(-\d+)?$")


def _structures_stream(value, param: str) -> DataStream:
    if isinstance(value, StandardizedOutput):
        return value.streams.structures
    if isinstance(value, DataStream):
        return value
    raise ValueError(f"{param} must be DataStream or StandardizedOutput, got {type(value)}")


class SelectionTransfer(BaseConfig):
    """
    SelectionTransfer: map a residue selection from reference structures onto other structures.

    Inputs:
        selection: a table column reference (keyed by reference ids) or a literal selection string
            applied to every reference, e.g. "A45-60+A72". Chainless selections need a single-chain reference.
        reference: structures the selection is numbered on. One structure is broadcast to every
            target; otherwise each target is paired with its reference by id matching.
        structures: structures to map the selection onto; they define the output ids.

    Outputs:
        Tables:
            selections: id | selection | n_residues | source_selection | unmapped | identity | reference.id
                selection is chain-aware in target numbering; unmapped lists reference residues with
                no counterpart; identity is the fraction of identical aligned pairs over the chains used.
            missing: id | removed_by | kind | cause
    """

    TOOL_NAME = "SelectionTransfer"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== SelectionTransfer ==="
echo "Uses biopipelines environment (no additional installation needed)."
touch "$INSTALL_SUCCESS"
echo "=== SelectionTransfer ready ==="
"""

    config_json = Path(lambda self: self.configuration_path("selection_transfer_config.json"))
    structures_json = Path(lambda self: self.configuration_path("structures.json"))
    reference_json = Path(lambda self: self.configuration_path("reference.json"))
    selections_csv = Path(lambda self: self.table_path("selections"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    local_missing_csv = Path(lambda self: self.execution_path("local_missing.csv"))
    transfer_py = Path(lambda self: self.pipe_script_path("pipe_selection_transfer.py"))

    def __init__(self,
                 selection: Union[TableReference, str],
                 reference: Union[DataStream, StandardizedOutput],
                 structures: Union[DataStream, StandardizedOutput],
                 **kwargs):
        self.selection = selection.strip() if isinstance(selection, str) else selection
        self.reference = reference
        self.structures = structures
        self.reference_stream = _structures_stream(reference, "reference")
        self.structures_stream = _structures_stream(structures, "structures")
        super().__init__(**kwargs)

    def validate_params(self):
        if isinstance(self.selection, str):
            bad = [p for p in self.selection.split("+") if not _SPAN.match(p.strip())]
            if not self.selection or bad:
                raise ValueError(f"selection is not a residue selection: {self.selection!r} "
                                 f"(expected e.g. 'A10-20+A45' or '10-20')")
        elif not isinstance(self.selection, TableReference):
            raise TypeError(f"selection must be a table column reference or a selection string, "
                            f"got {type(self.selection).__name__}")
        if len(self.reference_stream) == 0:
            raise ValueError("reference must not be empty")
        if len(self.structures_stream) == 0:
            raise ValueError("structures must not be empty")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        if isinstance(self.selection, str):
            lines.append(f"SELECTION: {self.selection}")
        else:
            lines.append(f"SELECTION: {os.path.basename(self.selection.path)}:{self.selection.column}")
        return lines

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        self.reference_stream.save_json(self.reference_json)
        if isinstance(self.selection, str):
            selection = {"literal": self.selection}
        else:
            selection = {"table": self.selection.path, "column": self.selection.column}
        with open(self.config_json, "w") as f:
            json.dump({"selection": selection}, f, indent=2)

        script = "#!/bin/bash\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += f"""echo "Transferring selection onto target structures"
python "{self.transfer_py}" \\
    --config-json "{self.config_json}" \\
    --structures-json "{self.structures_json}" \\
    --reference-json "{self.reference_json}" \\
    --selections-csv "{self.selections_csv}" \\
    --local-missing-csv "{self.local_missing_csv}"
"""
        script += self.generate_missing_propagation(
            self.structures, self.reference,
            local_missing=self.local_missing_csv, missing_csv=self.missing_csv)
        script += self.generate_completion_check_footer()
        return script

    def get_output_files(self) -> Dict[str, Any]:
        tables = {
            "selections": TableInfo(
                name="selections",
                path=self.selections_csv,
                columns=["id", "selection", "n_residues", "source_selection",
                         "unmapped", "identity", "reference.id"],
                description="Selection mapped onto each target structure",
            ),
            "missing": self.missing_table_info(self.missing_csv),
        }
        return {
            "tables": tables,
            "output_folder": self.output_folder,
        }
