# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""ProLIF tool: protein-ligand interaction fingerprints.

For each complex PDB, computes the ProLIF interaction-fingerprint matrix and
exports it in long form (one row per residue x interaction-type pair, with a
binary `present` flag) so it joins cleanly with the framework's tables.

Reference: https://github.com/chemosim-lab/ProLIF
"""

import os
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream


class ProLIF(BaseConfig):
    """
    ProLIF: compute interaction fingerprints for protein-ligand complexes.

    Inputs:
        structures: PDB complexes with explicit hydrogens. Crystal PDBs lack
                    Hs; chain an upstream protonation step (e.g. OpenBabel
                    with add_hydrogens=True, or Reduce) before this tool.
        ligand: compounds stream carrying the canonical SMILES and residue
                code for the ligand. The SMILES is used as a bond-order
                template for the PDB-extracted ligand (via
                AssignBondOrdersFromTemplate) and the residue `code` locates
                the ligand in the structures PDBs. Both are read from the
                stream's `smiles` / `code` columns at runtime. Supply as
                Ligand("MK1"), Ligand(code="MK1"), or any compounds-producing
                tool's output (the code must resolve to a single value).

    Outputs:
        Streams: (none)
        Tables:
            fingerprints: id | residue | resn | chain | resnum | resi | interaction_type | present
                          (resn duplicates residue, resi duplicates resnum;
                          PyMOL-style aliases for downstream compatibility)
    """

    TOOL_NAME = "ProLIF"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        biopipelines = folders.get("biopipelines", "")
        env_check = cls._env_exists_check("prolif", env_manager)
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_check}; then
    echo "ProLIF already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block("prolif", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("prolif", env_manager, biopipelines)
        return f"""echo "=== Installing ProLIF ==="
{skip}{remove_block}
{env_block}

if {env_manager} run -n prolif python -c "import prolif" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== ProLIF installation complete ==="
else
    echo "ERROR: ProLIF verification failed (cannot import prolif)"
    exit 1
fi
"""

    structures_json = Path(lambda self: self.configuration_path("structures.json"))
    ligand_json = Path(lambda self: self.configuration_path("ligand.json"))
    fingerprints_csv = Path(lambda self: self.table_path("fingerprints"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_prolif.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 ligand: Union[DataStream, StandardizedOutput],
                 **kwargs):
        self.structures = structures
        self.ligand = ligand
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        if isinstance(ligand, StandardizedOutput):
            self.ligand_stream: DataStream = ligand.streams.compounds
        elif isinstance(ligand, DataStream):
            self.ligand_stream = ligand
        else:
            raise ValueError(f"ligand must be DataStream or StandardizedOutput, got {type(ligand)}")

        super().__init__(**kwargs)

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        if not self.ligand_stream or len(self.ligand_stream) == 0:
            raise ValueError("ligand parameter is required and must not be empty")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"STRUCTURES: {len(self.structures_stream)} complexes")
        lines.append(f"LIGAND: {len(self.ligand_stream)} (code resolved from stream at runtime)")
        return lines

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        self.ligand_stream.save_json(self.ligand_json)
        script = "#!/bin/bash\n"
        script += "# ProLIF interaction-fingerprint script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        upstream_missing_flag = self.upstream_missing_flag(self.structures, self.ligand)
        script += f"""echo "Running ProLIF on {len(self.structures_stream)} complex(es)"
python "{self.helper_py}" \\
    --structures-json "{self.structures_json}" \\
    --ligand-json "{self.ligand_json}" \\
    --extras-dir "{self.extras_path()}" \\
    --fingerprints-csv "{self.fingerprints_csv}" \\
    --missing-csv "{self.missing_csv}"{upstream_missing_flag}
"""
        script += self.generate_completion_check_footer()
        return script

    def get_output_files(self) -> Dict[str, Any]:
        tables = {
            "fingerprints": TableInfo(
                name="fingerprints",
                path=self.fingerprints_csv,
                columns=["id", "residue", "resn", "chain", "resnum", "resi", "interaction_type", "present"],
                description="ProLIF interaction fingerprint (long form)",
            ),
            "missing": self.missing_table_info(self.missing_csv),
        }
        return {"tables": tables, "output_folder": self.output_folder}
