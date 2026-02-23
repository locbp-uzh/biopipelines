# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
PoseBusters validation for computationally generated molecule poses.

Validates docking/prediction outputs by checking bond lengths, angles,
steric clashes, protein interactions, stereochemistry, and more.
Supports ``dock`` mode (ligand + protein) and ``redock`` mode
(+ reference ligand for RMSD comparison).

Reference:
    Buttenschoen et al. (2024) PoseBusters: AI-based docking is not solved.
    Chem Sci 15, 3130. https://github.com/maabuu/posebusters
"""

import os
import json
from typing import Dict, List, Any, Optional, Union

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


class PoseBusters(BaseConfig):
    """
    Pipeline tool for validating computationally generated molecule poses.

    Takes protein-ligand complex structures and runs PoseBusters quality checks
    including bond lengths, bond angles, internal steric clashes, volume overlap
    with protein, and (in redock mode) RMSD to a reference ligand pose.

    Generates CSV with boolean pass/fail columns for each check and an
    ``all_pass`` summary column.

    Commonly used for:
    - Validating docking outputs from Gnina
    - Validating structure predictions from Boltz2
    - Quality-gating predicted poses before downstream analysis
    """

    TOOL_NAME = "PoseBusters"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        """Install PoseBusters in a dedicated conda environment.

        Creates a ``posebusters`` conda environment with Python 3.10 and
        installs the ``posebusters`` package via pip.
        """
        if env_manager == "pip":
            skip = "" if force_reinstall else """# Check if already installed
if python -c "import posebusters" 2>/dev/null; then
    echo "PoseBusters already installed, skipping. Use force_reinstall=True to reinstall."
    exit 0
fi
"""
            return f"""echo "=== Installing PoseBusters (pip) ==="
{skip}pip install posebusters

echo "=== PoseBusters installation complete ==="
"""
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_manager} env list 2>/dev/null | grep -q "posebusters"; then
    echo "PoseBusters already installed, skipping. Use force_reinstall=True to reinstall."
    exit 0
fi
"""
        return f"""echo "=== Installing PoseBusters ==="
{skip}{env_manager} create -n posebusters python=3.10 -y
{env_manager} activate posebusters
pip install posebusters gemmi

echo "=== PoseBusters installation complete ==="
"""

    # Lazy path descriptors
    analysis_csv = Path(lambda self: os.path.join(self.output_folder, "posebusters.csv"))
    config_file = Path(lambda self: os.path.join(self.output_folder, "posebusters_config.json"))
    structures_ds_json = Path(lambda self: os.path.join(self.output_folder, "structures.json"))
    posebusters_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_posebusters.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 ligand: str,
                 reference_ligand: Union[DataStream, StandardizedOutput, None] = None,
                 reference_ligand_code: Optional[str] = None,
                 mode: str = "dock",
                 id_map: Dict[str, str] = {"*": "*_<S>"},
                 **kwargs):
        """
        Initialize PoseBusters validation tool.

        Args:
            structures: Input protein-ligand complexes as DataStream or StandardizedOutput
            ligand: 3-letter residue code for the ligand (e.g., 'LIG', 'ATP')
            reference_ligand: Reference ligand structure for redock mode (SDF or PDB)
            reference_ligand_code: Residue code in reference structure (defaults to ligand)
            mode: 'dock' or 'redock' (auto-set to 'redock' if reference_ligand provided)
            id_map: ID mapping pattern for matching structure IDs to table IDs
            **kwargs: Additional parameters

        Output:
            Streams: (none)
            Tables:
                posebusters (dock): id | source_structure | mol_pred_loaded | mol_cond_loaded | sanitization | all_atoms_connected | bond_lengths | bond_angles | internal_steric_clash | aromatic_ring_flatness | double_bond_flatness | internal_energy | volume_overlap_with_protein | minimum_distance_to_protein | all_pass
                posebusters (redock): adds mol_true_loaded | rmsd before all_pass
        """
        # Resolve structures input to DataStream
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        self.ligand_name = ligand

        # Resolve reference_ligand input
        self.reference_ligand_stream: Optional[DataStream] = None
        if reference_ligand is not None:
            if isinstance(reference_ligand, StandardizedOutput):
                self.reference_ligand_stream = reference_ligand.streams.structures
            elif isinstance(reference_ligand, DataStream):
                self.reference_ligand_stream = reference_ligand
            else:
                raise ValueError(f"reference_ligand must be DataStream, StandardizedOutput, or None, got {type(reference_ligand)}")

        self.reference_ligand_code = reference_ligand_code
        self.id_map = id_map

        # Auto-set mode to redock if reference_ligand is provided
        if self.reference_ligand_stream is not None and mode == "dock":
            self.mode = "redock"
        else:
            self.mode = mode

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate PoseBusters parameters."""
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures cannot be empty")

        if not self.ligand_name or not isinstance(self.ligand_name, str):
            raise ValueError("ligand must be a non-empty string")

        if len(self.ligand_name) > 3:
            raise ValueError(f"ligand code must be 1-3 characters, got '{self.ligand_name}'")

        if self.mode not in ("dock", "redock"):
            raise ValueError(f"mode must be 'dock' or 'redock', got '{self.mode}'")

        if self.mode == "redock" and self.reference_ligand_stream is None:
            raise ValueError("redock mode requires reference_ligand to be provided")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input structures."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"STRUCTURES: {len(self.structures_stream)} files",
            f"LIGAND: {self.ligand_name}",
            f"MODE: {self.mode}",
        ])

        if self.reference_ligand_stream is not None:
            config_lines.append(f"REFERENCE LIGAND: {len(self.reference_ligand_stream)} files")
            if self.reference_ligand_code:
                config_lines.append(f"REFERENCE LIGAND CODE: {self.reference_ligand_code}")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate PoseBusters execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        script_content = "#!/bin/bash\n"
        script_content += "# PoseBusters execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self.generate_script_run_posebusters()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_posebusters(self) -> str:
        """Generate the PoseBusters validation part of the script."""
        # Serialize structures DataStream to JSON
        with open(self.structures_ds_json, 'w') as f:
            json.dump(self.structures_stream.to_dict(), f, indent=2)

        # Build reference_pdb path for redock mode
        reference_pdb = None
        if self.reference_ligand_stream is not None:
            # For redock, serialize reference DataStream and pass path
            ref_ds_json = os.path.join(self.output_folder, "reference_ligand.json")
            with open(ref_ds_json, 'w') as f:
                json.dump(self.reference_ligand_stream.to_dict(), f, indent=2)
            reference_pdb = ref_ds_json

        config_data = {
            "structures_json": self.structures_ds_json,
            "ligand_name": self.ligand_name,
            "mode": self.mode,
            "reference_pdb": reference_pdb,
            "reference_ligand_code": self.reference_ligand_code if self.reference_ligand_code else self.ligand_name,
            "id_map": self.id_map,
            "output_csv": self.analysis_csv
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        return f"""echo "Running PoseBusters validation"
echo "Structures: {len(self.structures_stream)}"
echo "Ligand: {self.ligand_name}"
echo "Mode: {self.mode}"
echo "Output: {self.analysis_csv}"

python "{self.posebusters_py}" --config "{self.config_file}"

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after PoseBusters validation."""
        columns = [
            "id", "source_structure",
            "mol_pred_loaded", "mol_cond_loaded",
            "sanitization", "all_atoms_connected",
            "bond_lengths", "bond_angles", "internal_steric_clash",
            "aromatic_ring_flatness", "double_bond_flatness",
            "internal_energy", "volume_overlap_with_protein",
            "minimum_distance_to_protein",
            "all_pass"
        ]

        # Add redock-specific columns
        if self.mode == "redock":
            columns.insert(-1, "mol_true_loaded")
            columns.insert(-1, "rmsd")

        tables = {
            "posebusters": TableInfo(
                name="posebusters",
                path=self.analysis_csv,
                columns=columns,
                description=f"PoseBusters validation ({self.mode} mode): {self.ligand_name}",
                count=len(self.structures_stream)
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
            "tool_params": {
                "ligand_name": self.ligand_name,
                "mode": self.mode,
                "reference_ligand_code": self.reference_ligand_code,
            }
        })
        return base_dict
