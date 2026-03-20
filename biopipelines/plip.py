# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
PLIP (Protein-Ligand Interaction Profiler) analysis tool.

Detects non-covalent interactions between proteins and ligands: hydrogen bonds,
hydrophobic contacts, salt bridges, pi-stacking, pi-cation, halogen bonds,
water bridges, and metal complexes.
Outputs detailed interaction CSV and per-ligand summary CSV.
"""

import os
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


class PLIP(BaseConfig):
    """
    Pipeline tool for profiling protein-ligand non-covalent interactions.

    Takes protein structures in PDB format and automatically detects all
    binding sites and their non-covalent interactions using PLIP.

    Generates two CSVs: a detailed interaction table (one row per interaction)
    and a per-ligand summary table with interaction counts by type.

    Commonly used for:
    - Characterizing binding site interactions
    - Comparing interaction profiles across designs
    - Identifying key stabilizing contacts
    - Binding mode analysis
    """

    # Tool identification
    TOOL_NAME = "PLIP"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        if env_manager == "pip":
            skip = "" if force_reinstall else """# Check if already installed
if python -c "from plip.structure.preparation import PDBComplex" 2>/dev/null; then
    echo "PLIP already installed, skipping. Use force_reinstall=True to reinstall."
    exit 0
fi
"""
            return f"""echo "=== Installing PLIP (pip) ==="
{skip}pip install plip
echo "=== PLIP installation complete ==="
"""
        env_file = os.path.join(folders.get("biopipelines", "."), "Environments", "PLIP.yaml")
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_manager} env list 2>/dev/null | grep -q "PLIPEnv"; then
    echo "PLIP already installed, skipping. Use force_reinstall=True to reinstall."
    exit 0
fi
"""
        return f"""echo "=== Installing PLIP ==="
{skip}{env_manager} env create -f "{env_file}" -y
echo "=== PLIP installation complete ==="
"""

    # Lazy path descriptors
    interactions_csv = Path(lambda self: os.path.join(self.output_folder, "interactions.csv"))
    summary_csv = Path(lambda self: os.path.join(self.output_folder, "summary.csv"))
    config_file = Path(lambda self: os.path.join(self.output_folder, "plip_config.json"))
    structures_ds_json = Path(lambda self: os.path.join(self.output_folder, "structures.json"))
    plip_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_plip.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 **kwargs):
        """
        Initialize PLIP interaction profiling tool.

        Args:
            structures: Input structures as DataStream or StandardizedOutput.
                        Must be PDB format. CIF-only input raises an error;
                        mixed pdb|cif input warns that CIF files will be skipped.
            **kwargs: Additional parameters

        Output:
            Streams: (none)
            Tables:
                interactions: id | ligand_name | ligand_chain | interaction_type | protein_residue | protein_chain | protein_residue_number | distance
                summary: id | ligand_name | ligand_chain | hydrophobic | hbond | water_bridge | salt_bridge | pi_stacking | pi_cation | halogen_bond | metal_complex | total_interactions
        """
        # Resolve input to DataStream
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate PLIP parameters."""
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures cannot be empty")

        fmt = getattr(self.structures_stream, 'format', None)
        if fmt == "cif":
            raise ValueError(
                "PLIP does not support CIF format. Input structures are CIF-only. "
                "Use convert='pdb' in PDB() or RCSB() to ensure PDB format output."
            )
        if fmt == "pdb|cif":
            print(
                "  WARNING: PLIP does not support CIF format. Input contains mixed pdb|cif structures. "
                "CIF files will be skipped at runtime — only PDB files will be analyzed."
            )

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input structures."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"STRUCTURES: {len(self.structures_stream)} files",
            f"FORMAT: {getattr(self.structures_stream, 'format', 'unknown')}",
            f"MODE: auto-detect all binding sites",
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate PLIP execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        script_content = "#!/bin/bash\n"
        script_content += "# PLIP execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self.generate_script_run_plip()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_plip(self) -> str:
        """Generate the PLIP interaction profiling part of the script."""
        import json

        # Serialize structures DataStream to JSON for HelpScript to load
        self.structures_stream.save_json(self.structures_ds_json)

        config_data = {
            "structures_json": self.structures_ds_json,
            "interactions_csv": self.interactions_csv,
            "summary_csv": self.summary_csv,
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        return f"""echo "Running PLIP interaction profiling"
echo "Structures: {len(self.structures_stream)}"
echo "Output interactions: {self.interactions_csv}"
echo "Output summary: {self.summary_csv}"

python "{self.plip_py}" --config "{self.config_file}"

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after PLIP analysis."""
        tables = {
            "interactions": TableInfo(
                name="interactions",
                path=self.interactions_csv,
                columns=["id", "ligand_name", "ligand_chain", "interaction_type",
                         "protein_residue", "protein_chain", "protein_residue_number",
                         "distance"],
                description="Detailed non-covalent interactions detected by PLIP (one row per interaction)",
                count=None  # variable number of interactions
            ),
            "summary": TableInfo(
                name="summary",
                path=self.summary_csv,
                columns=["id", "ligand_name", "ligand_chain",
                         "hydrophobic", "hbond", "water_bridge", "salt_bridge",
                         "pi_stacking", "pi_cation", "halogen_bond", "metal_complex",
                         "total_interactions"],
                description="Per-ligand interaction summary counts",
                count=None  # variable number of ligands
            ),
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
                "format": getattr(self.structures_stream, 'format', 'unknown'),
            }
        })
        return base_dict
