# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
SASA (Solvent Accessible Surface Area) tool for calculating delta SASA.

Calculates the change in solvent accessible surface area when the binder protein
is present vs absent. This metric indicates how much of the ligand surface is
buried by protein binding.

delta_SASA = SASA_ligand_alone - SASA_ligand_in_complex
"""

import os
from typing import Dict, List, Any, Union

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


class SASA(BaseConfig):
    """
    SASA tool for calculating delta solvent accessible surface area.

    Computes how much ligand surface area is buried when bound to the protein.
    Uses PyMOL's get_area command for SASA calculations.

    delta_SASA = SASA(ligand alone) - SASA(ligand in complex)

    A higher delta_SASA indicates more ligand surface is buried by the protein,
    suggesting tighter binding.
    """

    TOOL_NAME = "SASA"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== SASA ==="
echo "Requires ProteinEnv (installed with PyMOL.install())"
echo "No additional installation needed."
echo "=== SASA ready ==="
"""

    # Lazy path descriptors
    results_csv = Path(lambda self: os.path.join(self.output_folder, "sasa_analysis.csv"))
    structures_ds_json = Path(lambda self: os.path.join(self.output_folder, "structures.json"))
    helper_script = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_sasa.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 ligand: str,
                 dot_density: int = 4,
                 **kwargs):
        """
        Initialize SASA configuration.

        Args:
            structures: Input structures as DataStream or StandardizedOutput
            ligand: Ligand identifier/code (e.g., ":X:", "LIG", "AMX")
            dot_density: Dot density for SASA calculation (1-4, higher = more accurate)
            **kwargs: Additional parameters

        Output:
            Streams: (none)
            Tables:
                sasa: id | structure | sasa_ligand_alone | sasa_ligand_complex | delta_sasa
        """
        # Resolve input to DataStream
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        self.ligand = ligand
        self.dot_density = dot_density

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate SASA parameters."""
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        if not self.ligand:
            raise ValueError("ligand parameter is required")
        if self.dot_density < 1 or self.dot_density > 4:
            raise ValueError("dot_density must be between 1 and 4")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input structures."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()
        config_lines.extend([
            f"INPUT STRUCTURES: {len(self.structures_stream)} files",
            f"LIGAND: {self.ligand}",
            f"DOT DENSITY: {self.dot_density}"
        ])
        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate SASA execution script."""
        import json

        # Serialize structures DataStream to JSON for HelpScript to load
        with open(self.structures_ds_json, 'w') as f:
            json.dump(self.structures_stream.to_dict(), f, indent=2)

        script_content = "#!/bin/bash\n"
        script_content += "# SASA execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self._generate_script_run_sasa()
        script_content += self.generate_completion_check_footer()

        return script_content

    def _generate_script_run_sasa(self) -> str:
        """Generate the SASA execution part of the script."""
        # Normalize ligand code (remove colons for PyMOL selection)
        ligand_resn = self.ligand.replace(":", "")

        return f"""echo "Running SASA analysis"
echo "Structures: {len(self.structures_stream)} files"
echo "Ligand: {self.ligand}"
echo "Output: {self.results_csv}"

python "{self.helper_script}" \\
    --structures "{self.structures_ds_json}" \\
    --ligand "{ligand_resn}" \\
    --output_csv "{self.results_csv}" \\
    --dot_density {self.dot_density}

if [ $? -eq 0 ]; then
    echo "SASA analysis completed successfully"
    echo "Results written to: {self.results_csv}"
else
    echo "Error: SASA analysis failed"
    exit 1
fi

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after SASA execution."""
        tables = {
            "sasa": TableInfo(
                name="sasa",
                path=self.results_csv,
                columns=["id", "structure", "sasa_ligand_alone", "sasa_ligand_complex", "delta_sasa"],
                description="Solvent accessible surface area analysis",
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
            "sasa_params": {
                "ligand": self.ligand,
                "dot_density": self.dot_density
            }
        })
        return base_dict
