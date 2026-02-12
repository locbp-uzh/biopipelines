# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Distance analysis for calculating distances between atoms and residues.

Analyzes protein structures to calculate distances between specific atoms and residues,
commonly used for ligand-protein binding site analysis and interaction validation.
Outputs CSV with distance metrics for all structures.
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


class Distance(BaseConfig):
    """
    Pipeline tool for analyzing structures to calculate distances between atoms and residues.

    Takes structures as input and outputs CSV with distance metrics for all structures.

    Commonly used for:
    - Ligand-protein binding site analysis
    - Metal coordination analysis
    - Catalytic site geometry verification
    - Interaction distance measurements
    """

    TOOL_NAME = "Distance"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== Distance ==="
echo "Uses biopipelines environment (no additional installation needed)."
echo "=== Distance ready ==="
"""

    # Lazy path descriptors
    analysis_csv = Path(lambda self: os.path.join(self.output_folder, "analysis.csv"))
    config_file = Path(lambda self: os.path.join(self.output_folder, "distance_config.json"))
    structures_ds_json = Path(lambda self: os.path.join(self.output_folder, "structures.json"))
    helper_script = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_distance.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 atom: Union[str, List[str], None] = None,
                 residue: Union[str, List[str], None] = None,
                 method: str = "min",
                 metric_name: str = None,
                 **kwargs):
        """
        Initialize atom-residue distance analysis tool.

        Args:
            structures: Input structures as DataStream or StandardizedOutput
            atom: Atom selection string, list of two selections, or None
            residue: Residue selection string, list of two selections, or None
            method: How to calculate distance ("min", "max", "mean", "closest")
            metric_name: Custom name for the distance column (default: "distance")

        Selection Modes:
            1. Atom-Residue mode: atom=str, residue=str
            2. Atom-Atom mode: atom=[str, str], residue=None
            3. Residue-Residue mode: atom=None, residue=[str, str]

        Selection Syntax:
            Atom selections:
            - 'LIG.Cl' -> ligand chlorine atoms
            - 'HAL.Br' -> halogen bromine atoms
            - 'name CA' -> all alpha carbon atoms

            Residue selections:
            - 'D in IGDWG' -> aspartic acid in sequence context
            - '145' -> residue number 145
            - '145-150' -> residue range 145 to 150
            - '145+147+150' -> specific residues 145, 147, and 150
            - '-1' -> last residue (C-terminus)
            - '-2' -> second-to-last residue
            - '1' -> first residue (N-terminus)
        """
        # Resolve input to DataStream
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        self.atom_selection = atom
        self.residue_selection = residue
        self.distance_metric = method
        self.custom_metric_name = metric_name

        # Validate selection mode
        if atom is None and residue is None:
            raise ValueError("Both atom and residue cannot be None")

        # Validate that if one is a list, the other must be None
        if isinstance(atom, list) and residue is not None:
            raise ValueError("If atom is a list, residue must be None")
        if isinstance(residue, list) and atom is not None:
            raise ValueError("If residue is a list, atom must be None")

        # Validate list lengths
        if isinstance(atom, list) and len(atom) != 2:
            raise ValueError("atom list must contain exactly 2 selections")
        if isinstance(residue, list) and len(residue) != 2:
            raise ValueError("residue list must contain exactly 2 selections")

        # Validate distance metric
        if method not in ["min", "max", "mean", "closest"]:
            raise ValueError(f"Invalid method: {method}. Options: min, max, mean, closest")

        super().__init__(**kwargs)

    def get_metric_name(self) -> str:
        """Get the default metric name."""
        if self.custom_metric_name:
            return self.custom_metric_name
        return "distance"

    def validate_params(self):
        """Validate Distance parameters."""
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")

        # At least one of atom or residue must be specified
        if self.atom_selection is None and self.residue_selection is None:
            raise ValueError("At least one of atom or residue must be specified")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input structures."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        # Determine mode
        if isinstance(self.atom_selection, list):
            mode = "Atom-Atom"
            config_lines.append(f"MODE: {mode}")
            config_lines.append(f"ATOM 1: {self.atom_selection[0]}")
            config_lines.append(f"ATOM 2: {self.atom_selection[1]}")
        elif isinstance(self.residue_selection, list):
            mode = "Residue-Residue"
            config_lines.append(f"MODE: {mode}")
            config_lines.append(f"RESIDUE 1: {self.residue_selection[0]}")
            config_lines.append(f"RESIDUE 2: {self.residue_selection[1]}")
        else:
            mode = "Atom-Residue"
            config_lines.append(f"MODE: {mode}")
            config_lines.append(f"ATOM SELECTION: {self.atom_selection}")
            config_lines.append(f"RESIDUE SELECTION: {self.residue_selection}")

        config_lines.extend([
            f"DISTANCE METRIC: {self.distance_metric}",
            f"OUTPUT METRIC: {self.get_metric_name()}"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate Distance execution script."""
        script_content = "#!/bin/bash\n"
        script_content += "# Distance execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self._generate_script_run_distance_analysis()
        script_content += self.generate_completion_check_footer()

        return script_content

    def _generate_script_run_distance_analysis(self) -> str:
        """Generate the distance analysis execution part of the script."""
        # Write config file at pipeline time
        os.makedirs(self.output_folder, exist_ok=True)

        # Serialize structures DataStream to JSON for HelpScript to load
        with open(self.structures_ds_json, 'w') as f:
            json.dump(self.structures_stream.to_dict(), f, indent=2)

        # Create config data
        config_data = {
            "structures_json": self.structures_ds_json,
            "atom_selection": self.atom_selection,
            "residue_selection": self.residue_selection,
            "distance_metric": self.distance_metric,
            "metric_name": self.get_metric_name(),
            "output_csv": self.analysis_csv
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        return f"""echo "Running distance analysis"
echo "Atom selection: {self.atom_selection}"
echo "Residue selection: {self.residue_selection}"
echo "Distance metric: {self.distance_metric}"
echo "Output: {self.analysis_csv}"

# Run Python analysis script
python "{self.helper_script}" \\
  --config "{self.config_file}"

if [ $? -eq 0 ]; then
    echo "Distance analysis completed successfully"
    echo "Results written to: {self.analysis_csv}"
else
    echo "Error: Distance analysis failed"
    exit 1
fi

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after distance analysis."""
        tables = {
            "distances": TableInfo(
                name="distances",
                path=self.analysis_csv,
                columns=["id", "source_structure", self.get_metric_name()],
                description=f"Distance analysis: {self.atom_selection} to {self.residue_selection}",
                count=len(self.structures_stream)
            )
        }

        return {
            "structures": DataStream.empty("structures", "pdb"),
            "sequences": DataStream.empty("sequences", "fasta"),
            "compounds": DataStream.empty("compounds", "sdf"),
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()

        # Determine mode
        if isinstance(self.atom_selection, list):
            mode = "atom-atom"
        elif isinstance(self.residue_selection, list):
            mode = "residue-residue"
        else:
            mode = "atom-residue"

        base_dict.update({
            "distance_params": {
                "mode": mode,
                "atom_selection": self.atom_selection,
                "residue_selection": self.residue_selection,
                "distance_metric": self.distance_metric,
                "metric_name": self.get_metric_name()
            }
        })
        return base_dict
