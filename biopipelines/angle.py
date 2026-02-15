# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Angle analysis for calculating angles and torsional angles between atoms.

Analyzes protein structures to calculate angles between 3 atoms (bond angle on middle atom)
or torsional angles between 4 atoms (dihedral angle). Outputs CSV with angle metrics
for all structures.
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


class Angle(BaseConfig):
    """
    Pipeline tool for analyzing structures to calculate angles between atoms.

    Takes structures and atom selections as input, outputs CSV with angle metrics.

    Supports two modes:
    - 3 atoms (A-B-C): Calculate bond angle at B
    - 4 atoms (A-B-C-D): Calculate torsional/dihedral angle

    Output unit is configurable: "degrees" (default) or "radians".

    Commonly used for:
    - Backbone phi/psi angle analysis
    - Side chain rotamer analysis
    - Metal coordination geometry
    - Ligand binding geometry verification
    """

    TOOL_NAME = "Angle"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== Angle ==="
echo "Uses biopipelines environment (no additional installation needed)."
echo "=== Angle ready ==="
"""

    # Lazy path descriptors
    analysis_csv = Path(lambda self: os.path.join(self.output_folder, "analysis.csv"))
    config_file = Path(lambda self: os.path.join(self.output_folder, "angle_config.json"))
    structures_ds_json = Path(lambda self: os.path.join(self.output_folder, "structures.json"))
    helper_script = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_angle.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 atoms: List[str],
                 metric_name: str = None,
                 unit: str = "degrees",
                 **kwargs):
        """
        Initialize angle analysis tool.

        Args:
            structures: Input structures as DataStream or StandardizedOutput
            atoms: List of 3 or 4 atom selection strings
            metric_name: Custom name for the angle column (default: "angle" or "torsion")
            unit: Output unit, "degrees" (default) or "radians"

        Selection Syntax (same as Distance):
            Atom selections:
            - 'LIG.Cl' -> ligand chlorine atoms
            - 'LIG.C1' -> ligand atom named C1
            - 'name CA' -> all alpha carbon atoms

            Residue selections:
            - 'D in IGDWG' -> aspartic acid in sequence context
            - '145' -> residue number 145
            - '145.CA' -> alpha carbon of residue 145
            - '-1' -> last residue (C-terminus)
            - '-1.C' -> carbonyl carbon of last residue

        Examples:
            # Bond angle at CA of residue 10 (N-CA-C angle)
            Angle(structures=boltz, atoms=['10.N', '10.CA', '10.C'])

            # Phi angle (C-N-CA-C)
            Angle(structures=boltz, atoms=['9.C', '10.N', '10.CA', '10.C'])

            # Psi angle (N-CA-C-N)
            Angle(structures=boltz, atoms=['10.N', '10.CA', '10.C', '11.N'])

            # Chi1 angle for residue 50
            Angle(structures=boltz, atoms=['50.N', '50.CA', '50.CB', '50.CG'])

            # Ligand geometry
            Angle(structures=boltz, atoms=['LIG.C1', 'LIG.C2', 'LIG.C3'])
        """
        # Resolve input to DataStream
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        self.atom_selections = atoms
        self.custom_metric_name = metric_name
        self.unit = unit

        # Validate atom count
        if not isinstance(atoms, list):
            raise ValueError("atoms must be a list of selection strings")
        if len(atoms) not in [3, 4]:
            raise ValueError(f"atoms must contain exactly 3 or 4 selections, got {len(atoms)}")

        # Validate unit
        if unit not in ("degrees", "radians"):
            raise ValueError(f"unit must be 'degrees' or 'radians', got '{unit}'")

        # Determine angle type
        self.is_torsion = len(atoms) == 4

        super().__init__(**kwargs)

    def get_metric_name(self) -> str:
        """Get the default metric name."""
        if self.custom_metric_name:
            return self.custom_metric_name
        return "torsion" if self.is_torsion else "angle"

    def validate_params(self):
        """Validate Angle parameters."""
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")

        if not self.atom_selections:
            raise ValueError("atoms parameter is required")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input structures."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        angle_type = "Torsional (4 atoms)" if self.is_torsion else "Bond (3 atoms)"
        config_lines.append(f"ANGLE TYPE: {angle_type}")

        for i, selection in enumerate(self.atom_selections):
            config_lines.append(f"ATOM {i+1}: {selection}")

        config_lines.append(f"OUTPUT METRIC: {self.get_metric_name()}")
        config_lines.append(f"UNIT: {self.unit}")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate Angle execution script."""
        script_content = "#!/bin/bash\n"
        script_content += "# Angle analysis execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self._generate_script_run_angle_analysis()
        script_content += self.generate_completion_check_footer()

        return script_content

    def _generate_script_run_angle_analysis(self) -> str:
        """Generate the angle analysis execution part of the script."""
        # Write config file at pipeline time
        os.makedirs(self.output_folder, exist_ok=True)

        # Serialize structures DataStream to JSON for HelpScript to load
        with open(self.structures_ds_json, 'w') as f:
            json.dump(self.structures_stream.to_dict(), f, indent=2)

        # Create config data
        config_data = {
            "structures_json": self.structures_ds_json,
            "atom_selections": self.atom_selections,
            "is_torsion": self.is_torsion,
            "metric_name": self.get_metric_name(),
            "unit": self.unit,
            "output_csv": self.analysis_csv
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        angle_type = "torsional" if self.is_torsion else "bond"
        atoms_str = " -> ".join(self.atom_selections)

        return f"""echo "Running {angle_type} angle analysis"
echo "Atoms: {atoms_str}"
echo "Output: {self.analysis_csv}"

# Run Python analysis script
python "{self.helper_script}" \\
  --config "{self.config_file}"

if [ $? -eq 0 ]; then
    echo "Angle analysis completed successfully"
    echo "Results written to: {self.analysis_csv}"
else
    echo "Error: Angle analysis failed"
    exit 1
fi

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after angle analysis."""
        angle_type = "Torsional angle" if self.is_torsion else "Bond angle"
        atoms_desc = " -> ".join(self.atom_selections)

        tables = {
            "angles": TableInfo(
                name="angles",
                path=self.analysis_csv,
                columns=["id", "source_structure", self.get_metric_name(), "unit"],
                description=f"{angle_type} ({self.unit}): {atoms_desc}",
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

        base_dict.update({
            "angle_params": {
                "is_torsion": self.is_torsion,
                "atom_selections": self.atom_selections,
                "metric_name": self.get_metric_name(),
                "unit": self.unit
            }
        })
        return base_dict
