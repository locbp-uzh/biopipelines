# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Angle analysis for calculating angles and torsional angles between atoms.

Analyzes protein structures to calculate angles between atoms. Three modes are supported,
selected by the shape of the `atoms` argument:

- (a, o, b)           : bond angle at o (a-o-b)
- (a, x1, x2, b)      : dihedral angle along axis x1-x2
- ((a1, a2), (b1, b2)): angle between vectors a1->a2 and b1->b2

Outputs CSV with angle metrics for all structures.
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

    # Angle modes
    MODE_BOND = "bond"
    MODE_TORSION = "torsion"
    MODE_VECTOR = "vector"

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
                 atoms: tuple,
                 metric_name: str = None,
                 unit: str = "degrees",
                 **kwargs):
        """
        Initialize angle analysis tool.

        Args:
            structures: Input structures as DataStream or StandardizedOutput
            atoms: Tuple specifying which angle to compute. Three forms are accepted:

                (a, o, b)
                    Bond angle at o: the angle formed by selections a, o, b with
                    vertex at o.

                (a, x1, x2, b)
                    Dihedral angle along axis x1-x2 (standard IUPAC torsion).

                ((a1, a2), (b1, b2))
                    Angle between vectors a1->a2 and b1->b2 (each selection may
                    resolve to multiple atoms; centroids are used).

            metric_name: Custom name for the angle column (default: "angle",
                         "torsion", or "vector_angle" depending on mode)
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

        Output:
            Streams: (none)
            Tables:
                angles: id | source_structure | <metric_name> | unit
        """
        # Resolve input to DataStream
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        self.custom_metric_name = metric_name
        self.unit = unit

        # Validate unit
        if unit not in ("degrees", "radians"):
            raise ValueError(f"unit must be 'degrees' or 'radians', got '{unit}'")

        # Detect mode from atoms shape
        self.mode, self.atom_selections = self._parse_atoms(atoms)

        super().__init__(**kwargs)

    @staticmethod
    def _parse_atoms(atoms):
        """
        Parse the atoms argument and return (mode, flat_selections).

        flat_selections is always a plain list of strings (or nested list for
        vector mode) that can be serialized to JSON and understood by the
        HelpScript.
        """
        if not isinstance(atoms, tuple):
            raise ValueError(
                "atoms must be a tuple: (a, o, b), (a, x1, x2, b), "
                "or ((a1, a2), (b1, b2))"
            )

        # Vector mode: tuple of two 2-tuples
        if (len(atoms) == 2
                and isinstance(atoms[0], tuple) and len(atoms[0]) == 2
                and isinstance(atoms[1], tuple) and len(atoms[1]) == 2):
            (a1, a2), (b1, b2) = atoms
            return Angle.MODE_VECTOR, [[a1, a2], [b1, b2]]

        # Bond angle: 3-tuple of strings
        if len(atoms) == 3:
            a, o, b = atoms
            return Angle.MODE_BOND, [a, o, b]

        # Dihedral: 4-tuple of strings
        if len(atoms) == 4:
            a, x1, x2, b = atoms
            return Angle.MODE_TORSION, [a, x1, x2, b]

        raise ValueError(
            "atoms must be a 3-tuple (a, o, b), a 4-tuple (a, x1, x2, b), "
            "or a pair of 2-tuples ((a1, a2), (b1, b2))"
        )

    def get_metric_name(self) -> str:
        """Get the default metric name."""
        if self.custom_metric_name:
            return self.custom_metric_name
        return {"bond": "angle", "torsion": "torsion", "vector": "vector_angle"}[self.mode]

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

        angle_type = {
            self.MODE_BOND: "Bond (3 atoms)",
            self.MODE_TORSION: "Torsional (4 atoms)",
            self.MODE_VECTOR: "Vector-vector",
        }[self.mode]
        config_lines.append(f"ANGLE TYPE: {angle_type}")

        if self.mode == self.MODE_VECTOR:
            (a1, a2), (b1, b2) = self.atom_selections
            config_lines.append(f"VECTOR 1: {a1} -> {a2}")
            config_lines.append(f"VECTOR 2: {b1} -> {b2}")
        else:
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
        # Write config file at configuration time
        os.makedirs(self.output_folder, exist_ok=True)

        # Serialize structures DataStream to JSON for HelpScript to load
        self.structures_stream.save_json(self.structures_ds_json)

        # Create config data
        config_data = {
            "structures_json": self.structures_ds_json,
            "atom_selections": self.atom_selections,
            "mode": self.mode,
            "metric_name": self.get_metric_name(),
            "unit": self.unit,
            "output_csv": self.analysis_csv
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        angle_type = {"bond": "bond", "torsion": "torsional", "vector": "vector-vector"}[self.mode]
        if self.mode == self.MODE_VECTOR:
            (a1, a2), (b1, b2) = self.atom_selections
            atoms_str = f"({a1}->{a2}) ^ ({b1}->{b2})"
        else:
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
        angle_type = {"bond": "Bond angle", "torsion": "Torsional angle",
                      "vector": "Vector-vector angle"}[self.mode]
        if self.mode == self.MODE_VECTOR:
            (a1, a2), (b1, b2) = self.atom_selections
            atoms_desc = f"({a1}->{a2}) ^ ({b1}->{b2})"
        else:
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
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()

        base_dict.update({
            "angle_params": {
                "mode": self.mode,
                "atom_selections": self.atom_selections,
                "metric_name": self.get_metric_name(),
                "unit": self.unit
            }
        })
        return base_dict
