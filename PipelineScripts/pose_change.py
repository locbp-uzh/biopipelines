# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
PoseChange - Measure ligand pose distance between reference and sample structures.

Calculates RMSD and distance metrics for ligand poses between a reference
holo structure (e.g., from X-ray crystallography) and designed holo structures.
Useful for validating whether design tools reproduce known binding poses.
"""

import os
from typing import Dict, List, Any, Union, Optional

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


class PoseChange(BaseConfig):
    """
    Measures ligand pose distance between reference and sample structures.

    Compares ligand conformations between a reference holo structure (e.g., XRC)
    and designed/predicted holo structures. Calculates RMSD, centroid distance,
    and orientation metrics for ligand poses.

    Commonly used for:
    - Validating binding pose prediction accuracy
    - Comparing designed structures to experimental references
    - Analyzing ligand pose consistency across designs
    """

    TOOL_NAME = "PoseChange"

    # Lazy path descriptors
    analysis_csv = Path(lambda self: os.path.join(self.output_folder, "pose_analysis.csv"))
    config_file = Path(lambda self: os.path.join(self.output_folder, "pose_config.json"))
    samples_ds_json = Path(lambda self: os.path.join(self.output_folder, "samples_structures.json"))
    pose_change_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_pose_change.py"))

    def __init__(self,
                 reference_structure: Union[DataStream, StandardizedOutput],
                 sample_structures: Union[DataStream, StandardizedOutput],
                 reference_ligand: str,
                 sample_ligand: Optional[str] = None,
                 alignment_selection: str = "protein",
                 calculate_centroid: bool = True,
                 calculate_orientation: bool = False,
                 **kwargs):
        """
        Initialize PoseChange tool.

        Args:
            reference_structure: Reference holo structure as DataStream or StandardizedOutput
            sample_structures: Sample structures to compare against reference as DataStream or StandardizedOutput
            reference_ligand: Ligand residue name in reference structure (e.g., 'LIG', 'ATP')
            sample_ligand: Ligand residue name in sample structures (default: same as reference_ligand)
            alignment_selection: Selection for protein alignment before pose comparison
                               (default: "protein", can be PyMOL selection like "chain A")
            calculate_centroid: Calculate ligand centroid distance (default: True)
            calculate_orientation: Calculate orientation angle difference (default: False)
            **kwargs: Additional parameters passed to BaseConfig
        """
        # Resolve reference structure to DataStream
        if isinstance(reference_structure, StandardizedOutput):
            self.reference_stream: DataStream = reference_structure.streams.structures
        elif isinstance(reference_structure, DataStream):
            self.reference_stream = reference_structure
        else:
            raise ValueError(f"reference_structure must be DataStream or StandardizedOutput, got {type(reference_structure)}")

        # Resolve sample structures to DataStream
        if isinstance(sample_structures, StandardizedOutput):
            self.samples_stream: DataStream = sample_structures.streams.structures
        elif isinstance(sample_structures, DataStream):
            self.samples_stream = sample_structures
        else:
            raise ValueError(f"sample_structures must be DataStream or StandardizedOutput, got {type(sample_structures)}")

        self.reference_ligand = reference_ligand
        self.sample_ligand = sample_ligand or reference_ligand
        self.alignment_selection = alignment_selection
        self.calculate_centroid = calculate_centroid
        self.calculate_orientation = calculate_orientation

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate tool parameters."""
        if not self.reference_stream or len(self.reference_stream) == 0:
            raise ValueError("reference_structure cannot be empty")

        if len(self.reference_stream) != 1:
            raise ValueError(f"Expected single reference structure, got {len(self.reference_stream)}")

        if not self.samples_stream or len(self.samples_stream) == 0:
            raise ValueError("sample_structures cannot be empty")

        if not self.reference_ligand:
            raise ValueError("reference_ligand parameter is required")

        if not isinstance(self.reference_ligand, str):
            raise ValueError(f"reference_ligand must be a string, got {type(self.reference_ligand)}")

        if len(self.reference_ligand) < 1 or len(self.reference_ligand) > 3:
            raise ValueError(f"reference_ligand must be 1-3 characters (residue name), got '{self.reference_ligand}'")

        if self.sample_ligand and (len(self.sample_ligand) < 1 or len(self.sample_ligand) > 3):
            raise ValueError(f"sample_ligand must be 1-3 characters, got '{self.sample_ligand}'")

        if not self.alignment_selection:
            raise ValueError("alignment_selection cannot be empty")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input structures."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"REFERENCE: {len(self.reference_stream)} structure",
            f"SAMPLE STRUCTURES: {len(self.samples_stream)}",
            f"REFERENCE LIGAND: {self.reference_ligand}",
            f"SAMPLE LIGAND: {self.sample_ligand}",
            f"ALIGNMENT: {self.alignment_selection}",
            f"METRICS: RMSD" + (", centroid" if self.calculate_centroid else "") + (", orientation" if self.calculate_orientation else "")
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate PoseChange execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        script_content = "#!/bin/bash\n"
        script_content += "# PoseChange execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self.generate_script_run_pose_change()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_pose_change(self) -> str:
        """Generate the pose distance analysis part of the script."""
        import json

        # Serialize sample structures DataStream to JSON for HelpScript to load
        with open(self.samples_ds_json, 'w') as f:
            json.dump(self.samples_stream.to_dict(), f, indent=2)

        config_data = {
            "reference_pdb": self.reference_stream.files[0],
            "reference_ligand": self.reference_ligand,
            "samples_json": self.samples_ds_json,
            "ligand": self.sample_ligand,
            "alignment_selection": self.alignment_selection,
            "calculate_centroid": self.calculate_centroid,
            "calculate_orientation": self.calculate_orientation,
            "output_csv": self.analysis_csv
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        return f"""echo "Running pose distance analysis"
echo "Reference: {os.path.basename(self.reference_stream.files[0])}"
echo "Reference ligand: {self.reference_ligand}"
echo "Sample structures: {len(self.samples_stream)}"
echo "Sample ligand: {self.sample_ligand}"
echo "Alignment: {self.alignment_selection}"
echo "Output: {self.analysis_csv}"

python "{self.pose_change_py}" --config "{self.config_file}"

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files."""
        # Determine output columns based on metrics
        columns = [
            "id",
            "target_structure",
            "reference_structure",
            "ligand_rmsd"
        ]

        if self.calculate_centroid:
            columns.append("centroid_distance")

        if self.calculate_orientation:
            columns.extend(["orientation_angle", "orientation_axis"])

        columns.extend([
            "alignment_rmsd",
            "num_ligand_atoms",
            "alignment_method"
        ])

        tables = {
            "changes": TableInfo(
                name="changes",
                path=self.analysis_csv,
                columns=columns,
                description=f"Ligand pose distance analysis comparing {self.sample_ligand} poses to reference",
                count=len(self.samples_stream)
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
            "tool_params": {
                "sample_ligand": self.sample_ligand,
                "reference_ligand": self.reference_ligand,
                "alignment_selection": self.alignment_selection,
                "calculate_centroid": self.calculate_centroid,
                "calculate_orientation": self.calculate_orientation
            }
        })
        return base_dict
