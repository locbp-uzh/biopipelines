"""
PoseDistance - Measure ligand pose distance between reference and sample structures.

Calculates RMSD and distance metrics for ligand poses between a reference
holo structure (e.g., from X-ray crystallography) and designed holo structures.
Useful for validating whether design tools reproduce known binding poses.
"""

import os
from typing import Dict, List, Any, Union, Optional

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo


class PoseDistance(BaseConfig):
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

    TOOL_NAME = "PoseDistance"
    

    def __init__(self,
                 reference_structure: Union[str, ToolOutput, StandardizedOutput],
                 sample_structures: Union[str, List[str], ToolOutput, StandardizedOutput],
                 reference_ligand: str,
                 sample_ligand: Optional[str] = None,
                 alignment_selection: str = "protein",
                 calculate_centroid: bool = True,
                 calculate_orientation: bool = False,
                 **kwargs):
        """
        Initialize PoseDistance tool.

        Args:
            reference_structure: Reference holo structure (XRC or designed reference)
            sample_structures: Sample structures to compare against reference
            reference_ligand: Ligand residue name in reference structure (e.g., 'LIG', 'ATP')
            sample_ligand: Ligand residue name in sample structures (default: same as reference_ligand)
            alignment_selection: Selection for protein alignment before pose comparison
                               (default: "protein", can be PyMOL selection like "chain A")
            calculate_centroid: Calculate ligand centroid distance (default: True)
            calculate_orientation: Calculate orientation angle difference (default: False)
            **kwargs: Additional parameters passed to BaseConfig
        """
        # Store parameters
        self.reference_structure_input = reference_structure
        self.sample_structures_input = sample_structures
        self.reference_ligand = reference_ligand
        self.sample_ligand = sample_ligand or reference_ligand
        self.alignment_selection = alignment_selection
        self.calculate_centroid = calculate_centroid
        self.calculate_orientation = calculate_orientation

        # Initialize base class
        super().__init__(**kwargs)

        # Set up dependencies
        if hasattr(reference_structure, 'config'):
            self.dependencies.append(reference_structure.config)
        if hasattr(sample_structures, 'config'):
            self.dependencies.append(sample_structures.config)

    def get_analysis_csv_path(self) -> str:
        """Get the path for the analysis CSV file - defined once, used everywhere."""
        return os.path.join(self.output_folder, "pose_analysis.csv")

    def validate_params(self):
        """Validate tool parameters."""
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
        """
        Configure inputs from previous tools.
        """
        self.folders = pipeline_folders

        # Resolve reference structure
        self.reference_structure_path = None
        if isinstance(self.reference_structure_input, str):
            # Direct file path
            self.reference_structure_path = self.reference_structure_input
        elif hasattr(self.reference_structure_input, 'structures'):
            # From tool output
            if isinstance(self.reference_structure_input.structures, list):
                if len(self.reference_structure_input.structures) == 0:
                    raise ValueError("No reference structure found in tool output")
                if len(self.reference_structure_input.structures) > 1:
                    raise ValueError(
                        f"Expected single reference structure, got {len(self.reference_structure_input.structures)}. "
                        f"Please provide a single PDB file."
                    )
                self.reference_structure_path = self.reference_structure_input.structures[0]
            else:
                self.reference_structure_path = self.reference_structure_input.structures

        if not self.reference_structure_path:
            raise ValueError("Could not resolve reference structure path")

        # Resolve sample structures
        self.sample_structure_paths = []
        self.sample_structure_ids = []

        if isinstance(self.sample_structures_input, str):
            # Direct file path
            self.sample_structure_paths = [self.sample_structures_input]
            self.sample_structure_ids = [os.path.splitext(os.path.basename(self.sample_structures_input))[0]]
        elif isinstance(self.sample_structures_input, list):
            # List of file paths
            self.sample_structure_paths = self.sample_structures_input
            self.sample_structure_ids = [os.path.splitext(os.path.basename(p))[0] for p in self.sample_structures_input]
        elif hasattr(self.sample_structures_input, 'structures'):
            # From tool output
            if isinstance(self.sample_structures_input.structures, list):
                self.sample_structure_paths = self.sample_structures_input.structures
            else:
                self.sample_structure_paths = [self.sample_structures_input.structures]

            # Get IDs if available
            if hasattr(self.sample_structures_input, 'structure_ids'):
                self.sample_structure_ids = self.sample_structures_input.structure_ids
            else:
                self.sample_structure_ids = [os.path.splitext(os.path.basename(p))[0] for p in self.sample_structure_paths]

        if not self.sample_structure_paths:
            raise ValueError("Could not resolve sample structure paths")

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"REFERENCE: {os.path.basename(self.reference_structure_path) if hasattr(self, 'reference_structure_path') else 'not configured'}",
            f"SAMPLE STRUCTURES: {len(self.sample_structure_paths) if hasattr(self, 'sample_structure_paths') else 0}",
            f"REFERENCE LIGAND: {self.reference_ligand}",
            f"SAMPLE LIGAND: {self.sample_ligand}",
            f"ALIGNMENT: {self.alignment_selection}",
            f"METRICS: RMSD" + (", centroid" if self.calculate_centroid else "") + (", orientation" if self.calculate_orientation else "")
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate PoseDistance execution script."""
        output_folder = self.output_folder
        os.makedirs(output_folder, exist_ok=True)

        # Output CSV path - defined once in get_analysis_csv_path()
        analysis_csv = self.get_analysis_csv_path()

        # Create config for helper script
        config_file = os.path.join(output_folder, "pose_config.json")
        config_data = {
            "reference_pdb": self.reference_structure_path,
            "reference_ligand": self.reference_ligand,
            "target_pdbs": self.sample_structure_paths,
            "target_ids": self.sample_structure_ids,
            "ligand": self.sample_ligand,
            "alignment_selection": self.alignment_selection,
            "calculate_centroid": self.calculate_centroid,
            "calculate_orientation": self.calculate_orientation,
            "output_csv": analysis_csv
        }

        import json
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        # Generate bash script
        script_content = f"""#!/bin/bash
# PoseDistance execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Running pose distance analysis"
echo "Reference: {os.path.basename(self.reference_structure_path)}"
echo "Reference ligand: {self.reference_ligand}"
echo "Sample structures: {len(self.sample_structure_paths)}"
echo "Sample ligand: {self.sample_ligand}"
echo "Alignment: {self.alignment_selection}"
echo "Output: {analysis_csv}"

# Run pose distance analysis
python "{os.path.join(self.folders['HelpScripts'], 'pipe_pose_distance.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Pose analysis completed successfully"
    echo "Results written to: {analysis_csv}"
else
    echo "Error: Pose analysis failed"
    exit 1
fi

{self.generate_completion_check_footer()}
"""
        return script_content

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files."""
        analysis_csv = self.get_analysis_csv_path()

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
            "analysis": TableInfo(
                name="analysis",
                path=analysis_csv,
                columns=columns,
                description=f"Ligand pose distance analysis comparing {self.sample_ligand} poses to reference",
                count=len(self.sample_structure_paths) if hasattr(self, 'sample_structure_paths') else 0
            )
        }

        return {
            "structures": [],
            "structure_ids": [],
            "compounds": [],
            "compound_ids": [],
            "sequences": [],
            "sequence_ids": [],
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
