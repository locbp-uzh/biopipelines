"""
PoseDistance - Measure ligand pose distance between reference and target structures.

Calculates RMSD and distance metrics for ligand poses between a reference
holo structure (e.g., from X-ray crystallography) and designed holo structures.
Useful for validating whether design tools reproduce known binding poses.
"""

import os
from typing import Dict, List, Any, Union, Optional

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo
    from .mixins import InputHandlerMixin, DatasheetNavigatorMixin, FilePathDescriptor
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo
    from mixins import InputHandlerMixin, DatasheetNavigatorMixin, FilePathDescriptor


class PoseDistance(InputHandlerMixin, DatasheetNavigatorMixin, BaseConfig):
    """
    Measures ligand pose distance between reference and target structures.

    Compares ligand conformations between a reference holo structure (e.g., XRC)
    and designed/predicted holo structures. Calculates RMSD, centroid distance,
    and orientation metrics for ligand poses.

    Commonly used for:
    - Validating binding pose prediction accuracy
    - Comparing designed structures to experimental references
    - Analyzing ligand pose consistency across designs
    """

    TOOL_NAME = "PoseDistance"
    DEFAULT_ENV = "ProteinEnv"

    # Automatic file path management with descriptors
    results_csv = FilePathDescriptor("pose_analysis.csv")
    config_file = FilePathDescriptor("pose_config.json")
    pose_script = FilePathDescriptor("pipe_pose_distance.py", folder_key="HelpScripts")

    def __init__(self,
                 reference: Union[str, ToolOutput, StandardizedOutput],
                 target_structures: Union[str, List[str], ToolOutput, StandardizedOutput],
                 ligand: str,
                 reference_ligand: Optional[str] = None,
                 alignment_selection: str = "protein",
                 calculate_centroid: bool = True,
                 calculate_orientation: bool = False,
                 **kwargs):
        """
        Initialize PoseDistance tool.

        Args:
            reference: Reference holo structure (XRC or designed reference)
            target_structures: Target structures to compare against reference
            ligand: Ligand residue name in target structures (e.g., 'LIG', 'ATP')
            reference_ligand: Ligand residue name in reference (default: same as ligand)
            alignment_selection: Selection for protein alignment before pose comparison
                               (default: "protein", can be PyMOL selection like "chain A")
            calculate_centroid: Calculate ligand centroid distance (default: True)
            calculate_orientation: Calculate orientation angle difference (default: False)
            **kwargs: Additional parameters passed to BaseConfig
        """
        # Store parameters
        self.reference_input = reference
        self.target_structures_input = target_structures
        self.ligand = ligand
        self.reference_ligand = reference_ligand or ligand
        self.alignment_selection = alignment_selection
        self.calculate_centroid = calculate_centroid
        self.calculate_orientation = calculate_orientation

        # Initialize base class
        super().__init__(**kwargs)

        # Set up dependencies
        if hasattr(reference, 'config'):
            self.dependencies.append(reference.config)
        if hasattr(target_structures, 'config'):
            self.dependencies.append(target_structures.config)

    def validate_params(self):
        """Validate tool parameters."""
        if not self.ligand:
            raise ValueError("ligand parameter is required")

        if not isinstance(self.ligand, str):
            raise ValueError(f"ligand must be a string, got {type(self.ligand)}")

        if len(self.ligand) < 1 or len(self.ligand) > 3:
            raise ValueError(f"ligand must be 1-3 characters (residue name), got '{self.ligand}'")

        if self.reference_ligand and (len(self.reference_ligand) < 1 or len(self.reference_ligand) > 3):
            raise ValueError(f"reference_ligand must be 1-3 characters, got '{self.reference_ligand}'")

        if not self.alignment_selection:
            raise ValueError("alignment_selection cannot be empty")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """
        Configure inputs using mixins - demonstrates new elegant patterns.
        """
        self.folders = pipeline_folders

        # ========================================================================
        # ELEGANT INPUT HANDLING - One line with resolve_input()!
        # ========================================================================
        # Resolve reference structure
        reference_resolved = self.resolve_input(self.reference_input, 'structures')
        if not reference_resolved.files:
            raise ValueError("No reference structure found")
        if len(reference_resolved.files) > 1:
            raise ValueError(
                f"Expected single reference structure, got {len(reference_resolved.files)}. "
                f"Please provide a single PDB file."
            )
        self.reference_pdb = reference_resolved.files[0]

        # Resolve target structures
        target_resolved = self.resolve_input(self.target_structures_input, 'structures')
        if not target_resolved.files:
            raise ValueError("No target structures found")
        self.target_pdbs = target_resolved.files
        self.target_ids = target_resolved.ids

        # Store datasheets if available for merging results
        self.target_datasheets = target_resolved.datasheets

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"REFERENCE: {os.path.basename(self.reference_pdb) if hasattr(self, 'reference_pdb') else 'not configured'}",
            f"TARGET STRUCTURES: {len(self.target_pdbs) if hasattr(self, 'target_pdbs') else 0}",
            f"LIGAND: {self.ligand}",
            f"REFERENCE LIGAND: {self.reference_ligand}",
            f"ALIGNMENT: {self.alignment_selection}",
            f"METRICS: RMSD" + (", centroid" if self.calculate_centroid else "") + (", orientation" if self.calculate_orientation else "")
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate PoseDistance execution script."""
        # Create config for helper script
        config_data = {
            "reference_pdb": self.reference_pdb,
            "reference_ligand": self.reference_ligand,
            "target_pdbs": self.target_pdbs,
            "target_ids": self.target_ids,
            "ligand": self.ligand,
            "alignment_selection": self.alignment_selection,
            "calculate_centroid": self.calculate_centroid,
            "calculate_orientation": self.calculate_orientation,
            "output_csv": self.results_csv
        }

        import json
        os.makedirs(self.output_folder, exist_ok=True)
        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        # Generate bash script
        script_content = f"""#!/bin/bash
# PoseDistance execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "=========================================="
echo "Pose Distance Analysis"
echo "=========================================="
echo "Reference: {os.path.basename(self.reference_pdb)}"
echo "Reference ligand: {self.reference_ligand}"
echo "Target structures: {len(self.target_pdbs)}"
echo "Target ligand: {self.ligand}"
echo "Alignment: {self.alignment_selection}"
echo "Output: {self.results_csv}"
echo "=========================================="

# Run pose distance analysis
python "{self.pose_script}" \\
  --config "{self.config_file}"

if [ $? -eq 0 ]; then
    echo "Pose analysis completed successfully"
    echo "Results written to: {self.results_csv}"
else
    echo "Error: Pose analysis failed"
    exit 1
fi

{self.generate_completion_check_footer()}
"""
        return script_content

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
            "alignment_selection"
        ])

        datasheets = {
            "analysis": DatasheetInfo(
                name="analysis",
                path=self.results_csv,
                columns=columns,
                description=f"Ligand pose distance analysis comparing {self.ligand} poses to reference",
                count=len(self.target_pdbs) if hasattr(self, 'target_pdbs') else 0
            )
        }

        return {
            "structures": [],
            "structure_ids": [],
            "sequences": [],
            "sequence_ids": [],
            "datasheets": datasheets,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                "ligand": self.ligand,
                "reference_ligand": self.reference_ligand,
                "alignment_selection": self.alignment_selection,
                "calculate_centroid": self.calculate_centroid,
                "calculate_orientation": self.calculate_orientation
            }
        })
        return base_dict
