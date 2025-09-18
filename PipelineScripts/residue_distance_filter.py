"""
ResidueDistanceFilter configuration for distance-based residue selection.

Analyzes protein structures to identify residues within or beyond a specified
distance from a reference ligand, generating PyMOL-formatted selections.
"""

import os
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo


class ResidueDistanceFilter(BaseConfig):
    """
    ResidueDistanceFilter configuration for distance-based residue selection.

    Analyzes protein structures to identify residues within or beyond a specified
    distance from a reference ligand, outputting PyMOL-formatted selections.
    """

    TOOL_NAME = "ResidueDistanceFilter"
    DEFAULT_ENV = "ProteinEnv"
    DEFAULT_RESOURCES = {"gpu": "V100", "memory": "8GB", "time": "6:00:00"}

    def __init__(self,
                 structures: Union[str, List[str], ToolOutput] = None,
                 ligand: str = "",
                 distance: float = 5.0,
                 reference_type: str = "ligand",
                 reference_selection: str = "",
                 **kwargs):
        """
        Initialize ResidueDistanceFilter configuration.

        Args:
            structures: Input structures (PDB files or ToolOutput from previous tool)
            ligand: Ligand identifier for distance reference (e.g., "LIG", "ATP")
            distance: Distance cutoff in Angstroms (default: 5.0)
            reference_type: Type of reference for distance calculation (default: "ligand")
                          Options: "ligand", "atoms", "residues"
            reference_selection: Specific selection if not using ligand
                               (e.g., "resname ATP", "resi 100-105")
            **kwargs: Additional parameters
        """
        # Store ResidueDistanceFilter-specific parameters
        self.structures = structures
        self.ligand = ligand
        self.distance = distance
        self.reference_type = reference_type
        self.reference_selection = reference_selection

        # Track input source type
        self.input_is_tool_output = isinstance(structures, ToolOutput)
        self.input_structures = structures

        # Initialize base class
        super().__init__(**kwargs)

        # Initialize file paths (will be set in configure_inputs)
        self._initialize_file_paths()

    def validate_params(self):
        """Validate ResidueDistanceFilter-specific parameters."""
        if not self.structures:
            raise ValueError("structures parameter is required")

        if not self.ligand and self.reference_type == "ligand":
            raise ValueError("ligand parameter is required when reference_type is 'ligand'")

        if self.distance <= 0:
            raise ValueError("distance must be positive")

        if self.reference_type not in ["ligand", "atoms", "residues"]:
            raise ValueError("reference_type must be 'ligand', 'atoms', or 'residues'")

        if self.reference_type != "ligand" and not self.reference_selection:
            raise ValueError("reference_selection is required when reference_type is not 'ligand'")

    def _initialize_file_paths(self):
        """Initialize common file paths used throughout the class."""
        self.selections_csv = None
        self.distance_filter_py = None

    def _setup_file_paths(self):
        """Set up all file paths after output_folder is known."""
        # Output files
        self.selections_csv = os.path.join(self.output_folder, "selections.csv")

        # Helper script paths (only set if folders are available)
        if hasattr(self, 'folders') and self.folders:
            self.distance_filter_py = os.path.join(self.folders["HelpScripts"], "pipe_residue_distance_filter.py")
        else:
            # Temporary placeholder when folders aren't available yet
            self.distance_filter_py = None

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input structures from various sources."""
        self.folders = pipeline_folders
        self._setup_file_paths()  # Set up all file paths now that we have folders

        if self.input_is_tool_output:
            # Input from previous tool (e.g., RFdiffusion, Boltz2)
            tool_output: ToolOutput = self.input_structures

            # Try to get structures - could be in various formats
            source_structures = []

            # Try different output types
            for struct_type in ["structures", "pdbs"]:
                struct_files = tool_output.get_output_files(struct_type)
                if struct_files:
                    source_structures = struct_files
                    break

            if not source_structures:
                raise ValueError(f"No structure outputs found from {tool_output.tool_type}")

            # Store source for script generation
            self.input_sources = {"structures": source_structures}

            # Add dependency
            self.dependencies.append(tool_output.config)

        elif isinstance(self.input_structures, list):
            # Direct list of structure file paths (from StandardizedOutput)
            if self.input_structures:
                self.input_sources = {"structures": self.input_structures}
            else:
                raise ValueError("Empty structure list provided")

        elif isinstance(self.input_structures, StandardizedOutput):
            # StandardizedOutput object (from tool.output)
            if self.input_structures.structures:
                self.input_sources = {"structures": self.input_structures.structures}
            else:
                raise ValueError("No structures found in StandardizedOutput")

        elif isinstance(self.input_structures, str):
            # String input - single PDB file
            if self.input_structures.endswith('.pdb'):
                pdb_source = os.path.join(pipeline_folders["PDBs"], self.input_structures)
                if os.path.exists(pdb_source):
                    self.input_sources = {"structures": [pdb_source]}
                else:
                    raise ValueError(f"PDB file not found: {pdb_source}")
            else:
                raise ValueError("String input must be a PDB file path")
        else:
            raise ValueError(f"Unsupported input type: {type(self.input_structures)}")

    def get_config_display(self) -> List[str]:
        """Get ResidueDistanceFilter configuration display lines."""
        config_lines = super().get_config_display()

        # Input information - show only structure count, not full details
        if self.input_is_tool_output:
            # Count structures from tool output
            structure_count = 0
            if hasattr(self.input_structures, 'structures') and self.input_structures.structures:
                structure_count = len(self.input_structures.structures)
            config_lines.append(f"INPUT: {self.input_structures.tool_type} output ({structure_count} structures)")
        else:
            config_lines.append(f"INPUT: {self.input_structures}")

        config_lines.extend([
            f"LIGAND: {self.ligand}",
            f"DISTANCE: {self.distance}Å",
            f"REFERENCE TYPE: {self.reference_type}",
        ])

        if self.reference_selection:
            config_lines.append(f"REFERENCE SELECTION: {self.reference_selection}")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """
        Generate ResidueDistanceFilter execution script.

        Args:
            script_path: Path where script should be written

        Returns:
            Script content as string
        """
        # Generate script content
        script_content = "#!/bin/bash\n"
        script_content += "# ResidueDistanceFilter execution script\n"
        script_content += "# Generated by ProteinNotebooks pipeline system\n\n"
        script_content += self.generate_completion_check_header()

        # Get specific structure files
        if hasattr(self, 'input_sources') and "structures" in self.input_sources:
            structure_files = self.input_sources["structures"]
            # Convert to comma-separated string for passing to script
            structure_files_str = ",".join(structure_files)
        else:
            raise ValueError("No structure sources found")

        # Determine reference specification
        if self.reference_type == "ligand":
            reference_spec = f"ligand:{self.ligand}"
        else:
            reference_spec = f"{self.reference_type}:{self.reference_selection}"

        script_content += f"""echo "Analyzing residue distances for {len(structure_files)} structures"
echo "Reference: {reference_spec}"
echo "Distance cutoff: {self.distance}Å"

# Run distance analysis
python {self.distance_filter_py} \\
    "{structure_files_str}" \\
    "{reference_spec}" \\
    {self.distance} \\
    "{self.selections_csv}"

echo "Distance analysis completed"
echo "Selections saved to: {self.selections_csv}"

"""
        script_content += self.generate_completion_check_footer()

        return script_content

    def _predict_structure_ids(self) -> List[str]:
        """
        Predict the structure IDs that will be analyzed.

        Returns:
            List of structure identifiers
        """
        structure_ids = []

        # Get structure files from input sources
        if hasattr(self, 'input_sources') and "structures" in self.input_sources:
            structure_files = self.input_sources["structures"]
            for struct_file in structure_files:
                # Extract ID from filename (remove extension)
                filename = os.path.basename(struct_file)
                if filename.endswith('.pdb'):
                    structure_ids.append(filename[:-4])
                elif filename.endswith('.cif'):
                    structure_ids.append(filename[:-4])
                else:
                    structure_ids.append(filename)

        return structure_ids

    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after ResidueDistanceFilter execution.

        Returns:
            Dictionary mapping output types to file paths with standard keys:
            - structures: Empty (no structures from ResidueDistanceFilter)
            - compounds: Empty (no compounds from ResidueDistanceFilter)
            - sequences: Empty (no sequences from ResidueDistanceFilter)
            - datasheets: Selections datasheet with distance-based residue selections
            - output_folder: Tool's output directory
        """
        # Ensure file paths are set up
        if not hasattr(self, 'selections_csv') or self.selections_csv is None:
            # Fallback if configure_inputs hasn't been called yet
            self._setup_file_paths()

        # Predict structure IDs for analysis count
        structure_ids = self._predict_structure_ids()

        # Organize datasheets by content type with detailed metadata
        datasheets = {
            "selections": DatasheetInfo(
                name="selections",
                path=self.selections_csv,
                columns=["id", "within", "beyond", "distance_cutoff", "reference_ligand"],
                description="PyMOL-formatted residue selections based on distance to ligand",
                count=len(structure_ids)
            )
        }

        return {
            "structures": [],
            "structure_ids": [],
            "compounds": [],
            "compound_ids": [],
            "sequences": [],
            "sequence_ids": [],
            "datasheets": datasheets,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including ResidueDistanceFilter-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "residue_distance_filter_params": {
                "ligand": self.ligand,
                "distance": self.distance,
                "reference_type": self.reference_type,
                "reference_selection": self.reference_selection,
                "input_type": "tool_output" if self.input_is_tool_output else "direct"
            }
        })
        return base_dict