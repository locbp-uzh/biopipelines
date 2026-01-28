"""
DistanceSelector configuration for distance-based residue selection.

Analyzes protein structures to identify residues within or beyond a specified
distance from a reference ligand, generating PyMOL-formatted selections.
"""

import os
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo


class DistanceSelector(BaseConfig):
    """
    DistanceSelector configuration for distance-based residue selection.

    Analyzes protein structures to identify residues within or beyond a specified
    distance from a reference ligand, outputting PyMOL-formatted selections.
    """

    TOOL_NAME = "DistanceSelector"
    

    def __init__(self,
                 structures: Union[str, List[str], ToolOutput] = None,
                 ligand: str = "",
                 distance: float = 5.0,
                 reference: str = "ligand",
                 residues: str = "",
                 restrict_to: Union[str, tuple, None] = None,
                 id_map: Dict[str, str] = {"*": "*_<N>"},
                 include_reference: bool = True,
                 # Deprecated parameters for backward compatibility
                 reference_type: str = None,
                 reference_selection: str = None,
                 **kwargs):
        """
        Initialize DistanceSelector configuration.

        Args:
            structures: Input structures (PDB files or ToolOutput from previous tool)
            ligand: Ligand identifier for distance reference (e.g., "LIG", "ATP")
            distance: Distance cutoff in Angstroms (default: 5.0)
            reference: Type of reference for distance calculation (default: "ligand")
                      Options: "ligand", "residues"
            residues: PyMOL-style residue selection when reference="residues"
                     (e.g., "87-100", "10+15+20-25")
            restrict_to: Optional selection to restrict distance search to.
                        Accepts:
                        - Table reference tuple: (table, "column")
                        - Direct selection string: "10-20+30-40"
                        - None: Consider all protein residues (default)
            id_map: ID mapping pattern for matching structure IDs to table IDs (default: {"*": "*_<N>"})
                   - Used when table IDs don't match structure IDs
                   - Example: structure ID "rifampicin_1_2" maps to table ID "rifampicin_1"
                   - Pattern {"*": "*_<N>"} strips last "_<number>" from structure ID
                   - Set to {"*": "*"} for no mapping (1:1 ID match)
            include_reference: Whether to include reference residues in "within" selection (default: True)
                             - True: "within" includes reference + nearby residues
                             - False: "within" includes only nearby residues (excludes reference itself)
            reference_type: DEPRECATED - use 'reference' instead
            reference_selection: DEPRECATED - use 'residues' instead
            **kwargs: Additional parameters
        """
        # Handle deprecated parameters
        if reference_type is not None:
            reference = reference_type
        if reference_selection is not None:
            residues = reference_selection

        # Store DistanceSelector-specific parameters
        self.structures = structures
        self.ligand = ligand
        self.distance = distance
        self.reference = reference
        self.residues = residues
        self.restrict_to_selection = restrict_to
        self.id_map = id_map
        self.include_reference = include_reference

        # Track input source type
        self.input_is_tool_output = isinstance(structures, ToolOutput)
        self.input_structures = structures

        # Initialize base class
        super().__init__(**kwargs)

        # Track dependency if restrict_to_selection is a table reference
        if isinstance(restrict_to, tuple) and len(restrict_to) == 2:
            table_obj, _ = restrict_to
            if hasattr(table_obj, 'config'):
                self.dependencies.append(table_obj.config)

        # Initialize file paths (will be set in configure_inputs)
        self._initialize_file_paths()

    def validate_params(self):
        """Validate DistanceSelector-specific parameters."""
        if not self.structures:
            raise ValueError("structures parameter is required")

        if not self.ligand and self.reference == "ligand":
            raise ValueError("ligand parameter is required when reference is 'ligand'")

        if self.distance <= 0:
            raise ValueError("distance must be positive")

        if self.reference not in ["ligand", "residues"]:
            raise ValueError("reference must be 'ligand' or 'residues'")

        if self.reference == "residues" and not self.residues:
            raise ValueError("residues parameter is required when reference is 'residues'")

        # Validate restrict_to_selection if provided
        if self.restrict_to_selection is not None:
            if isinstance(self.restrict_to_selection, tuple):
                if len(self.restrict_to_selection) != 2:
                    raise ValueError("Table reference must be a tuple of (table, column_name)")
            elif not isinstance(self.restrict_to_selection, str):
                raise ValueError("restrict_to_selection must be a string, tuple, or None")

    def _initialize_file_paths(self):
        """Initialize common file paths used throughout the class."""
        self.selections_csv = None
        self.distance_filter_py = None

    def _setup_file_paths(self):
        """Set up all file paths after output_folder is known."""
        # Output files
        self.selections_csv = os.path.join(self.output_folder, "selections.csv")
        # Input list file (avoids "Argument list too long" with many structures)
        self.structures_list_file = os.path.join(self.output_folder, ".input_structures.txt")

        # Helper script paths (only set if folders are available)
        if hasattr(self, 'folders') and self.folders:
            self.distance_selector_py = os.path.join(self.folders["HelpScripts"], "pipe_distance_selector.py")
        else:
            # Temporary placeholder when folders aren't available yet
            self.distance_selector_py = None

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
            # StandardizedOutput object (from tool)
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
        """Get DistanceSelector configuration display lines."""
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
            f"REFERENCE: {self.reference}",
            f"DISTANCE: {self.distance}Å",
        ])

        if self.reference == "ligand":
            config_lines.append(f"LIGAND: {self.ligand}")
        elif self.reference == "residues":
            config_lines.append(f"RESIDUES: {self.residues}")

        config_lines.append(f"INCLUDE REFERENCE: {self.include_reference}")

        if self.restrict_to_selection is not None:
            if isinstance(self.restrict_to_selection, tuple):
                table_obj, column = self.restrict_to_selection
                config_lines.append(f"RESTRICT TO: {column} from table")
            else:
                config_lines.append(f"RESTRICT TO: {self.restrict_to_selection}")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """
        Generate DistanceSelector execution script.

        Args:
            script_path: Path where script should be written

        Returns:
            Script content as string
        """
        # Generate script content
        script_content = "#!/bin/bash\n"
        script_content += "# DistanceSelector execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()

        # Get specific structure files
        if hasattr(self, 'input_sources') and "structures" in self.input_sources:
            structure_files = self.input_sources["structures"]
        else:
            raise ValueError("No structure sources found")

        # Write structure paths to list file (avoids "Argument list too long" error)
        with open(self.structures_list_file, 'w') as f:
            for struct in structure_files:
                f.write(f"{struct}\n")

        # Determine reference specification
        if self.reference == "ligand":
            reference_spec = f"ligand:{self.ligand}"
        else:
            reference_spec = f"{self.reference}:{self.residues}"

        # Resolve restrict_to_selection using base class method
        if self.restrict_to_selection is not None:
            restrict_spec = self.resolve_table_reference(self.restrict_to_selection)
        else:
            restrict_spec = ""  # Empty means no restriction

        restrict_echo = f'echo "Restricting to selection: {restrict_spec}"' if restrict_spec else ""

        # Serialize id_map as JSON string for passing to script
        import json
        id_map_json = json.dumps(self.id_map)

        # Convert include_reference to string for bash
        include_reference_str = "true" if self.include_reference else "false"

        script_content += f"""echo "Analyzing residue distances for {len(structure_files)} structures"
echo "Reference: {reference_spec}"
echo "Distance cutoff: {self.distance}Å"
echo "Include reference: {self.include_reference}"
{restrict_echo}

# Run distance analysis
python {self.distance_selector_py} "{self.structures_list_file}" "{reference_spec}" {self.distance} "{restrict_spec}" "{self.selections_csv}" '{id_map_json}' {include_reference_str}

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
        Get expected output files after DistanceSelector execution.

        Returns:
            Dictionary mapping output types to file paths with standard keys:
            - structures: Empty (no structures from DistanceSelector)
            - compounds: Empty (no compounds from DistanceSelector)
            - sequences: Empty (no sequences from DistanceSelector)
            - tables: Selections table with distance-based residue selections
            - output_folder: Tool's output directory
        """
        # Ensure file paths are set up
        if not hasattr(self, 'selections_csv') or self.selections_csv is None:
            # Fallback if configure_inputs hasn't been called yet
            self._setup_file_paths()

        # Predict structure IDs for analysis count
        structure_ids = self._predict_structure_ids()

        # Organize tables by content type with detailed metadata
        tables = {
            "selections": TableInfo(
                name="selections",
                path=self.selections_csv,
                columns=["id", "pdb", "within", "beyond", "distance_cutoff", "reference_ligand"],
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
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including DistanceSelector-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "distance_selector_params": {
                "ligand": self.ligand,
                "distance": self.distance,
                "reference": self.reference,
                "residues": self.residues,
                "restrict_to_selection": str(self.restrict_to_selection) if self.restrict_to_selection else None,
                "include_reference": self.include_reference,
                "input_type": "tool_output" if self.input_is_tool_output else "direct"
            }
        })
        return base_dict