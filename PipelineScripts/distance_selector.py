"""
DistanceSelector configuration for distance-based residue selection.

Analyzes protein structures to identify residues within or beyond a specified
distance from a reference ligand, generating PyMOL-formatted selections.
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


class DistanceSelector(BaseConfig):
    """
    DistanceSelector configuration for distance-based residue selection.

    Analyzes protein structures to identify residues within or beyond a specified
    distance from a reference ligand, outputting PyMOL-formatted selections.
    """

    TOOL_NAME = "DistanceSelector"

    # Lazy path descriptors
    selections_csv = Path(lambda self: os.path.join(self.output_folder, "selections.csv"))
    structures_list_file = Path(lambda self: os.path.join(self.output_folder, ".input_structures.txt"))
    distance_selector_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_distance_selector.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 ligand: str = "",
                 distance: float = 5.0,
                 reference: str = "ligand",
                 residues: str = "",
                 restrict_to: Union[str, tuple, None] = None,
                 id_map: Dict[str, str] = {"*": "*_<N>"},
                 include_reference: bool = True,
                 **kwargs):
        """
        Initialize DistanceSelector configuration.

        Args:
            structures: Input structures as DataStream or StandardizedOutput
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
            **kwargs: Additional parameters
        """
        # Resolve input to DataStream
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        # Store DistanceSelector-specific parameters
        self.ligand = ligand
        self.distance = distance
        self.reference = reference
        self.residues = residues
        self.restrict_to_selection = restrict_to
        self.id_map = id_map
        self.include_reference = include_reference

        # Initialize base class
        super().__init__(**kwargs)

        # Track dependency if restrict_to_selection is a table reference
        if isinstance(restrict_to, tuple) and len(restrict_to) == 2:
            table_obj, _ = restrict_to
            if hasattr(table_obj, 'config'):
                self.dependencies.append(table_obj.config)

    def validate_params(self):
        """Validate DistanceSelector-specific parameters."""
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")

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

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input structures."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get DistanceSelector configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"INPUT STRUCTURES: {len(self.structures_stream)} files",
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
        # Write structure paths to list file (avoids "Argument list too long" error)
        with open(self.structures_list_file, 'w') as f:
            for struct_path in self.structures_stream.files:
                f.write(f"{struct_path}\n")

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

        script_content = "#!/bin/bash\n"
        script_content += "# DistanceSelector execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += f"""echo "Analyzing residue distances for {len(self.structures_stream)} structures"
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

    def get_output_files(self) -> Dict[str, Any]:
        """
        Get expected output files after DistanceSelector execution.

        Returns:
            Dictionary with DataStream objects and tables
        """
        # Organize tables by content type with detailed metadata
        tables = {
            "selections": TableInfo(
                name="selections",
                path=self.selections_csv,
                columns=["id", "pdb", "within", "beyond", "distance_cutoff", "reference_ligand"],
                description="PyMOL-formatted residue selections based on distance to ligand",
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
        """Serialize configuration including DistanceSelector-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "distance_selector_params": {
                "ligand": self.ligand,
                "distance": self.distance,
                "reference": self.reference,
                "residues": self.residues,
                "restrict_to_selection": str(self.restrict_to_selection) if self.restrict_to_selection else None,
                "include_reference": self.include_reference
            }
        })
        return base_dict
