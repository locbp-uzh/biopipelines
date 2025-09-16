"""
FetchStructure tool for downloading protein structures from RCSB PDB.

Downloads protein structures by PDB ID in either PDB or CIF format,
with support for batch downloading and proper error handling.
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


class FetchStructure(BaseConfig):
    """
    Pipeline tool for downloading protein structures from RCSB PDB.
    
    Downloads structures by PDB ID in either PDB or CIF format, creating
    structure files and metadata datasheets for downstream processing.
    """
    
    # Tool identification
    TOOL_NAME = "FetchStructure"
    DEFAULT_ENV = "ProteinEnv"
    COMPATIBLE_ENVS = ["ProteinEnv", "Boltz2Env"]
    DEFAULT_RESOURCES = {"gpu": None, "memory": "2GB", "time": "1:00:00"}
    
    def __init__(self,
                 pdbs: Union[str, List[str]],
                 ids: Optional[Union[str, List[str]]] = None,
                 format: str = "pdb",
                 include_biological_assembly: bool = False,
                 remove_waters: bool = True,
                 **kwargs):
        """
        Initialize FetchStructure tool.

        Args:
            pdbs: PDB ID(s) to download. Can be single string or list of strings (e.g. "4ufc" or ["4ufc","1abc"])
            ids: Custom IDs for renaming. Can be single string or list of strings (e.g. "POI" or ["POI1","POI2"]). If None, uses pdbs as ids.
            format: File format ("pdb" or "cif", default: "pdb")
            include_biological_assembly: Whether to download biological assembly (default: False)
            remove_waters: Whether to remove water molecules from structures (default: True)
            **kwargs: Additional parameters

        Examples:
            # Download single structure with default naming
            fetch = pipeline.add(FetchStructure(
                pdbs="4ufc"
            ))

            # Download multiple structures with custom IDs
            fetch = pipeline.add(FetchStructure(
                pdbs=["4ufc", "1abc"],
                ids=["POI1", "POI2"],
                format="pdb"
            ))

            # Download with biological assembly and custom naming, keep waters
            fetch = pipeline.add(FetchStructure(
                pdbs="4ufc",
                ids="MyProtein",
                format="pdb",
                include_biological_assembly=True,
                remove_waters=False
            ))
        """
        # Normalize pdbs to list
        if isinstance(pdbs, str):
            self.pdb_ids = [pdbs.upper()]
        else:
            self.pdb_ids = [pdb_id.upper() for pdb_id in pdbs]

        # Handle custom IDs - default to pdb_ids if not provided
        if ids is None:
            self.custom_ids = self.pdb_ids.copy()
        else:
            if isinstance(ids, str):
                self.custom_ids = [ids]
            else:
                self.custom_ids = list(ids)

        # Validate that pdbs and ids have same length if ids provided
        if len(self.custom_ids) != len(self.pdb_ids):
            raise ValueError(f"Length mismatch: pdbs has {len(self.pdb_ids)} items but ids has {len(self.custom_ids)} items")

        self.format = format.lower()
        self.include_biological_assembly = include_biological_assembly
        self.remove_waters = remove_waters
        
        # Validate format
        if self.format not in ["pdb", "cif"]:
            raise ValueError(f"Invalid format: {self.format}. Must be 'pdb' or 'cif'")
        
        # Validate PDB IDs (basic format check)
        for pdb_id in self.pdb_ids:
            if not self._is_valid_pdb_id(pdb_id):
                raise ValueError(f"Invalid PDB ID format: {pdb_id}. Must be 4 characters (e.g., '1ABC')")
        
        # Initialize base class
        super().__init__(**kwargs)
    
    def _is_valid_pdb_id(self, pdb_id: str) -> bool:
        """Validate PDB ID format (4 characters, alphanumeric)."""
        return len(pdb_id) == 4 and pdb_id.isalnum()
    
    def validate_params(self):
        """Validate FetchStructure parameters."""
        if not self.pdb_ids:
            raise ValueError("pdbs cannot be empty")

        if not self.custom_ids:
            raise ValueError("ids cannot be empty")

        if len(self.pdb_ids) != len(self.custom_ids):
            raise ValueError("pdbs and ids must have same length")

        if self.format not in ["pdb", "cif"]:
            raise ValueError("format must be 'pdb' or 'cif'")
    
    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input parameters (no inputs for FetchStructure)."""
        self.folders = pipeline_folders
    
    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()
        
        config_lines.extend([
            f"PDB_IDS: {', '.join(self.pdb_ids)} ({len(self.pdb_ids)} structures)",
            f"CUSTOM_IDS: {', '.join(self.custom_ids)}",
            f"FORMAT: {self.format.upper()}",
            f"BIOLOGICAL_ASSEMBLY: {self.include_biological_assembly}",
            f"REMOVE_WATERS: {self.remove_waters}"
        ])
        
        return config_lines
    
    def generate_script(self, script_path: str) -> str:
        """
        Generate script to fetch structures from RCSB PDB.
        
        Args:
            script_path: Path where script should be written
            
        Returns:
            Script content as string
        """
        runtime_folder = os.path.dirname(script_path)
        output_folder = self.output_folder
        os.makedirs(output_folder, exist_ok=True)
        
        # Output files
        structures_datasheet = os.path.join(output_folder, "structures.csv")
        sequences_datasheet = os.path.join(output_folder, "sequences.csv")
        failed_datasheet = os.path.join(output_folder, "failed_downloads.csv")
        
        # Create config file for fetching
        config_file = os.path.join(output_folder, "fetch_config.json")
        config_data = {
            "pdb_ids": self.pdb_ids,
            "custom_ids": self.custom_ids,
            "format": self.format,
            "include_biological_assembly": self.include_biological_assembly,
            "remove_waters": self.remove_waters,
            "output_folder": output_folder,
            "structures_datasheet": structures_datasheet,
            "sequences_datasheet": sequences_datasheet,
            "failed_datasheet": failed_datasheet
        }
        
        import json
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)
        
        # Generate script content
        script_content = f"""#!/bin/bash
# FetchStructure execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Fetching {len(self.pdb_ids)} structures from RCSB PDB"
echo "Format: {self.format.upper()}"
echo "PDB IDs: {', '.join(self.pdb_ids)}"
echo "Custom IDs: {', '.join(self.custom_ids)}"
echo "Output folder: {output_folder}"

# Run Python structure fetching script
python "{os.path.join(self.folders['HelpScripts'], 'pipe_fetch_structure.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Successfully fetched structures"
    echo "Structures datasheet: {structures_datasheet}"
    echo "Sequences datasheet: {sequences_datasheet}"
    if [ -f "{failed_datasheet}" ]; then
        echo "Failed downloads logged: {failed_datasheet}"
    fi
else
    echo "Error: Failed to fetch structures"
    exit 1
fi

{self.generate_completion_check_footer()}
"""
        
        return script_content
    
    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after structure fetching.
        
        Returns:
            Dictionary with output file paths and datasheet information
        """
        # Generate structure file paths using custom IDs
        extension = ".pdb" if self.format == "pdb" else ".cif"
        structure_files = []
        for custom_id in self.custom_ids:
            structure_files.append(os.path.join(self.output_folder, f"{custom_id}{extension}"))

        # Structure IDs are the custom IDs
        structure_ids = self.custom_ids.copy()

        # Sequence IDs are also the custom IDs
        sequence_ids = self.custom_ids.copy()
        
        # Output datasheets
        structures_csv = os.path.join(self.output_folder, "structures.csv")
        sequences_csv = os.path.join(self.output_folder, "sequences.csv")
        failed_csv = os.path.join(self.output_folder, "failed_downloads.csv")
        
        # Define datasheets that will be created
        datasheets = {
            "structures": DatasheetInfo(
                name="structures",
                path=structures_csv,
                columns=["id", "pdb_id", "file_path", "format", "file_size", "download_date"],
                description="Successfully downloaded structure files",
                count=len(self.pdb_ids)  # May be fewer if some downloads fail
            ),
            "sequences": DatasheetInfo(
                name="sequences",
                path=sequences_csv,
                columns=["id", "sequence"],
                description="Protein sequences extracted from structures",
                count=len(self.pdb_ids)  # May be fewer if some downloads fail
            ),
            "failed": DatasheetInfo(
                name="failed",
                path=failed_csv,
                columns=["pdb_id", "error_message", "http_status", "attempted_url"],
                description="Failed structure downloads with error details",
                count="variable"
            )
        }
        
        return {
            "structures": structure_files,
            "structure_ids": structure_ids,
            "compounds": [],
            "compound_ids": [],
            "sequences": [sequences_csv],
            "sequence_ids": sequence_ids,
            "datasheets": datasheets,
            "output_folder": self.output_folder
        }
    
    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                "pdb_ids": self.pdb_ids,
                "custom_ids": self.custom_ids,
                "format": self.format,
                "include_biological_assembly": self.include_biological_assembly,
                "remove_waters": self.remove_waters
            }
        })
        return base_dict