"""
Ligand tool for fetching small molecule ligands from RCSB PDB.

Downloads ideal SDF files from RCSB and converts them to PDB format.
Fetches ligands with priority-based lookup: local_folder -> Ligands/ -> RCSB download.
Downloads are saved both in Ligands/ folder and tool output folder for reuse.
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


class Ligand(BaseConfig):
    """
    Pipeline tool for fetching small molecule ligands from RCSB PDB.

    Downloads ideal SDF files and converts to PDB format. Implements priority-based
    lookup: checks local_folder (if provided), then Ligands/ folder, then downloads
    from RCSB. Downloaded ligands are saved to both Ligands/ folder (for reuse) and
    tool output folder.
    """

    # Tool identification
    TOOL_NAME = "Ligand"
    DEFAULT_ENV = None  # Loaded from config.yaml

    def __init__(self,
                 codes: Union[str, List[str]],
                 ids: Optional[Union[str, List[str]]] = None,
                 local_folder: Optional[str] = None,
                 **kwargs):
        """
        Initialize Ligand tool.

        Args:
            codes: 3-letter ligand code(s) to fetch. Can be single string or list of strings (e.g. "ATP" or ["ATP","GDP"])
            ids: Custom IDs for renaming. Can be single string or list of strings (e.g. "ligand1" or ["ligand1","ligand2"]). If None, uses codes as ids.
            local_folder: Custom local folder to check first (before Ligands/). Default: None
            **kwargs: Additional parameters

        Fetch Priority:
            For each ligand code, searches in order:
            1. local_folder (if parameter provided)
            2. ./Ligands/ folder in repository
            3. Download from RCSB PDB (saved to both Ligands/ and output folder)

        Examples:
            # Fetch with automatic fallback (checks Ligands/, then downloads)
            lig = pipeline.add(Ligand(
                codes="ATP"
            ))

            # Fetch multiple ligands with custom IDs
            lig = pipeline.add(Ligand(
                codes=["ATP", "GDP"],
                ids=["ligand1", "ligand2"]
            ))

            # Check custom folder first, then Ligands/, then download
            lig = pipeline.add(Ligand(
                codes="PVY",
                local_folder="/path/to/my/ligands"
            ))
        """
        # Normalize codes to list - convert to uppercase for RCSB
        if isinstance(codes, str):
            self.ligand_codes = [codes.upper()]
        else:
            self.ligand_codes = [code.upper() for code in codes]

        # Handle custom IDs - default to ligand_codes if not provided
        if ids is None:
            self.custom_ids = self.ligand_codes.copy()
        else:
            if isinstance(ids, str):
                self.custom_ids = [ids]
            else:
                self.custom_ids = list(ids)

        # Validate that codes and ids have same length if ids provided
        if len(self.custom_ids) != len(self.ligand_codes):
            raise ValueError(f"Length mismatch: codes has {len(self.ligand_codes)} items but ids has {len(self.custom_ids)} items")

        self.local_folder = local_folder

        # Validate ligand code format (3-letter codes)
        for code in self.ligand_codes:
            if not code or len(code) > 3 or not code.isalnum():
                raise ValueError(f"Invalid ligand code: {code}. Must be 1-3 alphanumeric characters (e.g., 'ATP', 'GDP')")

        # Initialize base class
        super().__init__(**kwargs)

    def validate_params(self):
        """Validate Ligand parameters."""
        if not self.ligand_codes:
            raise ValueError("codes cannot be empty")

        if not self.custom_ids:
            raise ValueError("ids cannot be empty")

        if len(self.ligand_codes) != len(self.custom_ids):
            raise ValueError("codes and ids must have same length")

        for code in self.ligand_codes:
            if not code or len(code) > 3 or not code.isalnum():
                raise ValueError(f"Invalid ligand code: {code}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input parameters and check for local files."""
        self.folders = pipeline_folders

        # Check which files exist locally and which will need to be downloaded
        repo_ligands_folder = pipeline_folders.get('Ligands', '')
        self.found_locally = []
        self.needs_download = []

        for ligand_code, custom_id in zip(self.ligand_codes, self.custom_ids):
            found = False
            local_path = None

            # Check local_folder if specified
            if self.local_folder:
                local_path = os.path.join(self.local_folder, f"{ligand_code}.pdb")
                if os.path.exists(local_path):
                    self.found_locally.append((ligand_code, local_path))
                    found = True
                    print(f"  Found ligand {ligand_code} locally: {local_path}")

            # Check Ligands/ folder
            if not found and repo_ligands_folder:
                local_path = os.path.join(repo_ligands_folder, f"{ligand_code}.pdb")
                if os.path.exists(local_path):
                    self.found_locally.append((ligand_code, local_path))
                    found = True
                    print(f"  Found ligand {ligand_code} locally: {local_path}")

            if not found:
                # Not found locally - will download from RCSB
                self.needs_download.append(ligand_code)
                print(f"  Ligand {ligand_code} not found locally, will download from RCSB")

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"LIGAND_CODES: {', '.join(self.ligand_codes)} ({len(self.ligand_codes)} ligands)",
            f"CUSTOM_IDS: {', '.join(self.custom_ids)}",
            f"LOCAL_FOLDER: {self.local_folder if self.local_folder else 'None (uses Ligands/)'}"
        ])

        # Add status of files found/not found
        if hasattr(self, 'found_locally') and self.found_locally:
            config_lines.append(f"FOUND_LOCALLY: {', '.join([code for code, _ in self.found_locally])}")
        if hasattr(self, 'needs_download') and self.needs_download:
            config_lines.append(f"NEEDS_DOWNLOAD: {', '.join(self.needs_download)}")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """
        Generate script to fetch ligands from local folders or RCSB PDB.

        Args:
            script_path: Path where script should be written

        Returns:
            Script content as string
        """
        runtime_folder = os.path.dirname(script_path)
        output_folder = self.output_folder
        os.makedirs(output_folder, exist_ok=True)

        # Output files
        compounds_table = os.path.join(output_folder, "compounds.csv")
        failed_table = os.path.join(output_folder, "failed_downloads.csv")

        # Get Ligands folder path from folder manager
        repo_ligands_folder = self.folders['Ligands']

        # Create config file for fetching
        config_file = os.path.join(output_folder, "fetch_config.json")
        config_data = {
            "ligand_codes": self.ligand_codes,
            "custom_ids": self.custom_ids,
            "local_folder": self.local_folder,
            "repo_ligands_folder": repo_ligands_folder,
            "output_folder": output_folder,
            "compounds_table": compounds_table,
            "failed_table": failed_table
        }

        import json
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        # Generate script content
        script_content = f"""#!/bin/bash
# Ligand execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Fetching {len(self.ligand_codes)} ligands"
echo "Ligand codes: {', '.join(self.ligand_codes)}"
echo "Custom IDs: {', '.join(self.custom_ids)}"
echo "Priority: {'local_folder -> ' if self.local_folder else ''}Ligands/ -> RCSB download"
echo "Output folder: {output_folder}"

# Run Python ligand fetching script
python "{os.path.join(self.folders['HelpScripts'], 'pipe_ligand.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Successfully fetched ligands"
    echo "Compounds table: {compounds_table}"
    if [ -f "{failed_table}" ]; then
        echo "Failed downloads logged: {failed_table}"
    fi
else
    echo "Error: Failed to fetch ligands"
    exit 1
fi

{self.generate_completion_check_footer()}
"""

        return script_content

    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after ligand fetching.

        Returns:
            Dictionary with output file paths and table information
        """
        # Generate structure file paths using custom IDs
        structure_files = []
        for custom_id in self.custom_ids:
            structure_files.append(os.path.join(self.output_folder, f"{custom_id}.pdb"))

        # Structure IDs are the custom IDs
        structure_ids = self.custom_ids.copy()

        # Output tables
        compounds_csv = os.path.join(self.output_folder, "compounds.csv")
        failed_csv = os.path.join(self.output_folder, "failed_downloads.csv")

        # Define tables that will be created
        tables = {
            "compounds": TableInfo(
                name="compounds",
                path=compounds_csv,
                columns=["id", "code", "file_path", "format", "smiles", "ccd"],
                description="Successfully fetched ligand files with SMILES from RCSB",
                count=len(self.ligand_codes)  # May be fewer if some downloads fail
            ),
            "failed": TableInfo(
                name="failed",
                path=failed_csv,
                columns=["ligand_code", "error_message", "source", "attempted_path"],
                description="Failed ligand fetches with error details",
                count="variable"
            )
        }

        return {
            "structures": structure_files,
            "structure_ids": structure_ids,
            "compounds": [compounds_csv],
            "compound_ids": structure_ids,  # Same as structure_ids since each structure is a compound
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                "ligand_codes": self.ligand_codes,
                "custom_ids": self.custom_ids,
                "local_folder": self.local_folder
            }
        })
        return base_dict
