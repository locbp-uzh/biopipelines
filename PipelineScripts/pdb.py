"""
PDB tool for fetching protein structures from local folders or RCSB PDB.

Fetches structures with priority-based lookup: local_folder -> PDBs/ -> RCSB download.
Downloads are saved both in PDBs/ folder and tool output folder for reuse.
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


class PDB(BaseConfig):
    """
    Pipeline tool for fetching protein structures from local folders or RCSB PDB.

    Implements priority-based lookup: checks local_folder (if provided), then PDBs/
    folder, then downloads from RCSB. Downloaded structures are saved to both PDBs/
    folder (for reuse) and tool output folder.
    """

    # Tool identification
    TOOL_NAME = "PDB"
    DEFAULT_ENV = None  # Loaded from config.yaml
    
    def __init__(self,
                 pdbs: Union[str, List[str]],
                 ids: Optional[Union[str, List[str]]] = None,
                 format: str = "pdb",
                 local_folder: Optional[str] = None,
                 biological_assembly: bool = False,
                 remove_waters: bool = True,
                 **kwargs):
        """
        Initialize PDB tool.

        Args:
            pdbs: PDB ID(s) to fetch. Can be single string or list of strings (e.g. "4ufc" or ["4ufc","1abc"])
            ids: Custom IDs for renaming. Can be single string or list of strings (e.g. "POI" or ["POI1","POI2"]). If None, uses pdbs as ids.
            format: File format ("pdb" or "cif", default: "pdb")
            local_folder: Custom local folder to check first (before PDBs/). Default: None
            biological_assembly: Whether to download biological assembly from RCSB (default: False)
            remove_waters: Whether to remove water molecules from structures (default: True)
            **kwargs: Additional parameters

        Fetch Priority:
            For each PDB ID, searches in order:
            1. local_folder (if parameter provided)
            2. ./PDBs/ folder in repository
            3. Download from RCSB PDB (saved to both PDBs/ and output folder)

        Examples:
            # Fetch with automatic fallback (checks PDBs/, then downloads)
            pdb = pipeline.add(PDB(
                pdbs="4ufc"
            ))

            # Fetch multiple structures with custom IDs
            pdb = pipeline.add(PDB(
                pdbs=["4ufc", "1abc"],
                ids=["POI1", "POI2"],
                format="pdb"
            ))

            # Check custom folder first, then PDBs/, then download
            pdb = pipeline.add(PDB(
                pdbs="my_structure",
                local_folder="/path/to/my/pdbs"
            ))

            # Download with biological assembly and custom naming, keep waters
            pdb = pipeline.add(PDB(
                pdbs="4ufc",
                ids="MyProtein",
                format="pdb",
                biological_assembly=True,
                remove_waters=False
            ))
        """
        # Normalize pdbs to list - preserve original case for local file lookups
        if isinstance(pdbs, str):
            self.pdb_ids = [pdbs]
        else:
            self.pdb_ids = list(pdbs)

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
        self.local_folder = local_folder
        self.biological_assembly = biological_assembly
        self.remove_waters = remove_waters

        # Validate format
        if self.format not in ["pdb", "cif"]:
            raise ValueError(f"Invalid format: {self.format}. Must be 'pdb' or 'cif'")

        # Note: PDB ID format validation is skipped at init time because:
        # 1. local_folder may contain custom-named files
        # 2. PDBs/ folder may contain custom-named files
        # 3. Runtime script will validate RCSB format only if download is needed

        # Initialize base class
        super().__init__(**kwargs)
    
    def validate_params(self):
        """Validate PDB parameters."""
        if not self.pdb_ids:
            raise ValueError("pdbs cannot be empty")

        if not self.custom_ids:
            raise ValueError("ids cannot be empty")

        if len(self.pdb_ids) != len(self.custom_ids):
            raise ValueError("pdbs and ids must have same length")

        if self.format not in ["pdb", "cif"]:
            raise ValueError("format must be 'pdb' or 'cif'")
    
    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input parameters and check for local files."""
        self.folders = pipeline_folders

        # Check which files exist locally and which will need to be downloaded
        repo_pdbs_folder = pipeline_folders.get('PDBs', '')
        self.found_locally = []
        self.needs_download = []

        for pdb_id, custom_id in zip(self.pdb_ids, self.custom_ids):
            extension = ".pdb" if self.format == "pdb" else ".cif"
            found = False
            local_path = None

            # Check local_folder if specified
            if self.local_folder:
                local_path = os.path.join(self.local_folder, f"{pdb_id}{extension}")
                if os.path.exists(local_path):
                    self.found_locally.append((pdb_id, local_path))
                    found = True

            # Check PDBs/ folder
            if not found and repo_pdbs_folder:
                local_path = os.path.join(repo_pdbs_folder, f"{pdb_id}{extension}")
                if os.path.exists(local_path):
                    self.found_locally.append((pdb_id, local_path))
                    found = True

            if found:
                # Check if it's a valid RCSB ID and query for ligands
                rcsb_id = pdb_id.upper()
                if len(rcsb_id) == 4 and rcsb_id.isalnum():
                    # Valid RCSB format - check for ligands
                    try:
                        has_ligands = self._check_rcsb_exists_silent(rcsb_id, custom_id)
                        if has_ligands:
                            print(f"  Found PDB {pdb_id} locally: {local_path} (contains ligands)")
                        else:
                            print(f"  Found PDB {pdb_id} locally: {local_path}")
                    except:
                        # RCSB query failed, just show found locally
                        print(f"  Found PDB {pdb_id} locally: {local_path}")
                else:
                    # Custom file, not an RCSB ID
                    print(f"  Found PDB {pdb_id} locally: {local_path}")
            else:
                # Not found locally - check if valid RCSB ID
                rcsb_id = pdb_id.upper()
                if len(rcsb_id) != 4 or not rcsb_id.isalnum():
                    raise ValueError(f"PDB '{pdb_id}' not found locally and is not a valid RCSB PDB ID (must be 4 alphanumeric characters)")

                # Check if exists on RCSB
                has_ligands = self._check_rcsb_exists(rcsb_id, custom_id)
                self.needs_download.append(pdb_id)

                if has_ligands:
                    print(f"  PDB {pdb_id} not found locally, will download from RCSB (contains ligands)")
                else:
                    print(f"  PDB {pdb_id} not found locally, will download from RCSB")

        # Query RCSB API to get ligand information for all structures
        self._fetch_ligand_info_from_rcsb()

    def _check_ligands_in_rcsb(self, rcsb_id: str, custom_id: str = None) -> tuple:
        """
        Check RCSB for ligands (shared logic).

        Args:
            rcsb_id: RCSB PDB ID (uppercase, 4 characters)
            custom_id: Custom ID for compound naming (optional)

        Returns:
            Tuple of (has_ligands: bool, ligand_codes: list)
        """
        try:
            import requests
        except ImportError:
            return False, []

        try:
            url = f"https://data.rcsb.org/rest/v1/core/entry/{rcsb_id}"
            response = requests.get(url, timeout=10)
            response.raise_for_status()

            data = response.json()

            # Check for ligands
            if 'rcsb_entry_info' in data:
                entry_info = data['rcsb_entry_info']
                if 'nonpolymer_bound_components' in entry_info:
                    ligands = entry_info['nonpolymer_bound_components']
                    # Filter out common solvents/ions
                    common_solvents = {'HOH', 'WAT', 'H2O', 'NA', 'CL', 'CA', 'MG', 'K', 'ZN', 'SO4', 'PO4', 'GOL', 'EDO'}
                    real_ligands = [lig for lig in ligands if lig not in common_solvents]

                    # Store predicted compound IDs if custom_id provided
                    if custom_id and real_ligands:
                        for ligand_code in real_ligands:
                            self.predicted_compound_ids.append(f"{custom_id}_{ligand_code}")

                    return len(real_ligands) > 0, real_ligands

            return False, []

        except Exception:
            return False, []

    def _check_rcsb_exists_silent(self, rcsb_id: str, custom_id: str = None) -> bool:
        """
        Check if PDB has ligands on RCSB (silent, no exceptions).

        Args:
            rcsb_id: RCSB PDB ID (uppercase, 4 characters)
            custom_id: Custom ID for compound naming (optional)

        Returns:
            True if has ligands, False otherwise
        """
        has_ligands, _ = self._check_ligands_in_rcsb(rcsb_id, custom_id)
        return has_ligands

    def _check_rcsb_exists(self, rcsb_id: str, custom_id: str = None) -> bool:
        """
        Check if PDB exists on RCSB and return if it has ligands.

        Args:
            rcsb_id: RCSB PDB ID (uppercase, 4 characters)
            custom_id: Custom ID for compound naming (optional)

        Returns:
            True if has ligands, False if no ligands

        Raises:
            ValueError: If PDB doesn't exist on RCSB
        """
        try:
            import requests
        except ImportError:
            # Can't check, assume it exists
            return False

        try:
            url = f"https://data.rcsb.org/rest/v1/core/entry/{rcsb_id}"
            response = requests.get(url, timeout=10)
            response.raise_for_status()

            # Use shared logic for ligand checking
            has_ligands, _ = self._check_ligands_in_rcsb(rcsb_id, custom_id)
            return has_ligands

        except requests.exceptions.HTTPError as e:
            if e.response.status_code == 404:
                raise ValueError(f"PDB '{rcsb_id}' not found on RCSB (URL: {url})")
            else:
                raise ValueError(f"Error checking RCSB for '{rcsb_id}': {e}")
        except Exception as e:
            raise ValueError(f"Error checking RCSB for '{rcsb_id}': {e}")

    def _fetch_ligand_info_from_rcsb(self):
        """Fetch ligand information from RCSB API at pipeline runtime (silently for local files)."""
        self.predicted_compound_ids = []

        # Skip - ligand info already handled in _check_rcsb_exists for downloads
        # For local files, don't query RCSB (they may be custom files)
        return

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"PDB_IDS: {', '.join(self.pdb_ids)} ({len(self.pdb_ids)} structures)",
            f"CUSTOM_IDS: {', '.join(self.custom_ids)}",
            f"FORMAT: {self.format.upper()}",
            f"LOCAL_FOLDER: {self.local_folder if self.local_folder else 'None (uses PDBs/)'}",
            f"BIOLOGICAL_ASSEMBLY: {self.biological_assembly}",
            f"REMOVE_WATERS: {self.remove_waters}"
        ])

        # Add status of files found/not found
        if hasattr(self, 'found_locally') and self.found_locally:
            config_lines.append(f"FOUND_LOCALLY: {', '.join([pdb_id for pdb_id, _ in self.found_locally])}")
        if hasattr(self, 'needs_download') and self.needs_download:
            config_lines.append(f"NEEDS_DOWNLOAD: {', '.join(self.needs_download)}")

        return config_lines
    
    def generate_script(self, script_path: str) -> str:
        """
        Generate script to fetch structures from local folders or RCSB PDB.

        Args:
            script_path: Path where script should be written

        Returns:
            Script content as string
        """
        runtime_folder = os.path.dirname(script_path)
        output_folder = self.output_folder
        os.makedirs(output_folder, exist_ok=True)

        # Output files
        structures_table = os.path.join(output_folder, "structures.csv")
        sequences_table = os.path.join(output_folder, "sequences.csv")
        failed_table = os.path.join(output_folder, "failed_downloads.csv")
        compounds_table = os.path.join(output_folder, "compounds.csv")

        # Get PDBs folder path from folder manager
        repo_pdbs_folder = self.folders['PDBs']

        # Create config file for fetching
        config_file = os.path.join(output_folder, "fetch_config.json")
        config_data = {
            "pdb_ids": self.pdb_ids,
            "custom_ids": self.custom_ids,
            "format": self.format,
            "local_folder": self.local_folder,
            "repo_pdbs_folder": repo_pdbs_folder,
            "biological_assembly": self.biological_assembly,
            "remove_waters": self.remove_waters,
            "output_folder": output_folder,
            "structures_table": structures_table,
            "sequences_table": sequences_table,
            "failed_table": failed_table,
            "compounds_table": compounds_table
        }

        import json
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        # Generate script content
        script_content = f"""#!/bin/bash
# PDB execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Fetching {len(self.pdb_ids)} structures"
echo "Format: {self.format.upper()}"
echo "PDB IDs: {', '.join(self.pdb_ids)}"
echo "Custom IDs: {', '.join(self.custom_ids)}"
echo "Priority: {'local_folder -> ' if self.local_folder else ''}PDBs/ -> RCSB download"
echo "Output folder: {output_folder}"

# Run Python structure fetching script
python "{os.path.join(self.folders['HelpScripts'], 'pipe_pdb.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Successfully fetched structures"
    echo "Structures table: {structures_table}"
    echo "Sequences table: {sequences_table}"
    if [ -f "{failed_table}" ]; then
        echo "Failed downloads logged: {failed_table}"
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
            Dictionary with output file paths and table information
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
        
        # Output tables
        structures_csv = os.path.join(self.output_folder, "structures.csv")
        sequences_csv = os.path.join(self.output_folder, "sequences.csv")
        failed_csv = os.path.join(self.output_folder, "failed_downloads.csv")
        compounds_csv = os.path.join(self.output_folder, "compounds.csv")

        # Define tables that will be created
        tables = {
            "structures": TableInfo(
                name="structures",
                path=structures_csv,
                columns=["id", "pdb_id", "file_path", "format", "file_size", "source"],
                description="Successfully fetched structure files",
                count=len(self.pdb_ids)  # May be fewer if some downloads fail
            ),
            "sequences": TableInfo(
                name="sequences",
                path=sequences_csv,
                columns=["id", "sequence"],
                description="Protein sequences extracted from structures",
                count=len(self.pdb_ids)  # May be fewer if some downloads fail
            ),
            "compounds": TableInfo(
                name="compounds",
                path=compounds_csv,
                columns=["id", "code", "format", "smiles", "ccd"],
                description="Ligands extracted from PDB structures (SMILES from RCSB)",
                count="variable"
            ),
            "failed": TableInfo(
                name="failed",
                path=failed_csv,
                columns=["pdb_id", "error_message", "source", "attempted_path"],
                description="Failed structure fetches with error details",
                count="variable"
            )
        }

        # Get predicted compound IDs if available (set during configure_inputs)
        compound_ids = getattr(self, 'predicted_compound_ids', [])

        return {
            "structures": structure_files,
            "structure_ids": structure_ids,
            "compounds": [compounds_csv],
            "compound_ids": compound_ids,  # Populated from RCSB API query at pipeline runtime
            "sequences": [sequences_csv],
            "sequence_ids": sequence_ids,
            "tables": tables,
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
                "local_folder": self.local_folder,
                "biological_assembly": self.biological_assembly,
                "remove_waters": self.remove_waters
            }
        })
        return base_dict