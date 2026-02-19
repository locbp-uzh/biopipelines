# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
PDB tool for fetching protein structures from local folders or RCSB PDB.

Fetches structures with priority-based lookup: local_folder -> PDBs/ -> RCSB download.
Downloads are saved both in PDBs/ folder and tool output folder for reuse.
"""

import os
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


class PDBOperation:
    """Base class for PDB operations applied after loading structures."""

    def __init__(self, op_type: str, **kwargs):
        self.op_type = op_type
        self.params = kwargs

    def to_dict(self) -> Dict[str, Any]:
        """Convert operation to dictionary for serialization."""
        return {"op": self.op_type, **self.params}


class PDB(BaseConfig):
    """
    Pipeline tool for fetching protein structures from local folders or RCSB PDB.

    Implements priority-based lookup: checks local_folder (if provided), then PDBs/
    folder, then downloads from RCSB. Downloaded structures are saved to both PDBs/
    folder (for reuse) and tool output folder.

    Supports operations that are applied to structures after loading, similar to
    PyMOL and Plot tools. Operations are passed as positional arguments.

    Example:
        # Simple fetch
        pdb = PDB(pdbs="4ufc")

        # Fetch with ligand renaming (for RFdiffusion3 compatibility)
        pdb = PDB(
            pdbs="rifampicin.pdb",
            PDB.Rename("LIG", ":L:")
        )

        # Multiple operations
        pdb = PDB(
            pdbs="structure.pdb",
            PDB.Rename("LIG", ":L:"),
            PDB.Rename("HOH", ":W:")
        )
    """

    TOOL_NAME = "PDB"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== PDB ==="
echo "Uses biopipelines environment (no additional installation needed)."
echo "=== PDB ready ==="
"""

    # Lazy path descriptors
    structures_csv = Path(lambda self: os.path.join(self.output_folder, "structures.csv"))
    sequences_csv = Path(lambda self: os.path.join(self.output_folder, "sequences.csv"))
    compounds_csv = Path(lambda self: os.path.join(self.output_folder, "compounds.csv"))
    failed_csv = Path(lambda self: os.path.join(self.output_folder, "failed_downloads.csv"))
    config_file = Path(lambda self: os.path.join(self.output_folder, "fetch_config.json"))
    pdb_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_pdb.py"))

    # --- Static methods for creating operations ---

    @staticmethod
    def Rename(old: str, new: str) -> PDBOperation:
        """
        Rename a residue/ligand in the structure.

        Useful for renaming CCD ligand codes to non-CCD names for compatibility
        with tools like RFdiffusion3 that have issues with certain CCD codes.

        Args:
            old: Current residue name (e.g., "LIG", "ATP")
            new: New residue name (e.g., ":L:", "ATP1")

        Returns:
            PDBOperation for renaming

        Example:
            PDB(pdbs="structure.pdb", PDB.Rename("LIG", ":L:"))
        """
        return PDBOperation("rename", old=old, new=new)

    # --- Instance methods ---

    def __init__(self,
                 pdbs: Union[str, List[str], 'StandardizedOutput', 'DataStream'],
                 *args,
                 ids: Optional[Union[str, List[str]]] = None,
                 format: str = "pdb",
                 local_folder: Optional[str] = None,
                 biological_assembly: bool = False,
                 remove_waters: bool = True,
                 **kwargs):
        """
        Initialize PDB tool.

        Args:
            pdbs: PDB ID(s) to fetch, a folder path containing PDB files,
                  or a StandardizedOutput/DataStream from an upstream tool.
                  Can be single string, list of strings (e.g. "4ufc" or ["4ufc","1abc"]),
                  a folder path (absolute or relative to PDBs folder),
                  or a tool output whose structures will be used at execution time.
            *args: Operations to apply after loading (e.g., PDB.Rename("LIG", ":L:"))
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

            # Load all PDB files from a folder (absolute path)
            pdb = pipeline.add(PDB(
                pdbs="/path/to/my/structures"
            ))

            # Load all PDB files from a folder (relative to PDBs folder)
            pdb = pipeline.add(PDB(
                pdbs="my_structures_subfolder"
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

            # Rename ligand for RFdiffusion3 compatibility
            pdb = pipeline.add(PDB(
                pdbs="structure.pdb",
                PDB.Rename("LIG", ":L:")
            ))

            # Pass upstream tool output (files resolved at execution time)
            pdb = PDB(boltz_output, PDB.Rename("LIG", ":L:"))
        """
        # Extract operations from args
        self.operations = []
        for arg in args:
            if isinstance(arg, PDBOperation):
                self.operations.append(arg)
            else:
                raise ValueError(f"Unexpected positional argument: {arg}. Expected PDBOperation (e.g., PDB.Rename(...))")

        # Handle StandardizedOutput / DataStream input from upstream tools
        self.from_upstream = False
        if isinstance(pdbs, StandardizedOutput):
            self.structures_stream = pdbs.streams.structures
            self.from_upstream = True
        elif isinstance(pdbs, DataStream):
            self.structures_stream = pdbs
            self.from_upstream = True

        if self.from_upstream:
            self.pdb_ids = list(self.structures_stream.ids)
            self.format = self.structures_stream.format if self.structures_stream.format in ("pdb", "cif") else format.lower()
        else:
            # Check if pdbs is a folder path and load all files from it
            if isinstance(pdbs, str) and self._is_folder_path(pdbs):
                self.pdb_ids = self._load_files_from_folder(pdbs, format)
            # Normalize pdbs to list - preserve original case for local file lookups
            elif isinstance(pdbs, str):
                self.pdb_ids = [pdbs]
            elif isinstance(pdbs, list):
                self.pdb_ids = list(pdbs)
            else:
                raise ValueError(f"pdbs must be a string, list of strings, StandardizedOutput, or DataStream, got {type(pdbs)}")
            self.format = format.lower()

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

        if not self.from_upstream:
            self.local_folder = local_folder
        else:
            self.local_folder = None
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
    
    def _is_folder_path(self, path: str) -> bool:
        """
        Check if a string is a folder path.

        Args:
            path: Path to check

        Returns:
            True if path is a folder, False otherwise
        """
        # Check if absolute path exists and is a directory
        if os.path.isabs(path) and os.path.isdir(path):
            return True

        # Check if it's a relative path to PDBs folder (we'll check this at configure_inputs time)
        # For now, if it contains path separators or exists as a directory, treat it as a folder
        if os.path.isdir(path):
            return True

        return False

    def _load_files_from_folder(self, folder_path: str, format: str) -> List[str]:
        """
        Load all PDB/CIF files from a folder.

        Args:
            folder_path: Path to folder (absolute or relative to current directory)
            format: File format to look for ('pdb' or 'cif')

        Returns:
            List of file basenames without extension

        Note:
            This method is called at __init__ time. It tries to find the folder as an absolute path
            or relative to current directory. The folder will also be checked relative to PDBs folder
            later in configure_inputs if needed.
        """
        extension = f".{format}"

        # Try absolute path first
        if os.path.isabs(folder_path) and os.path.isdir(folder_path):
            target_folder = folder_path
        # Try relative path as-is
        elif os.path.isdir(folder_path):
            target_folder = folder_path
        else:
            # Store for later resolution relative to PDBs folder
            # This will be checked in configure_inputs when we have access to pipeline_folders
            self.folder_source = folder_path
            self.folder_needs_resolution = True
            return []  # Will be populated in configure_inputs

        # List all files with the specified extension
        files = [f for f in os.listdir(target_folder) if f.endswith(extension)]

        if not files:
            raise ValueError(f"No {format.upper()} files found in folder '{folder_path}'")

        # Extract basenames without extension
        basenames = [os.path.splitext(f)[0] for f in sorted(files)]

        print(f"  Found {len(basenames)} {format.upper()} files in folder: {folder_path}")

        # Store the folder path for use in configure_inputs
        self.folder_source = target_folder
        self.folder_needs_resolution = False

        return basenames

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

        # When input comes from an upstream tool, files will exist at execution time only
        if self.from_upstream:
            self.found_locally = []
            self.needs_download = []
            print(f"  PDB: using {len(self.pdb_ids)} structures from upstream tool (validated at execution time)")
            return

        # Check if folder needs resolution relative to PDBs
        if hasattr(self, 'folder_needs_resolution') and self.folder_needs_resolution:
            repo_pdbs_folder = pipeline_folders.get('PDBs', '')
            candidate_path = os.path.join(repo_pdbs_folder, self.folder_source)

            if os.path.isdir(candidate_path):
                # Found the folder relative to PDBs
                extension = f".{self.format}"
                files = [f for f in os.listdir(candidate_path) if f.endswith(extension)]

                if not files:
                    raise ValueError(f"No {self.format.upper()} files found in folder '{candidate_path}'")

                # Extract basenames without extension
                basenames = [os.path.splitext(f)[0] for f in sorted(files)]

                print(f"  Found {len(basenames)} {self.format.upper()} files in folder (relative to PDBs): {self.folder_source}")

                # Update pdb_ids and custom_ids
                self.pdb_ids = basenames
                self.custom_ids = basenames.copy()
                self.folder_source = candidate_path
                self.folder_needs_resolution = False
            else:
                raise ValueError(f"Folder '{self.folder_source}' not found (tried absolute, relative, and relative to PDBs folder)")

        # Check which files exist locally and which will need to be downloaded
        repo_pdbs_folder = pipeline_folders.get('PDBs', '')
        self.found_locally = []
        self.needs_download = []

        for pdb_id, custom_id in zip(self.pdb_ids, self.custom_ids):
            extension = ".pdb" if self.format == "pdb" else ".cif"
            found = False
            local_path = None

            # Check folder_source first if it was set (from folder loading)
            if hasattr(self, 'folder_source') and self.folder_source:
                local_path = os.path.join(self.folder_source, f"{pdb_id}{extension}")
                if os.path.exists(local_path):
                    self.found_locally.append((pdb_id, local_path))
                    found = True

            # Check local_folder if specified
            if not found and self.local_folder:
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
                        format_label = self.format.upper()
                        if has_ligands:
                            print(f"  Found {format_label} {pdb_id} locally: {local_path} (contains ligands)")
                        else:
                            print(f"  Found {format_label} {pdb_id} locally: {local_path}")
                    except:
                        # RCSB query failed, just show found locally
                        print(f"  Found {self.format.upper()} {pdb_id} locally: {local_path}")
                else:
                    # Custom file, not an RCSB ID
                    print(f"  Found {self.format.upper()} {pdb_id} locally: {local_path}")
            else:
                # Not found locally - check if valid RCSB ID
                rcsb_id = pdb_id.upper()
                if len(rcsb_id) != 4 or not rcsb_id.isalnum():
                    raise ValueError(f"PDB '{pdb_id}' not found locally and is not a valid RCSB PDB ID (must be 4 alphanumeric characters)")

                # Check if exists on RCSB
                has_ligands = self._check_rcsb_exists(rcsb_id, custom_id)
                self.needs_download.append(pdb_id)

                if has_ligands:
                    print(f"  {self.format.upper()} {pdb_id} not found locally, will download from RCSB (contains ligands)")
                else:
                    print(f"  {self.format.upper()} {pdb_id} not found locally, will download from RCSB")

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
                    # Filter out common solvents/ions/crystallization agents
                    common_solvents = {
                        'HOH', 'WAT', 'H2O',  # Water
                        'NA', 'CL', 'CA', 'MG', 'K', 'ZN', 'MN', 'FE', 'CU', 'NI', 'CO',  # Common ions
                        'SO4', 'PO4', 'NO3',  # Anions
                        'GOL', 'EDO', 'PEG', 'PGE', 'PE4', 'PE3', 'P6G', 'PG4', '1PE',  # Glycols and PEGs
                        'ACT', 'ACE', 'ACY',  # Acetate
                        'PYR', 'PYO',  # Pyruvate
                        'DMS', 'DMSO', 'BME', 'MPD', 'TRS', 'EPE'  # Common solvents
                    }
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
        """Fetch ligand information from RCSB API at configuration time (silently for local files)."""
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

        # Add operations if any
        if self.operations:
            op_summaries = []
            for op in self.operations:
                if op.op_type == "rename":
                    op_summaries.append(f"Rename({op.params['old']} → {op.params['new']})")
                else:
                    op_summaries.append(op.op_type)
            config_lines.append(f"OPERATIONS: {', '.join(op_summaries)}")

        # Add status of files found/not found
        if hasattr(self, 'found_locally') and self.found_locally:
            config_lines.append(f"FOUND_LOCALLY: {', '.join([pdb_id for pdb_id, _ in self.found_locally])}")
        if hasattr(self, 'needs_download') and self.needs_download:
            config_lines.append(f"NEEDS_DOWNLOAD: {', '.join(self.needs_download)}")

        return config_lines
    
    def generate_script(self, script_path: str) -> str:
        """Generate PDB execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        script_content = "#!/bin/bash\n"
        script_content += "# PDB execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self.generate_script_run_pdb()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_pdb(self) -> str:
        """Generate the PDB fetching execution part of the script."""
        import json

        # Get PDBs folder path from folder manager
        repo_pdbs_folder = self.folders['PDBs']

        # If folder_source is set (from folder loading), use it as local_folder for the runtime script
        effective_local_folder = getattr(self, 'folder_source', None) or self.local_folder

        config_data = {
            "pdb_ids": self.pdb_ids,
            "custom_ids": self.custom_ids,
            "format": self.format,
            "local_folder": effective_local_folder,
            "repo_pdbs_folder": repo_pdbs_folder,
            "biological_assembly": self.biological_assembly,
            "remove_waters": self.remove_waters,
            "output_folder": self.output_folder,
            "structures_table": self.structures_csv,
            "sequences_table": self.sequences_csv,
            "failed_table": self.failed_csv,
            "compounds_table": self.compounds_csv,
            "operations": [op.to_dict() for op in self.operations]
        }

        # When input comes from upstream, provide the source files for the runtime script
        if self.from_upstream:
            config_data["from_upstream"] = True
            config_data["upstream_files"] = list(self.structures_stream.files)
            config_data["upstream_files_contain_wildcards"] = self.structures_stream.files_contain_wildcards
            config_data["upstream_map_table"] = self.structures_stream.map_table

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        return f"""echo "Fetching {len(self.pdb_ids)} structures"
echo "Format: {self.format.upper()}"
echo "PDB IDs: {', '.join(self.pdb_ids)}"
echo "Custom IDs: {', '.join(self.custom_ids)}"
echo "Priority: {'local_folder -> ' if self.local_folder else ''}PDBs/ -> RCSB download"
echo "Output folder: {self.output_folder}"

python "{self.pdb_py}" --config "{self.config_file}"

"""
    
    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after structure fetching."""
        # Generate structure file paths using custom IDs
        extension = ".pdb" if self.format == "pdb" else ".cif"
        structure_files = [os.path.join(self.output_folder, f"{custom_id}{extension}")
                          for custom_id in self.custom_ids]

        tables = {
            "structures": TableInfo(
                name="structures",
                path=self.structures_csv,
                columns=["id", "pdb_id", "file_path", "format", "file_size", "source"],
                description="Successfully fetched structure files",
                count=len(self.pdb_ids)
            ),
            "sequences": TableInfo(
                name="sequences",
                path=self.sequences_csv,
                columns=["id", "sequence"],
                description="Protein sequences extracted from structures",
                count=len(self.pdb_ids)
            ),
            "compounds": TableInfo(
                name="compounds",
                path=self.compounds_csv,
                columns=["id", "code", "format", "smiles", "ccd"],
                description="Ligands extracted from PDB structures (SMILES from RCSB)",
                count="variable"
            ),
            "failed": TableInfo(
                name="failed",
                path=self.failed_csv,
                columns=["pdb_id", "error_message", "source", "attempted_path"],
                description="Failed structure fetches with error details",
                count="variable"
            )
        }

        # Create DataStreams
        structures = DataStream(
            name="structures",
            ids=self.custom_ids.copy(),
            files=structure_files,
            map_table=self.structures_csv,
            format=self.format
        )

        sequences = DataStream(
            name="sequences",
            ids=self.custom_ids.copy(),
            files=[self.sequences_csv],
            map_table=self.sequences_csv,
            format="csv"
        )

        # Get predicted compound IDs if available (set during configure_inputs)
        compound_ids = getattr(self, 'predicted_compound_ids', [])
        compounds = DataStream(
            name="compounds",
            ids=compound_ids,
            files=[],  # Value-based format - data is in map_table, not individual files
            map_table=self.compounds_csv,
            format="csv"
        )

        return {
            "structures": structures,
            "sequences": sequences,
            "compounds": compounds,
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
                "remove_waters": self.remove_waters,
                "operations": [op.to_dict() for op in self.operations]
            }
        })
        return base_dict