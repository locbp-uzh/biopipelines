# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
PDB tool for fetching protein structures from local folders or RCSB PDB.

Fetches structures with priority-based lookup: local_folder -> pdbs/ -> RCSB download.
Downloads are saved both in pdbs/ folder and tool output folder for reuse.
"""

import os
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .biopipelines_io import TableReference
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from biopipelines_io import TableReference


def _normalize_selection(selection):
    """Normalize a remove/operation selection to a str or TableReference.

    Accepts the documented column-reference forms and converts them to a single
    canonical TableReference, which serializes to the TABLE_REFERENCE:path:column
    string the pipe script resolves per-structure:

    - ``str``                       — a literal PyMOL-style selection, passed through.
    - ``TableReference``            — passed through (e.g. ``tbl.tables.x.col``).
    - ``(TableInfo, "col")`` tuple  — the user-manual column-reference form.
    - ``(path_str, "col")`` tuple   — a raw (path, column) pair.
    """
    if isinstance(selection, (str, TableReference)):
        return selection
    if isinstance(selection, tuple) and len(selection) == 2:
        table, column = selection
        if not isinstance(column, str):
            raise ValueError(
                f"remove selection tuple must be (TableInfo|path, column_name), "
                f"got column of type {type(column).__name__}"
            )
        if isinstance(table, TableInfo):
            # TableInfo routes bare attribute access to TableReference; the real
            # path is only on .info.
            return TableReference(table.info.path, column)
        if isinstance(table, str):
            return TableReference(table, column)
        raise ValueError(
            f"remove selection tuple's first element must be a TableInfo or path "
            f"string, got {type(table).__name__}"
        )
    raise ValueError(
        f"remove selection must be a PyMOL-style string, a TableReference "
        f"(e.g. tool.tables.structures.designed), or a (TableInfo, column) "
        f"tuple; got {type(selection).__name__}"
    )


class PDBOperation:
    """Base class for PDB operations applied after loading structures."""

    def __init__(self, op_type: str, **kwargs):
        self.op_type = op_type
        self.params = kwargs

    def to_dict(self) -> Dict[str, Any]:
        """Convert operation to dictionary for serialization."""
        # Stringify any non-JSON-native params (e.g. a TableReference selection
        # serializes to its TABLE_REFERENCE:path:column form, parsed downstream).
        params = {
            k: (v if isinstance(v, (str, int, float, bool, type(None))) else str(v))
            for k, v in self.params.items()
        }
        return {"op": self.op_type, **params}


class PDB(BaseConfig):
    """
    Pipeline tool for fetching protein structures from local folders or RCSB PDB.

    Implements priority-based lookup: checks local_folder (if provided), then pdbs/
    folder, then downloads from RCSB. Downloaded structures are saved to both pdbs/
    folder (for reuse) and tool output folder.

    Supports operations that are applied to structures after loading, similar to
    PyMOL and Plot tools. Operations are passed as positional arguments.

    Example:
        # Simple fetch
        pdb = PDB(pdbs="4ufc")

        # Fetch with ligand renaming (for RFdiffusion3 compatibility)
        pdb = PDB(
            pdbs="rifampicin.pdb",
            PDB.rename("LIG", ":L:")
        )

        # Multiple operations
        pdb = PDB(
            pdbs="structure.pdb",
            PDB.rename("LIG", ":L:"),
            PDB.rename("HOH", ":W:")
        )
    """

    TOOL_NAME = "PDB"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== PDB ==="
echo "Uses biopipelines environment (no additional installation needed)."
touch "$INSTALL_SUCCESS"
echo "=== PDB ready ==="
"""

    # Lazy path descriptors — all three declared streams (structures,
    # sequences, compounds) are content-bearing: their map_table CSV IS
    # the content. Each lives inside its own stream folder; the "failed"
    # table is a standalone TableInfo under tables/.
    structures_csv = Path(lambda self: self.stream_path("structures", "structures.csv"))
    sequences_csv = Path(lambda self: self.stream_path("sequences", "sequences.csv"))
    compounds_csv = Path(lambda self: self.stream_path("compounds", "compounds.csv"))
    failed_csv = Path(lambda self: self.table_path("failed"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    config_file = Path(lambda self: self.configuration_path("fetch_config.json"))
    pdb_py = Path(lambda self: self.pipe_script_path("pipe_pdb.py"))

    # --- Static methods for creating operations ---

    @staticmethod
    def rename(old: str, new: str) -> PDBOperation:
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
            PDB(pdbs="structure.pdb", PDB.rename("LIG", ":L:"))
        """
        return PDBOperation("rename", old=old, new=new)

    @staticmethod
    def remove(selection: Union[str, tuple, 'TableReference'], remove_hetatm: bool = True) -> PDBOperation:
        """
        Remove residues from the structure by selection.

        Useful for truncating a structure after docking/posing — e.g. dock a
        ligand against the full-length protein (which has the real pocket), then
        drop an N-terminal segment, carrying the bound ligand (HETATM) along.

        Args:
            selection: Residues to remove, as a PyMOL-style selection string or a
                table column reference resolved per-structure at runtime.
                - `"1-83"`              — residues 1-83 (any chain)
                - `"1-83+90-95"`        — multiple ranges
                - `"A1-83"`             — chain-qualified range
                - `tool.tables.x.col`   — a TableReference column (per-structure)
                - `(TableInfo, "col")`  — equivalent (TableInfo, column) tuple
                A column reference is matched to each structure by ID at runtime
                (a single-row table broadcasts to all structures; an unmatched
                structure is left unchanged).
            remove_hetatm: When True (default), a HETATM (ligand/ion) whose residue
                number falls inside the selected range is removed along with the
                protein residues. Set False to remove only protein ATOM residues and
                keep all HETATM records — useful when a ligand shares a residue
                number with a removed protein span but should be preserved. A ligand
                sitting *outside* the removed range is kept either way; the cleanest
                way to spare a ligand is to chain-qualify the selection (e.g.
                "A1-83") so only the protein chain is touched.

        Returns:
            PDBOperation for residue removal

        Example:
            # Truncate to residues >=84; a ligand outside 1-83 is kept
            PDB(docked, PDB.remove("1-83"))

            # Drop chain A residues 1-83, leaving a ligand on another chain intact
            PDB(docked, PDB.remove("A1-83"))

            # Per-structure truncation from an upstream table column
            PDB(docked, PDB.remove(rfd.tables.structures.flank))
        """
        return PDBOperation("remove", selection=_normalize_selection(selection),
                            remove_hetatm=remove_hetatm)

    @staticmethod
    def break_bond(atom1: str, atom2: str) -> PDBOperation:
        """
        Break a bond between two atoms (PyMOL `unbond`-style).

        Removes any CONECT/LINK records joining the two named atoms; coordinates
        are left untouched. Use it to sever a covalent attachment after it has
        served its purpose — e.g. after PLACER resamples a dye covalently tethered
        to a catalytic cysteine, break the bond so downstream design tools
        (RFdiffusion) see the ligand as a separate non-covalent HETATM.

        Atom selections use the standard BioPipelines `<residue>.<atom>` syntax
        (same as Distance/Angle):
            - `"A145.SG"`  — atom SG of residue 145 on chain A
            - `"LIG.C12"`  — atom C12 of residue LIG
            - `"145.SG"`   — chain-agnostic residue number

        Args:
            atom1: First atom of the bond (e.g. `"A145.SG"`).
            atom2: Second atom of the bond (e.g. `"LIG.C12"`).

        Returns:
            PDBOperation for bond removal

        Example:
            # Sever the BG-Cys145 thioether on all PLACER conformers before design
            PDB(placer_poses, PDB.break_bond("A145.SG", "LIG.C12"))
        """
        return PDBOperation("break_bond", atom1=atom1, atom2=atom2)

    @staticmethod
    def rotate_bond(atom1: str, atom2: str, angle: float) -> PDBOperation:
        """
        Rotate a fragment about the ``atom1``–``atom2`` bond (torsion rotation).

        Everything on ``atom2``'s side of the bond is rotated by ``angle`` degrees
        about the bond axis; ``atom1``'s side stays fixed. The moving fragment is
        found by intra-residue connectivity (atoms reachable from ``atom2`` without
        crossing the ``atom1``–``atom2`` bond), so only the chosen rotatable bond's
        downstream atoms move. Useful for re-orienting a flexible ligand moiety
        while keeping a covalent anchor in place — e.g. swing a dye's fluorophore
        toward a target region by rotating ~180° about a linker bond, keeping the
        benzyl/Cys-attached end fixed.

        Atom selections use the standard `<residue>.<atom>` syntax (as break_bond).

        Args:
            atom1: Fixed end of the bond axis (e.g. `"LIG.C60"`).
            atom2: Moving end; its side of the bond rotates (e.g. `"LIG.C47"`).
            angle: Rotation in degrees (e.g. `180`).

        Returns:
            PDBOperation for the bond rotation.

        Example:
            # Swing the fluorophore 180 deg about a linker bond, anchor end fixed
            PDB(complex, PDB.rotate_bond("LIG.C60", "LIG.C47", 180))
        """
        return PDBOperation("rotate_bond", atom1=atom1, atom2=atom2, angle=float(angle))

    # --- Instance methods ---

    def __init__(self,
                 pdbs: Union[str, List[str], Dict[str, str], 'StandardizedOutput', 'DataStream'],
                 *args,
                 ids: Optional[Union[str, List[str]]] = None,
                 convert: Optional[str] = None,
                 local_folder: Optional[str] = None,
                 biological_assembly: bool = False,
                 remove_waters: bool = True,
                 chain: Union[str, List[str]] = "auto",
                 split_chains: bool = False,
                 fetch_compounds: bool = True,
                 **kwargs):
        """
        Initialize PDB tool.

        Args:
            pdbs: PDB ID(s) to fetch, a folder path containing PDB files,
                  or a StandardizedOutput/DataStream from an upstream tool.
                  Can be single string, list of strings (e.g. "4ufc" or ["4ufc","1abc"]),
                  a dictionary mapping IDs to PDB codes (e.g. {"POI": "4ufc", "POI2": "1abc"}),
                  a folder path (absolute or relative to PDBs folder),
                  or a tool output whose structures will be used at execution time.
            *args: Operations to apply after loading (e.g., PDB.rename("LIG", ":L:"))
            ids: Custom IDs for renaming. Can be single string or list of strings (e.g. "POI" or ["POI1","POI2"]). If None, uses pdbs as ids.
                 Ignored when pdbs is a dictionary (ids come from dict keys).
            convert: Target format to convert structures to - "pdb", "cif", or None (default).
                     When None, no conversion is performed: structures are kept in whatever
                     format they are found locally or downloaded as from RCSB.
                     The structures DataStream format will be "pdb|cif" when None.
            local_folder: Custom local folder to check first (before pdbs/). Default: None
            biological_assembly: Whether to download biological assembly from RCSB (default: False)
            remove_waters: Whether to remove water molecules from structures (default: True)
            chain: Which chain(s) to extract the sequence from. Note that
                only an explicit single chain letter filters the structure
                file on disk; every other form leaves all chains intact in
                the .pdb/.cif and only affects the sequences stream.
                * "auto" (default) — structure file unchanged (all chains
                  kept on disk). Sequences stream emits one <id>,<sequence>
                  row per input, carrying the longest chain.
                * "all" — structure file unchanged. Sequences stream emits
                  one <id>_<chain_letter>,<sequence> row per chain present
                  in the structure (sequence pulled from RCSB FASTA when
                  the input is a PDB code, otherwise from the structure
                  file). No aggregate longest row.
                * List[str], e.g. ["A","C"] — structure file unchanged.
                  Sequences stream emits one <id>_<chain_letter>,<sequence>
                  row per listed chain. Cardinality is known at config
                  time (literal IDs).
                * Single chain letter "A"/"B"/... — structure file on
                  disk is filtered down to just that chain. Sequences
                  stream emits a single <id>,<sequence_for_that_chain>
                  row.
            split_chains: When True (and chain is "all" or a list), also split
                the structure file into one .pdb (or .cif) per chain letter at
                <output>/<custom_id>_<letter>.<ext>. The structures stream
                IDs become <custom_id>_<letter>. Mutually exclusive with the
                single-chain forms ("auto" / explicit chain letter).
            **kwargs: Additional parameters

        Fetch Priority:
            For each PDB ID, searches in order:
            1. local_folder (if parameter provided)
            2. ./pdbs/ folder in repository
            3. Download from RCSB PDB (saved to both pdbs/ and output folder)

        Output:
            Streams: structures (.pdb/.cif), sequences (.csv), compounds (.csv)
            Tables:
                structures: id | pdb_id | file_path | format | file_size | source
                sequences: id | sequence
                compounds: id | code | format | smiles | ccd
                failed: pdb_id | error_message | source | attempted_path
        """
        # Extract operations from args
        self.operations = []
        for arg in args:
            if isinstance(arg, PDBOperation):
                self.operations.append(arg)
            else:
                raise ValueError(f"Unexpected positional argument: {arg}. Expected PDBOperation (e.g., PDB.rename(...))")

        # Collected at RCSB-lookup time; initialize up front so any code path
        # that appends (e.g. _check_ligands_in_rcsb) cannot hit AttributeError.
        self.predicted_compound_ids = []

        # Dict input: {id: pdb_code} -> extract ids from keys, pdb codes from values
        if isinstance(pdbs, dict):
            if ids is not None:
                print("  Warning: 'ids' parameter ignored when pdbs is a dictionary (using dict keys)")
            if not pdbs:
                raise ValueError("Must provide at least one PDB ID")
            ids = list(pdbs.keys())
            pdbs = list(pdbs.values())

        # Handle StandardizedOutput / DataStream input from upstream tools
        self.from_upstream = False
        self.pdbs_input = pdbs  # kept for upstream missing-table detection
        if isinstance(pdbs, StandardizedOutput):
            self.structures_stream = pdbs.streams.structures
            self.from_upstream = True
        elif isinstance(pdbs, DataStream):
            self.structures_stream = pdbs
            self.from_upstream = True

        if self.from_upstream:
            self.pdb_ids = list(self.structures_stream.ids)
            upstream_fmt = self.structures_stream.format
            self.convert = upstream_fmt if upstream_fmt in ("pdb", "cif") else convert.lower() if convert else None
        else:
            # Check if pdbs is a folder path and load all files from it
            if isinstance(pdbs, str) and self._is_folder_path(pdbs):
                self.pdb_ids = self._load_files_from_folder(pdbs)
            # Normalize pdbs to list - preserve original case for local file lookups
            elif isinstance(pdbs, str):
                self.pdb_ids = [pdbs]
            elif isinstance(pdbs, list):
                self.pdb_ids = list(pdbs)
            else:
                raise ValueError(f"pdbs must be a string, list of strings, StandardizedOutput, or DataStream, got {type(pdbs)}")
            self.convert = convert.lower() if convert else None

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
        # Chain may be one of:
        #   - "auto" (default longest-only)
        #   - "all"  (one row per chain in the structure)
        #   - "A"/"B"/... (single explicit chain letter)
        #   - List[str] of explicit chain letters (subset of "all")
        if isinstance(chain, list):
            if not chain:
                raise ValueError("chain list cannot be empty")
            if not all(isinstance(c, str) and c for c in chain):
                raise ValueError("chain list must contain non-empty strings")
            self.chain = list(chain)
        else:
            self.chain = chain
        self.split_chains = split_chains
        self.fetch_compounds = fetch_compounds

        chain_is_multi = isinstance(self.chain, list) or self.chain == "all"
        if self.split_chains and not chain_is_multi:
            raise ValueError("split_chains=True requires chain=\"all\" or chain=[...]")

        # Validate convert
        if self.convert is not None and self.convert not in ["pdb", "cif"]:
            raise ValueError(f"Invalid convert: {self.convert}. Must be 'pdb', 'cif', or None")

        if isinstance(self.chain, list):
            for c in self.chain:
                _validate_freeform_string("chain", c)
        else:
            _validate_freeform_string("chain", self.chain)
        _validate_freeform_string("local_folder", self.local_folder)

        # Note: PDB ID format validation is skipped at init time because:
        # 1. local_folder may contain custom-named files
        # 2. pdbs/ folder may contain custom-named files
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

    def _load_files_from_folder(self, folder_path: str) -> List[str]:
        """
        Load all PDB/CIF files from a folder.

        Args:
            folder_path: Path to folder (absolute or relative to current directory)

        Returns:
            List of file basenames without extension

        Note:
            This method is called at __init__ time. It tries to find the folder as an absolute path
            or relative to current directory. The folder will also be checked relative to PDBs folder
            later in configure_inputs if needed.
        """
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

        # List all PDB and CIF files
        files = [f for f in os.listdir(target_folder) if f.endswith(".pdb") or f.endswith(".cif")]

        if not files:
            raise ValueError(f"No PDB or CIF files found in folder '{folder_path}'")

        # Extract basenames without extension
        basenames = [os.path.splitext(f)[0] for f in sorted(files)]

        print(f"  Found {len(basenames)} structure file(s) in folder: {folder_path}")

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

        if self.convert is not None and self.convert not in ["pdb", "cif"]:
            raise ValueError("convert must be 'pdb', 'cif', or None")
    
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
            repo_pdbs_folder = pipeline_folders.get('pdbs', '')
            candidate_path = os.path.join(repo_pdbs_folder, self.folder_source)

            if os.path.isdir(candidate_path):
                # Found the folder relative to PDBs; scan for both pdb and cif files
                files = [f for f in os.listdir(candidate_path) if f.endswith(".pdb") or f.endswith(".cif")]

                if not files:
                    raise ValueError(f"No PDB or CIF files found in folder '{candidate_path}'")

                # Extract basenames without extension
                basenames = [os.path.splitext(f)[0] for f in sorted(files)]

                print(f"  Found {len(basenames)} structure file(s) in folder (relative to PDBs): {self.folder_source}")

                # Update pdb_ids and custom_ids
                self.pdb_ids = basenames
                self.custom_ids = basenames.copy()
                self.folder_source = candidate_path
                self.folder_needs_resolution = False
            else:
                raise ValueError(f"Folder '{self.folder_source}' not found (tried absolute, relative, and relative to PDBs folder)")

        # Check which files exist locally and which will need to be downloaded
        repo_pdbs_folder = pipeline_folders.get('pdbs', '')
        self.found_locally = []
        self.needs_download = []
        # Track which formats are present locally (only relevant when convert=None)
        self.local_formats = {"pdb": False, "cif": False}

        for pdb_id, custom_id in zip(self.pdb_ids, self.custom_ids):
            found = False
            local_path = None

            # Check if pdb_id is already an absolute path to an existing file
            if os.path.isabs(pdb_id) and os.path.isfile(pdb_id):
                self.found_locally.append((pdb_id, pdb_id))
                local_path = pdb_id
                found = True

            if not found:
                if self.convert is not None:
                    # Specific convert target: prefer that extension locally (may convert if only other found)
                    extensions_to_try = [".pdb" if self.convert == "pdb" else ".cif", ".cif" if self.convert == "pdb" else ".pdb"]
                else:
                    # No conversion: accept both, pdb first then cif
                    extensions_to_try = [".pdb", ".cif"]

                for extension in extensions_to_try:
                    # Check folder_source first if it was set (from folder loading)
                    if not found and hasattr(self, 'folder_source') and self.folder_source:
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

                    # Check pdbs/ folder
                    if not found and repo_pdbs_folder:
                        local_path = os.path.join(repo_pdbs_folder, f"{pdb_id}{extension}")
                        if os.path.exists(local_path):
                            self.found_locally.append((pdb_id, local_path))
                            found = True

                    if found:
                        break

            if found and self.convert is None and local_path:
                fmt = "cif" if local_path.endswith(".cif") else "pdb"
                self.local_formats[fmt] = True

            if found:
                # Check if it's a valid RCSB ID and query for ligands
                rcsb_id = pdb_id.upper()
                if len(rcsb_id) == 4 and rcsb_id.isalnum():
                    # Valid RCSB format - check for ligands
                    try:
                        has_ligands = self._check_rcsb_exists_silent(rcsb_id, custom_id)
                        if has_ligands:
                            print(f"  Found {pdb_id} locally: {local_path} (contains ligands)")
                        else:
                            print(f"  Found {pdb_id} locally: {local_path}")
                    except:
                        # RCSB query failed, just show found locally
                        print(f"  Found {pdb_id} locally: {local_path}")
                else:
                    # Custom file, not an RCSB ID
                    print(f"  Found {pdb_id} locally: {local_path}")
            else:
                # Not found locally - check if valid RCSB ID.
                # wwPDB Format 3.3: 4 chars, first is a digit, rest alphanumeric.
                rcsb_id = pdb_id.upper()
                if len(rcsb_id) != 4 or not rcsb_id.isalnum() or not rcsb_id[0].isdigit():
                    raise ValueError(f"PDB '{pdb_id}' not found locally and is not a valid RCSB PDB ID (must be 4 characters, first a digit, rest alphanumeric)")

                # Check if exists on RCSB
                has_ligands = self._check_rcsb_exists(rcsb_id, custom_id)
                self.needs_download.append(pdb_id)

                if has_ligands:
                    print(f"  {pdb_id} not found locally, will download from RCSB (contains ligands)")
                else:
                    print(f"  {pdb_id} not found locally, will download from RCSB")

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
        except (requests.RequestException, ValueError):
            # Network failure, HTTP error, or non-JSON response — treat as "no
            # ligand info available" rather than crashing pipeline construction.
            return False, []

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

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        convert_display = f"convert to {self.convert.upper()}" if self.convert else "keep as-is (pdb|cif)"
        config_lines.extend([
            f"PDB_IDS: {', '.join(self.pdb_ids)} ({len(self.pdb_ids)} structures)",
            f"CUSTOM_IDS: {', '.join(self.custom_ids)}",
            f"CONVERT: {convert_display}",
            f"LOCAL_FOLDER: {self.local_folder if self.local_folder else 'None (uses pdbs/)'}",
            f"BIOLOGICAL_ASSEMBLY: {self.biological_assembly}",
            f"REMOVE_WATERS: {self.remove_waters}",
            f"CHAIN: {','.join(self.chain) if isinstance(self.chain, list) else self.chain}",
            f"SPLIT_CHAINS: {self.split_chains}"
        ])

        # Add operations if any
        if self.operations:
            op_summaries = []
            for op in self.operations:
                if op.op_type == "rename":
                    op_summaries.append(f"Rename({op.params['old']} -> {op.params['new']})")
                elif op.op_type == "remove":
                    op_summaries.append(f"remove({op.params['selection']})")
                elif op.op_type == "break_bond":
                    op_summaries.append(f"break_bond({op.params['atom1']}, {op.params['atom2']})")
                elif op.op_type == "rotate_bond":
                    op_summaries.append(f"rotate_bond({op.params['atom1']}, {op.params['atom2']}, {op.params['angle']})")
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
        # Layout sub-dirs already created by the pipeline.
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
        repo_pdbs_folder = self.folders['pdbs']

        # If folder_source is set (from folder loading), use it as local_folder for the runtime script
        effective_local_folder = getattr(self, 'folder_source', None) or self.local_folder

        # When convert is None, pipe_pdb.py keeps whatever format is found locally or downloaded.
        # New layout: PDBs land inside the structures/ stream folder;
        # each table CSV has its own canonical home.
        config_data = {
            "pdb_ids": self.pdb_ids,
            "custom_ids": self.custom_ids,
            "convert": self.convert,  # None means no conversion
            "local_folder": effective_local_folder,
            "repo_pdbs_folder": repo_pdbs_folder,
            "biological_assembly": self.biological_assembly,
            "remove_waters": self.remove_waters,
            "chain": self.chain,
            "split_chains": self.split_chains,
            "output_folder": self.stream_folder("structures"),
            "structures_table": self.structures_csv,
            "sequences_table": self.sequences_csv,
            "failed_table": self.failed_csv,
            "compounds_table": self.compounds_csv,
            "fetch_compounds": self.fetch_compounds,
            "operations": [op.to_dict() for op in self.operations]
        }

        # When input comes from upstream, provide the source files for the runtime script
        if self.from_upstream:
            config_data["from_upstream"] = True
            config_data["upstream_files"] = list(self.structures_stream.files)
            config_data["upstream_map_table"] = self.structures_stream.map_table
            # Ids the upstream already filtered out are excused (not failures).
            # The pipe script owns tables/missing.csv: it skips these ids and
            # writes them out under THIS tool's id (custom_id), so a rename via
            # ids= still excuses the renamed output at completion time.
            upstream_missing = self._collect_upstream_missing_paths(self.pdbs_input)
            if upstream_missing:
                config_data["upstream_missing"] = upstream_missing
                config_data["missing_out"] = self.missing_csv

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        convert_display = f"convert to {self.convert.upper()}" if self.convert else "keep as-is (pdb|cif)"
        return f"""echo "Fetching {len(self.pdb_ids)} structures"
echo "Convert: {convert_display}"
echo "PDB IDs: {', '.join(self.pdb_ids)}"
echo "Custom IDs: {', '.join(self.custom_ids)}"
echo "Priority: {'local_folder -> ' if self.local_folder else ''}pdbs/ -> RCSB download"
echo "Output folder: {self.output_folder}"

python "{self.pdb_py}" --config "{self.config_file}"

"""
    
    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after structure fetching."""
        structures_dir = self.stream_folder("structures")
        # When convert is set, output will always be that format; otherwise it can be either
        if self.convert is not None:
            extension = ".pdb" if self.convert == "pdb" else ".cif"
            stream_format = self.convert
        else:
            # No conversion: predict extension where possible
            local_formats = getattr(self, 'local_formats', {"pdb": False, "cif": False})
            has_downloads = bool(getattr(self, 'needs_download', []))
            only_pdb = local_formats["pdb"] and not local_formats["cif"]
            only_cif = local_formats["cif"] and not local_formats["pdb"]

            if not has_downloads and only_pdb:
                extension = ".pdb"
                stream_format = "pdb"
            elif not has_downloads and only_cif:
                extension = ".cif"
                stream_format = "cif"
            else:
                # Mixed formats or downloads present: extension unknown at config time
                extension = ".*"
                stream_format = "pdb|cif"

        # One file template, expanded per resolved id by the framework — works for
        # literal ids and lazy patterns (e.g. "_<1..10>" or "[_<chain>]") alike.
        structure_files = [os.path.join(structures_dir, f"<id>{extension}")]
        if self.split_chains and isinstance(self.chain, list):
            # Explicit chain list -> deterministic per-chain ids.
            structure_id_patterns = [f"{cid}_{ch}"
                                     for cid in self.custom_ids
                                     for ch in self.chain]
        elif self.split_chains:
            # chain="all" -> chain letters resolved at runtime; lazy id pattern.
            structure_id_patterns = [f"{cid}[_<chain>]" for cid in self.custom_ids]
        else:
            structure_id_patterns = list(self.custom_ids)

        tables = {
            "structures": TableInfo(
                name="structures",
                path=self.structures_csv,
                columns=["id", "pdb_id", "file_path", "format", "file_size", "source"],
                description="Successfully fetched structure files"
            ),
            "sequences": TableInfo(
                name="sequences",
                path=self.sequences_csv,
                columns=["id", "sequence"],
                description="Protein sequences extracted from structures"
            ),
            "compounds": TableInfo(
                name="compounds",
                path=self.compounds_csv,
                columns=["id", "code", "format", "smiles", "ccd"],
                description="Ligands extracted from PDB structures (SMILES from RCSB)"
            ),
            "failed": TableInfo(
                name="failed",
                path=self.failed_csv,
                columns=["pdb_id", "error_message", "source", "attempted_path"],
                description="Failed structure fetches with error details"
            )
        }

        # Excuse upstream-filtered ids: when the upstream input carries a
        # `missing` manifest, declare + own one here so the completion check
        # treats those ids' absent structures as expected, not failures.
        if self._collect_upstream_missing_paths(self.pdbs_input):
            tables["missing"] = self.missing_table_info(self.missing_csv)

        # Create DataStreams
        structures = DataStream(
            name="structures",
            ids=structure_id_patterns,
            files=structure_files,
            map_table=self.structures_csv,
            format=stream_format
        )

        # Sequences stream cardinality follows the chain parameter:
        #   "auto" / explicit chain letter -> exactly one row per input id
        #                                     (literal IDs).
        #   List[str]                      -> one row per (input id, chain
        #                                     letter) pair, count known at
        #                                     config time -> literal IDs.
        #   "all"                          -> one row per chain in the
        #                                     structure, count resolved at
        #                                     runtime via the RCSB FASTA
        #                                     fetch -> lazy IDs.
        if isinstance(self.chain, list):
            sequence_id_patterns = [f"{cid}_{ch}"
                                    for cid in self.custom_ids
                                    for ch in self.chain]
        elif self.chain == "all":
            sequence_id_patterns = [f"{cid}[_<chain>]" for cid in self.custom_ids]
        else:
            sequence_id_patterns = list(self.custom_ids)
        sequences = DataStream(
            name="sequences",
            ids=sequence_id_patterns,
            files=[],
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
                "convert": self.convert,
                "local_folder": self.local_folder,
                "biological_assembly": self.biological_assembly,
                "remove_waters": self.remove_waters,
                "chain": self.chain,
                "split_chains": self.split_chains,
                "operations": [op.to_dict() for op in self.operations]
            }
        })
        return base_dict