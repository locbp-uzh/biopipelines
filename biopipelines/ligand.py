# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Ligand tool for fetching small molecule ligands from RCSB PDB, PubChem, or SMILES strings.

Downloads SDF files and converts them to PDB format with proper atom numbering.
Supports lookup by CCD code (RCSB), compound name, CID, or CAS number (PubChem).
Also supports direct SMILES input for custom molecules.
Fetches ligands with priority-based lookup: local_folder -> ligands/ -> RCSB/PubChem download.
"""

import os
import re
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream


# Extended CCD codes are 1-5 alphanumeric; 4-5 char codes only exist in mmCIF.
_CCD_CODE_RE = re.compile(r'^[A-Za-z0-9]{1,5}$')


def _validate_ccd_code(code: str) -> str:
    """Validate and uppercase a CCD residue code (1-5 alphanumeric, extended format)."""
    if not code:
        raise ValueError("Residue code cannot be empty")
    if not _CCD_CODE_RE.match(code):
        raise ValueError(
            f"Invalid ligand code {code!r}: must be 1-5 alphanumeric characters")
    return code.upper()


class Ligand(BaseConfig):
    """
    Pipeline tool for fetching small molecule ligands from RCSB PDB, PubChem, or SMILES strings.

    Downloads SDF files and converts to PDB format with proper atom numbering.
    Implements priority-based lookup: checks local_folder (if provided), then
    ligands/ folder, then downloads from RCSB or PubChem based on lookup type.
    Also supports direct SMILES input for custom molecules.
    """

    TOOL_NAME = "Ligand"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== Ligand ==="
echo "Uses biopipelines environment (no additional installation needed)."
touch "$INSTALL_SUCCESS"
echo "=== Ligand ready ==="
"""

    # Lazy path descriptors
    # compounds is a value-based stream: its compounds.csv holds the chemistry
    # metadata (SMILES, code, …) and doubles as the TableInfo content. The
    # physical coordinate files (.pdb / .cif) live on the structures stream,
    # under structures/, with their own map_table.
    compounds_csv = Path(lambda self: self.stream_path("compounds", "compounds.csv"))
    structures_map = Path(lambda self: self.stream_map_path("structures"))
    failed_csv = Path(lambda self: self.table_path("failed"))
    config_file = Path(lambda self: self.configuration_path("fetch_config.json"))
    extract_structures_json = Path(lambda self: self.configuration_path("extract_structures.json"))
    images_folder = Path(lambda self: self.stream_folder("images") if self.generate_images else None)
    compound_images_py = Path(lambda self: self.pipe_script_path("pipe_compound_images.py"))
    ligand_py = Path(lambda self: self.pipe_script_path("pipe_ligand.py"))

    def __init__(self,
                 lookup: Optional[Union[str, List[str], Dict[str, str]]] = None,
                 ids: Optional[Union[str, List[str]]] = None,
                 codes: Optional[Union[str, List[str]]] = None,
                 code: Optional[Union[str, List[str]]] = None,
                 source: Optional[str] = None,
                 local_folder: Optional[str] = None,
                 smiles: Optional[Union[str, List[str], Dict[str, str]]] = None,
                 structures: Optional[Union['DataStream', 'StandardizedOutput']] = None,
                 generate_images: bool = False,
                 **kwargs):
        """
        Initialize Ligand tool.

        Args:
            lookup: Lookup value(s) for fetching. Can be:
                    - RCSB CCD codes: "ATP", "GDP", "HEM"
                    - PubChem CID: "2244"
                    - PubChem CAS: "50-78-2"
                    - PubChem name: "aspirin", "caffeine"
                    - Path to a .txt file: one SMILES per line
                    - Path to a .cdxml file: ChemDraw molecules
                    - Dictionary mapping IDs to lookup values: {"lig1": "ATP", "lig2": "GDP"}
                      Equivalent to Ligand(lookup=list(values), ids=list(keys))
                    Can be None if using smiles instead.
            ids: Output identifier(s) for filenames (e.g., "my_ligand" -> my_ligand.pdb).
                 If not provided, defaults to lookup values (for lookup), "smilesN" (for smiles),
                 or names/indices from CDXML (for cdxml).
                 Ignored when lookup or smiles is a dictionary (ids come from dict keys).
            codes: residue code(s) to carry on the compounds stream (e.g., "LIG").
                   1-5 alphanumeric (extended CCD). If not provided, defaults to the
                   lookup value (for lookup) or "LIG" (for smiles/cdxml).
            code: Code-only construction. `Ligand(code="ZIT")` builds a compounds
                  stream that merely names an existing HETATM residue code — no
                  download, no SMILES, no structures stream. The result is a
                  value-based compounds csv (format="csv", code set, smiles empty),
                  used to hand a residue code to HETATM-selector tools. Mutually
                  exclusive with lookup / smiles / codes. Accepts a list.
            source: Force source ("rcsb" or "pubchem"). If None, auto-detects.
                    Ignored when using smiles or cdxml.
            local_folder: Custom local folder to check first (before ligands/). Default: None
            smiles: SMILES string(s) for direct molecule input. Bypasses lookup entirely.
                    Can also be a dictionary mapping IDs to SMILES: {"lig1": "CCO", "lig2": "CC"}
                    Equivalent to Ligand(smiles=list(values), ids=list(keys))
            structures: a DataStream / StandardizedOutput of complex structures
                    (PDB/CIF) to carve a bound ligand out of. Requires `codes` (the
                    HETATM residue code(s) to extract). For each input structure the
                    matching HETATM block is written to a coordinate file KEEPING the
                    bound coordinates — no download, no SMILES, no bond-order
                    templating (run OpenBabel(structures=..., convert_3d="sdf") after
                    if a tool needs an SDF). A code absent from a structure is routed
                    to the `failed` table. Mutually exclusive with lookup/smiles/code.
                    Example: Ligand(structures=complex, codes="STI").
            generate_images: Generate PNG images for each ligand using RDKit. Default: False
            **kwargs: Additional parameters

        Output:
            Streams: compounds (value-based csv: chemistry/metadata), structures
                     (SDF for download/SMILES; the input structure's own format for
                     structures=... extract; absent in code-only mode), images
                     (.png, if generate_images=True). For other coordinate formats
                     run OpenBabel(compounds=lig, convert_3d="pdb"|"cif").
            Tables:
                compounds: id | format | code | lookup | source | ccd | cid | cas | smiles | name | formula | file_path
                failed: lookup | error_message | source | attempted_path
        """
        if "output_format" in kwargs:
            raise ValueError(
                "output_format was removed from Ligand: the structures stream is "
                "always SDF (download/SMILES) or the input structure's own format "
                "(structures=... extract). For PDB/CIF coordinates run "
                "OpenBabel(compounds=lig, convert_3d=\"pdb\"|\"cif\").")

        # Code-only construction: Ligand(code="ZIT"). Names an existing HETATM
        # residue code with no chemistry — produces a value-based compounds csv
        # (smiles empty) and no structures stream. Mutually exclusive with the
        # download/generation paths.
        self.code_only = False
        self._structures_only = False
        if code is not None:
            if lookup is not None or smiles is not None or codes is not None:
                raise ValueError("code=... is mutually exclusive with lookup, smiles, and codes")
            self.code_only = True
            if isinstance(code, str):
                code_list = [code]
            else:
                code_list = list(code)
            if not code_list:
                raise ValueError("code cannot be empty")
            # Ligand is the sole validator of the code; enforce 1-5 alphanumeric.
            self.residue_codes = [_validate_ccd_code(c) for c in code_list]
            # ids default to the codes themselves
            if ids is not None:
                self.custom_ids = [ids] if isinstance(ids, str) else list(ids)
            else:
                self.custom_ids = list(self.residue_codes)
            if len(self.custom_ids) != len(self.residue_codes):
                raise ValueError(
                    f"Length mismatch: ids has {len(self.custom_ids)} items but code has {len(self.residue_codes)}")
            # No download/generation inputs in code-only mode
            self.lookup_values = []
            self.smiles_values = []
            self.file_smiles_ids = []
            self.extract_structures_stream = None
            self.source = None
            self.local_folder = None
            # code-only names an existing HETATM: no coordinate file is written.
            self.structures_format = None
            self.generate_images = generate_images
            super().__init__(**kwargs)
            return

        # Posed-ligand modifier: Ligand("STI", structures=complex). `structures=`
        # does NOT replace the lookup/smiles path — it changes only WHERE the 3-D
        # coordinates come from. The chemistry (SMILES, code, metadata) still comes
        # from the normal lookup/smiles flow below; at runtime the coordinate file
        # for each id is carved out of the matching bound HETATM in the input
        # structure (keeping crystal coords) instead of being downloaded/embedded.
        # The compounds stream therefore still carries SMILES, so a downstream
        # OpenBabel(compounds=this, convert_3d="sdf") produces a posed,
        # bond-order-correct SDF. A code absent from every input structure is
        # routed to the failed table at runtime.
        self.extract_structures_stream = None
        self._structures_arg = structures
        if structures is not None:
            if code is not None:
                raise ValueError("structures=... is not compatible with code=... (code-only names a HETATM, it has no coordinates to source)")
            if isinstance(structures, StandardizedOutput):
                self.extract_structures_stream = structures.streams.structures
            elif isinstance(structures, DataStream):
                self.extract_structures_stream = structures
            else:
                raise ValueError(
                    f"structures must be DataStream or StandardizedOutput, got {type(structures).__name__}")

            # structures= with only codes (no lookup/smiles): carve the bound
            # HETATM, keeping coords. No chemistry source → compounds stream has
            # the code but no SMILES, so downstream bond-order templating is not
            # possible (perception is used instead). When lookup/smiles IS also
            # given, fall through to the normal path: chemistry comes from there
            # and the coordinates are still carved (config carries the stream).
            if lookup is None and smiles is None:
                if codes is None:
                    raise ValueError("structures=... requires codes=... (the HETATM residue code(s) to extract) when no lookup/smiles is given")
                self._structures_only = True
                code_src = [codes] if isinstance(codes, str) else list(codes)
                self.residue_codes = [_validate_ccd_code(c) for c in code_src]
                self.custom_ids = ([ids] if isinstance(ids, str) else list(ids)) if ids is not None else list(self.residue_codes)
                if len(self.custom_ids) != len(self.residue_codes):
                    raise ValueError(
                        f"Length mismatch: ids has {len(self.custom_ids)} items but codes has {len(self.residue_codes)}")
                self.lookup_values = []
                self.smiles_values = []
                self.file_smiles_ids = []
                self.source = None
                self.local_folder = None
                # Extract carves the bound HETATM block as-is, keeping crystal
                # coords. The carver dispatches on each input file's extension
                # (pdb/cif), so the output stream mirrors the input's format.
                self.structures_format = self.extract_structures_stream.format or "pdb"
                self.generate_images = generate_images
                super().__init__(**kwargs)
                if isinstance(structures, StandardizedOutput) and hasattr(structures, "config"):
                    self.dependencies.append(structures.config)
                return

        # Dict input for smiles: {id: smiles_string} -> extract ids and smiles
        _dict_ids = None
        if isinstance(smiles, dict):
            if ids is not None:
                print("  Warning: 'ids' parameter ignored when smiles is a dictionary (using dict keys)")
            if not smiles:
                raise ValueError("smiles dictionary cannot be empty")
            _dict_ids = list(smiles.keys())
            smiles = list(smiles.values())

        # Dict input for lookup: {id: lookup_value} -> extract ids and lookups
        if isinstance(lookup, dict):
            if ids is not None and _dict_ids is None:
                print("  Warning: 'ids' parameter ignored when lookup is a dictionary (using dict keys)")
            if not lookup:
                raise ValueError("lookup dictionary cannot be empty")
            if _dict_ids is not None:
                # Both smiles and lookup are dicts - combine ids
                _dict_ids = list(lookup.keys()) + _dict_ids
            else:
                _dict_ids = list(lookup.keys())
            lookup = list(lookup.values())

        if _dict_ids is not None:
            ids = _dict_ids

        # Handle smiles input
        if smiles is not None:
            if isinstance(smiles, str):
                self.smiles_values = [smiles]
            else:
                self.smiles_values = list(smiles)
        else:
            self.smiles_values = []

        # Handle lookup - detect file paths (.txt, .cdxml) and extract SMILES from them
        if lookup is not None:
            if isinstance(lookup, str):
                raw_lookups = [lookup]
            else:
                raw_lookups = list(lookup)
        else:
            raw_lookups = []

        self.lookup_values = []
        self.file_smiles_ids = []  # (smiles_list, ids_list) pairs from parsed files
        for val in raw_lookups:
            lower = val.lower()
            if lower.endswith('.txt'):
                file_smiles, file_ids = self._parse_txt_smiles(val)
                self.smiles_values.extend(file_smiles)
                self.file_smiles_ids.append((val, file_ids))
            elif lower.endswith('.cdxml'):
                file_smiles, file_ids = self._parse_cdxml_ligands(val)
                self.smiles_values.extend(file_smiles)
                self.file_smiles_ids.append((val, file_ids))
            else:
                self.lookup_values.append(val)

        # Validate: must have at least one of lookup or smiles
        if not self.lookup_values and not self.smiles_values:
            raise ValueError("Must provide at least one of 'lookup' or 'smiles'")

        # Total number of ligands
        total_count = len(self.lookup_values) + len(self.smiles_values)

        # Handle ids - default based on input type
        if ids is not None:
            if isinstance(ids, str):
                self.custom_ids = [ids]
            else:
                self.custom_ids = list(ids)
        else:
            # Default ids: lookup values, then user-supplied smiles as "smilesN", then file-derived IDs
            self.custom_ids = self.lookup_values.copy()
            # Count SMILES that came from files
            n_file_smiles = sum(len(ids) for _, ids in self.file_smiles_ids)
            n_plain_smiles = len(self.smiles_values) - n_file_smiles
            for i in range(n_plain_smiles):
                self.custom_ids.append(f"smiles{i + 1}")
            for _, file_ids in self.file_smiles_ids:
                self.custom_ids.extend(file_ids)

        # Handle codes - default based on input type
        if codes is not None:
            if isinstance(codes, str):
                self.residue_codes = [codes.upper()]
            else:
                self.residue_codes = [c.upper() for c in codes]
        else:
            # Default codes: the lookup value if it is itself a valid CCD code,
            # else "LIG" (a SMILES/name lookup has no meaningful residue code).
            self.residue_codes = [
                lv.upper() if _CCD_CODE_RE.match(lv) else "LIG"
                for lv in self.lookup_values
            ]
            for _ in range(len(self.smiles_values)):
                self.residue_codes.append("LIG")

        # Validate lengths
        if len(self.custom_ids) != total_count:
            raise ValueError(f"Length mismatch: ids has {len(self.custom_ids)} items but total ligands is {total_count}")
        if len(self.residue_codes) != total_count:
            raise ValueError(f"Length mismatch: codes has {len(self.residue_codes)} items but total ligands is {total_count}")

        # Validate source
        if source is not None and source not in ["rcsb", "pubchem"]:
            raise ValueError(f"Invalid source: {source}. Must be 'rcsb', 'pubchem', or None")
        self.source = source

        self.local_folder = local_folder

        # Download/SMILES paths emit SDF coordinates (native RCSB/RDKit form, no
        # fixed-column width limit, so an extended code is never a problem).
        self.structures_format = "sdf"
        self.generate_images = generate_images

        for code in self.residue_codes:
            _validate_ccd_code(code)

        # Initialize base class
        super().__init__(**kwargs)

        # structures= modifier: depend on the upstream tool producing the coords.
        if (self._structures_arg is not None
                and isinstance(self._structures_arg, StandardizedOutput)
                and hasattr(self._structures_arg, "config")):
            self.dependencies.append(self._structures_arg.config)

    @staticmethod
    def _parse_txt_smiles(txt_path: str):
        """
        Parse SMILES strings from a text file (one SMILES per line).

        The file name (without extension) is used as the root for IDs:
        e.g., "myligands.txt" -> myligands1, myligands2, ...

        Returns:
            (smiles_list, ids_list) - parallel lists
        """
        if not os.path.exists(txt_path):
            raise ValueError(f"TXT file not found: {txt_path}")

        root = os.path.splitext(os.path.basename(txt_path))[0]

        with open(txt_path, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]

        if not lines:
            raise ValueError(f"No SMILES found in TXT file: {txt_path}")

        smiles_list = lines
        ids_list = [f"{root}{i + 1}" for i in range(len(lines))]

        return smiles_list, ids_list

    @staticmethod
    def _parse_cdxml_ligands(cdxml_path: str):
        """
        Parse individual molecule fragments from a CDXML file.

        Each top-level <fragment> on the page is converted to a SMILES string.
        Names come from <chemicalproperty> elements whose BasisObjects reference the fragment;
        falls back to "ligandN" (1-based) if no name is defined.

        Returns:
            (smiles_list, ids_list) - parallel lists, one entry per molecule
        """
        import xml.etree.ElementTree as ET
        try:
            from rdkit import Chem, RDLogger
        except ImportError:
            raise ImportError("RDKit is required for CDXML parsing. Install with: conda install -c conda-forge rdkit")

        if not os.path.exists(cdxml_path):
            raise ValueError(f"CDXML file not found: {cdxml_path}")

        RDLogger.DisableLog('rdApp.warning')

        tree = ET.parse(cdxml_path)
        root = tree.getroot()
        page = root.find('page')
        if page is None:
            raise ValueError(f"No <page> element found in CDXML file: {cdxml_path}")

        CDXML_HEADER = '<?xml version="1.0" encoding="UTF-8" ?><CDXML><page>'
        CDXML_FOOTER = '</page></CDXML>'

        # Build id -> text for all <t> elements on the page
        t_id_to_text = {}
        for child in page:
            if child.tag == 't':
                tid = child.attrib.get('id', '')
                s_elem = child.find('s')
                if s_elem is not None and s_elem.text:
                    t_id_to_text[tid] = s_elem.text.strip()

        # Map fragment_id -> name via chemicalproperty (same as CompoundLibrary)
        frag_id_to_name = {}
        for child in page:
            if child.tag == 'chemicalproperty':
                display_id = child.attrib.get('ChemicalPropertyDisplayID', '')
                basis = child.attrib.get('BasisObjects', '').split()
                name = t_id_to_text.get(display_id, '')
                if name and basis:
                    frag_id_to_name[basis[0]] = name

        # Parse each top-level fragment as a separate molecule
        smiles_list = []
        ids_list = []
        frag_index = 0
        for child in page:
            if child.tag != 'fragment':
                continue
            frag_id = child.attrib.get('id', '')
            frag_str = ET.tostring(child, encoding='unicode')
            mols = Chem.MolsFromCDXML(CDXML_HEADER + frag_str + CDXML_FOOTER)
            if not mols or mols[0] is None:
                continue
            mol = mols[0]
            smiles = Chem.MolToSmiles(mol)
            if not smiles:
                continue
            frag_index += 1
            name = frag_id_to_name.get(frag_id, f"ligand{frag_index}")
            smiles_list.append(smiles)
            ids_list.append(name)

        RDLogger.EnableLog('rdApp.warning')

        if not smiles_list:
            raise ValueError(f"No valid molecules found in CDXML file: {cdxml_path}")

        return smiles_list, ids_list

    def _detect_lookup_type(self, lookup: str) -> str:
        """
        Detect the type of lookup value.

        Returns: "ccd" (RCSB), "cid" (PubChem), "cas" (PubChem), or "name" (PubChem)
        """
        # CCD codes: 1-5 alphanumeric, canonically uppercase. A value carrying a
        # lowercase letter (water, urea, aspirin) is a PubChem name, not a CCD.
        if _CCD_CODE_RE.match(lookup) and not lookup.isdigit() and lookup == lookup.upper():
            return "ccd"

        # PubChem CID: purely numeric
        if lookup.isdigit():
            return "cid"

        # CAS number: XX-XX-X format (digits-digits-digit)
        if re.match(r'^\d+-\d+-\d$', lookup):
            return "cas"

        # Default: compound name
        return "name"

    def validate_params(self):
        """Validate Ligand parameters."""
        if not self.custom_ids:
            raise ValueError("ids cannot be empty")

        if not self.residue_codes:
            raise ValueError("codes cannot be empty")

        # structures= modifier: the coordinate source must be non-empty. The
        # chemistry counts are validated by the lookup/smiles path below.
        if self.extract_structures_stream is not None and len(self.extract_structures_stream) == 0:
            raise ValueError("structures input must not be empty")

        if self.code_only or self._structures_only:
            # No chemistry source: per-id codes only (code-only names a HETATM;
            # structures-only carves the bound HETATM keeping coords).
            if len(self.custom_ids) != len(self.residue_codes):
                raise ValueError(
                    f"ids length ({len(self.custom_ids)}) must match code count ({len(self.residue_codes)})")
        else:
            if not self.lookup_values and not self.smiles_values:
                raise ValueError("Must have at least one of 'lookup', 'smiles', or 'code'")

            total_count = len(self.lookup_values) + len(self.smiles_values)
            if len(self.custom_ids) != total_count:
                raise ValueError(f"ids length ({len(self.custom_ids)}) must match total ligands ({total_count})")

            if len(self.residue_codes) != total_count:
                raise ValueError(f"codes length ({len(self.residue_codes)}) must match total ligands ({total_count})")

        for code in self.residue_codes:
            if not code:
                raise ValueError("Residue code cannot be empty")

        _validate_freeform_string("source", self.source)
        _validate_freeform_string("local_folder", self.local_folder)
        for i, code in enumerate(self.residue_codes):
            _validate_freeform_string(f"codes[{i}]", code)

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input parameters and check for local files."""
        self.folders = pipeline_folders


        # Check which files exist locally and which will need to be downloaded
        repo_ligands_folder = pipeline_folders.get('ligands', '')
        self.found_locally = []
        self.needs_download = []

        # Process lookup values
        for lookup, custom_id in zip(self.lookup_values, self.custom_ids[:len(self.lookup_values)]):
            found = False
            local_path = None

            # Check local_folder if specified
            if self.local_folder:
                local_path = os.path.join(self.local_folder, f"{lookup}.pdb")
                if os.path.exists(local_path):
                    self.found_locally.append((lookup, local_path))
                    found = True
                    print(f"  Found ligand {lookup} locally: {local_path}")

            # Check ligands/ folder
            if not found and repo_ligands_folder:
                local_path = os.path.join(repo_ligands_folder, f"{lookup}.pdb")
                if os.path.exists(local_path):
                    self.found_locally.append((lookup, local_path))
                    found = True
                    print(f"  Found ligand {lookup} locally: {local_path}")

            if not found:
                # Determine source for display
                lookup_type = self._detect_lookup_type(lookup)
                effective_source = self.source if self.source else ("rcsb" if lookup_type == "ccd" else "pubchem")

                self.needs_download.append(lookup)
                print(f"  Ligand {lookup} not found locally, will download from {effective_source}")

        # Process SMILES values (these are always "generated" at runtime)
        smiles_start_idx = len(self.lookup_values)
        for i, smiles in enumerate(self.smiles_values):
            custom_id = self.custom_ids[smiles_start_idx + i]
            smiles_preview = smiles[:30] + "..." if len(smiles) > 30 else smiles
            print(f"  SMILES ligand {custom_id}: {smiles_preview} (will be generated with RDKit)")

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"IDS: {', '.join(self.custom_ids)} ({len(self.custom_ids)} ligands)",
            f"CODES: {', '.join(self.residue_codes)}",
        ])
        if self.code_only:
            config_lines.append("MODE: code-only (names existing HETATM; no structures)")
            return config_lines
        config_lines.append(f"FORMAT: {self.structures_format.upper()}")
        if self.extract_structures_stream is not None:
            config_lines.append(
                f"COORDS: carved from {len(self.extract_structures_stream)} bound structure(s) "
                f"(chemistry from lookup/smiles, coordinates from HETATM)")

        if self.lookup_values:
            config_lines.append(f"LOOKUP: {', '.join(self.lookup_values)}")
            config_lines.append(f"SOURCE: {self.source if self.source else 'auto-detect'}")
            config_lines.append(f"LOCAL_FOLDER: {self.local_folder if self.local_folder else 'None (uses ligands/)'}")

        if self.file_smiles_ids:
            for file_path, file_ids in self.file_smiles_ids:
                config_lines.append(f"FILE: {os.path.basename(file_path)} ({len(file_ids)} molecules)")
        elif self.smiles_values:
            smiles_preview = [s[:20] + "..." if len(s) > 20 else s for s in self.smiles_values]
            config_lines.append(f"SMILES: {', '.join(smiles_preview)} ({len(self.smiles_values)} molecules)")

        # Add status of files found/not found
        if hasattr(self, 'found_locally') and self.found_locally:
            config_lines.append(f"FOUND_LOCALLY: {', '.join([lookup for lookup, _ in self.found_locally])}")
        if hasattr(self, 'needs_download') and self.needs_download:
            config_lines.append(f"NEEDS_DOWNLOAD: {', '.join(self.needs_download)}")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate script to fetch ligands from local folders or RCSB/PubChem."""
        # Layout sub-dirs already created by the pipeline.
        script_content = "#!/bin/bash\n"
        script_content += "# Ligand execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self.generate_script_run_ligand()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_ligand(self) -> str:
        """Generate the ligand fetching part of the script."""
        import json

        repo_ligands_folder = self.folders['ligands']

        config_data = {
            "custom_ids": self.custom_ids,
            "residue_codes": self.residue_codes,
            "lookup_values": self.lookup_values,
            "smiles_values": self.smiles_values,
            "source": self.source,
            "local_folder": self.local_folder,
            "output_format": self.structures_format,
            "code_only": self.code_only,
            "repo_ligands_folder": repo_ligands_folder,
            # Coordinate files live on the structures stream; compounds.csv (the
            # value-based chemistry table) lives in the compounds stream folder.
            "output_folder": self.stream_folder("structures"),
            "compounds_table": self.compounds_csv,
            "structures_table": self.structures_map,
            "failed_table": self.failed_csv
        }

        # structures= modifier: hand the pipe script the input complexes (id ->
        # file) it carves the per-id HETATM `code` out of, keeping bound coords —
        # the chemistry (smiles/metadata) still comes from the lookup/smiles path.
        if self.extract_structures_stream is not None:
            self.extract_structures_stream.save_json(self.extract_structures_json)
            config_data["extract_structures_json"] = self.extract_structures_json

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        total_count = len(self.lookup_values) + len(self.smiles_values)
        lookup_info = f"Lookup: {', '.join(self.lookup_values)}" if self.lookup_values else ""
        smiles_info = f"SMILES: {len(self.smiles_values)} molecule(s)" if self.smiles_values else ""

        image_script = ""
        if self.generate_images:
            image_script = f"""
echo "Generating ligand images"
python3 "{self.compound_images_py}" "{self.compounds_csv}" "{self.images_folder}"
"""

        return f"""echo "Processing {total_count} ligands"
echo "Custom IDs: {', '.join(self.custom_ids)}"
echo "Residue codes: {', '.join(self.residue_codes)}"
{f'echo "{lookup_info}"' if lookup_info else ''}
{f'echo "{smiles_info}"' if smiles_info else ''}
echo "Output folder: {self.output_folder}"

python "{self.ligand_py}" --config "{self.config_file}"
{image_script}
"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after ligand fetching."""
        # "sdf" (download/smiles) or the input's format (extract, possibly mixed).
        # A mixed "pdb|cif" extract can't predict a single extension -> glob.
        ext = "*" if "|" in (self.structures_format or "") else self.structures_format

        tables = {
            "compounds": TableInfo(
                name="compounds",
                path=self.compounds_csv,
                columns=["id", "format", "code", "lookup", "source", "ccd", "cid", "cas", "smiles", "name", "formula", "file_path"],
                description="Successfully fetched/generated ligand files with metadata"
            ),
            "failed": TableInfo(
                name="failed",
                path=self.failed_csv,
                columns=["lookup", "error_message", "source", "attempted_path"],
                description="Failed ligand fetches with error details"
            )
        }

        # compounds is always a value-based csv (chemistry/metadata in map_table).
        compounds = DataStream(
            name="compounds",
            ids=self.custom_ids.copy(),
            files=[],
            map_table=self.compounds_csv,
            format="csv"
        )

        result = {
            "compounds": compounds,
            "tables": tables,
            "output_folder": self.output_folder
        }

        # Code-only mode names an existing HETATM — no coordinate files, so
        # no structures stream. Otherwise the coordinate files live on the
        # structures stream (templated, with its own map_table).
        if self.code_only:
            result["structures"] = DataStream.empty("structures", "sdf")
        else:
            structures = DataStream(
                name="structures",
                ids=self.custom_ids.copy(),
                files=[self.stream_path("structures", f"<id>.{ext}")],
                map_table=self.structures_map,
                format=self.structures_format
            )
            result["structures"] = structures

        if self.generate_images and not self.code_only:
            image_files = [os.path.join(self.images_folder, f"{cid}.png")
                           for cid in self.custom_ids]
            result["images"] = DataStream(
                name="images",
                ids=self.custom_ids.copy(),
                files=image_files,
                format="png"
            )

        return result

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                "custom_ids": self.custom_ids,
                "residue_codes": self.residue_codes,
                "lookup_values": self.lookup_values,
                "smiles_values": self.smiles_values,
                "code_only": self.code_only,
                "source": self.source,
                "local_folder": self.local_folder,
                "structures_format": self.structures_format,
                "parsed_files": [fp for fp, _ in self.file_smiles_ids],
                "generate_images": self.generate_images
            }
        })
        return base_dict
