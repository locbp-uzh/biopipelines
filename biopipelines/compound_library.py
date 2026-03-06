# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
CompoundLibrary configuration for generating and processing ligand libraries.

Handles dictionary-based SMILES library generation, CSV output generation,
and optional covalent ligand CCD/PKL file preparation for Boltz2.
"""

import os
import json
import csv
import itertools
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


class CompoundLibrary(BaseConfig):
    """
    CompoundLibrary configuration for dictionary-based SMILES library generation.

    Expands dictionaries with substitution keys into complete compound libraries,
    generates CSV files with standardized format, and optionally prepares
    covalent ligand files for Boltz2.
    """

    TOOL_NAME = "CompoundLibrary"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== CompoundLibrary ==="
echo "Uses biopipelines environment (no additional installation needed)."
echo "=== CompoundLibrary ready ==="
"""

    # Lazy path descriptors
    compounds_csv = Path(lambda self: os.path.join(self.output_folder, "compounds.csv"))
    compound_properties_csv = Path(lambda self: os.path.join(self.output_folder, "compound_properties.csv"))
    summary_file = Path(lambda self: os.path.join(self.output_folder, "summary.txt"))
    library_dict_json = Path(lambda self: os.path.join(self.output_folder, "library_dict.json"))
    covalent_folder = Path(lambda self: os.path.join(self.output_folder, "covalent_library") if self.covalent else None)
    covalent_compounds_csv = Path(lambda self: os.path.join(self.output_folder, "covalent_library", "compounds.csv") if self.covalent else None)
    compound_expansion_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_compound_library.py"))
    smiles_properties_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_smiles_properties.py"))
    covalent_generation_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_compound_library.py"))

    def __init__(self,
                 library: Union[str, Dict[str, Union[str, List[str]]]],
                 primary_key: Optional[str] = None,
                 covalent: bool = False,
                 validate_smiles: bool = True,
                 conformer_method: str = "UFF",
                 **kwargs):
        """
        Initialize CompoundLibrary configuration.

        Args:
            library: Dictionary with expansion keys or path to existing CSV library
            primary_key: Root key for expansion when library is a dictionary
            covalent: Generate CCD/PKL files for covalent ligand binding (calls runtime script)
            validate_smiles: Validate SMILES strings during expansion
            conformer_method: Method for conformer generation ("UFF", "OpenFF", "DFT")
            **kwargs: Additional parameters

        Output:
            Streams: compounds (.csv)
            Tables:
                compounds: id | format | smiles | ccd | ...branching_keys
                covalent_compounds: id | format | smiles | ccd (if covalent=True)
        """
        # Store CompoundLibrary-specific parameters
        self.library = library
        self.primary_key = primary_key
        self.covalent = covalent
        self.validate_smiles = validate_smiles
        self.conformer_method = conformer_method

        # Track library source type
        self.library_dict = None
        self.library_csv = None
        self.library_cdxml = None
        self.expanded_compounds = []
        self.compound_ids = []

        # Initialize base class
        super().__init__(**kwargs)

        # For dictionary libraries, do expansion immediately to get compound IDs
        if isinstance(self.library, dict):
            self.library_dict = self.library
            self._expand_library()
        elif isinstance(self.library, str) and self.library.endswith('.cdxml'):
            self.library_cdxml = self.library
            self._expand_cdxml()

    def validate_params(self):
        """Validate CompoundLibrary-specific parameters."""
        if not self.library:
            raise ValueError("library parameter is required")

        # Validate library format
        if isinstance(self.library, dict):
            # Dictionary-based library
            if self.primary_key and self.primary_key not in self.library:
                raise ValueError(f"primary_key '{self.primary_key}' not found in library dictionary")
        elif isinstance(self.library, str):
            # CSV file path
            if not (self.library.endswith('.csv') or self.library.endswith('.cdxml')):
                raise ValueError("library file must have .csv or .cdxml extension")
        else:
            raise ValueError("library must be a dictionary or CSV file path")

        # Validate conformer method for covalent ligands
        if self.covalent:
            valid_methods = ["UFF", "OpenFF", "DFT"]
            if self.conformer_method not in valid_methods:
                raise ValueError(f"conformer_method must be one of: {valid_methods}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input library sources."""
        self.folders = pipeline_folders

        if isinstance(self.library, str) and self.library.endswith('.cdxml'):
            # CDXML file - resolve path
            if os.path.exists(self.library):
                self.library_cdxml = self.library
            else:
                project_path = os.path.join(pipeline_folders["biopipelines"], self.library)
                if os.path.exists(project_path):
                    self.library_cdxml = project_path
                else:
                    raise ValueError(f"CDXML file not found: {self.library}")
            # Expand if not already done (file resolved via project path)
            if not self.expanded_compounds:
                self._expand_cdxml()
        elif isinstance(self.library, str):
            # CSV file - check if it exists
            if os.path.exists(self.library):
                self.library_csv = self.library
            else:
                # Try in project directory
                project_path = os.path.join(pipeline_folders["biopipelines"], self.library)
                if os.path.exists(project_path):
                    self.library_csv = project_path
                else:
                    raise ValueError(f"Library CSV file not found: {self.library}")
        elif isinstance(self.library, dict):
            # Dictionary-based library (expansion already done in __init__)
            if not self.library_dict:
                self.library_dict = self.library
                self._expand_library()

    def _expand_library(self):
        """Expand dictionary library into individual compounds and generate IDs."""
        if not self.library_dict:
            return

        # Use the same expansion logic as boltz_compound_library.py but with <key> format
        library = self.library_dict.copy()
        library_keys = list(library.keys())
        primary_lib_key = self.primary_key

        # Find the primary key in the base configuration
        if primary_lib_key and primary_lib_key not in library:
            raise ValueError(f"Primary key '{primary_lib_key}' not found in library")

        # Expand each base compound from the primary key
        if primary_lib_key:
            final_compounds = []
            primary_value = library[primary_lib_key]

            # Handle both string and list values for primary key
            if isinstance(primary_value, str):
                # Single SMILES string
                final_compounds.append({'smiles': primary_value, 'branching': {}})
            elif isinstance(primary_value, list):
                # List of SMILES strings
                for base in primary_value:
                    final_compounds.append({'smiles': base, 'branching': {}})
            else:
                raise ValueError(f"Primary key '{primary_lib_key}' must be a string or list of strings")

            no_new_branching = False
            while not no_new_branching:
                no_new_branching = True
                updated_compounds = []

                # Process every compound in current list
                for compound in final_compounds:
                    key_found = False
                    # Check for every library key in the current compound's SMILES
                    for key in library_keys:
                        # Use <key> format instead of *key*
                        key_pattern = f"<{key}>"
                        if key_pattern in compound['smiles']:
                            key_found = True
                            no_new_branching = False
                            # For every possible substitution for the key, create a new compound
                            key_options = library[key]
                            # Handle both string and list values for expansion keys
                            if isinstance(key_options, str):
                                key_options = [key_options]

                            for option in key_options:
                                new_smiles = compound['smiles'].replace(key_pattern, option, 1)
                                new_branching = compound['branching'].copy()
                                # Strip <> from option for display (e.g. "<o-hydroxyphenyl>" -> "o-hydroxyphenyl")
                                clean_option = option.strip("<>") if option.startswith("<") and option.endswith(">") else option
                                new_branching[key] = clean_option
                                updated_compounds.append({'smiles': new_smiles, 'branching': new_branching})
                            # Process one key per compound per iteration
                            break

                    if not key_found:
                        updated_compounds.append(compound)

                final_compounds = updated_compounds

            # Generate compound IDs based on the pattern from boltz_compound_library.py
            num_compounds = len(final_compounds)
            characters = 4
            if num_compounds > 9: characters = 3
            if num_compounds > 99: characters = 2
            if num_compounds > 999: characters = 1
            if num_compounds > 99999: characters = 0

            compound_ids = []
            for u_l_n in range(num_compounds):
                u_l_n_str = str(u_l_n)
                n0 = 5 - characters - len(u_l_n_str)
                zeros_str = '0' * n0
                compound_name = primary_lib_key if num_compounds == 1 else primary_lib_key[:characters] + zeros_str + u_l_n_str
                compound_ids.append(compound_name)

            # Store expanded compounds
            self.expanded_compounds = final_compounds
            self.compound_ids = compound_ids
        else:
            # There is no expansion, every key corresponds to an item
            compound_ids = []
            final_compounds = []
            for name, smiles in library.items():
                compound_ids.append(name)
                final_compounds.append({'smiles': smiles, 'branching': {}})
            self.expanded_compounds = final_compounds
            self.compound_ids = compound_ids

    def _expand_cdxml(self):
        """
        Expand CDXML file using ChemDraw's native R-group substitution table format.

        Parses <altgroup> elements (R-group tables) and <bracketedgroup> elements
        (repeat units / polymer brackets) directly from the CDXML XML, then enumerates
        all combinations via RDKit molzip.
        """
        try:
            import warnings
            import xml.etree.ElementTree as ET
            from rdkit import Chem, RDLogger
            from rdkit.Chem import rdmolops, RWMol
        except ImportError:
            raise ImportError(
                "RDKit is required for CDXML R-group enumeration. "
                "Install with: conda install -c conda-forge rdkit"
            )

        # Suppress RDKit "Incomplete atom labelling" warnings that fire while parsing
        # scaffold/fragment CDXMLs that have unmatched ExternalConnectionPoint atoms.
        RDLogger.DisableLog('rdApp.warning')

        cdxml_path = self.library_cdxml
        if not cdxml_path or not os.path.exists(cdxml_path):
            raise ValueError(f"CDXML file not found: {cdxml_path}")

        # --- Parse XML ---
        tree = ET.parse(cdxml_path)
        root = tree.getroot()
        page = root.find('page')
        if page is None:
            raise ValueError("No <page> element found in CDXML file.")

        # Collect all existing numeric IDs to assign collision-free synthetic IDs
        all_ids = set()
        for elem in root.iter():
            for attr in ('id', 'B', 'E'):
                val = elem.attrib.get(attr, '')
                if val.isdigit():
                    all_ids.add(int(val))
        _next_id = [max(all_ids, default=0) + 100000]

        def _fresh_id():
            nid = _next_id[0]
            _next_id[0] += 1
            return str(nid)

        CDXML_HEADER = '<?xml version="1.0" encoding="UTF-8" ?><CDXML><page>'
        CDXML_FOOTER = '</page></CDXML>'

        def _make_dummy(rwmol, atom_idx, mapnum):
            atom = rwmol.GetAtomWithIdx(atom_idx)
            atom.SetAtomicNum(0)
            atom.SetAtomMapNum(mapnum)
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(0)

        def _parse_fragment_mol(frag_elem, conn_to_mapnum):
            """
            Parse a <fragment> element via MolsFromCDXML.
            ExternalConnectionPoint atoms (identified by CDX_NODE_ID) are converted to
            dummy * atoms with atomMapNum = conn_to_mapnum[ExternalConnectionNum].

            conn_to_mapnum: {ExternalConnectionNum (int) -> global_mapnum (int)}
            """
            # Build node_id -> ExternalConnectionNum from XML
            frag_conn = {}
            for n in frag_elem.findall('n'):
                if n.attrib.get('NodeType') == 'ExternalConnectionPoint':
                    frag_conn[n.attrib['id']] = int(n.attrib.get('ExternalConnectionNum', '1'))

            frag_str = ET.tostring(frag_elem, encoding='unicode')
            mols = Chem.MolsFromCDXML(CDXML_HEADER + frag_str + CDXML_FOOTER)
            if not mols or mols[0] is None:
                return None
            rwmol = RWMol(mols[0])
            for atom in rwmol.GetAtoms():
                cdx_id = atom.GetPropsAsDict().get('CDX_NODE_ID', None)
                if cdx_id is not None:
                    cdx_str = str(cdx_id)
                    if cdx_str in frag_conn:
                        conn_num = frag_conn[cdx_str]
                        mapnum = conn_to_mapnum.get(conn_num, conn_num)
                        _make_dummy(rwmol, atom.GetIdx(), mapnum)
            return rwmol.GetMol()

        def _frag_display_smiles(mol):
            """Return SMILES for a fragment with dummy attachment atoms removed."""
            try:
                rwmol = RWMol(mol)
                dummy_idxs = sorted(
                    [a.GetIdx() for a in rwmol.GetAtoms() if a.GetAtomicNum() == 0],
                    reverse=True
                )
                # Before removing each dummy, clear noImplicit + explicit Hs on its
                # neighbors so RDKit re-computes implicit Hs correctly after removal
                for idx in dummy_idxs:
                    for nb in rwmol.GetAtomWithIdx(idx).GetNeighbors():
                        nb.SetNoImplicit(False)
                        nb.SetNumExplicitHs(0)
                for idx in dummy_idxs:
                    rwmol.RemoveAtom(idx)
                Chem.SanitizeMol(rwmol)
                return Chem.MolToSmiles(rwmol)
            except Exception:
                return Chem.MolToSmiles(mol)

        def _assemble(scaffold_mol, frag_mols):
            """Combine scaffold + fragments and zip by atom map number."""
            params = rdmolops.MolzipParams()
            params.label = rdmolops.MolzipLabel.AtomMapNumber
            combined = scaffold_mol
            for frag in frag_mols:
                combined = Chem.CombineMols(combined, frag)
            result = Chem.molzip(combined, params)
            Chem.SanitizeMol(result)
            return result

        # --- Identify scaffold fragment and altgroups ---
        scaffold_frag_elem = None
        altgroups = {}
        for child in page:
            if child.tag == 'fragment' and scaffold_frag_elem is None:
                scaffold_frag_elem = child
            elif child.tag == 'altgroup':
                altgroups[child.attrib['id']] = child

        if scaffold_frag_elem is None:
            raise ValueError(f"No scaffold fragment found in CDXML file: {cdxml_path}")

        # --- Find R-group nodes (NamedAlternativeGroup) on scaffold ---
        rgroup_nodes = []  # (node_id, label, ag_id, bond_ordering, valence)
        for n in scaffold_frag_elem.findall('n'):
            if n.attrib.get('NodeType') == 'NamedAlternativeGroup':
                node_id = n.attrib['id']
                ag_id = n.attrib.get('AltGroupID', '')
                bond_ordering = n.attrib.get('BondOrdering', '').split()
                s_elem = n.find('./t/s')
                label = s_elem.text.strip() if s_elem is not None and s_elem.text else f'R{node_id}'
                valence = int(altgroups[ag_id].attrib.get('Valence', '1')) if ag_id in altgroups else 1
                rgroup_nodes.append((node_id, label, ag_id, bond_ordering, valence))

        # --- Find polymer/repeat unit brackets ---
        # bracketedgroup: BracketedObjectIDs (node ids inside bracket),
        #   bracketattachment/crossingbond: BondID + InnerAtomID
        # parameterizedBracketLabel on the graphic gives the repeat range (e.g. "1-3")
        bracket_groups = []  # list of dicts with parsed bracket info
        graphic_id_to_param_label = {}  # graphic id -> repeat range string
        for child in page:
            if child.tag == 'graphic' and child.attrib.get('BracketUsage') == 'MultipleGroup':
                gid = child.attrib.get('id', '')
                for tag in child.findall('objecttag'):
                    if tag.attrib.get('Name') == 'parameterizedBracketLabel':
                        s = tag.find('./t/s')
                        if s is not None and s.text:
                            graphic_id_to_param_label[gid] = s.text.strip()
            elif child.tag == 'bracketedgroup' and child.attrib.get('BracketUsage') == 'MultipleGroup':
                inner_node_ids = set(child.attrib.get('BracketedObjectIDs', '').split())
                default_repeat = int(child.attrib.get('RepeatCount', '1'))
                crossings = []  # list of (bond_id, inner_atom_node_id, graphic_id)
                for att in child.findall('bracketattachment'):
                    graphic_id = att.attrib.get('GraphicID', '')
                    for cb in att.findall('crossingbond'):
                        crossings.append((
                            cb.attrib.get('BondID', ''),
                            cb.attrib.get('InnerAtomID', ''),
                            graphic_id
                        ))
                bracket_groups.append({
                    'inner_node_ids': inner_node_ids,
                    'default_repeat': default_repeat,
                    'crossings': crossings,
                })

        # Resolve parameterized repeat ranges from graphic labels
        # Match each bracket_group's crossings to graphic_id_to_param_label
        for bg in bracket_groups:
            repeat_range = None
            for bond_id, inner_node_id, graphic_id in bg['crossings']:
                if graphic_id in graphic_id_to_param_label:
                    label_str = graphic_id_to_param_label[graphic_id]
                    if '-' in label_str:
                        parts = label_str.split('-')
                        try:
                            repeat_range = list(range(int(parts[0]), int(parts[1]) + 1))
                        except ValueError:
                            pass
                    else:
                        try:
                            repeat_range = [int(label_str)]
                        except ValueError:
                            pass
                    break
            if repeat_range is None:
                repeat_range = [bg['default_repeat']]
            bg['repeat_range'] = repeat_range

        # --- Assign global map numbers for R-groups ---
        # node_map_info: node_id -> [(global_mapnum, synthetic_id_str, ExternalConnectionNum)]
        map_counter = [1]
        node_map_info = {}
        label_order = []

        for node_id, label, ag_id, bond_ordering, valence in rgroup_nodes:
            label_order.append(label)
            if valence == 1:
                node_map_info[node_id] = [(map_counter[0], node_id, 1)]
                map_counter[0] += 1
            else:
                entries = []
                for i, bond_id in enumerate(bond_ordering):
                    syn_id = _fresh_id()
                    entries.append((map_counter[0], syn_id, i + 1))
                    map_counter[0] += 1
                node_map_info[node_id] = entries

        # --- Assign map numbers for repeat unit attachment points ---
        # Each bracket group needs 2 map numbers (left/right attachment)
        for bg in bracket_groups:
            crossings = bg['crossings']
            bg['map_left'] = map_counter[0]
            bg['map_right'] = map_counter[0] + 1
            map_counter[0] += 2
            # Identify which crossing is "left" (inner atom connects to outer-left)
            # and which is "right" — order determined by crossings list order
            if len(crossings) >= 2:
                bg['crossing_left'] = crossings[0]   # (bond_id, inner_node_id, graphic_id)
                bg['crossing_right'] = crossings[1]
            elif len(crossings) == 1:
                bg['crossing_left'] = crossings[0]
                bg['crossing_right'] = None

        # --- Build modified scaffold XML ---
        import copy
        scaffold_copy = copy.deepcopy(scaffold_frag_elem)
        id_to_node = {n.attrib['id']: n for n in scaffold_copy.findall('n')}
        id_to_bond = {b.attrib['id']: b for b in scaffold_copy.findall('b')}

        nodes_to_remove = []
        nodes_to_add = []

        # Replace NamedAlternativeGroup nodes with ExternalConnectionPoint dummies
        for node_id, label, ag_id, bond_ordering, valence in rgroup_nodes:
            n_elem = id_to_node[node_id]
            entries = node_map_info[node_id]
            if valence == 1:
                n_elem.attrib.clear()
                n_elem.attrib['id'] = node_id
                n_elem.attrib['NodeType'] = 'ExternalConnectionPoint'
                n_elem.attrib['ExternalConnectionNum'] = '1'
                for child in list(n_elem):
                    n_elem.remove(child)
            else:
                for i, bond_id in enumerate(bond_ordering):
                    global_mapnum, syn_id, conn_num = entries[i]
                    d = ET.Element('n')
                    d.attrib['id'] = syn_id
                    d.attrib['NodeType'] = 'ExternalConnectionPoint'
                    d.attrib['ExternalConnectionNum'] = str(conn_num)
                    nodes_to_add.append(d)
                    if bond_id in id_to_bond:
                        b = id_to_bond[bond_id]
                        if b.attrib.get('B') == node_id:
                            b.attrib['B'] = syn_id
                        elif b.attrib.get('E') == node_id:
                            b.attrib['E'] = syn_id
                nodes_to_remove.append(node_id)

        # Replace bracket inner atoms with ExternalConnectionPoint dummies
        for bg in bracket_groups:
            crossing_left = bg.get('crossing_left')
            crossing_right = bg.get('crossing_right')
            map_left = bg['map_left']
            map_right = bg['map_right']

            for crossing, mapnum in [(crossing_left, map_left), (crossing_right, map_right)]:
                if crossing is None:
                    continue
                bond_id, inner_node_id, _ = crossing
                syn_id = _fresh_id()
                # Track: inner_node_id -> syn_id (for repeat unit mol building later)
                # Add a new ExternalConnectionPoint node to the scaffold
                d = ET.Element('n')
                d.attrib['id'] = syn_id
                d.attrib['NodeType'] = 'ExternalConnectionPoint'
                # Store mapnum in a custom attrib for later CDX_NODE_ID lookup
                d.attrib['ExternalConnectionNum'] = '1'  # placeholder, we use CDX_NODE_ID to find it
                nodes_to_add.append(d)

                # Rewire the crossing bond: replace the inner atom end with the dummy
                if bond_id in id_to_bond:
                    b = id_to_bond[bond_id]
                    if b.attrib.get('B') == inner_node_id:
                        b.attrib['B'] = syn_id
                    elif b.attrib.get('E') == inner_node_id:
                        b.attrib['E'] = syn_id

                # Store the syn_id -> mapnum for post-parse labelling
                bg.setdefault('syn_id_to_mapnum', {})[syn_id] = mapnum

            # Remove ALL inner bracket atoms from scaffold (they belong to the repeat unit only)
            inner_ids = bg['inner_node_ids']
            for inner_node_id in inner_ids:
                if inner_node_id in id_to_node:
                    nodes_to_remove.append(inner_node_id)

            # Remove bonds that connect only to inner bracket atoms
            for b_elem in list(scaffold_copy.findall('b')):
                b_b = b_elem.attrib.get('B', '')
                b_e = b_elem.attrib.get('E', '')
                if b_b in inner_ids and b_e in inner_ids:
                    scaffold_copy.remove(b_elem)

        for node_id in nodes_to_remove:
            if node_id in id_to_node:
                scaffold_copy.remove(id_to_node[node_id])
        for n_elem in nodes_to_add:
            scaffold_copy.append(n_elem)

        # Build CDX_NODE_ID -> global_mapnum mapping for scaffold post-parse
        cdxid_to_mapnum = {}
        for node_id, entries in node_map_info.items():
            for global_mapnum, syn_id, conn_num in entries:
                if syn_id.isdigit():
                    cdxid_to_mapnum[int(syn_id)] = global_mapnum
                else:
                    # valence=1: syn_id == node_id (original)
                    if node_id.isdigit():
                        cdxid_to_mapnum[int(node_id)] = global_mapnum
        for bg in bracket_groups:
            for syn_id, mapnum in bg.get('syn_id_to_mapnum', {}).items():
                if syn_id.isdigit():
                    cdxid_to_mapnum[int(syn_id)] = mapnum

        # Parse scaffold mol, convert attachment atoms to * dummies
        scaffold_cdxml = CDXML_HEADER + ET.tostring(scaffold_copy, encoding='unicode') + CDXML_FOOTER
        scaffold_mols = Chem.MolsFromCDXML(scaffold_cdxml)
        if not scaffold_mols or scaffold_mols[0] is None:
            RDLogger.EnableLog('rdApp.warning')
            raise ValueError(f"Failed to parse modified scaffold from: {cdxml_path}")
        scaffold_mol = scaffold_mols[0]
        rwmol = RWMol(scaffold_mol)
        for atom in rwmol.GetAtoms():
            cdx_id = atom.GetPropsAsDict().get('CDX_NODE_ID', None)
            if cdx_id is not None and cdx_id in cdxid_to_mapnum:
                _make_dummy(rwmol, atom.GetIdx(), cdxid_to_mapnum[cdx_id])
        scaffold_mol = rwmol.GetMol()

        # --- Parse altgroup fragments ---
        rgroup_fragments = {}  # label -> [mol, ...]
        for node_id, label, ag_id, bond_ordering, valence in rgroup_nodes:
            if ag_id not in altgroups:
                continue
            conn_to_mapnum = {conn_num: global_mapnum
                              for global_mapnum, syn_id, conn_num in node_map_info[node_id]}
            frags = []
            for frag_elem in altgroups[ag_id].findall('fragment'):
                mol = _parse_fragment_mol(frag_elem, conn_to_mapnum)
                if mol is not None:
                    frags.append(mol)
            rgroup_fragments[label] = frags

        if not rgroup_fragments and not bracket_groups:
            raise ValueError(
                f"No R-group fragments or repeat units found in CDXML file: {cdxml_path}. "
                "Use ChemDraw's R-group substitution table (Rn labels with an R-group table)."
            )

        # --- Parse repeat unit mols ---
        # For each bracket group, build a mini-fragment with the inner atoms + 2 ECP nodes
        repeat_unit_mols = []  # list of (mol_with_two_dummy_atoms, map_left, map_right, repeat_range)
        for bg in bracket_groups:
            inner_ids = bg['inner_node_ids']
            map_left = bg['map_left']
            map_right = bg['map_right']
            crossing_left = bg.get('crossing_left')
            crossing_right = bg.get('crossing_right')

            if crossing_left is None:
                continue

            # Build mini-fragment XML: inner atoms + bonds between them + 2 ECP attachment nodes
            syn_left = _fresh_id()
            syn_right = _fresh_id()

            frag_elem = ET.Element('fragment')
            frag_elem.attrib['id'] = _fresh_id()

            # Add inner nodes from original scaffold
            for n in scaffold_frag_elem.findall('n'):
                if n.attrib.get('id', '') in inner_ids:
                    frag_elem.append(copy.deepcopy(n))

            # Add ECP attachment nodes
            d_left = ET.SubElement(frag_elem, 'n')
            d_left.attrib['id'] = syn_left
            d_left.attrib['NodeType'] = 'ExternalConnectionPoint'
            d_left.attrib['ExternalConnectionNum'] = '1'

            d_right = ET.SubElement(frag_elem, 'n')
            d_right.attrib['id'] = syn_right
            d_right.attrib['NodeType'] = 'ExternalConnectionPoint'
            d_right.attrib['ExternalConnectionNum'] = '2'

            # Add bonds between inner atoms
            for b in scaffold_frag_elem.findall('b'):
                b_b = b.attrib.get('B', '')
                b_e = b.attrib.get('E', '')
                if b_b in inner_ids and b_e in inner_ids:
                    frag_elem.append(copy.deepcopy(b))

            # Add crossing bonds: ECP_left -- inner_left_atom, inner_right_atom -- ECP_right
            inner_left_node_id = crossing_left[1]   # InnerAtomID for left crossing
            inner_right_node_id = crossing_right[1] if crossing_right else None

            b_left = ET.SubElement(frag_elem, 'b')
            b_left.attrib['id'] = _fresh_id()
            b_left.attrib['B'] = syn_left
            b_left.attrib['E'] = inner_left_node_id

            if inner_right_node_id:
                b_right = ET.SubElement(frag_elem, 'b')
                b_right.attrib['id'] = _fresh_id()
                b_right.attrib['B'] = inner_right_node_id
                b_right.attrib['E'] = syn_right

            frag_str = ET.tostring(frag_elem, encoding='unicode')
            mols = Chem.MolsFromCDXML(CDXML_HEADER + frag_str + CDXML_FOOTER)
            if not mols or mols[0] is None:
                continue

            rwmol = RWMol(mols[0])
            for atom in rwmol.GetAtoms():
                cdx_id = atom.GetPropsAsDict().get('CDX_NODE_ID', None)
                if cdx_id is not None:
                    if cdx_id == int(syn_left):
                        _make_dummy(rwmol, atom.GetIdx(), map_left)
                    elif cdx_id == int(syn_right):
                        _make_dummy(rwmol, atom.GetIdx(), map_right)
            repeat_mol = rwmol.GetMol()
            repeat_unit_mols.append((repeat_mol, map_left, map_right, bg['repeat_range']))

        # --- Build expanded scaffold variants (one per repeat count combo) ---
        # If there are repeat units, we need to enumerate repeat counts and
        # chain the repeat unit n times, then zip with the scaffold.
        # Strategy: produce one scaffold_mol variant per repeat count combination.
        if repeat_unit_mols:
            scaffold_variants = []
            repeat_count_lists = [ru[3] for ru in repeat_unit_mols]
            for count_combo in itertools.product(*repeat_count_lists):
                variant_scaffold = scaffold_mol
                for (repeat_mol, map_left, map_right, _), count in zip(repeat_unit_mols, count_combo):
                    # Chain `count` copies: left_dummy -- [unit]^n -- right_dummy
                    # Build by iteratively zipping copies via intermediate map numbers
                    chained = self._chain_repeat_unit(
                        repeat_mol, count, map_left, map_right,
                        map_counter, _fresh_id
                    )
                    if chained is None:
                        variant_scaffold = None
                        break
                    variant_scaffold = Chem.CombineMols(variant_scaffold, chained)
                if variant_scaffold is not None:
                    scaffold_variants.append((variant_scaffold, count_combo))
        else:
            scaffold_variants = [(scaffold_mol, ())]

        # --- Enumerate all R-group + repeat count combinations ---
        frag_lists = [rgroup_fragments[lbl] for lbl in label_order]
        final_compounds = []
        failed_count = 0

        for scaffold_variant, count_combo in scaffold_variants:
            for rgroup_combo in (itertools.product(*frag_lists) if frag_lists else [()]):
                try:
                    result = _assemble(scaffold_variant, list(rgroup_combo))
                    smiles = Chem.MolToSmiles(result)
                    branching = {lbl: _frag_display_smiles(frag)
                                 for lbl, frag in zip(label_order, rgroup_combo)}
                    if count_combo:
                        for i, (_, _, _, repeat_range) in enumerate(repeat_unit_mols):
                            branching[f'repeat_unit_{i+1}_count'] = str(count_combo[i])
                    final_compounds.append({'smiles': smiles, 'branching': branching})
                except Exception:
                    failed_count += 1

        if not final_compounds:
            raise ValueError(
                f"All {failed_count} enumeration combinations failed to assemble. "
                "Check that R-group labels and repeat unit brackets are correctly drawn in ChemDraw."
            )

        if failed_count > 0:
            import warnings
            warnings.warn(f"{failed_count} combination(s) failed to assemble and were skipped.")

        # Generate compound IDs
        base_name = os.path.splitext(os.path.basename(self.library_cdxml))[0]
        num_compounds = len(final_compounds)
        characters = 4
        if num_compounds > 9: characters = 3
        if num_compounds > 99: characters = 2
        if num_compounds > 999: characters = 1
        if num_compounds > 99999: characters = 0

        compound_ids = []
        for u_l_n in range(num_compounds):
            u_l_n_str = str(u_l_n)
            n0 = 5 - characters - len(u_l_n_str)
            zeros_str = '0' * n0
            compound_name = base_name if num_compounds == 1 else base_name[:characters] + zeros_str + u_l_n_str
            compound_ids.append(compound_name)

        RDLogger.EnableLog('rdApp.warning')
        self.expanded_compounds = final_compounds
        self.compound_ids = compound_ids

    @staticmethod
    def _chain_repeat_unit(repeat_mol, count, map_left, map_right, map_counter, fresh_id_fn):
        """
        Build a molecule representing `count` copies of repeat_mol chained together.

        repeat_mol has [*:map_left] and [*:map_right] attachment points.
        Returns a mol with [*:map_left] on the left end and [*:map_right] on the right end,
        ready to be zipped into the scaffold.

        For count=1: returns repeat_mol as-is.
        For count>1: chains copies by zipping right end of copy i to left end of copy i+1,
                     using fresh intermediate map numbers.
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import rdmolops, RWMol

            if count <= 0:
                return None
            if count == 1:
                return repeat_mol

            params = rdmolops.MolzipParams()
            params.label = rdmolops.MolzipLabel.AtomMapNumber

            # Start with first copy, keep map_left on left, use a temp map on right
            result = repeat_mol

            for i in range(1, count):
                # Assign a fresh intermediate map number for the junction
                junction_map = map_counter[0]
                map_counter[0] += 1

                # On current result: relabel map_right -> junction_map
                rwresult = RWMol(result)
                for atom in rwresult.GetAtoms():
                    if atom.GetAtomMapNum() == map_right:
                        atom.SetAtomMapNum(junction_map)
                result = rwresult.GetMol()

                # On next copy: relabel map_left -> junction_map
                rwcopy = RWMol(repeat_mol)
                for atom in rwcopy.GetAtoms():
                    if atom.GetAtomMapNum() == map_left:
                        atom.SetAtomMapNum(junction_map)
                next_copy = rwcopy.GetMol()

                # Zip them together
                combined = Chem.CombineMols(result, next_copy)
                result = Chem.molzip(combined, params)
                Chem.SanitizeMol(result)

            # The final result has map_left on left and map_right on right
            return result

        except Exception:
            return None

    @staticmethod
    def _assemble_molecule(core, fragments):
        """
        Assemble a molecule from core scaffold and R-group fragments using molzip.

        Args:
            core: RDKit Mol - core scaffold with * dummy atoms (atomMapNum > 0)
            fragments: tuple of RDKit Mol - one fragment per R-group position

        Returns:
            Assembled RDKit Mol, or None on failure
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import rdmolops

            combined = core
            for frag in fragments:
                combined = Chem.CombineMols(combined, frag)

            params = rdmolops.MolzipParams()
            params.label = rdmolops.MolzipLabel.AtomMapNumber
            assembled = Chem.molzip(combined, params)
            Chem.SanitizeMol(assembled)
            return assembled
        except Exception:
            return None

    def _generate_csv_script(self, compounds_data_json: str, source_label: str) -> str:
        """Generate inline Python script to write compounds CSV from JSON data."""
        return f"""
echo "Generating compound library from {source_label} ({len(self.expanded_compounds)} compounds)"

# Generate CSV with expanded compounds
python3 -c "
import csv
import json

# Load compound data from JSON file
with open('{compounds_data_json}', 'r') as f:
    data = json.load(f)
expanded_compounds = data['expanded_compounds']
compound_ids = data['compound_ids']

# Write CSV file with standardized format
with open('{self.compounds_csv}', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)

    # Write header with standardized columns
    header = ['id', 'format', 'smiles', 'ccd']

    # Add branching columns
    all_branch_keys = set()
    for comp in expanded_compounds:
        all_branch_keys.update(comp['branching'].keys())
    all_branch_keys = sorted(list(all_branch_keys))
    header.extend(all_branch_keys)

    writer.writerow(header)

    # Write compound data
    for i, compound_data in enumerate(expanded_compounds):
        compound_id = compound_ids[i]
        row = [
            compound_id,
            'smiles',  # format
            compound_data['smiles'],
            ''  # ccd (empty for non-covalent)
        ]

        # Add branching information
        for key in all_branch_keys:
            row.append(compound_data['branching'].get(key, ''))

        writer.writerow(row)

print(f'Generated compound library: {{len(expanded_compounds)}} compounds')
"
"""

    def generate_script(self, script_path: str) -> str:
        """
        Generate bash script for CompoundLibrary processing.

        Args:
            script_path: Path where script should be written

        Returns:
            Script content as string
        """
        script_content = "#!/bin/bash\n"
        script_content += "# CompoundLibrary processing script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += "echo \"Processing compound library\"\n"

        # Create output directories
        script_content += f"""
# Create output directories
mkdir -p "{self.output_folder}"
"""
        if self.covalent:
            script_content += f'mkdir -p "{self.covalent_folder}"\n'

        if self.library_csv:
            # Process existing CSV file
            script_content += f"""
echo "Loading compound library from CSV: {os.path.basename(self.library_csv)}"
cp "{self.library_csv}" "{self.compounds_csv}"
"""
        elif self.library_dict or self.library_cdxml:
            # Generate from dictionary or CDXML (both use pre-expanded compounds)
            os.makedirs(self.output_folder, exist_ok=True)

            if self.library_dict:
                with open(self.library_dict_json, 'w') as f:
                    json.dump(self.library_dict, f, indent=2)
                source_label = "dictionary"
            else:
                source_label = f"CDXML ({os.path.basename(self.library_cdxml)})"

            compounds_data_json = os.path.join(self.output_folder, "compounds_data.json")
            with open(compounds_data_json, 'w') as f:
                json.dump({'expanded_compounds': self.expanded_compounds, 'compound_ids': self.compound_ids}, f, indent=2)

            script_content += self._generate_csv_script(compounds_data_json, source_label)

        # Generate covalent ligand files if requested
        if self.covalent:
            script_content += f"""
echo "Generating covalent ligand CCD/PKL files"
# Note: This requires Boltz2 cache folder to be available
# The script will generate CCD/PKL files based on the compounds CSV
if [ -n "$BOLTZ_CACHE_FOLDER" ] && [ -d "$BOLTZ_CACHE_FOLDER" ]; then
    python3 "{self.covalent_generation_py}" \\
        "$BOLTZ_CACHE_FOLDER" \\
        "base_config.txt" \\
        "{self.compounds_csv}" \\
        "{self.covalent_folder}" \\
        "{self.covalent_folder}" \\
        "{self.covalent_compounds_csv}" \\
        "{self.conformer_method}"
else
    echo "Warning: BOLTZ_CACHE_FOLDER not set or not found. Skipping covalent file generation."
    echo "To generate covalent files, ensure Boltz2 environment is available and BOLTZ_CACHE_FOLDER is set."
fi
"""

        # Generate summary
        library_type = "Dictionary" if self.library_dict else ("CDXML" if self.library_cdxml else "CSV file")
        primary_key_str = self.primary_key if self.primary_key else "None"
        covalent_str = str(self.covalent)
        conformer_method_str = self.conformer_method
        compounds_csv_basename = os.path.basename(self.compounds_csv)
        is_covalent = self.covalent

        script_content += f"""
echo "Generating library summary"
python3 -c "
import pandas as pd
import os

# Read library file
df = pd.read_csv('{self.compounds_csv}')
compound_count = len(df)

# Write summary
with open('{self.summary_file}', 'w') as f:
    f.write('Compound Library Summary\\n')
    f.write('========================\\n')
    f.write(f'Library type: {library_type}\\n')
    if '{primary_key_str}' != 'None':
        f.write(f'Primary key: {primary_key_str}\\n')
    f.write(f'Total compounds: {{compound_count}}\\n')
    f.write(f'Covalent ligands: {covalent_str}\\n')
    f.write(f'Conformer method: {conformer_method_str}\\n')
    f.write(f'Output file: {compounds_csv_basename}\\n')
    if {is_covalent}:
        f.write(f'Covalent library folder: covalent_library/\\n')

print(f'Library processed: {{compound_count}} compounds')
print(f'Output: {self.compounds_csv}')
"
"""

        script_content += self.generate_completion_check_footer()

        return script_content

    def get_output_files(self) -> Dict[str, Any]:
        """
        Get expected output files after CompoundLibrary processing.

        Returns:
            Dictionary with DataStream objects and tables
        """
        # Generate predicted compound IDs if not already done
        if not self.compound_ids and self.library_dict:
            if self.covalent:
                # For covalent: use library keys as-is (they're already expanded)
                self.compound_ids = list(self.library_dict.keys())
            else:
                # For simple library: use library keys as-is
                self.compound_ids = list(self.library_dict.keys())

        # Build tables with rich metadata
        tables = {}
        columns = ["id", "format", "smiles", "ccd"]
        if self.library_dict or self.library_cdxml:
            # Add branching columns
            all_branch_keys = set()
            for comp in self.expanded_compounds:
                all_branch_keys.update(comp['branching'].keys())
            columns.extend(sorted(list(all_branch_keys)))

        tables["compounds"] = TableInfo(
            name="compounds",
            path=self.compounds_csv,
            columns=columns,
            description="Generated compound library with SMILES and metadata",
            count=len(self.expanded_compounds) if self.expanded_compounds else 0
        )

        if self.covalent and self.covalent_compounds_csv:
            tables["covalent_compounds"] = TableInfo(
                name="covalent_compounds",
                path=self.covalent_compounds_csv,
                columns=["id", "format", "smiles", "ccd"],
                description="Covalent compound library with CCD identifiers"
            )

        # Create compounds DataStream
        compounds = DataStream(
            name="compounds",
            ids=self.compound_ids,
            files=[],  # Value-based format - data is in map_table, not individual files
            map_table=self.compounds_csv,
            format="csv"
        )

        return {
            "compounds": compounds,
            "tables": tables,
            "output_folder": self.output_folder
        }

    def get_config_display(self) -> List[str]:
        """Get configuration display lines for pipeline output."""
        config_lines = super().get_config_display()

        if isinstance(self.library, dict):
            config_lines.append(f"LIBRARY: Dictionary ({len(self.library)} keys)")
            if self.primary_key:
                config_lines.append(f"PRIMARY KEY: {self.primary_key}")
        elif self.library_cdxml:
            config_lines.append(f"LIBRARY: CDXML ({os.path.basename(self.library_cdxml)})")
            # Show R-group positions
            if self.expanded_compounds:
                all_positions = set()
                for comp in self.expanded_compounds:
                    all_positions.update(comp['branching'].keys())
                if all_positions:
                    config_lines.append(f"R-GROUP POSITIONS: {', '.join(sorted(all_positions))}")
        else:
            config_lines.append(f"LIBRARY: {os.path.basename(self.library)}")

        config_lines.extend([
            f"COVALENT LIGANDS: {self.covalent}",
            f"VALIDATE SMILES: {self.validate_smiles}"
        ])

        if self.covalent:
            config_lines.append(f"CONFORMER METHOD: {self.conformer_method}")

        return config_lines

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including CompoundLibrary-specific parameters."""
        base_dict = super().to_dict()
        if self.library_cdxml:
            library_type = "cdxml"
        elif self.library_dict:
            library_type = "dictionary"
        else:
            library_type = "csv"

        base_dict.update({
            "compound_library_params": {
                "library": self.library if isinstance(self.library, str) else "<dictionary>",
                "library_type": library_type,
                "primary_key": self.primary_key,
                "covalent": self.covalent,
                "validate_smiles": self.validate_smiles,
                "conformer_method": self.conformer_method,
                "num_compounds": len(self.expanded_compounds) if self.expanded_compounds else 0
            }
        })
        return base_dict
