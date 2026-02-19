#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Runtime helper script for Ligand tool.

Fetches small molecule ligands with priority-based lookup: local_folder -> Ligands/ -> RCSB/PubChem download.
Downloads SDF files and converts them to PDB format with proper atom numbering.
Supports RCSB (CCD codes), PubChem (name, CID, CAS), and direct SMILES input.
"""

import os
import sys
import argparse
import json
import re
import pandas as pd
from datetime import datetime
from typing import Dict, List, Any, Tuple, Optional
from pathlib import Path


def detect_lookup_type(lookup: str) -> str:
    """
    Detect the type of lookup value.

    Args:
        lookup: The lookup value to analyze

    Returns:
        "ccd" (RCSB), "cid" (PubChem), "cas" (PubChem), or "name" (PubChem)
    """
    # CCD codes: 1-3 uppercase alphanumeric characters
    if re.match(r'^[A-Z0-9]{1,3}$', lookup.upper()) and not lookup.isdigit():
        return "ccd"

    # PubChem CID: purely numeric
    if lookup.isdigit():
        return "cid"

    # CAS number: XX-XX-X format (digits-digits-digit)
    if re.match(r'^\d+-\d+-\d$', lookup):
        return "cas"

    # Default: compound name
    return "name"


def find_local_ligand(lookup: str, local_folder: str,
                      repo_ligands_folder: str, output_format: str = "pdb") -> Tuple[Optional[str], Optional[str]]:
    """
    Find ligand file locally.

    Priority: local_folder (if given) -> repo_ligands_folder -> None

    Args:
        lookup: Lookup value (used as filename stem)
        local_folder: Custom local folder (can be None)
        repo_ligands_folder: Repository Ligands folder
        output_format: "pdb" or "cif" - the format to look for

    Returns:
        Tuple of (path to local structure file, path to local CSV file) or (None, None) if not found
    """
    search_locations = []
    ext = output_format  # "pdb" or "cif"

    if local_folder:
        search_locations.append(local_folder)
    search_locations.append(repo_ligands_folder)

    for location in search_locations:
        file_candidate = os.path.join(location, f"{lookup}.{ext}")
        csv_candidate = os.path.join(location, f"{lookup}.csv")
        if os.path.exists(file_candidate):
            print(f"Found {lookup} locally: {file_candidate}")
            csv_path = csv_candidate if os.path.exists(csv_candidate) else None
            return file_candidate, csv_path

    return None, None


def load_local_metadata(csv_path: str) -> Dict[str, Any]:
    """
    Load metadata from local CSV file.

    Args:
        csv_path: Path to the CSV metadata file

    Returns:
        Dictionary with metadata fields
    """
    try:
        df = pd.read_csv(csv_path)
        if len(df) > 0:
            row = df.iloc[0].to_dict()
            # Convert NaN to empty strings
            return {k: ('' if pd.isna(v) else v) for k, v in row.items()}
    except Exception as e:
        print(f"  Warning: Could not load metadata from {csv_path}: {e}")
    return {}


def copy_local_ligand(lookup: str, custom_id: str, residue_code: str,
                      source_path: str, source_csv_path: Optional[str],
                      output_folder: str, output_format: str = "pdb") -> Tuple[bool, str, Dict[str, Any]]:
    """
    Copy local ligand file to output folder with residue renaming.

    Args:
        lookup: Original lookup value
        custom_id: Custom ID for output filename
        residue_code: Residue code to use
        source_path: Path to local ligand file (PDB or CIF)
        source_csv_path: Path to local ligand CSV metadata file (can be None)
        output_folder: Directory to save the ligand
        output_format: "pdb" or "cif"

    Returns:
        Tuple of (success: bool, file_path: str, metadata: dict)
    """
    try:
        ext = output_format

        # Read local file
        with open(source_path, 'r') as f:
            content = f.read()

        # Rename residue code in the file content
        if output_format == "pdb":
            content = rename_residue_chain_A(content, residue_code)
        elif output_format == "cif":
            content = rename_residue_in_cif(content, residue_code)

        filename = f"{custom_id}.{ext}"
        output_path = os.path.join(output_folder, filename)

        with open(output_path, 'w') as f:
            f.write(content)

        file_size = os.path.getsize(output_path)

        # Load metadata from CSV if available
        metadata = {}
        if source_csv_path:
            metadata = load_local_metadata(source_csv_path)

        # Update/override certain fields
        metadata.update({
            "file_size": file_size,
            "source": metadata.get('source', 'local'),
            "source_path": source_path,
        })

        print(f"Successfully copied {lookup} as {custom_id}.{ext}: {file_size} bytes (from local)")
        return True, output_path, metadata

    except Exception as e:
        error_msg = f"Error copying local file {lookup}: {str(e)}"
        print(f"Error: {error_msg}")
        metadata = {
            "error_message": error_msg,
            "source": "local_copy_failed",
            "attempted_path": source_path
        }
        return False, "", metadata


def convert_smiles_to_pdb_rdkit(smiles: str, residue_code: str) -> Optional[str]:
    """
    Convert SMILES to PDB format using RDKit.

    IMPORTANT: This preserves RDKit's atom naming convention, which is required
    for compatibility with tools like RFdiffusion3 that internally use RDKit
    to regenerate conformers and expect matching atom names.

    Args:
        smiles: SMILES string of the molecule
        residue_code: 3-letter residue code to use in PDB (e.g., "AMX")

    Returns:
        PDB content as string, or None if conversion fails
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"  Error: Could not parse SMILES: {smiles[:50]}...")
            return None

        # Add hydrogens for proper 3D geometry
        mol = Chem.AddHs(mol)

        # Generate 3D coordinates using ETKDG (recommended method)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42  # For reproducibility
        result = AllChem.EmbedMolecule(mol, params)

        if result == -1:
            # ETKDG failed, try with random coordinates
            print(f"  Warning: ETKDG embedding failed, trying random coordinates...")
            params.useRandomCoords = True
            result = AllChem.EmbedMolecule(mol, params)

            if result == -1:
                print(f"  Error: Could not generate 3D coordinates")
                return None

        # Optimize geometry with MMFF force field
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
        except Exception as e:
            print(f"  Warning: MMFF optimization failed: {e}, using unoptimized coordinates")

        # Remove hydrogens for cleaner output (RFD3 will add them back)
        mol = Chem.RemoveHs(mol)

        # Generate PDB block - RDKit uses its own atom naming
        pdb_block = Chem.MolToPDBBlock(mol)

        if not pdb_block:
            print(f"  Error: Could not generate PDB block")
            return None

        # Only rename the residue code, preserve RDKit's atom naming
        pdb_content = rename_residue_chain_A(pdb_block, residue_code)

        return pdb_content

    except ImportError:
        print("  Error: RDKit not available. Install with: pip install rdkit")
        return None
    except Exception as e:
        print(f"  Error converting SMILES to PDB with RDKit: {str(e)}")
        return None


def convert_smiles_to_cif_rdkit(smiles: str, residue_code: str) -> Optional[str]:
    """
    Convert SMILES to mmCIF format using RDKit with explicit bond orders.

    Generates a proper mmCIF structure file with:
    - _atom_site category for atomic coordinates (required by biotite/foundry)
    - _chem_comp_bond category for explicit bond orders (for RFdiffusion3)

    This format is required for tools like RFdiffusion3 that use biotite to parse
    structures and need explicit bond order information.

    Args:
        smiles: SMILES string of the molecule
        residue_code: Residue code to use (e.g., "LIG", "L_G")

    Returns:
        mmCIF content as string, or None if conversion fails
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"  Error: Could not parse SMILES: {smiles[:50]}...")
            return None

        # Add hydrogens for proper 3D geometry
        mol = Chem.AddHs(mol)

        # Generate 3D coordinates using ETKDG (recommended method)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42  # For reproducibility
        result = AllChem.EmbedMolecule(mol, params)

        if result == -1:
            # ETKDG failed, try with random coordinates
            print(f"  Warning: ETKDG embedding failed, trying random coordinates...")
            params.useRandomCoords = True
            result = AllChem.EmbedMolecule(mol, params)

            if result == -1:
                print(f"  Error: Could not generate 3D coordinates")
                return None

        # Optimize geometry with MMFF force field
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
        except Exception as e:
            print(f"  Warning: MMFF optimization failed: {e}, using unoptimized coordinates")

        # Remove hydrogens for cleaner output
        mol = Chem.RemoveHs(mol)

        # Get conformer for coordinates
        conf = mol.GetConformer()

        # Build mmCIF content with proper _atom_site category
        cif_lines = []
        cif_lines.append("data_" + residue_code)
        cif_lines.append("#")

        # _entry section
        cif_lines.append(f"_entry.id {residue_code}")
        cif_lines.append("#")

        # _chem_comp section (component metadata)
        cif_lines.append(f"_chem_comp.id {residue_code}")
        cif_lines.append("_chem_comp.name 'Small molecule ligand'")
        cif_lines.append("_chem_comp.type NON-POLYMER")
        cif_lines.append("#")

        # Track atom names for bond section
        atom_names = []
        element_counts = {}

        # Pre-generate atom names
        for atom in mol.GetAtoms():
            element = atom.GetSymbol()
            element_counts[element] = element_counts.get(element, 0) + 1
            atom_name = f"{element}{element_counts[element]}"
            atom_names.append(atom_name)

        # _atom_site section (required by biotite for structure parsing)
        # Include all standard mmCIF fields that biotite/foundry may require
        cif_lines.append("loop_")
        cif_lines.append("_atom_site.group_PDB")
        cif_lines.append("_atom_site.id")
        cif_lines.append("_atom_site.type_symbol")
        cif_lines.append("_atom_site.label_atom_id")
        cif_lines.append("_atom_site.label_alt_id")
        cif_lines.append("_atom_site.label_comp_id")
        cif_lines.append("_atom_site.label_asym_id")
        cif_lines.append("_atom_site.label_entity_id")
        cif_lines.append("_atom_site.label_seq_id")
        cif_lines.append("_atom_site.pdbx_PDB_ins_code")
        cif_lines.append("_atom_site.pdbx_PDB_model_num")
        cif_lines.append("_atom_site.Cartn_x")
        cif_lines.append("_atom_site.Cartn_y")
        cif_lines.append("_atom_site.Cartn_z")
        cif_lines.append("_atom_site.occupancy")
        cif_lines.append("_atom_site.B_iso_or_equiv")
        cif_lines.append("_atom_site.auth_atom_id")
        cif_lines.append("_atom_site.auth_comp_id")
        cif_lines.append("_atom_site.auth_asym_id")
        cif_lines.append("_atom_site.auth_seq_id")
        cif_lines.append("_atom_site.pdbx_formal_charge")

        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            element = atom.GetSymbol()
            charge = atom.GetFormalCharge()
            pos = conf.GetAtomPosition(idx)
            atom_name = atom_names[idx]
            atom_id = idx + 1  # 1-indexed

            # Format: group_PDB id type_symbol label_atom_id label_alt_id label_comp_id
            #         label_asym_id label_entity_id label_seq_id pdbx_PDB_ins_code
            #         pdbx_PDB_model_num Cartn_x Cartn_y Cartn_z occupancy B_iso
            #         auth_atom_id auth_comp_id auth_asym_id auth_seq_id pdbx_formal_charge
            # Use . for empty/null values, HETATM for ligand atoms
            cif_lines.append(
                f"HETATM {atom_id} {element} {atom_name} . {residue_code} A 1 . . 1 "
                f"{pos.x:.3f} {pos.y:.3f} {pos.z:.3f} 1.00 0.00 "
                f"{atom_name} {residue_code} A 1 {charge}"
            )

        cif_lines.append("#")

        # Note: _chem_comp_bond is intentionally omitted as it's for CCD files,
        # not structure files. Foundry/biotite infers bonds from coordinates.

        return "\n".join(cif_lines)

    except ImportError:
        print("  Error: RDKit not available. Install with: pip install rdkit")
        return None
    except Exception as e:
        print(f"  Error converting SMILES to CIF with RDKit: {str(e)}")
        import traceback
        traceback.print_exc()
        return None


def rename_residue_in_cif(cif_content: str, new_residue_code: str) -> str:
    """
    Rename the residue code in CIF content.

    Updates all occurrences of the residue code in the CIF file:
    - data_<code> header
    - _entry.id
    - _chem_comp.id
    - _atom_site.label_comp_id column
    - _atom_site.auth_comp_id column

    Args:
        cif_content: CIF file content as string
        new_residue_code: New residue code to use (e.g., "LIG", "G")

    Returns:
        CIF content with residue code renamed
    """
    lines = cif_content.split('\n')
    output_lines = []

    # First pass: find the original residue code from data_ line or _entry.id
    original_code = None
    for line in lines:
        if line.startswith('data_'):
            original_code = line[5:].strip()
            break
        if line.startswith('_entry.id '):
            original_code = line.split()[1].strip()
            break

    if original_code is None:
        # Can't determine original code, return as-is
        return cif_content

    # If codes are the same, no changes needed
    if original_code == new_residue_code:
        return cif_content

    # Track if we're in the _atom_site loop and which columns have comp_id
    in_atom_site_loop = False
    atom_site_columns = []
    label_comp_id_idx = None
    auth_comp_id_idx = None

    for line in lines:
        # Check for entering _atom_site loop
        if line.strip() == 'loop_':
            in_atom_site_loop = False
            atom_site_columns = []
            label_comp_id_idx = None
            auth_comp_id_idx = None
            output_lines.append(line)
            continue

        # Track _atom_site column definitions
        if line.startswith('_atom_site.'):
            in_atom_site_loop = True
            col_name = line.strip()
            atom_site_columns.append(col_name)
            if col_name == '_atom_site.label_comp_id':
                label_comp_id_idx = len(atom_site_columns) - 1
            elif col_name == '_atom_site.auth_comp_id':
                auth_comp_id_idx = len(atom_site_columns) - 1
            output_lines.append(line)
            continue

        # Handle data_ line
        if line.startswith('data_'):
            output_lines.append(f'data_{new_residue_code}')
            continue

        # Handle _entry.id line
        if line.startswith('_entry.id '):
            output_lines.append(f'_entry.id {new_residue_code}')
            continue

        # Handle _chem_comp.id line
        if line.startswith('_chem_comp.id '):
            output_lines.append(f'_chem_comp.id {new_residue_code}')
            continue

        # Handle _atom_site data rows (HETATM/ATOM lines)
        if in_atom_site_loop and (line.startswith('HETATM') or line.startswith('ATOM')):
            # Split the line into fields
            parts = line.split()
            if len(parts) >= len(atom_site_columns):
                # Replace residue codes in the appropriate columns
                if label_comp_id_idx is not None and label_comp_id_idx < len(parts):
                    if parts[label_comp_id_idx] == original_code:
                        parts[label_comp_id_idx] = new_residue_code
                if auth_comp_id_idx is not None and auth_comp_id_idx < len(parts):
                    if parts[auth_comp_id_idx] == original_code:
                        parts[auth_comp_id_idx] = new_residue_code
                output_lines.append(' '.join(parts))
                continue

        # Handle # separator - exit atom_site loop
        if line.strip() == '#':
            in_atom_site_loop = False

        output_lines.append(line)

    return '\n'.join(output_lines)


def rename_residue_chain_A(pdb_content: str, residue_code: str) -> str:
    """
    Rename only the residue code in PDB content, preserving atom names.
    Also sets chain to chain A.

    This is critical for RFdiffusion3 compatibility - it internally uses RDKit
    to regenerate conformers and expects the original RDKit atom naming.

    Args:
        pdb_content: PDB content (e.g., from RDKit)
        residue_code: 3-letter residue code to use (e.g., "AMX")

    Returns:
        PDB content with residue renamed but atom names preserved
    """
    lines = pdb_content.split('\n')
    output_lines = []

    # Ensure residue code is uppercase and exactly 3 chars
    res_code = residue_code.upper()[:3].ljust(3)

    for line in lines:
        if line.startswith(('HETATM', 'ATOM')):
            # PDB format: columns 17-19 are residue name (0-indexed: 17:20)
            # Preserve everything else including atom names (columns 12-15)
            if len(line) >= 22:
                new_line = line[:17] + res_code + line[20] + "A" + line[22:]
                output_lines.append(new_line)
            else:
                output_lines.append(line)
        elif line.startswith('CONECT') or line.startswith('END'):
            output_lines.append(line)
        # Skip other lines (HEADER, COMPND, etc. from RDKit)

    # Ensure END record
    if output_lines and not output_lines[-1].startswith('END'):
        output_lines.append('END')

    return '\n'.join(output_lines)


def _extract_element_from_line(line: str) -> str:
    """
    Extract element symbol from a PDB ATOM/HETATM line.

    Tries multiple sources: element column (76-77), atom name, or first alpha chars.

    Args:
        line: PDB ATOM/HETATM line

    Returns:
        Element symbol (1-2 chars, uppercase)
    """
    # Try element column (76-77) first - most reliable
    if len(line) >= 78:
        element = line[76:78].strip()
        if element and element[0].isalpha():
            return element.upper()

    # Try to extract from atom name (columns 12-15)
    if len(line) >= 16:
        atom_name = line[12:16].strip()
        if atom_name:
            # Extract leading alpha characters
            element = ''
            for char in atom_name:
                if char.isalpha():
                    element += char
                else:
                    break
            if element:
                # Handle 2-char elements (e.g., CL, BR, FE)
                element = element.upper()
                if len(element) >= 2 and element[:2] in ['CL', 'BR', 'FE', 'ZN', 'MG', 'CA', 'NA', 'MN', 'CU', 'CO', 'NI', 'SE']:
                    return element[:2]
                return element[0]

    return 'X'  # Unknown element


def normalize_pdb_content(pdb_content: str, residue_code: str) -> str:
    """
    Normalize PDB content with proper atom numbering.

    Args:
        pdb_content: Raw PDB content from OpenBabel conversion or local file
        residue_code: 3-letter residue code to use (e.g., "LIG")

    Returns:
        Normalized PDB content with:
        - Sequential atom serial numbers (1, 2, 3, ...)
        - Properly numbered atom names (C1, C2, N1, O1, ...)
        - Chain ID = 'A'
        - Residue number = 1
        - Residue name = residue_code (3-letter, uppercase)
        - Element symbol in columns 76-77
        - Renumbered CONECT records
    """
    lines = pdb_content.split('\n')
    atom_lines = []
    conect_lines = []
    other_lines = []
    end_line = None

    # Ensure residue code is uppercase and max 3 chars
    res_code = residue_code.upper()[:3].ljust(3)

    # First pass: categorize lines and build serial number mapping
    old_to_new_serial = {}
    new_serial = 0

    for line in lines:
        if line.startswith(('HETATM', 'ATOM')):
            new_serial += 1
            # Extract old serial number (columns 6-10)
            try:
                old_serial = int(line[6:11].strip())
                old_to_new_serial[old_serial] = new_serial
            except (ValueError, IndexError):
                pass
            atom_lines.append(line)
        elif line.startswith('CONECT'):
            conect_lines.append(line)
        elif line.startswith('END'):
            end_line = line
        else:
            other_lines.append(line)

    # Second pass: process atoms with new serial numbers and atom names
    output_lines = []
    element_counts = {}
    atom_serial = 0

    for line in atom_lines:
        atom_serial += 1

        # Extract element symbol
        element = _extract_element_from_line(line)

        # Generate numbered atom name
        element_counts[element] = element_counts.get(element, 0) + 1
        atom_num = element_counts[element]

        # Format atom name (4 chars, columns 12-15)
        atom_name_str = f"{element}{atom_num}"
        if len(element) == 1:
            atom_name = f" {atom_name_str:<3}"
        else:
            atom_name = f"{atom_name_str:<4}"

        record = 'HETATM'
        serial_str = f"{atom_serial:5d}"
        alt_loc = ' '
        chain = 'A'
        res_seq = '   1'
        icode = ' '

        # Get coordinates and other data (columns 27-66)
        if len(line) >= 54:
            coords = line[27:66]
        elif len(line) >= 27:
            coords = line[27:].ljust(39)
        else:
            coords = ' ' * 39

        coords = coords.ljust(39)
        element_col = f"{element:>2}"

        new_line = f"{record}{serial_str} {atom_name}{alt_loc}{res_code} {chain}{res_seq}{icode}{coords}          {element_col}"
        output_lines.append(new_line)

    # Third pass: renumber CONECT records
    for line in conect_lines:
        # CONECT records have atom serial numbers in columns 6-10, 11-15, 16-20, 21-25, 26-30
        # Each field is 5 characters wide
        try:
            new_conect = "CONECT"
            # Parse all serial numbers from the CONECT line
            i = 6
            while i + 5 <= len(line):
                field = line[i:i+5].strip()
                if field:
                    try:
                        old_num = int(field)
                        new_num = old_to_new_serial.get(old_num, old_num)
                        new_conect += f"{new_num:5d}"
                    except ValueError:
                        break
                else:
                    break
                i += 5

            if len(new_conect) > 6:  # Has at least one connection
                output_lines.append(new_conect)
        except Exception:
            # If parsing fails, skip this CONECT record
            continue

    # Add END record
    if end_line:
        output_lines.append(end_line)
    else:
        output_lines.append("END")

    return '\n'.join(output_lines)


def fetch_smiles_from_rcsb(ligand_code: str) -> Optional[str]:
    """
    Fetch SMILES string for a ligand from RCSB REST API.

    Args:
        ligand_code: 3-letter ligand code (e.g., 'ATP', 'GDP')

    Returns:
        SMILES string or None if not found
    """
    try:
        import requests

        # RCSB REST API endpoint for ligand info
        url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_code}"

        headers = {
            'User-Agent': 'BioPipelines-Ligand/1.0 (https://github.com/locbp-uzh/biopipelines)'
        }

        response = requests.get(url, headers=headers, timeout=10)
        response.raise_for_status()

        data = response.json()

        # SMILES is in pdbx_chem_comp_descriptor array
        # Each entry has: comp_id, descriptor, program, program_version, type
        # Types include: SMILES, SMILES_CANONICAL, InChI, InChIKey
        if 'pdbx_chem_comp_descriptor' in data:
            descriptors = data['pdbx_chem_comp_descriptor']
            if isinstance(descriptors, list):
                # Prefer SMILES_CANONICAL, fall back to SMILES
                for preferred_type in ['SMILES_CANONICAL', 'SMILES']:
                    for desc in descriptors:
                        if desc.get('type') == preferred_type and desc.get('descriptor'):
                            return desc['descriptor']

        return None

    except Exception as e:
        print(f"  Warning: Could not fetch SMILES for {ligand_code} from RCSB: {str(e)}")
        return None


def fetch_smiles_from_pubchem(lookup: str, lookup_type: str) -> Optional[str]:
    """
    Fetch SMILES string for a compound from PubChem.

    Args:
        lookup: The lookup value (name, CID, or CAS)
        lookup_type: "cid", "cas", or "name"

    Returns:
        SMILES string or None if not found
    """
    try:
        import requests

        base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        headers = {'User-Agent': 'BioPipelines-Ligand/1.0 (https://github.com/locbp-uzh/biopipelines)'}

        # Resolve to CID first
        if lookup_type == "cid":
            cid = int(lookup)
        else:
            search_url = f"{base_url}/compound/name/{lookup}/cids/JSON"
            response = requests.get(search_url, headers=headers, timeout=30)
            if response.status_code == 404:
                return None
            response.raise_for_status()
            data = response.json()
            if 'IdentifierList' not in data or not data['IdentifierList'].get('CID'):
                return None
            cid = data['IdentifierList']['CID'][0]

        # Get SMILES - PubChem may return different field names
        props_url = f"{base_url}/compound/cid/{cid}/property/CanonicalSMILES/JSON"
        response = requests.get(props_url, headers=headers, timeout=30)
        response.raise_for_status()
        props = response.json()['PropertyTable']['Properties'][0]
        # PubChem returns ConnectivitySMILES when requesting CanonicalSMILES
        return (props.get('CanonicalSMILES') or
                props.get('ConnectivitySMILES') or
                props.get('SMILES'))

    except Exception as e:
        print(f"  Warning: Could not fetch SMILES for {lookup} from PubChem: {str(e)}")
        return None


def fetch_properties_from_rcsb(ligand_code: str) -> Dict[str, Any]:
    """
    Fetch compound properties from RCSB REST API.

    Args:
        ligand_code: 3-letter ligand code (e.g., 'ATP', 'GDP')

    Returns:
        Dictionary with smiles, name, formula fields
    """
    try:
        import requests

        url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_code}"
        headers = {
            'User-Agent': 'BioPipelines-Ligand/1.0 (https://github.com/locbp-uzh/biopipelines)'
        }

        response = requests.get(url, headers=headers, timeout=10)
        response.raise_for_status()

        data = response.json()
        result = {'ccd': ligand_code}

        if 'chem_comp' in data:
            chem_comp = data['chem_comp']
            result['name'] = chem_comp.get('name', '')
            result['formula'] = chem_comp.get('formula', '')

        # SMILES is in pdbx_chem_comp_descriptor array
        if 'pdbx_chem_comp_descriptor' in data:
            descriptors = data['pdbx_chem_comp_descriptor']
            if isinstance(descriptors, list):
                for preferred_type in ['SMILES_CANONICAL', 'SMILES']:
                    for desc in descriptors:
                        if desc.get('type') == preferred_type and desc.get('descriptor'):
                            result['smiles'] = desc['descriptor']
                            break
                    if 'smiles' in result:
                        break

        return result

    except Exception as e:
        print(f"  Warning: Could not fetch properties for {ligand_code} from RCSB: {str(e)}")
        return {}


def convert_sdf_to_pdb(sdf_content: str, ligand_code: str) -> Optional[str]:
    """
    Convert SDF content to PDB format using OpenBabel.

    Args:
        sdf_content: SDF file content as string
        ligand_code: 3-letter ligand code for residue naming

    Returns:
        PDB content as string, or None if conversion fails
    """
    try:
        import tempfile
        import subprocess

        # Create temporary files for conversion
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sdf', delete=False) as sdf_file:
            sdf_file.write(sdf_content)
            sdf_path = sdf_file.name

        with tempfile.NamedTemporaryFile(mode='r', suffix='.pdb', delete=False) as pdb_file:
            pdb_path = pdb_file.name

        try:
            # Use obabel command line tool
            result = subprocess.run(
                ['obabel', sdf_path, '-O', pdb_path],
                capture_output=True,
                text=True,
                timeout=30
            )

            if result.returncode != 0:
                print(f"  Error: OpenBabel conversion failed: {result.stderr}")
                return None

            with open(pdb_path, 'r') as f:
                pdb_content = f.read()

            # Normalize the PDB content (sets residue name, renumbers atoms, etc.)
            return normalize_pdb_content(pdb_content, ligand_code)

        finally:
            # Clean up temporary files
            if os.path.exists(sdf_path):
                os.unlink(sdf_path)
            if os.path.exists(pdb_path):
                os.unlink(pdb_path)

    except Exception as e:
        print(f"  Error converting SDF to PDB with OpenBabel: {str(e)}")
        return None


def fetch_properties_from_pubchem(lookup: str, lookup_type: str) -> Dict[str, Any]:
    """
    Fetch compound properties from PubChem (without downloading SDF).

    Args:
        lookup: The lookup value (name, CID, or CAS)
        lookup_type: "cid", "cas", or "name"

    Returns:
        Dictionary with: cid, smiles, name, formula

    Raises:
        ValueError: If compound not found
    """
    import requests

    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    headers = {'User-Agent': 'BioPipelines-Ligand/1.0 (https://github.com/locbp-uzh/biopipelines)'}

    # Step 1: Resolve to CID
    if lookup_type == "cid":
        cid = int(lookup)
        print(f"  Using PubChem CID: {cid}")
    else:
        # Search by name or CAS
        search_url = f"{base_url}/compound/name/{lookup}/cids/JSON"
        print(f"  Searching PubChem for: {lookup}")
        response = requests.get(search_url, headers=headers, timeout=30)

        if response.status_code == 404:
            raise ValueError(f"Compound '{lookup}' not found on PubChem")

        response.raise_for_status()
        data = response.json()

        if 'IdentifierList' not in data or not data['IdentifierList'].get('CID'):
            raise ValueError(f"Compound '{lookup}' not found on PubChem")

        cid = data['IdentifierList']['CID'][0]  # Take first match
        print(f"  Resolved to PubChem CID: {cid}")

    # Step 2: Get properties (SMILES, name, formula)
    props_url = f"{base_url}/compound/cid/{cid}/property/CanonicalSMILES,IsomericSMILES,IUPACName,MolecularFormula/JSON"
    response = requests.get(props_url, headers=headers, timeout=30)
    response.raise_for_status()
    props = response.json()['PropertyTable']['Properties'][0]

    # PubChem returns different SMILES field names depending on query
    smiles = (props.get('CanonicalSMILES') or
              props.get('ConnectivitySMILES') or
              props.get('SMILES') or
              props.get('IsomericSMILES') or '')
    name = props.get('IUPACName', '')
    formula = props.get('MolecularFormula', '')

    print(f"  SMILES: {smiles[:50]}..." if len(smiles) > 50 else f"  SMILES: {smiles}")
    print(f"  Name: {name[:50]}..." if len(name) > 50 else f"  Name: {name}")

    return {
        'cid': str(cid),
        'smiles': smiles,
        'name': name,
        'formula': formula,
        'source': 'pubchem'
    }


def fetch_sdf_from_pubchem(cid: str) -> str:
    """
    Download SDF content from PubChem for a given CID.

    Args:
        cid: PubChem compound ID

    Returns:
        SDF content as string

    Raises:
        ValueError: If download fails
    """
    import requests

    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    headers = {'User-Agent': 'BioPipelines-Ligand/1.0 (https://github.com/locbp-uzh/biopipelines)'}

    # Try 3D SDF first, fall back to 2D
    sdf_url = f"{base_url}/compound/cid/{cid}/SDF?record_type=3d"
    print(f"  Downloading 3D SDF (fallback)...")
    response = requests.get(sdf_url, headers=headers, timeout=30)

    if response.status_code == 404:
        print(f"  No 3D structure available, using 2D...")
        sdf_url = f"{base_url}/compound/cid/{cid}/SDF"
        response = requests.get(sdf_url, headers=headers, timeout=30)

    response.raise_for_status()
    sdf_content = response.text

    if not sdf_content.strip():
        raise ValueError(f"Downloaded SDF file for CID {cid} is empty")

    return sdf_content


def download_from_rcsb(ligand_code: str, custom_id: str, residue_code: str,
                       output_folder: str, repo_ligands_folder: str,
                       output_format: str = "pdb") -> Tuple[bool, str, Dict[str, Any]]:
    """
    Download a single ligand from RCSB PDB.

    Uses RDKit to generate PDB or CIF from SMILES. CIF format includes explicit
    bond orders which is recommended for tools like RFdiffusion3.

    Args:
        ligand_code: 3-letter ligand code (e.g., 'ATP')
        custom_id: Custom ID for renaming the ligand
        residue_code: Code to use in output file
        output_folder: Directory to save the ligand
        repo_ligands_folder: Repository Ligands folder for caching
        output_format: "pdb" or "cif" (default: "pdb")

    Returns:
        Tuple of (success: bool, file_path: str, metadata: dict)
    """
    ligand_code = ligand_code.upper()

    # Validate ligand code format
    if not ligand_code or len(ligand_code) > 3 or not ligand_code.isalnum():
        error_msg = f"Invalid ligand code format: {ligand_code}. Must be 1-3 alphanumeric characters."
        print(f"Error: {error_msg}")
        metadata = {
            "error_message": error_msg,
            "source": "rcsb_invalid_code",
            "attempted_path": "N/A"
        }
        return False, "", metadata

    try:
        import requests

        print(f"Fetching {ligand_code} from RCSB...")

        # Get properties (including SMILES) from RCSB
        props = fetch_properties_from_rcsb(ligand_code)
        smiles = props.get('smiles', '')

        content = None
        ext = output_format  # "pdb" or "cif"

        # Convert SMILES to requested format using RDKit
        if smiles:
            print(f"  SMILES: {smiles[:50]}..." if len(smiles) > 50 else f"  SMILES: {smiles}")
            if output_format == "cif":
                print(f"  Converting SMILES to CIF using RDKit (with bond orders)...")
                content = convert_smiles_to_cif_rdkit(smiles, residue_code)
            else:
                print(f"  Converting SMILES to PDB using RDKit...")
                content = convert_smiles_to_pdb_rdkit(smiles, residue_code)

        # Fallback for PDB: OpenBabel SDF conversion (only if RDKit failed)
        if content is None and output_format == "pdb":
            print(f"  RDKit conversion failed or no SMILES, falling back to OpenBabel SDF...")
            sdf_url = f"https://files.rcsb.org/ligands/download/{ligand_code}_ideal.sdf"
            headers = {
                'User-Agent': 'BioPipelines-Ligand/1.0 (https://github.com/locbp-uzh/biopipelines)'
            }
            response = requests.get(sdf_url, headers=headers, timeout=30)
            response.raise_for_status()
            sdf_content = response.text

            if not sdf_content.strip():
                raise ValueError(f"Downloaded SDF file is empty")

            content = convert_sdf_to_pdb(sdf_content, residue_code)

        if content is None:
            raise ValueError(f"Failed to convert ligand to {output_format.upper()} format")

        # Save to Ligands/ folder for caching (using ligand_code as filename)
        os.makedirs(repo_ligands_folder, exist_ok=True)
        cache_path = os.path.join(repo_ligands_folder, f"{ligand_code}.{ext}")
        with open(cache_path, 'w') as f:
            f.write(content)
        print(f"  Cached to Ligands/ folder: {cache_path}")

        # Save companion CSV with metadata
        cache_csv_path = os.path.join(repo_ligands_folder, f"{ligand_code}.csv")
        cache_metadata = {
            'id': ligand_code,
            'format': 'ccd',
            'code': residue_code,
            'lookup': ligand_code,
            'source': 'rcsb',
            'ccd': ligand_code,
            'cid': '',
            'cas': '',
            'smiles': smiles,
            'name': props.get('name', ''),
            'formula': props.get('formula', ''),
        }
        pd.DataFrame([cache_metadata]).to_csv(cache_csv_path, index=False)
        print(f"  Cached metadata: {cache_csv_path}")

        # Save to output folder (using custom_id)
        output_filename = f"{custom_id}.{ext}"
        output_path = os.path.join(output_folder, output_filename)
        with open(output_path, 'w') as f:
            f.write(content)

        file_size = os.path.getsize(output_path)

        metadata = {
            "file_size": file_size,
            "source": "rcsb",
            "ccd": ligand_code,
            "cid": "",
            "cas": "",
            "smiles": smiles,
            "name": props.get('name', ''),
            "formula": props.get('formula', ''),
        }

        print(f"Successfully downloaded {ligand_code} as {custom_id}.{ext}: {file_size} bytes")
        return True, output_path, metadata

    except Exception as e:
        error_msg = f"Error downloading {ligand_code} from RCSB: {str(e)}"
        print(f"Error: {error_msg}")

        metadata = {
            "error_message": error_msg,
            "source": "rcsb_download_failed",
            "attempted_path": f"rcsb:{ligand_code}"
        }

        return False, "", metadata


def download_from_pubchem(lookup: str, lookup_type: str, custom_id: str, residue_code: str,
                          output_folder: str, repo_ligands_folder: str,
                          output_format: str = "pdb") -> Tuple[bool, str, Dict[str, Any]]:
    """
    Download a ligand from PubChem.

    Uses RDKit to generate PDB or CIF from SMILES. CIF format includes explicit
    bond orders which is recommended for tools like RFdiffusion3.

    Args:
        lookup: The lookup value (name, CID, or CAS)
        lookup_type: "cid", "cas", or "name"
        custom_id: Custom ID for renaming the ligand
        residue_code: Code to use in output file
        output_folder: Directory to save the ligand
        repo_ligands_folder: Repository Ligands folder for caching
        output_format: "pdb" or "cif" (default: "pdb")

    Returns:
        Tuple of (success: bool, file_path: str, metadata: dict)
    """
    try:
        print(f"Fetching {lookup} from PubChem ({lookup_type})...")

        # Fetch properties from PubChem (SMILES, name, formula - no SDF yet)
        pubchem_data = fetch_properties_from_pubchem(lookup, lookup_type)
        smiles = pubchem_data.get('smiles', '')
        cid = pubchem_data.get('cid', '')

        content = None
        ext = output_format  # "pdb" or "cif"

        # Convert SMILES to requested format using RDKit
        if smiles:
            if output_format == "cif":
                print(f"  Converting SMILES to CIF using RDKit (with bond orders)...")
                content = convert_smiles_to_cif_rdkit(smiles, residue_code)
            else:
                print(f"  Converting SMILES to PDB using RDKit...")
                content = convert_smiles_to_pdb_rdkit(smiles, residue_code)

        # Fallback for PDB: Download SDF and use OpenBabel (only if RDKit failed)
        if content is None and output_format == "pdb":
            print(f"  RDKit conversion failed or no SMILES, falling back to OpenBabel SDF...")
            sdf_content = fetch_sdf_from_pubchem(cid)
            content = convert_sdf_to_pdb(sdf_content, residue_code)

        if content is None:
            raise ValueError(f"Failed to convert ligand to {output_format.upper()} format")

        # Save to Ligands/ folder for caching (always cache both formats if possible)
        os.makedirs(repo_ligands_folder, exist_ok=True)
        cache_path = os.path.join(repo_ligands_folder, f"{lookup}.{ext}")
        with open(cache_path, 'w') as f:
            f.write(content)
        print(f"  Cached to Ligands/ folder: {cache_path}")

        # Save companion CSV with metadata
        cache_csv_path = os.path.join(repo_ligands_folder, f"{lookup}.csv")
        cache_metadata = {
            'id': lookup,
            'format': 'smiles',
            'code': residue_code,
            'lookup': lookup,
            'source': 'pubchem',
            'ccd': '',
            'cid': cid,
            'cas': lookup if lookup_type == 'cas' else '',
            'smiles': smiles,
            'name': pubchem_data.get('name', ''),
            'formula': pubchem_data.get('formula', ''),
        }
        pd.DataFrame([cache_metadata]).to_csv(cache_csv_path, index=False)
        print(f"  Cached metadata: {cache_csv_path}")

        # Save to output folder (using custom_id)
        output_filename = f"{custom_id}.{ext}"
        output_path = os.path.join(output_folder, output_filename)
        with open(output_path, 'w') as f:
            f.write(content)

        file_size = os.path.getsize(output_path)

        metadata = {
            "file_size": file_size,
            "source": "pubchem",
            "ccd": "",
            "cid": cid,
            "cas": lookup if lookup_type == 'cas' else '',
            "smiles": smiles,
            "name": pubchem_data.get('name', ''),
            "formula": pubchem_data.get('formula', ''),
        }

        print(f"Successfully downloaded {lookup} as {custom_id}.{ext}: {file_size} bytes")
        return True, output_path, metadata

    except Exception as e:
        error_msg = f"Error downloading {lookup} from PubChem: {str(e)}"
        print(f"Error: {error_msg}")

        metadata = {
            "error_message": error_msg,
            "source": "pubchem_download_failed",
            "attempted_path": f"pubchem:{lookup}"
        }

        return False, "", metadata


def generate_from_smiles(smiles: str, custom_id: str, residue_code: str,
                          output_folder: str, output_format: str = "pdb") -> Tuple[bool, str, Dict[str, Any]]:
    """
    Generate a ligand structure file from a SMILES string using RDKit.

    Args:
        smiles: SMILES string of the molecule
        custom_id: Custom ID for output filename
        residue_code: Residue code to use in output file
        output_folder: Directory to save the ligand
        output_format: "pdb" or "cif" (default: "pdb")

    Returns:
        Tuple of (success: bool, file_path: str, metadata: dict)
    """
    try:
        print(f"Generating structure from SMILES: {smiles[:50]}..." if len(smiles) > 50 else f"Generating structure from SMILES: {smiles}")

        content = None
        ext = output_format

        # Convert SMILES to requested format using RDKit
        if output_format == "cif":
            print(f"  Converting SMILES to CIF using RDKit (with bond orders)...")
            content = convert_smiles_to_cif_rdkit(smiles, residue_code)
        else:
            print(f"  Converting SMILES to PDB using RDKit...")
            content = convert_smiles_to_pdb_rdkit(smiles, residue_code)

        if content is None:
            raise ValueError(f"Failed to convert SMILES to {output_format.upper()} format")

        # Save to output folder
        output_filename = f"{custom_id}.{ext}"
        output_path = os.path.join(output_folder, output_filename)
        with open(output_path, 'w') as f:
            f.write(content)

        file_size = os.path.getsize(output_path)

        metadata = {
            "file_size": file_size,
            "source": "smiles",
            "ccd": "",
            "cid": "",
            "cas": "",
            "smiles": smiles,
            "name": "",
            "formula": "",
        }

        print(f"Successfully generated {custom_id}.{ext}: {file_size} bytes")
        return True, output_path, metadata

    except Exception as e:
        error_msg = f"Error generating from SMILES: {str(e)}"
        print(f"Error: {error_msg}")

        metadata = {
            "error_message": error_msg,
            "source": "smiles_generation_failed",
            "attempted_path": f"smiles:{smiles[:30]}..."
        }

        return False, "", metadata


def fetch_ligands(config_data: Dict[str, Any]) -> int:
    """
    Fetch multiple ligands with priority-based lookup and/or generate from SMILES.

    Priority for lookup: local_folder -> Ligands/ -> RCSB/PubChem download
    SMILES are always generated directly using RDKit.

    Args:
        config_data: Configuration dictionary with fetch parameters

    Returns:
        Number of failed fetches
    """
    custom_ids = config_data['custom_ids']
    residue_codes = config_data['residue_codes']
    lookup_values = config_data.get('lookup_values', [])
    smiles_values = config_data.get('smiles_values', [])
    source = config_data.get('source')  # "rcsb", "pubchem", or None (auto-detect)
    local_folder = config_data.get('local_folder')
    output_format = config_data.get('output_format', 'pdb')  # "pdb" or "cif"
    repo_ligands_folder = config_data['repo_ligands_folder']
    output_folder = config_data['output_folder']
    compounds_table = config_data['compounds_table']
    failed_table = config_data['failed_table']

    total_count = len(lookup_values) + len(smiles_values)
    print(f"Processing {total_count} ligands ({len(lookup_values)} lookup, {len(smiles_values)} SMILES)")
    print(f"Output format: {output_format.upper()}")
    if lookup_values:
        if source:
            print(f"Forced source: {source}")
        else:
            print(f"Source: auto-detect")
        print(f"Priority: {'local_folder -> ' if local_folder else ''}Ligands/ -> download")

    # Create output directory
    os.makedirs(output_folder, exist_ok=True)

    # Track results
    successful_downloads = []
    failed_downloads = []

    # Process lookup values first
    for i, (custom_id, residue_code, lookup) in enumerate(zip(
            custom_ids[:len(lookup_values)],
            residue_codes[:len(lookup_values)],
            lookup_values), 1):
        print(f"\n[{i}/{total_count}] Processing lookup {lookup} -> {custom_id} (code: {residue_code})")

        # Detect lookup type
        lookup_type = detect_lookup_type(lookup)
        print(f"  Detected lookup type: {lookup_type}")

        # Determine effective source
        if source:
            effective_source = source
        else:
            # Auto-detect
            effective_source = "rcsb" if lookup_type == "ccd" else "pubchem"
        print(f"  Source: {effective_source}")

        # Try to find locally first (check for requested format)
        local_file_path, local_csv_path = find_local_ligand(lookup, local_folder, repo_ligands_folder, output_format)

        if local_file_path:
            # Copy from local
            success, file_path, metadata = copy_local_ligand(
                lookup, custom_id, residue_code, local_file_path, local_csv_path, output_folder, output_format
            )
        else:
            # Download from appropriate source
            print(f"  {lookup} not found locally, downloading from {effective_source}")

            if effective_source == "rcsb":
                success, file_path, metadata = download_from_rcsb(
                    lookup, custom_id, residue_code, output_folder, repo_ligands_folder, output_format
                )
            else:
                success, file_path, metadata = download_from_pubchem(
                    lookup, lookup_type, custom_id, residue_code, output_folder, repo_ligands_folder, output_format
                )

        if success:
            # Determine format: 'ccd' if has CCD code, otherwise 'smiles'
            ccd_val = metadata.get('ccd', '')
            fmt = 'ccd' if ccd_val else 'smiles'
            successful_downloads.append({
                'id': custom_id,
                'format': fmt,
                'code': residue_code,
                'lookup': lookup,
                'source': metadata.get('source', effective_source),
                'ccd': ccd_val,
                'cid': metadata.get('cid', ''),
                'cas': metadata.get('cas', ''),
                'smiles': metadata.get('smiles', ''),
                'name': metadata.get('name', ''),
                'formula': metadata.get('formula', ''),
                'file_path': file_path
            })
        else:
            failed_downloads.append({
                'lookup': lookup,
                'error_message': metadata.get('error_message', 'Unknown error'),
                'source': metadata.get('source', 'unknown'),
                'attempted_path': metadata.get('attempted_path', '')
            })

    # Process SMILES values
    smiles_start_idx = len(lookup_values)
    for i, (custom_id, residue_code, smiles) in enumerate(zip(
            custom_ids[smiles_start_idx:],
            residue_codes[smiles_start_idx:],
            smiles_values), smiles_start_idx + 1):
        print(f"\n[{i}/{total_count}] Processing SMILES -> {custom_id} (code: {residue_code})")

        success, file_path, metadata = generate_from_smiles(
            smiles, custom_id, residue_code, output_folder, output_format
        )

        if success:
            successful_downloads.append({
                'id': custom_id,
                'format': 'smiles',
                'code': residue_code,
                'lookup': '',  # No lookup for direct SMILES
                'source': 'smiles',
                'ccd': '',
                'cid': '',
                'cas': '',
                'smiles': smiles,
                'name': '',
                'formula': '',
                'file_path': file_path
            })
        else:
            failed_downloads.append({
                'lookup': f"SMILES:{custom_id}",
                'error_message': metadata.get('error_message', 'Unknown error'),
                'source': 'smiles_generation_failed',
                'attempted_path': metadata.get('attempted_path', '')
            })

    # Save successful downloads table
    if successful_downloads:
        df_success = pd.DataFrame(successful_downloads)
        df_success.to_csv(compounds_table, index=False)
        print(f"\nSuccessful fetches saved: {compounds_table} ({len(successful_downloads)} ligands)")
    else:
        # Create empty table with proper columns
        columns = ["id", "format", "code", "lookup", "source", "ccd", "cid", "cas", "smiles", "name", "formula", "file_path"]
        empty_df = pd.DataFrame(columns=columns)
        empty_df.to_csv(compounds_table, index=False)
        print(f"No successful fetches - created empty table: {compounds_table}")

    # Save failed downloads table
    if failed_downloads:
        df_failed = pd.DataFrame(failed_downloads)
        df_failed.to_csv(failed_table, index=False)
        print(f"Failed fetches saved: {failed_table} ({len(failed_downloads)} failures)")
    else:
        # Create empty failed downloads table
        empty_failed_df = pd.DataFrame(columns=["lookup", "error_message", "source", "attempted_path"])
        empty_failed_df.to_csv(failed_table, index=False)
        print("No failed fetches")

    # Summary
    print(f"\n=== FETCH SUMMARY ===")
    print(f"Requested: {total_count} ligands")
    print(f"Successful: {len(successful_downloads)}")
    print(f"Failed: {len(failed_downloads)}")
    if total_count > 0:
        print(f"Success rate: {len(successful_downloads)/total_count*100:.1f}%")

    if successful_downloads:
        # Count how many have SMILES
        with_smiles = sum(1 for item in successful_downloads if item['smiles'])
        print(f"SMILES available: {with_smiles}/{len(successful_downloads)}")

    # Log any failed fetches
    if failed_downloads:
        print(f"\nFailed fetches:")
        for failure in failed_downloads:
            print(f"  - {failure['lookup']}: {failure['error_message']}")

    # Return the number of failures
    return len(failed_downloads)


def main():
    parser = argparse.ArgumentParser(description='Fetch ligands with priority-based lookup')
    parser.add_argument('--config', required=True, help='JSON config file with fetch parameters')

    args = parser.parse_args()

    # Load configuration
    if not os.path.exists(args.config):
        print(f"Error: Config file not found: {args.config}")
        sys.exit(1)

    try:
        with open(args.config, 'r') as f:
            config_data = json.load(f)
    except Exception as e:
        print(f"Error loading config: {e}")
        sys.exit(1)

    # Validate required parameters
    required_params = ['custom_ids', 'residue_codes', 'repo_ligands_folder',
                       'output_folder', 'compounds_table', 'failed_table']
    for param in required_params:
        if param not in config_data:
            print(f"Error: Missing required parameter: {param}")
            sys.exit(1)

    # Ensure at least one of lookup_values or smiles_values is present
    lookup_values = config_data.get('lookup_values', [])
    smiles_values = config_data.get('smiles_values', [])
    if not lookup_values and not smiles_values:
        print("Error: Must provide at least one of 'lookup_values' or 'smiles_values'")
        sys.exit(1)

    try:
        failed_count = fetch_ligands(config_data)

        # Fail if ANY fetches failed
        if failed_count > 0:
            print(f"\nERROR: {failed_count} ligand fetch(es) failed")
            print("Pipeline cannot continue with incomplete ligand set")
            print("Check failed_downloads.csv for details")
            sys.exit(1)

        print("\nAll ligands fetched successfully")

    except Exception as e:
        print(f"Error fetching ligands: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
