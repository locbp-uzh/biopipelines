#!/usr/bin/env python3
"""
Runtime helper script for Ligand tool.

Fetches small molecule ligands with priority-based lookup: local_folder -> Ligands/ -> RCSB/PubChem download.
Downloads SDF files and converts them to PDB format with proper atom numbering.
Supports both RCSB (CCD codes) and PubChem (name, CID, CAS) as sources.
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
                      repo_ligands_folder: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Find ligand file locally.

    Priority: local_folder (if given) -> repo_ligands_folder -> None

    Args:
        lookup: Lookup value (used as filename stem)
        local_folder: Custom local folder (can be None)
        repo_ligands_folder: Repository Ligands folder

    Returns:
        Tuple of (path to local PDB file, path to local CSV file) or (None, None) if not found
    """
    search_locations = []

    if local_folder:
        search_locations.append(local_folder)
    search_locations.append(repo_ligands_folder)

    for location in search_locations:
        pdb_candidate = os.path.join(location, f"{lookup}.pdb")
        csv_candidate = os.path.join(location, f"{lookup}.csv")
        if os.path.exists(pdb_candidate):
            print(f"Found {lookup} locally: {pdb_candidate}")
            csv_path = csv_candidate if os.path.exists(csv_candidate) else None
            return pdb_candidate, csv_path

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
                      source_pdb_path: str, source_csv_path: Optional[str],
                      output_folder: str) -> Tuple[bool, str, Dict[str, Any]]:
    """
    Copy local ligand file to output folder with normalization.

    Args:
        lookup: Original lookup value
        custom_id: Custom ID for output filename
        residue_code: 3-letter residue code to use in PDB
        source_pdb_path: Path to local ligand PDB file
        source_csv_path: Path to local ligand CSV metadata file (can be None)
        output_folder: Directory to save the ligand

    Returns:
        Tuple of (success: bool, file_path: str, metadata: dict)
    """
    try:
        with open(source_pdb_path, 'r') as f:
            pdb_content = f.read()

        # Normalize the PDB content
        normalized_pdb = normalize_pdb_content(pdb_content, residue_code)

        filename = f"{custom_id}.pdb"
        output_path = os.path.join(output_folder, filename)

        with open(output_path, 'w') as f:
            f.write(normalized_pdb)

        file_size = os.path.getsize(output_path)

        # Load metadata from CSV if available
        metadata = {}
        if source_csv_path:
            metadata = load_local_metadata(source_csv_path)

        # Update/override certain fields
        metadata.update({
            "file_size": file_size,
            "source": metadata.get('source', 'local'),
            "source_path": source_pdb_path,
        })

        # Try to fetch SMILES if not in metadata
        if not metadata.get('smiles'):
            lookup_type = detect_lookup_type(lookup)
            if lookup_type == "ccd":
                smiles = fetch_smiles_from_rcsb(lookup.upper())
                if smiles:
                    metadata['smiles'] = smiles
                    metadata['ccd'] = lookup.upper()
            else:
                # Try PubChem for name/CID/CAS lookups
                smiles = fetch_smiles_from_pubchem(lookup, lookup_type)
                if smiles:
                    metadata['smiles'] = smiles

        print(f"Successfully copied {lookup} as {custom_id}: {file_size} bytes (from local)")
        return True, output_path, metadata

    except Exception as e:
        error_msg = f"Error copying local file {lookup}: {str(e)}"
        print(f"Error: {error_msg}")
        metadata = {
            "error_message": error_msg,
            "source": "local_copy_failed",
            "attempted_path": source_pdb_path
        }
        return False, "", metadata


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
    """
    lines = pdb_content.split('\n')
    output_lines = []
    atom_serial = 0

    # Track element counts for generating unique atom names
    element_counts = {}

    # Ensure residue code is uppercase and max 3 chars
    res_code = residue_code.upper()[:3].ljust(3)

    for line in lines:
        if line.startswith(('HETATM', 'ATOM')):
            atom_serial += 1

            # Extract element symbol
            element = _extract_element_from_line(line)

            # Generate numbered atom name
            element_counts[element] = element_counts.get(element, 0) + 1
            atom_num = element_counts[element]

            # Format atom name (4 chars, columns 12-15)
            # PDB convention: 1-char elements right-justify in col 13-14
            # For ligands with numbers: typically "C1", "C12", etc.
            atom_name_str = f"{element}{atom_num}"
            if len(element) == 1:
                # 1-char element: " C1 " or " C12"
                atom_name = f" {atom_name_str:<3}"
            else:
                # 2-char element: "CL1 " or "CL12"
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

            # Pad coords to proper length if needed
            coords = coords.ljust(39)

            # Format element for columns 76-77 (right-justified)
            element_col = f"{element:>2}"

            # Build the new line (80 chars total for proper PDB format)
            # Columns: 0-5(record), 6-10(serial), 11(space), 12-15(atom), 16(altloc),
            #          17-19(resname), 20(space), 21(chain), 22-25(resseq), 26(icode),
            #          27-65(coords+occupancy+temp), 66-75(spaces), 76-77(element)
            new_line = f"{record}{serial_str} {atom_name}{alt_loc}{res_code} {chain}{res_seq}{icode}{coords}          {element_col}"
            output_lines.append(new_line)

        elif line.startswith('CONECT'):
            # Skip CONECT records - they reference old serial numbers
            continue
        elif line.startswith('END'):
            output_lines.append(line)
        else:
            # Keep other records (HEADER, REMARK, etc.)
            output_lines.append(line)

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
            'User-Agent': 'BioPipelines-Ligand/1.0 (https://gitlab.uzh.ch/locbp/public/biopipelines)'
        }

        response = requests.get(url, headers=headers, timeout=10)
        response.raise_for_status()

        data = response.json()

        # Extract SMILES from the response
        if 'chem_comp' in data:
            chem_comp = data['chem_comp']
            # Try different possible fields for SMILES
            for field in ['pdbx_smiles_canonical', 'smiles', 'smiles_canonical']:
                if field in chem_comp and chem_comp[field]:
                    return chem_comp[field]

        # Alternative location in descriptors
        if 'rcsb_chem_comp_descriptor' in data:
            descriptors = data['rcsb_chem_comp_descriptor']
            if isinstance(descriptors, dict):
                if 'smiles_canonical' in descriptors:
                    return descriptors['smiles_canonical']
                elif 'smiles' in descriptors:
                    return descriptors['smiles']

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
        headers = {'User-Agent': 'BioPipelines-Ligand/1.0 (https://gitlab.uzh.ch/locbp/public/biopipelines)'}

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

        # Get SMILES
        props_url = f"{base_url}/compound/cid/{cid}/property/CanonicalSMILES/JSON"
        response = requests.get(props_url, headers=headers, timeout=30)
        response.raise_for_status()
        props = response.json()['PropertyTable']['Properties'][0]
        return props.get('CanonicalSMILES', None)

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
            'User-Agent': 'BioPipelines-Ligand/1.0 (https://gitlab.uzh.ch/locbp/public/biopipelines)'
        }

        response = requests.get(url, headers=headers, timeout=10)
        response.raise_for_status()

        data = response.json()
        result = {'ccd': ligand_code}

        if 'chem_comp' in data:
            chem_comp = data['chem_comp']
            result['name'] = chem_comp.get('name', '')
            result['formula'] = chem_comp.get('formula', '')

            # Get SMILES
            for field in ['pdbx_smiles_canonical', 'smiles', 'smiles_canonical']:
                if field in chem_comp and chem_comp[field]:
                    result['smiles'] = chem_comp[field]
                    break

        # Alternative location for SMILES
        if 'smiles' not in result and 'rcsb_chem_comp_descriptor' in data:
            descriptors = data['rcsb_chem_comp_descriptor']
            if isinstance(descriptors, dict):
                result['smiles'] = descriptors.get('smiles_canonical') or descriptors.get('smiles', '')

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


def fetch_from_pubchem(lookup: str, lookup_type: str) -> Dict[str, Any]:
    """
    Fetch compound data from PubChem.

    Args:
        lookup: The lookup value (name, CID, or CAS)
        lookup_type: "cid", "cas", or "name"

    Returns:
        Dictionary with: cid, smiles, name, formula, sdf_content

    Raises:
        ValueError: If compound not found or download fails
    """
    import requests

    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    headers = {'User-Agent': 'BioPipelines-Ligand/1.0 (https://gitlab.uzh.ch/locbp/public/biopipelines)'}

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
    props_url = f"{base_url}/compound/cid/{cid}/property/CanonicalSMILES,IUPACName,MolecularFormula/JSON"
    response = requests.get(props_url, headers=headers, timeout=30)
    response.raise_for_status()
    props = response.json()['PropertyTable']['Properties'][0]

    smiles = props.get('CanonicalSMILES', '')
    name = props.get('IUPACName', '')
    formula = props.get('MolecularFormula', '')

    print(f"  SMILES: {smiles[:50]}..." if len(smiles) > 50 else f"  SMILES: {smiles}")
    print(f"  Name: {name[:50]}..." if len(name) > 50 else f"  Name: {name}")

    # Step 3: Get 3D SDF (try 3D first, fall back to 2D)
    sdf_url = f"{base_url}/compound/cid/{cid}/SDF?record_type=3d"
    print(f"  Downloading 3D SDF...")
    response = requests.get(sdf_url, headers=headers, timeout=30)

    if response.status_code == 404:
        # No 3D structure available, try 2D
        print(f"  No 3D structure available, using 2D...")
        sdf_url = f"{base_url}/compound/cid/{cid}/SDF"
        response = requests.get(sdf_url, headers=headers, timeout=30)

    response.raise_for_status()
    sdf_content = response.text

    if not sdf_content.strip():
        raise ValueError(f"Downloaded SDF file for CID {cid} is empty")

    return {
        'cid': str(cid),
        'smiles': smiles,
        'name': name,
        'formula': formula,
        'sdf_content': sdf_content,
        'source': 'pubchem'
    }


def download_from_rcsb(ligand_code: str, custom_id: str, residue_code: str,
                       output_folder: str, repo_ligands_folder: str) -> Tuple[bool, str, Dict[str, Any]]:
    """
    Download a single ligand from RCSB PDB.

    Downloads the ideal SDF file and converts it to PDB format.

    Args:
        ligand_code: 3-letter ligand code (e.g., 'ATP')
        custom_id: Custom ID for renaming the ligand
        residue_code: 3-letter code to use in PDB file
        output_folder: Directory to save the ligand
        repo_ligands_folder: Repository Ligands folder for caching

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

    # RCSB ideal SDF download URL
    sdf_url = f"https://files.rcsb.org/ligands/download/{ligand_code}_ideal.sdf"

    try:
        print(f"Downloading {ligand_code} ideal SDF from RCSB: {sdf_url}")

        import requests

        # Download with proper headers
        headers = {
            'User-Agent': 'BioPipelines-Ligand/1.0 (https://gitlab.uzh.ch/locbp/public/biopipelines)'
        }

        response = requests.get(sdf_url, headers=headers, timeout=30)
        response.raise_for_status()

        sdf_content = response.text

        # Validate SDF content
        if not sdf_content.strip():
            raise ValueError(f"Downloaded SDF file is empty")

        # Get additional properties from RCSB
        props = fetch_properties_from_rcsb(ligand_code)

        # Convert SDF to PDB (normalization happens inside)
        pdb_content = convert_sdf_to_pdb(sdf_content, residue_code)
        if pdb_content is None:
            raise ValueError(f"Failed to convert SDF to PDB format")

        # Save to Ligands/ folder for caching (using ligand_code as filename)
        os.makedirs(repo_ligands_folder, exist_ok=True)
        cache_pdb_path = os.path.join(repo_ligands_folder, f"{ligand_code}.pdb")
        with open(cache_pdb_path, 'w') as f:
            f.write(pdb_content)
        print(f"  Cached to Ligands/ folder: {cache_pdb_path}")

        # Save companion CSV with metadata
        cache_csv_path = os.path.join(repo_ligands_folder, f"{ligand_code}.csv")
        cache_metadata = {
            'id': ligand_code,
            'code': residue_code,
            'lookup': ligand_code,
            'source': 'rcsb',
            'ccd': ligand_code,
            'cid': '',
            'cas': '',
            'smiles': props.get('smiles', ''),
            'name': props.get('name', ''),
            'formula': props.get('formula', ''),
        }
        pd.DataFrame([cache_metadata]).to_csv(cache_csv_path, index=False)
        print(f"  Cached metadata: {cache_csv_path}")

        # Save to output folder (using custom_id)
        output_filename = f"{custom_id}.pdb"
        output_path = os.path.join(output_folder, output_filename)
        with open(output_path, 'w') as f:
            f.write(pdb_content)

        file_size = os.path.getsize(output_path)

        metadata = {
            "file_size": file_size,
            "source": "rcsb",
            "url": sdf_url,
            "ccd": ligand_code,
            "cid": "",
            "cas": "",
            "smiles": props.get('smiles', ''),
            "name": props.get('name', ''),
            "formula": props.get('formula', ''),
        }

        print(f"Successfully downloaded {ligand_code} as {custom_id}: {file_size} bytes")
        return True, output_path, metadata

    except Exception as e:
        error_type = type(e).__name__
        error_msg = f"Error downloading {ligand_code} from RCSB: {str(e)}"
        print(f"Error: {error_msg}")

        metadata = {
            "error_message": error_msg,
            "source": "rcsb_download_failed",
            "attempted_path": sdf_url
        }

        return False, "", metadata


def download_from_pubchem(lookup: str, lookup_type: str, custom_id: str, residue_code: str,
                          output_folder: str, repo_ligands_folder: str) -> Tuple[bool, str, Dict[str, Any]]:
    """
    Download a ligand from PubChem.

    Args:
        lookup: The lookup value (name, CID, or CAS)
        lookup_type: "cid", "cas", or "name"
        custom_id: Custom ID for renaming the ligand
        residue_code: 3-letter code to use in PDB file
        output_folder: Directory to save the ligand
        repo_ligands_folder: Repository Ligands folder for caching

    Returns:
        Tuple of (success: bool, file_path: str, metadata: dict)
    """
    try:
        print(f"Fetching {lookup} from PubChem ({lookup_type})...")

        # Fetch from PubChem
        pubchem_data = fetch_from_pubchem(lookup, lookup_type)

        # Convert SDF to PDB
        pdb_content = convert_sdf_to_pdb(pubchem_data['sdf_content'], residue_code)
        if pdb_content is None:
            raise ValueError(f"Failed to convert SDF to PDB format")

        # Save to Ligands/ folder for caching
        os.makedirs(repo_ligands_folder, exist_ok=True)
        cache_pdb_path = os.path.join(repo_ligands_folder, f"{lookup}.pdb")
        with open(cache_pdb_path, 'w') as f:
            f.write(pdb_content)
        print(f"  Cached to Ligands/ folder: {cache_pdb_path}")

        # Save companion CSV with metadata
        cache_csv_path = os.path.join(repo_ligands_folder, f"{lookup}.csv")
        cache_metadata = {
            'id': lookup,
            'code': residue_code,
            'lookup': lookup,
            'source': 'pubchem',
            'ccd': '',
            'cid': pubchem_data.get('cid', ''),
            'cas': lookup if lookup_type == 'cas' else '',
            'smiles': pubchem_data.get('smiles', ''),
            'name': pubchem_data.get('name', ''),
            'formula': pubchem_data.get('formula', ''),
        }
        pd.DataFrame([cache_metadata]).to_csv(cache_csv_path, index=False)
        print(f"  Cached metadata: {cache_csv_path}")

        # Save to output folder (using custom_id)
        output_filename = f"{custom_id}.pdb"
        output_path = os.path.join(output_folder, output_filename)
        with open(output_path, 'w') as f:
            f.write(pdb_content)

        file_size = os.path.getsize(output_path)

        metadata = {
            "file_size": file_size,
            "source": "pubchem",
            "ccd": "",
            "cid": pubchem_data.get('cid', ''),
            "cas": lookup if lookup_type == 'cas' else '',
            "smiles": pubchem_data.get('smiles', ''),
            "name": pubchem_data.get('name', ''),
            "formula": pubchem_data.get('formula', ''),
        }

        print(f"Successfully downloaded {lookup} as {custom_id}: {file_size} bytes")
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


def fetch_ligands(config_data: Dict[str, Any]) -> int:
    """
    Fetch multiple ligands with priority-based lookup.

    Priority: local_folder -> Ligands/ -> RCSB/PubChem download

    Args:
        config_data: Configuration dictionary with fetch parameters

    Returns:
        Number of failed fetches
    """
    custom_ids = config_data['custom_ids']
    residue_codes = config_data['residue_codes']
    lookup_values = config_data['lookup_values']
    source = config_data.get('source')  # "rcsb", "pubchem", or None (auto-detect)
    local_folder = config_data.get('local_folder')
    repo_ligands_folder = config_data['repo_ligands_folder']
    output_folder = config_data['output_folder']
    compounds_table = config_data['compounds_table']
    failed_table = config_data['failed_table']

    print(f"Fetching {len(lookup_values)} ligands")
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

    # Fetch each ligand
    for i, (custom_id, residue_code, lookup) in enumerate(zip(custom_ids, residue_codes, lookup_values), 1):
        print(f"\n[{i}/{len(lookup_values)}] Processing {lookup} -> {custom_id} (code: {residue_code})")

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

        # Try to find locally first
        local_pdb_path, local_csv_path = find_local_ligand(lookup, local_folder, repo_ligands_folder)

        if local_pdb_path:
            # Copy from local
            success, file_path, metadata = copy_local_ligand(
                lookup, custom_id, residue_code, local_pdb_path, local_csv_path, output_folder
            )
        else:
            # Download from appropriate source
            print(f"  {lookup} not found locally, downloading from {effective_source}")

            if effective_source == "rcsb":
                success, file_path, metadata = download_from_rcsb(
                    lookup, custom_id, residue_code, output_folder, repo_ligands_folder
                )
            else:
                success, file_path, metadata = download_from_pubchem(
                    lookup, lookup_type, custom_id, residue_code, output_folder, repo_ligands_folder
                )

        if success:
            successful_downloads.append({
                'id': custom_id,
                'code': residue_code,
                'lookup': lookup,
                'source': metadata.get('source', effective_source),
                'ccd': metadata.get('ccd', ''),
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

    # Save successful downloads table
    if successful_downloads:
        df_success = pd.DataFrame(successful_downloads)
        df_success.to_csv(compounds_table, index=False)
        print(f"\nSuccessful fetches saved: {compounds_table} ({len(successful_downloads)} ligands)")
    else:
        # Create empty table with proper columns
        columns = ["id", "code", "lookup", "source", "ccd", "cid", "cas", "smiles", "name", "formula", "file_path"]
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
    print(f"Requested: {len(lookup_values)} ligands")
    print(f"Successful: {len(successful_downloads)}")
    print(f"Failed: {len(failed_downloads)}")
    print(f"Success rate: {len(successful_downloads)/len(lookup_values)*100:.1f}%")

    if successful_downloads:
        # Count how many have SMILES
        with_smiles = sum(1 for item in successful_downloads if item['smiles'])
        print(f"SMILES fetched: {with_smiles}/{len(successful_downloads)}")

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
    required_params = ['custom_ids', 'residue_codes', 'lookup_values', 'repo_ligands_folder',
                       'output_folder', 'compounds_table', 'failed_table']
    for param in required_params:
        if param not in config_data:
            print(f"Error: Missing required parameter: {param}")
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
