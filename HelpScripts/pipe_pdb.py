#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Runtime helper script for PDB tool.

Fetches protein structures with priority-based lookup: local_folder -> PDBs/ -> RCSB download.
Downloads are saved to both PDBs/ folder (for reuse) and tool output folder.
"""

import os
import sys
import argparse
import json
import pandas as pd
from datetime import datetime
from typing import Dict, List, Any, Tuple, Optional
from pathlib import Path

# Import the local PDB parser
sys.path.append(os.path.dirname(__file__))
from pdb_parser import get_protein_sequence, parse_pdb_file


def convert_cif_to_pdb(cif_content: str) -> str:
    """
    Convert CIF format to PDB format using BioPython.

    Args:
        cif_content: CIF file content

    Returns:
        PDB format content

    Raises:
        Exception: If conversion fails
    """
    import tempfile
    from Bio.PDB import MMCIFParser, PDBIO

    # Write CIF content to temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.cif', delete=False) as tmp_cif:
        tmp_cif.write(cif_content)
        tmp_cif.flush()
        cif_path = tmp_cif.name

    try:
        # Parse CIF file
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("structure", cif_path)

        # Write to PDB format
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp_pdb:
            pdb_path = tmp_pdb.name

        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_path)

        # Read PDB content
        with open(pdb_path, 'r') as f:
            pdb_content = f.read()

        # Clean up temporary files
        os.unlink(cif_path)
        os.unlink(pdb_path)

        return pdb_content

    except Exception as e:
        # Clean up on error
        if os.path.exists(cif_path):
            os.unlink(cif_path)
        raise Exception(f"CIF to PDB conversion failed: {str(e)}")


def build_rename_mapping(operations: List[Dict[str, Any]]) -> Dict[str, str]:
    """
    Build a mapping from old residue names to new names based on rename operations.

    Args:
        operations: List of operation dictionaries

    Returns:
        Dict mapping old_name -> new_name for rename operations
    """
    mapping = {}
    for op in operations:
        if op.get("op") == "rename":
            mapping[op["old"]] = op["new"]
    return mapping


def apply_operations(content: str, format: str, operations: List[Dict[str, Any]]) -> str:
    """
    Apply a sequence of operations to structure content.

    Args:
        content: Structure file content
        format: File format ("pdb" or "cif")
        operations: List of operation dictionaries with 'op' key and parameters

    Returns:
        Modified structure content
    """
    for op in operations:
        op_type = op.get("op")
        if op_type == "rename":
            content = apply_rename_operation(content, format, op["old"], op["new"])
        else:
            print(f"Warning: Unknown operation type '{op_type}', skipping")

    return content


def apply_rename_operation(content: str, format: str, old_name: str, new_name: str) -> str:
    """
    Rename a residue/ligand in the structure.

    For PDB format, modifies columns 18-20 (residue name) in ATOM/HETATM records.
    The new name is padded/truncated to fit the 3-character field, unless it's
    a special format like ":L:" which uses the full field width.

    Args:
        content: Structure file content
        format: File format ("pdb" or "cif")
        old_name: Current residue name (e.g., "LIG")
        new_name: New residue name (e.g., ":L:")

    Returns:
        Structure content with renamed residues
    """
    if format != "pdb":
        print(f"Warning: Rename operation not yet implemented for CIF format")
        return content

    lines = content.split('\n')
    modified_lines = []
    rename_count = 0

    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            # PDB format: columns 18-20 (0-indexed: 17-20) contain residue name
            # The field is 3 characters, right-justified with spaces
            current_res_name = line[17:20].strip()

            if current_res_name == old_name:
                # Format new name to fit in 3-character field
                # Right-justify, but for special names like ":L:" use as-is
                if len(new_name) <= 3:
                    formatted_name = new_name.rjust(3)
                else:
                    # For longer names, truncate (though this shouldn't normally happen)
                    formatted_name = new_name[:3]
                    print(f"Warning: Truncating residue name '{new_name}' to '{formatted_name}'")

                # Reconstruct the line with new residue name
                modified_line = line[:17] + formatted_name + line[20:]
                modified_lines.append(modified_line)
                rename_count += 1
            else:
                modified_lines.append(line)
        else:
            modified_lines.append(line)

    if rename_count > 0:
        print(f"  Renamed {rename_count} atoms: {old_name} → {new_name}")
    else:
        print(f"  Warning: No atoms found with residue name '{old_name}'")

    return '\n'.join(modified_lines)


def remove_waters_from_content(content: str, format: str) -> str:
    """
    Remove water molecules from structure content.

    Args:
        content: Structure file content
        format: File format ("pdb" or "cif")

    Returns:
        Structure content with waters removed
    """
    if format != "pdb":
        # For CIF format, just return as-is for now
        return content

    lines = content.split('\n')
    filtered_lines = []

    for line in lines:
        # Skip water molecules (HOH) and other common solvent molecules
        if line.startswith(('ATOM', 'HETATM')):
            res_name = line[17:20].strip()
            if res_name in ['HOH', 'WAT', 'H2O', 'SOL', 'TIP3', 'TIP4', 'SPC']:
                continue  # Skip water lines
        filtered_lines.append(line)

    return '\n'.join(filtered_lines)


def filter_chain_from_content(content: str, format: str, chain: str) -> str:
    """
    Filter structure content to keep only the specified chain.

    Args:
        content: Structure file content
        format: File format ("pdb" or "cif")
        chain: Chain identifier to keep (e.g., "A", "E")

    Returns:
        Structure content with only the specified chain
    """
    if format == "pdb":
        lines = content.split('\n')
        filtered_lines = []

        for line in lines:
            if line.startswith(('ATOM', 'HETATM')):
                if len(line) > 21 and line[21] != chain:
                    continue
            filtered_lines.append(line)

        return '\n'.join(filtered_lines)

    elif format == "cif":
        # CIF _atom_site lines are whitespace-delimited.
        # auth_asym_id (chain ID) column index is determined from the header.
        lines = content.split('\n')
        filtered_lines = []
        in_atom_site = False
        column_names = []
        chain_col_idx = None

        for line in lines:
            stripped = line.strip()

            # Detect start of _atom_site loop
            if stripped == 'loop_' and not in_atom_site:
                filtered_lines.append(line)
                in_atom_site = True
                column_names = []
                chain_col_idx = None
                continue

            if in_atom_site:
                if stripped.startswith('_atom_site.'):
                    column_names.append(stripped)
                    if stripped == '_atom_site.auth_asym_id':
                        chain_col_idx = len(column_names) - 1
                    filtered_lines.append(line)
                    continue

                # Data line within _atom_site block
                if column_names and (stripped.startswith('ATOM') or stripped.startswith('HETATM')):
                    if chain_col_idx is not None:
                        fields = stripped.split()
                        if len(fields) > chain_col_idx and fields[chain_col_idx] != chain:
                            continue
                    filtered_lines.append(line)
                    continue

                # End of _atom_site block
                in_atom_site = False
                column_names = []
                chain_col_idx = None

            filtered_lines.append(line)

        return '\n'.join(filtered_lines)

    return content


def extract_sequence_from_structure(content: str, format: str, chain: str = "longest") -> str:
    """
    Extract protein sequence from structure content.

    By default returns the sequence of the longest chain. Structures often contain
    multiple identical or non-identical chains (e.g. dimers, complexes); concatenating
    them would produce an artificially long sequence that inflates memory usage in
    downstream tools such as LigandMPNN.

    Args:
        content: Structure file content
        format: File format ("pdb" or "cif")
        chain: Which chain to extract. "longest" (default) picks the longest chain.
            Specify a chain letter (e.g. "A", "B") to select that chain.

    Returns:
        Protein sequence from the selected chain
    """
    import tempfile

    if format == "cif":
        return extract_sequence_from_cif(content, chain=chain)

    # PDB format
    try:
        # Write content to temporary file for parsing
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
            tmp.write(content)
            tmp.flush()

            # Parse using the existing PDB parser
            atoms = parse_pdb_file(tmp.name)
            sequences_dict = get_protein_sequence(atoms)

            # Clean up temporary file
            os.unlink(tmp.name)

            if not sequences_dict:
                return ""

            # If a specific chain is requested, return it
            if chain != "longest":
                if chain in sequences_dict:
                    return sequences_dict[chain]
                else:
                    available = ', '.join(sorted(sequences_dict.keys()))
                    print(f"Warning: Chain '{chain}' not found in structure. Available chains: {available}. Falling back to longest chain.")

            # Return the longest chain sequence (sequences_dict already contains only
            # standard amino-acid chains; water/ions are excluded by get_protein_sequence)
            candidates = [seq for seq in sequences_dict.values() if len(seq) >= 10]
            if not candidates:
                return max(sequences_dict.values(), key=len)
            return max(candidates, key=len)

    except Exception as e:
        print(f"Warning: Could not extract sequence from PDB - {str(e)}")
        return ""


def extract_sequence_from_cif(content: str, chain: str = "longest") -> str:
    """
    Extract protein sequence from CIF content using BioPython.

    By default returns the sequence of the longest chain. Structures often contain
    multiple identical or non-identical chains (e.g. dimers, complexes); concatenating
    them would produce an artificially long sequence that inflates memory usage in
    downstream tools such as LigandMPNN.

    Args:
        content: CIF file content
        chain: Which chain to extract. "longest" (default) picks the longest chain.
            Specify a chain letter (e.g. "A", "B") to select that chain.

    Returns:
        Protein sequence from the selected chain
    """
    import tempfile

    # Standard amino acid mapping (3-letter to 1-letter)
    aa_map = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }

    try:
        from Bio.PDB import MMCIFParser

        # Write content to temporary file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.cif', delete=False) as tmp:
            tmp.write(content)
            tmp.flush()
            cif_path = tmp.name

        try:
            # Parse CIF file
            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure("structure", cif_path)

            # Extract sequences from all chains
            sequences = {}
            for model in structure:
                for chain in model:
                    chain_id = chain.id
                    residues = []
                    for residue in chain:
                        res_name = residue.get_resname()
                        if res_name in aa_map:
                            res_id = residue.get_id()[1]  # Residue number
                            residues.append((res_id, aa_map[res_name]))

                    if residues:
                        # Sort by residue number and concatenate
                        residues.sort(key=lambda x: x[0])
                        sequences[chain_id] = ''.join([r[1] for r in residues])

            # Clean up
            os.unlink(cif_path)

            if not sequences:
                return ""

            # If a specific chain is requested, return it
            if chain != "longest":
                if chain in sequences:
                    print(f"  Extracted sequence from CIF: {len(sequences[chain])} residues (chain {chain})")
                    return sequences[chain]
                else:
                    available = ', '.join(sorted(sequences.keys()))
                    print(f"Warning: Chain '{chain}' not found in CIF. Available chains: {available}. Falling back to longest chain.")

            # Return the longest chain sequence (sequences already contains only
            # standard amino-acid chains; water/ions are excluded by the aa_map filter above)
            candidates = [seq for seq in sequences.values() if len(seq) >= 10]
            if not candidates:
                longest = max(sequences.values(), key=len)
            else:
                longest = max(candidates, key=len)
            print(f"  Extracted sequence from CIF: {len(longest)} residues (longest of {len(sequences)} chain(s))")
            return longest

        except Exception as e:
            # Clean up on error
            if os.path.exists(cif_path):
                os.unlink(cif_path)
            raise e

    except ImportError:
        print("Warning: BioPython not available for CIF sequence extraction")
        return ""
    except Exception as e:
        print(f"Warning: Could not extract sequence from CIF - {str(e)}")
        return ""


def extract_ligands_from_structure(content: str, format: str) -> List[str]:
    """
    Extract ligand identifiers (3-letter codes) from structure content.

    Args:
        content: Structure file content
        format: File format ("pdb" or "cif")

    Returns:
        List of unique ligand 3-letter codes (e.g., ['ATP', 'GDP'])
    """
    if format != "pdb":
        print("Warning: Ligand extraction not implemented for CIF format")
        return []

    ligands = set()
    lines = content.split('\n')

    # Common non-ligand residues to exclude (solvents, ions, crystallization agents, caps)
    exclude_residues = {
        'HOH', 'WAT', 'H2O', 'SOL', 'TIP3', 'TIP4', 'SPC',  # Water variants
        'NA', 'CL', 'K', 'CA', 'MG', 'ZN', 'MN', 'FE', 'CU', 'NI', 'CO',  # Common ions
        'SO4', 'PO4', 'NO3',  # Anions
        'GOL', 'EDO', 'PEG', 'PGE', 'PE4', 'PE3', 'P6G', 'PG4', '1PE',  # Glycols and PEGs
        'ACT', 'ACE', 'ACY',  # Acetate
        'PYR', 'PYO',  # Pyruvate
        'DMS', 'BME', 'MPD', 'TRS', 'EPE',  # Common solvents
        'NME', 'NH2'  # Caps and common modifications
    }

    for line in lines:
        if line.startswith('HETATM'):
            res_name = line[17:20].strip()
            # Only include if not in exclusion list and is 3 characters
            if res_name and len(res_name) <= 3 and res_name not in exclude_residues:
                ligands.add(res_name)

    return sorted(list(ligands))


def fetch_ligand_smiles_from_rcsb(ligand_code: str) -> Optional[str]:
    """
    Fetch SMILES string for a ligand from RCSB REST API.

    Args:
        ligand_code: 3-letter ligand code (e.g., 'ATP', 'PVY')

    Returns:
        SMILES string or None if not found
    """
    try:
        import requests

        # RCSB REST API endpoint for ligand info
        url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_code}"

        headers = {
            'User-Agent': 'BioPipelines-PDB/1.0 (https://github.com/locbp-uzh/biopipelines)'
        }

        response = requests.get(url, headers=headers, timeout=10)
        response.raise_for_status()

        data = response.json()

        # Extract SMILES from the response
        # RCSB provides SMILES in the 'chem_comp' section
        if 'chem_comp' in data:
            chem_comp = data['chem_comp']
            # Try different possible fields for SMILES
            for field in ['pdbx_smiles_canonical', 'smiles', 'smiles_canonical']:
                if field in chem_comp and chem_comp[field]:
                    print(f"  Found SMILES for {ligand_code}: {chem_comp[field][:50]}...")
                    return chem_comp[field]

        # Alternative location in descriptors
        if 'rcsb_chem_comp_descriptor' in data:
            descriptors = data['rcsb_chem_comp_descriptor']
            if isinstance(descriptors, dict):
                if 'smiles_canonical' in descriptors:
                    print(f"  Found SMILES for {ligand_code}: {descriptors['smiles_canonical'][:50]}...")
                    return descriptors['smiles_canonical']
                elif 'smiles' in descriptors:
                    print(f"  Found SMILES for {ligand_code}: {descriptors['smiles'][:50]}...")
                    return descriptors['smiles']

        print(f"  No SMILES found for {ligand_code} in RCSB response")
        return None

    except Exception as e:
        print(f"  Warning: Could not fetch SMILES for {ligand_code}: {str(e)}")
        return None


def find_local_structure(pdb_id: str, format: str, local_folder: str,
                        repo_pdbs_folder: str) -> Optional[Tuple[str, str]]:
    """
    Find structure file locally with fallback to CIF if PDB not found.

    Priority for PDB format: X.pdb -> X.cif (will convert)
    Priority for CIF format: X.cif -> X.pdb (will convert, though unusual)
    Search locations: local_folder (if given) -> repo_pdbs_folder

    Args:
        pdb_id: PDB identifier (with or without extension)
        format: Requested format ("pdb" or "cif")
        local_folder: Custom local folder (can be None)
        repo_pdbs_folder: Repository PDBs folder

    Returns:
        Tuple of (path to local file, actual format) or None if not found
    """
    # Remove extension if provided
    pdb_id_base = pdb_id.replace('.pdb', '').replace('.cif', '')

    # Determine primary and fallback extensions
    if format == "pdb":
        primary_ext = ".pdb"
        fallback_ext = ".cif"
    else:
        primary_ext = ".cif"
        fallback_ext = ".pdb"

    search_locations = []
    if local_folder:
        search_locations.append(local_folder)
    search_locations.append(repo_pdbs_folder)

    # Try primary format first
    for location in search_locations:
        candidate = os.path.join(location, f"{pdb_id_base}{primary_ext}")
        if os.path.exists(candidate):
            print(f"Found {pdb_id_base} locally: {candidate}")
            return candidate, format

    # Try fallback format
    for location in search_locations:
        candidate = os.path.join(location, f"{pdb_id_base}{fallback_ext}")
        if os.path.exists(candidate):
            fallback_format = "cif" if fallback_ext == ".cif" else "pdb"
            print(f"Found {pdb_id_base} as {fallback_format.upper()} (will convert to {format.upper()}): {candidate}")
            return candidate, fallback_format

    return None


def copy_local_structure(pdb_id: str, custom_id: str, source_path: str,
                        source_format: str, target_format: str, remove_waters: bool,
                        output_folder: str, operations: List[Dict[str, Any]] = None,
                        chain: str = "longest") -> Tuple[bool, str, str, List[Dict[str, str]], Dict[str, Any]]:
    """
    Copy local structure file to output folder with optional format conversion.

    Args:
        pdb_id: PDB identifier
        custom_id: Custom ID for output filename
        source_path: Path to local structure file
        source_format: Format of source file ("pdb" or "cif")
        target_format: Desired output format ("pdb" or "cif")
        remove_waters: Whether to remove water molecules
        output_folder: Directory to save the structure
        operations: List of operations to apply (e.g., rename)

    Returns:
        Tuple of (success: bool, file_path: str, sequence: str, ligands: List[Dict], metadata: dict)
    """
    if operations is None:
        operations = []

    try:
        with open(source_path, 'r') as f:
            content = f.read()

        # Convert format if needed
        needs_conversion = (source_format != target_format)
        if needs_conversion:
            if source_format == "cif" and target_format == "pdb":
                print(f"  Converting CIF to PDB format...")
                content = convert_cif_to_pdb(content)
            elif source_format == "pdb" and target_format == "cif":
                raise NotImplementedError("PDB to CIF conversion not implemented")

        if remove_waters:
            content = remove_waters_from_content(content, target_format)

        # Filter to specific chain if requested
        if chain != "longest":
            content = filter_chain_from_content(content, target_format, chain)

        # Extract ligands BEFORE applying rename operations (to get original CCD codes for SMILES lookup)
        original_ligand_codes = extract_ligands_from_structure(content, target_format)
        rename_mapping = build_rename_mapping(operations) if operations else {}

        # Apply operations (e.g., rename)
        if operations:
            content = apply_operations(content, target_format, operations)

        extension = ".pdb" if target_format == "pdb" else ".cif"
        filename = f"{custom_id}{extension}"
        output_path = os.path.join(output_folder, filename)

        with open(output_path, 'w') as f:
            f.write(content)

        file_size = os.path.getsize(output_path)
        sequence = extract_sequence_from_structure(content, target_format, chain=chain)

        # Build ligands list using original codes for SMILES lookup, renamed codes for output
        ligands = []
        if original_ligand_codes:
            print(f"  Found {len(original_ligand_codes)} ligand(s) in structure: {', '.join(original_ligand_codes)}")
            for original_code in original_ligand_codes:
                # Use original code for SMILES fetch from RCSB
                smiles = fetch_ligand_smiles_from_rcsb(original_code)
                # Use renamed code (if any) for output
                output_code = rename_mapping.get(original_code, original_code)
                ligands.append({
                    'id': f"{custom_id}_{output_code}",
                    'code': output_code,
                    'format': 'smiles' if smiles else '',
                    'smiles': smiles if smiles else '',
                    'ccd': original_code  # Keep original CCD code for reference
                })

        source_info = f"local ({source_format.upper()})"
        if needs_conversion:
            source_info += f" -> converted to {target_format.upper()}"

        metadata = {
            "file_size": file_size,
            "source": source_info,
            "source_path": source_path
        }

        print(f"Successfully processed {pdb_id} as {custom_id}: {file_size} bytes ({source_info})")
        return True, output_path, sequence, ligands, metadata

    except Exception as e:
        error_msg = f"Error processing local file {pdb_id}: {str(e)}"
        print(f"Error: {error_msg}")
        metadata = {
            "error_message": error_msg,
            "source": "local_copy_failed",
            "attempted_path": source_path
        }
        return False, "", "", [], metadata


def download_from_rcsb(pdb_id: str, custom_id: str, format: str, biological_assembly: bool,
                   remove_waters: bool, output_folder: str, repo_pdbs_folder: str,
                   operations: List[Dict[str, Any]] = None,
                   chain: str = "longest") -> Tuple[bool, str, str, List[Dict[str, str]], Dict[str, Any]]:
    """
    Download a single structure from RCSB PDB and save to both PDBs/ and output folder.

    If PDB format is requested but not available, automatically falls back to CIF and converts.

    Args:
        pdb_id: PDB identifier (4 characters)
        custom_id: Custom ID for renaming the structure
        format: File format ("pdb" or "cif")
        biological_assembly: Whether to download biological assembly
        remove_waters: Whether to remove water molecules
        output_folder: Directory to save the structure
        repo_pdbs_folder: Repository PDBs folder for caching
        operations: List of operations to apply (e.g., rename)

    Returns:
        Tuple of (success: bool, file_path: str, sequence: str, ligands: List[Dict], metadata: dict)
    """
    if operations is None:
        operations = []

    pdb_id = pdb_id.upper()

    # Validate RCSB PDB ID format before attempting download
    if len(pdb_id) != 4 or not pdb_id.isalnum():
        error_msg = f"Invalid RCSB PDB ID format: {pdb_id}. RCSB IDs must be 4 alphanumeric characters (e.g., '4UFC'). File not found locally."
        print(f"Error: {error_msg}")
        metadata = {
            "error_message": error_msg,
            "source": "rcsb_invalid_id",
            "attempted_path": "N/A"
        }
        return False, "", "", [], metadata

    # Try primary format first, then fallback if needed
    download_format = format
    fallback_attempted = False

    def attempt_download(download_fmt: str):
        """Helper to attempt download in specified format."""
        extension = ".pdb" if download_fmt == "pdb" else ".cif"

        if download_fmt == "pdb":
            if biological_assembly:
                url = f"https://files.rcsb.org/download/{pdb_id}.pdb1.gz"
            else:
                url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        else:  # cif
            if biological_assembly:
                url = f"https://files.rcsb.org/download/{pdb_id}-assembly1.cif.gz"
            else:
                url = f"https://files.rcsb.org/download/{pdb_id}.cif"

        return url, extension

    url, extension = attempt_download(download_format)

    try:
        print(f"Downloading {pdb_id} from RCSB: {url}")

        # Import requests only when needed for download
        try:
            import requests
        except ImportError:
            error_msg = f"Cannot download from RCSB: 'requests' module not available. Please install with: pip install requests"
            print(f"Error: {error_msg}")
            metadata = {
                "error_message": error_msg,
                "source": "rcsb_missing_dependency",
                "attempted_path": url
            }
            return False, "", "", metadata

        # Download with proper headers
        headers = {
            'User-Agent': 'BioPipelines-PDB/1.0 (https://github.com/locbp-uzh/biopipelines)'
        }

        response = requests.get(url, headers=headers, timeout=30)
        response.raise_for_status()

        # Handle gzipped content
        if url.endswith('.gz'):
            import gzip
            downloaded_content = gzip.decompress(response.content).decode('utf-8')
        else:
            downloaded_content = response.text

        # Save downloaded file to PDBs/ folder for caching BEFORE conversion (using pdb_id and download format)
        os.makedirs(repo_pdbs_folder, exist_ok=True)
        cache_extension = ".pdb" if download_format == "pdb" else ".cif"
        cache_filename = f"{pdb_id}{cache_extension}"
        cache_path = os.path.join(repo_pdbs_folder, cache_filename)
        with open(cache_path, 'w') as f:
            f.write(downloaded_content)
        print(f"Cached to PDBs/ folder: {cache_path}")

        # Convert format if needed
        needs_conversion = (download_format != format)
        if needs_conversion:
            if download_format == "cif" and format == "pdb":
                print(f"  Converting CIF to PDB format...")
                content = convert_cif_to_pdb(downloaded_content)
            elif download_format == "pdb" and format == "cif":
                raise NotImplementedError("PDB to CIF conversion not implemented")
        else:
            content = downloaded_content

        # Remove waters if requested (use target format after conversion)
        if remove_waters:
            content = remove_waters_from_content(content, format)

        # Filter to specific chain if requested
        if chain != "longest":
            content = filter_chain_from_content(content, format, chain)

        # Extract ligands BEFORE applying rename operations (to get original CCD codes for SMILES lookup)
        original_ligand_codes = extract_ligands_from_structure(content, format)
        rename_mapping = build_rename_mapping(operations) if operations else {}

        # Apply operations (e.g., rename)
        if operations:
            content = apply_operations(content, format, operations)

        # Validate file content (basic check - use target format)
        if format == "pdb":
            if not (content.startswith("HEADER") or content.startswith("ATOM") or content.startswith("MODEL") or content.startswith("REMARK")):
                raise ValueError(f"Downloaded file does not appear to be valid PDB format")
        else:  # cif
            if not ("data_" in content or "_entry.id" in content):
                raise ValueError(f"Downloaded file does not appear to be valid CIF format")

        # Save to output folder (using custom_id and target format)
        output_extension = ".pdb" if format == "pdb" else ".cif"
        output_filename = f"{custom_id}{output_extension}"
        output_path = os.path.join(output_folder, output_filename)
        with open(output_path, 'w') as f:
            f.write(content)

        # Get file size
        file_size = os.path.getsize(output_path)

        # Extract sequence from structure
        sequence = extract_sequence_from_structure(content, format, chain=chain)

        # Build ligands list using original codes for SMILES lookup, renamed codes for output
        ligands = []
        if original_ligand_codes:
            print(f"  Found {len(original_ligand_codes)} ligand(s) in structure: {', '.join(original_ligand_codes)}")
            for original_code in original_ligand_codes:
                # Use original code for SMILES fetch from RCSB
                smiles = fetch_ligand_smiles_from_rcsb(original_code)
                # Use renamed code (if any) for output
                output_code = rename_mapping.get(original_code, original_code)
                ligands.append({
                    'id': f"{custom_id}_{output_code}",
                    'code': output_code,
                    'format': 'smiles' if smiles else '',
                    'smiles': smiles if smiles else '',
                    'ccd': original_code  # Keep original CCD code for reference
                })

        source_info = f"rcsb_download ({download_format.upper()})"
        if needs_conversion:
            source_info += f" -> converted to {format.upper()}"

        metadata = {
            "file_size": file_size,
            "source": source_info,
            "url": url
        }

        print(f"Successfully downloaded {pdb_id} as {custom_id}: {file_size} bytes ({source_info})")
        return True, output_path, sequence, ligands, metadata

    except Exception as e:
        # Try fallback to CIF if PDB download failed and we haven't tried fallback yet
        if format == "pdb" and not fallback_attempted and download_format == "pdb":
            print(f"  PDB download failed, trying CIF format as fallback...")
            fallback_attempted = True
            download_format = "cif"
            url, extension = attempt_download(download_format)

            try:
                print(f"  Downloading {pdb_id} from RCSB: {url}")
                import requests
                headers = {'User-Agent': 'BioPipelines-PDB/1.0 (https://github.com/locbp-uzh/biopipelines)'}
                response = requests.get(url, headers=headers, timeout=30)
                response.raise_for_status()

                # Handle gzipped content
                if url.endswith('.gz'):
                    import gzip
                    cif_content = gzip.decompress(response.content).decode('utf-8')
                else:
                    cif_content = response.text

                # Convert CIF to PDB
                print(f"  Converting CIF to PDB format...")
                content = convert_cif_to_pdb(cif_content)

                # Remove waters if requested
                if remove_waters:
                    content = remove_waters_from_content(content, "pdb")

                # Filter to specific chain if requested
                if chain != "longest":
                    content = filter_chain_from_content(content, "pdb", chain)

                # Extract ligands BEFORE applying rename operations (to get original CCD codes for SMILES lookup)
                original_ligand_codes = extract_ligands_from_structure(content, "pdb")
                rename_mapping = build_rename_mapping(operations) if operations else {}

                # Apply operations (e.g., rename)
                if operations:
                    content = apply_operations(content, "pdb", operations)

                # Cache the CIF file
                os.makedirs(repo_pdbs_folder, exist_ok=True)
                cache_path = os.path.join(repo_pdbs_folder, f"{pdb_id}.cif")
                with open(cache_path, 'w') as f:
                    f.write(cif_content)
                print(f"  Cached CIF to PDBs/ folder: {cache_path}")

                # Save PDB to output folder
                output_path = os.path.join(output_folder, f"{custom_id}.pdb")
                with open(output_path, 'w') as f:
                    f.write(content)

                file_size = os.path.getsize(output_path)
                sequence = extract_sequence_from_structure(content, "pdb", chain=chain)

                # Build ligands list using original codes for SMILES lookup, renamed codes for output
                ligands = []
                if original_ligand_codes:
                    print(f"  Found {len(original_ligand_codes)} ligand(s) in structure: {', '.join(original_ligand_codes)}")
                    for original_code in original_ligand_codes:
                        # Use original code for SMILES fetch from RCSB
                        smiles = fetch_ligand_smiles_from_rcsb(original_code)
                        # Use renamed code (if any) for output
                        output_code = rename_mapping.get(original_code, original_code)
                        ligands.append({
                            'id': f"{custom_id}_{output_code}",
                            'code': output_code,
                            'format': 'smiles' if smiles else '',
                            'smiles': smiles if smiles else '',
                            'ccd': original_code  # Keep original CCD code for reference
                        })

                metadata = {
                    "file_size": file_size,
                    "source": "rcsb_download (CIF) -> converted to PDB",
                    "url": url
                }

                print(f"Successfully downloaded {pdb_id} as {custom_id}: {file_size} bytes (CIF fallback -> PDB)")
                return True, output_path, sequence, ligands, metadata

            except Exception as fallback_e:
                print(f"  CIF fallback also failed: {str(fallback_e)}")
                # Fall through to original error handling

        # Handle both requests exceptions and other errors
        error_type = type(e).__name__

        # Try to extract HTTP status if it's a requests exception
        http_status = 'unknown'
        if 'requests' in sys.modules and hasattr(e, 'response') and e.response:
            http_status = getattr(e.response, 'status_code', 'unknown')

        if 'RequestException' in error_type or 'HTTPError' in error_type:
            error_msg = f"HTTP error downloading {pdb_id}: {str(e)}"
            source = f"rcsb_download_failed_{http_status}"
        else:
            error_msg = f"Unexpected error downloading {pdb_id}: {str(e)}"
            source = "rcsb_processing_error"

        print(f"Error: {error_msg}")

        metadata = {
            "error_message": error_msg,
            "source": source,
            "attempted_path": url
        }

        return False, "", "", [], metadata


def resolve_upstream_file(pdb_id: str, upstream_files: List[str],
                         files_contain_wildcards: bool) -> Optional[str]:
    """
    Resolve the source file for a structure from upstream tool output.

    Args:
        pdb_id: The structure ID to resolve
        upstream_files: List of file paths from the upstream tool
        files_contain_wildcards: Whether the paths contain glob patterns

    Returns:
        Resolved file path, or None if not found
    """
    import glob as glob_module

    if not upstream_files:
        return None

    if files_contain_wildcards:
        # Try each pattern
        for pattern in upstream_files:
            expanded = glob_module.glob(pattern)
            if not expanded:
                continue
            # Try to match by ID
            for fp in expanded:
                basename = os.path.splitext(os.path.basename(fp))[0]
                if basename == pdb_id or basename.startswith(f"{pdb_id}_") or basename.startswith(f"{pdb_id}-"):
                    return fp
            # If single pattern and single expansion, use it
            if len(upstream_files) == 1 and len(expanded) == 1:
                return expanded[0]
        return None

    # Direct file list - match by index or by name
    if len(upstream_files) == 1:
        # Single file for all IDs
        if os.path.exists(upstream_files[0]):
            return upstream_files[0]
        return None

    # Try to match by name
    for fp in upstream_files:
        basename = os.path.splitext(os.path.basename(fp))[0]
        if basename == pdb_id:
            return fp

    return None


def fetch_structures(config_data: Dict[str, Any]) -> int:
    """
    Fetch multiple structures with priority-based lookup: local_folder -> PDBs/ -> RCSB download.
    Also handles upstream tool outputs where files already exist at known paths.

    Args:
        config_data: Configuration dictionary with fetch parameters

    Returns:
        Number of failed fetches
    """
    pdb_ids = config_data['pdb_ids']
    custom_ids = config_data.get('custom_ids', pdb_ids)
    format = config_data['format']
    local_folder = config_data.get('local_folder')
    repo_pdbs_folder = config_data['repo_pdbs_folder']
    biological_assembly = config_data.get('biological_assembly', False)
    remove_waters = config_data.get('remove_waters', True)
    output_folder = config_data['output_folder']
    structures_table = config_data['structures_table']
    sequences_table = config_data['sequences_table']
    failed_table = config_data['failed_table']
    compounds_table = config_data.get('compounds_table', os.path.join(output_folder, 'compounds.csv'))
    operations = config_data.get('operations', [])
    chain = config_data.get('chain', 'longest')
    from_upstream = config_data.get('from_upstream', False)
    upstream_files = config_data.get('upstream_files', [])
    upstream_wildcards = config_data.get('upstream_files_contain_wildcards', False)

    if from_upstream:
        print(f"Processing {len(pdb_ids)} structures from upstream tool")
    else:
        print(f"Fetching {len(pdb_ids)} structures in {format.upper()} format")
        print(f"Priority: {'local_folder -> ' if local_folder else ''}PDBs/ -> RCSB download")
    if biological_assembly:
        print("Including biological assemblies")
    if remove_waters:
        print("Water molecules will be removed")
    if operations:
        op_summaries = []
        for op in operations:
            if op.get("op") == "rename":
                op_summaries.append(f"Rename({op['old']} → {op['new']})")
            else:
                op_summaries.append(op.get("op", "unknown"))
        print(f"Operations: {', '.join(op_summaries)}")

    # Create output directory
    os.makedirs(output_folder, exist_ok=True)

    # Track results
    successful_downloads = []
    successful_sequences = []
    all_ligands = []
    failed_downloads = []

    # Fetch each structure
    for i, (pdb_id, custom_id) in enumerate(zip(pdb_ids, custom_ids), 1):
        print(f"\n[{i}/{len(pdb_ids)}] Processing {pdb_id} -> {custom_id}")

        if from_upstream:
            # Resolve file from upstream tool output
            if len(upstream_files) == len(pdb_ids):
                source_path = upstream_files[i - 1]
            else:
                source_path = resolve_upstream_file(pdb_id, upstream_files, upstream_wildcards)

            if source_path and os.path.exists(source_path):
                # Detect source format from extension
                source_format = "cif" if source_path.endswith(".cif") else "pdb"
                success, file_path, sequence, ligands, metadata = copy_local_structure(
                    pdb_id, custom_id, source_path, source_format, format, remove_waters, output_folder, operations, chain=chain
                )
            else:
                error_msg = f"Upstream file not found for '{pdb_id}': {source_path}"
                print(f"Error: {error_msg}")
                success = False
                file_path = ""
                sequence = ""
                ligands = []
                metadata = {
                    "error_message": error_msg,
                    "source": "upstream_not_found",
                    "attempted_path": str(source_path)
                }
        else:
            # Try to find locally first (returns tuple of (path, source_format) or None)
            local_result = find_local_structure(pdb_id, format, local_folder, repo_pdbs_folder)

            if local_result:
                # Copy from local (may need conversion)
                local_path, source_format = local_result
                success, file_path, sequence, ligands, metadata = copy_local_structure(
                    pdb_id, custom_id, local_path, source_format, format, remove_waters, output_folder, operations, chain=chain
                )
            else:
                # Download from RCSB (with automatic CIF fallback if PDB fails)
                print(f"{pdb_id} not found locally, downloading from RCSB")
                success, file_path, sequence, ligands, metadata = download_from_rcsb(
                    pdb_id, custom_id, format, biological_assembly, remove_waters,
                    output_folder, repo_pdbs_folder, operations, chain=chain
                )

        if success:
            successful_downloads.append({
                'id': custom_id,
                'pdb_id': pdb_id,
                'file_path': file_path,
                'format': format,
                'file_size': metadata['file_size'],
                'source': metadata['source']
            })
            if sequence:  # Only add if sequence extraction was successful
                successful_sequences.append({
                    'id': custom_id,
                    'sequence': sequence
                })
            if ligands:  # Add ligands to the collection
                all_ligands.extend(ligands)
        else:
            failed_downloads.append({
                'pdb_id': pdb_id,
                'error_message': metadata['error_message'],
                'source': metadata['source'],
                'attempted_path': metadata['attempted_path']
            })
    
    # Save successful downloads table
    if successful_downloads:
        df_success = pd.DataFrame(successful_downloads)
        df_success.to_csv(structures_table, index=False)
        print(f"\nSuccessful fetches saved: {structures_table} ({len(successful_downloads)} structures)")
    else:
        # Create empty table with proper columns
        empty_df = pd.DataFrame(columns=["id", "pdb_id", "file_path", "format", "file_size", "source"])
        empty_df.to_csv(structures_table, index=False)
        print(f"No successful fetches - created empty table: {structures_table}")

    # Save sequences table
    if successful_sequences:
        df_sequences = pd.DataFrame(successful_sequences)
        df_sequences.to_csv(sequences_table, index=False)
        print(f"Sequences saved: {sequences_table} ({len(successful_sequences)} sequences)")
    else:
        # Create empty sequences table with proper columns
        empty_seq_df = pd.DataFrame(columns=["id", "sequence"])
        empty_seq_df.to_csv(sequences_table, index=False)
        print(f"No sequences extracted - created empty table: {sequences_table}")

    # Save compounds table
    if all_ligands:
        df_compounds = pd.DataFrame(all_ligands)
        df_compounds.to_csv(compounds_table, index=False)
        print(f"Compounds saved: {compounds_table} ({len(all_ligands)} ligands)")
    else:
        # Create empty compounds table with proper columns
        empty_compounds_df = pd.DataFrame(columns=["id", "code", "format", "smiles", "ccd"])
        empty_compounds_df.to_csv(compounds_table, index=False)
        print(f"No ligands found - created empty table: {compounds_table}")

    # Save failed downloads table (always create, even if empty)
    if failed_downloads:
        df_failed = pd.DataFrame(failed_downloads)
        df_failed.to_csv(failed_table, index=False)
        print(f"Failed fetches saved: {failed_table} ({len(failed_downloads)} failures)")
    else:
        # Create empty failed downloads table
        empty_failed_df = pd.DataFrame(columns=["pdb_id", "error_message", "source", "attempted_path"])
        empty_failed_df.to_csv(failed_table, index=False)
        print("No failed fetches")
    
    # Summary
    print(f"\n=== FETCH SUMMARY ===")
    print(f"Requested: {len(pdb_ids)} structures")
    print(f"Successful: {len(successful_downloads)}")
    print(f"Failed: {len(failed_downloads)}")
    print(f"Success rate: {len(successful_downloads)/len(pdb_ids)*100:.1f}%")

    if successful_downloads:
        total_size = sum(item['file_size'] for item in successful_downloads)
        print(f"Total downloaded: {total_size:,} bytes ({total_size/1024/1024:.2f} MB)")

    print(f"Sequences extracted: {len(successful_sequences)}")
    print(f"Ligands extracted: {len(all_ligands)}")
    
    # Log any failed fetches
    if failed_downloads:
        print(f"\nFailed fetches:")
        for failure in failed_downloads:
            print(f"  - {failure['pdb_id']}: {failure['error_message']}")

    # Return the number of failures
    return len(failed_downloads)


def main():
    parser = argparse.ArgumentParser(description='Fetch protein structures with priority-based lookup')
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
    required_params = ['pdb_ids', 'format', 'repo_pdbs_folder', 'output_folder', 'structures_table', 'sequences_table', 'failed_table']
    for param in required_params:
        if param not in config_data:
            print(f"Error: Missing required parameter: {param}")
            sys.exit(1)

    try:
        failed_count = fetch_structures(config_data)

        # Fail if ANY fetches failed
        if failed_count > 0:
            print(f"\nERROR: {failed_count} structure fetch(es) failed")
            print("Pipeline cannot continue with incomplete structure set")
            print("Check failed_downloads.csv for details")
            sys.exit(1)

        print("\nAll structures fetched successfully")

    except Exception as e:
        print(f"Error fetching structures: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()