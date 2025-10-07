#!/usr/bin/env python3
"""
Runtime helper script for PDB tool.

Fetches protein structures with priority-based lookup: local_folder -> PDBs/ -> RCSB download.
Downloads are saved to both PDBs/ folder (for reuse) and tool output folder.
"""

import os
import sys
import argparse
import json
import requests
import pandas as pd
from datetime import datetime
from typing import Dict, List, Any, Tuple
from pathlib import Path

# Import the local PDB parser
sys.path.append(os.path.dirname(__file__))
from pdb_parser import get_protein_sequence, parse_pdb_file


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


def extract_sequence_from_structure(content: str, format: str) -> str:
    """
    Extract protein sequence from structure content using text parsing.

    Args:
        content: Structure file content
        format: File format ("pdb" or "cif")

    Returns:
        Concatenated protein sequence from all chains
    """
    if format != "pdb":
        # For CIF format, return empty sequence for now
        print("Warning: Sequence extraction not implemented for CIF format")
        return ""

    try:
        # Write content to temporary file for parsing
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
            tmp.write(content)
            tmp.flush()

            # Parse using the existing PDB parser
            atoms = parse_pdb_file(tmp.name)
            sequences_dict = get_protein_sequence(atoms)

            # Clean up temporary file
            os.unlink(tmp.name)

            # Concatenate all chain sequences
            return "".join(sequences_dict.values())

    except Exception as e:
        print(f"Warning: Could not extract sequence - {str(e)}")
        return ""


def find_local_structure(pdb_id: str, format: str, local_folder: str,
                        repo_pdbs_folder: str) -> Optional[str]:
    """
    Find structure file locally.

    Priority: local_folder (if given) -> repo_pdbs_folder -> None

    Args:
        pdb_id: PDB identifier (4 characters)
        format: File format ("pdb" or "cif")
        local_folder: Custom local folder (can be None)
        repo_pdbs_folder: Repository PDBs folder

    Returns:
        Path to local file or None if not found
    """
    extension = ".pdb" if format == "pdb" else ".cif"
    search_locations = []

    if local_folder:
        search_locations.append(local_folder)
    search_locations.append(repo_pdbs_folder)

    for location in search_locations:
        candidate = os.path.join(location, f"{pdb_id}{extension}")
        if os.path.exists(candidate):
            print(f"Found {pdb_id} locally: {candidate}")
            return candidate

    return None


def copy_local_structure(pdb_id: str, custom_id: str, source_path: str,
                        format: str, remove_waters: bool,
                        output_folder: str) -> Tuple[bool, str, str, Dict[str, Any]]:
    """
    Copy local structure file to output folder.

    Args:
        pdb_id: PDB identifier
        custom_id: Custom ID for output filename
        source_path: Path to local structure file
        format: File format ("pdb" or "cif")
        remove_waters: Whether to remove water molecules
        output_folder: Directory to save the structure

    Returns:
        Tuple of (success: bool, file_path: str, sequence: str, metadata: dict)
    """
    try:
        with open(source_path, 'r') as f:
            content = f.read()

        if remove_waters:
            content = remove_waters_from_content(content, format)

        extension = ".pdb" if format == "pdb" else ".cif"
        filename = f"{custom_id}{extension}"
        output_path = os.path.join(output_folder, filename)

        with open(output_path, 'w') as f:
            f.write(content)

        file_size = os.path.getsize(output_path)
        sequence = extract_sequence_from_structure(content, format)

        metadata = {
            "file_size": file_size,
            "source": "local",
            "source_path": source_path
        }

        print(f"Successfully copied {pdb_id} as {custom_id}: {file_size} bytes (from local)")
        return True, output_path, sequence, metadata

    except Exception as e:
        error_msg = f"Error copying local file {pdb_id}: {str(e)}"
        print(f"Error: {error_msg}")
        metadata = {
            "error_message": error_msg,
            "source": "local_copy_failed",
            "attempted_path": source_path
        }
        return False, "", "", metadata


def download_from_rcsb(pdb_id: str, custom_id: str, format: str, include_biological_assembly: bool,
                   remove_waters: bool, output_folder: str, repo_pdbs_folder: str) -> Tuple[bool, str, str, Dict[str, Any]]:
    """
    Download a single structure from RCSB PDB and save to both PDBs/ and output folder.

    Args:
        pdb_id: PDB identifier (4 characters)
        custom_id: Custom ID for renaming the structure
        format: File format ("pdb" or "cif")
        include_biological_assembly: Whether to download biological assembly
        remove_waters: Whether to remove water molecules
        output_folder: Directory to save the structure
        repo_pdbs_folder: Repository PDBs folder for caching

    Returns:
        Tuple of (success: bool, file_path: str, sequence: str, metadata: dict)
    """
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
        return False, "", "", metadata

    # Determine URL based on format and assembly
    extension = ".pdb" if format == "pdb" else ".cif"

    if format == "pdb":
        if include_biological_assembly:
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb1.gz"
        else:
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    else:  # cif
        if include_biological_assembly:
            url = f"https://files.rcsb.org/download/{pdb_id}-assembly1.cif.gz"
        else:
            url = f"https://files.rcsb.org/download/{pdb_id}.cif"

    try:
        print(f"Downloading {pdb_id} from RCSB: {url}")

        # Download with proper headers
        headers = {
            'User-Agent': 'BioPipelines-PDB/1.0 (https://gitlab.uzh.ch/locbp/public/biopipelines)'
        }

        response = requests.get(url, headers=headers, timeout=30)
        response.raise_for_status()

        # Handle gzipped content
        if url.endswith('.gz'):
            import gzip
            content = gzip.decompress(response.content).decode('utf-8')
        else:
            content = response.text

        # Remove waters if requested
        if remove_waters:
            content = remove_waters_from_content(content, format)

        # Validate file content (basic check)
        if format == "pdb":
            if not (content.startswith("HEADER") or content.startswith("ATOM") or content.startswith("MODEL")):
                raise ValueError(f"Downloaded file does not appear to be valid PDB format")
        else:  # cif
            if not ("data_" in content or "_entry.id" in content):
                raise ValueError(f"Downloaded file does not appear to be valid CIF format")

        # Save to PDBs/ folder for caching (using pdb_id, not custom_id)
        os.makedirs(repo_pdbs_folder, exist_ok=True)
        cache_filename = f"{pdb_id}{extension}"
        cache_path = os.path.join(repo_pdbs_folder, cache_filename)
        with open(cache_path, 'w') as f:
            f.write(content)
        print(f"Cached to PDBs/ folder: {cache_path}")

        # Save to output folder (using custom_id)
        output_filename = f"{custom_id}{extension}"
        output_path = os.path.join(output_folder, output_filename)
        with open(output_path, 'w') as f:
            f.write(content)

        # Get file size
        file_size = os.path.getsize(output_path)

        metadata = {
            "file_size": file_size,
            "source": "rcsb_download",
            "url": url
        }

        # Extract sequence from structure
        sequence = extract_sequence_from_structure(content, format)

        print(f"Successfully downloaded {pdb_id} as {custom_id}: {file_size} bytes")
        return True, output_path, sequence, metadata
        
    except requests.exceptions.RequestException as e:
        error_msg = f"HTTP error downloading {pdb_id}: {str(e)}"
        print(f"Error: {error_msg}")

        # Try to get HTTP status code
        http_status = getattr(e.response, 'status_code', 'unknown') if hasattr(e, 'response') and e.response else 'network_error'

        metadata = {
            "error_message": error_msg,
            "source": f"rcsb_download_failed_{http_status}",
            "attempted_path": url
        }

        return False, "", "", metadata

    except Exception as e:
        error_msg = f"Unexpected error downloading {pdb_id}: {str(e)}"
        print(f"Error: {error_msg}")

        metadata = {
            "error_message": error_msg,
            "source": "rcsb_processing_error",
            "attempted_path": url
        }

        return False, "", "", metadata


def fetch_structures(config_data: Dict[str, Any]) -> int:
    """
    Fetch multiple structures with priority-based lookup: local_folder -> PDBs/ -> RCSB download.

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
    include_biological_assembly = config_data.get('include_biological_assembly', False)
    remove_waters = config_data.get('remove_waters', True)
    output_folder = config_data['output_folder']
    structures_datasheet = config_data['structures_datasheet']
    sequences_datasheet = config_data['sequences_datasheet']
    failed_datasheet = config_data['failed_datasheet']

    print(f"Fetching {len(pdb_ids)} structures in {format.upper()} format")
    print(f"Priority: {'local_folder -> ' if local_folder else ''}PDBs/ -> RCSB download")
    if include_biological_assembly:
        print("Including biological assemblies")
    if remove_waters:
        print("Water molecules will be removed")

    # Create output directory
    os.makedirs(output_folder, exist_ok=True)

    # Track results
    successful_downloads = []
    successful_sequences = []
    failed_downloads = []

    # Fetch each structure
    for i, (pdb_id, custom_id) in enumerate(zip(pdb_ids, custom_ids), 1):
        print(f"\n[{i}/{len(pdb_ids)}] Processing {pdb_id} -> {custom_id}")

        # Try to find locally first
        local_path = find_local_structure(pdb_id, format, local_folder, repo_pdbs_folder)

        if local_path:
            # Copy from local
            success, file_path, sequence, metadata = copy_local_structure(
                pdb_id, custom_id, local_path, format, remove_waters, output_folder
            )
        else:
            # Download from RCSB
            print(f"{pdb_id} not found locally, downloading from RCSB")
            success, file_path, sequence, metadata = download_from_rcsb(
                pdb_id, custom_id, format, include_biological_assembly, remove_waters,
                output_folder, repo_pdbs_folder
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
        else:
            failed_downloads.append({
                'pdb_id': pdb_id,
                'error_message': metadata['error_message'],
                'source': metadata['source'],
                'attempted_path': metadata['attempted_path']
            })
    
    # Save successful downloads datasheet
    if successful_downloads:
        df_success = pd.DataFrame(successful_downloads)
        df_success.to_csv(structures_datasheet, index=False)
        print(f"\nSuccessful fetches saved: {structures_datasheet} ({len(successful_downloads)} structures)")
    else:
        # Create empty datasheet with proper columns
        empty_df = pd.DataFrame(columns=["id", "pdb_id", "file_path", "format", "file_size", "source"])
        empty_df.to_csv(structures_datasheet, index=False)
        print(f"No successful fetches - created empty datasheet: {structures_datasheet}")

    # Save sequences datasheet
    if successful_sequences:
        df_sequences = pd.DataFrame(successful_sequences)
        df_sequences.to_csv(sequences_datasheet, index=False)
        print(f"Sequences saved: {sequences_datasheet} ({len(successful_sequences)} sequences)")
    else:
        # Create empty sequences datasheet with proper columns
        empty_seq_df = pd.DataFrame(columns=["id", "sequence"])
        empty_seq_df.to_csv(sequences_datasheet, index=False)
        print(f"No sequences extracted - created empty datasheet: {sequences_datasheet}")

    # Save failed downloads datasheet (always create, even if empty)
    if failed_downloads:
        df_failed = pd.DataFrame(failed_downloads)
        df_failed.to_csv(failed_datasheet, index=False)
        print(f"Failed fetches saved: {failed_datasheet} ({len(failed_downloads)} failures)")
    else:
        # Create empty failed downloads datasheet
        empty_failed_df = pd.DataFrame(columns=["pdb_id", "error_message", "source", "attempted_path"])
        empty_failed_df.to_csv(failed_datasheet, index=False)
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
    required_params = ['pdb_ids', 'format', 'repo_pdbs_folder', 'output_folder', 'structures_datasheet', 'sequences_datasheet', 'failed_datasheet']
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