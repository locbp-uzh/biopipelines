#!/usr/bin/env python3
"""
Runtime helper script for Ligand tool.

Fetches small molecule ligands with priority-based lookup: local_folder -> Ligands/ -> RCSB download.
Downloads ideal SDF files from RCSB and converts them to PDB format.
Downloads are saved to both Ligands/ folder (for reuse) and tool output folder.
"""

import os
import sys
import argparse
import json
import pandas as pd
from datetime import datetime
from typing import Dict, List, Any, Tuple, Optional
from pathlib import Path


def find_local_ligand(ligand_code: str, local_folder: str,
                     repo_ligands_folder: str) -> Optional[str]:
    """
    Find ligand file locally.

    Priority: local_folder (if given) -> repo_ligands_folder -> None

    Args:
        ligand_code: 3-letter ligand code
        local_folder: Custom local folder (can be None)
        repo_ligands_folder: Repository Ligands folder

    Returns:
        Path to local file or None if not found
    """
    search_locations = []

    if local_folder:
        search_locations.append(local_folder)
    search_locations.append(repo_ligands_folder)

    for location in search_locations:
        candidate = os.path.join(location, f"{ligand_code}.pdb")
        if os.path.exists(candidate):
            print(f"Found {ligand_code} locally: {candidate}")
            return candidate

    return None


def copy_local_ligand(ligand_code: str, custom_id: str, source_path: str,
                     output_folder: str) -> Tuple[bool, str, Dict[str, Any]]:
    """
    Copy local ligand file to output folder.

    Args:
        ligand_code: 3-letter ligand code
        custom_id: Custom ID for output filename
        source_path: Path to local ligand file
        output_folder: Directory to save the ligand

    Returns:
        Tuple of (success: bool, file_path: str, metadata: dict)
    """
    try:
        with open(source_path, 'r') as f:
            content = f.read()

        filename = f"{custom_id}.pdb"
        output_path = os.path.join(output_folder, filename)

        with open(output_path, 'w') as f:
            f.write(content)

        file_size = os.path.getsize(output_path)

        # Try to fetch SMILES from RCSB even for local files
        smiles = fetch_ligand_smiles_from_rcsb(ligand_code)

        metadata = {
            "file_size": file_size,
            "source": "local",
            "source_path": source_path,
            "smiles": smiles if smiles else "",
            "format": "smiles" if smiles else ""
        }

        print(f"Successfully copied {ligand_code} as {custom_id}: {file_size} bytes (from local)")
        return True, output_path, metadata

    except Exception as e:
        error_msg = f"Error copying local file {ligand_code}: {str(e)}"
        print(f"Error: {error_msg}")
        metadata = {
            "error_message": error_msg,
            "source": "local_copy_failed",
            "attempted_path": source_path
        }
        return False, "", metadata


def fetch_ligand_smiles_from_rcsb(ligand_code: str) -> Optional[str]:
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

            # Replace residue name with the ligand code
            lines = pdb_content.split('\n')
            pdb_lines = []
            for line in lines:
                if line.startswith(('HETATM', 'ATOM')):
                    # Replace residue name (columns 18-20, 1-indexed) with ligand code
                    # PDB format: columns are 1-indexed, so residue name is at positions 17-19 (0-indexed: 17:20)
                    line = line[:17] + ligand_code.ljust(3) + line[20:]
                pdb_lines.append(line)

            return '\n'.join(pdb_lines)

        finally:
            # Clean up temporary files
            if os.path.exists(sdf_path):
                os.unlink(sdf_path)
            if os.path.exists(pdb_path):
                os.unlink(pdb_path)

    except Exception as e:
        print(f"  Error converting SDF to PDB with OpenBabel: {str(e)}")
        return None


def download_from_rcsb(ligand_code: str, custom_id: str, output_folder: str,
                       repo_ligands_folder: str) -> Tuple[bool, str, Dict[str, Any]]:
    """
    Download a single ligand from RCSB PDB and save to both Ligands/ and output folder.

    Downloads the ideal SDF file and converts it to PDB format.

    Args:
        ligand_code: 3-letter ligand code (e.g., 'ATP')
        custom_id: Custom ID for renaming the ligand
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

        # Import requests only when needed for download
        try:
            import requests
        except ImportError:
            error_msg = f"Cannot download from RCSB: 'requests' module not available. Please install with: pip install requests"
            print(f"Error: {error_msg}")
            metadata = {
                "error_message": error_msg,
                "source": "rcsb_missing_dependency",
                "attempted_path": sdf_url
            }
            return False, "", metadata

        # Download with proper headers
        headers = {
            'User-Agent': 'BioPipelines-Ligand/1.0 (https://gitlab.uzh.ch/locbp/public/biopipelines)'
        }

        response = requests.get(sdf_url, headers=headers, timeout=30)
        response.raise_for_status()

        sdf_content = response.text

        # Validate SDF content (basic check)
        if not sdf_content.strip():
            raise ValueError(f"Downloaded SDF file is empty")

        # Convert SDF to PDB
        pdb_content = convert_sdf_to_pdb(sdf_content, ligand_code)
        if pdb_content is None:
            raise ValueError(f"Failed to convert SDF to PDB format")

        # Save to Ligands/ folder for caching (using ligand_code, not custom_id)
        os.makedirs(repo_ligands_folder, exist_ok=True)
        cache_filename = f"{ligand_code}.pdb"
        cache_path = os.path.join(repo_ligands_folder, cache_filename)
        with open(cache_path, 'w') as f:
            f.write(pdb_content)
        print(f"Cached to Ligands/ folder: {cache_path}")

        # Save to output folder (using custom_id)
        output_filename = f"{custom_id}.pdb"
        output_path = os.path.join(output_folder, output_filename)
        with open(output_path, 'w') as f:
            f.write(pdb_content)

        # Get file size
        file_size = os.path.getsize(output_path)

        # Fetch SMILES from RCSB API
        smiles = fetch_ligand_smiles_from_rcsb(ligand_code)

        metadata = {
            "file_size": file_size,
            "source": "rcsb_download",
            "url": sdf_url,
            "smiles": smiles if smiles else "",
            "format": "smiles" if smiles else ""
        }

        print(f"Successfully downloaded {ligand_code} as {custom_id}: {file_size} bytes")
        return True, output_path, metadata

    except Exception as e:
        # Handle both requests exceptions and other errors
        error_type = type(e).__name__

        # Try to extract HTTP status if it's a requests exception
        http_status = 'unknown'
        if 'requests' in sys.modules and hasattr(e, 'response') and e.response:
            http_status = getattr(e.response, 'status_code', 'unknown')

        if 'RequestException' in error_type or 'HTTPError' in error_type:
            error_msg = f"HTTP error downloading {ligand_code}: {str(e)}"
            source = f"rcsb_download_failed_{http_status}"
        else:
            error_msg = f"Unexpected error downloading {ligand_code}: {str(e)}"
            source = "rcsb_processing_error"

        print(f"Error: {error_msg}")

        metadata = {
            "error_message": error_msg,
            "source": source,
            "attempted_path": sdf_url
        }

        return False, "", metadata


def fetch_ligands(config_data: Dict[str, Any]) -> int:
    """
    Fetch multiple ligands with priority-based lookup: local_folder -> Ligands/ -> RCSB download.

    Args:
        config_data: Configuration dictionary with fetch parameters

    Returns:
        Number of failed fetches
    """
    ligand_codes = config_data['ligand_codes']
    custom_ids = config_data.get('custom_ids', ligand_codes)
    local_folder = config_data.get('local_folder')
    repo_ligands_folder = config_data['repo_ligands_folder']
    output_folder = config_data['output_folder']
    compounds_table = config_data['compounds_table']
    failed_table = config_data['failed_table']

    print(f"Fetching {len(ligand_codes)} ligands")
    print(f"Priority: {'local_folder -> ' if local_folder else ''}Ligands/ -> RCSB download")

    # Create output directory
    os.makedirs(output_folder, exist_ok=True)

    # Track results
    successful_downloads = []
    failed_downloads = []

    # Fetch each ligand
    for i, (ligand_code, custom_id) in enumerate(zip(ligand_codes, custom_ids), 1):
        print(f"\n[{i}/{len(ligand_codes)}] Processing {ligand_code} -> {custom_id}")

        # Try to find locally first
        local_path = find_local_ligand(ligand_code, local_folder, repo_ligands_folder)

        if local_path:
            # Copy from local
            success, file_path, metadata = copy_local_ligand(
                ligand_code, custom_id, local_path, output_folder
            )
        else:
            # Download from RCSB
            print(f"{ligand_code} not found locally, downloading from RCSB")
            success, file_path, metadata = download_from_rcsb(
                ligand_code, custom_id, output_folder, repo_ligands_folder
            )

        if success:
            successful_downloads.append({
                'id': custom_id,
                'code': ligand_code,
                'file_path': file_path,
                'format': metadata.get('format', ''),
                'smiles': metadata.get('smiles', ''),
                'ccd': ligand_code  # CCD code is the 3-letter code
            })
        else:
            failed_downloads.append({
                'ligand_code': ligand_code,
                'error_message': metadata['error_message'],
                'source': metadata['source'],
                'attempted_path': metadata['attempted_path']
            })

    # Save successful downloads table
    if successful_downloads:
        df_success = pd.DataFrame(successful_downloads)
        df_success.to_csv(compounds_table, index=False)
        print(f"\nSuccessful fetches saved: {compounds_table} ({len(successful_downloads)} ligands)")
    else:
        # Create empty table with proper columns
        empty_df = pd.DataFrame(columns=["id", "code", "file_path", "format", "smiles", "ccd"])
        empty_df.to_csv(compounds_table, index=False)
        print(f"No successful fetches - created empty table: {compounds_table}")

    # Save failed downloads table (always create, even if empty)
    if failed_downloads:
        df_failed = pd.DataFrame(failed_downloads)
        df_failed.to_csv(failed_table, index=False)
        print(f"Failed fetches saved: {failed_table} ({len(failed_downloads)} failures)")
    else:
        # Create empty failed downloads table
        empty_failed_df = pd.DataFrame(columns=["ligand_code", "error_message", "source", "attempted_path"])
        empty_failed_df.to_csv(failed_table, index=False)
        print("No failed fetches")

    # Summary
    print(f"\n=== FETCH SUMMARY ===")
    print(f"Requested: {len(ligand_codes)} ligands")
    print(f"Successful: {len(successful_downloads)}")
    print(f"Failed: {len(failed_downloads)}")
    print(f"Success rate: {len(successful_downloads)/len(ligand_codes)*100:.1f}%")

    if successful_downloads:
        # Count how many have SMILES
        with_smiles = sum(1 for item in successful_downloads if item['smiles'])
        print(f"SMILES fetched: {with_smiles}/{len(successful_downloads)}")

    # Log any failed fetches
    if failed_downloads:
        print(f"\nFailed fetches:")
        for failure in failed_downloads:
            print(f"  - {failure['ligand_code']}: {failure['error_message']}")

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
    required_params = ['ligand_codes', 'repo_ligands_folder', 'output_folder', 'compounds_table', 'failed_table']
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
