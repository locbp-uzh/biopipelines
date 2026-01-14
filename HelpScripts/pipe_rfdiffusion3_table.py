#!/usr/bin/env python3
"""
Create RFdiffusion3 results table.

This script processes RFdiffusion3 output to create a standardized CSV table
that tracks design information for pipeline integration.
"""

import argparse
import pandas as pd
import json
import os
import sys


def get_pdb_length(pdb_file):
    """
    Extract actual length from PDB file by counting unique residue numbers.

    Args:
        pdb_file: Path to PDB file

    Returns:
        Integer length or None if file missing/unreadable
    """
    if not os.path.exists(pdb_file):
        return None

    try:
        residues = set()
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith("ATOM"):
                    # Extract residue number from PDB ATOM line
                    # Format: ATOM serial name altLoc resName chainID resSeq iCode x y z ...
                    #         0-3   7-10  13-15 17     18-19   21     23-26  27
                    res_num = int(line[22:26].strip())
                    chain = line[21].strip()
                    residues.add((chain, res_num))

        return len(residues)
    except Exception as e:
        print(f"Warning: Could not parse PDB length from {pdb_file}: {e}")
        return None


def extract_metadata_from_json(json_file):
    """
    Extract metadata from JSON input file.

    Args:
        json_file: Path to JSON configuration file

    Returns:
        Dictionary with metadata (contig, length, hotspots, input)
    """
    if not os.path.exists(json_file):
        print(f"Warning: JSON file not found: {json_file}")
        return {}

    try:
        with open(json_file, 'r') as f:
            config = json.load(f)

        # Extract first (and typically only) entry
        if not config:
            return {}

        # Get first design entry
        design_key = list(config.keys())[0]
        design_config = config[design_key]

        metadata = {
            "contig": design_config.get("contig", ""),
            "length": design_config.get("length", ""),
            "hotspots": json.dumps(design_config.get("select_hotspots", {})),
            "input": design_config.get("input", "")
        }

        return metadata

    except Exception as e:
        print(f"Warning: Could not parse JSON metadata: {e}")
        return {}


def parse_rfd3_log(log_file):
    """
    Parse RFdiffusion3 log file to extract generation times if available.

    Args:
        log_file: Path to log file

    Returns:
        Dictionary mapping design indices to timing info
    """
    if not os.path.exists(log_file):
        return {}

    timing_info = {}

    try:
        with open(log_file, 'r') as f:
            current_design = None
            for line in f:
                # Look for design markers (implementation-specific)
                # This is a placeholder - actual log format may vary
                if "Processing design" in line or "Design" in line:
                    # Try to extract design index
                    pass

                # Look for timing information
                if "time" in line.lower() or "duration" in line.lower():
                    # Try to extract timing
                    pass

    except Exception as e:
        print(f"Warning: Could not parse log file: {e}")

    return timing_info


def main():
    parser = argparse.ArgumentParser(
        description='Create RFdiffusion3 results table'
    )
    parser.add_argument(
        '--output_folder',
        required=True,
        help="RFdiffusion3 output folder"
    )
    parser.add_argument(
        '--json_file',
        required=True,
        help="Path to JSON input file"
    )
    parser.add_argument(
        '--pipeline_name',
        required=True,
        help="Pipeline name for ID generation"
    )
    parser.add_argument(
        '--num_designs',
        type=int,
        required=True,
        help="Number of designs"
    )
    parser.add_argument(
        '--num_models',
        type=int,
        required=True,
        help="Number of models per design"
    )
    parser.add_argument(
        '--table_path',
        required=True,
        help="Output CSV table path"
    )
    parser.add_argument(
        '--design_startnum',
        type=int,
        default=1,
        help="Starting design number"
    )

    args = parser.parse_args()

    # Extract metadata from JSON
    metadata = extract_metadata_from_json(args.json_file)

    # Parse log file for timing (optional)
    timing_info = parse_rfd3_log(
        os.path.join(os.path.dirname(args.output_folder), "Logs",
                     f"{os.path.basename(args.output_folder).split('_')[0]}_RFdiffusion3.log")
    )

    # Build table
    designs = []
    for i in range(args.num_designs):
        for j in range(args.num_models):
            design_num = args.design_startnum + i
            model_num = args.design_startnum + j

            # Conditional naming: include model suffix only if num_models > 1
            if args.num_models > 1:
                structure_id = f"{args.pipeline_name}_d{design_num}_m{model_num}"
            else:
                structure_id = f"{args.pipeline_name}_{design_num}"

            pdb_file = f"{structure_id}.pdb"
            pdb_path = os.path.join(args.output_folder, pdb_file)

            # Get actual length from PDB
            actual_length = get_pdb_length(pdb_path)

            # Determine status
            if os.path.exists(pdb_path):
                status = "success"
            else:
                status = "missing"

            design_entry = {
                "id": structure_id,
                "design": design_num,
                "model": model_num,
                "pdb": pdb_file,
                "contig": metadata.get("contig", ""),
                "length": actual_length if actual_length else metadata.get("length", ""),
                "time": timing_info.get(i, None),
                "status": status
            }

            designs.append(design_entry)

    # Create DataFrame and save
    df = pd.DataFrame(designs)

    # Ensure output directory exists
    os.makedirs(os.path.dirname(args.table_path), exist_ok=True)

    df.to_csv(args.table_path, index=False)

    # Print summary
    print(f"Created results table: {args.table_path}")
    print(f"Total designs: {len(designs)}")
    successful = sum(1 for d in designs if d['status'] == 'success')
    print(f"Successful: {successful}/{len(designs)}")

    if successful < len(designs):
        missing = [d['id'] for d in designs if d['status'] == 'missing']
        print(f"Missing designs: {', '.join(missing)}")


if __name__ == "__main__":
    main()
