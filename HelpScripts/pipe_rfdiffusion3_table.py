#!/usr/bin/env python3
# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# This software is freely available for use, modification, and redistribution.
# If you use this software or any derivative of it in published work,
# you must cite the original author and this repository.
"""
Create RFdiffusion3 results table.

This script processes RFdiffusion3 output to create a standardized CSV table
that tracks design information for pipeline integration, including fixed/designed
regions parsed from the sampled_contig in specifications.
"""

import argparse
import pandas as pd
import json
import os
import sys
import re


def list_to_sele(a):
    """Convert list of residue numbers to PyMOL selection string."""
    if not a:
        return ""

    s = ""
    i = 0
    while i < len(a):
        if i > 0:
            s += "+"
        s += f"{a[i]}"
        # Represent consecutive indices with a dash
        if i < len(a) - 1:
            if int(a[i])+1 == int(a[i+1]):
                s += "-"
                j = i + 2
                while j < len(a):
                    if int(a[j]) != int(a[j-1])+1:
                        break
                    j += 1
                i = j - 1
                s += f"{a[i]}"
        i += 1
    return s


def parse_sampled_contig(sampled_contig, total_length):
    """
    Parse sampled_contig to determine fixed and designed residue regions.

    RFdiffusion3 sampled_contig format examples:
    - "17,A9-140" -> 17 designed residues (1-17), then fixed A9-140 mapped to 18-149
    - "A1-50,20,A60-100" -> fixed A1-50 (1-50), 20 designed (51-70), fixed A60-100 (71-111)
    - "160" -> 160 designed residues (de novo design with only length specified)

    Args:
        sampled_contig: The sampled contig string from specifications
        total_length: Total length of the protein (for validation)

    Returns:
        Tuple of (fixed_selection, designed_selection) as PyMOL selection strings
    """
    if not sampled_contig:
        return "", ""

    # Ensure sampled_contig is a string (may be int/float from pandas)
    sampled_contig = str(sampled_contig)

    fixed_residues = []
    designed_residues = []
    current_pos = 1

    # Split by comma to get segments
    segments = [s.strip() for s in sampled_contig.split(',')]

    for segment in segments:
        # Check if segment is a chain reference (e.g., "A9-140" or "A1-50")
        chain_match = re.match(r'^([A-Z])(\d+)-(\d+)$', segment)
        if chain_match:
            # Fixed region from input structure
            chain = chain_match.group(1)
            start = int(chain_match.group(2))
            end = int(chain_match.group(3))
            num_residues = end - start + 1

            # These residues are fixed
            for i in range(num_residues):
                fixed_residues.append(current_pos + i)
            current_pos += num_residues

        # Check if segment is just a number (de novo designed)
        elif segment.isdigit():
            num_residues = int(segment)
            # These residues are designed
            for i in range(num_residues):
                designed_residues.append(current_pos + i)
            current_pos += num_residues

        # Check for range without chain (e.g., "15-25" in original contig means variable length)
        # In sampled_contig this should be resolved to a specific number
        else:
            # Try to parse as a single chain reference without range (e.g., "A50")
            single_match = re.match(r'^([A-Z])(\d+)$', segment)
            if single_match:
                # Single residue from chain - fixed
                fixed_residues.append(current_pos)
                current_pos += 1
            else:
                print(f"Warning: Could not parse segment '{segment}' in sampled_contig")

    return list_to_sele(fixed_residues), list_to_sele(designed_residues)


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


def load_specifications(specs_csv):
    """
    Load specifications CSV to get sampled_contig per structure.

    Args:
        specs_csv: Path to specifications CSV file

    Returns:
        Dictionary mapping structure_id to sampled_contig
    """
    if not os.path.exists(specs_csv):
        print(f"Warning: Specifications CSV not found: {specs_csv}")
        return {}

    try:
        df = pd.read_csv(specs_csv)
        specs_dict = {}
        for _, row in df.iterrows():
            structure_id = row.get('id', '')
            sampled_contig = row.get('sampled_contig', '')
            if structure_id:
                # Ensure sampled_contig is a string (pandas may parse pure numbers as int/float)
                specs_dict[structure_id] = str(sampled_contig) if pd.notna(sampled_contig) else ''
        return specs_dict
    except Exception as e:
        print(f"Warning: Could not load specifications CSV: {e}")
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
    parser.add_argument(
        '--specifications_csv',
        required=False,
        default=None,
        help="Path to specifications CSV for sampled_contig data"
    )

    args = parser.parse_args()

    # Extract metadata from JSON
    metadata = extract_metadata_from_json(args.json_file)

    # Load specifications for sampled_contig data
    specs_csv = args.specifications_csv
    if specs_csv is None:
        specs_csv = os.path.join(args.output_folder, "rfdiffusion3_specifications.csv")
    specs_dict = load_specifications(specs_csv)

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

            # Get sampled_contig for this structure and parse fixed/designed
            sampled_contig = specs_dict.get(structure_id, "")
            fixed_sele, designed_sele = parse_sampled_contig(
                sampled_contig,
                actual_length if actual_length else 0
            )

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
                "fixed": fixed_sele,
                "designed": designed_sele,
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
