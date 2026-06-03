#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

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

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.sele_utils import list_to_sele


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
    Extract per-design-key metadata from the JSON input file.

    The JSON holds one entry per design key (the input pdb id in the multi-PDB
    case, or a single design prefix otherwise). Returns a dict mapping each key
    to its {contig, length, hotspots, input}, so a structure id can be matched
    back to the entry that produced it.

    Args:
        json_file: Path to JSON configuration file

    Returns:
        Dict mapping design_key -> {contig, length, hotspots, input}
    """
    if not os.path.exists(json_file):
        print(f"Warning: JSON file not found: {json_file}")
        return {}

    try:
        with open(json_file, 'r') as f:
            config = json.load(f)

        per_key = {}
        for design_key, design_config in (config or {}).items():
            per_key[design_key] = {
                "contig": design_config.get("contig", ""),
                "length": design_config.get("length", ""),
                "hotspots": json.dumps(design_config.get("select_hotspots", {})),
                "input": design_config.get("input", ""),
            }
        return per_key

    except Exception as e:
        print(f"Warning: Could not parse JSON metadata: {e}")
        return {}


def design_key_for_id(structure_id, known_keys):
    """Match an output structure id back to its JSON design key.

    Output ids are "<key>_d<D>_m<M>" or "<key>_<D>". Pick the longest known key
    that prefixes the id (keys can contain underscores, so prefer the longest
    match to avoid a shorter key shadowing a longer one).
    """
    candidates = [k for k in known_keys if structure_id.startswith(f"{k}_")]
    return max(candidates, key=len) if candidates else None


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
        '--table_path',
        required=True,
        help="Output CSV table path"
    )
    parser.add_argument(
        '--specifications_csv',
        required=False,
        default=None,
        help="Path to specifications CSV for sampled_contig data"
    )

    args = parser.parse_args()

    # Per-design-key metadata (contig/length) from the JSON input.
    metadata_by_key = extract_metadata_from_json(args.json_file)

    # Load specifications for sampled_contig data
    specs_csv = args.specifications_csv
    if specs_csv is None:
        specs_csv = os.path.join(args.output_folder, "tables", "specifications.csv")
    specs_dict = load_specifications(specs_csv)

    # Scan the structures folder for every produced PDB — one prefix per input
    # PDB in the multi-PDB case, one for the single/de-novo case. Each id is
    # "<key>_d<D>_m<M>" or "<key>_<D>", matched back to its JSON entry by key.
    import glob
    structures_dir = os.path.join(args.output_folder, "structures")
    pdb_paths = sorted(glob.glob(os.path.join(structures_dir, "*.pdb")))

    designs = []
    for pdb_path in pdb_paths:
        pdb_file = os.path.basename(pdb_path)
        structure_id = os.path.splitext(pdb_file)[0]

        key = design_key_for_id(structure_id, metadata_by_key.keys())
        meta = metadata_by_key.get(key, {}) if key else {}

        # Parse the design/model numbers off the id suffix.
        m = re.search(r'_d(\d+)_m(\d+)$', structure_id)
        if m:
            design_num, model_num = int(m.group(1)), int(m.group(2))
        else:
            m = re.search(r'_(\d+)$', structure_id)
            design_num = int(m.group(1)) if m else ""
            model_num = ""

        actual_length = get_pdb_length(pdb_path)

        sampled_contig = specs_dict.get(structure_id, "")
        fixed_sele, designed_sele = parse_sampled_contig(
            sampled_contig, actual_length if actual_length else 0
        )

        designs.append({
            "id": structure_id,
            "design": design_num,
            "model": model_num,
            "pdb": pdb_file,
            "fixed": fixed_sele,
            "designed": designed_sele,
            "contig": meta.get("contig", ""),
            "length": actual_length if actual_length else meta.get("length", ""),
            "time": None,
            "status": "success",
        })

    # Create DataFrame and save
    df = pd.DataFrame(designs, columns=[
        "id", "design", "model", "pdb", "fixed", "designed",
        "contig", "length", "time", "status"])

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
