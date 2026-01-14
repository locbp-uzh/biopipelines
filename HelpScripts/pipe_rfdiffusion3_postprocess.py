#!/usr/bin/env python
"""
Post-processing script for RFdiffusion3 outputs.

Processes raw RFdiffusion3 outputs:
1. Decompresses .cif.gz files
2. Converts CIF format to PDB using BioPython
3. Renames to BioPipelines convention ({prefix}_N.pdb)
4. Extracts metrics and specifications from JSON files
5. Creates metrics and specifications CSV files

Usage:
    python pipe_rfdiffusion3_postprocess.py \\
        --raw_folder /path/to/raw_output \\
        --output_folder /path/to/output \\
        --prefix MyPrefix \\
        --num_designs 5 \\
        --design_startnum 1 \\
        --metrics_csv /path/to/metrics.csv \\
        --specifications_csv /path/to/specs.csv
"""

import argparse
import gzip
import json
import os
import glob
import tempfile
import pandas as pd
from typing import Dict, List, Tuple, Optional
from Bio.PDB import MMCIFParser, PDBIO, Select


class AllAtoms(Select):
    """Select all atoms for PDB output."""
    def accept_atom(self, atom):
        return True


def find_cif_gz_files(raw_folder: str) -> List[Tuple[str, str, int, int]]:
    """
    Find all CIF.gz files in raw folder and extract design and model numbers.

    Args:
        raw_folder: Path to raw output folder

    Returns:
        List of tuples: (cif_gz_path, json_path, design_number, model_number)
    """
    # Pattern: *_design_*_model_*.cif.gz
    pattern = os.path.join(raw_folder, "*_design_*_model_*.cif.gz")
    cif_files = glob.glob(pattern)

    results = []
    for cif_path in cif_files:
        # Extract design and model numbers from filename
        basename = os.path.basename(cif_path)
        # Remove .cif.gz extension
        name_without_ext = basename.replace(".cif.gz", "")

        # Find design_N and model_M
        # Pattern: *_design_N_model_M
        parts = name_without_ext.split("_design_")
        if len(parts) == 2:
            # parts[1] should be "N_model_M"
            design_and_model = parts[1].split("_model_")
            if len(design_and_model) == 2:
                try:
                    design_num = int(design_and_model[0])
                    model_num = int(design_and_model[1])

                    # Find corresponding JSON file
                    json_path = cif_path.replace(".cif.gz", ".json")

                    if os.path.exists(json_path):
                        results.append((cif_path, json_path, design_num, model_num))
                    else:
                        print(f"Warning: JSON file not found for {basename}")
                        results.append((cif_path, None, design_num, model_num))
                except ValueError:
                    print(f"Warning: Could not parse design/model numbers from {basename}")
                    continue
            else:
                print(f"Warning: Unexpected filename format (model): {basename}")
                continue
        else:
            print(f"Warning: Unexpected filename format (design): {basename}")
            continue

    # Sort by design number, then model number
    results.sort(key=lambda x: (x[2], x[3]))
    return results


def decompress_cif(cif_gz_path: str, temp_dir: str) -> str:
    """
    Decompress CIF.gz file to temporary directory.

    Args:
        cif_gz_path: Path to compressed CIF file
        temp_dir: Temporary directory for decompression

    Returns:
        Path to decompressed CIF file
    """
    basename = os.path.basename(cif_gz_path).replace(".gz", "")
    cif_path = os.path.join(temp_dir, basename)

    with gzip.open(cif_gz_path, 'rb') as f_in:
        with open(cif_path, 'wb') as f_out:
            f_out.write(f_in.read())

    return cif_path


def convert_cif_to_pdb(cif_path: str, pdb_path: str) -> bool:
    """
    Convert CIF file to PDB format using BioPython.

    Args:
        cif_path: Path to input CIF file
        pdb_path: Path to output PDB file

    Returns:
        True if successful, False otherwise
    """
    try:
        # Parse CIF file
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("structure", cif_path)

        # Write as PDB
        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_path, select=AllAtoms())

        return True
    except Exception as e:
        print(f"Error converting {cif_path} to PDB: {e}")
        return False


def extract_metrics_from_json(json_path: Optional[str]) -> Dict:
    """
    Extract metrics from RFdiffusion3 JSON output.

    Args:
        json_path: Path to JSON file (or None)

    Returns:
        Dictionary with metrics data
    """
    if json_path is None or not os.path.exists(json_path):
        return {}

    try:
        with open(json_path, 'r') as f:
            data = json.load(f)

        metrics = data.get("metrics", {})
        return {
            "max_ca_deviation": metrics.get("max_ca_deviation"),
            "n_chainbreaks": metrics.get("n_chainbreaks"),
            "n_clashing_interresidue_w_sidechain": metrics.get("n_clashing.interresidue_clashes_w_sidechain"),
            "n_clashing_interresidue_w_backbone": metrics.get("n_clashing.interresidue_clashes_w_backbone"),
            "non_loop_fraction": metrics.get("non_loop_fraction"),
            "loop_fraction": metrics.get("loop_fraction"),
            "helix_fraction": metrics.get("helix_fraction"),
            "sheet_fraction": metrics.get("sheet_fraction"),
            "num_ss_elements": metrics.get("num_ss_elements"),
            "radius_of_gyration": metrics.get("radius_of_gyration"),
            "alanine_content": metrics.get("alanine_content"),
            "glycine_content": metrics.get("glycine_content"),
            "num_residues": metrics.get("num_residues")
        }
    except Exception as e:
        print(f"Error reading JSON file {json_path}: {e}")
        return {}


def extract_specifications_from_json(json_path: Optional[str]) -> Dict:
    """
    Extract specifications from RFdiffusion3 JSON output.

    Args:
        json_path: Path to JSON file (or None)

    Returns:
        Dictionary with specifications data
    """
    if json_path is None or not os.path.exists(json_path):
        return {}

    try:
        with open(json_path, 'r') as f:
            data = json.load(f)

        spec = data.get("specification", {}).get("extra", {})
        return {
            "sampled_contig": spec.get("sampled_contig"),
            "num_tokens_in": spec.get("num_tokens_in"),
            "num_residues_in": spec.get("num_residues_in"),
            "num_chains": spec.get("num_chains"),
            "num_atoms": spec.get("num_atoms"),
            "num_residues": spec.get("num_residues")
        }
    except Exception as e:
        print(f"Error reading JSON file {json_path}: {e}")
        return {}


def main():
    parser = argparse.ArgumentParser(description="Post-process RFdiffusion3 CIF.gz outputs")
    parser.add_argument("--raw_folder", required=True, help="Raw output folder with CIF.gz files")
    parser.add_argument("--output_folder", required=True, help="Final output folder for PDB files")
    parser.add_argument("--prefix", required=True, help="Prefix for output filenames")
    parser.add_argument("--num_designs", type=int, required=True, help="Expected number of designs")
    parser.add_argument("--num_models", type=int, required=True, help="Expected number of models per design")
    parser.add_argument("--design_startnum", type=int, required=True, help="Starting design number")
    parser.add_argument("--metrics_csv", required=True, help="Path to metrics CSV output")
    parser.add_argument("--specifications_csv", required=True, help="Path to specifications CSV output")

    args = parser.parse_args()

    print(f"Processing RFdiffusion3 outputs from: {args.raw_folder}")
    print(f"Output folder: {args.output_folder}")
    print(f"Prefix: {args.prefix}")
    print(f"Expected designs: {args.num_designs}")
    print(f"Expected models per design: {args.num_models}")

    # Find all CIF.gz files
    cif_files = find_cif_gz_files(args.raw_folder)

    if not cif_files:
        print(f"ERROR: No CIF.gz files found matching pattern in {args.raw_folder}")
        exit(1)

    print(f"Found {len(cif_files)} CIF.gz files")

    expected_total = args.num_designs * args.num_models
    if len(cif_files) != expected_total:
        print(f"WARNING: Found {len(cif_files)} files but expected {expected_total} ({args.num_designs} designs × {args.num_models} models)")

    # Create temporary directory for decompressed CIF files
    temp_dir = tempfile.mkdtemp(prefix="rfd3_postprocess_")

    try:
        # Process each CIF file
        metrics_data = []
        specs_data = []
        success_count = 0

        for idx, (cif_gz_path, json_path, design_num, model_num) in enumerate(cif_files):
            # RFdiffusion3 numbers designs as 0, 10, 20, ... (divide by 10)
            # Models are numbered as 0, 1, 2, ... (keep as is)
            output_design = args.design_startnum + (design_num // 10)
            output_model = args.design_startnum + model_num
            structure_id = f"{args.prefix}_d{output_design}_m{output_model}"

            print(f"\nProcessing design_{design_num}_model_{model_num} -> {structure_id}")

            # Decompress CIF.gz
            print(f"  Decompressing CIF.gz...")
            cif_path = decompress_cif(cif_gz_path, temp_dir)

            # Convert to PDB
            pdb_path = os.path.join(args.output_folder, f"{structure_id}.pdb")
            print(f"  Converting to PDB: {pdb_path}")

            if convert_cif_to_pdb(cif_path, pdb_path):
                success_count += 1
                print(f"  ✓ Successfully created {structure_id}.pdb")
            else:
                print(f"  ✗ Failed to convert {structure_id}")
                continue

            # Extract metrics
            print(f"  Extracting metrics from JSON...")
            metrics = extract_metrics_from_json(json_path)
            metrics["id"] = structure_id
            metrics["design"] = output_design
            metrics["model"] = output_model
            metrics_data.append(metrics)

            # Extract specifications
            specs = extract_specifications_from_json(json_path)
            specs["id"] = structure_id
            specs["design"] = output_design
            specs["model"] = output_model
            specs_data.append(specs)

        print(f"\n{'='*60}")
        print(f"Successfully processed {success_count}/{len(cif_files)} designs")

        # Create DataFrames and save CSVs
        if metrics_data:
            print(f"\nCreating metrics CSV: {args.metrics_csv}")
            metrics_df = pd.DataFrame(metrics_data)
            # Reorder columns with id, design, model first
            cols = ["id", "design", "model"] + [c for c in metrics_df.columns if c not in ["id", "design", "model"]]
            metrics_df = metrics_df[cols]
            metrics_df.to_csv(args.metrics_csv, index=False)
            print(f"  ✓ Saved {len(metrics_df)} rows")

        if specs_data:
            print(f"\nCreating specifications CSV: {args.specifications_csv}")
            specs_df = pd.DataFrame(specs_data)
            # Reorder columns with id, design, model first
            cols = ["id", "design", "model"] + [c for c in specs_df.columns if c not in ["id", "design", "model"]]
            specs_df = specs_df[cols]
            specs_df.to_csv(args.specifications_csv, index=False)
            print(f"  ✓ Saved {len(specs_df)} rows")

        print(f"\n{'='*60}")
        print("Post-processing complete!")

    finally:
        # Clean up temporary directory
        import shutil
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
            print(f"\nCleaned up temporary files from {temp_dir}")


if __name__ == "__main__":
    main()
