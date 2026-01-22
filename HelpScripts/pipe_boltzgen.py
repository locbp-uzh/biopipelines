#!/usr/bin/env python3
"""
Enhanced helper script for BoltzGen postprocessing.

This script handles:
1. Structure collection from BoltzGen output
2. Sequence extraction from designed structures
3. Table format conversion for BioPipelines compatibility
4. Output validation
"""

import os
import sys
import argparse
import pandas as pd
from pathlib import Path


def extract_sequence_from_structure(structure_file):
    """
    Extract protein sequence from structure file (PDB or CIF).

    Args:
        structure_file: Path to structure file

    Returns:
        str: Protein sequence in single-letter code
    """
    try:
        from Bio import PDB
        from Bio.PDB import MMCIFParser, PDBParser
        from Bio.PDB.Polypeptide import protein_letters_3to1
    except ImportError:
        print("Error: BioPython is required for sequence extraction")
        print("Install with: pip install biopython")
        sys.exit(1)

    # Determine parser based on file extension
    if str(structure_file).endswith('.cif'):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    try:
        structure = parser.get_structure('structure', structure_file)
    except Exception as e:
        print(f"Error parsing structure {structure_file}: {e}")
        return None

    # Extract sequence from first model, first chain
    sequences = []
    for model in structure:
        for chain in model:
            seq = ""
            for residue in chain:
                if PDB.is_aa(residue, standard=True):
                    resname = residue.get_resname()
                    seq += protein_letters_3to1[resname]
            if seq:
                sequences.append(seq)
                break  # Only take first chain with sequence
        if sequences:
            break  # Only take first model

    return sequences[0] if sequences else None


def collect_and_process_structures(final_designs_folder, output_folder, budget):
    """
    Collect structure files and extract sequences.

    Args:
        final_designs_folder: Path to final_ranked_designs/final_<budget>_designs/
        output_folder: BoltzGen output folder
        budget: Budget parameter

    Returns:
        tuple: (structures_list, sequences_df)
    """
    if not os.path.exists(final_designs_folder):
        print(f"Error: Final designs folder not found: {final_designs_folder}")
        return [], None

    # Collect all .cif and .pdb files
    structures = []
    for ext in ['*.cif', '*.pdb']:
        structures.extend(Path(final_designs_folder).glob(ext))

    # Sort by name for consistent ordering
    structures.sort()

    if not structures:
        print("Warning: No structures found in final designs folder")
        return [], None

    print(f"Found {len(structures)} final structures")

    # Extract sequences from structures
    sequences_data = []
    structures_data = []

    for struct_path in structures:
        # Get structure ID from filename (remove extension)
        struct_id = struct_path.stem

        # Extract sequence
        sequence = extract_sequence_from_structure(struct_path)

        if sequence:
            sequences_data.append({
                "id": struct_id,
                "sequence": sequence
            })

            structures_data.append({
                "id": struct_id,
                "pdb": str(struct_path.absolute())
            })
        else:
            print(f"Warning: Could not extract sequence from {struct_path.name}")

    # Create sequences DataFrame
    sequences_df = pd.DataFrame(sequences_data)

    # Create structures DataFrame
    structures_df = pd.DataFrame(structures_data)

    # Save sequences CSV
    sequences_csv = os.path.join(output_folder, "sequences.csv")
    sequences_df.to_csv(sequences_csv, index=False)
    print(f"Saved sequences to: {sequences_csv}")

    # Save structures list
    structures_csv = os.path.join(output_folder, "structures.csv")
    structures_df.to_csv(structures_csv, index=False)
    print(f"Saved structures list to: {structures_csv}")

    # Also save simple text list for backwards compatibility
    structures_list_file = os.path.join(output_folder, 'final_structures.txt')
    with open(structures_list_file, 'w') as f:
        for struct in structures:
            f.write(f"{struct.absolute()}\n")

    return [str(s.absolute()) for s in structures], sequences_df


def process_metrics(metrics_csv, output_folder):
    """
    Process BoltzGen metrics CSV for BioPipelines compatibility.

    Args:
        metrics_csv: Path to BoltzGen metrics CSV
        output_folder: Output folder path

    Returns:
        bool: True if successful
    """
    if not os.path.exists(metrics_csv):
        print(f"Warning: Metrics CSV not found: {metrics_csv}")
        return False

    try:
        # Read BoltzGen metrics
        df = pd.read_csv(metrics_csv)

        if len(df) == 0:
            print(f"Warning: Metrics CSV is empty")
            return False

        print(f"Processed {len(df)} metric entries")

        # BoltzGen metrics are already in good format, just validate columns
        expected_columns = ["id"]
        missing_columns = [col for col in expected_columns if col not in df.columns]

        if missing_columns:
            print(f"Warning: Missing expected columns: {missing_columns}")

        return True

    except Exception as e:
        print(f"Error processing metrics CSV: {e}")
        return False


def validate_outputs(output_folder, budget):
    """
    Validate that all expected outputs were generated.

    Args:
        output_folder: BoltzGen output folder
        budget: Budget parameter

    Returns:
        bool: True if all validations pass
    """
    success = True

    # Check for final designs folder
    final_designs_folder = os.path.join(
        output_folder,
        'final_ranked_designs',
        f'final_{budget}_designs'
    )
    if not os.path.exists(final_designs_folder):
        print(f"Error: Final designs folder not found: {final_designs_folder}")
        success = False

    # Check for metrics CSV
    metrics_csv = os.path.join(
        output_folder,
        'final_ranked_designs',
        f'final_designs_metrics_{budget}.csv'
    )
    if not os.path.exists(metrics_csv):
        print(f"Error: Metrics CSV not found: {metrics_csv}")
        success = False

    # Check for sequences CSV (created by this script)
    sequences_csv = os.path.join(output_folder, "sequences.csv")
    if not os.path.exists(sequences_csv):
        print(f"Warning: Sequences CSV not found (will be created): {sequences_csv}")

    return success


def main():
    """Main entry point for BoltzGen postprocessing."""
    parser = argparse.ArgumentParser(
        description='Enhanced helper script for BoltzGen postprocessing'
    )
    parser.add_argument(
        '--output_folder',
        type=str,
        required=True,
        help='BoltzGen output folder'
    )
    parser.add_argument(
        '--budget',
        type=int,
        required=True,
        help='Budget parameter used in BoltzGen run'
    )
    parser.add_argument(
        '--extract_sequences',
        action='store_true',
        default=True,
        help='Extract sequences from designed structures (default: True)'
    )
    parser.add_argument(
        '--validate',
        action='store_true',
        default=True,
        help='Validate output files were generated correctly (default: True)'
    )

    args = parser.parse_args()

    print("="*60)
    print("BoltzGen Postprocessing")
    print("="*60)

    # Define key paths
    final_designs_folder = os.path.join(
        args.output_folder,
        'final_ranked_designs',
        f'final_{args.budget}_designs'
    )
    metrics_csv = os.path.join(
        args.output_folder,
        'final_ranked_designs',
        f'final_designs_metrics_{args.budget}.csv'
    )

    success = True

    # Extract sequences and collect structures
    if args.extract_sequences:
        print("\n1. Extracting sequences from structures...")
        structures, sequences_df = collect_and_process_structures(
            final_designs_folder,
            args.output_folder,
            args.budget
        )

        if not structures:
            print("Error: No structures found or processed")
            success = False
        else:
            print(f"  ✓ Processed {len(structures)} structures")

    # Process metrics
    print("\n2. Processing metrics...")
    if not process_metrics(metrics_csv, args.output_folder):
        print("  ✗ Metrics processing failed")
        success = False
    else:
        print("  ✓ Metrics processed successfully")

    # Validate outputs
    if args.validate:
        print("\n3. Validating outputs...")
        if not validate_outputs(args.output_folder, args.budget):
            print("  ✗ Validation failed")
            success = False
        else:
            print("  ✓ All outputs validated")

    print("\n" + "="*60)
    if success:
        print("BoltzGen postprocessing completed successfully ✓")
        print("="*60)
        return 0
    else:
        print("BoltzGen postprocessing completed with errors ✗")
        print("="*60)
        return 1


if __name__ == "__main__":
    sys.exit(main())
