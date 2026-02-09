#!/usr/bin/env python3
# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Helper script for ESMFold - executed during SLURM runtime.

Loads ESMFold model via torch.hub and predicts structures from sequences.
Model is cached in shared folder to avoid repeated downloads.
"""

import argparse
import os
import sys
import pandas as pd
import torch


def setup_torch_hub_cache(cache_dir):
    """
    Configure torch.hub to use custom cache directory.

    Args:
        cache_dir: Path to cache directory
    """
    os.environ['TORCH_HOME'] = cache_dir
    torch.hub.set_dir(cache_dir)
    print(f"Torch hub cache set to: {cache_dir}")


def load_esmfold_model(cache_dir):
    """
    Load ESMFold model via torch.hub.

    Args:
        cache_dir: Directory for model cache

    Returns:
        ESMFold model ready for inference
    """
    print("Loading ESMFold model...")

    # Set cache directory
    setup_torch_hub_cache(cache_dir)

    # Load model via torch.hub
    try:
        model = torch.hub.load("facebookresearch/esm:main", "esmfold_v1")
        model = model.eval()

        # Move to GPU if available
        if torch.cuda.is_available():
            model = model.cuda()
            print("Model loaded on GPU")
        else:
            print("Model loaded on CPU (no GPU available)")

        return model
    except Exception as e:
        print(f"ERROR: Failed to load ESMFold model: {e}")
        sys.exit(1)


def predict_structure(model, sequence, seq_id, output_folder, num_recycles=4, chunk_size=None):
    """
    Predict structure for a single sequence.

    Args:
        model: ESMFold model
        sequence: Protein sequence string
        seq_id: Sequence identifier
        output_folder: Output directory
        num_recycles: Number of recycling iterations
        chunk_size: Chunk size for long sequences

    Returns:
        Path to output PDB file
    """
    print(f"Predicting structure for {seq_id} (length: {len(sequence)})")

    try:
        # Configure model parameters
        if hasattr(model, 'set_chunk_size') and chunk_size:
            model.set_chunk_size(chunk_size)

        # Run inference
        with torch.no_grad():
            output = model.infer_pdb(sequence)

        # Save PDB file
        output_file = os.path.join(output_folder, f"{seq_id}.pdb")
        with open(output_file, 'w') as f:
            f.write(output)

        print(f"Structure saved: {output_file}")
        return output_file

    except Exception as e:
        print(f"ERROR: Failed to predict structure for {seq_id}: {e}")
        return None


def main():
    parser = argparse.ArgumentParser(description="Run ESMFold structure prediction")
    parser.add_argument('--input', required=True, help='Input CSV file with id,sequence columns')
    parser.add_argument('--output', required=True, help='Output directory for PDB files')
    parser.add_argument('--model-cache', required=True, help='Model cache directory')
    parser.add_argument('--num-recycles', type=int, default=4, help='Number of recycling iterations')
    parser.add_argument('--chunk-size', type=int, default=None, help='Chunk size for long sequences')

    args = parser.parse_args()

    # Validate input file
    if not os.path.exists(args.input):
        print(f"ERROR: Input file not found: {args.input}")
        sys.exit(1)

    # Create output directory
    os.makedirs(args.output, exist_ok=True)

    # Create model cache directory
    os.makedirs(args.model_cache, exist_ok=True)

    # Load sequences
    print(f"Reading sequences from {args.input}")
    try:
        df = pd.read_csv(args.input)
        if 'id' not in df.columns or 'sequence' not in df.columns:
            print("ERROR: Input CSV must have 'id' and 'sequence' columns")
            sys.exit(1)
    except Exception as e:
        print(f"ERROR: Failed to read input CSV: {e}")
        sys.exit(1)

    print(f"Found {len(df)} sequences to predict")

    # Load ESMFold model
    model = load_esmfold_model(args.model_cache)

    # Predict structures
    success_count = 0
    failed_sequences = []

    for idx, row in df.iterrows():
        seq_id = row['id']
        sequence = row['sequence']

        result = predict_structure(
            model=model,
            sequence=sequence,
            seq_id=seq_id,
            output_folder=args.output,
            num_recycles=args.num_recycles,
            chunk_size=args.chunk_size
        )

        if result:
            success_count += 1
        else:
            failed_sequences.append(seq_id)

    # Report results
    print(f"\nESMFold prediction completed:")
    print(f"  Successful: {success_count}/{len(df)}")

    if failed_sequences:
        print(f"  Failed: {len(failed_sequences)}")
        print(f"  Failed IDs: {', '.join(failed_sequences)}")
        sys.exit(1)
    else:
        print("All structures predicted successfully")
        sys.exit(0)


if __name__ == "__main__":
    main()
