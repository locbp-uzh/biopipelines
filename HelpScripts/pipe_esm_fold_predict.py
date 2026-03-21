# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
ESMFold prediction via the official ESM library (sokrypton/ColabFold approach).

Used for Colab/pip installations. Downloads the ESMFold model weights directly
and runs inference using model.infer() + model.output_to_pdb().
Reads a FASTA file and writes one PDB per sequence to the output directory.
"""

import argparse
import os
import sys
import torch
import gc


def parse_fasta(fasta_path):
    """Parse a FASTA file into (id, sequence) pairs."""
    sequences = []
    current_id = None
    current_seq = []

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id is not None:
                    sequences.append((current_id, ''.join(current_seq)))
                current_id = line[1:].split()[0]  # Take first word after >
                current_seq = []
            elif line:
                current_seq.append(line)

    if current_id is not None:
        sequences.append((current_id, ''.join(current_seq)))

    return sequences


def main():
    parser = argparse.ArgumentParser(description='Run ESMFold prediction via official ESM library')
    parser.add_argument('--fasta', type=str, required=True, help='Input FASTA file')
    parser.add_argument('--output', type=str, required=True, help='Output directory for PDB files')
    parser.add_argument('--num-recycles', type=int, default=3, help='Number of recycles')
    parser.add_argument('--chunk-size', type=int, default=None, help='Chunk size for memory reduction')
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    # Parse sequences
    sequences = parse_fasta(args.fasta)
    if not sequences:
        print("ERROR: No sequences found in FASTA file", file=sys.stderr)
        sys.exit(1)

    # Load model (downloaded during install step)
    model_name = "esmfold.model"
    print(f"Loading ESMFold model from {model_name}...")
    model = torch.load(model_name, weights_only=False)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    model.eval().to(device).requires_grad_(False)

    if device == "cuda":
        torch.backends.cuda.matmul.allow_tf32 = True

    print(f"Using {args.num_recycles} recycles")

    # Predict structures
    failed = []
    completed = 0

    for seq_id, sequence in sequences:
        try:
            print(f"Predicting: {seq_id} ({len(sequence)} residues)")

            # Set chunk size per sequence (optimized for T4 GPUs)
            chunk_size = args.chunk_size
            if chunk_size is None:
                if len(sequence) > 700:
                    chunk_size = 64
                else:
                    chunk_size = 128
            model.set_chunk_size(chunk_size)

            if device == "cuda":
                torch.cuda.empty_cache()

            with torch.no_grad():
                output = model.infer(sequence, num_recycles=args.num_recycles)

            pdb_str = model.output_to_pdb(output)[0]

            output_path = os.path.join(args.output, f"{seq_id}.pdb")
            with open(output_path, 'w') as f:
                f.write(pdb_str)

            # Report mean pLDDT
            mean_plddt = output["plddt"][0, :, 1].mean().item()
            print(f"  -> {seq_id}.pdb (mean pLDDT: {mean_plddt:.1f})")
            completed += 1

            # Free memory between predictions
            del output
            gc.collect()
            if device == "cuda":
                torch.cuda.empty_cache()

        except Exception as e:
            print(f"WARNING: {seq_id} failed: {e}", file=sys.stderr)
            failed.append(seq_id)

    print(f"\nCompleted: {completed}/{len(sequences)}")
    if failed:
        print(f"Failed: {failed}", file=sys.stderr)
    if completed == 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
