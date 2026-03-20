# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
ESMFold prediction via HuggingFace Transformers API.

Used for Colab/pip installations where the fair-esm CLI is not available.
Reads a FASTA file and writes one PDB per sequence to the output directory.
"""

import argparse
import os
import sys
import torch


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


def convert_outputs_to_pdb(output, sequence_length):
    """Convert HuggingFace ESMFold output to PDB string."""
    from transformers.models.esm.openfold_utils.protein import Protein, to_pdb
    from transformers.models.esm.openfold_utils.feats import atom14_to_atom37

    final_atom_positions = atom14_to_atom37(output["positions"][-1], output)
    final_atom_mask = output["atom37_atom_exists"]

    # Build Protein object for PDB conversion
    protein = Protein(
        aatype=output["aatype"][0, :sequence_length].cpu().numpy(),
        atom_positions=final_atom_positions[0, :sequence_length].cpu().numpy(),
        atom_mask=final_atom_mask[0, :sequence_length].cpu().numpy(),
        residue_index=(output["residue_index"][0, :sequence_length] + 1).cpu().numpy(),
        b_factors=output["plddt"][0, :sequence_length, None].repeat(1, 37).cpu().numpy(),
        chain_index=None,
    )

    return to_pdb(protein)


def main():
    parser = argparse.ArgumentParser(description='Run ESMFold prediction via HuggingFace Transformers')
    parser.add_argument('--fasta', type=str, required=True, help='Input FASTA file')
    parser.add_argument('--output', type=str, required=True, help='Output directory for PDB files')
    parser.add_argument('--num-recycles', type=int, default=4, help='Number of recycles')
    parser.add_argument('--chunk-size', type=int, default=None, help='Chunk size for memory reduction')
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    # Parse sequences
    sequences = parse_fasta(args.fasta)
    if not sequences:
        print("ERROR: No sequences found in FASTA file", file=sys.stderr)
        sys.exit(1)

    print(f"Loading ESMFold model from HuggingFace...")
    from transformers import AutoTokenizer, EsmForProteinFolding

    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = model.to(device)
    model.eval()

    # Memory optimizations
    model.esm = model.esm.half()
    if device == "cuda":
        torch.backends.cuda.matmul.allow_tf32 = True

    # Set chunk size for memory reduction
    chunk_size = args.chunk_size
    if chunk_size is None:
        # Auto-detect: use chunk_size=64 for Colab T4 (16GB), None for larger GPUs
        if device == "cuda":
            gpu_mem_gb = torch.cuda.get_device_properties(0).total_mem / 1e9
            if gpu_mem_gb < 20:
                chunk_size = 64
    if chunk_size is not None:
        model.trunk.set_chunk_size(chunk_size)
        print(f"Using chunk size: {chunk_size}")

    # Set number of recycles
    model.config.num_recycles = args.num_recycles
    print(f"Using {args.num_recycles} recycles")

    # Predict structures
    failed = []
    completed = 0

    for seq_id, sequence in sequences:
        try:
            print(f"Predicting: {seq_id} ({len(sequence)} residues)")
            tokenized = tokenizer([sequence], return_tensors="pt", add_special_tokens=False)["input_ids"].to(device)

            with torch.no_grad():
                output = model(tokenized)

            pdb_str = convert_outputs_to_pdb(output, len(sequence))

            output_path = os.path.join(args.output, f"{seq_id}.pdb")
            with open(output_path, 'w') as f:
                f.write(pdb_str)

            # Report mean pLDDT
            mean_plddt = output["plddt"][0, :len(sequence)].mean().item()
            print(f"  -> {seq_id}.pdb (mean pLDDT: {mean_plddt:.1f})")
            completed += 1

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
