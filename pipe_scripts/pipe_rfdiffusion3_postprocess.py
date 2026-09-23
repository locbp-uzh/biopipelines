#!/usr/bin/env python
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

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
import math
import os
import glob
import tempfile
import pandas as pd
from typing import Dict, List, Tuple, Optional
from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Select
from Bio.SeqUtils import seq1


class AllAtoms(Select):
    """Select all atoms for PDB output."""
    def accept_atom(self, atom):
        return True


def find_cif_gz_files(raw_folder: str) -> List[Tuple[str, str, str, int, int]]:
    """
    Find all CIF.gz files in raw folder and extract key, design and model numbers.

    Foundry names each output ``<key>_<design>_model_<model>.cif.gz`` where
    ``<key>`` is the design key from the input JSON. In the multi-PDB case the
    key is the input pdb id; otherwise it is the single design prefix. The key
    is preserved here so the output structure id is ``<key>_d<D>_m<M>`` and
    multi-PDB outputs stay distinct and joinable to their parent.

    Args:
        raw_folder: Path to raw output folder

    Returns:
        List of tuples: (cif_gz_path, json_path, key, design_number, model_number)
    """
    # Pattern: *_*_model_*.cif.gz (no "design" word in newer RFD3 versions)
    pattern = os.path.join(raw_folder, "*_model_*.cif.gz")
    cif_files = glob.glob(pattern)

    results = []
    for cif_path in cif_files:
        # Extract key, design and model numbers from filename
        basename = os.path.basename(cif_path)
        # Remove .cif.gz extension
        name_without_ext = basename.replace(".cif.gz", "")

        # Find <key>_<design>_model_<model> pattern.
        # Split by _model_ to get the model number.
        parts = name_without_ext.split("_model_")
        if len(parts) == 2:
            try:
                model_num = int(parts[1])

                # The design number is the last underscore-separated part before
                # _model_; everything before it is the design key.
                prefix_and_design = parts[0]
                design_parts = prefix_and_design.rsplit("_", 1)
                if len(design_parts) == 2:
                    key, design_num = design_parts[0], int(design_parts[1])
                else:
                    print(f"Warning: Could not extract key/design from {basename}")
                    continue

                # Find corresponding JSON file
                json_path = cif_path.replace(".cif.gz", ".json")

                if os.path.exists(json_path):
                    results.append((cif_path, json_path, key, design_num, model_num))
                else:
                    print(f"Warning: JSON file not found for {basename}")
                    results.append((cif_path, None, key, design_num, model_num))
            except ValueError:
                print(f"Warning: Could not parse design/model numbers from {basename}")
                continue
        else:
            print(f"Warning: Unexpected filename format: {basename}")
            continue

    # Sort by key, then design number, then model number
    results.sort(key=lambda x: (x[2], x[3], x[4]))
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


def _pdbio_writes_shifted_columns() -> Optional[bool]:
    """Does this Biopython's PDBIO put segID/element/charge one column early? None if it cannot tell."""
    try:
        import importlib
        import sys
        importlib.import_module("Bio.PDB.PDBIO")
        template = sys.modules["Bio.PDB.PDBIO"]._ATOM_FORMAT_STRING
    except Exception:
        return None
    if "%s%s      %4s" in template:
        return False
    if "%s%s     %4s" in template:
        return True
    return None


def _element_from_atom_name(name_field: str) -> str:
    """The PDB convention: a blank or digit in column 13 means a one-letter element in column 14."""
    name_field = name_field.ljust(4)
    if name_field[0] == " " or name_field[0].isdigit():
        return name_field[1].strip().upper()
    return "".join(c for c in name_field[:2] if c.isalpha()).upper()


def _fix_element_columns(pdb_path: str) -> None:
    """Put segID, element and charge in the columns the PDB spec gives them (73-76, 77-78, 79-80).

    Biopython through 1.86 formats ATOM/HETATM with five spaces after the B-factor where the spec has six, so a two-letter symbol like "SI" lands in 76-77 and a spec-compliant reader takes "I " -> iodine; 1.87 writes the spec layout. The layout is read off the installed PDBIO's format string, falling back to the line length (79 shifted, 80 compliant), so a compliant file is never shifted a second time. A blank element is filled from the atom name by the PDB convention, so " CA " is carbon, not calcium.
    """
    shifted_writer = _pdbio_writes_shifted_columns()

    with open(pdb_path) as fh:
        lines = fh.readlines()

    fixed = []
    for line in lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            fixed.append(line)
            continue
        raw = line.rstrip("\r\n")
        shifted = shifted_writer if shifted_writer is not None else len(raw) == 79
        body = raw.ljust(80)
        if shifted:
            segid, element, charge = body[71:75], body[75:77].strip(), body[77:79].strip()
        else:
            segid, element, charge = body[72:76], body[76:78].strip(), body[78:80].strip()
        if not element:
            element = _element_from_atom_name(body[12:16])
        fixed.append(
            f"{body[:66]}{'':>6}{segid}{element.upper():>2}{charge:>2}".rstrip() + "\n"
        )

    with open(pdb_path, "w") as fh:
        fh.writelines(fixed)


# Covalent radii (Angstrom) for the elements a design PDB can contain.
# Cordero et al. 2008; the 0.45 A slack is the usual perception tolerance.
_COVALENT_RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07,
    "S": 1.05, "CL": 1.02, "BR": 1.20, "I": 1.39, "SI": 1.11, "SE": 1.20,
    "B": 0.84, "FE": 1.32, "ZN": 1.22, "MG": 1.41, "CA": 1.76, "NA": 1.66,
    "K": 2.03, "MN": 1.39, "CU": 1.32, "NI": 1.24, "CO": 1.26,
}
_BOND_TOLERANCE = 0.45


def _add_ligand_conect(pdb_path: str) -> int:
    """Append CONECT records for intra-ligand (HETATM) bonds.

    Without explicit connectivity a viewer must perceive ligand bonds by
    distance, and the cutoff it uses is element-dependent and version-
    dependent. For si-rhodamine the Si-C bonds span 1.785-1.892 A — all
    chemically normal — but the two Si-methyls are the longest of the four,
    so a viewer whose cutoff lands near 1.85 A drops exactly those two and
    keeps the ring bonds. CONECT records remove that ambiguity: the bonds
    are stated rather than guessed.

    Only HETATM-HETATM bonds are emitted. Protein connectivity is implied by
    residue templates and needs no CONECT; cross-linking a ligand to the
    protein would assert a covalent bond that may not exist.

    Returns the number of CONECT records written.
    """
    lines = open(pdb_path).read().splitlines()
    het = []
    for i, l in enumerate(lines):
        if l.startswith("HETATM"):
            b = l.ljust(80)
            el = b[76:78].strip().upper() or "".join(
                c for c in b[12:16].strip() if c.isalpha()
            )[:2].upper()
            het.append({
                "serial": int(b[6:11]), "el": el,
                "xyz": (float(b[30:38]), float(b[38:46]), float(b[46:54])),
            })
    if len(het) < 2:
        return 0

    bonds = {a["serial"]: [] for a in het}
    for i in range(len(het)):
        for j in range(i + 1, len(het)):
            a, b = het[i], het[j]
            ra = _COVALENT_RADII.get(a["el"])
            rb = _COVALENT_RADII.get(b["el"])
            if ra is None or rb is None:
                continue
            d = math.dist(a["xyz"], b["xyz"])
            if 0.4 < d <= ra + rb + _BOND_TOLERANCE:
                bonds[a["serial"]].append(b["serial"])
                bonds[b["serial"]].append(a["serial"])

    records = []
    for serial in sorted(bonds):
        partners = sorted(bonds[serial])
        # PDB allows at most four partners per CONECT line; wrap the rest.
        for k in range(0, len(partners), 4):
            chunk = partners[k:k + 4]
            records.append(
                "CONECT" + f"{serial:5d}" + "".join(f"{p:5d}" for p in chunk)
            )
    if not records:
        return 0

    out = [l for l in lines if not l.startswith(("CONECT", "END", "MASTER"))]
    out += records + ["END"]
    with open(pdb_path, "w") as fh:
        fh.write("\n".join(out) + "\n")
    return len(records)


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

        _fix_element_columns(pdb_path)
        _add_ligand_conect(pdb_path)

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
            "ligand_clashes": metrics.get("n_clashing.ligand_clashes"),
            "ligand_min_distance": metrics.get("n_clashing.ligand_min_distance"),
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


def extract_sequence_from_pdb(pdb_path: str) -> Dict[str, str]:
    """
    Extract amino acid sequences from a PDB file.

    Args:
        pdb_path: Path to PDB file

    Returns:
        Dictionary mapping chain IDs to sequences
    """
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", pdb_path)

        sequences = {}
        for model in structure:
            for chain in model:
                chain_id = chain.get_id()
                residues = []
                for residue in chain:
                    # Skip hetero atoms (water, ligands, etc.)
                    if residue.get_id()[0] != ' ':
                        continue
                    resname = residue.get_resname()
                    try:
                        one_letter = seq1(resname)
                        residues.append(one_letter)
                    except KeyError:
                        # Unknown residue, skip
                        continue
                if residues:
                    sequences[chain_id] = ''.join(residues)
            break  # Only process first model

        return sequences
    except Exception as e:
        print(f"Error extracting sequence from {pdb_path}: {e}")
        return {}


def main():
    parser = argparse.ArgumentParser(description="Post-process RFdiffusion3 CIF.gz outputs")
    parser.add_argument("--raw_folder", required=True, help="Raw output folder with CIF.gz files")
    parser.add_argument("--output_folder", required=True, help="Final output folder for PDB files")
    parser.add_argument("--num_designs", type=int, required=True, help="Expected number of designs")
    parser.add_argument("--num_models", type=int, required=True, help="Expected number of models per design")
    parser.add_argument("--design_startnum", type=int, required=True, help="Starting design number")
    parser.add_argument("--metrics_csv", required=True, help="Path to metrics CSV output")
    parser.add_argument("--specifications_csv", required=True, help="Path to specifications CSV output")
    parser.add_argument("--sequences_csv", required=True, help="Path to sequences CSV output")

    args = parser.parse_args()

    print(f"Processing RFdiffusion3 outputs from: {args.raw_folder}")
    print(f"Output folder: {args.output_folder}")
    print(f"Expected designs: {args.num_designs}")
    print(f"Expected models per design: {args.num_models}")

    # Find all CIF.gz files
    cif_files = find_cif_gz_files(args.raw_folder)

    if not cif_files:
        print(f"ERROR: No CIF.gz files found matching pattern in {args.raw_folder}")
        exit(1)

    print(f"Found {len(cif_files)} CIF.gz files")

    # Expected count scales with the number of design keys (one per input PDB
    # in the multi-PDB case, or one for the single/de-novo case), each
    # producing num_designs x num_models outputs.
    num_keys = len({c[2] for c in cif_files})
    expected_total = num_keys * args.num_designs * args.num_models
    if len(cif_files) != expected_total:
        print(f"WARNING: Found {len(cif_files)} files but expected {expected_total} "
              f"({num_keys} key(s) x {args.num_designs} designs x {args.num_models} models)")

    # Create temporary directory for decompressed CIF files
    temp_dir = tempfile.mkdtemp(prefix="rfd3_postprocess_")

    try:
        # Process each CIF file
        metrics_data = []
        specs_data = []
        sequences_data = []
        success_count = 0

        for idx, (cif_gz_path, json_path, key, design_num, model_num) in enumerate(cif_files):
            # With n_batches parameter, designs are numbered sequentially: 0, 1, 2, ...
            # Models within each design are also numbered sequentially: 0, 1, 2, ...
            output_design = args.design_startnum + design_num
            output_model = args.design_startnum + model_num

            # Name from the foundry design key (the input pdb id in the multi-PDB
            # case, or the single design prefix). Include the model suffix only
            # when num_models > 1.
            if args.num_models > 1:
                structure_id = f"{key}_d{output_design}_m{output_model}"
            else:
                structure_id = f"{key}_{output_design}"

            print(f"\nProcessing {design_num}_model_{model_num} -> {structure_id}")

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

            # Extract sequence from PDB
            print(f"  Extracting sequence from PDB...")
            chain_sequences = extract_sequence_from_pdb(pdb_path)
            for chain_id, sequence in chain_sequences.items():
                seq_id = f"{structure_id}_{chain_id}" if len(chain_sequences) > 1 else structure_id
                sequences_data.append({
                    "id": seq_id,
                    "source_id": structure_id,
                    "source_pdb": pdb_path,
                    "chain": chain_id,
                    "sequence": sequence,
                    "length": len(sequence)
                })
                print(f"    Chain {chain_id}: {len(sequence)} residues")

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

        if sequences_data:
            print(f"\nCreating sequences CSV: {args.sequences_csv}")
            sequences_df = pd.DataFrame(sequences_data)
            # Reorder columns
            cols = ["id", "source_id", "source_pdb", "chain", "sequence", "length"]
            sequences_df = sequences_df[cols]
            sequences_df.to_csv(args.sequences_csv, index=False)
            print(f"  ✓ Saved {len(sequences_df)} sequences")

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
