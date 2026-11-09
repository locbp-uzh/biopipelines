#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Run Foundry's `mpnn` once per structure and collect its designs into the `sequences` and `designs` streams.

`mpnn --config_json <file>` takes a JSON with engine-level keys (model_type, checkpoint_path, is_legacy_weights,
out_directory, write_fasta, write_structures) and an `inputs` list of per-structure settings. It writes
`<out>/<name>.fa`, appending, with one record per design headed `<name>_b<batch>_d<design>, key=value, ...`, and
one `<out>/<name>_b<batch>_d<design>.cif` per design. One call per structure keeps a bad backbone from taking
the rest down with it, since upstream's loop has no per-input recovery.

The `.fa` holds every chain of a design as one unseparated string, so the chains are recovered from the written
CIF: its standard amino-acid residues, per chain, in order. The two must agree in length, or the design is a
failure rather than a guess.

Three artifacts, kept separate on purpose: the `designs` and `sequences` tables carry only designs that were
produced; the local missing table says why a declared design is absent. A structure that failed is in the missing
table once per design it should have produced.
"""

import argparse
import csv
import json
import os
import re
import subprocess
import sys
from typing import Dict, List, Optional, Tuple

import yaml

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

SEQUENCE_COLUMNS = ["id", "design", "structures.id", "source_pdb", "chain", "sequence", "sample", "T", "seed",
                    "overall_confidence", "ligand_confidence", "seq_rec", "ligand_seq_rec"]
DESIGN_COLUMNS = ["id", "structures.id", "source_pdb", "n_chains", "overall_confidence",
                  "ligand_confidence", "seq_rec", "ligand_seq_rec"]

# FASTA header key -> table column. Keys a release does not write simply stay empty.
HEADER_KEYS = {
    "confidence": "overall_confidence",
    "ligand_interface_confidence": "ligand_confidence",
    "sequence_recovery": "seq_rec",
    "ligand_interface_sequence_recovery": "ligand_seq_rec",
}

HEADER_NAME = re.compile(r"^(?P<name>.+)_b(?P<batch>\d+)_d(?P<design>\d+)$")

# The residues a sequence record can hold. Anything else in a chain is context (ligand, nucleic acid, ion).
AMINO_ACIDS = frozenset((
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
    "PRO", "SER", "THR", "TRP", "TYR", "VAL", "UNK"))

# Per-input keys passed straight through from the wrapper's config, under upstream's own spelling.
PASSTHROUGH_KEYS = ("seed", "temperature", "structure_noise", "designed_residues", "fixed_residues",
                    "designed_chains", "fixed_chains", "omit", "bias")

CIF_TOKEN = re.compile(r"""'[^']*'(?=\s|$)|"[^"]*"(?=\s|$)|\S+""")


def build_config(config: Dict, structure_path: str, name: str, out_dir: str) -> Dict:
    """The JSON handed to `mpnn --config_json` for one structure."""
    entry = {
        "structure_path": structure_path,
        "name": name,
        "batch_size": config["num_sequences"],
        "number_of_batches": config["num_batches"],
    }
    entry.update({key: config[key] for key in PASSTHROUGH_KEYS if key in config})
    engine = {
        "model_type": config["model_type"],
        "is_legacy_weights": config["is_legacy_weights"],
        "out_directory": out_dir,
        "write_fasta": True,
        "write_structures": True,
        "inputs": [entry],
    }
    if config.get("checkpoint_path"):
        engine["checkpoint_path"] = config["checkpoint_path"]
    return engine


def parse_header(header: str) -> Tuple[Optional[Dict], Dict[str, float]]:
    """`<name>_b0_d1, confidence=0.5, ...` -> ({"name", "batch", "design"}, {column: value}).

    The first field is the decorated name; the rest are `key=value` pairs. An unknown key is ignored,
    an unparsable value is left out.
    """
    fields = [field.strip() for field in header.lstrip(">").split(",")]
    match = HEADER_NAME.match(fields[0]) if fields else None
    if not match:
        return None, {}
    ident = {"name": match.group("name"), "batch": int(match.group("batch")),
             "design": int(match.group("design"))}
    metrics: Dict[str, float] = {}
    for field in fields[1:]:
        if "=" not in field:
            continue
        key, value = field.split("=", 1)
        if key.strip() in HEADER_KEYS:
            try:
                metrics[HEADER_KEYS[key.strip()]] = float(value)
            except ValueError:
                continue
    return ident, metrics


def read_fasta(path: str) -> List[Tuple[str, str]]:
    """(header, sequence) records of a FASTA file; a later record with the same header replaces an earlier one."""
    records: Dict[str, str] = {}
    header = None
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                header = line
                records[header] = ""
            elif header is not None and line:
                records[header] += line
    return list(records.items())


def parse_atom_site(text: str) -> Dict[str, List[str]]:
    """The `_atom_site` loop of an mmCIF as {tag: column values}, tags without the category prefix."""
    lines = text.splitlines()
    columns: Dict[str, List[str]] = {}
    index = 0
    while index < len(lines):
        if lines[index].strip() != "loop_":
            index += 1
            continue
        index += 1
        tags = []
        while index < len(lines) and lines[index].startswith("_"):
            tags.append(lines[index].strip())
            index += 1
        if not tags or not tags[0].startswith("_atom_site."):
            continue
        columns = {tag[len("_atom_site."):]: [] for tag in tags}
        names = list(columns)
        while index < len(lines):
            line = lines[index]
            if not line.strip() or line.startswith(("_", "#", "loop_", "data_")):
                break
            values = CIF_TOKEN.findall(line)
            if len(values) == len(names):
                for name, value in zip(names, values):
                    columns[name].append(value.strip("'\""))
            index += 1
        break
    return columns


def chain_lengths(columns: Dict[str, List[str]]) -> List[Tuple[str, int]]:
    """(chain, residue count) in order of first appearance, counting standard amino-acid residues of model 1."""
    chain_key = "auth_asym_id" if "auth_asym_id" in columns else "label_asym_id"
    seq_key = "auth_seq_id" if "auth_seq_id" in columns else "label_seq_id"
    if chain_key not in columns or seq_key not in columns or "label_comp_id" not in columns:
        raise ValueError("the written structure has no _atom_site chain, residue number and component columns")
    models = columns.get("pdbx_PDB_model_num")
    insertion = columns.get("pdbx_PDB_ins_code")
    order: List[str] = []
    residues: Dict[str, set] = {}
    for row, (chain, number, comp) in enumerate(zip(columns[chain_key], columns[seq_key],
                                                    columns["label_comp_id"])):
        if models is not None and models[row] != models[0]:
            continue
        if comp not in AMINO_ACIDS:
            continue
        if chain not in residues:
            order.append(chain)
            residues[chain] = set()
        residues[chain].add((number, insertion[row] if insertion is not None else "?"))
    return [(chain, len(residues[chain])) for chain in order]


def split_chains(sequence: str, lengths: List[Tuple[str, int]]) -> Dict[str, str]:
    """{chain: sequence}, refusing a record whose length is not the sum of the chain lengths."""
    if sum(count for _, count in lengths) != len(sequence):
        raise ValueError(
            f"the designed sequence has {len(sequence)} residues but the written structure has "
            f"{sum(count for _, count in lengths)} standard amino-acid residues over chains "
            f"{'+'.join(chain for chain, _ in lengths) or '(none)'}")
    chains, start = {}, 0
    for chain, count in lengths:
        chains[chain] = sequence[start:start + count]
        start += count
    return chains


def select_chains(chains: Dict[str, str], requested) -> List[Tuple[Optional[str], str]]:
    """The (chain, sequence) rows a design contributes under the wrapper's `chains` setting.

    None expects a single chain and refuses more than one; "all" keeps every chain; a list keeps those.
    A chain is None for the single-chain default, whose ids carry no chain suffix.
    """
    if requested is None:
        if len(chains) != 1:
            raise ValueError(
                f"backbone has {len(chains)} chains ({'+'.join(chains)}); pass chains=\"all\" or "
                f"chains=[...] to say which reach the sequences stream")
        return [(None, next(iter(chains.values())))]
    if requested == "all":
        return list(chains.items())
    missing = [chain for chain in requested if chain not in chains]
    if missing:
        raise ValueError(f"requested chain(s) {missing} are not in this backbone ({'+'.join(chains)})")
    return [(chain, chains[chain]) for chain in requested]


def suffix_chains(requested) -> bool:
    """Only a setting that can emit several chains per design adds a `_<chain>` suffix."""
    return requested == "all" or (isinstance(requested, list) and len(requested) > 1)


def design_ids(name: str, total: int) -> List[str]:
    return [f"{name}_{k}" for k in range(1, total + 1)]


def collect_structure(name: str, source_path: str, out_dir: str, config: Dict) -> Tuple[List[Dict], List[Dict], List[Dict]]:
    """The sequence rows, design rows and missing rows one structure's `mpnn` run produced.

    Reads `<out_dir>/<name>.fa` and the per-design CIFs. A design that cannot be resolved is a failure row
    keyed by its declared id; a design the FASTA lacks altogether likewise.
    """
    batch_size = config["num_sequences"]
    total = batch_size * config["num_batches"]
    requested = config.get("chains")
    source = os.path.splitext(os.path.basename(source_path))[0]

    fasta = os.path.join(out_dir, f"{name}.fa")
    if not os.path.exists(fasta):
        return [], [], [{"id": did, "kind": "failure", "cause": "mpnn wrote no FASTA for this structure"}
                        for did in design_ids(name, total)]

    sequence_rows, design_rows, missing_rows = [], [], []
    produced = set()
    for header, sequence in read_fasta(fasta):
        ident, metrics = parse_header(header)
        if ident is None or ident["name"] != name:
            continue
        sample = ident["batch"] * batch_size + ident["design"] + 1
        did = f"{name}_{sample}"
        if sample > total or did in produced:
            continue
        try:
            cif = os.path.join(out_dir, f"{name}_b{ident['batch']}_d{ident['design']}.cif")
            if not os.path.exists(cif):
                raise ValueError("mpnn wrote no structure for this design, so its chains cannot be told apart")
            with open(cif) as handle:
                lengths = chain_lengths(parse_atom_site(handle.read()))
            chains = split_chains(sequence, lengths)
            picked = select_chains(chains, requested)
        except ValueError as error:
            missing_rows.append({"id": did, "kind": "failure", "cause": str(error)[:200]})
            produced.add(did)
            continue
        produced.add(did)
        design = {"id": did, "structures.id": name, "source_pdb": source, "n_chains": len(chains)}
        design.update(metrics)
        design_rows.append(design)
        for chain, chain_sequence in picked:
            row = {"id": f"{did}_{chain}" if suffix_chains(requested) else did,
                   "design": did, "structures.id": name, "source_pdb": source,
                   "chain": chain or "", "sequence": chain_sequence,
                   "sample": sample, "T": config["temperature"], "seed": config.get("seed")}
            row.update(metrics)
            sequence_rows.append(row)

    for did in design_ids(name, total):
        if did not in produced:
            missing_rows.append({"id": did, "kind": "failure", "cause": "mpnn wrote no record for this design"})
    return sequence_rows, design_rows, missing_rows


def run_mpnn(prefix_argv: List[str], config_path: str) -> Optional[str]:
    """Run `mpnn`; None on success, else the last line it said as the reason."""
    result = subprocess.run(prefix_argv + ["mpnn", "--config_json", config_path],
                            capture_output=True, text=True)
    if result.stdout:
        print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, end="", file=sys.stderr)
    if result.returncode == 0:
        return None
    lines = [line.strip() for line in result.stderr.splitlines() if line.strip()]
    return (lines[-1] if lines else f"mpnn exited with status {result.returncode}")[:200]


def write_rows(path: str, fieldnames: List[str], rows: List[Dict]) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: ("" if row.get(column) is None else row[column]) for column in fieldnames})


def main():
    from biopipelines.biopipelines_io import (load_datastream, iterate_files, step_id_from_table_path,
                                              container_argv_prefix)

    parser = argparse.ArgumentParser(description="Run FoundryMPNN and collect its designs")
    parser.add_argument("--config-yaml", required=True)
    parser.add_argument("--structures-json", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--sequences-csv", required=True)
    parser.add_argument("--designs-csv", required=True)
    parser.add_argument("--local-missing-csv", required=True)
    parser.add_argument("--container-prefix", default="")
    args = parser.parse_args()

    with open(args.config_yaml) as handle:
        config = yaml.safe_load(handle) or {}

    step_id = step_id_from_table_path(args.local_missing_csv)
    prefix_argv = container_argv_prefix(args.container_prefix)
    os.makedirs(args.out_dir, exist_ok=True)
    config_dir = os.path.join(args.out_dir, "_configs")
    os.makedirs(config_dir, exist_ok=True)

    sequence_rows, design_rows, missing_rows = [], [], []
    attempted = failures = 0
    total = config["num_sequences"] * config["num_batches"]

    for name, path in iterate_files(load_datastream(args.structures_json)):
        attempted += 1
        # Upstream appends to <name>.fa, so a retry would otherwise stack its records on the last attempt's.
        stale = os.path.join(args.out_dir, f"{name}.fa")
        if os.path.exists(stale):
            os.remove(stale)
        config_path = os.path.join(config_dir, f"{name}.json")
        with open(config_path, "w") as handle:
            json.dump(build_config(config, path, name, args.out_dir), handle, indent=2)

        reason = run_mpnn(prefix_argv, config_path)
        if reason is not None:
            print(f"WARNING: {name} failed: {reason}", file=sys.stderr)
            failures += 1
            missing_rows.extend({"id": did, "removed_by": step_id, "kind": "failure", "cause": reason}
                                for did in design_ids(name, total))
            continue

        sequences, designs, missing = collect_structure(name, path, args.out_dir, config)
        sequence_rows.extend(sequences)
        design_rows.extend(designs)
        missing_rows.extend({**row, "removed_by": step_id} for row in missing)
        if not designs:
            failures += 1

    write_rows(args.sequences_csv, SEQUENCE_COLUMNS, sequence_rows)
    write_rows(args.designs_csv, DESIGN_COLUMNS, design_rows)
    write_rows(args.local_missing_csv, ["id", "removed_by", "kind", "cause"], missing_rows)

    print(f"Collected {len(design_rows)} design(s) for {attempted} structure(s); "
          f"{len(missing_rows)} design(s) missing")

    # Nothing attempted, or every structure failed: exiting 0 would rubber-stamp a broken step.
    if not attempted:
        print("ERROR: no input entities to process (upstream stream is empty or its files are absent)",
              file=sys.stderr)
        sys.exit(1)
    if failures == attempted:
        print(f"ERROR: all {attempted} attempted structure(s) failed", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
