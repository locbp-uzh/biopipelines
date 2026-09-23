# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Convert ProteinMPNN / LigandMPNN .fa output into the canonical sequences stream.

A design over a multi-chain backbone is emitted by both models as one record whose chains are joined by a separator ('/' for ProteinMPNN, ':' for LigandMPNN). One row of a `sequences` stream is one polymer chain, so such a record is split into one row per chain here. The design itself keeps its identity in the `designs` stream, which is what a downstream consumer groups the chain rows by.
"""

import argparse

parser = argparse.ArgumentParser(description='Makes csv files from ProteinMPNN/LigandMPNN .fa output')
parser.add_argument('FA_FOLDER', type=str)
parser.add_argument('queries_csv_file', type=str)
parser.add_argument('queries_fasta_file', type=str)
parser.add_argument('-d', '--duplicates', action='store_true',
                    help='Allow duplicate sequences in output')
parser.add_argument('--id-map', type=str, default=None,
                    help='Path to JSON file mapping PDB basenames to stream IDs')
parser.add_argument('--missing-csv', type=str, default=None,
                    help='Path to write missing.csv for removed duplicates')
parser.add_argument('--step-tool-name', type=str, default=None,
                    help='Step and tool name for missing.csv removed_by column (e.g. 005_ProteinMPNN)')
parser.add_argument('--upstream-missing', type=str, default=None,
                    help='Path to upstream missing.csv to propagate')
parser.add_argument('--fill-gaps', type=str, default=None,
                    help='Replace X (unknown/gap residues) with this amino acid (e.g., G for glycine)')
parser.add_argument('--ds-json', type=str, default=None,
                    help='DataStream JSON for runtime id_map generation (alternative to --id-map)')
parser.add_argument('--chains', type=str, default='',
                    help='Which backbone chains reach the sequences stream: "" (expect a single '
                         'chain), "all", or a comma-separated list such as "A,B"')
parser.add_argument('--designs-csv', type=str, default=None,
                    help='Path to write the designs table (one row per design, before chain split)')

# Parse the arguments
args = parser.parse_args()

import os
import sys
import json
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.pdb_parser import parse_pdb_file, STANDARD_RESIDUES

# ProteinMPNN joins chains with '/', LigandMPNN with ':'. Neither is a residue letter, so a record carrying either is multi-chain.
CHAIN_SEPARATORS = "/:"

# protein_mpnn_run.py concatenates chains in this order; a PDB's own record order is only consulted when alphabetical does not fit.
CHAIN_ALPHABET = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"


def split_chains(sequence):
    """The record's chain segments, in the order the model wrote them."""
    for sep in CHAIN_SEPARATORS:
        if sep in sequence:
            return [seg for seg in sequence.split(sep)]
    return [sequence]


_ONE_LETTER = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
               "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
               "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def chain_residue_counts(pdb_path, with_sequences=False):
    """(appearance-ordered chain ids, {chain: residue count}) for the polymer chains of a PDB."""
    atoms = parse_pdb_file(pdb_path)
    order = []
    residues = {}
    for atom in atoms:
        if atom.res_name not in STANDARD_RESIDUES:
            continue
        if atom.chain not in residues:
            residues[atom.chain] = {}
            order.append(atom.chain)
        residues[atom.chain].setdefault(atom.res_num, atom.res_name)
    counts = {c: len(r) for c, r in residues.items()}
    if not with_sequences:
        return order, counts
    sequences = {c: "".join(_ONE_LETTER.get(r[n], "X") for n in sorted(r)) for c, r in residues.items()}
    return order, counts, sequences


def _identity(segments, candidate, sequences):
    return sum(a == b for seg, chain in zip(segments, candidate)
               for a, b in zip(seg, sequences.get(chain, "")))


def _alphabetical(chain_ids):
    return sorted(chain_ids, key=lambda c: (CHAIN_ALPHABET.index(c) if c in CHAIN_ALPHABET else len(CHAIN_ALPHABET), c))


def resolve_chain_order(segments, pdb_path):
    """Chain ids aligned to `segments`, or (None, reason) when the two cannot be reconciled.

    Segment lengths are matched against each chain's residue count rather than its sequence, because the model's three-to-one letter mapping and this parser's disagree on non-standard residues while the counts do not.
    """
    appearance, counts, sequences = chain_residue_counts(pdb_path, with_sequences=True)
    if not counts:
        return None, "no polymer chains found in the input structure"
    lengths = [len(seg) for seg in segments]
    matching = []
    for candidate in (_alphabetical(counts.keys()), appearance):
        if [counts[c] for c in candidate] == lengths and candidate not in matching:
            matching.append(candidate)
    if len(matching) > 1:
        # Equal lengths cannot tell the orders apart; the residues the model kept fixed can.
        scores = [_identity(segments, c, sequences) for c in matching]
        if scores[0] == scores[1]:
            print(f"Warning: {os.path.basename(pdb_path)}: chains {matching[0]} and {matching[1]} "
                  f"are indistinguishable by length and sequence; assuming {matching[0]}",
                  file=sys.stderr)
            return matching[0], None
        return matching[scores.index(max(scores))], None
    if matching:
        return matching[0], None
    return None, (f"cannot match {len(segments)} sequence segment(s) of length {lengths} to chains "
                  f"{ {c: counts[c] for c in _alphabetical(counts.keys())} }")


# Load optional ID map. When built from a DataStream JSON, its keys also gate
# which .fa files are converted, so a filtered upstream stream (ids restricted
# at runtime via ids_expanded) does not pull in filtered-out structures.
id_map = None
pdb_paths = {}
allowed_bases = None
if args.ds_json and os.path.exists(args.ds_json):
    from biopipelines.biopipelines_io import load_datastream, iterate_files
    ds = load_datastream(args.ds_json)
    id_map = {}
    for struct_id, pdb_path in iterate_files(ds):
        pdb_base = os.path.splitext(os.path.basename(pdb_path))[0]
        id_map[pdb_base] = struct_id
        pdb_paths[pdb_base] = pdb_path
    allowed_bases = set(id_map.keys())
elif args.id_map and os.path.exists(args.id_map):
    with open(args.id_map, 'r') as f:
        id_map = json.load(f)

requested_chains = None
if args.chains.strip().lower() == "all":
    requested_chains = "all"
elif args.chains.strip():
    requested_chains = [c.strip() for c in args.chains.split(",") if c.strip()]

# A chain suffix is appended only when the step can emit more than one chain per design, so a single-chain run keeps the ids it has always had.
suffix_chains = requested_chains == "all" or (isinstance(requested_chains, list) and len(requested_chains) > 1)

fa_files = sorted(os.listdir(args.FA_FOLDER))
seen_records = {}       # full multi-chain record -> first design id that carried it
sequence_rows = []      # one row per chain
design_rows = []        # one row per design
dropped = []            # {id, removed_by, kind, cause}


def drop(design_id, cause, kind='failure'):
    print(f"WARNING: {design_id}: {cause}", file=sys.stderr)
    if args.step_tool_name:
        dropped.append({'id': design_id, 'removed_by': args.step_tool_name,
                        'kind': kind, 'cause': cause})


for fa in fa_files:
    if not fa.endswith(".fa"):
        continue
    pdb_base = fa[:-3]
    if allowed_bases is not None and pdb_base not in allowed_bases:
        continue
    mapped_base = id_map.get(pdb_base, pdb_base) if id_map else pdb_base

    with open(os.path.join(args.FA_FOLDER, fa), "r") as handle:
        lines = [line.strip() for line in handle.readlines()]

    chain_order = None
    chain_order_error = None

    # Starts from 2 to skip the native sequence the model echoes back as its first record.
    for i in range(2, len(lines), 2):
        header = lines[i][1:]
        record = lines[i + 1]

        params = {}
        for p in header.split(", "):
            if '=' not in p:
                continue
            key, value = p.split("=", 1)
            # LigandMPNN spells the sample index `id`; ProteinMPNN spells it `sample`.
            params['sample' if key == 'id' else key] = value

        if 'sample' not in params:
            drop(f"{mapped_base}_?", f"FASTA header carries no sample index: {header!r}")
            continue
        design_id = f"{mapped_base}_{params['sample']}"

        if record in seen_records and not args.duplicates:
            drop(design_id, f"Duplicate of {seen_records[record]}", kind='filter')
            continue
        seen_records[record] = design_id

        segments = split_chains(record)

        if len(segments) == 1 and requested_chains is None:
            chains_here = [None]
        else:
            if chain_order is None and chain_order_error is None:
                if pdb_base not in pdb_paths:
                    chain_order_error = ("cannot resolve chain letters without the input structure; "
                                         "the step must pass --ds-json")
                else:
                    chain_order, chain_order_error = resolve_chain_order(segments, pdb_paths[pdb_base])
            if chain_order_error:
                drop(design_id, chain_order_error)
                continue
            if len(segments) != len(chain_order):
                drop(design_id, f"{len(segments)} sequence segment(s) against {len(chain_order)} chain(s)")
                continue
            if requested_chains is None:
                drop(design_id,
                     f"backbone has {len(chain_order)} chains ({'+'.join(chain_order)}); "
                     f"pass chains=\"all\" or chains=[...] to say which reach the sequences stream")
                continue
            chains_here = list(chain_order)

        design_row = {'id': design_id, 'structures.id': mapped_base, 'source_pdb': pdb_base,
                      'n_chains': len(segments)}
        design_row.update({k: v for k, v in params.items() if k != 'sample'})
        design_row['sample'] = params['sample']
        design_rows.append(design_row)

        emitted = 0
        for segment, chain in zip(segments, chains_here):
            if isinstance(requested_chains, list) and chain not in requested_chains:
                continue
            gap_indices = [pos + 1 for pos, aa in enumerate(segment) if aa == "X"]
            sequence = segment
            if gap_indices and args.fill_gaps:
                sequence = sequence.replace("X", args.fill_gaps)
                print(f"  Filled {len(gap_indices)} gap(s) in {design_id} chain {chain} with {args.fill_gaps}")
            row = {
                'id': f"{design_id}_{chain}" if suffix_chains else design_id,
                'design': design_id,
                'structures.id': mapped_base,
                'source_pdb': pdb_base,
                'chain': chain or "",
                'sequence': sequence,
                'gaps': "+".join(str(p) for p in gap_indices) if gap_indices else "",
            }
            row.update({k: v for k, v in params.items() if k != 'sample'})
            row['sample'] = params['sample']
            sequence_rows.append(row)
            emitted += 1

        if emitted == 0:
            drop(design_id, f"none of the requested chains {requested_chains} are in this backbone")
            design_rows.pop()

if sequence_rows:
    df = pd.DataFrame(sequence_rows)
    # `id` and `sequence` lead so the file reads as the sequences table it is.
    lead = [c for c in ('id', 'design', 'structures.id', 'source_pdb', 'chain', 'sequence') if c in df.columns]
    df = df[lead + [c for c in df.columns if c not in lead]]
    os.makedirs(os.path.dirname(args.queries_csv_file), exist_ok=True)
    df.to_csv(args.queries_csv_file, index=False)

    with open(args.queries_fasta_file, "w") as fasta:
        fasta.write("\n".join(f">{r['id']}\n{r['sequence']}" for r in sequence_rows))
    print(f"Wrote {len(sequence_rows)} sequence row(s) from {len(design_rows)} design(s)")
else:
    print("No sequences produced", file=sys.stderr)

if args.designs_csv:
    os.makedirs(os.path.dirname(args.designs_csv), exist_ok=True)
    if design_rows:
        designs_df = pd.DataFrame(design_rows)
        lead = [c for c in ('id', 'structures.id', 'source_pdb', 'n_chains') if c in designs_df.columns]
        designs_df = designs_df[lead + [c for c in designs_df.columns if c not in lead]]
    else:
        designs_df = pd.DataFrame(columns=['id', 'structures.id', 'source_pdb', 'n_chains'])
    designs_df.to_csv(args.designs_csv, index=False)
    print(f"Wrote {len(design_rows)} design row(s) to {args.designs_csv}")

# Write missing.csv if requested
if args.missing_csv:
    upstream_rows = []
    if args.upstream_missing and os.path.exists(args.upstream_missing):
        try:
            upstream_df = pd.read_csv(args.upstream_missing)
            if not upstream_df.empty:
                upstream_rows = upstream_df.to_dict('records')
                print(f"Loaded {len(upstream_rows)} upstream missing entries")
        except Exception as e:
            print(f"Warning: Could not read upstream missing.csv: {e}")

    # The per-PDB loop in the tool script may already have written failure rows here.
    local_rows = []
    if os.path.exists(args.missing_csv):
        try:
            existing = pd.read_csv(args.missing_csv)
            if not existing.empty:
                local_rows = existing.to_dict('records')
        except Exception as e:
            print(f"Warning: Could not read existing missing.csv: {e}")

    all_missing = upstream_rows + local_rows + dropped
    if all_missing:
        missing_df = pd.DataFrame(all_missing).drop_duplicates(subset=['id'], keep='first')
    else:
        missing_df = pd.DataFrame(columns=['id', 'removed_by', 'kind', 'cause'])
    missing_df.to_csv(args.missing_csv, index=False)
    n_fail = sum(1 for r in dropped if r['kind'] == 'failure')
    print(f"Created missing.csv with {len(dropped) - n_fail} filtered, {n_fail} failed, "
          f"{len(upstream_rows)} upstream entries")
