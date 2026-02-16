# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.


import argparse

parser = argparse.ArgumentParser(description='Makes csv files from ProteinMPNN.fa output ')
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

# Parse the arguments
args = parser.parse_args()

import os
import json
import pandas as pd

# Load optional ID map
id_map = None
if args.id_map and os.path.exists(args.id_map):
    with open(args.id_map, 'r') as f:
        id_map = json.load(f)

fa_files = os.listdir(args.FA_FOLDER)
seen_sequences = {}  # sequence -> first_seen_id
duplicate_entries = []  # list of {id, removed_by, cause} dicts
for fa in fa_files:
    if fa.endswith(".fa"):
        data = []
        with open(os.path.join(args.FA_FOLDER,fa),"r") as file:
            lines = [line.strip() for line in file.readlines()]
            #Starts from 2 to skip the original protein
            for i in range(2,len(lines),2):
                seq_data = dict()
                seq_data["id"] = "" #put here just to preserve the order
                raw_sequence = lines[i+1]
                # Detect gap positions (X = unknown/missing residues)
                gap_indices = [pos + 1 for pos, aa in enumerate(raw_sequence) if aa == "X"]
                seq_data["gaps"] = "+".join(str(p) for p in gap_indices) if gap_indices else ""
                if gap_indices and args.fill_gaps:
                    raw_sequence = raw_sequence.replace("X", args.fill_gaps)
                    print(f"  Filled {len(gap_indices)} gap(s) at position(s) {seq_data['gaps']} with {args.fill_gaps}")
                seq_data["sequence"] = raw_sequence
                params = lines[i][1:].split(", ")
                for p in params:
                    if not '=' in p:
                        continue
                    key,value = p.split("=")
                    if key == 'id': key = 'sample' #ligandmpnn is structured differently from proteinmpnn
                    seq_data[key] = value
                pdb_base = fa[:-3]
                mapped_base = id_map.get(pdb_base, pdb_base) if id_map else pdb_base
                seq_data["id"] = mapped_base + "_" + seq_data["sample"]
                seq_data["structures.id"] = mapped_base
                if seq_data["sequence"] in seen_sequences:
                    if args.duplicates:
                        data.append(seq_data)
                        # Don't update seen_sequences - keep first
                    else:
                        kept_id = seen_sequences[seq_data["sequence"]]
                        print(f"Skipped duplicate: {seq_data['id']}")
                        if args.step_tool_name:
                            duplicate_entries.append({
                                'id': seq_data['id'],
                                'removed_by': args.step_tool_name,
                                'cause': f"Duplicate of {kept_id}"
                            })
                else:
                    seen_sequences[seq_data["sequence"]] = seq_data["id"]
                    data.append(seq_data)
        if len(data) == 0: continue
        #Update CSV
        old_data = []
        if os.path.exists(args.queries_csv_file):
            with open(args.queries_csv_file,"r") as old_queries:
                old_data = old_queries.readlines()
        with open(args.queries_csv_file,"w") as queries:
            keys = list(data[0].keys())
            if len(old_data) != 0:
                queries.writelines(old_data)
            #Header is already part of old_data
            else:
                for i,k in enumerate(keys):
                    if i>0: queries.write(',')
                    queries.write(k)
            for seq_data in data:
                queries.write("\n")
                for i,k in enumerate(keys):
                    if i>0: queries.write(',')
                    queries.write(seq_data[k])
        #Update FASTA
        old_data = []
        if os.path.exists(args.queries_fasta_file):
            with open(args.queries_fasta_file,"r") as old_queries:
                old_data = old_queries.readlines()
        with open(args.queries_fasta_file,"w") as queries:
            if len(old_data) != 0:
                queries.writelines(old_data)
                queries.write("\n")
            for s_d,seq_data in enumerate(data):
                if s_d>0: queries.write('\n')
                queries.write(">")
                queries.write(seq_data["id"])
                queries.write('\n')
                queries.write(seq_data["sequence"])

# Write missing.csv if requested
if args.missing_csv:
    # Load upstream missing entries if provided
    upstream_rows = []
    if args.upstream_missing and os.path.exists(args.upstream_missing):
        try:
            upstream_df = pd.read_csv(args.upstream_missing)
            if not upstream_df.empty:
                upstream_rows = upstream_df.to_dict('records')
                print(f"Loaded {len(upstream_rows)} upstream missing entries")
        except Exception as e:
            print(f"Warning: Could not read upstream missing.csv: {e}")

    all_missing = upstream_rows + duplicate_entries
    if all_missing:
        missing_df = pd.DataFrame(all_missing)
    else:
        missing_df = pd.DataFrame(columns=['id', 'removed_by', 'cause'])
    missing_df.to_csv(args.missing_csv, index=False)
    print(f"Created missing.csv with {len(duplicate_entries)} duplicates, {len(upstream_rows)} upstream entries")
