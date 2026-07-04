#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Run LASErMPNN batch inference and collect designs into streams.

Runs ``python -m LASErMPNN.run_batch_inference`` over all input structures (fed
as a path list so each gets its own output subdir), then:
  * renames each ``<exec>/<stem>/design_{n}.pdb`` into the structures stream as
    ``<parent>_<n+1>.pdb`` and writes structures_map.csv (with a structures.id
    provenance column);
  * parses ``designs.fasta`` for the per-design sequence and log-prob score and
    writes the content-bearing sequences.csv;
  * writes a missing table for inputs/ranks that produced no design.

The repo has no setup.py, so the module is invoked from the clone's parent dir
with it on PYTHONPATH — matching the upstream Colab notebook.

Usage: see argparse below.
"""

import sys
import os
import re
import glob
import json
import shlex
import shutil
import argparse
import subprocess

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import load_datastream, iterate_files

BFACTOR_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pipe_lasermpnn_bfactor.py")


def _input_stem(pdb_path):
    """Reproduce LASErMPNN's subdir stem: filename with all suffixes stripped."""
    base = os.path.basename(pdb_path)
    # LASErMPNN: file.stem if '.' not in stem else stem.split('.')[0]
    stem = base.split(".")[0]
    return stem


def prepare_inputs(entries, prepared_dir, positions_json, fix_beta):
    """Return (list_txt_path, stem_to_id).

    Every input is staged into ``<prepared>/<design_id>.pdb`` so the LASErMPNN
    output subdir stem is always the unique design id (avoids stem collisions
    across different ids, whatever the source filenames were). When fix_beta is
    set the staged copy is B-factor-stamped; otherwise it is a plain copy.
    """
    os.makedirs(prepared_dir, exist_ok=True)
    input_paths = []
    stem_to_id = {}
    for design_id, pdb_path in entries:
        out_pdb = os.path.join(prepared_dir, f"{design_id}.pdb")
        if fix_beta:
            subprocess.run(
                [sys.executable, BFACTOR_SCRIPT, "stamp", positions_json,
                 design_id, pdb_path, out_pdb],
                check=True,
            )
        else:
            shutil.copyfile(pdb_path, out_pdb)
        stem_to_id[_input_stem(out_pdb)] = design_id
        input_paths.append(os.path.abspath(out_pdb))

    list_txt = os.path.join(prepared_dir, "inputs.txt")
    with open(list_txt, "w") as f:
        f.write("\n".join(input_paths) + "\n")
    return list_txt, stem_to_id


def run_inference(repo_parent, list_txt, out_dir, num_sequences, device, run_options):
    """Invoke run_batch_inference as a module from the clone's parent dir."""
    os.makedirs(out_dir, exist_ok=True)
    env = os.environ.copy()
    env["PYTHONPATH"] = repo_parent + os.pathsep + env.get("PYTHONPATH", "")
    # LASErMPNN imports logomaker/matplotlib at module load; force a headless
    # backend so Colab's inline backend (exported via MPLBACKEND) can't crash it.
    env["MPLBACKEND"] = "Agg"
    cmd = (
        [sys.executable, "-m", "LASErMPNN.run_batch_inference",
         list_txt, out_dir, str(num_sequences), "-d", device]
        + run_options
    )
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=repo_parent, env=env)


def parse_fasta_scores(fasta_path):
    """Parse designs.fasta -> {(pdb_code, design_idx): (sequence, score)}.

    Header: >{pdb_code}_design_{idx}_segment_{seg}_chain_{ch} score={score}
    Keeps the first chain record per design (single-chain designs are the norm).
    """
    result = {}
    if not os.path.exists(fasta_path):
        return result
    header_re = re.compile(r"^>(?P<code>.+)_design_(?P<idx>\d+)_segment_.*?(?:score=(?P<score>[-\d.eE]+))?\s*$")
    cur = None
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                m = header_re.match(line)
                if m:
                    code = m.group("code")
                    idx = int(m.group("idx"))
                    score = m.group("score")
                    cur = (code, idx)
                    if cur not in result:
                        result[cur] = {"sequence": "", "score": score if score is not None else ""}
                    else:
                        cur = None  # already captured first chain for this design
                else:
                    cur = None
            elif cur is not None and line:
                result[cur]["sequence"] += line
    return result


def collect(out_dir, structures_dir, structures_map, sequences_csv, missing_csv,
            entries, stem_to_id, num_sequences):
    """Rename designs into the structures stream and build the tables."""
    import pandas as pd

    os.makedirs(structures_dir, exist_ok=True)
    fasta_scores = parse_fasta_scores(os.path.join(out_dir, "designs.fasta"))

    struct_rows = []
    seq_rows = []
    missing_rows = []

    for design_id, _ in entries:
        # Find this input's LASErMPNN subdir by any stem that maps back to it.
        subdirs = [s for s, sid in stem_to_id.items() if sid == design_id]
        produced = {}
        for stem in subdirs:
            sub = os.path.join(out_dir, stem)
            for pdb in glob.glob(os.path.join(sub, "design_*.pdb")):
                m = re.search(r"design_(\d+)\.pdb$", os.path.basename(pdb))
                if m:
                    produced[int(m.group(1))] = pdb

        for n in range(num_sequences):
            out_id = f"{design_id}_{n + 1}"
            if n in produced:
                dest = os.path.join(structures_dir, f"{out_id}.pdb")
                shutil.copyfile(produced[n], dest)
                struct_rows.append({"id": out_id, "file": dest, "structures.id": design_id})

                # Sequence/score keyed by (input stem, design idx).
                seq, score = "", ""
                for stem in subdirs:
                    if (stem, n) in fasta_scores:
                        seq = fasta_scores[(stem, n)]["sequence"]
                        score = fasta_scores[(stem, n)]["score"]
                        break
                seq_rows.append({"id": out_id, "structures.id": design_id,
                                 "sequence": seq, "score": score})
            else:
                missing_rows.append({"id": out_id, "removed_by": "LASErMPNN",
                                     "kind": "failure", "cause": "design not produced"})

    os.makedirs(os.path.dirname(structures_map), exist_ok=True)
    pd.DataFrame(struct_rows, columns=["id", "file", "structures.id"]).to_csv(structures_map, index=False)

    os.makedirs(os.path.dirname(sequences_csv), exist_ok=True)
    pd.DataFrame(seq_rows, columns=["id", "structures.id", "sequence", "score"]).to_csv(sequences_csv, index=False)

    os.makedirs(os.path.dirname(missing_csv), exist_ok=True)
    pd.DataFrame(missing_rows, columns=["id", "removed_by", "kind", "cause"]).to_csv(missing_csv, index=False)

    print(f"Collected {len(struct_rows)} designs; {len(missing_rows)} missing.")
    if not struct_rows:
        print("ERROR: LASErMPNN produced no designs", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Run LASErMPNN and collect designs")
    parser.add_argument("--structures-json", required=True)
    parser.add_argument("--repo-parent", required=True)
    parser.add_argument("--exec-root", required=True)
    parser.add_argument("--structures-dir", required=True)
    parser.add_argument("--structures-map", required=True)
    parser.add_argument("--sequences-csv", required=True)
    parser.add_argument("--missing-csv", required=True)
    parser.add_argument("--num-sequences", type=int, required=True)
    parser.add_argument("--device", required=True)
    parser.add_argument("--run-options", required=True,
                        help="Shared run_batch_inference flags as one string")
    parser.add_argument("--positions-json", default=None)
    parser.add_argument("--fix-beta", action="store_true")
    args = parser.parse_args()

    ds = load_datastream(args.structures_json)
    entries = list(iterate_files(ds))
    if not entries:
        raise ValueError(f"No structures in DataStream: {args.structures_json}")

    prepared_dir = os.path.join(args.exec_root, "_inputs")
    list_txt, stem_to_id = prepare_inputs(entries, prepared_dir, args.positions_json, args.fix_beta)

    out_dir = os.path.join(args.exec_root, "designs")
    run_options = shlex.split(args.run_options)
    if args.fix_beta:
        run_options = run_options + ["--fix_beta"]
    run_inference(args.repo_parent, list_txt, out_dir, args.num_sequences,
                  args.device, run_options)

    collect(out_dir, args.structures_dir, args.structures_map, args.sequences_csv,
            args.missing_csv, entries, stem_to_id, args.num_sequences)


if __name__ == "__main__":
    main()
