#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
PLACER post-processing.

The driver runs PLACER.protocol.dump_output per input unit, which writes
upstream's canonical artefacts into the runs folder:

    {runs}/{unit_id}_model.pdb   — all N predicted models concatenated, each
                                   framed as `MODEL <n>` ... `ENDMDL`
    {runs}/{unit_id}.csv         — one row per model with score columns
                                   (model_idx, plddt, plddt_pde, fape, and in
                                   ligand mode also prmsd, rmsd, kabsch, lddt)

`unit_id` is `<structure>+<ligand>` in ligand mode and `<structure>` in
sidechain/apo mode (no ligand axis). The driver's `units.json` gives the parents
of every unit, which are read from it by exact id and never recovered from the id.

This script flattens those into BioPipelines streams:

  - splits the multi-model PDB into one file per model, named
    `<unit_id>_<n>.pdb` where <n> is the MODEL number (matches the CSV's
    model_idx, so file<->score alignment survives reranking);
  - writes the structures map_table (id | file | provenance columns);
  - emits the scores table (one row per sample, with provenance columns);
  - (ligand mode) passes the input ligand chemistry through unchanged into the
    compounds map_table — PLACER refines the pose, not the identity;
  - records every run of `units.json` that gave no model in the missing table, with
    the error `failures.json` holds for it, or a note that the driver never reached it.
"""

import argparse
import csv
import json
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import load_datastream, iterate_values, step_id_from_table_path  # noqa: E402
from biopipelines.combinatorics import generate_multiplied_ids  # noqa: E402

# Score columns we surface, per mode, in this order. dump_output may emit a
# subset; any missing column is left blank. Apo/sidechain runs do not predict a
# ligand, so the ligand-accuracy terms (prmsd/rmsd/kabsch) are absent.
SCORE_COLUMNS_LIGAND = ["prmsd", "plddt", "plddt_pde", "fape", "rmsd", "kabsch"]
SCORE_COLUMNS_SIDECHAIN = ["plddt", "plddt_pde", "fape"]

MODEL_RE = re.compile(r"^MODEL\s+(\d+)\s*$")


def split_multimodel_pdb(path):
    """Yield (model_number:int, pdb_text:str) for each MODEL...ENDMDL block.

    The per-model text keeps everything between the MODEL and ENDMDL lines
    (exclusive of those framing lines) and is wrapped with a fresh END so each
    file is a standalone single-model PDB.
    """
    with open(path) as f:
        lines = f.readlines()

    current_n = None
    buf = []
    for line in lines:
        m = MODEL_RE.match(line.rstrip("\n"))
        if m:
            current_n = int(m.group(1))
            buf = []
            continue
        if line.startswith("ENDMDL"):
            if current_n is not None:
                yield current_n, "".join(buf) + "END\n"
            current_n = None
            buf = []
            continue
        if current_n is not None:
            buf.append(line)


def load_units(runs_folder):
    """{unit id: {"structures.id": ..., "compounds.id": ...}} as the driver wrote it."""
    with open(os.path.join(runs_folder, "units.json")) as handle:
        return json.load(handle)


def load_failures(runs_folder):
    """{run id: error} the driver recorded for the runs that raised; empty when it wrote none."""
    path = os.path.join(runs_folder, "failures.json")
    if not os.path.isfile(path):
        return {}
    with open(path) as handle:
        return json.load(handle)


def load_scores(csv_path, score_columns):
    """Return {model_idx:int -> {score_col: value}} from a dump_output CSV.

    Keyed on the 1-based model number so it aligns with the split PDB files.
    """
    by_idx = {}
    df = pd.read_csv(csv_path)
    if "model_idx" not in df.columns:
        raise ValueError(f"{csv_path} has no 'model_idx' column; columns={list(df.columns)}")
    for _, row in df.iterrows():
        idx = int(row["model_idx"])
        by_idx[idx] = {col: row[col] for col in score_columns if col in df.columns}
    return by_idx


def write_compounds_passthrough(ligand_json, compounds_map):
    """Copy the input ligand chemistry into our compounds map_table verbatim
    (ids, code, smiles, ... unchanged). PLACER does not rename the ligand."""
    ds = load_datastream(ligand_json)
    rows = []
    fieldnames = None
    for cid, values in iterate_values(ds):
        row = {"id": cid}
        row.update(values)
        rows.append(row)
        if fieldnames is None:
            fieldnames = ["id"] + [k for k in values.keys() if k != "id"]
    if fieldnames is None:
        fieldnames = ["id"]
    with open(compounds_map, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=["ligand", "sidechain"])
    ap.add_argument("--runs-folder", required=True)
    ap.add_argument("--structures-folder", required=True)
    ap.add_argument("--structures-map", required=True)
    ap.add_argument("--scores-csv", required=True)
    ap.add_argument("--missing-csv", required=True)
    # Ligand-mode only: the compounds passthrough.
    ap.add_argument("--ligand-json", default=None)
    ap.add_argument("--compounds-map", default=None)
    args = ap.parse_args()

    if args.mode == "ligand" and not (args.ligand_json and args.compounds_map):
        print("ERROR: --ligand-json and --compounds-map are required in ligand mode", file=sys.stderr)
        sys.exit(1)

    is_ligand = args.mode == "ligand"
    score_columns = SCORE_COLUMNS_LIGAND if is_ligand else SCORE_COLUMNS_SIDECHAIN
    # Provenance columns on the maps/scores: ligand mode tracks both axes.
    prov_cols = ["structures.id", "compounds.id"] if is_ligand else ["structures.id"]

    dirs = [
        args.structures_folder,
        os.path.dirname(args.structures_map),
        os.path.dirname(args.scores_csv),
        os.path.dirname(args.missing_csv),
    ]
    if is_ligand:
        dirs.append(os.path.dirname(args.compounds_map))
    for d in dirs:
        os.makedirs(d, exist_ok=True)

    if not os.path.isdir(args.runs_folder):
        print(f"ERROR: runs folder missing: {args.runs_folder}", file=sys.stderr)
        sys.exit(1)

    # Discover units from the {unit_id}_model.pdb files the driver wrote.
    model_pdbs = sorted(
        f for f in os.listdir(args.runs_folder) if f.endswith("_model.pdb")
    )

    map_rows = []
    score_rows = []
    missing_rows = []
    step_id = step_id_from_table_path(args.missing_csv)

    units = load_units(args.runs_folder)
    unknown = [fname[: -len("_model.pdb")] for fname in model_pdbs if fname[: -len("_model.pdb")] not in units]
    if unknown:
        print(f"ERROR: run(s) that units.json does not list: {unknown}", file=sys.stderr)
        sys.exit(1)

    for fname in model_pdbs:
        unit_id = fname[: -len("_model.pdb")]
        model_path = os.path.join(args.runs_folder, fname)
        scores_path = os.path.join(args.runs_folder, f"{unit_id}.csv")

        scores_by_idx = load_scores(scores_path, score_columns) if os.path.isfile(scores_path) else {}

        prov = {col: units[unit_id][col] for col in prov_cols}

        n_models = 0
        for model_n, pdb_text in split_multimodel_pdb(model_path):
            out_id = generate_multiplied_ids([unit_id], [str(model_n)])[0][0]
            dst_pdb = os.path.join(args.structures_folder, f"{out_id}.pdb")
            with open(dst_pdb, "w") as f:
                f.write(pdb_text)

            map_rows.append({"id": out_id, "file": dst_pdb, **prov})

            srow = {"id": out_id, **prov, "sample": str(model_n)}
            model_scores = scores_by_idx.get(model_n, {})
            for col in score_columns:
                srow[col] = model_scores.get(col, "")
            score_rows.append(srow)
            n_models += 1

        if n_models == 0:
            missing_rows.append({
                "id": unit_id,
                "removed_by": step_id,
                "kind": "failure",
                "cause": "no models in _model.pdb",
            })

    failures = load_failures(args.runs_folder)
    modeled = {fname[: -len("_model.pdb")] for fname in model_pdbs}
    for unit_id in units:
        if unit_id not in modeled:
            cause = failures.get(unit_id, "PLACER wrote no model for this run")
            missing_rows.append({
                "id": unit_id,
                "removed_by": step_id,
                "kind": "failure",
                "cause": " ".join(str(cause).split()),
            })

    with open(args.structures_map, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["id", "file"] + prov_cols)
        w.writeheader()
        w.writerows(map_rows)

    with open(args.scores_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["id"] + prov_cols + ["sample"] + score_columns)
        w.writeheader()
        w.writerows(score_rows)

    if is_ligand:
        write_compounds_passthrough(args.ligand_json, args.compounds_map)

    with open(args.missing_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["id", "removed_by", "kind", "cause"])
        w.writeheader()
        w.writerows(missing_rows)

    with_models = len(modeled) - sum(1 for row in missing_rows if row["id"] in modeled)
    print(
        f"PLACER post-process ({args.mode}): {len(map_rows)} model(s) from "
        f"{with_models}/{len(units)} run(s); {len(missing_rows)} run(s) without a model"
    )

    if not map_rows:
        print("ERROR: PLACER produced no usable models", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
