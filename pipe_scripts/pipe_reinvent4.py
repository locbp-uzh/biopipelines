#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Runtime helper for REINVENT4.

Runs `reinvent` in its sampling mode once per input compound (once in total for the de novo generator), then turns each CSV it wrote into the compounds stream, the samples table and the inputs table.
One process per input is what lets every sampled molecule keep the id of the compound it came from: REINVENT's output carries the input only in the form it gave the model (randomized or standardized), never as the line that was read, and it aborts the whole run on one input it cannot read.
"""

import argparse
import json
import os
import subprocess
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from biopipelines.biopipelines_io import (  # noqa: E402
    load_datastream,
    iterate_values,
    step_id_from_table_path,
    container_argv_prefix,
)
from biopipelines.reinvent4 import (  # noqa: E402
    CONDITIONAL_GENERATORS,
    PROVENANCE_COLUMN,
    INPUT_COLUMNS,
    UPSTREAM_COLUMNS,
    GENERATED_CODE,
    REINVENT4,
)

MISSING_COLUMNS = ["id", "removed_by", "kind", "cause"]
NUMERIC_COLUMNS = ("nll", "tanimoto")
DEVICE_PROBE = "import torch; print(torch.cuda.is_available())"


def toml_value(value):
    """A Python scalar as a TOML literal; strings take JSON escaping, a subset of TOML's basic strings."""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return repr(value)
    return json.dumps(str(value))


def build_sampling_toml(cfg, smiles_file, output_file, device):
    """The REINVENT sampling configuration for one run."""
    parameters = {"model_file": cfg["model_file"]}
    if smiles_file is not None:
        parameters["smiles_file"] = smiles_file
    parameters["output_file"] = output_file
    parameters["num_smiles"] = cfg["num_compounds"]
    # Without it REINVENT also writes the invalid molecules.
    parameters["unique_molecules"] = True
    if cfg.get("sample_strategy") is not None:
        parameters["sample_strategy"] = cfg["sample_strategy"]
    if cfg.get("temperature") is not None:
        parameters["temperature"] = float(cfg["temperature"])

    lines = ['run_type = "sampling"', f"device = {toml_value(device)}", "", "[parameters]"]
    lines += [f"{key} = {toml_value(value)}" for key, value in parameters.items()]
    return "\n".join(lines) + "\n"


def read_samples(path, generator):
    """The rows of REINVENT's sampling CSV under this tool's column names."""
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    columns = UPSTREAM_COLUMNS[generator]
    absent = [column for column in columns if column not in frame.columns]
    if absent:
        raise ValueError(f"REINVENT's output lacks the column(s) {absent}; it has {list(frame.columns)}")
    rows = []
    for _, record in frame.iterrows():
        row = {ours: record[theirs] for theirs, ours in columns.items()}
        for column in NUMERIC_COLUMNS:
            if column in row:
                row[column] = float(row[column])
        rows.append(row)
    return rows


def input_text(smiles):
    """The input SMILES as one line of REINVENT's input file, or None if it cannot be one."""
    if smiles is None or (isinstance(smiles, float) and pd.isna(smiles)):
        return None
    text = str(smiles).strip()
    return text or None


def cuda_available(container_prefix):
    """Ask the interpreter REINVENT will run under, so a CPU build or a CPU node selects the CPU."""
    cmd = container_argv_prefix(container_prefix) + ["python", "-c", DEVICE_PROBE]
    answer = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout.split()
    if not answer or answer[-1] not in ("True", "False"):
        raise RuntimeError(f"could not tell whether CUDA is available, the probe printed {answer!r}")
    return answer[-1] == "True"


def last_line(path):
    """The last non-empty line of a log, for a failure's cause."""
    try:
        with open(path) as handle:
            lines = [line.strip() for line in handle if line.strip()]
    except OSError:
        return "no log"
    return lines[-1][:200] if lines else "empty log"


def run_reinvent(cfg, label, smiles, device, container_prefix):
    """Run one sampling and return its rows, or raise with the reason it produced none."""
    work = cfg["sampling_dir"]
    stem = os.path.join(work, label)
    smiles_file = f"{stem}.smi" if smiles is not None else None
    output_file = f"{stem}.csv"
    toml_file = f"{stem}.toml"
    log_file = f"{stem}.log"

    if smiles_file is not None:
        with open(smiles_file, "w") as handle:
            handle.write(smiles + "\n")
    with open(toml_file, "w") as handle:
        handle.write(build_sampling_toml(cfg, smiles_file, output_file, device))

    cmd = container_argv_prefix(container_prefix) + ["reinvent"]
    if cfg.get("seed") is not None:
        cmd += ["--seed", str(cfg["seed"])]
    cmd.append(toml_file)
    print(f"Running: {' '.join(cmd)}", flush=True)
    with open(log_file, "w") as log:
        result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, cwd=work)
    if result.returncode != 0:
        raise RuntimeError(f"reinvent exited with status {result.returncode}: {last_line(log_file)}")
    return read_samples(output_file, cfg["generator"])


def main():
    p = argparse.ArgumentParser(description="Sample molecules with REINVENT4")
    p.add_argument("--config", required=True)
    p.add_argument("--container-prefix", default="")
    args = p.parse_args()

    with open(args.config) as f:
        cfg = json.load(f)

    generator = cfg["generator"]
    conditional = generator in CONDITIONAL_GENERATORS
    step_id = step_id_from_table_path(cfg["local_missing_csv"])
    os.makedirs(cfg["sampling_dir"], exist_ok=True)

    # An input with no SMILES is a deliberate drop and one REINVENT cannot run is a failure; both keep a row in inputs.
    jobs, input_rows, local_missing_rows = [], [], []
    if conditional:
        compounds_ds = load_datastream(cfg["compounds_json"])
        for input_id, values in iterate_values(compounds_ds, columns=["smiles"]):
            smiles = input_text(values.get("smiles"))
            if smiles is None:
                print(f"WARNING: skipping {input_id} (no smiles)", file=sys.stderr)
                input_rows.append({"id": input_id, "smiles": "", "num_compounds": None})
                local_missing_rows.append({
                    "id": input_id, "removed_by": step_id, "kind": "filter",
                    "cause": "no smiles in compounds stream",
                })
            elif any(ch.isspace() for ch in smiles):
                print(f"WARNING: skipping {input_id} (whitespace in smiles)", file=sys.stderr)
                input_rows.append({"id": input_id, "smiles": smiles, "num_compounds": None})
                local_missing_rows.append({
                    "id": input_id, "removed_by": step_id, "kind": "failure",
                    "cause": "smiles contains whitespace, which REINVENT's input file cannot carry",
                })
            else:
                jobs.append((input_id, smiles))
    else:
        jobs.append((None, None))

    device = "cuda:0" if jobs and cuda_available(args.container_prefix) else "cpu"
    print(f"Device: {device}", flush=True)

    compound_rows, sample_rows = [], []
    failures = 0
    for input_id, smiles in jobs:
        label = input_id if conditional else "denovo"
        id_prefix = input_id if conditional else cfg["id_prefix"]
        try:
            rows = run_reinvent(cfg, label, smiles, device, args.container_prefix)
        except Exception as e:
            print(f"WARNING: {label} failed: {e}", file=sys.stderr)
            failures += 1
            if conditional:
                input_rows.append({"id": input_id, "smiles": smiles, "num_compounds": None})
                local_missing_rows.append({
                    "id": input_id, "removed_by": step_id, "kind": "failure", "cause": str(e)[:200],
                })
            continue

        for index, row in enumerate(rows, start=1):
            compound_id = f"{id_prefix}_{index}"
            compound = {"id": compound_id, "format": "smiles", "code": GENERATED_CODE,
                        "smiles": row["smiles"], "ccd": ""}
            sample = {"id": compound_id, **row}
            if conditional:
                compound[PROVENANCE_COLUMN] = input_id
                sample[PROVENANCE_COLUMN] = input_id
            compound_rows.append(compound)
            sample_rows.append(sample)

        if conditional:
            input_rows.append({"id": input_id, "smiles": smiles, "num_compounds": len(rows)})
            if not rows:
                local_missing_rows.append({
                    "id": input_id, "removed_by": step_id, "kind": "filter",
                    "cause": "REINVENT returned no valid molecule for this input",
                })

    outputs = [cfg["compounds_csv"], cfg["samples_csv"], cfg["local_missing_csv"]]
    if conditional:
        outputs.append(cfg["inputs_csv"])
    for path in outputs:
        os.makedirs(os.path.dirname(path), exist_ok=True)

    pd.DataFrame(compound_rows, columns=REINVENT4.compound_columns(generator)).to_csv(cfg["compounds_csv"], index=False)
    pd.DataFrame(sample_rows, columns=REINVENT4.sample_columns(generator)).to_csv(cfg["samples_csv"], index=False)
    if conditional:
        order = {input_id: position for position, input_id in enumerate(compounds_ds.ids_expanded)}
        input_rows.sort(key=lambda row: order.get(row["id"], len(order)))
        inputs = pd.DataFrame(input_rows, columns=INPUT_COLUMNS)
        inputs["num_compounds"] = inputs["num_compounds"].astype("Int64")
        inputs.to_csv(cfg["inputs_csv"], index=False)
    pd.DataFrame(local_missing_rows, columns=MISSING_COLUMNS).to_csv(cfg["local_missing_csv"], index=False)

    print(f"Compounds: {cfg['compounds_csv']} ({len(compound_rows)} rows)")
    print(f"Samples: {cfg['samples_csv']} ({len(sample_rows)} rows)")
    print(f"Failed runs: {failures} of {len(jobs)}")

    # Nothing attempted or every attempt failed is a broken step; inputs that yielded nothing are not.
    if not jobs:
        print("ERROR: no input compound with a usable SMILES (upstream stream is empty or carries none)",
              file=sys.stderr)
        sys.exit(1)
    if failures == len(jobs):
        print(f"ERROR: all {len(jobs)} REINVENT run(s) failed", file=sys.stderr)
        sys.exit(1)
    if not conditional and not compound_rows:
        print("ERROR: REINVENT returned no valid molecule", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
