#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Collect RF3's predictions into the framework's streams and tables.

`rf3 fold` writes, per example `<name>`, a folder `<out>/<name>/` holding the top-ranked sample as
`<name>_model.cif`, `<name>_summary_confidences.json`, `<name>_confidences.json` and
`<name>_ranking_scores.csv`, plus one `seed-<S>_sample-<N>/` folder per sample with the same three files
under the name `<name>_seed-<S>_sample-<N>`. An example that early-stopped has no structure.
This turns that into one `<id>` per complex (or `<id>_1..N` when the user asked for every sample), a
structures map_table, a pae stream (the confidences JSON, which carries the full PAE matrix) and a
confidence table.

Three artifacts, kept separate on purpose: the confidence table has a row for every id that ENTERED, with
NaN where a number is unavailable; the stream maps have only the ids whose file was really written; the
local missing table says why an id has no structure. An example that failed is in all three.
"""

import argparse
import csv
import glob
import gzip
import json
import math
import os
import re
import shutil
import sys
from typing import Dict, List, Optional

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'biopipelines'))
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

CONFIDENCE_COLUMNS = ["id", "plddt", "ptm", "iptm", "pae", "pde", "has_clash",
                      "ranking_score", "early_stopped", "seed", "sample"]

# Column -> key of the summary JSON, as upstream's source spells them (overall_* are means).
SUMMARY_KEYS = {
    "plddt": "overall_plddt",
    "ptm": "ptm",
    "iptm": "iptm",
    "pae": "overall_pae",
    "pde": "overall_pde",
    "has_clash": "has_clash",
    "ranking_score": "ranking_score",
}

SAMPLE_DIR_RE = re.compile(r"^seed-(?P<seed>-?\d+)_sample-(?P<sample>\d+)$")
MODEL_EXTENSIONS = ("cif", "cif.gz")


def parse_bool(text) -> bool:
    return str(text).strip().lower() in ("1", "true", "yes", "on")


def example_names(queries_json: str) -> List[str]:
    with open(queries_json) as handle:
        return [str(example["name"]) for example in json.load(handle)]


def read_json(path: str) -> Dict:
    """A JSON object, or empty when the file is absent or unreadable.

    Absent is not an error: the structure is the result, and a confidence file that failed to write
    should leave NaN in a column rather than discard a prediction that exists.
    """
    if not os.path.exists(path):
        return {}
    try:
        with open(path) as handle:
            data = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return {}
    return data if isinstance(data, dict) else {}


def read_summary(path: str) -> Dict[str, Optional[float]]:
    data = read_json(path)
    return {column: data.get(key) for column, key in SUMMARY_KEYS.items()}


def read_csv_rows(path: str) -> List[Dict[str, str]]:
    if not os.path.exists(path):
        return []
    try:
        with open(path, newline="") as handle:
            return list(csv.DictReader(handle))
    except OSError:
        return []


def find_model(directory: str, stem: str) -> Optional[str]:
    """`<stem>_model.cif`, or its gzipped spelling when the run was given compress_outputs."""
    for extension in MODEL_EXTENSIONS:
        path = os.path.join(directory, f"{stem}_model.{extension}")
        if os.path.exists(path):
            return path
    return None


def best_sample(example_dir: str, name: str):
    """(seed, sample) of the sample upstream ranked first, from its ranking_scores.csv, else (None, None).

    The top-level files are a copy of that sample, so this only labels them. Ties go to the first row,
    as upstream's own max() does.
    """
    best = None
    for row in read_csv_rows(os.path.join(example_dir, f"{name}_ranking_scores.csv")):
        try:
            score = float(row["ranking_score"])
        except (KeyError, TypeError, ValueError):
            continue
        if math.isnan(score):
            continue
        if best is None or score > best[0]:
            best = (score, row.get("seed"), row.get("sample"))
    if best is None:
        return None, None
    return _as_int(best[1]), _as_int(best[2])


def _as_int(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def early_stop_info(example_dir: str, name: str) -> Optional[Dict]:
    """Evidence that the example was abandoned by the pLDDT early stop, or None.

    Upstream records it differently across releases: a `<name>_ranking_scores.csv` row with an
    `early_stopped` column in the newer ones, a `score.csv` in the older ones. Either carries `mean_plddt`.
    """
    for path in (os.path.join(example_dir, f"{name}_ranking_scores.csv"),
                 os.path.join(example_dir, "score.csv")):
        for row in read_csv_rows(path):
            if parse_bool(row.get("early_stopped", "")):
                return {"mean_plddt": _as_float(row.get("mean_plddt"))}
    return None


def _as_float(value) -> Optional[float]:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(number) else number


def sample_folders(example_dir: str):
    """(seed, sample, folder) for every sample folder, ordered by seed then sample."""
    found = []
    for path in glob.glob(os.path.join(example_dir, "seed-*_sample-*")):
        match = SAMPLE_DIR_RE.match(os.path.basename(path))
        if match and os.path.isdir(path):
            found.append((int(match.group("seed")), int(match.group("sample")), path))
    return sorted(found)


def collect_example(prediction_folder: str, name: str, top_only: bool, n_samples: int,
                    threshold: Optional[float] = None) -> Dict:
    """What one example produced.

    Returns {"outputs": [...], "dropped": [...]}. An output is {"id", "structure", "confidences",
    "summary", "seed", "sample"}; a dropped entry is {"id", "kind", "cause", "early_stopped"} and stands for
    an id the wrapper declared that has no structure.
    """
    example_dir = os.path.join(prediction_folder, name)
    outputs: List[Dict] = []

    if top_only:
        model = find_model(example_dir, name)
        if model:
            seed, sample = best_sample(example_dir, name)
            outputs.append({
                "id": name, "structure": model,
                "confidences": os.path.join(example_dir, f"{name}_confidences.json"),
                "summary": read_summary(os.path.join(example_dir, f"{name}_summary_confidences.json")),
                "seed": seed, "sample": sample,
            })
        wanted = [name]
    else:
        for index, (seed, sample, folder) in enumerate(sample_folders(example_dir), start=1):
            stem = f"{name}_seed-{seed}_sample-{sample}"
            model = find_model(folder, stem)
            if not model:
                continue
            outputs.append({
                "id": f"{name}_{index}", "structure": model,
                "confidences": os.path.join(folder, f"{stem}_confidences.json"),
                "summary": read_summary(os.path.join(folder, f"{stem}_summary_confidences.json")),
                "seed": seed, "sample": sample,
            })
        wanted = [f"{name}_{k}" for k in range(1, n_samples + 1)]

    produced = {output["id"] for output in outputs}
    dropped: List[Dict] = []
    stopped = None if outputs else early_stop_info(example_dir, name)
    for wanted_id in wanted:
        if wanted_id in produced:
            continue
        if stopped is not None:
            plddt = stopped.get("mean_plddt")
            detail = f"mean pLDDT {plddt:.3f} after the first recycle" if plddt is not None \
                else "mean pLDDT after the first recycle"
            limit = f" below early_stopping_plddt_threshold {threshold}" if threshold is not None else \
                " below early_stopping_plddt_threshold"
            dropped.append({"id": wanted_id, "kind": "filter", "early_stopped": True,
                            "cause": f"early stopped: {detail}{limit}"})
        elif not os.path.isdir(example_dir):
            dropped.append({"id": wanted_id, "kind": "failure", "early_stopped": False,
                            "cause": "rf3 fold wrote no output folder for this input"})
        else:
            dropped.append({"id": wanted_id, "kind": "failure", "early_stopped": False,
                            "cause": "rf3 fold wrote no structure for this input"})
    return {"outputs": outputs, "dropped": dropped}


def skipped_example(name: str, cause: str, top_only: bool, n_samples: int) -> Dict:
    """What an example the config step could not write produced: nothing, under every id it declared."""
    wanted = [name] if top_only else [f"{name}_{k}" for k in range(1, n_samples + 1)]
    return {"outputs": [],
            "dropped": [{"id": wanted_id, "kind": "failure", "early_stopped": False, "cause": cause}
                        for wanted_id in wanted]}


def stage_structure(source: str, destination: str) -> None:
    """Copy the model to `<id>.cif`, decompressing when upstream gzipped it."""
    if source.endswith(".gz"):
        with gzip.open(source, "rb") as handle, open(destination, "wb") as target:
            shutil.copyfileobj(handle, target)
    else:
        shutil.copyfile(source, destination)


def confidence_row(out_id: str, summary: Dict, seed, sample, early_stopped: bool) -> Dict:
    row = {"id": out_id, "early_stopped": early_stopped, "seed": seed, "sample": sample}
    row.update({column: value for column, value in summary.items() if value is not None})
    return row


def write_rows(path: str, fieldnames: List[str], rows: List[Dict]) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in fieldnames})


def main():
    from biopipelines.biopipelines_io import read_skipped_queries, step_id_from_table_path

    parser = argparse.ArgumentParser(description="Collect RF3 predictions")
    parser.add_argument("--prediction-folder", required=True)
    parser.add_argument("--queries-json", required=True)
    parser.add_argument("--skipped-csv", default=None)
    parser.add_argument("--structures-dir", required=True)
    parser.add_argument("--structures-map-csv", required=True)
    parser.add_argument("--confidence-csv", required=True)
    parser.add_argument("--pae-dir", required=True)
    parser.add_argument("--pae-map-csv", required=True)
    parser.add_argument("--local-missing-csv", required=True)
    parser.add_argument("--top-only", default="True")
    parser.add_argument("--diffusion-batch-size", type=int, required=True)
    parser.add_argument("--early-stopping-plddt-threshold", type=float, default=None)
    args = parser.parse_args()

    names = example_names(args.queries_json)
    skipped = read_skipped_queries(args.skipped_csv)
    if not names and not skipped:
        print("ERROR: no examples to collect — the inputs JSON is empty")
        sys.exit(1)

    step_id = step_id_from_table_path(args.local_missing_csv)
    os.makedirs(args.structures_dir, exist_ok=True)
    os.makedirs(args.pae_dir, exist_ok=True)

    map_rows, pae_rows, confidence_rows, missing_rows = [], [], [], []
    failed_examples = 0

    top_only = parse_bool(args.top_only)
    results = [collect_example(args.prediction_folder, name, top_only,
                               args.diffusion_batch_size, args.early_stopping_plddt_threshold)
               for name in names]
    results += [skipped_example(name, cause, top_only, args.diffusion_batch_size)
                for name, cause in skipped.items()]

    for result in results:
        for output in result["outputs"]:
            destination = os.path.join(args.structures_dir, f"{output['id']}.cif")
            stage_structure(output["structure"], destination)
            map_rows.append({"id": output["id"], "file": destination})
            if os.path.exists(output["confidences"]):
                pae_destination = os.path.join(args.pae_dir, f"{output['id']}.json")
                shutil.copyfile(output["confidences"], pae_destination)
                pae_rows.append({"id": output["id"], "file": pae_destination})
            else:
                print(f"WARNING: no {os.path.basename(output['confidences'])} for {output['id']}, "
                      f"left out of the pae stream")
            confidence_rows.append(confidence_row(
                output["id"], output["summary"], output["seed"], output["sample"], False))
        for dropped in result["dropped"]:
            confidence_rows.append(confidence_row(dropped["id"], {}, None, None, dropped["early_stopped"]))
            missing_rows.append({"id": dropped["id"], "removed_by": step_id,
                                 "kind": dropped["kind"], "cause": dropped["cause"]})
        if any(dropped["kind"] == "failure" for dropped in result["dropped"]):
            failed_examples += 1

    write_rows(args.structures_map_csv, ["id", "file"], map_rows)
    write_rows(args.pae_map_csv, ["id", "file"], pae_rows)
    write_rows(args.confidence_csv, CONFIDENCE_COLUMNS, confidence_rows)
    write_rows(args.local_missing_csv, ["id", "removed_by", "kind", "cause"], missing_rows)

    print(f"Collected {len(map_rows)} structure(s) for {len(results)} example(s); "
          f"{len(missing_rows)} id(s) without one")

    # Every example failed means the step did not work; early-stopped examples were filtered, not failed.
    if failed_examples == len(results):
        print("ERROR: no RF3 prediction produced a structure")
        sys.exit(1)


if __name__ == "__main__":
    main()
