# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
A pipe script contains python code to execute at runtime, that is, it assumes that upstream tools have run and their outputs are present. They are invoked within the generate_script() of a tool, and can for example preprocess the inputs, run the model, standardize the outputs, and so on. A tool can call more than one pipe script if needed. It executes in the declared env or container (it has to follow activate_environment() in the generated script), and not necessarily inside biopipelines, so it imports biopipelines only for the shared IO helpers and stays dependency-light otherwise.

Contract with the wrapper (_TEMPLATE.py):
- Every value the wrapper sends must be consumed: a --flag needs a matching add_argument; a config.yaml key needs to be read. This file mirrors the config.yaml style the wrapper uses (params in the yaml, paths as flags).
- It reloads the input stream from the JSON the wrapper serialized (--structures-json).
- It writes each declared output to the path the wrapper passed, plus that stream's id>entity map CSV.
- Computation, failure-reporting and selection are separate operations, one artifact each: a result table carries one row per id that entered (unavailable values are NaN, the id is never dropped); a stream map carries only ids whose file was actually produced; the missing table records why an id has no entity. A failed id therefore appears in all three - NaN in the table, absent from the map, a kind="failure" row in missing. Do not filter a result table on a threshold either; emit it complete and let the user apply Panda.filter downstream.
- Missing entities can be managed in one of two ways:
    A. (matching the template) it writes only its own per-id failures to --local-missing-csv. It does not read the upstream missing table; the wrapper emits a separate propagation step that merges upstream + this local file into tables/missing.csv.
    B. accepts --upstream-missing, reads it (biopipelines_io.read_upstream_missing), and writes the merged tables/missing.csv here.
To adapt: copy to pipe_scripts/pipe_<yourtool>.py, edit.
"""

import argparse
import os
import sys
from typing import Dict

import pandas as pd
import yaml

# Add repo root to sys.path so the biopipelines package is importable when this runs out-of-env
# (in the tool's own conda env / container, where biopipelines may not be installed).
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

""" Dependencies
Shared IO helpers, imported after the sys.path shim above.
- load_datastream: rebuild the DataStream from the wrapper's JSON.
- iterate_files(ds)  -> (id, resolved_file_path)   for FILE-BASED streams (one file per id). Used below.
- iterate_values(ds, columns=[...]) -> (id, {col: val})  for VALUE-BASED streams (data inline in map_table, files=[]), e.g. a sequences stream's `sequence` or a compounds stream's `smiles`.
- step_id_from_table_path: derive this step's "<order>_<Tool>" id, stamped into missing.removed_by.
- container_argv_prefix: turn --container-prefix into an argv list to prepend to a subprocess call.     For example: apptainer ... run ...
"""
from biopipelines.biopipelines_io import (  # noqa: E402
    load_datastream,
    iterate_files,
    iterate_values,
    step_id_from_table_path,
    container_argv_prefix,
)

def main():
    p = argparse.ArgumentParser()
    # Params come from the yaml; only paths/plumbing are flags. (Add one add_argument per wrapper --flag.)
    p.add_argument("--config-yaml", required=True)
    p.add_argument("--structures-json", required=True)
    p.add_argument("--annotated-dir", required=True)
    p.add_argument("--annotated-map-csv", required=True)
    p.add_argument("--counts-csv", required=True)
    p.add_argument("--local-missing-csv", required=True)
    p.add_argument("--container-prefix", default="")
    args = p.parse_args()

    with open(args.config_yaml) as f:
        config = yaml.safe_load(f) or {}
    mode = config.get("mode", "default")  # read every param the wrapper wrote

    os.makedirs(args.annotated_dir, exist_ok=True)

    ds = load_datastream(args.structures_json)
    step_id = step_id_from_table_path(args.local_missing_csv)

    """ If this tool launched an external binary, build its argv prefix once:
        - cprefix = container_argv_prefix(args.container_prefix)
        - subprocess.run(cprefix + ["mkdssp", pdb_path, out_path], check=True)
        Returns [] in env mode, so the same call works either way."""

    """Computation, failure-reporting and selection are separate operations, and each output answers exactly one of them.
    - RESULT TABLE (counts.csv): one row per id that entered. An unavailable measurement is NaN; the id is never dropped. A table is a matrix of results, not a filter.
    - STREAM MAP (annotated_map.csv): only ids whose FILE was actually produced. This one IS selective.
    - MISSING (local_missing.csv): why an id has no stream entity - kind="failure" (raised) or kind="filter" (deliberately dropped).
    """
    annotated_rows, counts_rows, local_missing_rows = [], [], []
    attempted = 0   # ids this tool tried to process
    failures = 0    # of those, the ones that raised (a filter is NOT a failure)

    for sid, pdb_path in iterate_files(ds):
        attempted += 1
        try:
            n_residues = _count_residues(pdb_path)

            """ Filtering IDs
            Deliberate drop: if THIS tool chooses to exclude an id from its STREAM, record it as kind="filter" (excused; reason goes in cause) and skip writing the file.
            if n_residues < min_residues:
                counts_rows.append({"id": sid, "n_residues": n_residues})
                local_missing_rows.append({"id": sid,
                                           "removed_by": step_id,
                                           "kind": "filter",
                                           "cause": "too few residues"})
                continue"""

            out_path = os.path.join(args.annotated_dir, f"{sid}.txt")
            with open(out_path, "w") as f:
                f.write(f"{sid}\t{mode}\t{n_residues}\n")

            annotated_rows.append({"id": sid, "file": out_path})   
            counts_rows.append({"id": sid, "n_residues": n_residues})
        except Exception as e:
            print(f"WARNING: {sid} failed: {e}", file=sys.stderr)
            failures += 1
            counts_rows.append({"id": sid, "n_residues": None})
            local_missing_rows.append({
                "id": sid,
                "removed_by": step_id,
                "kind": "failure",
                "cause": str(e)[:200],
            })

    # Write every declared output.
    for path in (args.counts_csv, args.annotated_map_csv, args.local_missing_csv):
        os.makedirs(os.path.dirname(path), exist_ok=True)
    pd.DataFrame(counts_rows, columns=["id", "n_residues"]).to_csv(args.counts_csv, index=False)
    pd.DataFrame(annotated_rows, columns=["id", "file"]).to_csv(args.annotated_map_csv, index=False)
    # Local failures/drops only. The wrapper's propagation step merges these with upstream into tables/missing.csv.
    pd.DataFrame(local_missing_rows, columns=["id", "removed_by", "kind", "cause"]).to_csv(args.local_missing_csv, index=False)

    print(f"Counts: {args.counts_csv} ({len(counts_rows)} rows, {failures} NaN)")
    print(f"Annotated: {args.annotated_map_csv} ({len(annotated_rows)} rows)")
    print(f"Local failures: {failures}")

    # Two step-level errors, both non-zero:
    # - nothing was attempted: iterate_files skipped every declared id, so the upstream stream is empty or none of its files exist. Exiting 0 here would rubber-stamp a broken pipeline and hide the upstream tool that actually failed.
    # - every id that was attempted raised.
    # A run that attempted ids and deliberately FILTERED all of them still exits 0 - the tool worked, it just selected nothing. Never fail merely because the result table carries NaN.
    if not attempted:
        print("ERROR: no input entities to process (upstream stream is empty or its files are absent)",
              file=sys.stderr)
        sys.exit(1)
    if failures == attempted:
        print(f"ERROR: all {attempted} attempted id(s) failed", file=sys.stderr)
        sys.exit(1)


def _count_residues(pdb_path: str) -> int:
    """Trivial stand-in for the tool's real per-id work: count ATOM/HETATM residues."""
    seen = set()
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                seen.add((line[21], line[22:27]))  # (chain, resSeq+iCode)
    return len(seen)


"""SHARED-FILE reader/writer (delete if your tool has no shared-file stream).
A shared-file stream is one artifact whose records correspond to all ids. iterate_files yields (id, same_path) for every id, so parse it once here and index by id. 
When a Panda/Pool filters a shared-file stream it re-emits it with only the kept ids via biopipelines.stream_slicers.get_slicer(<format>), and get_merger when gathering. For this, the stream's declared `format` must have slicer registered there."""
def _read_shared_fasta(fasta_path: str) -> Dict[str, str]:
    """Parse a multi-record FASTA into {id: sequence}, reading the single file once."""
    records: Dict[str, str] = {}
    current_id = None
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                current_id = line[1:].split()[0]
                records[current_id] = ""
            elif current_id is not None:
                records[current_id] += line.strip()
    return records


def _write_shared_fasta(fasta_path: str, rows) -> None:
    """Write one multi-record FASTA covering all ids (rows: iterable of {'id','sequence'})."""
    os.makedirs(os.path.dirname(fasta_path), exist_ok=True)
    with open(fasta_path, "w") as f:
        for row in rows:
            f.write(f">{row['id']}\n{row['sequence']}\n")


if __name__ == "__main__":
    main()
