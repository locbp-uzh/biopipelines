#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""BioEmu2 runner.

Each combination of the sequences and ligands axes is one run. For every run this writes a system YAML, calls `bioemu2 sample-equilibrium`, splits the sampled trajectory into one PDB per conformer and stages the compact `.xtc` and its topology.

Three artifacts, kept separate on purpose.
The summary table has a row for every run that ENTERED, with empty cells where a run failed.
The stream maps have only what was really written.
The local missing table says why a run has no entity; a failed run is in all three.
The wrapper's propagation step merges that file with the upstream manifests into tables/missing.csv.
"""

import argparse
import collections
import copy
import hashlib
import io
import itertools
import json
import os
import shutil
import string
import subprocess
import sys
from typing import Dict, List, Optional, Tuple

_biopipelines_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'biopipelines')
sys.path.insert(0, _biopipelines_dir)
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pandas as pd
import yaml

from combinatorics import CombinatoricsConfig, predict_single_output_id
from ligand_utils import auth_ligand_field, ligand_chemistry
from pipe_axis_data import load_axis_data
from biopipelines.biopipelines_io import container_argv_prefix, step_id_from_table_path  # noqa: E402
from biopipelines.id_map_utils import MsaLookup  # noqa: E402


# Entity order is chain order: protein chains first, then ligands.
AXIS_ORDER = ("sequences", "ligands")

SUMMARY_COLUMNS = ["id", "sequences.id", "ligands.id", "n_samples", "n_residues"]
PHYSICALITY_COLUMNS = ["id", "run.id", "frame", "status", "reason"]
MISSING_COLUMNS = ["id", "removed_by", "kind", "cause"]

CAUSE_LIMIT = 300
# Lines of upstream's output kept to read the cause from; a traceback ends well inside this.
OUTPUT_TAIL = 200
# Count columns of the summary table: integers, blank where a run failed.
INTEGER_COLUMNS = ("n_samples", "n_residues")

# The name a bond names a residue by; X is upstream's letter for an unidentified residue.
THREE_LETTER = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY", "H": "HIS",
    "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN",
    "R": "ARG", "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR", "X": "UNK",
}


class RunError(ValueError):
    """A run that cannot be completed; its message is the cause recorded in the missing table."""


def axis_elements(config: CombinatoricsConfig) -> Dict[str, Dict]:
    """Per axis: the elements it contributes, iterated and static."""
    elements = {}
    for name, axis in config.axes.items():
        iterated, static, static_first = load_axis_data(axis.to_dict())
        elements[name] = {"mode": axis.mode, "iterated": iterated, "static": static,
                          "static_first": static_first}
    return elements


def plan_iteration(elements: Dict[str, Dict]) -> Tuple[List, List]:
    """Split the axes into those that iterate and those that are bundled whole."""
    iteration, static_only = [], []
    for name in [n for n in AXIS_ORDER if n in elements]:
        data = elements[name]
        if data["mode"] == "each":
            items = data["iterated"] or data["static"]
            extra = data["static"] if data["iterated"] else []
            iteration.append((name, items, extra, data["static_first"]))
        elif data["iterated"]:
            iteration.append((name, data["iterated"], data["static"], data["static_first"]))
        else:
            static_only.append((name, data["static"]))
    return iteration, static_only


def build_runs(elements: Dict[str, Dict]) -> List[Dict]:
    """One dict per run: its id (the framework's own, as the wrapper predicted it), its records per axis and its provenance."""
    iteration, static_only = plan_iteration(elements)
    combinations = (list(itertools.product(*[list(enumerate(items)) for _, items, _, _ in iteration]))
                    if iteration else [()])
    runs = []
    for combination in combinations:
        records, selections = {}, {}
        for index, (name, items, extra, static_first) in enumerate(iteration):
            position, record = combination[index]
            records[name] = (list(extra) if static_first else []) + [record] + ([] if static_first else list(extra))
            selections[name] = ("each", [str(item.get("id", "?")) for item in items], position,
                                [str(item.get("id", "?")) for item in extra], bool(static_first))
        for name, static in static_only:
            records[name] = list(static)
            selections[name] = ("bundle", [str(item.get("id", "?")) for item in static], None, [], False)
        if not selections:
            continue
        runs.append({
            "id": predict_single_output_id(**selections),
            "records": records,
            "provenance": {name: "+".join(str(r.get("id", "?")) for r in recs) for name, recs in records.items()},
        })
    return runs


def polymer_sequence(record: Dict) -> str:
    """The residues of one protein record, refusing an empty one by its id (an empty cell arrives as the float nan)."""
    sequence = record.get("sequence", "")
    if not isinstance(sequence, str) or not sequence.strip():
        raise RunError(f"{record.get('id', '?')}: empty sequence ({sequence!r})")
    return sequence.strip()


def ligand_entity(record: Dict, asym_id: str, prefer_ccd: bool) -> Dict:
    """The upstream ligand entity: a CCD code or SMILES, chosen as Boltz2 and Protenix choose.

    prefer_ccd serves a covalent bond, which names the ligand's residue and so needs its CCD code whenever the row has one.
    """
    try:
        field, value = auth_ligand_field(record)
        if prefer_ccd and ligand_chemistry(record)["ccd"]:
            field, value = "ccd", ligand_chemistry(record)["ccd"]
    except ValueError as exc:
        raise RunError(str(exc))
    return {"entity_type": "ligand", ("ccd_code" if field == "ccd" else "smiles"): value,
            "asym_ids": [asym_id]}


def covalent_bond(link: Dict, entities: List[Dict]) -> Dict:
    """An upstream `bonds` entry from the wrapper's covalent_linkage dict."""
    chain = str(link["chain"])
    protein = next((e for e in entities if e["entity_type"] == "protein" and chain in e["asym_ids"]), None)
    if protein is None:
        raise RunError(f"covalent_linkage names chain {chain}, which is not a protein chain of this system")
    position = int(link["position"])
    sequence = protein["sequence"]
    if not 1 <= position <= len(sequence):
        raise RunError(f"covalent_linkage position {position} is outside chain {chain} ({len(sequence)} residues)")
    residue = THREE_LETTER.get(sequence[position - 1].upper())
    if residue is None:
        raise RunError(f"covalent_linkage: residue {position} of chain {chain} is {sequence[position - 1]!r}, "
                       "which has no three-letter name here")
    ligand = next((e for e in entities if e["entity_type"] == "ligand"), None)
    if ligand is None:
        raise RunError("covalent_linkage needs a ligand in the system")
    if "ccd_code" not in ligand:
        raise RunError("covalent_linkage needs the ligand as a CCD code, and this ligand has only a SMILES")
    return {
        "atom1": {"name": link["protein_atom"], "res_name": residue, "res_id": position, "asym_id": chain},
        "atom2": {"name": link["ligand_atom"], "res_name": ligand["ccd_code"], "res_id": 1,
                  "asym_id": ligand["asym_ids"][0]},
    }


def build_system(records: Dict[str, List[Dict]], contacts: Optional[List[Dict]],
                 covalent: Optional[Dict]) -> Tuple[Dict, List[Tuple[Dict, str]]]:
    """(system dict, [(entity, record id)]) for one run; chains are A, B, ... in axis order."""
    chains = iter(string.ascii_uppercase)
    entities, sources = [], []
    for axis in AXIS_ORDER:
        for record in records.get(axis, []):
            for member in record.get("__members__", [record]):
                asym_id = next(chains, None)
                if asym_id is None:
                    raise RunError(f"a system holds at most {len(string.ascii_uppercase)} chains")
                if axis == "sequences":
                    entity = {"entity_type": "protein", "sequence": polymer_sequence(member),
                              "asym_ids": [asym_id]}
                else:
                    entity = ligand_entity(member, asym_id, prefer_ccd=bool(covalent))
                entities.append(entity)
                sources.append((entity, str(member.get("id", "?"))))
    system = {"entities": entities}
    if contacts:
        system["contacts"] = copy.deepcopy(contacts)
    if covalent:
        system["bonds"] = [covalent_bond(covalent, entities)]
    return system, sources


def write_system_yaml(path: str, system: Dict) -> None:
    with open(path, "w") as handle:
        yaml.safe_dump(system, handle, sort_keys=False)


def build_command(prefix: str, system_yaml: str, out_dir: str, config: Dict) -> List[str]:
    command = container_argv_prefix(prefix) + [
        "bioemu2", "sample-equilibrium",
        "--system", system_yaml,
        "--out", out_dir,
        "--num-samples", str(config["num_samples"]),
        "--batch-size", str(config["batch_size"]),
        "--sampler", str(config["sampler"]),
        "--physicality-filter", str(config["physicality_filter"]),
    ]
    if config.get("seed") is not None:
        command += ["--seed", str(config["seed"])]
    return command


def msa_cache_entry(sequence: str, a3m: str) -> Tuple[str, Dict]:
    """(file name, JSON content) of the cache entry bioemu2's ColabFold client reads for one protein chain.

    The key is the sha256 of the sorted-keys JSON of the request, the checksum is that of the JSON list of alignments, and a hit is refused when the first record of the alignment is not the sequence.
    """
    request = {"sequences": [sequence], "use_pairing": False, "pairing_strategy": "greedy", "use_env": True}
    key = hashlib.sha256(json.dumps(request, sort_keys=True).encode("utf-8")).hexdigest()
    entry = {"request": request, "a3ms": [a3m],
             "sha256": hashlib.sha256(json.dumps([a3m]).encode("utf-8")).hexdigest()}
    return f"{key}.json", entry


def a3m_query(a3m: str) -> str:
    """The first record of an alignment with insertions (lower case) and gaps removed."""
    residues = []
    for line in a3m.splitlines():
        if line.startswith(">"):
            if residues:
                break
            continue
        residues.append(line.strip())
    return "".join(c for c in "".join(residues) if c.isupper())


def a3m_from_boltz_csv(path: str) -> str:
    """A Boltz MSA csv (`key,sequence`, the query first) as a3m text."""
    frame = pd.read_csv(path, dtype=str)
    lines = []
    for index, sequence in enumerate(frame["sequence"].tolist()):
        lines += [">101" if index == 0 else f">hit_{index}", str(sequence)]
    return "\n".join(lines) + "\n"


def load_alignment(path: str) -> str:
    if not os.path.isfile(path):
        raise RunError(f"the alignment {path} does not exist")
    extension = os.path.splitext(path)[1].lower()
    if extension == ".a3m":
        with open(path) as handle:
            return handle.read()
    if extension == ".csv":
        return a3m_from_boltz_csv(path)
    raise RunError(f"the alignment {path} is neither .a3m nor .csv")


def stage_msa(sources: List[Tuple[Dict, str]], msas: MsaLookup, cache_dir: str) -> List[str]:
    """Write the supplied alignment of a one-protein system into the step's MSA cache; returns the ids that got one.

    A protein with no alignment among the msas is left to bioemu2's own ColabFold search, as documented, and the log says so.
    """
    proteins = [(entity, rid) for entity, rid in sources if entity["entity_type"] == "protein"]
    supplied = []
    for entity, rid in proteins:
        path = msas.path(rid, entity["sequence"])
        if path is None:
            print(f"  WARNING: {msas.cause(rid)}; it is searched through ColabFold")
        else:
            supplied.append((entity, rid, path))
    if not supplied:
        return []
    if len(proteins) > 1:
        raise RunError(f"msas serves a system with one protein chain, and this one has {len(proteins)}")
    (entity, rid, path), = supplied
    a3m = load_alignment(os.path.abspath(path))
    if a3m_query(a3m) != entity["sequence"].upper():
        raise RunError(f"the first record of the alignment for {rid} is not its sequence")
    name, entry = msa_cache_entry(entity["sequence"], a3m)
    os.makedirs(cache_dir, exist_ok=True)
    with open(os.path.join(cache_dir, name), "w") as handle:
        json.dump(entry, handle)
    return [rid]


def last_informative_line(lines) -> str:
    """The last line of upstream's output that carries words, whitespace collapsed (blank lines, bar redraws and rules do not)."""
    for line in reversed(list(lines)):
        text = " ".join(line.split())
        if any(character.isalnum() for character in text):
            return text
    return ""


def sample(command: List[str], env: Dict[str, str]) -> None:
    print(f"  running: {' '.join(command)}")
    tail = collections.deque(maxlen=OUTPUT_TAIL)
    # Unbuffered, so a traceback is not followed by stdout that python flushed at exit; newline="" keeps upstream's own line endings.
    with subprocess.Popen(command, env={**env, "PYTHONUNBUFFERED": "1"}, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT) as process:
        for line in io.TextIOWrapper(process.stdout, encoding="utf-8", errors="replace", newline=""):
            sys.stdout.write(line)
            sys.stdout.flush()
            tail.append(line)
    if process.returncode != 0:
        cause = f"bioemu2 sample-equilibrium exited with status {process.returncode}"
        detail = last_informative_line(tail)
        if detail:
            cause = f"{cause}: {detail[:max(CAUSE_LIMIT - len(cause) - 2, 0)]}"
        raise RunError(cause)


def number_atoms(topology) -> None:
    """Give every atom an integer serial: mdtraj reads them from a CIF as text, and its PDB writer takes only integers."""
    for atom in topology.atoms:
        atom.serial = atom.index + 1


def split_trajectory(topology_cif: str, samples_xtc: str, run_id: str,
                     structures_dir: str) -> Tuple[List[Tuple[str, str]], int]:
    """One PDB per sampled frame as `<run>_<k>.pdb`; returns ([(structure id, path)], residues in the system)."""
    import mdtraj as md

    trajectory = md.load(samples_xtc, top=topology_cif)
    number_atoms(trajectory.topology)
    os.makedirs(structures_dir, exist_ok=True)
    frames = []
    for index in range(trajectory.n_frames):
        frame_id = f"{run_id}_{index + 1}"
        path = os.path.join(structures_dir, f"{frame_id}.pdb")
        trajectory[index].save_pdb(path)
        frames.append((frame_id, path))
    return frames, int(trajectory.n_residues)


def read_physicality(path: str) -> Dict[int, Dict]:
    """{frame: {status, reason}} from bioemu2's physicality.json, empty when the file is absent or unreadable."""
    if not os.path.isfile(path):
        return {}
    try:
        with open(path) as handle:
            records = json.load(handle)
        return {int(record["frame"]): record for record in records}
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print(f"WARNING: cannot read {path}: {exc}", file=sys.stderr)
        return {}


def process_run(run: Dict, config: Dict, paths: argparse.Namespace,
                msas: Optional[MsaLookup]) -> Dict:
    """Sample one run and stage its outputs; the rows each output gets are in the returned dict."""
    run_id = run["id"]
    system, sources = build_system(run["records"], config.get("contacts"), config.get("covalent_linkage"))

    env = dict(os.environ)
    if msas is not None:
        staged = stage_msa(sources, msas, paths.msa_cache_dir)
        if staged:
            env["BIOEMU2_MSA_CACHE"] = paths.msa_cache_dir
            print(f"  alignment supplied for {', '.join(staged)}")

    work = os.path.join(paths.runs_dir, run_id)
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    system_yaml = os.path.join(work, "system.yaml")
    write_system_yaml(system_yaml, system)
    sample(build_command(paths.container_prefix, system_yaml, work, config), env)

    topology_cif = os.path.join(work, "topology.cif")
    samples_xtc = os.path.join(work, "samples.xtc")
    if not (os.path.isfile(topology_cif) and os.path.isfile(samples_xtc)):
        raise RunError("bioemu2 wrote no samples.xtc and topology.cif")

    frames, n_residues = split_trajectory(topology_cif, samples_xtc, run_id, paths.structures_dir)
    if len(frames) != config["num_samples"]:
        print(f"WARNING: {run_id}: {len(frames)} conformers written, {config['num_samples']} requested",
              file=sys.stderr)

    os.makedirs(paths.trajectories_dir, exist_ok=True)
    destination_xtc = os.path.join(paths.trajectories_dir, f"{run_id}.xtc")
    destination_top = os.path.join(paths.trajectories_dir, f"{run_id}_topology.cif")
    shutil.copyfile(samples_xtc, destination_xtc)
    shutil.copyfile(topology_cif, destination_top)

    report = read_physicality(os.path.join(work, "physicality.json"))
    provenance = {"sequences.id": run["provenance"].get("sequences", ""),
                  "ligands.id": run["provenance"].get("ligands", "")}
    return {
        "structures": [{"id": fid, "file": path, "run.id": run_id, **provenance} for fid, path in frames],
        "trajectory": {"id": run_id, "file": destination_xtc, "topology": destination_top},
        "physicality": [{"id": fid, "run.id": run_id, "frame": index,
                         "status": report.get(index, {}).get("status"),
                         "reason": report.get(index, {}).get("reason")}
                        for index, (fid, _) in enumerate(frames)],
        "summary": {"id": run_id, **provenance, "n_samples": len(frames), "n_residues": n_residues},
    }


def write_table(path: str, columns: List[str], rows: List[Dict]) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    frame = pd.DataFrame(rows, columns=columns)
    # A blank next to integers would make pandas write 3.0; the nullable type keeps 3 and writes the blank empty.
    for column in INTEGER_COLUMNS:
        if column in frame.columns:
            frame[column] = frame[column].astype("Int64")
    frame.to_csv(path, index=False)


def main():
    parser = argparse.ArgumentParser(description="Sample BioEmu2 equilibrium ensembles")
    parser.add_argument("--config-yaml", required=True)
    parser.add_argument("--combinatorics-config", required=True)
    parser.add_argument("--runs-dir", required=True)
    parser.add_argument("--structures-dir", required=True)
    parser.add_argument("--structures-map-csv", required=True)
    parser.add_argument("--trajectories-dir", required=True)
    parser.add_argument("--trajectories-map-csv", required=True)
    parser.add_argument("--summary-csv", required=True)
    parser.add_argument("--physicality-csv", required=True)
    parser.add_argument("--local-missing-csv", required=True)
    parser.add_argument("--container-prefix", default="")
    parser.add_argument("--msas-json", default=None)
    parser.add_argument("--msa-cache-dir", default=None)
    args = parser.parse_args()

    with open(args.config_yaml) as handle:
        config = yaml.safe_load(handle) or {}

    runs = build_runs(axis_elements(CombinatoricsConfig.load(args.combinatorics_config)))
    msas = MsaLookup.load(args.msas_json) if args.msas_json else None
    if msas is not None and not args.msa_cache_dir:
        parser.error("--msas-json needs --msa-cache-dir")
    step_id = step_id_from_table_path(args.local_missing_csv)

    structure_rows, trajectory_rows, physicality_rows, summary_rows, missing_rows = [], [], [], [], []
    attempted = failures = 0
    for run in runs:
        attempted += 1
        print(f"Run {run['id']}")
        try:
            outcome = process_run(run, config, args, msas)
        except Exception as exc:
            cause = str(exc) if isinstance(exc, RunError) else f"{type(exc).__name__}: {exc}"
            print(f"WARNING: {run['id']} failed: {cause}", file=sys.stderr)
            failures += 1
            summary_rows.append({"id": run["id"], "sequences.id": run["provenance"].get("sequences", ""),
                                 "ligands.id": run["provenance"].get("ligands", "")})
            missing_rows.append({"id": run["id"], "removed_by": step_id, "kind": "failure",
                                 "cause": cause[:CAUSE_LIMIT]})
            continue
        structure_rows += outcome["structures"]
        trajectory_rows.append(outcome["trajectory"])
        physicality_rows += outcome["physicality"]
        summary_rows.append(outcome["summary"])
        print(f"  {run['id']}: {outcome['summary']['n_samples']} conformers, "
              f"{outcome['summary']['n_residues']} residues")

    write_table(args.structures_map_csv, ["id", "file", "run.id", "sequences.id", "ligands.id"], structure_rows)
    write_table(args.trajectories_map_csv, ["id", "file", "topology"], trajectory_rows)
    write_table(args.summary_csv, SUMMARY_COLUMNS, summary_rows)
    write_table(args.physicality_csv, PHYSICALITY_COLUMNS, physicality_rows)
    write_table(args.local_missing_csv, MISSING_COLUMNS, missing_rows)
    print(f"Structures map: {args.structures_map_csv} ({len(structure_rows)} rows)")
    print(f"Trajectories map: {args.trajectories_map_csv} ({len(trajectory_rows)} rows)")
    print(f"Summary: {args.summary_csv} ({len(summary_rows)} rows)")
    print(f"Physicality: {args.physicality_csv} ({len(physicality_rows)} rows)")
    print(f"Local failures: {failures}")

    # Nothing attempted means the upstream stream was empty; every run failing means the step did not work.
    if not attempted:
        print("ERROR: no runs to sample (an input axis resolved to zero records)", file=sys.stderr)
        sys.exit(1)
    if failures == attempted:
        print(f"ERROR: all {attempted} attempted run(s) failed", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
