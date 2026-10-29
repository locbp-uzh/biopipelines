#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Build one Protenix input JSON per query from the combinatorics config.

This runs at execution time because the sequences do not exist before it: they arrive in the map_tables the previous step wrote.
The wrapper only declares the pairing; this resolves it.

Each query is its own file in one folder, because `protenix pred -i <folder>` loads the model once and keeps a failure to the file it came from.
A query that cannot be written (a bad residue letter, a protein chain with no alignment, a bond naming an entity it does not have) is left out and recorded in --skipped-csv with the reason.
The postprocess turns those into rows of the confidence and missing tables.
"""

import argparse
import csv
import itertools
import json
import os
import shutil
import sys
from typing import Dict, List, Optional, Tuple

_biopipelines_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'biopipelines')
sys.path.insert(0, _biopipelines_dir)
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pandas as pd

from combinatorics import CombinatoricsConfig, predict_single_output_id
from ligand_utils import auth_ligand_field
from nucleic_acids import DOUBLE_STRANDED_ENTITY_TYPES, reverse_complement
from pipe_axis_data import load_axis_data

# Entity order is chain order, so it is fixed: covalent_bonds names entities by their position in it.
AXIS_ORDER = ("proteins", "ssDNA", "dsDNA", "ssRNA", "dsRNA", "ligands")

ENTITY_KEYS = {
    "protein": "proteinChain",
    "ssdna": "dnaSequence",
    "dsdna": "dnaSequence",
    "ssrna": "rnaSequence",
    "dsrna": "rnaSequence",
}

# The letters upstream's input format allows per entity type.
ALLOWED_LETTERS = {
    "proteinChain": set("ACDEFGHIKLMNPQRSTVWYX"),
    "dnaSequence": set("ATGCN"),
    "rnaSequence": set("AUGCN"),
}


class QueryError(ValueError):
    """A query that cannot be written; its message is the cause recorded in the missing table."""


def polymer_sequence(record: Dict, entity_key: str, record_id: str) -> str:
    """The residues of one polymer record, refusing an empty or out-of-alphabet one by its id.

    An empty FASTA read through pandas arrives as the float nan, which would otherwise become a query with no residues.
    """
    sequence = record.get("sequence", "")
    if not isinstance(sequence, str) or not sequence.strip():
        raise QueryError(f"{record_id}: empty sequence ({sequence!r})")
    sequence = sequence.strip()
    invalid = sorted(set(sequence) - ALLOWED_LETTERS[entity_key])
    if invalid:
        raise QueryError(
            f"{record_id}: letter(s) {''.join(invalid)} are not allowed in a {entity_key} "
            f"(allowed: {''.join(sorted(ALLOWED_LETTERS[entity_key]))})")
    return sequence


def load_msa_lookup(msas_json: Optional[str]) -> Dict[str, str]:
    """{protein id: absolute alignment path} from the msas stream, when one was passed.

    A row is keyed by its `id` and by the `sequences.id` it carries: MMseqs2, Boltz2 and AlphaFold name their rows `<id>_msa` and keep the protein id in `sequences.id`.
    """
    if not msas_json:
        return {}
    with open(msas_json) as handle:
        stream = json.load(handle)
    table = stream.get("map_table")
    if not table or not os.path.exists(table):
        return {}
    frame = pd.read_csv(table, dtype={"id": str, "sequences.id": str})
    column = next((c for c in ("file", "path") if c in frame.columns), None)
    if column is None:
        return {}
    lookup = {}
    for key in ("id", "sequences.id"):
        if key in frame.columns:
            lookup.update({str(row[key]): os.path.abspath(str(row[column])) for _, row in frame.iterrows()
                           if isinstance(row.get(column), str) and not pd.isna(row[key])})
    unsupported = sorted({os.path.splitext(p)[1].lower() for p in lookup.values()} - {".a3m"})
    if unsupported:
        raise ValueError(
            f"Protenix reads precomputed alignments as .a3m; the msas stream carries "
            f"{', '.join(unsupported)}. Generate the MSA with MMseqs2(output_format=\"a3m\").")
    return lookup


def build_entity(entity_type: str, record: Dict, msa_lookup: Optional[Dict[str, str]],
                 rev_comp: bool = False) -> Dict:
    """One upstream entity (count 1) for one record.

    msa_lookup is None when no alignments were supplied, and a dict (possibly empty) when they were: then a protein without one is an error, since a path that is absent would make upstream search the server.
    """
    record_id = str(record.get("id", "?"))

    if entity_type == "ligand":
        # Same CCD-versus-SMILES decision as Boltz2, so both models are handed the same chemistry.
        field, value = auth_ligand_field(record)
        return {"ligand": {"ligand": f"CCD_{value}" if field == "ccd" else value, "count": 1}}

    entity_key = ENTITY_KEYS[entity_type]
    sequence = polymer_sequence(record, entity_key, record_id)
    if rev_comp:
        sequence = reverse_complement(sequence, entity_type)
    entity = {"sequence": sequence, "count": 1}

    if entity_key == "proteinChain" and msa_lookup is not None:
        path = msa_lookup.get(record_id)
        if not path:
            raise QueryError(f"{record_id}: no alignment for this protein chain in the msas stream")
        if not os.path.exists(path):
            raise QueryError(f"{record_id}: the alignment {path} does not exist")
        entity["unpairedMsaPath"] = path
    return {entity_key: entity}


def axis_elements(config: CombinatoricsConfig) -> Dict[str, Dict]:
    """Per axis: its entity type and the elements it contributes, iterated and static."""
    elements = {}
    for name, axis in config.axes.items():
        iterated, static, static_first = load_axis_data(axis.to_dict())
        elements[name] = {
            "entity_type": axis.entity_type or "protein",
            "mode": axis.mode,
            "iterated": iterated,
            "static": static,
            "static_first": static_first,
        }
    return elements


def plan_iteration(elements: Dict[str, Dict]) -> Tuple[List, List]:
    """Split the axes into those that iterate and those that are bundled whole."""
    ordered = [name for name in AXIS_ORDER if name in elements]
    ordered += [name for name in elements if name not in ordered]

    iteration, static_only = [], []
    for name in ordered:
        data = elements[name]
        if data["mode"] == "each":
            items = data["iterated"] or data["static"]
            extra = data["static"] if data["iterated"] else []
            iteration.append((name, data["entity_type"], items, extra, data["static_first"]))
        elif data["iterated"]:
            iteration.append((name, data["entity_type"], data["iterated"], data["static"],
                              data["static_first"]))
        else:
            static_only.append((name, data["entity_type"], data["static"]))
    return iteration, static_only


def collect_items(combination, iteration, static_only) -> Tuple[List, Dict]:
    """The (entity type, record) pairs of one query in AXIS_ORDER, and its id selections.

    An iterating axis and a bundled one still land in axis order, so the entity numbers `covalent_bonds` uses do not depend on which axes iterate.
    """
    per_axis, selections = {}, {}
    for index, (name, entity_type, iterated, extra, static_first) in enumerate(iteration):
        position, record = combination[index]
        before = [(entity_type, item) for item in extra] if static_first else []
        after = [] if static_first else [(entity_type, item) for item in extra]
        per_axis[name] = before + [(entity_type, record)] + after
        selections[name] = (
            "each" if iterated else "bundle",
            [str(item.get("id", "?")) for item in iterated],
            position,
            [str(item.get("id", "?")) for item in extra],
            bool(static_first),
        )

    for name, entity_type, static in static_only:
        per_axis[name] = [(entity_type, item) for item in static]
        selections[name] = ("bundle", [str(item.get("id", "?")) for item in static], None, [], False)

    names = [name for name in AXIS_ORDER if name in per_axis]
    names += [name for name in per_axis if name not in names]
    return [item for name in names for item in per_axis[name]], selections


def check_bonds(bonds: List[Dict], n_entities: int) -> None:
    """Refuse a bond naming an entity the query does not have; the rest is upstream's to check."""
    for index, bond in enumerate(bonds):
        for key in ("entity1", "entity2"):
            number = int(bond[key])
            if not 1 <= number <= n_entities:
                raise QueryError(
                    f"covalent_bonds[{index}]: {key}={bond[key]} but this query has {n_entities} entities")


def build_job(query_id: str, items: List, msa_lookup: Optional[Dict[str, str]],
              bonds: Optional[List[Dict]]) -> Dict:
    entities = []
    for entity_type, record in items:
        for member in record.get("__members__", [record]):
            entities.append(build_entity(entity_type, member, msa_lookup))
            if entity_type in DOUBLE_STRANDED_ENTITY_TYPES:
                entities.append(build_entity(entity_type, member, msa_lookup, rev_comp=True))
    job = {"name": query_id, "sequences": entities}
    if bonds:
        check_bonds(bonds, len(entities))
        job["covalent_bonds"] = [dict(bond) for bond in bonds]
    return job


def build_queries(elements: Dict[str, Dict], msa_lookup: Optional[Dict[str, str]],
                  bonds: Optional[List[Dict]]) -> Tuple[Dict[str, Dict], Dict[str, str]]:
    """({query id: upstream job}, {query id: cause}) for the queries written and the ones skipped."""
    iteration, static_only = plan_iteration(elements)
    # One empty combination only when nothing iterates; an iterating axis with no records yields no query.
    combinations = (list(itertools.product(*[list(enumerate(iterated))
                                             for _, _, iterated, _, _ in iteration]))
                    if iteration else [()])

    jobs: Dict[str, Dict] = {}
    skipped: Dict[str, str] = {}
    for combination in combinations:
        items, selections = collect_items(combination, iteration, static_only)
        if not selections:
            continue
        # The id comes from the framework's own helper: it is what the wrapper predicted at configuration time.
        query_id = predict_single_output_id(**selections)
        try:
            jobs[query_id] = build_job(query_id, items, msa_lookup, bonds)
        except ValueError as exc:
            skipped[query_id] = str(exc)
    return jobs, skipped


def write_csv(path: str, header: List[str], rows: List[List]) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description="Build the Protenix query JSONs")
    parser.add_argument("--combinatorics-config", required=True)
    parser.add_argument("--queries-dir", required=True)
    parser.add_argument("--query-ids-csv", required=True)
    parser.add_argument("--skipped-csv", required=True)
    parser.add_argument("--msas-json")
    parser.add_argument("--covalent-bonds-json")
    args = parser.parse_args()

    config = CombinatoricsConfig.load(args.combinatorics_config)
    msa_lookup = load_msa_lookup(args.msas_json) if args.msas_json else None
    bonds = None
    if args.covalent_bonds_json:
        with open(args.covalent_bonds_json) as handle:
            bonds = json.load(handle)

    jobs, skipped = build_queries(axis_elements(config), msa_lookup, bonds)
    if not jobs and not skipped:
        print("ERROR: no queries were built; every input axis resolved to zero records")
        sys.exit(1)

    write_csv(args.query_ids_csv, ["id"], [[query_id] for query_id in [*jobs, *skipped]])
    write_csv(args.skipped_csv, ["id", "cause"], [[query_id, cause] for query_id, cause in skipped.items()])
    for query_id, cause in skipped.items():
        print(f"WARNING: query {query_id} skipped: {cause}")

    # A stale folder would hand upstream last run's queries, and its own -update-msa files, as inputs.
    if os.path.isdir(args.queries_dir):
        shutil.rmtree(args.queries_dir)
    os.makedirs(args.queries_dir)
    for query_id, job in jobs.items():
        with open(os.path.join(args.queries_dir, f"{query_id}.json"), "w") as handle:
            json.dump([job], handle, indent=2)

    print(f"Wrote {len(jobs)} quer(y/ies) to {args.queries_dir}; {len(skipped)} skipped")
    # Nothing left to predict is a step-level error: exiting 0 would hide it behind an empty prediction folder.
    if not jobs:
        print("ERROR: every query was skipped")
        sys.exit(1)


if __name__ == "__main__":
    main()
