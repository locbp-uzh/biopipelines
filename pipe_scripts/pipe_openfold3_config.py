#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Build OpenFold3's query JSON from the combinatorics config.

This runs at execution time because the sequences do not exist before it: they arrive in the
map_tables the previous step wrote. The wrapper only declares the pairing; this resolves it.

One JSON holds every query — `run_openfold predict` folds a whole batch in one process, so a
campaign costs one model load rather than one per complex.
"""

import argparse
import json
import os
import sys
from typing import Dict, List, Optional

_biopipelines_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'biopipelines')
sys.path.insert(0, _biopipelines_dir)
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pandas as pd

from combinatorics import CombinatoricsConfig, predict_single_output_id
from ligand_utils import auth_ligand_field
from nucleic_acids import DOUBLE_STRANDED_ENTITY_TYPES, RNA_ENTITY_TYPES, reverse_complement
from pipe_axis_data import load_axis_data

# The order chains are written in, so a query's chain ids are stable across runs. An unstable
# chain order would change `chains_ptm:0` from one run to the next and make two runs of the same
# campaign look different when nothing was.
AXIS_ORDER = ("proteins", "ssDNA", "dsDNA", "ssRNA", "dsRNA", "ligands")

MOLECULE_TYPES = {
    "protein": "protein",
    "ssdna": "dna",
    "dsdna": "dna",
    "ssrna": "rna",
    "dsrna": "rna",
    "ligand": "ligand",
}


def chain_ids(count: int, start: int) -> List[str]:
    """A, B, ... Z, AA, AB — enough labels for a complex with more than 26 chains."""
    labels = []
    for index in range(start, start + count):
        label, n = "", index
        while True:
            label = chr(65 + n % 26) + label
            n = n // 26 - 1
            if n < 0:
                break
        labels.append(label)
    return labels


def polymer_sequence(record: Dict, entity_type: str, record_id: str) -> str:
    """The residues of one polymer record, refusing an empty one loudly.

    An empty FASTA read through pandas arrives here as the float nan. Left alone it becomes a
    query with no residues, and the failure surfaces deep inside the model with nothing naming
    the input that caused it.
    """
    sequence = record.get("sequence", "")
    if not isinstance(sequence, str) or not sequence.strip():
        raise ValueError(
            f"{entity_type} '{record_id}' has an empty sequence ({sequence!r}) — "
            f"check the table or FASTA it came from")
    return sequence.strip()


def msa_paths_for(record_id: str, msa_lookup: Dict[str, str]) -> Optional[List[str]]:
    """The precomputed alignment for one chain, if one was supplied for it."""
    path = msa_lookup.get(str(record_id))
    return [path] if path else None


def build_chain(entity_type: str, record: Dict, label: str, msa_lookup: Dict[str, str],
                template: Optional[str], template_chains: Optional[List[str]],
                rev_comp: bool = False) -> Dict:
    record_id = str(record.get("id", "?"))
    chain: Dict = {"molecule_type": MOLECULE_TYPES[entity_type], "chain_ids": [label]}

    if entity_type == "ligand":
        # auth_ligand_field decides CCD-vs-SMILES from the row's own provenance, the same way
        # Boltz2 does, so the two models are handed the same chemistry for the same input.
        field, value = auth_ligand_field(record)
        if field == "ccd":
            chain["ccd_codes"] = [value]
        else:
            chain["smiles"] = value
        return chain

    sequence = polymer_sequence(record, entity_type, record_id)
    if rev_comp:
        sequence = reverse_complement(sequence, entity_type)
    chain["sequence"] = sequence

    if entity_type == "protein":
        paths = msa_paths_for(record_id, msa_lookup)
        if paths:
            chain["main_msa_file_paths"] = paths
        if template and (not template_chains or label in template_chains):
            chain["template_cif_paths"] = [template]
            chain["template_cif_chain_ids"] = [None]
    return chain


def load_msa_lookup(msas_json: Optional[str]) -> Dict[str, str]:
    """{id: alignment path} from the msas stream, when one was passed."""
    if not msas_json:
        return {}
    with open(msas_json) as handle:
        stream = json.load(handle)
    table = stream.get("map_table")
    if not table or not os.path.exists(table):
        return {}
    frame = pd.read_csv(table, dtype={"id": str})
    column = next((c for c in ("file", "path") if c in frame.columns), None)
    if column is None:
        return {}
    lookup = {str(row["id"]): str(row[column]) for _, row in frame.iterrows()
              if isinstance(row.get(column), str)}
    unsupported = sorted({os.path.splitext(p)[1].lower() for p in lookup.values()}
                         - {".a3m", ".sto", ".npz"})
    if unsupported:
        raise ValueError(
            f"OpenFold3 reads precomputed alignments as .a3m, .sto or .npz; the msas stream "
            f"carries {', '.join(unsupported)}. Generate the MSA with MMseqs2(output_format="
            f"\"a3m\") or leave use_msa_server on.")
    return lookup


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


def build_queries(elements: Dict[str, Dict], msa_lookup: Dict[str, str],
                  template: Optional[str], template_chains: Optional[List[str]]) -> Dict:
    import itertools

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

    queries: Dict[str, Dict] = {}
    # One empty combination only when nothing iterates; an iterating axis with no records yields no query, not an IndexError.
    combinations = (list(itertools.product(*[list(enumerate(items))
                                             for _, _, items, _, _ in iteration]))
                    if iteration else [()])

    for combination in combinations:
        chains: List[Dict] = []
        counter = 0

        def add(entity_type, record):
            nonlocal counter
            for member in record.get("__members__", [record]):
                label = chain_ids(1, counter)[0]
                counter += 1
                chains.append(build_chain(entity_type, member, label, msa_lookup,
                                          template, template_chains))
                if entity_type in DOUBLE_STRANDED_ENTITY_TYPES:
                    label2 = chain_ids(1, counter)[0]
                    counter += 1
                    chains.append(build_chain(entity_type, member, label2, msa_lookup,
                                              template, template_chains, rev_comp=True))

        # The id is composed by the framework's own helper, not by joining strings here: it is
        # what the wrapper predicted at configuration time, and a second spelling of the same
        # rule is a second thing to keep in step.
        axis_selections = {}
        for index, (name, entity_type, items, extra, static_first) in enumerate(iteration):
            position, record = combination[index]
            if static_first:
                for item in extra:
                    add(entity_type, item)
            add(entity_type, record)
            if not static_first:
                for item in extra:
                    add(entity_type, item)
            axis_selections[name] = (
                "each" if items else "bundle",
                [str(item.get("id", "?")) for item in items],
                position,
                [str(item.get("id", "?")) for item in extra],
                bool(static_first),
            )

        for name, entity_type, items in static_only:
            for item in items:
                add(entity_type, item)
            axis_selections[name] = (
                "bundle", [str(item.get("id", "?")) for item in items], None, [], False)

        if not axis_selections:
            continue
        query_id = predict_single_output_id(**axis_selections)
        queries[query_id] = {"chains": chains}

    return {"queries": queries}


def main():
    parser = argparse.ArgumentParser(description="Build the OpenFold3 query JSON")
    parser.add_argument("--combinatorics-config", required=True)
    parser.add_argument("--queries-json", required=True)
    parser.add_argument("--msas-json")
    parser.add_argument("--template")
    parser.add_argument("--template-chains")
    args = parser.parse_args()

    config = CombinatoricsConfig.load(args.combinatorics_config)
    msa_lookup = load_msa_lookup(args.msas_json)
    template_chains = args.template_chains.split(",") if args.template_chains else None

    payload = build_queries(axis_elements(config), msa_lookup, args.template, template_chains)
    if not payload["queries"]:
        print("ERROR: no queries were built — every input axis resolved to zero records")
        sys.exit(1)

    os.makedirs(os.path.dirname(os.path.abspath(args.queries_json)), exist_ok=True)
    with open(args.queries_json, "w") as handle:
        json.dump(payload, handle, indent=2)
    print(f"Wrote {len(payload['queries'])} quer(y/ies) to {args.queries_json}")


if __name__ == "__main__":
    main()
