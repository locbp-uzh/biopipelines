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
import re
import shutil
import sys
from typing import Dict, List, Optional, Tuple

_biopipelines_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'biopipelines')
sys.path.insert(0, _biopipelines_dir)
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from combinatorics import CombinatoricsConfig, predict_single_output_id
from biopipelines.biopipelines_io import write_skipped_queries
from biopipelines.id_map_utils import MissingAlignment, MsaLookup
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


ALIGNMENT_EXTENSIONS = (".a3m", ".sto", ".npz")


def stage_alignment(path: str, key: str, directory: str) -> str:
    """Copy an .a3m or .sto alignment to <directory>/<key>/colabfold_main.<ext>, the name upstream parses.

    Upstream reads an alignment file only when its stem is one of the databases it names (`colabfold_main`, `uniref90_hits`, ...) and skips any other, which leaves the chain with no alignment and raises an IndexError.
    A `.npz` is already parsed and is passed through unchanged.
    """
    extension = os.path.splitext(path)[1].lower()
    if extension not in ALIGNMENT_EXTENSIONS:
        raise ValueError(
            f"OpenFold3 reads precomputed alignments as .a3m, .sto or .npz; the msas stream "
            f"carries {extension or 'a file with no extension'} ({path}). Generate the MSA with "
            f"MMseqs2(output_format=\"a3m\") or leave use_msa_server on.")
    if extension == ".npz":
        return path
    folder = os.path.join(directory, re.sub(r"[^A-Za-z0-9_.+-]", "_", key))
    os.makedirs(folder, exist_ok=True)
    staged = os.path.join(folder, f"colabfold_main{extension}")
    shutil.copyfile(path, staged)
    return staged


class Alignments:
    """The supplied alignments, each staged once under the name upstream reads when its protein first asks."""

    def __init__(self, lookup: MsaLookup, directory: str):
        self.lookup = lookup
        self.directory = directory
        self._staged: Dict[str, str] = {}

    def paths_for(self, protein_id: str, sequence: Optional[str]) -> List[str]:
        """The staged alignment of one protein chain; raises MissingAlignment when it has none."""
        if protein_id not in self._staged:
            path = self.lookup.require(protein_id, sequence)
            if not os.path.exists(path):
                raise MissingAlignment(f"the alignment {path} of {protein_id} does not exist")
            self._staged[protein_id] = stage_alignment(path, protein_id, self.directory)
        return [self._staged[protein_id]]


def build_chain(entity_type: str, record: Dict, label: str, alignments: Optional[Alignments],
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
        if alignments is not None:
            chain["main_msa_file_paths"] = alignments.paths_for(record_id, sequence)
        if template and (not template_chains or label in template_chains):
            chain["template_cif_paths"] = [template]
            chain["template_cif_chain_ids"] = [None]
    return chain


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


def build_queries(elements: Dict[str, Dict], alignments: Optional[Alignments],
                  template: Optional[str], template_chains: Optional[List[str]]
                  ) -> Tuple[Dict, Dict[str, str]]:
    """({"queries": {query id: query}}, {query id: cause}) for the queries written and the ones skipped.

    A query with a protein chain that has no supplied alignment is skipped with that cause, so it folds nowhere instead of folding single-sequence.
    """
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
    skipped: Dict[str, str] = {}
    # One empty combination only when nothing iterates; an iterating axis with no records yields no query, not an IndexError.
    combinations = (list(itertools.product(*[list(enumerate(items))
                                             for _, _, items, _, _ in iteration]))
                    if iteration else [()])

    for combination in combinations:
        chains: List[Dict] = []
        causes: List[str] = []
        counter = 0

        def add(entity_type, record):
            nonlocal counter
            for member in record.get("__members__", [record]):
                label = chain_ids(1, counter)[0]
                counter += 1
                try:
                    chains.append(build_chain(entity_type, member, label, alignments,
                                              template, template_chains))
                except MissingAlignment as exc:
                    causes.append(str(exc))
                if entity_type in DOUBLE_STRANDED_ENTITY_TYPES:
                    label2 = chain_ids(1, counter)[0]
                    counter += 1
                    chains.append(build_chain(entity_type, member, label2, alignments,
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
        if causes:
            skipped[query_id] = "; ".join(dict.fromkeys(causes))
        else:
            queries[query_id] = {"chains": chains}

    return {"queries": queries}, skipped


def main():
    parser = argparse.ArgumentParser(description="Build the OpenFold3 query JSON")
    parser.add_argument("--combinatorics-config", required=True)
    parser.add_argument("--queries-json", required=True)
    parser.add_argument("--skipped-csv", required=True)
    parser.add_argument("--msas-json")
    parser.add_argument("--template")
    parser.add_argument("--template-chains")
    args = parser.parse_args()

    config = CombinatoricsConfig.load(args.combinatorics_config)
    alignments = None
    if args.msas_json:
        alignments = Alignments(MsaLookup.load(args.msas_json),
                                os.path.join(os.path.dirname(os.path.abspath(args.queries_json)), "msas"))
    template_chains = args.template_chains.split(",") if args.template_chains else None

    payload, skipped = build_queries(axis_elements(config), alignments, args.template, template_chains)
    write_skipped_queries(args.skipped_csv, skipped)
    for query_id, cause in skipped.items():
        print(f"WARNING: query {query_id} skipped: {cause}")
    # Written even when empty, so a stale file from an earlier run is never mistaken for this one's queries.
    os.makedirs(os.path.dirname(os.path.abspath(args.queries_json)), exist_ok=True)
    with open(args.queries_json, "w") as handle:
        json.dump(payload, handle, indent=2)
    print(f"Wrote {len(payload['queries'])} quer(y/ies) to {args.queries_json}; {len(skipped)} skipped")
    if not payload["queries"]:
        print("ERROR: every query was skipped" if skipped
              else "ERROR: no queries were built — every input axis resolved to zero records")
        sys.exit(1)


if __name__ == "__main__":
    main()
