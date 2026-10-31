#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Build RF3's inputs JSON from the combinatorics config.

This runs at execution time because the sequences do not exist before it: they arrive in the
map_tables the previous step wrote. The wrapper only declares the pairing; this resolves it.

One JSON list holds every example, because `rf3 fold` predicts a whole list in one process and so
loads the model once rather than once per complex. Each example is `{"name": <id>, "components": [...]}`,
and `name` is the id the postprocess finds the outputs under.
"""

import argparse
import json
import os
import sys
from typing import Dict, List, Optional, Tuple

_biopipelines_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'biopipelines')
sys.path.insert(0, _biopipelines_dir)
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from combinatorics import CombinatoricsConfig, predict_single_output_id, predict_single_provenance
from biopipelines.biopipelines_io import write_skipped_queries
from biopipelines.id_map_utils import MissingAlignment, MsaLookup
from ligand_utils import auth_ligand_field
from nucleic_acids import DOUBLE_STRANDED_ENTITY_TYPES, reverse_complement
from pipe_axis_data import load_axis_data
from pipe_query_provenance import write_query_provenance

# The order chains are written in, so an example's chain ids are stable across runs.
AXIS_ORDER = ("proteins", "ssDNA", "dsDNA", "ssRNA", "dsRNA", "ligands")

# mmCIF entity_poly.type strings that AtomWorks' ChainType accepts; given explicitly because "ACGT" also reads as a protein.
NUCLEIC_ACID_CHAIN_TYPES = {
    "ssdna": "polydeoxyribonucleotide",
    "dsdna": "polydeoxyribonucleotide",
    "ssrna": "polyribonucleotide",
    "dsrna": "polyribonucleotide",
}

MSA_EXTENSIONS = (".a3m", ".a3m.gz", ".fasta")


class UnreadableAlignment(ValueError):
    """An alignment RF3 cannot read; the example it belongs to is skipped with this message as its cause."""


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

    An empty FASTA read through pandas arrives here as the float nan. Left alone it becomes an
    example with no residues, and the failure surfaces deep inside the model with nothing naming
    the input that caused it.
    """
    sequence = record.get("sequence", "")
    if not isinstance(sequence, str) or not sequence.strip():
        raise ValueError(
            f"{entity_type} '{record_id}' has an empty sequence ({sequence!r}) — "
            f"check the table or FASTA it came from")
    return sequence.strip()


def build_component(entity_type: str, record: Dict, label: str, msas: Optional[MsaLookup],
                    rev_comp: bool = False) -> Dict:
    record_id = str(record.get("id", "?"))

    if entity_type == "ligand":
        # CCD-vs-SMILES is decided from the row's provenance, as Boltz2 does, so both models get the same chemistry.
        field, value = auth_ligand_field(record)
        return {"ccd_code": value, "chain_id": label} if field == "ccd" \
            else {"smiles": value, "chain_id": label}

    sequence = polymer_sequence(record, entity_type, record_id)
    if rev_comp:
        sequence = reverse_complement(sequence, entity_type)
    component: Dict = {"seq": sequence, "chain_id": label}

    if entity_type in NUCLEIC_ACID_CHAIN_TYPES:
        component["chain_type"] = NUCLEIC_ACID_CHAIN_TYPES[entity_type]
    elif entity_type == "protein" and msas is not None:
        path = msas.require(record_id, sequence)
        if not path.lower().endswith(MSA_EXTENSIONS):
            raise UnreadableAlignment(
                f"RF3 reads precomputed alignments as {', '.join(MSA_EXTENSIONS)}; the msas stream "
                f"carries {path}. Generate the MSA with MMseqs2(output_format=\"a3m\").")
        if not os.path.exists(path):
            raise MissingAlignment(f"the alignment {path} of {record_id} does not exist")
        component["msa_path"] = path
    return component


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


def build_examples(elements: Dict[str, Dict], msas: Optional[MsaLookup],
                   provenance: Optional[Dict[str, Dict[str, str]]] = None
                   ) -> Tuple[List[Dict], Dict[str, str]]:
    """([example, ...], {example id: cause}) for the examples written and the ones skipped.

    An example with a protein chain whose alignment is absent, not on disk or not readable is skipped with that cause, so it is predicted nowhere instead of single-sequence.
    When a dict is passed as `provenance`, every example that entered, skipped or not, is recorded in it with the input each axis contributed.
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

    examples: List[Dict] = []
    skipped: Dict[str, str] = {}
    # One empty combination only when nothing iterates; an iterating axis with no records yields no example, not an IndexError.
    combinations = (list(itertools.product(*[list(enumerate(items))
                                             for _, _, items, _, _ in iteration]))
                    if iteration else [()])

    for combination in combinations:
        components: List[Dict] = []
        causes: List[str] = []
        counter = 0

        def add(entity_type, record):
            nonlocal counter
            for member in record.get("__members__", [record]):
                label = chain_ids(1, counter)[0]
                counter += 1
                try:
                    components.append(build_component(entity_type, member, label, msas))
                except (MissingAlignment, UnreadableAlignment) as exc:
                    causes.append(str(exc))
                if entity_type in DOUBLE_STRANDED_ENTITY_TYPES:
                    label2 = chain_ids(1, counter)[0]
                    counter += 1
                    components.append(build_component(entity_type, member, label2, msas,
                                                      rev_comp=True))

        # The id comes from the framework's helper so it matches what the wrapper predicted at configuration time.
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
        name = predict_single_output_id(**axis_selections)
        if provenance is not None:
            provenance[name] = predict_single_provenance(**axis_selections)
        if causes:
            skipped[name] = "; ".join(dict.fromkeys(causes))
        else:
            examples.append({"name": name, "components": components})

    return examples, skipped


def main():
    parser = argparse.ArgumentParser(description="Build the RF3 inputs JSON")
    parser.add_argument("--combinatorics-config", required=True)
    parser.add_argument("--queries-json", required=True)
    parser.add_argument("--skipped-csv", required=True)
    parser.add_argument("--provenance-csv")
    parser.add_argument("--msas-json")
    args = parser.parse_args()

    config = CombinatoricsConfig.load(args.combinatorics_config)
    msas = MsaLookup.load(args.msas_json) if args.msas_json else None
    provenance: Dict[str, Dict[str, str]] = {}
    examples, skipped = build_examples(axis_elements(config), msas, provenance)
    write_skipped_queries(args.skipped_csv, skipped)
    if args.provenance_csv:
        write_query_provenance(args.provenance_csv, provenance)
    for name, cause in skipped.items():
        print(f"WARNING: example {name} skipped: {cause}")

    # Written even when empty, so a stale file from an earlier run is never mistaken for this one's examples.
    os.makedirs(os.path.dirname(os.path.abspath(args.queries_json)), exist_ok=True)
    with open(args.queries_json, "w") as handle:
        json.dump(examples, handle, indent=2)
    print(f"Wrote {len(examples)} example(s) to {args.queries_json}; {len(skipped)} skipped")
    if not examples:
        print("ERROR: every example was skipped" if skipped
              else "ERROR: no examples were built — every input axis resolved to zero records")
        sys.exit(1)


if __name__ == "__main__":
    main()
