#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
ESMFold2 inference script.

Reads the combinatorics config written by the ESMFold2 tool, assembles one
complex per output id (proteins / nucleic acids / ligands as chains), folds it
with ESMFold2 over `num_seeds` independent seeds, keeps the best sample by ipTM
(complex) / pLDDT (monomer), and writes one mmCIF per id plus a per-id scores
JSON consumed by the post-processing script.

Runs under the 'esmfold2' conda environment.
"""

import argparse
import json
import os
import sys
import traceback

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.combinatorics import predict_single_output_id

# Double-stranded axes add a reverse-complement chain; ligand axes drive the
# affinity/first-ligand bookkeeping (none here, but kept for symmetry).
DOUBLE_STRANDED = {"dsdna", "dsrna"}
LIGAND_TYPES = {"ligand"}
_DNA_COMPLEMENT = str.maketrans("ACGT", "TGCA")
_RNA_COMPLEMENT = str.maketrans("ACGU", "UGCA")


def reverse_complement(seq, entity_type):
    table = _RNA_COMPLEMENT if "rna" in entity_type else _DNA_COMPLEMENT
    return seq.translate(table)[::-1]


def load_axis_records(axis):
    """Split an axis's source CSV rows into (iterated, static, static_first)."""
    iterated, static = [], []
    min_iter_order, min_static_order = float("inf"), float("inf")
    for source in axis.get("sources", []):
        if isinstance(source, dict):
            path, is_iter, order = source.get("path"), source.get("iterate", True), source.get("order", 0)
        else:
            path, is_iter, order = source, True, 0
        if not path or not os.path.exists(path):
            print(f"WARNING: source file not found: {path}", file=sys.stderr)
            continue
        records = pd.read_csv(path).to_dict("records")
        if is_iter:
            iterated.extend(records)
            min_iter_order = min(min_iter_order, order)
        else:
            static.extend(records)
            min_static_order = min(min_static_order, order)
    return iterated, static, min_static_order < min_iter_order


def build_complexes(config, msa_by_seq):
    """Yield (complex_id, [chain dicts]) for every output id.

    Each chain dict is {entity_type, id (chain letter), sequence|ccd, msa_path}.
    Mirrors the framework cartesian-product assignment so ids match
    predict_output_ids_with_provenance exactly.
    """
    import itertools

    axes = config.get("axes", {})
    iteration_axes, static_only_axes = [], []
    for name, axis in axes.items():
        entity_type = axis.get("entity_type")
        mode = axis.get("mode", "each")
        iterated, static, static_first = load_axis_records(axis)
        if mode == "each":
            items = iterated if iterated else static
            statics = static if iterated else []
            iteration_axes.append(dict(name=name, entity_type=entity_type, mode=mode,
                                       items=items, static=statics, static_first=static_first,
                                       all_ids=[r["id"] for r in items]))
        elif mode == "bundle":
            if iterated:
                iteration_axes.append(dict(name=name, entity_type=entity_type, mode=mode,
                                           items=iterated, static=static, static_first=static_first,
                                           all_ids=[r["id"] for r in iterated]))
            else:
                static_only_axes.append(dict(name=name, entity_type=entity_type, items=static))

    def chain_letter(counter):
        c = chr(65 + counter[0])
        counter[0] += 1
        return c

    def add_item(chains, entity_type, item, counter):
        if entity_type in LIGAND_TYPES:
            chains.append(dict(entity_type=entity_type, id=chain_letter(counter),
                               ccd=item.get("code") or "", smiles=item.get("smiles") or ""))
            return
        seq = str(item.get("sequence", ""))
        msa_path = msa_by_seq.get(seq) if entity_type == "protein" else None
        chains.append(dict(entity_type=entity_type, id=chain_letter(counter),
                           sequence=seq, msa_path=msa_path))
        if entity_type in DOUBLE_STRANDED:
            chains.append(dict(entity_type=entity_type, id=chain_letter(counter),
                               sequence=reverse_complement(seq, entity_type), msa_path=None))

    if not iteration_axes:
        chains, counter, selections = [], [0], {}
        for sa in static_only_axes:
            selections[sa["name"]] = ("bundle", [r["id"] for r in sa["items"]], None, [], False)
            for item in sa["items"]:
                add_item(chains, sa["entity_type"], item, counter)
        yield predict_single_output_id(**selections), chains
        return

    item_lists = [list(enumerate(ia["items"])) for ia in iteration_axes]
    for combo in itertools.product(*item_lists):
        chains, counter, selections = [], [0], {}
        for i, ia in enumerate(iteration_axes):
            idx, item = combo[i]
            static_ids = [r["id"] for r in ia["static"]]
            selections[ia["name"]] = (ia["mode"], ia["all_ids"], idx, static_ids, ia["static_first"])
            if ia["mode"] == "bundle" and ia["static_first"] and ia["static"]:
                for s in ia["static"]:
                    add_item(chains, ia["entity_type"], s, counter)
                add_item(chains, ia["entity_type"], item, counter)
            else:
                add_item(chains, ia["entity_type"], item, counter)
                for s in ia["static"]:
                    add_item(chains, ia["entity_type"], s, counter)
        for sa in static_only_axes:
            selections[sa["name"]] = ("bundle", [r["id"] for r in sa["items"]], None, [], False)
            for item in sa["items"]:
                add_item(chains, sa["entity_type"], item, counter)
        yield predict_single_output_id(**selections), chains


def load_msa_by_sequence(msas_table):
    """Map protein sequence -> MSA file path from the recycled msas table."""
    if not msas_table or not os.path.exists(msas_table):
        return {}
    df = pd.read_csv(msas_table)
    out = {}
    for _, row in df.iterrows():
        seq = str(row.get("sequence", "") or "")
        msa_file = row.get("msa_file", "") or ""
        if seq and msa_file and os.path.exists(str(msa_file)):
            out[seq] = str(msa_file)
    return out


def to_structure_input(chains, input_builder, MSA, max_depth):
    """Translate generic chain dicts into a StructurePredictionInput."""
    sequences = []
    for ch in chains:
        et = ch["entity_type"]
        if et == "protein":
            msa = None
            if ch.get("msa_path"):
                if not ch["msa_path"].endswith((".a3m", ".a3m.gz")):
                    raise ValueError(
                        f"ESMFold2 recycles MSAs in a3m format; got {ch['msa_path']}. "
                        "Convert a Boltz2 CSV msas output with MSA(source, convert='a3m').")
                msa = MSA.from_a3m(ch["msa_path"], max_sequences=max_depth)
            sequences.append(input_builder.ProteinInput(id=ch["id"], sequence=ch["sequence"], msa=msa))
        elif et in ("ssdna", "dsdna"):
            sequences.append(input_builder.DNAInput(id=ch["id"], sequence=ch["sequence"]))
        elif et in ("ssrna", "dsrna"):
            sequences.append(input_builder.RNAInput(id=ch["id"], sequence=ch["sequence"]))
        elif et == "ligand":
            if ch.get("ccd"):
                sequences.append(input_builder.LigandInput(id=ch["id"], ccd=[ch["ccd"]]))
            elif ch.get("smiles"):
                sequences.append(input_builder.LigandInput(id=ch["id"], smiles=ch["smiles"]))
            else:
                raise ValueError(f"ligand chain {ch['id']} has neither ccd nor smiles")
        else:
            raise ValueError(f"unknown entity_type: {et}")
    return input_builder.StructurePredictionInput(sequences=sequences)


def main():
    p = argparse.ArgumentParser(description="ESMFold2 complex inference")
    p.add_argument("--combinatorics-config", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--model-name", default="biohub/ESMFold2")
    p.add_argument("--num-loops", type=int, default=10)
    p.add_argument("--num-sampling-steps", type=int, default=100)
    p.add_argument("--num-diffusion-samples", type=int, default=1)
    p.add_argument("--num-seeds", type=int, default=1)
    p.add_argument("--msa-max-depth", type=int, default=1024)
    p.add_argument("--include-pae", action="store_true")
    p.add_argument("--all-samples", action="store_true")
    p.add_argument("--msas-table", default=None)
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    import torch
    from esm.utils.structure import input_builder
    from esm.utils.msa import MSA
    from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model
    from esm.models.esmfold2 import ESMFold2InputBuilder

    with open(args.combinatorics_config) as f:
        config = json.load(f)
    msa_by_seq = load_msa_by_sequence(args.msas_table)

    print(f"Loading ESMFold2 model: {args.model_name}", flush=True)
    model = ESMFold2Model.from_pretrained(args.model_name).cuda().eval()
    folder = ESMFold2InputBuilder()

    complexes = list(build_complexes(config, msa_by_seq))
    print(f"Assembled {len(complexes)} complex(es)", flush=True)

    scores = []
    failed = []
    for cid, chains in complexes:
        try:
            spi = to_structure_input(chains, input_builder, MSA, args.msa_max_depth)
            is_monomer = sum(1 for c in chains if c["entity_type"] != "ligand") <= 1

            # fold() returns a list of per-diffusion-sample results (a single
            # MolecularComplexResult when num_diffusion_samples==1); normalize to a list.
            samples = []  # (rank, metrics, result) across all seeds x diffusion samples
            for seed in range(args.num_seeds):
                with torch.no_grad():
                    out = folder.fold(
                        model, spi,
                        num_loops=args.num_loops,
                        num_sampling_steps=args.num_sampling_steps,
                        num_diffusion_samples=args.num_diffusion_samples,
                        seed=seed,
                    )
                results = out if isinstance(out, (list, tuple)) else [out]
                for res in results:
                    # res.plddt is a 0-1 fraction; framework convention (B-factors, Boltz2/ESMFold tables) is 0-100.
                    m = dict(plddt=float(res.plddt.mean()) * 100.0, ptm=float(res.ptm),
                             iptm=float(getattr(res, "interface_ptm", getattr(res, "iptm", float("nan")))))
                    if args.include_pae and getattr(res, "pae", None) is not None:
                        m["max_pae"] = float(res.pae.max())
                    rank = m["plddt"] if is_monomer else m["iptm"]
                    samples.append((rank, m, res))

            def write(out_id, metrics, res):
                out_path = os.path.join(args.output_dir, f"{out_id}.cif")
                with open(out_path, "w") as fh:
                    fh.write(res.complex.to_mmcif())
                scores.append(dict(id=out_id, **metrics))

            if args.all_samples:
                # surface every diffusion sample (best seed if multiple seeds) as <id>_k
                samples.sort(key=lambda t: t[0], reverse=True)
                keep = samples[:args.num_diffusion_samples]
                for k, (_r, m, res) in enumerate(keep, start=1):
                    write(f"{cid}_{k}", m, res)
                best_m = keep[0][1]
            else:
                _r, best_m, best_res = max(samples, key=lambda t: t[0])
                write(cid, best_m, best_res)

            print(f"Folded {cid} | pLDDT {best_m['plddt']:.1f} | pTM {best_m['ptm']:.3f} | "
                  f"ipTM {best_m['iptm']:.3f}", flush=True)
        except Exception as e:
            print(f"WARNING: {cid} failed: {e}", file=sys.stderr)
            traceback.print_exc()
            failed.append(cid)

    with open(os.path.join(args.output_dir, "ESMFold2_scores.json"), "w") as f:
        json.dump(scores, f)

    if failed:
        print(f"Failed {len(failed)}/{len(complexes)}: {failed}", file=sys.stderr)
    if not scores:
        print("ERROR: no complexes folded successfully", file=sys.stderr)
        sys.exit(1)
    print(f"ESMFold2 complete: {len(scores)} structure(s) written", flush=True)


if __name__ == "__main__":
    main()
