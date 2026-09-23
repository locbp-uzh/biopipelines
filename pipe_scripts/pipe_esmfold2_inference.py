#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
ESMFold2 inference script.

Reads the combinatorics config written by the ESMFold2 tool, assembles one
complex per output id (proteins / nucleic acids / ligands as chains, plus any
modification / covalent-bond constraints from --constraints), folds it
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

# Let the caching allocator grow a segment instead of demanding one contiguous
# block. Folding many complexes of DIFFERENT sizes in one process fragments the
# pool badly (a 718-residue complex at 20 diffusion samples leaves holes that a
# 644-residue one cannot use), and PyTorch's own OOM message recommends exactly
# this setting for that pattern. Must be set before torch initialises CUDA,
# hence here rather than in the batch script. Complements the explicit
# empty_cache() in the per-complex loop; neither alone was sufficient.
os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "biopipelines"))
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
        # keep_default_na=False: a blank cell (e.g. the `ccd` column of a SMILES
        # ligand row) must read back as "", not NaN. Every consumer below tests
        # truthiness with `item.get(x) or ""`, and NaN is truthy in Python
        # (`bool(float("nan"))` is True), so the default read left a blank ccd
        # as NaN, which then satisfied `if ch.get("ccd"):` at the branch below
        # and got wrapped as `ccd=[nan]` -- crashing downstream in the ESM
        # library with "TypeError: sequence item 0: expected str instance,
        # float found" when it tried to `','.join(item.ccd)`. This surfaced
        # only after fixing the ccd/code column mix-up two commits back: `code`
        # always holds a real string ("LIG" by default) so it never hit the
        # NaN branch, but the correct `ccd` column is genuinely blank for every
        # SMILES ligand.
        records = pd.read_csv(path, keep_default_na=False, dtype={"id": str}).to_dict("records")
        if isinstance(source, dict) and source.get("ids"):
            from id_patterns import select_ids
            by_id = {str(r["id"]): r for r in records if "id" in r}
            records = [by_id[i] for i in select_ids([str(p) for p in source["ids"]], list(by_id))]
        if is_iter:
            if isinstance(source, dict) and source.get("group_by"):
                from pipe_axis_data import group_records
                records = group_records(records, source)
            iterated.extend(records)
            min_iter_order = min(min_iter_order, order)
        else:
            static.extend(records)
            min_static_order = min(min_static_order, order)
    return iterated, static, min_static_order < min_iter_order


def build_complexes(config, msa_by_seq, glycans=None):
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

    def add_glycan_chains(chains):
        """One ligand chain per glycan, each an ordered list of CCD codes.

        Glycans ride on the constraints payload rather than the compounds axis:
        they are not combinatorial, and Ligand cannot express a multi-residue
        chain anyway (_validate_ccd_code caps a code at 5 characters). Appended
        after the axis chains, then lettered with everything else, so they sort
        after the polymers like any other ligand.
        """
        for codes in (glycans or []):
            chains.append(dict(entity_type="ligand", id=None,
                               ccd=list(codes), smiles=""))
        return chains

    def assign_chain_ids(chains):
        """Letter the chains in entity order: polymers first, ligands last.

        Letters used to be handed out in call order, which followed whether an
        axis ITERATED, not what it contained. A bundled protein axis is static
        and a bare ligand axis iterates, so `ESMFold2(proteins=Bundle(a, b),
        ligands=lig)` gave the ligand chain A and pushed the proteins to B/C --
        contradicting this tool's own documented order and, worse, moving the
        protein letters whenever a ligand was added or removed. Any
        covalent_bonds/modifications constraint written against those letters
        then silently addressed a different chain.

        Ordering by entity type instead keeps polymer letters stable no matter
        how many ligands are attached, which is what makes multi-ligand work
        (glycans especially) tractable. The sort is stable, so chains within a
        group keep their assembly order and a double-stranded pair stays
        adjacent.
        """
        rank = {"protein": 0,
                "ssdna": 1, "dsdna": 1, "ssrna": 1, "dsrna": 1}
        chains.sort(key=lambda c: rank.get(c["entity_type"], 2))
        counter = [0]
        for c in chains:
            c["id"] = chain_letter(counter)
        return chains

    def add_item(chains, entity_type, item, counter):
        if "__members__" in item:
            # A Grouped element: its member rows are the chains of one complex, in row order.
            for member in item["__members__"]:
                add_item(chains, entity_type, member, counter)
            return
        if entity_type in LIGAND_TYPES:
            # `ccd` comes from the compounds table's `ccd` column, NOT its `code`
            # column. `ccd` is populated only for a genuine RCSB CCD lookup
            # (pipe_ligand.py:1058), while `code` is the residue label every
            # ligand carries and which defaults to "LIG" (pipe_ligand.py:337
            # leaves `ccd` empty for a SMILES ligand).
            #
            # Reading `code` here sent every SMILES ligand down the CCD branch
            # below, because that branch prefers ccd whenever it is non-empty.
            # LIG is itself a real CCD entry (C15H11N3, an indazole-pyridine),
            # so a SMILES-defined dye folded as that heterocycle with the
            # supplied SMILES sitting unused in the row. The protein folded
            # fine, so the failure is silent: confidence, ipTM and PAE are all
            # reported for a complex containing the wrong molecule.
            #
            # `format` is what actually distinguishes the two, so gate on it and
            # only fall back to `code` once the compound is known to be a CCD
            # entry -- `ccd` is blank for anything defined by SMILES.
            fmt = str(item.get("format") or "").strip().lower()
            smiles = str(item.get("smiles") or "").strip()
            if fmt == "ccd":
                ccd = str(item.get("ccd") or item.get("code") or "").strip()
            else:
                ccd = ""
            chains.append(dict(entity_type=entity_type, id=None,
                               ccd=[ccd] if ccd else [], smiles=smiles))
            return
        seq = str(item.get("sequence", ""))
        msa_path = msa_by_seq.get(seq) if entity_type == "protein" else None
        chains.append(dict(entity_type=entity_type, id=None,
                           sequence=seq, msa_path=msa_path))
        if entity_type in DOUBLE_STRANDED:
            chains.append(dict(entity_type=entity_type, id=None,
                               sequence=reverse_complement(seq, entity_type), msa_path=None))

    if not iteration_axes:
        chains, counter, selections = [], [0], {}
        for sa in static_only_axes:
            selections[sa["name"]] = ("bundle", [r["id"] for r in sa["items"]], None, [], False)
            for item in sa["items"]:
                add_item(chains, sa["entity_type"], item, counter)
        yield predict_single_output_id(**selections), assign_chain_ids(add_glycan_chains(chains))
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
        yield predict_single_output_id(**selections), assign_chain_ids(add_glycan_chains(chains))


def load_msa_by_sequence(msas_table):
    """Map protein sequence -> MSA file path from the recycled msas table."""
    if not msas_table or not os.path.exists(msas_table):
        return {}
    df = pd.read_csv(msas_table)
    out = {}
    for _, row in df.iterrows():
        seq = str(row.get("sequence", "") or "")
        msa_file = row.get("file", "") or ""
        if seq and msa_file and os.path.exists(str(msa_file)):
            out[seq] = str(msa_file)
    return out


def load_constraints(path):
    """Read the modification / covalent-bond config, or {} if absent."""
    if not path:
        return {}
    with open(path) as f:
        return json.load(f)


def group_modifications(constraints, chain_ids):
    """Group modification specs by chain, as (label, position, ccd) triples.

    Positions stay 1-indexed here; the conversion to upstream's 0-indexed
    Modification happens once, at the point of construction.
    """
    grouped = {}
    for i, mod in enumerate(constraints.get("modifications", [])):
        label = f"modifications[{i}]"
        chain = mod["chain"]
        if chain not in chain_ids:
            raise ValueError(
                f"{label}: chain {chain!r} is not in this complex (chains: {chain_ids})")
        grouped.setdefault(chain, []).append((label, int(mod["position"]), mod["ccd"]))
    return grouped


def build_modifications(input_builder, chain_id, sequence, grouped):
    """Turn this chain's specs into upstream Modification objects, or None."""
    specs = grouped.get(chain_id)
    if not specs:
        return None
    out = []
    for label, position, ccd in specs:
        if position > len(sequence):
            raise ValueError(
                f"{label}: position {position} is past the end of chain {chain_id} "
                f"({len(sequence)} residues)")
        out.append(input_builder.Modification(position=position - 1, ccd=ccd))
    return out


def build_covalent_bonds(input_builder, constraints, chain_ids):
    """Build the bonds with placeholder atom indices; resolve_covalent_atom_indices fills them."""
    bonds = []
    for i, bond in enumerate(constraints.get("covalent_bonds", [])):
        chain1, res1, _name1 = bond["atom1"]
        chain2, res2, _name2 = bond["atom2"]
        for key, chain in (("atom1", chain1), ("atom2", chain2)):
            if chain not in chain_ids:
                raise ValueError(
                    f"covalent_bonds[{i}].{key}: chain {chain!r} is not in this complex "
                    f"(chains: {chain_ids})")
        bonds.append(input_builder.CovalentBond(
            chain_id1=chain1, res_idx1=int(res1) - 1, atom_idx1=0,
            chain_id2=chain2, res_idx2=int(res2) - 1, atom_idx2=0))
    return bonds


def resolve_covalent_atom_indices(spi, specs):
    """Replace the placeholder atom indices with the ones upstream will index by.

    Upstream identifies a bonded atom by its position in the residue's atom list,
    and that list depends on whether the chain carries a covalent bond at all —
    leaving atoms are stripped only then. So the bonds must already be declared
    before the atoms can be enumerated, which is why they are built with
    placeholders and patched here rather than resolved up front.

    Upstream silently drops a bond whose atom index is out of range; every lookup
    failure below raises instead.
    """
    from esm.models.esmfold2.prepare_input import build_chains_from_input

    chains, tokens, atoms = build_chains_from_input(spi)
    asym_by_chain = {c.chain_id: c.asym_id for c in chains}

    residue_atoms = {}
    for atom in atoms:
        if not atom.is_valid or atom.token_index >= len(tokens):
            continue
        token = tokens[atom.token_index]
        residue_atoms.setdefault((token.asym_id, token.residue_index), []).append(atom)

    def atom_index(label, chain_id, res_idx, atom_name):
        found = residue_atoms.get((asym_by_chain[chain_id], res_idx))
        if not found:
            raise ValueError(
                f"{label}: chain {chain_id} has no residue at position {res_idx + 1}")
        names = [a.name for a in found]
        if atom_name not in names:
            residue_name = tokens[found[0].token_index].residue_name
            raise ValueError(
                f"{label}: atom {atom_name!r} is not in {residue_name} at "
                f"{chain_id}{res_idx + 1}; available: {names}")
        return names.index(atom_name)

    for i, (bond, spec) in enumerate(zip(spi.covalent_bonds, specs)):
        bond.atom_idx1 = atom_index(
            f"covalent_bonds[{i}].atom1", bond.chain_id1, bond.res_idx1, spec["atom1"][2])
        bond.atom_idx2 = atom_index(
            f"covalent_bonds[{i}].atom2", bond.chain_id2, bond.res_idx2, spec["atom2"][2])


def to_structure_input(chains, input_builder, MSA, max_depth, constraints):
    """Translate generic chain dicts into a StructurePredictionInput."""
    chain_ids = [ch["id"] for ch in chains]
    grouped_mods = group_modifications(constraints, chain_ids)

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
            sequences.append(input_builder.ProteinInput(
                id=ch["id"], sequence=ch["sequence"], msa=msa,
                modifications=build_modifications(
                    input_builder, ch["id"], ch["sequence"], grouped_mods)))
        elif et in ("ssdna", "dsdna"):
            sequences.append(input_builder.DNAInput(
                id=ch["id"], sequence=ch["sequence"],
                modifications=build_modifications(
                    input_builder, ch["id"], ch["sequence"], grouped_mods)))
        elif et in ("ssrna", "dsrna"):
            sequences.append(input_builder.RNAInput(
                id=ch["id"], sequence=ch["sequence"],
                modifications=build_modifications(
                    input_builder, ch["id"], ch["sequence"], grouped_mods)))
        elif et == "ligand":
            if ch["id"] in grouped_mods:
                raise ValueError(
                    f"{grouped_mods[ch['id']][0][0]}: chain {ch['id']} is a ligand; "
                    "modifications apply to polymer chains only")
            # ccd is only set when the compound's format is ccd (see add_item),
            # so these two are mutually exclusive and the order is not a policy.
            if ch.get("ccd"):
                print(f"ligand chain {ch['id']}: CCD {'-'.join(ch['ccd'])}", flush=True)
                sequences.append(input_builder.LigandInput(id=ch["id"], ccd=list(ch["ccd"])))
            elif ch.get("smiles"):
                print(f"ligand chain {ch['id']}: SMILES {ch['smiles']}", flush=True)
                sequences.append(input_builder.LigandInput(id=ch["id"], smiles=ch["smiles"]))
            else:
                raise ValueError(
                    f"ligand chain {ch['id']} has neither ccd nor smiles. A SMILES "
                    "compound needs a non-empty smiles column; a CCD compound needs "
                    "format=ccd and a code.")
        else:
            raise ValueError(f"unknown entity_type: {et}")

    bonds = build_covalent_bonds(input_builder, constraints, chain_ids)
    spi = input_builder.StructurePredictionInput(
        sequences=sequences,
        covalent_bonds=bonds or None,
    )
    if bonds:
        resolve_covalent_atom_indices(spi, constraints["covalent_bonds"])
    return spi


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
    p.add_argument("--constraints", default=None)
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
    constraints = load_constraints(args.constraints)

    print(f"Loading ESMFold2 model: {args.model_name}", flush=True)
    model = ESMFold2Model.from_pretrained(args.model_name).cuda().eval()
    folder = ESMFold2InputBuilder()

    complexes = list(build_complexes(config, msa_by_seq,
                                     glycans=constraints.get("glycans")))
    print(f"Assembled {len(complexes)} complex(es)", flush=True)

    scores = []
    failed = []
    for cid, chains in complexes:
        try:
            spi = to_structure_input(chains, input_builder, MSA, args.msa_max_depth, constraints)
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
                    # PER-CHAIN-PAIR ipTM, e.g. iptm_chain_0_2 for chain A vs C.
                    #
                    # The scalar `iptm` above is max-over-tokens of each token's
                    # mean inter-chain TM, pooled across ALL chain pairs, so in a
                    # protein+protein+ligand complex it is decided by the large
                    # protein-protein interface and carries essentially no
                    # information about the ligand (34 atoms cannot move a max
                    # taken over ~700 tokens). The model computes the full
                    # chain x chain matrix on every forward pass anyway; without
                    # this it was discarded, and recovering it needs a refold
                    # because the PAE logits it derives from are not saved.
                    #
                    # Boltz2 already reports the equivalent (ligand_iptm,
                    # pair_chains_iptm-i-j) and there it correlates +0.52 with
                    # dye pLDDT against +0.15 for the scalar -- i.e. the chain
                    # pair is the part that knows about the ligand.
                    #
                    # NOT comparable in absolute value to protein-protein ipTM:
                    # d0 is computed complex-wide, then averaged over only
                    # (n_ligand_atoms x n_chain_tokens) pairs, so ligand pairs
                    # sit systematically lower. Use it to rank designs, not as
                    # an absolute score.
                    pci = getattr(res, "pair_chains_iptm", None)
                    if pci is not None:
                        try:
                            n = int(pci.shape[0])
                            for a in range(n):
                                for b in range(a + 1, n):
                                    m[f"iptm_chain_{a}_{b}"] = float(pci[a][b])
                        except Exception:
                            pass
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
        finally:
            # RELEASE THE PER-DESIGN GPU TENSORS BEFORE THE NEXT COMPLEX.
            #
            # `samples` holds one result per diffusion sample, and each result
            # carries CUDA tensors (coordinates, pLDDT, PAE) plus the pair
            # representation they were decoded from. Rebinding it on the next
            # iteration drops the references, but the caching allocator keeps
            # the freed blocks — and they are large and badly shaped, so the
            # next complex fails on an allocation far smaller than the memory
            # nominally available.
            #
            # Invisible at 1-5 samples, which is why it survived this long.
            # Measured at 20 samples on an H200 (139.80 GiB): after a few
            # designs the process held 136.66 GiB, of which 17-26 GiB was
            # reserved-but-unallocated fragmentation, and a 10.93 GiB request
            # died. Ten of seventeen runs lost 34-85% of their designs this
            # way; the survivors were the ones whose binders happened to be
            # short enough to stay under the ceiling.
            #
            # Order-dependent by nature: a long binder early fragments the pool
            # for everything after it, so the same job can succeed or fail on
            # the same inputs in a different order.
            samples = None
            spi = None
            try:
                import torch as _torch
                if _torch.cuda.is_available():
                    _torch.cuda.empty_cache()
            except Exception:
                pass

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
