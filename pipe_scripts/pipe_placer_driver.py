#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
PLACER batch driver.

Runs in the `placer` env with the cloned baker-laboratory/PLACER repo on
sys.path. Mirrors what upstream `run_PLACER.py` does per input, but iterates
the BioPipelines structures stream and configures PLACERinput from the wrapper
inputs rather than parsing CLI flags.

Two modes (selected by --mode):

  ligand    — resolve the single `code` from the compounds stream and predict
              every copy of that ligand (predict_ligand([code]) + predict_multi).
              Output IDs are <structure>+<ligand>. PLACER errors (and the pair is
              recorded as failed) if the structure has no HETATM for `code`.

  sidechain — apo / sidechain repacking. No ligand is predicted; `--target-res`
              (chain-resno) is the crop center (target_res). Output IDs are just
              <structure>.

Common flow per input::

    placer = PLACER.PLACER()                 # loads default weights, once
    pl = PLACER.PLACERinput()
    pl.pdb(struct_path); pl.name(out_id)
    <mode-specific selectors>
    if exclude_sm: pl.exclude_sm(True)
    outputs = placer.run(pl, nsamples)
    PLACER.protocol.dump_output(outputs, f"{runs}/{out_id}", rerank=rerank)

`dump_output` writes the canonical artefacts upstream produces:
    {runs}/{out_id}_model.pdb   — all N models concatenated (MODEL/ENDMDL)
    {runs}/{out_id}.csv         — per-model scores

The post-processing step (biopipelines env) splits the multi-model PDB into
per-sample files and builds the BioPipelines streams.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import load_datastream, iterate_files  # noqa: E402
from biopipelines.ligand_utils import resolve_ligand_code  # noqa: E402
from biopipelines.pdb_parser import field_res_name, field_res_seq  # noqa: E402


def _count_ligand_copies(struct_path, code):
    """Count distinct copies of HETATM residue `code` (by chain+resseq) in a
    PDB/mmCIF, to decide whether PLACER's predict_multi path applies. mmCIF is
    handled by gemmi when available; PDB by fixed-column parsing."""
    code = code.strip()
    ext = os.path.splitext(struct_path)[1].lower()
    seen = set()
    if ext in (".cif", ".cif.gz"):
        try:
            import gemmi
            st = gemmi.read_structure(struct_path)
            for model in st:
                for chain in model:
                    for res in chain:
                        if res.name.strip() == code and res.het_flag == "H":
                            seen.add((chain.name, str(res.seqid)))
            return len(seen)
        except Exception:
            return 1  # be permissive: fall back to single-ligand path
    with open(struct_path) as fh:
        for ln in fh:
            if ln.startswith("HETATM") and field_res_name(ln) == code:
                seen.add((ln[20:22].strip(), ln[22:27].strip()))
    return len(seen)


def _atom_sel_to_placer_tuple(sel, struct_path):
    """Translate a '<residue>.<atom>' selection into PLACER's (chain, resno,
    name3, atom) tuple, resolving the residue name (name3) from the structure.

    The residue part may be chain-prefixed ("A145"), a bare number ("145"), or a
    residue name ("LIG"). Reads ATOM/HETATM lines from a PDB to find name3 and,
    when only a name was given, the residue number/chain.
    """
    res_part, atom_name = sel.rsplit(".", 1)
    chain = resnum = resname = None
    i = 0
    while i < len(res_part) and not (res_part[i].isdigit() or res_part[i] == "-"):
        i += 1
    if i > 0 and i < len(res_part):
        chain = res_part[:i]; res_part = res_part[i:]
    if res_part.lstrip("-").isdigit():
        resnum = int(res_part)
    else:
        resname = res_part
    # Scan the structure to fill in the missing fields (name3 always; resnum/chain
    # if a residue name was given).
    with open(struct_path) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            ln_chain = line[21] if len(line) > 21 else " "
            ln_resname = field_res_name(line)
            try:
                ln_resnum = int(field_res_seq(line))
            except ValueError:
                continue
            if chain is not None and ln_chain != chain:
                continue
            if resnum is not None and ln_resnum != resnum:
                continue
            if resname is not None and ln_resname != resname:
                continue
            # First matching residue wins.
            return (ln_chain, ln_resnum, ln_resname, atom_name.strip())
    raise ValueError(f"bond atom selection {sel!r} matched no residue in {struct_path}")


def parse_target_res(s):
    """Parse a "chain-resno" or "chain-resno-name3" crop-center string into the
    tuple PLACERinput.target_res expects, matching upstream run_PLACER.py."""
    parts = s.split("-")
    if len(parts) == 2:
        return (parts[0], int(parts[1]))
    if len(parts) == 3:
        return (parts[0], int(parts[1]), parts[2])
    raise ValueError(
        f"target_res must be 'chain-resno' or 'chain-resno-name3', got {s!r}"
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--placer-repo", required=True, help="Cloned PLACER repo root")
    ap.add_argument("--structures-json", required=True)
    ap.add_argument("--runs-folder", required=True, help="One {out_id}_model.pdb + {out_id}.csv per input")
    ap.add_argument("--nsamples", type=int, required=True)
    ap.add_argument("--rerank", default=None, choices=["prmsd", "plddt", "plddt_pde"])
    ap.add_argument("--exclude-sm", action="store_true", help="Drop all small molecules (true apo)")
    ap.add_argument("--mode", required=True, choices=["ligand", "sidechain"])
    ap.add_argument("--ligand-json", default=None, help="compounds stream JSON (ligand mode)")
    ap.add_argument("--target-res", default=None, help="chain-resno crop center (sidechain mode)")
    ap.add_argument("--bonds-json", default=None,
                    help="JSON list of [atom1, atom2, length] covalent bonds to enforce; "
                         "atoms in '<residue>.<atom>' syntax (e.g. 'A145.SG').")
    args = ap.parse_args()

    if args.mode == "ligand" and not args.ligand_json:
        print("ERROR: --ligand-json is required in ligand mode", file=sys.stderr)
        sys.exit(1)
    if args.mode == "sidechain" and not args.target_res:
        print("ERROR: --target-res is required in sidechain mode", file=sys.stderr)
        sys.exit(1)

    if not os.path.isfile(os.path.join(args.placer_repo, "PLACER.py")):
        print(f"ERROR: PLACER repo not found at {args.placer_repo}", file=sys.stderr)
        sys.exit(1)

    # Import the upstream package from the clone (its PLACER.py sits at the
    # repo root and re-exports PLACERinput / protocol).
    sys.path.insert(0, args.placer_repo)
    import PLACER  # noqa: E402

    structures = list(iterate_files(load_datastream(args.structures_json)))
    if not structures:
        print("ERROR: no input structures", file=sys.stderr)
        sys.exit(1)

    # In ligand mode, each structure pairs with each ligand id; in sidechain
    # mode there is no ligand axis (one "unit" per structure).
    if args.mode == "ligand":
        ligand_code = resolve_ligand_code(args.ligand_json)
        print(f"Ligand code (predict_ligand name3): {ligand_code}")
        ligand_ids = list(load_datastream(args.ligand_json).ids_expanded)
        if not ligand_ids:
            print("ERROR: ligand stream is empty", file=sys.stderr)
            sys.exit(1)
    else:
        ligand_code = None
        ligand_ids = [None]  # single pass per structure
        print(f"Sidechain/apo mode; target_res = {args.target_res}")

    # Optional covalent bonds to enforce (list of [atom1, atom2, length]).
    bond_specs = []
    if args.bonds_json:
        with open(args.bonds_json) as bf:
            bond_specs = json.load(bf)
        print(f"Covalent bonds to enforce: {bond_specs}")

    os.makedirs(args.runs_folder, exist_ok=True)

    # Load the model once; reuse across inputs.
    placer = PLACER.PLACER()

    failed = []
    n_ok = 0
    for struct_id, struct_path in structures:
        ext = os.path.splitext(struct_path)[1].lower()
        for lig_id in ligand_ids:
            out_id = f"{struct_id}+{lig_id}" if lig_id is not None else struct_id
            try:
                pl = PLACER.PLACERinput()
                if ext in (".cif", ".cif.gz"):
                    pl.cif(struct_path)
                else:
                    pl.pdb(struct_path)
                pl.name(out_id)

                if args.mode == "ligand":
                    # name3 selector: predict every copy of `ligand_code`.
                    pl.predict_ligand([ligand_code])
                    # predict_multi(True) scores all matched copies, but upstream
                    # PLACER asserts N_ligands > 1 inside it — so it must only be
                    # enabled when the structure actually holds >1 copy of the code.
                    # A single docked ligand (the common case here) uses the plain
                    # single-ligand path.
                    if _count_ligand_copies(struct_path, ligand_code) > 1:
                        pl.predict_multi(True)
                    # Enforce any covalent bonds (resolve name3/resnum per structure).
                    if bond_specs:
                        resolved = []
                        for a1, a2, blen in bond_specs:
                            t1 = _atom_sel_to_placer_tuple(a1, struct_path)
                            t2 = _atom_sel_to_placer_tuple(a2, struct_path)
                            resolved.append([t1, t2, float(blen)])
                        pl.bonds(resolved)
                        print(f"  {out_id}: enforcing {len(resolved)} covalent bond(s): {resolved}")
                else:
                    # Crop center on a protein residue; no ligand predicted.
                    # PLACERinput.target_res wants a tuple; parse "chain-resno"
                    # or "chain-resno-name3" exactly as upstream run_PLACER.py.
                    pl.target_res(parse_target_res(args.target_res))

                if args.exclude_sm:
                    pl.exclude_sm(True)

                outputs = placer.run(pl, args.nsamples)
                prefix = os.path.join(args.runs_folder, out_id)
                PLACER.protocol.dump_output(outputs, prefix, rerank=args.rerank)
                print(f"  [ok] {out_id}: {len(outputs)} model(s) -> {prefix}_model.pdb")
                n_ok += 1
            except Exception as exc:
                print(f"  [fail] {out_id}: {type(exc).__name__}: {exc}", file=sys.stderr)
                failed.append((out_id, f"{type(exc).__name__}: {exc}"))

    print(f"PLACER: {n_ok} input(s) succeeded, {len(failed)} failed")
    if failed:
        for out_id, msg in failed:
            print(f"  - {out_id}: {msg}", file=sys.stderr)

    # Exit non-zero only if every input failed — partial success is still useful.
    if n_ok == 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
