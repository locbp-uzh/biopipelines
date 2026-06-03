#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
PocketGen batch driver.

Replaces upstream `generate_new.py`'s `if __name__ == "__main__"` block. We
import the same `utils.*`/`models.*` symbols (running inside the `pocketgen`
env with the cloned repo on sys.path) but iterate `names = sorted(...)` from
the staging folder, honour an env-var checkpoint path, and skip Vina (it's a
post-hoc scoring step we don't need and pulls in qvina / pybel deps that
break Colab installs).

CLI:
    pipe_pocketgen_driver.py \\
        --pocketgen-repo /path/to/PocketGen \\
        --config /path/to/configs/train_model.yml \\
        --checkpoint /path/to/checkpoints/pocketgen.pt \\
        --target /path/to/staging \\
        --device cuda:0
"""

import argparse
import os
import sys
from functools import partial

import torch
from torch.utils.data import DataLoader
from torch_geometric.transforms import Compose
from tqdm import tqdm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pocketgen-repo", required=True, help="Cloned PocketGen repo root")
    ap.add_argument("--config", required=True)
    ap.add_argument("--checkpoint", required=True)
    ap.add_argument("--target", required=True, help="Staging folder; one subdir per pair")
    ap.add_argument("--device", default="cuda:0")
    args = ap.parse_args()

    if not os.path.isdir(args.pocketgen_repo):
        print(f"ERROR: PocketGen repo not found at {args.pocketgen_repo}", file=sys.stderr)
        sys.exit(1)

    # Run from the repo root so relative imports (utils.*, models.*) resolve
    # exactly as they do for upstream's generate_new.py.
    os.chdir(args.pocketgen_repo)
    sys.path.insert(0, args.pocketgen_repo)

    # Upstream bug: utils/relax.py uses `unit.kilojoules_per_mole` etc.
    # without importing `unit` as a bare name (it does
    # `from simtk.unit import kilocalories_per_mole, angstroms`, but the
    # body uses `unit.<...>`). Inject `unit` into the module's globals
    # before anything imports utils.relax. With openmm 8.x the unit module
    # lives at openmm.unit; simtk.unit is also present as a back-compat
    # shim, but openmm.unit is the canonical path.
    import utils.relax as _relax  # noqa: E402
    if not hasattr(_relax, "unit"):
        from openmm import unit as _openmm_unit
        _relax.unit = _openmm_unit

    # Upstream module imports (same set as generate_new.py, minus Vina).
    import esm  # noqa: E402
    from models.PD import Pocket_Design_new  # noqa: E402
    from utils.misc import load_config, seed_all  # noqa: E402
    from utils.transforms import FeaturizeProteinAtom, FeaturizeLigandAtom  # noqa: E402
    from utils.data import collate_mols_block  # noqa: E402
    from utils.protein_ligand import PDBProtein, parse_sdf_file  # noqa: E402
    from utils.data import torchify_dict  # noqa: E402

    config = load_config(args.config)
    seed_all(2089)

    protein_featurizer = FeaturizeProteinAtom()
    ligand_featurizer = FeaturizeLigandAtom()
    transform = Compose([protein_featurizer, ligand_featurizer])

    # ESM2 backbone
    pretrained_model, alphabet = esm.pretrained.load_model_and_alphabet_hub(
        "esm2_t33_650M_UR50D"
    )
    batch_converter = alphabet.get_batch_converter()
    ckpt = torch.load(args.checkpoint, map_location=args.device)
    del pretrained_model

    model = Pocket_Design_new(
        config.model,
        protein_atom_feature_dim=protein_featurizer.feature_dim,
        ligand_atom_feature_dim=ligand_featurizer.feature_dim,
        device=args.device,
    ).to(args.device)
    model.load_state_dict(ckpt["model"])

    # Discover pairs from the staging folder.
    names = sorted(
        d for d in os.listdir(args.target)
        if os.path.isdir(os.path.join(args.target, d))
    )
    if not names:
        print(f"ERROR: no pair folders under {args.target}", file=sys.stderr)
        sys.exit(1)

    # Reuse the helper functions exactly as upstream defines them in
    # generate_new.py. They're not exported as a module so we re-implement
    # the name2data piece inline (a copy of generate_new.py:name2data with
    # `transform` bound from this scope).
    def name2data(name: str):
        pdb_path = os.path.join(args.target, name, f"{name}.pdb")
        lig_path = os.path.join(args.target, name, f"{name}_ligand.sdf")
        pocket_path = os.path.join(args.target, name, f"{name}_pocket.pdb")

        with open(pdb_path) as f:
            pdb_block = f.read()
        protein = PDBProtein(pdb_block)
        seq = "".join(protein.to_dict_residue()["seq"])
        ligand = parse_sdf_file(lig_path, feat=False)
        r10_idx, r10_residues = protein.query_residues_ligand(
            ligand, radius=10, selected_residue=None, return_mask=False
        )
        full_seq_idx, _ = protein.query_residues_ligand(
            ligand, radius=3.5, selected_residue=r10_residues, return_mask=False
        )
        assert len(r10_idx) == len(r10_residues)

        with open(pocket_path, "w") as f:
            f.write(protein.residues_to_pdb_block(r10_residues))

        with open(pocket_path) as f:
            pocket = PDBProtein(f.read())

        pocket_dict = pocket.to_dict_atom()
        residue_dict = pocket.to_dict_residue()
        _, residue_dict["protein_edit_residue"] = pocket.query_residues_ligand(ligand)
        full_seq_idx.sort()
        r10_idx.sort()

        def _build(protein_dict, ligand_dict, residue_dict, seq, full_seq_idx, r10_idx):
            inst = {}
            for k, v in protein_dict.items():
                inst[f"protein_{k}"] = v
            for k, v in ligand_dict.items():
                inst[f"ligand_{k}"] = v
            for k, v in residue_dict.items():
                inst[k] = v
            inst["seq"] = seq
            inst["full_seq_idx"] = full_seq_idx
            inst["r10_idx"] = r10_idx
            return inst

        data = _build(
            torchify_dict(pocket_dict),
            torchify_dict(ligand),
            torchify_dict(residue_dict),
            seq,
            torch.tensor(full_seq_idx),
            torch.tensor(r10_idx),
        )
        data["protein_filename"] = pocket_path
        data["ligand_filename"] = lig_path
        data["whole_protein_name"] = pdb_path
        return transform(data)

    failed = []
    for name in tqdm(names, desc="PocketGen pairs"):
        try:
            data = name2data(name)
        except Exception as exc:
            print(f"  [fail] {name}: name2data: {type(exc).__name__}: {exc}", file=sys.stderr)
            failed.append(name)
            continue

        datalist = [data for _ in range(8)]
        dir_name = os.path.dirname(data["protein_filename"])
        os.makedirs(dir_name, exist_ok=True)

        model.generate_id = 0
        model.generate_id1 = 0
        loader = DataLoader(
            datalist,
            batch_size=4,
            shuffle=False,
            num_workers=getattr(config.train, "num_workers", 0),
            collate_fn=partial(collate_mols_block, batch_converter=batch_converter),
        )

        try:
            with torch.no_grad():
                model.eval()
                for batch in loader:
                    for key in batch:
                        if torch.is_tensor(batch[key]):
                            batch[key] = batch[key].to(args.device)
                    aar, rmsd, _ = model.generate(batch, dir_name)
                    print(f"  {name}: aar={aar} rmsd={rmsd}")
        except Exception as exc:
            print(f"  [fail] {name}: generate: {type(exc).__name__}: {exc}", file=sys.stderr)
            failed.append(name)
            continue

    if failed:
        print(f"WARNING: {len(failed)} pair(s) failed: {failed}", file=sys.stderr)

    # Exit non-zero only if every pair failed — partial success is still useful.
    if len(failed) == len(names):
        sys.exit(1)


if __name__ == "__main__":
    main()
