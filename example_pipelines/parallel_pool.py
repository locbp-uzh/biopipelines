# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP
#
# Licensed under the MIT License. See LICENSE file in the project root for details.
"""
Demo: ``with Parallel():`` + ``Pool(runs=runs)``.

CPU-only and using only basic entity types so the pipeline runs end-to-end
on the cluster without any GPU. Three independent demos, each pooling a
different entity type, all wrapped in a single Pipeline:

  1. **Sequence pool** — three Sequence runs in parallel; pool their
     ``sequences`` stream; sort the result.
  2. **Ligand pool** — three Ligand runs in parallel; pool their
     ``compounds`` stream.
  3. **PDB pool** — three PDB downloads in parallel; pool their
     ``structures`` + ``sequences`` + ``compounds`` streams (PDB emits
     all three).

After the run, inspect the gathered output of each demo: every map_table
should carry ids of the form ``<orig_id>_<pool_idx>`` and a ``pool.path``
column with values 1, 2, 3.
"""

from biopipelines.pipeline import *
from biopipelines import Sequence, Ligand, PDB, Pool, Panda


with Pipeline(project="Examples",
              job="ParallelPoolDemo",
              description="Parallel + Pool fan-out / fan-in demo (CPU-only)"):

    # ────────────────────────────────────────────────────────────────────────
    # Demo 1 — Sequence pool
    # ────────────────────────────────────────────────────────────────────────
    Suffix("sequences")
    seq_runs = []
    with Parallel():
        for label, payload in [
            ("alpha", "MKTAYIAKQRQISFVKSHFSRQLE"),
            ("beta",  "AETGFRYIQLDIDDFSEKQADYLG"),
            ("gamma", "GGGGAEPCKKLVVAGRHCAYAGC"),
        ]:
            Resources(time="00:30:00", memory="4GB")
            seq_runs.append(Sequence(seq=payload, ids=label, type="protein"))

    Resources(time="00:30:00", memory="4GB")
    seq_pool = Pool(runs=seq_runs)
    Panda(tables=seq_pool.tables.sequences,
          operations=[Panda.sort("id", ascending=True)])

    # ────────────────────────────────────────────────────────────────────────
    # Demo 2 — Ligand pool (compounds stream)
    # ────────────────────────────────────────────────────────────────────────
    Suffix("ligands")
    lig_runs = []
    with Parallel():
        for code in ["CCO", "CCC", "CC(=O)O"]:
            Resources(time="00:30:00", memory="4GB")
            lig_runs.append(Ligand(smiles=code, ids="lig"))

    Resources(time="00:30:00", memory="4GB")
    lig_pool = Pool(runs=lig_runs)
    Panda(tables=lig_pool.tables.compounds,
          operations=[Panda.sort("id", ascending=True)])

    # ────────────────────────────────────────────────────────────────────────
    # Demo 3 — PDB pool (structures + sequences + compounds streams)
    # ────────────────────────────────────────────────────────────────────────
    Suffix("pdbs")
    pdb_runs = []
    with Parallel():
        for code in ["1ubq", "1crn", "1l2y"]:
            # Same custom id "p" across runs → after Pool, ids become
            # p_1, p_2, p_3 — the canonical "same upstream tool, same
            # parameters, ids collide" case the framework is designed for.
            Resources(time="00:30:00", memory="4GB")
            pdb_runs.append(PDB(code, ids="p"))

    Resources(time="00:30:00", memory="4GB")
    pdb_pool = Pool(runs=pdb_runs)
    Panda(tables=pdb_pool.tables.structures,
          operations=[Panda.sort("id", ascending=True)])
