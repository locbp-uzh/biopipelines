"""
Quick smoke test for refactored tools - minimal resource usage
Tests core DataStream I/O flow through multiple tools
"""

import os, sys
sys.path.insert(0, os.getcwd())

from PipelineScripts.pipeline import *
from PipelineScripts.pdb import PDB
from PipelineScripts.rfdiffusion import RFdiffusion
from PipelineScripts.protein_mpnn import ProteinMPNN
from PipelineScripts.boltz2 import Boltz2
from PipelineScripts.filter import Filter
from PipelineScripts.select_best import SelectBest

with Pipeline(project="RefactorTest",
              job="QuickSmoke",
              description="Quick smoke test of refactored DataStream I/O"):

    Resources(gpu="V100",
              time="1:00:00",
              memory="16GB")

    # 1. RFdiffusion: de novo backbone (minimal)
    rfd = RFdiffusion(
        contigs="30-40",
        num_designs=2,
        steps=5  # Minimal steps
    )
    print(f"[1] RFdiffusion output: {rfd}")

    # 2. ProteinMPNN: sequence design
    pmpnn = ProteinMPNN(
        structures=rfd,
        num_sequences=2,
        redesigned=rfd.tables.structures.designed
    )
    print(f"[2] ProteinMPNN output: {pmpnn}")

    # 3. Boltz2: structure prediction with ligand (skips MMseqs2 for quick test)
    boltz = Boltz2(
        proteins=pmpnn,
        ligands="CCO",  # Ethanol - simple ligand to trigger affinity mode
        num_seeds=1,
        num_recycles=1
    )
    print(f"[3] Boltz2 output: {boltz}")

    # 4. Filter: expression-based filtering
    filtered = Filter(
        data=boltz,
        expression="confidence > 0.5"
    )
    print(f"[4] Filter output: {filtered}")

    # 5. SelectBest: top selection
    best = SelectBest(
        pool=boltz,
        tables=boltz.tables.structures,
        metric="confidence",
        mode="max"
    )
    print(f"[5] SelectBest output: {best}")

    print("\n=== Quick Smoke Test Complete ===")
    print("All DataStream I/O connections verified at pipeline time.")
