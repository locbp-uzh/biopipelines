"""
Test pipeline for refactored tools - RFdiffusion AllAtom
Tests: RFdiffusionAllAtom, LigandMPNN with designed regions
"""

import os, sys
sys.path.insert(0, os.getcwd())

from PipelineScripts.pipeline import *
from PipelineScripts.pdb import PDB
from PipelineScripts.ligand import Ligand
from PipelineScripts.rfdiffusion_allatom import RFdiffusionAllAtom
from PipelineScripts.ligand_mpnn import LigandMPNN
from PipelineScripts.boltz2 import Boltz2
from PipelineScripts.mmseqs2 import MMseqs2

with Pipeline(project="RefactorTest",
              job="RFdiffusionAA",
              description="Test RFdiffusion AllAtom with ligand"):

    Resources(gpu="A100",
              time="6:00:00",
              memory="8GB")

    # Load HaloTag with ligand
    halotag_tmr = PDB("6U32", "HT")

    # Test RFdiffusionAllAtom
    rfdaa = RFdiffusionAllAtom(
        pdb=halotag_tmr,
        ligand=halotag_tmr,
        contigs="20-20,A5-300,20-20",
        num_designs=2
    )

    # Test LigandMPNN with designed regions from RFDAA
    lmpnn = LigandMPNN(
        structures=rfdaa,
        ligand="TMR",
        num_sequences=2,
        redesigned=rfdaa.tables.structures.designed
    )

    # Fold with Boltz2
    boltz = Boltz2(
        proteins=lmpnn,
        ligands=halotag_tmr
    )

    print("=== RFdiffusionAA Test Summary ===")
    print(f"RFDAA structures: {len(rfdaa.structures.ids)}")
    print(f"LigandMPNN sequences: {len(lmpnn.sequences.ids)}")
    print(f"Boltz2 structures: {len(boltz.structures.ids)}")
