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
              memory="32GB")

    # Load HaloTag with ligand
    halotag = PDB("6U32", "HT")

    # Get ligand from the structure
    lig = Ligand("TMR", source=halotag)

    # Test RFdiffusionAllAtom
    rfdaa = RFdiffusionAllAtom(
        pdb=halotag,
        ligand=lig,
        contigs="A1-100,0 B1-50",
        num_designs=2,
        steps=10
    )

    # Test LigandMPNN with designed regions from RFDAA
    lmpnn = LigandMPNN(
        structures=rfdaa,
        ligand="LIG",
        num_sequences=2,
        redesigned=rfdaa.tables.structures.designed
    )

    # Get MSAs
    msas = MMseqs2(sequences=lmpnn)

    # Fold with Boltz2
    boltz = Boltz2(
        proteins=lmpnn,
        ligands=lig,
        msas=msas,
        num_seeds=1,
        num_recycles=1
    )

    print("=== RFdiffusionAA Test Summary ===")
    print(f"RFDAA structures: {len(rfdaa.structures.ids)}")
    print(f"LigandMPNN sequences: {len(lmpnn.sequences.ids)}")
    print(f"Boltz2 structures: {len(boltz.structures.ids)}")
