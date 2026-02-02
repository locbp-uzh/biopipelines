"""
Test pipeline for refactored tools - Boltz2 tests
Tests: PDB, CompoundLibrary, Boltz2, MMseqs2
"""

import os, sys
sys.path.insert(0, os.getcwd())

from PipelineScripts.pipeline import *
from PipelineScripts.entities import *
from PipelineScripts.boltz2 import Boltz2
from PipelineScripts.mmseqs2 import MMseqs2

with Pipeline(project="RefactorTest",
              job="Boltz2",
              description="Test Boltz2 with various inputs"):

    Resources(gpu="A100",
              time="4:00:00",
              memory="8GB")

    # Test 1: Boltz2 with direct sequence
    boltz_seq = Boltz2(
        proteins=Sequence("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH")
    )

    # Test 2: PDB fetch and Boltz2 with structure input
    lysozyme = PDB("1AKI", ids="LYZ")

    boltz_pdb = Boltz2(
        proteins=lysozyme
    )

    # Test 3: Boltz2 with ligand (tests compound input)
    boltz_ligand = Boltz2(
        proteins=Sequence("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH"),
        ligands=Ligand("ethanol")
    )

    # Test 4: CompoundLibrary with multiple ligands
    compounds = CompoundLibrary({
        'ethanol': 'CCO',
        'methanol': 'CO',
        'propanol': 'CCCO'
    })

    boltz_multi = Boltz2(
        proteins=lysozyme,
        ligands=compounds
    )

    # Test 5: Boltz2 with MSAs

    boltz_msa = Boltz2(
        proteins=lysozyme,
        msas=boltz_pdb
    )

    print("=== Boltz2 Test Summary ===")
    print(f"boltz_seq structures: {len(boltz_seq.structures.ids)}")
    print(f"boltz_pdb structures: {len(boltz_pdb.structures.ids)}")
    print(f"boltz_ligand structures: {len(boltz_ligand.structures.ids)}")
    print(f"boltz_multi structures: {len(boltz_multi.structures.ids)}")
    print(f"boltz_msa structures: {len(boltz_msa.structures.ids)}")
