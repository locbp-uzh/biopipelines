"""
Test pipeline for refactored tools - Boltz2 tests
Tests: PDB, CompoundLibrary, Boltz2, MMseqs2
"""

import os, sys
sys.path.insert(0, os.getcwd())

from PipelineScripts.pipeline import *
from PipelineScripts.pdb import PDB
from PipelineScripts.compound_library import CompoundLibrary
from PipelineScripts.boltz2 import Boltz2
from PipelineScripts.mmseqs2 import MMseqs2

with Pipeline(project="RefactorTest",
              job="Boltz2",
              description="Test Boltz2 with various inputs"):

    Resources(gpu="A100",
              time="4:00:00",
              memory="32GB")

    # Test 1: Boltz2 with direct sequence
    boltz_seq = Boltz2(
        proteins="MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH",
        num_seeds=1,
        num_recycles=1  # Low for quick test
    )

    # Test 2: PDB fetch and Boltz2 with structure input
    lysozyme = PDB("1AKI", "LYZ")

    boltz_pdb = Boltz2(
        proteins=lysozyme,
        num_seeds=1,
        num_recycles=1
    )

    # Test 3: Boltz2 with ligand (tests compound input)
    boltz_ligand = Boltz2(
        proteins="MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH",
        ligands="CCO",  # Ethanol - simple ligand
        affinity=True,
        num_seeds=1,
        num_recycles=1
    )

    # Test 4: CompoundLibrary with multiple ligands
    compounds = CompoundLibrary({
        'ethanol': 'CCO',
        'methanol': 'CO',
        'propanol': 'CCCO'
    })

    boltz_multi = Boltz2(
        proteins=lysozyme,
        ligands=compounds,
        num_seeds=1,
        num_recycles=1
    )

    # Test 5: Boltz2 with MSAs
    msas = MMseqs2(sequences=lysozyme)

    boltz_msa = Boltz2(
        proteins=lysozyme,
        msas=msas,
        num_seeds=1,
        num_recycles=1
    )

    print("=== Boltz2 Test Summary ===")
    print(f"boltz_seq structures: {len(boltz_seq.structures.ids)}")
    print(f"boltz_pdb structures: {len(boltz_pdb.structures.ids)}")
    print(f"boltz_ligand structures: {len(boltz_ligand.structures.ids)}")
    print(f"boltz_multi structures: {len(boltz_multi.structures.ids)}")
    print(f"boltz_msa structures: {len(boltz_msa.structures.ids)}")
