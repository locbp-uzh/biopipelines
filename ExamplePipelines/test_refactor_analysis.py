"""
Test pipeline for refactored tools - Analysis tools
Tests: MergeTables, ConcatenateTables, PLIP, DistanceSelector
"""

import os, sys
sys.path.insert(0, os.getcwd())

from PipelineScripts.pipeline import *
from PipelineScripts.pdb import PDB
from PipelineScripts.compound_library import CompoundLibrary
from PipelineScripts.boltz2 import Boltz2
from PipelineScripts.plip import PLIP
from PipelineScripts.merge_tables import MergeTables
from PipelineScripts.concatenate_tables import ConcatenateTables

with Pipeline(project="RefactorTest",
              job="Analysis",
              description="Test analysis tools: PLIP, MergeTables, ConcatenateTables"):

    Resources(gpu="A100",
              time="4:00:00",
              memory="32GB")

    # Create two separate Boltz2 runs to test table operations
    lysozyme = PDB("1AKI", "LYZ")

    compounds1 = CompoundLibrary({
        'ethanol': 'CCO',
        'methanol': 'CO'
    })

    compounds2 = CompoundLibrary({
        'propanol': 'CCCO',
        'butanol': 'CCCCO'
    })

    boltz1 = Boltz2(
        proteins=lysozyme,
        ligands=compounds1,
        num_seeds=1,
        num_recycles=1
    )

    boltz2 = Boltz2(
        proteins=lysozyme,
        ligands=compounds2,
        num_seeds=1,
        num_recycles=1
    )

    # Test PLIP for protein-ligand interactions
    plip1 = PLIP(structures=boltz1, ligand="LIG")
    plip2 = PLIP(structures=boltz2, ligand="LIG")

    # Test ConcatenateTables - stack rows from multiple sources
    concat = ConcatenateTables(
        tables=[plip1.tables.contacts, plip2.tables.contacts]
    )

    # Test MergeTables - join tables by id column
    merged = MergeTables(
        tables=[boltz1.tables.structures, plip1.tables.summary],
        on="id"
    )

    print("=== Analysis Test Summary ===")
    print(f"PLIP1 contacts table: {plip1.tables.contacts}")
    print(f"PLIP2 contacts table: {plip2.tables.contacts}")
    print(f"Concatenated table: {concat.tables.concatenated}")
    print(f"Merged table: {merged.tables.merged}")
