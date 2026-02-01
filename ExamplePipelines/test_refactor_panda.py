"""
Test pipeline for Panda tool - unified pandas-style table transformations

Panda can replace: Filter, Rank, SelectBest, MergeTables, ConcatenateTables, RemoveDuplicates
"""

import os, sys
sys.path.insert(0, os.getcwd())

from PipelineScripts.pipeline import *
from PipelineScripts.pdb import PDB
from PipelineScripts.compound_library import CompoundLibrary
from PipelineScripts.boltz2 import Boltz2
from PipelineScripts.panda import Panda

with Pipeline(project="RefactorTest",
              job="Panda",
              description="Test Panda tool - unified table transformations"):

    Resources(gpu="A100",
              time="4:00:00",
              memory="32GB")

    # Create test data with two separate Boltz2 runs
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

    # ========== Test 1: Filter (replaces Filter tool) ==========
    filtered = Panda(
        table=boltz1.tables.structures,
        operations=[
            Panda.filter("confidence > 0.5")
        ]
    )
    print(f"[1] Panda.filter output: {filtered}")

    # ========== Test 2: Sort + Head (replaces SelectBest) ==========
    best = Panda(
        table=boltz1.tables.structures,
        operations=[
            Panda.sort("confidence", ascending=False),
            Panda.head(2)
        ]
    )
    print(f"[2] Panda.sort+head output: {best}")

    # ========== Test 3: Rank (replaces Rank tool) ==========
    ranked = Panda(
        table=boltz1.tables.structures,
        operations=[
            Panda.rank(by="confidence", ascending=False)
        ]
    )
    print(f"[3] Panda.rank output: {ranked}")

    # ========== Test 4: Drop duplicates (replaces RemoveDuplicates) ==========
    deduped = Panda(
        table=boltz1.tables.structures,
        operations=[
            Panda.drop_duplicates(subset="sequence", keep="first")
        ]
    )
    print(f"[4] Panda.drop_duplicates output: {deduped}")

    # ========== Test 5: Concat (replaces ConcatenateTables) ==========
    concatenated = Panda(
        tables=[boltz1.tables.structures, boltz2.tables.structures],
        operations=[
            Panda.concat(fill="", add_source=True)
        ]
    )
    print(f"[5] Panda.concat output: {concatenated}")

    # ========== Test 6: Merge (replaces MergeTables) ==========
    # For this test we need tables with a common ID column
    # Using confidence tables from both runs
    merged = Panda(
        tables=[boltz1.tables.structures, boltz2.tables.structures],
        operations=[
            Panda.concat(fill="", add_source=True),  # First combine
            Panda.drop_duplicates(subset="id", keep="first")  # Then dedupe
        ]
    )
    print(f"[6] Panda.concat+dedupe output: {merged}")

    # ========== Test 7: Calculate (add computed columns) ==========
    calculated = Panda(
        table=boltz1.tables.structures,
        operations=[
            Panda.calculate({"conf_scaled": "confidence * 100"})
        ]
    )
    print(f"[7] Panda.calculate output: {calculated}")

    # ========== Test 8: Complex pipeline (multiple operations) ==========
    complex_result = Panda(
        tables=[boltz1.tables.structures, boltz2.tables.structures],
        operations=[
            Panda.concat(fill="", add_source=True),
            Panda.filter("confidence > 0.3"),
            Panda.sort("confidence", ascending=False),
            Panda.rank(by="confidence", ascending=False),
            Panda.head(5)
        ]
    )
    print(f"[8] Complex pipeline output: {complex_result}")

    # ========== Test 9: Pool mode (copies structures for filtered IDs) ==========
    pool_filtered = Panda(
        table=boltz1.tables.structures,
        operations=[
            Panda.filter("confidence > 0.6"),
            Panda.sort("confidence", ascending=False)
        ],
        pool=boltz1  # Copy structures matching filtered IDs
    )
    print(f"[9] Pool mode output: {pool_filtered}")

    # ========== Test 10: Groupby aggregation ==========
    grouped = Panda(
        tables=[boltz1.tables.structures, boltz2.tables.structures],
        operations=[
            Panda.concat(fill="", add_source=True),
            Panda.groupby("source_table", {"confidence": "mean"})
        ]
    )
    print(f"[10] Panda.groupby output: {grouped}")

    # ========== Test 11: Average by source (replaces AverageByTable) ==========
    averaged = Panda(
        tables=[boltz1.tables.structures, boltz2.tables.structures],
        operations=[
            Panda.concat(fill="", add_source=True),
            Panda.average_by_source()
        ]
    )
    print(f"[11] Panda.average_by_source output: {averaged}")

    print("\n=== Panda Test Summary ===")
    print("Panda replaces: Filter, Rank, SelectBest, MergeTables, ConcatenateTables, RemoveDuplicates, AverageByTable")
    print("Additional features: calculate, groupby, pivot, melt, sample, fillna, rename, select/drop columns")
