"""
Test pipeline for refactored tools - Basic tests
Tests: RFdiffusion, ProteinMPNN, MMseqs2, Filter, SelectBest
"""

import os, sys
sys.path.insert(0, os.getcwd())

from PipelineScripts.pipeline import *
from PipelineScripts.rfdiffusion import RFdiffusion
from PipelineScripts.protein_mpnn import ProteinMPNN
from PipelineScripts.mmseqs2 import MMseqs2
from PipelineScripts.filter import Filter
from PipelineScripts.select_best import SelectBest
from PipelineScripts.rank import Rank

with Pipeline(project="RefactorTest",
              job="Basic",
              description="Test basic refactored tools: RFdiffusion -> ProteinMPNN -> MMseqs2 -> Filter"):

    Resources(gpu="V100",
              time="2:00:00",
              memory="16GB")

    # Test 1: RFdiffusion de novo
    rfd = RFdiffusion(contigs="50-70",
                      num_designs=3,
                      steps=10)  # Low steps for quick test

    # Test 2: ProteinMPNN with DataStream input from RFdiffusion
    pmpnn = ProteinMPNN(structures=rfd,
                        num_sequences=2,
                        redesigned=rfd.tables.structures.designed)

    # Test 3: MMseqs2 for MSA generation
    msas = MMseqs2(sequences=pmpnn)

    # Test 4: Rank sequences by score
    ranked = Rank(data=pmpnn,
                  sort_by="score",
                  ascending=False)

    # Test 5: Filter by score threshold
    filtered = Filter(data=pmpnn,
                      expression="score > -2.0")

    # Test 6: SelectBest from filtered
    best = SelectBest(data=filtered,
                      n=2,
                      metric="score",
                      ascending=False)

    print("=== Output Summary ===")
    print(f"RFdiffusion structures: {len(rfd.structures.ids)}")
    print(f"ProteinMPNN sequences: {len(pmpnn.sequences.ids)}")
