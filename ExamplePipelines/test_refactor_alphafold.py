"""
Test pipeline for refactored tools - AlphaFold and ESMFold
Tests: AlphaFold, ESMFold with various inputs
"""

import os, sys
sys.path.insert(0, os.getcwd())

from PipelineScripts.pipeline import *
from PipelineScripts.pdb import PDB
from PipelineScripts.rfdiffusion import RFdiffusion
from PipelineScripts.protein_mpnn import ProteinMPNN
from PipelineScripts.mmseqs2 import MMseqs2
from PipelineScripts.alphafold import AlphaFold
from PipelineScripts.esmfold import ESMFold

with Pipeline(project="RefactorTest",
              job="Folding",
              description="Test AlphaFold and ESMFold"):

    Resources(gpu="any",
              time="6:00:00",
              memory="16GB")

    # Generate some sequences to fold
    rfd = RFdiffusion(
        contigs="40-60",
        num_designs=2
    )

    pmpnn = ProteinMPNN(
        structures=rfd,
        num_sequences=2,
        redesigned=rfd.tables.structures.designed
    )

    # Test AlphaFold with sequences and MSAs
    af = AlphaFold(
        sequences=pmpnn
    )

    print("=== Folding Test Summary ===")
    print(f"AlphaFold structures: {len(af.structures.ids)}")
