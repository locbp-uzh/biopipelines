"""
Test pipeline for refactored tools - Sequence design tools
Tests: LigandMPNN, MutationComposer, SDM, Fuse, StitchSequences
"""

import os, sys
sys.path.insert(0, os.getcwd())

from PipelineScripts.pipeline import *
from PipelineScripts.pdb import PDB
from PipelineScripts.rfdiffusion_allatom import RFdiffusionAllAtom
from PipelineScripts.ligand_mpnn import LigandMPNN
from PipelineScripts.mutation_composer import MutationComposer
from PipelineScripts.site_directed_mutagenesis import SDM
from PipelineScripts.fuse import Fuse
from PipelineScripts.stitch_sequences import StitchSequences

with Pipeline(project="RefactorTest",
              job="SequenceDesign",
              description="Test sequence design tools"):

    Resources(gpu="V100",
              time="4:00:00",
              memory="16GB")

    # Load a PDB with ligand for LigandMPNN test
    halotag = PDB("6U32", "HT")

    # Test LigandMPNN
    lmpnn = LigandMPNN(
        structures=halotag,
        ligand="TMR",  # Ligand code in 6U32
        num_sequences=3
    )

    # Test MutationComposer - generate mutation combinations
    mutations = MutationComposer(
        sequences=halotag,
        positions="A10,A20,A30",
        mutations="ALA,VAL"
    )

    # Test SDM (Site-Directed Mutagenesis)
    sdm = SDM(
        sequences=halotag,
        mutations="A10V,A20L"
    )

    # Test Fuse - combine chain sequences
    fused = Fuse(
        sequences=halotag,
        chains="A",
        linker="GSGSGS"
    )

    # Test StitchSequences - combine sequences from different sources
    stitched = StitchSequences(
        sequences=[lmpnn, sdm],
        mode="concatenate"
    )

    print("=== Sequence Design Test Summary ===")
    print(f"LigandMPNN sequences: {len(lmpnn.sequences.ids)}")
    print(f"MutationComposer sequences: {len(mutations.sequences.ids)}")
    print(f"SDM sequences: {len(sdm.sequences.ids)}")
