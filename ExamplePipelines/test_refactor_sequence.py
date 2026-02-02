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

    Resources(gpu="any",
              time="4:00:00",
              memory="16GB")

    # Load a PDB with ligand for LigandMPNN test
    halotag = PDB("6U32", ids= "HT")

    # Test LigandMPNN
    lmpnn = LigandMPNN(
        structures=halotag,
        ligand="TMR",  # Ligand code in 6U32
        num_sequences=3
    )

    # Test SDM (Site-Directed Mutagenesis)
    sdm = SDM(
        original=halotag,
        position=14,
        mode="charged"
    )

    # Test Fuse - combine chain sequences
    fused = Fuse(
        proteins=[halotag,halotag],
        name="HTdimer",
        linker_lengths=["4-5"]
)
