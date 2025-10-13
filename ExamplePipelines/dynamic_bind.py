#!/usr/bin/env python3
"""
Example pipeline using DynamicBind for ligand-specific protein-ligand complex prediction.

This example demonstrates how to:
1. Use DynamicBind to predict protein-ligand binding conformations
2. Handle different input formats (PDB files, SMILES strings, tool outputs)
3. Generate predicted complex structures with affinity scores
"""

import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import Pipeline
from PipelineScripts.dynamic_bind import DynamicBind
from PipelineScripts.rfdiffusion_allatom import RFdiffusionAllAtom
from PipelineScripts.ligand_mpnn import LigandMPNN
from PipelineScripts.compound_library import CompoundLibrary
from PipelineScripts.pdb import PDB

# Create pipeline
pipeline = Pipeline("Examples","DynamicBind","Test of several input formats for DynamicBind")

pipeline.resources(gpu="32GB", memory="16GB", time = "1:00:00")

HaloTag_TMR = pipeline.add(PDB("6U32",)) #HT7 bound to rhodamine

# Example 1: Single protein from PDB, SMILES string
dynamicbind = pipeline.add(DynamicBind(
    proteins=HaloTag_TMR,            # Input protein PDB file (from PDBs folder)
    ligands="ClCCCCCCOCCOCCNC(=O)C", # SMILES string (halotag linker acetamide)
))
print(dynamicbind)

# Example 2: Single protein from PDBs folder, SMILES string
dynamicbind = pipeline.add(DynamicBind(
    proteins=HaloTag_TMR,# HaloTag
    ligands=HaloTag_TMR, # Tetramethyl RHodamine
))
print(dynamicbind)

# Example 3: Multiple proteins, SMILES string
dehalogenases = pipeline.add(PDB(pdbs=["1BN6","4RAS"]))
dynamicbind = pipeline.add(DynamicBind(
     proteins=dehalogenases,
     ligands="ClCCCCCCOCCOCCNC(=O)C",  # Acetic acid
     samples_per_complex=10
))
print(dynamicbind)

# Example 4: Using tool outputs (chaining with other tools)
rfdaa = pipeline.add(RFdiffusionAllAtom(pdb=HaloTag_TMR, ligand="PVY", contigs="A6-140,100-140", num_designs=5))
lmpnn = pipeline.add(LigandMPNN(structures=rfdaa,ligand="PVY",redesigned=rfdaa.datasheets.structures.designed))
dynamicbind = DynamicBind(
    proteins=rfdaa,  # All structures from RFdiffusionAllAtom
    ligands="ClCCCCCCOCCOCCNC(=O)C"
)
print(dynamicbind)

# Example 5: One protein, multiple ligands
compound_library = pipeline.add(CompoundLibrary({"acetamide":"ClCCCCCCOCCOCCNC(=O)C",
                                                 "sulfonamide":"ClCCCCCCOCCOCCNS(=O)(=O)C"}))
dynamicbind=pipeline.add(DynamicBind(proteins=HaloTag_TMR,
                                      ligands=compound_library))
print(dynamicbind)

# Example 6: Multiple proteins, multiple ligands
compound_library = pipeline.add(CompoundLibrary({"acetamide":"ClCCCCCCOCCOCCNC(=O)C",
                                                 "sulfonamide":"ClCCCCCCOCCOCCNS(=O)(=O)C"}))
dynamicbind=pipeline.add(DynamicBind(proteins=rfdaa,
                                      ligands=compound_library))
print(dynamicbind)

pipeline.slurm()