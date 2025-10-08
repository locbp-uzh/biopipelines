#!/usr/bin/env python3
"""
Example pipeline using DynamicBind for ligand-specific protein-ligand complex prediction.

This example demonstrates how to:
1. Use DynamicBind to predict protein-ligand binding conformations
2. Handle different input formats (PDB files, SMILES strings, tool outputs)
3. Generate predicted complex structures with affinity scores
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from PipelineScripts import Pipeline, DynamicBind, RFdiffusionAllAtom, CompoundLibrary

# Create pipeline
pipeline = Pipeline(name="DynamicBindExample")

# Example 1: Single protein from PDBs folder, SMILES string
dynamicbind1 = DynamicBind(
    proteins="example_protein.pdb",  # Input protein PDB file (from PDBs folder)
    ligands="CCO",                   # SMILES string (ethanol)
    samples_per_complex=10,          # Number of samples generated per complex
    inference_steps=20,              # Number of diffusion steps
    no_relax=False,                  # Relax final structures with OpenMM
    seed=42                          # Random seed for reproducibility
)

# Example 2: Multiple proteins, SMILES string
# dynamicbind2 = DynamicBind(
#     proteins=["protein1.pdb", "protein2.pdb"],
#     ligands="CC(=O)O",  # Acetic acid
#     samples_per_complex=10
# )

# Example 3: Using tool outputs (chaining with other tools)
# rfdaa_output = RFdiffusionAllAtom(ligand="ZIT", contigs="A1-100", num_designs=5)
# compound_lib = CompoundLibrary(smiles_list=["CCO", "CC(=O)O", "c1ccccc1"])
# dynamicbind3 = DynamicBind(
#     proteins=rfdaa_output,  # All structures from RFdiffusionAllAtom
#     ligands=compound_lib,   # All compounds from library
#     samples_per_complex=10
# )

# Add to pipeline
pipeline.add_tool(dynamicbind1)

# Run pipeline
if __name__ == "__main__":
    pipeline.run()

    # Access outputs
    output = dynamicbind1.get_output()
    print(f"\nDynamicBind completed!")
    print(f"Output folder: {output.output_folder}")
    print(f"Structures datasheet: {output.datasheets['structures'].path}")
    print(f"Compounds datasheet: {output.datasheets['compounds'].path}")
