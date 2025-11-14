#!/usr/bin/env python3
"""
Test script for DNAEncoder tool.
"""

from PipelineScripts.pipeline import Pipeline
from PipelineScripts.dna_encoder import DNAEncoder
from PipelineScripts.load_output import LoadOutput

# Test with debug mode (doesn't actually submit to SLURM)
with Pipeline("TestProject", "DNAEncoderTest", "Testing DNA encoder tool", debug=True):

    # For testing, you would need a sequences output from a previous tool
    # Example: Load sequences from a previous LigandMPNN run
    # sequences = LoadOutput("path/to/ToolOutputs/XXX_LigandMPNN.json")

    # For now, just create the tool to test the structure
    # Uncomment below when you have actual sequences to test with

    # Test 1: E. coli optimization
    # dna_ec = DNAEncoder(
    #     sequences=sequences,
    #     organism="EC"
    # )

    # Test 2: Human optimization
    # dna_hs = DNAEncoder(
    #     sequences=sequences,
    #     organism="HS"
    # )

    # Test 3: Multi-organism optimization (E. coli & Human)
    # dna_mixed = DNAEncoder(
    #     sequences=sequences,
    #     organism="EC&HS"
    # )

    print("DNAEncoder tool loaded successfully!")
    print("To test with actual data, provide sequences from a previous tool output")

print("\nTest completed!")
