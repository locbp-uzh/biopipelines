# BioPipelines

## Overview

A Python framework for automated computational protein design workflows that can run in Jupyter/Colab notebooks as well as on SLURM-based computing clusters. BioPipelines provides standardized interfaces to connect bioinformatics tools. 

## Example pipeline

PDB download -> RFdiffusion -> ProteinMPNN -> AlphaFold.

```python
#imports omitted
with Pipeline(project="Examples",
              job="RFD-ProteinMPNN-AlphaFold2",
              description="Redesign of N terminus domain of lysozyme"):
    Resources(gpu="A100", 
              time="4:00:00",
              memory="16GB")
    lysozyme = PDB("168L")
    rfd = RFdiffusion(pdb=lysozyme,
                        contigs='50-70/A81-140', #redesign N terminus
                        num_designs=3)
    pmpnn = ProteinMPNN(structures=rfd, 
                      num_sequences=2)
    af = AlphaFold(proteins=pmpnn)
```

## More info

Check out the following resources for a detailed explanation of it.

- **[User Manual](Docs/UserManual.md)**
- **[Tool Reference](Docs/ToolReference.md)**
- **[Examples](ExamplePipelines/)**
- **[Developer Manual](DocsDeveloperManual.md)**

