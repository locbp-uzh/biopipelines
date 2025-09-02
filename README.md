# BioPipelines

A Python package for computational protein design workflows integrating multiple modeling tools.

## Overview

BioPipelines provides a unified interface for automated protein modeling workflows using:

- **RFdiffusion/RFdiffusion-AllAtom**: Structure generation and design
- **ProteinMPNN/LigandMPNN**: Sequence design and optimization  
- **AlphaFold2/ColabFold**: Structure prediction and validation
- **Boltz1/Boltz2**: Advanced structure prediction with ligand support
- **OmegaFold**: Alternative folding method

## Quick Start

```python
from PipelineScripts.pipeline import Pipeline
from PipelineScripts.ligand_mpnn import LigandMPNN
from PipelineScripts.boltz2 import Boltz2

# Create a pipeline
pipeline = Pipeline(
    pipeline_name="MyPipeline",
    job_name="TestJob",
    job_description="Example workflow"
)

# Configure resources
pipeline.resources(
    gpu="V100",
    time="24:00:00",
    memory="16GB"
)

# Add tools to pipeline
lmpnn = pipeline.add(
    LigandMPNN(
        structures="protein.pdb",
        num_sequences=10,
        design="145-180"
    )
)

boltz = pipeline.add(
    Boltz2(
        proteins=lmpnn.output,
        ligand="SMILES_STRING_HERE"
    )
)

# Execute pipeline
pipeline.save()
pipeline.slurm()
```



## Documentation

- [User Manual](Docs/user_manual.md) - Complete usage guide
- [Developer Guide](Docs/developer_guide.md) - Development and customization
- [Architecture Overview](Docs/ARCHITECTURE_OVERVIEW.md) - System design
- [Tool I/O Reference](Docs/tool_io_reference.md) - Input/output specifications

## Example Usage

See `pipeline_setup.py` for a complete example implementing:
1. LigandMPNN sequence design
2. Boltz2 structure prediction with ligand binding
3. Filtering based on distance criteria and confidence scores

## Hardware Requirements

- **GPU**: T4/V100/A100 depending on workload
- **Memory**: 8-16GB for typical workflows
- **Environment**: Designed for SLURM cluster systems

## License

MIT License - see LICENSE file for details.

## Contributing

This package is developed by LOCBP at the University of Zürich for computational protein design research.