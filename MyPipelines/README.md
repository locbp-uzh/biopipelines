# My Pipelines

This folder is for personal pipeline development and testing.

## Usage

Create your test pipelines here - they won't be committed to git:

```python
from PipelineScripts.pipeline import Pipeline
from PipelineScripts.rfdiffusion import RFdiffusion
from PipelineScripts.protein_mpnn import ProteinMPNN

pipeline = Pipeline("MyPipeline", "TestJob", "Testing new features")

rfd = pipeline.add(RFdiffusion(
    contigs="50-100",
    num_designs=5
))

pmpnn = pipeline.add(ProteinMPNN(
    structures=rfd,
    num_sequences=10
))

pipeline.slurm()
```

## Tips

- **Local testing**: Use `debug=True` when creating pipeline:
  ```python
  pipeline = Pipeline("Test", "Job", "desc", debug=True)
  ```

- **Pipeline outputs**: Go to `/shares/<user>/MyPipeline/TestJob_NNN/`

- **Keep examples**: Save working pipelines here for reference

- **Not committed**: Your pipeline files here are gitignored (except this README)

## File Organization

```
MyPipelines/
├── README.md              # This file (committed)
├── test_rfdiffusion.py    # Your pipelines (not committed)
├── binding_analysis.py    # Your pipelines (not committed)
└── ...                    # Your pipelines (not committed)
```
