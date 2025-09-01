# ProteinNotebooks User Manual

A comprehensive guide for using the computational protein design notebooks and the new pipeline system.

## Quick Start

### For Existing Jupyter Notebooks

If you're familiar with the existing notebook system, here's a quick reference:

1. **RFdiffusion.ipynb**: Generate protein backbones from scratch or with templates
2. **ProteinMPNN.ipynb**: Design sequences for protein structures  
3. **AlphaFold2.ipynb**: Fold sequences into 3D structures with confidence scores
4. **Boltz2.ipynb**: Advanced structure prediction with ligand support
5. **Combined workflows**: RFD-AA_MPNNs_Boltz2.ipynb, Boltz2_LigandMPNN_Boltz2.ipynb

### For New Python Pipeline System

```python
from ConfigScripts import RFdiffusion, ProteinMPNN, AlphaFold
from Utilities import Pipeline

# Create pipeline
pipeline = Pipeline.Pipeline("MyWorkflow")

# Add tools (automatic I/O chaining)
rfd = pipeline.add(RFdiffusion.RFdiffusion(contigs="A1-100,10-20"))
mpnn = pipeline.add(ProteinMPNN.ProteinMPNN(input_structures=rfd))
af = pipeline.add(AlphaFold.AlphaFold(sequences=mpnn))

# Generate and run
pipeline.save()  # Creates bash scripts
pipeline.slurm(gpu="V100")  # Creates SLURM submission
```

## Environment Setup

### Required Conda Environments

The system uses multiple conda environments for different tools:

```bash
# Primary environment (Python 3.9)
source activate ProteinEnv
# Tools: RFdiffusion, ProteinMPNN, AlphaFold2, OmegaFold

# Boltz2 environment (Python 3.11)  
source activate Boltz2Env
# Tools: Boltz2 predictions

# LigandMPNN environment (Python 3.11)
source activate ligandmpnn_env
# Tools: LigandMPNN (incompatible with Boltz2Env)
```

### Installation Files

Environment setup files are in `InstallationFiles/`:
- `ProteinEnv.yml`: Primary environment specification
- Installation guides for each tool

## Common Workflows

### 1. De Novo Protein Design

**Goal**: Create new proteins from scratch

**Jupyter Approach**:
```
RFdiffusion.ipynb → ProteinMPNN.ipynb → AlphaFold2.ipynb
```

**Pipeline Approach**:
```python
rfd = pipeline.add(RFdiffusion.RFdiffusion(
    contigs="A1-100",  # 100 residue protein
    num_designs=10
))

mpnn = pipeline.add(ProteinMPNN.ProteinMPNN(
    input_structures=rfd,
    num_sequences=5
))

af = pipeline.add(AlphaFold.AlphaFold(sequences=mpnn))
```

### 2. Motif Scaffolding  

**Goal**: Design proteins around functional motifs

**Key Parameters**:
- `contigs`: Specify motif and scaffold regions
- `inpaint`: Mask existing structure regions
- `active_site=True`: For small motifs

**Example**:
```python
rfd = pipeline.add(RFdiffusion.RFdiffusion(
    pdb="motif_template.pdb",
    contigs="A1-20,40-60",  # Motif + scaffold
    active_site=True,       # Small motif optimization
    inpaint="21-39"        # Mask connecting region
))
```

### 3. Ligand-Aware Design

**Goal**: Design proteins that bind specific ligands

**Jupyter Approach**:
```
Boltz2.ipynb (apo) → LigandMPNN → Boltz2.ipynb (holo)
```

**Key Features**:
- SMILES string ligand specification
- Binding affinity prediction
- Comparative apo/holo analysis

### 4. Sequence Optimization

**Goal**: Improve existing protein sequences

**Pipeline Approach**:
```python
mpnn = pipeline.add(ProteinMPNN.ProteinMPNN(
    input_structures="existing_structure.pdb",
    fixed_positions="10-15+20-25",  # Keep important residues
    plddt_threshold=70,             # Auto-fix confident regions
    num_sequences=20
))
```

## Parameter Reference

### RFdiffusion Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `pdb` | "" | Template PDB file (optional) |
| `contigs` | required | Region specification (e.g., "A1-100,10-20") |
| `inpaint` | "" | Regions to mask/redesign |
| `num_designs` | 1 | Number of backbone designs |
| `active_site` | False | Use active site model for small motifs |
| `steps` | 50 | Diffusion steps |
| `reproducible` | False | Deterministic sampling |

### ProteinMPNN Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `input_structures` | required | PDB files or previous tool output |
| `num_sequences` | 1 | Sequences per structure |
| `fixed_positions` | "" | PyMOL selection format |
| `plddt_threshold` | 100 | Auto-fix residues above threshold |
| `sampling_temp` | 0.1 | Sequence sampling temperature |
| `model_name` | "v_48_020" | ProteinMPNN model variant |
| `soluble_model` | True | Use soluble protein model |

### AlphaFold Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `sequences` | required | FASTA/CSV or previous tool output |
| `num_relax` | 0 | Models to relax with AMBER |
| `num_recycle` | 3 | Recycling iterations |
| `random_seed` | 0 | Reproducibility seed |
| `only_first` | True | Use only best model for ranking |

### Boltz2 Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `entities` | required | YAML protein/ligand specification |
| `affinity` | True | Calculate binding affinity |
| `msa_server` | "public" | MSA generation method |
| `output_format` | "pdb" | Structure output format |

## Resource Requirements

### GPU Requirements

| Tool | Minimum | Recommended | Notes |
|------|---------|-------------|--------|
| RFdiffusion | T4 | V100 | Memory scales with design size |
| ProteinMPNN | T4 | V100 | Fast, minimal GPU usage |
| AlphaFold2 | V100 | V100-32GB | Memory intensive |
| OmegaFold | V100-32GB | A100 | Model 2 requires high memory |
| Boltz2 | V100 | A100 | Complex structures need more memory |

### Time Estimates

| Workflow | Structures | Sequences | Time | GPU |
|----------|------------|-----------|------|-----|
| Basic RFD→MPNN→AF | 5 designs | 3 per design | 2-4h | V100 |
| Extensive design | 50 designs | 5 per design | 8-12h | V100 |
| Boltz2 analysis | 10 ligands | - | 4-6h | V100 |

## File Organization

### Input Files

Place input files in these locations:
- **PDB structures**: `PDBs/` folder
- **Ligand libraries**: CSV files with SMILES
- **Custom parameters**: Modify notebook cells or pipeline parameters

### Output Structure

All outputs follow this hierarchy:
```
{USER_FOLDER}/ProteinOutput/{ToolName}/{JobName}/
├── RunTime/              # Generated scripts
├── {ToolName}Output/     # Raw tool outputs  
├── Results/              # Processed results
├── *.log                # Execution logs
└── *_config.txt         # Job configuration
```

### Key Output Files

- **RFdiffusion**: `*_design_N.pdb` - Generated backbones
- **ProteinMPNN**: `seqs/*.fa` - Generated sequences  
- **AlphaFold**: Ranked PDB structures with confidence scores
- **Boltz2**: Complex predictions with binding analysis

## Troubleshooting

### Common Issues

**"PDB file not found"**
- Check file is in `PDBs/` folder
- Verify filename (with/without .pdb extension)

**"Environment not compatible"**
- Use correct environment for each tool
- Check `COMPATIBLE_ENVS` in error message

**"Out of memory"**
- Reduce `num_designs` or `num_sequences`
- Use smaller GPU memory requirements
- Switch to higher memory GPU (A100)

**"No sequences generated"**
- Check input PDB files are valid
- Verify fixed positions don't over-constrain design
- Lower `plddt_threshold` if too restrictive

### Performance Optimization

**Faster Execution**:
- Use batch processing for multiple similar jobs
- Group tools by environment to minimize switches
- Cache MSAs in Boltz2 workflows

**Memory Management**:
- Process designs in smaller batches
- Use `only_first=True` for AlphaFold to save memory
- Clean intermediate files periodically

## Migration Guide

### From Notebooks to Pipeline

**Old Notebook Approach**:
1. Edit parameters in notebook cells
2. Run cells sequentially  
3. Manually manage file paths
4. Generate SLURM scripts manually

**New Pipeline Approach**:
1. Configure tools with parameters
2. Chain tools automatically
3. Generate all scripts automatically
4. Handle environment switching

**Migration Steps**:
1. Identify your common workflows
2. Convert parameter settings to Python
3. Use `example_pipeline.py` as template
4. Test with small examples first

### Advantages of New System

- **No manual file management**: Automatic I/O chaining
- **Environment safety**: Prevents incompatible tool combinations
- **Reproducibility**: Complete parameter tracking
- **Scalability**: Easy batch processing and job arrays
- **Validation**: Parameter and dependency checking

## Support

For issues and feature requests:
- Check existing notebooks for parameter examples
- Review `example_pipeline.py` for usage patterns  
- Consult `tool_io_reference.md` for detailed I/O specs
- See `developer_guide.md` for customization