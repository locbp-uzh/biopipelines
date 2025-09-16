# LigandMPNN Tool Documentation

## Overview

LigandMPNN performs ligand-aware protein sequence design, generating sequences that are optimized for binding to specific ligands. The tool redesigns protein regions around ligand binding sites while maintaining structural compatibility.

## Basic Usage

```python
from PipelineScripts.ligand_mpnn import LigandMPNN

# Basic usage with structure file and position specification
lmpnn = pipeline.add(LigandMPNN(
    structures="protein_with_ligand.pdb",
    ligand="LIG",
    redesigned="145-180",  # PyMOL selection format
    num_sequences=10
))

# Using tool output as input (standardized format)
lmpnn = pipeline.add(LigandMPNN(
    input=rfdiffusion.output,  # Uses standardized input
    ligand="ATP", 
    num_sequences=20
))

# Legacy format (still supported)
lmpnn = pipeline.add(LigandMPNN(
    structures=boltz.output,
    ligand="HEM",
    num_sequences=10
))
```

## Parameters

### Input Parameters
- **`input`** (ToolOutput | StandardizedOutput | Dict): Standardized input with structures and datasheets
- **`structures`** (str | List[str] | ToolOutput): Input structures (legacy format)
  - File path to PDB file(s), or output from another tool

### Required Parameters  
- **`ligand`** (str): Ligand identifier/residue name in the structure
  - Must match the ligand name in the PDB file (e.g., "LIG", "ATP", "HEM")

### Core Parameters
- **`num_sequences`** (int, default=1): Total number of sequences to generate
- **`fixed`** (str): Fixed positions in LigandMPNN format ("A3 A4 A5") or datasheet reference
- **`redesigned`** (str): Positions to redesign in LigandMPNN format ("A3 A4 A5") or datasheet reference
- **`design_within`** (float, default=5.0): Distance (Å) from ligand for automatic position selection

### Processing Parameters
- **`model`** (str, default="v_32_010"): LigandMPNN model version
- **`batch_size`** (int, default=1): Batch size for processing
- **`name`** (str): Custom job name for output files
- **`datasheets`** (List[str]): Input datasheet files for metadata

### Advanced Parameters (via kwargs)
Additional LigandMPNN-specific parameters can be passed through kwargs for fine-tuning.

## Input Formats

### Standardized Input (Recommended)
```python
# Use .output from previous tools
lmpnn = LigandMPNN(input=previous_tool.output, ligand="ATP")
```

### Structure Input
- **PDB files**: Single file path or list of paths
- **Tool outputs**: ToolOutput objects from RFdiffusion, Boltz2, etc.

### Position Specification
- **PyMOL selection format**: Range notation ("145-180", "50-75+100-120") - **Recommended**
- **LigandMPNN format**: Space-separated residue identifiers ("A145 A146 A147 A148") 
- **Datasheet references**: Reference to datasheet columns (e.g., "input.datasheets.structures.fixed")
- **Automatic selection**: Uses `design_within` distance if positions not specified

## Output Structure

```
output_folder/
├── sequences.csv              # Main output - generated sequences with metadata
├── {name}_queries.csv         # Detailed sequence information
├── {name}_queries.fasta       # FASTA format sequences
├── parsed_pdbs.jsonl          # Processed structure information
├── lmpnn_commands.sh          # Generated LigandMPNN commands
├── lmpnn_positions_replacement.sh  # Position constraint script
└── logs/                      # Execution and debug logs
```

### Primary Output: sequences.csv
| Column | Description |
|--------|-------------|
| `id` | Unique sequence identifier |
| `sequence` | Generated amino acid sequence |
| `score` | LigandMPNN confidence score |
| `source_id` | Source structure identifier |
| `T` | Temperature parameter used |

## Position Selection Logic

LigandMPNN determines redesign positions through:

1. **Manual specification**: `fixed` and `redesigned` parameters override automatic selection
2. **Datasheet integration**: Reads position constraints from input datasheets
3. **Distance-based**: Uses `design_within` parameter for automatic selection around ligand
4. **Format conversion**: PyMOL ranges automatically converted to LigandMPNN format

### Manual Position Control

**PyMOL Selection Format (Recommended)**
```python
# Range-based specification (most common)
lmpnn = LigandMPNN(
    structures="complex.pdb",
    ligand="ATP",
    fixed="1-50+75-100",        # Keep residues 1-50 and 75-100 fixed
    redesigned="145-180"        # Redesign residues 145-180 only
)
```

**LigandMPNN Format (Also Supported)**
```python
# Explicit residue specification  
lmpnn = LigandMPNN(
    structures="complex.pdb",
    ligand="ATP",
    fixed="A1 A2 A3 B10 B11",           # Keep these specific residues fixed
    redesigned="A145 A146 A147 A148"    # Redesign only these residues
)
```

### Automatic Position Selection
```python
# Distance-based selection around ligand
lmpnn = LigandMPNN(
    structures="complex.pdb",
    ligand="ATP", 
    design_within=4.0  # Redesign residues within 4Å of ligand
)
```

## Integration with Other Tools

### Common Workflow Patterns

**Structure Generation → Sequence Design**
```python
# Generate structures, then design sequences
rfd = RFdiffusion(pdb="template.pdb", contigs="A1-50/B1-30")
lmpnn = LigandMPNN(input=rfd.output, ligand="HEM", num_sequences=50)
```

**Sequence Design → Structure Prediction**
```python
# Design sequences, then predict structures
lmpnn = LigandMPNN(structures="protein.pdb", ligand="ATP", num_sequences=20)
boltz = Boltz2(input=lmpnn.output, ligands="SMILES_STRING")
```

**Iterative Design with Analysis**
```python
# Multi-cycle optimization
for cycle in range(3):
    lmpnn = LigandMPNN(input=best_structures.output, ligand="LIG")
    boltz = Boltz2(input=lmpnn.output, ligands="CCO")
    analysis = ResidueAtomDistance(structures=boltz.output, ...)
    best_structures = SelectBest(pool=boltz.output, datasheets=analysis.output)
```

## Batch Processing

LigandMPNN calculates batches automatically:
```python
# Generate 100 sequences in batches of 10
lmpnn = LigandMPNN(
    structures="protein.pdb",
    ligand="ATP",
    num_sequences=100,  # Total sequences
    batch_size=10       # Process 10 at a time (10 batches)
)
```

## Model Versions

Available LigandMPNN models:
- **`v_32_010`** (default): Standard LigandMPNN model
- Other versions available based on LigandMPNN installation

## Environment Requirements

- **Conda environment**: `ligandmpnn_env`
- **Compatible environments**: Only `ligandmpnn_env` supported
- **Resources**: 
  - **GPU**: V100 recommended (default)
  - **Memory**: 16GB default
  - **Time**: 12:00:00 default

## Performance Considerations

### Resource Scaling
- **GPU usage**: V100 optimal for standard workflows
- **Memory scaling**: Increases with structure size and batch size
- **Time scaling**: Roughly linear with num_sequences

### Optimization Tips
- **Batch size**: Adjust based on GPU memory (1-10 typically)
- **Sequence count**: Start with 10-50 sequences for testing
- **Position constraints**: Specific constraints faster than distance-based selection

## Troubleshooting

### Common Issues

**"Ligand not found in structure"**
- Verify ligand name exactly matches PDB residue name
- Check all input structures contain the specified ligand
- Use molecular viewers to confirm ligand naming

**"Invalid position format"**
- Use LigandMPNN format: "A10 A11 B5" (space-separated)
- Check chain IDs match structure file
- Verify residue numbers exist in structure

**"Model not found"**
- Confirm LigandMPNN installation and model availability
- Check conda environment activation
- Verify model version spelling

### Input Validation
The tool validates:
- Input structures exist and are accessible
- Ligand parameter is non-empty
- num_sequences and batch_size are positive
- design_within is positive

## Related Tools

- **[ProteinMPNN](ProteinMPNN.md)**: General protein sequence design (no ligand awareness)
- **[MutationComposer](MutationComposer.md)**: Mutation-guided sequence generation
- **[Boltz2](Boltz2.md)**: Structure prediction with ligand support
- **[RFdiffusion](RFdiffusion.md)**: Structure generation for LigandMPNN input

## Examples

Complete workflow examples in [ExamplePipelines/](../ExamplePipelines/):
- `ligandmpnn_boltz2.py` - Basic LigandMPNN → Boltz2 workflow
- `ligandmpnn_boltz2_cycle.py` - Iterative optimization cycles
- `rfdaa_ligandmpnn_boltz2.py` - Full pipeline from structure generation
- `ligandmpnn_composer_cycle.py` - Mutation-guided iterative design