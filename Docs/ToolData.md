# Tool I/O Reference Guide

Detailed specification of input/output formats and folder structures for all modeling tools.

## Standard Output Structure

All tools follow this consistent folder hierarchy:
```
{USER_FOLDER}/ProteinOutput/{ToolName}/{JobName}/
├── RunTime/              # Generated scripts and intermediate files
│   ├── *.sh             # Execution scripts
│   ├── *.jsonl          # Parsed data files
│   └── *.csv            # Configuration files
├── {ToolName}Output/     # Raw tool outputs
├── Results/              # Processed and ranked results
├── *.log                # Execution logs
└── *_config.txt         # Job configuration record
```

## StandardizedOutput API

### Overview
All tools now return outputs through a unified `StandardizedOutput` object that provides:
- **Dot-notation access**: `tool.output.structures`
- **Dictionary access**: `tool.output['structures']` 
- **Rich datasheet metadata**: Organized by content type with descriptions
- **Data lineage tracking**: `source_id` columns for complete traceability

### Core Structure
```python
tool.output.structures      # List[str] - Structure file paths
tool.output.structure_ids    # List[str] - Structure identifiers
tool.output.compounds       # List[str] - Compound CSV file paths  
tool.output.compound_ids     # List[str] - Compound identifiers
tool.output.sequences        # List[str] - Sequence file paths
tool.output.sequence_ids     # List[str] - Sequence identifiers
tool.output.datasheets       # DatasheetContainer - Rich metadata
tool.output.output_folder    # str - Tool's output directory
```

### Datasheet Organization
Datasheets are organized by content type with rich metadata:

```python
tool.output.datasheets.structures    # Path to structures datasheet
tool.output.datasheets.sequences     # Path to sequences datasheet  
tool.output.datasheets.compounds     # Path to compounds datasheet

# Rich metadata access
info = tool.output.datasheets.info("sequences")
print(info.columns)      # ["id", "source_id", "sequence", "score"]
print(info.description)  # "ProteinMPNN sequence generation results..."
print(info.count)        # 25 entries
```

### Data Lineage with source_id
Every tool output includes `source_id` columns for complete traceability:

**RFdiffusion** → **ProteinMPNN** → **AlphaFold** → **Boltz2**
```
structure_001 → sequence_001_1 → folded_001_1 → complex_001_1
             → sequence_001_2 → folded_001_2 → complex_001_2
```

### Datasheet Column References
Tools can explicitly reference specific columns from upstream datasheets using dot notation:

```python
# LigandMPNN can reference RFdiffusion datasheet columns
lmpnn = LigandMPNN(
    input=rfd.output,
    ligand='RFP',
    fixed_positions='input.datasheets.structures.fixed',      # References RFdiffusion's fixed column
    designed_positions='input.datasheets.structures.designed' # References RFdiffusion's designed column
)

# Mixed approach: some from datasheet, some custom
lmpnn = LigandMPNN(
    input=rfd.output, 
    ligand='RFP',
    fixed_positions='input.datasheets.structures.fixed',  # From upstream datasheet
    designed_positions='90-95+100-105'                    # Custom PyMOL selection
)
```

**Syntax**: `input.datasheets.{datasheet_name}.{column_name}`
- **input**: References the input source
- **datasheets**: Access datasheet container
- **datasheet_name**: Name of specific datasheet (e.g., "structures", "sequences")  
- **column_name**: Specific column to reference (e.g., "fixed", "designed")

### Backward Compatibility
Legacy access patterns continue to work:
```python
# Legacy - still works
tool.output.datasheets["main"]
tool.output.get("pdbs", [])

# New - preferred
tool.output.datasheets.sequences  
tool.output.structures
```

## RFdiffusion

**Purpose**: Generate protein backbone structures from contigs specification

### Input Format
- **PDB file**: Optional template structure
- **Contigs**: `"A1-100,10-20"` format specifying chain regions
- **Inpaint**: Same format as contigs for masking regions

### Output Structure
```
RFdiffusion/{JobName}/
├── RunTime/
│   ├── rfd.sh                    # RFdiffusion execution script
│   ├── rfd_pmpnn.sh             # Combined RFD+MPNN script
│   └── {job}_parsed_pdbs.jsonl   # Parsed structures for MPNN
├── RFdiffusion/
│   ├── {job}_design_0.pdb       # Generated design 0
│   ├── {job}_design_1.pdb       # Generated design 1
│   └── ...                      # Additional designs
├── rfdiffusion_results.csv      # Structure metadata
└── _rfd.log                     # RFdiffusion execution log
```

### StandardizedOutput
```python
rfd.output.structures       # ["/path/to/design_0.pdb", "/path/to/design_1.pdb"]
rfd.output.structure_ids     # ["JobName_0", "JobName_1"] 
rfd.output.datasheets.structures  # "/path/to/rfdiffusion_results.csv"

# Datasheet columns: ["id", "source_id", "pdb", "fixed", "designed", "contigs", "time", "status"]
```

### Connection to Next Tool
→ **ProteinMPNN**: Uses `rfd.output.structures` as input

---

## ProteinMPNN

**Purpose**: Generate sequences for given protein backbones

### Input Format
- **Input**: `input=rfd.output` (StandardizedOutput from RFdiffusion)
- **Fixed positions**: PyMOL selection format `"10-15+17+20-54"`
- **pLDDT threshold**: Float value for confidence-based fixing

### Output Structure
```
ProteinMPNN/{JobName}/
├── RunTime/
│   ├── pmpnn.sh                 # MPNN execution script
│   ├── {job}_parsed_pdbs.jsonl  # Parsed input structures
│   ├── {job}_fixed_pos.jsonl    # Fixed position constraints
│   └── pmpnn_sele.csv          # Selection information
├── seqs/                        # Generated sequences
│   ├── {design_name}.fa         # FASTA sequences per design
│   └── ...
├── proteinmpnn_results.csv      # Sequence metadata with scores
├── {job}_queries.csv            # Formatted for folding tools
├── {job}_queries.fasta         # Combined FASTA sequences
└── _pmpnn.log                  # MPNN execution log
```

### StandardizedOutput
```python
pmpnn.output.sequences       # ["/path/to/JobName_queries.csv"]
pmpnn.output.sequence_ids     # ["design_0_1", "design_0_2", "design_1_1"]
pmpnn.output.datasheets.sequences  # "/path/to/proteinmpnn_results.csv"

# Datasheet columns: ["id", "source_id", "source_pdb", "sequence", "score", "seq_recovery", "rmsd"]
# source_id links sequences back to their RFdiffusion structure origins
```

### Connection to Next Tool
→ **AlphaFold**: Uses `pmpnn.output` as input

---

## AlphaFold2 (ColabFold)

**Purpose**: Fold sequences into 3D structures with confidence scores

### Input Format
- **Input**: `input=pmpnn.output` (StandardizedOutput from ProteinMPNN)

### Output Structure
```
AlphaFold/{JobName}/
├── Folding/                           # Raw ColabFold outputs
│   ├── {sequence_id}/
│   │   ├── {id}_relaxed_rank_001_*.pdb    # Best relaxed structure
│   │   ├── {id}_unrelaxed_rank_001_*.pdb  # Best unrelaxed structure
│   │   ├── {id}_scores_rank_001_*.json    # Confidence scores
│   │   ├── {id}_PAE_rank_001_*.json       # Predicted Aligned Error
│   │   └── ...                            # Additional ranked models
│   └── ...
├── {sequence_id}.pdb                  # Clean extracted best structures
├── {job}_queries.csv                  # Input sequences used
└── _alphafold.log                     # Folding execution log
```

### StandardizedOutput  
```python
af.output.structures         # ["/path/to/sequence_001_1.pdb", "/path/to/sequence_001_2.pdb"]
af.output.structure_ids      # ["sequence_001_1", "sequence_001_2"]
af.output.datasheets.structures  # "/path/to/JobName_queries.csv"

# Datasheet columns: ["id", "source_id", "sequence"] 
# source_id links structures back to their ProteinMPNN sequence origins
```

### Metrics Available
- **pLDDT**: Per-residue confidence (0-100)
- **pTM**: Template modeling score
- **PAE**: Predicted aligned error matrix

### Connection to Next Tool
→ **Boltz2**: Uses `af.output` for ligand binding predictions

---

## Boltz2

**Purpose**: Advanced structure prediction with ligand support and affinity calculation

### Input Format  
- **Input**: Sequences from ProteinMPNN or LigandMPNN
- **Ligand**: SMILES string or ligand library CSV

### Output Structure
```
Boltz2/{JobName}/
├── ApoPredictions/
│   └── boltz_results_ApoConfig/
│       ├── apo_config/
│       │   ├── sample_0.pdb           # Apo structure
│       │   ├── sample_0_confidence.json
│       │   └── sample_0_summary.json
│       └── msa/                       # MSA files
├── BoundPredictions/
│   └── boltz_results_BoundConfig/
│       └── bound_config_{ligand}/     # Per-ligand predictions
├── Library/
│   ├── library_scores.csv            # Affinity rankings
│   └── {job}.pse                     # PyMOL session
├── MSAs/                              # Cached MSA files
├── apo_msas.csv                      # Input sequences with MSAs
└── *.log                             # Execution logs
```

### StandardizedOutput
```python
boltz.output.structures         # ["/path/to/BoundPredictions/..."]
boltz.output.datasheets.compounds    # "/path/to/Library/library_scores.csv"
boltz.output.datasheets.sequences    # "/path/to/apo_msas.csv"
boltz.output.datasheets.apo_structures    # Apo predictions summary
boltz.output.datasheets.bound_structures  # Bound predictions summary

# Compounds datasheet: ["id", "source_id", "sequence", "ligand_smiles", "confidence_score", "complex_plddt", "iptm", "ptm"]
# source_id links complexes back to their sequence origins (from ProteinMPNN/LigandMPNN)
```

### Confidence Metrics
- **confidence_score**: 0.8 * complex_plddt + 0.2 * iptm
- **ptm/iptm**: Template modeling scores
- **ligand_iptm**: Protein-ligand interface confidence  
- **complex_plddt**: Overall structure confidence
- **complex_pde**: Predicted distance error

### Connection to Next Tool
→ **LigandMPNN**: Uses holo structures for sequence optimization
→ **Analysis**: Comparative binding affinity analysis

---

## LigandMPNN

**Purpose**: Design sequences optimized for ligand binding

### Input Format  
- **Input**: `input=rfd.output` or structures from previous tools
- **Ligand**: Ligand identifier from structure context
- **Position control**: Multiple flexible options for specifying fixed/designed regions

### Position Control Options

#### 1. Explicit Datasheet Column References (NEW)
```python
lmpnn = LigandMPNN(
    input=rfd.output,
    ligand='RFP',
    fixed_positions='input.datasheets.structures.fixed',
    designed_positions='input.datasheets.structures.designed'
)
```

#### 2. Mixed Approach (NEW)
```python
lmpnn = LigandMPNN(
    input=rfd.output,
    ligand='RFP',
    fixed_positions='input.datasheets.structures.fixed',  # From RFdiffusion
    designed_positions='60-75'  # Custom PyMOL selection
)
```

#### 3. Auto-Datasheet Mode
```python
lmpnn = LigandMPNN(input=rfd.output, ligand='RFP')
# Automatically uses upstream datasheet for positions
```

#### 4. Direct PyMOL Selections
```python
lmpnn = LigandMPNN(
    structures=['protein.pdb'],
    ligand='RFP',
    fixed_positions='10-50+80-120',
    designed_positions='60-75'
)
```

### Output Structure
```
LigandMPNN/{JobName}/
├── seqs/                            # Generated sequences
│   ├── {structure_name}_*.fa        # FASTA sequences per structure
│   └── ...
├── {job}_queries.csv                # Formatted sequences with scores
├── {job}_queries.fasta             # Combined FASTA sequences  
├── ligandmpnn_positions.json        # Fixed/designed positions used
└── _ligandmpnn.log                 # Execution log
```

### StandardizedOutput
```python
lmpnn.output.sequences       # ["/path/to/JobName_queries.csv"]
lmpnn.output.sequence_ids    # ["structure_0_1", "structure_0_2", ...]
lmpnn.output.datasheets.sequences  # "/path/to/JobName_queries.csv"

# Datasheet columns: ["id", "source_id", "sequence", "ligand_binding_score"]
# source_id links sequences back to their structure origins
```

### Connection to Next Tool
→ **Boltz2**: Uses `lmpnn.output` to validate designed sequences

---

## OmegaFold

**Purpose**: Alternative protein folding method

### Input Format
- **FASTA sequences**: Standard protein sequences

### Output Structure
```
OmegaFold/{JobName}/
├── {sequence_id}.pdb               # Folded structure
└── _omegafold.log                  # Execution log
```

### Output Files
- **PDB structures**: Folded protein structures
- **Limited metrics**: Fewer confidence scores than AlphaFold

---

## CompoundLibrary

**Purpose**: Generate and process compound libraries from dictionaries or CSV files

### Input Format
- **Library dictionary**: Python dictionary with expansion keys
- **CSV file**: Existing compound library file
- **Primary key**: Root key for expansion (e.g., "HT7_Cy5_Cl")
- **Covalent option**: Generate CCD/PKL files for covalent ligands

### Dictionary Format
```python
library = {
    "HT7_Cy5_Cl": "*Cl_Linker**Cy5*",
    "Cl_Linker": ["ClCC","ClCCC"],
    "Cy5": ["C[*STEREOCENTER*](C)(N)C(=O)O"],
    "STEREOCENTER": ["@","@@"]
}
```

**Key differences from boltz_compound_library.py:**
- Keys use `<key>` format in expansion (not `*<key>*`)
- CSV generation happens at runtime when added to pipeline
- Uses standardized compound output format

### Output Structure
```
{PipelineName}/{JobName}/{N}_CompoundLibrary/
├── {job}_compounds.csv           # Generated compound library
├── {job}_compound_properties.csv # Molecular properties (optional)
├── {job}_summary.txt            # Generation summary
├── {job}_library_dict.json      # Original dictionary for reference
└── covalent_library/            # CCD/PKL files (covalent=True only)
    ├── compounds.csv
    └── *.ccd/*.pkl files
```

### StandardizedOutput
```python
compounds.output.compounds       # ["/path/to/{N}_CompoundLibrary/JobName_compounds.csv"]
compounds.output.compound_ids    # ["HT7_0000", "HT7_0001", "HT7_0002", ...] 
compounds.output.datasheets.compounds  # "/path/to/{N}_CompoundLibrary/JobName_compounds.csv"

# Datasheet columns: ["id", "smiles", "format", "ccd", "properties..."]
# format: "smiles" or "ccd"
# ccd: CCD identifier for covalent compounds, empty for non-covalent
```

### Connection to Next Tool
→ **Boltz2**: Uses explicit parameter `ligands=compounds.output` for ligand library predictions

---

## Parameter Passing: Explicit vs Ambiguous

### Recommended: Explicit Parameter Passing
For tools that accept multiple input types (like Boltz2), use explicit parameters to avoid ambiguity:

```python
# RECOMMENDED: Explicit parameters
boltz = pipeline.add(Boltz2(
    proteins=lmpnn.output,      # Explicitly specify proteins
    ligands=compounds.output,   # Explicitly specify ligands  
    msas=previous_boltz.output  # Explicitly specify MSAs for recycling
))
```

### Deprecated: Ambiguous Input Parameter
The generic `input=` parameter is ambiguous and should be avoided:

```python
# DISCOURAGED: Ambiguous (what does input contain?)
boltz = pipeline.add(Boltz2(
    proteins="MGHHHMNDETFRTGS",
    input=compounds.output  # Unclear - compounds could be ligands, proteins, etc.
))
```

### Supported Explicit Parameters for Boltz2
- **proteins=tool.output**: Use sequences from previous tool (ProteinMPNN, etc.)  
- **ligands=tool.output**: Use compounds from previous tool (CompoundLibrary, etc.)
- **msas=tool.output**: Use MSAs from previous tool (previous Boltz2, etc.)

## Pipeline Integration Examples

### Complete De Novo Design Pipeline
```python
# 1. Generate compound library
compounds = pipeline.add(CompoundLibrary(
    library={
        "HT7_scaffold": ["<linker><fluorophore>"],
        "linker": ["CCOCC", "CCOCCOCC"], 
        "fluorophore": ["c1ccc(N)cc1"]
    },
    primary_key="HT7_scaffold"
))

# 2. Predict protein-ligand complexes
boltz = pipeline.add(Boltz2(
    proteins="MGHHHMNDETFRTGS",
    ligands=compounds.output,  # EXPLICIT: use compounds as ligands
    affinity=True
), env="Boltz2Env")

# 3. Analyze results
results = boltz.output
print(f"Predicted {len(results.structure_ids)} structures")
print(f"Used {len(results.compound_ids)} compounds")
```

### Protein Design + Ligand Binding Pipeline
```python
# 1. Design protein sequences  
sequences = pipeline.add(ProteinMPNN(
    pdb="target.pdb",
    fixed_positions="10-15+20-25"
))

# 2. Generate compound library
compounds = pipeline.add(CompoundLibrary(
    library="drug_library.csv"
))

# 3. Predict complexes with designed sequences and compounds
boltz = pipeline.add(Boltz2(
    proteins=sequences.output,    # EXPLICIT: use sequences from ProteinMPNN
    ligands=compounds.output,     # EXPLICIT: use compounds from CompoundLibrary
    affinity=True
), env="Boltz2Env")
```

### MSA Recycling Pipeline  
```python
# 1. Generate apo structures
apo_boltz = pipeline.add(Boltz2(
    proteins="MGHHHMNDETFRTGS"
))

# 2. Generate holo structures (reusing MSAs)
holo_boltz = pipeline.add(Boltz2(
    proteins="MGHHHMNDETFRTGS",
    ligands="CCO", 
    msas=apo_boltz.output        # EXPLICIT: recycle MSAs from apo prediction
), env="Boltz2Env")
```

### Legacy Examples (Deprecated)
```python
# 1. Generate backbone structures
rfd = RFdiffusion(contigs="A1-100,10-20", num_designs=5)
rfd.configure_inputs(folders)
rfd.generate_script("rfd.sh")

# 2. Design sequences - automatic input chaining  
pmpnn = ProteinMPNN(input=rfd.output, num_sequences=3)
pmpnn.configure_inputs(folders)
pmpnn.generate_script("pmpnn.sh")

# 3. Fold sequences - automatic input chaining
af = AlphaFold(input=pmpnn.output, num_relax=1) 
af.configure_inputs(folders)
af.generate_script("alphafold.sh")

# 4. Validate with ligand binding - automatic input chaining
boltz = Boltz2(sequences=af.output, ligand_smiles="CCO")
boltz.configure_inputs(folders)
boltz.generate_script("boltz.sh")
```

### Ligand-Aware Design with Explicit Position Control
```python
# 1. Generate structures with ligand context
rfd = RFdiffusionAllAtom(contigs="A1-100,10-20", ligand="RFP", num_designs=3)
rfd.configure_inputs(folders)
rfd.generate_script("rfd.sh")

# 2. Optimize sequences with explicit position control
lmpnn = LigandMPNN(
    input=rfd.output,
    ligand='RFP',
    fixed_positions='input.datasheets.structures.fixed',      # Use RFdiffusion fixed regions
    designed_positions='input.datasheets.structures.designed' # Use RFdiffusion designed regions
)
lmpnn.configure_inputs(folders)
lmpnn.generate_script("lmpnn.sh")

# 3. Alternative: Mixed position control  
lmpnn_mixed = LigandMPNN(
    input=rfd.output,
    ligand='RFP',
    fixed_positions='input.datasheets.structures.fixed',  # From RFdiffusion datasheet
    designed_positions='90-95+100-105'                    # Custom binding site focus
)
```

### Data Lineage Tracking
```python
# Complete traceability through source_id columns:
# RFdiffusion: structure_001 (no source_id - de novo)
# ↓
# ProteinMPNN: sequence_001_1 (source_id = "structure_001") 
# ↓  
# AlphaFold: folded_001_1 (source_id = "sequence_001_1")
# ↓
# Boltz2: complex_001_1 (source_id = "folded_001_1")

# Access lineage information
print(boltz.output.datasheets.compounds)  # Shows source_id column linking back
```

### Output Access Patterns
```python
# New standardized access (recommended)
structures = rfd.output.structures          # ["/path/to/design_0.pdb", ...]
structure_ids = rfd.output.structure_ids    # ["JobName_0", "JobName_1"]
sequences_csv = pmpnn.output.datasheets.sequences  # "/path/to/proteinmpnn_results.csv"

# Legacy access (still supported)
main_datasheet = pmpnn.output.datasheets["main"]  # Still works
pdbs = rfd.output.get("pdbs", [])                 # Still works
```

### Rich Metadata Display
```python
# Get detailed information about outputs
info = pmpnn.output.datasheets.info("sequences")
print(f"Datasheet: {info.path}")
print(f"Description: {info.description}")  
print(f"Columns: {info.columns}")
print(f"Count: {info.count} sequences")

# Output:
# Datasheet: /path/to/proteinmpnn_results.csv
# Description: ProteinMPNN sequence generation results with scores and structure recovery metrics
# Columns: ['id', 'source_id', 'source_pdb', 'sequence', 'score', 'seq_recovery', 'rmsd']
# Count: 15 sequences
```

This unified interface enables seamless tool chaining with complete data provenance and rich metadata throughout the protein design pipeline.