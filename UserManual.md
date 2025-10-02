# BioPipelines User Manual

## Getting Started

### Setup
**⚠️ Important**: Always run pipelines from the BioPipelines root directory.

```bash
cd /path/to/biopipelines
jupyter notebook my_pipeline.ipynb
```

### Basic Pipeline Structure
```python
from PipelineScripts.pipeline import Pipeline
from PipelineScripts.tool1 import Tool1
from PipelineScripts.tool2 import Tool2

# 1. Create pipeline
pipeline = Pipeline("PipelineName", "JobName", "Description")

# 2. Configure resources
pipeline.resources(gpu="V100", time="24:00:00", memory="16GB")

# 3. Add tools sequentially
tool1_output = pipeline.add(Tool1(parameters),env="YourEnv")
tool2_output = pipeline.add(Tool2(input=tool1_output, other_parameters))

# 4. Execute
pipeline.save()           # Generate scripts
pipeline.slurm()  # Submit to cluster
```

Q: requirements for each tool/environments in shares

## Tool Reference

### Structure Generation

**RFdiffusion** - Backbone structure generation
- Inputs: `pdb` (template), `contigs` (design spec), `num_designs` (default: 1)
- Outputs: `structures` [id, source_id, pdb, fixed, designed, contigs, time, status]

**RFdiffusionAllAtom** - All-atom generation with ligands
- Inputs: `ligand` (required), `pdb` (template), `contigs`, `num_designs` (default: 1)
- Outputs: `structures` [id, source_id, pdb, fixed, designed, contigs, time, status]

### Sequence Generation

**LigandMPNN** - Ligand-aware sequence design
- Inputs: `structures`, `ligand` (required), `num_sequences` (default: 1), `design_within` (default: 5.0Å)
- Outputs: `sequences` [id, sequence, sample, T, seed, overall_confidence, ligand_confidence, seq_rec]

**ProteinMPNN** - General protein sequence design
- Inputs: `structures`, `num_sequences` (default: 1), `fixed`, `redesigned`, `plddt_threshold` (default: 100.0)
- Outputs: `sequences` [id, source_id, source_pdb, sequence, score, seq_recovery, rmsd]
- Note: Sample numbering is 0-based (sample=0 is original/template, sample=1+ are designs)

**MutationComposer** - Mutation-guided sequence generation
- Inputs: `frequencies`, `num_sequences` (default: 10), `mode`, `combination_strategy`
- Outputs: `sequences` [id, sequence, mutations, mutation_positions]

**Fuse** - Protein-protein fusion design
- Inputs: `proteins/sequences`, `linker` (default: GGGGSGGGGSGGGGSGGGGS), `linker_lengths`
- Outputs: `sequences` [id, sequence, lengths]

**StitchSequences** - Combine sequences from multiple sequence generation tools
- Inputs: `sequences` (list of ToolOutputs), `selections` (position specifications)
- Outputs: `sequences` [id, sequence]

**SiteDirectedMutagenesis** - Introduce specific mutations into sequences
- Inputs: `sequences`, `mutations` (mutation specifications)
- Outputs: `sequences` [id, sequence, mutations_applied]

### Folding/Cofolding

**AlphaFold** - Structure prediction (ColabFold)
- Inputs: `sequences`, `num_relax` (default: 0), `num_recycle` (default: 3)
- Outputs: `structures` [id, source_id, sequence]

**Boltz2** - Advanced prediction with noncovalent interactions
- Inputs: `proteins`, `ligands`, `affinity` (default: True), `recycling_steps`, `diffusion_samples`
- Outputs: `confidence` [id, input_file, confidence_score, ptm, iptm, complex_plddt, complex_iplddt], `affinity` [id, input_file, affinity_pred_value, affinity_probability_binary]

### Analysis

**MutationProfiler** - Statistical analysis of sequence variations
- Inputs: `original`, `mutants`, `include_original` (default: True)
- Outputs: `profile` [position, original, count, frequency], `mutations` [position, original, A-Z], `absolute_frequencies`, `relative_frequencies`

**ResidueAtomDistance** - Distance measurements and contacts
- Inputs: `input`, `atom` (e.g., 'LIG.Cl'), `residue` (e.g., 'D in IGDWG'), `method` (default: "min")
- Outputs: `analysis` [id, source_structure, {metric_name}]

**PLIP** - Protein-ligand interaction analysis
- Inputs: `structures`, `ligand`, `output_format`, `create_pymol` (default: True)
- Outputs: `interactions` [id, ligand_id, interaction_type, residue, distance, angle, energy]

**DistanceSelector** - Distance-based residue selection
- Inputs: `structures`, `ligand`, `distance` (default: 5.0Å), `reference_type`
- Outputs: `selections` [id, within, beyond, distance_cutoff, reference_ligand]

**ProteinLigandContacts** - Analyze protein-ligand contact networks
- Inputs: `structures`, `selections` (protein regions), `ligand`, `contact_threshold` (default: 5.0Å)
- Outputs: `contact_analysis` [id, source_structure, selections, ligand, contacts, protein_ligand_distance]

**ConformationalChange** - Quantify structural changes between reference and target structures
- Inputs: `reference_structures` (single PDB or list), `target_structures`, `selection` (regions to analyze), `alignment` (align/super/cealign)
- Outputs: `conformational_analysis` [id, reference_structure, target_structure, selection, num_residues, RMSD, max_distance, mean_distance, sum_over_square_root]

**Confidence** - Extract confidence scores from structures
- Inputs: `structures`, `score_type`, `residue_range`
- Outputs: `confidence` [id, residue, confidence_score, score_type]

### Data Management

**Filter** - Expression-based result filtering
- Inputs: `data` (required), `pool`, `expression` (required), `max_items`, `sort_by`
- Outputs: `{input_name}` (filtered), `missing` [id, structure, msa]

**MergeDatasheets** - Combine analysis results
- Inputs: `datasheets` (list), `prefixes`, `key` (default: "id"), `calculate`, `id_map`
- Outputs: `merged` [key, source_structure, {prefixed_metrics}, {calculated_columns}]

**ConcatenateDatasheets** - Merge datasets across cycles
- Inputs: `datasheets` (list), `fill` (missing value filler)
- Outputs: `concatenated` [common columns, source_datasheet]

**RemoveDuplicates** - Sequence deduplication
- Inputs: `pool`, `history`, `compare` (default: "sequence")
- Outputs: Same as input (deduplicated), `missing` [id, structure, msa]

**SelectBest** - Select top-performing structures/sequences
- Inputs: `pool` (structures/sequences), `datasheets`, `metric`, `mode` ("min"/"max"), `max_items`
- Outputs: Same as input pool (filtered to best items)

**SliceDatasheet** - Extract subset of rows from datasheet
- Inputs: `datasheet`, `start_row`, `end_row`, `max_rows`
- Outputs: `sliced` [same columns as input, subset of rows]

**AverageByDatasheet** - Calculate average metrics across datasheets
- Inputs: `datasheets` (list), `group_by` (column), `metrics` (columns to average)
- Outputs: `averaged` [group_column, {averaged_metrics}]

**ExtractMetrics** - Extract specific metrics from multiple datasheets
- Inputs: `datasheets` (list), `metrics` (column names)
- Outputs: `metrics` [cycle, {extracted_metrics}]

### Visualization

**PyMOL** - Automated molecular visualization sessions
- Inputs: `structures`, `color_by`, `reference_structure`, `alignment`, `session_name`
- Outputs: .pse session files

### Sequence Analysis

**MMseqs2** - Multiple sequence alignment and homology search
- Inputs: `sequences`, `database`, `search_type`, `e_value`
- Outputs: `msas` [id, sequence_id, sequence, msa_file]

### Utilities

**LoadOutput** - Import results from previous pipelines
- Inputs: `output_json` (path), `filter`, `validate_files` (default: True)
- Outputs: Same as original tool

**FetchStructure** - Download structures from databases
- Inputs: `structure_ids`, `database` ("pdb"/"alphafold"), `format`
- Outputs: `structures` [id, structure_file, source_database]

**CompoundLibrary** - Generate compound libraries and properties
- Inputs: `smiles` (list), `properties`, `filters`
- Outputs: `compounds` [id, smiles, {calculated_properties}]

## Core Concepts

### Datasheets
All tools communicate through standardized CSV datasheets:
- **Structure tools**: `id, structure_file, [metrics]`
- **Sequence tools**: `id, sequence, source_id, [scores]`
- **Analysis tools**: `id, source_structure, [measurements]`

#### Accessing Datasheet Information

Tool outputs provide two ways to access datasheet information:

**String Paths (Default)**
```python
# Returns file path as string
path = tool.output.datasheets.sequences  # "/path/to/sequences.csv"
path = tool.o.datasheets.concatenated    # "/path/to/concatenated.csv"
```

**DatasheetInfo Objects (with Metadata)**
```python
# Returns DatasheetInfo object with metadata
info = tool.output.datasheets.info('sequences')     # DatasheetInfo object
info = tool.o.datasheets.info('concatenated')      # DatasheetInfo object

# Access properties
print(info.path)         # File path
print(info.columns)      # Column names
print(info.description)  # Description
print(info.count)        # Row count
```

**Usage in Tools**: Most tools accept both string paths and DatasheetInfo objects:
```python
# Both work for RemoveDuplicates
history=tool.output.datasheets.concatenated           # String path
history=tool.output.datasheets.info('concatenated')   # DatasheetInfo object
```

### Tool Chaining
Connect tools through their `.output` attributes:
```python
tool1 = pipeline.add(ToolA(input_params))
tool2 = pipeline.add(ToolB(input=tool1.output))  # Uses tool1's output
```

### Resource Management
Configure compute resources per pipeline:
```python
pipeline.resources(
    gpu="V100",        # GPU type: T4, V100, A100, H100, 32GB, 80GB
    memory="16GB",     # RAM allocation
    time="24:00:00",   # Wall time limit
    nodes=1            # Number of compute nodes
)
```

## Common Workflow Patterns

### Basic Design Cycle
```python
# Structure → Sequences → Prediction → Analysis
rfd = RFdiffusion(pdb="template.pdb", contigs="A1-100")
lmpnn = LigandMPNN(structures=rfd.output, ligand="ATP")
boltz = Boltz2(proteins=lmpnn.output, ligands="SMILES")
analysis = ResidueAtomDistance(structures=boltz.output, residues="A50", atoms="B101")
```

### Iterative Optimization
```python
best_structure = initial_structure
for cycle in range(5):
    sequences = LigandMPNN(structures=best_structure, num_sequences=20)
    structures = Boltz2(proteins=sequences)
    analysis = ResidueAtomDistance(structures=structures,...)
    best_structure = SelectBest(pool=structures, datasheets=analysis.distance, metric="distance",mode="min")
```

### Multi-Target Analysis
```python
# Analyze multiple conditions and merge results
apo_structures = Boltz2(proteins=sequences.output)
holo_structures = Boltz2(proteins=sequences.output, ligands="ATP")
apo_analysis = ResidueAtomDistance(structures=apo_structures.output, ...)
holo_analysis = ResidueAtomDistance(structures=holo_structures.output, ...)
combined = MergeDatasheets(
    datasheets=[apo_analysis.output, holo_analysis.output],
    prefixes=["apo_", "holo_"],
    calculate={"stability_change": "holo_energy - apo_energy"}
)
```

## Environment Configuration

### Conda Environments
Tools automatically switch between required environments:
- `ProteinEnv` - General protein tools
- `Boltz2Env` - Boltz prediction
- `ligandmpnn_env` - LigandMPNN/ProteinMPNN
- `RFDEnv` - RFdiffusion tools

### SLURM Integration
BioPipelines generates optimized SLURM submission scripts with:
- Automatic resource allocation
- Environment activation
- Job dependency management
- Log file organization

## Troubleshooting

### Common Issues
- **Path errors**: Ensure running from BioPipelines root directory
- **Environment issues**: Check conda environment availability
- **Resource limits**: Adjust GPU/memory requirements for your cluster
- **Tool failures**: Check individual tool logs in job output directories

### Debug Mode
Add debugging to pipelines:
```python
pipeline.debug = True  # Enable verbose logging
pipeline.dry_run = True  # Generate scripts without submission
```

### Log Locations
- **Pipeline logs**: `job_folder/pipeline.log`
- **Tool logs**: `job_folder/XXX_ToolName/tool.log`
- **SLURM logs**: `job_folder/slurm-JOBID.out`

## Examples

Complete pipeline examples are available in [ExamplePipelines/](ExamplePipelines/):

- **`ligandmpnn_boltz2.py`** - Basic sequence design and folding
- **`ligandmpnn_boltz2_cycle.py`** - Iterative optimization
- **`rfdaa_ligandmpnn_boltz2.py`** - Complete de novo design pipeline
- **`ligandmpnn_composer_cycle.py`** - Mutation-guided iterative design

## Advanced Features

### Custom Environments
Specify custom conda environments:
```python
tool = pipeline.add(SomeTool(env="custom_env", **params))
```

### Checkpoint Management
Load previous results:
```python
previous_run = pipeline.add(LoadOutput("/path/to/previous/job/ToolOutputs/<tool_output>.json"))
next_tool = pipeline.add(SomeTool(input=previous_run.output))
```