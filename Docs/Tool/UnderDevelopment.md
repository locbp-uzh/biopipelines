# Analysis

[← Back to Tool Reference](../ToolReference.md)

---


### PLIP (Protein-Ligand Interaction Profiler)

**UNDER DEVELOPMENT** - This tool is still being tested and refined.

Analyzes protein-ligand interactions using PLIP. Identifies hydrogen bonds, hydrophobic contacts, salt bridges, and other interaction types.

**Environment**: `biopipelines`

**Installation**: requires to download a singularity image in the containers folder.

**Note**: This tool is not fully debugged yet and may require adjustments.

**Parameters**:
- `structures`: Union[str, List[str], ToolOutput, StandardizedOutput] (required) - Input structures
- `ligand`: str = "" - Specific ligand identifier (if empty, analyzes all ligands)
- `output_format`: List[str] = None - Output formats (default: ['xml', 'txt', 'pymol'])
- `create_pymol`: bool = True - Generate PyMOL session files
- `create_images`: bool = False - Generate ray-traced images
- `analyze_peptides`: bool = False - Include protein-peptide interactions
- `analyze_intra`: bool = False - Include intra-chain interactions
- `analyze_dna`: bool = False - Include DNA/RNA interactions
- `max_threads`: int = 4 - Maximum threads for parallel processing
- `verbose`: bool = True - Enable verbose output

**Outputs**:
- `tables.interactions`:

  | id | ligand_id | interaction_type | residue | distance | angle | energy |
  |----|-----------|------------------|---------|----------|-------|--------|

**Example**:
```python
from PipelineScripts.plip import PLIP

plip = PLIP(
    structures=boltz,
    ligand="LIG",
    create_pymol=True
)
```


---

### DynamicBind

**UNDER DEVELOPMENT** - Installation failed.

Predicts ligand-specific protein-ligand complex structures using equivariant diffusion models. Recovers ligand-specific conformations from unbound protein structures.

**Key Features:**
- Predicts protein-ligand binding conformations
- Handles multiple ligands per protein
- Generates multiple poses with affinity predictions
- Outputs relaxed structures with confidence scores

**Installation:**

Requires two mamba environments.

```bash
cd /home/$USER/data
git clone https://github.com/luwei0917/DynamicBind.git
cd DynamicBind
wget https://zenodo.org/records/10183369/files/workdir.zip
unzip workdir.zip

# Request additional resources during installation
srun --mem=32G --cpus-per-task=8 --time=02:00:00 --pty bash

mamba create -n dynamicbind python=3.9 -y
mamba activate dynamicbind
pip install torch==1.13.1+cu117 torchvision==0.14.1+cu117 torchaudio==0.13.1 --extra-index-url https://download.pytorch.org/whl/cu117

# Setup LD_LIBRARY_PATH inside the environment
mkdir -p $CONDA_PREFIX/etc/mamba/activate.d
echo 'export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH' > \
     $CONDA_PREFIX/etc/mamba/activate.d/env_vars.sh
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

# Continue installation
mamba install -c conda-forge rdkit biopython scipy pandas pyyaml networkx -y
pip install numpy==1.26.4 e3nn==0.4.4 spyrmsd fair-esm tqdm
pip install \
    torch-scatter==2.1.1+pt113cu117 \
    torch-sparse==0.6.17+pt113cu117 \
    torch-cluster==1.6.1+pt113cu117 \
    torch-spline-conv \
    pyg_lib \
    -f https://data.pyg.org/whl/torch-1.13.1+cu117.html \
    --no-cache-dir
pip install torch-geometric

# Create relax environment
mamba create --name relax python=3.8
mamba activate relax
mamba install -c conda-forge openmm pdbfixer biopython openmmforcefields openff-toolkit ambertools=22 compilers -y
```

Config: `DynamicBind: null` (tool manages environments)

Model weights location: `/home/{username}/data/DynamicBind/workdir`

**Parameters:**
- `proteins`: Input protein(s) - PDB filename, list of PDBs, or tool output (all structures used)
- `ligands`: SMILES string or tool output with compounds table (must have 'smiles' column)
- `samples_per_complex`: Number of samples per complex (default: 10)
- `inference_steps`: Diffusion steps (default: 20)
- `no_relax`: Skip OpenMM relaxation (default: False)
- `hts`: High-throughput screening mode (default: False)
- `seed`: Random seed (default: 42)

**Example:**

```python
from PipelineScripts.dynamic_bind import DynamicBind
from PipelineScripts.rfdiffusion_allatom import RFdiffusionAllAtom
from PipelineScripts.compound_library import CompoundLibrary

# Basic usage - single protein, SMILES string
db = DynamicBind(
    proteins="target.pdb",
    ligands="CCO",  # SMILES string
    samples_per_complex=10
)

# Multiple proteins from PDBs folder
db = DynamicBind(
    proteins=["protein1.pdb", "protein2.pdb"],
    ligands="CC(=O)O"
)

# Chain with other tools - all structures from tool output
rfdaa = RFdiffusionAllAtom(ligand="ZIT", contigs="A1-100", num_designs=5)
compounds = CompoundLibrary(smiles_list=["CCO", "CC(=O)O"])
db = DynamicBind(proteins=rfdaa, ligands=compounds)
```

**Outputs:**
- Structures: SDF files with pattern `rank{N}_ligand_lddt{score}_affinity{score}_relaxed.sdf`
- Tables:
  - `dynamicbind_results.csv`: columns `id`, `ligand_id`, `structure`, `affinity`, `lddt`, `rank`
  - `compounds.csv`: compounds used (available as `output.compounds`)

---


---

### SequenceMetricCorrelation

Computes correlation signals between sequence mutations and performance metrics. Quantifies how specific mutations affect metrics using standardized effect sizes.

**Correlation Formulas**:

1D correlation (position-level):
```
c(i) = (m̄ᵢ₋ₘᵤₜ - m̄ᵢ₋wₜ) / √(s²ᵢ₋ₘᵤₜ + s²ᵢ₋wₜ)
```
Where:
- m̄ᵢ₋ₘᵤₜ = mean metric for sequences with mutation at position i
- m̄ᵢ₋wₜ = mean metric for sequences with wildtype at position i
- s²ᵢ₋ₘᵤₜ, s²ᵢ₋wₜ = unbiased variances (ddof=1)

2D correlation (amino acid-specific):
```
c(i,aa) = (m̄ᵢ,ₐₐ - m̄ᵢ,¬ₐₐ) / √(s²ᵢ,ₐₐ + s²ᵢ,¬ₐₐ)
```
Where:
- m̄ᵢ,ₐₐ = mean metric for sequences with amino acid 'aa' at position i
- m̄ᵢ,¬ₐₐ = mean metric for sequences with other amino acids at position i
- s²ᵢ,ₐₐ, s²ᵢ,¬ₐₐ = unbiased variances (ddof=1)

Note: When n ≤ 1 for either group, correlation is set to 0.

**Installation**: Same environment as MutationProfiler.

**Parameters**:
- `mutants`: Union[ToolOutput, StandardizedOutput, List[...]] (required) - Mutant sequences
- `data`: Union[ToolOutput, StandardizedOutput, TableInfo, str, List[...]] (required) - Table(s) with metric values
- `original`: Union[str, ToolOutput, StandardizedOutput] (required) - Reference sequence
- `metric`: str (required) - Metric column name to analyze
- `positions`: Optional[str] = None - PyMOL-style position filter (e.g., "141+143+145+147-149")

**Outputs**:
- `tables.correlation_1d`:

  | position | wt_aa | correlation | mean_mutated | mean_wt | var_mutated | var_wt | n_mutated | n_wt |
  |----------|-------|-------------|--------------|---------|-------------|--------|-----------|------|

- `tables.correlation_2d`:

  | position | wt_aa | A | C | D | E | F | G | H | I | K | L | M | N | P | Q | R | S | T | V | W | Y |
  |----------|-------|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|

**Example**:
```python
from PipelineScripts.sequence_metric_correlation import SequenceMetricCorrelation

# Single cycle
correlation = SequenceMetricCorrelation(
    mutants=filtered.tables.sequences,
    data=filtered.tables.merged,
    original=original_holo,
    metric="affinity_pred_value"
)

# Multi-cycle accumulation
correlation = SequenceMetricCorrelation(
    mutants=[cycle1.tables.sequences, cycle2.tables.sequences],
    data=[cycle1.tables.merged, cycle2.tables.merged],
    original=original_holo,
    metric="affinity_pred_value",
    positions="141+143+145+147-149"  # Filter plot to these positions
)
```

---

### BayesianAdjuster

Adjusts mutation frequencies using correlation signals via Bayesian log-odds updates. Boosts beneficial mutations and suppresses detrimental ones based on correlation evidence.

**Bayesian Update Formula**:
```
p(i,aa|c) = σ(σ⁻¹(p₀(i,aa)) + γ·c(i,aa))
```
Where:
- p₀(i,aa) = prior probability from MutationProfiler frequencies
- c(i,aa) = correlation signal from SequenceMetricCorrelation
- γ = strength hyperparameter (higher = more aggressive adjustment)
- σ(x) = sigmoid function = 1/(1+e⁻ˣ)
- σ⁻¹(p) = logit function = log(p/(1-p))

**Installation**: Same environment as MutationProfiler.

**Parameters**:
- `frequencies`: Union[TableInfo, str] (required) - Frequency table from MutationProfiler (typically absolute_frequencies)
- `correlations`: Union[TableInfo, str] (required) - Correlation table from SequenceMetricCorrelation (correlation_2d)
- `mode`: str = "min" - Optimization direction:
  - "min" = lower metric values are better (e.g., binding affinity)
  - "max" = higher metric values are better (e.g., activity)
- `gamma`: float = 3.0 - Strength hyperparameter for Bayesian update (higher = more aggressive)
- `kappa`: float = 10.0 - Pseudo-observations for sample size shrinkage (planned feature, not yet implemented). Intended to down-weight correlations from small sample sizes.
- `pseudocount`: float = 0.01 - Pseudocount added to all amino acids before adjustment. Ensures no amino acid has zero probability, allowing correlation signals to resurrect beneficial mutations. After adding, frequencies are normalized to preserve original sum.
- `positions`: Optional[str] = None - PyMOL-style position filter (e.g., "141+143+145+147-149")

**Outputs**:
- `tables.adjusted_probabilities`: Raw Bayesian-adjusted probabilities
- `tables.absolute_probabilities`: Normalized as absolute probabilities
- `tables.relative_probabilities`: Normalized as relative probabilities (original AA = 0)
- `tables.adjustment_log`: Detailed log of all adjustments

**Example**:
```python
from PipelineScripts.bayesian_adjuster import BayesianAdjuster
from PipelineScripts.mutation_composer import MutationComposer

# Basic usage
adjuster = BayesianAdjuster(
    frequencies=profiler.tables.absolute_frequencies,
    correlations=correlation.tables.correlation_2d,
    mode="min",  # Minimize affinity
    gamma=3.0,
    pseudocount=0.01
)

# More aggressive adjustment with position filter
adjuster = BayesianAdjuster(
    frequencies=profiler.tables.absolute_frequencies,
    correlations=correlation.tables.correlation_2d,
    mode="min",
    gamma=5.0,  # More aggressive
    pseudocount=0.005,  # Lower threshold for resurrected mutations
    positions="141+143+145+147-149"
)

# Use adjusted frequencies in MutationComposer
composer = MutationComposer(
    frequencies=adjuster.tables.absolute_probabilities,
    num_sequences=50,
    mode="weighted_random"
)
```

**Key Features**:
- **Pseudocount mechanism**: Allows resurrection of beneficial mutations that weren't initially present (probability = 0) but have strong positive correlations
- **Position filtering**: Consistent x-axis with MutationProfiler and SequenceMetricCorrelation for easy visual comparison
- **Dual normalization**: Provides both absolute (comparable to MutationProfiler) and relative (original AA excluded) probabilities

---

### SequenceMetricAnalysis

Analyzes correlations between sequence mutations and performance metrics across pools. Tracks position-specific mutation statistics and generates scored mutation tables for data-driven sequence optimization in iterative design cycles.

**Environment**: `ProteinEnv`

**Parameters**:
- `sequences`: Union[ToolOutput, StandardizedOutput] (required) - Sequence pool with 'id' and 'sequence' columns
- `metrics`: Union[ToolOutput, StandardizedOutput, TableInfo, str] (required) - Table with metric values (must have matching 'id' column)
- `reference_sequence`: Union[str, ToolOutput, StandardizedOutput] (required) - Reference sequence for mutation calling (can be string or tool output)
- `metric_columns`: Union[str, List[str]] (required) - Metric column name(s) to analyze
- `primary_metric`: str (required) - Which metric to use for scoring in mutation_deltas/mutation_zscores tables
- `mode`: str = "minimize" - Optimization direction: "minimize" (lower is better) or "maximize" (higher is better)
- `min_observations`: int = 3 - Minimum observation count to assign non-zero scores in delta/zscore tables
- `history`: Optional[Union[ToolOutput, StandardizedOutput]] = None - Previous analysis results to accumulate with

**Outputs**:
- `tables.mutation_statistics`:

  | position | wt_aa | mut_aa | count | affinity_mean | affinity_std | affinity_min | affinity_max | plddt_mean | plddt_std | ... |
  |----------|-------|--------|-------|---------------|--------------|--------------|--------------|------------|-----------|-----|

- `tables.mutation_deltas`: (MutationComposer-compatible)

  | position | wt_aa | A | C | D | E | F | G | H | I | K | L | M | N | P | Q | R | S | T | V | W | Y |
  |----------|-------|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|

  Values are delta scores (improvement vs WT): positive = beneficial, negative = detrimental, zero = neutral/insufficient data

- `tables.mutation_zscores`: (MutationComposer-compatible)

  | position | wt_aa | A | C | D | E | F | G | H | I | K | L | M | N | P | Q | R | S | T | V | W | Y |
  |----------|-------|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|

  Values are z-scores (standardized): positive = better than average, negative = worse than average, zero = neutral/insufficient data

- `tables.top_mutations`:

  | position | wt_aa | best_mutation | count | delta_score | zscore | affinity_mean | plddt_mean | ... |
  |----------|-------|---------------|-------|-------------|--------|---------------|------------|-----|

- `tables.coverage`:

  | position | wt_aa | n_observations | n_mutations_tested | coverage_fraction | max_count | min_count | mean_count |
  |----------|-------|----------------|-------------------|-------------------|-----------|-----------|------------|

**Example**:
```python
from PipelineScripts.sequence_metric_analysis import SequenceMetricAnalysis

# First cycle: initialize analysis
analysis = SequenceMetricAnalysis(
    sequences=filtered.tables.sequences,
    metrics=filtered.tables.merged,
    reference_sequence=original_holo,
    metric_columns=["affinity_pred_value", "complex_plddt", "contacts"],
    primary_metric="affinity_pred_value",
    mode="minimize"
)

# Subsequent cycles: accumulate with history
analysis = SequenceMetricAnalysis(
    sequences=filtered.tables.sequences,
    metrics=filtered.tables.merged,
    reference_sequence=original_holo,
    metric_columns=["affinity_pred_value", "complex_plddt"],
    primary_metric="affinity_pred_value",
    mode="minimize",
    min_observations=5,
    history=analysis  # Accumulates observations
)
```

**Key Features**:
- **Multi-metric tracking**: Analyze multiple metrics simultaneously (affinity, pLDDT, contacts, etc.)
- **Accumulation**: Pass previous analysis as history to accumulate statistics across cycles
- **Dual scoring**: Both delta scores (interpretable units) and z-scores (standardized) for flexibility
- **MutationComposer-compatible**: mutation_deltas and mutation_zscores tables use same format as MutationProfiler frequencies
- **Confidence filtering**: min_observations threshold ensures reliable mutation scores
---


---

### ESMFold

**UNDER DEVELOPMENT**

Predicts protein structures using Meta's ESM-2 with ESMFold. Fast single-sequence prediction without requiring MSAs. Models are cached in shared folder for reuse.

**References**: https://github.com/facebookresearch/esm

**Installation**: 
```bash
srun --mem=16GB --time=1:00:00 --pty bash
module load cuda/11.8.0
nvcc --version

mamba create -n esmfold python=3.9
mamba activate esmfold

mamba install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia

pip install "fair-esm[esmfold]"

pip install 'dllogger @ git+https://github.com/NVIDIA/dllogger.git'
pip install 'openfold @ git+https://github.com/aqlaboratory/openfold.git@4b41059694619831a7db195b7e0988fc4ff3a307'

```

**Environment**: `esmfold`

**Parameters**:
- `sequences`: Union[str, List[str], ToolOutput, StandardizedOutput] (required) - Input sequences
- `tables`: Optional[List[str]] = None - Input table files
- `name`: str = "" - Job name
- `num_recycle`: int = 4 - Number of recycling iterations
- `chunk_size`: Optional[int] = None - Chunk size for long sequences (auto if None)

**Outputs**:
- `structures`: List of predicted PDB files
- `tables.structures`:

  | id | sequence |
  |----|----------|

**Example**:
```python
from PipelineScripts.esmfold import ESMFold

esm = ESMFold(
    sequences=lmpnn,
    num_recycle=4
)
```
