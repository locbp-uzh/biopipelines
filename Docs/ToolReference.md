# BioPipelines Tool Reference

## Index

### Tool Reference

- [Structure Generation](#structure-generation)
  - [RFdiffusion](#rfdiffusion)
  - [RFdiffusionAllAtom](#rfdiffusionallatom)
  - [BoltzGen](#boltzgen)
- [Sequence Design](#sequence-design)
  - [ProteinMPNN](#proteinmpnn)
  - [LigandMPNN](#ligandmpnn)
  - [MutationComposer](#mutationcomposer)
  - [SDM (SiteDirectedMutagenesis)](#sdm-sitedirectedmutagenesis)
  - [Fuse](#fuse)
  - [StitchSequences](#stitchsequences)
  - [DNAEncoder](#dnaencoder)
- [Structure Prediction](#structure-prediction)
  - [AlphaFold](#alphafold)
  - [ESMFold](#esmfold)
  - [Boltz2](#boltz2)
  - [RF3](#rf3)
  - [OnionNet](#onionnet)
  - [OnionNet2](#onionnet2)
- [Analysis](#analysis)
  - [DynamicBind](#dynamicbind)
  - [ResidueAtomDistance](#residueatomdistance)
  - [PLIP (Protein-Ligand Interaction Profiler)](#plip-protein-ligand-interaction-profiler)
  - [DistanceSelector](#distanceselector)
  - [SelectionEditor](#selectioneditor)
  - [ConformationalChange](#conformationalchange)
  - [MutationProfiler](#mutationprofiler)
  - [SequenceMetricAnalysis](#sequencemetricanalysis)
  - [ProteinLigandContacts](#proteinligandcontacts)
- [Data Management](#data-management)
  - [Filter](#filter)
  - [Rank](#rank)
  - [SelectBest](#selectbest)
  - [RemoveDuplicates](#removeduplicates)
  - [MergeTables](#mergetables)
  - [ConcatenateTables](#concatenatetables)
  - [SliceTable](#slicetable)
  - [ExtractMetrics](#extractmetrics)
  - [AverageByTable](#averagebytable)
- [Utilities](#utilities)
  - [LoadOutput](#loadoutput)
  - [MMseqs2Server](#mmseqs2server)
  - [CompoundLibrary](#compoundlibrary)
  - [PDB](#pdb)
  - [PyMOL](#pymol)


# BioPipelines Tool Reference

This document contains the complete, verified tool reference with all parameters, default environments, and output specifications.

## Structure Generation

### RFdiffusion

Generates novel protein backbone structures using diffusion models. Designs de novo proteins or scaffolds functional motifs into new contexts.

**References**: https://www.nature.com/articles/s41586-023-06415-8

**Installation**: Go to the data directory, then clone the RFdiffusion git repository and download the weights as indicated in https://github.com/RosettaCommons/RFdiffusion. The environment they provide won't work on our cluster, so instead go to the biopipelines directory and run this:
```bash
conda env create -f Environments/ProteinEnv.yaml
conda activate ProteinEnv
pip install -r Environments/ProteinEnv_pip_requirements.txt
```
Followed by the original instructions:
```bash
cd env/SE3Transformer
pip install --no-cache-dir -r requirements.txt
python setup.py install
cd ../.. # change into the root directory of the repository
pip install -e . # install the rfdiffusion module from the root of the repository
```

**Parameters**:
- `pdb`: Union[str, ToolOutput, StandardizedOutput] = "" - Input PDB template (optional)
- `contigs`: str (required) - Contig specification (e.g., "A1-100", "10-20,A6-140")
- `inpaint`: str = "" - Inpainting specification
- `num_designs`: int = 1 - Number of backbone designs to generate
- `active_site`: bool = False - Use active site model for small motifs (<30 residues)
- `steps`: int = 50 - Number of diffusion steps
- `partial_steps`: int = 0 - Partial diffusion steps
- `reproducible`: bool = False - Use deterministic sampling
- `design_startnum`: int = 1 - Starting number for design IDs

**Outputs**:
- `structures`: List of generated PDB files
- `tables.structures`:

  | id | source_id | pdb | fixed | designed | contigs | time | status |
  |----|-----------|-----|-------|----------|---------|------|--------|

**Example**:
```python
rfd = RFdiffusion(
    contigs="50-100",
    num_designs=10

```

---

### RFdiffusionAllAtom

Generates protein structures with explicit modeling of ligands and small molecules. Diffusion model that handles all-atom representation including non-protein entities.

**References**: https://www.science.org/doi/10.1126/science.adl2528.

**Installation**: Works in ProteinEnv environment as installed in RFdiffusion.

**Parameters**:
- `ligand`: str (required) - Ligand identifier in PDB (e.g., 'LIG', 'ATP')
- `pdb`: Union[str, ToolOutput, StandardizedOutput] = "" - Input PDB template
- `contigs`: str (required) - Contig specification
- `inpaint`: str = "" - Inpainting specification
- `num_designs`: int = 1 - Number of designs
- `active_site`: bool = False - Use active site model
- `steps`: int = 200 - Diffusion steps (default higher for AllAtom)
- `partial_steps`: int = 0 - Partial diffusion steps
- `reproducible`: bool = False - Deterministic sampling
- `design_startnum`: int = 1 - Starting design number
- `ppi_design`: bool = False - Enable protein-protein interface design
- `ppi_hotspot_residues`: List[str] = None - Hotspot residues for PPI
- `ppi_binder_length`: int = None - Length of PPI binder
- `autogenerate_contigs`: bool = False - Auto-infer fixed contig segments
- `model_only_neighbors`: bool = False - Only remodel residues near ligand
- `num_recycles`: int = 1 - Number of diffusion recycles
- `scaffold_guided`: bool = False - Enforce strict adherence to scaffold
- `align_motif`: bool = True - Pre-align functional motif
- `deterministic`: bool = False - Fixed RNG seeds
- `inpaint_str`: str = None - Secondary structure pattern for inpainting
- `inpaint_seq`: str = None - Sequence pattern for inpainting
- `inpaint_length`: int = None - Target length for inpainted regions
- `guiding_potentials`: str = None - Custom external potentials

**Outputs**:
- `structures`: List of generated PDB files
- `tables.structures`:

  | id | source_id | pdb | fixed | designed | contigs | time | status |
  |----|-----------|-----|-------|----------|---------|------|--------|

**Example**:
```python
rfdaa = RFdiffusionAllAtom(
    pdb=template,
    ligand='LIG',
    contigs='10-20,A6-140',
    num_designs=5

```

---

### BoltzGen

Generates protein binders (proteins, peptides, or nanobodies) targeting specified molecules using an end-to-end pipeline that combines diffusion-based backbone generation, inverse folding, structure prediction, and multi-metric filtering.

**References**: https://github.com/HannesStark/boltzgen

**Installation**: Install BoltzGen via pip with Python ≥3.11:
```bash
conda create -n boltzgen python=3.11
conda activate boltzgen
pip install boltzgen
```

Models (~6GB) download automatically to `~/.cache` on first run, or specify custom cache with `cache_dir` parameter.

**Parameters**:
- `design_spec`: Union[str, Dict] (required) - YAML configuration string/dict or path to YAML file defining:
  - Target entities (proteins, ligands from .cif/.pdb files)
  - Binder specification (sequence ranges like "80..140")
  - Binding site constraints (binding_types, not_binding)
  - Secondary structure specifications (helix/sheet/loop)
  - Structural constraints (disulfide bonds, covalent connections)
- `protocol`: str = "protein-anything" - Design protocol:
  - "protein-anything": General protein binder design
  - "peptide-anything": Peptide binder design
  - "protein-small_molecule": Protein-small molecule binder design
  - "nanobody-anything": Nanobody design
- `num_designs`: int = 10000 - Total designs to generate
- `budget`: int = 100 - Final number after diversity filtering
- `design_checkpoints`: Optional[List[str]] = None - Model checkpoint paths
- `step_scale`: Optional[float] = None - Diffusion step scale
- `noise_scale`: Optional[float] = None - Diffusion noise scale
- `diffusion_batch_size`: Optional[int] = None - Samples per trunk run (auto if None)
- `inverse_fold_num_sequences`: int = 1 - Sequences per backbone
- `skip_inverse_folding`: bool = False - Skip sequence redesign
- `alpha`: float = 0.5 - Quality/diversity trade-off (0.0=quality only, 1.0=diversity only)
- `filter_biased`: bool = True - Remove amino acid composition outliers
- `additional_filters`: Optional[List[str]] = None - Hard threshold expressions (e.g., ["design_ALA>0.3"])
- `metrics_override`: Optional[Dict[str, float]] = None - Per-metric ranking weights
- `devices`: Optional[int] = None - Number of GPUs (auto-detect if None)
- `reuse`: bool = False - Resume interrupted runs
- `steps`: Optional[List[str]] = None - Run only specific pipeline steps
- `cache_dir`: Optional[str] = None - Model download location

**Outputs**:
- `structures`: Final designed structure files (in `final_ranked_designs/final_<budget>_designs/`)
- `tables.all_metrics`:

  | id | RMSD | hydrogen_bonds | packing_quality | interface_contacts | binding_energy | design_plddt |
  |----|------|----------------|-----------------|-------------------|----------------|--------------|

- `tables.final_metrics`:

  | id | RMSD | hydrogen_bonds | packing_quality | interface_contacts | binding_energy | design_plddt | rank |
  |----|------|----------------|-----------------|-------------------|----------------|--------------|------|

- `tables.aggregate_metrics`:

  | metric | mean | std | min | max |
  |--------|------|-----|-----|-----|

- `tables.per_target_metrics`:

  | target_id | num_designs | avg_rmsd | avg_plddt |
  |-----------|-------------|----------|-----------|

**Important Notes**:
- **Residue indexing**: All residue indices start at 1 and use canonical mmcif `label_asym_id`, not `auth_asym_id`
- **Sequence specification**: Use ranges like "80..140" for random length, "15..20AAAA" for random prefix + fixed tail
- **Output structure**: Pipeline generates intermediate designs, inverse-folded structures, and final ranked designs with comprehensive metrics

**Example**:
```python
# YAML design specification
design_yaml = """
entities:
  - protein:
      id: C
      sequence: 80..140
  - file:
      path: target.cif
      include:
        - chain:
            id: A
"""

boltzgen = BoltzGen(
    design_spec=design_yaml,
    protocol="protein-anything",
    num_designs=10000,
    budget=100,
    alpha=0.5


# Access final designs
best_designs = boltzgen.tables.final_metrics
```

**Example with File**:
```python
# Using external YAML file
boltzgen = BoltzGen(
    design_spec="designs/my_binder.yaml",
    protocol="peptide-anything",
    num_designs=5000,
    budget=50,
    additional_filters=["design_ALA<0.3", "filter_rmsd_design<2.5"]

```

---

## Sequence Design

### ProteinMPNN

Designs protein sequences for given backbone structures. Uses graph neural networks to optimize sequences for structure stability while respecting fixed/designed region constraints.

**References**: https://www.science.org/doi/10.1126/science.add2187.

**Installation**: Go to your data folder and clone the official repository (https://github.com/dauparas/ProteinMPNN). The model will then work in the same environment as RFdiffusion.
```bash
git clone https://github.com/dauparas/ProteinMPN
```

**Parameters**:
- `structures`: Union[str, List[str], ToolOutput] (required) - Input structures
- `tables`: Optional[List[str]] = None - Input table files
- `num_sequences`: int = 1 - Number of sequences per structure
- `fixed`: str = "" - Fixed positions (PyMOL selection or table reference)
- `redesigned`: str = "" - Redesigned positions (PyMOL selection or table reference)
- `fixed_chain`: str = "A" - Chain to apply fixed positions
- `plddt_threshold`: float = 100.0 - pLDDT threshold for automatic fixing (residues above threshold are fixed)
- `sampling_temp`: float = 0.1 - Sampling temperature
- `model_name`: str = "v_48_020" - ProteinMPNN model variant
- `soluble_model`: bool = True - Use soluble protein model

**Outputs**:
- `sequences`: CSV file with generated sequences
- `tables.sequences`:

  | id | source_id | source_pdb | sequence | score | seq_recovery | rmsd |
  |----|-----------|------------|----------|-------|--------------|------|

**Note**: Sample 0 is the original/template sequence, samples 1+ are designs.

**Example**:
```python
pmpnn = ProteinMPNN(
    structures=rfd,
    num_sequences=10,
    fixed="1-10+50-60",
    redesigned="20-40"

```

---

### LigandMPNN

Designs protein sequences optimized for ligand binding. Specialized version of ProteinMPNN that considers protein-ligand interactions during sequence design.

**References**: https://www.nature.com/articles/s41592-025-02626-1.

**Installation**: As from the official repository (https://github.com/dauparas/LigandMPNN), go to your data folder then run:
```bash
git clone https://github.com/dauparas/LigandMPNN.git
cd LigandMPNN
bash get_model_params.sh "./model_params"
conda create -n ligandmpnn_env python=3.11
pip3 install -r requirements.txt
```

**Parameters**:
- `structures`: Union[str, List[str], ToolOutput] (required) - Input structures
- `ligand`: str (required) - Ligand identifier for binding site focus
- `tables`: Optional[List[str]] = None - Input table files
- `name`: str = "" - Job name for output files
- `num_sequences`: int = 1 - Number of sequences to generate per structure
- `fixed`: str = "" - Fixed positions (LigandMPNN format "A3 A4 A5" or table reference)
- `redesigned`: str = "" - Designed positions (LigandMPNN format or table reference)
- `design_within`: float = 5.0 - Distance in Angstroms from ligand for post-generation analysis only (does not control design). For actually designing residues within a distance, use [DistanceSelector](#distanceselector) to select positions first.
- `model`: str = "v_32_010" - LigandMPNN model version (v_32_005, v_32_010, v_32_020, v_32_025)
- `batch_size`: int = 1 - Batch size for processing

**Outputs**:
- `sequences`: CSV file with generated sequences
- `tables.sequences`:

  | id | sequence | sample | T | seed | overall_confidence | ligand_confidence | seq_rec |
  |----|----------|--------|---|------|-------------------|-------------------|---------|

**Example**:
```python
lmpnn = LigandMPNN(
    structures=rfdaa,
    ligand="LIG",
    num_sequences=5,
    redesigned=rfdaa.tables.structures.designed

```

---

### MutationComposer

Generates new protein sequences by composing mutations based on frequency analysis. Creates combinatorial mutants from mutation profiles with different sampling strategies.

**Installation**: Same environment as MutationProfiler.

**Parameters**:
- `frequencies`: Union[List, ToolOutput, StandardizedOutput, TableInfo, str] (required) - Mutation frequency table(s) from MutationProfiler
- `num_sequences`: int = 10 - Number of sequences to generate
- `mode`: str = "single_point" - Generation strategy:
  - "single_point": One mutation per sequence
  - "weighted_random": Random mutations weighted by frequency
  - "hotspot_focused": Focus on high-frequency positions
  - "top_mutations": Use only top N mutations
- `min_frequency`: float = 0.01 - Minimum frequency threshold for mutations
- `max_mutations`: int = None - Maximum mutations per sequence
- `random_seed`: int = None - Random seed for reproducibility
- `prefix`: str = "" - Prefix for sequence IDs
- `hotspot_count`: int = 10 - Number of top hotspot positions (for hotspot_focused mode)
- `combination_strategy`: str = "average" - Strategy for combining multiple tables (average, maximum, stack, round_robin)

**Outputs**:
- `sequences`: CSV file with composed sequences
- `tables.sequences`:

  | id | sequence | mutations | mutation_positions |
  |----|----------|-----------|-------------------|

**Example**:
```python
profiler = MutationProfiler(original=ref, mutants=variants
composer = MutationComposer(
    frequencies=profiler.tables.relative_frequencies,
    num_sequences=50,
    mode="weighted_random",
    max_mutations=5

```

---

### SDM (SiteDirectedMutagenesis)

Performs site-directed mutagenesis at specified positions. Generates systematic amino acid substitutions for experimental library design or computational scanning.

**Environment**: `ProteinEnv`

**Parameters**:
- `original`: Union[str, ToolOutput, StandardizedOutput] (required) - Input structure/sequence
- `position`: int (required) - Target position for mutagenesis (1-indexed)
- `mode`: str = "saturation" - Mutagenesis strategy:
  - "saturation": All 20 amino acids
  - "hydrophobic": Hydrophobic residues only
  - "hydrophilic": Hydrophilic residues only
  - "charged": Charged residues only
  - "polar": Polar residues only
  - "nonpolar": Nonpolar residues only
  - "aromatic": Aromatic residues only
  - "aliphatic": Aliphatic residues only
  - "positive": Positively charged residues only
  - "negative": Negatively charged residues only
- `include_original`: bool = False - Include original amino acid in output
- `exclude`: str = "" - Amino acids to exclude (single letter codes as string, e.g., "CP")
- `prefix`: str = "" - Prefix for sequence IDs

**Outputs**:
- `sequences`: CSV file with mutant sequences
- `tables.sequences`:

  | id | sequence | mutation | position | original_aa | new_aa |
  |----|----------|----------|----------|-------------|--------|

- `tables.missing_sequences`:

  | id | sequence | reason |
  |----|----------|--------|

**Example**:
```python
sdm = SDM(
    original=template,
    position=42,
    mode="saturation",
    exclude="CP"

```

---

### Fuse

Concatenates multiple protein sequences with flexible linkers. Creates fusion proteins with customizable linker lengths for domain engineering.

**Environment**: `ProteinEnv`

**Parameters**:
- `proteins`: Union[List[str], str] (required) - List of protein sequences or PDB file paths
- `sequences`: Union[List[str], str] = None - Alias for proteins
- `name`: str = "" - Job name for output files
- `linker`: str = "GGGGSGGGGSGGGGSGGGGS" - Linker sequence that will be cut based on `linker_lengths` if specified
- `linker_lengths`: List[str] = None - List of length ranges for each junction to generate multiple variants by cutting the linker (e.g., ["1-6", "1-6"])

**Outputs**:
- `sequences`: CSV file with fused sequences
- `tables.sequences`:

  | id | sequence | lengths |
  |----|----------|---------|

**Example**:
```python
fused = Fuse(
    proteins=[domain1, domain2, domain3],
    linker="GGGGS"

```

---

### StitchSequences

Combines different regions from multiple sequence sources into chimeric proteins. Useful for creating hybrid sequences from different design outputs.

**Environment**: `ProteinEnv`

**Parameters**:
- `sequences`: List[Union[ToolOutput, StandardizedOutput]] (required) - List of sequence outputs to stitch
- `selections`: Union[List[Union[str, ToolOutput]], str] = None - Position specifications for each sequence (e.g., ["1-50", "51-100"])
- `id_map`: Dict[str, str] = None - ID mapping pattern (default: {"*": "*_N"})

**Outputs**:
- `sequences`: CSV file with stitched sequences
- `tables.sequences`:

  | id | sequence |
  |----|----------|

**Example**:
```python
stitched = StitchSequences(
    sequences=[lmpnn1, lmpnn2],
    selections=["1-50", "51-100"]

```

---

### DNAEncoder

Reverse-translates protein sequences to DNA with organism-specific codon optimization. Uses thresholded weighted codon sampling based on CoCoPUTs genome frequency tables.

**Environment**: `biopipelines`

**Parameters**:
- `sequences`: Union[ToolOutput, StandardizedOutput] (required) - Input protein sequences
- `organism`: str = "EC" - Target organism for codon optimization:
  - "EC" (Escherichia coli)
  - "SC" (Saccharomyces cerevisiae)
  - "HS" (Homo sapiens)
  - Combinations: "EC&HS", "EC&SC", "HS&SC", "EC&HS&SC"

**Outputs**:
- `tables.dna`:

  | id | protein_sequence | dna_sequence | organism | method |
  |----|------------------|--------------|----------|--------|

- Excel file with color-coded codons (red <5‰, orange 5-10‰, black ≥10‰)

**Example**:
```python
dna = DNAEncoder(
    sequences=lmpnn,
    organism="EC&HS"  # Conservative optimization for both E. coli and human
)
```

**Note**: Uses thresholded weighted sampling (codons ≥10‰, fallback to ≥5‰). For multi-organism optimization, uses minimum frequency across organisms. Please cite CoCoPUTs (HIVE) when using.

---

## Structure Prediction

### AlphaFold

Predicts protein structures from amino acid sequences using AlphaFold2. Generates high-confidence 3D models with optional relaxation.

**Environment**: `ProteinEnv`

**Parameters**:
- `sequences`: Union[str, List[str], ToolOutput, Dict[str, Any]] (required) - Input sequences or dict with sequences
- `tables`: Optional[List[str]] = None - Input table files
- `name`: str = "" - Job name
- `num_relax`: int = 0 - Number of best models to relax with AMBER
- `num_recycle`: int = 3 - Number of recycling iterations
- `rand_seed`: int = 0 - Random seed (0 = random)

**Outputs**:
- `structures`: List of predicted PDB files
- `tables.structures`:

  | id | source_id | sequence |
  |----|-----------|----------|

- `tables.confidence`:

  | id | structure | plddt | max_pae | ptm |
  |----|-----------|-------|---------|-----|

**Example**:
```python
af = AlphaFold(
    sequences=lmpnn,
    num_relax=1,
    num_recycle=5

```

---

### ESMFold

Predicts protein structures using Meta's ESM-2 with ESMFold. Fast single-sequence prediction without requiring MSAs. Models are cached in shared folder for reuse.

**References**: https://github.com/facebookresearch/esm

**Installation**: ESMFold works in the ProteinEnv environment. Install requirements:
```bash
conda activate ProteinEnv
pip install "fair-esm[esmfold]"
pip install 'dllogger @ git+https://github.com/NVIDIA/dllogger.git'
pip install 'openfold @ git+https://github.com/aqlaboratory/openfold.git@4b41059694619831a7db195b7e0988fc4ff3a307'
```

**Environment**: `ProteinEnv`

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
esm = ESMFold(
    sequences=lmpnn,
    num_recycle=4

```

---

### Boltz2

Predicts biomolecular complexes including proteins, nucleic acids, and small molecules. State-of-the-art model for protein-ligand and protein-protein complex prediction.

**Environment**: `Boltz2Env`

**Parameters**:
- `config`: Optional[str] = None - Direct YAML configuration string
- `proteins`: Union[str, List[str], ToolOutput] (required) - Protein sequences
- `ligands`: Union[str, ToolOutput, StandardizedOutput, None] = None - Ligand SMILES string or compound library ToolOutput
- `msas`: Optional[Union[str, ToolOutput]] = None - Pre-computed MSA files for recycling (pass entire ToolOutput, not .msas)
- `sequences`: Union[str, List[str], ToolOutput] = None - Legacy parameter (use proteins instead)
- `ligand_smiles`: Optional[str] = None - Legacy parameter (use ligands instead)
- `ligand_library`: Optional[str] = None - Path to CSV file with ligand library (deprecated, use CompoundLibrary tool)
- `primary_key`: Optional[str] = None - Key column in library to filter by
- `library_repr`: str = "SMILES" - Ligand representation (SMILES, CCD)
- `library_type`: str = "noncovalent" - Binding type (noncovalent, covalent)
- `affinity`: bool = True - Calculate binding affinity predictions
- `output_format`: str = "pdb" - Output format (pdb, mmcif)
- `msa_server`: str = "public" - MSA generation (public, local)
- `global_msas_cache`: bool = False - Enable global MSA caching across jobs
- `recycling_steps`: Optional[int] = None - Number of recycling steps (default: model-specific)
- `diffusion_samples`: Optional[int] = None - Number of diffusion samples (default: model-specific)
- `use_potentials`: bool = False - Enable external potentials

**Outputs**:
- `structures`: List of predicted complex PDB files
- `tables.confidence`:

  | id | input_file | confidence_score | ptm | iptm | complex_plddt | complex_iplddt |
  |----|------------|------------------|-----|------|---------------|----------------|

- `tables.affinity`:

  | id | input_file | affinity_pred_value | affinity_probability_binary |
  |----|------------|---------------------|----------------------------|

**Example**:
```python
boltz_apo = Boltz2(proteins=lmpnn
boltz_holo = Boltz2(
    proteins=lmpnn,
    ligands="CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin SMILES
    msas=boltz_apo,  # Pass entire ToolOutput
    affinity=True
), env="Boltz2Env")
```

---

### RF3

⚠️ **UNDER DEVELOPMENT** - This tool is still being tested and refined.

Predicts biomolecular structures using RoseTTAFold3. Supports protein-only and protein-ligand complex prediction with batch processing capabilities.

**Important**: RF3 requires MSAs to be provided. You can obtain MSAs either by:
1. Running Boltz2 first on your proteins (which uses public MMseqs2 server automatically)
2. Using our MMseqs2 implementation to generate MSAs separately

**Environment**: `modelforge`

**Installation**:
The official RF3 repository uses `uv` for installation, but for consistency with BioPipelines we use conda. Run the following in your data folder:
```bash
cd /home/$USER/data
git clone https://github.com/RosettaCommons/modelforge.git
cd modelforge
conda create -n modelforge python=3.12
conda activate modelforge
pip install -e .
```

**Parameters**:
- `proteins`: Union[str, List[str], ToolOutput, StandardizedOutput] (required) - Protein sequences
- `ligands`: Union[str, ToolOutput, StandardizedOutput, None] = None - Ligand SMILES string or compound library ToolOutput
- `msas`: Optional[Union[str, ToolOutput]] = None - Pre-computed MSA files (optional)
- `output_format`: str = "pdb" - Output format ("pdb" or "cif")
- `checkpoint_path`: Optional[str] = None - Path to RF3 checkpoint file
- `early_stopping_plddt`: Optional[float] = None - pLDDT threshold for early stopping
- `use_templates`: bool = False - Enable template-based prediction

**Note**: The `num_models` parameter is not currently supported due to configuration override limitations in RF3.

**Outputs**:
- `structures`: List of predicted structure files
- `tables.structures`:

  | id | model_id | file_path | plddt_score |
  |----|----------|-----------|-------------|

- `tables.confidence`:

  | id | model_id | plddt_score | ptm_score |
  |----|----------|-------------|-----------|

**Example**:
```python
# Apo prediction
rf3_apo = RF3(
    proteins=sequences


# Protein-ligand complex prediction
rf3_holo = RF3(
    proteins=sequences,
    ligands="CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin SMILES
    early_stopping_plddt=85.0


# Batch prediction with compound library
compounds = CompoundLibrary({
    "AspA": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "AspB": "CC(=O)OC1=CC=CC=C1C(=O)O"
}
rf3_batch = RF3(
    proteins=protein_sequences,
    ligands=compounds

```

---

### OnionNet

⚠️ **UNDER DEVELOPMENT** - This tool is still being tested and refined.

Predicts protein-ligand binding affinities from complex structures using OnionNet CNN-based model with rotation-free element-pair-specific contacts.

**Environment**: `OnionNetEnv`

**Installation**:
```bash
cd /home/$USER/data
git lfs clone https://github.com/zhenglz/onionnet.git
cd onionnet
conda env create -f onet_env.yaml
```

**Parameters**:
- `structures`: Union[str, ToolOutput, StandardizedOutput] (required) - Protein-ligand complex structures
- `model_weights`: Optional[str] = None - Path to model weights file (.h5)
- `scaler_model`: Optional[str] = None - Path to scaler model file
- `output_format`: str = "csv" - Output format ("csv" or "json")

**Outputs**:
- `tables.affinities`:

  | id | structure_path | predicted_affinity_pKa |
  |----|----------------|------------------------|

**Example**:
```python
# Predict affinities from Boltz2 structures
affinity = OnionNet(
    structures=boltz_output,
    model_weights="/path/to/weights.h5",
    scaler_model="/path/to/scaler.model"

```

---

### OnionNet2

⚠️ **UNDER DEVELOPMENT** - This tool is still being tested and refined.

Predicts protein-ligand binding affinities using OnionNet-2, an improved version with higher accuracy and lower computational cost. Uses residue-atom contacting shells in CNN architecture.

**Environment**: `OnionNet2Env`

**Installation**:
```bash
cd /home/$USER/data
git clone https://github.com/zchwang/OnionNet-2.git
cd OnionNet-2
conda create -n OnionNet2Env python=3.8
conda activate OnionNet2Env
pip install tensorflow==2.3 pandas==1.3.4 scikit-learn==0.22.1 numpy==1.18.5 scipy==1.4.1
```

**Parameters**:
- `structures`: Union[str, ToolOutput, StandardizedOutput] (required) - Protein-ligand complex structures
- `model_path`: Optional[str] = None - Path to trained model file
- `scaler_path`: Optional[str] = None - Path to scaler file
- `shells`: int = 62 - Number of contacting shells
- `output_format`: str = "csv" - Output format ("csv" or "json")

**Outputs**:
- `tables.affinities`:

  | id | structure_path | predicted_affinity_pKa |
  |----|----------------|------------------------|

**Example**:
```python
# Predict affinities with OnionNet-2
affinity = OnionNet2(
    structures=boltz_output,
    model_path="/path/to/model.h5",
    scaler_path="/path/to/scaler.pkl",
    shells=62

```

---

## Analysis

### DynamicBind

⚠️ **UNDER DEVELOPMENT** - Installation failed.

Predicts ligand-specific protein-ligand complex structures using equivariant diffusion models. Recovers ligand-specific conformations from unbound protein structures.

**Key Features:**
- Predicts protein-ligand binding conformations
- Handles multiple ligands per protein
- Generates multiple poses with affinity predictions
- Outputs relaxed structures with confidence scores

**Installation:**

Requires two conda environments.

```bash
cd /home/$USER/data
git clone https://github.com/luwei0917/DynamicBind.git
cd DynamicBind
wget https://zenodo.org/records/10183369/files/workdir.zip
unzip workdir.zip

# Request additional resources during installation
srun --mem=32G --cpus-per-task=8 --time=02:00:00 --pty bash

conda create -n dynamicbind python=3.9 -y
conda activate dynamicbind
pip install torch==1.13.1+cu117 torchvision==0.14.1+cu117 torchaudio==0.13.1 --extra-index-url https://download.pytorch.org/whl/cu117

# Setup LD_LIBRARY_PATH inside the environment
mkdir -p $CONDA_PREFIX/etc/conda/activate.d
echo 'export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH' > \
     $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

# Continue installation
conda install -c conda-forge rdkit biopython scipy pandas pyyaml networkx -y
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
conda create --name relax python=3.8
conda activate relax
conda install -c conda-forge openmm pdbfixer biopython openmmforcefields openff-toolkit ambertools=22 compilers -y
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
# Basic usage - single protein, SMILES string
db = DynamicBind(
    proteins="target.pdb",
    ligands="CCO",  # SMILES string
    samples_per_complex=10


# Multiple proteins from PDBs folder
db = DynamicBind(
    proteins=["protein1.pdb", "protein2.pdb"],
    ligands="CC(=O)O"


# Chain with other tools - all structures from tool output
rfdaa = RFdiffusionAllAtom(ligand="ZIT", contigs="A1-100", num_designs=5
compounds = CompoundLibrary(smiles_list=["CCO", "CC(=O)O"]
db = DynamicBind(proteins=rfdaa, ligands=compounds
```

**Outputs:**
- Structures: SDF files with pattern `rank{N}_ligand_lddt{score}_affinity{score}_relaxed.sdf`
- Tables:
  - `dynamicbind_results.csv`: columns `id`, `ligand_id`, `structure`, `affinity`, `lddt`, `rank`
  - `compounds.csv`: compounds used (available as `output.compounds`)

---

### ResidueAtomDistance

Measures distances between specific atoms and residues in structures. Useful for tracking ligand-protein interactions or structural features.

**Environment**: `ProteinEnv`

**Parameters**:
- `structures`: Union[ToolOutput, StandardizedOutput] (required) - Input structures
- `atom`: str (required) - Atom selection (e.g., 'LIG.Cl', 'name CA', 'A10.CA')
- `residue`: str (required) - Residue selection (e.g., 'D in IGDWG', '145', 'resn ALA')
- `method`: str = "min" - Distance calculation method (min, max, mean, closest)
- `metric_name`: str = None - Custom name for distance column in output (default: "distance")

**Outputs**:
- `tables.analysis`:

  | id | source_structure | {metric_name} |
  |----|------------------|---------------|

**Example**:
```python
distances = ResidueAtomDistance(
    structures=boltz,
    atom="LIG.Cl",
    residue="D in IGDWG",
    method="min",
    metric_name="chlorine_distance"

```

---

### PLIP (Protein-Ligand Interaction Profiler)

⚠️ **UNDER DEVELOPMENT** - This tool is still being tested and refined.

Analyzes protein-ligand interactions using PLIP. Identifies hydrogen bonds, hydrophobic contacts, salt bridges, and other interaction types.

**Environment**: `ProteinEnv`

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
plip = PLIP(
    structures=boltz,
    ligand="LIG",
    create_pymol=True

```

---

### DistanceSelector

Selects protein residues based on proximity to ligands or other reference points. Generates position specifications for downstream design tools. 

**Installation**: It requires an environment containing pandas (e.g. biopipelines).

**Parameters**:
- `structures`: Union[str, List[str], ToolOutput] (required) - Input structures
- `ligand`: str (required) - Ligand identifier for distance reference
- `distance`: float = 5.0 - Distance cutoff in Angstroms
- `reference_type`: str = "ligand" - Type of reference (ligand, atoms, residues)
- `reference_selection`: str = "" - Specific PyMOL selection if not using ligand

**Outputs**:
- `tables.selections`:

  | id | pdb | within | beyond | distance_cutoff | reference_ligand |
  |----|-----|--------|--------|-----------------|------------------|

**Example**:
```python
selector = DistanceSelector(
    structures=boltz,
    ligand="ATP",
    distance=8.0

```

---

### SelectionEditor

Modifies PyMOL-formatted selection strings (e.g., "3-45+58-60") with structure-aware operations. Validates all operations against actual PDB residue numbering and automatically merges overlapping/adjacent ranges.

**Installation**: Requires an environment containing pandas (e.g. biopipelines).

**Parameters**:
- `selection`: tuple (required) - Table column reference (e.g., `tool.tables.structures.designed`)
- `structures`: Optional[Union[ToolOutput, List[str]]] = None - Input structures (auto-detected from selection source if not provided)
- `expand`: int = 0 - Number of residues to add on each side of intervals
- `shrink`: int = 0 - Number of residues to remove from each side of intervals
- `shift`: int = 0 - Number of residues to shift all intervals (+/-)
- `invert`: bool = False - Whether to invert the selection (select complement)

**Outputs**:
- `tables.selections`:

  | id | pdb | {column_name} | original_{column_name} |
  |----|-----|---------------|------------------------|

**Example**:
```python
# Expand binding site selection by 2 residues
distances = DistanceSelector(
    structures=rfdaa,
    ligand="LIG",
    distance=5


expanded = SelectionEditor(
    selection=distances.tables.selections.within,
    expand=2


# Use expanded selection with LigandMPNN
lmpnn = LigandMPNN(
    structures=rfdaa,
    ligand="LIG",
    redesigned=expanded.tables.selections.within


# Invert selection to get everything except binding site
fixed_region = SelectionEditor(
    selection=distances.tables.selections.within,
    invert=True


# Shrink selection
tighter = SelectionEditor(
    selection=distances.tables.selections.within,
    shrink=1

```

**Key Features**:
- **Structure-aware**: Validates against actual PDB residue numbers (handles gaps, non-sequential numbering)
- **Range merging**: Automatically merges overlapping/adjacent intervals (e.g., "1-4+6-10" with expand=1 becomes "1-11")
- **Auto-detection**: Structures parameter is optional; can auto-detect from selection source

---

### ConformationalChange

⚠️ **UNDER DEVELOPMENT** - This tool is still being tested and refined.

Quantifies structural changes between reference and target structures. Calculates RMSD and distance metrics for specified regions.

**Environment**: `ProteinEnv`

**Note**: This tool is not fully debugged yet and may require adjustments.

**Parameters**:
- `reference_structures`: Union[str, ToolOutput, StandardizedOutput] (required) - Reference structures
- `target_structures`: Union[ToolOutput, StandardizedOutput] (required) - Target structures to compare
- `selection`: Union[str, ToolOutput] (required) - Region specification (PyMOL selection or table reference)
- `alignment`: str = "align" - Alignment method (align, super, cealign)

**Outputs**:
- `tables.conformational_analysis`:

  | id | reference_structure | target_structure | selection | num_residues | RMSD | max_distance | mean_distance | sum_over_square_root |
  |----|---------------------|------------------|-----------|--------------|------|--------------|---------------|---------------------|

**Example**:
```python
conf_change = ConformationalChange(
    reference_structures=apo_structures,
    target_structures=holo_structures,
    selection="resi 10-50",  # PyMOL selection
    alignment="super"

```

---

### MutationProfiler

Analyzes mutation patterns across sequence sets. Calculates position-specific amino acid frequencies for understanding sequence diversity.

**Installation**: This tool only needs the a small environment:
```bash
conda create -n MutationEnv seaborn matplotlib pandas logomaker scipy
```

**Parameters**:
- `original`: Union[ToolOutput, StandardizedOutput] (required) - Original/reference sequences
- `mutants`: Union[ToolOutput, StandardizedOutput] (required) - Mutant sequences to analyze
- `include_original`: bool = True - Include original sequence in frequency analysis

**Outputs**:
- `tables.profile`:

  | position | original | count | frequency |
  |----------|----------|-------|-----------|

- `tables.mutations`:

  | position | original | A | C | D | E | F | G | H | I | K | L | M | N | P | Q | R | S | T | V | W | Y |
  |----------|----------|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|

- `tables.absolute_frequencies`:

  | position | original | A | C | ... | Y |
  |----------|----------|---|---|-----|---|

- `tables.relative_frequencies`:

  | position | original | A | C | ... | Y |
  |----------|----------|---|---|-----|---|

**Example**:
```python
profiler = MutationProfiler(
    original=template,
    mutants=lmpnn

```

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


**Key Features**:
- **Multi-metric tracking**: Analyze multiple metrics simultaneously (affinity, pLDDT, contacts, etc.)
- **Accumulation**: Pass previous analysis as history to accumulate statistics across cycles
- **Dual scoring**: Both delta scores (interpretable units) and z-scores (standardized) for flexibility
- **MutationComposer-compatible**: mutation_deltas and mutation_zscores tables use same format as MutationProfiler frequencies
- **Confidence filtering**: min_observations threshold ensures reliable mutation scores
```
---

### ProteinLigandContacts
**Environment**: `ProteinEnv`

Analyzes contacts between selected protein regions and ligands. For each selected residue, calculates the minimum distance to any ligand atom. Returns contact counts and distance statistics.

**Parameters**:
- `structures`: Union[ToolOutput, StandardizedOutput] (required) - Input structures
- `selections`: Union[str, ToolOutput] = None - Protein region selections (string format: '10-20+30-40', table reference, or None for all protein)
- `ligand`: str (required) - Ligand residue name (3-letter code, e.g., 'LIG', 'ATP', 'GDP')
- `contact_threshold`: float = 5.0 - Distance threshold for counting contacts (Angstroms)
- `contact_metric_name`: str = None - Custom name for contact count column (default: "contacts")

**Outputs**:
- `tables.contact_analysis`:

  | id | source_structure | selections | ligand | contacts | min_distance | max_distance | mean_distance | sum_distances_sqrt_normalized |
  |----|------------------|------------|--------|----------|--------------|--------------|---------------|-------------------------------|

**Output Columns**:
- `id`: Structure identifier
- `source_structure`: Path to input structure file
- `selections`: Protein residues analyzed (e.g., '10-20+30-40' or 'all_protein')
- `ligand`: Ligand residue name
- `contacts`: Number of residues within contact_threshold distance
- `min_distance`: Minimum distance from any selected residue to ligand (Å)
- `max_distance`: Maximum distance from any selected residue to ligand (Å)
- `mean_distance`: Mean distance from selected residues to ligand (Å)
- `sum_distances_sqrt_normalized`: Sum of distances divided by √(number of residues)

**Example**:
```python
# Analyze contacts with specific protein regions
contacts = ProteinLigandContacts(
    structures=rfdaa,
    selections=rfdaa.tables.structures.designed,
    ligand="LIG",
    contact_threshold=5.0


# Use fixed selection for all structures
contacts = ProteinLigandContacts(
    structures=boltz,
    selections='10-20+30-40',
    ligand="ATP",
    contact_threshold=4.0


# Analyze all protein residues
contacts = ProteinLigandContacts(
    structures=boltz,
    ligand="GDP"

```

---

### PoseDistance

Measures ligand pose distance between reference holo structure and sample structures. Calculates RMSD and geometric metrics to quantify how well designed structures reproduce known binding poses.

**Environment**: `ProteinEnv`

**Parameters**:
- `reference_structure`: Union[str, ToolOutput, StandardizedOutput] (required) - Reference holo structure (e.g., XRC structure)
- `sample_structures`: Union[str, List[str], ToolOutput, StandardizedOutput] (required) - Designed/predicted structures to compare
- `reference_ligand`: str (required) - Ligand residue name in reference structure (e.g., 'LIG', 'ATP')
- `sample_ligand`: Optional[str] = None - Ligand residue name in sample structures (default: same as reference_ligand)
- `alignment_selection`: str = "protein" - PyMOL selection for protein alignment (e.g., "chain A", "backbone")
- `calculate_centroid`: bool = True - Calculate ligand centroid distance
- `calculate_orientation`: bool = False - Calculate orientation angle difference

**Outputs**:
- `tables.analysis`:

  | id | target_structure | reference_structure | ligand_rmsd | centroid_distance | alignment_rmsd | num_ligand_atoms | alignment_selection |
  |----|------------------|---------------------|-------------|-------------------|----------------|------------------|---------------------|

**Output Columns**:
- `ligand_rmsd`: RMSD between ligand poses after protein alignment (Å)
- `centroid_distance`: Distance between ligand centroids (Å)
- `alignment_rmsd`: RMSD of protein alignment (Å)
- `num_ligand_atoms`: Number of atoms in ligand

**Example**:
```python
# Compare designed structures to XRC reference
xrc = PDB(pdbs="4ufc", ids="reference"
designed = Boltz2(proteins=sequences, ligands="CCO"

pose_analysis = PoseDistance(
    reference_structure=xrc,
    sample_structures=designed,
    reference_ligand="ATP",
    sample_ligand="LIG",
    alignment_selection="chain A and backbone"


# Filter structures with RMSD < 2.0 Å
good_poses = Filter(
    data=pose_analysis.tables.analysis,
    pool=designed,
    expression="ligand_rmsd < 2.0"

```

---

## Data Management

### Filter

Filters structures or sequences based on metric criteria. Uses pandas query expressions to select items meeting specified thresholds.

**Environment**: `ProteinEnv`

**Parameters**:
- `data`: Union[ToolOutput, StandardizedOutput] (required) - Table input to filter
- `pool`: Union[ToolOutput, StandardizedOutput] = None - Structure/sequence pool for copying filtered items
- `expression`: str (required) - Pandas query-style filter expression (e.g., "distance < 3.5 and confidence > 0.8")
- `max_items`: Optional[int] = None - Maximum items to keep after filtering
- `sort_by`: Optional[str] = None - Column name to sort by before applying max_items
- `sort_ascending`: bool = True - Sort order (True = ascending, False = descending)

**Outputs**:
- Filtered pool with same structure as input
- All the tables of the upstream tool given as input will be copied and fitered based on the expression, and maintain the same table name. For example, after filtering an output of MergeTables, you can access filtered.tables.merged.
- `tables.missing`:

  | id | structure | msa |
  |----|-----------|-----|

**Example**:
```python
filtered = Filter(
    data=distances.tables.analysis,
    pool=boltz,
    expression="distance < 3.5 and confidence_score > 0.85",
    max_items=10,
    sort_by="distance"

```

---

### Rank

Ranks entries based on a metric (column or computed expression), renames IDs to sequential format, and optionally copies structures/compounds in ranked order. Useful for generating ranked lists with standardized naming.

**Environment**: `ProteinEnv`

**Parameters**:
- `data`: Union[ToolOutput, StandardizedOutput] (required) - Table input to rank
- `pool`: Union[ToolOutput, StandardizedOutput] = None - Structure/sequence pool for copying ranked items
- `metric`: str (required) - Column name or expression for ranking (e.g., "pLDDT" or "0.8*pLDDT+0.2*affinity")
- `ascending`: bool = False - Sort order (False = descending/best first, True = ascending)
- `prefix`: str = "rank" - Prefix for renamed IDs (e.g., "rank" produces rank_1, rank_2, ...)
- `top`: Optional[int] = None - Limit to top N entries after ranking

**Outputs**:
- Ranked pool with same structure as input (when pool provided)
- `tables.ranked`:

  | id | source_id | metric | {variable_columns} | {original_columns} |
  |----|-----------|--------|-------------------|-------------------|

**Output Columns**:
- `id`: Renamed IDs (e.g., rank_1, rank_2, ...)
- `source_id`: Original IDs
- `metric`: Computed metric column (if expression used)
- Individual variable columns if metric is an expression with multiple variables (e.g., pLDDT, affinity)
- All other original columns preserved

**Metric Types**:
- **Column reference**: Simple column name (e.g., "pLDDT", "affinity")
- **Expression**: Computed metric using pandas eval syntax (e.g., "0.8*pLDDT + 0.2*affinity", "pLDDT - 2*rmsd")

**Example**:
```python
# Rank by single column
ranked = Rank(
    data=analysis.tables.merged,
    pool=boltz,
    metric="pLDDT",
    ascending=False,  # Higher is better
    prefix="model",
    top=10
)

# Rank by computed expression with pool mode
ranked = Rank(
    data=merged.tables.merged,
    pool=boltz,
    metric="0.8*pLDDT + 0.2*binding_affinity",
    prefix="design",
    top=20
)

# Rank with ascending order (lower is better)
ranked = Rank(
    data=distances.tables.analysis,
    pool=structures,
    metric="distance",
    ascending=True,  # Lower is better
    prefix="candidate"
)
```

---

### SelectBest

Selects the single best structure or sequence based on optimization criteria. Supports single or multi-objective optimization with configurable weights.

**Environment**: `ProteinEnv`

**Parameters**:
- `pool`: Union[ToolOutput, StandardizedOutput, List[Union[ToolOutput, StandardizedOutput]]] (required) - Single or list of tool outputs to select from
- `tables`: Union[List[Union[ToolOutput, StandardizedOutput, TableInfo, str]], List[str]] (required) - Tables to evaluate for selection
- `metric`: str (required) - Primary metric to optimize
- `mode`: str = "max" - Optimization direction ("max" or "min")
- `weights`: Optional[Dict[str, float]] = None - Dictionary of {metric_name: weight} for multi-metric selection
- `tie_breaker`: str = "first" - How to break ties ("first", "random", or metric name)
- `composite_function`: str = "weighted_sum" - How to combine metrics (weighted_sum, product, min, max)
- `name`: str = "best" - Name for output structure file

**Outputs**:
- Single best structure/sequence with same format as input pool

**Example**:
```python
best = SelectBest(
    pool=boltz,
    tables=[distances.tables.analysis],
    metric="distance",
    mode="min"


# Multi-objective selection
best_multi = SelectBest(
    pool=boltz,
    tables=[analysis.tables.merged],
    metric="composite_score",
    weights={"binding_affinity": 0.6, "pLDDT": 0.4},
    mode="max"

```

---

### RemoveDuplicates

Removes duplicate structures or sequences from a pool. Supports deduplication by sequence, structure similarity, or ID matching.

**Environment**: `ProteinEnv`

**Parameters**:
- `pool`: Union[ToolOutput, StandardizedOutput] (required) - Items to deduplicate
- `history`: Optional[Union[ToolOutput, StandardizedOutput, List]] = None - Previous tables for cross-cycle deduplication
- `compare`: str = "sequence" - Comparison method (sequence, structure, id)
- `similarity_threshold`: float = 1.0 - Similarity threshold for structure comparison (1.0 = exact match)

**Outputs**:
- Deduplicated pool with same structure as input
- `tables.removed`:

  | id | reason |
  |----|--------|

**Example**:
```python
unique = RemoveDuplicates(
    pool=lmpnn,
    compare="sequence"

```

---

### MergeTables

Combines multiple tables by joining on a common key column. Enables integration of metrics from different analysis tools.

**Environment**: `ProteinEnv`

**Parameters**:
- `tables`: List[Union[ToolOutput, StandardizedOutput, TableInfo, str]] (required) - List of tables to merge
- `key`: str = "id" - Join column name
- `prefixes`: Optional[List[str]] = None - Prefixes for columns from each table
- `suffixes`: Optional[List[str]] = None - Suffixes for columns from each table
- `how`: str = "inner" - Join type (inner, outer, left, right)
- `calculate`: Optional[Dict[str, str]] = None - Derived column expressions {new_col: expression}

**Outputs**:
- `tables.merged`: Combined table with columns from all inputs

**Example**:
```python
merged = MergeTables(
    tables=[distances.tables.analysis, plip.tables.interactions],
    prefixes=["dist_", "plip_"],
    key="id",
    calculate={"score": "dist_distance + plip_energy"}

```

---

### ConcatenateTables

Stacks multiple tables vertically (row-wise). Useful for combining results from multiple cycles or parallel runs.

**Environment**: `ProteinEnv`

**Parameters**:
- `tables`: List[Union[ToolOutput, StandardizedOutput, TableInfo, str]] (required) - List of tables to concatenate
- `fill`: str = "N/A" - Value for missing columns
- `ignore_index`: bool = True - Reset index in concatenated output

**Outputs**:
- `tables.concatenated`: Row-wise concatenation of all input tables

**Example**:
```python
concat = ConcatenateTables(
    tables=[cycle1_results, cycle2_results, cycle3_results],
    fill="N/A"

```

---

### SliceTable

Extracts a subset of rows and/or columns from a table. Enables data sampling and column selection.

**Environment**: `ProteinEnv`

**Parameters**:
- `table`: Union[ToolOutput, StandardizedOutput, TableInfo, str] (required) - Input table to slice
- `start`: int = 0 - Starting row index
- `end`: Optional[int] = None - Ending row index (None = to end)
- `step`: int = 1 - Step size for slicing
- `columns`: Optional[List[str]] = None - Specific columns to keep (None = all columns)

**Outputs**:
- `tables.sliced`: Sliced table

**Example**:
```python
sliced = SliceTable(
    table=results.tables.analysis,
    start=0,
    end=100,
    columns=["id", "distance", "confidence"]

```

---

### ExtractMetrics

Extracts and aggregates specific metrics from tables. Supports grouping and various aggregation functions for data summarization.

**Environment**: `ProteinEnv`

**Parameters**:
- `tables`: List[Union[ToolOutput, StandardizedOutput, TableInfo, str]] (required) - Input tables
- `metrics`: List[str] (required) - Metric column names to extract
- `group_by`: Optional[str] = None - Column to group by for aggregation
- `aggregation`: str = "mean" - Aggregation function (mean, median, min, max, sum, std)
- `pivot`: bool = False - Pivot metrics to columns

**Outputs**:
- `tables.extracted`: Extracted metrics table

**Example**:
```python
metrics = ExtractMetrics(
    tables=[boltz.tables.confidence],
    metrics=["complex_plddt", "ptm"],
    group_by="input_file",
    aggregation="mean"

```

---

### AverageByTable

Computes averages of metrics grouped by a specified column. Useful for summarizing results across multiple structures or cycles.

**Environment**: `ProteinEnv`

**Parameters**:
- `tables`: List[Union[ToolOutput, StandardizedOutput, TableInfo, str]] (required) - Input tables
- `group_by`: str (required) - Column to group by
- `metrics`: List[str] (required) - Metric columns to average
- `weights`: Optional[Dict[str, float]] = None - Weights for each metric

**Outputs**:
- `tables.averaged`: Averaged metrics by group

**Example**:
```python
averaged = AverageByTable(
    tables=[cycle1.tables.analysis, cycle2.tables.analysis],
    group_by="structure_id",
    metrics=["distance", "confidence"]

```

---

## Utilities

### LoadOutput

Loads previously saved tool outputs for reuse in new pipelines. Enables incremental development and filtering of existing results at pipeline runtime.

**Environment**: `ProteinEnv`

**Parameters**:
- `output_json`: str (required) - Path to tool output JSON file (in ToolOutputs folder)
- `filter`: Optional[str] = None - Pandas query-style filter expression
- `validate_files`: bool = True - Check file existence when loading

**Outputs**:
- Same structure as original tool that created the output

**Example**:
```python
previous_boltz = LoadOutput(
    output_json="/path/to/job/ToolOutputs/003_Boltz2.json",
    filter="confidence_score > 0.8"

```

**PracticalTips**:
- Instead of copy-pasting ligand SMILES across pipelines, you can create a compound library, and then load the smiles passing an id filter:
```python
# Pipeline 1 to store the library
# imports, pipeline instantiation, ...
CompoundLibrary({
  "Compound1": "CCNCNNCC(=O)C",
  "Compound2": "CCNCNNCC(=O)C",
  ...
}
# submit

# Pipeline 2 running calculatios with one of the compounds
# imports, pipeline instantiation, ...
compound1 = LoadOutput(
    output_json="/path/to/job/ToolOutputs/001_CompoundLibrary.json",
    filter='id == "Compound1"' #quotes are important for proper pandas query here: x is a column name; "x" is a string.

boltz = Boltz2(proteins=HaloTag,
                            ligands=compound1
#submit

---

### MMseqs2

Generates multiple sequence alignments (MSAs) for protein sequences. Used for improving structure prediction quality by providing evolutionary information.

**Environment**: `ProteinEnv`

**Parameters**:
- `sequences`: Union[str, List[str], ToolOutput, StandardizedOutput] (required) - Input sequences
- `output_format`: str = "csv" - Output format (csv, a3m)
- `timeout`: int = 3600 - Timeout in seconds for server response

**Outputs**:
- `tables.msas`:

  | id | sequence_id | sequence | msa_file |
  |----|-------------|----------|----------|

**Example**:
```python
msas = MMseqs2(
    sequences=lmpnn,
    timeout=7200

```

---

### MMseqs2Server

Runs an MMseqs2 server for local MSA generation. Automatically started by MMseqs2 client when needed; manual setup typically not required.

**Environment**: None (doesn't require conda)

**Parameters**:
- `port`: int = 8000 - Server port
- `host`: str = "0.0.0.0" - Server host
- `workers`: int = 4 - Number of worker processes

**Note**: This server tool must be run separately to provide MMseqs2 as a service. However, the MMseqs2 client automatically runs it when needed, so manual server setup is typically not required.

---

### CompoundLibrary

Creates and manages chemical compound libraries. Supports combinatorial SMILES expansion and optional generation of covalent binding files.

**Environment**: `ProteinEnv`

**Parameters**:
- `library`: Union[str, Dict[str, Union[str, List[str]]]] (required) - Dictionary with expansion keys or path to CSV file
- `primary_key`: Optional[str] = None - Root key for expansion when library is dictionary
- `covalent`: bool = False - Generate CCD/PKL files for covalent ligand binding
- `validate_smiles`: bool = True - Validate SMILES strings during expansion
- `conformer_method`: str = "UFF" - Method for conformer generation (UFF, OpenFF, DFT)

**Outputs**:
- `compounds`: CSV file with compound library
- `tables.compounds`:

  | id | format | smiles | ccd | {branching_keys} |
  |----|--------|--------|-----|------------------|

**Examples**:
```python
# With primary_key for combinatorial expansion
library = CompoundLibrary(
    library={
        "scaffold": "<linker><fluorophore>",
        "linker": ["CCOCC", "CCOCCOCC"],
        "fluorophore": ["c1ccc(N)cc1"]
    },
    primary_key="scaffold"


# Without primary_key - direct SMILES list
library = CompoundLibrary(
    library={
        "compounds": ["CCO", "CCCO", "CCCCO", "c1ccccc1"]
    }

```

---

### PDB

Fetches protein structures with automatic fallback: checks local folders first, then downloads from RCSB PDB if not found locally. Downloads are automatically cached in PDBs/ folder for reuse.

**Environment**: `ProteinEnv`

**Parameters**:
- `pdbs`: Union[str, List[str]] (required) - PDB IDs to fetch (e.g., "4ufc" or ["4ufc","1abc"])
- `ids`: Optional[Union[str, List[str]]] - Custom IDs for renaming (defaults to pdbs)
- `format`: str = "pdb" - File format ("pdb" or "cif")
- `local_folder`: Optional[str] = None - Custom local folder to check first (before PDBs/)
- `biological_assembly`: bool = False - Download biological assembly from RCSB
- `remove_waters`: bool = True - Remove water molecules

**Fetch Priority**:
For each PDB ID, searches in order:
1. `local_folder` (if parameter provided)
2. `./PDBs/` folder in repository
3. Download from RCSB PDB (cached to PDBs/ for reuse)

**Outputs**:
- `structures`: List of structure files
- `tables.structures`:

  | id | pdb_id | file_path | format | file_size | source |
  |----|--------|-----------|--------|-----------|--------|

- `tables.sequences`:

  | id | sequence |
  |----|----------|

- `tables.failed`:

  | pdb_id | error_message | source | attempted_path |
  |--------|---------------|--------|----------------|

**Examples**:
```python
# Automatic fallback: checks PDBs/, then downloads
pdb = PDB(
    pdbs=["4ufc", "1abc"],
    ids=["POI1", "POI2"]


# Check custom folder first, then PDBs/, then download
pdb = PDB(
    pdbs=["my_structure"],
    local_folder="/path/to/my/pdbs"


# Download biological assembly
pdb = PDB(
    pdbs="4ufc",
    ids="MyProtein",
    biological_assembly=True,
    remove_waters=False

```

---

### PyMOL

Creates PyMOL visualizations and session files for structural analysis. Supports alignment, coloring schemes, and ray-traced image generation.

**Environment**: `ProteinEnv`

**Note**: This tool is not fully debugged yet and may require adjustments.

**Parameters**:
- `structures`: Union[str, List[str], ToolOutput, StandardizedOutput] (required) - Structures to visualize
- `reference_structure`: Optional[Union[str, ToolOutput]] = None - Reference for alignment
- `color_by`: str = "chain" - Coloring scheme (chain, b_factor, element, ss, spectrum)
- `alignment`: str = "align" - Alignment method (align, super, cealign)
- `ray_trace`: bool = False - Generate ray-traced images
- `image_size`: Tuple[int, int] = (1920, 1080) - Image dimensions if ray tracing

**Outputs**:
- PyMOL session files (.pse)
- Images (if ray_trace=True)

**Example**:
```python
pymol = PyMOL(
    structures=best,
    reference_structure=template,
    color_by="b_factor",
    alignment="super"

```

---
