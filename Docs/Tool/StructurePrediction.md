# Structure Prediction

[← Back to Tool Reference](../ToolReference.md)

---

### AlphaFold

Predicts protein structures from amino acid sequences using AlphaFold2. Generates high-confidence 3D models with optional relaxation.

**Resources**: GPU. H100 NVL is not compatible.

**Environment**: `biopipelines`

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
from PipelineScripts.alphafold import AlphaFold

af = AlphaFold(
    sequences=lmpnn,
    num_relax=1,
    num_recycle=5
)
```

---

### GhostFold

Database-free protein structure prediction using synthetic MSAs. GhostFold generates structure-aware multiple sequence alignments from single sequences using ProstT5, eliminating the need for large sequence databases while maintaining prediction accuracy. Uses ColabFold for the actual structure prediction.

**References**: https://github.com/brineylab/ghostfold

**Installation**:
```bash
# Clone GhostFold repository
cd ~/data
git clone https://github.com/brineylab/ghostfold
cd ghostfold

# Create conda environment
mamba create -n ghostfold python=3.10
mamba activate ghostfold

# Install dependencies
mamba install pytorch torchvision torchaudio -c pytorch
mamba install transformers sentencepiece -c conda-forge

pip install rich biopython
# Login to HuggingFace for ProstT5 model access
huggingface-cli login
```

Note: GhostFold requires LocalColabFold for structure prediction. Ensure ColabFold is installed and the `colabfold` environment is available.

**Resources**: GPU (A100 or better recommended for ProstT5 + ColabFold)

**Environment**: `ghostfold`

**Parameters**:
- `sequences`: Union[str, List[str], ToolOutput, Dict[str, Any]] (required) - Input sequences (FASTA, CSV, or upstream tool output)
- `name`: str = "" - Job name for output files
- `msa_only`: bool = False - Only generate synthetic MSAs, skip structure prediction
- `num_recycle`: int = 10 - Number of ColabFold recycling iterations
- `num_models`: int = 5 - Number of AlphaFold2 models to use (1-5)
- `num_seeds`: int = 5 - Number of random seeds for prediction diversity
- `subsample`: bool = False - Enable multi-level MSA subsampling (tests max-seq: 16, 32, 64, 128)
- `mask_msa`: Optional[float] = None - MSA masking fraction (0.0-1.0, e.g., 0.15 for 15%)

**Outputs**:

When `msa_only=False` (default - full prediction):
- `structures`: List of predicted PDB files
- `msas`: List of synthetic MSA files (.a3m)
- `tables.structures`:

  | id | sequence |
  |----|----------|

- `tables.confidence`:

  | id | structure | plddt | max_pae | ptm |
  |----|-----------|-------|---------|-----|

- `tables.msas`:

  | id | sequence_id | sequence | msa_file |
  |----|-------------|----------|----------|

When `msa_only=True` (MSA generation only):
- `msas`: List of synthetic MSA files (.a3m)
- `tables.msas`:

  | id | sequence_id | sequence | msa_file |
  |----|-------------|----------|----------|

**Example**:
```python
from PipelineScripts import Pipeline, ProteinMPNN, GhostFold, Resources

# Full prediction pipeline
with Pipeline("MyProject", "GhostFoldTest", "Database-free structure prediction"):
    Resources(gpu="A100", time="4:00:00", memory="32GB")

    pmpnn = ProteinMPNN(pdb="input.pdb", num_sequences=5)
    gf = GhostFold(sequences=pmpnn, num_recycle=10)

# MSA-only mode (for use with other folding tools like Boltz2)
with Pipeline("MyProject", "MSAOnly", "Generate synthetic MSAs"):
    Resources(gpu="A100", time="2:00:00", memory="32GB")

    pmpnn = ProteinMPNN(pdb="input.pdb", num_sequences=5)
    msas = GhostFold(sequences=pmpnn, msa_only=True)

    # Use MSAs with Boltz2
    boltz = Boltz2(proteins=pmpnn, msas=msas)
```

---

### ESMFold

Predicts protein structures using Meta's ESM-2 with ESMFold. Fast single-sequence prediction without requiring MSAs. Models are cached in shared folder for reuse.

**References**: https://github.com/facebookresearch/esm

**Installation**: 
```bash
mamba create -n esmfold python=3.12
mamba activate esmfold
pip install 
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
from PipelineScripts.esmfold import ESMFold

esm = ESMFold(
    sequences=lmpnn,
    num_recycle=4
)
```

---

### Boltz2

Predicts biomolecular complexes including proteins, nucleic acids, and small molecules. State-of-the-art model for protein-ligand and protein-protein complex prediction.

**Installation**:
```bash
mamba create -n Boltz2Env python=3.11
mamba activate Boltz2Env
pip install boltz[cuda] -U
```
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
from PipelineScripts.boltz2 import Boltz2

boltz_apo = Boltz2(proteins=lmpnn)
boltz_holo = Boltz2(
    proteins=lmpnn,
    ligands="CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin SMILES
    msas=boltz_apo,  # Pass entire ToolOutput
    affinity=True
)
```

---

### RF3

**UNDER DEVELOPMENT** - This tool is still being tested and refined.

Predicts biomolecular structures using RoseTTAFold3. Supports protein-only and protein-ligand complex prediction with batch processing capabilities.

**Important**: RF3 requires MSAs to be provided. You can obtain MSAs either by:
1. Running Boltz2 first on your proteins (which uses public MMseqs2 server automatically)
2. Using our MMseqs2 implementation to generate MSAs separately

**Environment**: `modelforge`

**Installation**:
The official RF3 repository uses `uv` for installation, but for consistency with BioPipelines we use mamba. Run the following in your data folder:
```bash
cd /home/$USER/data
git clone https://github.com/RosettaCommons/modelforge.git
cd modelforge
mamba create -n modelforge python=3.12
mamba activate modelforge
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
from PipelineScripts.rf3 import RF3
from PipelineScripts.compound_library import CompoundLibrary

# Apo prediction
rf3_apo = RF3(
    proteins=sequences
)

# Protein-ligand complex prediction
rf3_holo = RF3(
    proteins=sequences,
    ligands="CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin SMILES
    early_stopping_plddt=85.0
)

# Batch prediction with compound library
compounds = CompoundLibrary({
    "AspA": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "AspB": "CC(=O)OC1=CC=CC=C1C(=O)O"
})
rf3_batch = RF3(
    proteins=protein_sequences,
    ligands=compounds
)
```

---

### OnionNet

**UNDER DEVELOPMENT** - This tool is still being tested and refined.

Predicts protein-ligand binding affinities from complex structures using OnionNet CNN-based model with rotation-free element-pair-specific contacts.

**Environment**: `OnionNetEnv`

**Installation**:
```bash
cd /home/$USER/data
git lfs clone https://github.com/zhenglz/onionnet.git
cd onionnet
mamba env create -f onet_env.yaml
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
from PipelineScripts.onion_net import OnionNet

# Predict affinities from Boltz2 structures
affinity = OnionNet(
    structures=boltz_output,
    model_weights="/path/to/weights.h5",
    scaler_model="/path/to/scaler.model"
)
```

---

### OnionNet2

**UNDER DEVELOPMENT** - This tool is still being tested and refined.

Predicts protein-ligand binding affinities using OnionNet-2, an improved version with higher accuracy and lower computational cost. Uses residue-atom contacting shells in CNN architecture.

**Environment**: `OnionNet2Env`

**Installation**:
```bash
cd /home/$USER/data
git clone https://github.com/zchwang/OnionNet-2.git
cd OnionNet-2
mamba create -n OnionNet2Env python=3.8
mamba activate OnionNet2Env
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
from PipelineScripts.onion_net import OnionNet2

# Predict affinities with OnionNet-2
affinity = OnionNet2(
    structures=boltz_output,
    model_path="/path/to/model.h5",
    scaler_path="/path/to/scaler.pkl",
    shells=62
)
```

---
