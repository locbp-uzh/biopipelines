# Structure Prediction

[← Back to Tool Reference](../ToolReference.md)

---

### AlphaFold

Predicts protein structures from amino acid sequences using AlphaFold2. Generates high-confidence 3D models with optional relaxation. Supports single-sequence prediction mode for fast predictions without MSA generation.

**Resources**: GPU. H100 NVL is not compatible.

**Environment**: `localcolabfold` and `biopipelines`

**Installation**:
```bash
cd data
wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh
bash install_colabbatch_linux.sh
rm install_colabbatch_linux.sh
```

**Parameters**:
- `sequences`: Union[str, List[str], ToolOutput, Dict[str, Any]] (required) - Input sequences or dict with sequences
- `tables`: Optional[List[str]] = None - Input table files
- `name`: str = "" - Job name
- `num_relax`: int = 0 - Number of best models to relax with AMBER
- `num_recycle`: int = 3 - Number of recycling iterations
- `rand_seed`: int = 0 - Random seed (0 = random)
- `msa_mode`: Optional[str] = None - MSA generation mode:
  - `"mmseqs2_uniref_env"` - MMseqs2 with UniRef + environmental databases (default)
  - `"mmseqs2_uniref"` - MMseqs2 with UniRef only
  - `"single_sequence"` - No MSA generation (fast single-sequence prediction)

**Streams**: `structures`, `msas`

**Tables**:
- `structures`:

  | id | source_id | sequence |
  |----|-----------|----------|

- `confidence`:

  | id | structure | plddt | max_pae | ptm |
  |----|-----------|-------|---------|-----|

- `msas` (only when MSAs are generated):

  | id | sequence_id | sequence | msa_file |
  |----|-------------|----------|----------|

**Example**:
```python
from biopipelines.alphafold import AlphaFold

# Standard prediction with MSA
af = AlphaFold(
    sequences=lmpnn,
    num_relax=1,
    num_recycle=5
)

# Fast single-sequence prediction (no MSA)
af_fast = AlphaFold(
    sequences=lmpnn,
    msa_mode="single_sequence",
    num_recycle=3
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
- `recycling_steps`: Optional[int] = None - Number of recycling steps (default: model-specific)
- `diffusion_samples`: Optional[int] = None - Number of diffusion samples (default: model-specific)
- `use_potentials`: bool = False - Enable external potentials
- `template`: Optional[str] = None - Path to PDB template file for structure guidance
- `template_chain_ids`: Optional[List[str]] = None - Chain IDs to apply template to (e.g., ["A", "B"])
- `template_force`: bool = True - Force template usage
- `template_threshold`: float = 5.0 - RMSD threshold for template matching
- `pocket_residues`: Optional[List[int]] = None - Residue positions defining binding pocket (e.g., [50, 51, 52])
- `pocket_max_distance`: float = 7.0 - Maximum distance for pocket constraint
- `pocket_force`: bool = True - Force pocket constraint
- `glycosylation`: Optional[Dict[str, List[int]]] = None - N-glycosylation sites per chain (e.g., {"A": [164]})
- `covalent_linkage`: Optional[Dict[str, Any]] = None - Covalent attachment specification (see examples)

**Streams**: `structures`

**Tables**:
- `confidence`:

  | id | sequences.id | compounds.id | input_file | confidence_score | ptm | iptm | complex_plddt | complex_iplddt |
  |----|--------------|--------------|------------|------------------|-----|------|---------------|----------------|

- `affinity`:

  | id | sequences.id | compounds.id | input_file | affinity_pred_value | affinity_probability_binary |
  |----|--------------|--------------|------------|---------------------|----------------------------|

  Provenance columns (`sequences.id`, `compounds.id`) track which protein and ligand produced each row, enabling filtering and joins without parsing the ID string.

**Example**:
```python
from biopipelines.boltz2 import Boltz2

# Basic apo and holo prediction
boltz_apo = Boltz2(proteins=lmpnn)
boltz_holo = Boltz2(
    proteins=lmpnn,
    ligands="CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin SMILES
    msas=boltz_apo,  # Pass entire ToolOutput
    affinity=True
)

# With template guidance
boltz_template = Boltz2(
    proteins=lmpnn,
    ligands=compounds,
    template="reference.pdb",
    template_chain_ids=["A"]
)

# With pocket constraint
boltz_pocket = Boltz2(
    proteins=lmpnn,
    ligands=compounds,
    pocket_residues=[50, 51, 52, 80, 81, 82],
    pocket_max_distance=7.0
)

# With N-glycosylation (adds NAG at Asn-164)
boltz_glyco = Boltz2(
    proteins=lmpnn,
    ligands=compounds,
    glycosylation={"A": [164]}
)

# With covalent ligand attachment (e.g., to Cys-50)
boltz_covalent = Boltz2(
    proteins=lmpnn,
    ligands=covalent_inhibitor,
    covalent_linkage={
        "chain": "A",
        "position": 50,
        "protein_atom": "SG",  # Cysteine sulfur
        "ligand_atom": "C1"    # Ligand attachment atom
    }
)

# Multi-chain complex from SplitChains
split = SplitChains(sequences=fused, split_positions=[200])
boltz_complex = Boltz2(proteins=split, ligands=compounds)
```

---
