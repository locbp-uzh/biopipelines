# Structure Generation

[← Back to Tool Reference](../ToolReference.md)

---

### RFdiffusion

Generates novel protein backbone structures using diffusion models. Designs de novo proteins or scaffolds functional motifs into new contexts.

**References**: https://www.nature.com/articles/s41586-023-06415-8

**Installation**: Go to the data directory, then clone the RFdiffusion git repository and download the weights as indicated in https://github.com/RosettaCommons/RFdiffusion. The environment they provide won't work on our cluster, so instead go to the biopipelines directory and run this:
```bash
mamba env create -f Environments/ProteinEnv.yaml
mamba activate ProteinEnv
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
from PipelineScripts.rfdiffusion import RFdiffusion

rfd = RFdiffusion(
    contigs="50-100",
    num_designs=10
)
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
from PipelineScripts.rfdiffusion_allatom import RFdiffusionAllAtom

rfdaa = RFdiffusionAllAtom(
    pdb=template,
    ligand='LIG',
    contigs='10-20,A6-140',
    num_designs=5
)
```

---

### RFdiffusion3

Third-generation all-atom diffusion model for protein design. Operates at the atomic level (4 backbone + 10 sidechain atoms per residue) rather than residue-level, enabling precise design of sidechain interactions with ligands, catalytic residues, DNA/RNA, and symmetric assemblies. Approximately 10× faster than RFdiffusion2 with higher success rates on enzyme design benchmarks.

**References**:
- Paper: https://www.biorxiv.org/content/10.1101/2024.11.13.623358v1
- Github: https://github.com/RosettaCommons/foundry/blob/production/models/rfd3
- Input: https://github.com/RosettaCommons/foundry/blob/production/models/rfd3/docs/input.md

**Installation**:

RFdiffusion3 is part of the foundry framework, which requires Python ≥3.12. Install in a dedicated conda environment:

```bash
srun --pty --time=02:00:00 --mem=16GB bash
mamba create -n foundry python=3.12
mamba activate foundry
pip install "rc-foundry[all]"
#Download model weights. Default is ~/.foundry/checkpoints
mkdir /home/$USER/data/rfdiffusion3
foundry install rfd3 --checkpoint-dir /home/$USER/data/rfdiffusion3
```
**Parameters**:
- `input`: Union[str, ToolOutput, StandardizedOutput] = None - Input PDB structure (optional for de novo design)
- `contig`: str (required) - Contig specification defining protein chain ranges and design regions
  - Format: "A1-50,60-80,B1-100" (chain letters indicate motif from input PDB)
  - "/0" indicates chain break (adds 200aa jump)
  - Numbers without prefix = design new residues
- `length`: str = None - Length range for designed regions (e.g., "50-150")
- `select_unfixed_sequence`: str = None - Residue positions where sequence can change (e.g., "A69-76,A153-154")
- `select_fixed_atoms`: Dict[str, str] = None - Chain ranges to control atomic constraints
- `select_hotspots`: Dict[str, str] = None - Hotspot residues and/or specific atoms (e.g., {"A181": "CG2,CG1", "A407": "NE1,CZ2"})
- `infer_ori_strategy`: str = None - Strategy for orientation inference (e.g., "hotspots")
- `num_designs`: int = 1 - Number of designs to generate
- `num_designs`: int = 1 - Number of models to generate for each design. Roughly equivalent in practice to number of sequences per design in backbone-only design
- `design_startnum`: int = 1 - Starting number for design IDs
- `design_name`: str = None - Base name for output designs
- `out_dir`: str = None - Output directory path
- `symmetric`: bool = False - Enable symmetric design mode
- `symmetry_type`: str = None - Symmetry specification (e.g., "C3", "D2")

**Outputs**:
- `structures`: List of generated PDB files
- `tables.structures`:

  | id | source_id | pdb | contig | length | design_name | status |
  |----|-----------|-----|--------|--------|-------------|--------|

**Example**:
```python
from PipelineScripts.rfdiffusion3 import RFdiffusion3

# De novo protein design
rfd3 = RFdiffusion3(
    contig="80-120",
    num_designs=10
)

# Binder design with hotspots
# You can see the name of each atom in pymol picking mode
# If your ligand pdb doesn't have numbering of atoms (e.g. C, C, ...) you can assign individual ids in pymol with the command rename, and then save the pdb
rfd3 = RFdiffusion3(
    input=target_pdb,
    contig="50-80,A1-100",
    select_hotspots={"A45": "NE,CZ", "A67": "OG"},
    infer_ori_strategy="hotspots",
    num_designs=20
)

# Enzyme active site design
rfd3 = RFdiffusion3(
    input=template_pdb,
    contig="A1-50,40-60,A100-150",
    select_unfixed_sequence="A25-30,A125-130",
    length="90-110",
    num_designs=50
)
```

---

### BoltzGen

Generates protein binders (proteins, peptides, or nanobodies) targeting specified molecules using an end-to-end pipeline that combines diffusion-based backbone generation, inverse folding, structure prediction, and multi-metric filtering.

**References**:
- https://github.com/HannesStark/boltzgen
- https://www.biorxiv.org/content/10.1101/2025.11.20.689494v1.full.pdf

**Installation**: Install BoltzGen via pip with Python ≥3.11:
```bash
mamba create -n boltzgen python=3.11
mamba activate boltzgen
pip install boltzgen
```

**Resources**: GPU needed. Generation of 1000 designs for a small-molecule binder of 140-180 residues from design to affinity takes approx. 12-14 h (A100) or 6-7 h (H100). The analysis step requires 64 GB of memory.

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
Tables are only predicted if the appropriate step is chosen. There are many columns, please refer to the repository for more info.

After analysis you can access:

- `tables.aggregate_metrics`:
- `tables.per_target_metrics`:

After filtering you can access:

- `tables.all_designs_metrics`:
- `tables.final_metrics`:

**Notes**:
- **Residue indexing**: All residue indices start at 1 and use canonical mmcif `label_asym_id`, not `auth_asym_id`
- **Sequence specification**: Use ranges like "80..140" for random length, "15..20AAAA" for random prefix + fixed tail
- **Output structure**: Pipeline generates intermediate designs, inverse-folded structures, and final ranked designs with comprehensive metrics
- **Small molecule binders**: From the preprint:
```text
Our computational design pipeline generated on the order of 10 or 20 thousands designs of length 140 to 180 residues for each target. The initial set was filtered to a highest-confidence subset of 100 designs based on RMSD < 2.5 Å relative to the Boltz-2 refolded models. Within this filtered pool, designs were ranked using a composite metric (interaction score + Boltz score), which effectively integrates predicted structural fidelity with biophysical interaction quality. We slightly modify the weights of Boltz-2 metrics: design_iiptm: 1.1, design_ptm: 1.1, neg_min_design_to_target_pae: 1.1, plip_hbonds_refolded: 2, plip_saltbridge_refolded: 2, and delta_sasa_refolded: 2. To identify essential interactions, we fragmented the chemical groups of rucaparib and calculated the number of hydrogen bonds formed with each fragment. We prioritized designs forming hydrogen bonds with the carboxamide chemical groups, as this interaction is considered essential for specific binding.
```

**Example**:
```python
from PipelineScripts.boltzgen import BoltzGen

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
)

# Access final designs
best_designs = boltzgen.tables.final_metrics
```

**Example with File**:
```python
from PipelineScripts.boltzgen import BoltzGen

# Using external YAML file
boltzgen = BoltzGen(
    design_spec="designs/my_binder.yaml",
    protocol="peptide-anything",
    num_designs=5000,
    budget=50,
    additional_filters=["design_ALA<0.3", "filter_rmsd_design<2.5"]
)
```

---
