# Structure Generation

[← Back to Tool Reference](../tool_reference.md)

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

The design can be specified manually (Option 1, a full YAML/dict) or built automatically from BioPipelines inputs (Option 2). Provide either `design_spec`, or one of the automatic-build inputs.
- `design_spec`: Union[str, Dict] = None - Option 1: YAML configuration string/dict or path to YAML file defining:
  - Target entities (proteins, ligands from .cif/.pdb files)
  - Binder specification (sequence ranges like "80..140")
  - Binding site constraints (binding_types, not_binding)
  - Secondary structure specifications (helix/sheet/loop)
  - Structural constraints (disulfide bonds, covalent connections)
- `ligand`: Union[DataStream, StandardizedOutput] = None - Option 2 (for `protocol="protein-small_molecule"`): compounds stream to build a small-molecule target from.
- `target_structure`: Union[str, DataStream, StandardizedOutput] = None - Option 2 (for `protocol="protein-anything"`): target structure to build the design spec from.
- `binder_spec`: Union[str, Dict] = None - Option 2: binder specification (e.g. sequence range) for automatic building.
- `binding_region`: str = None - Option 2: binding-site region on the target.
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
- `alpha`: float = None - Quality/diversity trade-off (0.0=quality only, 1.0=diversity only; None uses the pipeline default)
- `filter_biased`: bool = True - Remove amino acid composition outliers
- `refolding_rmsd_threshold`: Optional[float] = None - RMSD cutoff against the refolded model for the highest-confidence subset
- `additional_filters`: Optional[List[str]] = None - Hard threshold expressions (e.g., ["design_ALA>0.3"])
- `metrics_override`: Optional[Dict[str, float]] = None - Per-metric ranking weights
- `devices`: Optional[int] = None - Number of GPUs (auto-detect if None)
- `reuse`: Union[DataStream, StandardizedOutput, str, None] = None - Resume an interrupted run (pass the previous run's output or path)
- `steps`: Optional[List[str]] = None - Run only specific pipeline steps
- `cache_dir`: Optional[str] = None - Model download location

**Tables** (only predicted if the appropriate step is chosen; there are many columns, please refer to the repository for more info):

After analysis:
- `aggregate_metrics`
- `per_target_metrics`

After filtering:
- `all_designs_metrics`
- `final_metrics`

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
from biopipelines.boltzgen import BoltzGen

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
from biopipelines.boltzgen import BoltzGen

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

### PocketGen

Designs a protein binding pocket tailored to a chosen target ligand. For each (scaffold, ligand) pair, PocketGen jointly redesigns the pocket residues' **sequence and side-chain atoms** so the pocket accommodates the ligand. One designed pocket per pair (8 samples per pair internally).

The scaffold must have the **ligand bound as HETATM** in the binding site — those coordinates define the pocket PocketGen edits. The ligand is supplied as a compounds stream carrying both the 3-letter `code` (present in the scaffold) and the canonical `smiles` (used as a bond-order template), typically `Ligand("CODE", smiles="...")`. See [The Ligand Contract](../developer_manual.md#the-ligand-contract-compounds--chemistry-structures--coordinates).

**References**: https://github.com/zaixizhang/PocketGen · https://www.nature.com/articles/s42256-024-00920-9

**Resources**: GPU.

**Environment**: `pocketgen`

**Installation**: `PocketGen.install()` clones the repo, creates the env, and downloads the pretrained checkpoint.

**Parameters**:
- `structures`: DataStream | StandardizedOutput (required) — Scaffold backbones with the ligand bound as HETATM.
- `ligand`: DataStream | StandardizedOutput (required) — Compounds stream with the bound ligand's `code` and `smiles`.

**Streams**:
- `structures` — relaxed PDB per pair (designed pocket grafted into the scaffold).
- `sequences` — designed pocket sequence per pair.

**Tables**:
- `sequences`: | id | structures.id | compounds.id | sequence |
- `missing`: | id | removed_by | cause |

**Example**:
```python
from biopipelines.pocketgen import PocketGen
from biopipelines.entities import PDB, Ligand

scaffold = PDB("2p16", convert="pdb")          # ligand bound as HETATM
lig = Ligand("GG2", smiles="...")               # CCD code + canonical SMILES
pg = PocketGen(structures=scaffold, ligand=lig)
```

---

### RFdiffusion

Generates novel protein backbone structures using diffusion models. Designs de novo proteins or scaffolds functional motifs into new contexts.

**References**: https://www.nature.com/articles/s41586-023-06415-8

**Installation**: Go to the data directory, then clone the RFdiffusion git repository and download the weights as indicated in https://github.com/RosettaCommons/RFdiffusion. Create the SE3nv environment from the biopipelines directory:
```bash
mamba env create -f environments/SE3nv.yaml
mamba activate SE3nv
pip install -r environments/SE3nv_pip_requirements.txt
```
If the BioPipelines SE3nv.yaml fails (e.g., due to CUDA version mismatch), try the official environment file from the cloned RFdiffusion repository instead: `mamba env create -f env/SE3nv.yml`. Then install SE3Transformer and RFdiffusion:
```bash
cd env/SE3Transformer
pip install --no-cache-dir -r requirements.txt
python setup.py install
cd ../.. # change into the root directory of the repository
pip install -e . # install the rfdiffusion module from the root of the repository
```

**Parameters**:
- `contigs`: str | (TableInfo, column) (required) - Contig specification (e.g., "A1-100", "10-20,A6-140"). May be a table column reference resolved per-input-PDB at runtime — useful when each input structure needs its own contig.
- `pdb`: Optional[Union[DataStream, StandardizedOutput]] = None - Input PDB template(s) (optional). **When the input stream holds multiple structures, RFdiffusion now runs once per input PDB** (a runtime loop over the input ids), producing `num_designs` designs for each — so an upstream pose ensemble (e.g. PLACER's models) propagates fully into design instead of only the first structure being used.
- `inpaint`: str | (TableInfo, column) = "" - Inpainting specification (also accepts a per-PDB table column reference)
- `inpaint_str`: str | (TableInfo, column) = "" - Secondary-structure inpainting specification (also accepts a per-PDB table column reference)
- `num_designs`: int = 1 - Number of backbone designs to generate per input PDB
- `active_site`: bool = False - Use active site model for small motifs (<30 residues)
- `steps`: int = 50 - Number of diffusion steps
- `partial_steps`: int = 0 - Partial diffusion steps
- `reproducible`: bool = False - Use deterministic sampling
- `design_startnum`: int = 1 - Starting number for design IDs

**Streams**: `structures`

**Tables**:
- `structures`:

  | id | source_id | pdb | fixed | designed | contigs | time | status |
  |----|-----------|-----|-------|----------|---------|------|--------|

**Example**:
```python
from biopipelines.rfdiffusion import RFdiffusion

rfd = RFdiffusion(
    contigs="50-100",
    num_designs=10
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
- `contig`: str | (TableInfo, column) = "" - Contig specification defining protein chain ranges and design regions
  - Format: "A1-50,60-80,B1-100" (chain letters indicate motif from input PDB)
  - "/0" indicates chain break (adds 200aa jump)
  - Numbers without prefix = design new residues
  - May be a table column reference resolved per-input-PDB at runtime (each input structure gets its own contig)
- `length`: Union[str, int] = None - Length range for designed regions (e.g., "50-150")
- `pdb`: Optional[Union[DataStream, StandardizedOutput]] = None - Input PDB structure(s), optional for de novo design. **Multiple input structures are iterated** — one design run per input PDB (`num_designs` each), so a pose ensemble propagates fully into design instead of only the first structure being used.
- `ligand`: Optional[Union[str, DataStream, StandardizedOutput]] = None - Ligand to design around (compounds stream or 3-letter code)
- `num_designs`: int = 1 - Number of designs to generate
- `num_models`: int = 1 - Number of models per design. Roughly equivalent to number of sequences per design in backbone-only design
- `prefix`: str = None - Base name for output designs
- `select_hotspots`: Union[str, Dict[str, str]] = None - Hotspot residues and/or specific atoms (e.g., {"A181": "CG2,CG1", "A407": "NE1,CZ2"})
- `select_fixed_atoms`: Union[bool, str, Dict[str, str]] = None - Chain ranges to control atomic constraints
- `select_buried`: Union[str, Dict[str, str]] = None - Residues/atoms to drive toward a buried environment
- `select_exposed`: Union[str, Dict[str, str]] = None - Residues/atoms to drive toward a solvent-exposed environment
- `select_hbond_donor`: Dict[str, str] = None - Atoms to constrain as hydrogen-bond donors
- `select_hbond_acceptor`: Dict[str, str] = None - Atoms to constrain as hydrogen-bond acceptors
- `json_config`: Union[str, Dict] = None - Raw foundry JSON config (string or dict) to pass through directly
- `design_startnum`: int = 1 - Starting number for design IDs

**Streams**: `structures`

**Tables**:
- `structures`:

  | id | source_id | pdb | contig | length | design_name | status |
  |----|-----------|-----|--------|--------|-------------|--------|

**Example**:
```python
from biopipelines.rfdiffusion3 import RFdiffusion3

# De novo protein design
rfd3 = RFdiffusion3(
    contig="80-120",
    num_designs=10
)

# Binder design with hotspots
# You can see the name of each atom in pymol picking mode
# If your ligand pdb doesn't have numbering of atoms (e.g. C, C, ...) you can assign individual ids in pymol with the command rename, and then save the pdb
rfd3 = RFdiffusion3(
    pdb=target_pdb,
    contig="50-80,A1-100",
    select_hotspots={"A45": "NE,CZ", "A67": "OG"},
    num_designs=20
)

# Enzyme active site design
rfd3 = RFdiffusion3(
    pdb=template_pdb,
    contig="A1-50,40-60,A100-150",
    length="90-110",
    num_designs=50
)
```

---

### RFdiffusionAllAtom

Generates protein structures with explicit modeling of ligands and small molecules. Diffusion model that handles all-atom representation including non-protein entities.

**References**: https://www.science.org/doi/10.1126/science.adl2528.

**Installation**: Works in SE3nv environment as installed with RFdiffusion.

**Parameters**:
- `ligand`: Union[str, DataStream, StandardizedOutput] (required) - Compounds stream (`Ligand(code="ATP")` or any compounds-producing tool) or a 3-letter code; the residue `code` is read from the stream at runtime and passed to `inference.ligand=`. One ligand is reused across all input PDBs.
- `pdb`: Optional[Union[DataStream, StandardizedOutput]] = None - Input PDB template(s). **Multiple input structures are iterated** — one design run per input PDB (`num_designs` each), so a pose ensemble propagates fully into design instead of only the first structure being used.
- `contigs`: str | (TableInfo, column) = "" - Contig specification. Accepts a table column reference resolved per-input-PDB at runtime (each pose can get its own contig).
- `inpaint`: str | (TableInfo, column) = "" - Inpainting specification (also accepts a per-PDB table column reference)
- `num_designs`: int = 1 - Number of designs per input PDB
- `active_site`: bool = False - Use active site model
- `steps`: int = 50 - Inference diffusion steps (Baker-lab default; trained at 200 but ~50 is equivalent quality at ~4x speed, ~20 also reported equivalent)
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
- `guiding_potentials`: str = None - Guiding-potential spec (`potentials.guiding_potentials`). See "Directing the design toward the ligand" below.
- `guide_scale`: float = None - Global multiplier on the guiding potential (`potentials.guide_scale`).
- `guide_decay`: str = None - How potential influence decays over the trajectory (`potentials.guide_decay`): "constant", "linear", "quadratic", or "cubic".

**Directing the design toward the ligand (pocket enclosure):**

A plain contig (e.g. `40-60,A84-182`) constrains only what is *fixed*; it does not tell the model the newly generated region should *contact* the ligand, so the new segment often drifts away from the bound molecule. To pull the diffusing backbone onto the ligand, add the `ligand_ncontacts` guiding potential (an attractive protein–ligand contact term).

Use the typed builder `RFdiffusionAllAtom.GuidingPotential.ligand_ncontacts(...)` rather than hand-writing the hydra string — you get named arguments, defaults, and validation:

```python
rfdaa = RFdiffusionAllAtom(
    pdb=poses,
    ligand=Ligand(code="LIG"),
    contigs="40-60,A84-182",
    num_designs=2,
    guiding_potentials=RFdiffusionAllAtom.GuidingPotential.ligand_ncontacts(
        weight=3, r_0=8, d_0=4),
    guide_decay="quadratic",
)
```

`ligand_ncontacts(weight, r_0, d_0)` — pull the design onto the ligand. `weight` 1–10 (higher = pulls harder), `r_0` ~8 Å (contact switching distance), `d_0` ~4 Å (always-in-contact). Start moderate and inspect; raise `weight`/`guide_scale` if the region still drifts off the ligand.

Notes:
- A raw hydra string is still accepted for `guiding_potentials` (back-compat).
- Alternative directing levers: **RFdiffusion3** conditions natively on ligand burial via `select_buried` (RASA) and on contacts via `select_hotspots` / `select_hbond_donor/acceptor`; **BoltzGen** has a `protein-small_molecule` protocol with `binding_region`. See the RFdiffusion3 and BoltzGen sections.

**Streams**: `structures`

**Tables**:
- `structures`:

  | id | source_id | pdb | fixed | designed | contigs | time | status |
  |----|-----------|-----|-------|----------|---------|------|--------|

**Example**:
```python
from biopipelines.rfdiffusion_allatom import RFdiffusionAllAtom
from biopipelines.ligand import Ligand

rfdaa = RFdiffusionAllAtom(
    pdb=template,
    ligand=Ligand(code="LIG"),
    contigs='10-20,A6-140',
    num_designs=5
)
```

---
