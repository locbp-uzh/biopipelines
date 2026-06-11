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

### RFdiffusion2

Enzyme active-site scaffolding and small-molecule binder design. Its defining capability is **atomic** motif specification: an active site is given as individual side-chain atoms (`contig_atoms`) that the model scaffolds either at fixed sequence positions (indexed) or anywhere along the chain (`contig_as_guidepost=True`, unindexed). Also supports RASA-conditioned buried-binder design. Distinct from RFdiffusionAllAtom (earlier all-atom model) and RFdiffusion3 (the foundry/rfd3 successor).

**References**:
- Github: https://github.com/RosettaCommons/RFdiffusion2
- Docs: https://rosettacommons.github.io/RFdiffusion2/

**Installation**:

`RFdiffusion2.install()` clones the repo and runs the upstream `setup.py`, which downloads the model weights and the Apptainer container image (`rf_diffusion/exec/bakerlab_rf_diffusion_aa.sif`). The tool runs through that container — point `containers.RFdiffusion2` at the `.sif` in `config.yaml`. The download is large (30+ minutes); `python setup.py overwrite` resumes if interrupted.

**Parameters**:
- `ligand`: Union[str, DataStream, StandardizedOutput] - Bound ligand(s) as a compounds stream (`Ligand(code="NAD")` or a tool's compounds output). Every id's 3-letter `code` is read at runtime and joined into `inference.ligand='NAD,OXM'`. A bare string is shorthand for `Ligand(code=...)`.
- `pdb`: Union[DataStream, StandardizedOutput] - Input active-site/theozyme (or ligand-only binder target) structure(s). **Multiple input structures are iterated** — one design run per input PDB (`num_designs` each).
- `contigs`: str | (TableInfo, column) - Contig specification (RFdiffusion2 uses `,` separators, e.g. `"46,A106-106,59,A166-166,2,A169-169,23,A193-193,46"`, or a bare length like `"150"` for binder design). A string is broadcast to every input PDB; a table column reference is resolved per input-PDB id at runtime.
- `contig_atoms`: Dict[str, str] = None - Per-residue side-chain atoms defining the atomic motif, e.g. `{"A106": "NE,CD,CZ", "A166": "OD1,CG"}` → `contigmap.contig_atoms`.
- `contig_as_guidepost`: bool = True - True scaffolds the motif unindexed (the model places it); False keeps the indexed contig positions. → `inference.contig_as_guidepost`.
- `only_guidepost_positions`: str = None - Restrict guidepost treatment to these positions (e.g. `"A106"`); the rest stay indexed.
- `partially_fixed_ligand`: Dict[str, List[str]] = None - Per-ligand subset of atoms to keep fixed, e.g. `{"NAD": ["O7N", "C7N"]}` → `++inference.partially_fixed_ligand` (appended; not in the base config struct).
- `length`: str = None - Total design length for binder mode (e.g. `"150-150"`) → `contigmap.length`.
- `relative_sasa`: float = None - Relative-SASA target for RASA-conditioned binder design (e.g. `0` to fully bury). Enables `inference.conditions.relative_sasa_v2`.
- `num_designs`: int = 1 - Number of designs per input PDB.
- `design_startnum`: int = 1 - Starting design number.
- `deterministic`: bool = False - Fixed RNG for reproducible outputs.
- `seed_offset`: int = None - Seed offset when deterministic.

The sweep/benchmark orchestrator (`benchmark/pipeline.py`, `sweep.*`), the `stop_step='end'` sequence-fitting/folding follow-on, and PyMOL visualization are **not** wrapped — the framework owns scheduling and downstream sequence/fold steps are separate tools.

**Outputs**:
- `structures` (`.pdb`, with `.trb` trajectory siblings)
- `tables.structures`: `id | pdb | fixed | designed | source_fixed | plddt_mean | status`

**Example** (active-site scaffolding, the open-source demo case):

```python
from biopipelines.rfdiffusion2 import RFdiffusion2
from biopipelines.entities import PDB, Ligand

site = PDB("/path/to/M0584_1ldm.pdb")
rfd2 = RFdiffusion2(
    pdb=site,
    ligand=Ligand(code=["NAD", "OXM"]),  # one compounds stream, two codes -> inference.ligand='NAD,OXM'
    contigs="46,A106-106,59,A166-166,2,A169-169,23,A193-193,46",
    contig_atoms={"A106": "NE,CD,CZ", "A166": "OD1,CG",
                  "A169": "NH2,CZ", "A193": "NE2,CD2,CE1"},
    contig_as_guidepost=True,
    num_designs=8,
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
- `select_buried`: Union[str, Dict[str, str]] = None - Residues/atoms to drive toward a buried environment (RASA). For a custom (non-CCD) ligand, code it `UNL` first (see the custom-ligand example) or the selection silently hits zero atoms.
- `select_partially_buried`: Union[str, Dict[str, str]] = None - Third RASA label between buried and exposed
- `select_exposed`: Union[str, Dict[str, str]] = None - Residues/atoms to drive toward a solvent-exposed environment
- `select_hbond_donor`: Dict[str, str] = None - Atoms to constrain as hydrogen-bond donors
- `select_hbond_acceptor`: Dict[str, str] = None - Atoms to constrain as hydrogen-bond acceptors
- `unindex`: Union[str, Dict[str, str]] = None - Unindexed motif components (no fixed sequence index) — the atomic-motif enzyme path (catalytic tip atoms) and diffused nucleic-acid motifs. Requires a contig or length; must not overlap the contig.
- `select_unfixed_sequence`: Union[bool, str, Dict[str, str]] = None - Components whose sequence is freed (diffused) rather than fixed. Excludes ligands/DNA (those keep fixed sequence).
- `redesign_motif_sidechains`: Union[bool, str] = None - Fixed-backbone sequence design over the motif when a contig is provided.
- `symmetry`: Union[str, Dict] = None - Symmetry group id ("C3", "D2", "T", "O", "I"); a string is shorthand for `{"id": ...}`. A dict may also set `is_unsym_motif` (comma list of contigs/ligands left unsymmetrized, e.g. DNA strands) and `is_symmetric_motif`. Setting this auto-selects the symmetry sampler.
- `ori_token`: Optional[List[float]] = None - Explicit origin/center-of-mass coordinates.
- `infer_ori_strategy`: str = None - Center-of-mass guidance strategy, "com" or "hotspots".
- `is_non_loopy`: bool = None - Non-loopy global conditioning.
- `plddt_enhanced`: bool = None - pLDDT enhancement (upstream default True).
- `partial_t`: float = None - Angstroms of noise for partial diffusion (≤15 recommended). Requires an input structure; `length` must not be set.
- `cfg`: bool = None - Enable classifier-free guidance (improves adherence to ligand burial / H-bond conditions).
- `cfg_scale`: float = None - CFG guidance strength (paper uses 2.0 for diffused-ligand binders).
- `step_scale`: float = None - Sampler step scale η (paper 1.5; higher = less diverse, more designable).
- `noise_scale`: float = None - Sampler noise γ₀ (paper 0.6; 0.0 for ODE sampling).
- `num_steps`: int = None - Number of denoising timesteps (paper 200).
- `center_option`: str = None - Centering, "all", "motif", or "diffuse".
- `seed`: int = None - Inference seed for reproducibility.
- `json_config`: Union[str, Dict] = None - Raw foundry JSON config (string or dict) to pass through directly
- `design_startnum`: int = 1 - Starting number for design IDs

The `cfg*`, `step_scale`, `noise_scale`, `num_steps`, `center_option`, and `seed` parameters are emitted as Hydra `inference_sampler.*` overrides on the `rfd3 design` command line (mapped to foundry keys: `noise_scale`→`gamma_0`, `num_steps`→`num_timesteps`). All other parameters above are written into the per-design inputs JSON.

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

# Enzyme atomic-motif scaffolding (unindexed catalytic tip atoms + substrate)
rfd3 = RFdiffusion3(
    pdb=active_site_pdb,
    contig="80-150",
    unindex="A50-52",            # catalytic residues, sequence position found during diffusion
    redesign_motif_sidechains=True,
    num_designs=50
)

# Symmetric oligomer (C3) — symmetric noise auto-selects the symmetry sampler
rfd3 = RFdiffusion3(length="100-120", symmetry="C3", num_designs=10)

# Diffused small-molecule binder: placed ligand + buried + classifier-free guidance
# (paper Fig 3c inference settings: η=1.5, γ₀=0.6, 200 steps, CFG scale 2).
# The ligand must carry coordinates — foundry appends it to an input atom array, so
# supply it bound via Ligand(code=, structures=) or pass a pdb with it as HETATM.
rfd3 = RFdiffusion3(
    length="80-120",
    ligand=Ligand(code="SAM", structures=posed_sam),
    select_buried="B1",
    cfg=True, cfg_scale=2.0,
    step_scale=1.5, noise_scale=0.6, num_steps=200,
    num_designs=400,
)

# CUSTOM (non-CCD) ligand — e.g. a SMILES-derived dye. If its residue code collides
# with a real CCD entry (LIG, SAM, …), atomworks loads that component's reference
# conformer and the atom names won't match the structure, so a RASA selector resolves
# to zero atoms and fails with "could not broadcast ... (0,3)". Code it UNL (atomworks'
# DO_NOT_MATCH_CCD sentinel) so its atoms are read from the structure:
prepped = PDB(my_pose, PDB.rename("LIG", "UNL"))
rfd3_custom = RFdiffusion3(
    pdb=prepped, ligand=Ligand(code="UNL", structures=prepped),
    contig="65-95,A84-182", select_buried="B1", cfg=True, cfg_scale=2.0,
)

# Protein-DNA binder: design against an input (e.g. AF3-predicted) DNA structure.
# Keep the DNA strands fixed via select_fixed_atoms; free their sequence to co-diffuse
# the DNA conformation with select_unfixed_sequence.
rfd3 = RFdiffusion3(
    pdb=dna_structure,
    contig="80-120",
    num_designs=100,
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
