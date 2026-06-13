# Analysis

[← Back to Tool Reference](../tool_reference.md)

### ADMETAI

ADMET endpoint predictions for compound libraries via the ADMET-AI Chemprop-RDKit model. Reads SMILES from a compounds stream and writes one row per input compound with all upstream-reported ADMET properties (~40+ regression and classification scores covering absorption, distribution, metabolism, excretion, and toxicity).

**Environment**: `admet_ai`

**Installation**: `ADMETAI.install()` creates a dedicated conda env (Python 3.10) and installs `admet-ai` via pip. Verification instantiates `ADMETModel`, which downloads the bundled weights and warms the cache.

**Parameters**:
- `compounds`: Union[DataStream, StandardizedOutput] (required) — Input compounds. SMILES are read from the `smiles` column of the compounds map_table.

**Tables**:
- `admet`:

  | id | smiles | <ADMET endpoint columns from upstream model> |
  |----|--------|----------------------------------------------|

  The endpoint column set is determined by the installed `admet-ai` version's `ADMETModel.predict()` output and is preserved verbatim in the CSV.

**Example**:
```python
from biopipelines.admet_ai import ADMETAI
from biopipelines.compound_library import CompoundLibrary

with Pipeline("Project", "ADMET", description="Library ADMET screen"):
    Resources(memory="8GB", time="1:00:00")
    library = CompoundLibrary("my_library.csv")
    admet = ADMETAI(compounds=library)

# Filter for high oral bioavailability candidates downstream
from biopipelines.panda import Panda
filtered = Panda(
    tables=admet.tables.admet,
    operations=[Panda.filter("HIA_Hou > 0.8")],
)
```

**Reference**: Swanson et al. (2024) ADMET-AI: a machine learning ADMET platform for evaluation of large-scale chemical libraries. *Bioinformatics* 40, btae416. https://github.com/swansonk14/admet_ai

---

### AF2BIND

Small-molecule binding-residue prediction from a single protein structure, leveraging AlphaFold2's internal pair representation through a trained linear head (sokrypton/af2bind). It probes the target with 20 "bait" amino acids and reads out a per-residue binding probability `p_bind` — no ligand required. GPU-bound.

**Environment**: AF2BIND **reuses LocalColabFold's conda env** on the cluster — it already ships a working `jaxlib+cuda` (the stack ColabFold runs on the cluster's GPUs), so `AF2BIND.install()` only pip-adds ColabDesign `@v1.1.1` rather than building a fragile `jax[cuda]` env. On Colab, ColabDesign is installed into the base Python (which has JAX+GPU). **AF2BIND therefore requires `AlphaFold.install()` (LocalColabFold) first.** It reads the AlphaFold2 network params from the **shared `AlphaFoldParams` cache** (the dir whose `params/` holds `params_model_*.npz`, configured in `folders.cache`) — on the cluster this is LocalColabFold's existing params, so nothing is re-downloaded. Only the AF2BIND linear-head weights are tool-specific and land under the `AF2BIND` folder.

**Parameters**:
- `structures`: Union[DataStream, StandardizedOutput] (required) — Input PDB structures.
- `chain`: str = `"A"` — Target chain to score.
- `mask_sidechains`: bool = True — Hide target side chains from AF2 (selects the side-chain-masked head, the more transferable default).
- `mask_sequence`: bool = False — Also hide target sequence identity.
- `top_k`: int = 15 — Number of highest-`p_bind` residues reported as a chain-aware selection string in `summary.binding_residues`.

**Streams**:
- `binding`: per-residue resi-csv (one `<id>.csv` per input), columns `id | chain | resi | resn | p_bind`. Consumable by the Selection tool, e.g. `Selection.add(af2bind.streams.binding, include="p_bind>0.5")`.

**Tables**:
- `binding`: `id | chain | resi | resn | p_bind` — one row per target-chain residue (same data as the `binding` stream, combined across inputs).
- `summary`: `id | n_residues | top_resi | top_p_bind | binding_residues` — `binding_residues` is the top-k selection (e.g. `"A45+A78-80"`), usable directly as a downstream selection column.
- `missing`: `id | cause`

**Example**:
```python
from biopipelines.af2bind import AF2BIND

target = PDB("6W70", convert="pdb")
bind = AF2BIND(structures=target, chain="A", top_k=20)
bind.tables.binding  # per-residue p_bind
bind.tables.summary  # top binding residues as a selection string
```

---

### Angle

Calculates angles between atoms in structures. Three modes are supported, selected by the shape of the `atoms` tuple. Useful for backbone phi/psi analysis, side chain rotamers, ligand geometry verification, and measuring relative orientations between structural elements.

**Environment**: `biopipelines`

**Parameters**:
- `structures`: Union[ToolOutput, StandardizedOutput] (required) - Input structures
- `atoms`: tuple (required) - Atom selections; form determines the angle mode (see below)
- `metric_name`: str = None - Custom name for angle column (default: `"angle"`, `"torsion"`, or `"vector_angle"` depending on mode)
- `unit`: str = "degrees" - Output unit: `"degrees"` or `"radians"`

**Selection Syntax**:
Same as Distance, with additional support for `residue.atom` format:
- `'10.CA'` - Alpha carbon of residue 10
- `'-1.C'` - Carbonyl carbon of last residue (C-terminus)
- `'LIG.C1'` - Atom C1 of ligand residue LIG
- `'D in IGDWG'` - Aspartic acid in sequence context (uses centroid if multiple atoms)

**Tables**:
- `angles`:

  | id | source_structure | {metric_name} | unit |
  |----|------------------|---------------|------|

**Angle Modes**:

| Form | Mode | Description | Range |
|------|------|-------------|-------|
| `(a, o, b)` | Bond | Angle at vertex o (a–o–b) | 0–180° |
| `(a, x1, x2, b)` | Torsional | Dihedral along axis x1–x2 | −180° to 180° |
| `((a1, a2), (b1, b2))` | Vector | Angle between vectors a1→a2 and b1→b2 | 0–180° |

Each selection may resolve to multiple atoms; centroids are used in that case.

**Example**:
```python
from biopipelines.angle import Angle

# Bond angle at CA (N-CA-C)
bond_angle = Angle(
    structures=boltz,
    atoms=('10.N', '10.CA', '10.C'),
    metric_name="nca_angle"
)

# Phi angle (C-N-CA-C)
phi = Angle(
    structures=boltz,
    atoms=('9.C', '10.N', '10.CA', '10.C'),
    metric_name="phi"
)

# Psi angle (N-CA-C-N)
psi = Angle(
    structures=boltz,
    atoms=('10.N', '10.CA', '10.C', '11.N'),
    metric_name="psi"
)

# Chi1 angle
chi1 = Angle(
    structures=boltz,
    atoms=('50.N', '50.CA', '50.CB', '50.CG'),
    metric_name="chi1"
)

# Ligand geometry
ligand_angle = Angle(
    structures=boltz,
    atoms=('LIG.C1', 'LIG.C2', 'LIG.C3')
)

# Angle between two structural vectors (e.g. helix axis proxies)
orientation = Angle(
    structures=boltz,
    atoms=(('66.NE1', '66.CA'), ('-173.OH', '-173.CA')),
    metric_name="orientation"
)

# Output in radians (for use with cos/sin in Panda.calculate)
orientation_rad = Angle(
    structures=boltz,
    atoms=(('66.NE1', '66.CA'), ('-173.OH', '-173.CA')),
    metric_name="orientation",
    unit="radians"
)
```

---

### APBS

Electrostatic surface potential. For each input PDB, runs `pdb2pqr` (PROPKA protonation at the given pH) then `apbs` to compute the Poisson-Boltzmann potential on a grid. Exposes the protonated/charged structure (PQR) and the potential grid (DX), plus per-structure charge metrics.

**References**: https://github.com/Electrostatics/apbs

**Environment**: `apbs`

**Parameters**:
- `structures`: DataStream | StandardizedOutput (required) — Input PDBs.
- `ph`: float = 7.4 — Protonation pH for pdb2pqr/PROPKA.
- `forcefield`: str = "AMBER" — pdb2pqr force field for charge/radius assignment.
- `ion_concentration`: float = 0.15 — Mobile-ion concentration (M).
- `grid_dim`: int = 65 — Grid points per dimension.
- `pdie`: float = 2.0 — Solute (protein) dielectric.
- `sdie`: float = 78.5 — Solvent dielectric.
- `solver`: str = "lpbe" — `"lpbe"` (linearized) or `"npbe"` (nonlinear) Poisson-Boltzmann.

**Streams**: `structures` (PQR), `grids` (DX potential grid)

**Tables**:
- `electrostatics`: | id | net_charge | n_basic | n_acidic | isoelectric_point |
- `missing`: | id | cause |

**Example**:
```python
from biopipelines.apbs import APBS

target = PDB("4UFC", convert="pdb")
elec = APBS(structures=target, ph=7.4)
elec.tables.electrostatics
```

---

### BioEmu

Emulates a protein's equilibrium structural ensemble from sequence with a generative model. For each input sequence, samples `num_samples` statistically independent conformers approximating the equilibrium distribution — a fast alternative to running MD. Emits per-conformer PDBs plus a compact trajectory.

**References**: https://github.com/microsoft/bioemu · https://www.science.org/doi/10.1126/science.adv9817

**Resources**: GPU.

**Environment**: `bioemu`

**Installation**: `BioEmu.install()` creates the env and installs the package (weights download on first run).

**Parameters**:
- `sequences`: DataStream | StandardizedOutput (required) — Input sequences.
- `num_samples`: int = 10 — Conformers to sample per sequence (filtering may drop a few).
- `batch_size`: int = 10 — Samples per forward pass.
- `reconstruct_sidechains`: bool = False — Rebuild all-atom side chains (default backbone-only).
- `filter_samples`: bool = True — Drop structurally implausible samples.
- `msa_host_url`: str = "" — Optional custom MSA server URL.

**Streams**:
- `structures` — per-conformer PDBs, ids `<seq>_<1..N>`.
- `trajectories` — one `samples.xtc` + topology per sequence.

**Tables**:
- `summary`: | id | sequence.id | n_samples | n_residues |
- `missing`: | id | cause |

**Example**:
```python
from biopipelines.bioemu import BioEmu

seqs = Sequence("YYDPETGTWY", ids="chignolin")
ens = BioEmu(sequences=seqs, num_samples=10)
```

---

### CABSflex

Fast coarse-grained Monte Carlo simulation of protein structure flexibility. Produces conformational ensembles and per-residue RMSF profiles. Optionally rebuilds models to all-atom representation using MODELLER.

**WARNING**: MODELLER requires a license key (free for academics). Get one at https://salilab.org/modeller/registration.html and set `KEY_MODELLER` before running.

**Environment**: `CABSflex` (Python 2.7 with `modeller` and `cabs`)

**Installation**: `conda create -n CABSflex python=2.7 && conda install -c salilab modeller -c lcbio cabs`. On Colab: `pip install modeller cabs`.

**Parameters**:
- `structures`: Union[DataStream, StandardizedOutput] (required) - Input protein structures
- `num_models`: int = 10 - Number of cluster medoids / final models per input structure
- `mc_cycles`: int = 50 - Monte Carlo cycles between trajectory frames
- `mc_steps`: int = 50 - Monte Carlo steps
- `mc_annealing`: int = 20 - Temperature annealing cycles
- `temperature`: Optional[str] = None - Temperature range as "TINIT TFINAL" (default: "1.4 1.4")
- `flexibility`: Optional[str] = None - Residue flexibility: float (0=flexible, 1=stiff), 'bf', 'bfi', 'bfg', or filename
- `filtering_count`: int = 1000 - Number of low-energy models for clustering
- `aa_rebuild`: bool = False - Rebuild to all-atom with MODELLER
- `restraints`: str = "ss2" - Restraint mode: 'all', 'ss1', 'ss2'
- `restraints_gap`: int = 3 - Min gap along chain for restraints
- `restraints_min`: float = 3.8 - Min distance in Å for restraints
- `restraints_max`: float = 8.0 - Max distance in Å for restraints
- `restraints_reduce`: Optional[float] = None - Randomly drop a fraction of restraints (0–1) to loosen the simulation
- `weighted_fit`: Optional[str] = None - Fit method: 'gauss', 'flex', 'ss', 'off', or filename (None = CABSflex default)
- `pdb_output`: str = "M" - Which PDB output to write ('M' = models, etc.)
- `max_parallel`: int = 1 - Number of input structures to simulate in parallel

**Streams**:
- `structures`: PDB ensemble models (`num_models` per input structure)
- `images`: SVG plots (RMSF, RMSD, energy) per input structure
- `rmsf`: Per-residue RMSF CSVs (`resi-csv` format, columns: id, chain, resi, rmsf), one file per input structure

**Tables**:
- `structures`:

  | id | file | structures.id |
  |----|------|---------------|

- `rmsf_all` (merged across all input structures):

  | id | chain | resi | rmsf |
  |----|-------|------|------|

**Example**:
```python
from biopipelines.cabsflex import CABSflex

# Basic flexibility simulation
protein = PDB("1AHN")
flex = CABSflex(structures=protein, num_models=10)

# Access merged RMSF table
flex.tables.rmsf_all  # all structures combined

# Access per-structure RMSF stream
for rmsf in flex.streams.rmsf:
    print(f"{rmsf.ids[0]}: {rmsf.files[0]}")

# Access ensemble models
flex.streams.structures  # 10 PDB models per input

# Custom simulation parameters
flex = CABSflex(
    structures=protein,
    num_models=5,
    mc_cycles=100,
    mc_annealing=30,
    temperature="1.4 2.0",
    flexibility=0.5,
    aa_rebuild=True
)
```

---

### Aggrescan3D

Structure-based prediction of protein aggregation propensity. A3D scores each residue using an experimentally-derived aggregation scale combined with structural context (solvent exposure), so it picks up aggregation-prone patches that are only apparent in the folded state. Negative scores indicate soluble residues, positive scores aggregation-prone ones; the per-structure average is a global solubility/aggregation readout. Complements the sequence-based solubility predictor PLM_Sol (which needs no structure). The `aggregation` stream uses the same `resi-csv` schema as CABSflex's `rmsf`, so `Selection` thresholds on it unchanged (e.g. to pick aggregation-prone residues for redesign).

This wrapper runs A3D in static mode only. The FoldX-backed repair / solubility-enhancing mutation modes and the CABS-flex-backed dynamic mode are not exposed.

**Environment**: `Aggrescan3D` (Python 2.7 with `aggrescan3d`)

**Installation**: `Aggrescan3D.install()` creates a dedicated env: `conda create -n Aggrescan3D python=2.7 && conda install -c lcbio -c conda-forge aggrescan3d`. The `aggrescan3d` package pulls its own SASA backend; pinning a standalone modern `freesasa` would conflict with the Python 2.7 requirement.

**Parameters**:
- `structures`: Union[DataStream, StandardizedOutput] (required) - Input protein structures
- `chains`: Optional[str] = None - Restrict the analysis to specific chain(s), e.g. `"A"` or `"AB"` (A3D `-C`). None analyses all chains.
- `distance`: float = 10.0 - Distance value (Å) used in the A3D score calculation (A3D `-D`)
- `max_parallel`: int = 1 - Number of input structures to score in parallel

**Streams**:
- `aggregation`: Per-residue A3D scores (`resi-csv` format, columns: id, chain, resi, score), one file per input structure
- `structures`: A3D output PDB per input, B-factor field replaced with the A3D score (paint it in PyMOL to visualise aggregation patches)
- `images`: Per-chain A3D score plots (PNG)

**Tables**:
- `scores` (per-structure global summary):

  | id | structures.id | avg_score | min_score | max_score | n_residues | n_aggregation_prone |
  |----|---------------|-----------|-----------|-----------|------------|---------------------|

- `aggregation_all` (merged across all input structures):

  | id | chain | resi | score |
  |----|-------|------|-------|

**Example**:
```python
from biopipelines.aggrescan3d import Aggrescan3D

# Score aggregation propensity of a predicted/designed structure
af = AlphaFold(proteins=pmpnn)
agg = Aggrescan3D(structures=af)

# Per-structure global solubility ranking
agg.tables.scores            # avg_score: more negative = more soluble

# Threshold on the per-residue profile, exactly as with CABSflex's rmsf stream
hotspots = Selection(
    Selection.add(agg.streams.aggregation, include="score>0"),
    structures=af
)
```

**Reference**: Kuriata et al. (2019) Aggrescan3D (A3D) 2.0: prediction and engineering of protein solubility. *Nucleic Acids Research* 47, W300–W307. https://bitbucket.org/lcbio/aggrescan3d

---

### EnsembleAnalysis

Per-residue RMSF and ensemble-level metrics from a conformer ensemble. Where CABSflex couples RMSF to its own coarse-grained sampling, EnsembleAnalysis analyzes *any* ensemble: it superposes the conformers (least-squares on CA or backbone) and reports per-residue fluctuation, so it overlays RMSF profiles from NMR ensembles, PLACER dumps, BioEmu samples, or any pool of conformers on the same footing. The `rmsf` stream uses the same `resi-csv` schema as CABSflex, so `Selection` thresholds on it unchanged.

**Environment**: `biopipelines` (numpy + the shared PDB parser; no external tool, no dedicated env)

**Installation**: none — `EnsembleAnalysis.install()` is a no-op marker; the tool runs in the base `biopipelines` environment.

**Inputs**: the ensemble partition is defined explicitly by a required `groups=` stream (the Consensus convention): its ids ARE the output ensemble ids, known at config time. Each `structures` id is matched to a group and all of its models pooled into that group's ensemble. The two input shapes compose freely under this rule:
- A **multi-model file** per id (several `MODEL` records, e.g. NMR/PLACER) — pass `groups=` a one-id stream to treat the file as a single ensemble.
- A **conformer set** — single-model files (e.g. BioEmu's `<seq>_1..N`); pass `groups=` the upstream sampler's input (e.g. the folded sequence) so the pooled ensemble id is the parent.

**Parameters**:
- `structures`: Union[DataStream, StandardizedOutput] (required) - The ensemble source.
- `groups`: Union[DataStream, StandardizedOutput] (required) - Stream whose ids define the ensemble partition; the output ensemble ids are these group ids. Pass the upstream sampler's input (e.g. the folded sequence), or a one-id stream for a single multi-model file.
- `selection`: str = "CA" - Atom set for superposition and RMSF: `"CA"` (alpha carbon, the standard RMSF convention) or `"backbone"` (N, CA, C, O).
- `reference`: str = "mean" - Superposition reference: `"mean"` (iterative align-to-mean) or `"first"` (align all frames to the first model).

**Streams**:
- `rmsf`: Per-residue `resi-csv`, one file per ensemble id. Columns: id, chain, resi, rmsf, rmsd_mean (RMSF = fluctuation about the mean; rmsd_mean = deviation about the reference).

**Tables**:
- `residues` (merged across all ensembles — the `rmsf` stream concatenated):

  | id | chain | resi | rmsf | rmsd_mean |
  |----|-------|------|------|-----------|

- `ensemble` (one row per ensemble):

  | id | n_frames | n_residues | mean_rmsf | max_rmsf | rg_mean | rg_std |
  |----|----------|------------|-----------|----------|---------|--------|

- `frames` (one row per conformer):

  | id | frame | rmsd_to_ref | rg |
  |----|-------|-------------|----|

**Example**:
```python
from biopipelines.ensemble_analysis import EnsembleAnalysis
from biopipelines.bioemu import BioEmu
from biopipelines.selection import Selection

# Overlay RMSF from a BioEmu ensemble (per-conformer PDBs pooled per sequence)
ens = BioEmu(sequences=seqs, num_samples=20)
rmsf = EnsembleAnalysis(structures=ens, groups=seqs)   # one profile per sequence

# Or from a multi-model NMR ensemble (one file, many MODELs): a one-id groups stream
nmr = PDB("2N5E")
rmsf = EnsembleAnalysis(structures=nmr, groups=nmr, selection="backbone")

# Threshold on flexibility, exactly as with CABSflex's rmsf stream
flexible = Selection.add(rmsf.streams.rmsf, include="rmsf>1.0")

# Merged per-residue table and per-ensemble summary
rmsf.tables.residues
rmsf.tables.ensemble
```

---

### ConformationalChange

Quantifies structural changes between reference and target structures using PyMOL's alignment RMSD.

**Environment**: `ProteinEnv`

**Parameters**:
- `reference_structures`: Union[DataStream, StandardizedOutput] (required) - Reference structures. Can be one or the same number as targets.
- `target_structures`: Union[DataStream, StandardizedOutput] (required) - Target structures to compare
- `selection`: Optional[str] = None - Residue range (e.g., '10-20+30-40'). None = whole structure.
- `alignment`: str = "align" - Alignment method (align, super, cealign). Rule of thumb: sequence similarity > 50% -> align; otherwise cealign.
- `atoms`: str = "all" - Which atoms to use for alignment:
  - `"all"` (default): all atoms
  - `"CA"`: alpha-carbon only
  - `"backbone"`: backbone atoms (CA+C+N+O)
  - Any `+`-separated atom names, e.g. `"CA+CB"`

**Tables**:
- `changes`:

  | id | reference_structure | target_structure | selection | num_aligned_atoms | RMSD |
  |----|---------------------|------------------|-----------|-------------------|------|

**Example**:
```python
from biopipelines.conformational_change import ConformationalChange

# All-atom RMSD on a specific region
conf_change = ConformationalChange(
    reference_structures=apo_structures,
    target_structures=holo_structures,
    selection="10-50",
    alignment="super"
)

# CA-only RMSD
conf_change = ConformationalChange(
    reference_structures=design,
    target_structures=refolded,
    atoms="CA"
)

# Backbone RMSD on redesigned region
conf_change = ConformationalChange(
    reference_structures=kinase,
    target_structures=refolded,
    selection=backbones.tables.structures.fixed,
    atoms="backbone"
)
```

---

### Contacts

**Environment**: `ProteinEnv`

Analyzes contacts between selected protein regions and ligands. For each selected residue, calculates the minimum distance to any ligand atom. Returns contact counts and distance statistics.

**Parameters**:
- `structures`: Union[ToolOutput, StandardizedOutput] (required) - Input structures
- `selections`: Union[str, Tuple[TableInfo, str], None] = None - Protein region selections. A string (`'10-20+30-40'`) applies to every structure; a `(table, "column")` reference resolves a per-structure selection by ID match at runtime (a single-row table broadcasts to all structures, otherwise structures absent from the table are skipped); `None` analyzes all protein residues.
- `ligand`: Union[DataStream, StandardizedOutput, None] = None - Compounds stream naming the ligand whose contacts are counted. Required only when `reference="ligand"`; the residue `code` is read from the stream's map_table at runtime (Ligand Contract).
- `reference`: Union[str, Tuple[TableInfo, str]] = `"ligand"` - What contacts are counted against. The literal `"ligand"` counts contacts against the ligand resolved from the `ligand` stream; a residue selection string accepted by `pdb_parser.resolve_selection` (`"84-182"`, `"A84-182"`, `"10+15+20"`) counts against a fixed protein region for all structures; a `(table, "column")` reference resolves per-structure residue references by ID match at runtime (single-row table broadcasts to all structures, otherwise structures absent from the table are skipped). Reference residues are excluded from the protein side so a residue is not counted as contacting itself.
- `contact_threshold`: float = 5.0 - Distance threshold for counting contacts (Angstroms)
- `contact_metric_name`: str = None - Custom name for contact count column (default: "contacts")

**Tables**:
- `contacts`:

  | id | source_structure | selections | ligand | contacts | min_distance | max_distance | mean_distance | sum_distances_sqrt_normalized |
  |----|------------------|------------|--------|----------|--------------|--------------|---------------|-------------------------------|

**Output Columns**:
- `id`: Structure identifier
- `source_structure`: Path to input structure file
- `selections`: Protein residues analyzed (e.g., '10-20+30-40' or 'all_protein')
- `ligand`: Reference identifier — the ligand residue name when `reference="ligand"`, otherwise the residue selection string (e.g. `"84-182"`).
- `contacts`: Number of residues within contact_threshold distance
- `min_distance`: Minimum distance from any selected residue to the reference (Å)
- `max_distance`: Maximum distance from any selected residue to the reference (Å)
- `mean_distance`: Mean distance from selected residues to the reference (Å)
- `sum_distances_sqrt_normalized`: Sum of distances divided by √(number of residues)

**Example**:
```python
from biopipelines.contacts import Contacts

# Ligand code read from the complex's compounds stream at runtime.
contacts = Contacts(
    structures=rfdaa,
    selections=rfdaa.tables.structures.designed,
    ligand=rfdaa,      # residue code resolved from the compounds stream
    contact_threshold=5.0
)

# A standalone Ligand works too.
contacts = Contacts(
    structures=boltz,
    selections='10-20+30-40',
    ligand=Ligand("ATP"),
    contact_threshold=4.0
)

# Analyze all protein residues
contacts = Contacts(
    structures=boltz,
    ligand=boltz
)

# Residue-vs-residue: count contacts between the designed N-segment ("1-50")
# and the fixed core ("84-182") to gate new-segment <-> core packing.
core_contacts = Contacts(
    structures=rfdaa,
    selections="1-50",
    reference="A84-182",
)
```

---

### Distance

Measures distances between specific atoms and residues in structures. Useful for tracking ligand-protein interactions or structural features.

**Environment**: `biopipelines`

**Parameters**:
- `structures`: Union[ToolOutput, StandardizedOutput] (required) - Input structures
- `atom`: str (required) - Atom selection (e.g., 'LIG.Cl', 'name CA', 'A10.CA')
- `residue`: str (required) - Residue selection (e.g., 'D in IGDWG', '145', 'resn ALA')
- `method`: str = "min" - Distance calculation method (min, max, mean, closest)
- `metric_name`: str = None - Custom name for distance column in output (default: "distance")
- `unit`: str = "angstrom" - Output unit: "angstrom" (default) or "nm"

**Tables**:
- `distances`:

  | id | source_structure | {metric_name} | unit |
  |----|------------------|---------------|------|

**Example**:
```python
from biopipelines.distance import Distance

distances = Distance(
    structures=boltz,
    atom="LIG.Cl",
    residue="D in IGDWG",
    method="min",
    metric_name="chlorine_distance"
)
```


---

### DistanceSelector

Selects protein residues based on proximity to a reference — a ligand, a residue range, or any atom/atom-set expressible in the framework's selection grammar (e.g. `LIG.O132`, `LIG.Cl+LIG.O3`, `87.CA`, `A141.CB`, `first.CA`, `D in IGDWG`). For a sequence-context selection like `D in IGDWG`, the context is matched against the chain sequence and the **first** occurrence of the target residue within the matched window is taken — pick a context unique in the sequence, and put the intended residue ahead of any other copy of the same letter in the window. Partitions residues into `within` / `beyond` by distance cutoff, a top-K cap, or both, and emits the per-residue distances as a `resi-csv` stream for downstream thresholding via [Selection](data_management.md#selection).

**Installation**: It requires an environment containing pandas (e.g. biopipelines).

**Parameters**:
- `structures`: DataStream | StandardizedOutput (required) — Input structures.
- `ligand`: Ligand | StandardizedOutput | None — Required only when `reference="ligand"`. The residue `code` is read from the compounds stream's map_table at runtime (Ligand Contract).
- `distance`: float | None = 5.0 — Distance cutoff in Angstroms. Optional when `top_k` is set.
- `reference`: str = "ligand" — Either the literal `"ligand"` (resolves the ligand's residue code from the compounds stream at runtime) or any selection expression accepted by `pdb_parser.resolve_selection`.
- `top_k`: int | None = None — Cap to the K nearest residues. Combined with `distance`, takes the K nearest among those within the cutoff. At least one of `distance` / `top_k` must be set.
- `mode`: `"min"` | `"centroid"` = `"min"` — `"min"` measures min over reference × residue atoms; `"centroid"` measures from each residue's atoms to the reference atom set's centroid.
- `atoms`: `"all"` | `"backbone"` | `"CA"` | `"sidechain"` = `"all"` — Restricts which atoms enter the distance calc on BOTH protein and reference sides (`backbone` = N, CA, C, O).
- `restrict_to`: str | (TableInfo, column) | None = None — Restrict the candidate residue set. Accepts a residue selection string, a bare chain spec (`"B"` or `"chain B"`), or a per-structure column reference.
- `include_reference`: bool = True — When the reference is itself protein residues, whether to include them in `within`. No-op for ligand/atom references.

**Streams**:
- `distances` — Per-residue `resi-csv`, one CSV per structure, columns: `id | chain | resi | distance`. Consumable by [Selection](data_management.md#selection), e.g. `Selection.add(ds.streams.distances, include="distance<5.0")`.

**Tables**:
- `selections`:

  | id | pdb | within | beyond | distance_cutoff | top_k | mode | reference |
  |----|-----|--------|--------|-----------------|-------|------|-----------|

**Example**:
```python
from biopipelines.distance_selector import DistanceSelector

# Ligand reference resolved from the complex's compounds stream.
ds = DistanceSelector(structures=boltz, ligand=boltz, distance=8.0)

# Atom-level reference: residues near a specific ligand atom.
ds = DistanceSelector(structures=boltz, reference="LIG.O132", distance=6.0)

# Multiple ligand atoms (union): residues near either boron atom. The
# 'LIG.B' prefix form selects every atom whose name starts with B (e.g. B41, B42).
ds = DistanceSelector(structures=boltz, reference="LIG.B41+LIG.B42", distance=5.0)

# Top-K nearest residues to a ligand atom, ignoring cutoff.
ds = DistanceSelector(structures=boltz, reference="LIG.O132", top_k=3, distance=None)

# Cα-only pocket within 8 Å of an active-site residue.
ds = DistanceSelector(structures=af, reference="87.CA", distance=8.0, atoms="CA")
```

---

### Consensus

General per-group aggregator for `resi-csv` streams. Collapses an N-file (per-id) `resi-csv` into per-`(chain, resi)` aggregate values *within each group*, and emits a new `resi-csv` stream carrying the aggregate columns. It is value-agnostic — it knows nothing about selections, distances, or any column's meaning; it groups the input ids by `groups=` (using the framework's id matching) and applies aggregation operations across each group, the way [Panda](data_management.md#panda) applies table operations.

The canonical use is to turn a per-pose proximity profile (e.g. [DistanceSelector](#distanceselector)'s `distances`) into a *consensus pocket* held fixed across an ensemble: aggregate the fraction of poses in which each residue is within a cutoff, then threshold that fraction with [Selection](data_management.md#selection).

**Installation**: It requires an environment containing pandas (e.g. biopipelines).

**Parameters**:
- `stream`: DataStream | StandardizedOutput (required) — A `resi-csv` stream (one CSV per id, rows keyed by `(chain, resi)` with value columns).
- `groups`: DataStream | StandardizedOutput (required) — The stream whose ids define the partition. Each input id is matched to a group id via the framework's id matching (`x3_2` → group `x3`); for a single group id, all input ids collapse into one consensus.
- `operations`: list of `ConsensusOp` (required) — One or more aggregations, each adding one output column.

**Operations** (factory methods on `Consensus`):
- `Consensus.fraction(predicate, name="frequency")` — fraction of a group's ids whose row at that `(chain, resi)` satisfies `predicate` (`"column op value"`, e.g. `"distance<=6"`). The denominator is the group size, so ids lacking that residue count as not-passing.
- `Consensus.mean(column, name=...)` / `Consensus.min(column, name=...)` / `Consensus.max(column, name=...)` — aggregate `column` across the group's rows present at that `(chain, resi)`.
- `Consensus.count(name="count")` — number of the group's ids with a row at that `(chain, resi)`.

**Streams**:
- The output stream **keeps the input stream's name** (Panda/Pool convention) — `Consensus(dsel.streams.distances, …)` is read back as `consensus.streams.distances`. It is a `resi-csv` with columns `id | chain | resi | <op columns…> | n_group`, broadcast **one file per input id** (each id carries its own group's aggregated rows), so a downstream stream-consuming tool stays id-keyed and matches per id exactly. Consumable by [Selection](data_management.md#selection), e.g. `Selection.add(c.streams.distances, include="frequency>=0.5")`.

**Example**:
```python
from biopipelines.consensus import Consensus
from biopipelines.selection import Selection

# Per-pose distances to a bound ligand across an ensemble of one protein.
dsel = DistanceSelector(structures=poses, ligand=lig, distance=6.0, restrict_to="A")

# Consensus pocket: residue within 6 A of the ligand in >= 50% of the poses.
consensus = Consensus(dsel.streams.distances, groups=protein,
                      operations=[Consensus.fraction("distance<=6", name="frequency")])

# Pocket = within; surface = its complement, held fixed across every design.
surface = Selection(Selection.add(consensus.streams.distances, include="frequency>=0.5"),
                    Selection.invert(), structures=poses)
ProteinMPNN(structures=poses, redesigned=surface.tables.selections.selection)
```

---

### DSSP

Per-residue secondary-structure assignment. For each input PDB, runs DSSP (`mkdssp`) and reports the 8-state code per residue plus a per-structure helix/sheet/coil summary.

**References**: https://github.com/PDB-REDO/dssp

**Environment**: `dssp`

**Parameters**:
- `structures`: DataStream | StandardizedOutput (required) — Input PDBs.

**Streams**:
- `dssp` — raw DSSP output per structure.
- `ss` — per-residue resi-csv (consumable by [Selection](data_management.md#selection)).

**Tables**:
- `secondary_structure`: | id | chain | resnum | resi | resname | resn | ss_code | ss_simple |
- `summary`: | id | n_residues | helix_frac | sheet_frac | coil_frac |

**Example**:
```python
from biopipelines.dssp import DSSP

target = PDB("4UFC", convert="pdb")
ss = DSSP(structures=target)
ss.tables.summary  # helix/sheet/coil fractions per structure
```

---

### FPocket

Geometry-based binding-pocket detection. For each input PDB, runs the `fpocket` binary and reports detected pockets with their alpha-sphere count, volume, druggability score, and constituent residues. No ligand or model inference required.

**Environment**: `fpocket` (bioconda binary).

**Parameters**:
- `structures`: Union[DataStream, StandardizedOutput] (required) — Input PDB structures.
- `min_alpha_spheres`: int = 35 — Minimum alpha-spheres per reported pocket (fpocket `-i`).
- `min_radius`: float = 3.0 — Minimum alpha-sphere radius in Å (fpocket `-m`).
- `max_radius`: float = 6.0 — Maximum alpha-sphere radius in Å (fpocket `-M`).
- `clustering_distance`: float = 2.4 — Single-linkage clustering distance for grouping alpha-spheres into pockets (fpocket `-D`).

**Streams**:
- `structures`: the original PDB annotated with pocket alpha-spheres as `HETATM STP` records (`<id>.pdb`; the `_out` suffix fpocket uses internally is dropped).

**Tables**:
- `pockets`: `id | pocket_idx | druggability | volume | n_alpha_spheres | n_residues | residues | pocket_file`
- `summary`: `id | n_pockets | top_druggability | top_volume | pymol_script`
- `missing`: `id | cause`

**Example**:
```python
from biopipelines.fpocket import FPocket

target = PDB("4UFC", convert="pdb")
pockets = FPocket(structures=target)
pockets.tables.pockets  # one row per detected pocket
```

---

### GEMS

Protein–ligand **binding-affinity** prediction via a graph neural network with protein/ligand language-model embeddings. For each id, scores a protein structure against a posed ligand and reports a predicted pKd.

> **The ligand must be posed inside the pocket.** GEMS reads the ligand's 3-D coordinates, so the SDF/PDB must hold the bound pose — a SMILES embedded from scratch lands at the origin and produces an empty graph. Pass a ligand-producing tool's output (e.g. a Boltz2 complex, or `Ligand("STI", structures=complex)` → `OpenBabel(convert_3d="sdf")`) as `ligands=`.

**References**: https://github.com/camlab-ethz/GEMS

**Resources**: GPU.

**Environment**: `gems`

**Parameters**:
- `structures`: DataStream | StandardizedOutput (required) — Per-id protein PDBs.
- `ligands`: StandardizedOutput = None — A ligand-producing tool output; coordinate files from its `structures` stream, SMILES from its `compounds` stream. Mutually exclusive with `ligands_3d`/`ligands_smiles`.
- `ligands_3d`: DataStream = None — Per-id ligand coordinate files passed directly.
- `ligands_smiles`: DataStream = None — Optional `smiles` column for bond-order templating, paired with `ligands_3d`.
- `skip_ligand_embedding`: bool = False — Use the faster GEMS18e variant (no ChemBERTa ligand embedding, less accurate).

**Tables**:
- `affinity`: | id | structures.id | ligands.id | pkd_pred |
- `missing`: | id | cause |

**Example**:
```python
from biopipelines.gems import GEMS

# complex already holds the bound ligand pose
aff = GEMS(structures=complex, ligands=complex)
aff.tables.affinity
```

---

### OpenMM

Energy-minimises protein structures (Amber14 + implicit GBn2 solvent) to relieve clashes and bad geometry before downstream metric calculation. Scope is intentionally narrow — minimisation only, no trajectory production.

**References**: https://github.com/openmm/openmm

**Environment**: `openmm`

**Parameters**:
- `structures`: DataStream | StandardizedOutput (required) — Input structures.
- `max_iterations`: int = 1000 — Minimizer iterations.
- `tolerance_kj_per_mol_nm`: float = 10.0 — Convergence tolerance.
- `forcefield`: str = "amber14-all" — Force field.
- `solvent`: str = "implicit-gbn2" — Solvent model.
- `platform`: str = "auto" — Compute platform (`"auto"`, `"CUDA"`, `"CPU"`, ...).
- `restraint_selection`: str = "" — Residues to position-restrain during minimization.
- `restraint_k`: float = 1000.0 — Restraint force constant (kJ/mol/nm²).

**Streams**: `structures` (minimized PDB per input)

**Tables**:
- `energies`: | id | energy_initial_kj_mol | energy_final_kj_mol | delta_kj_mol |

**Example**:
```python
from biopipelines.openmm import OpenMM

relaxed = OpenMM(structures=boltz, max_iterations=2000)
```

---

### P2Rank

Template-free, machine-learning ligand binding-site prediction. For each input structure, runs `prank predict` and reports the predicted pockets (ranked by ligandability) and per-residue scores. P2Rank scores and clusters points on the solvent-accessible surface; no ligand, template, or MSA is needed, and it runs on CPU.

**Environment**: `p2rank` (bioconda; bundles the `prank` launcher, needs Java 17+).

**Parameters**:
- `structures`: Union[DataStream, StandardizedOutput] (required) — Input PDB/CIF structures.
- `config`: str = `"default"` — P2Rank model preset (`prank predict -c <config>`). One of `"default"`, `"alphafold"` (tuned for AlphaFold models using pLDDT in the B-factor column), or `"conservation"` (needs precomputed conservation scores).
- `threads`: int = 1 — Number of CPU threads P2Rank uses (`prank predict -threads`).
- `visualizations`: bool = False — Emit P2Rank's PyMOL/visualization artefacts (`prank predict -visualizations`); off by default for speed.

**Streams**:
- `residues`: per-residue resi-csv (one `<id>.csv` per input), columns `id | chain | resi | resn | pocket_idx | score | probability`. Consumable by the Selection tool to retrieve pocket residues, e.g. `Selection.add(p2rank.streams.residues, include="probability>0.5")`.

**Tables**:
- `pockets`: `id | pocket_idx | rank | score | probability | n_residues | residues | center_x | center_y | center_z` — one row per predicted pocket; `residues` is a chain-aware selection string (e.g. `"A12+A45-47"`).
- `residues`: `id | chain | resi | resn | pocket_idx | score | probability` — one row per residue (`pocket_idx` 0 = unassigned). Same data as the `residues` stream, combined across all inputs.
- `summary`: `id | n_pockets | top_score | top_probability | top_residues` — `top_residues` is the rank-1 pocket's residue selection.
- `missing`: `id | cause`

**Example**:
```python
from biopipelines.p2rank import P2Rank
from biopipelines.selection import Selection

# Standard prediction
target = PDB("2SRC", convert="pdb")
pred = P2Rank(structures=target)
pred.tables.pockets   # ranked pockets
pred.tables.summary   # top pocket residues, usable as a selection downstream

# Threshold per-residue ligandability into a Selection
pocket_sel = Selection(structures=target, ops=[
    Selection.add(pred.streams.residues, include="probability>0.5"),
])

# On an AlphaFold model
af = AlphaFold(proteins=Sequence("MKT...", ids="poi"))
pred = P2Rank(structures=af, config="alphafold")
```

---

### PLM_Sol

Sequence-based protein solubility prediction. PLM_Sol embeds each sequence with ProtT5-XL and classifies it with a biLSTM/TextCNN head trained on the updated E. coli solubility dataset (UESolDS); it is among the strongest sequence-only solubility predictors. It needs no structure — only sequence — so it complements the structure-based aggregation scorer Aggrescan3D. The output `solubility` is a probability in [0,1] (higher = more soluble); `soluble` is the ≥ 0.5 binary call.

**Environment**: `plm_sol` (Python 3.8 + torch 2.0.1 + `bio-embeddings`)

**Installation**: `PLM_Sol.install()` clones the repo (for its inference code and committed model checkpoint) and creates a dedicated env. The ProtT5-XL (`prottrans_t5_xl_u50`) embedding weights download lazily on the first prediction.

**GPU note**: ProtT5-XL embedding generation is impractically slow on CPU — run PLM_Sol on a GPU runtime (Colab GPU / cluster GPU node).

**Parameters**:
- `sequences`: Union[DataStream, StandardizedOutput] (required) - Input protein sequences

**Tables**:
- `solubility` (sorted by solubility descending):

  | id | sequences.id | solubility | soluble |
  |----|--------------|------------|---------|

- `missing` (sequences PLM_Sol could not score):

  | id | removed_by | kind | cause |
  |----|------------|------|-------|

**Example**:
```python
from biopipelines.plm_sol import PLM_Sol

# Score solubility of designed sequences straight off the sequence axis
seqs = ProteinMPNN(structures=rfd, num_sequences=8)
sol = PLM_Sol(sequences=seqs)

sol.tables.solubility        # probability per sequence, most soluble first
```

**Reference**: Zhang et al. (2024) PLM_Sol: predicting protein solubility by benchmarking multiple protein language models with the updated *Escherichia coli* protein solubility dataset. *Briefings in Bioinformatics* 25, bbae404. https://github.com/Violet969/PLM_Sol

---

### PLIP

Protein interaction profiler. Reads complex PDBs and reports detected non-covalent interactions: hydrogen bonds, hydrophobic contacts, π-stacking, salt bridges, halogen bonds, water bridges, and metal complexes. Profiles protein–ligand, protein–peptide / protein–protein, and intra-chain interactions, selected by `mode`.

**References**: https://github.com/pharmai/plip

**Environment**: `plip` (on HPC runs via an apptainer container).

**Parameters**:
- `structures`: DataStream | StandardizedOutput (required) — Protein complexes (PDB).
- `ligand`: DataStream | StandardizedOutput = None — Optional compounds stream naming the ligand; if omitted PLIP detects ligands automatically. Only valid in `mode="ligand"`.
- `mode`: str = `"ligand"` — Analysis mode:
  - `"ligand"` (default): protein–ligand interactions. Each complex must contain a ligand in HETATM records.
  - `"peptide"`: protein–peptide / protein–protein interactions. The chain(s) named in `chains` (deposited as ATOM records) are treated as the peptide/partner "ligand" and profiled against the rest of the structure (PLIP's `--peptides`/`--inter`).
  - `"intra"`: intra-chain interactions within the single chain named in `chains` (PLIP's `--intra`). Slower than ligand mode on large structures.
- `chains`: str | List[str] = None — Chain ID(s) to treat as the peptide/partner (`mode="peptide"`) or the single chain to profile (`mode="intra"`). Required when `mode != "ligand"`; `"intra"` takes exactly one chain. Must be empty in `mode="ligand"`.
- `generate_pse`: bool = True — Emit a PyMOL session per complex (`sessions` stream).

**Streams**: `sessions` (PyMOL `.pse`, when `generate_pse=True`)

**Tables**:
- `interactions`: | id | ligand | interaction_type | residue | chain | ... |
- `summary`: | id | n_hbonds | n_hydrophobic | n_pi_stacking | ... |
- `missing`: | id | cause |

In `peptide`/`intra` mode the `ligand` column labels the peptide/partner binding site reported by PLIP rather than a HETATM code.

**Example**:
```python
from biopipelines.plip import PLIP

# Protein–ligand (default)
interactions = PLIP(structures=boltz_holo)
interactions.tables.summary  # interaction counts per complex

# Protein–peptide / protein–protein: chain I is the peptide partner
ppi = PLIP(structures=complex, mode="peptide", chains="I")

# Two partner chains profiled against the rest of the structure
ppi = PLIP(structures=complex, mode="peptide", chains=["B", "C"])

# Intra-chain interactions within chain A
intra = PLIP(structures=complex, mode="intra", chains="A")
```

---

### PoseBusters

Validates computationally generated molecule poses by checking bond lengths, bond angles, internal steric clashes, volume overlap with protein, and more. Supports `dock` mode (ligand + protein) and `redock` mode (+ reference ligand for RMSD comparison).

**Environment**: `posebusters`

**Parameters**:
- `structures`: Union[DataStream, StandardizedOutput] (required) - Protein-ligand complexes (PDB/CIF)
- `ligand`: Union[DataStream, StandardizedOutput] (required) - Compounds stream (`Ligand(code="LIG")` or any compounds-producing tool) naming the ligand. The residue `code` is read from the stream at runtime, and is reused as the reference-structure residue code in redock mode. **Provide a SMILES for non-trivial ligands.** The ligand is extracted from each complex and rebuilt as an SDF for validation. If the stream carries a `smiles`, it is used as a **bond-order template** (`AssignBondOrdersFromTemplate`); otherwise bond orders are perceived from coordinates alone (`rdDetermineBonds`). Coordinate-only perception **fails on charged/conjugated ligands** (dyes, Si-rhodamines, metallo-organics, zwitterions): the molecule won't build, so `mol_pred_loaded` fails and **every** check — including the coordinate-based distance/overlap checks — reports False, i.e. the whole step yields `all_pass=False`. For such ligands you **must** pass `Ligand(smiles=..., codes="LIG")` whose `codes` matches the residue code in the structures. A bare `Ligand(code="LIG")` (no SMILES) is fine only for simple ligands RDKit can perceive from coordinates, and additionally triggers a (failing) RCSB SMILES lookup unless the code is a real CCD code. Note: covalent-derived poses fail `minimum_distance_to_protein` by construction (the attachment atom sits at bonding distance); exclude that check via `exclude=` or `check=[...]` when validating them.
- `reference_ligand`: Union[DataStream, StandardizedOutput, None] = None - Reference ligand structure for redock mode
- `mode`: str = "dock" - Validation mode: 'dock' or 'redock' (auto-set to 'redock' if reference_ligand provided)
- `check`: Union[str, List[str]] = "pose" - Which checks contribute to `all_pass`. All check columns are still written to the CSV — excluded ones are placed to the right of `all_pass` so the audit trail is preserved.

  Presets:
  - `"pose"` (default): every check except `no_radicals`. That one check fails spuriously on predicted poses because RDKit's bond perception from coordinates flags unsaturated atoms as radicals.
  - `"all"`: every available check (PoseBusters default behaviour).
  - `"loading"`: only that the molecule loaded.
  - `"geometry"`: loading + bond lengths/angles/clashes/flatness/energy.
  - `"chemistry"`: loading + sanity + radicals + connectivity.

  List form: `check=["all_atoms_connected", "bond_lengths"]` to compute `all_pass` over a custom subset.

- `exclude`: str | List[str] = None - Check column(s) to drop from whatever `check` selected. Easier than re-listing every check to keep when you only need to omit one or two. E.g. `exclude="minimum_distance_to_protein"` keeps the `pose` preset but ignores the protein-distance check (which fails by construction on covalent-derived poses — the attachment atom sits at bonding distance).

**Tables**:
- `posebusters` (dock mode):

  | id | source_structure | mol_pred_loaded | mol_cond_loaded | sanitization | inchi_convertible | all_atoms_connected | bond_lengths | bond_angles | internal_steric_clash | aromatic_ring_flatness | non-aromatic_ring_non-flatness | double_bond_flatness | internal_energy | protein-ligand_maximum_distance | minimum_distance_to_protein | minimum_distance_to_organic_cofactors | minimum_distance_to_inorganic_cofactors | minimum_distance_to_waters | volume_overlap_with_protein | volume_overlap_with_organic_cofactors | volume_overlap_with_inorganic_cofactors | volume_overlap_with_waters | all_pass | no_radicals |
  |----|------------------|-----------------|-----------------|--------------|-------------------|---------------------|--------------|-------------|-----------------------|------------------------|--------------------------------|----------------------|-----------------|--------------------------------|-----------------------------|---------------------------------------|-----------------------------------------|----------------------------|-----------------------------|---------------------------------------|-----------------------------------------|----------------------------|----------|-------------|

  Columns to the left of `all_pass` are the checks selected by `check=` (default preset `"pose"` shown). Columns to the right are still computed but excluded from `all_pass`. In `redock` mode, additional columns `mol_true_loaded`, `molecular_formula`, `molecular_bonds`, `double_bond_stereochemistry`, `tetrahedral_chirality`, `stereochemistry_preserved`, and `rmsd_≤_2å` are also included.

**Output Columns**:
- `id`: Structure identifier (suffixed with `_lig1`, `_lig2` etc. when multiple ligand copies exist)
- `source_structure`: Path to input structure file
- Boolean check columns: Each PoseBusters test (True = pass, False = fail)
- `all_pass`: True only if all checks selected by `check=` pass

**Example**:
```python
from biopipelines.posebusters import PoseBusters

# Default — "pose" preset, drops electronics so all_pass reflects geometry/placement only
validation = PoseBusters(
    structures=boltz_holo,
    ligand=boltz_holo  # reads the LIG code Boltz2 assigned; or Ligand(code="LIG")
)

# Custom subset — all_pass over only these two columns
validation = PoseBusters(
    structures=boltz_holo,
    ligand=boltz_holo,
    check=["all_atoms_connected", "bond_lengths"]
)

# Reproduce the original PoseBusters strict behaviour
validation = PoseBusters(
    structures=boltz_holo,
    ligand=boltz_holo,  # reads the LIG code Boltz2 assigned
    check="all"
)

# Redock mode — validate against crystal reference
xrc = PDB("1ABC")
validation = PoseBusters(
    structures=boltz_holo,
    ligand=Ligand(code="ATP"),
    reference_ligand=xrc,
    mode="redock"
)
```

---

### PoseChange

Measures ligand pose distance between reference holo structure and sample structures. Calculates RMSD and geometric metrics to quantify how well designed structures reproduce known binding poses.

**Environment**: `ProteinEnv`

**Parameters**:
- `reference_structure`: Union[str, ToolOutput, StandardizedOutput] (required) - Reference holo structure (e.g., XRC structure)
- `sample_structures`: Union[str, List[str], ToolOutput, StandardizedOutput] (required) - Designed/predicted structures to compare
- `reference_ligand`: Union[Ligand, ToolOutput] (required) - The ligand in the reference structure. Pass a `Ligand` / any compounds-producing tool output; the residue `code` is read from the compounds stream's map_table at runtime (Ligand Contract).
- `sample_ligand`: Union[Ligand, ToolOutput, None] = None - The ligand in the sample structures (default: same as reference_ligand). Pass a `Ligand` / compounds stream.
- `reference_alignment`: Optional[str] = None - PyMOL selection for reference structure alignment (default: "not resn {reference_ligand}", built at runtime once the code is resolved)
- `target_alignment`: Optional[str] = None - PyMOL selection for target structure alignment (default: "not resn {sample_ligand}", built at runtime once the code is resolved)

**Tables**:
- `changes`:

  | id | target_structure | reference_structure | ligand_rmsd | centroid_distance | alignment_rmsd | num_ligand_atoms | reference_alignment | target_alignment |
  |----|------------------|---------------------|-------------|-------------------|----------------|------------------|---------------------|------------------|

**Output Columns**:
- `ligand_rmsd`: RMSD between ligand poses after protein alignment (Å)
- `centroid_distance`: Distance between ligand centroids (Å)
- `alignment_rmsd`: RMSD of protein alignment (Å)
- `num_ligand_atoms`: Number of atoms in ligand
- `reference_alignment`: PyMOL selection used for reference structure alignment
- `target_alignment`: PyMOL selection used for target structure alignment

**Example**:
```python
from biopipelines.pose_change import PoseChange
from biopipelines.pdb import PDB
from biopipelines.entities import Ligand
from biopipelines.boltz2 import Boltz2
from biopipelines.panda import Panda

# Compare designed structures to XRC reference
xrc = PDB(pdbs="4ufc", ids="reference")
atp = Ligand("ATP")
ethanol = Ligand(smiles="CCO", codes="LIG")
designed = Boltz2(proteins=sequences, ligands=ethanol)

pose_analysis = PoseChange(
    reference_structure=xrc,
    sample_structures=designed,
    reference_ligand=atp,        # code read from the compounds stream at runtime
    sample_ligand=designed,      # Boltz2's assigned code, read at runtime
    # By default aligns on "not resn <ref_code>" / "not resn <sample_code>".
    # Override individually if needed:
    # reference_alignment="chain A and not resn ATP",
    # target_alignment="chain A and not resn LIG"
)

# Filter structures with RMSD < 2.0 Å
good_poses = Panda(
    tables=pose_analysis.tables.changes,
    pool=designed,
    operations=[Panda.filter("ligand_rmsd < 2.0")]
)
```

---

### Prodigy

Binding-affinity prediction for **protein–protein** complexes. For each complex PDB, reports the predicted dissociation constant Kd and binding free energy ΔG between two chain groups, from interfacial contacts.

**References**: https://github.com/haddocking/prodigy

**Environment**: `prodigy`

**Parameters**:
- `structures`: DataStream | StandardizedOutput (required) — Protein–protein complex PDBs.
- `interface`: str = "A B" — The two chain groups forming the interface (e.g. `"A B"`, `"A,B C"`).
- `temperature`: float = 25.0 — Temperature (°C) for the Kd conversion.

**Tables**:
- `affinity`: | id | interface | delta_g_kcal_mol | kd_M |

**Example**:
```python
from biopipelines.prodigy import Prodigy

complex = PDB("1A22", convert="pdb")
aff = Prodigy(structures=complex, interface="A B")
```

---

### ProLIF

Protein–ligand interaction **fingerprints**. For each complex, computes the ProLIF interaction matrix and exports it in long form (one row per residue × interaction-type pair, with a binary `present` flag) so it joins cleanly with other tables.

**References**: https://github.com/chemosim-lab/ProLIF

**Environment**: `prolif`

**Parameters**:
- `structures`: DataStream | StandardizedOutput (required) — Complex PDBs.
- `ligand`: DataStream | StandardizedOutput (required) — Compounds stream naming the bound ligand (the `code` is read at runtime; the bound HETATM block is extracted and bond-orders applied from the SMILES template).

**Tables**:
- `fingerprints`: | id | residue | resn | chain | resnum | resi | interaction_type | present |

**Example**:
```python
from biopipelines.prolif import ProLIF
from biopipelines.entities import Ligand

fp = ProLIF(structures=boltz_holo, ligand=Ligand(code="LIG"))
```

---

### Reduce

Adds and optimises explicit hydrogens on protein and ligand atoms via the Richardson lab's `reduce` binary. Preserves protein residue topology and ligand HETATM records, so the output is suitable as input to interaction-fingerprint or MD-prep tools.

**References**: https://github.com/rlabduke/reduce

**Environment**: `reduce`

**Parameters**:
- `structures`: DataStream | StandardizedOutput (required) — Input structures.

**Streams**: `structures` (protonated PDB per input)

**Example**:
```python
from biopipelines.reduce import Reduce

protonated = Reduce(structures=boltz)
prolif = ProLIF(structures=protonated, ligand=Ligand(code="LIG"))
```

---

### RTMScore

Scores and ranks docking poses by residue–atom distance likelihood (graph transformers). Takes complex structures + a posed ligand; for each pose, auto-extracts the pocket within `cutoff` Å of the ligand and reports a score (higher = more likely a true binding pose).

> **The ligand must be posed inside the pocket** (same caveat as [GEMS](#gems)) — pass a ligand-producing tool output as `ligands=`, not a from-scratch SMILES.

**References**: https://github.com/sc8668/RTMScore

**Resources**: GPU.

**Environment**: `rtmscore`

**Parameters**:
- `structures`: DataStream | StandardizedOutput (required) — Complex PDBs (protein source; pocket auto-extracted around the ligand).
- `ligands`: StandardizedOutput = None — A ligand-producing tool output (SDF references to score + SMILES metadata). Mutually exclusive with `ligands_3d`/`ligands_smiles`.
- `ligands_3d`: DataStream = None — Per-id ligand coordinate files passed directly (each file may hold multiple poses → one row per pose).
- `ligands_smiles`: DataStream = None — Optional `smiles` column for bond-order templating, paired with `ligands_3d`.
- `cutoff`: float = 10.0 — Pocket cutoff (Å) around the reference ligand.
- `model`: str = "model1" — Checkpoint (`model1`..`model4`).

**Tables**:
- `scores`: | id | structures.id | ligands.id | pose | rtmscore |
- `missing`: | id | cause |

**Example**:
```python
from biopipelines.rtmscore import RTMScore

scored = RTMScore(structures=complex, ligands=complex)
scored.tables.scores
```

---

### SASA

Solvent-accessible surface area analysis. Computes the change in SASA of a ligand when bound to vs separated from its protein partner — `delta_SASA = SASA(ligand alone) − SASA(ligand in complex)`. A larger delta indicates more ligand surface buried by the binder, typically interpreted as tighter packing.

**Environment**: `ProteinEnv` (installed alongside `PyMOL`).

**Installation**: no extra step beyond `PyMOL.install()`; SASA reuses the same env.

**Parameters**:
- `structures`: Union[DataStream, StandardizedOutput] (required) — Protein-ligand complexes.
- `ligand`: Union[str, Ligand, ToolOutput] (required) — The bound ligand. Pass a `Ligand` / any compounds-producing tool output (e.g. a Boltz2 result); the residue `code` is read from the compounds stream's map_table at runtime (Ligand Contract). A bare code string (e.g. `"LIG"`, `"AMX"`, `":X:"`) is also accepted for back-compat.
- `dot_density`: int = 4 — PyMOL `get_area` dot density (1–4, higher = more accurate, slower).

**Tables**:
- `sasa`:

  | id | structure | sasa_ligand_alone | sasa_ligand_complex | delta_sasa |
  |----|-----------|-------------------|---------------------|------------|

**Example**:
```python
from biopipelines.sasa import SASA

complex = PDB("9RTM", convert="pdb")
sasa = SASA(structures=complex, ligand=complex)  # code read from the compounds stream
sasa.tables.sasa  # delta-SASA per structure
```

---

### ThermoMPNN

Predicted change in fold stability (ddG) upon point mutation, from structure. ThermoMPNN is a lightweight GNN built on the ProteinMPNN backbone, trained on the Megascale stability dataset. For each input structure it scores every position × 19 substitutions in a single forward pass — materially faster than physics-based ddG estimation. Negative `ddG_pred` = predicted stabilising (lower folding free energy). Complements VespaG: ThermoMPNN scores fold *stability* from structure, VespaG scores functional *fitness* from sequence.

**Environment**: `thermompnn`

**Installation**: `ThermoMPNN.install()` clones `Kuhlman-Lab/ThermoMPNN` (the default model checkpoint ships with the repo, so no weights download) and creates a torch env. Runs on CPU; a GPU speeds up large/many inputs.

**Parameters**:
- `structures`: Union[DataStream, StandardizedOutput] (required) — Input backbone/complex PDBs.
- `chain`: str = `"A"` — Chain to score.
- `mutations`: Union[str, (TableInfo, column)] = `""` — Optional. Empty (default) runs site-saturation over all positions of `chain`. Otherwise scores only the named point mutations, given as `+`-joined wildtype-position-mutant tokens (e.g. `"A42G+L50V"`, 1-indexed), or a table column reference whose per-structure cell holds such a string. In explicit mode the native saturation pass still runs and is then filtered to the requested mutants; any requested mutation absent from the output (e.g. a wildtype/position mismatch) is reported in `missing`.

**Tables**:
- `ddg`:

  | id | structures.id | chain | position | wildtype | mutation | ddG_pred |
  |----|---------------|-------|----------|----------|----------|----------|

  One row per scored mutation, sorted by `ddG_pred` ascending (most stabilising first).

- `missing`: `id | removed_by | cause`

**Example**:
```python
from biopipelines.thermompnn import ThermoMPNN
from biopipelines.panda import Panda

target = PDB("1AKI", convert="pdb")

# Full saturation landscape
ddg = ThermoMPNN(structures=target, chain="A")

# Keep only stabilising mutations
stabilising = Panda(
    tables=ddg.tables.ddg,
    operations=[Panda.filter("ddG_pred < 0")],
)

# Score a specific shortlist instead
ddg_subset = ThermoMPNN(structures=target, chain="A", mutations="A42G+L50V")
```

**Reference**: Dieckhaus et al. (2024) Transfer learning to leverage larger datasets for improved prediction of protein stability changes. *PNAS* 121, e2314853121. https://github.com/Kuhlman-Lab/ThermoMPNN

---

### VespaG

Zero-shot single-substitution fitness prediction from sequence. VespaG is a small expert-guided model distilled from ESM-2, leaderboard-competitive on ProteinGym while being MSA-free and fast. Higher score = predicted more tolerated/beneficial substitution. Complements ThermoMPNN (stability from structure) by scoring evolutionary/functional fitness from sequence alone.

**Environment**: `vespag`

**Installation**: `VespaG.install()` creates a dedicated env and pip-installs `vespag`. VespaG's own weights ship with the wheel; the ESM-2 (`esm2_t36_3B_UR50D`) embeddings it generates internally pull a **~11 GB weight download lazily on the first prediction**.

**GPU note**: embedding generation is impractically slow on CPU — run VespaG on a GPU runtime (Colab GPU / cluster GPU node).

**Parameters**:
- `sequences`: Union[DataStream, StandardizedOutput] (required) — Amino-acid sequences (uses the `sequences` stream / its content-bearing CSV).
- `mutations`: Union[str, (TableInfo, column)] = `""` — Optional. Empty (default) scores the full single-site saturation landscape of each sequence. Otherwise scores only the named point mutations, given as `+`-joined wildtype-position-mutant tokens (e.g. `"A42G+L50V"`, 1-indexed), or a table column reference whose per-sequence cell holds such a string. Explicit mutations are translated into a VespaG mutation file at runtime; any requested mutation absent from the output is reported in `missing`.

**Tables**:
- `fitness`:

  | id | sequences.id | position | wildtype | mutation | fitness |
  |----|--------------|----------|----------|----------|---------|

  One row per scored mutation, sorted by `fitness` descending (most tolerated/beneficial first).

- `missing`: `id | removed_by | cause`

**Example**:
```python
from biopipelines.vespag import VespaG
from biopipelines.sequence import Sequence
from biopipelines.panda import Panda

with Pipeline("Project", "Fitness", description="Zero-shot fitness scan"):
    Resources(gpu="1", memory="24GB", time="1:00:00")  # GPU for ESM-2 embeddings
    poi = Sequence("MKTAYIAKQR...", ids="poi")

    # Full saturation landscape
    fit = VespaG(sequences=poi)

    # Most-tolerated substitutions
    tolerated = Panda(
        tables=fit.tables.fitness,
        operations=[Panda.filter("fitness > 0"), Panda.sort("fitness", ascending=False)],
    )

# Score a specific shortlist instead
fit_subset = VespaG(sequences=poi, mutations="A42G+L50V")
```

**Reference**: Marquet et al. (2024) Expert-guided protein language models enable accurate and blazingly fast fitness prediction. *Bioinformatics* 40, btae621. https://github.com/JSchlensok/VespaG

---

### XTB

Semi-empirical GFN2-xTB interaction-energy scoring. For each complex, splits the structure into a ligand fragment (by 3-letter code) and the surrounding protein, runs single-point GFN2-xTB on protein, ligand, and complex, and reports `E_complex − E_protein − E_ligand` — a physics-grounded ranking signal complementary to ML scorers.

> A full GFN2-xTB single point on a whole protein (~1500 atoms) takes >1 h. For routine use, restrict to the binding pocket (a [DistanceSelector](#distanceselector) `within` selection) or use the cheaper `gfnff` method.

**References**: https://github.com/grimme-lab/xtb

**Environment**: `xtb`

**Parameters**:
- `structures`: DataStream | StandardizedOutput (required) — Complex PDBs.
- `ligand`: DataStream | StandardizedOutput (required) — Compounds stream naming the ligand fragment (3-letter `code`).
- `method`: str = "gfn2" — `"gfn2"` (accurate) or `"gfnff"` (fast force-field).
- `solvent`: str = None — Implicit solvent (e.g. `"water"`); `None` = gas phase.
- `charge`: int = 0 — Total system charge.
- `opt`: bool = False — Geometry-optimize before the single point.

**Tables**:
- `interaction_energies`: | id | e_complex_kj | e_protein_kj | e_ligand_kj | e_interaction_kj | ... |
- `missing`: | id | cause |

**Example**:
```python
from biopipelines.xtb import XTB
from biopipelines.entities import Ligand

# Restrict to the pocket for speed
e_int = XTB(structures=boltz_holo, ligand=Ligand(code="LIG"), method="gfn2")
```
