# Analysis

[← Back to Tool Reference](../tool_reference.md)

### ADMETAI

ADMET endpoint predictions for compound libraries via the ADMET-AI Chemprop-RDKit model. Reads SMILES from a compounds stream and writes one row per input compound with all upstream-reported ADMET properties (~40+ regression and classification scores covering absorption, distribution, metabolism, excretion, and toxicity).

**Tags**: predict-property, small-molecule

**Environment**: `admet_ai`

**Installation**: `ADMETAI.install()` creates a dedicated conda env (Python 3.10) and installs `admet-ai` via pip. Verification instantiates `ADMETModel`, which downloads the bundled weights and warms the cache. On Daint the env is a plain venv on the host Python 3.11 and `admet-ai` is pinned below 2.0 — see `skills/biopipelines/references/daint_backend.md`.

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

**Tags**: detect, predict-property, protein, pocket, residues, small-molecule

**Environment**: `af2bind` — its own conda env (`environments/af2bind.yaml`), built by `AF2BIND.install()`. The conda layer is only python + numpy<2 + pandas; the JAX/ColabDesign stack (jax 0.4.23 + dm-haiku 0.0.10 + chex 0.1.7 + optax 0.1.7 + colabdesign 1.1.1) is pip-installed as an explicit step, because that set is the last one where ColabDesign v1.1.1's haiku still finds `jax.linear_util`. On Colab, ColabDesign is installed into the base Python (which has JAX+GPU). AF2BIND does not depend on the AlphaFold tool's environment. The two share only the **`AlphaFoldParams` cache** (the dir whose `params/` holds `params_model_*.npz`, configured in `folders.cache`): AF2BIND downloads the AF2 network params only if that folder is empty, so a co-located AlphaFold install's weights are reused and nothing is re-downloaded. Only the AF2BIND linear-head weights are tool-specific and land under the `AF2BIND` folder.

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
- `missing`: `id | removed_by | kind | cause`

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

**Tags**: measure, protein, small-molecule, residues

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

**Tags**: measure, protein, energy

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
- `electrostatics`: | id | net_charge | n_basic | n_acidic | pI | mean_potential |
- `missing`: | id | removed_by | kind | cause |

**Example**:
```python
from biopipelines.apbs import APBS

target = PDB("4UFC", convert="pdb")
elec = APBS(structures=target, ph=7.4)
elec.tables.electrostatics
```

---

### BFactor

Reads the B-factor column per residue and summarizes it over any number of named selections.

The B-factor column is where a structure predictor puts its per-residue confidence: Boltz2, AlphaFold and ESMFold all write pLDDT there in [0,100], while an experimental structure carries a real temperature factor in A^2. Every tool in this catalog reports confidence as one number for the whole structure (`complex_plddt`, `plddt_mean`), which averages the part under scrutiny together with everything around it. On a 476-residue fusion whose designed linker is four residues, a connector the model cannot place at all moves the structure-level score by less than a percent. This reads the vector the predictor already wrote, so the region that differs between designs can be scored on its own.

**Tags**: measure, protein, residues

**Environment**: `biopipelines`

**Values are never rescaled.** A 0-100 pLDDT and a 0-100 temperature factor are indistinguishable by range, so any auto-normalization would corrupt one of them; the caller knows which produced the file. One value per residue: the CA atom's B-factor when the residue has a CA, otherwise the mean over that residue's atoms.

**Parameters**:
- `structures`: DataStream | StandardizedOutput (required) - PDB/mmCIF structures.
- `selections`: Dict[str, str | TableReference] = None - `{name: selection}`. Each selection is a framework selection string (`"A75-77+A274"`; an unqualified `"75-77"` when only one chain has residues 75-77; if several do, that id is recorded as a failure naming them; a part that is not a span is refused) or a table column reference (`rfd.tables.structures.designed`) resolved per structure id. One group of summary columns is emitted per name; a single-row table broadcasts to every structure. Names must be alphanumeric with underscores, since they prefix the summary columns.

**Streams**:
- `bfactors`: per-residue `resi-csv`, one `<id>.csv` per input, columns `id | chain | resi | icode | bfactor`. Consumable by [Selection](data_management.md#selection) and [Consensus](#consensus), e.g. `Selection.add(bf.streams.bfactors, include="bfactor>=70")` to keep the confidently placed residues, or `Consensus.mean("bfactor")` to average a position across an ensemble.

**Tables**:
- `summary`:

  | id | n_residues | all_mean | all_sd | <name>_mean | <name>_sd | <name>_min | <name>_n | <name>_delta |
  |----|------------|----------|--------|-------------|-----------|------------|----------|--------------|

  One group of five columns per named selection. `<name>_delta` is the selection's mean minus the whole structure's: negative means the region is less certain than the model around it, which is the comparison a structure-level average cannot express. `<name>_n` is how many residues the selection actually matched - read it before trusting the mean, since a selection that matches nothing reports an empty mean rather than a zero.

- `missing`: | id | removed_by | kind | cause |

**Units.** Whatever the input file carries. Boltz2, AlphaFold and ESMFold write pLDDT in 0-100, where higher is more confident (verified against Boltz2 output: CA values span 26.4-97.3). A crystallographic temperature factor occupies the same field and the same rough range, but *higher means more disordered*, so the sign of "good" flips with the source and the range cannot tell you which you have. There is no cross-source comparison to be made without knowing which produced the file.

**Example**:
```python
from biopipelines.bfactor import BFactor

# Confidence over a diffused linker, scored from the column that defined it.
rfd = RFdiffusion(pdb=posed, contigs="A3-76/0-4/B2-197/0-4/A324-523", length="471-476")
mpnn = ProteinMPNN(structures=rfd, redesigned=rfd.tables.structures.designed)
folded = Boltz2(proteins=mpnn)

bf = BFactor(structures=folded,
             selections={"linker": rfd.tables.structures.designed,
                         "motif": rfd.tables.structures.fixed})

# Rank designs by how well the predictor places the designed span itself.
best = Panda(tables=bf.tables.summary,
             operations=[Panda.sort("linker_mean", ascending=False), Panda.head(10)])
```

---

### BindingData

Experimentally measured protein–ligand binding affinities for a compounds stream, retrieved from the ChEMBL REST API and the BindingDB RESTful service. For each input compound the tool resolves the molecule in each configured source and collects its Ki/Kd/IC50/EC50 records against protein targets. No model, no GPU — two HTTP APIs.

Sources are queried independently and their rows concatenated, so a measurement curated by both databases appears twice, distinguished by the `source` column; deduplicate downstream with `Panda` if needed. Affinity values are reported in nM. ChEMBL records carry `pchembl_value` (-log10 of the molar value, comparable across measurement types); BindingDB glues its relation onto the value (`">30000"`), which the tool splits into the separate `relation` and `affinity_nm` columns.

Every input id yields at least one row: a compound queried successfully with no matching record keeps an all-NaN row so the `affinities` table stays a complete matrix over the input stream, and is additionally recorded in `missing` with `kind="filter"`.

**Tags**: fetch, small-molecule, protein, complex, binding

**Environment**: `biopipelines` (no installation needed).

**Parameters**:
- `compounds`: Union[DataStream, StandardizedOutput] (required) — Input compounds. SMILES are read from the `smiles` column of the compounds map_table.
- `sources`: Union[str, List[str]] = `"chembl"` — Which databases to query: `"chembl"`, `"bindingdb"`, or a list of both.
- `affinity_types`: List[str] = None — Measurement types to keep, any of `"Ki"`, `"Kd"`, `"IC50"`, `"EC50"`. None keeps all four.
- `match`: str = `"exact"` — How a compound is matched to database molecules: `"exact"` (ChEMBL flexmatch / BindingDB similarity 1.0), `"similarity"` (Tanimoto above `similarity`), or `"substructure"` (ChEMBL only; rejected at construction when `bindingdb` is in `sources`).
- `similarity`: float = 0.85 — Tanimoto cutoff in [0, 1] for `match="similarity"`. Ignored otherwise.
- `max_affinity`: float = None — Keep only records at or below this value in nM. None keeps every record. A record whose value cannot be parsed is excluded when a ceiling is set.
- `organism`: str = None — Restrict to targets whose organism field contains this string (e.g. `"Homo sapiens"`). None keeps every organism.
- `max_records`: int = 1000 — Cap applied **twice** on the ChEMBL path: to the number of database molecules a query resolves to, and to the activities fetched per resolved molecule. The worst case is therefore `max_records²` records and `max_records + 1` requests per compound, so leave it low for a `match="similarity"` run over a large library. There is no pagination, so anything past ChEMBL's first page (its own cap is 1000) is not fetched.

**Tables**:
- `affinities`:

  | id | smiles | source | source_molecule_id | target_id | target_name | organism | affinity_type | relation | affinity_nm | pchembl_value | assay_id | assay_description | assay_type | document_id | match_similarity | queried_utc | chembl_release |
  |----|--------|--------|--------------------|-----------|-------------|----------|---------------|----------|-------------|---------------|----------|-------------------|------------|-------------|------------------|-------------|----------------|

  One row per compound–target measurement. `target_id`, `pchembl_value`, `assay_*` and `document_id` are ChEMBL-only (BindingDB's API returns no identifiers for them); `match_similarity` is populated for `match="similarity"`.
  `queried_utc` and `chembl_release` record **when** the lookup ran and **which ChEMBL release answered**. Both tables are a snapshot of databases that grow with every release, so two runs months apart legitimately disagree; without these columns there is nothing in the output to tell them apart, and a number cannot be traced back to what produced it. `chembl_release` is blank when ChEMBL was not queried or its status endpoint declined.


- `targets`:

  | target_id | source | target_name | organism | n_compounds | n_records | best_affinity_nm | queried_utc | chembl_release |
  |-----------|--------|-------------|----------|-------------|-----------|------------------|-------------|----------------|

  One row per distinct target and source. Grouping keys on `target_id` where the source supplies one — distinct ChEMBL targets share a `pref_name` (human and sheep COX-1 both read "Prostaglandin G/H synthase 1"), so keying on the name would merge them. BindingDB returns no target id, leaving the name as its only key.

- `missing`: `id | removed_by | kind | cause` — compounds with no SMILES or no matching record, plus upstream-filtered ids.

**Example**:
```python
from biopipelines import BindingData, CompoundLibrary, Panda

with Pipeline("Project", "Affinities", description="Known binders of a library"):
    Resources(time="1:00:00")
    library = CompoundLibrary("my_library.csv")
    affinities = BindingData(
        compounds=library,
        sources=["chembl", "bindingdb"],
        max_affinity=1000,          # nM
        organism="Homo sapiens",
    )

    # Keep only sub-micromolar dissociation constants
    Panda(tables=affinities.tables.affinities,
          operations=[Panda.filter("affinity_type == 'Kd'")])
```

**References**: Zdrazil et al. (2024) The ChEMBL Database in 2023. *Nucleic Acids Res* 52, D1180. https://www.ebi.ac.uk/chembl/api/data — Gilson et al. (2016) BindingDB in 2015. *Nucleic Acids Res* 44, D1045. https://bindingdb.org/rwd/bind/BindingDBRESTfulAPI.jsp

---

### BioEmu

Emulates a protein's equilibrium structural ensemble from sequence with a generative model. For each input sequence, samples `num_samples` statistically independent conformers approximating the equilibrium distribution — a fast alternative to running MD. Emits per-conformer PDBs plus a compact trajectory.

**Tags**: sample-ensemble, protein, ensemble, flexibility, all-atom, sidechain

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
- `missing`: | id | removed_by | kind | cause |

**Example**:
```python
from biopipelines.bioemu import BioEmu

seqs = Sequence("YYDPETGTWY", ids="chignolin")
ens = BioEmu(sequences=seqs, num_samples=10)
```

---

### BioEmu2

Ligand-capable successor of BioEmu (microsoft/bioemu2). Samples `num_samples` conformers of a system that holds protein chains and, optionally, ligands, by calling `bioemu2 sample-equilibrium`, and splits the sampled trajectory into one PDB per conformer. Optional contact and covalent-bond constraints go into the system, so the same wrapper covers free-pose, restrained and covalent-anchor sampling.

**Runs and ids.** A run is one system. The `sequences` and `ligands` axes combine as in `Boltz2` and `Protenix`: by default every sequence is paired with every ligand (ids `<sequence>+<ligand>`), and `Bundle` puts several rows in one system (`Bundle(ligands)` is one system with all of the ligands, id `<ligand>+<ligand>+<sequence>`). Without `ligands` the ensemble is apo and the run takes the id of its sequence. Protein chains are A, B, ... in the order the sequences are given, and the ligands take the next letters, which is what `contacts` and `covalent_linkage` refer to. The conformers of a run are `<run>_1` to `<run>_<num_samples>`.

**Ligands.** A ligand follows the Ligand Contract: its CCD code is passed to bioemu2 as `ccd_code` when the compounds row has one (a `Ligand("ATP")` fetched from RCSB does), and its SMILES otherwise, exactly as `Boltz2` and `Protenix` choose. bioemu2 places the ligand from scratch. A row with neither (a `Ligand(codes="LIG")` that only names a residue) fails its run and is listed in `missing`, instead of the protein being sampled alone.

**Constraints.** `contacts` entries are upstream's own and are written into the system unchanged. `covalent_linkage` is the wrapper's short form of one upstream `bonds` entry: it names the protein atom by chain, residue and atom name and the atom of the first ligand, and the residue names come from the sequence and the ligand's CCD code. Because the bond names the ligand's residue, that ligand has to be given by CCD code.

**Alignments.** By default bioemu2 queries the ColabFold server for each new sequence and caches the result. With `msas` the supplied `.a3m` or Boltz-style `.csv` is written, for that step only, into the cache directory bioemu2 reads (`BIOEMU2_MSA_CACHE`), so no search is made for that sequence. The alignment is matched to the sequence by its `sequences.id`, and its first record has to be the sequence; a sequence with no alignment is searched as usual. A system with several protein chains cannot use `msas`: bioemu2 keys the cache by the whole set of chains.

The cache of the searched sequences is `~/.cache/bioemu2/colabfold_msa` unless `BIOEMU2_MSA_CACHE` is exported, which bioemu2 reads itself. Pointing it at storage that is not purged (on Daint the group store, since scratch is) keeps a search from being paid for again by the next job. A run with a supplied alignment uses the step's own directory instead, so a stale entry of the same sequence cannot take its place.

**Physicality.** bioemu2 checks every sampled frame (bond lengths, clashes, chirality, chain breaks, for the protein only) and reports it, but never removes a frame. All conformers are therefore in `structures` and `trajectories`, and the `physicality` table says which ones failed and the first check they failed. Filter on it with `Panda` to drop them.

**Failures.** A run that bioemu2 does not complete, or that the wrapper cannot build (an empty sequence, a ligand with no chemistry, a bond to a residue the chain does not have), keeps its row in `summary` with empty `n_samples` and `n_residues` (the counts of the other rows stay integers), has no entry in the stream maps, and gets a `failure` row in `missing` with the reason. When bioemu2 itself exits non-zero the reason is its exit status followed by the last line of its output (for a Python traceback, the exception), and the whole output is in the step log. The step exits non-zero only when no run was attempted or every run failed.

**Tags**: sample-ensemble, protein, small-molecule, complex, ensemble, flexibility, all-atom, covalent, msa

**References**: https://github.com/microsoft/bioemu2

**Resources**: GPU of compute capability 7.5 or newer with the default install, which is upstream's `cuda` extra (CUDA 13 kernels). A V100 needs upstream's `cuda12` extra instead; see its README.

**Environment**: `bioemu2` (dedicated env, Python 3.12, which upstream requires)

**Installation**: `BioEmu2.install()` builds the env from `environments/bioemu2.yaml` and installs `bioemu2[cuda]>=0.1.0rc0` from `environments/bioemu2.pip.txt`; the specifier names a pre-release because the release candidate is the only version PyPI has, and it admits that for bioemu2 alone. It then downloads the `bioemu2-preview` checkpoint of the `microsoft/bioemu2` Hugging Face repository into the Hugging Face cache, so the first run needs no download on the compute node; if that fails the install still succeeds and the first run downloads it. The install is verified by the Python version, importing `bioemu2`, `mdtraj` and `torch`, and `bioemu2 sample-equilibrium --help`. An env left by an earlier spec that ran an older Python is rebuilt.

**Daint**: the `daint` variant names `environments: BioEmu2: "bioemu2"` and an image under `edf:` (`<group>/edf/bioemu2.toml`), so the whole step runs in the image and activates a venv built on top of it, as `Boltz2` does with `Boltz2Env`. `BioEmu2.install()` builds that venv inside the image, and `environments/daint/provision_venv.sh bioemu2` builds it on a bare node from the repo's own spec. An image that already carries bioemu2 is used by setting the entry to `""`: nothing is activated, the helper runs on the image's Python (which needs `pandas`, `pyyaml` and `mdtraj`), and `BioEmu2.install()` is refused by name, as for `BioEmu`. A run on a GH200 was reported to work under BioEmu2 1.0 with the image supplying the interpreter; neither route has been run since the rewrite, and the recipe of the image is not in the repo. See `skills/biopipelines/references/daint_backend.md`.

**Runtime environment**:
- `BIOEMU2_SRUN_ENVIRONMENT`: an EDF path. Where the helper runs, the `bioemu2 sample-equilibrium` call becomes `srun --environment=<EDF> bioemu2 ...`, and the helper's environment is cut to what SLURM, enroot and bioemu2 read, because the host's CUDA, Triton and `PYTHONPATH` variables reach the image and break bioemu2's torch and triton init. It is for a helper that runs outside the image, whose Python then needs `mdtraj`, `pandas` and `pyyaml`. Unset (the default) nothing changes, and nothing of it is written into the generated script. Leave it unset when the whole step already runs in the image, which the `edf:` entry of the `daint` variant does on either route.
- `BIOEMU2_MSA_CACHE`: where bioemu2 keeps ColabFold alignments, see Alignments above.

**Parameters**:
- `sequences`: DataStream | StandardizedOutput (required) - protein sequences, one chain per row. `Bundle` puts several in one system
- `ligands`: DataStream | StandardizedOutput | None = None - compounds, by CCD code or SMILES. None for an apo ensemble. `Bundle` puts several in one system
- `num_samples`: int = 100 - conformers per run
- `batch_size`: int = 16 - samples denoised together per forward pass; it caps GPU memory at long sequence times large batch. With a `seed`, the first batch is the same at any batch size and the later ones are not, so changing it changes the conformers after the first batch
- `seed`: int | None = None - RNG seed; None gives every run its own. The same seed and batch size gave identical conformers on repeat
- `sampler`: str = "edm" - `edm` (unsteered), `edm_gd` (gradient guidance), `edm_fk` (Feynman-Kac resampling) or `edm_fkc` (Feynman-Kac corrector). The steered samplers use upstream's default physicality potentials
- `contacts`: list[dict] | None = None - upstream contact dicts, `{"token1": {"asym_id": "A", "res_id": 3}, "token2": {"asym_id": "B", "atom_name": "PG"}, "max_distance": 6.0}`. A token is a chain with a 1-based residue index or with a ligand atom name
- `covalent_linkage`: dict | None = None - `{"chain": "A", "position": 3, "protein_atom": "SG", "ligand_atom": "C1"}`. Requires `ligands`, given by CCD code
- `msas`: DataStream | StandardizedOutput | None = None - precomputed `.a3m` or Boltz-style `.csv` alignments, for single-chain systems. Formerly `precomputed_msa`, which still works and is deprecated. An alignment is matched to its sequence by id (maps written before the id rule included), else by sequence; a sequence with none is searched through ColabFold as usual and the log names it
- `physicality_filter`: str = "high_quality" - `high_quality`, `high_quality_per_residue` (violation counts per residue), `no_violations` or `no_filter` (no check, every frame `not_checked`). It only changes what the `physicality` table reports

**Streams**: `structures`, `trajectories`

`structures` holds one `<run>_<k>.pdb` per conformer. `trajectories` holds one `<run>.xtc` per run, and the `topology` column of its map names the `<run>_topology.cif` that goes with it. In a structure the protein chains come first and each ligand follows as the next chain letter, one residue numbered 1. A ligand is written as ATOM records, not HETATM, under its CCD code as the residue name, or `LIG` when it was given by SMILES, with atom names that bioemu2 assigns (`PG` of ATP, `C1`, `O1` of a SMILES ethanol). `Contacts` and `PoseBusters` read the ATP ligand of such a structure when given `Ligand("ATP")`. The physicality check does not cover the ligand, and `PoseBusters` failed three of four ATP conformers (bond lengths in two, sanitization in one).

**Tables**:
- `summary`:

  | id | sequences.id | ligands.id | n_samples | n_residues |
  |----|--------------|------------|-----------|------------|

  One row per run that entered, with `n_samples` and `n_residues` empty where the run failed. `n_residues` counts every residue of the sampled system, ligands included. `ligands.id` is empty for an apo run.

- `physicality`:

  | id | run.id | frame | status | reason |
  |----|--------|-------|--------|--------|

  One row per conformer, with `id` the structure id and `frame` its zero-based index in the trajectory. `status` is `passed`, `rejected` or `not_checked`, and `reason` is the first failed check, empty when none. Both are empty if bioemu2 wrote no report for the run.

- `missing`:

  | id | removed_by | kind | cause |
  |----|------------|------|-------|

**Example**:
```python
from biopipelines import BioEmu2, Ligand, Sequence

protein = Sequence("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ", ids="design_1")
atp = Ligand("ATP")

# Apo and holo ensembles, 100 conformers each
apo = BioEmu2(sequences=protein)
holo = BioEmu2(sequences=protein, ligands=atp, seed=1)

# A restraint between residue 3 of the protein and the gamma phosphorus of the ligand, under a steered sampler
steered = BioEmu2(sequences=protein, ligands=atp, sampler="edm_gd",
                  contacts=[{"token1": {"asym_id": "A", "res_id": 3},
                             "token2": {"asym_id": "B", "atom_name": "PG"}, "max_distance": 6.0}])
```

---

### CABSflex

Fast coarse-grained Monte Carlo simulation of protein structure flexibility. Produces conformational ensembles and per-residue RMSF profiles. Optionally rebuilds models to all-atom representation using MODELLER.

**WARNING**: MODELLER requires a license key (free for academics). Get one at https://salilab.org/modeller/registration.html and set `KEY_MODELLER` before running.

**Tags**: sample-ensemble, measure, protein, ensemble, residues, flexibility, all-atom

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

**Tags**: predict-property, protein, residues, solubility

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

- `missing`: `id | removed_by | kind | cause` — structures the scorer failed on, plus upstream-filtered ids.

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

### StructureCluster

Groups a structures stream into fold families and ranks the families, so a campaign of thousands of designs can be read as "which topologies came out, and which of them fold well" instead of a directory of models. Each structure is reduced to a length-normalized CA trace (resampled to `n_points` positions along the residue index) and compared by distance-matrix RMSD — no superposition is ever computed, which is what makes 10,000 inputs tractable in the base environment. dRMSD is mapped onto a bounded score with the TM-score normalization `1/(1+(dRMSD/d0)^2)`, so `threshold` reads on the familiar 0..1 scale and 0.5 is a sensible start.

This is a proxy for TM-score, not TM-align. dRMSD on a resampled trace is cheaper and stricter: it is sensitive to overall size and internal register in a way an optimal superposition is not, and it will not recognize two folds related by a large rigid-body domain motion. Tune `threshold` against the run's own `similarity_to_representative` distribution rather than transferring a value from the TM-score literature.

Clustering is sphere exclusion (leader / Taylor-Butina) in *ranking* order: structures are visited best-first, and each unassigned structure becomes a representative that absorbs everything within `threshold` of it. Two consequences worth knowing — each cluster's representative is its best-scoring member (the model you actually want to open), and a structure joins the *first* representative that captures it rather than its nearest one. Cost is O(N·n_clusters), not O(N²).

For clustering with the real TM-align aligner at fixed TM-score thresholds, or of whole complexes, use [Foldseek](#foldseek); it counts folds but does not rank them.

**Tags**: data, protein, ensemble

**Environment**: `biopipelines` (numpy + pandas + the shared PDB parser, which already reads mmCIF; no external tool, no dedicated env)

**Installation**: none — the tool runs in the base `biopipelines` environment.

**Parameters**:
- `structures`: DataStream | StandardizedOutput — predicted structures (PDB or mmCIF). Analyzes `chain`, or the longest protein chain when unset — which in a de novo pipeline is the designed chain, and keeps a HETATM ligand chain out of the trace automatically.
- `metrics`: TableInfo | StandardizedOutput | str | list = None — per-id table(s) carrying the `rank_by` columns, joined on `id` (e.g. `fold.tables.confidence`, `sasa.tables.sasa`). Without it, clusters are ranked by size.
- `rank_by`: list[str] = ["plddt", "iptm"] — columns to average per cluster and rank on, in priority order.
- `ascending`: bool | list[bool] = False — sort direction; the default (higher is better) is right for pLDDT, ipTM and buried surface alike.
- `threshold`: float = 0.5 — similarity cutoff in (0, 1] for joining a cluster. Higher means tighter, more numerous clusters.
- `n_points`: int = 64 — CA trace resampling length. Higher sharpens discrimination and costs O(n_points²) memory per structure.
- `chain`: str = "" — chain to analyze; unset takes the longest protein chain.
- `min_residues`: int = 20 — shorter chains are dropped from the clustering with `kind="filter"` (they keep their `assignments` row).
- `max_structures`: int = 0 — input ceiling; 0 disables it.

**Tables**:
- `assignments`: | id | cluster | cluster_rank | is_representative | is_medoid | similarity_to_representative | n_residues | chain | radius_of_gyration | helix_frac | strand_frac | coil_frac | relative_contact_order |
- `clusters`: | cluster | cluster_rank | size | fraction | representative | medoid | mean_n_residues | mean_helix_frac | mean_strand_frac | mean_coil_frac | mean_relative_contact_order | mean_radius_of_gyration |

Both tables carry further columns named after `rank_by`: `assignments` appends each ranked metric as-is, `clusters` appends the `mean_`, `median_`, `min_` and `max_` of each.
- `missing`: | id | removed_by | kind | cause |

Cluster labels are assigned *after* ranking, so `cluster_001` is the best family, and both tables keep a row for every entity that entered (unparseable structures carry NaN rather than vanishing).

**Topology descriptors**: `helix_frac` / `strand_frac` / `coil_frac` come from CA geometry alone (P-SEA-style i→i+2/3/4 distance criteria) and are an estimate, not DSSP — on crystallographic input they track it closely (86% helix for bacteriorhodopsin, 68% strand for a beta-sandwich design), but ~1 Å of per-atom backbone noise pulls them toward coil by 20–30 points, so read them as a comparison between models from the same predictor. Run `DSSP` on the representatives when a real hydrogen-bond assignment matters. `relative_contact_order` (mean sequence separation of CA contacts / chain length) is the noise-robust companion: low for local helical topologies, high for folds stitched together by long-range beta contacts.

**Example**:
```python
from biopipelines import ESMFold2, SASA, StructureCluster, Panda

folds = ESMFold2(proteins=seqs, ligands=dye)
burial = SASA(structures=folds, ligand=dye, mode="ligand")

families = StructureCluster(
    structures=folds,
    metrics=[folds.tables.confidence, burial.tables.sasa],
    rank_by=["plddt", "iptm", "delta_sasa"],
    threshold=0.5,
)

# one structure per family, best family first
reps = Panda(tables=families.tables.assignments,
             operations=[Panda.filter("is_representative == True"),
                         Panda.sort("cluster_rank")],
             pool=folds)
```

---

### Foldseek

Clusters structures at fixed TM-score thresholds with Foldseek and reports each structure's cluster and the cluster count of the set. It answers "how many distinct folds came out, and which designs share one" under a protocol you can quote: `foldseek easy-cluster` for single chains, `foldseek easy-multimercluster` for complexes, run once per threshold in `tm_thresholds`. Cluster counts describe diversity under that protocol; they are not a per-design quality score.

**Not StructureCluster.** [StructureCluster](#structurecluster) is a base-environment triage tool: a superposition-free distance-matrix proxy for TM-score, built for thousands of models, that ranks the families it finds on metrics you hand it. Foldseek runs the real aligner (TM-align by default) at published TM-score cutoffs, can cluster whole complexes, and ranks nothing. Use StructureCluster to decide which families to open; use Foldseek to report diversity.

**Tags**: measure, protein, complex

**Environment**: `foldseek` (bioconda binary; `Foldseek.install()` builds it from `environments/foldseek.yaml`). A BioContainers image is registered in `environments/_containers.yaml` and used when a `containers:` destination is set.

**TM-score threshold and coverage.** Foldseek keeps an alignment when its TM-score exceeds the threshold (`--tmscore-threshold`), normalized by the aligned region by default, and when it passes the E-value and coverage filters. Left unset, `evalue` and `coverage` are not passed and Foldseek applies its own defaults, which differ by mode (read from its source): single-chain `easy-cluster` uses `-e 0.01` and `-c 0.8` with cov-mode 0, so the alignment must cover 80% of both structures; `easy-multimercluster` uses `-e 10` and `-c 0.0`, no coverage requirement. Counts shift with these filters, so state them when you compare studies. Each threshold is a separate Foldseek run, so the cost grows with the number of thresholds.

**Below TM 0.5.** Below TM 0.5 the thresholds do nothing under the defaults, because the E-value and coverage filters act before the TM-score threshold does and structures of different folds never reach it. In a test set of ubiquitin (two copies), three T4 lysozyme variants, GFP and a 12-residue ubiquitin fragment, thresholds of 0.05, 0.3, 0.5 to 0.8, 0.95 and 0.99 all gave the same 4 clusters. With `evalue=10, coverage=0.0` the threshold 0.05 merged ubiquitin with GFP (3 clusters), while 0.3 and 0.5 still gave 4, so those two filters, not the TM-score, set the limit; the 12-residue fragment stayed alone in every run. `cov_mode` changed which structure represents a cluster (`cov_mode=1` chose the first copy and the first lysozyme) but not its members. In multimer mode, `coverage=0.8, cov_mode=2` merged a trypsin-ubiquitin complex with the ubiquitin monomer (5 clusters became 4), which the default, with no coverage requirement, kept apart.

**Parameters**:
- `structures`: Union[DataStream, StandardizedOutput] (required) — PDB or mmCIF structures.
- `tm_thresholds`: list[float] = [0.5, 0.6, 0.7, 0.8] — TM-score thresholds, each in (0, 1) with at most two decimals.
- `multimer`: bool = False — cluster complexes with `easy-multimercluster` instead of single chains with `easy-cluster`. In this mode `tm_thresholds` sets `--chain-tm-threshold` (every chain must align above it), and the whole complex must also exceed `multimer_tm_threshold` with an interface above `interface_lddt_threshold`. Chain names must not contain underscores. Without it every structure must hold one protein chain; a file with several chains is recorded as a failure instead of being assigned to one of them by guess.
- `alignment_type`: int = 1 — Foldseek `--alignment-type`: 1 is TM-align, 2 is the 3Di+AA local alignment followed by the TM-score filter.
- `min_seq_id`: float = 0.0 — Foldseek `--min-seq-id`, the minimum sequence identity to cluster. Single-chain mode only.
- `evalue`: float = None — Foldseek `-e`, the E-value above which an alignment is dropped. Must be positive. `None` passes no flag, so Foldseek's own default applies (0.01 for single chains, 10 for complexes).
- `coverage`: float = None — Foldseek `-c`, the fraction of residues the alignment must cover, read through `cov_mode`. In [0, 1]; 0 switches the requirement off. `None` passes no flag, so Foldseek's own default applies (0.8 for single chains, 0.0 for complexes).
- `cov_mode`: int = 0 — Foldseek `--cov-mode`, always passed. One of 0 to 5, as Foldseek defines them: 0 coverage of query and target, 1 coverage of target, 2 coverage of query, 3 target length at least `coverage` of query length, 4 query length at least `coverage` of target length, 5 the shorter sequence at least `coverage` of the other length. With `multimer=True` it also picks which complex TM-score must exceed `multimer_tm_threshold` (0 both, 1 target, 2 query), and only 0, 1 and 2 are accepted.
- `multimer_tm_threshold`: float = 0.65 — Foldseek `--multimer-tm-threshold`. `multimer=True` only.
- `interface_lddt_threshold`: float = 0.65 — Foldseek `--interface-lddt-threshold`. `multimer=True` only.
- `threads`: int = 1 — CPU threads Foldseek uses (`--threads`); match the CPUs requested with `Resources(cpus=...)`.

**Tables**:
- `tm_clusters`: `id | cluster_<t> | representative_<t> | n_clusters_<t>` for every threshold, where `<t>` is the percent form of the threshold (`tm50`, `tm60`, `tm70`, `tm80` by default). `cluster_<t>` is the cluster number (1 is the largest), `representative_<t>` the id of the cluster's representative structure, and `n_clusters_<t>` the number of clusters the whole set splits into (the same value on every row). One row per structure that entered; a structure Foldseek did not report carries NaN.
- `missing`: `id | removed_by | kind | cause`

The raw Foldseek cluster TSV of every threshold stays under the step's `_execution/foldseek_tm<NN>/` folder.

**Example**:
```python
from biopipelines import ESMFold2, Foldseek, Panda

folds = ESMFold2(proteins=seqs)
diversity = Foldseek(structures=folds, threads=8)

# one structure per cluster at TM 0.5
reps = Panda(tables=diversity.tables.tm_clusters,
             operations=[Panda.filter("id == representative_tm50")],
             pool=folds)

# stricter filters: alignments must cover 90% of both structures, E-value below 1e-3
strict = Foldseek(structures=folds, evalue=1e-3, coverage=0.9, cov_mode=0)

# no coverage requirement at all (single chains default to 0.8)
loose = Foldseek(structures=folds, coverage=0.0)

# complexes: every chain above TM 0.6, whole complex above 0.65
complexes = Foldseek(structures=folds, tm_thresholds=[0.6], multimer=True)

# complexes: at least 80% of the query complex must be covered by the alignment
covered = Foldseek(structures=folds, multimer=True, coverage=0.8, cov_mode=2)
```

The parameters are recorded in the step's `_configuration/config.yaml` (`evalue` and `coverage` are `null` when unset) and in the config display.

**Reference**: van Kempen et al. Fast and accurate protein structure search with Foldseek. *Nature Biotechnology*, doi:10.1038/s41587-023-01773-0; Barrio-Hernandez et al. Clustering predicted structures at the scale of the known protein universe. *Nature*, doi:10.1038/s41586-023-06510-w; Kim et al. Rapid and sensitive protein complex alignment with Foldseek-Multimer. *Nature Methods*, doi:10.1038/s41592-025-02593-7. https://github.com/steineggerlab/foldseek

---

### USalign

Aligns every pair of a table of structure pairs three ways with US-align: in sequence order, with a circular permutation, and non-sequentially. Each pair gives one row per way, with the TM-score normalized by either structure's length, the RMSD, the aligned length, the residue alignment as two gapped sequences and the superposition as a rotation matrix and a translation. It answers "does this design keep the fold of its parent, and does it keep it only once a circular permutation is allowed?": the sequential alignment of a circular permutant scores low although the fold is the same, while the permuted and non-sequential alignments do not.

**Not Foldseek or StructureCluster.** [Foldseek](#foldseek) and [StructureCluster](#structurecluster) group a whole set into families. USalign scores the pairs you name, however few, and keeps the alignment itself. For the sequential and permuted alignments alone, with the TM-score difference per pair, see [TMalign](#tmalign).

**Tags**: measure, protein

**Environment**: `usalign` (bioconda `usalign`, which ships both `USalign` and `TMalign`; `USalign.install()` builds it from `environments/usalign.yaml` and TMalign runs in the same env). bioconda also builds it for aarch64, so a `config.daint.yaml` entry exists, but the env has not been built on Daint yet.

**Pairs and ids.** The pair table is any table with one row per pair, such as a CSV loaded with `Table`. `structure_1` and `structure_2` name its two columns as references, `pairs.tables.data.<column>` or `(pairs.tables.data, "<column>")`, and must come from the same table. Both ids of a pair must be ids of `structures`. A pair with an id that is not there, or whose structure file an upstream tool removed, keeps its rows with NaN values and one `kind="failure"` row per alignment in `missing`. A pair listed twice is aligned once. `(a, b)` and `(b, a)` are different pairs: structure_1 is the one that is rotated onto structure_2, and the circular permutation is searched in structure_1.

**Alignments.** One row per pair and mode:

| `mode` | US-align flag | What it is |
|---|---|---|
| `sequential` | `-mm 0` | The standard alignment: aligned residues keep their order in both structures. |
| `circular_permutation` | `-mm 3` | Allows a circular permutation of structure_1. When none is found it is the sequential alignment, with identical numbers. |
| `non_sequential` | `-mm 5` (fNS) or `-mm 6` (sNS) | Pairs residues without regard to sequence order. `fNS` is fully non-sequential; `sNS` is semi-non-sequential and keeps the order within secondary-structure elements. |

Only the first chain of each structure is aligned, which is US-align's default for a pair of monomers. A more elaborate alignment is not a better one: the permuted or non-sequential TM-score can fall slightly below the sequential one.

**Reading the columns.** `tm_score_1` is normalized by `length_1` and `tm_score_2` by `length_2`; US-align recommends the one normalized by the reference structure, `tm_score_2`. `aligned_1` and `aligned_2` are the two aligned sequences, `-` marking a gap. In a permuted alignment a `*` in `aligned_1` marks the circular permutation point, and its position gives the cut. The rotation is `rotation_ij` (the matrix u) and `translation_i` (the vector t) of US-align: a coordinate x of structure_1 lands at `t + u x`.

**Splitting a large pair table.** The framework has no SLURM job arrays and nothing that splits a table. What it has: `with Parallel():` submits one job per iteration, `DataStream.chunks()` splits a stream, and `Pool` merges the siblings' tables. A table is split by row with `Panda.filter("index % N == k")`, where `index` is the row number of the pandas query; it works on a pair table that has no `id` column. Each slice then gets its own USalign, in the same job as the slice:

```python
pairs = Table("pairs.csv")                           # columns: parent, variant
runs = []
with Parallel():
    for k in range(8):
        Resources(time="04:00:00", memory="4GB")     # one sibling job per slice
        part = Panda(tables=pairs.tables.data, operations=[Panda.filter(f"index % 8 == {k}")])
        runs.append(USalign(structures=folds,
                            structure_1=part.tables.result.variant,
                            structure_2=part.tables.result.parent))
Resources(time="00:30:00", memory="4GB")
cp = Pool(runs=runs, keep_ids=True)                  # cp.tables.alignments, cp.tables.missing
```

`keep_ids=True` leaves the ids as USalign wrote them, which stay unique because a pair is in one slice; a pair repeated in the table would put it in two, and Pool then raises naming it. Without `keep_ids` each id gains `_<k>` and `original.id` holds the id before the suffix. Use at most as many slices as the table has rows: a slice with no pair fails its step.

**Parameters**:
- `structures`: Union[DataStream, StandardizedOutput] (required) — PDB or mmCIF structures; every id used in the pair table must be one of them.
- `structure_1`: TableReference | (TableInfo, str) (required) — column of the pair table with the id of the structure that is rotated.
- `structure_2`: TableReference | (TableInfo, str) (required) — column of the same table with the id of the structure it is rotated onto.
- `non_sequential_variant`: str = "fNS" — "fNS" (`-mm 5`, fully non-sequential) or "sNS" (`-mm 6`, semi-non-sequential).

**Tables**:
- `alignments`: | id | structure_1 | structure_2 | mode | length_1 | length_2 | aligned_length | rmsd | tm_score_1 | tm_score_2 | rotation_11 | rotation_12 | rotation_13 | rotation_21 | rotation_22 | rotation_23 | rotation_31 | rotation_32 | rotation_33 | translation_1 | translation_2 | translation_3 | aligned_1 | aligned_2 |
- `missing`: | id | removed_by | kind | cause |

`id` is `<structure_1>+<structure_2>_<mode>`. One row per distinct pair and mode; a pair or alignment that could not be run carries NaN in every number and empty aligned sequences, and a row of the same `id` in `missing`. The step fails only when every alignment failed. The tool filters nothing: apply `Panda.filter` to the table.

**Example**:
```python
from biopipelines import ESMFold2, Table, USalign, Panda

pairs = Table("pairs.csv")                           # one row per pair: parent, variant
folds = ESMFold2(proteins=seqs)

cp = USalign(structures=folds,
             structure_1=pairs.tables.data.variant,
             structure_2=pairs.tables.data.parent)

# variants whose permuted alignment to the parent is good
permuted = Panda(tables=cp.tables.alignments,
                 operations=[Panda.filter("mode == 'circular_permutation' and tm_score_2 > 0.7")])

# the semi-non-sequential variant instead of the fully non-sequential one
sns = USalign(structures=folds, structure_1=pairs.tables.data.variant,
              structure_2=pairs.tables.data.parent, non_sequential_variant="sNS")
```

**Reference**: Zhang et al. US-align: universal structure alignments of proteins, nucleic acids, and macromolecular complexes. *Nature Methods* 19, 1109-1115 (2022), doi:10.1038/s41592-022-01585-1. https://github.com/pylelab/USalign

---

### TMalign

Aligns every pair of a table of structure pairs in sequence order and with a circular permutation using TM-align, and reports for each pair how much the permutation gained. It is USalign without the non-sequential alignment, plus a per-pair `delta` table: use it to ask "which of these pairs are circular permutants of each other?" at the lowest cost, and [USalign](#usalign) when the non-sequential alignment is wanted too. The pair table, the id rules, the columns of `alignments` and the splitting of a large table across jobs are the same as for USalign; read them there.

**Tags**: measure, protein

**Environment**: `usalign` (the same env as USalign: the bioconda `usalign` package ships `TMalign` too, and the separate bioconda `tmalign` package claims the same `bin/TMalign`, so it is not installed beside). `TMalign.install()` builds it through `USalign`.

**Alignments.** One row per pair and mode, `mode` being `sequential` (plain `TMalign`) or `circular_permutation` (`TMalign -cp`, a circular permutation of structure_1). When TM-align finds no circular permutation the two rows hold identical numbers. Only the first chain of each structure is aligned.

**The TM-score difference.** `delta` holds, per pair, the TM-score of the permuted alignment minus the sequential one: `delta_tm_1` for the TM-scores normalized by structure_1, `delta_tm_2` for those normalized by structure_2. A large positive value says the pair is a circular permutation of one fold; 0 says TM-align found none. It is not always non-negative: TM-align adopts the permuted alignment on its own criterion, and the TM-score of that alignment can be slightly lower for one normalization (1UBQ against 1CRN gave -0.001 and -0.031).

**Parameters**:
- `structures`: Union[DataStream, StandardizedOutput] (required) — PDB or mmCIF structures; every id used in the pair table must be one of them.
- `structure_1`: TableReference | (TableInfo, str) (required) — column of the pair table with the id of the structure that is rotated.
- `structure_2`: TableReference | (TableInfo, str) (required) — column of the same table with the id of the structure it is rotated onto.

**Tables**:
- `alignments`: | id | structure_1 | structure_2 | mode | length_1 | length_2 | aligned_length | rmsd | tm_score_1 | tm_score_2 | rotation_11 | rotation_12 | rotation_13 | rotation_21 | rotation_22 | rotation_23 | rotation_31 | rotation_32 | rotation_33 | translation_1 | translation_2 | translation_3 | aligned_1 | aligned_2 |
- `delta`: | id | structure_1 | structure_2 | delta_tm_1 | delta_tm_2 |
- `missing`: | id | removed_by | kind | cause |

`alignments` ids are `<structure_1>+<structure_2>_<mode>`; `delta` has one row per distinct pair, with `id` `<structure_1>+<structure_2>`, and NaN where either alignment of the pair gave no TM-score. The ids in `missing` are those of the `alignments` rows that failed.

**Example**:
```python
from biopipelines import ESMFold2, Table, TMalign, Panda

pairs = Table("pairs.csv")                           # one row per pair: parent, variant
folds = ESMFold2(proteins=seqs)

cp = TMalign(structures=folds,
             structure_1=pairs.tables.data.variant,
             structure_2=pairs.tables.data.parent)

# pairs the circular permutation explains
permuted = Panda(tables=cp.tables.delta,
                 operations=[Panda.filter("delta_tm_2 > 0.1")])
```

**Reference**: Zhang and Skolnick. TM-align: a protein structure alignment algorithm based on the TM-score. *Nucleic Acids Research* 33, 2302-2309 (2005), doi:10.1093/nar/gki524; circular permutation as implemented in US-align, Zhang et al. *Nature Methods* 19, 1109-1115 (2022), doi:10.1038/s41592-022-01585-1. https://github.com/pylelab/USalign

---

### EnsembleAnalysis

Per-residue RMSF and ensemble-level metrics from a conformer ensemble. Where CABSflex couples RMSF to its own coarse-grained sampling, EnsembleAnalysis analyzes *any* ensemble: it superposes the conformers (least-squares on CA or backbone) and reports per-residue fluctuation, so it overlays RMSF profiles from NMR ensembles, PLACER dumps, BioEmu samples, or any pool of conformers on the same footing. The `rmsf` stream uses the same `resi-csv` schema as CABSflex, so `Selection` thresholds on it unchanged.

**Tags**: data, measure, protein, ensemble, residues, flexibility

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

- `missing` (only when an input axis carries an upstream manifest):

  | id | removed_by | kind | cause |
  |----|------------|------|-------|

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

**Tags**: measure, protein

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
- `pairing`: str = "sequence" - How atoms are put into correspondence:
  - `"sequence"` (default): align/super/cealign. Pairs by sequence similarity, so residues it cannot match are **excluded from the RMSD**. Right for homologues.
  - `"ordered"`: `cmd.fit(matchmaker=-1)` — Nth atom to Nth atom. Right when the structures share numbering and atom order but **not** sequence, e.g. a design and the refold of an inverse-folded sequence.
  - `"identifier"`: `cmd.fit(matchmaker=0)` — pairs on chain/resi/**resn**/name, so mutated positions are dropped.
- `cycles`: int = 5 - PyMOL's outlier-rejection cycles. Each discards the worst-fitting pairs and refits, so the reported RMSD describes only the survivors. **Set 0 to measure every atom.**
- `cutoff`: float = 2.0 - Rejection threshold (Å) for those cycles.
- `frame`: Optional[str | (TableInfo, "column")] = None - Superpose on this selection, then measure `selection` in that frame without refitting. Answers "given the cores are aligned, how far is the designed part from where it was designed?" A `(TableInfo, "column")` reference resolves per input structure at runtime, so each comparison can be framed on its own region.

  With a `frame`, the measured `selection` is compared where the fit left it. Under `pairing="sequence"` or `"ordered"` that comparison pairs atoms **by position**, so the selection must contain the same number of atoms in both structures — it raises if not, rather than reporting the `0.000 Å` that positional pairing over mismatched counts would otherwise produce. Use `pairing="identifier"` to pair by chain / residue / atom name instead when the two selections differ.

> **Defaults measure similarity, not fidelity.** For "did this fold as designed", use `pairing="ordered", cycles=0`. Measured on one real design whose segment is 11.43 Å from its design over all 200 backbone atoms: the defaults report **1.89 Å** (105 atoms, 29% dropped), `cycles=0` gives 5.69 Å (148 atoms — sequence pairing still excluded 52), `pairing="ordered"` gives the true 11.43 Å, and adding `frame="1-96"` gives **18.89 Å** — where the segment actually sits once the cores are aligned. `cealign` is not a way out: it reports the best common fragment path, 32 of 200 atoms, 3.75 Å.

**Tables**:
- `changes`:

  | id | reference_structure | target_structure | selection | num_aligned_atoms | RMSD | RMSD_before | num_atoms_before | num_residues_aligned | atoms_dropped_pct |
  |----|---------------------|------------------|-----------|-------------------|------|-------------|------------------|----------------------|-------------------|

  `RMSD` covers the atoms that survived refinement; `RMSD_before` / `num_atoms_before` are the same before rejection. A large `atoms_dropped_pct` means the RMSD is not describing the whole selection.

  **`atoms_dropped_pct` is only meaningful for `pairing="sequence"`.** That path goes through `cmd.align` / `cmd.super`, which report how many atoms their outlier-rejection cycles kept. `cmd.fit`, used by `pairing="ordered"` and `"identifier"`, reports no such count and does not narrow the selection, so `num_atoms_before` and `num_aligned_atoms` are necessarily equal there and `atoms_dropped_pct` is `0` by construction — not evidence that nothing was rejected. `cycles` still applies in those modes, and `RMSD` vs `RMSD_before` still shows its effect.

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

**Tags**: measure, protein, small-molecule, complex, pocket, interactions

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

  | id | source_structure | selections | ligand | {contact_metric_name} | min_distance | max_distance | mean_distance | sum_distances_sqrt_normalized |
  |----|------------------|------------|--------|-----------------------|--------------|--------------|---------------|-------------------------------|

  The contact-count column is named `contacts` unless `contact_metric_name=` renames it.

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

**Tags**: measure, protein, small-molecule, residues

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

**Tags**: measure, protein, small-molecule, pocket, residues

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

- `missing` (only when an input axis carries an upstream manifest):

  | id | removed_by | kind | cause |
  |----|------------|------|-------|

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

**Tags**: data, residues, pocket, ensemble, protein

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

**Tables**:
- `missing` (only when an input axis carries an upstream manifest): | id | removed_by | kind | cause |

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

**Tags**: measure, protein, residues

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
- `missing`: | id | removed_by | kind | cause |

**Example**:
```python
from biopipelines.dssp import DSSP

target = PDB("4UFC", convert="cif")
ss = DSSP(structures=target)
ss.tables.summary  # helix/sheet/coil fractions per structure
```

**Malformed refinement records are handled for you.** `mkdssp` 4.6.1 converts PDB input to mmCIF internally, and for some entries that conversion emits two `refine` rows with the same key, which its validator rejects: `Duplicate Key violation, cat: refine`. When that happens the runner retries with `REMARK 3` stripped and notes it in the log. `REMARK 3` is refinement metadata and contributes nothing to a secondary-structure assignment — on 4UFC the stripped result is byte-for-byte identical to feeding the same entry as mmCIF. So `convert="pdb"` and `convert="cif"` both work; neither needs special handling.

Note that `parse error at line 1: This file does not seem to be an mmCIF file` appears on **every** PDB input, including successful ones — it is `mkdssp` probing formats before falling back to its PDB reader. Verified on 1UBQ, 168L and 4HHB, which all print it and produce full output. It is never the cause of a failure.

---

### Embedding

One fixed-length vector per input from a registered embedding model, and the cosine similarity of the pairs you ask for. The tool is model-agnostic: every model is one entry of a registry in `biopipelines/embedding.py`, and today the two entries read structures with the Progres family of networks (Progres, CIRPIN), which embed the C-alpha graph of a protein into 128 numbers so that two structures of the same fold land close together, and CIRPIN's training makes that hold across circular permutations. A model of another input kind, for instance a sequence embedder, is one more registry entry plus an adapter class, not a new tool. It sits in Analysis next to ESM2, whose embeddings it complements, rather than in Sequence Statistics, whose tools work on frequency and correlation tables of mutations.

**Tags**: predict-property, protein, residues

**Environment**: `Embedding` (Python 3.12, `torch` 2.8.0+cu128 with the matching `torch-scatter` wheel pinned as a pair, `torch-geometric`, `progres`, `scipy`; the rest unpinned). It needs no GPU: the network has 12 MB of weights and the default `device="cpu"` is upstream's own. Not configured on Daint: the repository's config marks PyG wheels as unavailable on aarch64, and no run there exists.

**Installation**: `Embedding.install()` builds the env and checks that `torch`, `torch_scatter`, `torch_geometric`, `einops` and `progres.progres` import. It downloads no weights, and importing `progres.progres` downloads none either. `force_reinstall=True` rebuilds the env.

**Parameters**:
- `inputs`: DataStream | StandardizedOutput (required) - the structures to embed, PDB or mmCIF. From a standardized output the `structures` stream is read (the stream the registry entry names as its `input_kind`). Only the first chain of the first model is read, from the C-alpha atoms of its `ATOM` records; `HETATM` records, further chains and further models are ignored without a warning, so give one chain per id. A structure with fewer than 3 C-alpha atoms fails its id. The reader guesses the format from the extension (`.cif` or `.mmcif` is read as mmCIF, anything else as PDB), so a structure file must carry its real extension.
- `model`: str (required) - `"progres"` or `"cirpin"`. Required, with no default, because the two share one network and the name is what records which weights you meant.
- `weights`: str (required) - path to the checkpoint file of that model: a torch file holding a dict whose `model` key is the state dict. It must exist when the pipeline is configured, and a URL or a directory is refused. The step exports it as `TRAINED_MODEL_FP` and prints the first 12 digits of its SHA-256 in the log, so the log names the weights that ran. No default is downloaded. The file is read with `weights_only=True`.
- `first`, `second`: column references (default `None`) - the two members of each pair to compare, as columns of one table: `first=pairs.tables.data.query, second=pairs.tables.data.target` (or `(table, "column")` tuples). Both or neither, and both from the same table. Ids are compared as text with the ids of `inputs`; a pair is one table row, in the table's order, repeats included. Omit both for no `similarities` table.
- `node_embeddings`: bool = False - also write the `node_embeddings` stream. For the three registered models it is the output of the per-node decoder (`node_dec`) before the sum over residues, not normalized.
- `device`: str = `"cpu"` - `"cpu"`, `"cuda"` (fails the step when no GPU is visible rather than falling back) or `"auto"`.

**Registry contract.** The registry is `EMBEDDING_MODELS` in `biopipelines/embedding.py`; nothing else in the wrapper or the pipe script names a model. Each entry is a frozen `EmbeddingModel` declaring:

| Field | Meaning | The three registered entries |
|---|---|---|
| `name` | value of `model` | `progres`, `cirpin` |
| `input_kind` | stream the model reads | `structures` |
| `formats` | file formats accepted on that stream | `pdb`, `cif`, `mmcif` |
| `adapter` | class in `pipe_scripts/pipe_embedding.py` (`ADAPTERS`) that builds the model from the weights and embeds one input | `progres_graph` |
| `weights_env` | environment variable the script sets to the weights path; the adapter reads it | `TRAINED_MODEL_FP` |
| `upstream_default` | where upstream ships default weights, shown when `weights` is wrong; empty when none is named | Progres: Zenodo download and a copy in the CIRPIN repository; CIRPIN: its repository |
| `dimension` | length of one embedding vector; another length fails that id | 128 |
| `node_module` | submodule whose output is the per-node embedding; empty means `node_embeddings=True` is refused | `node_dec` |
| `extra_args` | keyword arguments the adapter passes to the upstream reader | `fileformat="guess"` |

Adding a model is one `EmbeddingModel` in the registry and, when its network or reader differs, one adapter class (`embed(path) -> (vector, nodes)`) in the pipe script. No adapter exists yet for a model that reads sequences, and none is built until one can be tested.

**How each model is selected.**
- `progres` and `cirpin` are selected by the weights file alone. They are the same `Model` of the `progres` package: the CIRPIN repository's copy of the Progres checkpoint (`Progres_trained_model.pt`) and its CIRPIN checkpoint (`CIRPIN_model_5k_cp_epoch301.pt`) both hold a dict with the keys `epoch` and `model`, and both load strictly into that `Model`. The Zenodo `trained_model.pt` that the pip package downloads was not fetched, so its identity with the repository's copy is not checked. CIRPIN differs in how it was trained (synthetic circular permutations of the training chains), not in the network or the graph features.
- The CIRPIN repository also ships a modified copy of `progres.py`, which reads `TRAINED_MODEL_FP` at import, lets `load_trained_model` take a path, skips the data download, and returns per-node embeddings through `embed_nodes`. This tool does not use that file. It uses the pip package's `Model`, `read_coords` and `coords_to_graph` and loads the checkpoint itself, because the pip package's `load_trained_model(device)` takes no path and first downloads its default weights and databases (about 850 MB, per its README), which a compute node without internet cannot do. The step still exports `TRAINED_MODEL_FP`, the variable the modified copy reads, so both agree on the weights.

**Streams**:
- `embeddings` - `<id>.npy` per input (format `npy`), float32, shape (128,), L2-normalized by the network, so the dot product of two vectors is their cosine.
- `node_embeddings` - `<id>.npy` per input (format `npy`), float32, shape (residues, dimension). Empty unless `node_embeddings=True`.

**Tables**:
- `similarities` (absent without `first` and `second`):

  | id | first | second | cosine_similarity |
  |----|-------|--------|-------------------|

  One row per row of the pair table; `id` is `<first>+<second>`. `cosine_similarity` is NaN when either member has no embedding, and the pair is then a `kind="failure"` row in `missing` under its pair id, with the cause naming the member and whether it is not in the input stream, was dropped upstream, or failed to embed.
- `missing`: | id | removed_by | kind | cause | - an input that failed to embed (no C-alpha atoms, fewer than 3, a non-finite vector) under its own id, and an uncomparable pair under its pair id; both `kind="failure"`.

**Cosine versus the Progres score.** `cosine_similarity` is the plain cosine in [-1, 1]. Progres's own score is `(1 + cosine) / 2`, running 0 to 1, and its search threshold (0.8 by default) is on that scale; convert with `Panda.calculate({"progres_score": "(1 + cosine_similarity) / 2"})`.

**Example**:
```python
from biopipelines.embedding import Embedding

with Pipeline("Project", "FoldSearch", description="CIRPIN embeddings of designs"):
    Resources(gpu="none", memory="8GB", time="1:00:00")
    designs = Boltz2(proteins=sequences)
    pairs = Table("pairs.csv")   # columns: query, target
    vectors = Embedding(
        designs, model="cirpin",
        weights="/data/weights/CIRPIN_model_5k_cp_epoch301.pt",
        first=pairs.tables.data.query, second=pairs.tables.data.target,
        node_embeddings=True,
    )
    scored = Panda(tables=vectors.tables.similarities,
                   operations=[Panda.calculate({"progres_score": "(1 + cosine_similarity) / 2"})])
```

**Verified and not verified.** Run on a Windows CPU with the real `progres` 1.1.0 package (torch 2.11, torch_geometric 2.8.1, einops 0.8.2) and the real Progres and CIRPIN checkpoints from the CIRPIN repository, on synthetic C-alpha traces: both checkpoints load strictly and under `weights_only=True`; each embedding is (128,) with norm 1; the first chain alone is embedded for a two-chain file; the node embeddings are (residues, 128) and `graph_dec(sum of nodes)`, normalized, reproduces the embedding to 1e-7 (so the hooked tensor is the one the network pools); the whole pipe script runs on the real weights, pairs included. On a synthetic helix and the same helix circularly permuted by 15 residues the cosine is 0.994 with CIRPIN and 0.680 with Progres, as the CIRPIN training predicts, but that is one synthetic example, not a benchmark. `torch_scatter` had no wheel on that machine and was replaced by a ten-line sum-scatter, which is what the network uses it for; the compiled `torch_scatter`, a GPU, the cluster, Colab and container installs, and real protein structures have not been run. The embeddings of real structures were not compared against `progres embed` output.

**References**: Progres: Greener & Jamali (2025) Fast protein structure searching using structure graph embeddings. *Bioinformatics Advances* 5, vbaf042. https://github.com/greener-group/progres. CIRPIN: Kolodziej, Abulnaga & Ovchinnikov (2025) CIRPIN: Learning Circular Permutation-Invariant Representations to Uncover Putative Protein Homologs. *bioRxiv* 10.1101/2025.11.18.689110. https://github.com/aidenkoloj/CIRPIN

---

### ESM2

Pseudo-log-likelihood and mean-pooled embeddings of protein sequences from the ESM-2 language model. It runs through Hugging Face `transformers` on the released `facebook/esm2_t*_UR50D` checkpoints (MIT license). Per sequence the tool writes how probable the model finds it (`pll`, `pll_per_residue`) and a mean-pooled embedding vector, so a pipeline can rank designs by naturalness or feed the vectors to its own downstream analysis. VespaG and PLM_Sol already use ESM-2 and ProtT5 internally but expose neither; this tool exposes the model itself. It sits with them in Analysis rather than in Sequence Statistics, whose tools work on frequency and correlation tables of mutations, not on sequences.

**Tags**: predict-property, protein, sequence-only

**Environment**: `esm2` (Python 3.12, `torch`, `transformers`, unpinned). GPU is optional: small checkpoints on short sequences are practical on CPU, while the 3B checkpoint and long sequences want a GPU (the cost of `pll` grows with the square of the length times the model size). Only the 8M checkpoint has been run, on CPU.

**Installation**: `ESM2.install()` creates the env and downloads the default checkpoint; `ESM2.install(models=["esm2_t6_8M_UR50D", "esm2_t36_3B_UR50D"])` downloads the ones named. Weights go to `folders.repositories.ESM2`, one subfolder per checkpoint. A Hub id follows its `main` branch, which can change between submitting and running, so each checkpoint is resolved to a commit hash at install time, downloaded as that commit, and `<repo>@<commit>` is written to `REVISION` beside the files. Runs load that local folder only (offline), print the revision in the step log, and a checkpoint that already has a `REVISION` is not downloaded again; delete its folder to take a newer commit. The 8M to 650M repositories ship safetensors; the 3B one ships two pickle shards (about 11 GB), whose loading needs torch >= 2.6 (transformers refuses `torch.load` below it). `force_reinstall=True` rebuilds the env and keeps the weights.

**Parameters**:
- `sequences`: DataStream | StandardizedOutput (required) - protein sequences, one chain per id. Upper-case letters of the ESM-2 alphabet only (the 20 standard amino acids and `X`, `B`, `U`, `Z`, `O`); any other character fails that id and names the characters. Sequences are never truncated.
- `model`: str = `"esm2_t33_650M_UR50D"` - one of `esm2_t6_8M_UR50D`, `esm2_t12_35M_UR50D`, `esm2_t30_150M_UR50D`, `esm2_t33_650M_UR50D`, `esm2_t36_3B_UR50D` (hidden sizes 320, 480, 640, 1280, 2560). The default is the size of the model card's own example: 2.6 GB of weights, which fits any current GPU, and the largest checkpoint published as safetensors. `esm2_t48_15B_UR50D` is not offered: its 60 GB of fp32 weights do not fit one GPU, and this tool does not shard or run in half precision.
- `compute_pll`: bool = True - write the `scores` table. One forward pass per residue, so on a large library turn it off when only embeddings are wanted.
- `compute_embeddings`: bool = True - write the `embeddings` stream. One forward pass per sequence. At least one of the two must be on.
- `batch_size`: int = 16 - masked copies of one sequence scored per forward pass. Changes memory and speed, not the values.
- `device`: str = `"auto"` - `"auto"` uses a GPU when one is visible, otherwise the CPU; `"cpu"`; `"cuda"`, which fails the step when no GPU is visible instead of falling back.

**Streams**:
- `embeddings` - `<id>.npy` per sequence (format `npy`), the mean-pooled vector described below. Empty when `compute_embeddings=False`.

**Tables**:
- `scores` (absent when `compute_pll=False`):

  | id | length | pll | pll_per_residue |
  |----|--------|-----|-----------------|

  Exactly one row per id that entered; `length`, `pll` and `pll_per_residue` are NaN for an id that failed.

- `missing`: | id | removed_by | kind | cause |

**Pseudo-log-likelihood.** For a sequence x_1..x_L tokenized as `<cls> x_1 .. x_L <eos>`, residue i is replaced by `<mask>` while every other token stays, the model runs once per i, and log p(x_i | all other tokens) is read from the log-softmax over the model's full 33-token output vocabulary at that position, at the true token. `pll` is the sum over the L residues and `pll_per_residue` their mean; both are at most 0, and values closer to 0 mean the model finds the sequence more probable. The sum grows with length, so compare sequences of different lengths on `pll_per_residue`, which still depends on composition (a low-complexity repeat can score well). This is the pseudo-log-likelihood of Salazar et al. (2020), and the same sum the `pseudo-ppl` strategy of upstream's `examples/variant-prediction/predict.py` returns (it is named for a perplexity but sums log-probabilities and does not normalize by length). Computed in fp32 with the model in eval mode, so it is deterministic.

**Embeddings.** One forward pass on the unmasked sequence; the encoder's last hidden state (after its final layer norm) at the L residue positions, `<cls>` and `<eos>` excluded, averaged over positions into one float32 vector of the hidden size. fair-esm's `esm-extract --include mean` for the last layer computes the same average, except that it truncates at 1022 residues by default and this tool never truncates.

**Not offered: variant-effect scores.** Wild-type and masked marginals (Meier et al. 2021) score a substitution against a reference sequence, which a single input sequence does not carry, so there is nothing for them to be computed against. VespaG scores the single-substitution landscape of a sequence.

**Example**:
```python
from biopipelines.esm2 import ESM2

with Pipeline("Project", "Naturalness", description="ESM-2 scoring of designs"):
    Resources(gpu="1", memory="16GB", time="1:00:00")
    designs = ProteinMPNN(structures=backbones, num_sequences=8)
    esm = ESM2(sequences=designs)

    # Most natural-looking designs first
    ranked = Panda(tables=esm.tables.scores, operations=[Panda.sort("pll_per_residue", ascending=False)])

# Embeddings only, from a small checkpoint on CPU
vectors = ESM2(sequences=designs, model="esm2_t12_35M_UR50D", compute_pll=False, device="cpu")
```

**Verified and not verified.** Checked against the real `esm2_t6_8M_UR50D` checkpoint on CPU (Windows, torch 2.11, transformers 5.14): the batched pseudo-log-likelihood equals a naive loop of single-masked forward passes to 1e-6, the result does not depend on `batch_size`, and the pooled vector matches the encoder's last hidden state. Not yet run: the larger checkpoints, GPU, the 3B sharded pickle load, and the install on a cluster, Colab or container.

**Reference**: Lin et al. (2023) Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science* 379, 1123-1130. https://www.science.org/doi/10.1126/science.ade2574 - code: https://github.com/facebookresearch/esm - weights: https://huggingface.co/facebook/esm2_t33_650M_UR50D. Salazar et al. (2020) Masked Language Model Scoring. *ACL*, 2699-2712. Meier et al. (2021) Language models enable zero-shot prediction of the effects of mutations on protein function. *NeurIPS*.

---

### FPocket

Geometry-based binding-pocket detection. For each input PDB, runs the `fpocket` binary and reports detected pockets with their alpha-sphere count, volume, druggability score, and constituent residues. No ligand or model inference required.

**Tags**: detect, predict-property, protein, pocket

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
- `missing`: `id | removed_by | kind | cause`

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

**Tags**: predict-property, complex, small-molecule, binding

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

  One row per structure and ligand pair; `id` is `<structure>+<ligand>` (for example `prot1+lig1`), and `structures.id` and `ligands.id` name the two parents. A pair GEMS gave no prediction for keeps its row with a NaN `pkd_pred` and a `failure` row of `missing`; the step exits non-zero only when no pair got a prediction.

  **Units.** `pkd_pred` is a predicted pKd — that is −log10(Kd) with Kd in molar — so **higher is stronger**: 9 is nanomolar, 6 is micromolar. The opposite direction to Boltz2's `affinity_pred_value`.
- `missing`: | id | removed_by | kind | cause |

**Example**:
```python
from biopipelines.gems import GEMS

# complex already holds the bound ligand pose
aff = GEMS(structures=complex, ligands=complex)
aff.tables.affinity
```

---

### IPSAE

Scores the interface of a predicted complex with ipSAE, by running the upstream `ipsae.py` (Dunbrack, MIT) once per structure. ipSAE reads the predictor's PAE matrix and keeps only the residue pairs whose predicted aligned error is below a cutoff, which makes the score of an interface less sensitive to disordered regions and to how far a sequence was trimmed than AlphaFold's whole-complex ipTM. Upstream also reports ipTM recomputed from the PAE, pDockQ, pDockQ2 and LIS. The formulas are upstream's; this wrapper stages the inputs, runs the unmodified script and parses the score table it writes.

Each chain pair is scored in both directions (`asym`: chain 1's residues as the alignment reference, PAE of chain 2's residues) and as their maximum (`max`). The `ipsae` table summarizes the interface of one chain of interest against its partners from the `max` rows; `ipsae_pairs` keeps every row upstream printed for every chain of the structure.

**Tags**: measure, protein, complex

**Environment**: `biopipelines` (upstream needs numpy only; CPU).

**Installation**: `IPSAE.install()` downloads `ipsae.py` from the `main` branch of https://github.com/DunbrackLab/IPSAE, unpinned, into the folder `folders.repositories.IPSAE` of the active config, and checks that it parses and that numpy imports. Every run prints the sha256 of the script it used at the top of the step log.

**Parameters**:
- `structures`: DataStream | StandardizedOutput (required) - predicted complexes, PDB or mmCIF.
- `pae`: DataStream | StandardizedOutput (required) - per-id PAE matrices: a predictor output carrying a `pae` stream (`predictor.streams.pae`) or a DataStream of PAE files. Format `npz` (Boltz layout, key `pae`) or `json` (ColabFold or AlphaFold3 layout, key `pae` or `predicted_aligned_error`). A predictor that exposes no `pae` stream is refused at configuration time, naming the streams it does have. The file is paired with the structure by id.
- `binder_chain`: str = "A" - the chain whose interface the `ipsae` table summarizes. Must be a chain of the structures, otherwise that id is recorded as a failure listing the chains it has.
- `target_chains`: str | List[str] = None - chains scored against `binder_chain`. `None` means every other chain of the structure.
- `pae_cutoff`: float = 10.0 - residue pairs with PAE below this (A) enter the score.
- `dist_cutoff`: float = 10.0 - distance cutoff (A) upstream uses when counting interface residues (`dist1`, `dist2`).

The defaults of `binder_chain`, `target_chains` and `pae_cutoff` are those of the `ipsae` block of RosettaCommons Chaperone (https://github.com/RosettaCommons/chaperone); `dist_cutoff` is the 10 A of upstream's AlphaFold3 and Boltz usage examples. Chaperone's headline `ipsae` is the lowest per-target score (its `ipsae_max` is the highest) where this tool reports the highest; the per-target rows are in `ipsae_pairs` for either reading.

**Tables**:
- `ipsae`:

  | id | binder_chain | target_chain | ipsae | ipsae_d0chn | ipsae_d0dom | iptm_af | iptm_d0chn | pdockq | pdockq2 | lis |
  |----|--------------|--------------|-------|-------------|-------------|---------|------------|--------|---------|-----|

  Exactly one row per id that entered. The row is the `max`-type result of `binder_chain` against the target chain with the highest `ipsae` (the first listed wins a tie), and `target_chain` names that chain. All score columns are NaN for an id that failed.

- `ipsae_pairs`:

  | id | chain1 | chain2 | type | pae_cutoff | dist_cutoff | ipsae | ipsae_d0chn | ipsae_d0dom | iptm_af | iptm_d0chn | pdockq | pdockq2 | lis | n0res | n0chn | n0dom | d0res | d0chn | d0dom | nres1 | nres2 | dist1 | dist2 |
  |----|--------|--------|------|------------|-------------|-------|-------------|-------------|---------|------------|--------|---------|-----|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|

  Upstream's columns, lowercased: for each chain pair, one `asym` row per direction and one `max` row, so a structure of N chains has 3 N (N-1) / 2 rows. `ipsae` uses d0 from the number of chain 2 residues with PAE below the cutoff for each chain 1 residue; `ipsae_d0chn` and `ipsae_d0dom` use d0 from the summed chain lengths and from the residues of the pair with any PAE below the cutoff. `n0*` and `d0*` are the residue counts and d0 values behind them, `nres1` and `nres2` the residues of each chain with PAE below the cutoff to the other, and `dist1` and `dist2` those that are also within `dist_cutoff`. An id that failed has no rows here and its row in `ipsae` is NaN.

- `missing`: | id | removed_by | kind | cause |

An id is recorded in `missing` as a failure when its PAE file is absent, when `binder_chain` or a listed target chain is not in the structure, when the structure has no chain besides `binder_chain`, or when upstream raises or writes no score table (it exits 0 after printing a message when an input file is unusable, so the table, not the exit code, decides). Ids already dropped upstream are carried over.

**Values upstream cannot compute from the supplied PAE file are NaN.** `pdockq` and `pdockq2` need pLDDT, and `iptm_af` the predictor's own ipTM. An `npz` PAE file carries neither (upstream looks for `plddt_*.npz` and `confidence_*.json` beside it, which are not staged), so those three columns are NaN for every `npz`. A ColabFold `json` supplies them through its `plddt` and `iptm` keys; an AlphaFold3 `json` supplies pLDDT through `atom_plddts` but its ipTM lives in a separate summary file, so `iptm_af` is NaN there. Left alone, upstream prints a number computed from zeros in all these cases. `ipsae`, `ipsae_d0chn`, `ipsae_d0dom`, `iptm_d0chn` and `lis` come from the PAE matrix alone and are always present.

**What each score sees.** `ipsae`, `ipsae_d0chn`, `ipsae_d0dom`, `iptm_d0chn` and `lis` are functions of the PAE matrix and the chain lengths, not of the coordinates, so they do not know where the chains are. Chain B of a ColabFold complex moved 300 A away, with its json unchanged, keeps the ipSAE of the intact complex (0.925), and `dist1` and `dist2` fall to 0; a Boltz2 `npz` behaves the same (0.962 before and after). Only `pdockq` and `pdockq2`, which count contacts in the structure, fall (0.58 and 0.94 to 0.0 and 0.0 in that test), and they exist only for a `json` that carries pLDDT. A predictor's own PAE and coordinates agree, so this matters for a structure that was edited or docked after the PAE was written.

**OpenFold3 `json`.** `pdockq` and `pdockq2` are NaN for an OpenFold3 output even though its `*_confidences.json` holds a per-atom `plddt`. For an mmCIF structure upstream reads pLDDT only from `atom_plddts`, the AlphaFold3 key, and the wrapper leaves NaN what upstream could not read. `iptm_af` is NaN there too.

**Staging.** Upstream writes next to its input and builds file names by substituting substrings of the whole path, so each id is copied to `model.pdb` or `model.cif` plus `pae_model.<ext>` in `_execution/runs/<id>/`, the script runs there, and the copies are removed afterwards. Upstream's own outputs (`model_<pae>_<dist>.txt`, `_byres.txt`, `.pml`) stay in that folder; the `_byres.txt` file holds the per-residue scores. The structure atoms must match the PAE tokens the way upstream expects (one CA per residue; ligand atoms are skipped), otherwise upstream fails and the id is recorded as a failure with its message.

**Example**:
```python
from biopipelines.ipsae import IPSAE

with Pipeline("Project", "Interface", description="ipSAE of a predicted complex"):
    folded = Boltz2(proteins=designs, msas=msas)
    scores = IPSAE(structures=folded, pae=folded, binder_chain="A", target_chains="B")

# Keep the complexes whose interface is confidently placed
good = Panda(tables=scores.tables.ipsae, operations=[Panda.filter("ipsae >= 0.6")])
```

**Reference**: Dunbrack (2025) *Rēs ipSAE loquuntur*: What's wrong with AlphaFold's ipTM score and how to fix it. *bioRxiv*. https://www.biorxiv.org/content/10.1101/2025.02.10.637595 - code: https://github.com/DunbrackLab/IPSAE

---

### LigandAtomSelector

Distance-based residue selection referenced to a **subset of atoms inside one ligand residue**, rather than to the whole residue. [DistanceSelector](#distanceselector) measures from every atom of its reference; this measures only from the atom names you name.

The case it exists for is a chimeric ligand that occupies a single residue — a dye conjugated to a peptide, say, all written as one `LIG` — where you want the residues lining the peptide half and not those near the dye. Its output schema is identical to `DistanceSelector`, so the result drops into `LigandMPNN(redesigned=...)` unchanged.

**Tags**: measure, residues, small-molecule, pocket, protein

**Environment**: `biopipelines`

**Parameters**:
- `structures`: DataStream | StandardizedOutput (required) — Input structures.
- `ligand`: str (required) — Ligand residue name, e.g. `"LIG"`.
- `atoms`: str (required) — `+`-joined atom names within that residue to measure from, e.g. `"C61+C62+S57+O49"`. **Every name must exist in the ligand residue of every input structure**, or that structure is skipped — so a stream whose ligand atom naming is not uniform loses members silently apart from the `missing` row.
- `distance`: float = 5.0 — Cutoff in Å.
- `restrict_to`: str | (TableInfo, "column") | None = None — Restrict the search to a selection: a chain-aware string (`"A10-20+A30-40"`), a per-structure table column, or `None` for all protein residues.
- `include_reference`: bool = True — Kept for API symmetry with `DistanceSelector`. A ligand is not a protein residue, so this has no effect on the protein-residue output.

**Tables**:
- `selections`:

  | id | pdb | within | beyond | distance_cutoff | reference_ligand |
  |----|-----|--------|--------|-----------------|------------------|

  `within` and `beyond` are chain-aware selection strings (`"A12+A15-18"`), ready to hand to a design tool.

Atom names come from whatever wrote the structure. A ligand carved from a crystal file, one generated by RDKit, and one round-tripped through PDBQT can all name the same atoms differently — check the names in an actual input before relying on them.

**Example**:
```python
from biopipelines import LigandAtomSelector, LigandMPNN

# Redesign only the pocket around the glutathione half of a dye conjugate
pocket = LigandAtomSelector(structures=complexes, ligand="LIG",
                            atoms="C61+C62+S57+O49", distance=6.0)

designs = LigandMPNN(structures=complexes, ligand="LIG",
                     redesigned=pocket.tables.selections.within)
```

---

### OpenMM

Energy-minimises protein structures (Amber14 + implicit GBn2 solvent) to relieve clashes and bad geometry before downstream metric calculation. Scope is intentionally narrow — minimisation only, no trajectory production.

**Tags**: refine-structure, protein, small-molecule, complex, energy, covalent, sidechain, all-atom

**References**: https://github.com/openmm/openmm

**Environment**: `openmm` (includes the OpenFF/AmberTools stack the `ligand=` path needs, ~740 MB)

**Parameters**:
- `structures`: DataStream | StandardizedOutput (required) — Input structures.
- `max_iterations`: int = 1000 — Minimizer iterations.
- `tolerance_kj_per_mol_nm`: float = 10.0 — Convergence tolerance.
- `forcefield`: str = "amber14-all" — Force field.
- `solvent`: str = "implicit-gbn2" — Solvent model.
- `platform`: str = "auto" — Compute platform (`"auto"`, `"CUDA"`, `"CPU"`, ...).
- `restraint_selection`: str = "" — Residues to position-restrain during minimization.
- `restraint_k`: float = 1000.0 — Restraint force constant (kJ/mol/nm²).
- `ligand`: compounds stream = None — A small molecule bound in the input structures, given as a compounds stream (a `Ligand`, or any tool's `compounds` output). The ligand's HETATM block is carved out by residue name and its topology is built from the stream's `smiles` column via OpenFF — **not** perceived from the PDB coordinates, so bond orders are correct even in a structure with no CONECT records. Without this parameter a bound ligand is simply dropped from the system and the protein minimises as if the pocket were empty.
- `ligand_charge_method`: str = "am1bcc" — Partial-charge model, `"am1bcc"` or `"nagl"`. AM1-BCC is the reference method and takes minutes for a ~45-atom molecule; NAGL is a graph neural network that is far faster but covers a narrower element set (no Si, P, or metals). Charges are cached across structures on the canonical graph, so a many-pose run over one molecule computes them once. See the note below on mixing methods.
- `ligand_forcefield`: str = "gaff-2.11" — Small-molecule force field, `"gaff-2.11"` or `"openff-2.0.0"`.
- `covalent_anchor`: str = None — Protein atom name holding the ligand covalently, e.g. `"SG"`, or `"SG62"` to pin the residue number as well. The closest anchor/ligand-heavy-atom pair within `covalent_max_distance` is held by a stiff harmonic bond during minimisation. Read the caveat below before using this.
- `covalent_k`: float = 300000.0 — Restraint force constant, kJ/mol/nm².
- `covalent_length`: float = 0.18 — Restraint equilibrium length in **nm** (0.18 nm ≈ a C–S single bond).
- `covalent_max_distance`: float = 3.0 — How far the anchor may sit from the ligand and still count as bonded, in **Å**. Note the unit differs from `covalent_length`.
- `mobile_selection`: str | (TableInfo, "column") = "" — Chain-aware selection (`A10-20+B5`; a bare chain letter such as `B` selects the whole chain) whose side chains are free to move; every other protein atom is frozen and the ligand stays mobile. Backbone atoms of the selected residues are frozen too, so only rotamers relax. Empty (the default) means the whole structure is mobile.
- `frozen_selection`: str | (TableInfo, "column") = "" — Chain-aware selection in the same syntax (a bare chain letter freezes the whole chain), frozen outright, all atoms. Empty means nothing is frozen. Mutually exclusive with `mobile_selection`.

Frozen atoms have their mass set to zero, which makes them **immovable**, unlike `restraint_selection` which only applies a harmonic penalty and lets the whole structure drift. Both selections accept a per-structure `(TableInfo, "column")` reference, so a different region can be held for each input. A selection that cannot be parsed, an empty per-structure cell, or one that matches no protein atom is an error rather than a silent no-op: a literal is rejected at configuration time, and a structure whose selection fails is reported as failed instead of being minimised unconstrained.

**Streams**: `structures` (minimized PDB per input), `compounds` (only when `ligand=` was given: the input compounds stream, passed straight through)

The `compounds` stream is the *input* stream, unchanged — the minimiser reuses the input's residue codes rather than assigning its own, so the chemistry is identical and downstream tools keep it. Without `ligand=` the key is an empty placeholder and there is no compounds output at all.

**Tables**:
- `energies`: | id | energy_initial_kj_mol | energy_final_kj_mol | delta_kj_mol | n_mobile_atoms | n_frozen_atoms | charge_method |
- `missing` (only when an input axis carries an upstream manifest): | id | removed_by | kind | cause |

**Reading the energies.** `energy_final_kj_mol` and `delta_kj_mol` include any restraint and covalent-bond terms that were active, so they are **not comparable across structures minimised under different selections**. Compare them within one selection regime, or use them as a relative clash indicator rather than an absolute energy.

**Mixing charge methods.** With `ligand_charge_method="nagl"`, a molecule outside NAGL's element coverage falls back to AM1-BCC with a warning on stderr. That is deliberate — a silicon rhodamine or a phosphonate cannot get NAGL charges at all — but it means one run can produce structures parameterised two different ways, whose energies are not on the same scale. If that matters, set `"am1bcc"` explicitly rather than relying on the fallback.

**The covalent restraint is a tether, not a bond.** `covalent_anchor` adds a stiff harmonic bond to an already-built system; it does not rebuild the topology. The linked pair *is* excluded from the nonbonded terms, as a real bonded pair would be — without that the two atoms would be held a bond length apart while still repelling each other through full Lennard-Jones and Coulomb, and that clash would dominate the reported energies. What the tether still does **not** give you is the rest of a real bond: no angle or torsion terms cross the junction, so the ligand can pivot freely about the attachment point, and both molecules remain chemically saturated — the ligand keeps its hydrogens and the anchor cysteine keeps its HG, so the junction carries two more hydrogens and different formal charges than the real adduct. Use it to hold a covalent ligand in its pocket during minimisation; do not read the resulting energies as those of the adduct itself.

**Example**:
```python
from biopipelines import OpenMM, Ligand

relaxed = OpenMM(structures=boltz, max_iterations=2000)

# With a bound ligand, freezing everything but the pocket side chains
lig = Ligand(smiles="CC(=O)Oc1ccccc1C(=O)O", codes="LIG")
relaxed = OpenMM(structures=boltz, ligand=lig,
                 mobile_selection="A45+A48+A102",
                 ligand_charge_method="am1bcc")

# Per-structure pocket, read from an upstream selector's table
relaxed = OpenMM(structures=boltz, ligand=lig,
                 mobile_selection=(pocket.tables.selection, "within"))
```

---

### P2Rank

Template-free, machine-learning ligand binding-site prediction. For each input structure, runs `prank predict` and reports the predicted pockets (ranked by ligandability) and per-residue scores. P2Rank scores and clusters points on the solvent-accessible surface; no ligand, template, or MSA is needed, and it runs on CPU.

**Tags**: detect, pocket, residues, protein

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
- `missing`: `id | removed_by | kind | cause`

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

### Peppr

Scores predicted structures against reference structures with [peppr](https://github.com/aivant/peppr). The atoms of each prediction are matched to its reference, chain permutations and symmetric ligands included, and the CASP/CAPRI metrics are computed: backbone RMSD, TM-score, lDDT, lDDT-PLI, fnat, iRMSD, LRMSD and DockQ. Every value is a plain number per id; nothing is filtered on it, so apply `Panda.filter` downstream.

**Tags**: measure, protein, small-molecule, complex

**Environment**: `peppr`

**Installation**: `Peppr.install()` creates a dedicated conda env (Python 3.11 or newer, which `ProteinEnv` does not guarantee) and installs `peppr` through pip. `biotite` is kept below 2.0 in `environments/peppr.yaml`: peppr 0.13.1 imports a biotite class that 2.0 renamed, so an unbounded install fails at `import peppr`. The verification step imports peppr, so it catches this.

**Parameters**:
- `structures`: Union[DataStream, StandardizedOutput] (required) - Predicted structures (PDB or mmCIF).
- `references`: Union[DataStream, StandardizedOutput] (required) - Reference structures (PDB or mmCIF). A stream with one id is the reference of every prediction. A stream with several is matched to the predictions by id: the same id, or a parent or child of it (a refold `design_3_1` finds `design_3`). A sibling id is not a match, so a prediction whose reference is absent is not scored against another molecule's.
- `metrics`: Optional[List[str]] = None - Keys to compute, in column order. Default: all eight.

| key | peppr class | what the value is |
|-----|-------------|-------------------|
| `rmsd` | `MonomerRMSD` | Backbone RMSD (Cα, and C3' for nucleic acids) of each polymer chain after superposing it alone, averaged over chains weighted by atom count. Å. |
| `tm_score` | `MonomerTMScore` | TM-score of each protein chain, averaged the same way. |
| `lddt` | `MonomerLDDTScore` | lDDT of the contacts inside each polymer chain, averaged the same way. Contacts between chains are not included. |
| `lddt_pli` | `LDDTPLIScore` | CASP lDDT-PLI, the lDDT of the polymer-ligand interactions. NaN when no ligand touches the polymer in the reference. |
| `fnat` | `ContactFraction` | Fraction of the reference's contacts between two polymer chains that the prediction reproduces, averaged over chain pairs. |
| `irmsd` | `InterfaceRMSD` | Interface RMSD between two polymer chains, averaged over pairs. Å. |
| `lrmsd` | `LigandRMSD` | Backbone RMSD of the shorter chain of each polymer pair in contact, after superposing the longer one, averaged over pairs. Å. |
| `dockq` | `DockQScore` | DockQ of each chain pair, polymer-ligand pairs included, averaged over pairs. 0 to 1, higher is better. |

`rmsd`, `tm_score` and `lddt` score each chain on its own superposition, so they say nothing about how the chains sit relative to each other; `fnat`, `irmsd`, `lrmsd` and `dockq` do. The four pair metrics are NaN for a structure with one chain. peppr's own `MonomerRMSD` takes an RMSD threshold, but only to bin the summary it can print over many systems; it never changes the value, so there is no `rmsd_threshold` parameter here.

**Input requirements** (peppr's, not this wrapper's): each molecule, polymer chain or ligand, needs its own chain id, and bonds must be known. biotite reads them: polymer bonds come from the Chemical Component Dictionary, ligand bonds from `CONECT` records in a PDB file or from `chem_comp_bond` (else the dictionary) in an mmCIF file. A ligand is matched to its counterpart through its bonds and their orders, so a ligand whose file gives no usable bonds cannot be matched: its id lands in `missing`, with a `cause` that names the ligand and what its bonds lack. Tested with Boltz2 predictions of benzamidine against the 3PTB crystal: an mmCIF whose ligand was given as SMILES (residue `LIG1`) has no `chem_comp_bond` and a name the dictionary does not know, so the ligand has no bonds and the id fails; an mmCIF whose ligand was given by its CCD code (`Ligand("BEN")`, residue `BEN`) takes its bonds from the dictionary and is scored (`lddt_pli` 0.74, `dockq` 0.95); the PDB output of the same prediction carries `CONECT` records, which give bonds of unknown order, and fails too. A ligand in a PDB file without `CONECT` records has no bonds either. To score a predicted ligand, give it its CCD code or supply a file with `chem_comp_bond`. Hydrogens, solvent and monoatomic ions are dropped. Prediction and reference must contain the same molecules (sequence identity of at least 0.95 per chain); otherwise the id fails with peppr's message in `missing`.

**Tables**:
- `peppr`: one row per prediction that entered, NaN where a metric could not be computed.

  | id | rmsd | tm_score | lddt | lddt_pli | fnat | irmsd | lrmsd | dockq |
  |----|------|----------|------|----------|------|-------|-------|-------|

  The columns are `id` plus the requested `metrics`, in that order. A metric peppr cannot compute for a system, such as an interface metric on a monomer, is NaN and is not an error.
- `missing`: ids that could not be evaluated at all, with `kind` = `failure` and peppr's reason in `cause`. These are a prediction with no matching reference, a file that cannot be read, atoms that cannot be matched to the reference, or every requested metric failing. Their row in `peppr` is all NaN. Ids dropped by an upstream filter are carried over from the upstream `missing` table.

  | id | removed_by | kind | cause |
  |----|------------|------|-------|

**Example**:
```python
from biopipelines.peppr import Peppr

# Refolds against one crystal structure
scores = Peppr(structures=boltz, references=crystal)

# Each refold against the design it derives from, interface metrics only
scores = Peppr(structures=refolds, references=designs, metrics=["fnat", "irmsd", "lrmsd", "dockq"])

scores.tables.peppr          # id | fnat | irmsd | lrmsd | dockq
```

**Reference**: peppr, https://github.com/aivant/peppr. The pair metrics and DockQ follow Mirabello & Wallner (2024), https://doi.org/10.1093/bioinformatics/btae586.

---

### PLM_Sol

Sequence-based protein solubility prediction. PLM_Sol embeds each sequence with ProtT5-XL and classifies it with a biLSTM/TextCNN head trained on the updated E. coli solubility dataset (UESolDS); it is among the strongest sequence-only solubility predictors. It needs no structure — only sequence — so it complements the structure-based aggregation scorer Aggrescan3D. The output `solubility` is a probability in [0,1] (higher = more soluble); `soluble` is the ≥ 0.5 binary call.

**Tags**: predict-property, protein, solubility, sequence-only

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

### SequenceIdentityGate

Novelty gate: the highest global sequence identity of each designed sequence to a set of natural references, with the nearest reference. A design that copies a natural protein can score well on every structural metric because the predictor may have memorized it too; this table says how close each design sits to the nearest natural sequence. The tool never filters: `passed` states whether `max_identity` is strictly below `threshold`, and every design keeps its row, so the gate is applied downstream with `Panda.filter`.

Identity comes from a global (Needleman-Wunsch) alignment of the design against every reference: identical aligned positions divided by the alignment length, gaps included, so a design that shares a long stretch with a longer natural protein is not scored as if it were that protein.

**The alignment.** BLOSUM62 with gap-open -10 and gap-extend -0.5, which are the gap scores of EMBOSS Needle (its protein defaults of 10 and 0.5), but with end gaps penalized like any other gap. Needle does not penalize end gaps unless asked to (`-endweight`); Biopython's `PairwiseAligner`, which runs here, does by default, and the tool keeps that. For a design shorter than its reference this changes which positions are paired and so the identity. In a synthetic test of 60 random 300-residue proteins, each with an 80-residue design cut from it and half of its residues replaced, free end gaps (`end_gap_score=0`) and penalized ones gave a different identity for 25 designs, by up to about 0.02; with real sequences the difference depends on how conserved the stretch is. Free end gaps are not offered as a parameter. Confirm that the penalized choice suits your references. Among alignments tied on score Biopython's first one decides.

**The threshold.** The default 0.4 is the value of the `sequence_identity_gate` block of RosettaCommons Chaperone. That block describes it as "Maximum scaffold identity a design may have and still pass. Set per campaign; 0.40 is the RNA motif-scaffolding bench...", scores the scaffold region only (`exclude_motif=True`, because "motif identity is 100% by construction") and compares RNA designs with Rfam and RNAcentral homologs. Here the same number gates protein sequences over their whole length, so it was not measured for proteins, and a motif copied unchanged from a natural protein counts toward the identity and raises `max_identity` for every design that carries it. Set `threshold` for each campaign, and confirm the value you choose rather than keeping the default by inertia. Excluding the motif is not offered, because it needs the motif residues mapped from structure numbering onto sequence positions and no stream in the framework carries that mapping.

**Tags**: measure, protein

**Environment**: `biopipelines` (Biopython, already a base dependency; nothing to install).

**Parameters**:
- `sequences`: Union[DataStream, StandardizedOutput] (required) — protein sequences to score. One row is one polymer chain, so a multi-chain design is scored chain by chain.
- `references`: Union[DataStream, StandardizedOutput] (required) — natural protein sequences as a `sequences` stream; `Sequence("references.fasta")` turns a FASTA file into one. A reference that cannot be scored (empty, or a letter outside the BLOSUM62 alphabet such as `U`) stops the step instead of silently weakening the gate.
- `threshold`: float = 0.4 — maximum identity, as a fraction in [0, 1], a design may have and still pass. Set it per campaign: 0.4 is Chaperone's RNA scaffold-only value, not a protein threshold (see **The threshold**).

**Tables**:
- `identity_gate`: `id | max_identity | nearest_reference | passed`. `max_identity` is a fraction in [0, 1], `nearest_reference` the id of the reference that gave it (the first one on a tie), `passed` is `max_identity < threshold`. A design that could not be scored carries NaN in all three.
- `missing`: `id | removed_by | kind | cause`

Cost is one alignment per design and unique reference, so thousands of designs against hundreds of references take minutes to hours on one CPU.

**Example**:
```python
from biopipelines import ProteinMPNN, Sequence, SequenceIdentityGate, Panda

designs = ProteinMPNN(structures=rfd, num_sequences=8)
naturals = Sequence("naturals.fasta")
gate = SequenceIdentityGate(sequences=designs, references=naturals, threshold=0.4)

# keep the designs that are not close to any natural sequence
novel = Panda(tables=gate.tables.identity_gate,
              operations=[Panda.filter("passed == True")],
              pool=designs)
```

---

### PLIP

Protein interaction profiler. Reads complex PDBs and reports detected non-covalent interactions: hydrogen bonds, hydrophobic contacts, π-stacking, salt bridges, halogen bonds, water bridges, and metal complexes. Profiles protein–ligand, protein–peptide / protein–protein, and intra-chain interactions, selected by `mode`.

**Tags**: detect, complex, small-molecule, nucleic-acid, pocket, residues, interactions

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
- `interactions`: | id | ligand | interaction_type | residue | chain | resnum | distance | details |
- `summary`: | id | n_hbonds | n_hydrophobic | n_pi_stacking | n_salt_bridges | n_halogen | n_water_bridges | session_files |
- `missing`: | id | removed_by | kind | cause |

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

**Tags**: validate, small-molecule, complex, covalent

**Environment**: `posebusters`

**Parameters**:
- `structures`: Union[DataStream, StandardizedOutput] (required) - Protein-ligand complexes (PDB/CIF)
- `ligand`: Union[DataStream, StandardizedOutput] (required) - Compounds stream (`Ligand(codes="LIG")` or any compounds-producing tool) naming the ligand. The residue `code` is read from the stream at runtime, and is reused as the reference-structure residue code in redock mode. **Provide a SMILES for non-trivial ligands.** The ligand is extracted from each complex and rebuilt as an SDF for validation. If the stream carries a `smiles`, it is used as a **bond-order template** (`AssignBondOrdersFromTemplate`); otherwise bond orders are perceived from coordinates alone (`rdDetermineBonds`). Coordinate-only perception **fails on charged/conjugated ligands** (dyes, Si-rhodamines, metallo-organics, zwitterions): the molecule won't build, so `mol_pred_loaded` fails and **every** check — including the coordinate-based distance/overlap checks — reports False, i.e. the whole step yields `all_pass=False`. For such ligands you **must** pass `Ligand(smiles=..., codes="LIG")` whose `codes` matches the residue code in the structures. A bare `Ligand(codes="LIG")` (no SMILES) is fine only for simple ligands RDKit can perceive from coordinates, and additionally triggers a (failing) RCSB SMILES lookup unless the code is a real CCD code. Note: covalent-derived poses fail `minimum_distance_to_protein` by construction (the attachment atom sits at bonding distance); exclude that check via `exclude=` or `check=[...]` when validating them.
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

  A structure the tool could not validate, either because no ligand residue with the stream's code could be extracted from it or because PoseBusters raised on it, keeps its row with `all_pass` False and no check columns, and has a `failure` row in `missing` giving the reason. A pose that was validated and failed a check is a result, not a failure, and has no `missing` row.
- `missing`: | id | removed_by | kind | cause |

  The ids dropped upstream from `structures` or `ligand`, and this step's own failures (the structures above) as `failure` rows with `removed_by` the step id and the reason in `cause`.

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
    ligand=boltz_holo  # reads the LIG code Boltz2 assigned; or Ligand(codes="LIG")
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
    ligand=Ligand(codes="ATP"),
    reference_ligand=xrc,
    mode="redock"
)
```

---

### PoseChange

Measures ligand pose distance between reference holo structure and sample structures. Calculates RMSD and geometric metrics to quantify how well designed structures reproduce known binding poses.

**Tags**: measure, small-molecule, complex

**Environment**: `ProteinEnv`

**Parameters**:
- `reference_structure`: Union[str, ToolOutput, StandardizedOutput] (required) - Reference holo structure (e.g., XRC structure)
- `sample_structures`: Union[str, List[str], ToolOutput, StandardizedOutput] (required) - Designed/predicted structures to compare
- `reference_ligand`: Union[Ligand, ToolOutput] (required) - The ligand in the reference structure. Pass a `Ligand` / any compounds-producing tool output; the residue `code` is read from the compounds stream's map_table at runtime (Ligand Contract).
- `sample_ligand`: Union[Ligand, ToolOutput, None] = None - The ligand in the sample structures (default: same as reference_ligand). Pass a `Ligand` / compounds stream.
- `reference_alignment`: Optional[str] = None - PyMOL selection for reference structure alignment (default: "not resn {reference_ligand}", built at runtime once the code is resolved)
- `target_alignment`: Optional[str] = None - PyMOL selection for target structure alignment (default: "not resn {sample_ligand}", built at runtime once the code is resolved)

- `heavy_only`: bool = False - Compare heavy atoms only. Set this whenever one side has crossed a format that drops hydrogens — a PDBQT round trip through Vina or GNINA keeps only the polar ones — because the two ligands must carry the same atom count for the comparison to run at all, and a mismatch fails every row. A pose RMSD over hydrogens is noise anyway.

**Tables**:
- `changes`:

  | id | target_structure | reference_structure | ligand_rmsd | centroid_distance | orientation_angle | orientation_axis | alignment_rmsd | num_ligand_atoms | rmsd_pairing | reference_alignment | target_alignment |
  |----|------------------|---------------------|-------------|-------------------|-------------------|------------------|----------------|------------------|--------------|---------------------|------------------|

**Output Columns**:
- `ligand_rmsd`: RMSD between ligand poses after protein alignment (Å)
- `centroid_distance`: Distance between ligand centroids (Å)
- `alignment_rmsd`: RMSD of protein alignment (Å)
- `num_ligand_atoms`: Number of atoms in ligand
- `reference_alignment`: PyMOL selection used for reference structure alignment
- `target_alignment`: PyMOL selection used for target structure alignment

- `rmsd_pairing`: how the ligand atoms were put into correspondence — `name` (matched by atom name) or `graph` (matched by RDKit substructure, used when names are absent or duplicated)



**Ligand atoms are paired by identity, not by order.** `ligand_rmsd` matches the two ligands atom-name to atom-name, falling back to an RDKit graph match when names are unusable, and **raises** rather than pairing by position if neither works. Two files listing the same atoms in a different order therefore give ~0 Å instead of a large meaningless number. This changed the values every earlier run reported, so do not compare `ligand_rmsd` across the two behaviours. A reference and sample whose ligands differ in atom count are refused up front — usually a protonation difference, which is what `heavy_only=True` is for.

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

### LigandDdG

Protein-ligand binding energy with PyRosetta, with the ligand parameterized on the fly. ddG is the Rosetta score of the complex minus the score with the protein chain translated 10 000 Å away from the ligand; more negative means stronger predicted binding. Two scoring schemes, both ported from [Chaperone](https://github.com/RosettaCommons/chaperone) blocks, selected by `scoring`:

- `"soft"` (`ligand_ddg_soft`, the default): the interface within 10 Å of the ligand is repacked with beta_nov16, then the complex is scored with the soft-repulsive beta_nov16_soft against the separated state. This is repeated `replicates` times on the same pose; with at least three values the best and worst are dropped and the rest averaged. The soft score forgives the small clashes of predicted structures that would otherwise dominate the hard one.
- `"hard"` (`ligand_ddg`): backbone-constrained soft-then-hard minimization of the interface, then Rosetta's Ddg filter (beta_nov16, `repack=0`) across the ligand's jump. A single run. The protocol is the active part of Chaperone's `ddg_protocol.xml`, run through PyRosetta's RosettaScripts parser.

The numbers are Rosetta energy units, not kcal/mol and not measured affinities, and a hard value must not be compared with a soft one.

**Tags**: measure, complex, protein, small-molecule, binding, energy

**References**: https://github.com/RosettaCommons/chaperone, https://doi.org/10.1093/bioinformatics/btq007

**License**: PyRosetta is free for academic and non-profit use and needs a paid commercial license from the University of Washington (license@uw.edu). Installing it accepts the PyRosetta license; the Rosetta scripts fetched at install time fall under the Rosetta license.

**Environment**: `PyRosetta`, shared with [PPIDdG](#ppiddg) (see there for what the install builds); CPU only, no Daint entry.

**Ligand parameters.** Rosetta needs a `.params` file for a non-standard residue, and the PyRosetta wheel does not ship the tool that writes it. For each structure the ligand is cut out of the complex and rebuilt: bond orders and charges come from the SMILES of the `ligand` compounds stream, hydrogens are added at `ph` with Dimorphite-DL when the structure carries none (a ligand that already has hydrogens is used as written), atoms are renamed to unique element-first names, and Rosetta's `molfile_to_params.py` writes the params. Parameters are generated per structure, so atom names never have to agree across structures. If that cannot work for a ligand (no SMILES, or `molfile_to_params.py` cannot type it), the structure is recorded as a failure with the message, and `params_file` supplies pre-made parameters instead. Only the protein chain and the ligand are kept; other chains are dropped. Everything else on the protein chain stays as it is, hetero residues and waters included, except hydrogens, which are discarded for Rosetta to rebuild.

**Parameters**:
- `structures`: DataStream | StandardizedOutput (required) — PDB complexes; mmCIF is not read.
- `ligand`: str | DataStream | StandardizedOutput (required) — Compounds stream naming the ligand (`Ligand(smiles=..., codes="LIG")` or any compounds-producing tool); a bare string is shorthand for `Ligand(codes=...)`. The residue code is read from the stream's `code` column and must be the residue name of the ligand in the structures, at most 3 characters in a PDB. Unless `params_file` is given, the stream must carry a single SMILES; a code-only ligand cannot be parameterized.
- `scoring`: str = "soft" — `"soft"` or `"hard"`, as above.
- `binder_chain`: str = "A" — The protein chain, translated away to give the unbound state.
- `ligand_chain`: str = "B" — The chain carrying the ligand. Exactly one residue named like the ligand code must sit on it, so a structure with several copies of the ligand on that chain fails. It has to be a chain of its own, as predictors write it: the ligand of a crystal file usually shares the protein chain and has to be moved to its own first. A ligand with no atom within 6 Å of `binder_chain` (the ligand of another copy in a crystal, say) is a failure for that structure, since no pair of atoms is in range and the ddG would read as exactly zero.
- `replicates`: int = None — Pack + score cycles for `scoring="soft"`, default 5. Refused with `scoring="hard"`, which runs once.
- `ph`: float = 7.4 — pH for Dimorphite-DL protonation of a ligand the structure carries no hydrogens for. Ignored with `params_file` and for a ligand that already has hydrogens.
- `params_file`: str = None — Path to a pre-made Rosetta `.params` file. The ligand is then taken from the structure unchanged, so its atom names and hydrogens must already match the params; no SMILES is needed and nothing is generated.
- `seed`: int = 42 — Base random seed. Each replicate of each structure is seeded from it with a distinct, fixed value, so a run repeats exactly; PyRosetta otherwise seeds itself from the system. With `scoring="soft"` the seed does not change the value on the five streptavidin complexes tried (ten seeds each gave one ddG to four decimals, so `ddg_soft_std` is 0): the repack of the interface reaches the same optimum every time, and `replicates=1` gives the same number as 5.

**Tables**:
- `ligand_ddg`: the columns follow `scoring`.

  | id | ddg_norepack |
  |----|--------------|

  `scoring="hard"`: the Ddg filter value.

  | id | ddg_soft | ddg_soft_std | n_replicates_used |
  |----|----------|--------------|-------------------|

  `scoring="soft"`: the trimmed mean of the replicate ddG values, the population standard deviation of the retained ones, and how many were retained (`replicates`, minus two when there are at least three). One row per input structure; a structure that fails keeps its row with NaN values, and the reason is in `missing`. No threshold is applied; filter downstream with `Panda.filter`.

- `missing`: | id | removed_by | kind | cause |

**Example**:
```python
from biopipelines.ligand_ddg import LigandDdG

# Rank co-folded designs by how favorably Rosetta scores the bound ligand.
ligand = Ligand(smiles="CC(=O)Nc1ccc(O)cc1", codes="LIG")
folded = Boltz2(proteins=designs, ligands=ligand)
ddg = LigandDdG(structures=folded, ligand=ligand, scoring="soft", replicates=5)

best = Panda(tables=ddg.tables.ligand_ddg, operations=[Panda.sort("ddg_soft"), Panda.head(10)])
```

---

### PPIDdG

Protein-protein binding energy with PyRosetta. Each complex is relaxed with FastRelax (optionally constrained to the input coordinates), then ddG = score(relaxed complex) - score(`binder_chain` translated 10 000 Å away) is computed with beta_nov16 and no repacking of the separated state. Mirrors the `ppi_ddg` block of [Chaperone](https://github.com/RosettaCommons/chaperone).

The numbers are Rosetta energy units, not kcal/mol and not measured affinities; more negative means stronger predicted binding. They compare only across runs with identical chains, relaxation settings and score function, and the cheap setting (`relax_rounds=1`, `replicates=1`) is noisy: FastRelax is stochastic, so raise `replicates` before ranking close designs.

**Tags**: measure, complex, protein, binding, energy

**References**: https://github.com/RosettaCommons/chaperone, https://doi.org/10.1093/bioinformatics/btq007

**License**: PyRosetta is free for academic and non-profit use and needs a paid commercial license from the University of Washington (license@uw.edu). Installing it accepts the PyRosetta license.

**Environment**: `PyRosetta` (shared with [LigandDdG](#ligandddg)); CPU only, x86-64 Linux or Colab. PyRosetta has no aarch64 Linux wheel, so there is no Daint entry.

**Installation**: `PPIDdG.install()` (or `LigandDdG.install()`) creates the `PyRosetta` environment (Python 3.10, which `molfile_to_params.py` needs because it opens files in a mode Python 3.11 removed) from `environments/PyRosetta.yaml` and `PyRosetta.pip.txt`: `pyrosetta` from the RosettaCommons index, unpinned, plus RDKit and Dimorphite-DL. It then sparse-clones `source/scripts/python/public` of [RosettaCommons/rosetta](https://github.com/RosettaCommons/rosetta) into `<env prefix>/share/rosetta_public`, because the PyRosetta wheel does not ship `molfile_to_params.py`, which LigandDdG needs.

**Parameters**:
- `structures`: DataStream | StandardizedOutput (required) — PDB complexes; mmCIF is not read. Intended for the 20 standard amino acids: a ligand or other hetero residue makes Rosetta fail on that structure.
- `binder_chain`: str = "A" — The chain translated away to give the unbound state.
- `target_chain`: str = "B" — The partner that stays in place. It must exist in the structure. As in Chaperone's block, every chain other than `binder_chain` stays in place, so extra chains count as part of the target.
- `relax_rounds`: int = 1 — FastRelax repeats per replicate. One round takes roughly 1-3 minutes per structure on a CPU core; 3-5 gives lower, more converged energies.
- `replicates`: int = 1 — Independent relax + ddG cycles, each starting from the input. Each adds one full relax run; with 1 the standard deviation is 0.
- `constrain_to_start_coords`: bool = True — Harmonic coordinate constraints to the input during FastRelax, keeping the backbone near the supplied geometry. The reported scores use plain beta_nov16, so the constraint weight never enters them.
- `constraint_weight`: float = 1.0 — Weight of those constraints; higher stays closer to the input.
- `seed`: int = 42 — Base random seed. Each replicate of each structure is seeded from it with a distinct, fixed value, so a run repeats exactly (PyRosetta otherwise seeds itself from the system and a rerun gives other numbers). A different seed gives a different relaxation, with ddG differences of a few energy units at one replicate.

**Tables**:
- `ppi_ddg`:

  | id | ddg | ddg_std | n_replicates | relax_score |
  |----|-----|---------|--------------|-------------|

  `ddg` is the mean over replicates, `ddg_std` their population standard deviation, `relax_score` the total beta_nov16 score of the relaxed complex from the last replicate. One row per input structure: a structure that fails keeps its row with NaN values, and the reason is in `missing`. No threshold is applied; filter downstream with `Panda.filter`.

- `missing`: | id | removed_by | kind | cause |

**Example**:
```python
from biopipelines.ppi_ddg import PPIDdG

# Chain A of each predicted complex is the designed binder, chain B the target.
folded = Boltz2(proteins=designs)
ddg = PPIDdG(structures=folded, binder_chain="A", target_chain="B", relax_rounds=3, replicates=3)

# Most favorable first, with the replicate spread alongside for judging whether two scores differ.
best = Panda(tables=ddg.tables.ppi_ddg, operations=[Panda.sort("ddg"), Panda.head(10)])
```

---

### Prodigy

Binding-affinity prediction for **protein–protein** complexes. For each complex PDB, reports the predicted dissociation constant Kd and binding free energy ΔG between two chain groups, from interfacial contacts.

**Tags**: predict-property, complex, protein, binding

**References**: https://github.com/haddocking/prodigy

**Environment**: `prodigy`

**Parameters**:
- `structures`: DataStream | StandardizedOutput (required) — Protein–protein complex PDBs.
- `interface`: str = "A B" — The two chain groups forming the interface (e.g. `"A B"`, `"A,B C"`).
- `temperature`: float = 25.0 — Temperature (°C) for the Kd conversion.

**Tables**:
- `affinity`: | id | interface | delta_g_kcal_mol | kd_M | n_intermol_contacts | percent_charged_nis | percent_apolar_nis |
- `missing`: | id | removed_by | kind | cause |

**Example**:
```python
from biopipelines.prodigy import Prodigy

complex = PDB("1A22", convert="pdb")
aff = Prodigy(structures=complex, interface="A B")
```

---

### ProLIF

Protein–ligand interaction **fingerprints**. For each complex, computes the ProLIF interaction matrix and exports it in long form (one row per residue × interaction-type pair, with a binary `present` flag) so it joins cleanly with other tables.

**Tags**: detect, complex, small-molecule, residues, interactions

**References**: https://github.com/chemosim-lab/ProLIF

**Environment**: `prolif`

**Parameters**:
- `structures`: DataStream | StandardizedOutput (required) — Complex PDBs.
- `ligand`: DataStream | StandardizedOutput (required) — Compounds stream naming the bound ligand (the `code` is read at runtime; the bound HETATM block is extracted and bond-orders applied from the SMILES template).

**Tables**:
- `fingerprints`: | id | residue | resn | chain | resnum | resi | interaction_type | present |

  Long format: several rows per complex, none for a complex whose ligand makes no interaction ProLIF detects. Such a complex succeeded, so the step exits non-zero only when no complex was attempted or every one failed.
- `missing`: | id | removed_by | kind | cause |

**Example**:
```python
from biopipelines.prolif import ProLIF
from biopipelines.entities import Ligand

fp = ProLIF(structures=boltz_holo, ligand=Ligand(codes="LIG"))
```

---

### Reduce

Adds and optimises explicit hydrogens on protein and ligand atoms via the Richardson lab's `reduce` binary. Preserves protein residue topology and ligand HETATM records, so the output is suitable as input to interaction-fingerprint or MD-prep tools.

**Tags**: refine-structure, protein, small-molecule, all-atom, sidechain

**References**: https://github.com/rlabduke/reduce

**Environment**: `reduce`

**Parameters**:
- `structures`: DataStream | StandardizedOutput (required) — Input structures.

**Streams**: `structures` (protonated PDB per input)

**Example**:
```python
from biopipelines.reduce import Reduce

protonated = Reduce(structures=boltz)
prolif = ProLIF(structures=protonated, ligand=Ligand(codes="LIG"))
```

---

### RTMScore

Scores and ranks docking poses by residue–atom distance likelihood (graph transformers). Takes complex structures + a posed ligand; for each pose, auto-extracts the pocket within `cutoff` Å of the ligand and reports a score (higher = more likely a true binding pose).

> **The ligand must be posed inside the pocket** (same caveat as [GEMS](#gems)) — pass a ligand-producing tool output as `ligands=`, not a from-scratch SMILES.

**Tags**: predict-property, complex, small-molecule, pocket, binding

**References**: https://github.com/sc8668/RTMScore

**Resources**: GPU.

**Environment**: `rtmscore`

**Parameters**:
- `structures`: DataStream | StandardizedOutput (required) — Complex PDBs (protein source; pocket auto-extracted around the ligand).
- `ligands`: StandardizedOutput = None — A ligand-producing tool output (SDF references to score + SMILES metadata). Mutually exclusive with `ligands_3d`/`ligands_smiles`.
- `ligands_3d`: DataStream = None — Per-id ligand coordinate files passed directly (each file may hold multiple poses → one row per pose).
- `ligands_smiles`: DataStream = None — Optional `smiles` column for bond-order templating, paired with `ligands_3d`.
- `cutoff`: float = 10.0 — Pocket cutoff (Å) around the reference ligand.
- `model`: str = "model1" — Checkpoint (`model1`, `model2` or `model3`).

**Tables**:
- `scores`: | id | structures.id | ligands.id | pose | rtmscore |

  **Units.** `rtmscore` is dimensionless — an aggregated statistical potential summed from negative log-likelihoods. Upstream does not state which direction is better, so verify against the paper before ranking on it. It is a **pose-ranking score, not an affinity**: it is not calibrated to a binding constant and is not comparable to Prodigy's ΔG, GEMS's pKd or Boltz2's log-IC50.

  Long format: one row per pose of a pair, none for a pair RTMScore answered with no pose. Such a pair succeeded, so the step exits non-zero only when every pair failed.
- `missing`: | id | removed_by | kind | cause |



**Which pairs are scored.** By default every protein is scored against every ligand — the full cross product. When the protein and ligand id lists are **identical**, RTMScore instead scores only the diagonal, pairing each protein with the ligand of the same id. That is the shape produced by carving each ligand out of its own complex, where 375 of 376 cross-product cells would be scoring a pose against a protein it was never in.

The switch is inferred from the id lists, which has two consequences worth knowing. The lists are compared **in order**, so the same ids arriving in a different order fall back to the cross product. And a genuine N×N screen whose two streams happen to share ids collapses to N diagonal scores — nothing records that the off-diagonal pairs were not computed, so rename one side's ids if you want them.

**Example**:
```python
from biopipelines.rtmscore import RTMScore

scored = RTMScore(structures=complex, ligands=complex)
scored.tables.scores
```

---

### SASA

Solvent-accessible surface area analysis. Computes the change in SASA of a ligand when bound to vs separated from its protein partner — `delta_SASA = SASA(ligand alone) − SASA(ligand in complex)`. A larger delta indicates more ligand surface buried by the binder, typically interpreted as tighter packing.

**Tags**: measure, complex, small-molecule, residues, sasa

**Environment**: `ProteinEnv` (installed alongside `PyMOL`).

**Installation**: no extra step beyond `PyMOL.install()`; SASA reuses the same env.

**Parameters**:
- `structures`: Union[DataStream, StandardizedOutput] (required) — Protein-ligand complexes.
- `ligand`: Union[str, Ligand, ToolOutput] (required) — The bound ligand. Pass a `Ligand` / any compounds-producing tool output (e.g. a Boltz2 result); the residue `code` is read from the compounds stream's map_table at runtime (Ligand Contract). A bare code string (e.g. `"LIG"`, `"AMX"`, `":X:"`) is also accepted for back-compat.
- `dot_density`: int = 4 — PyMOL `get_area` dot density (1–4, higher = more accurate, slower).

**Streams** (`mode="residues"`):
- `accessibility`: per-residue relative accessibility (`resi-csv`), one file per input structure. In this mode SASA declares no tables.

**Tables** (`mode="ligand"`, the default):
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

**Single point mutations only.** ThermoMPNN predicts the ddG of *one* substitution at a time. Every value in the `ddg` table is an independent single-mutant ddG measured against the wildtype. It does **not** model the combined stability of several simultaneous mutations: summing per-mutation ddGs ignores epistasis and will mis-estimate a real multi-mutant. In the `mutations=` filter, each `+`-joined token is scored separately (it is a selection of single mutants, not a combined variant). For genuine double-mutant ddG (additive and epistatic) use the separate [ThermoMPNN-D](https://github.com/Kuhlman-Lab/ThermoMPNN-D) model; to judge the overall stability of a full multi-mutation design, fold it and compare a global metric to the wildtype.

**Tags**: predict-property, protein, residues, stability

**Environment**: `thermompnn`

**Installation**: `ThermoMPNN.install()` clones `Kuhlman-Lab/ThermoMPNN` (the default model checkpoint ships with the repo, so no weights download) and creates a torch env. Runs on CPU; a GPU speeds up large/many inputs.

**Parameters**:
- `structures`: Union[DataStream, StandardizedOutput] (required) — Wildtype backbone/complex PDBs. ThermoMPNN runs its site-saturation pass on these; they define `wildtype`/`resi` per position.
- `chain`: str = `"A"` — Chain to score.
- `mutations`: Union[str, (TableInfo, column)] = `""` — Optional. Empty (default) runs site-saturation over all positions of `chain`. Otherwise scores only the named point mutations, given as `+`-joined wildtype-position-mutant tokens (e.g. `"A42G+L50V"`, 1-indexed), or a table column reference whose per-structure cell holds such a string. In explicit mode the native saturation pass still runs and is then filtered to the requested mutants; any requested mutation absent from the output (e.g. a wildtype/position mismatch) is reported in `missing`. Each token is an independent single mutant.

**Tables**:
- `ddg`:

  | id | structures.id | chain | position | wildtype | mutation | ddG_pred |
  |----|---------------|-------|----------|----------|----------|----------|

  One row per scored mutation, sorted by `ddG_pred` ascending (most stabilising first).

- `missing`: `id | removed_by | kind | cause`

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

# Score a specific shortlist instead (each token is an independent single mutant)
ddg_subset = ThermoMPNN(structures=target, chain="A", mutations="A42G+L50V")
```

**Reference**: Dieckhaus et al. (2024) Transfer learning to leverage larger datasets for improved prediction of protein stability changes. *PNAS* 121, e2314853121. https://github.com/Kuhlman-Lab/ThermoMPNN

---

### VespaG

Zero-shot single-substitution fitness prediction from sequence. VespaG is a small expert-guided model distilled from ESM-2, leaderboard-competitive on ProteinGym while being MSA-free and fast. Higher score = predicted more tolerated/beneficial substitution. Complements ThermoMPNN (stability from structure) by scoring evolutionary/functional fitness from sequence alone.

**Tags**: predict-property, protein, residues, fitness, sequence-only

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

- `missing`: `id | removed_by | kind | cause`

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

**Tags**: measure, complex, small-molecule, energy

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
- `interaction_energies`: | id | e_complex_kj | e_protein_kj | e_ligand_kj | e_interaction_kj | e_interaction_kcal | charge_complex |
- `missing`: | id | removed_by | kind | cause |

**Example**:
```python
from biopipelines.xtb import XTB
from biopipelines.entities import Ligand

# Restrict to the pocket for speed
e_int = XTB(structures=boltz_holo, ligand=Ligand(codes="LIG"), method="gfn2")
```
