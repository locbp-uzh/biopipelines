# Structure Prediction & Docking

[← Back to Tool Reference](../tool_reference.md)

---

### AlphaFold

Predicts protein structures from amino acid sequences using AlphaFold2. Generates high-confidence 3D models with optional relaxation. Supports single-sequence prediction mode for fast predictions without MSA generation.

Single-sequence mode has no parameter of its own: forward ColabFold's flag, which a Python keyword cannot spell because it holds a dash, as `AlphaFold(proteins=seqs, **{"msa-mode": "single_sequence"})`. Written `msa_mode=` it reaches `colabfold_batch` with an underscore, which that parser rejects (`unrecognized arguments: --msa_mode`).

By default each input sequence is folded as a separate monomer. Wrap `proteins` in `Bundle(...)` to fold the bundled sequences together as one multi-chain complex: ColabFold receives a single colon-joined query (`SEQ_A:SEQ_B`), auto-selects the multimer model, and runs its default paired+unpaired MSA pipeline. `Bundle(static, Each(...))` holds `static` fixed and folds it against each iterated sequence (one complex per element). A bundled complex's output id is the chain ids joined with `+` (e.g. `p1+p2`). Pre-computed `msas=` cannot be combined with a `Bundle` (a bundled complex always uses ColabFold's own MSA pipeline).

**Tags**: predict-structure, protein, complex

**Resources**: GPU. Verified on A100 and on H100 (sm_90); the CUDA 12 jax stack ships sm_90 kernels, so Hopper cards no longer need avoiding.

**Environment**: `biopipelines`, plus a conda **prefix** env at `<AlphaFold>/colabfold-conda` that `AlphaFold.install()` builds. It is addressed by path rather than by an `environments:` name because `MMseqs2Server` resolves `colabfold_search` and `mmseqs` inside that same `bin/`, and `AF2BIND` and `BioEmu` read the AF2 weights from `<AlphaFoldParams>` beside it.

**Installation**: `AlphaFold.install()`. It follows the [official ColabFold instructions](https://github.com/sokrypton/ColabFold#installation) — conda supplies python and the pinned `mmseqs2` build `colabfold_search` is tested against, pip supplies ColabFold itself and the CUDA 12 wheels:

```bash
mamba create -y -p <AlphaFold>/colabfold-conda -c conda-forge -c bioconda python=3.13 pip mmseqs2=18.8cc5c
<AlphaFold>/colabfold-conda/bin/pip install "colabfold[alphafold,openmm]==1.6.3" "jax[cuda12]" "openmm[cuda12]"
```

The versions are pinned in the wrapper as `AlphaFold.COLABFOLD_VERSION` and `AlphaFold.MMSEQS2_VERSION`; bump them there, not here.

Installation is verified by running `colabfold_batch --help` and importing `jax_cuda12_plugin`. The second check is what distinguishes a GPU-capable install from a CPU-only one, and it passes without a GPU present, so it is meaningful on a login node.

`AlphaFold.install(force_reinstall=True)` removes and rebuilds `colabfold-conda` but keeps `<AlphaFoldParams>`, so the weights are not re-downloaded.

**Parameters**:
- `proteins`: DataStream | StandardizedOutput | Bundle | Each (required) - Input protein sequences. A bare input folds one monomer per sequence; `Bundle(...)` folds the bundled sequences as one multi-chain complex.
- `msas`: StandardizedOutput = None - Pre-computed MSAs in A3M format. When provided,
  ColabFold uses these instead of generating new MSAs, and every query needs one, matched to its id:
  a query with none is not predicted, keeps a blank row of `confidence` and lands in `missing` as a `failure`.
  **Note**: Prefer an A3M the search produced — `MMseqs2(output_format="a3m")`,
  AlphaFold or ColabFold. `MSA(source, convert="a3m")` rebuilds a usable A3M only when
  the CSV has its `.a3m_meta` sidecar; without one (any CSV written before the sidecar
  existed) the `>` headers are synthesised indices, which fold as a monomer but cannot
  pair complex chains.
- `num_relax`: int = 0 - Number of best models to relax with AMBER
- `num_recycle`: int = 3 - Number of recycling iterations
- `rand_seed`: int = 0 - Random seed (0 = random)

**Streams**: `structures`, `pae`, `msas`

`pae` holds one `<id>.json` per structure id (the `structures` ids, so `p1+p2` for a bundled complex): ColabFold's scores JSON of the best model, copied as ColabFold wrote it. Its `pae` key is the full predicted aligned error matrix (a list of lists, residues x residues, in angstroms), next to `plddt` (per residue) and the scalars `max_pae` and `ptm` (plus `iptm` for a complex). A model without a PAE head (plain `alphafold2`, not the ptm or multimer models) writes no `pae` key; the file is still handed over, so a consumer reads the missing key as an error of its own. The stream's map lists only the ids whose JSON exists.

**Tables**:
- `structures`:

  | id | file |
  |----|------|

- `confidence`:

  | id | structure | plddt | max_pae | ptm |
  |----|-----------|-------|---------|-----|

- `msas`:

  | id | sequence | file |
  |----|----------|----------|

- `missing` (when `msas` is given, or an input axis carries an upstream manifest):

  | id | removed_by | kind | cause |
  |----|------------|------|-------|

**Example**:
```python
from biopipelines.alphafold import AlphaFold
from biopipelines.msa import MSA

# Standard prediction
af = AlphaFold(
    proteins=lmpnn,
    num_relax=1,
    num_recycle=5
)

# With pre-computed MSAs from a previous AlphaFold run
af = AlphaFold(proteins=seq, msas=previous_af_result)

# Fold a single A:B complex (one prediction, id "chainA+chainB")
from biopipelines.combinatorics import Bundle, Each
af = AlphaFold(proteins=Bundle(seqs_a, seqs_b))

# One complex per binder: fixed receptor folded against each binder
af = AlphaFold(proteins=Bundle(receptor, Each(binders)))
```

---

### Boltz2

Predicts biomolecular complexes including proteins, nucleic acids, and small molecules. State-of-the-art model for protein-ligand and protein-protein complex prediction.

**Tags**: predict-structure, predict-property, protein, small-molecule, nucleic-acid, complex, binding, covalent, all-atom, sequence-only

**Installation**:
```bash
mamba create -n Boltz2Env python=3.11
mamba activate Boltz2Env
pip install boltz[cuda] -U
```
**Environment**: `Boltz2Env`

**Parameters**:
- `config`: Optional[str] = None - Direct YAML configuration string
- `proteins`: Optional[Union[DataStream, StandardizedOutput]] = None - Protein sequences
- `ssDNA`: Optional[Union[DataStream, StandardizedOutput]] = None - Single-stranded DNA sequences
- `dsDNA`: Optional[Union[DataStream, StandardizedOutput]] = None - Double-stranded DNA sequences
- `ssRNA`: Optional[Union[DataStream, StandardizedOutput]] = None - Single-stranded RNA sequences
- `dsRNA`: Optional[Union[DataStream, StandardizedOutput]] = None - Double-stranded RNA sequences
- `ligands`: Optional[Union[DataStream, StandardizedOutput]] = None - Ligand compounds stream (e.g. from `Ligand` / `CompoundLibrary`). May be the **only** input: a ligands-only job folds or scores compounds on their own. Pass `affinity=False` when you do — `affinity` defaults to `True` and a ligands-only prediction has no receptor to bind to, so `Boltz2(ligands=lib)` raises rather than silently scoring nothing.
- `msas`: Optional[StandardizedOutput] = None - Pre-computed MSA files for recycling (pass entire tool output, not .msas). Every protein chain needs an alignment matched to its protein id (or, failing that, its sequence): a complex with a chain that has none is not predicted, keeps a blank row of `confidence` and lands in `missing` as a `failure`; Boltz2 no longer folds that chain single-sequence
  Supports recycling from AlphaFold via `MSA(af_result, convert="csv")`.
- `single_sequence`: bool = False - Fold without an MSA. Mutually exclusive with `msas`. The `msas` stream comes back empty, and accuracy drops for anything that is not already well represented in the model — use it for speed or when no homologs exist, not as a default.
- `affinity`: bool = True - Calculate binding affinity predictions
- `output_format`: str = "pdb" - Output format (pdb, mmcif)
- `recycling_steps`: Optional[int] = None - Number of recycling steps (default: model-specific)
- `diffusion_samples`: Optional[int] = None - Number of diffusion samples (default: model-specific)
- `sampling_steps`: Optional[int] = None - Diffusion sampling steps (boltz default 200)
- `step_scale`: Optional[float] = None - Diffusion temperature; lower samples more diverse poses (boltz-2 default 1.5)
- `max_msa_seqs`: Optional[int] = None - Cap on MSA depth (boltz default 8192). The first thing to lower when a large complex runs out of VRAM, since peak memory scales with MSA depth
- `subsample_msa`: bool = False - Sample the MSA down rather than truncating it, keeping its diversity at a lower depth
- `num_subsampled_msa`: Optional[int] = None - How many sequences to subsample (boltz default 1024). Requires `subsample_msa=True`; boltz ignores the count without the flag, so the wrapper refuses the pair rather than letting the run proceed at full depth
- `max_parallel_samples`: Optional[int] = None - Diffusion samples held on the GPU at once (boltz default 5). Lower it to trade wall time for peak memory when `diffusion_samples` is high
- `no_kernels`: bool = False - Disable the trifast/cuequivariance triangular kernels. Needed on GPUs those kernels do not support; slower where they do
- `top_only`: bool = True - When True (default), keep only the top model as `<id>`. When False, surface every diffusion sample as a separate structure `<id>_1..N` — e.g. to feed a Boltz2 pose ensemble (with native covalent linkage / full scaffold) into downstream design. Pair with `diffusion_samples=N` to control how many.
- `use_potentials`: bool = False - Enable external potentials
- `template`: Optional[str] = None - Path to PDB template file for structure guidance
- `template_chain_ids`: Optional[List[str]] = None - Chain IDs to apply template to (e.g., ["A", "B"])
- `template_force`: bool = True - Force template usage
- `template_threshold`: float = 5.0 - RMSD threshold for template matching
- `pocket_residues`: Optional[List[int]] = None - Residue positions defining binding pocket (e.g., [50, 51, 52])
- `pocket_max_distance`: float = 7.0 - Maximum distance for pocket constraint
- `pocket_force`: bool = True - Force pocket constraint
- `glycosylation`: Optional[Dict[str, List[int]]] = None - N-glycosylation sites per chain (e.g., {"A": [164]})
- `covalent_linkage`: Optional[Dict[str, Any]] = None - Covalent attachment specification (see examples). The `ligand_atom` is a Boltz atom name — get it from `Boltz2.predict_atom_names(...)` (see below).
- `contacts`: Optional[List[Dict[str, Any]]] = None - Contact constraints between atoms/residues. A ligand token is `[chain_id, atom_name]` (e.g. `["B", "C8"]`); get the atom name from `Boltz2.predict_atom_names(...)` (see below).
- `disulfide_bonds`: Optional[List[Dict[str, Any]]] = None - Disulfide bond constraints
- `metal_coord`: Optional[List[Dict[str, Any]]] = None - Metal-coordination bond constraints

**Streams**: `structures`, `pae`, `sequences`, `compounds` (only when ligands are folded), `msas` (empty when `single_sequence=True`)

`pae` holds one `<id>.npz` per structure id: the full predicted aligned error matrix of that model, as a numpy array under the key `pae` with shape tokens x tokens (`numpy.load(path)["pae"]`). Its ids are exactly the `structures` ids, so with `top_only=False` each diffusion sample `<id>_<k>` has its own matrix, taken from the same model as its structure. The wrapper passes `--write_full_pae` to `boltz predict` and copies `pae_<input>_model_<n>.npz` from the raw `execution/` dump; the stream's map lists only the ids whose file Boltz wrote. The file is a compressed npz, but the matrix grows with the square of the token count, so a large library adds up.

**Tables**:
- `structures`:

  | id | file |
  |----|------|

- `confidence`:

  | id | input_file | confidence_score | ptm | iptm | complex_plddt | complex_iplddt |
  |----|------------|------------------|-----|------|---------------|----------------|

- `affinity` (only when `affinity=True`):

  | id | input_file | affinity_pred_value | affinity_probability_binary |
  |----|------------|---------------------|-----------------------------|

  **Units.** `affinity_pred_value` is `log10(IC50)` with IC50 in µM, so **lower is stronger** — IC50 10⁻⁹ M gives −3, 10⁻⁶ M gives 0, 10⁻⁴ M gives 2. `affinity_probability_binary` is the probability the ligand binds at all, 0–1, **higher is more likely**. The two come from different training data and answer different questions: the probability for hit-finding, the value for hit-to-lead optimization. Neither is comparable to another tool's affinity number.

  Prefer `affinity_probability_binary` (binder probability, higher = more likely a binder) for ranking — it is more reliable than the `affinity_pred_value` regression score in most cases.

- `sequences`:

  | id | sequence |
  |----|----------|

- `msas` (only when `single_sequence=False`):

  | id | sequence | file |
  |----|----------|----------|

- `compounds` (only when ligands are folded; `code` carries the residue code Boltz assigned):

  | id | format | code | smiles | ccd |
  |----|--------|------|--------|-----|

- `missing` (when `msas` is given, or an input axis carries an upstream manifest):

  | id | removed_by | kind | cause |
  |----|------------|------|-------|

  Beyond the declared columns above, `structures`, `confidence` and `affinity` each gain one `<axis>.id` provenance column per input axis, named after the parameter (`proteins.id`, `ligands.id`, …): the row `ubq+lig1+lig_2` has `proteins.id` = `ubq+lig1` and `ligands.id` = `lig_2`, appended at runtime by the postprocessing step (`pipe_boltz_postprocessing.py`); the `pae` map carries the same ones. They track which protein and ligand produced each row, so filtering and joins need no ID-string parsing. The same step also writes an undeclared `value` column onto `structures` and `pae`; treat it as an implementation detail rather than part of the contract.

**Example**:
```python
from biopipelines.boltz2 import Boltz2

# Basic apo and holo prediction
boltz_apo = Boltz2(proteins=lmpnn)
boltz_holo = Boltz2(
    proteins=lmpnn,
    ligands=Ligand(smiles="CC(=O)OC1=CC=CC=C1C(=O)O"),  # Aspirin SMILES
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

```

**Predicting ligand atom names** (`Boltz2.predict_atom_names`)

Constraints that target a ligand atom — `contacts`, `covalent_linkage`, `metal_coord` — need the atom name Boltz assigns to the ligand (e.g. `["B", "C8"]`). Boltz names ligand atoms deterministically from chemistry, so you can look them up *before* running anything. `Boltz2.predict_atom_names(ligand)` returns a 2-D depiction with every heavy atom labelled by its Boltz name (renders inline in Jupyter and is written to a PNG), plus a `.names` lookup. It is a configuration-time helper: no GPU, no Pipeline context required, and it works in a plain Python session.

```python
from biopipelines.boltz2 import Boltz2
from biopipelines import Ligand

res = Boltz2.predict_atom_names(Ligand("ampicillin"))
res          # labelled image renders inline in a notebook
res.path     # PNG written to ./ampicillin_atoms.png (override with path=)
res.names    # {"ampicillin": {"S24": ..., "C25": ..., ...}}
```

It accepts a SMILES string, a `Ligand` / `CompoundLibrary` (one labelled panel per compound), or the output one returns inside a pipeline. Two naming regimes, matching Boltz exactly:

- **SMILES ligands** (raw SMILES, `Ligand(smiles=...)`, a PubChem name/CID/CAS lookup, `CompoundLibrary`): `element.upper()` + canonical rank index, with hydrogens included in the ranking — so heavy-atom indices are not `1..N` (ethanol → `C8`, `C9`, `O7`).
- **CCD ligands** (`Ligand("ATP")` and any RCSB CCD code): the CCD's own atom names from the RCSB definition (`PG`, `O5'`, `C1'`, …).

A network call is made only to resolve a name/CID/CAS or CCD lookup; a raw SMILES or `Ligand(smiles=...)` is fully offline.

---

### DiffDock

Blind diffusion-based docking: samples candidate ligand poses over translation, rotation, and torsion, then re-ranks them with a confidence model. Works directly from a protein PDB + a ligand (SMILES or SDF) **without a binding-box hint** — useful when the pocket is unknown. Output IDs follow the multi-axis pattern `<protein>+<ligand>_rank<N>`; the rank-1 pose is the `structures` stream, and the full ranked list is in the table.

**Tags**: dock, protein, small-molecule, complex

**References**: https://github.com/gcorso/DiffDock · https://arxiv.org/abs/2210.01776

**Resources**: GPU (CPU fallback exists but is much slower).

**Environment**: `diffdock`

**Installation**: `DiffDock.install()` clones the repo + the facebookresearch/esm side-clone, builds the env, and installs the matching PyG wheels.

**Parameters**:
- `structures`: DataStream | StandardizedOutput (required) — Protein backbones (PDB).
- `compounds`: DataStream | StandardizedOutput (required) — Ligands as SMILES or SDF.
- `samples_per_complex`: int = 10 — Poses sampled per (protein, ligand) pair.
- `inference_steps`: int = 20 — Reverse-diffusion steps.
- `actual_steps`: int = None — Steps actually executed (≤ `inference_steps`; default `inference_steps - 1`).
- `no_final_step_noise`: bool = True — Disable noise on the last step.
- `batch_size`: int = 32 — DiffDock internal batch size.

**Streams**: `structures` (rank-N pose SDF per pair)

**Tables**:
- `confidence`: | id | structures.id | compounds.id | rank | confidence |
- `missing`: | id | removed_by | kind | cause |

**Example**:
```python
from biopipelines.diffdock import DiffDock

target = PDB("4ufc", convert="pdb")
lig = Ligand(smiles="CC(=O)Oc1ccccc1C(=O)O", ids="ASA")
dock = DiffDock(structures=target, compounds=lig)
```

---

### DynamicBind

Flexible-backbone docking. DynamicBind is an equivariant generative model that predicts a **ligand-specific protein conformation**, letting the receptor backbone flex toward its bound state rather than docking into a rigid pocket. Reports per-pose lDDT and a predicted affinity; can optionally render a transition movie.

**Tags**: dock, sample-ensemble, protein, small-molecule, complex, binding, flexibility, all-atom, sidechain

**References**: https://github.com/luwei0917/DynamicBind · https://www.nature.com/articles/s41467-024-45461-2

**Resources**: GPU. Pin `gpu="A100"` — DynamicBind's torch build is cu11x and crashes ("no kernel image") on newer cards like the H100.

**Environment**: `dynamicbind` (inference) + `dynamicbind_relax` (OpenMM pose relaxation); both created at install.

**Installation**: `DynamicBind.install()` clones the repo, creates both envs, and fetches the weights bundle from Zenodo.

**Parameters**:
- `structures`: DataStream | StandardizedOutput (required) — Protein structures (PDB).
- `compounds`: DataStream | StandardizedOutput (required) — Ligands (SMILES or SDF).
- `num_samples`: int = 40 — Poses sampled per pair.
- `num_saved`: int = None — How many top poses to keep (default: all).
- `inference_steps`: int = 20 — Reverse-diffusion steps.
- `num_workers`: int = 1 — Data-loading workers.
- `rigid_protein`: bool = False — Keep the protein rigid (disables the flexing that distinguishes DynamicBind).
- `make_movie`: bool = False — Render a conformational-transition movie per pose (`movies` stream).
- `seed`: int = 42 — Random seed.

**Streams**: `structures` (predicted complex per pose), `movies` (only when `make_movie=True`)

**Tables**:
- `affinity`: | id | structures.id | compounds.id | rank | lddt | affinity |

  **Units.** `affinity` is passed through from upstream, which does not document its scale or direction — the DynamicBind README states only that it writes an "Affinity Score". Treat it as an ordinal score within one run and do not compare it to another tool's number, or convert it to a Kd, without checking the paper first.
- `missing`: | id | removed_by | kind | cause |

**Example**:
```python
from biopipelines.dynamicbind import DynamicBind

with Pipeline("Examples", "DynamicBind-demo"):
    Resources(gpu="A100", memory="32GB", time="6:00:00")
    target = PDB("4ufc", convert="pdb")
    lig = Ligand(smiles="CC(=O)Oc1ccccc1C(=O)O", ids="ASA")
    dock = DynamicBind(structures=target, compounds=lig)
```

---

### ESMFold

Predicts a protein's 3-D structure directly from its amino-acid sequence using a protein language model — **no MSA or templates needed**, so it is much faster to set up than AlphaFold for single-sequence predictions. One PDB per input sequence, with per-residue pLDDT (in the B-factor column) and a global pTM.

**Tags**: predict-structure, protein, sequence-only, all-atom

**References**: https://github.com/facebookresearch/esm

**Resources**: GPU. Long sequences may need `chunk_size` / `max_tokens_per_batch` lowered to fit memory; set `memory="32GB"`.

**Environment**: `esmfold`

**Installation**: `ESMFold.install()` creates the env (mirrors the One-command-install-ESMfold recipe on the cluster; the ColabFold notebook recipe on Colab).

**Parameters**:
- `sequences`: DataStream | StandardizedOutput (required) — Input sequences.
- `chunk_size`: int = None — Axial-attention chunk size to reduce GPU memory (try 128/64/32 for long sequences that OOM). `None` = no chunking.
- `num_recycles`: int = 4 — Recycling iterations.
- `max_tokens_per_batch`: int = 1024 — Max tokens per forward pass; lower to avoid OOM.
- `cpu_offload`: bool = False — Offload weights to CPU to save GPU memory.
- `cpu_only`: bool = False — Run entirely on CPU (very slow).

**Streams**: `structures`

**Tables**:
- `structures`: | id | file | sequences.id |
- `confidence`: | id | file | plddt | ptm |

**Example**:
```python
from biopipelines.esmfold import ESMFold

seqs = Sequence(["MKTVRQERLKSIVRILERSKEPVSGAQ"], ids=["p1"])
esm = ESMFold(sequences=seqs)
```

---

### ESMFold2

Predicts all-atom **biomolecular complexes** (proteins, DNA, RNA, ligands; antibody scFvs as multi-chain proteins) from sequence, building on ESMC representations through a recurrent folding stack and a diffusion module. Folds MSA-free in the single-sequence regime, or recycles a precomputed MSA per protein chain. Reports per-atom pLDDT, pTM, interface pTM (ipTM), and optionally PAE. Inference-time scaling is a first-class control: more folding loops, diffusion samples, or seeds raise quality, with best-of-N selection by ipTM (complex) / pLDDT (monomer).

**Tags**: predict-structure, protein, small-molecule, nucleic-acid, complex, all-atom, covalent, sequence-only

**References**: https://github.com/evolutionaryscale/esm — "Language Modeling Materializes a World Model of Protein Biology" (bioRxiv 2026). MIT-licensed; weights `biohub/ESMFold2` on HuggingFace.

**Resources**: GPU, large memory (built on ESMC-scale representations); set e.g. `memory="64GB"` and an explicit `gpu="A100"` / `gpu="H100"` rather than `gpu="any"`.

**Environment**: `esmfold2`

**Installation**: `ESMFold2.install()` creates the env and pip-installs the `esm` package from its `main` branch (esm >= 3.4 is required: it carries the model as `esm.models.esmfold2.EsmFold2Model`). Verifies by importing `EsmFold2Model`. An environment built from esm 3.3.0 cannot load the current `biohub/ESMFold2` checkpoint; rebuild it with `ESMFold2.install(force_reinstall=True)`. The Daint environment pins esm to a 3.3.0 commit (`environments/esmfold2.pip.daint.txt`), which this version cannot use.

**Parameters**:
- `proteins`: DataStream | StandardizedOutput = None — Protein sequences.
- `ssDNA` / `dsDNA` / `ssRNA` / `dsRNA`: DataStream | StandardizedOutput = None — Nucleic-acid sequences (double-stranded axes add the reverse-complement chain automatically).
- `ligands`: DataStream | StandardizedOutput = None — Compounds stream (`Ligand` / `CompoundLibrary`); folded by CCD code (SMILES ligands supported too).
- `msas`: StandardizedOutput = None — Precomputed MSAs to recycle per protein chain (e.g. from MMseqs2 or a previous run). Omit for single-sequence folding. With it, every protein chain needs an alignment matched to its protein id (or, failing that, its sequence): a complex with a chain that has none is not folded, keeps a blank row of `confidence` and lands in `missing` as a `failure`.
- `modifications`: Optional[List[Dict]] = None — Per-residue CCD substitutions on a polymer chain, e.g. `[{"chain": "A", "position": 122, "ccd": "SEP"}]`. The residue is replaced by its CCD entry and atom-tokenized, so `SEP` / `TPO` / `PTR` give real phospho-residues with their phosphate rather than the unmodified parent. Works on DNA/RNA chains too (`PSU`, …). CCD codes only — the upstream `Modification.smiles` field is not implemented, and a `smiles` key is rejected rather than ignored.
- `covalent_bonds`: Optional[List[Dict]] = None — Inter-entity covalent links, e.g. `[{"atom1": ["A", 122, "SG"], "atom2": ["C", 1, "C1"]}]`. Atom names are resolved against the residue's CCD atoms at runtime; an unknown chain, residue, or atom name raises instead of being silently dropped.
- `num_loops`: int = 10 — Recurrent folding loops (inference-time scaling).
- `num_sampling_steps`: int = 100 — Diffusion sampling steps.
- `num_diffusion_samples`: int = 1 — Diffusion samples per seed.
- `num_seeds`: int = 1 — Independent seeds; the best by ipTM (complex) / pLDDT (monomer) is kept.
- `msa_max_depth`: int = 1024 — MSA subsampling depth per loop.
- `include_pae`: bool = False — Also tabulate PAE (adds `max_pae` to the confidence table) and write each structure's full PAE matrix to the `pae` stream.
- `top_only`: bool = True — Keep only the best sample per complex. When False, surface every diffusion sample as `<id>_1..N` (pair with `num_diffusion_samples=N`).
- `model_name`: str = "biohub/ESMFold2" — HuggingFace checkpoint id (e.g. a Fast variant).

**Streams**: `structures` (mmCIF), `pae` (only with `include_pae=True`), `compounds` (ligand passthrough when ligands are folded)

`pae` holds one `<id>.npz` per structure, under the same ids as `structures`. The array under key `pae` is the N x N predicted aligned error in Å (0 to 32), float32, over all chains concatenated in token order, where a token is one residue, or one atom for ligands and modified residues. This is the layout Boltz writes, so one reader serves both. The stream's map_table lists only the ids whose matrix was written; a matrix that is not square is skipped with a warning and the structure is kept. Read one with `np.load(path)["pae"]`.

**Tables**:
- `structures`: | id | file |
- `confidence`: | id | file | plddt | ptm | iptm | max_pae | — `max_pae` is present only with `include_pae=True`.
- `compounds` (only when ligands are folded): | id | format | code | smiles | ccd |
- `missing`: | id | removed_by | kind | cause | — always declared. A complex that raises while folding (a modification past the end of a chain, an unknown atom name, an out-of-memory error) keeps its id as a row of `confidence` with every score blank and gets a `failure` row here carrying the first line of the error; the other ids are still folded. Ids dropped upstream are carried over.

`structures` gains one `<axis>.id` provenance column per input axis at runtime, appended after the declared columns (`pipe_esmfold2_postprocessing.py:79`). `confidence` does not — it is written without them.

**Example**:
```python
from biopipelines.esmfold2 import ESMFold2
from biopipelines.combinatorics import Each

# Fold each designed binder against a fixed receptor, 4 seeds, keep best by ipTM
rec = Sequence(["MIEIKDKQLTGLRFIDLFAGLGGFRLALES..."], ids=["EGFR"])
cplx = ESMFold2(proteins=Each(rec, binders), num_seeds=4)

# Protein + ligand complex with recycled MSAs and PAE
holo = ESMFold2(proteins=rec, ligands=Ligand("ATP"), msas=prev_run, include_pae=True)

# Phospho-substrate folded against a phosphatase
pair = ESMFold2(
    proteins=Each(substrate, phosphatase),
    modifications=[{"chain": "A", "position": 122, "ccd": "SEP"}],
    include_pae=True,
)
```

Chain ids in `modifications` and `covalent_bonds` are assigned in axis order — proteins, then ssDNA/dsDNA/ssRNA/dsRNA, then ligands — starting at `A`, with a double-stranded axis taking two consecutive letters. They name a position in the complex rather than a sequence id, so when an axis iterates, the same constraint applies to every item of that axis.

There is no `pocket` parameter. Upstream defines `PocketConditioning` and stores it on `StructurePredictionInput`, but no released checkpoint consumes it: esm 3.3.0 zeroed `pocket_feature` and discarded the constraint, and esm 3.4 raises `NotImplementedError`. Verify that path before adding one.

---

### Gnina

Molecular docking and pose scoring with a CNN. Combines AutoDock Vina search with a convolutional neural network for more accurate binding pose prediction. Supports multi-run docking with statistical analysis across independent runs, optional conformer generation, and pose consistency analysis.

`mode` selects the search regime:
- `"docking"` (default) — full Vina search + CNN rescore. The ligand 3-D pose is built from `compounds` chemistry; the binding box is determined in order: explicit `center`+`size` > `autobox_ligand` > crystal ligand HETATM records in the input PDB.
- `"score"` — no search. Score the ligand exactly where it already sits in the input complex (`gnina --score_only`). Use this to get GNINA's Vina/CNN scores for an already-posed complex (a Boltz2 co-fold, a crystal pose, or a Gnina docking output) without re-docking.
- `"minimize"` — local energy minimization of the in-pocket pose, then score (`gnina --minimize`). Same inputs as `"score"`; emits the refined pose.

In `"score"`/`"minimize"` mode the ligand is extracted from each complex's HETATM records (identified by the `compounds` residue code), bond-order-templated against its SMILES, and autoboxed on its own location (`autobox_add` padding). The `compounds` stream **must carry SMILES** in these modes — a code-only `Ligand`/`ligand="LIG"` string is rejected (no chemistry to template bond orders from). The docking-only parameters (`center`, `size`, `autobox_ligand`, `exhaustiveness`, `num_modes`, `num_runs`, `generate_conformers` and the conformer knobs, `rmsd_threshold`, `top_only`, `rank_conformers_by`) must stay at their defaults outside `"docking"` mode — setting one raises.

**Tags**: dock, predict-property, protein, small-molecule, pocket, binding, energy

**Resources**: GPU recommended (CPU fallback available but slow).

**Environment**: `biopipelines` (plus CUDA modules configured under `gnina:` in your `config.<variant>.yaml`)

**Installation**: Downloads pre-built `gnina.1.3.2` binary from GitHub.

**Parameters**:
- `structures`: Union[DataStream, StandardizedOutput] (required) - Protein structures. In `score`/`minimize` mode, complexes that already carry the ligand as HETATM.
- `compounds`: Union[DataStream, StandardizedOutput] (required) - Ligands (from `Ligand()`, `CompoundLibrary()`, etc.). In `score`/`minimize` mode must carry SMILES.
- `mode`: str = "docking" - One of `"docking"`, `"score"`, `"minimize"`.
- `scoring`: str = "vina" - Vina force field (`"vina"`, `"vinardo"`, `"ad4"`); only valid on [Vina](#vina), which selects the engine without a CNN.
- `autobox_ligand`: Union[DataStream, StandardizedOutput, str, None] = None - Reference ligand for automatic box (docking only)
- `center`: Optional[str] = None - Explicit box center as `"x,y,z"`
- `size`: Union[float, str, None] = None - Box dimensions in Angstroms (single float = cubic, `"x,y,z"` = asymmetric)
- `autobox_add`: float = 4.0 - Padding around autobox ligand (Angstroms)
- `exhaustiveness`: int = 8 - Search exhaustiveness
- `num_modes`: int = 9 - Docked poses per run
- `num_runs`: int = 1 - Independent docking runs per conformer
- `seed`: int = 42 - Base random seed (each run uses `seed + run_index`)
- `cnn_scoring`: str = "rescore" - CNN scoring mode (`"rescore"`, `"refinement"`, `"none"`, `"rescore_only"`)
- `generate_conformers`: bool = False - Generate RDKit ETKDGv3 conformers per ligand
- `num_conformers`: int = 50 - Conformers to generate before filtering
- `energy_window`: float = 2.0 - Max relative MMFF94 energy (kcal/mol) for conformer filtering
- `conformer_rmsd`: float = 1.0 - Heavy-atom RMSD cutoff (Angstroms) for Butina clustering
- `conformer_energies`: Optional[tuple] = None - Pre-computed energies as `(TableInfo, "column_name")`
- `cnn_score_threshold`: float = 0.0 - Min CNNscore to accept a pose (0-1); 0 keeps every pose. GNINA has no such filter of its own — it sorts by CNNscore and returns all poses — so any non-zero value is a choice this wrapper adds. CNNscore predicts whether a pose is near the *crystal* pose, which is uncorrelated with binding strength (rho +0.025 over 578 de novo designs) and is trained on natural complexes, so designed pockets with synthetic ligands score low throughout. Raise it to screen for pose confidence; leave it at 0 when every input must be scored.
- `rmsd_threshold`: float = 2.0 - RMSD cutoff (Angstroms) for pose consistency clustering
- `protonate`: bool = True - Add hydrogens with OpenBabel before docking
- `pH`: float = 7.4 - Protonation pH
- `top_only`: bool = True - keep only the best conformer of each (protein, ligand) pair as `<protein>+<ligand>`. False surfaces every conformer that has an accepted pose as `<protein>+<ligand>_conf<k>`, the id of its `docking_summary` row. Docking mode only.
- `rank_conformers_by`: str = "auto" - the `docking_summary` column that picks the best conformer when `top_only` is True. `"auto"` is `pseudo_binding_energy` when every posed conformer of the pair has one and `best_vina` otherwise. Or name one: `pseudo_binding_energy`, `best_vina` or `mean_vina` (lowest wins), `best_cnn_score` or `mean_cnn_affinity` (highest wins; not with `cnn_scoring="none"`). A named metric is never replaced by another: if a posed conformer of a pair lacks it, that pair is a `failure` row of `missing` and the other pairs go on. An unknown name, `pseudo_binding_energy` without `generate_conformers` or `conformer_energies`, and any value other than `"auto"` together with `top_only=False` are refused at construction. Docking mode only.

**Streams**: `structures` (combined protein+ligand PDB files). In `docking` mode with `top_only=True` one entry per (protein, ligand) pair under the id `<protein>+<ligand>`, the best pose of the pair's best conformer, written as `<id>_best.pdb`. With `top_only=False` one entry per conformer that has an accepted pose, under `<protein>+<ligand>_conf<k>` (the `docking_summary` id of that conformer), written as `<id>_best.pdb`; how many conformers dock is known only at run time (a generated set is filtered by energy and clustered by RMSD), so these ids are declared lazily as `<protein>+<ligand>[_conf<#>]`, and a conformer without an accepted pose has no entry. In `score`/`minimize` mode the scored pose under the complex's own id, written as `<id>_score.pdb` or `<id>_minimize.pdb`.

**Tables** (`mode="docking"`):
- `docking_results` - One row per accepted pose across all runs:

  | id | structures.id | compounds.id | conformer_id | run | pose | vina_score | cnn_score | cnn_affinity |
  |----|---------------|--------------|--------------|-----|------|------------|-----------|--------------|

  IDs are unique per pose: `<protein>+<ligand>_conf<k>_r<run>_p<pose>`, the summary id of its conformer with the run and the pose appended. `structures.id`, `compounds.id`, `conformer_id`, `run` and `pose` carry the parts, so no consumer needs to split an id.

- `docking_summary` - One row per (protein, ligand, conformer) group, aggregated across runs:

  | id | structures.id | compounds.id | conformer_id | best_vina | mean_vina | std_vina | best_cnn_score | mean_cnn_affinity | std_cnn_affinity | pose_consistency | conformer_energy | pseudo_binding_energy | best_pose_file |
  |----|---------------|--------------|--------------|-----------|-----------|---------|----------------|-------------------|------------------|-----------------|------------------|-----------------------|----------------|

  IDs are unique per (protein, ligand, conformer) group: `<protein>+<ligand>_conf<k>`, the pair's stream id with the conformer appended. A ligand id may contain `_` (`a_b`+`c` and `a`+`b_c` are `a_b+c` and `a+b_c`). A `+` inside an id is accepted, but two pairs that would join to the same id (`a`+`b+c` and `a+b`+`c` are both `a+b+c`) are refused, at configuration time and again before the first dock.

  - `mean_vina` / `std_vina`: mean and std of the best Vina score per run (not all poses — just the top pose from each independent run). Robust metric for ranking.
  - `mean_cnn_affinity` / `std_cnn_affinity`: same logic for CNN affinity.
  - `pose_consistency`: fraction of runs whose best pose clusters together by RMSD. Measures reproducibility (1.0 = all runs converge).
  - `conformer_energy`: MMFF94 relative strain energy (only with `generate_conformers=True`).
  - `pseudo_binding_energy`: `best_vina + conformer_energy`. Penalizes Vina score by conformer strain. Only populated when conformer generation is enabled.
  - `best_pose_file`: the complex PDB of the group's best pose (lowest Vina score over its runs). With `top_only=True` the structures stream holds one file per pair, so among a pair's conformers the one `rank_conformers_by` picks (by default the lowest `pseudo_binding_energy`, the lowest `best_vina` when a conformer has no energy) is written to `structures/<protein>+<ligand>_best.pdb`; the other conformers' files go to `execution/conformer_poses/<summary id>_best.pdb`. With `top_only=False` every conformer's file is `structures/<summary id>_best.pdb` and `conformer_poses/` is not written.

**Tables** (`mode="score"` / `"minimize"`):
- `scores` - One row per complex that entered, whether or not it could be scored; a complex that failed keeps its row with every score empty (NaN) and its reason in `missing`. Empty cells also where GNINA did not emit that term for the chosen `cnn_scoring`:

  | id | structures.id | compounds.id | vina_affinity | cnn_score | cnn_affinity | cnn_vs | cnn_affinity_variance |
  |----|---------------|--------------|---------------|-----------|--------------|--------|-----------------------|

  - `vina_affinity`: GNINA's `minimizedAffinity` (kcal/mol; lower = stronger).
  - `cnn_score`: CNN pose score (0-1).
  - `cnn_affinity` / `cnn_affinity_variance`: CNN predicted affinity and its variance across the ensemble.
  - `cnn_vs`: CNN virtual-screening score (`cnn_score × cnn_affinity`).

**Tables** (all modes):
- `missing` - IDs removed upstream, and this step's own failures as `failure` rows with `removed_by` the step id and the reason in `cause`. In `docking` mode one row per (protein, ligand) pair that has no entry in the structures stream: a receptor with unparseable coordinates or a ligand that could not be prepared, a binary that failed or timed out, a run that gave no accepted pose, or a `rank_conformers_by` metric a posed conformer lacks. In `score`/`minimize` mode one row per complex (no HETATM with the ligand code, SMILES that doesn't match the extracted skeleton, a binary failure). In `docking` mode the step itself exits non-zero only when no pair gave a structure. | id | removed_by | kind | cause |

**Example**:
```python
from biopipelines.gnina import Gnina

# Basic docking (autobox from crystal ligand in PDB)
docked = Gnina(structures=boltz, compounds=ligand)

# Explicit box
docked = Gnina(
    structures=boltz,
    compounds=ligand,
    center="10.5,20.3,15.0",
    size=25.0
)

# With conformer generation for flexible ligand docking
docked = Gnina(
    structures=boltz,
    compounds=compound_library,
    generate_conformers=True,
    num_conformers=100,
    energy_window=3.0,
    num_runs=10,
    exhaustiveness=64
)

# Pick each pair's conformer by its CNN score instead of the default, or keep every conformer
best_by_cnn = Gnina(structures=boltz, compounds=compound_library, generate_conformers=True,
                    rank_conformers_by="best_cnn_score")
every_conformer = Gnina(structures=boltz, compounds=compound_library, generate_conformers=True, top_only=False)

# Score an already-posed complex without docking (e.g. a Boltz2 co-fold)
boltz = Boltz2(proteins=poi, ligands=Ligand("STI"))
scored = Gnina(structures=boltz, compounds=Ligand("STI"), mode="score")
# scored.tables.scores -> vina_affinity | cnn_score | cnn_affinity | ...

# Local minimization then score
refined = Gnina(structures=boltz, compounds=Ligand("STI"), mode="minimize")
```

---

### Vina

AutoDock Vina docking. `Vina` is the AutoDock Vina engine as its own class: it keeps the box definition, conformer generation, multi-run pose-consistency analysis and the three `mode` values, and drops everything CNN-specific. Use it where the GNINA binary cannot run — GNINA publishes only prebuilt x86-64 releases, while conda-forge builds `vina` for linux-aarch64, so this is the docking route on ARM clusters such as CSCS Daint (GH200).

Two behavioural differences follow from Vina having no CNN:
- `cnn_scoring` and `cnn_score_threshold` are rejected (setting either raises), and `scoring` selects the force field instead: `"vina"` (default), `"vinardo"`, or `"ad4"`.
- **Units.** `vina_affinity` / `best_vina` / `mean_vina` are the Vina score in **kcal/mol**, so **lower is stronger** (same column semantics as Gnina's `vina_affinity`). It is a docking score, not a measured or calibrated affinity — do not compare it to Prodigy's ΔG or convert it to a Kd.
- The output tables carry **no `cnn_*` columns** — they are absent, not blank. `docking_results` is `id | structures.id | compounds.id | conformer_id | run | pose | vina_score`, and `docking_summary` drops the three CNN aggregates.

Vina reads and writes PDBQT only, so the wrapper converts the receptor and ligand in and the poses back to SDF with OpenBabel, re-attaching each affinity so the shared parsing and aggregation stages are unchanged. Vina's own `--autobox` has no reference-file form, so an `autobox_ligand` (or an auto-detected crystal ligand) is turned into an explicit centre and size from its coordinates. Receptor PDBQT is prepared with OpenBabel rather than the Meeko/ADFR route — adequate for ranking and pose selection, but treat absolute affinities with caution.

**Tags**: dock, predict-property, protein, small-molecule, pocket, binding, energy

**Resources**: CPU only. On a whole-node-billed cluster, pack several runs with `Parallel(pack=N)` rather than giving one dock an idle-GPU node. The thread pool is sized from `SLURM_CPUS_PER_TASK`.

**Environment**: `vina` (conda-forge: `vina`, `rdkit`, `openbabel`)

**Installation**: Creates the `vina` env from conda-forge; no binary download.

**Parameters**: as `Gnina`, except `cnn_scoring` / `cnn_score_threshold` are absent and `scoring: str = "vina"` is added. `top_only` and `rank_conformers_by` are as in Gnina, but `rank_conformers_by` cannot name `best_cnn_score` or `mean_cnn_affinity`: Vina's `docking_summary` has no such column, and naming one raises.

**Streams**: `structures` (combined protein+ligand PDB files), under the same ids as [Gnina](#gnina): `<protein>+<ligand>` in `docking` mode (`<protein>+<ligand>_conf<k>` with `top_only=False`), the complex's own id in `score`/`minimize` mode.

**Tables** (`mode="docking"`):
- `docking_results` — one row per accepted pose, id `<protein>+<ligand>_conf<k>_r<run>_p<pose>`: | id | structures.id | compounds.id | conformer_id | run | pose | vina_score |
- `docking_summary` — one row per (protein, ligand, conformer) group, id `<protein>+<ligand>_conf<k>`: | id | structures.id | compounds.id | conformer_id | best_vina | mean_vina | std_vina | pose_consistency | conformer_energy | pseudo_binding_energy | best_pose_file |

**Tables** (`mode="score"` / `"minimize"`):
- `scores` — one row per complex that entered (a complex that failed keeps its row, every score NaN): | id | structures.id | compounds.id | vina_affinity |

**Tables** (all modes):
- `missing` — IDs removed upstream or that failed locally: | id | removed_by | kind | cause |

Column meanings are as `Gnina`'s tables of the same name, minus every `cnn_*` term.

```python
from biopipelines import Vina

protein = PDB("9RTM", ids="rhotag", convert="pdb")
tmr = Ligand(smiles="CN(C)c1ccc2c(c1)OC1=CC(=[N+](C)C)C=CC1=C2c1ccccc1C(=O)O",
             ids="TMR", codes="LIG")

# Box recovered from the crystal ligand already in the PDB
docked = Vina(structures=protein, compounds=tmr, exhaustiveness=16)

# Explicit box, vinardo scoring, repeated runs for pose consistency
docked = Vina(structures=protein, compounds=tmr,
              center="22.87,8.65,55.72", size=24.0,
              scoring="vinardo", num_runs=2)
```

---

### NeuralPLexer

Predicts protein–ligand complex structures from a protein sequence/structure plus a ligand (SMILES/SDF), using a physics-inspired flow-based generative model — no binding-box hint required. Output IDs follow `<protein>+<ligand>_rank<N>`.

**Tags**: predict-structure, sample-ensemble, protein, small-molecule, complex, all-atom

**References**: https://github.com/zrqiao/NeuralPLexer · https://www.nature.com/articles/s42256-024-00792-z

**Resources**: GPU. On Colab it needs a **high-RAM runtime** — the openfold attention + complex sampling OOM the free T4's 12 GB system RAM (verified end-to-end on an A100 high-RAM runtime).

**Environment**: `neuralplexer`

**Installation**: `NeuralPLexer.install()` clones the repo, creates the env, builds the openfold CUDA extension, and fetches the (~8.7 GB) weights bundle from Zenodo.

**Parameters**:
- `structures`: DataStream | StandardizedOutput (required) — Protein structures.
- `compounds`: DataStream | StandardizedOutput (required) — Ligands (SMILES or SDF).
- `n_samples`: int = 16 — Complex samples per pair.
- `num_steps`: int = 40 — Sampling steps.
- `chunk_size`: int = 4 — Attention chunk size (lower to reduce memory).
- `sampler`: str = "langevin_simulated_annealing" — Sampling scheme.
- `cuda`: bool = True — Run on GPU.

**Streams**: `structures` (predicted complex PDB per sample)

**Tables**:
- `confidence`: | id | structures.id | compounds.id | rank | confidence |
- `missing`: | id | removed_by | kind | cause |

**Example**:
```python
from biopipelines.neuralplexer import NeuralPLexer

target = PDB("4ufc", convert="pdb")
lig = Ligand(smiles="CC(=O)Oc1ccccc1C(=O)O", ids="ASA")
cplx = NeuralPLexer(structures=target, compounds=lig, n_samples=16)
```

---

### OpenFold3

Open-source co-folding of proteins, nucleic acids and ligands: the OpenFold Consortium / AlQuraishi lab reimplementation of AlphaFold3, Apache 2.0 with weights included and no access request. It takes the same input axes as `Boltz2`, so the two are interchangeable in a pipeline and a prediction can be cross-checked against a second model.

It is **not** a drop-in replacement for `Boltz2` in what it predicts: there is no affinity head and ligands are non-covalent. For a binding-affinity readout or a covalent probe, use `Boltz2`.

The whole batch folds in one process — `run_openfold predict` takes a single JSON holding every query — so a campaign pays one model load rather than one per complex.

**Tags**: predict-structure, protein, small-molecule, nucleic-acid, complex, all-atom

**References**: https://github.com/aqlaboratory/openfold-3 · method: AlphaFold3 (Abramson et al. 2024, https://www.nature.com/articles/s41586-024-07487-w), which OpenFold3 reimplements; no OpenFold3 paper is cited here

**Resources**: GPU, CUDA 12.1 or newer, at least 32 GB of GPU memory. Weights (~2 GB) are downloaded by `OpenFold3.install()` (`setup_openfold`) into `$HOME/.openfold3`; they are not fetched on first prediction, and `run_openfold` refuses to start without them.

**Parameters**

- `proteins`: Optional[Union[DataStream, StandardizedOutput]] = None - protein sequences
- `ssDNA` / `dsDNA` / `ssRNA` / `dsRNA`: Optional[Union[DataStream, StandardizedOutput]] = None - nucleic-acid sequences. A double-stranded axis becomes two chains, the second the reverse complement
- `ligands`: Optional[Union[DataStream, StandardizedOutput]] = None - compounds, expressed by CCD code or SMILES exactly as `Boltz2` expresses them, so both models receive the same chemistry for the same input
- `msas`: Optional[Union[DataStream, StandardizedOutput]] = None - precomputed alignments (`.a3m`, `.sto` or `.npz`). Supplying them turns the MSA server off; passing `use_msa_server=True` as well is refused. Every protein chain of a query needs one, matched to its protein id: a query with a chain that has none is not predicted, keeps a blank row of `confidence` and lands in `missing` as a `failure`
- `use_msa_server`: Optional[bool] = None - generate alignments through the ColabFold server; unset means on, unless `msas` is given
- `num_diffusion_samples`: Optional[int] = None - structures sampled per query (upstream default 5)
- `num_model_seeds`: Optional[int] = None - random seeds per query (upstream default 1)
- `seeds`: Optional[List[int]] = None - the explicit seed values, when a run has to be reproducible by seed. Contradicting `num_model_seeds` is refused rather than silently resolved
- `output_format`: str = "cif" - "cif" or "pdb". A `.pdb` carries per-atom pLDDT in its B-factor column
- `top_only`: bool = True - keep only the best sample as `<id>`. False surfaces every sample as `<id>_1..N`
- `inference_ckpt_name`: Optional[str] = None - a checkpoint from the upstream list (default `openfold3_p2_v1`)
- `inference_ckpt_path`: Optional[str] = None - an explicit checkpoint file. This is how a privately fine-tuned model, such as a federated checkpoint, runs through the same pipeline
- `template`: Optional[str] = None - a CIF template applied to the protein chains
- `template_chain_ids`: Optional[List[str]] = None - which chains the template applies to
- `low_mem`: bool = False - add upstream's `low_mem` preset, for a GPU that cannot hold the default
- `devices`: Optional[int] = None - GPUs to distribute over (upstream default 1)
- `runner_yaml`: Optional[str] = None - a user-supplied runner YAML merged over the wrapper's. The escape hatch for upstream settings this wrapper does not name; the user's keys win

**Streams**: `structures`, `pae`

`pae` holds one `<id>.json` per structure, under the same ids as `structures`. It is upstream's full confidence file for the kept sample (`*_confidences.json`), copied unchanged: a JSON object with the keys `plddt` (N_atom floats), `pde` and `pae` (N_token x N_token nested lists, so `np.array(json.load(f)["pae"])` is the matrix). `pae` is the predicted aligned error in Å (0 to 32). Rows and columns are OpenFold3's tokens over all chains, following AlphaFold3 tokenization: one per standard residue, one per atom for ligands and modified residues. The stream's map_table lists only the ids whose file exists. It is not declared when a user `runner_yaml` sets `write_full_confidence_scores: false` or a `full_confidence_output_format` other than `json`. The file grows with the square of the token count, because upstream writes both `pde` and `pae` with one indented line per value, so a campaign of many large designs pays for it in disk.

**Tables**:
- `confidence`:

  | id | plddt | ptm | iptm | gpde | has_clash | disorder | ranking_score | seed | sample |
  |----|-------|-----|------|------|-----------|----------|---------------|------|--------|

  One row per query that entered, with blanks where a score is unavailable — a query that
  produced no structure keeps its row, and `missing` says why. `plddt`, `disorder` and
  `ranking_score` are upstream's `avg_plddt`, `disorder` and `sample_ranking_score`; `seed` and
  `sample` name which of the generated samples this row is, so a `top_only` result says which
  one was chosen.

- `missing`:

  | id | removed_by | kind | cause |
  |----|------------|------|-------|

With `top_only=True` (the default) the sample with the highest `ranking_score` is kept as
`<id>` — the best, not the first. With `top_only=False` every sample is surfaced as
`<id>_1..N`.

### PLACER

Atomic-level graph neural network that stochastically regenerates coordinates from a partially corrupted input structure, producing a scored conformational ensemble. The wrapper exposes PLACER's two task modes, **selected by the inputs given** (no explicit mode flag):

- **Ligand-pose mode** (`ligand` provided) — the input structure has the ligand bound as HETATM; PLACER resamples that ligand's pose (and pocket sidechains), keeping any other ligands fixed. The ligand `code` is resolved from the compounds stream at runtime and passed to PLACER's `predict_ligand` selector by `name3` (every copy of that residue is predicted via `predict_multi`). Returns the full score set including the ligand-accuracy terms.
- **Sidechain/apo mode** (`ligand` omitted, `target_res` provided) — a protein residue is the crop center and PLACER resamples the surrounding sidechains; no ligand is predicted. Returns only the lDDT-style confidences (no `prmsd`/`rmsd`/`kabsch`).

PLACER reads PDB and RCSB mmCIF directly — no SDF staging. Each input produces `nsamples` models; the output IDs multiply by a sample index.

> **Prefer mmCIF for multi-ligand / multi-chain inputs.** PLACER's PDB parser asserts that ligand chains and protein chains don't share a chain letter, and raises *"One or more of ligand chains already exist in parsed protein chains"* on structures where they collide (e.g. RCSB `4dtz` as a `.pdb`). Its mmCIF parser doesn't have this limitation, so for crystal complexes feed an mmCIF (`PDB("/path/4dtz.cif")` or any RCSB CIF) rather than the converted PDB. The driver dispatches on file extension automatically (`.cif`/`.cif.gz` → CIF parser, else PDB parser). Note RCSB-sourced mmCIF only — Rosetta/AF3 CIFs are not parsed correctly upstream.

**Tags**: refine-structure, sample-ensemble, protein, small-molecule, pocket, all-atom, sidechain, covalent

**Resources**: GPU (CUDA 12.1 / torch 2.3 stack).

**Environment**: `placer` (dedicated conda env; model weights ship in the cloned repo).

**Installation**: Clones `baker-laboratory/PLACER` and creates the `placer` env. Weights (`weights/PLACER_model_1.pt`) are bundled — no separate download.

**Parameters**:
- `structures`: Union[DataStream, StandardizedOutput] (required) — protein structures (PDB/mmCIF), ligand bound as HETATM for ligand-pose mode.
- `ligand`: Union[DataStream, StandardizedOutput, None] = None — compounds stream with the bound ligand's 3-letter `code`. Presence selects ligand-pose mode.
- `target_res`: Optional[str] = None — crop-center residue `"chain-resno"` (e.g. `"A-149"`). Required in sidechain/apo mode; rejected in ligand-pose mode.
- `exclude_sm`: bool = False — drop all small molecules from the prediction (true apo).
- `nsamples`: int = 10 — samples per input (upstream recommends 50–100).
- `rerank`: Optional[str] = None — rank models/CSV by `"prmsd"`, `"plddt"`, or `"plddt_pde"` (use `plddt`/`plddt_pde` in apo mode; `prmsd` is ligand-only).
- `bonds`: Optional[tuple | List[tuple]] = None — covalent bonds to enforce during resampling (ligand-pose mode only). Each is `(atom1, atom2, length)` with atoms in `<residue>.<atom>` syntax, e.g. `("A145.SG", "LIG.C12", 1.8)` to tether a ligand to a catalytic cysteine (true SNAP-Tag/AGT chemistry). The residue name (`name3`) is resolved per-structure from the input. Pair with [`PDB.break_bond`](inputs_io.md#pdb) to detach the ligand afterward for non-covalent downstream design.

**Streams**:
- `structures` — predicted PDB per sample (prmsd in the B-factor column).
- `compounds` (ligand mode only) — chemistry passthrough of the input ligand (code/SMILES unchanged).

**Tables**:
- `scores` — one row per sample.
  - ligand mode: `id | structures.id | compounds.id | sample | prmsd | plddt | plddt_pde | fape | rmsd | kabsch`
  - apo mode: `id | structures.id | sample | plddt | plddt_pde | fape`
- `missing` — `id | removed_by | kind | cause` for inputs PLACER could not process.

Output IDs: `<structure>+<ligand>_<sample>` (ligand mode) or `<structure>_<sample>` (apo mode).

**Example**:
```python
from biopipelines.placer import PLACER

# Ligand-pose mode: refine a bound ligand's pose
holo = PDB("4dtz")              # ligand bound as HETATM
lig = Ligand("LDP")             # residue code present in 4dtz
poses = PLACER(structures=holo, ligand=lig, nsamples=50, rerank="prmsd")

# Sidechain/apo mode: repack around a residue, no ligand
apo = PDB("dnHEM1_apo")
repacked = PLACER(structures=apo, target_res="A-149", exclude_sm=True,
                  nsamples=50, rerank="plddt")
```

---

### RF3

All-atom co-folding of proteins, nucleic acids and ligands with RoseTTAFold3, the structure-prediction network of RosettaCommons Foundry. It is run through its `rf3 fold` command and takes the same input axes as `Boltz2` and `OpenFold3`, so the three are interchangeable in a pipeline and a prediction can be cross-checked against another model.

It is **not** a drop-in replacement for `Boltz2` in what it predicts: there is no affinity head. Upstream also supports covalent bonds, templates, cyclic chains and ground-truth conformers; none of those is a named parameter here (see the forwarding note below).

The whole batch folds in one process: `rf3 fold` takes a single JSON list holding every complex, so a campaign pays one model load rather than one per complex. A single input that makes upstream raise stops the process, and the inputs after it are not predicted; they appear in `missing` as failures.

RF3 shares the `foundry` environment with `RFdiffusion3`, but not its install: `RF3.install()` creates the environment if it is absent and downloads the RF3 weights into the `RF3` folder of `repositories:`, which `RFdiffusion3.install()` does not do (it fetches only the rfd3 weights, and skips itself once those are present).

**Tags**: predict-structure, protein, small-molecule, nucleic-acid, complex, all-atom

**References**: https://github.com/RosettaCommons/foundry (`models/rf3`) · paper: Accelerating Biomolecular Modeling with AtomWorks and RF3, https://doi.org/10.1101/2025.08.14.670328

**Resources**: GPU. The weights are downloaded by `RF3.install()` (`foundry install rf3`) into `repositories.RF3`; `rf3 fold` is pointed at that folder through `FOUNDRY_CHECKPOINT_DIRS`. Not verified: GPU memory needs, and aarch64 or Colab operation (the wrapper is registered for cluster, Colab and container configs but only the x86-64 cluster is claimed).

**Parameters**

- `proteins`: Optional[Union[DataStream, StandardizedOutput]] = None - protein sequences
- `ssDNA` / `dsDNA` / `ssRNA` / `dsRNA`: Optional[Union[DataStream, StandardizedOutput]] = None - nucleic-acid sequences. A double-stranded axis becomes two chains, the second the reverse complement. The chain type is passed to upstream explicitly, because a string such as "ACGT" is also a valid protein sequence
- `ligands`: Optional[Union[DataStream, StandardizedOutput]] = None - compounds, expressed by CCD code or SMILES exactly as `Boltz2` expresses them. Upstream's own guidance is that a CCD code folds better than a SMILES string
- `msas`: Optional[Union[DataStream, StandardizedOutput]] = None - precomputed alignments for the protein axis, one per protein id (`.a3m`, `.a3m.gz` or `.fasta`). Upstream cannot compute an MSA itself, so without this the proteins fold single-sequence, and then the default early-stopping threshold often abandons them. With `msas`, an input with a protein chain that has no alignment is not predicted and lands in `missing` as a `failure`. Upstream pairs alignments across chains on the fly and does not accept pre-paired ones
- `n_recycles`: Optional[int] = None - trunk recycles (upstream default 10). A single recycle is refused while early stopping is on, because upstream spends the first recycle on the early-stopping check and has none left to fold with
- `diffusion_batch_size`: int = 5 - structures sampled per input from one forward pass of the trunk (upstream default 5). Always passed, so the ids declared for `top_only=False` are exact
- `num_steps`: Optional[int] = None - diffusion sampling steps. Upstream's config default is 50; its README states 200 and says 50 performs similarly and is more than twice as fast. Left to upstream unless given
- `seed`: Optional[int] = None - model seed; left to upstream when None
- `early_stopping_plddt_threshold`: float = 0.5 - an input whose mean all-atom pLDDT after the first recycle falls below this (0 to 1) is abandoned: no structure is written. 0.5 is upstream's default; 0 turns it off
- `top_only`: bool = True - keep only the best sample of each input as `<id>` (upstream's ranking score, 0.8 * ipTM + 0.2 * pTM minus a penalty for clashing chains). False surfaces every sample as `<id>_1..N`, N being `diffusion_batch_size`
- `checkpoint_path`: Optional[str] = None - an explicit checkpoint file (upstream `ckpt_path`), for a model other than the one `RF3.install()` fetched
- any other keyword is forwarded to `rf3 fold` as a Hydra `key=value` override, for example `template_selection`, `cyclic_chains` or `residue_renaming_dict`. Only the parameters above are validated or reflected in the declared outputs

**Streams**: `structures`, `pae`

`structures` holds one `<id>.cif` per complex. The B-factor column carries per-atom pLDDT on a 0 to 1 scale (upstream's `annotate_b_factor_with_plddt`, which its config turns on).

`pae` holds one `<id>.json` per structure, under the same ids as `structures`. It is upstream's full confidence file for that sample (`*_confidences.json`), copied unchanged: a JSON object with the keys `atom_chain_ids`, `atom_plddts` (N_atom floats, 0 to 1), `pae` (N_token x N_token nested lists, so `np.array(json.load(f)["pae"])` is the matrix, in Å from 0 to 32), `token_chain_ids` and `token_res_ids`. Rows and columns are RF3's tokens over all chains, one per standard residue and one per atom for ligands and modified residues. `token_res_ids` is a running token index in upstream's source, not a residue number. The stream's map_table lists only the ids whose file exists.

**Tables**:
- `confidence`:

  | id | plddt | ptm | iptm | pae | pde | has_clash | ranking_score | early_stopped | seed | sample |
  |----|-------|-----|------|-----|-----|-----------|---------------|---------------|------|--------|

  One row per id that entered, with blanks where RF3 reported nothing. `plddt`, `pae` and `pde` are upstream's `overall_*` means; `ptm`, `iptm`, `has_clash` and `ranking_score` are its own. Upstream reports an `iptm` of 0 for a single chain. `seed` and `sample` name which of the generated samples the row is, so a `top_only` result says which one was chosen. An input stopped early has `early_stopped` True and blanks elsewhere.

- `missing`:

  | id | removed_by | kind | cause |
  |----|------------|------|-------|

  An input stopped by the early-stopping threshold is a `filter` row: you set the threshold, so its absence is expected. An input for which `rf3 fold` wrote no structure for another reason is a `failure` row.

With `top_only=True` (the default) the sample with the highest `ranking_score` is kept as `<id>`. Newer upstream releases record an early stop as a row in `<id>_ranking_scores.csv`, older ones as a `score.csv`; both are read.

```python
complexes = RF3(proteins=binders, ligands=Ligand("HEM"), msas=binder_msas)
```

---

### Protenix

All-atom co-folding of proteins, nucleic acids and ligands with ByteDance's Protenix, a trainable reproduction of AlphaFold3, run as `protenix-v2` (about 464 M parameters) unless another checkpoint is named. It takes the same input axes as `Boltz2` and `OpenFold3`, so a prediction can be cross-checked against another model by changing the tool name.

**License.** The upstream code and the other checkpoints are Apache-2.0. The upstream README states that the weights of `protenix-v2` are proprietary and confidential, are not released under an open-source license, and may not be reproduced, distributed, sublicensed or transferred to a third party without the rights holder's prior written consent. Upstream's own download URL for `protenix-v2.pt` answered HTTP 403 when checked on 2026-10-28, while the checkpoints of the other models downloaded; whether it is reachable from the cluster has not been checked. Obtain the weights from the rights holder and put them at `<ProtenixCache>/checkpoint/protenix-v2.pt`. The wrapper does not fetch from third-party mirrors. Do not copy or share the checkpoint, and read the README before publishing anything built on it.

It is **not** a drop-in replacement for `Boltz2` in what it predicts: there is no affinity head, and covalent bonds are the explicit atom pairs of `covalent_bonds`, not a selection resolved against the sequence. Two upstream features are deliberately not wired. Contact and pocket constraints exist only for the older `protenix_base_constraint_v0.5.0` checkpoint (upstream's model table lists the constraint column as unsupported for `protenix-v2`). Templates need kalign plus per-chain search results or a sequence database.

**Alignments.** Three modes. With `msas`, each protein chain gets its alignment as an unpaired `.a3m`. Upstream's ColabFold note says Protenix pairs chains through species information that ColabFold-style alignments lack, so such a complex is predicted without cross-chain pairing. By default upstream's MSA service at protenix-server.com searches every protein chain, which means the sequences leave the cluster; supply `msas` for anything unpublished. With `use_msa_server=False` and no `msas`, every chain is folded from its single sequence, which upstream warns can degrade accuracy.

**Failures.** `protenix pred` records a failed query in `ERR/<id>.txt` and still exits 0, so the wrapper reads those files: a query with no structure gets a `failure` row in `missing` carrying the first line of upstream's error. The best-known case is `protenix-v2` refusing a complex of more than 2560 tokens. A query the wrapper cannot even write (a residue letter upstream does not accept, a protein chain with no alignment in `msas`, a covalent bond naming an entity the query does not have) is skipped the same way, with the reason.

**Tags**: predict-structure, protein, small-molecule, nucleic-acid, complex, all-atom, covalent, sequence-only

**References**: https://github.com/bytedance/Protenix · Protenix-v2: Zhang et al., bioRxiv 2026, https://doi.org/10.64898/2026.04.10.717613 · reproduces AlphaFold3 (Abramson et al. 2024, https://www.nature.com/articles/s41586-024-07487-w)

**Resources**: GPU. Upstream's cost table lists peak memory of about 6 GB at 500 tokens, 18 GB at 1000 and 67 GB at 2000; request `gpu="A100"` or larger for complexes beyond a few hundred tokens. x86-64 only: upstream pins `triton==3.3.1`, which PyPI ships for x86-64 only, so it is not set up on Daint (aarch64).

**Environment**: `protenix` (dedicated env, Python 3.11, `pip install protenix` from the official PyPI index).

**Installation**: `Protenix.install()` builds the env and downloads upstream's CCD files and the `protenix-v2` checkpoint into `ProtenixCache` (`$PROTENIX_ROOT_DIR`, the home directory when the config sets no cache), so the first prediction needs no network. When the checkpoint cannot be downloaded the install still succeeds, with a warning that names where the file goes; any other `model_name` runs, and upstream downloads that checkpoint on its first use.

**Parameters**

- `proteins`: Optional[Union[DataStream, StandardizedOutput]] = None - protein sequences. Letters other than the 20 standard residues and `X` fail that query
- `ssDNA` / `dsDNA` / `ssRNA` / `dsRNA`: Optional[Union[DataStream, StandardizedOutput]] = None - nucleic-acid sequences in upper case (A/C/G/T/N for DNA, A/C/G/U/N for RNA). A double-stranded axis becomes two chains, the second the reverse complement
- `ligands`: Optional[Union[DataStream, StandardizedOutput]] = None - compounds, expressed by CCD code or SMILES exactly as `Boltz2` expresses them, so both models receive the same chemistry for the same input. A single-atom ion is a ligand with its CCD code (`MG`), which upstream builds the same way as its `ion` entity
- `msas`: Optional[Union[DataStream, StandardizedOutput]] = None - precomputed `.a3m` alignments, from `MMseqs2(output_format="a3m")` or `MSA`. Every protein chain of a query needs one; a query with a chain that has none is skipped and listed in `missing`. A csv stream is refused at configuration. `use_msa_server=True` together with `msas` is refused, since upstream then searches every chain of the query and overwrites the supplied alignments
- `use_msa_server`: Optional[bool] = None - search alignments through upstream's service; unset means on, unless `msas` is given. False without `msas` folds every chain single-sequence
- `model_name`: str = "protenix-v2" - one of `protenix-v2`, `protenix_base_default_v1.0.0`, `protenix_base_20250630_v1.0.0`, `protenix_base_default_v0.5.0`, `protenix_mini_default_v0.5.0`, `protenix_tiny_default_v0.5.0`. Upstream's `--use_default_params` is always on, so cycles and diffusion steps follow the model (10 and 200; 4 and 5 for mini and tiny). The constraint and ESM checkpoints are not offered
- `num_diffusion_samples`: Optional[int] = None - structures sampled per seed (upstream default 5)
- `seeds`: Optional[List[int]] = None - random seeds (upstream default 101). A query yields seeds times samples structures
- `top_only`: bool = True - keep only the structure with the best `ranking_score` as `<id>`. False surfaces every structure as `<id>_1..N`, ordered by seed and then by upstream's sample index
- `covalent_bonds`: Optional[List[Dict]] = None - upstream's covalent bonds, each a dict with `entity1`, `position1`, `atom1`, `entity2`, `position2`, `atom2` and optionally `copy1` and `copy2` (both or neither). Entities are numbered from 1 in the order of the query: proteins, ssDNA, dsDNA (strand, then reverse complement), ssRNA, dsRNA, ligands, each axis in the order of its records. The same bonds apply to every query. For a ligand given as SMILES, `atom` is an index into the SMILES or an element with its occurrence (`C3`); for a CCD ligand it is the CCD atom name. Polymer-polymer bonds are not reliably handled upstream
- `include_pae`: bool = True - write the PAE matrix of every kept structure to the `pae` stream. This makes upstream dump a full-data JSON per sample (`--need_atom_confidence`), which is large; the wrapper reads the matrix out of it and deletes the JSONs
- `no_kernels`: bool = False - run the triangle operations with plain PyTorch (`--trimul_kernel torch --triatt_kernel torch`) instead of the cuEquivariance kernels, for a GPU those do not support

The layer norm always runs in plain PyTorch (`LAYERNORM_TYPE=torch`, exported by the generated script): upstream's default fused layer norm is JIT-compiled at the first import and needs nvcc and `CUDA_HOME`, which an environment holding only the pip CUDA wheels lacks. An exported `LAYERNORM_TYPE=fast_layernorm` selects the compiled one again.

**Streams**: `structures`, `pae`

`structures` holds one `<id>.cif` per predicted complex. `pae` holds one `<id>.npz` per structure under the same ids, with the full N_token x N_token predicted aligned error under the array key `pae` (float32, Å), the layout Boltz2 writes, which `IPSAE` and other readers of a `pae` stream accept. That Protenix's token order lines up with the residues of its CIF, as Boltz2's does, rests on both tokenizing like AlphaFold3 and has not been run. Rows and columns are Protenix's tokens over all chains: one per standard residue, one per atom for ligands and modified residues. The stream's map_table lists only the ids whose matrix could be read, and it is not declared with `include_pae=False`.

**Tables**:
- `confidence`:

  | id | plddt | ptm | iptm | gpde | has_clash | disorder | ranking_score | seed | sample |
  |----|-------|-----|------|------|-----------|----------|---------------|------|--------|

  One row per query that entered, with blanks where a score is unavailable: a query that produced no structure keeps its row, and `missing` says why. The scores are those of upstream's `summary_confidence` JSON; `plddt` is on a 0 to 100 scale, `ptm` and `iptm` from 0 to 1, and `ranking_score` is upstream's `0.8 iptm + 0.2 ptm + 0.5 disorder - 100 has_clash`, so a structure with a clash scores below -99. `disorder` is a constant 0 in upstream's current code. `seed` and `sample` name the structure the row describes, `sample` being upstream's index within the seed (0 is the best of that seed).

- `missing`:

  | id | removed_by | kind | cause |
  |----|------------|------|-------|

With `top_only=True` (the default) the structure with the highest `ranking_score` over all seeds and samples is kept as `<id>`. With `top_only=False` every structure is surfaced as `<id>_1..N`.

**Example**:
```python
from biopipelines import Protenix, MMseqs2, IPSAE

protein = Sequence("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH", ids="Hb")
msa = MMseqs2(sequences=protein, output_format="a3m")

# Supplied alignments, five seeds of five samples each, best structure kept
fold = Protenix(proteins=protein, msas=msa, seeds=[1, 2, 3, 4, 5], ligands=Ligand("ATP"))

# The pae stream has the Boltz2 layout, which IPSAE reads
iface = IPSAE(structures=fold, pae=fold)
```
---
