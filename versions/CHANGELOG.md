# Changelog

All notable changes to biopipelines and its tool wrappers are recorded here. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/) for the framework version (`biopipelines.__version__` in `biopipelines/__init__.py`).

Each release section has two subsections:

- **Framework** — changes to `biopipelines/`-internal infrastructure (`Pipeline`, `BaseConfig`, debug capture, CLI, shared `pipe_scripts/` helpers). Tracked by the framework version.
- **Tools** — wrapper changes that bumped a tool's `TOOL_VERSION`. Each bullet names the tool, the old → new version, and a one-line note. Mirrors the `current` field in `tool_changelog.yaml`.

The pre-commit hook (`versions/check_tool_edits.py`) refuses any commit that modifies a wrapper or its mapped pipe scripts without bumping the tool's `current` in `versions/tool_changelog.yaml` AND adding a matching bullet under `[unreleased] → Tools` here.

## [unreleased]

A patch release on top of `1.3.0`: opt-in container/environment execution switching for external tools, a JSON-config refactor of the resolver-script argument passing, and a `grain` parameter for `Panda.merge`. Framework version `1.3.1`.

### Framework

**Container ⇆ environment execution switching**

- New `biopipelines_io.container_argv_prefix(prefix)` splits a tool's `container_prefix()` string into argv tokens for a helper to prepend to a binary subprocess, so the external binary runs inside the configured `.sif` while the helper's own (biopipelines-importing) Python keeps running in the activated host env.
- New `BaseConfig.warn_container_unsupported()` emits a runtime warning when a container is configured for a tool that runs its model in-process or dispatches host-env sub-processes and therefore cannot use it; such a tool falls back to environment mode instead of silently ignoring the setting.
- The developer manual gains a "Containerization" section: the `container_prefix()` seam and the two tool architectures (heavy command emitted directly in bash vs. binary invoked inside a host helper), that environment activation still happens in container mode (the prefix wraps only the binary), the ESMFold-style split / warn fallback for in-process tools, and the per-iteration container-startup caveat with the persistent-instance (`apptainer instance start/exec/stop`) escape hatch.

### Tools

- **LASErMPNN** (new, 1.0) — ligand-conditioned inverse folding with all-atom sidechain packing (polizzilab/LASErMPNN). Reads the ligand from the input PDB's HETATM (no ligand-code argument); `fixed`/`redesigned` positions (string or `(TableInfo, column)` reference) are stamped into input B-factors and run with `--fix_beta`; outputs full-atom designed complexes (`structures`) and a content-bearing `sequences` stream (`id, structures.id, sequence, score`). Ships gpu/cpu/colab env variants; Colab-verified on a T4 across plain design, fixed/redesigned, `repack_only`, `ignore_ligand`, the soluble model, and the ALA/GLY budget.
- **Container-mode execution** added (heavy binary runs in the `.sif` when a container image is configured, host env otherwise): **ESMFold**, **NeuralPLexer**, **DiffDock**, **PocketGen** (heavy command emitted directly in bash); **XTB**, **APBS**, **DSSP**, **FPocket**, **P2Rank**, **Reduce**, **VespaG**, **RTMScore**, **ThermoMPNN** (prefix threaded into the helper and prepended to the binary subprocess). All 1.0 → 1.1.
- **Container-configured warning** added for tools that run their model in-process or dispatch host-env sub-processes and do not support container execution: **AF2BIND**, **BioEmu**, **OpenMM**, **ProLIF**, **Prodigy**, **Frame2Seq**, **PLM_Sol**, **PLACER**, **DynamicBind**, **Gnina**, **GEMS**. All 1.0 → 1.1.
- **LigandMPNN**, **ProteinMPNN**, **RFdiffusion**, **RFdiffusionAllAtom**, **RFDAA_PrepareLigand**, **RFdiffusion2**, **HBDesigner** (1.0 → 1.1) — resolver-script arguments moved from positional CLI into a JSON config file; generated bash changed, no user-facing parameter change. RFdiffusion2 and HBDesigner are newly tracked in `tool_changelog.yaml`.
- **Panda** (1.0 → 1.1) — `Panda.merge` gains a `grain` parameter (`finest`/`coarsest`) selecting which input's provenance level drives the output rows; `how` now defaults to `None` (derived from `grain`).

## [1.3.0] — 2026-06-26

A feature release on top of `1.2.0`: support for LSF and PBS/Torque batch schedulers alongside SLURM, a reworked stitch/combinatorics provenance scheme, and assorted stream-resolution and tool fixes. Contains one backward-incompatible API change (see below).

### Framework

**Batch schedulers (LSF, PBS/Torque)**

- The submission path is no longer SLURM-only. A new `biopipelines/schedulers.py` introduces a `SchedulerBackend` abstraction with `slurm`, `lsf`, and `pbs` backends; `Pipeline` generates native `#SBATCH` / `#BSUB` / `#PBS` directives for the configured scheduler, and the `submit` wrapper drives `sbatch` / `bsub` / `qsub`, capturing job ids to chain dependent batches. The same `Resources(...)` call is portable across all three.
- `machine.scheduler.name` now accepts `lsf` and `pbs` in addition to `slurm`, `colab`, and `none`; `bp-config auto` detects the scheduler from `sinfo` / `bsub` / `qsub` on PATH. GPU and memory specs are translated best-effort per scheduler (LSF memory to MB plus a matching `rusage`, walltime to `[HH:]MM`; PBS as separate `-l` requests); a constraint-style GPU spec with no native equivalent falls back to a plain GPU-count request with a printed warning. Scheduler-native knobs pass through `Resources(**options)` under the active scheduler's options key (`slurm_options` / `lsf_options` / `pbs_options`). Generated scripts are stemmed by scheduler (`slurm_batch*.sh`, `lsf_batch*.sh`, `pbs_batch*.sh`).
- The `submit` and `run` wrappers now export `BIOPIPELINES_CONFIG_VARIANT` and deep-merge the user overlay (`.config.<variant>.yaml`) when reading the scheduler block, so the bash and Python sides always agree on the active variant and on the resolved scheduler/modules/init.
- **Breaking:** `Pipeline.slurm()` is renamed to `Pipeline.generate_job_scripts()`, the `slurm_script` attribute to `job_script`, and `ConfigManager.get_slurm_modules()` to `get_scheduler_modules()`. The old names are removed (no aliases). Pipeline scripts that called `p.slurm()` directly must switch to `p.generate_job_scripts()`; the common `with Pipeline(...)` auto-submit path is unaffected.

**Combinatorics and stitching**

- `assemble_provenance_id()` is the single source of truth for `+`/`_` id assembly. StitchSequences over tool axes now emits `+`-joined provenance ids (e.g. `pdb_3_5+pdb_3_2`) plus `sequences_k.id` columns rather than positional `_`-joined ids, while raw/marker axes keep `_<i>` fallbacks; the config-time predictor mirrors the runtime exactly. Combinatorics accepts a bare filtered `DataStream` as an axis source.

**Stream resolution and id filtering**

- Id filtering supports pattern selection over `map_table` rows with bracket-preserving expansion. `map_table` is treated strictly as a runtime artifact: config-time `create_map_table` is dropped from `get_output_files` (with a documented exception for `Load`).

### Tools

- **DNAEncoder** — its stream is the content table (`sequences/sequences.csv`) and emits a canonical `sequence` column.
- **ProteinMPNN** — adds SolubleMPNN alias for soluble_model=True; switched parameter default back to False.

## [1.2.0] — 2026-06-11

A large release on top of `1.1.1`: a config-system overhaul, a reworked output/datastream contract, ~30 new tool wrappers, full Google Colab support alongside the cluster path, and a self-hosted MMseqs2 CPU search server.

### Framework

**Configuration system**

- New variant-aware config layout: `config.<variant>.yaml` ships committed defaults; user edits land in a gitignored `.config.<variant>.yaml` that is deep-merged on top at load time, so pulling repo-config changes never clobbers local settings. `machine.env_manager` and `machine.scheduler` are now dicts; folder layout moved from CamelCase to snake_case; config and execution output folders are separated.
- New `bp-config` CLI: `path` (prints the active config path — pipe it into an editor with `$EDITOR "$(bp-config path)"`, `--base` for the committed file), `edit`, `list`, `show` (faithful YAML dump with an arrow-key variant picker), `set`, and `auto`. `set`/`auto`/`edit` write only the overlay (`set` validates against the base, `edit` shows the merged tree and saves the diff, `auto` host-probes the machine block and diffs against base+overlay). `get`/`set` work non-interactively for scripting.
- `bp-config edit` opens an interactive arrow-key TUI (`biopipelines/config_editor.py`, prompt_toolkit-based): expand/collapse, in-line value editing with type-preserving coercion, search, and atomic save. Round-trip uses `ruamel.yaml` so comments and key ordering survive. Every config write takes a timestamped `<name>.<YYYYMMDD-HHMMSS>.bak`, so successive edits never clobber the previous backup.
- `prompt_toolkit>=3` and `ruamel.yaml>=0.17` are now required dependencies (not optional extras), so `bp-config edit` is always available. `rdkit` is also promoted to a core dependency. `pyproject.toml` declares minimum version floors for the framework's Python deps.

**Output and datastream contract**

- Per-tool outputs are standardized as typed streams. `DataStream.files` widens to `Union[str, List[str]]`: a single non-empty string now denotes one shared artifact covering all ids (e.g. a multi-record FASTA), in addition to the empty-list (value-based) and per-id-list forms. Every consumer handles the shared-file form explicitly — no silent copy-whole-file fallback.
- New `biopipelines/stream_slicers.py`: a registry of format-keyed slicer and merger functions (ships FASTA and CSV). Unknown formats raise rather than copying a whole file, so a filter step can never leak stale records downstream. FASTA header tokenizer handles rich upstream headers (`>4LCD, score=...`) while preserving NCBI-style `sp|P12345|FOO` ids.
- Canonical per-tool output layout unified as `_configuration/`, `_execution/`, `<stream>/`, `tables/`, `_extras/`.
- Missing-manifest cascade owned by `BaseConfig`: multi-axis propagation and provenance-aware completion remap; `missing.csv` carries a `kind` column (`FAILED` only on local failures, not upstream-filtered ids). Excusal matches on each path's owner id (threaded through the `<id>` substitution), not the filename basename, so templates that place `<id>` off the stem are excused correctly, including lazy fan-outs the config-time remap can't pre-expand.
- `get_mapped_ids` provenance overhaul: O(N²) → O(N) id-matching; closest-siblings-only design-group matching (stops a cross-design Cartesian explosion in StitchSequences); table lookup unified on `get_mapped_ids`. `pipe_check_completion` runs under Python 3.8 tool envs (`typing.Tuple`).

**Pipeline orchestration**

- `Parallel()` context manager + DAG batch-dependency machinery for splitting a run across SLURM jobs.
- `Pool` tool: gather N parallel runs into one `StandardizedOutput` (`runs=` list, `streams=` subset, `recount_prefix` opt-in renumber).
- `LoadMultiple` nests its `Load` steps under `Folder("LoadMultiple")` by default (`folder=None` opts out); accepts any-length step prefix.
- Internal-tool support: `_internal` flag hides a tool under `.internal/` with split public/exec numbering; `Folder()` grouping.
- Panda provenance: chained stream-map writes rename `<stream>.id` → `<stream>.-1.id` (negative integers reserved for chained-generation tracking, can't collide with Bundle's non-negative `<stream>.N.id`), so `Mock → Panda → Panda` chains stay resolvable. Implicit source tagging on concat; `zscore` op (population std, optional by-group + sign flip).

**Parsing and selections**

- mmCIF parsing in the shared structure reader (CIF-correct loop tokenization and null handling, prefer `auth_atom_id`); centralized fixed-column field accessors as the single source of truth for ATOM/HETATM layout, routed through every per-pipe-script reader.
- Insertion codes treated as residue identity throughout.

**Shell-safety**

- Systematic shell-safety pass across all tools: free-form string params and shell-unsafe input validated in `Pipeline`, config, and the tool base; tightened machine-enum validation; description escaping counted as its own layer. Layers documented in the developer manual.

**Colab support**

- Per-tool `.colab.yaml` env variants alongside cluster envs; install paths split by scheduler. `Pipeline` does not auto-enable `local_output` on Colab; Drive mount has no exec bit (executables kept local, weights ok on Drive). `MPLBACKEND=Agg` forced before matplotlib imports in vendored stacks.
- `machine.email` is now a single string (default `""`, no SLURM-notification spam), replacing the legacy `machine.emails` dict (still honoured for back-compat).

**MMseqs2 CPU search server (cluster)**

- Self-hosted CPU `colabfold_search` server (LocalColabFold's bundled mmseqs, not the GPU build): per-mode DB build, page-cache warm-up with mlock, right-sized resources (~900 GB / 64 cpu for the ~746 GB index), NFS-atomic single-server lock via `mkdir`, batched many-against-many queue, orphaned-job recovery, idle shutdown, and a multi-server pool. The client auto-starts the CPU server by default.

### Tools

**Added (~30 new wrappers)**

- Folding / structure: ESMFold, Frame2Seq, NeuralPLexer.
- Docking / pose: DiffDock, DynamicBind, PocketGen, P2Rank, PLACER, AF2BIND.
- Scoring: XTB, RTMScore, GEMS, ADMETAI, Prodigy.
- Cheminformatics / IO: RDKit, OpenBabel, UniProt, Reduce.
- Interactions / structure analysis: PLIP, ProLIF, FPocket, DSSP, APBS, OpenMM.
- Solubility / aggregation: Aggrescan3D, PLM_Sol.
- Ensembles: BioEmu, EnsembleAnalysis, Consensus.
- Mutation-effect: ThermoMPNN, VespaG.
- Structure generation: RFdiffusion2 (atomic active-site scaffolding via the bakerlab Apptainer image, HPC-only), RFdiffusion3.
- Inputs & I/O: Scripting (run a user-authored two-phase configuration/execution script as a typed pipeline step that declares its own streams/tables); Mock (feature/results testing).
- Every new tool ships cluster + Colab env variants and a per-tool install/smoke pipeline; config-time parameter coverage asserts each public constructor parameter reaches the emitted artifacts.

**Changed**

- RFdiffusion3 1.0 → 2.0: cover the full RFD3 paper example set. New input-spec keys (`symmetry`, `unindex`, `select_unfixed_sequence`, `select_partially_buried`, `redesign_motif_sidechains`, `ori_token`, `infer_ori_strategy`, `is_non_loopy`, `plddt_enhanced`, `partial_t`) and inference-sampler overrides (`cfg`, `cfg_scale`, `step_scale`, `noise_scale`, `num_steps`, `center_option`, `seed`). Symmetry auto-selects the symmetry sampler. Adds symmetric-oligomer, diffused-ligand+CFG, protein-DNA, and atomic-motif enzyme design to the existing de-novo / binder / rigid-ligand paths.
- DistanceSelector 1.0 → 2.0: selection-based `within`/`beyond` references and per-structure `(TableInfo, column)` references; excuses upstream-filtered ids in the missing-manifest cascade.
- PDB: `PDB.remove(selection)` drops residues by number range *or* residue name (e.g. `NAP`, `EDO`, digit-bearing CCD codes), mixing the two and rewriting dependent ANISOU/CONECT/LINK records to stay self-consistent; `PDB.rotate_bond` (torsion rotation about a bond) and `PDB.break_bond` (covalent-ligand support); `PDB.chain` accepts a chain list / `"all"` / `split_chains`; chain-qualified atom selections (`A141.CB`).
- Ligand: extended (1–5 char) CCD codes via one validator everywhere; `Ligand(structures=, codes=)` extract mode (posed coords from bound HETATM + chemistry from lookup); bare `ligand="LIG"` sugar; prefers the RCSB/PubChem ideal SDF; `output_format` removed (structures stream is SDF, redirects to OpenBabel). Code-only `Ligand(code=...)` now fails loud with an actionable error instead of emitting `smiles: nan` into a downstream YAML.
- Selection grammar: `+` union for atom-set references (`LIG.B41+LIG.B42`, `LIG.O3+LIG.Cl`), `LIG.B`-prefix name matching.
- Boltz2: `top_only` surfaces all diffusion samples with compact lazy ids and suffix-aware provenance; `predict_atom_names`; covalent-linkage selection strings; Mutagenesis derives synthetic per-mutant MSAs so Boltz2 skips per-variant server queries.
- PoseBusters: SMILES-templated ligand bond orders (fixes charged/conjugated dyes); `exclude=` drops checks from `all_pass`; fixed `all_pass` always reading False (`numpy.bool_` not counted).
- OpenBabel 1.0 → 1.1: read SMILES from the map_table when the compounds stream has no files; accept structures input.

## [1.1.1] — 2026-04-30

### Framework

- Add `__version__` to `biopipelines/__init__.py` and switch `pyproject.toml` to a dynamic version that reads from it.
- Add `TOOL_VERSION` class attribute to `BaseConfig` (default `"1.0"`) and explicit `TOOL_VERSION = "1.0"` on every concrete tool wrapper.
- Record `biopipelines_version` and per-tool `tool_version` in `_debug_capture/tools.json` (Reviewer #2 §6, B6 outstanding sub-point on tool versioning).
- Introduce `tool_changelog.yaml` as the source of truth for tool wrapper versions and the file map the pre-commit hook checks against.
- Add `scripts/check_tool_version_bumps.py` and `.pre-commit-config.yaml` to enforce that any wrapper change is accompanied by a `TOOL_VERSION` bump and a `CHANGELOG.md` entry.

### Tools

- All wrappers initialised at `TOOL_VERSION = "1.0"`.
