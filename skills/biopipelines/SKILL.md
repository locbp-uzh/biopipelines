---
name: biopipelines
description: >-
  Author and run BioPipelines (locbp-uzh, CSBJ 2026) computational protein &
  ligand design pipelines on a GPU. BioPipelines gives 86 comp-bio tools one
  common Pipeline/Resources API: structure prediction and docking (AlphaFold,
  Boltz2, ESMFold2, DiffDock, GNINA, NeuralPLexer, DynamicBind), de novo
  generation (RFdiffusion 1/2/3, RFdiffusionAllAtom, BoltzGen, PocketGen,
  HBDesigner), inverse folding / sequence design (ProteinMPNN, LigandMPNN,
  LASErMPNN, Frame2Seq, ThermoMPNN), MD (OpenMM), pocket detection (FPocket,
  P2Rank, AF2BIND), cheminformatics (RDKit, OpenBabel, ADMET-AI), and
  interaction/stability analysis (PLIP, ProLIF, PoseBusters, Prodigy). Use for
  protein binder/enzyme design, inverse folding, ligand docking, compound-library
  (incl. covalent) screening, and codon optimization. Load when the user wants a
  BioPipelines workflow, names one of these tools, or wants to run such a
  design/screening pipeline on GPU.
---

# BioPipelines

[BioPipelines](https://github.com/locbp-uzh/biopipelines) (Quargnali &
Rivera-Fuentes, LOC-BP UZH; CSBJ [10.34133/csbj.0129](https://spj.science.org/doi/10.34133/csbj.0129))
is a Python framework that puts 86 protein/ligand-modeling tools behind one declarative `Pipeline` API — the authoritative count lives in `docs/tool_index.md` and is CI-enforced against the sources. You describe a workflow as a chain of tools; the framework generates and runs the per-tool scripts, tracks IDs through typed data streams, and materializes outputs.

This skill ships *inside the repo* so that adding it to any agent is one step
(import the skill from this repository) and it tracks the framework as it
evolves. It covers **using** the framework and running it on **any single-node
GPU backend** (a managed-container provider like Modal/RunPod, a plain Docker
GPU box, or an interactive Slurm GPU shell) via the generic `container` config
variant.

## Absorb, don't recite

The `llm/` folder in this repo is the authoritative, author-maintained agent
contract. **Read those first; they override anything paraphrased here.** Do not
duplicate their content into your own notes — point at them:

- `llm/pipelines.md` — how to author a `Pipeline`, the tool/data-stream model, ID tracking.
- `llm/development.md` — conventions, gotchas, how tools are structured.
- `llm/cluster.md`, `llm/colab.md` — two of the three backend references (see `llm/daint.md` for the third). The repo ships five config variants: cluster, colab, daint, container and local.
- `llm/daint.md` — CSCS Alps/Daint: access, storage, and which tools have actually been verified by running there. `config.daint.yaml`'s `environments:` block is a routing table, not a verification record.
- `docs/tool_index.md`, `docs/tool_reference.md` — the full tool catalog and per-tool signatures.
- `references/container_backend.md` (in this skill) — the generic single-node GPU backend.
- `references/daint_backend.md` (in this skill) — the CSCS Alps/Daint backend (aarch64, venv, Container Engine).

## The API in one screen

```python
from biopipelines import Pipeline, Resources, Sequence, Ligand, UniProt, Boltz2

with Pipeline("Project", "job_name", description="..."):
    Resources(gpu="A100", memory="64GB", time="6:00:00", cpus=8)
    prot = Sequence("MSEQ...", type="protein")           # or UniProt("Q15436") to fetch it
    lig  = Ligand(smiles="C[N+]1=C(...)...")
    Boltz2(proteins=prot, ligands=lig, output_format="mmcif")
```

Every tool takes its inputs **by keyword**. `Boltz2`'s first positional parameter is `config` (a raw YAML string), so `Boltz2(prot, lig, ...)` binds the protein to `config` and the ligand to `proteins` — always write `proteins=` / `ligands=`. Check the real signature in `docs/tool_reference.md` before passing anything positionally; an unknown keyword is swallowed by `**kwargs` rather than rejected.

Each `Resources(...)` call **opens a new batch** (one scheduler job), inheriting whatever it does not specify from the previous batch. So call it once per resource profile — a CPU-only stage after a GPU stage gets its own `Resources(gpu="none", ...)` rather than inheriting the GPU and walltime of the heaviest stage. Inside a plain `Parallel()` block every sibling iteration **must** call `Resources()` to open its own batch; inside `Parallel(pack=N)` it is called **once** to describe the whole node allocation and `Run(...)` delimits each task.

`Tool.install()` is called **inside** the `with Pipeline(...)` block (a bare `Boltz2.install()` at module scope is a silent no-op) — or just use `bp-warm` (below), which wraps it for you.

## Running on any single-node GPU host (the `container` backend)

The repo ships a generic single-node variant so you do not hand-author a config
per provider. Set three env vars and point one config line at your persistent
mount:

```bash
export BIOPIPELINES_CONFIG_VARIANT=container   # select config.container.yaml
export BIOPIPELINES_OTF=1                       # run tools inline (no scheduler)
export BIOPIPELINES_LOCAL_OUTPUT=0              # honor configured output dir, NOT cwd
```

In `config.container.yaml` edit `folders.base.root:` to your persistent mount (default `/workspace`); the config's own paths — `home`, `data`, `scratch`, weight caches, `biopipelines_output` — all derive from it. The one path that does **not** is the micromamba env root: that is `MAMBA_ROOT_PREFIX`, hardcoded to `/workspace/micromamba` in `Dockerfile.container`, so a different `root:` needs `MAMBA_ROOT_PREFIX=<root>/micromamba` exported alongside it (see `references/container_backend.md`). First, warm the tools you need onto that mount once:

```bash
bp-warm Boltz2 ProteinMPNN        # builds per-tool micromamba envs + downloads weights
python my_pipeline.py             # subsequent runs reuse the warm env + cached weights
```

`Dockerfile.container` in the repo root builds a ready image (CUDA 12.4 +
micromamba + `pip install -e '.[colab]'`, with `build-essential`/`gcc` present
— several tools JIT-compile a CUDA/C helper at import and fail without it).

**Why `BIOPIPELINES_LOCAL_OUTPUT=0` matters:** with OTF on a non-Colab
scheduler the framework otherwise diverts output to the ephemeral `./outputs`
(cwd), silently overriding your configured `biopipelines_output`. On a container
that directory is lost at teardown. Setting it to `0` routes results to the
persistent mount.

## Reporting back

Check the completion markers before reporting anything as successful: they are empty files named `<NNN>_<ToolName>_COMPLETED` / `_FAILED` / `_WARNING`, written one level above each tool's output folder, and the per-tool log is `<Job>/Logs/<NNN>_<ToolName>.log`.

Save the structure (`.cif`/`.pdb`), the confidence/affinity JSON, and a summary
figure as artifacts. For a co-fold, report pTM, ipTM/ligand-ipTM, complex pLDDT,
and `affinity_probability_binary` (`affinity` defaults to `True`). Note that Boltz2 has no
covalent-mechanism knowledge: an unconstrained co-fold of a covalent ligand
finds a non-covalent pocket, not the reactive residue — use Boltz2's
`covalent_linkage` constraint when the mechanism is covalent.
