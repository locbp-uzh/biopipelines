---
name: biopipelines
description: >-
  Design and run computational protein and ligand workflows on a GPU: binder and enzyme design,
  de novo backbone generation, inverse folding and sequence redesign, structure prediction,
  protein-ligand docking, compound-library and covalent screening, saturation mutagenesis,
  stability and solubility engineering, codon optimization. BioPipelines (locbp-uzh, CSBJ 2026)
  puts 86 tools behind one declarative Pipeline API - RFdiffusion 1/2/3, BoltzGen, ProteinMPNN,
  LigandMPNN, LASErMPNN, Frame2Seq, ThermoMPNN, AlphaFold, Boltz2, ESMFold2, DiffDock, GNINA,
  Vina, NeuralPLexer, OpenMM, FPocket, P2Rank, PLIP, ProLIF, PoseBusters, Prodigy, RDKit,
  ADMET-AI - and tracks IDs through typed data streams so a multi-stage campaign stays
  traceable end to end. Load when the user wants to design, predict, dock, screen, redesign or
  analyze proteins or ligands, names one of these tools, or has a BioPipelines pipeline to
  write, run or debug.
---

# BioPipelines

[BioPipelines](https://github.com/locbp-uzh/biopipelines) (Quargnali & Rivera-Fuentes, LOC-BP UZH; CSBJ [10.34133/csbj.0129](https://spj.science.org/doi/10.34133/csbj.0129)) is a Python framework that puts 86 protein/ligand-modeling tools behind one declarative `Pipeline` API. You describe a workflow as a chain of tools; the framework generates and runs the per-tool scripts, tracks IDs through typed data streams, and materializes outputs.

This skill ships *inside the repo*, so adding it to any agent is one step and it tracks the framework as it evolves. It covers **using** the framework, and running it on **any single-node GPU backend** (Modal/RunPod, a Docker GPU box, or an interactive Slurm GPU shell) via the generic `container` config variant.

## First: is `bp-mcp` registered?

Look at your own tool list. If it holds `bp_tools`, `bp_setup`, `bp_runs`, `bp_status`, `bp_logs`, `bp_lineage`, `bp_table`, `bp_submit`, `bp_resubmit`, `bp_cancel`, `bp_visualize`, `bp_fetch`, `bp_project`, `bp_provenance` and `bp_reproduce`, **those tools are the interface to BioPipelines** and the file-reading and shell recipes in this skill are what a session without them falls back on.

| What you want | Call | Instead of |
| --- | --- | --- |
| Find a tool | `bp_tools()`, `bp_tools(tags=[...])`, `bp_tools(outputs="rmsd")` | reading `references/tool_index.md` |
| A tool's parameters | `bp_tools(name="LigandMPNN")` | opening `docs/tool/*.md` |
| Run a pipeline | `bp_submit(script=...)` | `ssh s3it "./submit ..."` |
| Status, logs, tables, lineage | `bp_status`, `bp_logs`, `bp_table`, `bp_lineage` | `ssh` + `ls` + `tail` + `scp` |
| Resume, stop, or show a run | `bp_resubmit`, `bp_cancel`, `bp_visualize` | `./resubmit`, `scancel`, `bp-visualize` over ssh |
| Project documents | `bp_project(project_dir=..., action=...)` | writing PROJECT.md by hand |

Do not run a script in this repo, or a shell command, to do what one of these does. `skills/biopipelines/build_tool_index.py` in particular is a **maintainer** script that regenerates the catalog file; it is not how you look a tool up. If no `bp_*` tool is present, follow the sections below, and see `references/mcp_server.md` to set the server up.

## Finding tools without the MCP server — read the index, not the prose docs

**`references/tool_index.md` is your catalog.** 88 callable entries, one line each: name, version, hardware (CPU/GPU), the platforms it is verified on, what it does, and a pointer to its full entry. It costs ~4k tokens; the prose docs it points into cost ~78k, and the README's badge table another ~32k.

The loop is: **read the index → pick your tools → read only those sections.** Each entry ends with `docs/tool/<file>.md#<anchor>` — read that section for the real signature before you write the stage.

Each index line ends with the tool's tags — a closed 38-term vocabulary in four facets, defined in `docs/tool_tags.md`. Grep them: `covalent`, `binder-design` and `motif-scaffolding` are the sort of capability the one-line summary never mentions.

Two rules, because both failure modes are silent:

- **Never read a whole `docs/tool/*.md` file to find a tool.** `analysis.md` alone is 24k tokens for 33 tools, and a campaign typically uses three of them.
- **Never guess a tool name or a parameter.** An unknown keyword is swallowed by `**kwargs` rather than rejected, so a typo becomes a silently ignored argument, not an error.

The index lists 88 entries against 86 registered tools: `SolubleMPNN` and `LoadMultiple` are public classes that share a parent's identity (`ProteinMPNN` and `Load`), so they are callable but carry no separate registration.

## The framework contract

`llm/` is the authoritative, author-maintained agent contract — it overrides anything paraphrased here. Read the file that matches what you are doing; do not read them all up front.

- `llm/pipelines.md` — authoring a `Pipeline`: the tool/data-stream model, ID tracking, combinatorics. Read this before writing a multi-stage workflow.
- `llm/development.md` — only when changing framework code, not when using it.
- `references/cluster_backend.md`, `references/colab_backend.md`, `references/daint_backend.md`, `references/container_backend.md` (in this skill) — one backend reference each. Read the one you are on. They live here rather than in `llm/` so any host can use them, not only a session that was told to read `llm/`.
- `docs/user_manual.md` — the long-form manual. Consult sections as needed; it is ~17k tokens end to end.
- `references/mcp_server.md` (in this skill) — installing and registering `bp-mcp`, which serves `bp_tools`.

## The API in one screen

```python
from biopipelines import Pipeline, Resources, Sequence, Ligand, UniProt, Boltz2

with Pipeline("Project", "job_name", description="..."):
    Resources(gpu="A100", memory="64GB", time="6:00:00", cpus=8)
    prot = Sequence("MSEQ...", type="protein")           # or UniProt("Q15436") to fetch it
    lig  = Ligand(smiles="C[N+]1=C(...)...")
    Boltz2(proteins=prot, ligands=lig, output_format="mmcif")
```

Every tool takes its inputs **by keyword**. `Boltz2`'s first positional parameter is `config` (a raw YAML string), so `Boltz2(prot, lig, ...)` binds the protein to `config` and the ligand to `proteins` — always write `proteins=` / `ligands=`.

Each `Resources(...)` call **opens a new batch** (one scheduler job), inheriting whatever it does not specify from the previous batch. Call it once per resource profile — a CPU-only stage after a GPU stage gets its own `Resources(gpu="none", ...)` rather than inheriting the GPU and walltime of the heaviest stage. Inside a plain `Parallel()` block every sibling iteration **must** call `Resources()` to open its own batch; inside `Parallel(pack=N)` it is called **once** to describe the whole node allocation and `Run(...)` delimits each task.

`Tool.install()` is called **inside** the `with Pipeline(...)` block (a bare `Boltz2.install()` at module scope is a silent no-op) — or just use `bp-warm`, which wraps it for you.

## Running on any single-node GPU host (the `container` backend)

The repo ships a generic single-node variant so you do not hand-author a config per provider. Set three env vars and point one config line at your persistent mount:

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

`Dockerfile.container` in the repo root builds a ready image (CUDA 12.4 + micromamba + `pip install -e '.[colab]'`, with `build-essential`/`gcc` present — several tools JIT-compile a CUDA/C helper at import and fail without it).

**Why `BIOPIPELINES_LOCAL_OUTPUT=0` matters:** with OTF on a non-Colab scheduler the framework otherwise diverts output to the ephemeral `./outputs` (cwd), silently overriding your configured `biopipelines_output`. On a container that directory is lost at teardown. Setting it to `0` routes results to the persistent mount.

## Reporting back

Check the completion markers before reporting anything as successful: they are empty files named `<NNN>_<ToolName>_COMPLETED` / `_FAILED` / `_WARNING`, written one level above each tool's output folder, and the per-tool log is `<Job>/Logs/<NNN>_<ToolName>.log`.

Save the structure (`.cif`/`.pdb`), the confidence/affinity JSON, and a summary figure as artifacts. For a co-fold, report pTM, ipTM/ligand-ipTM, complex pLDDT, and `affinity_probability_binary` (`affinity` defaults to `True`). Note that Boltz2 has no covalent-mechanism knowledge: an unconstrained co-fold of a covalent ligand finds a non-covalent pocket, not the reactive residue — use Boltz2's `covalent_linkage` constraint when the mechanism is covalent.
