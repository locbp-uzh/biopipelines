# The generic single-node GPU backend (`container` variant)

BioPipelines ships two native backends: SLURM (`cluster`) and Google Colab
(`colab`). The `container` variant generalizes the Colab (single-node,
micromamba, no scheduler) model to *any* single GPU node with a persistent
mount — a managed-container provider (Modal, RunPod, Vast.ai), a plain Docker
GPU host, or an interactive Slurm GPU shell.

## The three env vars

| Variable | Value | Effect |
|---|---|---|
| `BIOPIPELINES_CONFIG_VARIANT` | `container` | loads `config.container.yaml` (checked before the colab autodetect) |
| `BIOPIPELINES_OTF` | `1` | run each tool inline as it is added; no job submission |
| `BIOPIPELINES_LOCAL_OUTPUT` | `0` | write to the configured `biopipelines_output`, NOT the ephemeral cwd |

## One line to edit

`config.container.yaml` -> `folders.base.root:` (default `/workspace`). Point it
at the persistent mount. Everything derives from `<root>`:

- `home`, `data`, `scratch` under `<root>`
- micromamba env root under `<root>` (set `MAMBA_ROOT_PREFIX=<root>/micromamba`)
- weight caches (`BoltzCache`, `ColabFoldDatabases`, ...) under `<root>/cache`
- `biopipelines_output` = `<root>/outputs`

Because envs and weights live on the persistent mount, the FIRST job pays the
install/download cost (`bp-warm`) and every later job reuses it.

## Warm-up

`Tool.install()` only registers an install step inside an active `Pipeline`
context. Use the `bp-warm` CLI, which wraps the install in a context for you:

```bash
bp-warm Boltz2                     # one tool
bp-warm Boltz2 ProteinMPNN RFdiffusion
bp-warm --gpu A100 --time 2:00:00 Boltz2
bp-warm --force-reinstall Boltz2   # rebuild even if present
```

Installs are idempotent; re-running skips tools already present.

## Compiler requirement

Build the image from a CUDA base and install `build-essential gcc g++`, then
`export CC=gcc CXX=g++`. Several tools (Boltz-2 through Triton, cuda ops)
JIT-compile a CUDA/C helper at import; a runtime-only base image fails at import
with a compiler-not-found error. See `Dockerfile.container`.
