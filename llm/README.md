# `llm/` — using BioPipelines with an AI coding assistant

This folder holds the two session prompts and a template for cluster-specific resource notes. The backend references moved to `skills/biopipelines/references/` — they describe a cluster, not a prompt, so any host can read them there, while a file under `llm/` is only reachable by a session that was told to read `llm/`.

**Register `bp-mcp` first if you can.** With the MCP server in place the assistant calls tools — `bp_tools`, `bp_submit`, `bp_status`, `bp_logs`, `bp_table` — instead of composing `ssh` commands, and each run records what was done to it as a side effect. The prompts here still describe the manual path, and mark it as the fallback. Setup: `pip install -e ".[mcp]"` then `claude mcp add --scope user biopipelines -- bp-mcp`; details in `skills/biopipelines/references/mcp_server.md`.

## Quick start

1. **Pick a prompt** depending on what you want to do (see below) and load
   it into your assistant at the start of the session. A first message like

   > Read and follow `llm/pipelines.md`. <then your actual request>

   works for most assistants. For framework work, swap in `llm/development.md`.

2. **If you'll run on a cluster:** also do the one-time cluster setup —
   follow `skills/biopipelines/references/cluster_backend.md` (add an ssh alias, smoke-test with
   `ssh cluster echo ok`), then copy `resources.md.template` to
   `resources.md` and fill in your cluster's partitions, GPU types, and walltime policy by running the probe commands in `pipelines.md`. This gives the assistant honest defaults instead of generic guesses.

   **If you'll run on Google Colab:** no cluster setup needed — `pipelines.md` covers the in-notebook setup and `skills/biopipelines/references/colab_backend.md` covers the optional one-time Colab MCP server registration that lets the assistant execute cells itself. Skip the cluster reference and `resources.md`.

   **If you'll run on CSCS Alps/Daint:** read `skills/biopipelines/references/daint_backend.md` as well — it is a SLURM cluster, but aarch64, whole-node billing and the CSCS Container Engine make its defaults different.

## Which prompt to use

Two session prompts live here — pick one to load first:

- **`pipelines.md`** — when you want to *use* the framework: design a
  pipeline for a specific biological problem and run it. Covers both execution modes (SLURM cluster and Google Colab) — the prompt asks you which one applies up front and adapts its defaults accordingly.
- **`development.md`** — when you want to *change* the framework itself:
  add a tool wrapper, fix a bug, refactor internals, update docs.

If your task crosses both (e.g. you need a pipeline but also hit a bug in an existing tool), handle them in two separate sessions. The two prompts give different defaults and pull in different reference docs; mixing them tends to produce muddled answers.

The backend references the prompts point at are in `skills/biopipelines/references/` — load only the one you are on:

- **`cluster_backend.md`** — ssh alias setup and the cluster idioms, for a session without `bp-mcp`.
- **`colab_backend.md`** — Colab MCP server setup, Drive persistence, Colab gotchas.
- **`daint_backend.md`** — CSCS Alps/Daint: the variant's venv/EDF mechanism, then access, storage, node packing and which tools are verified there.
- **`container_backend.md`** — any single-node GPU box via the generic `container` variant.

## Recording what was done

A run's operational record lives **with the run**: `<Job>_NNN/_operations.jsonl`, next to its outputs and completion markers. It answers "what was done to this campaign" — submitted, resubmitted, cancelled, fetched — and it travels with the results when they are copied off the cluster. The same file is written on every backend, so Colab and container runs are covered identically.

This replaces `llm/log.sh`, which wrapped every `ssh`/`scp` into `logs/YYYY-MM-DD.log`. That shape logged the transport rather than the work, filed it by date rather than by run, covered only the cluster, and rested on the assistant remembering to use it. Anything already under `llm/logs/` is historical and untracked; nothing writes there now.

Where `bp-mcp` is registered, its tools reach the cluster themselves and write the record as they go — see `skills/biopipelines/references/mcp_server.md`.

## File inventory

| File                    | Purpose                                                  | Tracked? |
| ----------------------- | -------------------------------------------------------- | -------- |
| `pipelines.md`          | Prompt for pipeline-author sessions                      | yes      |
| `development.md`        | Prompt for framework-developer sessions                  | yes      |
| `resources.md.template` | Schema for the user-local `resources.md`                 | yes      |
| `resources.md`          | Your cluster-specific resource notes                     | no       |
| `logs/`                 | Historical logs from the retired `log.sh`                | no       |
