# `llm/` — using BioPipelines with an AI coding assistant

This folder contains everything you need to drive BioPipelines from an LLM
coding assistant (Claude Code, Cursor, Copilot Chat, etc.): five prompt files
that set up the assistant, a small ssh/scp wrapper that logs every cluster
command you run, and a template for cluster-specific resource notes.

## Quick start

1. **Pick a prompt** depending on what you want to do (see below) and load
   it into your assistant at the start of the session. A first message like

   > Read and follow `llm/pipelines.md`. <then your actual request>

   works for most assistants. For framework work, swap in `llm/development.md`.

2. **If you'll run on a cluster:** also do the one-time cluster setup —
   follow `cluster.md` (add an ssh alias, smoke-test with
   `ssh cluster echo ok`), then copy `resources.md.template` to
   `resources.md` and fill in your cluster's partitions, GPU types, and walltime policy by running the probe commands in `pipelines.md`. This gives the assistant honest defaults instead of generic guesses.

   **If you'll run on Google Colab:** no cluster setup needed — `pipelines.md`
   covers the in-notebook setup and `colab.md` covers the optional one-time Colab MCP server registration that lets the assistant execute cells itself. Skip `cluster.md` and `resources.md`.

   **If you'll run on CSCS Alps/Daint:** read `daint.md` as well — it is a SLURM
   cluster, but aarch64, whole-node billing and the CSCS Container Engine make its defaults different.

## Which prompt to use

Two of the five are session prompts — pick one of these to load first:

- **`pipelines.md`** — when you want to *use* the framework: design a
  pipeline for a specific biological problem and run it. Covers both execution modes (SLURM cluster and Google Colab) — the prompt asks you which one applies up front and adapts its defaults accordingly.
- **`development.md`** — when you want to *change* the framework itself:
  add a tool wrapper, fix a bug, refactor internals, update docs.

If your task crosses both (e.g. you need a pipeline but also hit a bug in an existing tool), handle them in two separate sessions. The two prompts give different defaults and pull in different reference docs; mixing them tends to produce muddled answers.

The other three are backend references that the two session prompts point at as needed — load them only for the backend you're on:

- **`cluster.md`** — ssh alias setup and the `log.sh`-wrapped cluster idioms.
- **`colab.md`** — Colab MCP server setup, Drive persistence, Colab gotchas.
- **`daint.md`** — CSCS Alps/Daint: access, storage, node packing, the Container
  Engine, and which tools are verified there.

## Logging cluster activity

(Cluster mode only — Colab notebooks log themselves via cell outputs.)

**All cluster commands must go through `llm/log.sh`.** This is not optional.
Every `ssh` and `scp` invocation the assistant makes has to be wrapped, so that
the command and its full output land in `logs/YYYY-MM-DD.log`.
The assistant should refuse to run a raw `ssh`/`scp` outside of `log.sh`.

```bash
# Wrong — leaves no trace.
ssh cluster 'squeue -u $(whoami)'

# Right — logged to logs/YYYY-MM-DD.log.
llm/log.sh ssh cluster 'squeue -u $(whoami)'
```

## File inventory

| File                    | Purpose                                                  | Tracked? |
| ----------------------- | -------------------------------------------------------- | -------- |
| `pipelines.md`          | Prompt for pipeline-author sessions                      | yes      |
| `development.md`        | Prompt for framework-developer sessions                  | yes      |
| `cluster.md`            | One-time ssh alias setup + `log.sh` usage                | yes      |
| `colab.md`              | Colab MCP server setup + Colab-specific gotchas          | yes      |
| `daint.md`              | CSCS Alps/Daint backend reference                        | yes      |
| `log.sh`                | Generic command logger (wraps ssh / scp / anything)      | yes      |
| `resources.md.template` | Schema for the user-local `resources.md`                 | yes      |
| `resources.md`          | Your cluster-specific resource notes                     | no       |
| `logs/`                 | Dated logs from `log.sh`                                 | no       |
