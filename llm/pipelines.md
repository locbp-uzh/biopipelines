# Working agreement — writing BioPipelines for a specific problem

This file is the working agreement for **pipeline-author sessions** on this
repo: designing and running a BioPipelines pipeline for a concrete biological
problem. Do not modify framework code in this mode — for that, see
`llm/development.md`.

**How to use this file.** Absorb the conventions below and apply them as
relevant — do not recite them back to the user, do not announce that you've
read them, do not run the interview as a checklist before the user has
described a task. Wait for the actual request, then act.

## Read at session start

Before responding to the user's first task, read these files end-to-end so
you have the framework's contract loaded — the `Pipeline`/`Resources` API,
typed streams, the cluster-vs-Colab differences, and which tools exist:

- `docs/user_manual.md`
- `docs/tool_reference.md`
- `docs/tool/*.md`
- `README.md` — the tool table records, per tool, the hardware (CPU/GPU), the environments it has been tested on (HPC/Colab), and the upstream references (repo + paper). Read all of this off the badge markup on each row, never the base64 `<img src=...>` blobs (those carry no information): the `alt="..."` attributes give the status (`alt="CPU"`, `alt="GPU"`, `alt="HPC ok"`, `alt="Colab ok"`), and the `<a href="...">` wrapping a badge whose `alt` is `repo` or `paper` gives the reference URL. A BP-native tool shows a `BP` badge and has no repo/paper link. The HPC/Colab badges are presence-means-supported: a tool runs only in the modes whose badge it carries, so a missing badge means that mode is unsupported (not merely unverified). Only MMseqs2/MMseqs2Server is currently HPC-only — don't propose it in a Colab pipeline.

Do not announce that you've read them, do not summarise
them, do not list what you found. Just have the context loaded so your
suggestions and questions are grounded in actual API shapes and tool names
instead of inferred guesses.

## Execution mode (cluster vs Colab)

BioPipelines runs in two modes that share the API but differ in
infrastructure:

- **Cluster** — `.py` scripts submitted via `./submit` to SLURM, conda envs
  via `mamba`/`conda`/`micromamba`, outputs persist on shared storage.
  Read `cluster.md` for more information.
- **Colab** — `.ipynb` cells executed inline, `micromamba` envs installed
  per session, runtime ephemeral (~12 h cap, anything not saved to Drive or
  downloaded is lost).

The auto-loaded config differs (`config.cluster.yaml` vs `config.colab.yaml`).
If the user hasn't made it clear which mode applies, infer from context
(notebook file mentioned → Colab; ssh/`./submit`/SLURM mentioned → cluster);
ask only if it remains genuinely ambiguous.

## Common concerns (both modes)

### Environment assumption
This session runs on a **local machine, not the execution environment**.
Project folders, datasets, conda/micromamba envs, and GPUs are not
accessible. GPU-bound code cannot be executed locally — verify it via review,
and if available, CI.

### Clarifying questions (only when the request leaves choices open)
Run this protocol only when the user's request leaves real choices unresolved
(a new pipeline, an open-ended exploration, an unfamiliar problem). Skip it
entirely for directive requests like "rerun the FRET pipeline with `n=200`"
— just do the thing. Skip categories the user has already pinned, and group
remaining questions into one message. The goal is to avoid hallucinating
defaults, not to interrogate the user.

**Do not narrate procedural steps.** No "Let me check X first", no "I'll
look at Y before answering". Reading docs, searching the codebase, and
finding similar example pipelines are internal work — the user shouldn't
see them. The only user-visible output during the interview is the
substantive question(s).

**Question format.** Use **closed-option questions** over open-ended
prose.

**First check whether `AskUserQuestion` is in your tool list for this
session.** If it isn't, you don't have it — go straight to the markdown
fallback below; do not pretend to invoke it. If it is available, **you
must use it** for any question with 2–4 discrete options. Markdown
numbered lists are the fallback for hosts without the tool, not a
stylistic preference when the tool is present.

`AskUserQuestion` constraints: 1–4 questions per call, each with 2–4
options; option `label` is 1–5 words (longer text goes in `description`);
recommended option first with "(Recommended)" suffix; "Other" is
auto-added.

Markdown fallback (use only when `AskUserQuestion` is unavailable):

> **Scope?** (default: a)
> 1. **(a) minimal** — wrapper + `pipe_*` script + install + one `tool_reference.md` entry
> 2. **(b) standard** — minimal + tests entry + `user_manual.md` mention
> 3. **(c) full** — standard + example pipeline

Either way, the question text and option set should be identical across
renderers — the format is a presentation choice, not a content choice.

**Do not use `I` and `you` when asking questions** Refer to yourself as
"The coding agent" and to the user as "The user".

#### 1. Problem framing
- What is the biological question or objective?
- What is the input (sequences, structures, ligands, tabular data) and where
  is it found? (Cluster: local files require scp. Colab: local files need to
  be uploaded to the runtime or mounted from Drive.)
- What does success look like — a ranked list, a binder, a fitted curve?

#### 2. Scale
- How many inputs (sequences / structures / conditions)?
- One-shot run, parameter sweep, or many parallel jobs? Or one-shot run,
  automatic evaluation, then scale-up?
- *Cluster only:* for sweeps, should this become N parallel SLURM jobs (see
  `example_pipelines/multiple_submission.py`) or one job with internal
  looping?
- *Colab only:* sweeps run sequentially on the single assigned GPU. Discuss
  whether the total runtime fits within the 12 h kernel limit, and split the
  notebook across sessions if not.

#### 3. Outputs
- Which artifacts must be kept long-term vs left in the run directory?
- Which plots / summary tables should the pipeline produce?
- *Cluster:* which artifacts must be pulled back locally vs left on the
  cluster? Pull only a representative subset.
- *Colab:* which artifacts must be saved to Drive or downloaded before the
  runtime disconnects? Outputs in `/content/BioPipelines/` evaporate when
  the kernel resets.

#### 4. Iteration plan
- Run once with **minimal parameters** first to validate the approach, then
  scale up. Confirm the user agrees before launching a large job.

### Code rules (both modes)
- Don't pass tool parameters that aren't relevant to the problem and don't
  restate defaults (e.g. `num_recycles` on AlphaFold).
- Recycle MSAs across runs whenever the inputs allow it.

### Reporting back
- Summarise what was run, the resource choices, and why.
- Surface any failures from success markers (`SUCCESS` / `FAILURE` /
  `WARNING`) rather than claiming success blind, and read the relevant tool
  `_log` files.
- *Cluster:* pull only a subset of artifacts locally — enough to inspect,
  not the full output tree.
- *Colab:* remind the user to save what they want to keep before the runtime
  ends.

---

## Cluster mode

### Resources

- **Site-specific guidance lives in `llm/resources.md`.** If that file
  doesn't exist yet, populate it first (see *Probing the cluster* below).
  Treat its values as the source of truth for time / GPU / partition
  defaults; the rules below are fallbacks for when no probe has been done.
- **Time**: default to `Resources(time="24:00:00")` and matching
  `#SBATCH --time=24:00:00`, then scale up in 24 h multiples for longer
  runs. Why: many SLURM sites bill or schedule by *requested* wall-time, so
  under-requesting risks job kill while over-requesting can cost queue
  priority — but on sites whose backfill window is exactly 24 h
  (e.g. UZH S3IT), asking for 6 h vs 24 h gets you the same queue position
  with more failure margin, which makes 24 h the safer default.
  `llm/resources.md` overrides this if the local site's policy differs
  (shorter partition limits, walltime-based charging, etc.).
- **GPU**: use `"any"` for small / short jobs. Request a specific class
  (`a100`, `h100`) only when memory or throughput demands it; explain the
  tradeoff to the user before deciding. Available GPU classes vary by site
  — check `llm/resources.md`.
- **CPU-only stages**: identify which tools don't need a GPU and avoid
  reserving one for them.

### Probing the cluster (one-time)

If `llm/resources.md` is missing or older than ~3 months, regenerate it
before suggesting concrete resource values. On the cluster, run:

```bash
llm/log.sh ssh cluster 'sinfo -o "%P %l %G %D %t" | sort -u'   # partitions, time limits, GPUs
llm/log.sh ssh cluster 'sinfo -o "%G" | sort -u'                # available GPU types
llm/log.sh ssh cluster 'scontrol show config | grep -E "SchedulerType|DefaultTime|MaxArraySize"'
llm/log.sh ssh cluster 'sacctmgr show qos format=Name,MaxWall,Priority 2>/dev/null | head'
```

Record the answers in `llm/resources.md` (use `llm/resources.md.template` as
the starting structure). The file is gitignored so each user keeps their own.
Ask the user about the default time for a job, in case billing is based on
time requested.

### Writing the pipeline

- Produce a `.py` script. Place it in `my_pipelines/`. The folder is
  gitignored, so `git push` will not propagate it — scp it to the cluster
  instead.
- Keep input data on the cluster. Do **not** stage CSVs / FASTA / PDB inputs
  inside the local working copy (also gitignored, but more importantly:
  large inputs don't belong in a code repo); scp them directly to the
  cluster.

### Running the pipeline

All cluster interaction is over plain `ssh` and `scp`, and **every call must
be wrapped in `llm/log.sh`** — no raw `ssh`/`scp` invocations. The wrapper
appends each command and its full output to `llm/logs/YYYY-MM-DD.log`, which
is the audit trail for the session. If the user runs a raw cluster command,
re-issue it through `log.sh` rather than relying on the unlogged result.

Typical idioms (assume an ssh alias `cluster` and a remote repo at
`~/biopipelines`):

| Step                              | Command                                                                                |
| --------------------------------- | -------------------------------------------------------------------------------------- |
| Sync remote repo to a branch      | `llm/log.sh ssh cluster "cd ~/biopipelines && git fetch && git reset --hard origin/<branch>"` |
| Copy a personal pipeline / inputs | `llm/log.sh scp my_pipelines/foo.py cluster:~/biopipelines/my_pipelines/`              |
| Submit                            | `llm/log.sh ssh cluster "cd ~/biopipelines && ./submit my_pipelines/foo.py"`           |
| Resume after cancel/fail          | `llm/log.sh ssh cluster "cd ~/biopipelines && ./resubmit <slurm.sh>"`                  |
| Watch the queue                   | `llm/log.sh ssh cluster 'squeue -u $(whoami)'`                                         |
| Inspect logs                      | `llm/log.sh ssh cluster "tail -n 100 <RunTime>/slurm.out"`                             |
| Pull a result artifact            | `llm/log.sh scp cluster:<RunTime>/slurm.out ./`                                        |
| Cancel a job                      | `llm/log.sh ssh cluster "scancel <jobid>"`                                             |

Outputs land at the path configured on the cluster
(typically `/shares/<group>/<user>/BioPipelines/<Project>/<Job>_NNN/`).

---

## Colab mode

### Resources

- **GPU**: Colab assigns one GPU per runtime. The user picks the *runtime
  type* (T4 free tier; A100/L4/V100 on paid tiers) from the Colab UI — your
  pipeline cannot request a class. Don't write `Resources(gpu=...)`-style
  hints into the notebook; they have no effect.
- **Time**: no SLURM walltime to set. The hard limit is the Colab kernel
  budget (~12 h on paid; shorter on free, with idle disconnects). Plan
  pipelines to either fit in one session or checkpoint to Drive between
  sessions.
- **CPU-only stages**: still worth identifying — they're faster on a
  cheaper runtime, and switching runtime type costs the user time
  reinstalling envs.

### Writing the pipeline

- Produce a Jupyter notebook (`.ipynb`), not a `.py` script. Example notebooks
  live in `examples/notebooks`, and personal development is in `my_pipelines/`
- First cell: clone the repo, install BioPipelines, and call `.install()` on
  the tools the pipeline needs. See the "Google Colab" section in
  `docs/user_manual.md` for the canonical setup snippet.
- The pipeline runs **inline as cells execute** — there is no `./submit`,
  no SLURM, no two-phase configure-then-execute. Treat the pipeline as
  ordinary notebook code that happens to use the BioPipelines API.
- Inputs: small files can be uploaded directly to the runtime
  (`/content/`); larger or persistent inputs should be mounted from Drive.
- The auto-detected config is `config.colab.yaml` — confirm it's selected
  by checking `Pipeline(...).config` in an early cell if anything looks off.

### Running the pipeline

There is nothing to ssh into. The "running" step is just executing the
notebook cells in order. `llm/log.sh` does not apply in Colab mode — the
notebook itself is the audit trail (cell outputs are saved with the
`.ipynb`).

Things to remind the user about during execution:

| Concern                | What to do                                                                  |
| ---------------------- | --------------------------------------------------------------------------- |
| Save the notebook      | File → Save in Drive, periodically. The autosave is *to Colab*, not Drive. |
| Save outputs           | All tools implement `.download(...)`.    |
| Long-running tools     | Run the long cell, then check on it — Colab will idle-disconnect a tab.    |
| Reinstalling each session | Tools install once per kernel; budget ~5–15 min of setup at the start.   |

Outputs land in `/content/BioPipelines/<Project>/<Job>_NNN/` by default,
which is **ephemeral**. Anything the user wants to keep must be copied to
Drive or downloaded.
