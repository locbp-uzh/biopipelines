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

**First, look at your own tool list.** If it contains tools named `bp_tools`, `bp_setup`, `bp_status` and the rest, the `bp-mcp` server is registered and those tools — not shell commands, not scripts in this repo — are how you interact with BioPipelines: finding a tool, submitting a run, reading its status, its logs, its tables. The sections below describe the same operations for a session without them. Where the two differ, the tool wins.

Before responding to the user's first task, load the framework's contract — the
`Pipeline`/`Resources` API, typed streams, and the cluster-vs-Colab differences:

- `docs/user_manual.md`
- `skills/biopipelines/references/tool_index.md` — the tool catalog, one line per tool
  (~4k tokens). **This replaces reading `docs/tool_reference.md` and `docs/tool/*.md` up
  front**, which together run to ~78k tokens. Read the index, pick the tools the task needs,
  then read only those sections: every entry ends with a `docs/tool/<file>.md#<anchor>`
  pointer to its authoritative entry, and that is where the real signature lives. Read a
  section before you write the stage that uses it — an unknown keyword is swallowed by
  `**kwargs` rather than rejected, so a guessed parameter is silently ignored instead of
  erroring. Never read a whole `docs/tool/*.md` file to find a tool: `analysis.md` alone is
  ~24k tokens for 33 tools, and a campaign typically uses three of them.

The index also carries, per tool, the hardware it needs (CPU/GPU) and the platforms it has been verified on, lifted from the README's badge table — so you never need to open `README.md`, which is ~32k tokens and written for people, not for you. Each line ends with the tool's tags, so a capability the one-line summary omits — `covalent`, `binder-design`, `symmetry` — is still greppable in the index alone.

Do not announce that you've read them, do not summarise
them, do not list what you found. Just have the context loaded so your
suggestions and questions are grounded in actual API shapes and tool names
instead of inferred guesses.

**With `bp_tools` in your tool list, none of the reading above is how you find a tool.** Call `bp_tools()` for the catalog and `bp_tools(name="Boltz2")` for one tool's section — same content, same files, no path to get wrong. It also answers the two questions the index file cannot: `bp_tools(outputs="rmsd")` lists every tool whose streams, tables or columns match that fragment, read from the tools' own source, and `bp_tools(tags=[...], exclude=[...])` filters by intent over the vocabulary in `docs/tool_tags.md`. Opening `docs/tool/*.md` by hand when the tool is available costs tokens and finds nothing extra.

## Execution mode (cluster vs Colab)

BioPipelines runs in two modes that share the API but differ in
infrastructure:

- **Cluster** — `.py` scripts submitted via `./submit` to SLURM, conda envs
  via `mamba`/`conda`/`micromamba`, outputs persist on shared storage.
  Read `skills/biopipelines/references/cluster_backend.md` for more information.
- **Colab** — `.ipynb` cells executed inline, `micromamba` envs installed
  per session, runtime ephemeral (~12 h cap, anything not saved to Drive or
  downloaded is lost). Read `skills/biopipelines/references/colab_backend.md` for more information, including whether this session can drive the runtime itself.

A third variant, `daint` (CSCS Alps, aarch64), is also a SLURM cluster but differs enough to have its own file — read `skills/biopipelines/references/daint_backend.md` if the user is on Daint.

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

**Show results, do not only describe them.** When a step has produced anything renderable — structures, images, plots, a scores table — do not answer with prose alone. Render it and open it. `bp-visualize <tool folder>` writes a self-contained HTML page for that one step; it runs on a login node, needs no scheduler, and works while the rest of the job is still queued, so a mid-run step can be shown the moment it finishes.

- **Mid-run, per step:** `bp-visualize` on the step that just completed. This is the one to reach for while a job is still going, and for "show me the best N" questions, which it answers directly: `--descending <table>.<column>` orders by any column of any of the step's tables and `--max-items N` caps it, overriding the structure viewer's default 5-item sample.
- **Whole run, at the end:** `pipeline.html` (below). `bp-visualize` does not replace it — one is a step, the other is the run.
- **Then actually open it. Do not just report a path.** Opening it is part of reporting the result, not an optional extra; a path in the transcript is something the user has to go and act on, which is the work you were asked to do.
  - Windows: `Start-Process -FilePath <browser.exe> -ArgumentList '"<page>"'` from PowerShell. Do **not** use `cmd.exe /c start` — it blocks while cold-starting a browser that is not already running (measured: one call sat 120 s and was backgrounded), and it exits 0 as soon as the shell accepts the file, so its exit code says nothing about whether a window appeared. Find the real browser from the `ProgId` under `HKCU:\SOFTWARE\Microsoft\Windows\Shell\Associations\UrlAssociations\http\UserChoice`; the default handler is not necessarily the browser that happens to be running.
  - macOS: `open <page>`. Linux: `xdg-open <page>`.
- **Verify the window, not the launcher's exit code.** A launcher returning 0 is not evidence the page opened. On Windows confirm with `Get-Process <browser> | Where-Object MainWindowTitle | Select-Object MainWindowTitle` — the title is the page's own `<title>`, so it proves both that a window exists and that it holds the right page. Only then tell the user it is open.
- **Put it somewhere stable before you open it.** A session scratchpad or `/tmp` is wiped and is not findable later. Copy pulled pages into `outputs/_views/` in the repo (gitignored), keep the step's own filename, and tell the user that folder. On the cluster the page already lives at `<step>/_extras/<step>_view.html`, which is its permanent home.
- **If you cannot render or cannot open it — no renderable output, the page fails to build, or a backend with no browser — say so and give the exact command the user can run.** Silently falling back to a prose summary is the failure mode this exists to remove.

The sort key must be qualified as `<table>.<column>` (`--descending confidence.plddt`, not `--descending plddt`); a bare column name is refused, naming the qualified spellings that would have worked. Look the columns up in `docs/tool_reference.md` or in the step's own `.expected_outputs.json`.

- **Always hand over the run page after any new result.** Every run writes `<Job>/RunTime/pipeline.html`, a self-contained record of the pipeline: one lane per batch with its resources, one card per step carrying its paths, tool version, environment and completion marker, arrows for the dataflow that was actually recovered, a completed-against-failed donut, and the machine the run used. It embeds the rendered outputs — table values and, when the 3D viewers are kept, the coordinates themselves — so it is the artifact to show, not a summary of it.
  - After a run finishes, refresh it: `regenerate_pipeline_page("<Job>/RunTime")`. `save()` writes the page at configuration time, when nothing has run and every node reads "pending", so an unrefreshed page is a plan rather than a result.
  - The 3D structure and grid viewers work offline: the py3Dmol library is vendored in the repo and inlined once per page, so a normal page is self-contained *with* working viewers (~695 KB against 28 KB without). `allow_external=True` is only needed if that vendored copy is absent.
  - *Cluster:* copy just that one file back. It needs no sibling files and no network, which is the whole point of it being self-contained — do not try to serve the output tree or link into it over SSH.
  - Tell the user the path and that it opens in a browser. Do not paraphrase the page in prose instead of providing it; the point is that they can look.
- Summarise what was run, the resource choices, and why.
- Surface any failures from the completion markers rather than claiming success blind. The markers are empty files written **one level above** each tool's output folder, named `<NNN>_<ToolName>_COMPLETED`, `<NNN>_<ToolName>_FAILED` or `<NNN>_<ToolName>_WARNING` (`<NNN>` is the tool's zero-padded step number, e.g. `003_Boltz2_COMPLETED`; see `pipe_scripts/pipe_check_completion.py`). So glob `*_FAILED` and `*_WARNING` in the job folder — there are no files named `SUCCESS`, `FAILURE` or `WARNING`.
- Then read the per-tool log. There are two copies of the same stream: `<Job>/Logs/<NNN>_<ToolName>.log` collects them job-wide, and `<tool output folder>/_log` sits next to the tool's own outputs (`pipeline.py:507` on the fly, and the `| tee` at `:764`). Reach for `_log` when you already have a tool folder open, and `Logs/` when you are scanning a whole job. There are no files named `SUCCESS`, `FAILURE` or `WARNING`.
- *Cluster:* pull only a subset of artifacts locally — enough to inspect,
  not the full output tree.
- *Colab:* remind the user to save what they want to keep before the runtime
  ends.

### Output recycling
- If a stage of the pipeline fails and has to be repeated, prefer loading
  results from upstream tools instead of running them again. This can be done
  with the *Load* tool: `Load("/path/to/tool_output_folder")`
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
ssh s3it 'sinfo -o "%P %l %G %D %t" | sort -u'   # partitions, time limits, GPUs
ssh s3it 'sinfo -o "%G" | sort -u'                # available GPU types
ssh s3it 'scontrol show config | grep -E "SchedulerType|DefaultTime|MaxArraySize"'
ssh s3it 'sacctmgr show qos format=Name,MaxWall,Priority 2>/dev/null | head'
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

**If `bp_submit`, `bp_status` or `bp_logs` appear in your tool list, they are how you reach the cluster. Use them.** Do not open a shell to do what one of them does: an `ssh` command that submits a run bypasses the operations record the tool would have written, guesses which `.out` file exists, and reads markers the tool already interprets. Shell out only for something no tool covers, and say which tool you checked first.

| What you want | Call | Not |
| --- | --- | --- |
| Run a pipeline on the cluster | `bp_submit(script=...)` | `ssh s3it "./submit ..."` |
| Did it work / where is it | `bp_status(job=...)` | `ssh s3it "ls ... && squeue"` |
| Why did a step fail | `bp_logs(job=..., step="007_Boltz2")` | `ssh s3it "tail Logs/..."` |
| Read a result table | `bp_table(job=..., step=..., table=...)` | `scp` then open the CSV |
| Which ids were dropped, and where | `bp_lineage(job=...)` | reading `missing.csv` by hand |
| A page a supervisor or referee can read | `bp_lineage(job=..., page=True)` | pasting numbers into a document |
| What produced this run | `bp_provenance(job=...)` | guessing from the commit date |
| How do two runs differ | `bp_provenance(job=..., against=...)` | diffing two output folders |
| Run this campaign again, elsewhere | `bp_reproduce(job=..., target=...)` | re-authoring the pipeline |
| Show me a finished step | `bp_visualize(job=..., step=...)` | describing the numbers in prose |
| Resume a run that failed partway | `bp_resubmit(job=...)` | resubmitting the whole thing |
| Stop queued jobs | `bp_cancel(job=...)`, then `confirm=True` | `ssh s3it "scancel ..."` |
| What runs exist | `bp_runs()` | `ssh s3it "ls <output root>"` |
| Bring a file back | `bp_fetch(remote_path=...)` | `scp s3it:...` |
| Which tool does X | `bp_tools(tags=[...])`, `bp_tools(outputs="rmsd")`, `bp_tools(name=...)` | grepping `docs/tool/` |
| Connect the first time | `bp_setup(host=..., repo=...)` | probing by hand |

The tools take `host=` and default to the saved one, so a call needs no ssh alias in it. Each run still records what was done to it in `<Job>_NNN/_operations.jsonl`, next to its outputs, so the record survives being copied off the cluster; the MCP tools write that record as they go, and a hand-run `ssh` does not.

#### Without the MCP server

Only when those tools are **not** in your tool list. Then all cluster interaction is plain `ssh` and `scp`, with the idioms below (assume an ssh alias named after the site — `s3it` here — and a remote repo at `~/biopipelines`):

| Step                              | Command                                                                                |
| --------------------------------- | -------------------------------------------------------------------------------------- |
| Sync remote repo to a branch      | `ssh s3it "cd ~/biopipelines && git fetch && git reset --hard origin/<branch>"` |
| Copy a personal pipeline / inputs | `scp my_pipelines/foo.py s3it:~/biopipelines/my_pipelines/`              |
| Submit                            | `ssh s3it "cd ~/biopipelines && ./submit my_pipelines/foo.py"`           |
| Resume after cancel/fail          | `ssh s3it "cd ~/biopipelines && ./resubmit <RunTime>/<job script>"`      |
| Watch the queue                   | `ssh s3it 'squeue -u $(whoami)'`                                         |
| List the scheduler output files   | `ssh s3it "ls <RunTime>/*.out"`                                          |
| Inspect a scheduler log           | `ssh s3it "tail -n 100 <RunTime>/job_batch1.out"`                        |
| Inspect a tool log                | `ssh s3it "tail -n 100 <Job>/Logs/<NNN>_<ToolName>.log"`                 |
| Cancel a job                      | `ssh s3it "scancel <jobid>"`                                             |
| Refresh the run page after a run  | `ssh s3it "cd ~/biopipelines && python -c 'from biopipelines.pipeline import regenerate_pipeline_page as r; r(\"<RunTime>\")'"` |
| Fetch the run page                | `scp s3it:<RunTime>/pipeline.html ./`                                    |
| Render one step (top 5 by a score) | `ssh s3it "cd ~/biopipelines && bp-visualize <Job>/<NNN>_<Tool> --descending <table>.<column> --max-items 5"` |
| Fetch a step's page               | `scp s3it:<Job>/<NNN>_<Tool>/_extras/<NNN>_<Tool>_view.html ./`           |

**Which `.out` file to tail.** Every `Resources()` call in the pipeline opens a new batch, and each batch is submitted as its own job. `submit` writes one scheduler output per batch: `<RunTime>/job_batch<N>.out` (from `<RunTime>/slurm_batch<N>.sh`). Only a pipeline that ended up with a *single* batch gets `<RunTime>/slurm.out` (from `<RunTime>/slurm.sh`). Since most pipelines call `Resources()` more than once, `ls <RunTime>/*.out` first rather than assuming `slurm.out` exists — and pass the matching `slurm_batch<N>.sh` to `resubmit` when resuming.

The scheduler `.out` files carry the submission and batch-driver output; the per-tool stdout/stderr is in `<Job>/Logs/<NNN>_<ToolName>.log`.

**Render and hand over individual steps as they finish.** `bp-visualize <Job>/<NNN>_<Tool>` writes `<...>/_extras/<NNN>_<Tool>_view.html` — one self-contained page for that step, built from its exported `ToolOutputs/` manifest and whatever is on disk. It needs no scheduler and no pipeline script, so it runs on the login node against a job that is still in progress. Copy that one file back and open it; do not wait for the whole run to describe a finished step in prose.

**Hand over `pipeline.html` once the run is done.** Refresh it on the cluster first, then copy that single file back and give the user its local path — it is self-contained, so it needs nothing else beside it and no network. This is the one artifact to pull whole; everything else should be a subset chosen for inspection.

Outputs land at the path configured on the cluster (typically `/shares/<group>/<user>/BioPipelines/<Project>/<Job>_NNN/`). Inside that job folder: one `<NNN>_<ToolName>/` directory per step, the completion markers, `Logs/` and `RunTime/` — so `<RunTime>` above is `<...>/<Job>_NNN/RunTime`.

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

- Produce a Jupyter notebook (`.ipynb`), not a `.py` script. Example notebooks live in `example_pipelines/notebooks/`, and personal development is in `my_pipelines/`
- First cell: clone the repo, install BioPipelines, and call `.install()` on
  the tools the pipeline needs. See the "Google Colab" section in
  `docs/user_manual.md` for the canonical setup snippet.
- The pipeline runs **inline as cells execute** — there is no `./submit`,
  no SLURM, no two-phase configure-then-execute. Treat the pipeline as
  ordinary notebook code that happens to use the BioPipelines API.
- Inputs: small files can be uploaded directly to the runtime
  (`/content/`); larger or persistent inputs should be mounted from Drive.
- The auto-detected config is `config.colab.yaml`. `Pipeline` takes a `config=` argument but never stores it, so `Pipeline(...).config` raises `AttributeError` — to confirm which variant is active, run `from biopipelines.config_manager import ConfigManager; ConfigManager().get_variant()` in an early cell, or `!bp-config show`.

### Running the pipeline

There is nothing to ssh into. The "running" step is just executing the
notebook cells in order. The notebook carries its own record too — cell
outputs are saved with the `.ipynb` — on top of the run's own
`<Job>_NNN/_operations.jsonl`.

**Who executes the cells depends on whether the Colab MCP server is registered.** If `mcp__colab-mcp__*` tools are in this session's tool list, you create the notebook, add cells and run them yourself, and read the outputs back directly — see `skills/biopipelines/references/colab_backend.md` for the full workflow. If they are not, you cannot touch the runtime: hand the user the notebook, ask them to run it, and work from what they paste back. Check the tool list; do not assume either.

Things to remind the user about during execution:

| Concern                | What to do                                                                  |
| ---------------------- | --------------------------------------------------------------------------- |
| Save the notebook      | File → Save in Drive, periodically. The autosave is *to Colab*, not Drive. |
| Save outputs           | All tools implement `.download()` — no arguments; it zips the tool's output folder and triggers a browser download on Colab. |
| Long-running tools     | Run the long cell, then check on it — Colab will idle-disconnect a tab.    |
| Reinstalling each session | Tools install once per kernel; budget ~5–15 min of setup at the start.   |
| See the run page       | `regenerate_pipeline_page("<Job>/RunTime")`, then display it inline or download it. The 3D viewers are embedded, so no network is needed. |

Outputs land in `/content/BioPipelines/<Project>/<Job>_NNN/` by default,
which is **ephemeral**. Anything the user wants to keep must be copied to
Drive or downloaded.
