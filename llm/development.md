# Working agreement — developing the BioPipelines framework

This file is the working agreement for **framework-development sessions** on
this repo: extending tools, fixing bugs, refactoring, improving
infrastructure. (For *using* the framework on a specific biological problem,
see `llm/pipelines.md`.)

**How to use this file.** Absorb the conventions below and apply them as
relevant — do not recite them back to the user, do not announce that you've
read them, do not perform the protocol as a checklist before the user has
described a task. Wait for the actual request, then act.

## Read at session start

Before responding to the user's first task, read these files end-to-end so
you have the framework's contract loaded — what streams are, how typed I/O
works, what `pipe_*` scripts are, which shared libraries handle which
concerns:

- `docs/user_manual.md`
- `docs/developer_manual.md`

Then read `docs/tool_reference.md` when the task involves writing or modifying a tool — use it to find, per design axis, the existing tools that exhibit the pattern you need. The skeleton to start from is `biopipelines/_TEMPLATE.py`, never a copy of one existing tool's file; see step 3 of *Order of operations for a new tool* below for why.

Do not announce that you've read them, do not summarise
them, do not list what you found. Just have the context loaded so your
proposals are grounded in the framework's actual contract instead of
inferred guesses.

## Aligning with the user (only when the request leaves choices open)

Most one-line directives ("rename `Pool.runs` to `Pool.jobs`", "fix the
typo in `tool_reference.md`") need no interview — just do them. Run the
interview only when the request leaves real choices open (a new tool, a
non-trivial refactor, an interface change).

**Never ask the user to invent things the framework already decides**
(typed streams, output shape, env-manager choice, …) — those are in
`developer_manual.md`, which you've already read. Skip anything the user
has already pinned.

**Do not narrate procedural steps.** No "I'll check the analogue first",
no "Let me look at X", no "Now I have a good template", no "want me to
read Y first?". Reading files, searching the codebase, picking an
analogue — those are internal work, the user shouldn't see them. The only
user-visible output during the interview is the substantive question(s)
or the proposed sketch.

**Do not write narrative comments in the codebase** Keep comments to a
minimum to understand and use the codebase.

**Do not use `I` and `you` when asking questions** Refer to yourself as
"The coding agent" and to the user as "The user".

### Order of operations for a new tool

**Read these three before writing any code — they are the contract, not background reading:**
- `biopipelines/_TEMPLATE.py` — the skeleton to copy to `biopipelines/<yourtool>.py`.
- `biopipelines/_TEMPLATE_CHECKLIST.md` — the step-by-step checklist. Work it top to bottom; §9 (Registration) is what CI and the pre-commit hook enforce, so skipping it means a red pipeline or a rejected commit, not a stylistic lapse.
- `pipe_scripts/_pipe_template.py` — the annotated execution-time counterpart.

1. **Settle the branch decision first** — see "Branch and PR
   conventions" below. Default to proposing a dedicated `tool-<name>`
   branch; do not touch any file until the user has chosen.

2. **Confirm the upstream reference.** Until the user confirms a URL
   or uploads files, you cannot meaningfully sketch anything: the
   tool's I/O, parameters, and env all depend on what the upstream
   actually is. Ask one focused question — "Is this `<your-best-guess>`
   at `<URL>`, or a different one?" — and **stop until they answer**. Do
   not search for analogues, do not propose an interface, do not write
   anything else. **Never clone a same-named repo from the internet on
   your own** (malicious lookalikes are a real risk for ML tooling).

3. After the reference is confirmed, do the rest internally: read the
   upstream docs the user shared, then **identify the patterns the new
   tool needs**, picking each axis independently from existing tools:

   - **Inputs** — which typed streams the tool consumes
     (`compounds`, `sequences`, `structures`, MSAs, multiple inputs, …).
     Tables? Table references?
   - **Configuration** — how the wrapper hands inputs and parameters to
     the upstream binary: bash variable injection, CLI flags, JSON/YAML
     config file written at config time, environment variables, or a mix.
   - **Execution pattern** — single-shot CLI, batched CLI, server-backed,
     in-process Python, or multi-stage.
   - **Outputs** — streams, tables.
   - **Env** — shared env, dedicated env, colab case.
   - **ID handling** — passes IDs through, generates new IDs, fans out to
     derived IDs.

   For each axis, find one or two existing tools in `tool_reference.md` that exhibit that pattern cleanly and use them as the reference **for that aspect only**. Start the file from `biopipelines/_TEMPLATE.py` and do not copy a single tool's file wholesale — that imports its quirks along with its shape. The new tool is an assembly of patterns picked per-axis, not a clone of one neighbour.

   **Most importantly, get the configuration-vs-execution split right.**
   This is the deepest mistake to avoid. Every line of code in the tool
   wrapper class runs at *configuration time* and must be bash-emitting
   only — no Python computation, no I/O on real data, no calls that
   depend on upstream outputs existing. Anything that runs against actual
   data (parsing, validating, transforming, building map tables) belongs
   in a `pipe_*.py` script that runs at *execution time*. See the
   "Two-Phase Execution" and "IDs: Configuration Time vs Execution Time"
   sections of `developer_manual.md`.

   **Output-folder layout.** Never `mkdir`; route every path through a
   `BaseConfig` helper (`configuration_path` / `execution_path` /
   `stream_folder` / `stream_map_path` / `table_path` / `extras_path`).
   See `developer_manual.md` → "Path Descriptors and the Canonical
   Layout" for the full table.

   **Don't copy canonical artifacts.** Stream `map_table`s and other
   upstream files are authoritative — pass their paths directly to the
   binary. Materialise a derived view under `_configuration/` only if
   the binary genuinely refuses the original schema, and say so in a
   comment.

   **One exception is mandatory, not optional: a filtered stream.** A
   `map_table` is the full warehouse; a filtered stream narrows `ids` but
   leaves the table intact. A consumer that reads the raw `map_table`
   (`cp`, `read_csv`, `os.listdir`) ignores the filter and reprocesses
   every row. So: if the consumer resolves through the DataStream at
   runtime (`load_datastream` + `iterate_files` / `iterate_values` /
   `Resolve.*`) it is filtered by construction — prefer that. Only if the
   binary needs a concrete CSV does the tool bypass datastream resolution,
   and then it **must** emit `self.generate_filtered_map_table_block(...)`
   before the consumer step instead of feeding the raw table. See
   `developer_manual.md` → "Filtered streams vs. the raw map_table".

4. Then present **a single concrete sketch** in one message — class name,
   input streams, output streams, key parameters with each marked
   **implemented** or **deferred** — and ask the user to confirm or
   correct it. Half-finished parameters in a public API are worse than
   missing ones, so the implemented/deferred split must be explicit.

   **Before adding any parameter that selects fields, columns, or layouts
   of a typed stream — stop.** Stream schemas are framework-fixed: a
   `compounds` stream always has its SMILES in the `smiles` column, a
   `sequences` stream always has its sequence in the canonical column,
   etc. Tools read those columns by their canonical names. If you find
   yourself sketching `smiles_column=...`, `sequence_column=...`,
   `id_field=...`, `output_layout=...`, or any similar selector over an
   invariant schema, the parameter is almost certainly wrong — either
   the framework already pins the answer or you should ask the user
   instead of inventing a knob. When in doubt, check `developer_manual.md`
   or ask.

5. Ask scope as one closed question (see *Asking questions* below).

6. **Write the wrapper and its `pipe_*` script**, following `_TEMPLATE_CHECKLIST.md` §2–§10.

7. **Register the tool.** This is not optional polish: `tests/test_registry_consistency.py` runs in both CI pipelines and fails the build on any of these, and `versions/check_tool_edits.py` runs as a pre-commit hook and rejects the commit. Work `_TEMPLATE_CHECKLIST.md` §9, which lists all of them — the ones most often missed:
   - `versions/tool_changelog.yaml`: a new entry with `files:`, `current:` equal to the class's `TOOL_VERSION` literal, and a `history:` entry.
   - `versions/CHANGELOG.md`: a bullet naming the tool under `[unreleased]` → `### Tools`. Note the path — there is no `CHANGELOG.md` at the repo root.
   - `docs/tool_index.md`: a table row **and** a bump of the `**Public-API count: N**` line; the row's version cell must equal the source `TOOL_VERSION`.
   - `docs/tool_reference.md`: an entry (the test greps for the tool name).
   - `docs/tool/<category>.md`: a real `##` or `###` heading whose text starts with the tool name. Prose mentioning the tool is not enough.
   - `README.md`: a tool-table row using the exact `<td><sub><b>Name</b>` markup the test greps for.
   - `environments:` entries in **three** config variants — `config.cluster.yaml`, `config.colab.yaml`, `config.container.yaml`. The test compares colab and container against cluster and fails on any gap; `config.daint.yaml` is exempt (aarch64).
   - `biopipelines/__init__.py`: export the class (not CI-enforced, but the tool is unreachable without it).

8. **Run the gates locally before pushing** — see *Tests and gates* below.

### Tests and gates

- **The local gate is `python -m pytest tests/ -q`.** Both CI pipelines (`.gitlab-ci.yml` and `.github/workflows/tests.yml`) run `pytest tests/` on Python 3.10–3.13, and on GitLab a red test stage blocks the mirror to the public GitHub repo. Run it before pushing.
- **The registration gate is `tests/test_registry_consistency.py`.** It is the test that fails when step 7 above is incomplete. Run it alone with `python -m pytest tests/test_registry_consistency.py -q` for a fast check while wiring a new tool up.
- **The per-tool parameter suites are opt-in.** `pyproject.toml` sets `addopts = "-ra -m 'not tool_parameters'"`, so they are excluded by default; they are GPU-bound and run with `pytest tests/tool_parameters -m tool_parameters -v`. Do not expect them to have run just because `pytest tests/` was green.
- **`pre-commit install` is required once per clone.** The `TOOL_VERSION` hook lives in `.pre-commit-config.yaml` and runs `python versions/check_tool_edits.py`; on a fresh clone that hook is not installed, so commits that skip the version bump go through silently and CI catches them later:
  ```bash
  pip install pre-commit
  pre-commit install
  pre-commit run --all-files   # on demand
  ```
  The hook fires whenever a staged file is listed under some tool's `files:` in `versions/tool_changelog.yaml`, and then demands the bumped `current`, a matching `TOOL_VERSION` literal in the wrapper class, and a staged `versions/CHANGELOG.md` mentioning the tool under `[unreleased]` → `Tools`.

### Asking questions: format

When you need user input, use **closed-option questions** over open-ended
prose. The user should be able to answer in 30 seconds by selecting an
option (or saying "other" with a short note).

**First check whether `AskUserQuestion` is in your tool list for this
session.** If it isn't, you don't have it ask questions in markdown style.
If it is available, **you must use it** for any question that has 2–4
discrete options. Markdown numbered lists are the fallback for hosts
without the tool, not a stylistic preference when the tool is present.

`AskUserQuestion` constraints: 1–4 questions per call, each with 2–4
options. Each option's `label` is 1–5 words; longer text goes in
`description`. Put the recommended option first and append "(Recommended)"
to its label. The tool auto-adds an "Other" option, so don't include one.

Either way, the question text and the option set should be identical
across renderers — the format is a presentation choice, not a content
choice.

## Environment assumptions

- This session runs on a **local machine, not on the cluster**. Project
  folders, datasets, and conda/mamba envs are not accessible. GPU-bound code
  cannot be executed locally — verify it via review and CI, not by running it.
- Never install packages into the user's environments without explicit
  consent. If possible rewrite a script to avoid the missing dep.
- When writing install specifications, always retrieve the official
  instructions from the given repository.

## Branch and PR conventions

**As soon as the task is clear, ask the user where the work should
land** — before reading files, before sketching, before any edit. The
question has two options; switching to a dedicated branch is the
recommendation, staying on the current branch is the alternative.

- **(a) New dedicated branch from `main`** *(Recommended)* — short
  `<kind>-<slug>` name. Examples: `tool-esmfold`, `fix-mpnn-paths`,
  `docs-userman-rewrite`, `refactor-pool-runs`. Wait for the user to
  confirm the proposed name before creating it.
- **(b) Stay on the current branch** — only when the user explicitly
  prefers it (e.g. quick fix, exploratory tweak, branch already
  scoped to this work).

Once the branch decision is made, proceed with the rest of the work.
Do not touch any file until this is settled.

Open PRs against `main` unless told otherwise. CI runs on push — that
is the test signal for code correctness.

## Code rules

- When implementing or refactoring, never introduce unnecessary fallbacks
  or guessed defaults. Validate at boundaries; trust internals. (Other
  framework-level rules — typed streams, bash-only config, shared
  libraries, etc. — are in `developer_manual.md`; comply with them.)
- **No selector parameters over invariant stream schemas.** Don't expose
  `<thing>_column`, `id_field`, `output_layout`, etc. — typed-stream
  schemas are fixed by the framework, and a parameter that lets the user
  override them is almost always a bug, not a feature. See the
  parameter-sketch rule under *Order of operations for a new tool*.

For new tools, the per-axis pattern identification in step 2 of the
order of operations above is internal work, not a user-visible
deliverable.

## Testing the change (when the change warrants it)

Trivial edits (typos, doc tweaks) need no end-to-end test. For any code or
config change, the section below applies.

Framework code should behave identically on cluster and Colab — they share
the same Python API and runtime logic. The two genuine differences are:

- **Where you can drive the test from.** On the cluster, you can operate end to end yourself (push → sync → submit → inspect → iterate) over `log.sh ssh`. On Colab it depends on whether the Colab MCP server is registered: if `mcp__colab-mcp__*` tools are in this session's tool list you can create the notebook, run its cells and read the outputs back yourself (see `colab.md`); if they are not, you cannot touch the runtime and the user has to execute the cells and paste back outputs. Check the tool list rather than assuming. Either way the cluster stays the default verification venue whenever it is available, since it needs no notebook round-trip at all.
- **Install paths can diverge.** `.install()` scripts, `config.colab.yaml`
  env entries, and micromamba behavior are the one area where Colab needs
  separate testing even when the runtime logic is unchanged.

**Never regress an already-debugged Colab path to add cluster support.** If a tool already has a verified `<tool>.colab.yaml`, do not edit it for the cluster — keep it and add a generic `<tool>.yaml` (picked up off Colab). Likewise in code: don't replace a working Colab step; put it under `if scheduler == "colab":` and add the cluster logic under `else:`.

### Cluster verification (default)

If your change touches Resources defaults, partition selection, or anything
that depends on what the cluster actually provides, consult `llm/resources.md`
(populate it first if absent — see `llm/pipelines.md` → "Probing the cluster").
Do not bake site-specific assumptions into framework code without checking.

You cannot run GPU code or SLURM submissions locally. To exercise a change
end to end, ssh into the cluster. **Every `ssh` and `scp` call must be
wrapped in `llm/log.sh`** — no raw cluster commands. The wrapper writes the
command and its output to `llm/logs/YYYY-MM-DD.log`, which is the audit
trail for the session. Read `cluster.md` for more information.

1. Push the branch.
2. On the cluster: `cd <repo> && git fetch && git checkout <branch> && git reset --hard origin/<branch>`. Confirm the branch and the hard-reset target with the user first — this is destructive on the remote checkout.
3. On the cluster: `./submit <pipeline.py>` with a minimal test pipeline.
4. Inspect the logs. `ls <RunTime>/*.out` first: a pipeline with more than one `Resources()` call is submitted as one job per batch and writes `<RunTime>/job_batch<N>.out`; only a single-batch pipeline writes `<RunTime>/slurm.out`. Then read the per-tool log at `<Job>/Logs/<NNN>_<ToolName>.log`, and check the completion markers (`<NNN>_<ToolName>_COMPLETED` / `_FAILED` / `_WARNING`) sitting directly in the job folder.
5. **Inspect the actual output files, not just the success marker.** "completed successfully" only means the completion-check found the expected paths — not that the result is right. Open the produced files: a structure PDB has plausible coordinates, a scores/affinity table has finite non-empty values, an N-sample run wrote N files. Also sanity-check timing — a GPU job that took minutes for a tiny input likely fell back to CPU.

Two recurring traps:

- **Probe before assuming an install bug.** When a tool errors on real data, check the *inputs* first (coordinates, pose, format) — a scorer's empty-pocket crash was a ligand placed at the origin, not a broken env. A one-shot diagnostic (print the relevant state) beats re-running the whole job on a guess.
- **`gpu="any"` can land on an H100 (sm_90).** Tools whose env pins torch ≤cu117 have no sm_90 kernel and crash ("no kernel image") or silently run on CPU. Use `gpu="A100"` for cu11x-pinned tools. See `llm/resources.md` → Gotchas.

To test pipelines, write them under `my_pipelines` (gitignored folder) and scp, instead of going through the push-sync cycle.

### Colab verification (only if install paths are touched)

If the change modifies any `.install()` script, `config.colab.yaml`, or the
micromamba env setup, also verify on Colab — cluster verification will not
catch install-side regressions.

**If `mcp__colab-mcp__*` tools are in this session's tool list**, drive it yourself: follow the workflow in `colab.md` → "Workflow — running a test pipeline on Colab", which builds the same notebook and executes the cells over MCP.

**Otherwise** you cannot drive Colab from this session; hand the user a notebook with:
1. Replaced git clone line:
   ```bash
   !git clone -b <branch> https://github.com/<org>/biopipelines.git
   ```
2. the affected `.install()` call
3. a minimal pipeline cell

and ask them to:

1. Disconnect and delete the current runtime, then reconnect (a fresh
   kernel is required — re-running on a kernel that already has the env
   masks broken `.install()` logic).
2. Run all the cells.
3. Paste back the cell outputs, plus the tool's log — `<Job>/Logs/<NNN>_<ToolName>.log`, or the same stream as `_log` inside the tool's own output folder — and which completion marker was written (`<NNN>_<ToolName>_COMPLETED` / `_FAILED` / `_WARNING`, one level above the tool folder). Or download the artifacts with `<tool>.download()`.

## Reporting back

- State what was changed and why, in terms of the framework's invariants.
- If you skipped a verification step (e.g. couldn't run GPU code), say so
  explicitly rather than implying success.
