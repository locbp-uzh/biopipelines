# `bp-mcp` — the BioPipelines MCP server

An agent with the skill alone can *read about* BioPipelines. With `bp-mcp` registered it can *query* it: ask what tools exist and pull any tool's documentation without opening a file.

The server runs on **your machine**, not the cluster. It needs no install rights there, leaves no long-lived process on a shared login node, and keeps ssh credentials in your own agent. A cluster-side deployment stays possible later; from the server's point of view that is just the case where the ssh hop is local.

## Install

The MCP SDK is an optional extra — authoring and running pipelines does not need it.

```bash
pip install -e ".[mcp]"
```

## Register

Claude Code, at user scope so it is available in every project and not committed:

```bash
claude mcp add --scope user biopipelines -- bp-mcp
```

Verify with `claude mcp list` — it should show `biopipelines … ✓ Connected`.

For any other host, the command is `bp-mcp` over **stdio**. If the entry point is not on `PATH` (an uninstalled checkout), `python -m biopipelines.mcp_server` is equivalent.

## Getting started

Three things have to exist before the tools can reach a cluster, and only the first is outside their reach:

1. **An ssh alias** in your `~/.ssh/config`, named after the site rather than after what kind of machine it is (`Host s3it`, `Host daint` — `cluster` stops telling you which machine once there are two). These tools cannot create it; ssh config is yours.
2. **A BioPipelines checkout on the cluster**, wherever you keep it.
3. **A config variant on the cluster that claims your username there** — `bp-config auto --variant cluster`, run once on the cluster.

Then, from the agent:

```
bp_setup(host="s3it", repo="~/biopipelines-locbp")
```

It walks the whole chain and stops at the first thing missing, so the report names what to fix rather than surfacing it as a config error three steps later:

```
Ready: s3it (~/biopipelines-locbp)
  [ok] ssh alias reachable: ok
  [ok] repository at ~/biopipelines-locbp: found
  [ok] biopipelines importable: 1.4.0
  [ok] config variant resolves: /shares/<group>/<user>/BioPipelines
  [ok] output root readable: 12 projects
```

On success it saves that host's settings to a gitignored `.bp-mcp.json`, and every later call can omit them. Pass `host=""` explicitly to read this machine instead while a cluster is saved.

## More than one cluster

Settings are kept **per host**, so a lab with access to two clusters configures both and nothing carries over between them:

```
bp_setup(host="s3it", repo="~/biopipelines-locbp")
bp_setup(host="daint", repo="$SCRATCH/biopipelines", variant="daint",
         python="$SCRATCH/venvs/biopipelines/bin/python", make_default=False)
```

`bp_setup()` with no arguments lists what is configured and marks which one an unnamed call uses; `make_default=False` adds a cluster without changing that, and `forget=True` drops one after a renamed alias or a typo. Every tool then takes `host=` to reach either.

Three things are per host rather than global, and each has a reason:

- **`repo`** — the checkout is in a different place on each machine (on Daint it has to live on `$SCRATCH`, whose inode quota is the only one that fits the environments). An earlier version kept one repo for whichever host was named, which sent one cluster's path to the other: it resolves to nothing, or worse to somebody else's tree.
- **`python`** — Daint's login nodes have no `python` at all; it is a venv on `$SCRATCH`.
- **`variant`** — the resolver and `./submit` both call the framework, which otherwise takes whatever the remote auto-detects. Generating a Daint run against the `cluster` variant produces scripts for a scheduler layout that is not there. Pinning it exports `BIOPIPELINES_CONFIG_VARIANT` on every remote command for that host.

## Tools

| Tool | Arguments | Returns |
|---|---|---|
| `bp_tools` | none | The full catalog: one line per tool with version, CPU/GPU, verified platforms, and purpose |
| `bp_tools` | `name` | That tool's complete documentation section, parameters included |
| `bp_tools` | `outputs` / `inputs` | Tools found by what they RETURN or CONSUME — a stream, table or column name, case-insensitive substring, read from the tools' own source |
| `bp_tools` | `tags`, `exclude` | Tools found by intent, over the closed vocabulary in `docs/tool_tags.md` |
| `bp_setup` | `host`, `repo`, `python`, `variant`, `make_default`, `forget` | Checks the connection chain and saves a working one. A `prelude` (shell run before every command) is set in `.bp-mcp.json` by hand, never through the tool |
| `bp_runs` | `host`, `project` | Projects, or the runs inside one, under that machine's output root |
| `bp_status` | `job`, `host`, `project` | Every step with its completion marker, the counts, the run header, and the first failure |
| `bp_logs` | `job`, `step`, `tail`, `host` | The tail of one step's log. `step` is `001_Tool`, or the tool name alone when unambiguous |
| `bp_ancestry` | `job`, `id`, `csv`, `host`, `project` | Which design came from which input. With `id`, walks one design back to its inputs showing every parent; without, summarizes the whole graph. `csv=True` writes `ancestry.csv`, one row per link |
| `bp_lineage` | `job`, `page`, `host`, `project` | Per step: IDs produced per stream and IDs dropped, for a whole run. `page=True` writes a self-contained HTML page instead — identity, attrition, resolved parameters and every dropped ID with its reason |
| `bp_table` | `job`, `step`, `table`, `limit` | One step's CSV as rows; omit `table` to list what that step wrote |
| `bp_submit` | `script`, `upload`, `verbose` | Runs a pipeline on the cluster; returns the run directories and job ids. A timeout is reported as an unknown outcome, not a failure: check `bp_runs` before submitting again |
| `bp_fetch` | `remote_path`, `local`, `overwrite` | Copies a file or a whole step folder back from the cluster; refuses to overwrite an existing local file unless `overwrite=True` |
| `bp_resubmit` | `job`, `script`, `keep_dependencies` | Resumes a run that failed partway. SPENDS COMPUTE |
| `bp_cancel` | `job`, `job_ids`, `confirm` | Without `confirm=True`, names the jobs still queued or running and cancels nothing; with it, cancels only those. DESTROYS WORK |
| `bp_visualize` | `job`, `step`, `descending` | One step's images inline, plus its interactive page fetched locally |
| `bp_project` | `project_dir`, `action` | `survey`, `scaffold`, or `history` on a project folder's documents |
| `bp_provenance` | `job`, `against` | what produced a run — versions, resolved parameters, environment digests; with `against`, how two runs differ |
| `bp_reproduce` | `job`, `target` | plan a rerun: the preserved script plus every difference the target would introduce. Submits nothing |

Three ways to find a tool, for three different questions. `name` when you know it. `outputs="rmsd"` when you know what the number is called — it matches stream, table and column names extracted from `get_output_files()`, so it is exact rather than a keyword guess over prose, and it prints each tool's `**Units.**` note alongside, which matters because Boltz2's affinity is log10(IC50) where lower is stronger while GEMS's is a pKd where higher is. `tags=[...]` when you only know what you want to achieve: terms in one facet are OR-ed and different facets AND-ed, and `exclude=[...]` drops tools carrying a tag — `tags=["protein", "binding"], exclude=["small-molecule"]` is how you reach Prodigy, and `exclude=["fetch"]` is how you ask what runs with no network. An unknown tag prints the whole vocabulary rather than an empty list.

An unrecognized `name` returns the closest documented names rather than an error, so a near-miss like `Boltz` answers with `Boltz2, BoltzGen` instead of leaving the model to invent a signature. `bp_logs` with a wrong step lists the steps that do have logs.

`bp_status` distinguishes the two failure shapes, because they call for different responses: everything after the failure pending means the run stopped there; later steps completed means the failure did not halt the pipeline. Guessing either way sends an agent looking for a cascade that may not exist.

## Project documents

A project folder holds its `<Job>_NNN/` directories beside two documents:

* **`PROJECT.md`** — goal, inputs, and the validated protocol: what we run and why we believe it. Rewritten as understanding improves, never appended to.
* **`HISTORY.md`** — append-only and dated: runs, findings, restructures.

`bp_project` surveys them, scaffolds what is missing, and appends history entries. It **never overwrites an existing file** — running `scaffold` on an established project adds only the document it lacks. `PROJECT.md` carries the science: an agent drafts and proposes, the scientist owns the content.

Both views read `docs/tool/*.md` through `biopipelines/tool_docs.py` — the same module that generates `references/tool_index.md`. The server and the shipped index cannot disagree about what a tool is or where its documentation lives.

## These tools replace the shell, they do not supplement it

Once the server is registered, reaching for `ssh`, `scp` or a script in the repo to do what a tool does is a regression, not a matter of taste. `ssh s3it "./submit ..."` skips the `_operations.jsonl` entry that `bp_submit` writes beside the run; `tail`ing a log by hand means first guessing whether the run produced `slurm.out` or `job_batch3.out`, which `bp_status` already resolved; opening `docs/tool/analysis.md` to read one signature costs 24k tokens for what `bp_tools(name=...)` returns in a few hundred.

**`bp-submit` and `bp_submit` are different things.** The hyphenated names (`bp-submit`, `bp-config`, `bp-run`, `bp-visualize`) are console scripts a person runs in a shell. The underscored names (`bp_submit`, `bp_status`, …) are the MCP tools. With the server registered, an agent calls the underscored ones; it does not shell out to the hyphenated ones for the same job.

Scripts under `skills/` and `versions/` are maintainer utilities. `build_tool_index.py` regenerates the catalog file after someone edits the prose docs — it is not how a tool is looked up, and running it during a pipeline session changes a tracked file for no reason.

Shell out for what genuinely has no tool: `sinfo` and `squeue`. Cancelling, resubmitting and visualizing are `bp_cancel`, `bp_resubmit` and `bp_visualize`.

## Why this matters more than it looks

An agent that guesses a parameter gets no error: BioPipelines tools absorb unknown keywords through `**kwargs`, so a typo becomes a silently ignored argument and the run proceeds with a default nobody chose. `bp_tools` exists so there is never a reason to guess.

## Reaching the cluster

The cluster is the main platform, so every run-reading tool takes `host` — the ssh alias from `~/.ssh/config`:

```
bp_runs(host="s3it")
bp_status(job="Binder_design_006", host="s3it")
bp_logs(job="Binder_design_006", step="009_Boltz2", host="s3it")
```

**The remote resolves its own output root.** `biopipelines_output` is defined per machine and expands `<username>` there, so the local config would answer for the wrong user. With a `host`, the tool runs the resolution on the far side with a self-contained snippet over `ConfigManager` (so an older cluster checkout still answers) and reads from what that returns. Pass `repo` if the checkout is not at `~/biopipelines`.

`job` may be a bare job name or a full path. A bare name is looked up under `<output_root>/<project>/`, or searched one level down when no `project` is given; a name found in more than one project is refused with the candidates listed, so pass `project=` then.

Every ssh call runs with `BatchMode=yes`, so a host key or passphrase prompt fails at once instead of hanging, and ssh's own failure (exit 255) is reported as an error rather than as a missing file. `host` must be an alias, never an option; remote paths handed to scp are limited to letters, digits, `._/+@=,-`, a leading `~` or `$VAR`, and no `..`.

Without a `host`, everything reads this machine's filesystem — the local, container and Colab variants, and a cluster share mounted locally.

**Read-only queries are not written to the run's operations log.** Recording every status poll is what made `llm/log.sh` useless: a `squeue` and a six-hour submission became the same kind of line. Only actions that change something earn a record.

## Submitting

```
bp_submit(script="my_pipelines/foo.py", upload="my_pipelines/foo.py")
```

`upload` copies the local file to `<repo>/<script>` first, then `./submit` runs there. The result gives the run directories created and the scheduler job ids.

**An accepted job is not a finished job.** `bp_submit` returns as soon as the scheduler takes it; nothing has executed. Poll `bp_status` rather than reporting success. The tool description says so too, because this is the easiest thing for an agent to get wrong.

It is also the first tool that **writes** to `<Job>_NNN/_operations.jsonl`: a run header (framework version, commit, config variant, host) and the submission itself. Those records are made before the run directory exists and flushed once `./submit` names it, which is what `RunLog`'s buffer is for. A scheduler rejection is reported as a failure even though a run directory was created — the directory existing is not the job being queued.

## Where the designs went

`bp_lineage` answers the question a campaign actually turns on — not "did it run" but "I started with 500 and have 12, which step removed the rest":

```
binder-campaign_001: lineage across 8 step(s) that report IDs
  004_Filter_stitched_001    [-26 dropped]
  022_Filter_stitched_004    [-40 dropped]
  ...
142 ID(s) dropped across 8 step(s); the largest single loss is 40 at 022_Filter_stitched_004.
```

Both numbers come from tables the framework already writes — `<stream>_map.csv` for the IDs a step produced, `missing.csv` for the ones that entered and did not come out — not from arithmetic on result tables. A whole run costs **one** ssh round trip, not one per step; a 253-step campaign was read in a single call.

Above 40 steps the listing keeps only the steps that lost IDs, because on a large campaign the attrition is the only part anyone reads.

`bp_table` then reads any one of those tables — `missing.csv` for the dropped ids, `analysis.csv` or `confidence.csv` for metrics — truncated to `limit` rows, since they run to thousands of lines.

## Pulling results back

`bp_fetch` copies a file or a whole step folder to this machine, deciding recursion by asking the remote rather than making the caller know which it is:

```
bp_fetch(remote_path=".../run_001/005_ProteinMPNN", local="./results/")
```

This is for the case where results exist only on the cluster — a common and easily-forgotten gap, since a run directory keeps growing after the interesting structures have been picked.

## Counts versus connections

`bp_lineage` and `bp_ancestry` answer different halves of the same question, and reaching for the wrong one wastes a call.

- **`bp_lineage`** — *how many*. Produced and dropped per step, from the map and missing tables. Works on every run, including old ones, because line counts need no provenance columns.
- **`bp_ancestry`** — *which came from which*. The parent graph, read from the `<axis>.id` columns of the same map tables and scoped by the wiring in `RunTime/pipeline_graph.json`.

Report a campaign's yield with the first; answer "prove this design came from that structure" with the second. A run whose tables predate the provenance columns reports its IDs as roots rather than inventing links, so an empty ancestry on an old run is a fact about the run, not a failure of the tool.
