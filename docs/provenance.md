# Provenance — what produced a run

Every run writes `manifest.json` into its job folder. Nothing has to be switched on for this: the question *"what exactly produced this result?"* is asked after a run, not before it, and a record you have to remember to enable is a record you will not have when it matters.

## What is in it

```json
{
  "schema": 1,
  "hash": "sha256:…",
  "created": "2026-09-22T…",
  "biopipelines": "1.4.1",
  "commit": "c81bad8",
  "variant": "cluster",
  "host": "…", "python": "3.10.14",
  "scheduler": "slurm", "env_manager": "mamba", "container_executor": "…",
  "project": "BinderDesign", "job": "designs",
  "tools": [
    {
      "order": 2,
      "tool": "ProteinMPNN", "class": "ProteinMPNN", "tool_version": "2.6",
      "environments": ["proteinmpnn"], "container_image": "",
      "parameters": {
        "passed":   {"num_sequences": 4, "structures": "<stream structures>"},
        "resolved": {"num_sequences": 4, "sampling_temp": 0.1, "model_name": "v_48_020", "...": "..."}
      }
    }
  ]
}
```

The distinction between `passed` and `resolved` is the point of the file. `passed` is the constructor call as the user wrote it, positional arguments included. `resolved` is that plus every default that applied — and the defaults are usually what determined the output. A manifest that recorded only `passed` could not tell you which ProteinMPNN weights produced a set of sequences, because nobody ever names `model_name`.

`forwarded` appears for tools that set `FORWARD_UNKNOWN_KWARGS`, holding the options passed straight through to the upstream command line.

## The hash

`hash` covers everything that would make a rerun differ, and deliberately excludes `created`, `host`, `project` and `job`. Running the same pipeline twice therefore produces the same hash; changing a sampling temperature, a tool version, a container image or the site variant produces a different one. A hash that drifted on its own could not flag a change that mattered.

## Comparing two runs

```python
from biopipelines import manifest

before = manifest.read("/path/to/job_001")
after  = manifest.read("/path/to/job_002")
for line in manifest.compare(before, after):
    print(line)
```

```
step 2 ProteinMPNN: sampling_temp 0.1 -> 0.3
step 3 Boltz2: tool_version '2.1' -> '2.2'
```

`compare()` returns an empty list only when the two runs agree on everything recorded, including their environments, and otherwise names each difference concretely rather than reporting that the two files differ — the latter is the answer that sends someone diffing JSON by hand.

## Relationship to `debug=True`

`Pipeline(debug=True)` writes a richer snapshot into `_debug_capture/`: full environment exports, `pip freeze`, `nvidia-smi`, the scheduler version, and a copy of the active config file. That capture is the forensic one and remains opt-in, because it costs a subprocess per environment on the executing node.

The manifest is the cheap half — pure Python at `save()` time — and is always written. Use the manifest to ask whether two runs differ; use `_debug_capture/` to find out why an environment behaved the way it did.

## The environment digest

The manifest is written before the job leaves your machine, so it cannot know which environment will execute it. The run records that itself: on the executing node, every pipeline exports each environment it uses *with build strings* and writes the export and its SHA-256 beside the run.

```
<job>/environments/proteinmpnn.txt      the export, plus pip freeze
<job>/environments/proteinmpnn.sha256   its digest
```

The export keeps build strings deliberately. `--no-builds` is right for an environment file you intend to recreate elsewhere and wrong for identity: two different builds of the same package versions export identically, which is exactly the distinction that now matters, because the same predictor ships in more than one build.

`manifest.read()` merges these digests under `environments_resolved`, and `compare()` reports them. **Note that an equal `hash` does not imply an equal environment** — the hash is sealed at `save()`, before the environment that runs the job is known. `compare()` therefore checks the digests before it short-circuits on the hash, so "same plan, rebuilt environment" is reported as a difference rather than as a match.

A run with no digests reports `environments: not recorded` rather than nothing, and comparing a run that has them with one that does not says the environments *cannot be compared*. Silence there would claim they match, which an absent record cannot support.

## Rerunning a campaign

The manifest records each tool's resolved parameters, but a parameter that was a DataStream is recorded as a name, so the wiring is not recoverable from it. What is recoverable is the campaign itself: `save()` copies the calling `.py` into `RunTime/`, so the authored pipeline survives next to its own record.

`bp_reproduce` pairs the two. It reads the manifest, finds the preserved script, asks the target host what it would bring to a rerun — BioPipelines commit, site variant, scheduler — and names every difference before anything is submitted:

```
BinderDesign/designs — 8 steps, recorded 2026-09-14T09:12:04
  BioPipelines 1.4.1 @ c81bad8, variant cluster, slurm
  steps: RFdiffusion -> ProteinMPNN -> Boltz2 -> PoseBusters -> Panda
  environments: boltz2 @ 4f1c9a0b2e77, proteinmpnn @ 9ab3c2e10d55
  script: /…/designs_003/RunTime/binder_campaign.py

A rerun on daint would differ:
  - config variant: recorded 'cluster', target has 'daint'
  - commit: recorded 'c81bad8', target has 'e3d2a37'

That may be fine — a second cluster differs by construction. It is stated so it is a
decision rather than a surprise.

To rerun: bp_submit with script='binder_campaign.py'. THIS SPENDS COMPUTE.
```

It submits nothing. A plan with differences is not a refusal — reproducing a campaign on a second cluster differs by construction — the point is that the differences are stated rather than discovered afterwards.

!!! note "What is still not covered"
    A target that cannot be probed is reported as *not probed*, never as matching. A run whose script was not preserved is reported as having nothing to rerun from. Neither is silently treated as success.

## Which design came from which input

`manifest.json` says what produced the run. It does not say which of the twelve designs you ordered came from which structure you started with, and that is the question a referee asks. `bp_ancestry` answers it from the run's own tables:

```
bp_ancestry(job="tagged_binder_001", id="9_Panda_1")

9_Panda_1   [009_Panda]
 └─ 005_Boltz2               3kzy_A_50_2+tag  [recorded]
    ├─ 002_Sequence             tag  [matched:parent]
    └─ 004_ProteinMPNN          3kzy_A_50_2  [matched:parent]
       └─ 003_RFdiffusion          3kzy_A_50  [recorded]
          └─ 001_PDB                  3kzy_A  [recorded]
```

Two parents, not one: the complex is the designed chain **and** the tag, and only one of them is the binder. Ancestry is a DAG for that reason, and a reading that collapses it to a chain hides half a binder campaign.

The bracket is the basis of each link, and it is shown because a link whose basis you cannot see is one you cannot audit:

- **`recorded`** — the step wrote the parent id into a `<axis>.id` column of its map table. Nothing is inferred.
- **`matched:<tier>`** — no column named a concrete parent, so the child id was resolved against the upstream ids using the same matcher the framework uses at runtime to wire one step into the next. An edge found this way is the wiring that actually happened.

Both tiers are needed because a map table does not always carry a usable parent. `Boltz2` above records `proteins.id = 3kzy_A_<1..100>_<1..2>` — a combinatorial declaration, not an identity — and elsewhere the framework deliberately drops a provenance column whose values the id already carries. A join answers neither case.

Parent search is scoped to the wiring recorded in `RunTime/pipeline_graph.json`. Without that scope a predicted complex matches its grandparent as readily as its parent.

An ID whose parent neither tier resolves is reported as a root and counted, never given an invented link:

```
bp_ancestry(job="tagged_binder_001")

  001_PDB                                1 ID(s)
  002_Sequence                           1 ID(s)
  003_RFdiffusion                      100 ID(s)  <- 001_PDB
  004_ProteinMPNN                      199 ID(s)  <- 003_RFdiffusion
  005_Boltz2                           199 ID(s)  <- 002_Sequence, 004_ProteinMPNN
  009_Panda                             24 ID(s)  <- 005_Boltz2

Established by: matched 398, recorded 323
```

`bp_ancestry(job, csv=True)` writes `ancestry.csv`, one row per link, for someone who will work on it outside these tools. The provenance page carries the same graph as a panel.

!!! note "Older runs"
    A run whose map tables predate the provenance columns has no ancestry to read. Where a column is absent **and** the id does not carry its parent, the ID is reported as a root. Saying so is the point; inventing an edge would make the record worthless for its one use.

## Where it sits

- `manifest.json` — the job folder, beside `RunTime/` and `Logs/`.
- `_operations.jsonl` — the operations record: what was submitted, cancelled, resubmitted and when.
- `_debug_capture/` — the opt-in environment snapshot.
- `RunTime/pipeline_graph.json` — the declared wiring, which scopes the ancestry search.

## A page for people who do not call tools

`bp_lineage(job, page=True)` writes `outputs/_views/<job>_provenance.html`: the run's recorded identity, the attrition step by step, the parameters as they resolved (explicitly-set ones in bold, defaults plain), and every dropped ID with the reason its own step recorded. Self-contained, no network, opens from `file://` after being copied off a cluster.

This is the same information `bp_provenance` and `bp_table` return, in the form a supervisor, a reviewer or a referee can actually read. The page reads through the same filesystem abstraction as the other tools, so a run that lives on a cluster renders locally without being fetched first.

!!! warning "Exports are redacted"
    `pip freeze` renders an editable install as the URL it was cloned from, so a clone made with a token in that URL prints the token — into an output root that is group-readable on a shared cluster. Every export therefore passes through a redaction filter (URL credentials, GitLab and GitHub token prefixes, `password=`/`token=`/`secret=`/`api_key=`) before it is written, and the unredacted capture is removed rather than kept. If redaction cannot run, the export is withheld and no digest is written.

    This does not make it safe to put credentials in a clone URL. Use an SSH key or a credential helper on the cluster; the filter is a backstop, not a permission.
