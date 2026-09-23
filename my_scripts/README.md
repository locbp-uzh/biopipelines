# Scripts

Default folder searched by the `Scripting` tool. With `folders.infrastructure.scripts` pointing here, `Scripting("foo.py", ...)` resolves to `my_scripts/foo.py` without needing an absolute path.

- **Not committed**: your scripts here are gitignored, so anything you drop in this folder stays on your machine. This README, `_template.py` and the three examples below are the exception and **are** committed. The `_` prefix marks the template as the folder's own furniture — copy it to start a step, do not edit it in place.

A script defines `configuration(inputs)` and `execution(inputs, outputs)`; see `biopipelines/scripting_api.py` for the contract.

## Examples (committed)

- `_template.py` — the starting point: passes structures through unchanged and records one number per structure, so it shows a file stream and a table together and runs against any pdb stream. `my_pipelines/_template.py` calls exactly this.
- `count_atoms.py` — table-only output (`outputs["count"].row(...)`).
- `filter_structures.py` — file-stream output (`outputs["structures"].file(id, name)`) with `outputs.drop(...)` for filtered ids.
- `align_protein_openings.py` — a real step rather than a teaching example: places a domain into the opening left by cutting another domain out of a scaffold, by facing the two openings across a common plane. Shows how to take parameters when `Scripting` accepts no free-form scalars (a table input, read per-id with `params.row(id)`), how to emit a file stream and a table together, and how to keep a script usable standalone via `argparse` under `__main__`. **Call it with `env="ProteinEnv"`** — its neighbour search is `scipy.spatial.cKDTree` and the default `biopipelines` env has no scipy, so the step fails only after the upstream folds have been paid for.

The first two use only the inputs/outputs proxies (no biopipelines import in the script). Note the runner still imports the framework to resolve inputs, so any `env` must carry biopipelines' deps — the default biopipelines env does, and a custom env does after `pip install -e ".[scripting]"` in the repo. That same install lets `execution` call framework helpers (`get_mapped_ids`, `DataStream`) directly.

## File Organization

```
my_scripts/
├── README.md              # This file (committed)
├── _template.py           # Starting point to copy (committed)
├── count_atoms.py         # Examples (committed)
├── filter_structures.py
├── align_protein_openings.py
├── my_step.py             # Your scripts (not committed)
└── ...                    # Your scripts (not committed)
```
