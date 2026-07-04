# Scripts

Default folder searched by the `Scripting` tool. With `folders.infrastructure.scripts` pointing here, `Scripting("foo.py", ...)` resolves to `my_scripts/foo.py` without needing an absolute path.

- **Not committed**: Your scripts here are gitignored (except this README and the committed examples).

A script defines `configuration(inputs)` and `execution(inputs, outputs)`; see `biopipelines/scripting_api.py` for the contract.

## Examples (committed)

- `count_atoms.py` — table-only output (`outputs["count"].row(...)`).
- `filter_structures.py` — file-stream output (`outputs["structures"].file(id, name)`) with `outputs.drop(...)` for filtered ids.

Both use only the inputs/outputs proxies (no biopipelines import in the script). Note the runner still imports the framework to resolve inputs, so any `env` must carry biopipelines' deps — the default biopipelines env does, and a custom env does after `pip install -e ".[scripting]"` in the repo. That same install lets `execution` call framework helpers (`get_mapped_ids`, `DataStream`) directly.

## File Organization

```
my_scripts/
├── README.md              # This file (committed)
├── count_atoms.py         # Examples (committed)
├── filter_structures.py
├── my_step.py             # Your scripts (not committed)
└── ...                    # Your scripts (not committed)
```
