# Scripts

Default folder searched by the `Scripting` tool. With
`folders.infrastructure.scripts` pointing here, `Scripting("foo.py", ...)`
resolves to `my_scripts/foo.py` without needing an absolute path.

- **Not committed**: Your scripts here are gitignored (except this README)

A script defines `configuration(inputs)` and `execution(inputs, outputs)`;
see `biopipelines/scripting_api.py` for the contract.

## File Organization

```
my_scripts/
├── README.md         # This file (committed)
├── my_step.py        # Your scripts (not committed)
└── ...               # Your scripts (not committed)
```
