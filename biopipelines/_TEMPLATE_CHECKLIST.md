# Checklist for addition of a new tool

Copy the template `_TEMPLATE.py` to `biopipelines/<yourtool>.py`, then work this list top to bottom. 

## 1. Tool's shape
Before editing, classify the tool based on the following:
- **Install:** 
    - base-env (biopipelines), no install (delete `_install_script`)
    - shares another tool's env (delegate, e.g. `return PyMOL._install_script(...)`)
    - dedicated env (adapt from `dssp.py`)
    - container (adapt from `gnina.py`: source in `environments/_containers.yaml` + `cls._container_pull_block(folders, force_reinstall)` in the install). 
    The last two require an `environments:` entry in **each** of `config.cluster.yaml`, `config.colab.yaml` and `config.container.yaml` — see §9. A container tool also gets a `containers:` destination line in `config.cluster.yaml` (image path) paired with its source in `environments/_containers.yaml`.
- **Inputs:** (determines output prediction)
    - single stream > check closest tool in function
    - multiple streams (combinatorial) > see Boltz2/ESMFold2/Gnina 

## 2. Identity & docstring
- [ ] `TOOL_NAME` unique; `TOOL_VERSION` set.
- [ ] `ENV_NAME` declared (the env this tool builds and runs in), unless the tool installs nothing and lives in the `biopipelines` env. **Nothing catches its absence** — not the pre-commit hook, not CI — the tool simply fails to install. Note it is the tool's *default* name, not the authoritative one: `_install_env` reads the `environments:` entry from the active config (see §3 and §9) and only falls back to `ENV_NAME` under `env_manager: pip`.
- [ ] Class docstring documents every input param and every output stream/table *with columns*. 

## 3. Install (`_install_script`)
- [ ] Kept / deleted / container per §1.
- [ ] If kept: gated on `_env_exists_check`; honors `force_reinstall`; verifies the binary/imports and `touch "$INSTALL_SUCCESS"` on success.
- [ ] The env name is **never** written into the install bash. Call `cls._install_env(env_manager)` and install into what it returns — that reads the same `environments:` entry `_load_environments` activates at run time, so install and runtime cannot diverge. With no entry under the active variant it raises, naming the tool, the variant and the config file; that is intended (29 tools have no `daint` entry because they do not work there), so do not add a fallback.
- [ ] If the installation diverges between cluster and Colab, shipped both `environments/<TOOL_NAME>.yaml` **and** `<TOOL_NAME>.colab.yaml`. The `.colab` one drops conda-only deps with no PyPI build. Fork the YAML rather than branching Python on scheduler. 

## 4. Path descriptors
- [ ] Every input JSON, output stream map, and table has a `Path(lambda self: ...)` descriptor - no path strings built inline.
- [ ] Every output stream returned in §4 has a matching `stream_map_path` descriptor here.

## 5. `__init__` (input contract)
- [ ] Normalizes each `StandardizedOutput`>`.streams.<name>` and passes `DataStream` through.
- [ ] Keeps the raw input handle (`self.structures`) for missing-propagation.
- [ ] Receives streams/standardized outputs, not bare values. Bare values (a sequence string, a PDB id, a SMILES) are turned into streams by the tools `Sequence`/`PDB`/`Ligand`/... - let the user pass those, so bare-value conversion lives in one place instead of every tool.
- [ ] `super().__init__(**kwargs)` is the last line.
- [ ] **Renaming a parameter?** Put the old spelling in `PARAMETER_ALIASES = {"<old>": "<new>"}` and, if it is retired rather than a synonym you intend to keep, list it in `DEPRECATED_ALIASES = ("<old>",)`. Without the alias the old key falls into `**kwargs`, never binds, and the pipeline runs with the parameter at its default — silently, or (on a `FORWARD_UNKNOWN_KWARGS` tool) as a bogus upstream flag. Two `ValueError`s at class creation catch the ways this goes wrong: an alias target no parameter names, and a `DEPRECATED_ALIASES` entry missing from `PARAMETER_ALIASES`. A reserved key (`name`) is *copied*, so the tool and the framework both see it; any other key passed under both spellings raises. Prefer a clean break when the old spelling only affected this tool's own output names, and a kept synonym when its value is a handle other code reaches for (`tool.tables.<name>`). See [Renaming a Parameter](../docs/developer_manual.md#renaming-a-parameter-parameter_aliases).

## 6. `validate_params`
- [ ] Fails fast with actionable `ValueError` for empty/invalid/mutually-exclusive inputs. This is the only guardrail before a cluster job queues.

## 7. `generate_script`
- [ ] Spine: `header()` > `activate_environment()` <-> command > `footer()`.
- [ ] Multi-phase? Split into `_generate_script_<phase>()` helpers with WHY docstrings; keep `generate_script` a short orchestrator (see `protein_mpnn.py`). Single-phase can inline (> `dssp.py`).
- [ ] Container switch, per architecture. `container_prefix()` returns `""` in env mode and the `apptainer exec … <image> ` prefix when `folders["container:<TOOL_NAME>"]` is set, so the *same* code runs either way. Wire whichever applies, even if the tool ships env-only:
    - heavy command emitted directly in the generated bash: prefix it, `f"{self.container_prefix()}python ..."` (see `dssp.py`, `rfdiffusion.py`).
    - a host helper (pipe script) dispatches the binary: pass `self.container_prefix()` into the pipe's config/flags, and in the pipe build the call as `container_argv_prefix(prefix) + cmd` so only the binary enters the `.sif` while the helper's biopipelines-importing Python stays on the host (see `gnina.py` / `pipe_gnina.py`).
- [ ] Only if the tool runs its model truly in-process (imports the model into the same interpreter as biopipelines) call `self.warn_container_unsupported()`, so a user who configures a container is told instead of silently getting env execution. This is a temporary state: such tools become containerizable once we build images that also contain biopipelines and run the whole helper inside the image.
- [ ] Serializes upstream streams to JSON (`save_json`) before invoking the pipe_script.
- [ ] If possible, serialize input parameters as well into a json inside the _configuration folder, and make the pipe scripts or the model resolve them at runtime. If this is not possible, make sure at least that every `--flag` emitted here has a matching `add_argument` in the paired pipe_script (§9).
- [ ] Missing-propagation: append `self.generate_missing_propagation(*inputs, local_missing=self.local_missing_csv, missing_csv=self.missing_csv)` (inputs can always be assumed upstream-filterable). It emits the standard block that merges upstream `missing` manifests and this tool's own local failures into `tables/missing.csv`; do not hand-thread `--upstream-missing` flags. ⚠ `local_missing=` is what carries the pipe's own rows into the merge. Drop it only if the tool can never skip an id. Alteratively, a tool that instead merges upstream itself in the pipe emits no propagation step. 
- [ ] Should support iteration of upstream stream using the declared IDs as a view over the entities. `iterate_files` already skips upstream-filtered ids (absent files); every id that enters but yields no output must land in `missing.csv`. `kind` is controlled: `"failure"` (a raise; NOT excused, flagged) vs `"filter"` (a deliberate drop; excused). The reason string goes in `cause`, not `kind`.

## 8. `get_output_files`
- [ ] Each output stream is classified as one of the three shapes: per-id files (`files=["<id>.ext"]`, consumed with `iterate_files`), shared file (`files="one.fasta"`, one artifact for all ids), or value-based (`files=[]`, datum inline in `map_table`, consumed with `iterate_values`). Content-bearing streams (sequences/compounds) are value-based; short strings/numbers per id should be value-based, not files.
- [ ] If a stream is shared-file, its declared `format` has a slicer registered in `biopipelines/stream_slicers.py` (`fasta`/`fa` and `csv` already do). Otherwise a downstream `Panda`/`Pool` filter raises `ValueError` — there is no copy-whole-file fallback. Add `@register("<fmt>")` (+ `@register_merger` if gatherable).
- [ ] Declares every stream (name, ids from input, `<id>` file template or `files=[]`, `map_table`, `format` tag) and every table (`TableInfo` with all columns).
- [ ] Includes the `missing` table (`id | removed_by | kind | cause`) if inputs can be upstream-filtered.
- [ ] Returns `"output_folder": self.output_folder`.

## 9. Registration
Two gates enforce this section, and both are mechanical. `versions/check_tool_edits.py` runs as a **pre-commit hook** and rejects the commit; `tests/test_registry_consistency.py` runs in **both CI pipelines** and fails the build. Neither is a style preference — a missing entry is a red pipeline.

⚠ The pre-commit hook only runs if it is installed. On a fresh clone: `pip install pre-commit && pre-commit install`. Without that, version-bump violations pass silently and surface later in CI.

**Version bookkeeping** (the pre-commit hook enforces the last two; `test_source_version_matches_tool_changelog` in CI enforces the first, and would otherwise never fire for this tool at all):
- [ ] Added the tool to `versions/tool_changelog.yaml` (`files:` = wrapper + pipe_script; `current:` = the `TOOL_VERSION` literal in the class body, character-for-character; a `history:` entry).
- [ ] Added a bullet naming the tool under `[unreleased]` → `### Tools` in **`versions/CHANGELOG.md`**. Note both parts: the file is `versions/CHANGELOG.md` (there is no `CHANGELOG.md` at the repo root), and the section is `[unreleased]`, not a released `[x.y.z]` one. The hook parses for exactly that heading pair.
- [ ] On any later edit to a file listed under the tool's `files:`: bump `TOOL_VERSION` in the class **and** `current` in the yaml together, and add another `[unreleased]` → `Tools` bullet.

**Registry consistency** (enforced by `tests/test_registry_consistency.py`):
- [ ] `docs/tool_index.md`: added a table row (`| <n> | <ToolName> | <Category> | <version> | ...`) **and** bumped the `**Public-API count: N**` line — the test asserts that number equals the count of `TOOL_NAME` definitions.
- [ ] The index row's version cell equals the source `TOOL_VERSION` exactly. A stale cell fails the build even though the row exists.
- [ ] `docs/tool_reference.md`: added an entry (the test greps for the tool name anywhere in the file, but write a real signature entry).
- [ ] `docs/tool/<category>.md`: added a section whose heading is `## <ToolName>` or `### <ToolName>` — a literal `##`/`###` heading starting with the tool name. Prose that merely mentions the tool does not satisfy the test and leaves the parameters undocumented.
- [ ] `README.md`: added a tool-table row using the exact markup the test greps for, `<td><sub><b><ToolName></b>` (with the name as its own `<b>` text). A row that documents two sibling tools together must be registered in the test's `README_COMBINED` map instead.
- [ ] `environments:` entries in **three** config variants — `config.cluster.yaml`, `config.colab.yaml` **and** `config.container.yaml`. The test diffs colab and container against cluster and fails on any tool present in one and absent from another. `config.daint.yaml` is exempt (aarch64: x86-64-only binaries and PyG wheels are unavailable there), so add it only if the tool actually works on aarch64.

**Not CI-enforced but still required:**
- [ ] Exported the class in `biopipelines/__init__.py` (both the `from .<module> import <ToolName>` line and the `__all__` entry). Without it `from biopipelines import <ToolName>` fails and the tool is unreachable.

## 10. The paired pipe_script(s) (`pipe_scripts/pipe_<yourtool>.py`)
See the annotated companion `_pipe_template.py` for the full pattern.
- [ ] Exists; its `add_argument` set matches the wrapper's `--flags` exactly.
- [ ] Reloads inputs from the JSON the wrapper serialized (works out-of-env via the import shim).
- [ ] Honors `--container-prefix` if it launches binaries itself.
- [ ] Writes each declared output to its declared path.
- [ ] Writes its own failures/drops to the **`--local-missing-csv`** path (`_execution/local_missing.csv`). `tables/missing.csv` is the canonical merged manifest, produced by the propagation step from upstream + this local file. The pipe never reads or writes upstream manifests. 
- [ ] Keeps computation, failure-reporting and selection separate — one artifact each. A failed id appears in all three:
    - *result table*: one row per id that ENTERED, always. An unavailable measurement is `NaN`; the id is never dropped, and required scalar fields (`id`) are never NaN. A table is a matrix of results, not a filter: dropping the row destroys the record that the id was processed and silently shortens a table a downstream `Panda`/merge joins on. 
    - *stream map*: only ids whose file was actually produced. 
    - *missing table*: why an id has no stream entity (`kind="failure"` raised / `kind="filter"` deliberately dropped, reason in `cause`).
- [ ] Does not filter a result table on a metric threshold. Emit it complete and let the user apply `Panda.filter` downstream, so selection stays an explicit pipeline step instead of a silent side effect of this tool. (`kind="filter"` is for narrowing a *stream*, not a table.)
- [ ] Exits non-zero on either step-level error: **nothing was attempted** (`iterate_files` skipped every declared id — the upstream stream is empty or its files are absent; exiting 0 would rubber-stamp a broken pipeline and hide the upstream tool that actually failed), or **every attempted id raised**. It exits 0 when ids were attempted and all were deliberately *filtered* — the tool worked, it just selected nothing. Never fail merely because the result table carries NaN.

## 11. Verify
- [ ] `python -c "import biopipelines.<yourtool>"`.
- [ ] `python -m pytest tests/ -q` is green. This is the local gate; both CI pipelines run `pytest tests/` on Python 3.10–3.13, and on GitLab a red test stage blocks the mirror to the public GitHub repo.
- [ ] `python -m pytest tests/test_registry_consistency.py -q` is green — the fast check for §9 while wiring the tool up.
- [ ] Minimal pipelines exploring different inputs and parameters produce the expected outputs.
- [ ] Installation works on the intended platform.
- [ ] *(Optional, GPU-bound.)* If the tool gets a per-tool parameter suite, it is excluded from the default run by `addopts = "-ra -m 'not tool_parameters'"` in `pyproject.toml`; run it explicitly with `pytest tests/tool_parameters -m tool_parameters -v`.

## 12. Documentation
The registry entries themselves are §9 — this section is the prose quality bar for them.
- [ ] Docstring in the class documents every parameter and every output stream/table with columns.
- [ ] The `docs/tool/<category>.md` section (added in §9) actually documents the parameters, not just the tool's existence.
- [ ] The README row (added in §9) carries a short description plus the repository and main reference-paper links, or the BP-native badge for a tool with no upstream.

## Gotchas (Python-level, not framework rules)
- ⚠ No backslash inside f-string `{...}` braces: py3.10/3.11 raise SyntaxError. Build such a fragment in a plain variable first, then interpolate it. One bad f-string breaks the whole tool's pipe script.