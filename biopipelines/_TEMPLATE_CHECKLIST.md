# Checklist for addition of a new tool

Copy the template `_TEMPLATE.py` to `biopipelines/<yourtool>.py`, then work this list top to bottom. 

## 1. Tool's shape
Before editing, classify the tool based on the following:
- **Install:** 
    - base-env (biopipelines), no install (delete `_install_script`)
    - shares another tool's env (delegate, e.g. `return PyMOL._install_script(...)`)
    - dedicated env (adapt from `dssp.py`)
    - container (adapt from `gnina.py`: source in `environments/_containers.yaml` + `cls._container_pull_block(folders, force_reinstall)` in the install). 
    The last two require entries in the config.<variant>.yaml file.
- **Inputs:** (determines output prediction)
    - single stream > check closest tool in function
    - multiple streams (combinatorial) > see Boltz2/ESMFold2/Gnina 

## 2. Identity & docstring
- [ ] `TOOL_NAME` unique; `TOOL_VERSION` set.
- [ ] Class docstring documents every input param and every output stream/table *with columns*. 

## 3. Install (`_install_script`)
- [ ] Kept / deleted / container per §1.
- [ ] If kept: gated on `_env_exists_check`; honors `force_reinstall`; verifies the binary/imports and `touch "$INSTALL_SUCCESS"` on success.
- [ ] If the installation diverges between cluster and Colab, shipped both `environments/<TOOL_NAME>.yaml` **and** `<TOOL_NAME>.colab.yaml`. The `.colab` one drops conda-only deps with no PyPI build. Fork the YAML rather than branching Python on scheduler. 

## 4. Path descriptors
- [ ] Every input JSON, output stream map, and table has a `Path(lambda self: ...)` descriptor - no path strings built inline.
- [ ] Every output stream returned in §4 has a matching `stream_map_path` descriptor here.

## 5. `__init__` (input contract)
- [ ] Normalizes each `StandardizedOutput`>`.streams.<name>` and passes `DataStream` through.
- [ ] Keeps the raw input handle (`self.structures`) for missing-propagation.
- [ ] Receives streams/standardized outputs, not bare values. Bare values (a sequence string, a PDB id, a SMILES) are turned into streams by the tools `Sequence`/`PDB`/`Ligand`/... - let the user pass those, so bare-value conversion lives in one place instead of every tool.
- [ ] `super().__init__(**kwargs)` is the last line.

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
Skipping any of these makes the commit rejected by `versions/check_tool_edits.py`:
- [ ] Added the tool to `versions/tool_changelog.yaml` (`files:` = wrapper + pipe_script; `current:` = the `TOOL_VERSION` literal; a `history:` entry).
- [ ] Added a `CHANGELOG.md` bullet under `[Version] > Tools` naming the tool.
- [ ] On any later edit to the files: bump `TOOL_VERSION` and `current` together.
- [ ] Exported the class in `biopipelines/__init__.py`.

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
- [ ] Minimal pipelines exploring different inputs and parameters produce the expected outputs.
- [ ] Installation works on the intended platform.

## 12. Documentation
- [ ] Docstring in the class.
- [ ] Entries under docs/tool_index.md, docs/tool_reference.md, docs/tool/<category>.md
- [ ] Short description and link to repository and main reference paper in the README.

## Gotchas (Python-level, not framework rules)
- ⚠ No backslash inside f-string `{...}` braces: py3.10/3.11 raise SyntaxError. Build such a fragment in a plain variable first, then interpolate it. One bad f-string breaks the whole tool's pipe script.