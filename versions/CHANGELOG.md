# Changelog

All notable changes to biopipelines and its tool wrappers are recorded here.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/) for the
framework version (`biopipelines.__version__` in `biopipelines/__init__.py`).

Each release section has two subsections:

- **Framework** — changes to `biopipelines/`-internal infrastructure
  (`Pipeline`, `BaseConfig`, debug capture, CLI, shared `pipe_scripts/`
  helpers). Tracked by the framework version.
- **Tools** — wrapper changes that bumped a tool's `TOOL_VERSION`. Each
  bullet names the tool, the old → new version, and a one-line note.
  Mirrors the `current` field in `tool_changelog.yaml`.

The pre-commit hook (`versions/check_tool_edits.py`) refuses any
commit that modifies a wrapper or its mapped pipe scripts without
bumping the tool's `current` in `versions/tool_changelog.yaml` AND adding a
matching bullet under `[unreleased] → Tools` here.

## [unreleased]

### Framework

- `bp-config path` prints the absolute path of the active
  `config.<variant>.yaml`, so users can pipe it into their editor
  (`$EDITOR "$(bp-config path)"`).
- `bp-config edit` opens the active config in an interactive arrow-key
  TUI (`biopipelines/config_editor.py`, prompt_toolkit-based). Supports
  expand/collapse, in-line value editing, type-preserving coercion,
  search, and atomic save with a `.bak` written alongside. Round-trip
  uses `ruamel.yaml` so comments and key ordering survive.
- `prompt_toolkit>=3` and `ruamel.yaml>=0.17` are now required
  dependencies (not optional extras), so `bp-config edit` is always
  available out of the box. Both are pure-Python and small.

### Tools

- RFdiffusion3 1.0 → 2.0: cover the full RFD3 paper example set. New
  input-spec keys (`symmetry`, `unindex`, `select_unfixed_sequence`,
  `select_partially_buried`, `redesign_motif_sidechains`, `ori_token`,
  `infer_ori_strategy`, `is_non_loopy`, `plddt_enhanced`, `partial_t`) and
  inference-sampler overrides (`cfg`, `cfg_scale`, `step_scale`,
  `noise_scale`, `num_steps`, `center_option`, `seed`). Symmetry auto-selects
  the symmetry sampler. Adds symmetric-oligomer, diffused-ligand+CFG,
  protein-DNA, and atomic-motif enzyme design to the existing de-novo / binder
  / rigid-ligand paths.

## [1.1.1] — 2026-04-30

### Framework

- Add `__version__` to `biopipelines/__init__.py` and switch
  `pyproject.toml` to a dynamic version that reads from it.
- Add `TOOL_VERSION` class attribute to `BaseConfig` (default `"1.0"`)
  and explicit `TOOL_VERSION = "1.0"` on every concrete tool wrapper.
- Record `biopipelines_version` and per-tool `tool_version` in
  `_debug_capture/tools.json` (Reviewer #2 §6, B6 outstanding sub-point
  on tool versioning).
- Introduce `tool_changelog.yaml` as the source of truth for tool
  wrapper versions and the file map the pre-commit hook checks against.
- Add `scripts/check_tool_version_bumps.py` and `.pre-commit-config.yaml`
  to enforce that any wrapper change is accompanied by a `TOOL_VERSION`
  bump and a `CHANGELOG.md` entry.

### Tools

- All wrappers initialised at `TOOL_VERSION = "1.0"`.
