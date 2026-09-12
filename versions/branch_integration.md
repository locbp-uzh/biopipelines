# Branch integration log

Record of merges between the long-lived branches of this repository, what conflicted, and how each conflict was resolved. One section per integration, newest first.

Branch roles:

- `main` — GitLab integration branch (`gitlab.uzh.ch/locbp/public/biopipelines-locbp`). Anyone in the lab can open an MR or push to it.
- `gianluca-main` — Gianluca's development branch.
- `pablo-main` — Pablo's development branch.
- `public` — mirrored to the public GitHub organization (`github.com/orgs/locbp-uzh`). Managed by pablo.riverafuentes and gianluca.quargnali only.

---

## 2026-08-31 — `main` → `gianluca-main`

Merge commit `3d2a772`. Parents: `b4eb4ba` (gianluca-main), `4bfbea7` (main). Pre-merge state preserved as tags `backup/gianluca-main-premerge` and `backup/main-premerge`.

### Divergence at the time of the merge

Merge base `bbc3152` (2026-07-28, "Resolve range-pattern ids from the map table at runtime").

| Direction | Commits | Files | Lines |
| --------- | ------- | ----- | ----- |
| `main` ahead of `gianluca-main` | 28 | 39 | +3462 / −36 |
| `gianluca-main` ahead of `main` | 66 | 73 | +5060 / −391 |

Only 13 files were touched by both sides, and only one of those was a code file (`biopipelines/esmfold2.py`). The rest were shared registries, docs and changelogs.

### What `main` contributed

| Item | Nature |
| ---- | ------ |
| `StructureCluster` 1.0 → 1.1 | New Analysis tool. Fold-similarity clustering of a structures stream by sphere exclusion over length-normalized C-alpha traces; 1.1 adds per-cluster medoids. |
| `CatalyticScoring` 1.0 | New tool. Per-protein catalytic metrics and oxygen-displacement consistency from a GNINA `conformer_ranking` table. |
| `LigandAtomSelector` 1.0 | New tool. `DistanceSelector` referenced to an atom-name subset inside one named ligand residue, for chimeric ligands. |
| `ESMFold2` 1.0 → 1.4 | `modifications`, `covalent_bonds`, `glycans`; entity-ordered chain ids; `format`-gated CCD handling; the missing-table propagation fix. |
| `Boltz2` 1.0 → 1.1 | `pipe_boltz_config_unified.py` now declares and honors `--single-sequence`. |
| `RDKit` strain | Minimizer iteration cap exposed as `max_iters`. |
| Daint tooling | `envs/daint/` container and venv provisioning, `*.daint.yaml` env specs, `esmfold2.pip.daint.txt`. |

### Conflict resolutions

**`biopipelines/esmfold2.py`** — both branches independently fixed the same bug: `get_output_files()` declares the `missing` table whenever an input carries an upstream missing manifest, but `generate_script` never emitted the matching propagation step, so `tables/missing.csv` was never written and the completion check stamped `_FAILED` on otherwise-correct steps.

Kept `main`'s form, which guards the emit on `_collect_upstream_missing_paths(*self._missing_input_sources())` so the declaration and the emit cannot drift apart. This is also the established pattern elsewhere in the package (`boltz2.py`, `consensus.py`, `distance_selector.py`, `ensemble_analysis.py`, `openmm.py`, `pdb.py`). The unconditional call from `gianluca-main` was dropped. The accompanying eight-line comment was collapsed to one line per `CLAUDE.md`.

**`versions/CHANGELOG.md`** — two pure insertion blocks at the same anchor in `[unreleased] → Tools`. Resolved as a union, `main`'s six entries above `gianluca-main`'s five.

**`docs/tool_index.md`** — the only conflict needing real reconciliation. `main` had inserted `StructureCluster` as a placeholder row numbered `40b` and left the header count at 81; `gianluca-main` had renumbered the whole table to 83. Resolution:

- Full renumber, placeholder row removed. Count 81 / 83 → **86**.
- Registered `CatalyticScoring` and `LigandAtomSelector`, which `main` added to `biopipelines/__init__.py` without ever indexing. Both classified Analysis: `CatalyticScoring` consumes GNINA tables, `LigandAtomSelector` mirrors `DistanceSelector`.
- Moved `BindingData` from the end of the table into the alphabetical Analysis block, where every other Analysis entry sits.
- Re-read every `Version` cell from the `TOOL_VERSION` literal in the corresponding source module, the file's own stated source of truth. This corrected eight stale entries: `LigandMPNN` 1.1 → 1.2, `Boltz2` 1.0 → 1.1, `ESMFold2` 1.0 → 1.4, `OpenMM` 1.1 → 1.2, `AiZynthFinder` 1.1 → 1.2, `MMseqs2` 1.0 → 1.1, `Ligand` 1.0 → 1.1, `RCSB` 1.2 → 1.3.

No tool defining a `TOOL_NAME` in `biopipelines/*.py` is now absent from the index, apart from the documented exclusions.

### Verification

`python -m pytest tests/ -q` → **826 passed, 2 skipped**. Both skips are pre-existing and environmental: `test_internal_and_folder.py:252` and `test_pipeline_generation.py:166` run `bash` on a generated `.sh` whose embedded Windows paths MSYS/Git-bash mangles; they run on POSIX CI. `import biopipelines` succeeds with 123 names in `__all__`.

### Known gaps left open

1. **`CatalyticScoring` and `LigandAtomSelector` are undocumented.** `main` added both tools and registered them in `biopipelines/__init__.py`, but neither appears in `docs/tool/*.md` nor in `docs/tool_reference.md`. They are now in `docs/tool_index.md`; the prose entries still need writing, which requires knowledge of the tools' intended use that the source alone does not supply. Both must be documented before `gianluca-main` goes anywhere near `public`.

2. **Five pre-existing `tool_changelog.yaml` mismatches on `gianluca-main`.** Present at `b4eb4ba`, absent from `main` at `4bfbea7`, and therefore not introduced by this merge. `versions/check_tool_edits.py` enforces agreement between a tool's `current` field and its `TOOL_VERSION` literal, so these will block any commit that stages the affected wrapper:

   | Tool | `tool_changelog.yaml` `current` | `TOOL_VERSION` in source |
   | ---- | ------------------------------- | ------------------------ |
   | `OpenMM` | 1.1 | 1.2 |
   | `ProteinMPNN` | 1.2 | 1.1 |
   | `RFDAA_PrepareLigand` | 1.2 | 1.1 |
   | `RFdiffusion` | 1.2 | 1.1 |
   | `RFdiffusionAllAtom` | 1.2 | 1.1 |

   `OpenMM` looks like a source bump whose changelog entry was never written. The other four look like changelog bumps whose source changes were reverted or never landed. Each needs deciding individually; do not resolve them by blanket-syncing one file to the other.

3. **`gianluca-main` has not been pushed.** The merge and the folder move below are local only.

---

## 2026-08-31 — `envs/daint` moved under `environments/`

Commit `883733a`. Not a branch integration, but recorded here because it rearranges what the merge above brought in.

A top-level `envs/` containing nothing but `daint/` was redundant next to `environments/`, which already holds every environment spec and an `environments/scripts/` for the generic container build. The six Daint scripts moved to `environments/daint/`, kept as a sibling of `scripts/` rather than merged into it because both directories carry a file named `Containerfile`.

Usage comments in all six files and the one real path (`build_base.sh`'s `HERE=`) were updated. Nothing outside the moved directory referenced the old location — verified by grep across `*.py`, `*.md`, `*.sh`, `*.yaml`, `*.yml`, `*.txt` and `Containerfile`, and `.gitignore` has no pattern that would swallow the new path.
