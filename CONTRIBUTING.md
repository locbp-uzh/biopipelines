# Contributing to BioPipelines

Thank you for your interest in contributing! Below are the guidelines for different types of contributors.

## Who can contribute

- **Core authors** (`@GQChem`, `@GianlucaQuargnali`, `@priveraf`): merge rights to `main` after self-review or peer review.
- **Lab members**: work in the private internal repository (`biopipelines-locbp`) and open PRs to `main` from there when a feature is ready.
- **External contributors**: fork this repository, develop on your fork, and open a PR to `main`.

## Workflow for external contributors

1. **Fork** this repository on GitHub.
2. **Clone** your fork locally and create a feature branch:
   ```bash
   git checkout -b my-feature
   ```
3. Make your changes. Keep commits focused and descriptive.
4. **Push** to your fork and open a **Pull Request against `main`**.
5. Fill in the PR template. A core author will review and may request changes before merging.

## Branch protection rules on `main`

- Direct pushes to `main` are not allowed.
- Every PR requires at least **1 approval from a core author** (enforced via CODEOWNERS).
- PRs must be up to date with `main` before merging.

## What we welcome

- New tool integrations (structural biology, ML models, docking, etc.)
- Enhancements to existing tools and helper scripts
- Bug fixes and robustness improvements
- Documentation improvements and new example pipelines

## What to avoid

- Do not include private project files, internal data, or unpublished sequences/structures. A file that must never be published carries a line beginning `#PRIVATE` in its header — a valid comment as-is in Python, YAML and shell. `python versions/check_private.py` lists any such file that a publish would carry into `public`, and exits non-zero if there is one; it walks the whole commit range, because merging publishes full ancestry and a file deleted later is still readable in history. It marks whole files only, so plumbing that merely names a private tool (configs, fixtures, generated indexes) still needs reading by eye.
- Do not add dependencies without discussion — open an issue first for significant new dependencies.
- Keep PRs focused. One feature or fix per PR is easier to review and merge.

## Code style

- Follow existing conventions in the codebase (class-based tools, `Pipeline` context manager pattern).
- Add docstrings to new public classes and functions.
- Test your pipeline with at least a minimal example before opening a PR.

## Questions?

Open a [GitHub Issue](../../issues) for bugs or feature requests, or reach out to the core authors directly.
