"""The vendored browser library must ship byte-identical and be packaged.

`renderers/vendor/3Dmol-min.js` is read as bytes and inlined into every page,
so the run page and `bp-visualize` work with no network. Two ways that breaks
silently: the bytes change (an end-of-line translation on a Windows checkout,
or an unrecorded upgrade), and the file does not ship at all, which only shows
up in a non-editable install -- CI installs with `-e`, so CI cannot see it.
"""
from __future__ import annotations

import hashlib
import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
VENDOR = REPO_ROOT / "biopipelines" / "renderers" / "vendor"
BLOB = VENDOR / "3Dmol-min.js"

# Recorded in renderers/vendor/README.md; the test reads it from there rather
# than duplicating it, so the note and the file cannot drift apart.
README = VENDOR / "README.md"


def _recorded_sha256() -> str:
    text = README.read_text(encoding="utf-8")
    match = re.search(r"\b([0-9a-f]{64})\b", text)
    assert match, "renderers/vendor/README.md records no sha256 for the vendored blob"
    return match.group(1)


def test_vendored_3dmol_matches_the_recorded_hash():
    """A changed blob must be a deliberate, recorded upgrade."""
    actual = hashlib.sha256(BLOB.read_bytes()).hexdigest()
    assert actual == _recorded_sha256(), (
        f"renderers/vendor/3Dmol-min.js is {actual}, README records "
        f"{_recorded_sha256()}. If this is an upgrade, update the README; if it "
        f"is an end-of-line translation, check .gitattributes."
    )


def test_vendored_blob_has_no_crlf():
    """`.gitattributes` exempts vendor/ from eol translation; prove it held."""
    assert b"\r\n" not in BLOB.read_bytes(), (
        "3Dmol-min.js contains CRLF, so a Windows checkout altered the bytes the "
        "page inlines"
    )


def test_the_license_ships_beside_the_blob():
    """Redistributing the minified bundle requires its license travel with it."""
    assert (VENDOR / "3Dmol-LICENSE.txt").is_file()
    # The webpack sidecar name is what a bundler-aware reader looks for.
    assert (VENDOR / "3Dmol-min.js.LICENSE.txt").is_file()


def test_renderers_are_declared_as_a_package():
    """Without this, a non-editable `pip install .` ships no renderers at all.

    They live inside the package and are resolved against its own directory, so
    the same `renderers/<name>.py` spelling in `config.*.yaml` works from a
    clone and from site-packages.
    """
    pyproject = (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert "[tool.setuptools.package-data]" in pyproject
    assert '"biopipelines.renderers" = ["vendor/*"]' in pyproject, (
        "the vendored blob is not package data"
    )
    assert (REPO_ROOT / "biopipelines" / "renderers" / "__init__.py").is_file(), (
        "packages.find cannot see a directory with no __init__.py"
    )


@pytest.mark.parametrize("module", [
    "pipeline_report.py", "structures.py", "streams.py", "tables.py",
    "grids.py", "images.py", "plots.py", "fasta.py",
])
def test_every_renderer_module_is_present(module):
    """The page renderer is loaded by path and raises FileNotFoundError if absent."""
    assert (REPO_ROOT / "biopipelines" / "renderers" / module).is_file()
