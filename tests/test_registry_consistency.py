"""Every public tool must be registered everywhere a reader would look for it.

These gaps are all mechanically detectable, and every one of them was found by
hand in an audit rather than by CI:

  * BindingData and LigandAtomSelector shipped with no README row, and
    LigandAtomSelector had no prose docs at all.
  * Scripting had no tool_changelog.yaml entry, so the pre-commit version gate
    could never fire for it.
  * OpenMM's ligand= path was mapped to an env carrying its dependencies on one
    config variant only, and failed with ModuleNotFoundError everywhere else --
    once per structure, inside the per-item except, rather than as one error.
  * MMseqs2 and MMseqs2Server were absent from config.daint.yaml's environment
    block, whose own comment claims every base-env tool is listed.

Exemptions are explicit and named, so adding one is a decision rather than an
oversight.
"""

import io
import pathlib
import re

import pytest

yaml = pytest.importorskip("yaml")

ROOT = pathlib.Path(__file__).resolve().parent.parent

# Not part of the advertised public API (mirrors docs/tool_index.md's own list).
INTERNAL = {"BoltzGenMerge", "BoltzGenImport", "RFDAA_PrepareLigand", "Mock",
            "TemplateTool", "base", "install"}

# Documented in the README under a combined row with a sibling tool.
README_COMBINED = {"MMseqs2": "MMseqs2 / MMseqs2Server",
                   "MMseqs2Server": "MMseqs2 / MMseqs2Server",
                   "Load": "Load / LoadMultiple"}

CONFIGS = ["config.cluster.yaml", "config.colab.yaml",
           "config.daint.yaml", "config.container.yaml"]

# Tools deliberately unavailable on a variant, with the reason. Daint is aarch64,
# so anything shipping only x86-64 binaries or PyG wheels cannot run there.
ENV_EXEMPT = {
    "config.daint.yaml": "aarch64: x86-64-only binaries and PyG wheels are unavailable",
}

_PAIR = re.compile(
    r'TOOL_NAME\s*=\s*["\']([^"\']+)["\'](?:(?!TOOL_NAME).)*?TOOL_VERSION\s*=\s*["\']([^"\']+)["\']',
    re.S)


def _read(rel):
    return io.open(ROOT / rel, encoding="utf-8").read()


@pytest.fixture(scope="module")
def tools():
    """{TOOL_NAME: (version, source file)} for every public tool."""
    found = {}
    for path in sorted((ROOT / "biopipelines").glob("*.py")):
        for m in _PAIR.finditer(io.open(path, encoding="utf-8").read()):
            if m.group(1) not in INTERNAL:
                found[m.group(1)] = (m.group(2), path.name)
    assert len(found) > 50, f"only found {len(found)} tools; the scan is broken"
    return found


@pytest.fixture(scope="module")
def changelog():
    data = yaml.safe_load(_read("versions/tool_changelog.yaml"))
    return data.get("tools", data)


def test_source_version_matches_tool_changelog(tools, changelog):
    """The gate the pre-commit hook enforces, checked for the whole tree.

    Four tools once carried a version bump copy-pasted from an unrelated tool,
    with that tool's note attached; nothing noticed until an audit.
    """
    drift = []
    for name, (version, src) in sorted(tools.items()):
        entry = changelog.get(name)
        if entry is None:
            drift.append(f"{name} ({src}): absent from tool_changelog.yaml")
        elif str(entry.get("current")) != str(version):
            drift.append(f"{name} ({src}): source {version}, yaml {entry.get('current')}")
    assert not drift, "tool_changelog.yaml disagrees with the sources:\n  " + "\n  ".join(drift)


def test_every_tool_is_in_the_index(tools):
    index = _read("docs/tool_index.md")
    missing = [n for n in sorted(tools) if not re.search(rf'\|\s*\d+\s*\|\s*{re.escape(n)}\s', index)]
    assert not missing, f"absent from docs/tool_index.md: {missing}"


def test_index_versions_match_the_sources(tools):
    index = _read("docs/tool_index.md")
    drift = []
    for name, (version, _src) in sorted(tools.items()):
        m = re.search(rf'\|\s*\d+\s*\|\s*{re.escape(name)}\s+\|[^|]+\|\s*([^\s|]+)\s*\|', index)
        if m and m.group(1) != version:
            drift.append(f"{name}: source {version}, index {m.group(1)}")
    assert not drift, "docs/tool_index.md is stale:\n  " + "\n  ".join(drift)


def test_index_count_matches_the_table(tools):
    index = _read("docs/tool_index.md")
    stated = int(re.search(r'\*\*Public-API count: (\d+)\*\*', index).group(1))
    assert stated == len(tools), (
        f"the index says {stated} tools; {len(tools)} define a TOOL_NAME")


def test_every_tool_has_a_reference_entry(tools):
    ref = _read("docs/tool_reference.md")
    missing = [n for n in sorted(tools) if n not in ref]
    assert not missing, (
        f"absent from docs/tool_reference.md, so unreachable from the docs index: {missing}")


def test_every_tool_has_prose_docs(tools):
    prose = "\n".join(_read(f"docs/tool/{p.name}")
                      for p in sorted((ROOT / "docs" / "tool").glob("*.md")))
    missing = [n for n in sorted(tools)
               if not re.search(rf'^#{{2,3}}\s+{re.escape(n)}\b', prose, re.M)]
    assert not missing, (
        f"no section in docs/tool/*.md, so the parameters are undocumented: {missing}")


def test_every_tool_has_a_readme_row(tools):
    readme = _read("README.md")
    rows = set(re.findall(r'<td><sub><b>([^<]+)</b>', readme))
    missing = []
    for name in sorted(tools):
        if name in rows or README_COMBINED.get(name) in rows:
            continue
        missing.append(name)
    assert not missing, f"no row in the README tool table: {missing}"


def test_env_mappings_agree_across_variants():
    """A tool mapped on one variant and absent on another is usually a bug.

    OpenMM's ligand stack lived on Daint alone for exactly this reason. Daint is
    exempt because it is aarch64 and genuinely cannot run some tools.
    """
    mapped = {}
    for config in CONFIGS:
        cfg = yaml.safe_load(_read(config))
        mapped[config] = set(((cfg or {}).get("environments", {}) or {}))

    reference = "config.cluster.yaml"
    gaps = []
    for config in CONFIGS:
        if config == reference or config in ENV_EXEMPT:
            continue
        for tool in sorted(mapped[reference] - mapped[config]):
            gaps.append(f"{tool}: in {reference}, absent from {config}")
    assert not gaps, (
        "environment mappings diverge between variants:\n  " + "\n  ".join(gaps))
