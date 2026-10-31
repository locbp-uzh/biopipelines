"""The table of version changes and the Tools bullets of [unreleased] agree with tool_changelog.yaml.

versions/release_table.py writes the table; a release entry is ready when every row, every bullet and the yaml name the same version of the same tools.
"""

import importlib.util
import re
from pathlib import Path

import pytest

yaml = pytest.importorskip("yaml")

ROOT = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="module")
def table():
    spec = importlib.util.spec_from_file_location("release_table", ROOT / "versions" / "release_table.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


BASE = 'Panda:\n  current: "1.5"\nGnina:\n  current: "2.2"\nDNAEncoder:\n  current: "1.2"\nVina:\n  current: "2.2"\n'
NOW = BASE.replace('"1.5"', '"1.8"').replace('"1.2"', '"1.3"') + 'FLOWR:\n  current: "1.3"\n'
CHANGELOG = """## [unreleased]

### Version changes

<!-- release-table:start -->
| Tool | 1.5.3 | 1.6.0 | summary |
|---|---|---|---|
| Panda | 1.5 | 1.8 | Kept by hand |
<!-- release-table:end -->

### Tools

- **FLOWR** new → 1.3. Pocket-conditioned 3D ligand generation. More text.

- **DNAEncoder** 1.2 → 1.3. An id dropped on purpose (`a.b` here) no longer aborts the step. More.

- **Panda** 1.5 → 1.8. First sentence of the bullet.

## [1.5.3]
"""


def test_changelog_table_lists_the_moved_tools_in_alphabetical_order(table):
    rows = table.build_table(BASE, NOW, CHANGELOG, "1.5.3", "1.6.0").splitlines()[2:]
    assert [row.split(" | ")[0].strip("| ") for row in rows] == ["DNAEncoder", "FLOWR", "Panda"]


def test_changelog_table_marks_a_new_tool_and_omits_an_unmoved_one(table):
    text = table.build_table(BASE, NOW, CHANGELOG, "1.5.3", "1.6.0")
    assert "| FLOWR | — | 1.3 |" in text and "| DNAEncoder | 1.2 | 1.3 |" in text
    assert "Gnina" not in text and "Vina" not in text


def test_changelog_table_keeps_a_summary_and_fills_a_new_row_from_the_bullet(table):
    text = table.build_table(BASE, NOW, CHANGELOG, "1.5.3", "1.6.0")
    assert "| Panda | 1.5 | 1.8 | Kept by hand |" in text
    assert "| FLOWR | — | 1.3 | Pocket-conditioned 3D ligand generation |" in text


def test_changelog_first_sentence_skips_a_period_inside_backticks(table):
    assert table.first_sentence("An id dropped on purpose (`a.b` here) no longer aborts the step. More.") == "An id dropped on purpose (`a.b` here) no longer aborts the step"
    assert table.first_sentence("Version 1.9 moves. More") == "Version 1.9 moves"


def test_changelog_table_replaces_the_old_one_and_inserts_one_when_absent(table):
    replaced = table.with_table(CHANGELOG, "| Tool |\n|---|\n| New |")
    assert "| New |" in replaced and "Kept by hand" not in replaced and replaced.endswith("## [1.5.3]\n")
    bare = "# Changelog\n\n## [unreleased]\n\n### Tools\n\n- **X** 1.0 → 1.1. Y.\n"
    inserted = table.with_table(bare, "| T |")
    assert inserted.index("### Version changes") < inserted.index("<!-- release-table:start -->") < inserted.index("### Tools")


def _release_text():
    text = (ROOT / "versions" / "CHANGELOG.md").read_text(encoding="utf-8").replace("\r\n", "\n")
    return text, text[text.index("## [unreleased]"):text.find("\n## [", text.index("## [unreleased]"))]


def test_changelog_table_rows_bullets_and_yaml_agree():
    """Skipped while no release table is being prepared (between a release and the next grouping)."""
    text, block = _release_text()
    if "<!-- release-table:start -->" not in block:
        pytest.skip("no table of version changes in [unreleased]")
    current = {tool: str(entry["current"]) for tool, entry in yaml.safe_load((ROOT / "versions" / "tool_changelog.yaml").read_text(encoding="utf-8")).items()}
    rows = dict(re.findall(r"^\| (\w+) \| [^|]+ \| ([^|]+) \|", block[block.index("<!-- release-table:start -->"):block.index("<!-- release-table:end -->")], re.M))
    rows.pop("Tool")
    tools = block[block.index("### Tools"):]
    bullets = re.findall(r"^- \*\*(\w+)\*\* (?:new|\d+\.\d+) → (\d+\.\d+)\. ", tools, re.M)
    problems = [f"{tool}: row {rows[tool].strip()}, yaml {current.get(tool)}" for tool in rows if rows[tool].strip() != current.get(tool)]
    problems += [f"{tool}: bullet ends at {to}, yaml {current.get(tool)}" for tool, to in bullets if to != current.get(tool)]
    problems += [f"{tool}: {n} bullets" for tool in rows if (n := sum(1 for t, _ in bullets if t == tool)) != 1]
    problems += [f"{tool}: bullet but no row" for tool, _ in bullets if tool not in rows]
    assert not problems, "the version table, the Tools bullets and tool_changelog.yaml disagree:\n  " + "\n  ".join(problems)


def test_changelog_has_no_tool_bullet_under_framework():
    _, block = _release_text()
    framework = block[block.index("### Framework"):block.index("### Tools")]
    stray = re.findall(r"^- (?:\*\*`?\w+`?\*\*|`\w+`) \d+\.\d+ → \d+\.\d+", framework, re.M)
    assert not stray, f"tool bullets under Framework: {stray}"
