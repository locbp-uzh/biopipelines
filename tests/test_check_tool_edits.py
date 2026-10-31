"""The pre-commit hook asks for a bullet headed by the tool, not for the tool's name anywhere in the Tools section.

Nine tool bullets sat under `### Framework` unnoticed because "Panda" occurs in the body of many Tools bullets, and the hook only looked for the substring.
"""

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="module")
def hook():
    pytest.importorskip("yaml")
    spec = importlib.util.spec_from_file_location("check_tool_edits", ROOT / "versions" / "check_tool_edits.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CHANGELOG = """# Changelog

## [unreleased]

### Framework

- `Panda` 1.7 → 1.8 - a pooled stream keeps its extra columns.

### Tools

- **Gnina** 2.8 → 2.9. A receptor with bad coordinates is refused, and `Panda.filter` on the table works again.

- **BoltzGenImport** 2.4 → 2.5. No behavioral change.

- **MMseqs2** 1.9 → 1.10, **MSA** 1.3 → 1.4. Both pass their sequence input on.

## [1.5.3]

- **Panda** 1.0 → 1.1. An older release.
"""


def test_a_tool_named_only_in_the_body_of_another_bullet_has_no_bullet(hook):
    section = hook.changelog_unreleased_tools_section(CHANGELOG)
    assert not hook.section_names_tool(section, "Panda")


def test_a_bullet_headed_by_the_tool_counts(hook):
    section = hook.changelog_unreleased_tools_section(CHANGELOG)
    assert hook.section_names_tool(section, "Gnina")


def test_every_tool_of_a_multi_tool_header_counts(hook):
    section = hook.changelog_unreleased_tools_section(CHANGELOG)
    assert hook.section_names_tool(section, "MMseqs2") and hook.section_names_tool(section, "MSA")


def test_a_tool_whose_name_starts_another_is_not_named_by_it(hook):
    section = hook.changelog_unreleased_tools_section(CHANGELOG)
    assert not hook.section_names_tool(section, "BoltzGen")
    assert not hook.section_names_tool(section, "MMseqs2Server")


def test_the_framework_section_and_released_entries_are_not_the_tools_section(hook):
    section = hook.changelog_unreleased_tools_section(CHANGELOG)
    assert "pooled stream" not in section and "An older release" not in section
