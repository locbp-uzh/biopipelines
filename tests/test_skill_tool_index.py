"""The skill's agent-facing tool index must not drift from the human docs it is generated from.

`skills/biopipelines/references/tool_index.md` is what an agent reads to find out which tools
exist; `docs/tool/*.md` is what a person reads. The first is generated from the second, so an
edit to the prose docs that is not followed by a regeneration leaves agents working from a
stale catalog — a tool they cannot see, or a signature pointer into a section that moved.

Regenerate with `python skills/biopipelines/build_tool_index.py`.
"""

import io
import pathlib
import re

from biopipelines import tool_docs

ROOT = pathlib.Path(__file__).resolve().parent.parent
INDEX = ROOT / "skills" / "biopipelines" / "references" / "tool_index.md"

# Public classes that deliberately share a parent's TOOL_NAME, so they never appear in the
# registry: SolubleMPNN is ProteinMPNN with the soluble model locked on, LoadMultiple is Load's
# many-inputs form. Both are exported and callable, so the agent index must list them.
SHARED_IDENTITY = {"SolubleMPNN", "LoadMultiple"}


def test_index_is_regenerated():
    """Byte-for-byte, so a doc edit cannot land without refreshing what agents read."""
    expected = tool_docs.render(tool_docs.collect())
    actual = io.open(INDEX, encoding="utf-8").read()
    assert actual == expected, (
        "skills/biopipelines/references/tool_index.md is stale -- run "
        "`python skills/biopipelines/build_tool_index.py`")


def test_every_registry_tool_is_in_the_agent_index():
    registry = set(re.findall(r"^\|\s*\d+\s*\|\s*([A-Za-z0-9_]+)\s*\|",
                              io.open(ROOT / "docs" / "tool_index.md", encoding="utf-8").read(),
                              re.M))
    listed = set(re.findall(r"^- \*\*([A-Za-z0-9_]+)\*\*",
                            io.open(INDEX, encoding="utf-8").read(), re.M))
    assert not registry - listed, f"absent from the agent index: {sorted(registry - listed)}"
    assert listed - registry == SHARED_IDENTITY, (
        f"unexpected entries with no registry row: {sorted(listed - registry - SHARED_IDENTITY)}")


def test_every_pointer_resolves():
    """A `docs/tool/file.md#anchor` the agent cannot follow is worse than no pointer."""
    broken = []
    for entry in tool_docs.collect():
        path, _, anchor = entry["path"].partition("#")
        text = io.open(ROOT / path, encoding="utf-8").read()
        headings = {tool_docs.slug(h) for h in re.findall(r"^#{2,3}\s+(.+?)\s*$", text, re.M)}
        headings |= set(re.findall(r"\{#([^}]+)\}", text))
        if anchor not in headings:
            broken.append(f"{entry['name']} -> {entry['path']}")
    assert not broken, "pointers into docs/tool/ that resolve to nothing:\n  " + "\n  ".join(broken)


def test_every_entry_carries_hardware_and_platforms():
    """The index exists so an agent never has to parse the README's badge markup.

    A tool with no hardware or no verified platform means the README table changed shape and
    the lift silently produced nothing -- the failure mode is an agent picking a GPU tool for
    a CPU-only run, which surfaces only at submission.
    """
    bare = [f"{e['name']} (hw={e['hardware']!r}, platforms={e['platforms']})"
            for e in tool_docs.collect() if not e["hardware"] or not e["platforms"]]
    assert not bare, "entries with no badges lifted from the README:\n  " + "\n  ".join(bare)


def test_index_stays_cheap_to_read():
    """The whole point is that an agent can afford this instead of the 78k-token prose docs."""
    size = INDEX.stat().st_size
    assert size < 40_000, f"the agent index has grown to {size} bytes (~{size // 4} tokens)"
