"""The agent contracts must offer the MCP tools before the shell commands they replace.

An agent reads these files top to bottom and acts on the first concrete recipe it finds. When
`cluster_backend.md` (then `llm/cluster.md`) opened with an `ssh cluster "./submit ..."` table and mentioned `bp_submit`
only in a later aside, that is what it ran -- losing the operations record the tool writes,
guessing which scheduler `.out` file exists, and re-deriving a status `bp_status` reports.

So ordering is part of the contract, not editorial taste: in every file that documents both
paths, the tool comes first and the shell fallback is explicitly marked as the case where the
tool is absent. These tests fail if an edit reverses that.
"""

import io
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent

CONTRACTS = {
    # The anchor deliberately omits the ssh alias: the examples were renamed `cluster` -> `s3it`
    # and this test failed on the alias rather than on the ordering it exists to protect.
    "skills/biopipelines/references/cluster_backend.md":
        ("bp_submit", '"cd ~/biopipelines && ./submit'),
    "llm/pipelines.md": ("bp_submit", "| Submit  "),
    "skills/biopipelines/SKILL.md": ("bp_tools", "references/tool_index.md` is your catalog"),
}


def read(rel):
    return io.open(ROOT / rel, encoding="utf-8").read()


@pytest.mark.parametrize("rel,pair", sorted(CONTRACTS.items()))
def test_the_tool_is_offered_before_the_command_it_replaces(rel, pair):
    tool, fallback = pair
    text = read(rel)
    at_tool, at_fallback = text.find(tool), text.find(fallback)
    assert at_tool != -1, f"{rel} never mentions {tool}"
    assert at_fallback != -1, f"{rel} no longer contains the fallback recipe {fallback!r}"
    assert at_tool < at_fallback, (
        f"{rel} shows {fallback!r} before {tool}; an agent runs the first recipe it reads")


@pytest.mark.parametrize("rel", sorted(CONTRACTS))
def test_the_fallback_says_it_is_the_fallback(rel):
    """A recipe with no condition on it reads as the recommended way."""
    text = read(rel).lower()
    assert "without the mcp" in text or "not in your tool list" in text, (
        f"{rel} documents the manual path without saying it applies only when bp-mcp is absent")


def test_the_server_tells_the_client_it_is_the_interface():
    """The instructions travel with the tools, so they reach an agent that reads no file here."""
    pytest.importorskip("mcp", reason="bp-mcp needs the MCP SDK")
    from biopipelines import mcp_server
    text = mcp_server.INSTRUCTIONS.lower()
    assert "do not `ssh`" in text or "do not `ssh" in text
    assert "docs/tool/*.md" in text, "opening the prose docs is the other thing to steer off"


def test_the_maintainer_script_is_named_as_such():
    """`build_tool_index.py` was being called as though it were the way to look a tool up."""
    assert "maintainer" in read("skills/biopipelines/SKILL.md").lower()
