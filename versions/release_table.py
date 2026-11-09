#!/usr/bin/env python3
"""Table of tool version changes for the head of the [unreleased] entry of versions/CHANGELOG.md.

    python versions/release_table.py                        # print the table
    python versions/release_table.py --write                # replace it in CHANGELOG.md, between its markers
    python versions/release_table.py --label 1.6.0 --write  # at release time, once the number is chosen

The previous release's column header is read from the newest released `## [x.y.z]` heading; this release's is
"unreleased" until --label names it, so a run without flags can never stamp a stale release number.

A row is a tool whose `current` in versions/tool_changelog.yaml differs from the one at the base ref ("—" before a tool that is new).
A summary already in the table is kept, and a new row takes the first sentence of the tool's bullet under [unreleased] -> Tools, to be shortened by hand.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
CHANGELOG = REPO_ROOT / "versions" / "CHANGELOG.md"
TOOL_CHANGELOG = REPO_ROOT / "versions" / "tool_changelog.yaml"
START = "<!-- release-table:start -->"
END = "<!-- release-table:end -->"
NEW = "—"


def versions(text: str) -> dict[str, str]:
    return {tool: str(entry.get("current")) for tool, entry in (yaml.safe_load(text) or {}).items()}


def first_sentence(text: str) -> str:
    """Up to the first period that ends a sentence, one inside backticks excepted."""
    inside = False
    for i, char in enumerate(text):
        if char == "`":
            inside = not inside
        elif char == "." and not inside and text[i + 1:i + 2] in ("", " "):
            return text[:i]
    return text


def unreleased(changelog: str) -> str:
    start = changelog.index("## [unreleased]")
    end = changelog.find("\n## [", start + 1)
    return changelog[start:end if end != -1 else len(changelog)]


def latest_released(changelog: str) -> str:
    """The version of the newest released entry, which is the column the table compares against."""
    match = re.search(r"^## \[(\d[\w.]*)\]", changelog, re.MULTILINE)
    if not match:
        raise SystemExit("no released '## [x.y.z]' heading in versions/CHANGELOG.md: pass --base-label")
    return match.group(1)


def bullet_bodies(changelog: str) -> dict[str, str]:
    """{tool: text after the header} for the bullets of [unreleased] -> Tools."""
    block = unreleased(changelog)
    tools = block[block.index("### Tools"):] if "### Tools" in block else ""
    return {m.group(1): m.group(2) for m in re.finditer(r"^- \*\*(\w+)\*\* (?:new|[\d.]+) → [\d.]+\. (.*)$", tools, re.MULTILINE)}


def kept_summaries(changelog: str) -> dict[str, str]:
    """{tool: summary} of the rows already between the markers."""
    if START not in changelog:
        return {}
    table = changelog[changelog.index(START):changelog.index(END)]
    rows = (line.strip("|").split("|", 3) for line in table.splitlines() if line.startswith("| ") and not line.startswith("| Tool |"))
    return {cells[0].strip(): cells[3].strip() for cells in rows if len(cells) == 4}


def build_table(base_yaml: str, now_yaml: str, changelog: str, base_label: str, label: str) -> str:
    base, now = versions(base_yaml), versions(now_yaml)
    kept, bodies = kept_summaries(changelog), bullet_bodies(changelog)
    rows = [f"| Tool | {base_label} | {label} | summary |", "|---|---|---|---|"]
    for tool in sorted(now, key=str.lower):
        if base.get(tool) == now[tool]:
            continue
        summary = kept.get(tool) or first_sentence(bodies.get(tool, ""))
        rows.append(f"| {tool} | {base.get(tool, NEW)} | {now[tool]} | {summary} |")
    return "\n".join(rows)


def with_table(changelog: str, table: str) -> str:
    block = f"{START}\n{table}\n{END}"
    if START in changelog:
        return changelog[:changelog.index(START)] + block + changelog[changelog.index(END) + len(END):]
    heading = "## [unreleased]\n"
    at = changelog.index(heading) + len(heading)
    return f"{changelog[:at]}\n### Version changes\n\n{block}\n{changelog[at:]}"


def git_show(ref: str, path: str) -> str:
    return subprocess.check_output(["git", "show", f"{ref}:{path}"], cwd=REPO_ROOT, text=True, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--base-ref", default="origin/main", help="ref whose tool_changelog.yaml is the previous release")
    parser.add_argument("--base-label", default=None, help="column header of the previous release (default: the newest released heading in the changelog)")
    parser.add_argument("--label", default="unreleased", help="column header of this release; the release number once it is chosen")
    parser.add_argument("--write", action="store_true", help="replace the table in versions/CHANGELOG.md")
    args = parser.parse_args()

    changelog = CHANGELOG.read_bytes().decode("utf-8")
    newline = "\r\n" if "\r\n" in changelog else "\n"
    changelog = changelog.replace("\r\n", "\n")
    base_label = args.base_label or latest_released(changelog)
    table = build_table(git_show(args.base_ref, "versions/tool_changelog.yaml"),
                        TOOL_CHANGELOG.read_text(encoding="utf-8"), changelog, base_label, args.label)
    if not args.write:
        sys.stdout.reconfigure(encoding="utf-8")
        print(table)
        return 0
    CHANGELOG.write_bytes(with_table(changelog, table).replace("\n", newline).encode("utf-8"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
