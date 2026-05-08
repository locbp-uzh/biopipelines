#!/usr/bin/env python3
"""Pre-commit hook: enforce TOOL_VERSION bumps for wrapper changes.

For every staged file, find which tools claim it in tool_changelog.yaml.
For each such tool, require BOTH:

  1. The tool's `current` field in tool_changelog.yaml is also staged
     (i.e. its version was bumped in this commit), and the bumped value
     matches the TOOL_VERSION literal in the wrapper module.
  2. CHANGELOG.md has a bullet mentioning the tool name in the
     [unreleased] -> Tools subsection (and CHANGELOG.md is staged).

If a wrapper-mapped file is staged without those bumps, exit 1 with a
message naming the tool and the file that triggered it.

Run order: invoked by pre-commit on every commit. Exits 0 quickly when
no wrapper files are staged.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import yaml


REPO_ROOT = Path(__file__).resolve().parent.parent
CHANGELOG = REPO_ROOT / "versions" / "CHANGELOG.md"
TOOL_CHANGELOG = REPO_ROOT / "versions" / "tool_changelog.yaml"


def staged_files() -> set[str]:
    out = subprocess.check_output(
        ["git", "diff", "--cached", "--name-only", "--diff-filter=ACMR"],
        cwd=REPO_ROOT,
        text=True,
    )
    return {line.strip().replace("\\", "/") for line in out.splitlines() if line.strip()}


def read_blob(ref: str, path: str) -> str | None:
    """Return the contents of `path` at `ref`, or None if missing."""
    try:
        return subprocess.check_output(
            ["git", "show", f"{ref}:{path}"],
            cwd=REPO_ROOT,
            text=True,
            stderr=subprocess.DEVNULL,
        )
    except subprocess.CalledProcessError:
        return None


def parse_tool_map(text: str) -> dict[str, dict]:
    return yaml.safe_load(text) or {}


def tool_version_in_wrapper(wrapper_path: Path, class_name: str) -> str | None:
    """Find TOOL_VERSION literal in the wrapper module's class body.

    Looks for ``class <class_name>(...):`` followed (within the class)
    by ``TOOL_VERSION = "x.y"``. Returns the version string or None.
    """
    if not wrapper_path.exists():
        return None
    src = wrapper_path.read_text(encoding="utf-8")
    # Find class block
    class_re = re.compile(rf"^class\s+{re.escape(class_name)}\s*\(", re.MULTILINE)
    m = class_re.search(src)
    if not m:
        return None
    # Scan from class start for the next class def (end of this class) or EOF
    rest = src[m.end():]
    end_m = re.search(r"\n^class\s+\w+\s*\(", rest, re.MULTILINE)
    block = rest[: end_m.start()] if end_m else rest
    v_m = re.search(r'^\s*TOOL_VERSION\s*=\s*"([^"]+)"', block, re.MULTILINE)
    return v_m.group(1) if v_m else None


def changelog_unreleased_tools_section(text: str) -> str:
    """Extract the body of the [unreleased] -> Tools subsection."""
    # Find [unreleased] header
    m = re.search(r"^##\s*\[unreleased\]\s*$", text, re.MULTILINE | re.IGNORECASE)
    if not m:
        return ""
    rest = text[m.end():]
    # Stop at next ## release header
    stop = re.search(r"^##\s+\[", rest, re.MULTILINE)
    unreleased = rest[: stop.start()] if stop else rest
    # Within unreleased, find ### Tools
    t = re.search(r"^###\s*Tools\s*$", unreleased, re.MULTILINE | re.IGNORECASE)
    if not t:
        return ""
    after = unreleased[t.end():]
    stop2 = re.search(r"^###\s+", after, re.MULTILINE)
    return after[: stop2.start()] if stop2 else after


def main() -> int:
    staged = staged_files()
    if not staged:
        return 0

    tool_changelog_text_now = TOOL_CHANGELOG.read_text(encoding="utf-8")
    tools_now = parse_tool_map(tool_changelog_text_now)

    # For diff: read HEAD's tool_changelog.yaml (None on first commit)
    head_text = read_blob("HEAD", "versions/tool_changelog.yaml")
    tools_head = parse_tool_map(head_text) if head_text else {}

    # Build reverse map: file -> [tools]
    file_to_tools: dict[str, list[str]] = {}
    for tool, entry in tools_now.items():
        for f in entry.get("files", []):
            file_to_tools.setdefault(f.replace("\\", "/"), []).append(tool)

    # Which tools have wrapper-mapped files staged?
    triggered: dict[str, list[str]] = {}
    for f in staged:
        for tool in file_to_tools.get(f, []):
            triggered.setdefault(tool, []).append(f)

    if not triggered:
        return 0

    # Need CHANGELOG.md staged + Tools section to mention each triggered tool
    changelog_staged = "versions/CHANGELOG.md" in staged
    changelog_text = CHANGELOG.read_text(encoding="utf-8") if CHANGELOG.exists() else ""
    tools_section = changelog_unreleased_tools_section(changelog_text)

    errors: list[str] = []
    for tool, files in sorted(triggered.items()):
        entry_now = tools_now.get(tool, {})
        entry_head = tools_head.get(tool, {})
        v_now = str(entry_now.get("current", ""))
        v_head = str(entry_head.get("current", "")) if entry_head else ""

        if v_now == v_head and v_head != "":
            errors.append(
                f"  - {tool}: TOOL_VERSION not bumped in tool_changelog.yaml "
                f"(still {v_now}). Triggered by: {', '.join(files)}"
            )
            continue

        # Wrapper module's TOOL_VERSION literal must match `current`
        wrapper_files = [
            Path(REPO_ROOT, f) for f in entry_now.get("files", [])
            if f.startswith("biopipelines/") and f.endswith(".py")
        ]
        wrapper_v = None
        for wf in wrapper_files:
            wrapper_v = tool_version_in_wrapper(wf, tool)
            if wrapper_v:
                break
        if wrapper_v is None:
            errors.append(
                f"  - {tool}: could not locate TOOL_VERSION in wrapper module "
                f"({', '.join(str(p.relative_to(REPO_ROOT)) for p in wrapper_files)})"
            )
        elif wrapper_v != v_now:
            errors.append(
                f"  - {tool}: tool_changelog.yaml has current={v_now!r} but "
                f"wrapper class has TOOL_VERSION = {wrapper_v!r}. They must match."
            )

        # CHANGELOG.md must mention the tool under [unreleased] -> Tools
        if not changelog_staged:
            errors.append(
                f"  - {tool}: CHANGELOG.md is not staged. Add a bullet under "
                f"[unreleased] -> Tools mentioning {tool}."
            )
        elif tool not in tools_section:
            errors.append(
                f"  - {tool}: no mention under [unreleased] -> Tools in "
                f"CHANGELOG.md. Add a bullet describing the change."
            )

    if errors:
        sys.stderr.write(
            "Wrapper changes require a TOOL_VERSION bump.\n"
            "Failures:\n" + "\n".join(errors) + "\n\n"
            "To fix:\n"
            "  1. Bump TOOL_VERSION in the wrapper class body.\n"
            "  2. Update `current` in tool_changelog.yaml (and append to history).\n"
            "  3. Add a bullet under [unreleased] -> Tools in CHANGELOG.md.\n"
            "  4. git add the three files and re-commit.\n"
            "Bypass (only when truly unrelated, e.g. typo in a comment): "
            "git commit --no-verify\n"
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
