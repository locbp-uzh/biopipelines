#!/usr/bin/env python3
"""Cherry-pick a commit from another branch and re-register its tools in the shared bookkeeping files.

Every tool commit edits the same files (CHANGELOG, tool_changelog.yaml, the tool indexes), so cherry-picking
between branches conflicts there every time. This resolves those conflicts mechanically and leaves real code
conflicts to you.

    python versions/integrate_commit.py <sha>                     # tools are detected from the commit
    python versions/integrate_commit.py <sha> --tool Foo:1.0      # force a tool and its version on this branch
    python versions/integrate_commit.py <sha> --no-commit         # stop before committing

For each tool the commit registers or bumps, in versions/tool_changelog.yaml: the new block (or the newest history
entry, renumbered to continue this branch's version); in versions/CHANGELOG.md: the commit's added bullets under
[unreleased] -> Tools; in docs/tool_index.md: the row, renumbered, with versions synced to the sources and the count
updated; then the agent-facing index is regenerated. A conflict in any other file stops the script with the
cherry-pick left open for you to resolve; rerun with --continue to finish the bookkeeping afterwards.
"""
import argparse
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GENERATED = ("docs/tool_index.md", "skills/biopipelines/references/tool_index.md",
             "versions/CHANGELOG.md", "versions/tool_changelog.yaml")
INTERNAL = {"BoltzGenMerge", "BoltzGenImport", "RFDAA_PrepareLigand", "Mock", "TemplateTool", "base", "install"}
PAIR = re.compile(r'TOOL_NAME\s*=\s*["\']([^"\']+)["\'](?:(?!TOOL_NAME).)*?TOOL_VERSION\s*=\s*["\']([^"\']+)["\']', re.S)
ROW = re.compile(r'^\|\s*(\d+)\s*\|\s*(\S+)\s*\|\s*([^|]+?)\s*\|\s*([^|]+?)\s*\|.*$', re.M)


def git(*args, check=False):
    out = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, encoding="utf-8")
    if check and out.returncode:
        sys.exit(out.stderr.strip())
    return out.stdout


def read(path):
    raw = (ROOT / path).read_bytes().decode("utf-8")
    return raw.replace("\r\n", "\n"), "\r\n" if "\r\n" in raw else "\n"


def write(path, text, nl):
    (ROOT / path).write_bytes(text.replace("\n", nl).encode("utf-8"))


def blocks(text):
    return {m.group(1): m.group(0) for m in re.finditer(r"^(\w[\w.-]*):\n(?:[ ].*\n|\n)*", text, re.M)}


def current(block):
    return re.search(r'  current: "([^"]+)"', block).group(1)


def bump_minor(version):
    major, minor = version.split(".")[:2]
    return f"{major}.{int(minor) + 1}"


def detect_tools(sha):
    """{tool: version on this branch} for each tool whose yaml entry the commit adds or advances."""
    theirs = blocks(git("show", f"{sha}:versions/tool_changelog.yaml"))
    parent = blocks(git("show", f"{sha}~1:versions/tool_changelog.yaml"))
    mine, _ = read("versions/tool_changelog.yaml")
    mine = blocks(mine)
    found = {}
    for tool, block in theirs.items():
        if tool in parent and current(parent[tool]) == current(block):
            continue
        found[tool] = current(block) if tool not in mine else bump_minor(current(mine[tool]))
    return found


def register_yaml(sha, tool, version):
    theirs = blocks(git("show", f"{sha}:versions/tool_changelog.yaml"))[tool]
    text, nl = read("versions/tool_changelog.yaml")
    mine = blocks(text).get(tool)
    if mine is None:
        text = text.rstrip("\n") + "\n\n" + theirs.rstrip("\n") + "\n"
    else:
        note = re.search(r"  history:\n  - version: '[^']+'\n    date: '([^']+)'\n    note: >-\n((?:      .*\n)+)", theirs)
        entry = f"  - version: '{version}'\n    date: '{note.group(1)}'\n    note: >-\n{note.group(2)}"
        new = re.sub(r'  current: "[^"]+"', f'  current: "{version}"', mine, count=1)
        empty = "  history: []\n"
        new = new.replace(empty, "  history:\n" + entry, 1) if empty in new else new.replace("  history:\n", "  history:\n" + entry, 1)
        text = text.replace(mine, new)
    write("versions/tool_changelog.yaml", text, nl)


def register_changelog(sha, renames):
    diff = git("show", sha, "--", "versions/CHANGELOG.md")
    added = [ln[1:] for ln in diff.splitlines() if ln.startswith("+") and not ln.startswith("+++") and ln[1:].strip()]
    bullets = [b for b in added if b.startswith("- ")]
    for tool, (old, new) in renames.items():
        bullets = [re.sub(rf"({re.escape(tool)}\W+){re.escape(old[0])} → {re.escape(old[1])}", rf"\g<1>{new[0]} → {new[1]}", b) for b in bullets]
    text, nl = read("versions/CHANGELOG.md")
    head = text.index("## [unreleased]")
    for section, wanted in (("### Tools\n\n", [b for b in bullets if b.startswith("- **")]),
                            ("### Framework\n\n", [b for b in bullets if not b.startswith("- **")])):
        anchor = text.index(section, head) + len(section)
        text = text[:anchor] + "".join(b + "\n\n" for b in wanted) + text[anchor:]
    write("versions/CHANGELOG.md", text, nl)


def sync_index(sha, tools):
    versions = {n: v for f in (ROOT / "biopipelines").glob("*.py") for n, v in PAIR.findall(f.read_text(encoding="utf-8"))}
    text, nl = read("docs/tool_index.md")
    theirs = git("show", f"{sha}:docs/tool_index.md")
    theirs_names = [m.group(2) for m in ROW.finditer(theirs)]
    for tool in tools:
        if re.search(rf"^\|\s*\d+\s*\|\s*{re.escape(tool)}\s*\|", text, re.M):
            continue
        row = next((m.group(0) for m in ROW.finditer(theirs) if m.group(2) == tool), None)
        if row is None:
            continue
        before = next((n for n in reversed(theirs_names[:theirs_names.index(tool)]) if re.search(rf"^\|\s*\d+\s*\|\s*{re.escape(n)}\s*\|", text, re.M)), None)
        at = re.search(rf"^\|\s*\d+\s*\|\s*{re.escape(before)}\s*\|.*$", text, re.M).end() if before else ROW.search(text).start()
        text = text[:at] + "\n" + row.rstrip("\r") + text[at:]
    count = 0
    lines = []
    for ln in text.split("\n"):
        m = ROW.match(ln)
        if m:
            count += 1
            if m.group(2) in versions:
                ln = re.sub(r"^(\|\s*\d+\s*\|\s*\S+\s*\|\s*[^|]+\|\s*)[\d.]+", lambda x: x.group(1) + versions[m.group(2)], ln, count=1)
            ln = re.sub(r"^\|\s*\d+\s*\|", f"| {count} |", ln, count=1)
        lines.append(ln)
    text = re.sub(r"Public-API count: \d+", f"Public-API count: {count}", "\n".join(lines))
    write("docs/tool_index.md", text, nl)
    missing = sorted(set(versions) - INTERNAL - {m.group(2) for m in ROW.finditer(text)})
    if missing:
        print("tools with no row in docs/tool_index.md:", ", ".join(missing))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sha")
    ap.add_argument("--tool", action="append", default=[], metavar="NAME[:VERSION]")
    ap.add_argument("--continue", dest="resume", action="store_true", help="finish after resolving real conflicts by hand")
    ap.add_argument("--no-commit", action="store_true")
    args = ap.parse_args()

    if not args.resume:
        pending = git("diff", "--name-only", "HEAD", "--", *GENERATED).split()
        if pending:
            sys.exit("uncommitted changes in " + ", ".join(pending) + ": a previous --no-commit run is not committed, and this run would overwrite its bookkeeping. Commit or reset first.")
        git("cherry-pick", "-x", "--no-commit", args.sha)
        stuck = [f for f in git("diff", "--name-only", "--diff-filter=U").split() if f not in GENERATED]
        if stuck:
            sys.exit("conflicts to resolve by hand: " + ", ".join(stuck) + "\nthen: python versions/integrate_commit.py " + args.sha + " --continue")
    git("checkout", "HEAD", "--", *GENERATED)

    tools = detect_tools(args.sha)
    for spec in args.tool:
        name, _, version = spec.partition(":")
        tools[name] = version or tools.get(name) or bump_minor(current(blocks(read("versions/tool_changelog.yaml")[0])[name]))
    renames = {}
    for tool, version in tools.items():
        theirs_version = current(blocks(git("show", f"{args.sha}:versions/tool_changelog.yaml"))[tool])
        parent = blocks(git("show", f"{args.sha}~1:versions/tool_changelog.yaml")).get(tool)
        mine = blocks(read("versions/tool_changelog.yaml")[0]).get(tool)
        if parent is not None and mine is not None:
            renames[tool] = ((current(parent), theirs_version), (current(mine), version))
        register_yaml(args.sha, tool, version)
    register_changelog(args.sha, renames)
    sync_index(args.sha, tools)
    subprocess.run([sys.executable, str(ROOT / "skills/biopipelines/build_tool_index.py")], cwd=ROOT, capture_output=True)
    print("integrated:", ", ".join(f"{t} {v}" for t, v in tools.items()) or "no tools registered")
    sources = {n: v for f in (ROOT / "biopipelines").glob("*.py") for n, v in PAIR.findall(f.read_text(encoding="utf-8"))}
    for tool, version in tools.items():
        if sources.get(tool) != version:
            print(f"set TOOL_VERSION = \"{version}\" in the {tool} wrapper (source says {sources.get(tool)})")
    if not args.no_commit:
        git("add", "-A", "--", *[str(p) for p in ("biopipelines", "pipe_scripts", "versions", "docs", "tests", "environments", "skills/biopipelines/references/tool_index.md")])
        message = git("log", "-1", "--format=%B", args.sha).rstrip() + f"\n\n(cherry picked from commit {args.sha})\n"
        out = subprocess.run(["git", "commit", "-q", "-m", message], cwd=ROOT, capture_output=True, text=True)
        sys.exit(out.stdout + out.stderr if out.returncode else 0)


if __name__ == "__main__":
    main()
