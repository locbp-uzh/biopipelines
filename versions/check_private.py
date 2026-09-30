#!/usr/bin/env python3
"""List files marked `#PRIVATE` that a publish would carry into a public branch.

A file with a line beginning `#PRIVATE` (a valid comment as-is in Python, YAML and
shell) must never reach `public`. Line-start, so a file that merely mentions the
marker -- this one, for instance -- is not itself flagged. Merging into `public` publishes full
ancestry, so a file that was private in one commit stays readable in history even
if a later commit deletes it -- which is why this checks the whole range rather
than only the tip.

    python versions/check_private.py                 # HEAD vs origin/public
    python versions/check_private.py origin/main     # HEAD vs that ref
    python versions/check_private.py --list          # every marked file at HEAD

Exit status is 1 when anything marked would be published, 0 otherwise.
"""
import subprocess
import sys

MARKER = "#PRIVATE"
DEFAULT_TARGET = "origin/public"


def _git(*args):
    return subprocess.run(["git", *args], capture_output=True, text=True).stdout


def _marked(ref, paths):
    """Paths whose blob at `ref` contains the marker. Binary files are skipped."""
    out = []
    for path in paths:
        blob = subprocess.run(["git", "show", f"{ref}:{path}"],
                              capture_output=True)
        if blob.returncode:
            continue
        try:
            text = blob.stdout.decode("utf-8")
        except UnicodeDecodeError:
            continue
        if any(ln.lstrip().startswith(MARKER) for ln in text.splitlines()):
            out.append(path)
    return out


def main(argv):
    if "--list" in argv:
        paths = _git("ls-tree", "-r", "--name-only", "HEAD").split()
        for p in sorted(_marked("HEAD", paths)):
            print(p)
        return 0

    target = argv[0] if argv else DEFAULT_TARGET
    if not _git("rev-parse", "--verify", "--quiet", target).strip():
        print(f"no such ref: {target}", file=sys.stderr)
        return 2

    # Every file touched anywhere in the range, not just those differing at the
    # tip: a file added and later deleted is still published with the history.
    paths = sorted(set(_git("log", "--name-only", "--pretty=format:",
                            f"{target}..HEAD").split("\n")) - {""})
    hits = sorted(set(_marked("HEAD", paths)) | set(_marked_in_history(target, paths)))
    if not hits:
        print(f"no {MARKER} files in {target}..HEAD")
        return 0

    print(f"{len(hits)} {MARKER} file(s) would be published to {target}:", file=sys.stderr)
    for p in hits:
        print(f"  {p}", file=sys.stderr)
    return 1


def _marked_in_history(target, paths):
    """Paths marked in any commit of the range, including ones deleted since."""
    out = set()
    for commit in _git("rev-list", f"{target}..HEAD").split():
        present = set(_git("ls-tree", "-r", "--name-only", commit).split())
        out |= set(_marked(commit, [p for p in paths if p in present]))
    return out


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
