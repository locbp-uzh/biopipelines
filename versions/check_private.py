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

What a merge publishes is exactly the objects reachable from HEAD and not from the
target, which `git rev-list --objects` names in one pass: every version of every
file, under the path it had. Reading them through a single `git cat-file --batch`
keeps the cost at one blob read per blob. Walking commit by commit and showing each
path at each commit re-read every unchanged file once per commit instead -- 211
commits over 385 paths asked git for 81,000 blobs to read 3,359, which on Windows,
where a process costs milliseconds, is twenty minutes against seconds.

Exit status is 1 when anything marked would be published, 0 otherwise.
"""
import subprocess
import sys
import tempfile

MARKER = b"#PRIVATE"
DEFAULT_TARGET = "origin/public"


def _git(*args):
    return subprocess.run(["git", *args], capture_output=True, text=True).stdout


def _is_marked(blob: bytes) -> bool:
    """True when some line of the blob starts with the marker, leading blanks aside."""
    return any(line.lstrip().startswith(MARKER) for line in blob.splitlines())


def _read_blobs(shas):
    """Yield (sha, contents) for each sha, through one `git cat-file --batch`.

    The request list goes in through a temporary file rather than a pipe we write
    ourselves: git answers as it reads, so a writer that is not yet reading fills
    git's stdout, git stops, and the two of us wait on each other forever.
    """
    if not shas:
        return
    with tempfile.TemporaryFile() as requests:
        requests.write(("\n".join(shas) + "\n").encode())
        requests.seek(0)
        proc = subprocess.Popen(["git", "cat-file", "--batch"],
                                stdin=requests, stdout=subprocess.PIPE)
    out = proc.stdout
    for _ in shas:
        header = out.readline()
        if not header:
            break
        parts = header.split()
        if len(parts) < 3:          # "<sha> missing"
            continue
        sha, kind, size = parts[0].decode(), parts[1], int(parts[2])
        body = out.read(size)
        out.read(1)                 # the newline git appends after the body
        if kind == b"blob":
            yield sha, body
    out.close()
    proc.wait()


def _marked_paths(pairs):
    """Paths of the (sha, path) pairs whose blob carries the marker."""
    by_sha = {}
    for sha, path in pairs:
        by_sha.setdefault(sha, set()).add(path)
    hits = set()
    for sha, blob in _read_blobs(list(by_sha)):
        if _is_marked(blob):
            hits |= by_sha[sha]
    return hits


def _pairs_at(ref):
    """(sha, path) of every blob in the tree of `ref`."""
    pairs = []
    for line in _git("ls-tree", "-r", ref).splitlines():
        meta, _, path = line.partition("\t")
        fields = meta.split()
        if len(fields) == 3 and fields[1] == "blob":
            pairs.append((fields[2], path))
    return pairs


def _pairs_published(target):
    """(sha, path) of every object a merge of HEAD into `target` would add.

    `rev-list --objects` lists them once each, whatever their history: a file added
    and later deleted is here, because its blob is still reachable from HEAD.
    """
    pairs = []
    for line in _git("rev-list", "--objects", f"{target}..HEAD").splitlines():
        sha, _, path = line.partition(" ")
        if path:                    # a commit has no path; trees do, and are skipped by type
            pairs.append((sha, path))
    return pairs


def main(argv):
    if "--list" in argv:
        for path in sorted(_marked_paths(_pairs_at("HEAD"))):
            print(path)
        return 0

    target = argv[0] if argv else DEFAULT_TARGET
    if not _git("rev-parse", "--verify", "--quiet", target).strip():
        print(f"no such ref: {target}", file=sys.stderr)
        return 2

    hits = sorted(_marked_paths(_pairs_published(target)))
    marker = MARKER.decode()
    if not hits:
        print(f"no {marker} files in {target}..HEAD")
        return 0

    print(f"{len(hits)} {marker} file(s) would be published to {target}:", file=sys.stderr)
    for path in hits:
        print(f"  {path}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
