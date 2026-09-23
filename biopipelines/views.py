"""Render one step's outputs and bring them back to the machine the agent runs on.

"Show results, do not only describe them" is already the instruction in `llm/pipelines.md`, and
until now the tool layer could not honor it: `bp-visualize` is a console script, so an agent
with only the MCP tools had to hand the user a shell command and a path.

`bp-visualize` builds a self-contained page on the login node — no scheduler, so it works while
the rest of the job is still queued. That page is the full artifact, and it has to be opened in
a browser: its structure views are 3Dmol.js running client-side, which no chat can host.

What a chat *can* host is an image, so the two are separated here. `images()` returns the PNGs
the step itself produced — plots, PyMOL renders, anything a tool wrote — which a host can show
inline; `render()` produces the interactive page and fetches it. An agent should offer both:
the picture in the conversation, the page for the parts a picture cannot carry.
"""

import re

from biopipelines.remote import quote

IMAGE_SUFFIXES = (".png", ".jpg", ".jpeg", ".gif", ".webp")
MAX_IMAGE_BYTES = 4_000_000  # a chat transport is not a file server; large renders stay on disk
MAX_IMAGES = 4

_WROTE = re.compile(r"(?:wrote|written to|Wrote)\s*:?\s*(\S+\.html)", re.I)


def _walk(fs, root, depth=2):
    """(path, name) for files under `root`, `depth` levels down. Small by construction."""
    walk_files = getattr(fs, "walk_files", None)
    if walk_files is not None:
        # Over ssh a listdir per folder is a round trip each; one find answers the whole walk.
        try:
            return walk_files(root, depth)
        except Exception:
            return []
    found = []
    stack = [(root, 0)]
    while stack:
        current, level = stack.pop()
        try:
            entries = fs.listdir(current)
        except Exception:
            continue
        for name, is_dir in entries:
            path = fs.join(current, name)
            if not is_dir:
                found.append((path, name))
                continue
            if level >= depth:
                continue
            # `_extras` is where a rendered page and its images land, so it is the one
            # underscore folder worth descending into; the rest is framework plumbing.
            if name.startswith("_") and name != "_extras":
                continue
            stack.append((path, level + 1))
    return found


def images(fs, job_dir, step, limit=MAX_IMAGES):
    """Image files this step wrote, as (path, name) — what a host can render inline."""
    found = [(path, name) for path, name in _walk(fs, fs.join(job_dir, step))
             if name.lower().endswith(IMAGE_SUFFIXES)]
    return sorted(found, key=lambda pair: pair[1])[:limit]


def page_path(output, fs, job_dir, step):
    """Where `bp-visualize` put its page: what it said, or the documented default."""
    said = _WROTE.search(output or "")
    if said:
        return said.group(1)
    return fs.join(job_dir, step, "_extras", f"{step}_view.html")


def render(ssh, job_dir, step, descending=None, ascending=None, max_items=5, ids=None):
    """Build the step's page on the cluster. Returns what happened and where the page is.

    It runs on the login node and needs no scheduler, so this works against a job that is still
    in progress — which is the point: a finished step should be shown the moment it finishes,
    not after the campaign ends.
    """
    folder = ssh.join(job_dir, step)
    if not ssh.is_dir(folder):
        return {"ok": False, "failures": [f"no step folder at {folder}"], "output": ""}

    # Not the bare name: a non-interactive ssh has only the system PATH, and bp-visualize
    # lives in the environment. `script()` resolves it beside the host's interpreter.
    binary = ssh.script("bp-visualize") if hasattr(ssh, "script") else "bp-visualize"
    parts = [f"cd {quote(ssh.repo)} && {getattr(ssh, 'env_prefix', '')}"
             f"{quote(binary)} {quote(folder)}"]
    if descending:
        parts.append(f"--descending {quote(descending)}")
    elif ascending:
        parts.append(f"--ascending {quote(ascending)}")
    if max_items:
        parts.append(f"--max-items {int(max_items)}")
    if ids:
        parts.append("--ids " + quote(",".join(ids)))

    command = " ".join(parts)
    code, out, err = ssh.run(command, timeout=max(ssh.timeout, 300))
    combined = out + ("\n" + err if err.strip() else "")
    if code != 0:
        return {"ok": False, "failures": [f"bp-visualize exited {code}"], "output": combined,
                "command": command}
    return {"ok": True, "failures": [], "output": combined, "command": command,
            "page": page_path(combined, ssh, job_dir, step)}


def summarize(result, page_local=None, shown=(), skipped=()):
    """What an agent reads back after rendering a step."""
    if not result["ok"]:
        lines = ["Could not render that step."]
        lines += [f"  - {f}" for f in result["failures"]]
        tail = "\n".join(result.get("output", "").splitlines()[-10:])
        if tail:
            lines.append("\nLast lines:\n" + tail)
        return "\n".join(lines)

    lines = [f"Rendered {result['page']}"]
    if page_local:
        lines.append(f"  fetched to {page_local} — open it for the interactive view")
    if shown:
        lines.append("  inline: " + ", ".join(shown))
    if skipped:
        lines.append("  too large to show inline, in the page instead: " + ", ".join(skipped))
    if not shown:
        lines.append("  This step wrote no image files, so there is nothing to show inline. "
                     "The page carries the structure views, which are interactive and only "
                     "render in a browser.")
    return "\n".join(lines)
