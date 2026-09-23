"""Keep a project folder's two documents: `PROJECT.md` and `HISTORY.md`.

A project folder holds the job directories a campaign produced plus the local analysis around
them. Two documents carry it:

* `PROJECT.md` — goal, inputs, and the validated protocol: what we do and why we believe it.
  Rewritten as understanding improves, never appended to. It carries judgment, so an agent
  drafts the skeleton and proposes edits, and the scientist owns the content.
* `HISTORY.md` — append-only and dated: runs, findings, restructures. Mechanical enough that
  an agent maintains it.

The rule the whole module is built around: **never overwrite what a human wrote.** `scaffold()`
creates only what is missing, and `append_history()` only ever adds at the end. An agent that
silently rewrites `PROJECT.md` destroys exactly the thing that makes these folders worth having.
"""

import datetime
import re

PROJECT = "PROJECT.md"
HISTORY = "HISTORY.md"
README = "README.md"

PROJECT_TEMPLATE = """# {name}

## Goal

<!-- What are we trying to make, and what counts as success? One paragraph. -->

## Inputs

<!-- Scaffolds, ligands, structures, datasets. Where each came from. -->

## Protocol

<!-- The validated protocol: what we run and why we believe it. For each step, state the
     choice, the evidence behind it, and what would falsify it. Rewrite this section as
     understanding improves — the chronology, including what failed, lives in HISTORY.md. -->

## Status

<!-- Where the campaign stands today, and what is live vs superseded. -->
"""

HISTORY_TEMPLATE = """# {name} — history

Append-only, newest at the bottom. Runs, findings and restructures, each dated. The current
protocol and the reasoning behind it live in `PROJECT.md`; this file is what happened, including
what did not work.
"""

README_TEMPLATE = """# {name}

<!-- One paragraph on what this project is, then where things are. The scientific content
     belongs in PROJECT.md. -->

| Path | Contents |
|---|---|
"""


def _today():
    return datetime.date.today().isoformat()


def _fs(fs=None):
    """The filesystem to act on — this machine unless a cluster's was passed.

    A project folder lives beside the runs, and the runs are on the cluster. Reading it through
    local `pathlib` reported an existing project as empty rather than as unreachable, which is
    the one answer worse than an error.
    """
    if fs is not None:
        return fs
    from biopipelines.remote import LocalFS
    return LocalFS()


def _name_of(project_dir):
    return str(project_dir).replace("\\", "/").rstrip("/").rsplit("/", 1)[-1]


def paths(project_dir, fs=None):
    fs = _fs(fs)
    return {"project": fs.join(project_dir, PROJECT),
            "history": fs.join(project_dir, HISTORY),
            "readme": fs.join(project_dir, README)}


def scaffold(project_dir, name=None, readme=False, fs=None):
    """Create the documents that are missing. Returns {name: 'created' | 'kept'}.

    Never touches a file that already exists — running this on an established project is safe
    and is how a project that predates the convention gets its missing half.
    """
    fs = _fs(fs)
    name = name or _name_of(project_dir)

    wanted = {"project": (PROJECT, PROJECT_TEMPLATE), "history": (HISTORY, HISTORY_TEMPLATE)}
    if readme:
        wanted["readme"] = (README, README_TEMPLATE)

    result = {}
    for key, (filename, template) in wanted.items():
        target = fs.join(project_dir, filename)
        if fs.is_file(target):
            result[key] = "kept"
            continue
        # `append_text` to an absent file writes it and its parent, on either filesystem, so
        # creating a document needs no separate primitive and no mkdir of its own.
        result[key] = "created" if fs.append_text(target, template.format(name=name)) else "failed"
    return result


def append_history(project_dir, title, body="", date=None, fs=None):
    """Add one dated entry to the end of HISTORY.md, creating the file if it is absent.

    Append-only by construction: there is no code path here that rewrites an existing entry.
    """
    fs = _fs(fs)
    target = fs.join(project_dir, HISTORY)
    if not fs.is_file(target):
        scaffold(project_dir, fs=fs)

    entry = f"\n## {date or _today()} — {title}\n"
    if body:
        entry += "\n" + body.strip() + "\n"
    if not fs.append_text(target, entry):
        raise OSError(f"could not append to {target}")
    return entry.strip()


def history_entries(project_dir, fs=None):
    """The dated entries already recorded, so an agent can avoid writing a duplicate."""
    fs = _fs(fs)
    target = fs.join(project_dir, HISTORY)
    if not fs.is_file(target):
        return []
    text = fs.read_text(target)
    return [{"date": m.group(1), "title": m.group(2).strip()}
            for m in re.finditer(r"^## (\d{4}-\d{2}-\d{2}) — (.+)$", text, re.M)]


def job_dirs(project_dir, fs=None):
    """Job folder names in the project, newest first — `<Job>_NNN/` beside the documents."""
    fs = _fs(fs)
    if not fs.is_dir(project_dir):
        return []
    return sorted((name for name, is_dir in fs.listdir(project_dir)
                   if is_dir and re.match(r"^.+_\d{3,}$", name)), reverse=True)


def survey(project_dir, fs=None):
    """What the project looks like now: which documents exist, which jobs, which entries.

    `exists` is reported separately from the rest: a folder that is not there is a different
    answer from a folder with nothing in it, and conflating them told a user their project was
    empty when the tool was simply looking at the wrong machine.
    """
    fs = _fs(fs)
    exists = fs.is_dir(project_dir)
    have = paths(project_dir, fs=fs)
    return {"project_dir": str(project_dir),
            "name": _name_of(project_dir),
            "exists": exists,
            "documents": {k: (exists and fs.is_file(v)) for k, v in have.items()},
            "jobs": job_dirs(project_dir, fs=fs) if exists else [],
            "entries": history_entries(project_dir, fs=fs) if exists else []}
