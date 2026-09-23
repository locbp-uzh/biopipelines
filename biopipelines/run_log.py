"""The operational record of what was done to a run, written next to the run's own outputs.

`<Job>_NNN/_operations.jsonl` answers "what happened to this campaign" — submitted, resubmitted,
cancelled, fetched, polled — and because it lives inside the job folder it travels with the
results when they are copied off the cluster, and is written identically on every backend.

This replaces `llm/log.sh`, which tee'd every `ssh`/`scp` into `llm/logs/YYYY-MM-DD.log`: that
logged the transport rather than the work, filed it by date rather than by run, covered only
the cluster, and depended on the agent remembering to use the wrapper.

It is deliberately *not* the provenance record. Which inputs, tool versions and parameters
produced which output IDs is already written per step into `_configuration/<tool>_config.json`
and the stream map tables. This file is the operational layer above that.

JSON Lines, because it is append-only, survives a truncated write, and stays greppable by a
human reading it over ssh.
"""

import datetime
import io
import json
import os
import pathlib

FILENAME = "_operations.jsonl"

# A record is a breadcrumb, not a transcript: a tool's full output belongs in
# `<Job>/Logs/<NNN>_<Tool>.log`, not inlined here.
OUTPUT_CHARS = 2000


def path_for(job_dir, fs=None):
    return (fs or _local()).join(str(job_dir), FILENAME)


def _local():
    from biopipelines.remote import LocalFS
    return LocalFS()


def _now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def record(job_dir, action, fs=None, **fields):
    """Append one operation to the run's log. Returns the record, or None if it could not be written.

    `fs` is where the run lives. It defaults to this machine, but a run on a cluster needs the
    ssh filesystem: writing a cluster path through local pathlib produces nothing, which is how
    the first real submission ended up with no record at all.

    Never raises: a run must not fail because its log could not be written. A caller that needs
    to know checks the return value — and should, because silence here once hid a total absence.
    """
    fs = fs or _local()
    entry = {"ts": _now(), "action": action}
    entry.update({k: v for k, v in fields.items() if v is not None})
    line = json.dumps(entry, default=str) + "\n"
    try:
        return entry if fs.append_text(path_for(job_dir, fs), line) else None
    except Exception:
        return None


def command(job_dir, argv, exit_code=None, output=None, host=None, fs=None, **fields):
    """Record a command run against the infrastructure — an ssh/scp call, a submission, a poll.

    `output` is truncated: this file is a breadcrumb trail, not a transcript.
    """
    if output is not None and len(output) > OUTPUT_CHARS:
        output = output[:OUTPUT_CHARS] + f"\n... [{len(output) - OUTPUT_CHARS} more chars]"
    return record(job_dir, "command", fs=fs,
                  argv=list(argv) if not isinstance(argv, str) else argv,
                  exit_code=exit_code, output=output, host=host, **fields)


def _git(*args):
    try:
        import subprocess
        out = subprocess.run(["git", *args],
                             cwd=pathlib.Path(__file__).resolve().parent.parent,
                             capture_output=True, text=True, timeout=5)
        return out.stdout.strip() if out.returncode == 0 else None
    except Exception:
        return None


def _git_commit():
    # Full sha: `--short` picks its length per repository, so two hosts can abbreviate the same commit differently.
    return _git("rev-parse", "HEAD") or None


def _git_dirty():
    """Uncommitted edits to tracked framework files; None when git is unavailable."""
    status = _git("status", "--porcelain", "--untracked-files=no")
    return None if status is None else bool(status)


def _header_fields():
    """What silently differs between two runs of the same script."""
    import platform
    from biopipelines import __version__
    return {"biopipelines": __version__,
            "commit": _git_commit(),
            "dirty": _git_dirty(),
            "config_variant": os.environ.get("BIOPIPELINES_CONFIG_VARIANT"),
            "host": platform.node() or None,
            "python": platform.python_version()}


def header(job_dir, fs=None, **fields):
    """The once-per-run record of what produced it — the first question asked when a rerun differs.

    Written as the run's first entry. Cheap enough to be unconditional.
    """
    return record(job_dir, "run", fs=fs, **{**_header_fields(), **fields})


class RunLog:
    """A run's log, usable before its directory exists.

    `bp_submit` knows the submission parameters before the job folder is created, and those are
    exactly what the record needs. Records made before `bind()` are held in memory and written
    in order once the directory is known; after `bind()` they go straight to disk.
    """

    def __init__(self, job_dir=None, fs=None):
        self._job_dir = None
        self._fs = fs
        self._pending = []
        if job_dir is not None:
            self.bind(job_dir)

    @property
    def job_dir(self):
        return self._job_dir

    @property
    def pending(self):
        return list(self._pending)

    def bind(self, job_dir):
        """Point at a directory and flush anything buffered, oldest first."""
        self._job_dir = str(job_dir)
        buffered, self._pending = self._pending, []
        for action, fields in buffered:
            record(self._job_dir, action, fs=self._fs, **fields)
        return self

    def record(self, action, **fields):
        if self._job_dir is None:
            # Stamp now, not at flush: the time the thing happened is the useful one.
            fields.setdefault("ts", _now())
            self._pending.append((action, fields))
            return None
        return record(self._job_dir, action, fs=self._fs, **fields)

    def command(self, argv, **fields):
        if self._job_dir is None:
            return self.record("command", argv=list(argv), **fields)
        return command(self._job_dir, argv, fs=self._fs, **fields)

    def header(self, submitted_from=False, **fields):
        # Collect the fields now even when buffering: they describe the machine and checkout
        # as they were at submission, which is the point.
        if submitted_from:
            # Written by the laptop for a run on a cluster: its identity is the submitter's, not the run's.
            return self.record("run", submitted_from=_header_fields(), **fields)
        return self.record("run", **{**_header_fields(), **fields})


def read(job_dir, fs=None):
    """Every operation recorded for this run, oldest first.

    A malformed line is skipped rather than fatal — a log truncated by a killed job should
    still be readable.
    """
    fs = fs or _local()
    target = path_for(job_dir, fs)
    if not fs.is_file(target):
        return []
    entries = []
    for line in fs.read_text(target).splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            entries.append(json.loads(line))
        except ValueError:
            continue
    return entries


def find_runs(root):
    """Job folders under `root` that carry an operations log, newest first."""
    root = pathlib.Path(root)
    if not root.exists():
        return []
    found = [p.parent for p in root.rglob(FILENAME)]
    return sorted(found, key=lambda p: os.path.getmtime(p / FILENAME), reverse=True)
