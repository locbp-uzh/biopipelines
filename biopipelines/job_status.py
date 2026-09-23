"""Read a run's state: which steps finished, which failed, and what their logs say.

A half-failed campaign is the case that matters. An agent that can only see stdout has to grep
its way through `Logs/` to find out that step 4 of 9 failed and the rest never ran; this module
answers that from the completion markers the framework already writes.

The unit is a job folder — `<output_root>/<Project>/<Job>_NNN/` — containing, per step, a
`<NNN>_<Tool>/` output folder, a `<NNN>_<Tool>_COMPLETED` / `_FAILED` / `_WARNING` marker
alongside it, and `Logs/<NNN>_<Tool>.log`.

Every read goes through a filesystem object (`remote.LocalFS` or `remote.Ssh`), so the same
code answers for a run on this disk and one on a cluster login node. The cluster is the main
platform; nothing here may assume local paths.
"""

import json
import re

from biopipelines.remote import LocalFS

MARKERS = {"COMPLETED": "completed", "FAILED": "failed", "WARNING": "warning"}

# `001_Mock_COMPLETED` -> index 001, tool Mock. The tool name itself may contain underscores.
_MARKER = re.compile(r"^(\d+)_(.+)_(" + "|".join(MARKERS) + r")$")
_STEP_DIR = re.compile(r"^(\d+)_(.+)$")

LOG_TAIL_LINES = 100


def steps(job_dir, fs=None):
    """Every step of the run, in execution order, with its status.

    A step folder with no marker is `pending`: either still running, or it never started
    because an earlier step failed. Both look the same on disk, which is itself worth saying.
    """
    fs = fs or LocalFS()
    if not fs.is_dir(job_dir):
        return []

    log_names = {name for name, is_dir in fs.listdir(fs.join(job_dir, "Logs")) if not is_dir}

    entries = fs.listdir(job_dir)
    found, directories = {}, []
    for name, is_dir in entries:
        if is_dir:
            if _STEP_DIR.match(name):
                directories.append(name)
        else:
            match = _MARKER.match(name)
            if match:
                index, tool, marker = match.groups()
                found[(index, tool)] = {"status": MARKERS[marker]}

    # A `Suffix` lands on the output folder and the log but not on the marker: the marker is
    # `001_PDB_COMPLETED` while the folder is `001_PDB_stitched_001/`. Matching them by prefix
    # is what stops one step being counted twice — once done, once invented as pending.
    for name in sorted(directories):
        key = next((k for k in found if name == f"{k[0]}_{k[1]}"
                    or name.startswith(f"{k[0]}_{k[1]}_")), None)
        if key is None:
            key = _STEP_DIR.match(name).groups()
            found.setdefault(key, {"status": "pending"})
        found[key].setdefault("folder", name)

    out = []
    for (index, tool), data in sorted(found.items()):
        step = f"{index}_{tool}"
        folder = data.get("folder", step)
        out.append({"index": index, "tool": tool, "step": step, "status": data["status"],
                    "has_log": f"{folder}.log" in log_names or f"{step}.log" in log_names,
                    "folder": folder})
    return out


def _header(job_dir, fs):
    """The run header, with the version, commit and variant taken from `manifest.json` when there is one.

    The manifest is written by the machine that generated the run; an operations header written by `bp_submit` describes the laptop that submitted it.
    """
    found = None
    path = fs.join(job_dir, "_operations.jsonl")
    if fs.is_file(path):
        for line in fs.read_text(path).splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                entry = json.loads(line)
            except ValueError:
                continue
            if entry.get("action") == "run":
                found = dict(entry)
                break
    try:
        from biopipelines import manifest
        record = manifest.read(job_dir, fs=fs, environments_too=False)
    except Exception:
        record = None
    if record:
        found = found or {}
        for key, source in (("biopipelines", "biopipelines"), ("commit", "commit"),
                            ("config_variant", "variant")):
            if record.get(source):
                found[key] = record[source]
    return found


def status(job_dir, fs=None):
    """A run's state in one object: per-step status, counts, and the header if one was written."""
    fs = fs or LocalFS()
    step_list = steps(job_dir, fs=fs)
    counts = {}
    for step in step_list:
        counts[step["status"]] = counts.get(step["status"], 0) + 1

    failed = [s["step"] for s in step_list if s["status"] == "failed"]
    return {"job": str(job_dir).replace("\\", "/").rstrip("/").rsplit("/", 1)[-1],
            "job_dir": str(job_dir),
            "exists": fs.is_dir(job_dir),
            "steps": step_list,
            "counts": counts,
            "failed": failed,
            # The first failure is the one to read; everything after it may be a consequence.
            "first_failure": failed[0] if failed else None,
            "header": _header(job_dir, fs) if step_list else None}


def log(job_dir, step, tail=LOG_TAIL_LINES, fs=None):
    """The tail of one step's log. `step` is `001_Mock`, or just the tool name if unambiguous."""
    fs = fs or LocalFS()
    logs = fs.join(job_dir, "Logs")
    names = [name for name, is_dir in fs.listdir(logs) if not is_dir and name.endswith(".log")]

    wanted = f"{step}.log"
    if wanted not in names:
        # `001_PDB` may be logged as `001_PDB_stitched_001.log` when the step carries a Suffix,
        # and a bare tool name as `006_AlphaFold.log`. Try both directions before giving up.
        matches = ([n for n in names if n.startswith(f"{step}_")]
                   or [n for n in names if n.endswith(f"_{step}.log")])
        if len(matches) != 1:
            return None
        wanted = matches[0]

    lines = fs.read_text(fs.join(logs, wanted)).splitlines()
    return {"step": wanted[:-4], "lines": len(lines), "shown": min(tail, len(lines)),
            "text": "\n".join(lines[-tail:])}


def summarize(state):
    """One screen an agent can act on, rather than a directory listing to interpret."""
    if not state["exists"]:
        return f"No run directory at {state['job_dir']}."
    if not state["steps"]:
        return f"{state['job']}: no steps on disk yet."

    counts = ", ".join(f"{n} {name}" for name, n in sorted(state["counts"].items()))
    lines = [f"{state['job']}: {counts}"]
    header = state["header"]
    if header:
        bits = [f"{k}={header[k]}" for k in ("biopipelines", "commit", "config_variant", "host")
                if header.get(k)]
        if bits:
            lines.append("  " + "  ".join(bits))
    for step in state["steps"]:
        mark = {"completed": "ok", "failed": "FAILED", "warning": "warning",
                "pending": "pending"}[step["status"]]
        lines.append(f"  {step['step']:<28} {mark}")

    if state["first_failure"]:
        note = (f"\nFirst failure: {state['first_failure']}. "
                f'Read it with bp_logs(step="{state["first_failure"]}").')
        # Only claim the run stopped when it actually did: a failed step does not always halt
        # the rest, and saying otherwise sends the agent looking for a cascade that never happened.
        after = state["steps"][[s["step"] for s in state["steps"]].index(state["first_failure"]) + 1:]
        if after and all(s["status"] == "pending" for s in after):
            note += " Everything after it is pending, so the run stopped there."
        elif any(s["status"] == "completed" for s in after):
            note += " Later steps still ran, so the failure did not halt the pipeline."
        lines.append(note)
    return "\n".join(lines)
