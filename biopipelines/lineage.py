"""Where the designs went: how many IDs each step produced, and how many it dropped.

This is the question a design campaign actually turns on — not "did it run" but "I started with
500 and ordered 12, and which step removed the rest". The framework already writes the answer
and nobody assembles it: each step's `<stream>_map.csv` lists the IDs it produced, and its
`missing.csv` lists the IDs that entered and did not come out.

Every other way of answering this involves reading result tables by hand and trusting the
arithmetic. It is also the piece an agent that merely generates code cannot reproduce, because
it is a property of having executed a typed graph rather than of the code that was written.

One `tally` call collects every count for a whole run, so a 36-step campaign costs one ssh
round trip rather than seventy.
"""

from biopipelines.remote import LocalFS

MAP_SUFFIX = "_map.csv"
MISSING = "missing.csv"


def _step_key(folder):
    """`005_ProteinMPNN_stitched_001` -> ('005', '005_ProteinMPNN_stitched_001')."""
    index = folder.split("_", 1)[0]
    return (index if index.isdigit() else "", folder)


def collect(job_dir, fs=None):
    """Per step: the IDs produced per stream, and how many were dropped.

    Steps with neither a map table nor a missing table are omitted — a tool that writes no
    stream has no lineage to report, and inventing a zero would read as attrition.
    """
    fs = fs or LocalFS()
    counts = fs.tally(job_dir, [f"*{MAP_SUFFIX}", MISSING])

    steps = {}
    for path, rows in counts.items():
        parts = path.replace("\\", "/").split("/")
        if len(parts) < 2:
            continue
        # `002_DSSP/dssp/dssp_map.csv` and `002_DSSP/tables/missing.csv` both belong to the
        # step, not to the stream folder they happen to sit in.
        folder, filename = parts[0], parts[-1]
        entry = steps.setdefault(folder, {"folder": folder, "streams": {}, "missing": 0})
        if filename == MISSING:
            entry["missing"] = rows
        elif filename.endswith(MAP_SUFFIX):
            entry["streams"][filename[: -len(MAP_SUFFIX)]] = rows

    ordered = sorted(steps.values(), key=lambda s: _step_key(s["folder"]))
    for step in ordered:
        step["produced"] = max(step["streams"].values()) if step["streams"] else 0
    _own_drops(job_dir, fs, ordered)
    return ordered


def _own_drops(job_dir, fs, ordered):
    """Count each dropped ID once, at the step that dropped it.

    Every step merges its upstream `missing.csv` into its own, so a raw row count charges one filter's drops again to every step after it. A row an earlier step already listed is inherited, not lost here. Left as the raw count when the rows cannot be read.
    """
    if not any(step["missing"] for step in ordered):
        return
    try:
        rows = fs.id_columns(job_dir, [MISSING], extra=("removed_by",))
    except Exception:
        return
    by_folder = {}
    for path, records in rows.items():
        by_folder.setdefault(path.replace("\\", "/").split("/")[0], []).extend(records)
    seen = set()
    for step in ordered:
        if step["folder"] not in by_folder:
            continue
        keys = {(r.get("id", ""), r.get("removed_by", "")) for r in by_folder[step["folder"]]}
        own = keys - seen
        step["inherited"] = len(keys) - len(own)
        step["missing"] = len(own)
        step["own_drops"] = own
        seen |= keys


def summarize(steps, job=None, max_steps=40):
    """The attrition of a campaign in one screen."""
    if not steps:
        return (f"No lineage tables in {job or 'this run'}. Steps write `<stream>_map.csv` and "
                "`missing.csv`; a run whose tools produce neither has nothing to trace.")

    lines = [f"{job or 'run'}: lineage across {len(steps)} step(s) that report IDs", ""]

    # A large campaign has hundreds of steps and listing them all buries the attrition, which is
    # the only part anyone reads. Above the cap, keep the steps that lost IDs and nothing else.
    shown, elided = steps, 0
    if len(steps) > max_steps:
        shown = [s for s in steps if s["missing"]] or steps[:max_steps]
        elided = len(steps) - len(shown)

    for step in shown:
        streams = ", ".join(f"{name} {n}" for name, n in sorted(step["streams"].items()))
        detail = streams or "no stream table"
        if step["missing"]:
            detail += f"  [-{step['missing']} dropped]"
        lines.append(f"  {step['folder']:<34} {detail}")
    if elided:
        lines.append(f"  … {elided} step(s) with no dropped IDs not listed")

    dropped = [s for s in steps if s["missing"]]
    total = sum(s["missing"] for s in dropped)
    if dropped:
        worst = max(dropped, key=lambda s: s["missing"])
        lines.append("")
        lines.append(f"{total} ID(s) dropped across {len(dropped)} step(s); "
                     f"the largest single loss is {worst['missing']} at {worst['folder']}.")
        lines.append("Read the dropped ids with bp_table(step=\"" + worst["folder"]
                     + "\", table=\"missing.csv\").")
    else:
        lines.append("")
        lines.append("No step reported dropped IDs.")
    return "\n".join(lines)


def tables(job_dir, step, fs=None):
    """The CSV tables one step wrote, as paths relative to the step folder.

    A step does not keep its tables beside itself: standalone tables land in `<step>/tables/`
    and a stream's map table in the stream's own folder, so listing `<step>/*.csv` reports a
    real run as having written nothing. Found live on a template run whose `002_Scripting`
    holds `tables/metrics.csv` and `structures/structures_map.csv`.

    Folders beginning with `_` are the framework's own (`_configuration`, `_extras`) and are
    left out — they are not results.
    """
    fs = fs or LocalFS()
    prefix = str(step).replace("\\", "/").strip("/") + "/"
    found = []
    for path in fs.tally(job_dir, ["*.csv"]):
        rel = path.replace("\\", "/")
        if not rel.startswith(prefix):
            continue
        inner = rel[len(prefix):]
        if inner.startswith("_") or "/_" in inner:
            continue
        found.append(inner)
    return sorted(found)


def _match_table(wanted, available):
    """Resolve what a caller typed against the paths a step actually wrote.

    `metrics`, `metrics.csv` and `tables/metrics.csv` are the same request; asking someone to
    know which folder the framework chose defeats the point of the tool.
    """
    text = str(wanted).replace("\\", "/").strip("/").lower()
    for rel in available:
        if rel.lower() == text:
            return rel
    for rel in available:
        base = rel.rsplit("/", 1)[-1].lower()
        if base == text or base[:-4] == text:
            return rel
    return None


def read_table(job_dir, step, table, limit=20, fs=None):
    """A step's table as rows. Truncated by default: these run to thousands of lines."""
    fs = fs or LocalFS()
    available = tables(job_dir, step, fs=fs)
    match = _match_table(table, available)
    path = fs.join(job_dir, step, *match.split("/")) if match else None
    if path is None or not fs.is_file(path):
        return {"error": f"no table {table!r} in {step}",
                "available": available}
    lines = fs.read_text(path).splitlines()
    if not lines:
        return {"step": step, "table": match, "rows": 0, "header": "", "lines": []}
    return {"step": step, "table": match, "header": lines[0],
            "rows": len(lines) - 1, "shown": min(limit, len(lines) - 1),
            "lines": lines[1:limit + 1]}
