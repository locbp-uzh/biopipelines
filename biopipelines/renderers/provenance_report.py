"""The run as a referee would want to read it: what produced it, and where the designs went.

`pipeline.html` shows the campaign's shape and its outputs. It does not say which software ran, and it does not say why a campaign that started with five hundred designs ended with twelve. Both facts exist already — the manifest holds the first, `<stream>_map.csv` and `missing.csv` the second — but they exist as an agent's API, reachable through `bp_provenance` and `bp_table`. A reviewer does not call MCP tools.

So this renders them as one page: the recorded identity of the run at the top, the attrition through the steps in the middle, and every dropped ID with the reason its own step gave, at the bottom. Self-contained, no network, opens from `file://` after being copied off a cluster.

It reads through `fs`, so a run that lives on a cluster renders without being fetched first.
"""

import csv
import html
import os
from typing import Any, Dict, List, Optional

from .. import ancestry, lineage, manifest

MAX_DROPPED_ROWS = 200
MAX_ANCESTRY_ROWS = 400



def _esc(value: Any) -> str:
    return html.escape("" if value is None else str(value), quote=True)


def build(job_dir: str, fs=None, max_dropped: int = MAX_DROPPED_ROWS) -> Dict[str, Any]:
    """Everything the page shows, gathered in one pass so rendering touches no filesystem."""
    steps = lineage.collect(job_dir, fs=fs)
    record = manifest.read(job_dir, fs=fs)
    by_folder = {}
    if record:
        for entry in record.get("tools") or []:
            by_folder[str(entry.get("order"))] = entry

    dropped = {}
    for step in steps:
        if not step.get("missing"):
            continue
        own = step.get("own_drops")
        try:
            # With own drops known, filter the whole table first: truncating first lost the rows past the limit.
            found = lineage.read_table(job_dir, step["folder"], "missing.csv",
                                       limit=10 ** 9 if own is not None else max_dropped, fs=fs)
        except Exception:
            continue
        if "error" not in found:
            if own is not None:
                # List only this step's own drops; the rows it inherited are shown where they were lost.
                import csv as _csv
                header = next(_csv.reader([found["header"]]))
                id_at = header.index("id") if "id" in header else None
                by_at = header.index("removed_by") if "removed_by" in header else None
                if id_at is not None:
                    kept = []
                    for line in found["lines"]:
                        cells = next(_csv.reader([line]), [])
                        key = (cells[id_at] if id_at < len(cells) else "",
                               cells[by_at] if by_at is not None and by_at < len(cells) else "")
                        if key in own:
                            kept.append(line)
                    found = {**found, "lines": kept[:max_dropped], "shown": min(len(kept), max_dropped),
                             "rows": len(own)}
            dropped[step["folder"]] = found
    # Ancestry is the one panel that can fail on a run whose tables predate provenance columns,
    # and a page missing its attrition because ancestry raised is worse than one missing ancestry.
    try:
        graph = ancestry.collect(job_dir, fs=fs)
    except Exception:
        graph = None
    return {"job_dir": job_dir, "steps": steps, "manifest": record, "dropped": dropped,
            "tools_by_order": by_folder, "ancestry": graph}


def _provenance_panel(record) -> str:
    if not record:
        return ('<div class="panel"><h2>Provenance</h2><p class="warn">This run wrote no '
                'manifest, so what produced it was never recorded. Runs from BioPipelines 1.5.0 '
                'onward record it automatically.</p></div>')
    rows = [("BioPipelines", f"{record.get('biopipelines')} @ {record.get('commit')}"),
            ("Config variant", record.get("variant")),
            ("Scheduler", record.get("scheduler")),
            ("Python", record.get("python")),
            ("Recorded", record.get("created")),
            ("Identity", record.get("hash"))]
    environments = record.get("environments_resolved") or {}
    if environments:
        rows.append(("Environments", "; ".join(
            f"{name} @ {digest[:16]}" for name, digest in sorted(environments.items()))))
    else:
        rows.append(("Environments", "not recorded by this run"))
    cells = "".join(f"<dt>{_esc(k)}</dt><dd>{_esc(v)}</dd>" for k, v in rows)
    return f'<div class="panel"><h2>Provenance</h2><dl class="kv">{cells}</dl></div>'


def _attrition_panel(steps) -> str:
    if not steps:
        return ('<div class="panel"><h2>Attrition</h2><p class="warn">No step in this run wrote '
                'a stream map or a missing table, so there is no ID-level lineage to show.</p>'
                '</div>')
    widest = max((step.get("produced") or 0) for step in steps) or 1
    rows = []
    for step in steps:
        produced = step.get("produced") or 0
        missing = step.get("missing") or 0
        width = max(1, round(100 * produced / widest))
        streams = ", ".join(f"{name} {count}"
                            for name, count in sorted((step.get("streams") or {}).items()))
        loss = (f'<span class="loss">&minus;{missing}</span>' if missing else
                '<span class="nil">&mdash;</span>')
        rows.append(
            f'<tr><td class="folder">{_esc(step["folder"])}</td>'
            f'<td class="barcell"><span class="bar" style="width:{width}%"></span>'
            f'<span class="n">{produced}</span></td>'
            f'<td>{loss}</td><td class="streams">{_esc(streams) or "&mdash;"}</td></tr>')
    total = sum(step.get("missing") or 0 for step in steps)
    note = (f"{total} ID(s) were dropped across this campaign; every one is listed below with "
            f"the reason its own step recorded.") if total else "No step dropped any ID."
    return ('<div class="panel"><h2>Attrition</h2>'
            '<table class="grid"><thead><tr><th>Step</th><th>Produced</th><th>Dropped</th>'
            f'<th>Streams</th></tr></thead><tbody>{"".join(rows)}</tbody></table>'
            f'<p class="note">{_esc(note)}</p></div>')


def _parameters_panel(steps, data) -> str:
    """Driven by the manifest, not by the lineage steps.

    A step appears in the lineage only once it has written a stream map, so keying off that would blank this panel for a pipeline that was saved but never ran — exactly when someone is checking what it is configured to do.
    """
    record = data.get("manifest") or {}
    blocks = []
    for tool in record.get("tools") or []:
        label = f"{tool.get('order'):03d}_{tool.get('tool')}" if isinstance(
            tool.get("order"), int) else str(tool.get("tool"))
        resolved = (tool.get("parameters") or {}).get("resolved") or {}
        passed = set((tool.get("parameters") or {}).get("passed") or {})
        if not resolved:
            continue
        cells = "".join(
            f'<tr class="{"set" if name in passed else ""}"><td>{_esc(name)}</td>'
            f'<td>{_esc(repr(value))}</td></tr>'
            for name, value in sorted(resolved.items()))
        image = tool.get("container_image")
        meta = (f"v{tool.get('tool_version')} · "
                f"env {', '.join(tool.get('environments') or []) or '—'}"
                + (f" · image {os.path.basename(str(image))}" if image else ""))
        blocks.append(
            f'<details class="params"><summary>{_esc(label)} '
            f'<span class="meta">{_esc(meta)}</span></summary>'
            f'<table class="grid params"><thead><tr><th>Parameter</th><th>Value</th></tr></thead>'
            f'<tbody>{cells}</tbody></table></details>')
    if not blocks:
        return ""
    return ('<div class="panel"><h2>Parameters as they resolved</h2>'
            '<p class="note">Bold rows were set explicitly; the rest are the defaults that '
            'applied, which are usually what determined the output.</p>'
            + "".join(blocks) + "</div>")


def _ancestry_panel(graph, max_rows: int = MAX_ANCESTRY_ROWS) -> str:
    """Which design came from which input, and on what basis each link was established."""
    if not graph or not graph.get("edges"):
        return ('<div class="panel"><h2>Ancestry</h2><p class="warn">No parent links could be '
                'built for this run. Steps record them in the <code>&lt;axis&gt;.id</code> '
                'columns of <code>&lt;stream&gt;_map.csv</code>; a run whose tables predate '
                'those columns has none to read.</p></div>')

    edges, steps = graph["edges"], graph.get("steps") or []
    wiring = []
    for step in steps:
        sources = sorted({e["parent_step"] for e in edges if e["child_step"] == step["folder"]})
        wiring.append(
            f'<tr><td class="folder">{_esc(step["folder"])}</td>'
            f'<td class="n">{len(step["produced"])}</td>'
            f'<td class="streams">{_esc(", ".join(sources)) or "&mdash; (source step)"}</td></tr>')

    tiers = {}
    for edge in edges:
        key = edge["tier"].split(":")[0]
        tiers[key] = tiers.get(key, 0) + 1
    tier_note = ", ".join(f"{count} {name}" for name, count in sorted(tiers.items()))

    shown = edges[:max_rows]
    rows = "".join(
        f'<tr><td>{_esc(e["child_step"])}</td><td>{_esc(e["child"])}</td>'
        f'<td>{_esc(e["parent_step"])}</td><td>{_esc(e["parent"])}</td>'
        f'<td class="streams">{_esc(e["tier"])}</td></tr>' for e in shown)
    more = ("" if len(shown) >= len(edges) else
            f'<p class="note">Showing {len(shown)} of {len(edges)} links.</p>')

    unresolved = graph.get("unresolved") or []
    orphan = ("" if not unresolved else
              f'<p class="note">{len(unresolved)} ID(s) have no traceable parent and are '
              f'reported as roots rather than given an invented link.</p>')

    return ('<div class="panel"><h2>Ancestry</h2>'
            '<p class="note">Which design came from which input. <code>recorded</code> means a '
            'step wrote the parent id down; <code>matched</code> means it was resolved with the '
            'same id ladder the framework uses to wire one step into the next, and the tier it '
            'answered at is shown.</p>'
            '<table class="grid"><thead><tr><th>Step</th><th>IDs</th><th>Parents from</th>'
            f'</tr></thead><tbody>{"".join(wiring)}</tbody></table>'
            f'<p class="note">Links established by: {_esc(tier_note)}.</p>{orphan}'
            '<details class="params"><summary>Every parent link '
            f'<span class="meta">{len(edges)} links</span></summary>'
            '<table class="grid"><thead><tr><th>Step</th><th>ID</th><th>Parent step</th>'
            f'<th>Parent ID</th><th>Basis</th></tr></thead><tbody>{rows}</tbody></table>'
            f'{more}</details></div>')


def _dropped_panel(dropped) -> str:
    if not dropped:
        return ""
    blocks = []
    for folder, found in sorted(dropped.items()):
        # csv.reader, not split(","): `cause` is often `str(e)[:200]`, pandas quotes any field
        # containing a comma, and splitting shifts the row so the reason truncates at the comma.
        header = "".join(f"<th>{_esc(c)}</th>"
                         for c in next(csv.reader([found["header"]]), []))
        body = "".join(
            "<tr>" + "".join(f"<td>{_esc(c)}</td>" for c in row) + "</tr>"
            for row in csv.reader(found["lines"]))
        more = ("" if found["shown"] >= found["rows"] else
                f'<p class="note">Showing {found["shown"]} of {found["rows"]} rows.</p>')
        blocks.append(f'<details class="params"><summary>{_esc(folder)} '
                      f'<span class="meta">{found["rows"]} dropped</span></summary>'
                      f'<table class="grid"><thead><tr>{header}</tr></thead>'
                      f'<tbody>{body}</tbody></table>{more}</details>')
    return ('<div class="panel"><h2>Every ID that dropped out</h2>'
            '<p class="note">Read from each step\'s own <code>missing.csv</code>. A design '
            'absent from the final table is here, with the step and the reason that removed '
            'it.</p>' + "".join(blocks) + "</div>")


_CSS = """
:root { --bg:#f7f7f8; --card:#fff; --ink:#1c1c1e; --muted:#6b6b70; --line:#d9d9de;
        --accent:#2b6cb0; --ok:#2f855a; --bad:#c53030; }
* { box-sizing:border-box; }
body { margin:0; background:var(--bg); color:var(--ink);
       font:13px/1.5 -apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif; }
header.page { padding:20px 24px 12px; border-bottom:1px solid var(--line); background:var(--card); }
header.page h1 { margin:0 0 2px; font-size:19px; }
header.page .sub { color:var(--muted); }
.wrap { max-width:1180px; margin:0 auto; padding:16px 24px 64px; }
.panel { background:var(--card); border:1px solid var(--line); border-radius:8px;
         padding:14px 16px; margin:14px 0; }
.panel h2 { margin:0 0 10px; font-size:12px; text-transform:uppercase;
            letter-spacing:.06em; color:var(--muted); }
.kv { display:grid; grid-template-columns:150px 1fr; gap:4px 14px; margin:0; }
.kv dt { color:var(--muted); }
.kv dd { margin:0; font-family:ui-monospace,Consolas,monospace; word-break:break-all; }
table.grid { width:100%; border-collapse:collapse; }
table.grid th { text-align:left; font-weight:600; color:var(--muted); font-size:11px;
                text-transform:uppercase; letter-spacing:.04em;
                border-bottom:1px solid var(--line); padding:4px 8px 4px 0; }
table.grid td { padding:3px 8px 3px 0; border-bottom:1px solid #f0f0f2;
                font-family:ui-monospace,Consolas,monospace; }
td.folder { white-space:nowrap; }
td.barcell { width:55%; }
.bar { display:inline-block; height:9px; background:var(--accent); border-radius:2px;
       vertical-align:middle; }
.bar + .n { margin-left:8px; color:var(--muted); }
.loss { color:var(--bad); font-weight:600; }
.nil { color:var(--muted); }
.note { color:var(--muted); margin:10px 0 0; }
.warn { color:var(--bad); margin:0; }
details.params { border-top:1px solid #f0f0f2; padding:6px 0; }
details.params summary { cursor:pointer; font-family:ui-monospace,Consolas,monospace; }
details.params .meta { color:var(--muted); font-family:inherit; }
table.params tr.set td { font-weight:700; }
"""

_PAGE = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title><style>{css}</style></head>
<body><header class="page"><h1>{heading}</h1>
<div class="sub">{sub}</div></header>
<div class="wrap">{panels}</div></body></html>
"""


def render(data: Dict[str, Any], out_path: str) -> str:
    record = data.get("manifest") or {}
    job = record.get("job") or os.path.basename(str(data["job_dir"]).rstrip("/\\"))
    project = record.get("project") or ""
    heading = f"{project} / {job}" if project else job
    steps = data.get("steps") or []
    # Only steps that wrote a stream map: a step reporting drops but no stream produces 0, and
    # calling that the campaign's yield reads as total loss rather than as a step with no stream.
    with_streams = [step for step in steps if step.get("streams")]
    if with_streams:
        sub = (f"{len(steps)} step(s) reporting IDs · {with_streams[0]['produced']} produced at "
               f"{with_streams[0]['folder']}, {with_streams[-1]['produced']} at "
               f"{with_streams[-1]['folder']}")
    else:
        sub = f"{len(steps)} step(s) reporting IDs · no step wrote a stream map"

    panels = (_provenance_panel(data.get("manifest"))
              + _attrition_panel(steps)
              + _parameters_panel(steps, data)
              + _ancestry_panel(data.get("ancestry"))
              + _dropped_panel(data.get("dropped") or {}))
    page = _PAGE.format(title=_esc(f"Provenance — {heading}"), css=_CSS,
                        heading=_esc(heading), sub=_esc(sub), panels=panels)
    directory = os.path.dirname(os.path.abspath(out_path))
    if directory:
        os.makedirs(directory, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as handle:
        handle.write(page)
    return out_path


def write(job_dir: str, out_path: str, fs=None, max_dropped: int = MAX_DROPPED_ROWS) -> str:
    """Build and render in one call. The page lands locally even when the run is on a cluster."""
    return render(build(job_dir, fs=fs, max_dropped=max_dropped), out_path)
