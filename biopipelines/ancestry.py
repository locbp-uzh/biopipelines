"""Which design came from which input: the per-ID parent graph of a campaign.

`lineage` counts what each step produced and dropped. This answers the other half — given
`9_Panda_1`, which Boltz2 complex, which MPNN sequence, which RFdiffusion backbone and which
PDB entry produced it. That is the question a referee asks, and the one an agent that generates
glue code cannot answer at all, because the answer is a property of having executed a typed
graph rather than of the code that was written.

The rules, and why each tier exists, are in `docs/design/id-ancestry.md`. In short: a
`<axis>.id` cell naming a concrete upstream id is taken as recorded; where the framework wrote
a combinatorial pattern instead, or pruned the column because the id already carries its
parent, the child is matched against the upstream ids with the same ladder the framework uses
at runtime to wire the streams together.

Parent search is scoped by `RunTime/pipeline_graph.json`, which records the declared wiring.
Without that scope a Boltz2 complex matches its grandparent as readily as its parent.
"""

import json

from biopipelines.remote import LocalFS

MAP_SUFFIX = "_map.csv"
GRAPH = "RunTime/pipeline_graph.json"
# Tables that carry stream content under their own name rather than `<stream>_map.csv`.
CONTENT_TABLES = ("sequences.csv", "confidence.csv")


def _step_key(folder):
    index = folder.split("_", 1)[0]
    return (int(index) if index.isdigit() else 0, folder)


def _is_pattern(value):
    """Combinatorial declarations like `3kzy_A_<1..100>_<1..2>` name no single id."""
    return "<" in value


def _read_graph(job_dir, fs):
    """Declared step wiring: {step folder: [upstream step folders]}. Empty when unrecorded."""
    try:
        raw = json.loads(fs.read_text(fs.join(job_dir, *GRAPH.split("/"))))
    except Exception:
        return {}, {}
    by_order = {}
    for step in raw.get("steps", []):
        folder = str(step.get("output_folder", "")).replace("\\", "/").rstrip("/").rsplit("/", 1)[-1]
        if folder:
            by_order[step.get("execution_order")] = folder
    upstream = {folder: [] for folder in by_order.values()}
    for edge in raw.get("edges", []):
        child, parent = by_order.get(edge.get("to_step")), by_order.get(edge.get("from_step"))
        if child and parent and parent not in upstream.get(child, []):
            upstream.setdefault(child, []).append(parent)
    tools = {by_order[s.get("execution_order")]: s.get("tool")
             for s in raw.get("steps", []) if s.get("execution_order") in by_order}
    return upstream, tools


def _read_steps(job_dir, fs):
    """Per step folder: the ids it produced, and the provenance cells recorded against each."""
    tables = fs.id_columns(job_dir, [f"*{MAP_SUFFIX}", *CONTENT_TABLES])
    steps = {}
    for path, rows in tables.items():
        parts = path.replace("\\", "/").split("/")
        if len(parts) < 2:
            continue
        folder = parts[0]
        entry = steps.setdefault(folder, {"folder": folder, "produced": [], "provenance": {},
                                          "streams": {}})
        stream = parts[-1][: -len(MAP_SUFFIX)] if parts[-1].endswith(MAP_SUFFIX) else parts[-1]
        entry["streams"][stream] = len({(r.get("id") or "").strip() for r in rows} - {""})
        for row in rows:
            rid = (row.get("id") or "").strip()
            if not rid:
                continue
            if rid not in entry["provenance"]:
                entry["produced"].append(rid)
                entry["provenance"][rid] = []
            for col, val in row.items():
                val = (val or "").strip()
                if col == "id" or not col.endswith(".id") or not val or val == rid:
                    continue
                if val not in entry["provenance"][rid]:
                    entry["provenance"][rid].append(val)
    return sorted(steps.values(), key=lambda s: _step_key(s["folder"]))


def _matched_edges(child_ids, parent_step):
    """The matcher's own answer for a whole step at once, as {child: (parent, tier)}."""
    from biopipelines.id_map_utils import get_mapped_ids_with_tiers
    if not child_ids or not parent_step["produced"]:
        return {}
    try:
        return get_mapped_ids_with_tiers(list(child_ids), list(parent_step["produced"]))
    except Exception:
        return {}


def collect(job_dir, fs=None):
    """The campaign's per-ID parent graph.

    Returns `{"steps": [...], "edges": [...], "unresolved": [...]}`. An edge carries the tier
    that established it, because an edge whose basis a reader cannot see is one they cannot
    audit.
    """
    fs = fs or LocalFS()
    steps = _read_steps(job_dir, fs)
    by_folder = {s["folder"]: s for s in steps}
    upstream, tools = _read_graph(job_dir, fs)
    order = {s["folder"]: i for i, s in enumerate(steps)}

    # Whether the RUN recorded its wiring, not whether this step has an upstream: a source step
    # legitimately has none, and treating that as "unrecorded" sends it scanning every earlier
    # step, which is how `002_Sequence` acquired `001_PDB` as a parent it never had.
    recorded_wiring = bool(upstream)

    edges, parents_of = [], {}
    for step in steps:
        folder = step["folder"]
        if recorded_wiring:
            candidates = [by_folder[p] for p in upstream.get(folder, []) if p in by_folder]
        else:
            # An older run: every earlier step is a candidate, and the transitive drop below
            # removes the grandparents that admits.
            candidates = [s for s in steps if order[s["folder"]] < order[folder]]

        accepted = {rid: [] for rid in step["produced"]}
        for rid in step["produced"]:
            for value in step["provenance"][rid]:
                if _is_pattern(value):
                    continue
                for cand in candidates:
                    if value in cand["provenance"]:
                        accepted[rid].append((cand["folder"], value, "recorded"))
                        break

        for cand in candidates:
            unanswered = [rid for rid in step["produced"]
                          if not any(e[0] == cand["folder"] for e in accepted[rid])]
            for rid, (parent, tier) in _matched_edges(unanswered, cand).items():
                # The matcher reports an unmatched id as (None, tier). Accepting that builds an
                # edge to nothing, which is precisely the invented link this module refuses.
                if not parent:
                    continue
                accepted[rid].append((cand["folder"], parent, f"matched:{tier}"))

        for rid, found in accepted.items():
            kept = found if recorded_wiring else _drop_transitive(found, parents_of)
            for parent_folder, parent_id, tier in kept:
                edges.append({"child": rid, "child_step": folder, "parent": parent_id,
                              "parent_step": parent_folder, "tier": tier})
            parents_of[(folder, rid)] = [(p, i) for p, i, _ in kept]

    produced_total = [(s["folder"], rid) for s in steps for rid in s["produced"]]
    first = steps[0]["folder"] if steps else None
    # A step the recorded wiring gives no upstream is a source; its ids are roots, not failures.
    sources = {f for f in by_folder if recorded_wiring and f in upstream and not upstream[f]}
    unresolved = [{"step": f, "id": i} for f, i in produced_total
                  if not parents_of.get((f, i)) and f != first and f not in sources]
    return {"steps": steps, "edges": edges, "unresolved": unresolved, "tools": tools}


def _drop_transitive(found, parents_of):
    """Remove a candidate parent that is already an ancestor of another candidate.

    Only needed where the wiring was not recorded: matching `3kzy_A_100_1+tag` against every
    earlier step finds the RFdiffusion backbone as readily as the MPNN sequence, and reporting
    both as parents makes the graph say the campaign had a shortcut it did not have.
    """
    if len(found) < 2:
        return found
    kept = []
    for cand in found:
        if not any(_reaches(other, cand, parents_of) for other in found if other is not cand):
            kept.append(cand)
    return kept


def _reaches(start, target, parents_of, depth=12):
    """Is `target` an ancestor of `start`?"""
    frontier, seen = [(start[0], start[1])], set()
    for _ in range(depth):
        nxt = []
        for node in frontier:
            if node in seen:
                continue
            seen.add(node)
            for parent in parents_of.get(node, []):
                if parent == (target[0], target[1]):
                    return True
                nxt.append(parent)
        if not nxt:
            return False
        frontier = nxt
    return False


def trace(graph, target, max_depth=20, step=None):
    """Every path from `target` back to the roots, as a nested structure.

    Multi-parent by construction: a complex has a designed chain and a tag, and a renderer that
    assumes one parent drops exactly the half of a binder campaign a reviewer wants to check.

    Nodes are (step, id): an id that passes through a step unchanged is a different node at each step, and keying on the id alone flattened that chain into a shortcut. Without `step`, the walk starts from the latest step that produced `target`.
    """
    by_child = {}
    for edge in graph["edges"]:
        by_child.setdefault((edge["child_step"], edge["child"]), []).append(edge)

    def walk(node, depth, seen):
        if depth >= max_depth or node in seen:
            return []
        branches = []
        for edge in by_child.get(node, []):
            parent = (edge["parent_step"], edge["parent"])
            branches.append({"id": edge["parent"], "step": edge["parent_step"],
                             "tier": edge["tier"],
                             "parents": walk(parent, depth + 1, seen | {node})})
        return branches

    producers = [s["folder"] for s in graph["steps"] if target in s["provenance"]]
    if step is not None:
        producers = [f for f in producers if f == step]
    if not producers:
        return None
    start = max(producers, key=_step_key)
    return {"id": target, "step": start, "tier": None, "parents": walk((start, target), 0, set())}


def render_trace(node, indent="", last=True, lines=None):
    """The walk-back as a tree a person reads without a viewer."""
    lines = [] if lines is None else lines
    if not indent:
        lines.append(f"{node['id']}   [{node['step'] or 'unknown step'}]")
    else:
        tier = f"  [{node['tier']}]" if node["tier"] else ""
        lines.append(f"{indent}{'└─ ' if last else '├─ '}"
                     f"{node['step'] or '?':<24} {node['id']}{tier}")
    children = node["parents"]
    for i, child in enumerate(children):
        deeper = indent + ("" if not indent else ("   " if last else "│  "))
        render_trace(child, deeper or " ", i == len(children) - 1, lines)
    return lines


def summarize(graph, job=None):
    """The graph in one screen: how it was established, and what it could not establish."""
    edges, steps = graph["edges"], graph["steps"]
    if not steps:
        return (f"No map tables in {job or 'this run'}, so there is no ancestry to build. "
                "Steps record parents in `<stream>_map.csv`.")

    tiers = {}
    for edge in edges:
        key = edge["tier"].split(":")[0]
        tiers[key] = tiers.get(key, 0) + 1
    produced = sum(max(s.get("streams", {}).values(), default=len(s["produced"])) for s in steps)

    lines = [f"{job or 'run'}: {len(edges)} parent link(s) over {produced} ID(s) "
             f"in {len(steps)} step(s)", ""]
    for step in steps:
        out = [e for e in edges if e["child_step"] == step["folder"]]
        sources = sorted({e["parent_step"] for e in out})
        # The largest stream, not the union: a predictor's msas/compounds streams echo its inputs' ids.
        count = max(step.get("streams", {}).values(), default=len(step["produced"]))
        detail = f"{count:>5} ID(s)"
        if sources:
            detail += f"  <- {', '.join(sources)}"
        lines.append(f"  {step['folder']:<34} {detail}")

    lines.append("")
    lines.append("Established by: " + ", ".join(f"{k} {v}" for k, v in sorted(tiers.items()))
                 if tiers else "No parent links were established.")
    if graph["unresolved"]:
        n = len(graph["unresolved"])
        first = graph["unresolved"][0]
        lines.append(f"{n} ID(s) have no traceable parent, the first being "
                     f"{first['id']} in {first['step']}. They are reported as roots rather "
                     f"than given an invented edge.")
    return "\n".join(lines)


def to_rows(graph):
    """The graph as one row per edge, for `ancestry.csv`."""
    return [{"child_step": e["child_step"], "child": e["child"],
             "parent_step": e["parent_step"], "parent": e["parent"], "tier": e["tier"]}
            for e in graph["edges"]]
