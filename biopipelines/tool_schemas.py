"""Find a tool by what it produces or consumes, rather than by what it is called.

"Which tool gives me an RMSD" is the question tags answer worst and this answers exactly. A tag
vocabulary would need one entry per metric, hand-applied and drifting; the answer is already in
the code — every tool declares its streams, tables and columns in `get_output_files()`, and
`versions/extract_output_schemas.py` reads them with an AST interpreter that
`tests/test_doc_schemas.py` already gates the docs against.

So this is derived from verified source, not from new metadata: 86 tools, 445 distinct column
names, 22 stream names, 67 table names, none of it written by hand.

Matching is case-insensitive substring, because the useful query is the concept: `rmsd` must
find `ligand_rmsd` and `rmsd_to_ref`, and `ddg` must find `ddG_pred`.

Inputs come from the documented parameter list instead — a parameter typed `DataStream` or
`StandardizedOutput` is a stream the tool consumes, which is how you ask "what takes structures".
"""

import importlib.util
import re

from biopipelines.tool_docs import ROOT, collect, section

_EXTRACTOR = ROOT / "versions" / "extract_output_schemas.py"

_outputs_cache = None
_inputs_cache = None


def _load_extractor():
    spec = importlib.util.spec_from_file_location("extract_output_schemas", _EXTRACTOR)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def outputs():
    """{tool: {"streams": set, "tables": set, "columns": set}} straight from the sources."""
    global _outputs_cache
    if _outputs_cache is not None:
        return _outputs_cache

    try:
        schemas = _load_extractor().extract_all()
    except Exception:
        # A tool whose source cannot be parsed must not take the search down with it; the
        # caller gets an empty index and says so rather than silently returning no matches.
        return {}

    found = {}
    for name, schema in schemas.items():
        tables = schema.get("tables") or {}
        columns = set()
        for table in tables.values():
            columns.update(str(c) for c in (table.get("columns") or []) if c)
        found[name] = {"streams": set((schema.get("streams") or {}).keys()),
                       "tables": set(tables.keys()),
                       "columns": columns}
    _outputs_cache = found
    return found


def inputs():
    """{tool: set of input stream parameters} — parameters typed as a DataStream."""
    global _inputs_cache
    if _inputs_cache is not None:
        return _inputs_cache

    found = {}
    for entry in collect():
        body = section(entry["name"]) or ""
        block = re.search(r"\*\*Parameters\*\*:(.*?)(?=\*\*[A-Z]|\Z)", body, re.S)
        names = set()
        if block:
            for line in block.group(1).splitlines():
                m = re.match(r"\s*-\s*`([A-Za-z_][A-Za-z_0-9]*)`\s*:(.*)", line)
                if m and re.search(r"DataStream|StandardizedOutput", m.group(2)):
                    names.add(m.group(1))
        found[entry["name"]] = names
    _inputs_cache = found
    return found


def units(name):
    """The tool's `**Units.**` note, if its docs carry one.

    A column name alone is not enough to use a number: Boltz2's `affinity_pred_value` is
    log10(IC50) where lower is stronger, while GEMS's `pkd_pred` is a pKd where higher is
    stronger. Returning the name without the scale invites exactly the comparison that is wrong.
    """
    body = section(name) or ""
    # Stop at a blank line OR the next markdown block — a note written without a trailing blank
    # line otherwise swallows the table row that follows it.
    m = re.search(r"\*\*Units\.\*\*\s*(.+?)(?=\n\s*\n|\n\s*[-|*]\s|\n\s*\*\*|\Z)", body, re.S)
    if not m:
        return ""
    return re.sub(r"\s+", " ", m.group(1)).strip()


def search(query, where="outputs"):
    """Tools whose outputs (or inputs) mention `query`, with what matched.

    Returns [(tool, {kind: [matched names]})], ordered by tool name. `where` is "outputs",
    "inputs" or "any".
    """
    needle = str(query).strip().lower()
    if not needle:
        return []

    out, ins = outputs(), inputs()
    names = sorted(set(out) | set(ins))
    hits = []
    for name in names:
        matched = {}
        if where in ("outputs", "any"):
            for kind, values in out.get(name, {}).items():
                found = sorted(v for v in values if needle in v.lower())
                if found:
                    matched[kind] = found
        if where in ("inputs", "any"):
            found = sorted(v for v in ins.get(name, set()) if needle in v.lower())
            if found:
                matched["inputs"] = found
        if matched:
            hits.append((name, matched))
    return hits


def summarize(query, where="outputs", limit=25):
    """One screen naming the tools and the exact stream, table or column that matched."""
    hits = search(query, where=where)
    if not outputs():
        return ("Could not read the output schemas from source "
                "(versions/extract_output_schemas.py). Search by tool name instead.")
    if not hits:
        label = {"outputs": "output", "inputs": "input", "any": "input or output"}[where]
        return (f"No tool has an {label} matching {query!r}. Matching is a case-insensitive "
                "substring over stream, table and column names, so try a shorter fragment "
                "— 'rmsd' rather than 'ligand_rmsd_to_ref'.")

    lines = [f"{len(hits)} tool(s) with an output matching {query!r}:"
             if where == "outputs" else f"{len(hits)} tool(s) matching {query!r} ({where}):", ""]
    for name, matched in hits[:limit]:
        detail = "; ".join(f"{kind}: {', '.join(values)}" for kind, values in sorted(matched.items()))
        lines.append(f"  {name:<24} {detail}")
        note = units(name)
        if note:
            lines.append(f"  {'':<24} units: {note[:300]}")
    if len(hits) > limit:
        lines.append(f"  … {len(hits) - limit} more")
    lines.append("")
    lines.append("These are declared in the tools' own source, not hand-maintained tags. Where a "
                 "`units:` line is shown, it is the tool's own documentation — read it before "
                 "comparing numbers across tools, because the scales differ and several run in "
                 "opposite directions. No units line means the scale is undocumented; open the "
                 "tool's section with bp_tools(name=...) rather than assuming.")
    return "\n".join(lines)
