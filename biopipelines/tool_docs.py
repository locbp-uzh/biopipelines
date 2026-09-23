"""Read the human tool docs the way an agent needs them: an index, then one tool at a time.

`docs/tool/*.md` is written for people — nine category files, ~78k tokens, one `###` (or `##`)
section per tool — and stays that way. An agent cannot afford to read it to find out what
exists, so this module serves two views over the same files: a one-line-per-tool index
(~4k tokens) and a single tool's section on demand.

Both the index generator (`skills/biopipelines/build_tool_index.py`) and the MCP server
(`biopipelines/mcp_server.py`) go through here, so the two can never disagree about what a
tool is or where its documentation lives.
"""

import io
import pathlib
import re

from biopipelines import tool_tags

ROOT = pathlib.Path(__file__).resolve().parent.parent
DOCS = ROOT / "docs" / "tool"

BLURB_CHARS = 150
PLATFORMS = {"HPC x86-64 ok": "x86-64", "HPC aarch64 ok": "aarch64", "Colab ok": "Colab"}

# Variants with no README row of their own; they run exactly what their parent runs.
BADGE_ALIAS = {"SolubleMPNN": "ProteinMPNN"}


def _read(rel):
    return io.open(ROOT / rel, encoding="utf-8").read()


def _tool_level(text):
    # `data_management.md` puts tools at `##`; every other file uses `###`. Read the level off
    # the file rather than normalizing the human docs, which stay as their authors wrote them.
    return "### " if "\n### " in text else "## "


def slug(name):
    # GitHub anchors keep underscores; only other punctuation becomes a hyphen.
    return re.sub(r"[^a-z0-9_]+", "-", name.lower()).strip("-")


def _names_and_anchor(heading):
    """`### Load / LoadMultiple {#load}` — one section, two callable names, an explicit anchor."""
    anchor = None
    m = re.search(r"\{#([^}]+)\}", heading)
    if m:
        anchor, heading = m.group(1), heading[: m.start()]
    names = [n.strip() for n in heading.split("/")]
    return [n for n in names if re.fullmatch(r"[A-Za-z0-9_]+", n)], anchor


def _blurb(body):
    """First prose sentence of a tool section, stripped of markdown noise."""
    for line in body.splitlines():
        line = line.strip()
        if not line or line.startswith(("#", "```", "|", "-", "*", ">", "!")):
            continue
        line = re.sub(r"`([^`]*)`", r"\1", line)
        line = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", line)
        line = re.sub(r"\*\*([^*]*)\*\*", r"\1", line)
        line = re.sub(r"\s+", " ", line)
        if line.startswith(("Streams", "Tables", "Example", "Parameters")):
            continue
        cut = line.find(". ")
        if cut != -1 and cut < BLURB_CHARS:
            line = line[: cut + 1]
        return line[:BLURB_CHARS].rstrip()
    return ""


def _versions():
    """{ToolName: version} from the generated registry index."""
    rows = re.findall(r"^\|\s*\d+\s*\|\s*([A-Za-z0-9_]+)\s*\|[^|]*\|\s*([0-9.]+)\s*\|",
                      _read("docs/tool_index.md"), re.M)
    return dict(rows)


def _badges():
    """{ToolName: (hardware, [platforms])} off the README table's `alt` attributes.

    The badges are presentation markup — an agent should never have to parse them, which is
    the whole reason they are lifted here. A row can carry both CPU and GPU.
    """
    found = {}
    for cell in re.split(r"<td><sub><b>", _read("README.md"))[1:]:
        alts = re.findall(r'alt="([^"]+)"', cell)
        hw = "/".join([a for a in ("GPU", "CPU") if a in alts])
        # Fixed order, not the README's badge order, so the output is stable.
        plat = [label for key, label in PLATFORMS.items() if key in alts]
        for name in cell.split("</b>")[0].split("/"):
            found[name.strip()] = (hw, plat)
    for variant, parent in BADGE_ALIAS.items():
        if parent in found:
            found[variant] = found[parent]
    return found


def _upstream():
    """{ToolName: {"repo": url, "paper": url}} from the README's badge links.

    An agent that knows a tool wraps `github.com/jwohlwend/boltz` can check a capability at its
    source; without it, verifying anything means guessing which project is meant. BP-native
    tools carry neither, which is itself the useful signal that there is no upstream.
    """
    found = {}
    for cell in re.split(r"<td><sub><b>", _read("README.md"))[1:]:
        links = {}
        for m in re.finditer(r'<a href="([^"]+)"[^>]*>.*?alt="(repo|paper)"', cell, re.S):
            links.setdefault(m.group(2), m.group(1))
        for name in cell.split("</b>")[0].split("/"):
            found[name.strip()] = links
    for variant, parent in BADGE_ALIAS.items():
        if parent in found:
            found[variant] = found[parent]
    return found


def collect():
    """Every callable tool, with the metadata the index line needs."""
    versions, badges, upstream = _versions(), _badges(), _upstream()
    entries = []
    for path in sorted(DOCS.glob("*.md")):
        text = io.open(path, encoding="utf-8").read()
        category = text.splitlines()[0].lstrip("# ").strip()
        level = _tool_level(text)
        for part in text.split("\n" + level)[1:]:
            heading, _, body = part.partition("\n")
            names, anchor = _names_and_anchor(heading)
            blurb = _blurb(body)
            tags = tool_tags.parse_tags_line(body) or []
            for name in names:
                hw, plat = badges.get(name, ("", []))
                links = upstream.get(name, {})
                entries.append({
                    "name": name,
                    "tags": tags,
                    "category": category,
                    "version": versions.get(name, ""),
                    "hardware": hw,
                    "platforms": plat,
                    "blurb": blurb,
                    "repo": links.get("repo", ""),
                    "paper": links.get("paper", ""),
                    "path": f"docs/tool/{path.name}#{anchor or slug(name)}",
                })
    return entries


def section(name):
    """The full markdown section for one tool, verbatim from the human docs.

    Returns None when the name is not a documented tool, so callers can offer suggestions
    rather than inventing a signature.
    """
    for path in sorted(DOCS.glob("*.md")):
        text = io.open(path, encoding="utf-8").read()
        level = _tool_level(text)
        for part in text.split("\n" + level)[1:]:
            heading, _, body = part.partition("\n")
            names, _anchor = _names_and_anchor(heading)
            if name in names:
                return (level + heading + "\n" + body).rstrip()
    return None


def suggest(name, limit=5):
    """Documented tool names closest to `name`, for an unknown-tool error worth reading."""
    import difflib
    known = [e["name"] for e in collect()]
    close = difflib.get_close_matches(name, known, n=limit, cutoff=0.5)
    lowered = name.lower()
    substr = [k for k in known if lowered in k.lower() and k not in close]
    return (close + substr)[:limit]


def render(entries):
    """The agent-facing index: one line per tool, grouped by the docs' own categories."""
    by_cat = {}
    for e in entries:
        by_cat.setdefault(e["category"], []).append(e)

    out = [
        "# BioPipelines tool index (agent-facing)",
        "",
        f"Generated by `skills/biopipelines/build_tool_index.py` — do not edit by hand. "
        f"{len(entries)} callable tools.",
        "",
        "Pick the tools you need here, then read **only** those sections from the linked file "
        "(or call `bp_tools` with the tool name). Never read a whole category file to find a tool, "
        "and never guess a name or a signature — every entry below points at the authoritative one.",
        "",
        "Each entry reads `Name vVERSION · HARDWARE · PLATFORMS — purpose  docs pointer` and, "
        "for a wrapped tool, its upstream repo and paper. A tool with neither is BioPipelines' "
        "own. Use the repo when you need to check a capability the docs here do not mention. "
        "Platforms are the ones the tool has been **verified** on (`x86-64` / `aarch64` = HPC, "
        "plus `Colab`). A platform missing from that list means *not verified*, which is not the "
        "same as *unsupported*: several tools run fine on a cluster without carrying the badge. "
        "Prefer a tool whose platforms cover the user's mode; when only the verification is "
        "missing, say so and ask rather than refusing the tool.",
        "",
        "Each line ends with the tool's tags — a closed 38-term vocabulary in four facets "
        "(action, subject, readout, capability) documented in `docs/tool_tags.md`. Grep them, "
        "or filter with `bp_tools(tags=[...], exclude=[...])`: within a facet the terms are OR-ed, across facets AND-ed. Tags are the intent layer; for a named metric ask "
        "`bp_tools(outputs=\"rmsd\")` instead, which reads stream, table and column names "
        "straight from the tools' source.",
        "",
    ]
    for cat in sorted(by_cat):
        out.append(f"## {cat}")
        out.append("")
        for e in sorted(by_cat[cat], key=lambda x: x["name"]):
            out.append(line(e))
        out.append("")
    return "\n".join(out)


def line(e):
    """One tool's index line — the same shape wherever a list of tools is printed."""
    meta = [f"v{e['version']}" if e["version"] else "", e["hardware"], ", ".join(e["platforms"])]
    head = " · ".join(m for m in meta if m)
    links = " ".join(f"[{kind}]({e[kind]})" for kind in ("repo", "paper") if e.get(kind))
    tail = f" `{e['path']}`" + (f" {links}" if links else "")
    tags = f" · tags: {', '.join(e['tags'])}" if e.get("tags") else ""
    return f"- **{e['name']}** {head} — {e['blurb']}{tail}{tags}"


def render_selection(entries, header):
    """A filtered list of tools, grouped by category, in the index's own line format."""
    by_cat = {}
    for e in entries:
        by_cat.setdefault(e["category"], []).append(e)
    out = [header, ""]
    for cat in sorted(by_cat):
        out.append(f"## {cat}")
        out.append("")
        for e in sorted(by_cat[cat], key=lambda x: x["name"]):
            out.append(line(e))
        out.append("")
    return "\n".join(out)
