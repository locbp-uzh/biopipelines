"""Renderer for a whole pipeline: one self-contained HTML page laid out as the pipeline's shape, each node carrying that step's own rendered tool output.

Node bodies come from the same renderers the notebook display uses, dispatched the same way (stream name -> format -> lowercased format -> _default). The page has to open from a file:// URL on a laptop after being copied off a cluster, so a renderer that reaches for a browser library on a CDN is satisfied in a three-link chain: the vendored copy under renderers/vendor/ is inlined once per page, else the remote tag is kept when allow_external=True, else the node falls back to the self-contained metadata table plus a note naming what was dropped.
"""

import html as html_module
import json
import os
import re


# A node body matching this renders blank off-network, so it is replaced rather than shipped.
_EXTERNAL_REF_RE = re.compile(r"""(?:src|href)\s*=\s*["']?\s*(?:https?:)?//""", re.I)

# Renderers mark a browser-library <script> with data-bp-lib so page assembly can swap it for a vendored copy.
_LIB_TAG_RE = re.compile(
    r"""<script\b[^>]*\bdata-bp-lib\s*=\s*["']([A-Za-z0-9_.-]+)["'][^>]*>\s*</script\s*>""",
    re.I,
)

_VENDOR_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "vendor")

_VENDORED_LIBS = {"3dmol": ("3Dmol-min.js", "3Dmol.js 2.5.2")}


# Without these a 40-step pipeline of base64 plots produces a page no browser opens.
_MAX_NODE_BODY_BYTES = 1_500_000
_MAX_TOTAL_BODY_BYTES = 12_000_000

_STATUS_LABELS = {
    "completed": "completed",
    "failed": "FAILED",
    "pending": "pending",
}


def _esc(value):
    return html_module.escape(str(value if value is not None else ""))


def _rel(path, root):
    """Display a path relative to the pipeline output folder when it is under it."""
    if not path:
        return ""
    if root and path.startswith(root):
        return os.path.relpath(path, root).replace("\\", "/")
    return path


def _fmt_resources(resources):
    if not resources:
        return "(none)"
    # The scheduler gates every GPU directive on the spec, so a count without one allocates nothing.
    declined = "gpu" in resources and str(resources["gpu"] or "none").strip().lower() in ("none", "")
    parts = []
    for key, value in resources.items():
        if declined and key in ("gpu", "gpus"):
            if key == "gpu":
                parts.append("gpu=none")
            continue
        if value is None or value == {} or value == "":
            continue
        if isinstance(value, dict):
            inner = ", ".join(f"{k}={v}" for k, v in value.items() if v is not None)
            if inner:
                parts.append(f"{key}({inner})")
        else:
            parts.append(f"{key}={value}")
    return ", ".join(parts) or "(none)"


# ── vendored browser libraries ──────────────────────────────────────────────

def _vendored_lib(name):
    """The on-disk path of the vendored copy of library ``name``, or None."""
    entry = _VENDORED_LIBS.get(name.lower())
    if not entry:
        return None
    path = os.path.join(_VENDOR_DIR, entry[0])
    return path if os.path.isfile(path) else None


def _claim_libraries(fragment, libs):
    """Strip every marked library tag that a vendored copy can satisfy, recording it in ``libs`` for one page-level embed.

    Tags left in place are the ones with no vendored copy, so the caller's external-reference check still decides their fate.
    """
    def _swap(match):
        name = match.group(1).lower()
        if _vendored_lib(name):
            libs.add(name)
            return ""
        return match.group(0)

    return _LIB_TAG_RE.sub(_swap, fragment)


def _library_scripts(libs):
    """One inline <script> per embedded library, so a page with many structure nodes carries the library once."""
    parts = []
    for name in sorted(libs):
        path = _vendored_lib(name)
        if not path:
            continue
        try:
            with open(path, "r", encoding="utf-8") as f:
                source = f.read()
        except Exception:
            continue
        # An inline script ends at the first "</script"; a library containing one would truncate the page.
        source = re.sub(r"</(?=\s*script)", lambda _m: "<\\/", source, flags=re.I)
        label = _VENDORED_LIBS[name][1]
        parts.append(
            f'<script data-bp-lib-embedded="{_esc(name)}">/* {label} '
            f"- vendored, see renderers/vendor/README.md */\n{source}\n</script>"
        )
    return "\n".join(parts)


# ── output reconstruction ────────────────────────────────────────────────────

def _looks_like_stream(value):
    return (isinstance(value, dict)
            and "format" in value
            and "ids" in value
            and "_type" not in value)


def _rebuild_output(step):
    """Rebuild a StandardizedOutput from the step's exported ToolOutputs JSON.

    Returns (output, exported_json, problem). The page is generated from the
    same on-disk artifacts before and after a run, so both cases go through
    here and the only difference is whether the files the streams name exist.
    """
    from biopipelines.base_config import StandardizedOutput
    from biopipelines.datastream import DataStream

    path = step.get("tool_outputs_json")
    if not path or not os.path.isfile(path):
        return None, None, f"no exported output inventory at {path or '(unset)'}"
    try:
        with open(path, encoding="utf-8") as f:
            exported = json.load(f)
    except Exception as e:
        return None, None, f"could not read {os.path.basename(path)}: {e}"

    structure = exported.get("output_structure") or {}
    rebuilt = {}
    for key, value in structure.items():
        if _looks_like_stream(value):
            try:
                stream = DataStream.from_dict(value)
            except Exception as e:
                rebuilt[key] = value
                continue
            # Only the runtime flag makes lazy ids read the map_table instead of expanding the config-time prefix.
            if stream.map_table and os.path.isfile(stream.map_table):
                stream._runtime_mode = True
            rebuilt[key] = stream
        else:
            rebuilt[key] = value

    try:
        return StandardizedOutput(rebuilt), exported, None
    except Exception as e:
        return None, exported, f"could not rebuild output object: {e}"


# ── node body: dispatch to the configured renderers ─────────────────────────

def _resolve_renderers_config(recorded_variant, notes):
    """The configured renderers, resolved once for the whole page.

    A page regenerated in a bare process resolves whatever config variant that process auto-detects, which need not be the one the pipeline ran under and may declare no renderers at all. Only then is the recorded variant activated, so a correct ambient config is never overridden.
    """
    from biopipelines.base_config import StandardizedOutput

    probe = StandardizedOutput({})
    try:
        config = probe._get_renderers_config()
    except Exception as e:
        notes.append(f"renderer config unavailable ({e}); nodes show metadata only")
        return {}
    if config.get("streams") or not recorded_variant:
        return config
    try:
        from biopipelines.config_manager import ConfigManager
        ConfigManager(variant=recorded_variant)
        config = probe._get_renderers_config()
    except Exception as e:
        notes.append(
            f"the active config declares no renderers and variant "
            f"{recorded_variant!r} could not be loaded ({e})"
        )
    if not config.get("streams"):
        notes.append("no stream renderers are configured; nodes show metadata only")
    return config


def _render_stream_bodies(output, config, notes, allow_external=False, libs=None):
    """Render every non-empty stream and table of ``output`` via the configured renderers.

    ``libs`` collects the browser libraries whose vendored copy satisfied a renderer, for the caller to embed once for the whole page.
    """
    from biopipelines.base_config import StandardizedOutput
    from biopipelines.datastream import DataStream

    if libs is None:
        libs = set()

    stream_renderers = config.get("streams", {})
    table_renderers = config.get("tables", {})
    metadata_script = stream_renderers.get("_metadata")
    default_script = stream_renderers.get("_default")

    def _run(script, payload):
        fn = StandardizedOutput._load_render_fn(script)
        return fn(payload, output) or ""

    parts = []

    for name, stream in output.streams.items():
        if not isinstance(stream, DataStream):
            continue
        if len(stream) == 0:
            notes.append(f"stream '{name}' declares no items")
            continue

        script = (stream_renderers.get(name)
                  or stream_renderers.get(stream.format)
                  or stream_renderers.get(stream.format.lower())
                  or default_script)

        rendered = ""
        if metadata_script:
            try:
                rendered += _run(metadata_script, stream)
            except Exception as e:
                notes.append(f"metadata renderer failed for '{name}': {e}")
        if not script:
            if not metadata_script:
                notes.append(f"no renderer configured for stream '{name}' ({stream.format})")
        elif script != metadata_script:
            try:
                specialized = _run(script, stream)
            except Exception as e:
                specialized = ""
                notes.append(
                    f"{os.path.basename(script)} failed for '{name}': {e}"
                )
            if specialized:
                specialized = _claim_libraries(specialized, libs)
            if (specialized and not allow_external
                    and _EXTERNAL_REF_RE.search(specialized)):
                notes.append(
                    f"{os.path.basename(script)} for '{name}' needs a CDN, which a "
                    f"file:// page cannot fetch, and no vendored copy is present in "
                    f"renderers/vendor — its metadata table is shown instead"
                )
            elif specialized:
                rendered += specialized
        parts.append(rendered)

    tables = getattr(output.tables, "_tables", {}) or {}
    for name, info in tables.items():
        script = table_renderers.get(name) or table_renderers.get("_default")
        if not script:
            notes.append(f"no renderer configured for table '{name}'")
            continue
        try:
            rendered = _run(script, info)
        except Exception as e:
            notes.append(f"{os.path.basename(script)} failed for table '{name}': {e}")
            continue
        if rendered:
            rendered = _claim_libraries(rendered, libs)
        if rendered and not allow_external and _EXTERNAL_REF_RE.search(rendered):
            notes.append(f"{os.path.basename(script)} for table '{name}' needs a CDN; skipped")
            continue
        parts.append(rendered)

    return "\n".join(p for p in parts if p)


# ── per-step status ─────────────────────────────────────────────────────────

def _step_status(step):
    if step.get("failed_marker") and os.path.isfile(step["failed_marker"]):
        return "failed"
    if step.get("completed_marker") and os.path.isfile(step["completed_marker"]):
        return "completed"
    return "pending"


# ── page assembly ───────────────────────────────────────────────────────────

_CSS = """
:root {
  --bg: #f7f7f8; --card: #ffffff; --ink: #1c1c1e; --muted: #6b6b70;
  --line: #d9d9de; --accent: #2b6cb0; --ok: #2f855a; --bad: #c53030;
  --pending: #8a8a90; --wire: #2b6cb0; --order: #b9b9c0;
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--bg); color: var(--ink);
       font: 13px/1.5 -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; }
a { color: var(--accent); }
header.page { padding: 20px 24px 12px; border-bottom: 1px solid var(--line); background: var(--card); }
header.page h1 { margin: 0 0 2px; font-size: 19px; }
header.page .sub { color: var(--muted); }
.wrap { max-width: 1180px; margin: 0 auto; padding: 16px 24px 64px; }
.panel { background: var(--card); border: 1px solid var(--line); border-radius: 8px;
         padding: 12px 14px; margin: 14px 0; }
.panel h2 { margin: 0 0 8px; font-size: 13px; text-transform: uppercase;
            letter-spacing: .06em; color: var(--muted); }
.kv { display: grid; grid-template-columns: max-content 1fr; gap: 2px 14px; }
.kv dt { color: var(--muted); }
.kv dd { margin: 0; font-family: ui-monospace, Consolas, monospace; word-break: break-all; }
.banner { border-left: 4px solid var(--accent); background: #eef4fb; padding: 10px 12px;
          border-radius: 4px; margin: 14px 0; }
.banner.bad { border-left-color: var(--bad); background: #fdeeee; }
.legend { display: flex; flex-wrap: wrap; gap: 16px; align-items: center; color: var(--muted); }
.legend .swatch { display: inline-block; width: 26px; height: 0; vertical-align: middle;
                  border-top: 2px solid var(--wire); margin-right: 6px; }
.legend .swatch.order { border-top: 2px dashed var(--order); }
.controls { display: flex; flex-wrap: wrap; gap: 10px; align-items: center; }
.controls button { font: inherit; padding: 4px 10px; border: 1px solid var(--line);
                   background: var(--card); border-radius: 5px; cursor: pointer; }
.controls input[type=search] { font: inherit; padding: 4px 8px; border: 1px solid var(--line);
                               border-radius: 5px; min-width: 200px; }
#flow { position: relative; }
#wires { position: absolute; inset: 0; width: 100%; height: 100%;
         pointer-events: none; overflow: visible; }
.lane { position: relative; margin: 18px 0; }
.lane > .lane-head { margin-left: 78px; padding: 6px 10px; border-radius: 6px 6px 0 0;
                     background: #ececf1; border: 1px solid var(--line); border-bottom: none;
                     font-weight: 600; }
.lane > .lane-head .lane-meta { font-weight: 400; color: var(--muted); }
.lane > .nodes { margin-left: 78px; border: 1px solid var(--line); border-radius: 0 0 6px 6px;
                 background: rgba(255,255,255,.55); padding: 8px; }
details.node { background: var(--card); border: 1px solid var(--line); border-radius: 6px;
               margin: 6px 0; }
details.node > summary { cursor: pointer; padding: 7px 10px; display: flex; gap: 10px;
                         align-items: center; flex-wrap: wrap; list-style: none; }
details.node > summary::-webkit-details-marker { display: none; }
details.node > summary::before { content: "\\25B8"; color: var(--muted); }
details.node[open] > summary::before { content: "\\25BE"; }
details.node.hot { outline: 2px solid var(--accent); }
.step-no { font-family: ui-monospace, Consolas, monospace; color: var(--muted); }
.tool { font-weight: 600; }
.ver, .suffix { color: var(--muted); }
.pill { margin-left: auto; font-size: 11px; padding: 1px 8px; border-radius: 999px;
        border: 1px solid var(--line); color: var(--pending); }
.pill.completed { color: var(--ok); border-color: var(--ok); }
.pill.failed { color: var(--bad); border-color: var(--bad); }
.node-body { border-top: 1px solid var(--line); padding: 10px; }
.chips { display: flex; flex-wrap: wrap; gap: 6px; margin: 2px 0 8px; }
.chip { border: 1px solid var(--line); border-radius: 999px; padding: 1px 9px;
        text-decoration: none; color: var(--accent); background: #f4f7fb; font-size: 12px; }
.chip.order { color: var(--muted); background: #f2f2f4; border-style: dashed; }
.notes { border-left: 3px solid #e0b34a; background: #fdf6e6; padding: 6px 10px;
         border-radius: 4px; margin: 8px 0; color: #6b5720; }
.notes ul { margin: 4px 0 0; padding-left: 18px; }
.render { overflow-x: auto; }
.internal { opacity: .72; }
footer.page { color: var(--muted); padding: 8px 24px 32px; }
.pipeline-panel { position: relative; }
.status-chart { position: absolute; top: 10px; right: 14px; text-align: center; }
.donut-num { font-size: 7px; font-weight: 700; fill: var(--fg); }
.donut-cap { font-size: 3.2px; fill: var(--muted); text-transform: uppercase; letter-spacing: .5px; }
.status-legend { font-size: 11px; color: var(--muted); margin-top: 2px; }
@media print { .controls, #wires { display: none; } details.node { break-inside: avoid; } }
"""

_JS = r"""
(function () {
  var flow = document.getElementById("flow");
  var svg = document.getElementById("wires");
  if (!flow || !svg) return;
  var edges = JSON.parse(document.getElementById("edge-data").textContent);
  var SVGNS = "http://www.w3.org/2000/svg";

  function nodeEl(step) { return document.getElementById("step-" + step); }

  var HEAD = 9;              // arrowhead length in px, in the same units as the paths
  var ANCHOR_GAP = 7;        // vertical spacing between edges meeting at one node
  var ANCHOR_INSET = 16;     // how far below a card's top edge the first anchor sits

  /* Only recovered dataflow is drawn. Execution order needs no arrow: it is the order the steps
     are listed in, and a dashed line for it invited being read as a dependency. */
  function visible(e) { return e.kind !== "order"; }

  function anchorPlan() {
    /* Several edges meeting at one node used to terminate on the identical point, which is what
       piled the arrowheads on top of each other. Each end gets its own slot on its own node. */
    var slots = {};
    function reserve(step, key) {
      var list = slots[step] || (slots[step] = []);
      if (list.indexOf(key) === -1) list.push(key);
      return list;
    }
    edges.forEach(function (e, i) {
      if (!visible(e)) return;
      reserve(e.from, "o" + i);
      reserve(e.to, "i" + i);
    });
    return slots;
  }

  function anchorY(step, key, rect, top, slots) {
    var list = slots[step] || [];
    var idx = Math.max(0, list.indexOf(key));
    var span = Math.max(1, list.length);
    /* Anchors stay in the card's header band so a tall expanded card does not drag them down. */
    var band = Math.min(rect.height - 6, ANCHOR_INSET + (span - 1) * ANCHOR_GAP);
    var offset = span === 1 ? Math.min(rect.height / 2, ANCHOR_INSET)
                            : ANCHOR_INSET * 0.6 + (band - ANCHOR_INSET * 0.6) * (idx / (span - 1));
    return top + rect.top + offset;
  }

  function draw() {
    while (svg.firstChild) svg.removeChild(svg.firstChild);
    var base = flow.getBoundingClientRect();
    /* Draw in CSS pixels: a viewBox stretched by preserveAspectRatio="none" scales x and y
       independently, and a marker inherits that, so the arrowhead skewed and drifted as soon as a
       card expanded and changed the height. */
    svg.removeAttribute("viewBox");
    svg.removeAttribute("preserveAspectRatio");
    svg.setAttribute("width", String(base.width));
    svg.setAttribute("height", String(base.height));

    var slots = anchorPlan();
    var corridor = {};

    edges.forEach(function (e, i) {
      if (!visible(e)) return;
      var a = nodeEl(e.from), b = nodeEl(e.to);
      if (!a || !b) return;
      var ra = a.getBoundingClientRect(), rb = b.getBoundingClientRect();
      var y1 = anchorY(e.from, "o" + i, ra, -base.top, slots);
      var y2 = anchorY(e.to, "i" + i, rb, -base.top, slots);
      var x = ra.left - base.left - 8;

      /* Two edges spanning the same distance would otherwise bow identically and overlap exactly. */
      var span = Math.abs(e.to - e.from);
      corridor[span] = (corridor[span] || 0) + 1;
      var bow = Math.min(70, 16 + 9 * span + 5 * (corridor[span] - 1));

      /* The curve stops where the head begins, and the head is a plain triangle drawn at that
         point. The final segment is horizontal by construction -- its control point shares y2 --
         so the head needs no rotation, and nothing depends on how a marker resolves refX. */
      var tipX = x, baseX = x - HEAD;
      var p = document.createElementNS(SVGNS, "path");
      p.setAttribute("d", "M" + x + "," + y1 +
                     " C" + (x - bow) + "," + y1 +
                     " " + (x - bow) + "," + y2 +
                     " " + baseX + "," + y2);
      p.setAttribute("fill", "none");
      p.setAttribute("stroke", "var(--wire)");
      p.setAttribute("stroke-width", "1.8");
      p.setAttribute("stroke-linecap", "round");
      p.setAttribute("data-from", e.from);
      p.setAttribute("data-to", e.to);
      svg.appendChild(p);

      var head = document.createElementNS(SVGNS, "path");
      head.setAttribute("d", "M" + baseX + "," + (y2 - HEAD * 0.45) +
                        " L" + tipX + "," + y2 +
                        " L" + baseX + "," + (y2 + HEAD * 0.45) + " z");
      head.setAttribute("fill", "var(--wire)");
      head.setAttribute("data-from", e.from);
      head.setAttribute("data-to", e.to);
      head.setAttribute("data-head", "1");
      svg.appendChild(head);
    });
  }

  var pending = null;
  function schedule() {
    if (pending) cancelAnimationFrame(pending);
    pending = requestAnimationFrame(function () { pending = null; draw(); });
  }

  flow.addEventListener("toggle", schedule, true);
  window.addEventListener("resize", schedule);
  window.addEventListener("load", schedule);

  document.getElementById("expand-all").addEventListener("click", function () {
    flow.querySelectorAll("details.node").forEach(function (d) { d.open = true; });
    schedule();
  });
  document.getElementById("collapse-all").addEventListener("click", function () {
    flow.querySelectorAll("details.node").forEach(function (d) { d.open = false; });
    schedule();
  });

  var filter = document.getElementById("filter");
  filter.addEventListener("input", function () {
    var q = filter.value.trim().toLowerCase();
    flow.querySelectorAll("details.node").forEach(function (d) {
      var hit = !q || (d.dataset.search || "").indexOf(q) !== -1;
      d.style.display = hit ? "" : "none";
    });
    schedule();
  });

  flow.addEventListener("mouseover", function (ev) {
    var node = ev.target.closest && ev.target.closest("details.node");
    if (!node) return;
    var step = node.id.replace("step-", "");
    svg.querySelectorAll("path").forEach(function (p) {
      var on = p.getAttribute("data-from") === step || p.getAttribute("data-to") === step;
      if (!p.getAttribute("data-head")) p.setAttribute("stroke-width", on ? "3" : "1.8");
      p.setAttribute("opacity", on ? "1" : "0.45");
    });
  });
  flow.addEventListener("mouseleave", function () {
    svg.querySelectorAll("path").forEach(function (p) { p.setAttribute("opacity", "1"); });
  });

  document.querySelectorAll("a.chip[data-jump]").forEach(function (a) {
    a.addEventListener("click", function (ev) {
      ev.preventDefault();
      var t = nodeEl(a.dataset.jump);
      if (!t) return;
      t.open = true;
      t.classList.add("hot");
      t.scrollIntoView({ behavior: "smooth", block: "center" });
      setTimeout(function () { t.classList.remove("hot"); }, 1600);
      schedule();
    });
  });

  draw();
})();
"""


def _node_metadata_html(step, root, exported):
    rows = [
        ("tool version", step.get("tool_version")),
        ("output folder", _rel(step.get("output_folder"), root)),
        ("step script", _rel(step.get("script_file"), root)),
        ("log", _rel(step.get("log_file"), root)),
        ("resources", _fmt_resources(step.get("resources"))),
        ("environment", ", ".join(step.get("environments") or []) or "(none)"),
        ("completion marker", _rel(step.get("completed_marker"), root)),
    ]
    if step.get("job_name"):
        rows.insert(0, ("job name", step["job_name"]))
    if exported:
        params = ((exported.get("configuration") or {}).get("tool_parameters") or {})
        tool_params = params.get("tool_params")
        if tool_params:
            rows.append(("tool parameters", json.dumps(tool_params, sort_keys=True)[:600]))
        rows.append(("inventory", _rel(step.get("tool_outputs_json"), root)))
    parts = ['<dl class="kv">']
    for key, value in rows:
        parts.append(f"<dt>{_esc(key)}</dt><dd>{_esc(value)}</dd>")
    parts.append("</dl>")
    return "".join(parts)


def _wiring_chips(step_no, in_edges, out_edges, order_in, order_out):
    chips = []
    for e in in_edges:
        label = f"&#8592; {e['from_step']:03d} {_esc(e['from_tool'])}"
        if e.get("stream"):
            label += f" · {_esc(e['stream'])}"
        args = ", ".join(e.get("arguments") or [])
        chips.append(
            f'<a class="chip" href="#step-{e["from_step"]}" data-jump="{e["from_step"]}" '
            f'title="recovered dataflow, stored as {_esc(e.get("shape"))} under {_esc(args)}">'
            f'{label}</a>'
        )
    for e in out_edges:
        label = f"&#8594; {e['to_step']:03d} {_esc(e['to_tool'])}"
        if e.get("stream"):
            label += f" · {_esc(e['stream'])}"
        chips.append(
            f'<a class="chip" href="#step-{e["to_step"]}" data-jump="{e["to_step"]}" '
            f'title="recovered dataflow">{label}</a>'
        )
    if not in_edges and order_in is not None:
        chips.append(
            f'<a class="chip order" href="#step-{order_in}" data-jump="{order_in}" '
            f'title="no dataflow recovered — this is the preceding step, not a dependency">'
            f'&#8592; {order_in:03d} (order only)</a>'
        )
    if not out_edges and order_out is not None:
        chips.append(
            f'<a class="chip order" href="#step-{order_out}" data-jump="{order_out}" '
            f'title="no dataflow recovered — this is the following step, not a dependency">'
            f'&#8594; {order_out:03d} (order only)</a>'
        )
    if not chips:
        return ""
    return '<div class="chips">' + "".join(chips) + "</div>"


def _content_summary(output):
    """One scannable line naming what the step declares, for the collapsed card."""
    if output is None:
        return ""
    from biopipelines.datastream import DataStream

    bits = []
    for name, stream in output.streams.items():
        if isinstance(stream, DataStream) and len(stream) > 0:
            bits.append(f"{name}×{len(stream)}")
    tables = list((getattr(output.tables, "_tables", {}) or {}).keys())
    if tables:
        bits.append("tables: " + ", ".join(tables[:4]) + ("…" if len(tables) > 4 else ""))
    line = " · ".join(bits[:6])
    return line if len(line) <= 90 else line[:87] + "…"


def _notes_html(notes):
    if not notes:
        return ""
    items = "".join(f"<li>{_esc(n)}</li>" for n in notes)
    return f'<div class="notes"><strong>Not rendered</strong><ul>{items}</ul></div>'


try:
    from biopipelines import __version__ as _framework_version
except Exception:
    _framework_version = None


def render_page(graph, out_path, allow_external=False):
    """Write the pipeline page for ``graph`` to ``out_path`` and return the path.

    A renderer that needs a browser library (the py3Dmol structure and grid viewers) is served from the vendored copy under renderers/vendor/, inlined once for the whole page, so those nodes render from a file:// URL with no network. ``allow_external=True`` matters only when that copy is missing: it then keeps the renderer's remote <script> tag instead of falling back to the self-contained metadata table, and the page says which renderer was dropped either way.
    """
    from biopipelines.base_config import StandardizedOutput

    meta = graph.get("pipeline", {})
    root = meta.get("output_folder", "")
    steps = graph.get("steps", [])
    batches = graph.get("batches", [])
    recovered = [e for e in graph.get("edges", []) if e.get("from_step") and e.get("to_step")]

    by_step = {s["execution_order"]: s for s in steps}
    ordered = sorted(steps, key=lambda s: s["execution_order"])
    in_edges = {s["execution_order"]: [] for s in steps}
    out_edges = {s["execution_order"]: [] for s in steps}
    for e in recovered:
        if e["to_step"] in in_edges:
            in_edges[e["to_step"]].append(e)
        if e["from_step"] in out_edges:
            out_edges[e["from_step"]].append(e)

    statuses = {s["execution_order"]: _step_status(s) for s in steps}
    any_ran = any(v != "pending" for v in statuses.values())

    # Consecutive steps only, and only without a recovered edge, so a sequence never overdraws a real dependency.
    wire_data = [
        {"from": e["from_step"], "to": e["to_step"], "kind": "flow"} for e in recovered
    ]
    recovered_pairs = {(e["from_step"], e["to_step"]) for e in recovered}
    for prev, cur in zip(ordered, ordered[1:]):
        pair = (prev["execution_order"], cur["execution_order"])
        if pair not in recovered_pairs:
            wire_data.append({"from": pair[0], "to": pair[1], "kind": "order"})

    page_notes = []
    renderers_cfg = _resolve_renderers_config(meta.get("config_variant"), page_notes)

    total_bytes = 0
    node_html = {}
    page_libs = set()
    for step in ordered:
        notes = list(page_notes)
        output, exported, problem = _rebuild_output(step)
        if problem:
            notes.append(problem)
        body = ""
        step_libs = set()
        if output is not None:
            try:
                body = _render_stream_bodies(output, renderers_cfg, notes, allow_external,
                                             step_libs)
            except Exception as e:
                notes.append(f"rendering this step's outputs failed: {e}")
        if len(body) > _MAX_NODE_BODY_BYTES:
            notes.append(
                f"rendered output truncated at {_MAX_NODE_BODY_BYTES // 1000} kB "
                f"(was {len(body) // 1000} kB) to keep the page openable"
            )
            body = ""
        if total_bytes + len(body) > _MAX_TOTAL_BODY_BYTES:
            notes.append(
                f"rendered output omitted: the page had already reached its "
                f"{_MAX_TOTAL_BODY_BYTES // 1000000} MB body budget"
            )
            body = ""
        if body:
            page_libs |= step_libs
        total_bytes += len(body)
        node_html[step["execution_order"]] = (
            _node_metadata_html(step, root, exported),
            _notes_html(notes),
            body,
            _content_summary(output),
        )

    order_index = {s["execution_order"]: i for i, s in enumerate(ordered)}

    lanes = []
    steps_by_batch = {}
    for step in ordered:
        steps_by_batch.setdefault(step.get("batch", -1), []).append(step)

    for batch in batches + ([{"index": -1, "resources": {}, "parents": [], "after_parents": [],
                              "packed": False}] if -1 in steps_by_batch else []):
        idx = batch["index"]
        members = steps_by_batch.get(idx, [])
        if not members:
            continue
        parents = batch.get("parents") or []
        after = batch.get("after_parents") or []
        dep_bits = []
        if parents:
            dep_bits.append("after success of batch " + ", ".join(str(p) for p in parents))
        if after:
            dep_bits.append("after start of batch " + ", ".join(str(p) for p in after))
        if batch.get("packed"):
            dep_bits.append("packed job steps")
        head = (
            f'<div class="lane-head">Batch {idx if idx >= 0 else "?"}'
            f'<span class="lane-meta"> &middot; {_esc(_fmt_resources(batch.get("resources")))}'
            + (f" &middot; {_esc('; '.join(dep_bits))}" if dep_bits else "")
            + f' &middot; {len(members)} step{"s" if len(members) != 1 else ""}</span></div>'
        )

        cards = []
        for step in members:
            no = step["execution_order"]
            metadata_html, notes_html, body, content = node_html[no]
            pos = order_index[no]
            order_in = ordered[pos - 1]["execution_order"] if pos > 0 else None
            order_out = ordered[pos + 1]["execution_order"] if pos + 1 < len(ordered) else None
            status = statuses[no]
            label = (f"{step.get('internal_order') or no:03d}"
                     if step.get("internal") else f"{step.get('public_step') or no:03d}")
            version = step.get("tool_version")
            search = " ".join(str(v) for v in (
                step.get("tool"), step.get("job_name"), step.get("suffix"), label
            ) if v).lower()
            cards.append(
                f'<details class="node{" internal" if step.get("internal") else ""}" '
                f'id="step-{no}" data-search="{_esc(search)}">'
                f'<summary>'
                f'<span class="step-no">{_esc(label)}</span>'
                f'<span class="tool">{_esc(step.get("tool"))}</span>'
                + (f'<span class="ver">v{_esc(version)}</span>' if version else "")
                + (f'<span class="suffix">{_esc(step.get("suffix"))}</span>'
                   if step.get("suffix") else "")
                + ('<span class="suffix">internal</span>' if step.get("internal") else "")
                + (f'<span class="ver">{_esc(content)}</span>' if content else "")
                + f'<span class="pill {status}">{_esc(_STATUS_LABELS[status])}</span>'
                f'</summary>'
                f'<div class="node-body">'
                + _wiring_chips(no, in_edges.get(no, []), out_edges.get(no, []),
                                order_in, order_out)
                + metadata_html + notes_html
                + (f'<div class="render">{body}</div>' if body else "")
                + '</div></details>'
            )
        lanes.append(f'<section class="lane">{head}<div class="nodes">'
                     + "".join(cards) + "</div></section>")

    n_no_input = sum(1 for s in ordered if not in_edges.get(s["execution_order"]))
    banner = (
        '<div class="banner">This page was written at save time: <strong>no step has run '
        'yet</strong>, so it shows structure, wiring and resources, and every node reports '
        'its outputs as pending. Regenerate it after the run '
        '(<code>regenerate_pipeline_page(&lt;RunTime folder&gt;)</code>) to fill in the '
        'rendered outputs.</div>'
        if not any_ran else ""
    )
    failed = [s for s in ordered if statuses[s["execution_order"]] == "failed"]
    if failed:
        banner += (
            '<div class="banner bad">'
            + _esc(f"{len(failed)} step(s) left a FAILED marker: "
                   + ", ".join(f"{s['execution_order']:03d} {s['tool']}" for s in failed))
            + "</div>"
        )

    # Counted here rather than in the graph so an older graph still gains the summary on regeneration.
    public = [st for st in steps if not st.get("internal")]
    environments = []
    for st in steps:
        for env in (st.get("environments") or []):
            if env not in environments:
                environments.append(env)
    pipeline_kv = "".join(
        f"<dt>{_esc(k)}</dt><dd>{_esc(v)}</dd>" for k, v in [
            ("project", meta.get("project")),
            ("job", meta.get("job")),
            ("description", meta.get("description")),
            ("steps", f"{len(public)} in {len(batches)} batch{'es' if len(batches) != 1 else ''}"),
            ("dataflow edges", str(len(recovered))),
            ("environments", ", ".join(environments)),
            ("failed", ", ".join(f"{st['execution_order']:03d} {st['tool']}" for st in failed) or None),
            ("config variant", meta.get("config_variant")),
            ("on the fly", meta.get("on_the_fly")),
            ("output folder", meta.get("output_folder")),
            ("biopipelines", meta.get("biopipelines_version") or _framework_version),
            ("python", meta.get("python_version")),
            ("host", meta.get("host")),
            ("graph saved", meta.get("saved_at")),
        ] if v not in (None, "")
    )

    tallies = [("completed", sum(1 for st in public if statuses[st["execution_order"]] == "completed"), "var(--ok)"),
               ("failed", sum(1 for st in public if statuses[st["execution_order"]] == "failed"), "var(--bad)"),
               ("pending", sum(1 for st in public
                               if statuses[st["execution_order"]] not in ("completed", "failed")), "var(--line)")]
    total = sum(count for _label, count, _colour in tallies) or 1
    # A donut out of dash offsets rather than arc maths: circumference 100 makes each share its own length.
    segments, offset = [], 0.0
    for label, count, colour in tallies:
        if not count:
            continue
        share = 100.0 * count / total
        segments.append(
            f'<circle class="slice" r="15.9155" cx="21" cy="21" fill="none" stroke="{colour}"'
            f' stroke-width="8" stroke-dasharray="{share:.3f} {100 - share:.3f}"'
            f' stroke-dashoffset="{-offset:.3f}" transform="rotate(-90 21 21)">'
            f'<title>{_esc(f"{count} {label}")}</title></circle>')
        offset += share
    legend = " &middot; ".join(f'{count} {label}' for label, count, _c in tallies if count)
    status_chart = (
        f'<div class="status-chart"><svg viewBox="0 0 42 42" width="86" height="86" role="img"'
        f' aria-label="{_esc(legend)}">{"".join(segments)}'
        f'<text x="21" y="22.5" text-anchor="middle" class="donut-num">{total}</text>'
        f'<text x="21" y="27.5" text-anchor="middle" class="donut-cap">steps</text></svg>'
        f'<div class="status-legend">{legend}</div></div>'
    )

    machine = graph.get("machine") or {}
    machine_panel = ""
    if machine:
        rows = "".join(
            f"<dt>{_esc(k.replace('_', ' '))}</dt><dd>{_esc(v if not isinstance(v, list) else ', '.join(map(str, v)))}</dd>"
            for k, v in machine.items())
        machine_panel = f'<div class="panel"><h2>Machine</h2><dl class="kv">{rows}</dl></div>'

    from datetime import datetime
    generated = datetime.now().isoformat(timespec="seconds")

    library_scripts = _library_scripts(page_libs)
    embedded = ", ".join(_esc(_VENDORED_LIBS[n][1]) for n in sorted(page_libs))
    # Measured from the bodies actually kept, not from allow_external, so the claim cannot drift from the page.
    if any(_EXTERNAL_REF_RE.search(b) for _m, _n, b, _c in node_html.values()):
        provenance = "CDN-dependent renderers kept: this page needs a network."
    elif page_libs:
        provenance = (f"self-contained, no network required &mdash; {embedded} "
                      f"embedded from biopipelines/renderers/vendor/.")
    else:
        provenance = "self-contained, no network required."

    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{_esc(meta.get("project"))} / {_esc(meta.get("job"))} &mdash; pipeline</title>
<style>{_CSS}</style>
{StandardizedOutput._CSS}
{library_scripts}
</head>
<body>
<header class="page">
  <h1>{_esc(meta.get("project"))} / {_esc(meta.get("job"))}</h1>
  <div class="sub">{len(ordered)} step{"s" if len(ordered) != 1 else ""} in
    {len([b for b in batches if steps_by_batch.get(b["index"])])} batch(es) &middot;
    {len(recovered)} recovered dataflow edge{"s" if len(recovered) != 1 else ""} &middot;
    {n_no_input} step(s) with no recovered input</div>
</header>
<div class="wrap">
  {banner}
  <div class="panel pipeline-panel"><h2>Pipeline</h2>{status_chart}<dl class="kv">{pipeline_kv}</dl></div>
  {machine_panel}

  <div class="panel controls">
    <button id="expand-all" type="button">Expand all</button>
    <button id="collapse-all" type="button">Collapse all</button>
    <input id="filter" type="search" placeholder="filter by tool or step&hellip;">
  </div>

  <div id="flow">
    <svg id="wires" aria-hidden="true"></svg>
    {"".join(lanes)}
  </div>

</div>
<footer class="page">An arrow is drawn only where the consuming tool holds an object stamped with its
  producer; a stream derived inside <code>DataStream</code> itself (<code>stream.chunks(...)</code>,
  iterating a bare stream) is a fresh object with no back-reference, so its edge is missing rather
  than guessed.<br>Generated {_esc(generated)} by renderers/pipeline_report.py &mdash;
  {provenance}</footer>
<script id="edge-data" type="application/json">{json.dumps(wire_data)}</script>
<script>{_JS}</script>
</body>
</html>
"""

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(document)
    return out_path
