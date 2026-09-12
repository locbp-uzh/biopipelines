# Renderer for pdb/cif structure streams: interactive 3D viewer (py3Dmol)

import glob
import hashlib
import os
import json
import random


_LIB_3DMOL_URL = "https://cdn.jsdelivr.net/npm/3dmol@2.5.2/build/3Dmol-min.js"

# The data-bp-lib marker lets a page assembler swap this tag for one inlined copy of renderers/vendor/3Dmol-min.js; standalone (notebook) output keeps the CDN URL.
_LIB_3DMOL_TAG = f'<script data-bp-lib="3dmol" src="{_LIB_3DMOL_URL}"></script>'

def _script_json(value):
    """``json.dumps`` for a value going inside an inline ``<script>``.

    An inline script ends at the first ``</script``, and ``json.dumps`` escapes quotes and backslashes but not that sequence -- so a PDB REMARK containing it would close the viewer's script early and execute whatever followed as markup. ``pipeline_report`` neutralizes the same sequence in the vendored library it inlines; file contents need it for the same reason.
    """
    return json.dumps(value).replace("</", "<\\/")


def _resolve_path(file_path):
    """Return file_path if it's a real file, or expand a wildcard to the first
    match on disk. Tools that don't know the final extension at declaration
    time store paths like ``<id>.*``; this resolves them at render."""
    if not file_path:
        return None
    if os.path.isfile(file_path):
        return file_path
    if any(ch in file_path for ch in "*?["):
        matches = sorted(glob.glob(file_path))
        if matches:
            return matches[0]
    return None


def _iter_id_file(stream):
    """Yield (id, file_path) pairs for a per-id structure stream.

    Prefer the map_table: it holds the fully-expanded ids and concrete paths
    the run actually produced. ``files_expanded`` only consults the map_table
    in runtime mode, so a post-run display of a lazy-id stream (Boltz2's
    ``<id>_<1..K>``, split-chain PDB, …) would otherwise expand to just the
    deterministic prefix and miss every real file. Fall back to the id/file
    zip for streams without a map on disk."""
    map_data = stream._get_map_data()
    if map_data is not None and len(map_data) > 0:
        file_col = next((c for c in ("file", "file_path") if c in map_data.columns), None)
        id_col = "id" if "id" in map_data.columns else None
        if file_col and id_col:
            for _, row in map_data.iterrows():
                yield str(row[id_col]), str(row[file_col])
            return
    for struct_id, file_path in zip(stream.ids_expanded, stream.files_expanded):
        yield struct_id, file_path


SAMPLE_HEAD = 3
SAMPLE_RANDOM = 2
MAX_EMBEDDED = SAMPLE_HEAD + SAMPLE_RANDOM


def _embed_budget(stream, output):
    """How many structures may be embedded: the default, or a caller's explicit cap.

    bp-visualize sets ``rendering_parameters[<stream>]["max_embedded"]`` when the user asked
    for a specific count. An explicit ask is a request, not a hint -- returning a sample of a
    requested top-20 would answer a different question than the one put.
    """
    params = getattr(output, "rendering_parameters", None) or {}
    requested = (params.get(stream.name) or {}).get("max_embedded")
    if isinstance(requested, int) and requested > 0:
        return requested
    return MAX_EMBEDDED


def sample_positions(total, seed="", budget=None):
    """Which indices of a structure stream to embed: the first few, plus a few from the rest.

    The page is meant to be copied off the cluster and opened locally, where a file:// link to a
    compute node's filesystem resolves to nothing -- so a structure is only inspectable if its
    contents are inline. That caps how many can go in: whole PDB files, embedded. The first few show
    what the run produced; sampling the rest is what shows whether quality holds across the run,
    which the head alone cannot.

    Deterministic in ``seed`` so regenerating a page from the same graph reproduces it exactly.
    """
    budget = MAX_EMBEDDED if budget is None else budget
    if total <= budget:
        return list(range(total))
    # An explicit budget is an ordered top-N: take the head and keep the caller's order.
    if budget != MAX_EMBEDDED:
        return list(range(budget))
    tail = random.Random(seed).sample(range(SAMPLE_HEAD, total), SAMPLE_RANDOM)
    return list(range(SAMPLE_HEAD)) + sorted(tail)


def render(stream, output):
    """Render an interactive 3D structure viewer for pdb/cif/pqr streams."""
    if not stream.has_only_formats("pdb", "cif", "pqr"):
        return ""

    pdb_data = []
    if stream.is_shared_file:
        # One shared structure file — render once. The stream's ids label
        # logical entries inside the file, not separate artifacts.
        path = _resolve_path(stream.files)
        if path:
            try:
                with open(path, "r") as f:
                    pdb_data.append((stream.name, f.read(), path))
            except Exception:
                pass
    else:
        # List first, then read only the sampled files: reading all of them to throw most away is
        # what made this slow on a large campaign.
        pairs = list(_iter_id_file(stream))
        sampled = sample_positions(len(pairs), seed=f"{stream.name}:{len(pairs)}",
                                   budget=_embed_budget(stream, output))
        for idx in sampled:
            struct_id, file_path = pairs[idx]
            resolved = _resolve_path(file_path)
            if resolved:
                try:
                    with open(resolved, "r") as f:
                        pdb_data.append((struct_id, f.read(), resolved))
                except Exception:
                    pass

    if not pdb_data:
        return ""

    # Detect format per-file from extension (pqr is parsed as pdb by 3Dmol).
    def _detect_fmt(file_path):
        if file_path.endswith(".cif"):
            return "cif"
        return "pdb"

    fmt = _detect_fmt(pdb_data[0][2])

    # Derived, not random: the page is regenerated after every step, and a fresh id each time made
    # two renders of an unchanged run differ.
    viewer_id = "bp3d_" + hashlib.sha1(
        f"{stream.name}:{[sid for sid, _, _fp in pdb_data]}".encode("utf-8")
    ).hexdigest()[:10]

    struct_ids_json = _script_json([sid for sid, _, _fp in pdb_data])
    struct_data_json = _script_json([content for _, content, _fp in pdb_data])

    truncated = len(stream) > len(pdb_data)
    total_label = f"{len(pdb_data)} structure{'s' if len(pdb_data) != 1 else ''}"
    if truncated:
        total_label += (
            f" (of {len(stream)} total — first {SAMPLE_HEAD}"
            f" plus {SAMPLE_RANDOM} sampled)"
        )

    colors = [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
        "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    ]
    colors_json = _script_json(colors)

    # Check for pLDDT coloring from rendering_parameters
    plddt_upper = None
    rendering_params = getattr(output, "rendering_parameters", None)
    if rendering_params:
        stream_params = rendering_params.get(stream.name, {})
        if stream_params.get("color_by") == "plddt":
            plddt_upper = stream_params.get("plddt_upper", 100)

    plddt_upper_json = _script_json(plddt_upper)

    # pLDDT color legend (shown only when pLDDT coloring is active)
    plddt_legend = ""
    if plddt_upper is not None:
        plddt_legend = f"""
  <div id="{viewer_id}_legend" style="display: flex; align-items: center; justify-content: center; gap: 8px; margin-top: 4px; font-family: monospace; font-size: 0.85em;">
    <span style="color: #888;">pLDDT:</span>
    <span style="display:inline-block;width:12px;height:12px;background:#126DFF;border-radius:2px;vertical-align:middle;"></span><span>High (&ge;90%)</span>
    <span style="display:inline-block;width:12px;height:12px;background:#0ECFF1;border-radius:2px;vertical-align:middle;"></span><span>Good (70-90%)</span>
    <span style="display:inline-block;width:12px;height:12px;background:#F6ED12;border-radius:2px;vertical-align:middle;"></span><span>Low (50-70%)</span>
    <span style="display:inline-block;width:12px;height:12px;background:#EE831D;border-radius:2px;vertical-align:middle;"></span><span>Very low (&lt;50%)</span>
  </div>"""

    btn = (
        'style="padding: 3px 10px; font-size: 0.82em; cursor: pointer; '
        'border: 1px solid #ccc; border-radius: 4px; background: #f5f5f5;"'
    )
    btn_active = (
        'style="padding: 3px 10px; font-size: 0.82em; cursor: pointer; '
        'border: 1px solid #888; border-radius: 4px; background: #ddeeff; font-weight: bold;"'
    )

    return f"""
{_LIB_3DMOL_TAG}
<div style="margin-top: 12px;">
  <strong>3D Structure Viewer</strong> ({total_label})
</div>
<div id="{viewer_id}_container" style="position: relative; width: 800px;">
  <div id="{viewer_id}_viewer" style="width: 800px; height: 500px; position: relative;"></div>
  <!-- Style toolbar -->
  <div id="{viewer_id}_toolbar" style="display: flex; align-items: center; justify-content: center; gap: 6px; margin-top: 6px; font-family: monospace; flex-wrap: wrap;">
    <span style="color: #888; font-size: 0.82em;">Style:</span>
    <button id="{viewer_id}_btn_cartoon" onclick="{viewer_id}_setStyle('cartoon')" {btn_active}>Cartoon</button>
    <button id="{viewer_id}_btn_stick" onclick="{viewer_id}_setStyle('stick')" {btn}>Sticks</button>
    <button id="{viewer_id}_btn_sphere" onclick="{viewer_id}_setStyle('sphere')" {btn}>Spheres</button>
    <button id="{viewer_id}_btn_surface" onclick="{viewer_id}_toggleSurface()" {btn}>Surface</button>
    <span style="color: #ccc;">|</span>
    <button id="{viewer_id}_btn_ligands" onclick="{viewer_id}_toggleLigands()" {btn_active}>Ligands: Element</button>
    <span style="color: #ccc;">|</span>
    <button id="{viewer_id}_btn_spin" onclick="{viewer_id}_toggleSpin()" {btn}>Spin</button>
  </div>
  <!-- Navigation -->
  <div style="display: flex; align-items: center; justify-content: center; gap: 12px; margin-top: 6px; font-family: monospace;">
    <button id="{viewer_id}_prev" onclick="{viewer_id}_navigate(-1)"
            style="padding: 4px 14px; font-size: 1.1em; cursor: pointer; border: 1px solid #ccc; border-radius: 4px; background: #f5f5f5;">&#9664;</button>
    <span id="{viewer_id}_label" style="min-width: 200px; text-align: center; font-size: 0.95em;"></span>
    <button id="{viewer_id}_next" onclick="{viewer_id}_navigate(1)"
            style="padding: 4px 14px; font-size: 1.1em; cursor: pointer; border: 1px solid #ccc; border-radius: 4px; background: #f5f5f5;">&#9654;</button>
    <button id="{viewer_id}_dl" onclick="{viewer_id}_download()"
            style="padding: 4px 10px; font-size: 0.85em; cursor: pointer; border: 1px solid #ccc; border-radius: 4px; background: #f5f5f5;" title="Download current structure">&#11015; .{fmt}</button>
  </div>{plddt_legend}
</div>
<script>
(function() {{
  var ids = {struct_ids_json};
  var data = {struct_data_json};
  var fmt = "{fmt}";
  var colors = {colors_json};
  var plddtUpper = {plddt_upper_json};
  var idx = 0;
  var viewer = null;

  // State
  var currentStyle = "cartoon";
  var showSurface = false;
  // Ligand display: "off", "element" (default CPK colors), "plddt" (pLDDT coloring)
  var ligandMode = "element";
  var spinning = false;

  var btnNormal = "padding: 3px 10px; font-size: 0.82em; cursor: pointer; border: 1px solid #ccc; border-radius: 4px; background: #f5f5f5; font-weight: normal;";
  var btnActive = "padding: 3px 10px; font-size: 0.82em; cursor: pointer; border: 1px solid #888; border-radius: 4px; background: #ddeeff; font-weight: bold;";

  function setBtn(name, active) {{
    var el = document.getElementById("{viewer_id}_btn_" + name);
    if (el) el.setAttribute("style", active ? btnActive : btnNormal);
  }}

  function getColor() {{
    return colors[idx % colors.length];
  }}

  function plddtColorfunc(atom) {{
    var b = atom.b;
    var upper = plddtUpper;
    if (b >= 0.9 * upper) return "#126DFF";
    if (b >= 0.7 * upper) return "#0ECFF1";
    if (b >= 0.5 * upper) return "#F6ED12";
    return "#EE831D";
  }}

  function applyStyles() {{
    if (!viewer) return;

    // Selectors
    var protSel = {{"hetflag": false}};
    var hetSel = {{"hetflag": true}};

    // Clear all styles
    viewer.setStyle({{}}, {{}});

    // Protein representation
    var protStyle = {{}};
    if (currentStyle === "cartoon") {{
      if (plddtUpper !== null) {{
        protStyle = {{"cartoon": {{"colorfunc": plddtColorfunc}}}};
      }} else {{
        protStyle = {{"cartoon": {{"color": getColor()}}}};
      }}
    }} else if (currentStyle === "stick") {{
      if (plddtUpper !== null) {{
        protStyle = {{"stick": {{"colorfunc": plddtColorfunc}}}};
      }} else {{
        protStyle = {{"stick": {{"color": getColor()}}}};
      }}
    }} else if (currentStyle === "sphere") {{
      if (plddtUpper !== null) {{
        protStyle = {{"sphere": {{"colorfunc": plddtColorfunc}}}};
      }} else {{
        protStyle = {{"sphere": {{"color": getColor()}}}};
      }}
    }}
    viewer.setStyle(protSel, protStyle);

    // Ligands (heteroatoms)
    if (ligandMode === "element") {{
      viewer.setStyle(hetSel, {{"stick": {{"colorscheme": "default"}}}});
    }} else if (ligandMode === "plddt" && plddtUpper !== null) {{
      viewer.setStyle(hetSel, {{"stick": {{"colorfunc": plddtColorfunc}}}});
    }} else if (ligandMode !== "off") {{
      // plddt requested but not available — fall back to element colors
      viewer.setStyle(hetSel, {{"stick": {{"colorscheme": "default"}}}});
    }} else {{
      viewer.setStyle(hetSel, {{}});
    }}

    // Surface
    viewer.removeAllSurfaces();
    if (showSurface) {{
      if (plddtUpper !== null) {{
        viewer.addSurface($3Dmol.SurfaceType.VDW, {{"opacity": 0.7, "colorfunc": plddtColorfunc}}, protSel);
      }} else {{
        viewer.addSurface($3Dmol.SurfaceType.VDW, {{"opacity": 0.7, "color": getColor()}}, protSel);
      }}
    }}

    viewer.render();
  }}

  function initViewer() {{
    if (typeof $3Dmol === "undefined") {{
      setTimeout(initViewer, 200);
      return;
    }}
    var el = document.getElementById("{viewer_id}_viewer");
    viewer = $3Dmol.createViewer(el, {{backgroundColor: "white"}});
    showStructure(0);
  }}

  function showStructure(i) {{
    if (!viewer) return;
    idx = i;
    if (idx < 0) idx = ids.length - 1;
    if (idx >= ids.length) idx = 0;
    viewer.removeAllModels();
    viewer.removeAllSurfaces();
    viewer.addModel(data[idx], fmt);
    applyStyles();
    viewer.zoomTo();
    viewer.render();
    // Update label. An id comes from a map_table, so it is not markup: escape before innerHTML.
    var esc = function(s) {{
      return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
                      .replace(/"/g, '&quot;').replace(/'/g, '&#39;');
    }};
    var labelHtml;
    if (plddtUpper !== null) {{
      labelHtml = esc(ids[idx]) + '  <span style="color:#888;">(' + (idx+1) + '/' + ids.length + ')</span>';
    }} else {{
      var color = colors[idx % colors.length];
      labelHtml = '<span style="display:inline-block;width:12px;height:12px;background:' + esc(color) +
        ';border-radius:2px;vertical-align:middle;margin-right:6px;"></span>' +
        esc(ids[idx]) + '  <span style="color:#888;">(' + (idx+1) + '/' + ids.length + ')</span>';
    }}
    document.getElementById("{viewer_id}_label").innerHTML = labelHtml;
  }}

  window.{viewer_id}_setStyle = function(style) {{
    currentStyle = style;
    setBtn("cartoon", style === "cartoon");
    setBtn("stick", style === "stick");
    setBtn("sphere", style === "sphere");
    applyStyles();
    viewer.render();
  }};

  window.{viewer_id}_toggleSurface = function() {{
    showSurface = !showSurface;
    setBtn("surface", showSurface);
    applyStyles();
  }};

  function updateLigandBtn() {{
    var el = document.getElementById("{viewer_id}_btn_ligands");
    if (!el) return;
    if (ligandMode === "off") {{
      el.setAttribute("style", btnNormal);
      el.textContent = "Ligands: Off";
    }} else if (ligandMode === "element") {{
      el.setAttribute("style", btnActive);
      el.textContent = "Ligands: Element";
    }} else {{
      el.setAttribute("style", btnActive);
      el.textContent = "Ligands: pLDDT";
    }}
  }}

  window.{viewer_id}_toggleLigands = function() {{
    if (plddtUpper !== null) {{
      // 3-state cycle: element → plddt → off → element ...
      if (ligandMode === "element") ligandMode = "plddt";
      else if (ligandMode === "plddt") ligandMode = "off";
      else ligandMode = "element";
    }} else {{
      // 2-state toggle: element ↔ off
      ligandMode = (ligandMode === "off") ? "element" : "off";
    }}
    updateLigandBtn();
    applyStyles();
  }};

  window.{viewer_id}_toggleSpin = function() {{
    spinning = !spinning;
    setBtn("spin", spinning);
    if (spinning) {{
      viewer.spin("y", 1);
    }} else {{
      viewer.spin(false);
    }}
  }};

  window.{viewer_id}_navigate = function(delta) {{
    showStructure(idx + delta);
  }};

  window.{viewer_id}_download = function() {{
    var blob = new Blob([data[idx]], {{type: "text/plain"}});
    var a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = ids[idx] + "." + fmt;
    a.click();
    URL.revokeObjectURL(a.href);
  }};

  initViewer();
}})();
</script>
"""
