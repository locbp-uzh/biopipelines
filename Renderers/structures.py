# Renderer for pdb/cif structure streams: interactive 3D viewer (py3Dmol)

import os
import json
import random


def render(stream, output):
    """Render an interactive 3D structure viewer for pdb/cif streams."""
    if stream.format not in ("pdb", "cif", "pdb|cif"):
        return ""

    max_structures = 50
    pdb_data = []
    for struct_id, file_path in zip(stream.ids, stream.files):
        if file_path and os.path.isfile(file_path):
            try:
                with open(file_path, "r") as f:
                    pdb_data.append((struct_id, f.read(), file_path))
            except Exception:
                pass
        if len(pdb_data) >= max_structures:
            break

    if not pdb_data:
        return ""

    # Detect format per-file from extension when format is "pdb|cif"
    def _detect_fmt(file_path, default):
        if file_path.endswith(".cif"):
            return "cif"
        if file_path.endswith(".pdb"):
            return "pdb"
        return default

    if stream.format == "pdb|cif":
        fmt = _detect_fmt(pdb_data[0][2], "pdb")
    else:
        fmt = "pdb" if stream.format == "pdb" else "cif"

    viewer_id = f"bp3d_{random.randint(100000, 999999)}"

    struct_ids_json = json.dumps([sid for sid, _, _fp in pdb_data])
    struct_data_json = json.dumps([content for _, content, _fp in pdb_data])

    truncated = len(stream) > max_structures
    total_label = f"{len(pdb_data)} structure{'s' if len(pdb_data) != 1 else ''}"
    if truncated:
        total_label += f" (of {len(stream)} total)"

    colors = [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
        "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    ]
    colors_json = json.dumps(colors)

    # Check for pLDDT coloring from rendering_parameters
    plddt_upper = None
    rendering_params = getattr(output, "rendering_parameters", None)
    if rendering_params:
        stream_params = rendering_params.get(stream.name, {})
        if stream_params.get("color_by") == "plddt":
            plddt_upper = stream_params.get("plddt_upper", 100)

    plddt_upper_json = json.dumps(plddt_upper)

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

    return f"""
<script src="https://cdn.jsdelivr.net/npm/3dmol@2.5.2/build/3Dmol-min.js"></script>
<div style="margin-top: 12px;">
  <strong>3D Structure Viewer</strong> ({total_label})
</div>
<div id="{viewer_id}_container" style="position: relative; width: 800px;">
  <div id="{viewer_id}_viewer" style="width: 800px; height: 500px; position: relative;"></div>
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
    viewer.addModel(data[idx], fmt);
    if (plddtUpper !== null) {{
      var upper = plddtUpper;
      viewer.setStyle({{}}, {{"cartoon": {{"colorfunc": function(atom) {{
        var b = atom.b;
        if (b >= 0.9 * upper) return "#126DFF";
        if (b >= 0.7 * upper) return "#0ECFF1";
        if (b >= 0.5 * upper) return "#F6ED12";
        return "#EE831D";
      }}}}}});
    }} else {{
      var color = colors[idx % colors.length];
      viewer.setStyle({{}}, {{"cartoon": {{"color": color}}}});
    }}
    viewer.zoomTo();
    viewer.render();
    var labelHtml;
    if (plddtUpper !== null) {{
      labelHtml = ids[idx] + '  <span style="color:#888;">(' + (idx+1) + '/' + ids.length + ')</span>';
    }} else {{
      var color = colors[idx % colors.length];
      labelHtml = '<span style="display:inline-block;width:12px;height:12px;background:' + color +
        ';border-radius:2px;vertical-align:middle;margin-right:6px;"></span>' +
        ids[idx] + '  <span style="color:#888;">(' + (idx+1) + '/' + ids.length + ')</span>';
    }}
    document.getElementById("{viewer_id}_label").innerHTML = labelHtml;
  }}

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
