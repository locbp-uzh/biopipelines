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
  </div>
</div>
<script>
(function() {{
  var ids = {struct_ids_json};
  var data = {struct_data_json};
  var fmt = "{fmt}";
  var colors = {colors_json};
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
    var color = colors[idx % colors.length];
    viewer.setStyle({{}}, {{"cartoon": {{"color": color}}}});
    viewer.zoomTo();
    viewer.render();
    document.getElementById("{viewer_id}_label").innerHTML =
      '<span style="display:inline-block;width:12px;height:12px;background:' + color +
      ';border-radius:2px;vertical-align:middle;margin-right:6px;"></span>' +
      ids[idx] + '  <span style="color:#888;">(' + (idx+1) + '/' + ids.length + ')</span>';
  }}

  window.{viewer_id}_navigate = function(delta) {{
    showStructure(idx + delta);
  }};

  initViewer();
}})();
</script>
"""
