# Renderer for dx volumetric streams (electrostatic potential, density maps).
# When the parent output also has a `structures` stream with the same ids,
# the grid is overlaid on the structure (positive isovalue blue, negative red).
# Otherwise the grid is rendered alone as a volumetric box.

import os
import json
import random


def _iter_id_file(stream):
    """Yield (id, file_path) pairs for a per-id stream.

    Prefer the map_table: it holds the fully-expanded ids and concrete paths the
    run actually produced. ``files_expanded`` only consults the map_table in
    runtime mode, so a post-run display of a lazy-id stream would otherwise
    expand to just the deterministic prefix and miss every real file. Fall back
    to the id/file zip for streams without a map on disk."""
    map_data = stream._get_map_data()
    if map_data is not None and len(map_data) > 0:
        file_col = next((c for c in ("file", "file_path") if c in map_data.columns), None)
        id_col = "id" if "id" in map_data.columns else None
        if file_col and id_col:
            for _, row in map_data.iterrows():
                yield str(row[id_col]), str(row[file_col])
            return
    for item_id, file_path in zip(stream.ids_expanded, stream.files_expanded):
        yield item_id, file_path


def render(stream, output):
    if stream.format not in ("dx",):
        return ""

    max_grids = 10
    grids = []
    for gid, file_path in _iter_id_file(stream):
        if file_path and os.path.isfile(file_path):
            try:
                with open(file_path, "r") as f:
                    grids.append((gid, f.read(), file_path))
            except Exception:
                pass
        if len(grids) >= max_grids:
            break
    if not grids:
        return ""

    # Pair each grid id with a structure from the parent StandardizedOutput.
    structures_by_id = {}
    structures_stream = None
    try:
        sibling = output.streams.structures if hasattr(output.streams, "structures") else None
        if sibling is not None and len(sibling) > 0 and sibling.format in ("pdb", "pqr"):
            structures_stream = sibling
            for sid, sfile in _iter_id_file(sibling):
                if sfile and os.path.isfile(sfile):
                    try:
                        with open(sfile, "r") as f:
                            structures_by_id[sid] = (f.read(), sfile)
                    except Exception:
                        pass
    except Exception:
        pass

    viewer_id = f"bpdx_{random.randint(100000, 999999)}"
    grid_ids = [g[0] for g in grids]
    grid_data = [g[1] for g in grids]
    struct_data = [structures_by_id.get(g[0], ("", ""))[0] for g in grids]
    has_structure = any(struct_data)

    truncated = len(stream) > max_grids
    total_label = f"{len(grids)} grid{'s' if len(grids) != 1 else ''}"
    if truncated:
        total_label += f" (of {len(stream)} total)"

    grid_ids_json = json.dumps(grid_ids)
    grid_data_json = json.dumps(grid_data)
    struct_data_json = json.dumps(struct_data)
    has_struct_json = json.dumps(has_structure)

    return f"""
<script src="https://cdn.jsdelivr.net/npm/3dmol@2.5.2/build/3Dmol-min.js"></script>
<div style="margin-top: 12px;">
  <strong>Electrostatic potential / volumetric grid</strong> ({total_label})
</div>
<div id="{viewer_id}_container" style="position: relative; width: 800px;">
  <div id="{viewer_id}_viewer" style="width: 800px; height: 500px; position: relative;"></div>
  <div id="{viewer_id}_toolbar" style="display: flex; align-items: center; justify-content: center; gap: 6px; margin-top: 6px; font-family: monospace; flex-wrap: wrap;">
    <span style="color: #888; font-size: 0.82em;">Isovalue (kT/e):</span>
    <input id="{viewer_id}_iso" type="range" min="0.5" max="5" step="0.25" value="1.0" oninput="bpdx_{viewer_id}_setIso(this.value)" style="width: 140px;">
    <span id="{viewer_id}_iso_label" style="font-size: 0.82em; min-width: 32px;">1.0</span>
    <span style="color: #ccc;">|</span>
    <button id="{viewer_id}_btn_spin" onclick="bpdx_{viewer_id}_toggleSpin()" style="padding: 3px 10px; font-size: 0.82em; cursor: pointer; border: 1px solid #ccc; border-radius: 4px; background: #f5f5f5;">Spin</button>
  </div>
  <div style="display: flex; align-items: center; justify-content: center; gap: 12px; margin-top: 6px; font-family: monospace;">
    <button onclick="bpdx_{viewer_id}_navigate(-1)" style="padding: 4px 14px; font-size: 1.1em; cursor: pointer; border: 1px solid #ccc; border-radius: 4px; background: #f5f5f5;">&#9664;</button>
    <span id="{viewer_id}_label" style="min-width: 200px; text-align: center; font-size: 0.95em;"></span>
    <button onclick="bpdx_{viewer_id}_navigate(1)" style="padding: 4px 14px; font-size: 1.1em; cursor: pointer; border: 1px solid #ccc; border-radius: 4px; background: #f5f5f5;">&#9654;</button>
  </div>
</div>
<script>
(function() {{
  var ids = {grid_ids_json};
  var dxData = {grid_data_json};
  var structData = {struct_data_json};
  var hasStruct = {has_struct_json};
  var idx = 0;
  var viewer = null;
  var isoVal = 1.0;
  var spinning = false;

  function showGrid(i) {{
    if (!viewer) return;
    idx = i;
    if (idx < 0) idx = ids.length - 1;
    if (idx >= ids.length) idx = 0;
    viewer.clear();
    if (hasStruct && structData[idx]) {{
      viewer.addModel(structData[idx], "pdb");
      viewer.setStyle({{}}, {{cartoon: {{color: "white"}}}});
    }}
    var voldata = new $3Dmol.VolumeData(dxData[idx], "dx");
    viewer.addIsosurface(voldata, {{
      isoval: isoVal, color: "blue", opacity: 0.5, smoothness: 5
    }});
    viewer.addIsosurface(voldata, {{
      isoval: -isoVal, color: "red", opacity: 0.5, smoothness: 5
    }});
    viewer.zoomTo();
    viewer.render();
    document.getElementById("{viewer_id}_label").textContent =
      ids[idx] + "  (" + (idx+1) + "/" + ids.length + ")";
  }}

  function initViewer() {{
    if (typeof $3Dmol === "undefined") {{
      setTimeout(initViewer, 200);
      return;
    }}
    var el = document.getElementById("{viewer_id}_viewer");
    viewer = $3Dmol.createViewer(el, {{backgroundColor: "white"}});
    showGrid(0);
  }}

  window.bpdx_{viewer_id}_setIso = function(v) {{
    isoVal = parseFloat(v);
    document.getElementById("{viewer_id}_iso_label").textContent = isoVal.toFixed(2);
    showGrid(idx);
  }};
  window.bpdx_{viewer_id}_navigate = function(delta) {{ showGrid(idx + delta); }};
  window.bpdx_{viewer_id}_toggleSpin = function() {{
    spinning = !spinning;
    var btn = document.getElementById("{viewer_id}_btn_spin");
    if (spinning) {{
      viewer.spin("y", 1);
      btn.style.background = "#ddeeff";
      btn.style.fontWeight = "bold";
    }} else {{
      viewer.spin(false);
      btn.style.background = "#f5f5f5";
      btn.style.fontWeight = "normal";
    }}
  }};

  initViewer();
}})();
</script>
"""
