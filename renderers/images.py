# Renderer for image streams: inline display with base64 encoding

import os
import base64
import html as html_module


def render(stream, output):
    """Render image stream items inline."""
    parts = []
    parts.append(
        f'<div class="bp-section">'
        f'<div class="bp-section-title">{html_module.escape(stream.name)} '
        f'<span style="font-weight: normal; color: #666;">({stream.format}, {len(stream)} items)</span></div>'
    )

    for item_id, file_path in zip(stream.ids, stream.files):
        if not file_path or not os.path.isfile(file_path):
            continue

        if stream.format.lower() == "svg":
            try:
                with open(file_path, "r") as f:
                    svg_content = f.read()
                parts.append(
                    f'<div style="margin: 4px 0;">'
                    f'<div style="color: #666; font-size: 0.85em;">{html_module.escape(str(item_id))}</div>'
                    f"{svg_content}</div>"
                )
            except Exception:
                pass
        else:
            try:
                with open(file_path, "rb") as f:
                    img_data = base64.b64encode(f.read()).decode("utf-8")
                mime = f"image/{stream.format.lower()}"
                if stream.format.lower() in ("jpg", "jpeg"):
                    mime = "image/jpeg"
                parts.append(
                    f'<div style="margin: 4px 0;">'
                    f'<div style="color: #666; font-size: 0.85em;">{html_module.escape(str(item_id))}</div>'
                    f'<img src="data:{mime};base64,{img_data}" '
                    f'style="max-width: 100%; height: auto;" /></div>'
                )
            except Exception:
                pass

    parts.append("</div>")
    return "\n".join(parts)
