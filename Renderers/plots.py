# Renderer for plot streams: side-by-side plot image + companion CSV table

import os
import base64
import html as html_module
import pandas as pd


def _render_compact_table(df):
    """Render a DataFrame as a compact HTML table with 2+...+2 row truncation."""
    columns = list(df.columns)
    rows = ['<table class="bp-table"><tr>']
    for col in columns:
        rows.append(f'<th>{html_module.escape(str(col))}</th>')
    rows.append('</tr>')

    n_rows = len(df)
    if n_rows <= 4:
        display_rows = list(range(n_rows))
        ellipsis_after = None
    else:
        display_rows = list(range(2)) + list(range(n_rows - 2, n_rows))
        ellipsis_after = 2

    for row_i, idx in enumerate(display_rows):
        if ellipsis_after is not None and row_i == ellipsis_after:
            rows.append(
                f'<tr class="bp-ellipsis"><td colspan="{len(columns)}">'
                f'... {n_rows - 4} more ...</td></tr>'
            )
        row = df.iloc[idx]
        rows.append('<tr>')
        for col in columns:
            val = str(row[col]) if pd.notna(row[col]) else ''
            display_val = val if len(val) <= 60 else val[:57] + '...'
            rows.append(f'<td>{html_module.escape(display_val)}</td>')
        rows.append('</tr>')

    rows.append('</table>')
    return "\n".join(rows)


def render(stream, output):
    """Render plots as side-by-side image + companion CSV table."""
    parts = []
    parts.append(
        '<div class="bp-section">'
        '<div class="bp-section-title">plots '
        f'<span style="font-weight: normal; color: #666;">'
        f'({stream.format}, {len(stream)} items)</span></div>'
    )

    for item_id, file_path in zip(stream.ids, stream.files):
        if not file_path or not os.path.isfile(file_path):
            continue

        # Read PNG as base64
        try:
            with open(file_path, "rb") as f:
                img_data = base64.b64encode(f.read()).decode("utf-8")
        except Exception:
            continue

        # Find companion CSV (same name, .csv extension)
        csv_path = os.path.splitext(file_path)[0] + ".csv"
        table_html = ""
        if os.path.isfile(csv_path):
            try:
                df = pd.read_csv(csv_path)
                if len(df) > 0:
                    table_html = _render_compact_table(df)
            except Exception:
                pass

        # Title from filename
        title = os.path.splitext(os.path.basename(file_path))[0].replace("_", " ")

        # Side-by-side: table left, image right
        parts.append(
            f'<div style="margin: 8px 0;">'
            f'<div style="color: #666; font-size: 0.85em; margin-bottom: 4px;">{html_module.escape(title)}</div>'
            f'<div style="display: flex; align-items: flex-start; gap: 16px; flex-wrap: wrap;">'
        )

        if table_html:
            parts.append(
                f'<div style="flex: 0 1 auto; max-width: 50%; overflow-x: auto;">'
                f'{table_html}</div>'
            )

        parts.append(
            f'<div style="flex: 1 1 auto;">'
            f'<img src="data:image/png;base64,{img_data}" '
            f'style="max-width: 100%; height: auto;" /></div>'
        )

        parts.append('</div></div>')

    parts.append('</div>')
    return "\n".join(parts)
