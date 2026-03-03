# Default renderer for tables: metadata + CSV content

import os
import html as html_module
import pandas as pd


def _rel(path, output_folder):
    """Make path relative to output_folder for display."""
    if output_folder and path.startswith(output_folder):
        return "<output_folder>/" + os.path.relpath(path, output_folder)
    return path


def _render_dataframe(df):
    """Render a DataFrame as an HTML table with 2+...+2 row truncation."""
    columns = list(df.columns)
    parts = ['<table class="bp-table"><tr>']
    for col in columns:
        parts.append(f'<th>{html_module.escape(str(col))}</th>')
    parts.append('</tr>')

    n_rows = len(df)
    if n_rows <= 4:
        display_rows = list(range(n_rows))
        ellipsis_after = None
    else:
        display_rows = list(range(2)) + list(range(n_rows - 2, n_rows))
        ellipsis_after = 2

    for row_i, idx in enumerate(display_rows):
        if ellipsis_after is not None and row_i == ellipsis_after:
            parts.append(
                f'<tr class="bp-ellipsis"><td colspan="{len(columns)}">'
                f'... {n_rows - 4} more ...</td></tr>'
            )
        row = df.iloc[idx]
        parts.append('<tr>')
        for col in columns:
            val = str(row[col]) if pd.notna(row[col]) else ''
            display_val = val if len(val) <= 60 else val[:57] + '...'
            parts.append(f'<td>{html_module.escape(display_val)}</td>')
        parts.append('</tr>')
    parts.append('</table>')
    return "\n".join(parts)


def render(table_info, output):
    """Render table metadata and CSV content."""
    rel = lambda p: _rel(p, output.output_folder)

    t_meta = table_info.info
    parts = []
    parts.append(
        f'<div class="bp-section">'
        f'<div class="bp-section-title">{html_module.escape(t_meta.name)}</div>'
    )

    # Metadata row
    t_headers = ['name', 'path', 'columns', 'description', 'count']
    t_values = [
        html_module.escape(t_meta.name),
        html_module.escape(rel(t_meta.path) if t_meta.path else ''),
        html_module.escape(", ".join(t_meta.columns) if t_meta.columns else ''),
        html_module.escape(t_meta.description),
        str(t_meta.count),
    ]
    parts.append('<table class="bp-table"><tr>')
    for h in t_headers:
        parts.append(f'<th>{h}</th>')
    parts.append('</tr><tr>')
    for v in t_values:
        parts.append(f'<td>{v}</td>')
    parts.append('</tr></table>')

    # CSV content
    if t_meta.path and os.path.exists(t_meta.path):
        try:
            table_df = pd.read_csv(t_meta.path)
            if len(table_df) > 0:
                parts.append(_render_dataframe(table_df))
        except Exception:
            pass

    parts.append('</div>')
    return "\n".join(parts)
