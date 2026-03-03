# Default renderer for streams: metadata table + data table

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


def render(stream, output):
    """Render stream metadata and data table."""
    rel = lambda p: _rel(p, output.output_folder)

    parts = []
    parts.append(
        f'<div class="bp-section">'
        f'<div class="bp-section-title">{html_module.escape(stream.name)}</div>'
    )

    # Metadata table (horizontal)
    meta_headers = ['name', 'format', 'items', 'map_table', 'files_contain_wildcards']
    meta_values = [
        html_module.escape(stream.name),
        html_module.escape(stream.format),
        str(len(stream)),
        html_module.escape(rel(stream.map_table) if stream.map_table else ''),
        str(stream.files_contain_wildcards),
    ]
    if stream.metadata:
        meta_headers.append('metadata')
        meta_values.append(html_module.escape(
            ", ".join(f"{k}={v}" for k, v in stream.metadata.items())
        ))
    parts.append('<table class="bp-table"><tr>')
    for h in meta_headers:
        parts.append(f'<th>{h}</th>')
    parts.append('</tr><tr>')
    for v in meta_values:
        parts.append(f'<td>{v}</td>')
    parts.append('</tr></table>')

    # Data: map_table or id/file fallback
    map_data = stream._get_map_data()
    if map_data is not None and len(map_data) > 0:
        parts.append(_render_dataframe(map_data))
    else:
        items = list(zip(stream.ids, stream.files)) if stream.files else [(iid, "") for iid in stream.ids]
        n_items = len(items)
        parts.append('<table class="bp-table"><tr><th>id</th><th>file</th></tr>')
        if n_items <= 4:
            display_items = items
            ellipsis_after = None
        else:
            display_items = items[:2] + items[-2:]
            ellipsis_after = 2

        for item_i, (item_id, item_file) in enumerate(display_items):
            if ellipsis_after is not None and item_i == ellipsis_after:
                parts.append(
                    f'<tr class="bp-ellipsis"><td colspan="2">'
                    f'... {n_items - 4} more ...</td></tr>'
                )
            parts.append(
                f"<tr><td>{html_module.escape(str(item_id))}</td>"
                f"<td>{html_module.escape(rel(item_file) if item_file else '')}</td></tr>"
            )
        parts.append('</table>')

    parts.append('</div>')
    return "\n".join(parts)
