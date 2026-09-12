# Default renderer for tables: metadata + CSV content

import os
import html as html_module
import pandas as pd


def _rel(path, output_folder):
    """Make path relative to output_folder for display."""
    if output_folder and path.startswith(output_folder):
        return "<output_folder>/" + os.path.relpath(path, output_folder)
    return path


_MAX_EXPANDED_ROWS = 50


def _render_rows(df, indices, columns, output_folder=None):
    parts = []
    for idx in indices:
        row = df.iloc[idx]
        parts.append('<tr>')
        for col in columns:
            val = str(row[col]) if pd.notna(row[col]) else ''
            if output_folder and val.startswith(output_folder):
                val = _rel(val, output_folder)
            display_val = val if len(val) <= 60 else val[:57] + '...'
            parts.append(f'<td>{html_module.escape(display_val)}</td>')
        parts.append('</tr>')
    return parts


def _render_dataframe(df, output_folder=None):
    """Render a DataFrame as a collapsible HTML table.

    Collapsed: 2 head + 2 tail rows. Expanded: up to 50 rows (2 head + 46
    middle + 2 tail with an inner ellipsis when the table is larger).
    """
    columns = list(df.columns)
    n_rows = len(df)
    n_cols = len(columns)

    header = ['<tr>']
    for col in columns:
        header.append(f'<th>{html_module.escape(str(col))}</th>')
    header.append('</tr>')
    header_html = "".join(header)

    if n_rows <= 4:
        parts = ['<table class="bp-table">', header_html]
        parts.extend(_render_rows(df, range(n_rows), columns, output_folder))
        parts.append('</table>')
        return "\n".join(parts)

    collapsed_indices = list(range(2)) + list(range(n_rows - 2, n_rows))
    collapsed_parts = ['<table class="bp-table">', header_html]
    collapsed_parts.extend(_render_rows(df, collapsed_indices[:2], columns, output_folder))
    collapsed_parts.append(
        f'<tr class="bp-ellipsis"><td colspan="{n_cols}">'
        f'... {n_rows - 4} more ...</td></tr>'
    )
    collapsed_parts.extend(_render_rows(df, collapsed_indices[2:], columns, output_folder))
    collapsed_parts.append('</table>')
    collapsed_html = "\n".join(collapsed_parts)

    if n_rows <= _MAX_EXPANDED_ROWS:
        expanded_indices = list(range(n_rows))
        expanded_ellipsis_at = None
    else:
        head_n = 2
        tail_n = 2
        middle_n = _MAX_EXPANDED_ROWS - head_n - tail_n
        expanded_indices = (
            list(range(head_n + middle_n))
            + list(range(n_rows - tail_n, n_rows))
        )
        expanded_ellipsis_at = head_n + middle_n

    expanded_parts = ['<table class="bp-table">', header_html]
    for row_i, idx in enumerate(expanded_indices):
        if expanded_ellipsis_at is not None and row_i == expanded_ellipsis_at:
            expanded_parts.append(
                f'<tr class="bp-ellipsis"><td colspan="{n_cols}">'
                f'... {n_rows - _MAX_EXPANDED_ROWS} more ...</td></tr>'
            )
        row = df.iloc[idx]
        expanded_parts.append('<tr>')
        for col in columns:
            val = str(row[col]) if pd.notna(row[col]) else ''
            if output_folder and val.startswith(output_folder):
                val = _rel(val, output_folder)
            display_val = val if len(val) <= 60 else val[:57] + '...'
            expanded_parts.append(f'<td>{html_module.escape(display_val)}</td>')
        expanded_parts.append('</tr>')
    expanded_parts.append('</table>')
    expanded_html = "\n".join(expanded_parts)

    summary_label = (
        f'▸ expand ({n_rows} rows, showing up to {_MAX_EXPANDED_ROWS})'
    )
    return (
        '<details class="bp-table-toggle">'
        f'<summary>{summary_label}</summary>'
        f'{expanded_html}'
        '</details>'
        f'<div class="bp-table-collapsed">{collapsed_html}</div>'
    )


def _project(df, output):
    """``df`` narrowed to the caller's selection, if one was requested.

    bp-visualize puts the ordered id list and the item cap in
    ``rendering_parameters["_selection"]``. Without this a page asking for the top 5 by
    pLDDT showed 5 structures beside a table of every row in file order -- the same page
    answering the question two different ways.
    """
    selection = (getattr(output, "rendering_parameters", None) or {}).get("_selection") or {}
    ids = selection.get("ids")
    cap = selection.get("max_items")
    if ids and "id" in df.columns:
        rank = {str(v): i for i, v in enumerate(ids)}
        keep = df[df["id"].astype(str).isin(rank)].copy()
        if len(keep) > 0:
            keep["_bp_rank"] = keep["id"].astype(str).map(rank)
            df = keep.sort_values("_bp_rank", kind="stable").drop(columns="_bp_rank")
    if isinstance(cap, int) and cap > 0:
        df = df.head(cap)
    return df


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
    t_headers = ['name', 'path', 'columns', 'description']
    t_values = [
        html_module.escape(t_meta.name),
        html_module.escape(rel(t_meta.path) if t_meta.path else ''),
        html_module.escape(", ".join(t_meta.columns) if t_meta.columns else ''),
        html_module.escape(t_meta.description),
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
            table_df = _project(pd.read_csv(t_meta.path), output)
            if len(table_df) > 0:
                parts.append(_render_dataframe(table_df, output_folder=output.output_folder))
        except Exception:
            pass

    parts.append('</div>')
    return "\n".join(parts)
