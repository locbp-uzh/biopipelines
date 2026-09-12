# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""bp-visualize — one self-contained HTML page for a single tool's output.

``RunTime/pipeline.html`` renders every step of a run. This renders one step, with explicit selection and ordering, so a finished step can be pulled off a cluster and opened in a browser while the rest of the job is still queued.

The step is reconstructed from its exported ``ToolOutputs/<step>.json``, which already carries every declared stream, every table, the output folder and the tool's ``rendering_parameters`` — so nothing here needs to know what the tool was.
"""

import argparse
import os
import sys
from typing import Any, Dict, List, Optional, Tuple

try:
    from .datastream import DataStream
    from .outputs import StandardizedOutput
except ImportError:
    sys.path.append(os.path.dirname(__file__))
    from datastream import DataStream
    from outputs import StandardizedOutput


TOOL_OUTPUTS_FOLDER = "ToolOutputs"
INTERNAL_FOLDER = ".internal"
VIEW_SUFFIX = "_view.html"
# Items embedded when --max-items is not given. Higher than the run page's own
# sampling budget (structures.MAX_EMBEDDED, 5) because a single-step page carries one
# step's structures rather than every node's, and because a threshold of 5 hides just
# one item of a 6-item stream -- the least useful place to stop.
DEFAULT_MAX_ITEMS = 10


class VisualizeError(Exception):
    """A bp-visualize failure that carries a message meant for the user."""


# ── locating the step ───────────────────────────────────────────────────────

def locate_tool_outputs_json(tool_folder: str) -> Tuple[str, str]:
    """``(job_root, tool_outputs_json)`` for a tool output folder.

    Walks up until an ancestor holds a ``ToolOutputs/`` directory, which is what makes this work through a ``Folder("group")`` nesting of any depth. An internal tool's manifest lives under ``ToolOutputs/.internal/``, mirroring its own placement.
    """
    tool_folder = os.path.abspath(tool_folder)
    if not os.path.isdir(tool_folder):
        raise VisualizeError(f"not a directory: {tool_folder}")

    step = os.path.basename(tool_folder)
    is_internal = os.path.basename(os.path.dirname(tool_folder)) == INTERNAL_FOLDER

    current = os.path.dirname(tool_folder)
    while True:
        candidate_root = os.path.join(current, TOOL_OUTPUTS_FOLDER)
        if os.path.isdir(candidate_root):
            names = [os.path.join(INTERNAL_FOLDER, f"{step}.json")] if is_internal else [f"{step}.json"]
            # An internal tool nests, but tolerate a flat manifest rather than failing on layout.
            names.append(f"{step}.json")
            for name in names:
                path = os.path.join(candidate_root, name)
                if os.path.isfile(path):
                    return current, path
            raise VisualizeError(
                f"no manifest for step {step!r} in {candidate_root}. "
                f"Present: {', '.join(sorted(os.listdir(candidate_root))) or '(empty)'}"
            )
        parent = os.path.dirname(current)
        if parent == current:
            raise VisualizeError(
                f"no {TOOL_OUTPUTS_FOLDER}/ directory above {tool_folder}. "
                f"Pass a tool output folder inside a job folder, e.g. "
                f"<project>/<job>_001/003_Boltz2"
            )
        current = parent


def load_step(tool_folder: str) -> Tuple[StandardizedOutput, Dict[str, Any], str]:
    """``(output, manifest, job_root)`` for a tool output folder.

    Reuses ``pipeline_report._rebuild_output`` so a step renders here exactly as it does in the pipeline page — including the ``_runtime_mode`` flag that makes lazy ids read the map_table instead of expanding a config-time prefix.
    """
    job_root, manifest_path = locate_tool_outputs_json(tool_folder)
    report = _page_renderer()
    output, manifest, problem = report._rebuild_output({"tool_outputs_json": manifest_path})
    if output is None:
        raise VisualizeError(f"could not rebuild {os.path.basename(tool_folder)}: {problem}")
    return output, manifest or {}, job_root


def _page_renderer():
    from .pipeline import _load_page_renderer
    return _load_page_renderer()


# ── selection and ordering ─────────────────────────────────────────────────

def parse_sort_key(spec: str, output: StandardizedOutput) -> Tuple[str, str]:
    """``("confidence", "plddt")`` from ``"confidence.plddt"``.

    A bare column name is refused rather than searched: two tables of one tool can carry the same column, and picking one silently is the kind of wrong answer that costs an afternoon. The error names the qualified spellings that would work.
    """
    tables = _tables_of(output)
    if "." in spec:
        table_name, _, column = spec.rpartition(".")
        if table_name not in tables:
            raise VisualizeError(
                f"no table {table_name!r} in this step. Available: "
                f"{', '.join(sorted(tables)) or '(none)'}"
            )
        return table_name, column

    candidates = [f"{name}.{spec}" for name, info in tables.items()
                  if spec in (getattr(_meta(info), "columns", None) or [])]
    hint = (f" Did you mean {' or '.join(candidates)}?" if candidates
            else f" Tables in this step: {', '.join(sorted(tables)) or '(none)'}.")
    raise VisualizeError(
        f"sort key {spec!r} must be qualified as <table>.<column>.{hint}"
    )


def _tables_of(output: StandardizedOutput) -> Dict[str, Any]:
    return getattr(output.tables, "_tables", {}) or {}


def _meta(info):
    """A table's metadata object.

    Reading a public attribute straight off a ``TableInfo`` yields a ``TableReference``
    for that column name, whether or not the column exists -- so ``info.columns`` is a
    reference to a column called "columns", not the schema. ``.info`` is the metadata.
    """
    return getattr(info, "info", info)


def ordered_ids(output: StandardizedOutput, table_name: str, column: str,
                descending: bool) -> List[str]:
    """The step's ids ordered by ``table.column``.

    Rows whose value does not parse as a number keep their table order at the end, so a text column still gives a stable, reproducible page instead of raising.
    """
    import pandas as pd

    info = _tables_of(output).get(table_name)
    if info is None:
        raise VisualizeError(f"no table {table_name!r} in this step")
    path = getattr(_meta(info), "path", None)
    if not path or not os.path.isfile(path):
        raise VisualizeError(
            f"table {table_name!r} has not been written yet ({path or 'no path'}); "
            f"the step may still be running"
        )
    df = pd.read_csv(path)
    if column not in df.columns:
        raise VisualizeError(
            f"no column {column!r} in {table_name} ({path}). Columns: "
            f"{', '.join(map(str, df.columns))}"
        )
    if "id" not in df.columns:
        raise VisualizeError(f"table {table_name!r} has no 'id' column to order ids by")

    values = pd.to_numeric(df[column], errors="coerce")
    df = df.assign(_bp_sort=values)
    ranked = df.sort_values("_bp_sort", ascending=not descending, kind="stable",
                            na_position="last")
    seen, out = set(), []
    for raw in ranked["id"].astype(str):
        if raw not in seen:
            seen.add(raw)
            out.append(raw)
    return out


def _reorder_stream(stream: DataStream, order: List[str]) -> DataStream:
    """``stream`` with its ids in ``order``; ids absent from ``order`` keep their relative position at the end."""
    ids = list(stream.ids_expanded)
    if not ids:
        return stream
    rank = {sid: i for i, sid in enumerate(order)}
    files = list(stream.files_expanded) if not stream.is_shared_file else []
    paired = list(zip(ids, files)) if files else [(sid, None) for sid in ids]
    paired.sort(key=lambda pair: (rank.get(pair[0], len(rank)),))
    return _rebuilt(stream, [p[0] for p in paired],
                    [p[1] for p in paired] if files else None)


def _is_ranked(stream: DataStream, order: List[str]) -> bool:
    """Whether the sort actually covers this stream's ids, so a cap has a criterion."""
    ranked = set(order)
    return any(sid in ranked for sid in stream.ids_expanded)


def _cap_stream(stream: DataStream, max_items: int) -> DataStream:
    ids = list(stream.ids_expanded)
    if max_items <= 0 or len(ids) <= max_items:
        return stream
    files = list(stream.files_expanded) if not stream.is_shared_file else []
    return _rebuilt(stream, ids[:max_items], files[:max_items] if files else None)


def _select_stream(stream: DataStream, keep: List[str]) -> DataStream:
    wanted = set(keep)
    ids = list(stream.ids_expanded)
    files = list(stream.files_expanded) if not stream.is_shared_file else []
    picked = [(sid, files[i] if files else None)
              for i, sid in enumerate(ids) if sid in wanted]
    return _rebuilt(stream, [p[0] for p in picked],
                    [p[1] for p in picked] if files else None)


def _rebuilt(stream: DataStream, ids: List[str], files: Optional[List[str]]) -> DataStream:
    """A concrete copy of ``stream`` carrying exactly ``ids``, in that order.

    The map_table is the full warehouse and every renderer prefers it — ``structures._iter_id_file`` and ``streams.render`` both read it before falling back to the id/file zip — so a copy that merely narrowed ``ids`` would still render every row. The path is kept for display, and the loaded frame is projected onto ``ids`` and seeded into the stream's cache, which is also what keeps a value-based stream (compounds, whose content lives only in the map) renderable at all.
    """
    copy = DataStream(
        name=stream.name,
        ids=ids,
        files=stream.files if stream.is_shared_file else (files or []),
        map_table=stream.map_table,
        format=stream.format,
        metadata=dict(stream.metadata),
        _runtime_mode=False,
    )
    copy._map_data = _projected_map(stream, ids)
    return copy


def _projected_map(stream: DataStream, ids: List[str]):
    """``stream``'s map rows for ``ids``, in that order, or None when it has no map."""
    try:
        frame = stream._get_map_data()
    except Exception:
        return None
    if frame is None or len(frame) == 0 or "id" not in frame.columns:
        return None
    indexed = frame.set_index(frame["id"].astype(str), drop=False)
    wanted = [sid for sid in ids if sid in indexed.index]
    if not wanted:
        return frame.iloc[0:0]
    return indexed.loc[wanted].reset_index(drop=True)


def apply_selection(output: StandardizedOutput,
                    ids: Optional[List[str]] = None,
                    sort: Optional[Tuple[str, str]] = None,
                    descending: bool = True,
                    max_items: int = 0,
                    streams: Optional[List[str]] = None,
                    tables: Optional[List[str]] = None) -> Tuple[StandardizedOutput, List[str]]:
    """``output`` narrowed and ordered, plus the notes describing what was done."""
    notes: List[str] = []
    data = dict(output._data)

    order: List[str] = []
    if sort is not None:
        order = ordered_ids(output, sort[0], sort[1], descending)
        notes.append(f"ordered by {sort[0]}.{sort[1]} "
                     f"{'descending' if descending else 'ascending'}")

    for name, stream in list(output.streams.items()):
        if not isinstance(stream, DataStream):
            continue
        if streams is not None and name not in streams:
            data.pop(name, None)
            continue
        working = stream
        if ids:
            working = _select_stream(working, ids)
            if len(working) == 0:
                notes.append(f"stream '{name}' has none of the requested ids")
        if order:
            working = _reorder_stream(working, order)
        # A cap alongside a sort means "the top N of what was ranked". Applying it to a
        # stream the sort table does not cover would drop items on no criterion at all
        # -- a Boltz2 step ranked by confidence cut its 6 sequences and 6 MSAs to 5.
        if max_items and (not order or _is_ranked(working, order)):
            before = len(working.ids_expanded)
            working = _cap_stream(working, max_items)
            if before > max_items:
                notes.append(f"stream '{name}': showing {max_items} of {before} items")
        elif max_items:
            notes.append(f"stream '{name}': shown whole, {sort[0]}.{sort[1]} does not rank it")
        data[name] = working

    if tables is not None:
        declared = _tables_of(output)
        data["tables"] = {k: v for k, v in declared.items() if k in tables}

    if max_items or order or ids:
        params = dict(getattr(output, "rendering_parameters", None) or {})
        if max_items:
            # The structure viewer samples 5 of a large stream by default; an explicit
            # cap is a request, not a hint, so raise its embed budget to match.
            for name in list(data):
                if isinstance(data.get(name), DataStream):
                    per = dict(params.get(name) or {})
                    per["max_embedded"] = max_items
                    params[name] = per
        # Tables read their CSV from disk and know nothing about the streams, so the
        # same selection has to reach them or the page answers one question two ways.
        # The order must be the effective one -- with both --ids and a sort, that is
        # the sort restricted to the requested ids, not the order they were typed in.
        if ids and order:
            wanted = set(ids)
            effective = [sid for sid in order if sid in wanted]
        else:
            effective = list(ids or []) or order
        params["_selection"] = {"ids": effective, "max_items": max_items}
        data["rendering_parameters"] = params

    return StandardizedOutput(data), notes


# ── page assembly ──────────────────────────────────────────────────────────

_SHELL = """<!doctype html>
<meta charset="utf-8">
<title>{title}</title>
{fragment_css}
<style>{css}
body {{ margin: 0; padding: 24px; }}
.bp-view {{ max-width: 1100px; margin: 0 auto; }}
.bp-view h1 {{ font-size: 1.3rem; margin: 0 0 4px; }}
.bp-view .bp-sub {{ color: var(--muted, #6b6b70); font-size: .85rem; margin-bottom: 18px;
                    font-family: ui-monospace, monospace; word-break: break-all; }}
.bp-view .bp-notes {{ border-left: 3px solid #d0d0d5; padding: 8px 12px; margin: 18px 0;
                      font-size: .82rem; color: var(--muted, #6b6b70); }}
.bp-view .bp-notes li {{ margin: 2px 0; }}
</style>
<div class="bp-view">
<h1>{heading}</h1>
<div class="bp-sub">{folder}</div>
{body}
{notes}
</div>
{libs}
"""


def _fragment_css() -> str:
    """The CSS for the classes the stream and table renderers emit.

    ``bp-table``, ``bp-section`` and ``bp-table-toggle`` come out of ``streams.py``,
    ``tables.py`` and ``structures.py``, but ``pipeline_report._CSS`` styles none of
    them -- it covers the page chrome only. ``render_page`` therefore emits two style
    blocks, its own and ``StandardizedOutput._CSS``, and a page carrying only the first
    shows every renderer table unstyled. Same source as the run page and the notebook,
    so all three agree.
    """
    return StandardizedOutput._CSS


def _notes_html(notes: List[str], esc) -> str:
    if not notes:
        return ""

    items = "".join(f"<li>{esc(n)}</li>" for n in notes)
    return f'<div class="bp-notes"><ul>{items}</ul></div>'


def build_page(output: StandardizedOutput, tool_folder: str, manifest: Dict[str, Any],
               extra_notes: Optional[List[str]] = None,
               allow_external: bool = False) -> str:
    """The complete HTML for one step. Self-contained: no network, no sibling files."""
    report = _page_renderer()
    notes = list(extra_notes or [])
    config = report._resolve_renderers_config(
        (manifest.get("export_metadata") or {}).get("config_variant"), notes)
    libs = set()
    body = report._render_stream_bodies(output, config, notes, allow_external, libs)
    if not body:
        notes.append("no renderer produced output for this step")

    tool = manifest.get("tool_name") or os.path.basename(tool_folder)
    step = os.path.basename(tool_folder)
    return _SHELL.format(
        title=report._esc(f"{step} — {tool}"),
        heading=report._esc(step),
        folder=report._esc(tool_folder),
        css=report._CSS,
        fragment_css=_fragment_css(),
        body=body,
        notes=_notes_html(notes, report._esc),
        libs=report._library_scripts(libs),
    )


def default_output_path(tool_folder: str) -> str:
    """``<tool folder>/_extras/<step>_view.html`` — beside the outputs it describes.

    ``_extras/`` is the canonical catch-all, and the completion check is declaration-driven (it never enumerates a folder), so a view written there cannot affect a step's verdict.
    """
    step = os.path.basename(os.path.abspath(tool_folder))
    return os.path.join(os.path.abspath(tool_folder), "_extras", f"{step}{VIEW_SUFFIX}")


def write_view(tool_folder: str, out_path: Optional[str] = None,
               sort_spec: Optional[str] = None, **selection) -> str:
    """Render one step to a self-contained page and return the path written.

    ``sort_spec`` is the unparsed ``<table>.<column>`` string; it is resolved here because naming the step's real tables in the error message needs the loaded output.

    ``max_items`` defaults to 0, meaning no cap -- each renderer keeps its own sampling. The CLI passes ``DEFAULT_MAX_ITEMS`` instead, so ``bp-visualize <folder>`` embeds 10 while ``write_view(folder)`` leaves the default sampling alone. The opinion belongs to the command, not to the function.
    """
    allow_external = selection.pop("allow_external", False)
    output, manifest, _job_root = load_step(tool_folder)
    if sort_spec:
        selection["sort"] = parse_sort_key(sort_spec, output)
    narrowed, notes = apply_selection(output, **selection)
    html = build_page(narrowed, os.path.abspath(tool_folder), manifest, notes,
                      allow_external=allow_external)
    target = os.path.abspath(out_path or default_output_path(tool_folder))
    os.makedirs(os.path.dirname(target), exist_ok=True)
    with open(target, "w", encoding="utf-8") as f:
        f.write(html)
    return target


# ── CLI ────────────────────────────────────────────────────────────────────

def _csv_list(value: str) -> List[str]:
    return [part.strip() for part in value.split(",") if part.strip()]


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="bp-visualize",
        description="Render one tool's output as a self-contained HTML page.",
        epilog=(
            "Examples:\n"
            "  bp-visualize <job>/003_Boltz2\n"
            "  bp-visualize <job>/003_Boltz2 --descending confidence.confidence_score --max-items 5\n"
            "  bp-visualize <job>/002_ESMFold --ascending confidence.plddt --streams structures\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("folder", help="the tool's output folder, e.g. <job>/003_Boltz2")
    parser.add_argument("-o", "--output", default=None,
                        help="where to write the page (default: <folder>/_extras/<step>_view.html)")
    parser.add_argument("--max-items", type=int, default=DEFAULT_MAX_ITEMS, metavar="N",
                        help=f"render at most N items (default {DEFAULT_MAX_ITEMS}); with an "
                             "ordering, caps only the streams that ordering ranks. 0 means no "
                             "cap, leaving each renderer its own sampling")
    order = parser.add_mutually_exclusive_group()
    order.add_argument("--descending", metavar="TABLE.COLUMN", default=None,
                       help="order items by a table column, highest first")
    order.add_argument("--ascending", metavar="TABLE.COLUMN", default=None,
                       help="order items by a table column, lowest first")
    parser.add_argument("--ids", type=_csv_list, default=None, metavar="ID,ID",
                        help="render only these ids")
    parser.add_argument("--streams", type=_csv_list, default=None, metavar="NAME,NAME",
                        help="render only these streams")
    parser.add_argument("--tables", type=_csv_list, default=None, metavar="NAME,NAME",
                        help="render only these tables")
    parser.add_argument("--allow-external", action="store_true",
                        help="keep a renderer's CDN <script> when no vendored copy exists")
    parser.add_argument("--open", action="store_true", dest="open_browser",
                        help="open the page in a browser (for local runs, not a login node)")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    args = _build_parser().parse_args(argv)

    try:
        path = write_view(
            args.folder,
            args.output,
            sort_spec=args.descending or args.ascending,
            ids=args.ids,
            descending=bool(args.descending),
            max_items=args.max_items,
            streams=args.streams,
            tables=args.tables,
            allow_external=args.allow_external,
        )
    except VisualizeError as e:
        print(f"bp-visualize: {e}", file=sys.stderr)
        return 1

    size_kb = os.path.getsize(path) // 1024
    print(f"Wrote {path} ({size_kb} kB)")
    if args.open_browser:
        import webbrowser
        webbrowser.open(f"file://{path.replace(os.sep, '/')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
