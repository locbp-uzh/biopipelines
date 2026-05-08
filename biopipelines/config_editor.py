# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Interactive TUI for editing a biopipelines config YAML.

Invoked by ``bp-config edit``. Loads the active ``config.<variant>.yaml``
with ruamel.yaml in round-trip mode (so comments and key ordering survive
a save), renders the tree with prompt_toolkit, and lets the user navigate
with arrow keys, expand/collapse branches, edit leaf values, search, and
save back to disk atomically (with a ``.bak`` alongside).

Bindings:
  ↑ / ↓        move cursor
  → / Enter    expand a collapsed branch, or edit a leaf
  ← / -        collapse the current branch (or move to parent)
  /            search by key name (case-insensitive substring)
  n / N        next / previous search hit
  s            save (writes ``<path>.bak`` first)
  q            quit (prompts if there are unsaved changes)
  ?            show help overlay
"""

from __future__ import annotations

import io
import os
import shutil
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, List, Optional, Tuple


# Sentinel used in row.value to mean "this row has no leaf value (it's a
# branch)". We can't use None because YAML legitimately stores None values.
_NO_VALUE = object()


@dataclass
class Row:
    """One visible line in the tree view.

    path is the sequence of keys from the document root to this node;
    e.g. ('folders', 'output') for a leaf, ('folders',) for a branch.
    For sequence elements the keys are integers.
    """
    path: Tuple[Any, ...]
    label: str                  # rendered key (str(key))
    value: Any = _NO_VALUE      # leaf scalar, or _NO_VALUE for branches
    has_children: bool = False
    expanded: bool = False
    depth: int = 0


@dataclass
class EditorState:
    file_path: Path
    doc: Any                    # ruamel.yaml round-trip CommentedMap/CommentedSeq
    rows: List[Row] = field(default_factory=list)
    cursor: int = 0
    dirty: bool = False
    status: str = ""
    search_term: str = ""
    show_help: bool = False


def _walk(node: Any, path: Tuple[Any, ...], depth: int,
          expanded_paths: set, out: List[Row]) -> None:
    """Recursively flatten the doc into rows, honouring expanded_paths."""
    from ruamel.yaml.comments import CommentedMap, CommentedSeq

    if isinstance(node, (dict, CommentedMap)):
        # Display dict entries alphabetically — most useful for sections
        # that grow (containers, environments, emails, tool_overrides);
        # the underlying YAML order on disk is preserved by ruamel and
        # not affected by this view-only sort.
        items = sorted(node.items(), key=lambda kv: str(kv[0]))
    elif isinstance(node, (list, CommentedSeq)):
        items = list(enumerate(node))
    else:
        return

    for key, value in items:
        child_path = path + (key,)
        is_branch = isinstance(value, (dict, list, CommentedMap, CommentedSeq))
        row = Row(
            path=child_path,
            label=str(key),
            value=_NO_VALUE if is_branch else value,
            has_children=is_branch and bool(value),
            expanded=is_branch and child_path in expanded_paths,
            depth=depth,
        )
        out.append(row)
        if is_branch and row.expanded:
            _walk(value, child_path, depth + 1, expanded_paths, out)


def _build_rows(state: EditorState, expanded_paths: set) -> None:
    state.rows = []
    _walk(state.doc, (), 0, expanded_paths, state.rows)


def _node_at_path(doc: Any, path: Tuple[Any, ...]) -> Any:
    node = doc
    for k in path:
        node = node[k]
    return node


def _set_leaf(doc: Any, path: Tuple[Any, ...], new_value: Any) -> None:
    parent = _node_at_path(doc, path[:-1])
    parent[path[-1]] = new_value


def _coerce(text: str, original: Any) -> Any:
    """Best-effort: keep the original Python type when the user edits a leaf.

    For str we just take text. For bool we accept true/false/yes/no
    (case-insensitive). For int/float we parse; failure raises ValueError.
    For None we accept empty string or 'null' as None, otherwise return the
    raw string (the user may be turning a None into a string deliberately).
    """
    if isinstance(original, bool):
        low = text.strip().lower()
        if low in ("true", "yes", "y", "1"):
            return True
        if low in ("false", "no", "n", "0"):
            return False
        raise ValueError(f"expected boolean, got {text!r}")
    if isinstance(original, int) and not isinstance(original, bool):
        return int(text)
    if isinstance(original, float):
        return float(text)
    if original is None:
        if text.strip() == "" or text.strip().lower() == "null":
            return None
        return text
    return text


def _render_value(v: Any) -> str:
    if v is None:
        return "null"
    if isinstance(v, bool):
        return "true" if v else "false"
    return str(v)


def run_editor(file_path: str) -> int:
    """Open file_path in the TUI editor. Returns shell exit code."""
    from prompt_toolkit import Application
    from prompt_toolkit.key_binding import KeyBindings
    from prompt_toolkit.layout import Layout, Window, HSplit, FormattedTextControl
    from prompt_toolkit.layout.dimension import Dimension
    from prompt_toolkit.styles import Style
    from prompt_toolkit.formatted_text import FormattedText
    from prompt_toolkit.widgets import TextArea
    from ruamel.yaml import YAML

    path = Path(file_path)
    if not path.exists():
        print(f"Config file not found: {path}", file=sys.stderr)
        return 1

    yaml = YAML()
    yaml.preserve_quotes = True
    with path.open("r", encoding="utf-8") as f:
        doc = yaml.load(f)

    if doc is None:
        print(f"Config file is empty: {path}", file=sys.stderr)
        return 1

    state = EditorState(file_path=path, doc=doc)
    # Start with every branch collapsed so the user sees just the
    # top-level section names; they can expand the one they want.
    expanded: set = set()
    from ruamel.yaml.comments import CommentedMap, CommentedSeq
    _build_rows(state, expanded)

    # ---- missing-path detection -------------------------------------
    # Filesystem paths in `folders.*` and `containers.*` get a red
    # marker if they don't resolve to an existing entry on the local
    # filesystem. Resolution is cached per session and refreshed on
    # save (and via R). Environment names (in `environments.*`) and
    # cluster modules are NOT checked — too cluster-specific to do
    # reliably from the login node and slow if done per render.
    import re as _re
    _PLACEHOLDER = _re.compile(r"<([^>]+)>")
    missing_paths: set = set()  # set of row.path tuples

    def _flatten_folders(folders_block: Any) -> Dict[str, str]:
        """Flat mapping of folder key -> raw path string from the YAML."""
        out: Dict[str, str] = {}
        if not isinstance(folders_block, (dict, CommentedMap)):
            return out
        for section, entries in folders_block.items():
            if isinstance(entries, (dict, CommentedMap)):
                for k, v in entries.items():
                    if isinstance(v, str):
                        out[str(k)] = v
        return out

    def _resolve_template(tmpl: str, lookup: Dict[str, str]) -> Optional[str]:
        """Substitute <key> placeholders against lookup. Returns None if any
        placeholder is unresolvable. Tolerates runtime placeholders like
        <username> / <cwd> / <user> by leaving them in place — the
        existence check below short-circuits to "unknown" for those."""
        result = tmpl
        for _ in range(10):
            matches = _PLACEHOLDER.findall(result)
            if not matches:
                break
            progressed = False
            for ph in matches:
                if ph in lookup and lookup[ph] != tmpl:
                    result = result.replace(f"<{ph}>", lookup[ph])
                    progressed = True
            if not progressed:
                break
        return result

    def refresh_missing() -> None:
        missing_paths.clear()
        if not isinstance(state.doc, (dict, CommentedMap)):
            return
        folders = state.doc.get("folders") or {}
        flat = _flatten_folders(folders)
        # Walk folders.<section>.<key> entries.
        for section_name, entries in folders.items():
            if not isinstance(entries, (dict, CommentedMap)):
                continue
            for k, v in entries.items():
                if not isinstance(v, str):
                    continue
                resolved = _resolve_template(v, flat) or ""
                # Skip rows that still have unresolved <runtime> placeholders
                # — we don't know the runtime values from inside the editor.
                if "<" in resolved:
                    continue
                if resolved and not os.path.exists(resolved):
                    missing_paths.add(("folders", str(section_name), str(k)))
        # Walk containers.<Tool> entries.
        containers = state.doc.get("containers") or {}
        if isinstance(containers, (dict, CommentedMap)):
            for k, v in containers.items():
                if not isinstance(v, str):
                    continue
                resolved = _resolve_template(v, flat) or ""
                if "<" in resolved:
                    continue
                if resolved and not os.path.exists(resolved):
                    missing_paths.add(("containers", str(k)))

    refresh_missing()

    # ---- rendering ----------------------------------------------------
    def header_text() -> FormattedText:
        dirty = " [modified]" if state.dirty else ""
        return FormattedText([("class:header",
                               f" biopipelines config: {state.file_path}{dirty} ")])

    # Vertical scroll offset: row index of the topmost visible line.
    # Updated on every render to keep the cursor row visible inside
    # whatever vertical space the body window has at runtime.
    view = {"top": 0}

    def _visible_rows() -> int:
        """How many rows fit in the body window right now."""
        from prompt_toolkit.application.current import get_app
        try:
            total = get_app().output.get_size().rows
        except Exception:
            total = 24
        # Subtract header (1) + footer (3 in help mode, 1 otherwise) +
        # the optional inline edit / search line. Be conservative: leave
        # at least 5 rows.
        chrome = 1 + (3 if state.show_help else 1)
        if editing_path or adding_path or searching[0]:
            chrome += 1
        return max(5, total - chrome)

    def body_text() -> FormattedText:
        # Clamp view_top so cursor is in view.
        n = len(state.rows)
        h = _visible_rows()
        if n == 0:
            view["top"] = 0
        else:
            if state.cursor < view["top"]:
                view["top"] = state.cursor
            elif state.cursor >= view["top"] + h:
                view["top"] = state.cursor - h + 1
            view["top"] = max(0, min(view["top"], max(0, n - h)))
        top = view["top"]
        bottom = min(n, top + h)
        out = []
        for i in range(top, bottom):
            row = state.rows[i]
            indent = "  " * row.depth
            if row.has_children:
                marker = "▾ " if row.expanded else "▸ "
            else:
                marker = "  "
            line_left = f"{indent}{marker}{row.label}"
            if row.value is _NO_VALUE:
                line = line_left
            else:
                line = f"{line_left}: {_render_value(row.value)}"
            style = "class:cursor" if i == state.cursor else ""
            if state.search_term and state.search_term.lower() in row.label.lower():
                style = (style + " class:match").strip()
            if row.path in missing_paths:
                style = (style + " class:missing").strip()
            out.append((style, line + "\n"))
        if n > bottom:
            out.append(("class:scrollbar", f"… {n - bottom} more below\n"))
        if top > 0:
            out.insert(0, ("class:scrollbar", f"… {top} more above\n"))
        if not out:
            out.append(("", "(empty document)\n"))
        return FormattedText(out)

    def footer_text() -> FormattedText:
        if state.show_help:
            help_lines = [
                "  ↑↓ navigate   →/Enter expand or edit value   ←/- collapse or parent",
                "  a add entry   d delete entry   r rename key   R refresh missing-path check",
                "  / search   n/N next/prev   s save   q quit   ? toggle help",
                "  red rows = folder/container path not found on local filesystem",
            ]
            return FormattedText([("class:footer", "\n".join(help_lines))])
        hint = (" ↑↓ move  →/Enter open  ←/- close  a add  d del  r rename  "
                "/ search  s save  q quit  ? help ")
        if state.status:
            hint += f"| {state.status}"
        return FormattedText([("class:footer", hint)])

    header = Window(content=FormattedTextControl(header_text), height=1,
                    style="class:header")
    body = Window(content=FormattedTextControl(body_text, focusable=True),
                  always_hide_cursor=True)
    # Footer height adapts to whether help is showing (4 lines) or not
    # (1 status line). Recomputed on every render via a lambda so the
    # body re-sizes when `?` toggles help on/off.
    from prompt_toolkit.layout.dimension import Dimension as _Dim
    footer = Window(
        content=FormattedTextControl(footer_text),
        height=lambda: _Dim.exact(4 if state.show_help else 1),
        style="class:footer",
    )
    root_container = HSplit([header, body, footer])
    layout = Layout(root_container)

    style = Style.from_dict({
        "header": "reverse",
        "footer": "reverse",
        "cursor": "reverse",
        "match": "bold underline",
        # Resolved path / container image not found on the local
        # filesystem. Marker only; the YAML still saves whatever the
        # user types.
        "missing": "fg:ansired",
        "scrollbar": "fg:ansigray italic",
    })

    # ---- helpers used by bindings ------------------------------------
    def current_row() -> Optional[Row]:
        if not state.rows:
            return None
        if not (0 <= state.cursor < len(state.rows)):
            state.cursor = max(0, min(state.cursor, len(state.rows) - 1))
        return state.rows[state.cursor]

    def rebuild() -> None:
        # Preserve cursor on the same path if possible after rebuild.
        target = state.rows[state.cursor].path if state.rows else None
        _build_rows(state, expanded)
        if target is not None:
            for i, r in enumerate(state.rows):
                if r.path == target:
                    state.cursor = i
                    break
            else:
                state.cursor = min(state.cursor, max(0, len(state.rows) - 1))

    def save_to_disk() -> None:
        bak = state.file_path.with_suffix(state.file_path.suffix + ".bak")
        shutil.copy2(state.file_path, bak)
        # Atomic write: dump to a temp in the same dir, then replace.
        buf = io.StringIO()
        yaml.dump(state.doc, buf)
        tmp = state.file_path.with_suffix(state.file_path.suffix + ".tmp")
        tmp.write_text(buf.getvalue(), encoding="utf-8")
        os.replace(tmp, state.file_path)
        state.dirty = False
        state.status = f"saved (backup: {bak.name})"

    # ---- inline edit prompt ------------------------------------------
    # Enter / Escape handlers for this textarea live on the app-level
    # KeyBindings (`kb`) below, gated by a has_focus filter. We can't
    # use TextArea(key_bindings=...) — that kwarg only exists on newer
    # prompt_toolkit releases — and TextArea.control.key_bindings is
    # None by default, so directly mutating it AttributeErrors.
    edit_textarea = TextArea(height=1, multiline=False, prompt="value: ")

    # Editor state machine. At most one of editing_path / adding_path is
    # active at a time. The associated "_mode" indicates which sub-state
    # the prompt is in.
    #
    #   editing_path + edit_mode="value"  -> edit a leaf's value
    #   editing_path + edit_mode="key"    -> rename a dict-leaf's key
    #   adding_path  + add_mode="key"     -> first prompt of dict-add
    #                                        (asks for the new key)
    #   adding_path  + add_mode="value"   -> second prompt of dict-add,
    #                                        or list-add's only prompt
    editing_path: List[Tuple[Any, ...]] = []
    edit_mode: List[str] = []
    adding_path: List[Tuple[Any, ...]] = []
    add_mode: List[str] = []
    # Stashed key from the first prompt of a dict-add, consumed by the
    # second prompt to assemble the final {key: value} entry.
    pending_key: List[str] = []

    def _is_dict_entry(path: Tuple[Any, ...]) -> bool:
        """True if path points at a leaf whose parent is a dict (vs a list)."""
        if not path:
            return False
        from ruamel.yaml.comments import CommentedMap
        parent = _node_at_path(state.doc, path[:-1])
        return isinstance(parent, (dict, CommentedMap))

    def _branch_kind(path: Tuple[Any, ...]) -> Optional[str]:
        """Return 'dict' / 'list' for a branch path, or None for leaves."""
        from ruamel.yaml.comments import CommentedMap, CommentedSeq
        try:
            node = _node_at_path(state.doc, path) if path else state.doc
        except (KeyError, IndexError, TypeError):
            return None
        if isinstance(node, (dict, CommentedMap)):
            return "dict"
        if isinstance(node, (list, CommentedSeq)):
            return "list"
        return None

    def _reset_prompt() -> None:
        editing_path.clear()
        edit_mode.clear()
        adding_path.clear()
        add_mode.clear()
        pending_key.clear()
        layout.container = root_container
        layout.focus(body)

    def _open_prompt(prompt_text: str, initial: str = "") -> None:
        """Show the inline edit textarea with the given prompt label and
        initial text. The textarea is reused across edit / rename / add
        prompts — its prompt is the BeforeInput processor's ``text``."""
        edit_textarea.text = initial
        for proc in edit_textarea.control.input_processors or []:
            if hasattr(proc, "text"):
                proc.text = prompt_text
                break
        edit_textarea.buffer.cursor_position = len(initial)
        layout.container = HSplit([header, body, edit_textarea, footer])
        layout.focus(edit_textarea)

    def begin_edit(row: Row) -> None:
        """Edit a leaf's value (Enter on a leaf row)."""
        _reset_prompt()
        editing_path.append(row.path)
        edit_mode.append("value")
        _open_prompt("value: ", _render_value(row.value))

    def begin_rename(row: Row) -> None:
        """Rename a dict-leaf's key (`r` on a dict-leaf row)."""
        if not _is_dict_entry(row.path):
            state.status = "rename only applies to dict entries"
            return
        _reset_prompt()
        editing_path.append(row.path)
        edit_mode.append("key")
        _open_prompt("new key: ", str(row.label))

    def begin_add(target_path: Tuple[Any, ...]) -> None:
        """Add a new entry under target_path.

        For dict targets: opens a two-prompt flow — first the key,
        committing it advances to a second prompt for the value.
        For list targets: opens a single value prompt; commit appends.

        If target_path points at a leaf whose value is None, replace it
        with an empty dict in place first (handles the ``containers:
        null`` case where the user wants to start adding entries).
        """
        # If target_path is a None-valued leaf, promote it to {}.
        if target_path:
            try:
                node = _node_at_path(state.doc, target_path)
            except (KeyError, IndexError, TypeError):
                node = "missing"
            if node is None:
                from ruamel.yaml.comments import CommentedMap
                _set_leaf(state.doc, target_path, CommentedMap())
                state.dirty = True
        kind = _branch_kind(target_path)
        if kind is None:
            state.status = "cannot add here: not a dict or list"
            return
        _reset_prompt()
        adding_path.append(target_path)
        if kind == "dict":
            add_mode.append("key")
            _open_prompt("new key: ", "")
        else:
            add_mode.append("value")
            _open_prompt("new value: ", "")

    def end_edit(commit: bool) -> None:
        # ---- adding flow -------------------------------------------------
        if adding_path:
            target = adding_path[0]
            mode = add_mode[0] if add_mode else "value"
            text = edit_textarea.text
            if not commit:
                _reset_prompt()
                state.status = "add cancelled"
                return
            if mode == "key":
                new_key = text.strip()
                if not new_key:
                    state.status = "key cannot be empty"
                    _reset_prompt()
                    return
                parent = _node_at_path(state.doc, target) if target else state.doc
                if new_key in parent:
                    state.status = f"key {new_key!r} already exists; edit it instead"
                    _reset_prompt()
                    return
                # Stash the key and re-prompt for the value.
                pending_key.append(new_key)
                add_mode.clear()
                add_mode.append("value")
                _open_prompt(f"value for {new_key!r}: ", "")
                return
            # mode == "value"
            try:
                new_val = _coerce(text, "")
            except ValueError as e:
                state.status = f"invalid value: {e}"
                _reset_prompt()
                return
            parent = _node_at_path(state.doc, target) if target else state.doc
            if pending_key:
                key = pending_key[0]
                parent[key] = new_val
                state.dirty = True
                state.status = f"added {'.'.join(map(str, target + (key,)))}"
            else:
                # list append
                parent.append(new_val)
                state.dirty = True
                state.status = f"appended to {'.'.join(map(str, target))}"
            rebuild()
            _reset_prompt()
            return

        # ---- editing flow ------------------------------------------------
        if not editing_path:
            _reset_prompt()
            return
        if not commit:
            _reset_prompt()
            state.status = "edit cancelled"
            return
        target_path = editing_path[0]
        row = next((r for r in state.rows if r.path == target_path), None)
        if row is None:
            _reset_prompt()
            return
        mode = edit_mode[0] if edit_mode else "value"
        text = edit_textarea.text
        if mode == "key":
            new_key = text.strip()
            if not new_key:
                state.status = "key cannot be empty"
                _reset_prompt()
                return
            old_key = target_path[-1]
            if new_key == old_key:
                _reset_prompt()
                state.status = "key unchanged"
                return
            parent = _node_at_path(state.doc, target_path[:-1])
            if new_key in parent:
                state.status = f"key {new_key!r} already exists"
                _reset_prompt()
                return
            items = list(parent.items())
            parent.clear()
            for k, v in items:
                parent[new_key if k == old_key else k] = v
            state.dirty = True
            state.status = f"renamed {'.'.join(map(str, target_path))} -> {new_key}"
            rebuild()
            _reset_prompt()
            return
        # mode == "value"
        try:
            new_val = _coerce(text, row.value)
        except ValueError as e:
            state.status = f"invalid value: {e}"
            _reset_prompt()
            return
        _set_leaf(state.doc, target_path, new_val)
        state.dirty = True
        state.status = f"set {'.'.join(map(str, target_path))}"
        rebuild()
        _reset_prompt()

    # ---- key bindings ------------------------------------------------
    # Letter / punctuation bindings (`-`, `s`, `?`, `q`, `Q`, `/`, `n`,
    # `N`) must be silenced while the user is typing inside the inline
    # edit or search textareas — otherwise typing 'q' inside an email
    # value triggers the quit handler. Gate every printable-character
    # binding on `_body_focused` so they only fire from the tree view.
    from prompt_toolkit.filters import has_focus
    kb = KeyBindings()
    _body_focused = has_focus(body)

    @kb.add("up", filter=_body_focused)
    def _(event):
        if state.cursor > 0:
            state.cursor -= 1
            state.status = ""

    @kb.add("down", filter=_body_focused)
    def _(event):
        if state.cursor < len(state.rows) - 1:
            state.cursor += 1
            state.status = ""

    @kb.add("pageup", filter=_body_focused)
    def _(event):
        state.cursor = max(0, state.cursor - 10)

    @kb.add("pagedown", filter=_body_focused)
    def _(event):
        state.cursor = min(len(state.rows) - 1, state.cursor + 10)

    @kb.add("home", filter=_body_focused)
    def _(event):
        state.cursor = 0

    @kb.add("end", filter=_body_focused)
    def _(event):
        state.cursor = max(0, len(state.rows) - 1)

    @kb.add("right", filter=_body_focused)
    @kb.add("enter", filter=_body_focused)
    def _(event):
        row = current_row()
        if row is None:
            return
        if row.has_children and not row.expanded:
            expanded.add(row.path)
            rebuild()
            state.status = ""
        elif not row.has_children:
            begin_edit(row)
        # If branch is already expanded, do nothing.

    @kb.add("left", filter=_body_focused)
    @kb.add("-", filter=_body_focused)
    def _(event):
        row = current_row()
        if row is None:
            return
        if row.expanded:
            expanded.discard(row.path)
            rebuild()
        elif row.depth > 0:
            # Move to parent.
            parent_path = row.path[:-1]
            for i, r in enumerate(state.rows):
                if r.path == parent_path:
                    state.cursor = i
                    break
        state.status = ""

    @kb.add("s", filter=_body_focused)
    def _(event):
        try:
            save_to_disk()
            refresh_missing()
        except Exception as e:
            state.status = f"save failed: {e}"

    @kb.add("R", filter=_body_focused)
    def _(event):
        """Re-check which folder / container paths are missing."""
        refresh_missing()
        n = len(missing_paths)
        state.status = f"refreshed: {n} missing path{'s' if n != 1 else ''}"

    @kb.add("a", filter=_body_focused)
    def _(event):
        """Add a new entry.

        * Branch row: add under that branch.
        * Leaf row whose value is ``None``: promote it in place to an
          empty dict and add inside (handles ``containers: null``).
        * Other leaf row: add a sibling under the leaf's parent.
        """
        row = current_row()
        if row is None:
            begin_add(())  # empty doc: add at root
            return
        if row.has_children:
            begin_add(row.path)
            return
        # Leaf row. If the value is None, treat it as a candidate empty
        # dict — begin_add() will promote it. Otherwise, add a sibling.
        if row.value is None:
            begin_add(row.path)
        else:
            begin_add(row.path[:-1])

    @kb.add("r", filter=_body_focused)
    def _(event):
        """Rename a dict-leaf's key in place (separate from value edit)."""
        row = current_row()
        if row is None:
            return
        if row.has_children:
            # Branch rename: same idea, walk up — rename the dict key
            # whose value is this branch. Skip if depth == 0 (top-level
            # sections shouldn't be casually renamed via this binding).
            if row.depth == 0:
                state.status = "cannot rename a top-level section"
                return
            begin_rename(row)
        else:
            begin_rename(row)

    @kb.add("d", filter=_body_focused)
    def _(event):
        """Delete the entry at the cursor. No confirmation: changes are
        only persisted on `s`, and the user can quit without saving."""
        row = current_row()
        if row is None:
            return
        if not row.path:
            state.status = "cannot delete the document root"
            return
        try:
            parent = _node_at_path(state.doc, row.path[:-1])
            key = row.path[-1]
            del parent[key]
        except (KeyError, IndexError, TypeError) as e:
            state.status = f"delete failed: {e}"
            return
        state.dirty = True
        state.status = f"deleted {'.'.join(map(str, row.path))}"
        # Cursor likely points past the end now; clamp.
        rebuild()
        if state.cursor >= len(state.rows):
            state.cursor = max(0, len(state.rows) - 1)

    @kb.add("?", filter=_body_focused)
    def _(event):
        state.show_help = not state.show_help

    @kb.add("q", filter=_body_focused)
    @kb.add("c-c")
    def _(event):
        if state.dirty:
            state.status = "unsaved changes — press Q to discard, s to save"
            return
        event.app.exit(result=0)

    @kb.add("Q", filter=_body_focused)
    def _(event):
        event.app.exit(result=0)

    # Inline-edit bindings: Enter / Escape, only while the edit textarea
    # has focus.
    _edit_focused = has_focus(edit_textarea)

    @kb.add("enter", filter=_edit_focused)
    def _(event):
        end_edit(commit=True)

    @kb.add("escape", filter=_edit_focused, eager=True)
    def _(event):
        end_edit(commit=False)

    # Search (Enter / Escape handlers live on the app-level `kb` below,
    # gated by a has_focus filter).
    search_textarea = TextArea(height=1, multiline=False, prompt="/")
    searching: List[bool] = [False]

    def begin_search() -> None:
        searching[0] = True
        search_textarea.text = ""
        layout.container = HSplit([header, body, search_textarea, footer])
        layout.focus(search_textarea)

    def end_search(commit: bool) -> None:
        if commit:
            state.search_term = search_textarea.text
            jump_to_match(forward=True)
        searching[0] = False
        layout.container = root_container
        layout.focus(body)

    def jump_to_match(forward: bool) -> None:
        if not state.search_term:
            return
        n = len(state.rows)
        if n == 0:
            return
        rng = range(state.cursor + 1, n) if forward else range(state.cursor - 1, -1, -1)
        for i in rng:
            if state.search_term.lower() in state.rows[i].label.lower():
                state.cursor = i
                return
        # Wrap around
        wrap = range(0, n) if forward else range(n - 1, -1, -1)
        for i in wrap:
            if state.search_term.lower() in state.rows[i].label.lower():
                state.cursor = i
                return
        state.status = f"no match for {state.search_term!r}"

    @kb.add("/", filter=_body_focused)
    def _(event):
        begin_search()

    @kb.add("n", filter=_body_focused)
    def _(event):
        jump_to_match(forward=True)

    @kb.add("N", filter=_body_focused)
    def _(event):
        jump_to_match(forward=False)

    _search_focused = has_focus(search_textarea)

    @kb.add("enter", filter=_search_focused)
    def _(event):
        end_search(commit=True)

    @kb.add("escape", filter=_search_focused, eager=True)
    def _(event):
        end_search(commit=False)

    app = Application(
        layout=layout,
        key_bindings=kb,
        style=style,
        full_screen=True,
        mouse_support=False,
    )
    return app.run() or 0
