#!/usr/bin/env python3
"""Extract every tool's declared output schema (streams, tables, columns) from source.

`docs/tool/*.md` documents each tool's `**Streams**` and `**Tables**` blocks, but the only
ground truth is what the tool's `get_output_files()` actually returns. This module reads that
method with a small abstract interpreter over the AST and reports, per tool, the stream names,
the table names and each table's column list. `tests/test_doc_schemas.py` compares the docs
against it, so the two cannot drift apart again.

Why an interpreter rather than a regex: `get_output_files()` builds its result dict
incrementally (`out["compounds"] = ...` under an `if`), pulls column lists from local
variables, module constants, class attributes and `self.missing_table_info()`, concatenates
them (`["id"] + provenance_columns`) and sometimes returns a different shape per `mode`. All of
that has to be followed to get the real answer, and anything that cannot be followed is
reported (`unresolved`, `exact=False`, `open_ended=True`) instead of silently guessed.

Three conventions the extractor encodes, all read off the tools that the docs already got right:

* ``DataStream.empty(name, fmt)`` is a placeholder, not an output. `Reduce` and `UniProt`
  return empty `sequences`/`compounds`/`structures` streams purely so downstream attribute
  access works, and neither documents them. A stream counts as an output only when some code
  path builds a real `DataStream(...)` for it.
* A key whose branches are an unresolvable expression and that same placeholder *is* an output,
  emitted conditionally. `OpenMM` returns ``self.ligand_stream if ... else
  DataStream.empty("compounds", "csv")``: the placeholder branch proves the key is a stream key,
  so the branch the interpreter cannot follow is a real stream on the path that takes it. Such a
  key is reported as a stream with ``always=False`` plus an ``unresolved`` note naming the
  expression; treating it as a placeholder would drop a documented output instead.
* A column whose *name* comes from a parameter (`Distance(metric_name=...)`) is unresolvable by
  definition. It is reported as ``None`` in the column list, and the docs spell it `{metric_name}`.

Run it directly for a human-readable dump, or with `--check` for the same comparison the test
makes:

    python versions/extract_output_schemas.py            # every tool's declared schema
    python versions/extract_output_schemas.py PLACER     # one tool
    python versions/extract_output_schemas.py --json     # machine-readable
    python versions/extract_output_schemas.py --check     # docs vs code, exit 1 on drift
"""

from __future__ import annotations

import ast
import io
import json
import os
import pathlib
import re
import sys
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

ROOT = pathlib.Path(__file__).resolve().parent.parent
PKG = ROOT / "biopipelines"

# Not part of the advertised public API, so not documented in docs/tool/ either. Mirrors tests/test_registry_consistency.py's own INTERNAL list.
INTERNAL = {"BoltzGenMerge", "BoltzGenImport", "RFDAA_PrepareLigand", "Mock",
            "TemplateTool", "base", "install"}

# The canonical `missing` schema, fixed once in BaseConfig.missing_table_info().
MISSING_COLUMNS = ("id", "removed_by", "kind", "cause")

# Keys of the get_output_files() dict that are neither streams nor the tables mapping. Streams are classified by value type; this only keeps the report honest about what it skipped.
NON_STREAM_KEYS = {"tables", "output_folder", "rendering_parameters"}

# How many alternatives to keep when a column list forks (mode/flag branches).
MAX_VARIANTS = 8


# --------------------------------------------------------------------------------------
# Abstract values
# --------------------------------------------------------------------------------------

class Unknown:
    """A value the interpreter could not resolve statically."""

    __slots__ = ("why",)

    def __init__(self, why: str = ""):
        self.why = why

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return f"Unknown({self.why!r})"


UNKNOWN = Unknown()


class Stream:
    """A `DataStream(...)` (real) or `DataStream.empty(...)` (placeholder) construction."""

    __slots__ = ("name", "empty")

    def __init__(self, name: Any, empty: bool = False):
        self.name = name if isinstance(name, str) else None
        self.empty = empty

    def __repr__(self) -> str:  # pragma: no cover
        return f"Stream({self.name!r}, empty={self.empty})"


class Table:
    """A `TableInfo(...)` construction (or `self.missing_table_info()`)."""

    __slots__ = ("name", "columns")

    def __init__(self, name: Any, columns: Any):
        self.name = name if isinstance(name, str) else None
        self.columns = columns

    def __repr__(self) -> str:  # pragma: no cover
        return f"Table({self.name!r}, {self.columns!r})"


class SymList:
    """A list whose items are known, plus `extra` for "and possibly more items"."""

    __slots__ = ("items", "extra")

    def __init__(self, items: Sequence[Any], extra: bool = False):
        self.items = list(items)
        self.extra = extra

    def __repr__(self) -> str:  # pragma: no cover
        return f"SymList({self.items!r}, extra={self.extra})"


class SymDict:
    """A dict being built up. `entries` maps a literal key to (value, always-present)."""

    __slots__ = ("entries", "dynamic")

    def __init__(self):
        self.entries: Dict[str, Tuple[Any, bool]] = {}
        self.dynamic: bool = False

    def set(self, key: Any, value: Any, always: bool = True) -> None:
        if isinstance(key, str):
            prev = self.entries.get(key)
            if prev is not None and not prev[1]:
                always = False
            self.entries[key] = (value, always)
        else:
            self.dynamic = True

    def clone(self) -> "SymDict":
        out = SymDict()
        out.entries = dict(self.entries)
        out.dynamic = self.dynamic
        return out

    def __repr__(self) -> str:  # pragma: no cover
        return f"SymDict({self.entries!r}, dynamic={self.dynamic})"


class Variants:
    """Two or more alternative values for one name (from an if/else or a ternary)."""

    __slots__ = ("options",)

    def __init__(self, options: Sequence[Any]):
        self.options = list(options)

    def __repr__(self) -> str:  # pragma: no cover
        return f"Variants({self.options!r})"


def _flatten(value: Any) -> List[Any]:
    """Every alternative a value can take."""
    if isinstance(value, Variants):
        out: List[Any] = []
        for opt in value.options:
            out.extend(_flatten(opt))
        return out
    return [value]


def _merge_values(*values: Any) -> Any:
    opts: List[Any] = []
    for value in values:
        for opt in _flatten(value):
            if not any(opt is kept for kept in opts):
                opts.append(opt)
    if not opts:
        return UNKNOWN
    return opts[0] if len(opts) == 1 else Variants(opts[:MAX_VARIANTS])


def _merge_dicts(a: SymDict, b: SymDict) -> SymDict:
    out = SymDict()
    out.dynamic = a.dynamic or b.dynamic
    keys = list(a.entries) + [k for k in b.entries if k not in a.entries]
    for key in keys:
        av, bv = a.entries.get(key), b.entries.get(key)
        if av is not None and bv is not None:
            out.entries[key] = (_merge_values(av[0], bv[0]), av[1] and bv[1])
        else:
            present = av if av is not None else bv
            out.entries[key] = (present[0], False)
    return out


def _deep_clone(value: Any, memo: Dict[int, Any]) -> Any:
    """Clone a value graph once per object, so aliases stay aliases in the copy.

    `result = {"tables": tables}` followed by `tables["msas"] = ...` under an `if` only
    works out if the branch's copy of `result` still points at the branch's copy of
    `tables`; cloning each name independently silently loses the second write.
    """
    if isinstance(value, (SymDict, SymList, Variants)):
        hit = memo.get(id(value))
        if hit is not None:
            return hit
    if isinstance(value, SymDict):
        out = SymDict()
        memo[id(value)] = out
        out.dynamic = value.dynamic
        out.entries = {k: (_deep_clone(v, memo), always)
                       for k, (v, always) in value.entries.items()}
        return out
    if isinstance(value, SymList):
        out_list = SymList([], value.extra)
        memo[id(value)] = out_list
        out_list.items = [_deep_clone(item, memo) for item in value.items]
        return out_list
    if isinstance(value, Variants):
        out_var = Variants([])
        memo[id(value)] = out_var
        out_var.options = [_deep_clone(option, memo) for option in value.options]
        return out_var
    return value


def _clone_env(env: Dict[str, Any]) -> Dict[str, Any]:
    memo: Dict[int, Any] = {}
    return {name: _deep_clone(value, memo) for name, value in env.items()}


_MISSING = object()


def _merge_env(base: Dict[str, Any], a: Dict[str, Any], b: Dict[str, Any]) -> Dict[str, Any]:
    out = dict(base)
    for key in set(a) | set(b):
        # A name bound to Python `None` in one branch is a real value, not an absent name.
        av = a.get(key, base.get(key, _MISSING))
        bv = b.get(key, base.get(key, _MISSING))
        if av is _MISSING and bv is _MISSING:
            continue
        if av is _MISSING:
            out[key] = bv
        elif bv is _MISSING:
            out[key] = av
        elif isinstance(av, SymDict) and isinstance(bv, SymDict):
            out[key] = _merge_dicts(av, bv)
        else:
            out[key] = _merge_values(av, bv)
    return out


# --------------------------------------------------------------------------------------
# The interpreter
# --------------------------------------------------------------------------------------

class _Return(Exception):
    """Raised to stop walking a straight-line path once it has returned."""


class SchemaInterpreter:
    """Abstractly execute one class's `get_output_files()` and collect its return dicts."""

    def __init__(self, tool: str, cls: ast.ClassDef, module_consts: Dict[str, Any],
                 self_attrs: Dict[str, Any], methods: Optional[Dict[str, ast.FunctionDef]] = None):
        self.tool = tool
        self.cls = cls
        self.consts = module_consts
        self.self_attrs = self_attrs
        self.methods: Dict[str, ast.FunctionDef] = methods if methods is not None else {
            n.name: n for n in cls.body if isinstance(n, ast.FunctionDef)
        }
        # (return value, the branch conditions that path was taken under)
        self.returns: List[Tuple[Any, Tuple[str, ...]]] = []
        self.notes: List[str] = []
        self._labels: List[str] = []
        self._depth = 0
        self._loop_depth = 0

    # -- entry point ------------------------------------------------------------------

    def run(self) -> List[Tuple[Any, Tuple[str, ...]]]:
        fn = self.methods.get("get_output_files")
        if fn is None:
            raise KeyError("get_output_files")
        try:
            self.exec_body(fn.body, {})
        except _Return:
            pass
        return self.returns

    def note(self, msg: str) -> None:
        if msg not in self.notes:
            self.notes.append(msg)

    # -- statements -------------------------------------------------------------------

    def exec_body(self, body: Sequence[ast.stmt], env: Dict[str, Any]) -> None:
        for stmt in body:
            self.exec_stmt(stmt, env)

    def exec_stmt(self, stmt: ast.stmt, env: Dict[str, Any]) -> None:
        if isinstance(stmt, (ast.Assign, ast.AnnAssign)):
            self._exec_assign(stmt, env)
        elif isinstance(stmt, ast.AugAssign):
            self._exec_augassign(stmt, env)
        elif isinstance(stmt, ast.Expr):
            self._exec_expr_stmt(stmt.value, env)
        elif isinstance(stmt, ast.If):
            self._exec_if(stmt, env)
        elif isinstance(stmt, (ast.For, ast.AsyncFor, ast.While)):
            self._exec_loop(stmt, env)
        elif isinstance(stmt, ast.Try):
            self._exec_try(stmt, env)
        elif isinstance(stmt, (ast.With, ast.AsyncWith)):
            self.exec_body(stmt.body, env)
        elif isinstance(stmt, ast.Return):
            if stmt.value is not None:
                self.returns.append((self.eval(stmt.value, env), tuple(self._labels)))
            raise _Return()
        # Everything else (raise / pass / import / assert / nested def / ...) teaches nothing about the output schema.

    def _exec_assign(self, stmt: ast.stmt, env: Dict[str, Any]) -> None:
        if isinstance(stmt, ast.AnnAssign):
            if stmt.value is None:
                return
            targets: List[ast.expr] = [stmt.target]
        else:
            targets = list(stmt.targets)  # type: ignore[union-attr]
        value = self.eval(stmt.value, env)
        for target in targets:
            self._bind(target, value, env)

    def _bind(self, target: ast.expr, value: Any, env: Dict[str, Any]) -> None:
        if isinstance(target, ast.Name):
            env[target.id] = value
        elif isinstance(target, ast.Subscript):
            container = self.eval(target.value, env)
            key = self.eval_slice(target.slice, env)
            for option in _flatten(container):
                if isinstance(option, SymDict):
                    option.set(key, value)
        elif isinstance(target, (ast.Tuple, ast.List)):
            for elt in target.elts:
                self._bind(elt, UNKNOWN, env)
        # self.X = ... inside get_output_files never changes the declared schema.

    def _exec_augassign(self, stmt: ast.AugAssign, env: Dict[str, Any]) -> None:
        if not isinstance(stmt.op, ast.Add) or not isinstance(stmt.target, ast.Name):
            self._bind(stmt.target, UNKNOWN, env)
            return
        current = env.get(stmt.target.id, UNKNOWN)
        grown = self._add(current, self.eval(stmt.value, env))
        if self._loop_depth and isinstance(grown, SymList):
            grown.extra = True
        env[stmt.target.id] = grown

    def _exec_expr_stmt(self, value: ast.expr, env: Dict[str, Any]) -> None:
        if not isinstance(value, ast.Call) or not isinstance(value.func, ast.Attribute):
            return
        method = value.func.attr
        target = self.eval(value.func.value, env)
        arg = self.eval(value.args[0], env) if value.args else UNKNOWN
        for option in _flatten(target):
            # A list grown inside a loop gains one item per iteration, and the interpreter only walks the body once, so whatever it collected is a lower bound.
            if self._loop_depth and isinstance(option, SymList) and method in (
                    "append", "add", "extend", "insert"):
                option.extra = True
            if method == "update" and isinstance(option, SymDict):
                if isinstance(arg, SymDict):
                    for key, (val, always) in arg.entries.items():
                        option.set(key, val, always)
                    option.dynamic = option.dynamic or arg.dynamic
                else:
                    option.dynamic = True
            elif method in ("append", "add") and isinstance(option, SymList):
                option.items.append(arg)
            elif method == "extend" and isinstance(option, SymList):
                if isinstance(arg, SymList):
                    option.items.extend(arg.items)
                    option.extra = option.extra or arg.extra
                else:
                    option.extra = True

    def _exec_if(self, stmt: ast.If, env: Dict[str, Any]) -> None:
        decided = self._truth(stmt.test, env)
        if decided is True:
            self.exec_body(stmt.body, env)
            return
        if decided is False:
            self.exec_body(stmt.orelse, env)
            return

        label = self._label(stmt.test, env)
        base = dict(env)
        env_true, env_false = _clone_env(env), _clone_env(env)

        returned_true = self._branch(stmt.body, env_true, label)
        returned_false = self._branch(stmt.orelse, env_false,
                                      f"not {label}" if label else None)

        if returned_true and returned_false:
            raise _Return()
        if returned_true:
            merged = env_false
        elif returned_false:
            merged = env_true
        else:
            merged = _merge_env(base, env_true, env_false)
        env.clear()
        env.update(merged)

    def _branch(self, body: Sequence[ast.stmt], env: Dict[str, Any],
                label: Optional[str]) -> bool:
        if label:
            self._labels.append(label)
        try:
            self.exec_body(body, env)
            return False
        except _Return:
            return True
        finally:
            if label:
                self._labels.pop()

    def _label(self, test: ast.expr, env: Dict[str, Any]) -> Optional[str]:
        """A short human label for a branch on a `self.<attr>` comparison."""
        if isinstance(test, ast.Compare) and len(test.ops) == 1:
            left, right = test.left, test.comparators[0]
            if (isinstance(left, ast.Attribute) and isinstance(left.value, ast.Name)
                    and left.value.id == "self" and isinstance(right, ast.Constant)):
                op = "==" if isinstance(test.ops[0], (ast.Eq, ast.Is)) else (
                    "!=" if isinstance(test.ops[0], (ast.NotEq, ast.IsNot)) else None)
                if op:
                    return f"{left.attr}{op}{right.value!r}"
        if (isinstance(test, ast.Attribute) and isinstance(test.value, ast.Name)
                and test.value.id == "self"):
            return f"{test.attr}"
        return None

    def _exec_loop(self, stmt: ast.stmt, env: Dict[str, Any]) -> None:
        target = getattr(stmt, "target", None)
        if target is not None:
            self._bind(target, UNKNOWN, env)
        # One pass: enough to see which names a loop contributes. Anything it appends to is marked open-ended, since the real iteration count is not known here.
        self._loop_depth += 1
        try:
            self.exec_body(stmt.body, env)
            self.exec_body(getattr(stmt, "orelse", None) or [], env)
        finally:
            self._loop_depth -= 1

    def _exec_try(self, stmt: ast.Try, env: Dict[str, Any]) -> None:
        try:
            self.exec_body(stmt.body, env)
        except _Return:
            self.exec_body(stmt.finalbody or [], env)
            raise
        for handler in stmt.handlers:
            try:
                self.exec_body(handler.body, dict(env))
            except _Return:
                pass
        self.exec_body(stmt.orelse or [], env)
        self.exec_body(stmt.finalbody or [], env)

    # -- expressions ------------------------------------------------------------------

    def eval_slice(self, node: ast.expr, env: Dict[str, Any]) -> Any:
        if node.__class__.__name__ == "Index":  # pragma: no cover - Python < 3.9
            node = node.value  # type: ignore[attr-defined]
        return self.eval(node, env)

    def eval(self, node: Optional[ast.expr], env: Dict[str, Any]) -> Any:
        if node is None:
            return UNKNOWN
        handler = getattr(self, "_eval_" + type(node).__name__, None)
        return UNKNOWN if handler is None else handler(node, env)

    def _truth(self, node: ast.expr, env: Dict[str, Any]) -> Optional[bool]:
        """True/False when the condition folds to a constant, else None."""
        value = self.eval(node, env)
        if isinstance(value, bool):
            return value
        return None

    def _eval_Constant(self, node: ast.Constant, env) -> Any:
        return node.value

    def _eval_Name(self, node: ast.Name, env) -> Any:
        if node.id in env:
            return env[node.id]
        if node.id in self.consts:
            return self.consts[node.id]
        return Unknown(f"name {node.id}")

    def _eval_List(self, node, env) -> Any:
        return self._eval_seq(node, env)

    _eval_Tuple = _eval_Set = _eval_List

    def _eval_seq(self, node, env) -> SymList:
        items: List[Any] = []
        extra = False
        for elt in node.elts:
            if isinstance(elt, ast.Starred):
                inner = self.eval(elt.value, env)
                if isinstance(inner, SymList):
                    items.extend(inner.items)
                    extra = extra or inner.extra
                else:
                    extra = True
            else:
                items.append(self.eval(elt, env))
        return SymList(items, extra)

    def _eval_Dict(self, node: ast.Dict, env) -> SymDict:
        out = SymDict()
        for key, value in zip(node.keys, node.values):
            if key is None:  # {**other}
                inner = self.eval(value, env)
                if isinstance(inner, SymDict):
                    for k, (v, always) in inner.entries.items():
                        out.set(k, v, always)
                    out.dynamic = out.dynamic or inner.dynamic
                else:
                    out.dynamic = True
                continue
            out.set(self.eval(key, env), self.eval(value, env))
        return out

    def _eval_JoinedStr(self, node: ast.JoinedStr, env) -> Any:
        parts = []
        for value in node.values:
            if isinstance(value, ast.Constant) and isinstance(value.value, str):
                parts.append(value.value)
            else:
                resolved = self.eval(getattr(value, "value", None), env)
                if isinstance(resolved, str):
                    parts.append(resolved)
                else:
                    return Unknown("f-string")
        return "".join(parts)

    def _eval_BinOp(self, node: ast.BinOp, env) -> Any:
        if isinstance(node.op, ast.Add):
            return self._add(self.eval(node.left, env), self.eval(node.right, env))
        return UNKNOWN

    def _add(self, left: Any, right: Any) -> Any:
        # Distribute over branches so `["a"] + ([] if x else ["b"])` keeps both shapes.
        lefts, rights = _flatten(left), _flatten(right)
        if len(lefts) > 1 or len(rights) > 1:
            return _merge_values(*[self._add(a, b) for a in lefts for b in rights])
        left, right = lefts[0], rights[0]
        if isinstance(left, str) and isinstance(right, str):
            return left + right
        if isinstance(left, SymList) or isinstance(right, SymList):
            items: List[Any] = []
            extra = False
            for side in (left, right):
                if isinstance(side, SymList):
                    items.extend(side.items)
                    extra = extra or side.extra
                else:
                    extra = True
            return SymList(items, extra)
        return UNKNOWN

    def _eval_IfExp(self, node: ast.IfExp, env) -> Any:
        decided = self._truth(node.test, env)
        if decided is True:
            return self.eval(node.body, env)
        if decided is False:
            return self.eval(node.orelse, env)
        return _merge_values(self.eval(node.body, env), self.eval(node.orelse, env))

    def _eval_BoolOp(self, node: ast.BoolOp, env) -> Any:
        return _merge_values(*[self.eval(v, env) for v in node.values])

    def _eval_UnaryOp(self, node: ast.UnaryOp, env) -> Any:
        if isinstance(node.op, ast.Not):
            inner = self.eval(node.operand, env)
            return (not inner) if isinstance(inner, bool) else UNKNOWN
        return UNKNOWN

    def _eval_Compare(self, node: ast.Compare, env) -> Any:
        if len(node.ops) != 1:
            return UNKNOWN
        left = self.eval(node.left, env)
        right = self.eval(node.comparators[0], env)
        if isinstance(left, (Unknown, SymList, SymDict, Stream, Table, Variants)):
            return UNKNOWN
        if isinstance(right, (Unknown, SymList, SymDict, Stream, Table, Variants)):
            return UNKNOWN
        op = node.ops[0]
        if isinstance(op, (ast.Eq, ast.Is)):
            return left == right
        if isinstance(op, (ast.NotEq, ast.IsNot)):
            return left != right
        return UNKNOWN

    def _eval_Attribute(self, node: ast.Attribute, env) -> Any:
        if isinstance(node.value, ast.Name) and node.value.id == "self":
            if node.attr in self.self_attrs:
                return self.self_attrs[node.attr]
            return Unknown(f"self.{node.attr}")
        return Unknown("attribute")

    def _eval_Subscript(self, node: ast.Subscript, env) -> Any:
        container = self.eval(node.value, env)
        key = self.eval_slice(node.slice, env)
        if isinstance(container, SymDict) and isinstance(key, str):
            entry = container.entries.get(key)
            return entry[0] if entry else Unknown(f"key {key}")
        if isinstance(container, SymList) and isinstance(key, int):
            try:
                return container.items[key]
            except IndexError:
                return UNKNOWN
        return UNKNOWN

    def _eval_ListComp(self, node, env) -> Any:
        return SymList([], extra=True)

    _eval_SetComp = _eval_GeneratorExp = _eval_ListComp

    def _eval_DictComp(self, node, env) -> Any:
        out = SymDict()
        out.dynamic = True
        return out

    def _eval_Call(self, node: ast.Call, env) -> Any:
        func = node.func
        kwargs = {kw.arg: kw.value for kw in node.keywords if kw.arg}

        def arg(name: str, pos: int) -> Any:
            if name in kwargs:
                return self.eval(kwargs[name], env)
            if len(node.args) > pos:
                return self.eval(node.args[pos], env)
            return UNKNOWN

        if isinstance(func, ast.Name):
            if func.id == "DataStream":
                return Stream(arg("name", 0))
            if func.id == "TableInfo":
                return Table(arg("name", 0), arg("columns", 2))
            if func.id in ("list", "tuple", "sorted", "set"):
                inner = self.eval(node.args[0], env) if node.args else UNKNOWN
                if isinstance(inner, SymList):
                    return SymList(inner.items, inner.extra)
                return SymList([], extra=True)
            if func.id == "dict":
                out = SymDict()
                if node.args:
                    inner = self.eval(node.args[0], env)
                    if isinstance(inner, SymDict):
                        for k, (v, always) in inner.entries.items():
                            out.set(k, v, always)
                    else:
                        out.dynamic = True
                for key, value in kwargs.items():
                    out.set(key, self.eval(value, env))
                return out
            return Unknown(f"{func.id}()")

        if isinstance(func, ast.Attribute):
            owner = func.value
            if isinstance(owner, ast.Name) and owner.id == "DataStream" and func.attr == "empty":
                return Stream(arg("name", 0), empty=True)
            if isinstance(owner, ast.Name) and owner.id == "self":
                return self._eval_self_call(func.attr, node, env)
            inner = self.eval(owner, env)
            if func.attr == "copy":
                if isinstance(inner, SymList):
                    return SymList(inner.items, inner.extra)
                if isinstance(inner, SymDict):
                    return inner.clone()
            if func.attr in ("get", "pop") and isinstance(inner, SymDict):
                key = self.eval(node.args[0], env) if node.args else UNKNOWN
                if isinstance(key, str):
                    entry = inner.entries.get(key)
                    if entry:
                        return entry[0]
            return Unknown(f".{func.attr}()")

        return UNKNOWN

    def _eval_self_call(self, name: str, node: ast.Call, env) -> Any:
        if name == "missing_table_info":
            return Table("missing", SymList(list(MISSING_COLUMNS)))
        fn = self.methods.get(name)
        if fn is None or self._depth >= 3:
            return Unknown(f"self.{name}()")
        # Inline the helper: only its return value matters, and only insofar as it resolves without knowing the arguments.
        sub = SchemaInterpreter(self.tool, self.cls, self.consts, self.self_attrs,
                                methods=self.methods)
        sub._depth = self._depth + 1
        sub_env: Dict[str, Any] = {a.arg: UNKNOWN for a in fn.args.args if a.arg != "self"}
        try:
            sub.exec_body(fn.body, sub_env)
        except _Return:
            pass
        for msg in sub.notes:
            self.note(msg)
        if not sub.returns:
            return Unknown(f"self.{name}() -> no static return")
        return _merge_values(*[value for value, _ in sub.returns])


# --------------------------------------------------------------------------------------
# Module / class scanning
# --------------------------------------------------------------------------------------

def _literal(value: Any) -> bool:
    return isinstance(value, (str, int, float, bool)) or isinstance(value, SymList)


def _module_constants(tree: ast.Module) -> Dict[str, Any]:
    """Module-level assignments that resolve to literals or literal lists."""
    consts: Dict[str, Any] = {}
    helper = SchemaInterpreter("<module>", ast.ClassDef(
        name="", bases=[], keywords=[], body=[], decorator_list=[]), consts, {}, methods={})
    for stmt in tree.body:
        if not isinstance(stmt, ast.Assign):
            continue
        value = helper.eval(stmt.value, {})
        if not _literal(value):
            continue
        for target in stmt.targets:
            if isinstance(target, ast.Name):
                consts[target.id] = value
    return consts


def _class_attrs(chain: Sequence[ast.ClassDef], consts: Dict[str, Any]) -> Dict[str, Any]:
    """Class-body assignments, with the most derived class in `chain` winning."""
    out: Dict[str, Any] = {}
    helper = SchemaInterpreter("<class>", chain[0], consts, {}, methods={})
    for cls in reversed(list(chain)):
        for stmt in cls.body:
            if not isinstance(stmt, ast.Assign):
                continue
            value = helper.eval(stmt.value, {})
            if not _literal(value):
                continue
            for target in stmt.targets:
                if isinstance(target, ast.Name):
                    out[target.id] = value
    return out


def _self_attrs(chain: Sequence[ast.ClassDef], consts: Dict[str, Any],
                class_attrs: Dict[str, Any]) -> Dict[str, Any]:
    """`self.X = <literal>` assignments in __init__ that resolve statically."""
    out: Dict[str, Any] = {}
    helper = SchemaInterpreter("<init>", chain[0], consts, class_attrs, methods={})
    seen: Set[str] = set()
    for cls in reversed(list(chain)):
        init = next((n for n in cls.body
                     if isinstance(n, ast.FunctionDef) and n.name == "__init__"), None)
        if init is None:
            continue
        for stmt in ast.walk(init):
            if not isinstance(stmt, ast.Assign):
                continue
            for target in stmt.targets:
                if not (isinstance(target, ast.Attribute)
                        and isinstance(target.value, ast.Name)
                        and target.value.id == "self"):
                    continue
                value = helper.eval(stmt.value, {})
                # A second, different assignment means the value depends on how the tool was configured, so it is not a static fact.
                if target.attr in seen and out.get(target.attr) != value:
                    out[target.attr] = UNKNOWN
                else:
                    out[target.attr] = value if _literal(value) else UNKNOWN
                seen.add(target.attr)
    merged = dict(class_attrs)
    merged.update({k: v for k, v in out.items() if not isinstance(v, Unknown)})
    return merged


def _tool_classes() -> Dict[str, Tuple[ast.ClassDef, ast.Module, pathlib.Path]]:
    """{TOOL_NAME: (class node, module tree, path)} for every public tool."""
    found: Dict[str, Tuple[ast.ClassDef, ast.Module, pathlib.Path]] = {}
    for path in sorted(PKG.glob("*.py")):
        src = io.open(path, encoding="utf-8", errors="replace").read()
        try:
            tree = ast.parse(src)
        except SyntaxError:  # pragma: no cover
            continue
        for node in tree.body:
            if not isinstance(node, ast.ClassDef):
                continue
            name = None
            for stmt in node.body:
                if isinstance(stmt, ast.Assign):
                    for target in stmt.targets:
                        if (isinstance(target, ast.Name) and target.id == "TOOL_NAME"
                                and isinstance(stmt.value, ast.Constant)):
                            name = stmt.value.value
            if name and name not in INTERNAL:
                found[name] = (node, tree, path)
    return found


def _base_chain(cls: ast.ClassDef, tree: ast.Module) -> List[ast.ClassDef]:
    """`cls` followed by every same-module base class, most derived first."""
    by_name = {n.name: n for n in tree.body if isinstance(n, ast.ClassDef)}
    chain: List[ast.ClassDef] = [cls]
    queue = [cls]
    while queue:
        current = queue.pop(0)
        for base in current.bases:
            name = base.id if isinstance(base, ast.Name) else (
                base.attr if isinstance(base, ast.Attribute) else None)
            parent = by_name.get(name) if name else None
            if parent is not None and parent not in chain:
                chain.append(parent)
                queue.append(parent)
    return chain


# --------------------------------------------------------------------------------------
# Column resolution
# --------------------------------------------------------------------------------------

def _column_variants(value: Any) -> Tuple[List[List[Optional[str]]], bool, List[str], List[str]]:
    """Resolve a `TableInfo(columns=...)` value.

    Returns (variants, open_ended, problems, defaults). Each variant is an ordered list of
    column names, with ``None`` where the *name* is a parameter (`Distance(metric_name=...)`)
    and so cannot be known here; `defaults` collects the literal fallbacks such a position
    can take. `open_ended` means a variant may carry further columns beyond the ones listed,
    which happens when a list is extended by a comprehension or inside a loop.
    """
    variants: List[List[Optional[str]]] = []
    open_ended = False
    problems: List[str] = []
    defaults: List[str] = []

    for option in _flatten(value):
        if isinstance(option, SymList):
            if option.extra:
                open_ended = True
                problems.append("column list is extended by a value resolved at runtime")
            expansions: List[List[Optional[str]]] = [[]]
            for item in option.items:
                choices = _flatten(item)
                literals = [c for c in choices if isinstance(c, str)]
                if len(literals) < len(choices):
                    # A configurable column name: one position, name chosen at runtime.
                    names: List[Optional[str]] = [None]
                    for literal in literals:
                        if literal not in defaults:
                            defaults.append(literal)
                else:
                    names = list(dict.fromkeys(literals))
                grown: List[List[Optional[str]]] = []
                for prefix in expansions:
                    for name in names:
                        grown.append(prefix + [name])
                expansions = grown[:MAX_VARIANTS]
            for expansion in expansions:
                if expansion not in variants:
                    variants.append(expansion)
        elif isinstance(option, Unknown):
            problems.append(f"column list unresolved ({option.why or 'dynamic'})")
        else:
            problems.append(f"column list is not a list ({type(option).__name__})")

    return variants[:MAX_VARIANTS], open_ended, problems, defaults


# --------------------------------------------------------------------------------------
# Public API
# --------------------------------------------------------------------------------------

def extract_tool_schema(tool: str, classes: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """The declared output schema of one tool.

    Keys:
      ``streams``  {name: {"always": bool}} for every real `DataStream` returned as a
                   top-level key. `DataStream.empty(...)`-only keys are excluded and listed
                   under ``placeholder_streams`` instead; a key that is a placeholder on one
                   branch and an unresolvable expression on the other counts as a stream with
                   ``always=False``.
      ``tables``   {name: {"always", "variants", "columns", "optional", "exact",
                   "open_ended", "problems"}}. ``columns`` is the ordered union across
                   variants; ``exact`` marks a single, fully-resolved, closed column list.
      ``dynamic_streams`` / ``dynamic_tables``  names are computed from the tool's inputs
      ``unresolved``  human-readable notes for anything not statically resolvable
    """
    classes = classes or _tool_classes()
    if tool not in classes:
        raise KeyError(tool)
    cls, tree, path = classes[tool]
    chain = _base_chain(cls, tree)
    owner = next((c for c in chain
                  if any(isinstance(n, ast.FunctionDef) and n.name == "get_output_files"
                         for n in c.body)), None)
    if owner is None:
        raise KeyError(f"{tool}: no get_output_files() in {path.name}")

    consts = _module_constants(tree)
    class_attrs = _class_attrs(chain, consts)
    attrs = _self_attrs(chain, consts, class_attrs)
    # Methods resolve most-derived-first so an override is what gets inlined.
    methods: Dict[str, ast.FunctionDef] = {}
    for candidate in reversed(chain):
        for node in candidate.body:
            if isinstance(node, ast.FunctionDef):
                methods[node.name] = node

    interp = SchemaInterpreter(tool, owner, consts, attrs, methods=methods)
    returns = interp.run()

    streams: Dict[str, Dict[str, Any]] = {}
    placeholders: Set[str] = set()
    tables: Dict[str, Dict[str, Any]] = {}
    unresolved: List[str] = list(interp.notes)
    dynamic_streams = False
    dynamic_tables = False
    paths: List[Dict[str, Any]] = []

    # One "shape" per (return path x alternative result dict x alternative tables dict). A name is unconditional only when every shape declares it, unconditionally.
    stream_hits: Dict[str, int] = {}
    table_hits: Dict[str, int] = {}
    n_shapes = 0

    for value, labels in returns:
        for option in _flatten(value):
            if not isinstance(option, SymDict):
                dynamic_streams = True
                unresolved.append("get_output_files() returns a value the extractor "
                                  f"cannot resolve ({type(option).__name__})")
                continue
            if option.dynamic:
                dynamic_streams = True

            shape_streams: Dict[str, bool] = {}
            for key, (val, always) in option.entries.items():
                if key == "tables" or key in NON_STREAM_KEYS:
                    continue
                options = _flatten(val)
                real = [o for o in options if isinstance(o, Stream) and not o.empty]
                empty = [o for o in options if isinstance(o, Stream) and o.empty]
                if not real and empty:
                    # An explicit placeholder on one branch proves the key is a stream key, so an unresolvable alternative beside it is a DataStream the tool really emits on that path.
                    unresolvable = [o for o in options if isinstance(o, Unknown)]
                    if unresolvable:
                        shape_streams[key] = False
                        why = ", ".join(sorted({o.why or "an unresolvable expression"
                                                for o in unresolvable}))
                        unresolved.append(
                            f"stream `{key}`: conditional — {why} on one branch and a "
                            f"DataStream.empty placeholder on the other, so it is emitted "
                            f"only when the tool's input supplies it")
                        continue
                    placeholders.add(empty[0].name or key)
                    continue
                for opt in real:
                    name = key
                    certain = always and len(options) == len(real)
                    shape_streams[name] = shape_streams.get(name, True) and certain

            table_dicts = _flatten(option.entries["tables"][0]) if "tables" in option.entries \
                else [SymDict()]
            for tval in table_dicts:
                n_shapes += 1
                for name, certain in shape_streams.items():
                    stream_hits[name] = stream_hits.get(name, 0) + 1
                    entry = streams.setdefault(name, {"always": True})
                    if not certain:
                        entry["always"] = False
                if not isinstance(tval, SymDict):
                    dynamic_tables = True
                    continue
                if tval.dynamic:
                    dynamic_tables = True
                for tname, (tinfo, talways) in tval.entries.items():
                    table_hits[tname] = table_hits.get(tname, 0) + 1
                    entry = tables.setdefault(tname, {
                        "always": True, "variants": [], "problems": [],
                        "open_ended": False, "defaults": []})
                    if not talways:
                        entry["always"] = False
                    for opt in _flatten(tinfo):
                        if not isinstance(opt, Table):
                            _add_problem(entry, "table value is not a TableInfo")
                            continue
                        variants, open_ended, problems, defaults = _column_variants(opt.columns)
                        for variant in variants:
                            if variant not in entry["variants"]:
                                entry["variants"].append(variant)
                        entry["open_ended"] = entry["open_ended"] or open_ended
                        for default in defaults:
                            if default not in entry["defaults"]:
                                entry["defaults"].append(default)
                        for problem in problems:
                            _add_problem(entry, problem)
                paths.append({"conditions": list(labels),
                              "streams": sorted(shape_streams),
                              "tables": sorted(tval.entries)})

    n_shapes = max(n_shapes, 1)
    for name, entry in streams.items():
        if stream_hits.get(name, 0) < n_shapes:
            entry["always"] = False
    for name, entry in tables.items():
        if table_hits.get(name, 0) < n_shapes:
            entry["always"] = False

    for name, entry in tables.items():
        columns: List[Optional[str]] = []
        for variant in entry["variants"]:
            for col in variant:
                if col not in columns:
                    columns.append(col)
        entry["columns"] = columns
        entry["optional"] = [c for c in columns
                             if c is not None and any(c not in v for v in entry["variants"])]
        entry["exact"] = (len(entry["variants"]) == 1 and not entry["open_ended"]
                          and not entry["problems"])
        if not entry["variants"]:
            unresolved.append(f"table `{name}`: "
                              + "; ".join(entry["problems"] or ["column list unresolved"]))
        elif entry["problems"]:
            unresolved.append(f"table `{name}`: " + "; ".join(entry["problems"]))

    if dynamic_streams:
        unresolved.append("stream names are built at runtime from the tool's inputs")
    if dynamic_tables:
        unresolved.append("table names are built at runtime from the tool's inputs")

    return {
        "tool": tool,
        "source": str(path.relative_to(ROOT)).replace(os.sep, "/"),
        "defined_in": owner.name,
        "streams": streams,
        "placeholder_streams": sorted(placeholders - set(streams)),
        "tables": tables,
        "dynamic_streams": dynamic_streams,
        "dynamic_tables": dynamic_tables,
        "unresolved": sorted(set(unresolved)),
        "paths": paths,
    }


def _add_problem(entry: Dict[str, Any], problem: str) -> None:
    if problem not in entry["problems"]:
        entry["problems"].append(problem)


def extract_all() -> Dict[str, Dict[str, Any]]:
    """{TOOL_NAME: schema} for every public tool that declares get_output_files()."""
    classes = _tool_classes()
    out: Dict[str, Dict[str, Any]] = {}
    for tool in sorted(classes):
        try:
            out[tool] = extract_tool_schema(tool, classes)
        except KeyError:
            continue
    return out


def render_columns(columns: Sequence[Optional[str]]) -> str:
    """A column list as a markdown-ready `a | b | c`, with `{...}` for runtime names."""
    return " | ".join(c if c is not None else "{...}" for c in columns)


# --------------------------------------------------------------------------------------
# Reading the docs back
# --------------------------------------------------------------------------------------

DOCS = ROOT / "docs" / "tool"

# Doc headings that are not a tool of their own.
DOC_ALIASES = {"Load / LoadMultiple": "Load"}

# Sections that deliberately carry no schema block because they document a thin variant of another tool, whose section owns the schema.
DOC_NO_SCHEMA = {"SolubleMPNN": "ProteinMPNN, whose section owns the schema"}

_HEADING = re.compile(r"^(#{2,3})\s+(.*?)\s*(?:\{#[^}]*\})?\s*$")
_LABEL = re.compile(r"^\*\*([A-Za-z][^*]*)\*\*")
_BACKTICKED = re.compile(r"`([^`]+)`")
_BULLET = re.compile(r"^\s*[-*]\s+(.*)$")
_PIPE_ROW = re.compile(r"\|(.+)\|")
_RUNTIME_CELL = re.compile(r"^\{.*\}$")


def _split_top_level(text: str) -> List[str]:
    """Split on commas/semicolons that are not inside parentheses or backticks."""
    parts, depth, tick, current = [], 0, False, []
    for char in text:
        if char == "`":
            tick = not tick
        elif not tick and char == "(":
            depth += 1
        elif not tick and char == ")":
            depth = max(0, depth - 1)
        if char in ",;" and depth == 0 and not tick:
            parts.append("".join(current))
            current = []
        else:
            current.append(char)
    parts.append("".join(current))
    return [p.strip() for p in parts if p.strip()]


def _cells(row: str) -> Tuple[List[Optional[str]], bool]:
    """The column names in one pipe-separated row, plus a "more columns follow" flag."""
    out: List[Optional[str]] = []
    open_ended = False
    for raw in row.split("|"):
        cell = raw.strip().strip("`*_:.,").strip()
        if not cell:
            continue
        if cell.startswith("...") or cell in ("…", "etc"):
            open_ended = True
            continue
        if _RUNTIME_CELL.match(cell):
            out.append(None)
            continue
        if " " in cell:  # an annotation or prose, not a column name
            open_ended = True
            continue
        out.append(cell)
    return out, open_ended


def _is_separator(row: str) -> bool:
    return not row.strip().strip("|").strip().strip("-:| ")


def _row_from_text(text: str) -> Optional[str]:
    """The pipe-separated column row inside one line of a Tables block, if any.

    Three spellings are in use across docs/tool/ and all three are accepted: the list
    inside a backtick span (`` `id | chain | resi` ``), a bare pipe row on the bullet line
    (`- \\`confidence\\`: | id | file |`), and a markdown table row on its own line.
    """
    for span in _BACKTICKED.findall(text):
        if "|" in span:
            return span
    bare = _BACKTICKED.sub(" ", text)
    if "|" in bare and not _is_separator(bare):
        return bare[bare.index("|"):bare.rindex("|") + 1]
    return None


def parse_doc_sections() -> Dict[str, Dict[str, Any]]:
    """{section title: parsed schema blocks} across every docs/tool/*.md file."""
    sections: Dict[str, Dict[str, Any]] = {}
    for path in sorted(DOCS.glob("*.md")):
        lines = io.open(path, encoding="utf-8").read().splitlines()
        title: Optional[str] = None
        start = 0
        for index, line in enumerate(lines + ["## <eof>"]):
            match = _HEADING.match(line)
            if not match:
                continue
            if title is not None:
                sections[title] = _parse_section(title, path, lines[start:index])
            title = DOC_ALIASES.get(match.group(2), match.group(2))
            start = index + 1
    return sections


def _parse_section(title: str, path: pathlib.Path, body: List[str]) -> Dict[str, Any]:
    """Pull the `**Streams**` / `**Tables**` blocks out of one doc section."""
    parsed: Dict[str, Any] = {
        "title": title,
        "file": str(path.relative_to(ROOT)).replace(os.sep, "/"),
        "has_streams_block": False,
        "has_tables_block": False,
        "streams": [],
        "prose_stream_bullets": 0,
        "tables": {},
        "table_order": [],
    }
    index, fenced = 0, False
    while index < len(body):
        line = body[index]
        if line.lstrip().startswith("```"):
            fenced = not fenced
        if fenced:
            index += 1
            continue
        label = _LABEL.match(line)
        if label and label.group(1).strip() in ("Streams", "Tables"):
            kind = label.group(1).strip()
            block, index = _collect_block(body, index)
            if kind == "Streams":
                parsed["has_streams_block"] = True
                names, prose = _parse_streams_block(block)
                for name in names:
                    if name not in parsed["streams"]:
                        parsed["streams"].append(name)
                parsed["prose_stream_bullets"] += prose
            else:
                parsed["has_tables_block"] = True
                for name, rows, open_ended in _parse_tables_block(block):
                    entry = parsed["tables"].setdefault(
                        name, {"rows": [], "open_ended": False})
                    entry["rows"].extend(rows)
                    entry["open_ended"] = entry["open_ended"] or open_ended
                    if name not in parsed["table_order"]:
                        parsed["table_order"].append(name)
            continue
        index += 1
    return parsed


def _collect_block(body: List[str], index: int) -> Tuple[List[str], int]:
    """The lines of one `**Label**` block: from the label to the next label / rule / fence."""
    block = [body[index]]
    cursor = index + 1
    while cursor < len(body):
        line = body[cursor]
        if line.lstrip().startswith("```") or line.strip() == "---":
            break
        if _LABEL.match(line):
            break
        block.append(line)
        cursor += 1
    return block, cursor


def _parse_streams_block(block: List[str]) -> Tuple[List[str], int]:
    names: List[str] = []
    prose = 0
    head = block[0]
    remainder = head.split("**", 2)[-1].lstrip(":").strip()
    for part in _split_top_level(remainder):
        found = _BACKTICKED.search(part)
        if found and part.startswith("`"):
            names.append(found.group(1))
    for line in block[1:]:
        bullet = _BULLET.match(line)
        if not bullet:
            continue
        text = bullet.group(1).strip()
        if text.startswith("`"):
            names.append(_BACKTICKED.match(text).group(1))
        else:
            prose += 1
    return names, prose


def _parse_tables_block(block: List[str]) -> List[Tuple[str, List[List[Optional[str]]], bool]]:
    """Every `- \\`name\\`` bullet in a Tables block, with the pipe rows under it."""
    found: List[Tuple[str, List[List[Optional[str]]], bool]] = []
    current: Optional[str] = None
    rows: List[List[Optional[str]]] = []
    open_ended = False

    def flush() -> None:
        if current is not None:
            found.append((current, rows, open_ended))

    def take(text: str) -> None:
        nonlocal open_ended
        row = _row_from_text(text)
        if row is None:
            return
        cells, more = _cells(row)
        if len(cells) >= 2 and cells not in rows:
            rows.append(cells)
            open_ended = open_ended or more

    for line in block:
        bullet = _BULLET.match(line)
        indent = len(line) - len(line.lstrip())
        if bullet and indent <= 1 and bullet.group(1).strip().startswith("`"):
            flush()
            text = bullet.group(1).strip()
            current = _BACKTICKED.match(text).group(1)
            rows, open_ended = [], False
            take(text[text.index("`", 1) + 1:])
            continue
        if current is None:
            continue
        take(line)
    flush()
    return found


# --------------------------------------------------------------------------------------
# Comparing docs against code
# --------------------------------------------------------------------------------------

def compare(tool: str, schema: Dict[str, Any], doc: Optional[Dict[str, Any]]) -> List[str]:
    """Every way `doc` disagrees with `schema`, as complete human-readable messages."""
    problems: List[str] = []
    where = f"{doc['file']} section `{tool}`" if doc else "docs/tool/"

    if doc is None:
        if schema["streams"] or schema["tables"]:
            problems.append(f"{tool}: no docs/tool/*.md section documents this tool")
        return problems

    code_streams = list(schema["streams"])
    doc_streams = list(doc["streams"])

    if schema["dynamic_streams"]:
        # The tool mirrors its inputs' stream names, so the docs can only describe the rule in prose. Still require every statically-declared stream to be listed.
        for name in code_streams:
            if name not in doc_streams:
                problems.append(
                    f"{tool}: {where} **Streams** omits `{name}`, which "
                    f"{schema['source']} declares (this tool also builds stream names at "
                    f"runtime, so extra prose entries are allowed)")
    else:
        if code_streams and not doc["has_streams_block"]:
            problems.append(f"{tool}: {where} has no **Streams** block, but "
                            f"{schema['source']} returns "
                            + ", ".join(f"`{n}`" for n in code_streams))
        elif set(doc_streams) != set(code_streams):
            missing = [n for n in code_streams if n not in doc_streams]
            extra = [n for n in doc_streams if n not in code_streams]
            detail = []
            if missing:
                detail.append("missing " + ", ".join(f"`{n}`" for n in missing))
            if extra:
                placeholder = [n for n in extra if n in schema["placeholder_streams"]]
                note = (" (a DataStream.empty placeholder, not an output)"
                        if placeholder == extra and extra else "")
                detail.append("documents nonexistent " + ", ".join(f"`{n}`" for n in extra)
                              + note)
            problems.append(f"{tool}: {where} **Streams** " + "; ".join(detail)
                            + f" — code declares [{', '.join(code_streams) or 'none'}]")

    code_tables = list(schema["tables"])
    doc_tables = list(doc["tables"])

    if code_tables and not doc["has_tables_block"]:
        problems.append(f"{tool}: {where} has no **Tables** block, but "
                        f"{schema['source']} returns "
                        + ", ".join(f"`{n}`" for n in code_tables))
        return problems

    for name in code_tables:
        if name not in doc_tables:
            entry = schema["tables"][name]
            problems.append(
                f"{tool}: {where} **Tables** omits `{name}` "
                f"({'conditional; ' if not entry['always'] else ''}"
                f"columns: {render_columns(entry['columns'])})")
    if not schema["dynamic_tables"]:
        for name in doc_tables:
            if name not in code_tables:
                problems.append(f"{tool}: {where} documents table `{name}`, which "
                                f"{schema['source']} never declares")

    for name in code_tables:
        if name not in doc_tables:
            continue
        entry = schema["tables"][name]
        if not entry["variants"]:
            # The column list is computed from the tool's inputs (Panda.result, a Table(...) the user defines): there is nothing static to check against.
            continue
        rows = doc["tables"][name]["rows"]
        doc_open = doc["tables"][name]["open_ended"]
        if not rows:
            if not [c for c in entry["columns"] if c is not None]:
                # Every column name is chosen at runtime, so there is nothing to spell out.
                continue
            problems.append(f"{tool}: {where} table `{name}` documents no column row; "
                            f"expected: {render_columns(entry['columns'])}")
            continue
        if entry["exact"] and len(rows) == 1:
            if rows[0] != entry["columns"]:
                problems.append(
                    f"{tool}: {where} table `{name}` columns disagree\n"
                    f"    docs: {render_columns(rows[0])}\n"
                    f"    code: {render_columns(entry['columns'])}")
            continue
        doc_union = [c for row in rows for c in row]
        code_union = entry["columns"]
        if entry["open_ended"]:
            # The extra columns are named at runtime, so the docs stand in for them with a concrete example instead of a placeholder cell.
            code_union = [c for c in code_union if c is not None]
        missing = [c for c in code_union if c not in doc_union]
        extra = [c for c in doc_union if c not in code_union]
        if missing:
            problems.append(
                f"{tool}: {where} table `{name}` omits column(s) "
                + ", ".join(f"`{c}`" if c else "`{...}`" for c in missing)
                + f"\n    docs: {render_columns(doc_union)}\n"
                f"    code: {render_columns(code_union)}")
        if extra and not entry["open_ended"] and not doc_open:
            problems.append(
                f"{tool}: {where} table `{name}` documents column(s) the code never "
                "declares: " + ", ".join(f"`{c}`" if c else "`{...}`" for c in extra)
                + f"\n    code: {render_columns(code_union)}")
    return problems


def compare_all() -> Tuple[Dict[str, List[str]], List[str]]:
    """({tool: problems}, orphan doc sections) across the whole repo."""
    schemas = extract_all()
    docs = parse_doc_sections()
    problems: Dict[str, List[str]] = {}
    for tool, schema in schemas.items():
        found = compare(tool, schema, docs.get(tool))
        if found:
            problems[tool] = found
    orphans = [title for title, doc in docs.items()
               if title not in schemas and title not in DOC_NO_SCHEMA
               and (doc["has_streams_block"] or doc["has_tables_block"])]
    return problems, orphans


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------

def _format(schema: Dict[str, Any]) -> str:
    lines = [f"### {schema['tool']}  ({schema['source']}::{schema['defined_in']})"]
    if schema["streams"]:
        lines.append("  streams: " + ", ".join(
            f"`{n}`" + ("" if e["always"] else " (conditional)")
            for n, e in schema["streams"].items()))
    else:
        lines.append("  streams: (none)")
    if schema["placeholder_streams"]:
        lines.append("  placeholders (not documented): "
                     + ", ".join(f"`{n}`" for n in schema["placeholder_streams"]))
    if schema["tables"]:
        for name, entry in schema["tables"].items():
            flags = []
            if not entry["always"]:
                flags.append("conditional")
            if entry["open_ended"]:
                flags.append("open-ended")
            if not entry["exact"]:
                flags.append("inexact")
            suffix = f" [{', '.join(flags)}]" if flags else ""
            lines.append(f"  table `{name}`{suffix}: {render_columns(entry['columns'])}")
            if len(entry["variants"]) > 1:
                for variant in entry["variants"]:
                    lines.append("      variant: " + render_columns(variant))
    else:
        lines.append("  tables: (none)")
    for path in schema["paths"]:
        if path["conditions"]:
            lines.append("  path " + " and ".join(path["conditions"])
                         + f": streams={path['streams']} tables={path['tables']}")
    for note in schema["unresolved"]:
        lines.append("  !! " + note)
    return "\n".join(lines)


def main(argv: Sequence[str]) -> int:
    args = [a for a in argv if not a.startswith("--")]
    schemas = extract_all()
    if args:
        unknown = [a for a in args if a not in schemas]
        if unknown:
            print("unknown tool(s): " + ", ".join(unknown), file=sys.stderr)
            return 2
        schemas = {k: v for k, v in schemas.items() if k in args}
    if "--check" in argv:
        problems, orphans = compare_all()
        total = 0
        for tool in sorted(problems):
            if args and tool not in args:
                continue
            print(f"### {tool}")
            for problem in problems[tool]:
                print("  - " + problem)
                total += 1
        for orphan in orphans:
            print(f"### <orphan doc section> {orphan}")
            total += 1
        correct = len(schemas) - len([t for t in problems if not args or t in args])
        print(f"\n{total} problem(s) across {len(problems)} tool(s); "
              f"{correct} tool(s) already correct")
        return 1 if total else 0
    if "--json" in argv:
        print(json.dumps(schemas, indent=2, sort_keys=True, default=str))
    else:
        for schema in schemas.values():
            print(_format(schema))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
