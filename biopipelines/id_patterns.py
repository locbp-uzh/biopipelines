# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Pattern-based ID and file path utilities for lazy DataStream expansion.

Provides pattern syntax for compact representation of ID lists:
    <0..49>       Numeric range (inclusive): 50 IDs
    <A B C>       Explicit set (space-separated): 3 IDs
    <0..2>_<A B>  Multi-slot: Cartesian product → 6 IDs

Bracket syntax for runtime-dependent (lazy) patterns:
    [_<?>]        Bracket segment — cannot be expanded at config time
    [_<#><A V>]   A slot may still constrain the shape the value will take
    prot_<0..4>[_<?>]  Mixed: deterministic prefix + lazy suffix

All functions are standalone (no DataStream imports).
"""

import re
from functools import lru_cache
from itertools import product
from typing import List, Optional, Tuple

# Optional so this module stays standalone-importable, which some pipe scripts rely on.
try:
    from . import contract_enforcement as _contract_enforcement
except ImportError:
    try:
        import contract_enforcement as _contract_enforcement
    except ImportError:
        _contract_enforcement = None

# Regex to find pattern slots: <...>, wherever they sit — a slot inside a [...] bracket matches too, which is what makes such an id a pattern that cannot yet be expanded.
_SLOT_RE = re.compile(r'<([^>\[\]]+)>')

# Regex to find bracket segments: [...]
_BRACKET_RE = re.compile(r'\[([^\]]+)\]')

# Separates a parent id from a child suffix, so no single suffix value contains it.
_SUFFIX_DELIMITER = '_'

# Public spelling of the delimiter and of "one segment", for id_map_utils, which needs the same rule
# to strip suffixes while its own slot notation stays deliberately different. See SEGMENT_REGEX.
SUFFIX_DELIMITER = _SUFFIX_DELIMITER
SEGMENT_REGEX = f'[^{re.escape(_SUFFIX_DELIMITER)}]+'
DIGIT_REGEX = r'\d+'


# `<...>` holds literal values and compacts a *predictable* id; `[...]` marks an *unpredictable* one.
# A class matches without enumerating, so it is only meaningful where matching happens: inside a bracket.
DIGIT_CLASS = '<#>'
SEGMENT_CLASS = '<?>'
MATCH_CLASSES = (DIGIT_CLASS, SEGMENT_CLASS)


class LazyPatternError(Exception):
    """Raised when trying to fully expand a pattern that contains brackets."""
    pass


def is_lazy(s: str) -> bool:
    """True if string contains a [...] bracket segment (runtime-dependent)."""
    return bool(_BRACKET_RE.search(s))


def has_brackets(s: str) -> bool:
    """Alias for is_lazy — True if string contains [...] brackets."""
    return is_lazy(s)


def strip_brackets(s: str) -> str:
    """Remove all [...] bracket segments from a pattern string.

    Returns the deterministic prefix (the part that can be expanded at config time).

    'prot_<0..4>[_<?>]'   → 'prot_<0..4>'
    'base[_<#><A V>]'     → 'base'
    'literal'             → 'literal'
    """
    return _BRACKET_RE.sub('', s)


def partial_expand(s: str) -> List[str]:
    """Expand the deterministic, out-of-bracket ``<..>`` slots, keep ``[...]`` verbatim.

    A pattern is just an id with unresolved slots; ``[...]`` segments stay
    unresolved (their values live in a runtime map_table). Expanding only the
    deterministic slots keeps the lazy bracket intact, so each result is still a
    valid lazy pattern, not a fabricated concrete id.

    'prot_<0..2>'                → ['prot_0', 'prot_1', 'prot_2']
    'a_<1 2>[_<N>]'              → ['a_1[_<N>]', 'a_2[_<N>]']
    'a_<1 2>[_<N>]_5<A S>'       → ['a_1[_<N>]_5A', 'a_1[_<N>]_5S',
                                    'a_2[_<N>]_5A', 'a_2[_<N>]_5S']
    'base[_<N>]'                 → ['base[_<N>]']
    """
    bracket_spans = [(m.start(), m.end()) for m in _BRACKET_RE.finditer(s)]
    slots = [m for m in _find_slots(s)
             if not any(b0 <= m.start() < b1 for b0, b1 in bracket_spans)]
    if not slots:
        return [s]
    slot_values = [_parse_slot(m.group(1)) for m in slots]
    results = []
    for combo in product(*slot_values):
        result = s
        for m, val in zip(reversed(slots), reversed(combo)):
            result = result[:m.start()] + val + result[m.end():]
        results.append(result)
    return results


def partial_expand_ids(ids: List[str]) -> List[str]:
    """Apply :func:`partial_expand` to each pattern and concatenate."""
    result = []
    for s in ids:
        result.extend(partial_expand(s))
    return result


def try_expand(s: str) -> Tuple[List[str], bool]:
    """Partially expand a pattern, reporting whether the result is concrete.

    Returns ``(ids, is_complete)``: ``is_complete`` is False when ``[...]``
    brackets remain (the ids are still lazy patterns, not concrete). Brackets
    are preserved, not stripped.

    'prot_<0..2>'                → (['prot_0', 'prot_1', 'prot_2'], True)
    'a_<1 2>[_<N>]'              → (['a_1[_<N>]', 'a_2[_<N>]'], False)
    """
    return partial_expand(s), not is_lazy(s)


def glob_from_lazy(s: str) -> str:
    """Replace [...] bracket segments with '*' to produce a glob pattern.

    Deliberately loose: a bracket is an optional, repeatable suffix group, which no glob can express, so this only *enumerates* candidates. Callers must then narrow with :func:`select_ids`, which is the authority on what a pattern covers — see ``pipe_check_completion._resolve_lazy_pairs``.

    'prot_<0..2>[_<N><A I L V>]+9DP'  → 'prot_<0..2>*+9DP'
    'prot_<0..2>[_<N><A I L V>]'       → 'prot_<0..2>*'
    'literal'                           → 'literal'
    """
    return _BRACKET_RE.sub('*', s)


def glob_from_lazy_ids(ids: List[str]) -> List[str]:
    """Expand deterministic slots and insert '*' where brackets were.

    For each ID pattern, replaces [...] with '*' then expands <..> slots.
    Returns glob-ready strings suitable for file matching.

    ['prot_<0..1>[_<N><A V>]+X']  → ['prot_0*+X', 'prot_1*+X']
    ['prot_<0..1>']               → ['prot_0', 'prot_1']
    """
    result = []
    for s in ids:
        globbed = glob_from_lazy(s)          # [...] → *
        result.extend(expand_pattern(globbed))  # expand <..> slots (no brackets left)
    return result


def _translate_literal(text: str) -> str:
    """Regex source for a deterministic chunk: matched exactly, apart from glob ``*``/``?`` that callers may still pass."""
    out = []
    for ch in text:
        if ch == '*':
            out.append('.*')
        elif ch == '?':
            out.append('.')
        else:
            out.append(re.escape(ch))
    return ''.join(out)


def _translate_slot_class(slot_content: str) -> str:
    """Regex for one slot inside a bracket: the shape it names, or a segment when it names none.

    A bracket defers to runtime, so a slot there cannot be enumerated -- but it can still say what shape the value will take, and saying so is what stops `design_1[_<#>]` from covering `design_1_ZZZ`.

    A slot holds a literal value set everywhere, bracket or not, so `<N>` is the one-element set {'N'} here too and `parent_[<N>]` repeats that letter. Only `<#>` and `<?>` are shapes. Write `[<?>]` for "any segment"; a bare word no longer means that.
    """
    slot_content = slot_content.strip()
    if f'<{slot_content}>' == DIGIT_CLASS:
        return DIGIT_REGEX
    if f'<{slot_content}>' == SEGMENT_CLASS:
        return SEGMENT_REGEX
    try:
        values = _parse_slot(slot_content)
    except ValueError:
        return SEGMENT_REGEX
    if len(values) == 1:
        return _translate_literal(values[0])
    # A single-character alternation is a character class, which the engine matches without branching.
    if all(len(v) == 1 for v in values):
        return '[' + ''.join(re.escape(v) for v in sorted(set(values))) + ']'
    return '(?:' + '|'.join(re.escape(v) for v in values) + ')'

# `\d+` and `[^_]+` both match arbitrarily long runs, and `\d` is a subset of `[^_]`.
_UNBOUNDED_CLASSES = (SEGMENT_REGEX, DIGIT_REGEX)


def _collapse_unbounded(parts: List[str]) -> str:
    r"""Linear regex for a bracket body that is nothing but unbounded classes.

    ``(?:\d+[^_]+)*`` has two unbounded quantifiers competing for the same
    characters -- ``<?>`` is ``[^_]+``, which includes digits -- so rejecting an
    id costs one attempt per ordered partition of its run: 395 ms at 22
    characters, doubling with each one. Committing to the first class would be
    linear but changes the language, so the shape is rewritten instead, to one
    that matches exactly the same strings.

    Two rewrites, both about who can absorb whose characters. A piece with a
    strictly wider neighbour gives up its quantifier, because that neighbour
    can take whatever it left: ``\d+[^_]+`` is ``\d[^_]+``. Within a run of
    equal-width pieces one quantifier covers the run: ``\d+\d+`` is
    ``\d+\d``, still "two or more digits", which no single class says. The
    outer ``*`` becomes ``?`` because two iterations are then always subsumed
    by one.

    ``[<?><#><?>]`` is why the rule is per-neighbour rather than "keep the
    widest": both ``[^_]+`` must keep their quantifiers, and only the digit
    between them gives one up.
    """
    if len(parts) == 1:
        return parts[0][:-1] + '*'
    width = {DIGIT_REGEX: 0, SEGMENT_REGEX: 1}
    keep = [True] * len(parts)
    for i, piece in enumerate(parts):
        neighbours = [parts[j] for j in (i - 1, i + 1) if 0 <= j < len(parts)]
        if any(width[n] > width[piece] for n in neighbours):
            keep[i] = False
    for i in range(len(parts) - 1):
        if keep[i] and keep[i + 1] and parts[i] == parts[i + 1]:
            keep[i] = False
    body = ''.join(p if keep[i] else p[:-1] for i, p in enumerate(parts))
    return '(?:' + body + ')?'


def _translate_bracket(content: str) -> str:
    """Regex source for one ``[...]`` segment: a zero-or-more repetition of its content, with each slot value bounded by the suffix delimiter.

    Two rules, and where the delimiter sits does the rest — no special cases:

    ``parent[_<N>]`` puts the delimiter *inside* the group, so each repetition carries its own and the group is genuinely optional and repeatable: ``parent``, ``parent_4E``, ``parent_4E_6U``.

    ``parent_[<N>]`` puts it *outside*, so the value cannot cross a delimiter and the group collapses to at most one suffix: ``parent_4E`` yes, ``parent_4E_6U`` no. Zero repetitions would give the bare ``parent_``, which is never a real id — so in practice this spelling means exactly one, which is what a flat renumber (``name_1``, ``name_2``) wants.

    Bounding the value is also what keeps selection exact: ``design_1[_<N>]`` cannot reach into ``design_10_5``.
    """
    slots = list(_SLOT_RE.finditer(content))
    if not slots:
        return '.*'  # no slot to bound, so keep the whole-segment wildcard
    parts: List[str] = []
    pos = 0
    for m in slots:
        lit = content[pos:m.start()]
        if lit:
            parts.append(_translate_literal(lit))
        parts.append(_translate_slot_class(m.group(1)))
        pos = m.end()
    lit = content[pos:]
    if lit:
        parts.append(_translate_literal(lit))
    if all(p in _UNBOUNDED_CLASSES for p in parts):
        return _collapse_unbounded(parts)
    return '(?:' + ''.join(parts) + ')*'


@lru_cache(maxsize=4096)
def _regex_from_lazy(s: str) -> re.Pattern:
    """Compile one bracket-bearing pattern into a regex for :func:`re.Pattern.fullmatch`.

    Everything outside ``[...]`` is literal (deterministic ``<..>`` slots are expected to be expanded already, by :func:`partial_expand`). Each bracket becomes an optional repeatable group — see :func:`_translate_bracket` — so ``design_1[_<?>]`` covers ``design_1``, ``design_1_5`` and ``design_1_5_7``, but never reaches into ``design_10_5``.

    A bracket that starts with a slot (``'design_1[<?>]'``) has no separator to delimit repetitions, so it collapses to a single unbounded run rather than a repetition.
    """
    parts = []
    pos = 0
    for m in _BRACKET_RE.finditer(s):
        parts.append(_translate_literal(s[pos:m.start()]))
        parts.append(_translate_bracket(m.group(1)))
        pos = m.end()
    parts.append(_translate_literal(s[pos:]))
    return re.compile(''.join(parts), re.DOTALL)


def select_ids(patterns: List[str], row_ids: List[str], where: str = "") -> List[str]:
    """Select the row ids a pattern set covers, in row order.

    The row ids (read from a runtime map_table) are the source of truth; the patterns only *select* among them. Deterministic ``<..>`` slots and literals match exactly, anchored at both ends; a ``[...]`` bracket keeps its literal text and wildcards only its ``<..>`` slots. Ids the patterns do not cover are dropped; patterns never fabricate ids absent from the rows, and never pull in ids they do not cover. This honors upstream filtering automatically — a filtered stream simply has fewer rows.

    Anchoring matters at ten or more ids: ``'design_<1..2>[_<N>]'`` selects ``design_1_*`` and ``design_2_*`` only, never ``design_10_*``.

    patterns=['a_<1 2>[_<N>]'], rows=['a_1_x','a_2_y','b_1'] → ['a_1_x','a_2_y']
    patterns=['a_<0..2>'],      rows=['a_0','a_2']           → ['a_0','a_2']

    Selecting none of a non-empty row set is reported via ``contract_enforcement.check_pattern_selection`` — an empty result is otherwise indistinguishable from having had nothing to select. ``where`` names the caller in that report.
    """
    regexes = [_regex_from_lazy(p) for p in partial_expand_ids(patterns)]
    selected = [rid for rid in row_ids if any(r.fullmatch(rid) for r in regexes)]
    if _contract_enforcement is not None:
        _contract_enforcement.report(
            _contract_enforcement.check_pattern_selection(patterns, row_ids, selected, where)
        )
    return selected


def resolve_pattern_ids(patterns: List[str], map_table: str) -> List[str]:
    """Runtime: select the ids of one map_table that a pattern set covers.

    Reads ``map_table``'s ``id`` column as the authoritative row set and returns
    :func:`select_ids` over it. Multi-source consumers call this once per source
    and concatenate.
    """
    import pandas as pd
    df = pd.read_csv(map_table, dtype={'id': str})
    if 'id' not in df.columns:
        raise KeyError(f"map_table {map_table} has no 'id' column")
    return select_ids(patterns, [str(v) for v in df['id'].tolist()], where=map_table)


def contains_pattern(s: str) -> bool:
    """True if the string carries a ``<..>`` slot, in or out of a ``[...]`` bracket.

    This answers *"is this a pattern?"*, not *"can I expand it now?"* — the two questions have different answers for a lazy id, and conflating them is how a caller ends up handing :func:`expand_ids` something that raises. A caller guarding a full expansion wants :func:`can_expand`; a caller asking whether an ids list is still symbolic (so a 1:1 files/ids length check does not apply) wants :func:`is_pattern`.
    """
    return bool(_SLOT_RE.search(s))


def is_pattern(s: str) -> bool:
    """True if the string is a pattern rather than a concrete id.

    Either kind of unresolved part counts: a ``<..>`` slot, or a ``[...]`` bracket whose content has no slot at all (``'[_x]'``).
    """
    return contains_pattern(s) or is_lazy(s)


def can_expand(s: str) -> bool:
    """True if :func:`expand_pattern` would turn this string into concrete ids now.

    A pattern with a ``[...]`` bracket cannot: its values live in a runtime map_table, so expansion is deferred to :func:`select_ids`, and :func:`expand_pattern` raises :class:`LazyPatternError`.
    """
    return contains_pattern(s) and not is_lazy(s)


def can_expand_ids(ids: List[str]) -> bool:
    """The list-level guard for :func:`expand_ids`: something to expand, and nothing lazy.

    A mixed list (one deterministic pattern, one lazy) answers False — expanding it would raise on the lazy element. Use :func:`partial_expand_ids` there, which expands what it can and leaves brackets intact.
    """
    return any(contains_pattern(s) for s in ids) and not any(is_lazy(s) for s in ids)


def is_literal(s: str) -> bool:
    """True if the string is a concrete id with nothing left to resolve."""
    return not is_pattern(s)


def _parse_slot(slot_content: str) -> List[str]:
    """Parse a single slot's content into its values.

    '0..49'  → ['0', '1', ..., '49']
    'A B C'  → ['A', 'B', 'C']

    A slot outside a bracket is a finite set of literals, so '<N>' is the one-element set {'N'} -- the letter, not a number. Write '<#>' for a number, and only inside a bracket, where nothing is enumerated.
    """
    slot_content = slot_content.strip()
    if f'<{slot_content}>' in MATCH_CLASSES:
        raise ValueError(
            f"'<{slot_content}>' is a match class and cannot be expanded, because a class has no "
            f"finite set of values. Use it inside a [...] bracket, which defers to runtime; outside "
            f"one, write the values you mean, e.g. '<0..9>'."
        )
    if '..' in slot_content:
        parts = slot_content.split('..')
        if len(parts) != 2:
            raise ValueError(f"Invalid range syntax: '{slot_content}'")
        start, end = int(parts[0]), int(parts[1])
        return [str(i) for i in range(start, end + 1)]
    else:
        values = slot_content.split()
        if not values:
            raise ValueError(f"Empty pattern slot: '<{slot_content}>'")
        return values


def _find_slots(s: str) -> List[re.Match]:
    """Find all pattern slots in a string."""
    return list(_SLOT_RE.finditer(s))


def count_pattern(s: str) -> int:
    """Count how many IDs a pattern expands to (without expanding).

    For lazy patterns (with brackets), counts only the deterministic prefix.

    'base_<0..49>'               → 50
    '<0..2>_<A B>'               → 6
    'literal'                    → 1
    'prot_<0..4>[_<N><A V>]'     → 5 (prefix only)
    """
    if is_lazy(s):
        s = strip_brackets(s)
    slots = _find_slots(s)
    if not slots:
        return 1
    total = 1
    for m in slots:
        total *= len(_parse_slot(m.group(1)))
    return total


def count_ids(ids: List[str]) -> int:
    """Sum of counts across all list elements."""
    return sum(count_pattern(s) for s in ids)


def expand_pattern(s: str) -> List[str]:
    """Expand a pattern string into all its IDs.

    'base_<0..2>'       → ['base_0', 'base_1', 'base_2']
    '<0..1>_<A B>'      → ['0_A', '0_B', '1_A', '1_B']
    'literal'           → ['literal']

    Raises LazyPatternError if the string contains [...] brackets.
    """
    if is_lazy(s):
        raise LazyPatternError(
            f"Cannot fully expand lazy pattern '{s}': contains bracket segments. "
            f"Use try_expand() for partial expansion or expand at runtime."
        )
    slots = _find_slots(s)
    if not slots:
        return [s]

    # Parse all slot values
    slot_values = [_parse_slot(m.group(1)) for m in slots]

    # Build results via cartesian product
    results = []
    for combo in product(*slot_values):
        result = s
        # Replace slots in reverse order to preserve positions
        for m, val in zip(reversed(slots), reversed(combo)):
            result = result[:m.start()] + val + result[m.end():]
        results.append(result)
    return results


def expand_ids(ids: List[str]) -> List[str]:
    """Expand each element in the list and concatenate results.

    Raises LazyPatternError if any element contains brackets.
    """
    result = []
    for s in ids:
        result.extend(expand_pattern(s))
    return result


def try_expand_ids(ids: List[str]) -> Tuple[List[str], bool]:
    """Partially expand all IDs, preserving ``[...]`` brackets.

    Returns:
        (expanded_ids, is_complete) — is_complete is False if any ID is still lazy.
    """
    result = partial_expand_ids(ids)
    complete = not any(is_lazy(s) for s in ids)
    return result, complete


def dedup_parent_children(ids: List[str]) -> List[str]:
    """Remove literal IDs that are already covered by a pattern in the list.

    ['prot_<0..2>', 'prot_0', 'prot_1'] -> ['prot_<0..2>']
    ['prot_0', 'prot_1', 'other']       -> ['prot_0', 'prot_1', 'other']
    """
    patterns = [s for s in ids if contains_pattern(s)]
    if not patterns:
        return ids

    covered = set()
    for p in patterns:
        try:
            covered.update(expand_pattern(p))
        except LazyPatternError:
            expanded, _ = try_expand(p)
            covered.update(expanded)

    return [s for s in ids if contains_pattern(s) or s not in covered]


def expand_at(s: str, index: int) -> str:
    """Get a single expanded element by index without full expansion.

    For single-slot patterns, this is O(1). For multi-slot, it computes
    the cartesian index decomposition.
    """
    slots = _find_slots(s)
    if not slots:
        if index != 0:
            raise IndexError(f"Index {index} out of range for literal '{s}'")
        return s

    slot_values = [_parse_slot(m.group(1)) for m in slots]
    slot_sizes = [len(v) for v in slot_values]

    # Decompose flat index into per-slot indices (row-major order)
    total = 1
    for sz in slot_sizes:
        total *= sz
    if index < 0 or index >= total:
        raise IndexError(f"Index {index} out of range for pattern '{s}' ({total} items)")

    indices = []
    remaining = index
    for sz in slot_sizes:
        total //= sz
        slot_idx = remaining // total
        remaining %= total
        indices.append(slot_idx)

    # Substitute
    result = s
    for m, slot_vals, si in zip(reversed(slots), reversed(slot_values), reversed(indices)):
        result = result[:m.start()] + slot_vals[si] + result[m.end():]
    return result


def expand_file_pattern(template: str, item_id: str) -> str:
    """Substitute <id> in a file template with an actual ID.

    '<id>.pdb'  + '5HG6_0'  → '5HG6_0.pdb'
    """
    return template.replace('<id>', item_id)


def file_has_glob(template: str) -> bool:
    """True if the file template contains a glob wildcard '*'."""
    return '*' in template


def glob_from_file_pattern(template: str) -> str:
    """Turn a '<id>' file template into a glob matching every id's file.

    '<id>_best.pdb' → '*_best.pdb'
    """
    return template.replace('<id>', '*')


def _basename(path: str) -> str:
    """Last path component, splitting on either separator.

    A path recorded on a cluster is read back on Windows and vice versa, so
    os.path.basename would keep the foreign separator inside the name.
    """
    return path.replace('\\', '/').rsplit('/', 1)[-1]


def id_from_file_pattern(template: str, path: str) -> Optional[str]:
    """Recover the id a path was built from, inverting :func:`expand_file_pattern`.

    '<id>_best.pdb' + 'design_1_best.pdb' → 'design_1'

    Only the file name is inverted; the template's directory is ignored, since a
    caller may hold the same file reached by a different route. Returns None when
    the template carries no '<id>' or the name does not fit it, which is also how
    a file that is not this stream's at all gets rejected.
    """
    tmpl = _basename(template)
    if '<id>' not in tmpl:
        return None
    parts = tmpl.split('<id>')
    source = _translate_literal(parts[0])
    for i, part in enumerate(parts[1:], start=1):
        source += r'(?P<id>.+)' if i == 1 else r'(?P=id)'
        source += _translate_literal(part)
    match = re.fullmatch(source, _basename(path), re.DOTALL)
    return match.group('id') if match else None


# ── Composition helpers ──

def make_range(base: str, start: int, end: int) -> str:
    """Build a range-pattern string.

    make_range('design', 0, 49) → 'design_<0..49>'
    """
    return f"{base}_<{start}..{end}>"


def make_set(base: str, values: List[str]) -> str:
    """Build a set-pattern string.

    make_set('pos', ['42A', '42V', '42W']) → 'pos_<42A 42V 42W>'
    """
    return f"{base}_<{' '.join(values)}>"


def append_suffix(parent_ids: List[str], suffix: str) -> List[str]:
    """Append a pattern suffix to each parent ID.

    append_suffix(['5HG6_<0..4>'], '<1..3>')
    → ['5HG6_<0..4>_<1..3>']
    """
    return [f"{pid}_{suffix}" for pid in parent_ids]
