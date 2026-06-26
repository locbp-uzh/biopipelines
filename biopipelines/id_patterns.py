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
    [_<N><S A L K>]  Bracket segment — cannot be expanded at config time
    prot_<0..4>[_<N><S A L K>]  Mixed: deterministic prefix + lazy suffix

All functions are standalone (no DataStream imports).
"""

import re
from itertools import product
from typing import List, Optional, Tuple

# Regex to find pattern slots: <...> but not inside [...] (runtime bracket syntax)
_SLOT_RE = re.compile(r'<([^>\[\]]+)>')

# Regex to find bracket segments: [...]
_BRACKET_RE = re.compile(r'\[([^\]]+)\]')


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

    'prot_<0..4>[_<N><S A L K>]' → 'prot_<0..4>'
    'base[_<N><A V>]'            → 'base'
    'literal'                    → 'literal'
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


def select_ids(patterns: List[str], row_ids: List[str]) -> List[str]:
    """Select the row ids a pattern set covers, in row order.

    The row ids (read from a runtime map_table) are the source of truth; the
    patterns only *select* among them. Deterministic ``<..>`` slots and literals
    match exactly; ``[...]`` brackets become ``*`` and match by glob. Ids the
    patterns do not cover are dropped; patterns never fabricate ids absent from
    the rows. This honors upstream filtering automatically — a filtered stream
    simply has fewer rows.

    patterns=['a_<1 2>[_<N>]'], rows=['a_1_x','a_2_y','b_1'] → ['a_1_x','a_2_y']
    patterns=['a_<0..2>'],      rows=['a_0','a_2']           → ['a_0','a_2']
    """
    from fnmatch import fnmatchcase
    globs = glob_from_lazy_ids(patterns)
    return [rid for rid in row_ids if any(fnmatchcase(rid, g) for g in globs)]


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
    return select_ids(patterns, [str(v) for v in df['id'].tolist()])


def contains_pattern(s: str) -> bool:
    """True if string contains any <..> pattern slot (outside brackets)."""
    return bool(_SLOT_RE.search(s))


def is_literal(s: str) -> bool:
    """True if string has no pattern slots at all."""
    return not contains_pattern(s)


def _parse_slot(slot_content: str) -> List[str]:
    """Parse a single slot's content into its values.

    '0..49'  → ['0', '1', ..., '49']
    'A B C'  → ['A', 'B', 'C']
    """
    slot_content = slot_content.strip()
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
