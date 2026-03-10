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


def try_expand(s: str) -> Tuple[List[str], bool]:
    """Expand the deterministic prefix of a pattern.

    Returns:
        (ids, is_complete) where:
        - If no brackets: ids = full expansion, is_complete = True
        - If brackets: ids = expansion of prefix only, is_complete = False

    'prot_<0..2>'                → (['prot_0', 'prot_1', 'prot_2'], True)
    'prot_<0..1>[_<N><A V>]'     → (['prot_0', 'prot_1'], False)
    'base[_<N><A V>]'            → (['base'], False)
    """
    if not is_lazy(s):
        return expand_pattern(s), True
    prefix = strip_brackets(s)
    return expand_pattern(prefix), False


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
    """Expand deterministic prefixes of all IDs.

    Returns:
        (expanded_ids, is_complete) — is_complete is False if any ID had brackets.
    """
    result = []
    complete = True
    for s in ids:
        expanded, is_complete = try_expand(s)
        result.extend(expanded)
        if not is_complete:
            complete = False
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
        result = result[:m.start()] + slot_vals[si] + result[m.end()]
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
