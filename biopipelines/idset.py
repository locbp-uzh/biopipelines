# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""An ordered set of ids, and the four operations the framework performs on ids.

An id was a `str` and an id set a `List[str]`, so the rules that govern them — a slot value never crosses the `_` suffix delimiter, `[...]` is a zero-or-more repetition, `+` composes axes while `_` separates a parent from its child — had nowhere to live except as free functions in `id_patterns`, to be remembered at each call site. `IdSet` gives them one home: a function that takes an `IdSet` cannot be handed a bare list that skipped the rule.

The operations are bundling, cartesian product, suffix multiplication and renaming. That is the whole surface, because those are the only ways ids are combined; everything else reads or selects.

This is a pure value type. It holds ids and nothing else — no map_table, no provenance, no file handles — so it stays comparable, hashable and testable in isolation, and equality is over the ids alone.
"""

from typing import Callable, Iterable, Iterator, List, Mapping, Tuple, Union

try:
    from . import id_patterns
except ImportError:
    import os
    import sys
    sys.path.append(os.path.dirname(__file__))
    import id_patterns


# `+` composes axes and `_` separates a parent from its child. Keeping them distinct is what lets
# `prot1+lig1` be told apart from a suffix pattern, so neither may appear inside the other's operand.
AXIS_SEPARATOR = "+"
SUFFIX_DELIMITER = "_"

# A bundled axis is one entity contributing one prefix; an iterated one contributes an id per row.
AXIS_MODES = ("bundle", "each")


class IdSet:
    """An ordered, immutable collection of ids, which may be patterns."""

    __slots__ = ("_ids",)

    def __init__(self, ids: Union[str, Iterable[str]] = ()):
        if isinstance(ids, str):
            ids = (ids,)
        collected: List[str] = []
        for value in ids:
            if not isinstance(value, str):
                raise TypeError(
                    f"An id must be a string, got {type(value).__name__}: {value!r}"
                )
            collected.append(value)
        self._ids: Tuple[str, ...] = tuple(collected)

    # ── identity ──────────────────────────────────────────────────────────────

    @property
    def ids(self) -> List[str]:
        """The ids as declared, patterns unexpanded."""
        return list(self._ids)

    def __len__(self) -> int:
        """How many ids are declared, which is not how many they expand to.

        A pattern counts once here. Use :meth:`count` for the expanded total when it is knowable, and :meth:`enumerated` for the ids themselves.
        """
        return len(self._ids)

    def __iter__(self) -> Iterator[str]:
        return iter(self._ids)

    def __contains__(self, value: object) -> bool:
        return value in self._ids

    def __getitem__(self, index):
        if isinstance(index, slice):
            return IdSet(self._ids[index])
        return self._ids[index]

    def __eq__(self, other: object) -> bool:
        if isinstance(other, IdSet):
            return self._ids == other._ids
        return NotImplemented

    def __hash__(self) -> int:
        return hash(self._ids)

    def __bool__(self) -> bool:
        return bool(self._ids)

    def __repr__(self) -> str:
        return f"IdSet({list(self._ids)!r})"

    # ── the rule ──────────────────────────────────────────────────────────────

    @property
    def is_lazy(self) -> bool:
        """True when any id carries a `[...]` bracket, so its values are not known yet."""
        return any(id_patterns.is_lazy(value) for value in self._ids)

    @property
    def can_expand(self) -> bool:
        """True when every id can be expanded to concrete ids right now."""
        return id_patterns.can_expand_ids(list(self._ids))

    def enumerated(self) -> "IdSet":
        """Every id this set stands for.

        An expandable set yields all its concrete ids. A set that is not fully expandable yields ids that still carry their unexpandable parts, brackets preserved — the honest answer rather than a guess, since a lazy bracket's values only exist once the upstream tool has run.
        """
        return IdSet(id_patterns.partial_expand_ids(list(self._ids)))

    def count(self) -> int:
        """How many ids this set expands to, counting a lazy bracket's prefix only."""
        return id_patterns.count_ids(list(self._ids))

    def select(self, row_ids: Iterable[str], where: str = "") -> "IdSet":
        """The subset of ``row_ids`` that these ids cover, in row order.

        Matching is exact outside the unresolved slots, so `design_1[_<#>]` covers `design_1_7` and never reaches into `design_10_7`.
        """
        return IdSet(id_patterns.select_ids(list(self._ids), list(row_ids), where))

    def deduplicated(self) -> "IdSet":
        """Drop literal ids already covered by a pattern in this set."""
        return IdSet(id_patterns.dedup_parent_children(list(self._ids)))

    # ── the four operations ───────────────────────────────────────────────────

    def bundled(self) -> "IdSet":
        """This whole set as one id: the members joined by `+`, duplicates dropped.

        A bundled axis is one entity however many ids it holds, so it contributes a single prefix to a product instead of multiplying it.
        """
        seen = set()
        unique = [value for value in self._ids if not (value in seen or seen.add(value))]
        return IdSet(AXIS_SEPARATOR.join(unique)) if unique else IdSet()

    def product(self, *others: "IdSet") -> "IdSet":
        """Cartesian product with the other sets, joined by `+`, row-major left to right.

        `IdSet(['p1','p2']).product(IdSet(['l1','l2']))` is `['p1+l1', 'p1+l2', 'p2+l1', 'p2+l2']`.
        """
        for other in others:
            if not isinstance(other, IdSet):
                raise TypeError(f"product expects IdSet, got {type(other).__name__}")
        combined = [(value,) for value in self._ids]
        for other in others:
            combined = [row + (value,) for row in combined for value in other._ids]
        return IdSet(AXIS_SEPARATOR.join(row) for row in combined)

    def multiplied_by_suffix(self, suffix: str) -> "IdSet":
        """Each id gains ``_<suffix>``, turning every parent into a child.

        The suffix may be a pattern, in which case the result is a pattern: appending `<1..3>` to `['5HG6_<0..4>']` gives `['5HG6_<0..4>_<1..3>']`.
        """
        if not isinstance(suffix, str):
            raise TypeError(f"A suffix must be a string, got {type(suffix).__name__}")
        if AXIS_SEPARATOR in suffix:
            raise ValueError(
                f"A suffix may not contain {AXIS_SEPARATOR!r} ({suffix!r}): that separator "
                "composes axes, and mixing the two makes an id impossible to decompose."
            )
        return IdSet(id_patterns.append_suffix(list(self._ids), suffix))

    def renamed(self, mapping: Union[Mapping[str, str], Callable[[str], str]]) -> "IdSet":
        """These ids under a substitution, cardinality preserved.

        A rename changes a cell, not an id count, which is what keeps provenance columns joinable. An id the mapping does not name is left as it is; a mapping that would collapse two ids into one raises, because the collapse would silently lose a row.
        """
        if callable(mapping):
            renamed = [mapping(value) for value in self._ids]
        else:
            renamed = [mapping.get(value, value) for value in self._ids]
        for value in renamed:
            if not isinstance(value, str):
                raise TypeError(f"A rename must produce strings, got {type(value).__name__}")
        if len(set(renamed)) != len(set(self._ids)):
            raise ValueError(
                f"Renaming would collapse distinct ids: {len(set(self._ids))} distinct in, "
                f"{len(set(renamed))} out. A rename must preserve cardinality so provenance "
                "columns stay joinable."
            )
        return IdSet(renamed)


def bundle(*sets: IdSet) -> IdSet:
    """Several bundled axes as one prefix, each axis deduplicated within itself and then joined.

    Deduplication is per axis, not across them: two axes may legitimately carry the same id, and each still contributes it, or the id would stop being reconstructible from the per-axis provenance columns.
    """
    parts: List[str] = []
    for one in sets:
        if not isinstance(one, IdSet):
            raise TypeError(f"bundle expects IdSet, got {type(one).__name__}")
        collapsed = one.bundled()
        if collapsed:
            parts.append(collapsed[0])
    return IdSet(AXIS_SEPARATOR.join(parts)) if parts else IdSet()


def product(*sets: IdSet) -> IdSet:
    """Cartesian product of the given sets, joined by `+`, row-major left to right."""
    if not sets:
        return IdSet()
    first, rest = sets[0], sets[1:]
    if not isinstance(first, IdSet):
        raise TypeError(f"product expects IdSet, got {type(first).__name__}")
    return first.product(*rest)


def decorate_with_static(ids: IdSet, static: IdSet, static_first: bool = False) -> IdSet:
    """Each id joined with the static companions of its own axis, on the side the declared order puts them.

    A `Bundle(Each(a), b)` axis iterates over `a` while `b` rides along on every row, so the static part decorates that axis's contribution rather than becoming an axis of its own.
    """
    if not isinstance(ids, IdSet) or not isinstance(static, IdSet):
        raise TypeError("decorate_with_static expects IdSet operands")
    if not static:
        return ids
    joined = AXIS_SEPARATOR.join(static)
    if static_first:
        return IdSet(f"{joined}{AXIS_SEPARATOR}{value}" for value in ids)
    return IdSet(f"{value}{AXIS_SEPARATOR}{joined}" for value in ids)


def compose_axes(axes: Iterable[Tuple]) -> IdSet:
    """The ids a set of axes produces, following the framework's naming convention.

    Each axis is ``(ids, mode)`` or ``(ids, mode, static, static_first)`` with mode ``"bundle"`` or ``"each"``. A bundled axis is one entity, so it contributes a single `+`-joined prefix instead of multiplying the product; an iterated axis contributes one id per row, decorated with its own static companions if it has any.

    Bundle prefixes come first, ahead of the iterated axes and regardless of declared order, which is why an id cannot be decomposed back by splitting on the separator and pairing positions with declared axes: a bundle occupies as many positions as it has members, and not the position it was declared in. Record what each axis contributed while composing instead.
    """
    bundles, iterated = [], []
    for axis in axes:
        ids, mode = axis[0], axis[1]
        static, static_first = (axis[2], axis[3]) if len(axis) > 2 else (IdSet(), False)
        if not isinstance(ids, IdSet):
            raise TypeError(f"compose_axes expects IdSet, got {type(ids).__name__}")
        if mode not in AXIS_MODES:
            # A misspelling used to fall through to the iterated branch and quietly pick one id.
            raise ValueError(f"Unknown axis mode {mode!r}; expected one of {AXIS_MODES}")
        if mode == "bundle":
            bundles.append(ids)
        else:
            iterated.append(decorate_with_static(ids, static, static_first))

    prefix = bundle(*bundles) if bundles else IdSet()
    if not iterated:
        return prefix
    product_of_iterated = iterated[0].product(*iterated[1:])
    if not prefix:
        return product_of_iterated
    return prefix.product(product_of_iterated)
