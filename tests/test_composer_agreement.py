# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Composing a whole id set must equal composing each row and collecting.

`predict_single_output_id`'s docstring claims it "mirrors predict_output_ids() logic but for a single output row, ensuring pipeline-time and SLURM-time ID generation always agree". Two implementations of one rule agree only until one of them is edited, and nothing checked it: the claim was aspirational. These cases make it a property, swept over every axis-mode combination, so a future divergence fails here instead of on a compute node as a FAILED completion marker.

The same property is asserted against `IdSet`'s bundling and product, which is the type the rule is meant to live in -- so all three routes to an id have to agree, not just the two functions.
"""

import itertools

import pytest

from biopipelines import combinatorics as cb
from biopipelines.idset import IdSet, compose_axes


AXIS_IDS = {
    "compounds": ["LIG1", "LIG2"],
    "structures": ["seedA_1", "seedA_2", "seedA_3"],
    "sequences": ["s1"],
}


def _whole_set(modes):
    """Compose every row at once, the configuration-time route."""
    selections = {name: (modes[name], AXIS_IDS[name], None, [], False) for name in modes}
    iterated = [n for n in modes if modes[n] != "bundle"]
    if not iterated:
        return [cb.predict_single_output_id(**selections)]
    rows = []
    for combination in itertools.product(*[range(len(AXIS_IDS[n])) for n in iterated]):
        picked = dict(zip(iterated, combination))
        rows.append(cb.predict_single_output_id(**{
            name: (modes[name], AXIS_IDS[name],
                   picked.get(name) if modes[name] != "bundle" else None, [], False)
            for name in modes
        }))
    return rows


# The composed ids for every mode combination, written out rather than recomputed.
#
# `_whole_set` and `_per_row` both called `predict_single_output_id`, differing only in how the test
# enumerated indices, and `_via_idset` called `compose_axes` -- which `predict_single_output_id` also
# calls. So the three "agreement" tests compared one implementation with itself: inverting the
# bundle-hoisting convention left all sixteen of them green. A literal table is the second opinion.
#
# The convention they pin: every bundled axis is joined and hoisted ahead of the iterated axes, both
# in declared axis order; the iterated axes then vary right-most-fastest.
EXPECTED_TWO_AXIS = {
    ("each", "each"): ["LIG1+seedA_1", "LIG1+seedA_2", "LIG1+seedA_3",
                       "LIG2+seedA_1", "LIG2+seedA_2", "LIG2+seedA_3"],
    ("each", "bundle"): ["seedA_1+seedA_2+seedA_3+LIG1", "seedA_1+seedA_2+seedA_3+LIG2"],
    ("bundle", "each"): ["LIG1+LIG2+seedA_1", "LIG1+LIG2+seedA_2", "LIG1+LIG2+seedA_3"],
    ("bundle", "bundle"): ["LIG1+LIG2+seedA_1+seedA_2+seedA_3"],
}

EXPECTED_THREE_AXIS = {
    ("each", "each", "each"): ["LIG1+seedA_1+s1", "LIG1+seedA_2+s1", "LIG1+seedA_3+s1",
                               "LIG2+seedA_1+s1", "LIG2+seedA_2+s1", "LIG2+seedA_3+s1"],
    ("each", "each", "bundle"): ["s1+LIG1+seedA_1", "s1+LIG1+seedA_2", "s1+LIG1+seedA_3",
                                 "s1+LIG2+seedA_1", "s1+LIG2+seedA_2", "s1+LIG2+seedA_3"],
    ("each", "bundle", "each"): ["seedA_1+seedA_2+seedA_3+LIG1+s1",
                                 "seedA_1+seedA_2+seedA_3+LIG2+s1"],
    ("each", "bundle", "bundle"): ["seedA_1+seedA_2+seedA_3+s1+LIG1",
                                   "seedA_1+seedA_2+seedA_3+s1+LIG2"],
    ("bundle", "each", "each"): ["LIG1+LIG2+seedA_1+s1", "LIG1+LIG2+seedA_2+s1",
                                 "LIG1+LIG2+seedA_3+s1"],
    ("bundle", "each", "bundle"): ["LIG1+LIG2+s1+seedA_1", "LIG1+LIG2+s1+seedA_2",
                                   "LIG1+LIG2+s1+seedA_3"],
    ("bundle", "bundle", "each"): ["LIG1+LIG2+seedA_1+seedA_2+seedA_3+s1"],
    ("bundle", "bundle", "bundle"): ["LIG1+LIG2+seedA_1+seedA_2+seedA_3+s1"],
}


def _via_idset(modes, order):
    """Compose the same ids through IdSet, using the axis-aware entry point.

    `product` alone cannot reproduce the convention: `bundled()` collapses an axis to a plain one-element set, so by the time a product runs there is nothing left to say which operand was a bundle, and the convention hoists bundle prefixes to the front. `compose_axes` is where that knowledge belongs.
    """
    return compose_axes([(IdSet(AXIS_IDS[name]), modes[name]) for name in order]).ids


ALL_MODES = ["each", "bundle"]

TWO_AXIS_NAMES = ["compounds", "structures"]
THREE_AXIS_NAMES = ["compounds", "structures", "sequences"]


@pytest.mark.parametrize("combination", sorted(EXPECTED_TWO_AXIS))
def test_the_configuration_time_composer_matches_the_written_table(combination):
    modes = dict(zip(TWO_AXIS_NAMES, combination))
    assert _whole_set(modes) == EXPECTED_TWO_AXIS[combination]


@pytest.mark.parametrize("combination", sorted(EXPECTED_THREE_AXIS))
def test_the_table_holds_for_three_axes_too(combination):
    modes = dict(zip(THREE_AXIS_NAMES, combination))
    assert _whole_set(modes) == EXPECTED_THREE_AXIS[combination]


@pytest.mark.parametrize("combination", sorted(EXPECTED_TWO_AXIS))
def test_idset_composes_the_same_ids_as_the_written_table(combination):
    """If IdSet is to own the rule, it must produce what the rule produces -- checked against the table, not against the other caller of the same function."""
    modes = dict(zip(TWO_AXIS_NAMES, combination))
    assert _via_idset(modes, TWO_AXIS_NAMES) == EXPECTED_TWO_AXIS[combination]


def test_a_bundled_axis_contributes_one_prefix_however_many_ids_it_holds():
    modes = {"compounds": "bundle", "structures": "each"}
    assert _whole_set(modes) == [
        "LIG1+LIG2+seedA_1", "LIG1+LIG2+seedA_2", "LIG1+LIG2+seedA_3"]
    assert _via_idset(modes, ["compounds", "structures"]) == _whole_set(modes)


def test_every_axis_bundled_yields_exactly_one_id():
    modes = {"compounds": "bundle", "structures": "bundle"}
    assert len(_whole_set(modes)) == 1


def test_a_bundle_is_hoisted_ahead_of_the_iterated_axes_whatever_the_declared_order():
    """This is why an id cannot be decomposed by splitting on `+` and pairing positions with declared axes: the bundled axis here was declared second and lands first, occupying as many positions as it has members."""
    modes = {"compounds": "each", "structures": "bundle"}
    assert _whole_set(modes) == [
        "seedA_1+seedA_2+seedA_3+LIG1", "seedA_1+seedA_2+seedA_3+LIG2"]
    assert _via_idset(modes, ["compounds", "structures"]) == _whole_set(modes)


def test_product_on_its_own_respects_declared_order():
    """`product` is the honest cartesian product and does not know about modes; the convention lives in compose_axes."""
    assert IdSet(["a", "b"]).product(IdSet(["c"])).ids == ["a+c", "b+c"]


# ── the static dimension, from Bundle(Each(a), b) ─────────────────────────────

@pytest.mark.parametrize("static, static_first, expected", [
    ([], False, ["LIG1", "LIG2"]),
    (["s1"], False, ["LIG1+s1", "LIG2+s1"]),
    (["s1"], True, ["s1+LIG1", "s1+LIG2"]),
    (["s1", "s2"], False, ["LIG1+s1+s2", "LIG2+s1+s2"]),
    (["s1", "s2"], True, ["s1+s2+LIG1", "s1+s2+LIG2"]),
])
def test_static_companions_decorate_their_own_axis(static, static_first, expected):
    """A `Bundle(Each(a), b)` axis iterates over `a` while `b` rides along on every row, so the static part decorates that axis's contribution rather than becoming an axis of its own."""
    ids = AXIS_IDS["compounds"]
    got = [cb.predict_single_output_id(
        compounds=("each", ids, i, static, static_first)) for i in range(len(ids))]
    assert got == expected


def test_a_static_decorated_axis_still_products_with_the_others():
    """The static part attaches to its own axis before the axes multiply, so it lands between the two axis contributions rather than at either end of the id."""
    got = [cb.predict_single_output_id(
        compounds=("each", ["LIG1", "LIG2"], i, ["sx"], False),
        sequences=("each", ["q1"], 0, [], False))
        for i in range(2)]
    assert got == ["LIG1+sx+q1", "LIG2+sx+q1"]


def test_the_runtime_composer_is_implemented_on_idset():
    """Not two implementations that happen to agree: predict_single_output_id calls compose_axes, so there is one rule for the runtime path. Verified byte-identical over a 5634-case sweep of axis counts, modes, static sets and orderings when it was moved."""
    import inspect
    source = inspect.getsource(cb.predict_single_output_id)
    assert "compose_axes" in source


# ── the whole-set predictor, which nothing above actually calls ───────────────

def _predict_whole_set(streams, modes, tmp_path):
    """The real configuration-time predictor, not a per-row loop standing in for it."""
    from biopipelines.combinatorics import Bundle, Each, predict_output_ids_with_provenance
    from biopipelines.datastream import DataStream, create_map_table
    from biopipelines.outputs import StandardizedOutput

    named = {}
    for index, (name, ids) in enumerate(streams):
        table = tmp_path / f"agree_{index}_{name}.csv"
        create_map_table(str(table), ids=ids)
        stream = DataStream(name=name, ids=list(ids), map_table=str(table), format="csv")
        wrapper = Bundle if modes[name] == "bundle" else Each
        named[name] = (wrapper(StandardizedOutput({name: stream})), name)
    return predict_output_ids_with_provenance(**named)[0]


def _produce_per_row(streams, modes):
    ids_by_name = dict(streams)
    iterated = [n for n, _ in streams if modes[n] != "bundle"]
    index_space = [range(len(ids_by_name[n])) for n in iterated] or [[None]]
    rows = []
    for combination in itertools.product(*index_space):
        picked = dict(zip(iterated, combination))
        rows.append(cb.predict_single_output_id(**{
            name: (modes[name], ids_by_name[name],
                   picked.get(name) if modes[name] != "bundle" else None, [], False)
            for name, _ in streams
        }))
    return rows


DISTINCT = [("sequences", ["X", "W"]), ("compounds", ["Y", "Z"]), ("structures", ["s1"])]
# Two bundled axes carrying the same id: the case a sweep over distinct ids never reaches.
SHARED = [("sequences", ["X", "Y"]), ("compounds", ["Y", "Z"]), ("structures", ["s1"])]


@pytest.mark.parametrize("streams", [DISTINCT, SHARED], ids=["distinct_ids", "shared_id"])
@pytest.mark.parametrize("first, second", list(itertools.product(ALL_MODES, repeat=2)))
def test_the_two_predictors_agree(tmp_path, streams, first, second):
    """The claim that per-row and whole-set generation agree was only ever checked between two per-row loops -- `predict_output_ids_with_provenance` was never called. It deduplicates a bundle per axis while `bundle()` deduplicated across all of them, so two bundled axes sharing an id declared one id and produced another, and the step got a FAILED marker."""
    modes = {streams[0][0]: first, streams[1][0]: second, streams[2][0]: "each"}
    assert _produce_per_row(streams, modes) == _predict_whole_set(streams, modes, tmp_path)


def test_two_bundled_axes_sharing_an_id_each_still_contribute_it(tmp_path):
    modes = {"sequences": "bundle", "compounds": "bundle", "structures": "each"}
    assert _predict_whole_set(SHARED, modes, tmp_path) == ["X+Y+Y+Z+s1"]
    assert _produce_per_row(SHARED, modes) == ["X+Y+Y+Z+s1"]
