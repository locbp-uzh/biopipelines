# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Ids and provenance that Fuse, BioEmu2 and Mock used to build by hand.

The repo rule is that an output id or a provenance value is produced only by `combinatorics.py`. Three sites wrote their own: Fuse joined its slot and junction parts with "+" at configuration time and again in its pipe script, BioEmu2 built each run's provenance with a "+" join of the records it had picked (static companions included, which the configuration-time predictor leaves out), and Mock took its bundled provenance from a join of its own and carried a positional split of the id that could not be reached.

Nothing here runs a model or a cluster job; the pipe scripts are driven on synthetic tables.
"""

import itertools
from types import SimpleNamespace

import pytest

from biopipelines import combinatorics as cb
from biopipelines.combinatorics import Bundle, Each, predict_output_ids_with_provenance


# ── the shared functions ──────────────────────────────────────────────────────

class TestJoinAxisIds:

    def test_parts_are_joined_by_the_axis_separator_and_numbers_become_text(self):
        assert cb.join_axis_ids(["donor", 3, "cam"]) == "donor+3+cam"

    def test_one_part_is_that_part(self):
        assert cb.join_axis_ids(["p1"]) == "p1"

    def test_no_parts_is_refused_rather_than_named_with_an_empty_id(self):
        with pytest.raises(ValueError):
            cb.join_axis_ids([])

    def test_it_is_the_composition_the_framework_uses_for_iterated_axes(self):
        picked = cb.predict_single_output_id(
            sequences=("each", ["a", "b"], 1, [], False), ligands=("each", ["x", "y"], 0, [], False))
        assert picked == cb.join_axis_ids(["b", "x"])


def _source(stream_name, ids):
    """A stand-in for a StandardizedOutput with one named stream, which is all the predictor reads."""
    stream = SimpleNamespace(ids=list(ids), map_table="/fake/map.csv")
    return SimpleNamespace(streams=SimpleNamespace(**{stream_name: stream}))


def _selections(shape):
    ids = {"sequences": ["s1", "s2"], "ligands": ["l1", "l2", "l3"], "cofactors": ["c1"]}
    return {name: (mode, ids[name], idx, static, static_first)
            for name, (mode, idx, static, static_first) in shape.items()}


PROVENANCE_SHAPES = {
    "each x each": {"sequences": ("each", 1, [], False), "ligands": ("each", 2, [], False)},
    "each x bundle": {"sequences": ("each", 0, [], False), "ligands": ("bundle", None, [], False)},
    "bundle x each": {"sequences": ("bundle", None, [], False), "ligands": ("each", 1, [], False)},
    "static after": {"ligands": ("each", 1, ["c1"], False)},
    "static first": {"ligands": ("each", 1, ["c1"], True)},
    "single bundle": {"cofactors": ("bundle", None, [], False)},
}


class TestPredictSingleProvenance:

    def test_an_iterated_axis_gives_the_id_it_selected_and_never_its_static_companions(self):
        provenance = cb.predict_single_provenance(**_selections(PROVENANCE_SHAPES["static first"]))
        assert provenance == {"ligands": "l2"}

    def test_a_bundled_axis_gives_its_ids_joined_and_each_member_when_there_are_several(self):
        provenance = cb.predict_single_provenance(**_selections(PROVENANCE_SHAPES["each x bundle"]))
        assert provenance == {"sequences": "s1", "ligands": "l1+l2+l3",
                              "ligands.1": "l1", "ligands.2": "l2", "ligands.3": "l3"}

    def test_a_bundle_of_one_id_has_no_member_column(self):
        assert cb.predict_single_provenance(**_selections(PROVENANCE_SHAPES["single bundle"])) == {
            "cofactors": "c1"}

    def test_an_unknown_mode_is_refused(self):
        with pytest.raises(ValueError, match="unknown mode"):
            cb.predict_single_provenance(sequences=("eech", ["s1"], 0, [], False))

    @pytest.mark.parametrize("shape", sorted(PROVENANCE_SHAPES))
    def test_it_agrees_with_the_whole_set_predictor_row_for_row(self, shape):
        """The configuration-time predictor writes the provenance columns; a pipe script that recomputes them per row must land on the same cells."""
        selections = _selections(PROVENANCE_SHAPES[shape])
        sources = {}
        for name, (mode, ids, idx, static, static_first) in selections.items():
            source = _source(name, ids)
            if static:
                companion = _source(name, static)
                pair = (companion, Each(source)) if static_first else (Each(source), companion)
                sources[name] = (Bundle(*pair), name)
            else:
                sources[name] = (Bundle(source) if mode == "bundle" else source, name)
        _ids, whole = predict_output_ids_with_provenance(**sources)

        iterated = [name for name, sel in selections.items() if sel[0] == "each"]
        rows = itertools.product(*[range(len(selections[name][1])) for name in iterated]) if iterated else [()]
        picked_rows = [dict(zip(iterated, row)) for row in rows]
        assert len(picked_rows) == len(next(iter(whole.values())))
        for position, picked in enumerate(picked_rows):
            row_selections = {name: (mode, ids, picked.get(name), static, static_first)
                              for name, (mode, ids, _idx, static, static_first) in selections.items()}
            single = cb.predict_single_provenance(**row_selections)
            assert single == {key: values[position] for key, values in whole.items()}
