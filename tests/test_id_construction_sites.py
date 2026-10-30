# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Ids and provenance that Fuse, BioEmu2 and Mock used to build by hand.

The repo rule is that an output id or a provenance value is produced only by `combinatorics.py`. Three sites wrote their own: Fuse joined its slot and junction parts with "+" at configuration time and again in its pipe script, BioEmu2 built each run's provenance with a "+" join of the records it had picked (static companions included, which the configuration-time predictor leaves out), and Mock took its bundled provenance from a join of its own and carried a positional split of the id that could not be reached.

Nothing here runs a model or a cluster job; the pipe scripts are driven on synthetic tables.
"""

import itertools
import os
import re
import sys
from types import SimpleNamespace

import pandas as pd
import pytest

from biopipelines import combinatorics as cb
from biopipelines import id_patterns
from biopipelines.combinatorics import (Bundle, CombinatoricsConfig, Each, generate_combinatorics_config,
                                        predict_output_ids_with_provenance)
from biopipelines.datastream import DataStream
from biopipelines.fuse import Fuse

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(REPO_ROOT, "pipe_scripts")
sys.path.insert(0, SCRIPTS)

import pipe_bioemu2 as bioemu2_pipe  # noqa: E402
import pipe_fuse_queries as fuse_pipe  # noqa: E402

PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"


def _script_text(name):
    with open(os.path.join(SCRIPTS, name)) as handle:
        return handle.read()


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


# ── Fuse ──────────────────────────────────────────────────────────────────────

def _fuse_slot(tmp_path, name, ids):
    table = tmp_path / f"{name}.csv"
    pd.DataFrame({"id": ids, "sequence": ["ACDE"[: 2 + i % 3] for i in range(len(ids))]}).to_csv(
        table, index=False)
    return {"ids": list(ids), "map_table": str(table), "files": []}


def _fuse_config_side(slots, linker_lengths):
    fuse = Fuse.__new__(Fuse)
    fuse.input_slots = slots
    fuse.linker_lengths = linker_lengths
    return fuse._predict_sequence_ids()


FUSE_CASES = {
    "no linkers": (None, [["a1", "a2"], ["b1"]]),
    "one junction": (["1-3"], [["a1", "a2"], ["b1", "b2"]]),
    "a junction without a linker": ([None, "2-3"], [["a1"], ["b1", "b2"], ["c1"]]),
    "a set of lengths": (["2+5"], [["a1"], ["b1"]]),
}


class TestFuse:

    @pytest.mark.parametrize("case", sorted(FUSE_CASES))
    def test_the_ids_the_pipe_script_writes_are_the_ids_the_wrapper_declares(self, tmp_path, case):
        linker_lengths, slot_ids = FUSE_CASES[case]
        slots = [_fuse_slot(tmp_path, f"slot{i}", ids) for i, ids in enumerate(slot_ids)]
        declared, declared_provenance = _fuse_config_side(slots, linker_lengths)
        rows = fuse_pipe.generate_fusion_sequences(slots, "GS", linker_lengths or [], "fuse")
        assert [row["id"] for row in rows] == declared
        assert len(set(declared)) == len(declared), "two distinct rows share an id"
        for index in range(len(slots)):
            assert [row[f"sequences_{index + 1}.id"] for row in rows] == declared_provenance[
                f"sequences_{index + 1}"]

    def test_the_documented_shape_is_slot_then_length_then_slot_joined_by_plus(self, tmp_path):
        slots = [_fuse_slot(tmp_path, "s0", ["a1", "a2"]), _fuse_slot(tmp_path, "s1", ["b1"])]
        declared, _ = _fuse_config_side(slots, ["1-2"])
        assert declared == ["a1+1+b1", "a1+2+b1", "a2+1+b1", "a2+2+b1"]

    def test_the_wrapper_and_its_pipe_script_build_the_id_with_the_one_shared_function(self):
        import biopipelines.fuse as fuse_module
        assert fuse_module.join_axis_ids is cb.join_axis_ids
        assert fuse_pipe.join_axis_ids is cb.join_axis_ids

    def test_neither_side_joins_the_parts_itself(self):
        assert not re.search(r"""['"]\+['"]\.join""", _script_text("pipe_fuse_queries.py"))
        with open(os.path.join(REPO_ROOT, "biopipelines", "fuse.py")) as handle:
            assert not re.search(r"""['"]\+['"]\.join""", handle.read())


# ── BioEmu2 ───────────────────────────────────────────────────────────────────

def _bioemu2_streams(tmp_path):
    pd.DataFrame({"id": ["p1", "p2"], "sequence": [PROTEIN, PROTEIN]}).to_csv(tmp_path / "s.csv", index=False)
    pd.DataFrame({"id": ["l1", "l2"], "smiles": ["CCO", "CCN"]}).to_csv(tmp_path / "l.csv", index=False)
    pd.DataFrame({"id": ["cof"], "smiles": ["CC"]}).to_csv(tmp_path / "c.csv", index=False)
    return (DataStream(name="sequences", ids=["p1", "p2"], map_table=str(tmp_path / "s.csv"), format="csv"),
            DataStream(name="compounds", ids=["l1", "l2"], map_table=str(tmp_path / "l.csv"), format="csv"),
            DataStream(name="compounds", ids=["cof"], map_table=str(tmp_path / "c.csv"), format="csv"))


class TestBioEmu2:

    @pytest.mark.parametrize("shape", ["each", "bundle", "static after", "static first"])
    def test_a_runs_provenance_is_what_the_wrapper_declares_for_that_row(self, tmp_path, shape):
        """A ligand axis with a static companion, `Bundle(Each(ligands), cofactor)`, was recorded as the companion joined to the ligand (`cof+l1`), which is no id of the compounds stream and so joins to nothing."""
        seqs, ligs, cof = _bioemu2_streams(tmp_path)
        ligands = {"each": ligs, "bundle": Bundle(ligs), "static after": Bundle(Each(ligs), cof),
                   "static first": Bundle(cof, Each(ligs))}[shape]
        path = str(tmp_path / "c.json")
        declared = generate_combinatorics_config(
            path, sequences=(seqs, "sequences", "protein"), ligands=(ligands, "compounds", "ligand"))
        runs = bioemu2_pipe.build_runs(bioemu2_pipe.axis_elements(CombinatoricsConfig.load(path)))
        assert [run["id"] for run in runs] == declared.predicted_ids
        for position, run in enumerate(runs):
            for column, values in declared.provenance.items():
                assert run["provenance"][column] == values[position], (run["id"], column)

    def test_a_conformer_is_named_by_its_run_and_its_one_based_index(self):
        assert bioemu2_pipe.frame_ids("p1+l1", 3) == ["p1+l1_1", "p1+l1_2", "p1+l1_3"]

    @pytest.mark.parametrize("count", [1, 4, 12])
    def test_the_conformer_ids_are_the_ones_the_wrapper_declares_as_a_pattern(self, count):
        declared = id_patterns.expand_ids(
            cb.generate_multiplied_ids_pattern(["p1+l1"], f"<1..{count}>"))
        assert bioemu2_pipe.frame_ids("p1+l1", count) == declared

    def test_the_pipe_script_joins_and_splits_no_id_itself(self):
        text = _script_text("pipe_bioemu2.py")
        assert not re.search(r"""['"]\+['"]\.join""", text)
        assert not re.search(r"""split\(['"]\+['"]""", text)
        assert not re.search(r"""f['"]\{run_id\}_\{""", text)
