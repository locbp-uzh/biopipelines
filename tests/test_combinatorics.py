"""Unit tests for biopipelines.combinatorics — Bundle/Each + ID prediction."""

import json
from types import SimpleNamespace

import pytest

from biopipelines import combinatorics as cb


def _mock_source(stream_name: str, ids, map_table: str = "/fake/map.csv"):
    """Build a minimal stand-in for a StandardizedOutput with one named stream.

    Combinatorics code paths access two things on the source:
      - ``source.streams.<stream_name>.ids``         (via _collect_iterated_ids_from_value)
      - ``source.streams.<stream_name>.map_table``   (via _extract_source_paths)

    The map_table string just needs to be truthy; the combinatorics ID-prediction
    paths exercised in these tests read IDs off ``.ids`` directly rather than
    parsing the CSV.
    """
    stream = SimpleNamespace(ids=list(ids), map_table=map_table)
    streams = SimpleNamespace(**{stream_name: stream})
    return SimpleNamespace(streams=streams)


# ── Each / Bundle wrappers ────────────────────────────────────────────────────

class TestWrappers:
    def test_each_requires_sources(self):
        with pytest.raises(ValueError):
            cb.Each()

    def test_bundle_requires_sources(self):
        with pytest.raises(ValueError):
            cb.Bundle()

    def test_each_stores_sources(self):
        e = cb.Each("a", "b")
        assert e.sources == ("a", "b")

    def test_bundle_stores_sources(self):
        b = cb.Bundle("x")
        assert b.sources == ("x",)

    def test_repr_shapes(self):
        assert repr(cb.Each("a")) == "Each('a')"
        assert repr(cb.Bundle("x", "y")) == "Bundle('x', 'y')"

    def test_is_combinatorics_wrapper(self):
        assert cb.is_combinatorics_wrapper(cb.Each("a")) is True
        assert cb.is_combinatorics_wrapper(cb.Bundle("a")) is True
        assert cb.is_combinatorics_wrapper("plain") is False

    def test_contains_combinatorics_wrapper_in_list(self):
        assert cb.contains_combinatorics_wrapper([cb.Each("a")]) is True
        assert cb.contains_combinatorics_wrapper(["a", "b"]) is False

    def test_get_mode(self):
        assert cb.get_mode(cb.Bundle("a")) == "bundle"
        assert cb.get_mode(cb.Each("a")) == "each"
        assert cb.get_mode("bare") == "each"


# ── AxisConfig / CombinatoricsConfig serialization ────────────────────────────

class TestAxisConfig:
    def test_to_dict_roundtrip(self):
        ax = cb.AxisConfig(
            name="proteins",
            mode="each",
            sources=[{"path": "/tmp/x.csv", "iterate": True}],
            entity_type="protein",
        )
        d = ax.to_dict()
        assert d["name"] == "proteins"
        assert d["mode"] == "each"
        assert d["entity_type"] == "protein"
        assert d["sources"] == [{"path": "/tmp/x.csv", "iterate": True}]


class TestCombinatoricsConfig:
    def test_save_and_load_roundtrip(self, tmp_path):
        ax = cb.AxisConfig(name="proteins", mode="each",
                           sources=[{"path": "p.csv", "iterate": True}],
                           entity_type="protein")
        cfg = cb.CombinatoricsConfig(
            axes={"proteins": ax},
            predicted_ids=["p1", "p2"],
            provenance={"proteins": ["p1", "p2"]},
        )
        path = tmp_path / "combo.json"
        cfg.save(str(path))

        # File is valid JSON
        data = json.loads(path.read_text())
        assert data["predicted_ids"] == ["p1", "p2"]
        assert data["provenance"] == {"proteins": ["p1", "p2"]}
        assert data["axes"]["proteins"]["mode"] == "each"

        # Round-trip through load()
        loaded = cb.CombinatoricsConfig.load(str(path))
        assert loaded.predicted_ids == ["p1", "p2"]
        assert loaded.get_mode("proteins") == "each"
        assert loaded.get_source_paths("proteins") == ["p.csv"]

    def test_get_mode_default_each_for_missing_axis(self):
        cfg = cb.CombinatoricsConfig()
        assert cfg.get_mode("unknown") == "each"
        assert cfg.get_sources("unknown") == []


# ── predict_output_ids / provenance ───────────────────────────────────────────

class TestPredictOutputIds:
    def test_single_axis_returns_input_ids(self):
        source = _mock_source("sequences", ["p1", "p2", "p3"])
        ids = cb.predict_output_ids(proteins=(source, "sequences"))
        assert ids == ["p1", "p2", "p3"]

    def test_two_axis_cartesian_product(self):
        prots = _mock_source("sequences", ["p1", "p2"])
        ligs = _mock_source("compounds", ["l1", "l2", "l3"])
        ids, prov = cb.predict_output_ids_with_provenance(
            proteins=(prots, "sequences"),
            ligands=(ligs, "compounds"),
        )
        assert ids == [
            "p1+l1", "p1+l2", "p1+l3",
            "p2+l1", "p2+l2", "p2+l3",
        ]
        assert prov["proteins"] == ["p1", "p1", "p1", "p2", "p2", "p2"]
        assert prov["ligands"] == ["l1", "l2", "l3", "l1", "l2", "l3"]

    def test_bundle_collapses_to_single(self):
        prots = _mock_source("sequences", ["p1", "p2"])
        ids = cb.predict_output_ids(
            proteins=(cb.Bundle(prots), "sequences"),
        )
        assert ids == ["p1+p2"]

    def test_empty_inputs_raises(self):
        with pytest.raises(ValueError):
            cb.predict_output_ids_with_provenance()


# ── generate_multiplied_ids ───────────────────────────────────────────────────

class TestGenerateMultipliedIds:
    def test_basic_multiply(self):
        out, prov = cb.generate_multiplied_ids(
            ["prot_1", "prot_2"],
            ["1", "2", "3"],
            input_stream_name="structures",
        )
        assert out == [
            "prot_1_1", "prot_1_2", "prot_1_3",
            "prot_2_1", "prot_2_2", "prot_2_3",
        ]
        assert prov == {
            "structures": [
                "prot_1", "prot_1", "prot_1",
                "prot_2", "prot_2", "prot_2",
            ]
        }

    def test_no_stream_name_empty_provenance(self):
        out, prov = cb.generate_multiplied_ids(["a"], ["1", "2"])
        assert out == ["a_1", "a_2"]
        assert prov == {}


class TestGenerateMultipliedIdsPattern:
    def test_compose_pattern(self):
        assert cb.generate_multiplied_ids_pattern(
            ["5HG6_<0..4>"], "<1..3>"
        ) == ["5HG6_<0..4>_<1..3>"]


# ── predict_single_output_id ──────────────────────────────────────────────────

class TestPredictSingleOutputId:
    def test_cartesian_selection(self):
        out = cb.predict_single_output_id(
            sequences=("each", ["prot1", "prot2"], 0, [], False),
            compounds=("each", ["lig1", "lig2", "lig3"], 1, [], False),
        )
        assert out == "prot1+lig2"

    def test_single_axis(self):
        out = cb.predict_single_output_id(
            sequences=("each", ["only"], 0, [], False),
        )
        assert out == "only"


# ── Mock-driven combinatorics end-to-end ──────────────────────────────────────

def test_combinatorics_each_each_via_mock(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """Two Mocks feeding a downstream Mock via Each/Each → cartesian IDs."""
    from biopipelines.mock import Mock
    from biopipelines.combinatorics import Each

    pipeline = new_pipeline("comb_each_each_mock")
    with pipeline:
        a = Mock(ids=["a1", "a2"],
                 streams={"sequences": {"format": "fa", "file": "<id>.fa"}})
        b = Mock(ids=["b1", "b2"],
                 streams={"compounds": {"format": "sdf", "file": "<id>.sdf"}})
        c = Mock(
            source=[Each(a.streams.sequences), Each(b.streams.compounds)],
            streams={"pairs": {"format": "pdb", "file": "<id>.pdb"}},
        )
        pipeline.save()

    ids = list(c.streams.pairs.ids)
    record_case(input="Each(Mock_a)+Each(Mock_b)",
                expected=["a1+b1", "a1+b2", "a2+b1", "a2+b2"],
                actual=ids)
    assert ids == ["a1+b1", "a1+b2", "a2+b1", "a2+b2"]


def test_combinatorics_bundle_via_mock(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """Bundle(Mock_a, Mock_b) → single bundled axis (no cartesian blow-up)."""
    from biopipelines.mock import Mock
    from biopipelines.combinatorics import Bundle

    pipeline = new_pipeline("comb_bundle_mock")
    with pipeline:
        a = Mock(ids=["x1", "x2"],
                 streams={"sequences": {"format": "fa", "file": "<id>.fa"}})
        b = Mock(ids=["y1", "y2"],
                 streams={"compounds": {"format": "sdf", "file": "<id>.sdf"}})
        c = Mock(
            source=Bundle(a.streams.sequences, b.streams.compounds),
            streams={"out": {"format": "pdb", "file": "<id>.pdb"}},
        )
        pipeline.save()

    ids = list(c.streams.out.ids)
    record_case(input="Bundle(Mock_a, Mock_b)",
                expected=("len < 4 (no cartesian)", True),
                actual=(f"len={len(ids)}", len(ids) < 4))
    assert len(ids) < 4  # Each×Each would give 4


def test_combinatorics_mixed_bundle_each_via_mock(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """Bundle(Each(a), b) — iterate over a, hold b constant per-combo."""
    from biopipelines.mock import Mock
    from biopipelines.combinatorics import Bundle, Each

    pipeline = new_pipeline("comb_mixed_mock")
    with pipeline:
        a = Mock(ids=["a1", "a2", "a3"],
                 streams={"sequences": {"format": "fa", "file": "<id>.fa"}})
        b = Mock(ids=["b1"],
                 streams={"compounds": {"format": "sdf", "file": "<id>.sdf"}})
        c = Mock(
            source=Bundle(Each(a.streams.sequences), b.streams.compounds),
            streams={"out": {"format": "pdb", "file": "<id>.pdb"}},
        )
        pipeline.save()

    ids = list(c.streams.out.ids)
    record_case(input="Bundle(Each(a[3]), b[1])",
                expected=("len == 3 (iterate a, hold b)", 3),
                actual=(f"len={len(ids)}", len(ids)))
    assert len(ids) == 3


# ── CombinatoricsConfig round-trip ────────────────────────────────────────────

def test_combinatorics_config_to_dict_roundtrip(record_case):
    """CombinatoricsConfig.to_dict() preserves axes, predicted_ids, provenance."""
    axis = cb.AxisConfig(
        name="sequences",
        mode="each",
        sources=[{"path": "/fake/a.csv", "iterate": True}],
        entity_type="protein",
    )
    cfg = cb.CombinatoricsConfig(
        axes={"sequences": axis},
        predicted_ids=["a", "b"],
        provenance={"sequences": ["a", "b"]},
    )
    d = cfg.to_dict()
    record_case(
        input="CombinatoricsConfig with 1 axis, 2 ids",
        expected=("predicted_ids and provenance present", True),
        actual=("keys",
                "predicted_ids" in d and "provenance" in d and "axes" in d),
    )
    assert d["predicted_ids"] == ["a", "b"]
    assert d["provenance"] == {"sequences": ["a", "b"]}
    assert "sequences" in d["axes"]
    assert d["axes"]["sequences"]["sources"] == [
        {"path": "/fake/a.csv", "iterate": True}
    ]
