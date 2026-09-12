# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""What Mock declares it will produce is what Mock produces.

Mock predicted ids at configuration time through `predict_output_ids_with_provenance`, which honors Bundle and Each, but `_collect_datastreams` was documented as flattening Bundle/Each into a plain DataStream list and did exactly that -- so the config JSON carried no mode and `pipe_mock.py` composed ids with a hand-rolled cartesian product over every axis. A bundled axis that should contribute one prefix iterated instead: two declared ids against four produced, and a FAILED completion marker.

Mock exists to exercise the framework, so a naming rule it cannot express is a hole in the test surface for every real tool. These cases pin that both sides go through one composer.
"""

import importlib.util
import itertools
import os

import pytest

from biopipelines.combinatorics import Bundle, Each
from biopipelines.datastream import DataStream, create_map_table
from biopipelines.mock import Mock


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _pipe_mock():
    spec = importlib.util.spec_from_file_location(
        "pipe_mock_under_test", os.path.join(REPO_ROOT, "pipe_scripts", "pipe_mock.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_table_counter = itertools.count()


def _stream(tmp_path, name, ids):
    # The table filename is counted, not named after the stream: two streams in one axis may share a name.
    table = tmp_path / f"map_{next(_table_counter)}_{name}.csv"
    create_map_table(str(table), ids=ids)
    return DataStream(name=name, ids=list(ids), map_table=str(table), format="csv")


def _predict_and_produce(tmp_path, source, streams):
    parents, _prov, _axis_names, axes = Mock._resolve_from_source(source)

    source_streams = []
    for stream in streams:
        json_path = tmp_path / f"{stream.name}.json"
        stream.save_json(str(json_path))
        source_streams.append({"name": stream.name, "path": str(json_path)})

    produced, _contributions = _pipe_mock()._resolve_parent_ids(
        {"parent_ids": parents, "axes": axes, "source_streams": source_streams})
    return parents, produced


# ── the reported failure ──────────────────────────────────────────────────────

def test_a_bundled_axis_contributes_a_prefix_not_an_iteration(tmp_path):
    compounds = _stream(tmp_path, "compounds", ["LIG1", "LIG2"])
    structures = _stream(tmp_path, "structures", ["seedA_1", "seedA_2"])

    predicted, produced = _predict_and_produce(
        tmp_path, [Bundle(compounds), Each(structures)], [compounds, structures])

    assert predicted == ["LIG1+LIG2+seedA_1", "LIG1+LIG2+seedA_2"]
    assert produced == predicted, "a bundled axis iterated instead of collapsing to one prefix"


# ── every axis-mode combination ───────────────────────────────────────────────

@pytest.mark.parametrize("first_mode, second_mode", [
    (Each, Each),
    (Bundle, Each),
    (Each, Bundle),
    (Bundle, Bundle),
])
def test_predict_equals_produce_for_every_axis_mode(tmp_path, first_mode, second_mode):
    compounds = _stream(tmp_path, "compounds", ["LIG1", "LIG2"])
    structures = _stream(tmp_path, "structures", ["seedA_1", "seedA_2", "seedA_3"])

    predicted, produced = _predict_and_produce(
        tmp_path,
        [first_mode(compounds), second_mode(structures)],
        [compounds, structures],
    )
    assert produced == predicted


def test_a_single_iterated_axis_is_unchanged(tmp_path):
    structures = _stream(tmp_path, "structures", ["p1", "p2", "p3"])
    predicted, produced = _predict_and_produce(tmp_path, Each(structures), [structures])
    assert predicted == ["p1", "p2", "p3"]
    assert produced == predicted


def test_a_config_without_axes_still_resolves(tmp_path):
    """A config written before Mock recorded its modes must keep working: one stream is one iterated axis."""
    structures = _stream(tmp_path, "structures", ["p1", "p2"])
    json_path = tmp_path / "structures.json"
    structures.save_json(str(json_path))

    produced, _contributions = _pipe_mock()._resolve_parent_ids({
        "parent_ids": ["p1", "p2"],
        "source_streams": [{"name": "structures", "path": str(json_path)}],
    })
    assert produced == ["p1", "p2"]


# ── provenance, which was recovered by splitting and could not be ─────────────

def _provenance_both_ways(tmp_path, source, streams):
    parents, prov, axis_names, axes = Mock._resolve_from_source(source)
    source_streams = []
    for stream in streams:
        json_path = tmp_path / f"{stream.name}.json"
        stream.save_json(str(json_path))
        source_streams.append({"name": stream.name, "path": str(json_path)})

    module = _pipe_mock()
    cfg = {"parent_ids": parents, "axes": axes, "provenance": prov,
           "axis_names": axis_names, "source_streams": source_streams}
    _ids, contributions = module._resolve_parent_ids(cfg)
    names = list(prov.keys())
    runtime = module._provenance_from_contributions(contributions, names)
    return {name: prov[name] for name in names}, runtime


def test_runtime_provenance_matches_what_was_declared(tmp_path):
    """`_rebuild_runtime_provenance` split the parent id on '+' and assigned part i to axis i. A bundled axis contributes several parts and the composer hoists it ahead of the iterated axes, so every column came out wrong: with Bundle(compounds) and Each(structures), `structures` was given a ligand id. Provenance is now kept from composition, where the contribution of each axis is still known."""
    compounds = _stream(tmp_path, "compounds", ["LIG1", "LIG2"])
    structures = _stream(tmp_path, "structures", ["seedA_1", "seedA_2"])
    declared, runtime = _provenance_both_ways(
        tmp_path, [Bundle(compounds), Each(structures)], [compounds, structures])
    assert runtime == declared


def test_the_positional_split_it_replaces_really_was_wrong(tmp_path):
    """Kept as evidence rather than as a guard: this is what the old path produced for the case above."""
    compounds = _stream(tmp_path, "compounds", ["LIG1", "LIG2"])
    structures = _stream(tmp_path, "structures", ["seedA_1", "seedA_2"])
    parents, prov, axis_names, axes = Mock._resolve_from_source(
        [Bundle(compounds), Each(structures)])
    module = _pipe_mock()
    cfg = {"parent_ids": parents, "axes": axes, "provenance": prov,
           "axis_names": axis_names,
           "source_streams": [{"name": "compounds", "path": "unused"},
                              {"name": "structures", "path": "unused"}]}
    positional = module._rebuild_runtime_provenance(cfg, parents)
    assert positional["structures"] == ["LIG1", "LIG1"], "a ligand id in the structures column"
    assert positional != {name: prov[name] for name in prov}


@pytest.mark.parametrize("first_mode, second_mode", [
    (Each, Each),
    (Bundle, Each),
    (Each, Bundle),
    (Bundle, Bundle),
])
def test_provenance_agrees_for_every_axis_mode(tmp_path, first_mode, second_mode):
    compounds = _stream(tmp_path, "compounds", ["LIG1", "LIG2"])
    structures = _stream(tmp_path, "structures", ["seedA_1", "seedA_2", "seedA_3"])
    declared, runtime = _provenance_both_ways(
        tmp_path, [first_mode(compounds), second_mode(structures)], [compounds, structures])
    assert runtime == declared


# ── Mock as the framework's test surface ──────────────────────────────────────

def _axis_case(tmp_path, build, streams):
    parents, prov, axis_names, axes = Mock._resolve_from_source(build)
    source_streams = []
    for position, stream in enumerate(streams):
        json_path = tmp_path / f"axis_{position}.json"
        stream.save_json(str(json_path))
        source_streams.append({"name": stream.name, "path": str(json_path)})
    produced, _contributions = _pipe_mock()._resolve_parent_ids(
        {"parent_ids": parents, "axes": axes, "provenance": prov,
         "axis_names": axis_names, "source_streams": source_streams})
    return parents, produced


def test_the_static_pattern_agrees(tmp_path):
    """`Bundle(Each(a), b)` iterates over `a` while `b` rides along on every row. `get_mode` reports the outer wrapper as a bundle and hides the inner `Each`, so recording only a mode was too coarse: predict gave `['a1+b1','a2+b1']` and produce gave `['a1+a2+b1']`. The axis metadata records which streams iterate and which are static."""
    iterated = _stream(tmp_path, "sequences", ["a1", "a2"])
    static = _stream(tmp_path, "compounds", ["b1"])
    predicted, produced = _axis_case(
        tmp_path, Bundle(Each(iterated), static), [iterated, static])
    assert predicted == ["a1+b1", "a2+b1"]
    assert produced == predicted


def test_a_static_source_declared_first_lands_first(tmp_path):
    iterated = _stream(tmp_path, "sequences", ["a1", "a2"])
    static = _stream(tmp_path, "compounds", ["b1"])
    predicted, produced = _axis_case(
        tmp_path, Bundle(static, Each(iterated)), [static, iterated])
    assert predicted == ["b1+a1", "b1+a2"]
    assert produced == predicted


def test_a_pure_bundle_of_one_axis_is_one_id(tmp_path):
    stream = _stream(tmp_path, "sequences", ["a1", "a2"])
    predicted, produced = _axis_case(tmp_path, Bundle(stream), [stream])
    assert predicted == ["a1+a2"]
    assert produced == predicted


def test_two_streams_in_one_axis_may_share_a_name(tmp_path):
    """Streams are recorded by index rather than name for exactly this case: keyed by name, the static source would overwrite the iterated one."""
    iterated = _stream(tmp_path, "sequences", ["a1", "a2"])
    static = _stream(tmp_path, "sequences", ["b1"])
    predicted, produced = _axis_case(
        tmp_path, Bundle(Each(iterated), static), [iterated, static])
    assert produced == predicted
    assert predicted == ["a1+b1", "a2+b1"]
