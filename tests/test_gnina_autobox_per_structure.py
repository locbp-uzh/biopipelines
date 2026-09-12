"""A per-structure autobox stream must resolve one reference per receptor.

`autobox_ligand` accepts a DataStream, but the wrapper resolved only `ids[0]` and
baked that single path into the config for the whole step. With a stream holding
one reference per receptor that is wrong twice over: every receptor gets boxed on
the first structure's ligand, and — because the emitted command carried the
unexpanded `<id>.pdb` template — GNINA failed outright with

    Error: could not open ".../boxref/structures/<id>.pdb" for reading.

on every pair. The failure was loud here, but the silent version (all receptors
boxed on one ligand) is the dangerous one.
"""

import inspect
import pathlib

import pytest


def test_single_item_stream_still_resolves_at_config_time():
    """A stream of one is a shared reference and needs no runtime lookup."""
    src = inspect.getsource(
        __import__("biopipelines.gnina", fromlist=["x"]).Gnina._generate_script_resolve_autobox)
    assert "AUTOBOX_LIGAND_ID" in src, "the single-reference path must be kept"


def test_multi_item_stream_defers_to_runtime():
    src = inspect.getsource(
        __import__("biopipelines.gnina", fromlist=["x"]).Gnina._generate_script_resolve_autobox)
    assert "autobox_ligand_ds" in src, \
        "a per-structure stream must hand the datastream to the runtime, not a single path"
    assert "len(self.autobox_ligand_stream) != 1" in src, \
        "the two cases must be distinguished by stream length"


def test_autobox_id_not_pinned_to_first_for_multi_streams():
    """ids[0] as *the* reference is exactly the bug; it may only stand for a lone item."""
    src = inspect.getsource(
        __import__("biopipelines.gnina", fromlist=["x"]).Gnina.__init__)
    assert "len(self.autobox_ligand_stream) == 1" in src, \
        "autobox_ligand_id may only be set when the stream holds a single reference"


@pytest.fixture(scope="module")
def pipe_gnina():
    """Load the pipe script directly; pipe_scripts is not an importable package."""
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "pipe_gnina", pathlib.Path(__file__).resolve().parent.parent
        / "pipe_scripts" / "pipe_gnina.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_runtime_resolves_reference_per_receptor(pipe_gnina):
    """Behavioural, not textual: each receptor must get ITS OWN reference.

    This replaced a source-grep asserting the literal string
    `autobox_by_id.get(protein_id)`. Hoisting the lookup out of the loop into one
    provenance-aware pass improved the behaviour and broke the grep -- the failure
    mode of pinning text rather than behaviour.
    """
    refs = {"recA": "/refs/a.sdf", "recB": "/refs/b.sdf"}
    out = pipe_gnina.resolve_autobox_references(["recA", "recB"], refs)

    assert out["recA"] == "/refs/a.sdf"
    assert out["recB"] == "/refs/b.sdf", (
        "every receptor got the same reference -- the bug this feature fixed")


def test_receptor_without_a_reference_resolves_to_none(pipe_gnina):
    """A missing reference must be falsy so the caller falls back to the HETATM."""
    out = pipe_gnina.resolve_autobox_references(
        ["recA", "recZ"], {"recA": "/refs/a.sdf"})
    assert out["recA"] == "/refs/a.sdf"
    assert not out["recZ"]


def test_no_references_is_an_empty_map(pipe_gnina):
    assert pipe_gnina.resolve_autobox_references(["recA"], {}) == {}


def test_the_map_table_is_consulted_for_renamed_ids(pipe_gnina, tmp_path, monkeypatch):
    """An upstream rename leaves no string relation, so provenance must be passed.

    Panda_1 has no suffix relation to design_3; without map_table_paths the lookup
    finds nothing and every receptor silently falls back to its own crystal HETATM.
    """
    import json

    map_csv = tmp_path / "structures_map.csv"
    map_csv.write_text("id,file\nrecA,/x/a.pdb\n", encoding="utf-8")
    sj = tmp_path / "structures.json"
    sj.write_text(json.dumps({"name": "structures", "ids": ["recA"],
                              "files": ["/x/a.pdb"], "map_table": str(map_csv),
                              "format": "pdb"}), encoding="utf-8")

    seen = {}

    def spy(**kwargs):
        seen.update(kwargs)
        return {}

    monkeypatch.setattr(pipe_gnina, "get_mapped_ids", spy)
    pipe_gnina.resolve_autobox_references(["recA"], {"other": "/refs/o.sdf"},
                                          structures_json=str(sj))

    assert seen.get("map_table_paths") == [str(map_csv)], (
        "the structures map_table was not passed to get_mapped_ids, so a renamed "
        "id cannot resolve through provenance")
