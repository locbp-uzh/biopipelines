"""Gnina and Vina: two (protein, ligand) pairs that would share one pair id are refused.

A `+` inside an id lets two different pairs join to the same string: protein `p1` with ligand `b+c` and protein
`p1+b` with ligand `c` are both `p1+b+c`. The two shared one docking scratch file, one `<id>_best.pdb` and one row
of every table, so one pair vanished from the outputs with no `missing` row and the step reported COMPLETED. The
wrapper refuses the inputs when the pipeline is configured, and the pipe script refuses them again before it docks
anything, because an upstream stream with lazy ids is only known at run time.
"""

import sys

import pytest

from tests._helpers import declared_stream


@pytest.mark.parametrize("cls", ["Gnina", "Vina"])
def test_the_wrapper_refuses_two_pairs_with_one_id(local_config, isolated_cwd, new_pipeline, cls):
    with pytest.raises(ValueError, match=r"p1\+b\+c"):
        declared_stream(new_pipeline, ["p1", "p1+b"], ["b+c", "c"], cls=cls)


@pytest.mark.parametrize("cls", ["Gnina", "Vina"])
def test_distinct_pair_ids_with_a_plus_inside_an_id_are_kept(local_config, isolated_cwd, new_pipeline, cls):
    declared, _, _ = declared_stream(new_pipeline, ["p1", "p2"], ["b+c", "c"], cls=cls)

    assert list(declared.ids) == ["p1+b+c", "p1+c", "p2+b+c", "p2+c"]


def test_the_pipe_script_refuses_two_pairs_with_one_id_before_it_docks(pipe_gnina, monkeypatch, tmp_path):
    class Stream:
        def __init__(self, ids):
            self.ids_expanded = list(ids)

    streams = {"structures.json": Stream(["p1", "p1+b"]), "compounds.json": Stream(["b+c", "c"])}
    monkeypatch.setattr(pipe_gnina, "load_datastream", lambda path: streams[path])
    monkeypatch.setattr(pipe_gnina, "iterate_files", lambda stream: [(i, f"{i}.pdb") for i in stream.ids_expanded])

    def docked(*args, **kwargs):
        raise AssertionError("the pipe script prepared a protein before it checked the pair ids")

    monkeypatch.setattr(pipe_gnina, "prepare_proteins", docked)
    config = tmp_path / "gnina_config.json"
    config.write_text('{"mode": "docking", "backend": "gnina", "structures_json": "structures.json", '
                      '"compounds_json": "compounds.json", "upstream_missing": null}')
    monkeypatch.setattr(sys, "argv", ["pipe_gnina.py", str(config)])

    with pytest.raises(ValueError, match=r"p1\+b\+c"):
        pipe_gnina.main()
