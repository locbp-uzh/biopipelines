"""`autobox_ligand` must be judged by its type, not by how many ids it holds.

The type check and the type error hung off a length test, so a stream that was
a perfectly good `DataStream` but happened to be empty fell past every branch
and raised `autobox_ligand must be DataStream, StandardizedOutput, or str, got
<class 'DataStream'>` -- a message that contradicts itself and names the wrong
problem. An empty upstream is an ordinary outcome of a filter.
"""
from __future__ import annotations

import pytest

from biopipelines.datastream import DataStream


def _stream(ids):
    return DataStream(
        name="structures",
        ids=list(ids),
        files=[f"/tmp/{i}.pdb" for i in ids],
        format="pdb",
    )


def _gnina(**kwargs):
    from biopipelines.gnina import Gnina
    return Gnina(
        structures=_stream(["rec"]),
        compounds=_stream(["lig"]),
        **kwargs,
    )


def test_empty_autobox_ligand_stream_is_accepted(local_config, isolated_cwd):
    """Zero ids is a filtered-to-nothing upstream, not a type error."""
    tool = _gnina(autobox_ligand=_stream([]))

    assert tool.autobox_ligand_stream is not None
    # Nothing to point at, so there is no single global reference.
    assert tool.autobox_ligand_id is None
    assert tool.autobox_ligand_path is None


def test_single_id_autobox_ligand_becomes_a_global_reference(local_config, isolated_cwd):
    tool = _gnina(autobox_ligand=_stream(["only"]))

    assert tool.autobox_ligand_id == "only"


def test_multi_id_autobox_ligand_is_resolved_per_structure(local_config, isolated_cwd):
    """More than one id means the box is chosen by id at runtime, not up front."""
    tool = _gnina(autobox_ligand=_stream(["a", "b"]))

    assert tool.autobox_ligand_stream is not None
    assert tool.autobox_ligand_id is None


def test_string_autobox_ligand_is_taken_as_a_path(local_config, isolated_cwd):
    tool = _gnina(autobox_ligand="/data/ref.sdf")

    assert tool.autobox_ligand_path == "/data/ref.sdf"
    assert tool.autobox_ligand_stream is None


def test_a_genuinely_wrong_type_still_raises(local_config, isolated_cwd):
    """The error must survive the restructure -- and name the real type."""
    with pytest.raises(ValueError, match="autobox_ligand must be"):
        _gnina(autobox_ligand=42)


def test_omitting_autobox_ligand_leaves_every_field_unset(local_config, isolated_cwd):
    tool = _gnina()

    assert tool.autobox_ligand_stream is None
    assert tool.autobox_ligand_id is None
    assert tool.autobox_ligand_path is None
