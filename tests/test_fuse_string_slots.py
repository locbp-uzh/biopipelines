"""A bare string in `Fuse(sequences=[...])` is a sequence, promoted to an internal `Sequence`.

The documentation spelled the slots as strings (`Fuse(sequences=[N, mid, C])` with `N="GNH..."`) while the code
accepted only streams and tool outputs, so the documented call raised "Slot 0: expected DataStream or
StandardizedOutput, got <class 'str'>". Other tools already take the shorthand through `resolve_basic_input`
(`ligand="LIG"` is an internal `Ligand(codes="LIG")`), and a string slot now does the same.
"""

import json
import pathlib

import pytest

from biopipelines.fuse import Fuse
from biopipelines.mock import Mock


def _configure(new_pipeline, slots):
    pipeline = new_pipeline("fuse_strings")
    with pipeline:
        Fuse(sequences=slots(), linker="GS", linker_lengths=["2"])
        pipeline.save()
    (config_path,) = pathlib.Path.cwd().rglob("fuse_config.json")
    with open(config_path) as handle:
        return json.load(handle)


def test_two_string_slots_are_two_internal_sequences(local_config, isolated_cwd, new_pipeline):
    config = _configure(new_pipeline, lambda: ["MKTAYIAKQR", "QISFVKSHFSRQ"])
    assert len(config["slots"]) == 2
    assert all(slot["ids"] == ["seq1"] for slot in config["slots"])


def test_a_string_slot_can_stand_beside_a_tool_output(local_config, isolated_cwd, new_pipeline):
    def slots():
        parent = Mock(ids=["a1", "a2"], streams={"sequences": {"format": "csv", "values": ["AC", "DE"]}},
                      map_table_strategy="config")
        return [parent, "QISFVKSHFSRQ"]
    config = _configure(new_pipeline, slots)
    assert config["slots"][0]["ids"] == ["a1", "a2"]
    assert config["slots"][1]["ids"] == ["seq1"]


def test_anything_that_is_not_a_sequence_source_is_still_refused(local_config, isolated_cwd, new_pipeline):
    pipeline = new_pipeline("fuse_refused")
    with pipeline:
        Mock(ids=["a1"], streams={"sequences": {"format": "csv", "values": ["AC"]}}, map_table_strategy="config")
        with pytest.raises(ValueError, match="Slot 0"):
            Fuse(sequences=[42, "QISFVKSHFSRQ"], linker="GS", linker_lengths=["2"])
