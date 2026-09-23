"""A step whose filter kept nothing is not a step that failed.

Found on a real campaign. `007_Panda` filtered on `plddt > 80 and RMSD < 2.0`, the four designs
came back at 8.6-11.8 A, and all four were dropped — correctly, and recorded in `missing.csv`.
Panda logged "operations completed successfully" and the step was then marked FAILED, because
its `fasta` output is one shared multi-record file for the whole step and, with no ids left,
nothing wrote it.

`_filter_expected_missing` excuses a path by its *owner id*, and a shared artifact has none, so
it could never be excused that way. The rule this pins is simpler than trying to give shared
files an owner: when every declared id of a stream is excused, that stream has nothing left to
write, artifacts included.

The distinction that matters is the last two cases. An empty result must stay distinguishable
from a broken one — if only some ids were filtered, the survivors' files are still required.
"""

import json
import os
import sys

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "pipe_scripts"))

import pipe_check_completion  # noqa: E402


def build(tmp_path, ids, dropped, write_files=(), shared=None):
    """A step folder declaring `ids`, with `dropped` recorded in its missing table."""
    step = tmp_path / "007_Panda"
    (step / "structures").mkdir(parents=True)
    (step / "tables").mkdir(parents=True)

    for name in write_files:
        (step / "structures" / name).write_text("ATOM\n", encoding="utf-8")
    if shared:
        (step / "fasta").mkdir(exist_ok=True)
        (step / "fasta" / shared).write_text(">x\nSEQ\n", encoding="utf-8")

    rows = ["id,removed_by,kind,cause"]
    rows += [f"{i},007_Panda,filter,Filtered by: plddt > 80 and RMSD < 2.0" for i in dropped]
    (step / "tables" / "missing.csv").write_text("\n".join(rows) + "\n", encoding="utf-8")

    expected = {
        "output_folder": str(step),
        "structures": {
            "name": "structures", "format": "pdb", "ids": list(ids),
            # A `<id>` template, the way a real stream serializes it: that is what
            # gives each path an owner id, which is what excusal matches on.
            "files": [str(step / "structures" / "<id>.pdb")],
            "map_table": str(step / "structures" / "structures_map.csv"),
        },
        # One shared artifact for the whole step, the shape that could never be excused.
        "fasta": {
            "name": "fasta", "format": "fasta", "ids": list(ids),
            "files": [str(step / "fasta" / "sequences.fasta")],
            "map_table": str(step / "fasta" / "fasta_map.csv"),
        },
        "tables": {"missing": {"path": str(step / "tables" / "missing.csv")}},
    }
    return str(step), expected


def test_every_id_filtered_is_not_a_failure(tmp_path):
    ids = ["3kzy_A_1_1", "3kzy_A_2_1", "3kzy_A_3_1", "3kzy_A_4_1"]
    step, expected = build(tmp_path, ids, dropped=ids)
    ok, missing = pipe_check_completion.check_expected_outputs(expected, "Panda", step)
    assert ok, f"a fully-filtered step was reported as failed: {missing}"


def test_the_shared_artifact_is_what_used_to_fail(tmp_path):
    """Pin the specific shape: per-id files excused, the one shared FASTA not."""
    ids = ["a", "b"]
    step, expected = build(tmp_path, ids, dropped=ids)
    ok, missing = pipe_check_completion.check_expected_outputs(expected, "Panda", step)
    assert "fasta" not in missing and ok


def test_a_partial_filter_still_requires_the_survivors(tmp_path):
    """The case that must NOT be excused — otherwise a broken step looks like a filtered one."""
    ids = ["a", "b", "c"]
    step, expected = build(tmp_path, ids, dropped=["a"])
    ok, missing = pipe_check_completion.check_expected_outputs(expected, "Panda", step)
    assert not ok, "only one id was filtered; b and c still owe their outputs"
    assert "structures" in missing


def test_a_partial_filter_with_its_survivors_present_passes(tmp_path):
    ids = ["a", "b", "c"]
    step, expected = build(tmp_path, ids, dropped=["a"],
                           write_files=["b.pdb", "c.pdb"], shared="sequences.fasta")
    ok, missing = pipe_check_completion.check_expected_outputs(expected, "Panda", step)
    assert ok, missing


def test_nothing_filtered_and_nothing_written_is_still_a_failure(tmp_path):
    """With no drops recorded, absent outputs mean the tool did not do its job."""
    step, expected = build(tmp_path, ["a", "b"], dropped=[])
    ok, missing = pipe_check_completion.check_expected_outputs(expected, "Panda", step)
    assert not ok and missing
