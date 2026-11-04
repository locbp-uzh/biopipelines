"""A pooled Panda promises the streams of all its pools, but a winner can come from a pool that has one fewer.

Found by running the iterative Gnina example on the cluster. `best_sequence` pools a Boltz2 output
(sequences, structures, pae, msas, compounds) with a ReMap output (sequences only). The winner of the
sort came from ReMap, so no structure, PAE or MSA of it exists anywhere, and Panda was marked FAILED
at every cycle for files its source never had, though its table and its sequences were right.

The winner's own pool is known at run time, so Panda records which streams that pool does not carry
at all, and the completion check does not demand them. A stream the pool does carry and is missing a
file of stays a failure, which is what the last case pins.
"""

import os

import pandas as pd

from pipe_scripts.pipe_panda import PANDA_SOURCE, streams_not_carried
from tests._helpers import sys_path_prepended

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
with sys_path_prepended(os.path.join(REPO_ROOT, "pipe_scripts")):
    import pipe_check_completion

BOLTZ_LIKE = {"sequences": {"a": "a.csv"}, "structures": {"a": "a.pdb"}, "pae": {"a": "a.npz"}}
REMAP_LIKE = {"sequences": {"b": "b.csv"}}


def winner_rows(source):
    return pd.DataFrame({"id": ["best_1"], "original_id": ["b"], PANDA_SOURCE: [source]})


class TestWhatIsNotCarried:

    def test_a_winner_from_the_pool_with_fewer_streams_lacks_the_others(self):
        absent = streams_not_carried(winner_rows(1), [BOLTZ_LIKE, REMAP_LIKE])
        assert absent == {"structures": ["best_1"], "pae": ["best_1"]}

    def test_a_winner_from_the_pool_with_every_stream_lacks_nothing(self):
        assert streams_not_carried(winner_rows(0), [BOLTZ_LIKE, REMAP_LIKE]) == {}


def build_step(tmp_path, carried=None, write_structure=False):
    """A step declaring `best_1` in two streams, the structure of which its source never had."""
    step = tmp_path / "014_Panda"
    (step / "sequences").mkdir(parents=True)
    (step / "structures").mkdir()
    (step / "tables").mkdir()
    (step / "_execution").mkdir()
    (step / "sequences" / "best_1.csv").write_text("id,sequence\nbest_1,MK\n", encoding="utf-8")
    if write_structure:
        (step / "structures" / "best_1.pdb").write_text("ATOM\n", encoding="utf-8")
    (step / "tables" / "missing.csv").write_text("id,removed_by,kind,cause\n", encoding="utf-8")
    if carried is not None:
        import json
        (step / "_execution" / "not_carried.json").write_text(json.dumps(carried), encoding="utf-8")
    expected = {
        "output_folder": str(step),
        "sequences": {"name": "sequences", "format": "csv", "ids": ["best_1"],
                      "files": [str(step / "sequences" / "<id>.csv")]},
        "structures": {"name": "structures", "format": "pdb", "ids": ["best_1"],
                       "files": [str(step / "structures" / "<id>.pdb")]},
        "tables": {"missing": {"path": str(step / "tables" / "missing.csv")}},
    }
    return str(step), expected


class TestTheCompletionCheck:

    def test_a_stream_the_source_pool_never_carried_is_not_demanded(self, tmp_path):
        step, expected = build_step(tmp_path, carried={"structures": ["best_1"]})
        ok, missing = pipe_check_completion.check_expected_outputs(expected, "Panda", step)
        assert ok, f"a file its source never had was demanded: {missing}"

    def test_without_the_record_the_missing_structure_is_still_a_failure(self, tmp_path):
        step, expected = build_step(tmp_path)
        ok, missing = pipe_check_completion.check_expected_outputs(expected, "Panda", step)
        assert not ok and "structures" in missing

    def test_a_stream_the_pool_carried_is_still_demanded_for_the_other_ids(self, tmp_path):
        step, expected = build_step(tmp_path, carried={"structures": ["best_2"]})
        ok, missing = pipe_check_completion.check_expected_outputs(expected, "Panda", step)
        assert not ok and "structures" in missing, "the record excuses only the ids it names"

    def test_a_file_that_exists_is_fine_whatever_the_record_says(self, tmp_path):
        step, expected = build_step(tmp_path, carried={"structures": ["best_1"]}, write_structure=True)
        ok, _ = pipe_check_completion.check_expected_outputs(expected, "Panda", step)
        assert ok
