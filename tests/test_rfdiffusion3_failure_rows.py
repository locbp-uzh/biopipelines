# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""RFdiffusion3: an input whose run crashed is a `failure` row of `missing.csv`, and it costs only itself.

`rfd3 design` takes every input in one call and stops at the first it cannot parse, so one bad input left the inputs after it without designs, and nothing named it. The inputs the batch produced nothing for now run one by one (`pipe_rfdiffusion3_pending.py` names them and cuts each into a JSON of its own), each with its stderr captured; the shared update-map script turns what is still absent into `failure` rows whose cause is the error line of that run. A batch that ran clean pays nothing: there is no pending input and the loop does not start.

The pending-input helper is driven on synthetic raw folders and the generated script is inspected.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
PENDING = REPO / "pipe_scripts" / "pipe_rfdiffusion3_pending.py"

CONFIG = {
    "a": {"contig": "30-40,/0,A1-50", "input": "a.pdb"},
    "a_1": {"contig": "30-40,/0,A1-50", "input": "a_1.pdb"},
    "crash_2": {"contig": "30-40,/0,A1-50", "input": "crash_2.pdb"},
    "late_3": {"contig": "30-40,/0,A1-50", "input": "late_3.pdb"},
}


def _pending(tmp_path, made):
    raw = tmp_path / "raw"
    raw.mkdir()
    for name in made:
        (raw / name).write_text("x")
    config = tmp_path / "inputs.json"
    config.write_text(json.dumps(CONFIG))
    done = subprocess.run([sys.executable, str(PENDING), "--json", str(config), "--raw-folder", str(raw)],
                          capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    return done.stdout.split()


class TestThePendingInputs:

    def test_an_input_the_batch_never_reached_is_pending_with_the_one_that_crashed(self, tmp_path):
        made = ["a_0_model_0.cif.gz", "a_0_model_0.json", "a_1_0_model_0.cif.gz", "a_1_0_model_0.json"]
        assert _pending(tmp_path, made) == ["crash_2", "late_3"]

    def test_a_key_that_continues_another_keys_name_does_not_stand_in_for_it(self, tmp_path):
        """`a_1_0_model_0` is the first design of key `a_1`, not the second batch of key `a`."""
        assert _pending(tmp_path, ["a_1_0_model_0.cif.gz"]) == ["a", "crash_2", "late_3"]

    def test_a_batch_that_ran_clean_has_nothing_pending(self, tmp_path):
        made = [f"{key}_0_model_0.cif.gz" for key in CONFIG]
        assert _pending(tmp_path, made) == []

    def test_one_entry_is_cut_into_a_json_of_its_own(self, tmp_path):
        config = tmp_path / "inputs.json"
        config.write_text(json.dumps(CONFIG))
        single = tmp_path / "single.json"
        done = subprocess.run([sys.executable, str(PENDING), "--json", str(config), "--key", "crash_2",
                               "--single-json", str(single)], capture_output=True, text=True)
        assert done.returncode == 0, done.stderr
        assert json.loads(single.read_text()) == {"crash_2": CONFIG["crash_2"]}


class TestTheScript:

    def _body(self, new_pipeline, tmp_path, **kwargs):
        import pathlib
        from biopipelines.datastream import DataStream, create_map_table
        from biopipelines.rfdiffusion3 import RFdiffusion3

        folder = tmp_path / "scaffolds"
        folder.mkdir()
        for pid in ("p_1", "q_2"):
            (folder / f"{pid}.pdb").write_text("END\n")
        template = str(folder / "<id>.pdb")
        table = tmp_path / "structures_map.csv"
        create_map_table(str(table), ids=["p_1", "q_2"], files=[template])
        stream = DataStream(name="structures", ids=["p_1", "q_2"], files=[template], map_table=str(table),
                            format="pdb")
        pipeline = new_pipeline("rfd3_failure_rows")
        with pipeline:
            built = RFdiffusion3(pdb=stream, contig="30-40,/0,A1-5", num_designs=2, num_models=1, **kwargs)
            script_path = pipeline.save()
        body = next(pathlib.Path(script_path).parent.glob("*_RFdiffusion3.sh")).read_text(
            encoding="utf-8", errors="replace")
        return built, body

    def test_a_crashed_input_is_rerun_alone_with_its_stderr_captured(self, local_config, isolated_cwd, new_pipeline,
                                                                     tmp_path):
        _, body = self._body(new_pipeline, tmp_path)
        assert "run_rfd3()" in body
        retry = body[body.index("for KEY in"):]
        retry = retry[:retry.index("\ndone")]
        assert "pipe_rfdiffusion3_pending.py" in retry
        assert 'tee "' in retry and '/$KEY.stderr" >&2' in retry
        assert body.index("run_rfd3 ") < body.index("for KEY in")
        assert "RERUN_FAILS" in retry and "break" in retry, "a failing environment must not reload the model for every input"

    def test_the_rerun_uses_the_command_of_the_batch(self, local_config, isolated_cwd, new_pipeline, tmp_path):
        _, body = self._body(new_pipeline, tmp_path, seed=3)
        function = body[body.index("run_rfd3()"):]
        function = function[:function.index("\n}")]
        assert "rfd3 design" in function and 'inputs="$1"' in function and "seed=3" in function
        assert body.count("rfd3 design") == 1

    @pytest.mark.parametrize("prefix", [None, "mot"])
    def test_the_failure_rows_are_written_and_merged_after_the_map(self, prefix, local_config, isolated_cwd,
                                                                    new_pipeline, tmp_path):
        extra = {"prefix": prefix} if prefix else {}
        built, body = self._body(new_pipeline, tmp_path, **extra)
        update = body[body.index("pipe_update_structures_map.py"):]
        update = update[:update.index("\n")]
        assert "--local-missing-csv" in update and "--stderr-dir" in update
        assert ('--key-prefix "mot"' in update) == bool(prefix)
        assert body.index("pipe_update_structures_map.py") < body.index("BP_MAIN_RC=$?") < body.index(
            "pipe_propagate_missing.py")
        assert "missing" in built.tool.get_output_files()["tables"]
