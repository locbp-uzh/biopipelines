# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""MSA: a listed alignment whose file is gone is a `failure` row of this step's `missing.csv`.

`pipe_msa.py` printed "MSA file not found" and moved on, so an id the producer's map listed whose file had since gone (a cleaned folder, a moved run) was declared by the step and absent from its outputs with no cause on record: the step was FAILED and the log line was the only trace. A blank `file` cell crashed the script instead, because a NaN reaches `os.path.exists`. Both are now `failure` rows of `execution/local_missing.csv`, which the propagation merges with the upstream manifests, as the other converters do.

Not covered, because it needs a cluster: a conversion of real ColabFold output.
"""

import json
import os
import pathlib
import subprocess
import sys

import pandas as pd

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PIPE = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_msa.py")
PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"


def _a3m(folder, name):
    path = folder / f"{name}.a3m"
    path.write_text(f"#{len(PROTEIN)}\t1\n>101\n{PROTEIN}\n>UniRef100_X\n{PROTEIN}\n")
    return str(path)


def _run(tmp_path, rows):
    step = tmp_path / "007_MSA"
    table = tmp_path / "in.csv"
    pd.DataFrame(rows).to_csv(table, index=False)
    config = tmp_path / "config.json"
    config.write_text(json.dumps({
        "input_msa_table": str(table), "convert": "csv", "output_folder": str(step / "msas"),
        "output_msas_csv": str(step / "msas" / "msas_map.csv"),
        "local_missing_csv": str(step / "_execution" / "local_missing.csv")}))
    done = subprocess.run([sys.executable, PIPE, "--config", str(config)], capture_output=True, text=True)
    return done, step


class TestAListedAlignmentThatIsGone:

    def test_each_one_is_a_failure_row_with_its_cause(self, tmp_path):
        present = _a3m(tmp_path, "p1")
        gone = str(tmp_path / "cleaned" / "p2.a3m")
        done, step = _run(tmp_path, [{"id": "p1", "file": present}, {"id": "p2", "file": gone},
                                     {"id": "p3", "file": ""}])
        assert done.returncode == 0, done.stdout + done.stderr
        converted = pd.read_csv(step / "msas" / "msas_map.csv")
        assert converted["id"].tolist() == ["p1"]
        local = pd.read_csv(step / "_execution" / "local_missing.csv", dtype=str, keep_default_na=False)
        assert local["id"].tolist() == ["p2", "p3"]
        assert set(local["kind"]) == {"failure"} and set(local["removed_by"]) == {"007_MSA"}
        assert gone in local["cause"][0]
        assert "no file" in local["cause"][1]

    def test_a_complete_table_leaves_no_row(self, tmp_path):
        done, step = _run(tmp_path, [{"id": "p1", "file": _a3m(tmp_path, "p1")}])
        assert done.returncode == 0, done.stdout + done.stderr
        local = pd.read_csv(step / "_execution" / "local_missing.csv")
        assert list(local.columns) == ["id", "removed_by", "kind", "cause"] and local.empty


class TestTheWrapper:

    @staticmethod
    def _step(new_pipeline, tmp_path):
        from biopipelines.mock import Mock
        from biopipelines.msa import MSA

        pipeline = new_pipeline("msa_missing")
        with pipeline:
            source = Mock(ids=["p1"], streams={"msas": {"format": "a3m", "file": "<id>.a3m"}},
                          tables={"msas": {"columns": ["id", "file"]},
                                  "missing": {"columns": ["id", "removed_by", "kind", "cause"]}},
                          map_table_strategy="config")
            built = MSA(source, convert="csv")
            script_path = pipeline.save()
        body = next(pathlib.Path(script_path).parent.glob("*_MSA.sh")).read_text(encoding="utf-8", errors="replace")
        return built, body

    def test_the_local_file_is_handed_to_the_pipe_and_merged_with_the_upstream_manifest(
            self, local_config, isolated_cwd, new_pipeline, tmp_path):
        built, body = self._step(new_pipeline, tmp_path)
        config = json.loads(open(built.producer.config_json).read())
        assert config["local_missing_csv"] == built.producer.local_missing_csv
        merge = body[body.index("pipe_propagate_missing.py"):]
        assert f'--local-missing "{built.producer.local_missing_csv}"' in merge
        assert body.index("pipe_msa.py") < body.index("BP_MAIN_RC=$?") < body.index("pipe_propagate_missing.py")

    def test_the_missing_table_is_declared(self, local_config, isolated_cwd, new_pipeline, tmp_path):
        built, _ = self._step(new_pipeline, tmp_path)
        assert built.tables.missing.info.columns == ["id", "removed_by", "kind", "cause"]
