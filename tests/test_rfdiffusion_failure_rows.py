# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""An input whose run crashed is a `failure` row of the step's `missing.csv`, with the error the run printed.

RFdiffusion, RFdiffusion2, RFdiffusion3 and RFdiffusionAllAtom run once per input structure and leave no trace of a run that crashed but the absence of its designs: no `missing.csv`, so the step stopped at FAILED with nothing naming the input or the reason. Each run's stderr is now captured next to its output, `pipe_update_structures_map.py` writes a `failure` row for every design a started run did not make (the cause is the real error line, not the hint printed after it), and `generate_missing_propagation(local_missing=...)` merges those rows with the upstream manifests into `tables/missing.csv`, as Gnina and the other tools do. The good ids and their parents are written as before; the script exits non-zero only when the runs that started made no file at all.

The shared script is driven on synthetic output folders, and the four wrappers' generated scripts are inspected.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "pipe_scripts" / "pipe_update_structures_map.py"

HYDRA_CRASH = """[2026-10-30 13:30:07][INFO] - Using contig: ['A1-50/10-20']
Error executing job with overrides: ["contigmap.contigs=['A1-50/10-20']"]
Traceback (most recent call last):
  File "run_inference.py", line 84, in main
    x_init, seq_init = sampler.sample_init()
  File "contigs.py", line 342, in get_idx0
    assert val in self.parsed_pdb["pdb_idx"], f"{val} is not in pdb file!"
AssertionError: ('A', 47) is not in pdb file!

Set the environment variable HYDRA_FULL_ERROR=1 for a complete stack trace.
"""

FOUNDRY_CRASH = """13:33:41 ERROR rfd3.inference.input_parsing: [rank: 0] Failed in validator: load_input
  Error: [component=A47] Residue A47 not found in atom array.
  Error type: ComponentValidationError
ValidationError: 1 validation error for DesignInputSpecification
  Value error, [component=A47] Residue A47 not found in atom array.
    For further information visit https://errors.pydantic.dev/2.14/v/value_error
"""


def _stream(tmp_path, ids, present):
    """An input stream of `ids`; only the `present` ones have a file, as after an upstream filter."""
    inputs = tmp_path / "inputs"
    inputs.mkdir(exist_ok=True)
    for pid in present:
        (inputs / f"{pid}.pdb").write_text("ATOM\n")
    path = tmp_path / "parents.json"
    path.write_text(json.dumps({"name": "structures", "ids": ids, "files": [str(inputs / "<id>.pdb")],
                                "map_table": "", "format": "pdb"}))
    return path


def _run(tmp_path, parents, made, stderr=None, key_prefix=None, suffix="<1..2>", extra=True):
    step = tmp_path / "005_RFdiffusion"
    structures = step / "structures"
    structures.mkdir(parents=True, exist_ok=True)
    for name in made:
        (structures / f"{name}.pdb").write_text("ATOM\n")
    errors = step / "_execution" / "stderr"
    errors.mkdir(parents=True, exist_ok=True)
    for key, text in (stderr or {}).items():
        (errors / f"{key}.stderr").write_text(text)
    cmd = [sys.executable, str(SCRIPT), "--structures-map", str(structures / "structures_map.csv"),
           "--output-folder", str(structures), "--parents-json", str(parents), "--suffix-pattern", suffix,
           "--parent-column", "structures.id"]
    if key_prefix:
        cmd += ["--key-prefix", key_prefix]
    local = step / "_execution" / "local_missing.csv"
    if extra:
        cmd += ["--local-missing-csv", str(local), "--stderr-dir", str(errors)]
    done = subprocess.run(cmd, capture_output=True, text=True)
    failures = pd.read_csv(local, dtype=str, keep_default_na=False) if local.exists() else None
    mapped = pd.read_csv(structures / "structures_map.csv", dtype=str, keep_default_na=False) \
        if (structures / "structures_map.csv").exists() else None
    return done, failures, mapped


class TestACrashedInputMidBatch:

    def test_it_is_a_failure_row_per_design_it_did_not_make_with_the_real_error_line(self, tmp_path):
        parents = _stream(tmp_path, ["a_1", "crash_2", "c_3"], present=["a_1", "crash_2", "c_3"])
        done, failures, _ = _run(tmp_path, parents, made=["a_1_1", "a_1_2", "c_3_1", "c_3_2"],
                                 stderr={"crash_2": HYDRA_CRASH})
        assert done.returncode == 0, done.stderr
        assert list(failures["id"]) == ["crash_2_1", "crash_2_2"]
        assert set(failures["kind"]) == {"failure"}
        assert set(failures["removed_by"]) == {"005_RFdiffusion"}
        assert set(failures["cause"]) == {"AssertionError: ('A', 47) is not in pdb file!"}

    def test_the_good_ids_and_their_parents_are_what_they_were(self, tmp_path):
        parents = _stream(tmp_path, ["a_1", "crash_2", "c_3"], present=["a_1", "crash_2", "c_3"])
        done, _, mapped = _run(tmp_path, parents, made=["a_1_1", "a_1_2", "c_3_1", "c_3_2"],
                               stderr={"crash_2": HYDRA_CRASH})
        assert done.returncode == 0, done.stderr
        assert dict(zip(mapped["id"], mapped["structures.id"])) == {
            "a_1_1": "a_1", "a_1_2": "a_1", "c_3_1": "c_3", "c_3_2": "c_3"}

    def test_a_run_that_made_some_designs_is_a_failure_for_the_rest(self, tmp_path):
        parents = _stream(tmp_path, ["a_1", "half_2"], present=["a_1", "half_2"])
        done, failures, _ = _run(tmp_path, parents, made=["a_1_1", "a_1_2", "half_2_1"],
                                 stderr={"half_2": HYDRA_CRASH})
        assert done.returncode == 0, done.stderr
        assert list(failures["id"]) == ["half_2_2"]

    def test_the_foundry_validator_line_is_the_cause_not_the_exception_that_wraps_it(self, tmp_path):
        parents = _stream(tmp_path, ["a_1", "crash_2"], present=["a_1", "crash_2"])
        done, failures, mapped = _run(tmp_path, parents, made=["mot_a_1_1", "mot_a_1_2"],
                                      stderr={"mot_crash_2": FOUNDRY_CRASH}, key_prefix="mot")
        assert done.returncode == 0, done.stderr
        assert list(failures["id"]) == ["mot_crash_2_1", "mot_crash_2_2"]
        assert set(failures["cause"]) == {"[component=A47] Residue A47 not found in atom array."}
        assert dict(zip(mapped["id"], mapped["structures.id"])) == {"mot_a_1_1": "a_1", "mot_a_1_2": "a_1"}

    def test_a_run_with_no_captured_error_still_gets_a_row(self, tmp_path):
        parents = _stream(tmp_path, ["a_1", "crash_2"], present=["a_1", "crash_2"])
        done, failures, _ = _run(tmp_path, parents, made=["a_1_1", "a_1_2"])
        assert done.returncode == 0, done.stderr
        assert list(failures["id"]) == ["crash_2_1", "crash_2_2"]
        assert all(failures["cause"])

    def test_a_parent_an_upstream_filter_dropped_never_ran_and_gets_no_row(self, tmp_path):
        parents = _stream(tmp_path, ["a_1", "dropped_2"], present=["a_1"])
        done, failures, _ = _run(tmp_path, parents, made=["a_1_1", "a_1_2"])
        assert done.returncode == 0, done.stderr
        assert len(failures) == 0


class TestNothingProduced:

    def test_the_script_exits_non_zero_when_every_run_that_started_made_nothing(self, tmp_path):
        parents = _stream(tmp_path, ["a_1", "b_2"], present=["a_1", "b_2"])
        done, failures, _ = _run(tmp_path, parents, made=[], stderr={"a_1": HYDRA_CRASH, "b_2": HYDRA_CRASH})
        assert done.returncode != 0
        assert list(failures["id"]) == ["a_1_1", "a_1_2", "b_2_1", "b_2_2"]

    def test_it_exits_zero_when_an_upstream_filter_left_no_run_to_start(self, tmp_path):
        parents = _stream(tmp_path, ["a_1", "b_2"], present=[])
        done, failures, _ = _run(tmp_path, parents, made=[])
        assert done.returncode == 0, done.stderr
        assert len(failures) == 0

    def test_a_caller_without_the_flag_is_left_alone(self, tmp_path):
        parents = _stream(tmp_path, ["a_1", "crash_2"], present=["a_1", "crash_2"])
        done, failures, mapped = _run(tmp_path, parents, made=["a_1_1"], extra=False)
        assert done.returncode == 0, done.stderr
        assert failures is None
        assert list(mapped["id"]) == ["a_1_1"]


def _body(new_pipeline, build, suffix):
    import pathlib

    pipeline = new_pipeline("rfd_failure_rows")
    with pipeline:
        build()
        script_path = pipeline.save()
    return next(pathlib.Path(script_path).parent.glob(f"*_{suffix}.sh")).read_text(encoding="utf-8", errors="replace")


def _structures(tmp_path, ids=("p_1", "q_2")):
    from biopipelines.datastream import DataStream, create_map_table

    folder = tmp_path / "scaffolds"
    folder.mkdir(exist_ok=True)
    for pid in ids:
        (folder / f"{pid}.pdb").write_text("END\n")
    template = str(folder / "<id>.pdb")
    table = tmp_path / "structures_map.csv"
    create_map_table(str(table), ids=list(ids), files=[template])
    return DataStream(name="structures", ids=list(ids), files=[template], map_table=str(table), format="pdb")


def _loop(body):
    start = body.index("for STRUCT_ID in")
    return body[start:body.index("\ndone", start)]


def _checks(body):
    loop = _loop(body)
    assert 'tee "' in loop and '/$STRUCT_ID.stderr" >&2' in loop, "each run's stderr is captured next to the output"
    update = body[body.index("pipe_update_structures_map.py"):]
    update = update[:update.index("\n")]
    assert "--local-missing-csv" in update and "--stderr-dir" in update
    merge = body[body.index("pipe_propagate_missing.py"):]
    assert "--local-missing" in merge[:merge.index("\n\n")] or "--local-missing" in merge
    assert body.index("pipe_update_structures_map.py") < body.index("BP_MAIN_RC=$?") < body.index("pipe_propagate_missing.py")


class TestTheWrappers:

    def test_rfdiffusion(self, local_config, isolated_cwd, new_pipeline, tmp_path):
        from biopipelines.rfdiffusion import RFdiffusion

        body = _body(new_pipeline,
                     lambda: RFdiffusion(contigs="A1-5/10-20", pdb=_structures(tmp_path), num_designs=2),
                     "RFdiffusion")
        _checks(body)

    def test_rfdiffusion2(self, local_config, isolated_cwd, new_pipeline, tmp_path):
        from biopipelines.rfdiffusion2 import RFdiffusion2

        body = _body(new_pipeline,
                     lambda: RFdiffusion2(ligand="LIG", pdb=_structures(tmp_path), contigs="10", num_designs=2),
                     "RFdiffusion2")
        _checks(body)

    def test_rfdiffusion_allatom(self, local_config, isolated_cwd, new_pipeline, tmp_path):
        from biopipelines.rfdiffusion_allatom import RFdiffusionAllAtom

        body = _body(new_pipeline,
                     lambda: RFdiffusionAllAtom(ligand="LIG", pdb=_structures(tmp_path), contigs="A1-5,10-20",
                                                num_designs=2),
                     "RFdiffusionAllAtom")
        _checks(body)

    @pytest.mark.parametrize("tool", ["RFdiffusion", "RFdiffusion2", "RFdiffusionAllAtom"])
    def test_the_manifest_is_declared_without_an_upstream_filter(self, tool, local_config, isolated_cwd,
                                                                  new_pipeline, tmp_path):
        from biopipelines.rfdiffusion import RFdiffusion
        from biopipelines.rfdiffusion2 import RFdiffusion2
        from biopipelines.rfdiffusion_allatom import RFdiffusionAllAtom

        build = {
            "RFdiffusion": lambda: RFdiffusion(contigs="A1-5/10-20", pdb=_structures(tmp_path)),
            "RFdiffusion2": lambda: RFdiffusion2(ligand="LIG", pdb=_structures(tmp_path), contigs="10"),
            "RFdiffusionAllAtom": lambda: RFdiffusionAllAtom(ligand="LIG", pdb=_structures(tmp_path),
                                                             contigs="A1-5,10-20"),
        }[tool]
        pipeline = new_pipeline("rfd_declares_missing")
        with pipeline:
            built = build()
        assert "missing" in built.tool.get_output_files()["tables"]
