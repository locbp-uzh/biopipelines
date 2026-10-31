"""Generated shell scripts are UTF-8 with LF endings on every host, so a Windows-generated pipeline does not carry cp1252 bytes or CRLF."""

import pathlib

import pytest


@pytest.mark.usefixtures("local_config", "isolated_cwd")
def test_every_generated_script_is_utf8_with_lf(new_pipeline):
    from biopipelines.rfdiffusion3 import RFdiffusion3

    pipeline = new_pipeline("script_encoding")
    with pipeline:
        RFdiffusion3(contigs="A50-100,80-100", num_designs=1)
        script_path = pipeline.save()

    scripts = sorted(pathlib.Path(script_path).parent.glob("*.sh"))
    assert scripts, "no generated shell script found"
    for script in scripts:
        data = script.read_bytes()
        data.decode("utf-8")
        assert b"\r" not in data, f"{script.name} has CRLF line endings"


@pytest.mark.usefixtures("slurm_local_config", "isolated_cwd")
@pytest.mark.parametrize("batches", [1, 2], ids=["job", "job_batchN"])
def test_the_scheduler_job_scripts_are_utf8_with_lf(batches):
    """`save()` does not write the job scripts; `generate_job_scripts()` does, for one batch and for several."""
    from biopipelines.mock import Mock
    from biopipelines.pipeline import Pipeline, Resources

    pipeline = Pipeline(project="TestSuite", job=f"job_encoding_{batches}", description="encoding",
                        on_the_fly=False, local_output=True, config="slurm_local")
    with pipeline:
        for index in range(batches):
            Resources(gpu="none", time="1:00:00")
            Mock(ids=[f"s{index}"], streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
        pipeline.save()
        pipeline.generate_job_scripts()

    job_scripts = [pathlib.Path(pipeline.job_script)]
    job_scripts += sorted(pathlib.Path(pipeline.folders["runtime"]).glob(f"*{pathlib.Path(pipeline.job_script).suffix}"))
    for script in dict.fromkeys(job_scripts):
        data = script.read_bytes()
        data.decode("utf-8")
        assert b"\r" not in data, f"{script.name} has CRLF line endings"
