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
