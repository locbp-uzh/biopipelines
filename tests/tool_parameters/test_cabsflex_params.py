"""Parameter coverage for biopipelines.cabsflex.CABSflex."""

import pytest

from ._helpers import assert_substrings_in, read_pipeline_sh


pytestmark = pytest.mark.tool_parameters


def _build(local_config, isolated_cwd, new_pipeline, **kwargs):
    from biopipelines.mock import Mock
    from biopipelines.cabsflex import CABSflex

    pipeline = new_pipeline("cabsflex_params")
    with pipeline:
        m = Mock(
            ids=["s1"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            map_table_strategy="config",
        )
        CABSflex(structures=m.streams.structures, **kwargs)
        script_path = pipeline.save()
    return read_pipeline_sh(script_path)


def test_num_models(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, num_models=20)
    assert "-k 20" in content


def test_mc_cycles(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, mc_cycles=100)
    assert "-y 100" in content


def test_mc_steps(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, mc_steps=75)
    assert "-s 75" in content


def test_mc_annealing(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, mc_annealing=15)
    assert "-a 15" in content


def test_temperature(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, temperature="2.0")
    assert "-t 2.0" in content


def test_filtering_count(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, filtering_count=500)
    assert "-n 500" in content


def test_aa_rebuild(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, aa_rebuild=True)
    assert " -A " in content or content.endswith(" -A")


def test_pdb_output(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, pdb_output="M")
    assert "-o M" in content


def test_restraints(local_config, isolated_cwd, new_pipeline):
    content = _build(
        local_config, isolated_cwd, new_pipeline,
        restraints="ss1", restraints_gap=4, restraints_min=4.0, restraints_max=9.0,
    )
    assert "-g ss1 4 4.0 9.0" in content


def test_smoke_all_params(local_config, isolated_cwd, new_pipeline):
    content = _build(
        local_config, isolated_cwd, new_pipeline,
        num_models=5,
        mc_cycles=20,
        mc_steps=30,
        mc_annealing=10,
        temperature="1.5",
        filtering_count=200,
        aa_rebuild=True,
        pdb_output="M",
    )
    assert_substrings_in(content, [
        "-k 5",
        "-y 20",
        "-s 30",
        "-a 10",
        "-t 1.5",
        "-n 200",
        "-A",
        "-o M",
    ])
