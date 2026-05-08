"""Parameter coverage for biopipelines.posebusters.PoseBusters.

PoseBusters writes its config (mode, ligand, reference) to a JSON config
file, then invokes pipe_posebusters.py with --config. The helper grepping
both the bash and the JSON catches both surfaces.
"""

import pytest

from ._helpers import assert_substrings_in, read_all_emitted_artifacts


pytestmark = pytest.mark.tool_parameters


def _build(local_config, isolated_cwd, new_pipeline, **kwargs):
    from biopipelines.mock import Mock
    from biopipelines.posebusters import PoseBusters

    kwargs.setdefault("ligand", "ATP")
    pipeline = new_pipeline("pb_params")
    with pipeline:
        m = Mock(
            ids=["s1"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            map_table_strategy="config",
        )
        PoseBusters(structures=m.streams.structures, **kwargs)
        script_path = pipeline.save()
    return read_all_emitted_artifacts(script_path)


def test_ligand(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, ligand="HEM")
    assert "HEM" in content


def test_mode_dock(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, mode="dock")
    assert '"mode": "dock"' in content


def test_smoke_all_params(local_config, isolated_cwd, new_pipeline):
    content = _build(
        local_config, isolated_cwd, new_pipeline,
        ligand="LIG",
        mode="dock",
    )
    assert_substrings_in(content, ['"ligand_name": "LIG"', '"mode": "dock"'])
