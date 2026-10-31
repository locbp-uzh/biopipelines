# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Config-time coverage for the two LigandMPNN options that map onto upstream flags.

``save_stats`` passes ``--save_stats 1`` (upstream then writes ``execution/stats/<input file name>.pt``) and ``use_atom_context=False`` passes ``--ligand_mpnn_use_atom_context 0``.
Both default to upstream's behavior and are emitted only when they change it, so a pipeline that sets neither is unchanged.
Same pattern as test_ligand_mpnn_params.py: a Mock input, then grep the emitted scripts; no GPU.
"""

import os

import pytest

from ._helpers import assert_substrings_in, read_all_emitted_artifacts

pytestmark = pytest.mark.tool_parameters


def _build_lmpnn(new_pipeline, **kw):
    """Returns the LigandMPNN step as the pipeline holds it and everything the run emitted."""
    from biopipelines.ligand import Ligand
    from biopipelines.ligand_mpnn import LigandMPNN
    from biopipelines.mock import Mock

    pipeline = new_pipeline("lmpnn_newflags")
    with pipeline:
        m = Mock(
            ids=["s1"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            map_table_strategy="config",
        )
        LigandMPNN(structures=m.streams.structures, ligand=Ligand(code="ATP"), **kw)
        script_path = pipeline.save()
    (tool,) = [t for t in pipeline.tools if t.TOOL_NAME == "LigandMPNN"]
    return tool, read_all_emitted_artifacts(script_path)


def test_save_stats_emits_flag(local_config, isolated_cwd, new_pipeline):
    _, content = _build_lmpnn(new_pipeline, save_stats=True)
    assert_substrings_in(content, ["--save_stats 1"], label="save_stats")


def test_save_stats_default_off(local_config, isolated_cwd, new_pipeline):
    _, content = _build_lmpnn(new_pipeline)
    assert "--save_stats" not in content and '"--save_stats"' not in content


def test_use_atom_context_false_emits_flag(local_config, isolated_cwd, new_pipeline):
    _, content = _build_lmpnn(new_pipeline, use_atom_context=False)
    assert_substrings_in(content, ["--ligand_mpnn_use_atom_context 0"], label="use_atom_context")


def test_use_atom_context_default_true_absent(local_config, isolated_cwd, new_pipeline):
    # True is upstream's default; the flag is emitted only to turn context OFF.
    _, content = _build_lmpnn(new_pipeline)
    assert "--ligand_mpnn_use_atom_context" not in content \
        and '"--ligand_mpnn_use_atom_context"' not in content


def test_the_stats_land_in_the_execution_folder(local_config, isolated_cwd, new_pipeline):
    """Upstream creates <out_folder>/stats itself, and `stats_folder` is where the docstring tells the user to read it."""
    tool, _ = _build_lmpnn(new_pipeline, save_stats=True)
    assert os.path.normpath(tool.stats_folder) == os.path.normpath(os.path.join(tool.execution_folder, "stats"))


def test_both_options_are_recorded_in_the_serialized_configuration(local_config, isolated_cwd, new_pipeline):
    tool, _ = _build_lmpnn(new_pipeline, save_stats=True, use_atom_context=False)
    params = tool.to_dict()["lmpnn_params"]
    assert params["save_stats"] is True and params["use_atom_context"] is False
    tool, _ = _build_lmpnn(new_pipeline)
    params = tool.to_dict()["lmpnn_params"]
    assert params["save_stats"] is False and params["use_atom_context"] is True
