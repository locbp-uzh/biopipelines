"""Parameter coverage for biopipelines.gnina.Gnina.

Gnina serialises most kwargs into a JSON config file at config-time
(``configuration/gnina_config.json``). The helper concatenates that JSON
with pipeline.sh so a substring grep covers both.
"""

import pytest

from ._helpers import assert_substrings_in, read_all_emitted_artifacts


pytestmark = pytest.mark.tool_parameters


def _build_with_outputs(local_config, isolated_cwd, new_pipeline, **gnina_kwargs):
    """Like _build but returns (content, StandardizedOutput) for output assertions."""
    from biopipelines.mock import Mock
    from biopipelines.gnina import Gnina

    pipeline = new_pipeline("gnina_params")
    with pipeline:
        struct = Mock(
            ids=["s1"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            map_table_strategy="config",
        )
        comp = Mock(
            ids=["c1"],
            streams={"compounds": {"format": "sdf", "file": "<id>.sdf"}},
            map_table_strategy="config",
        )
        out = Gnina(
            structures=struct.streams.structures,
            compounds=comp.streams.compounds,
            **gnina_kwargs,
        )
        script_path = pipeline.save()
    return read_all_emitted_artifacts(script_path), out


def _build(local_config, isolated_cwd, new_pipeline, **gnina_kwargs):
    from biopipelines.mock import Mock
    from biopipelines.gnina import Gnina

    pipeline = new_pipeline("gnina_params")
    with pipeline:
        struct = Mock(
            ids=["s1"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            map_table_strategy="config",
        )
        comp = Mock(
            ids=["c1"],
            streams={"compounds": {"format": "sdf", "file": "<id>.sdf"}},
            map_table_strategy="config",
        )
        Gnina(
            structures=struct.streams.structures,
            compounds=comp.streams.compounds,
            **gnina_kwargs,
        )
        script_path = pipeline.save()
    return read_all_emitted_artifacts(script_path)


def test_exhaustiveness(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, exhaustiveness=16)
    assert '"exhaustiveness": 16' in content


def test_num_modes(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, num_modes=20)
    assert '"num_modes": 20' in content


def test_num_runs(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, num_runs=3)
    assert '"num_runs": 3' in content


def test_seed(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, seed=7)
    assert '"seed": 7' in content


def test_cnn_scoring(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, cnn_scoring="refinement")
    assert '"cnn_scoring": "refinement"' in content


def test_generate_conformers(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, generate_conformers=True, num_conformers=12)
    assert '"generate_conformers": true' in content
    assert '"num_conformers": 12' in content


def test_pH(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, pH=6.0, protonate=True)
    assert '"pH": 6.0' in content
    assert '"protonate": true' in content


def test_thresholds(local_config, isolated_cwd, new_pipeline):
    content = _build(
        local_config, isolated_cwd, new_pipeline,
        cnn_score_threshold=0.7,
        rmsd_threshold=1.5,
        energy_window=3.0,
        conformer_rmsd=0.8,
    )
    assert_substrings_in(content, [
        '"cnn_score_threshold": 0.7',
        '"rmsd_threshold": 1.5',
        '"energy_window": 3.0',
        '"conformer_rmsd": 0.8',
    ])


def test_box_center_size(local_config, isolated_cwd, new_pipeline):
    content = _build(
        local_config, isolated_cwd, new_pipeline,
        center="(10, 20, 30)",
        size=18.0,
    )
    assert "10" in content and "20" in content and "30" in content
    assert "18" in content


def test_mode_default_is_docking(local_config, isolated_cwd, new_pipeline):
    content, out = _build_with_outputs(local_config, isolated_cwd, new_pipeline)
    assert '"mode": "docking"' in content
    assert sorted(out.tables._tables.keys()) == ["docking_results", "docking_summary"]


@pytest.mark.parametrize("mode", ["score", "minimize"])
def test_score_modes_emit_scores_table(local_config, isolated_cwd, new_pipeline, mode):
    content, out = _build_with_outputs(local_config, isolated_cwd, new_pipeline, mode=mode)
    assert f'"mode": "{mode}"' in content
    assert "scores" in out.tables._tables
    assert "missing" in out.tables._tables
    assert "docking_results" not in out.tables._tables
    assert out.tables.scores.info.columns == [
        "id", "structures.id", "compounds.id", "vina_affinity",
        "cnn_score", "cnn_affinity", "cnn_vs", "cnn_affinity_variance",
    ]
    assert out.streams.structures.files[0].endswith(f"<id>_{mode}.pdb")


def test_score_mode_rejects_docking_only_params(local_config, isolated_cwd, new_pipeline):
    with pytest.raises(ValueError, match="docking-only parameters"):
        _build_with_outputs(local_config, isolated_cwd, new_pipeline,
                            mode="score", exhaustiveness=16)


def test_invalid_mode_raises(local_config, isolated_cwd, new_pipeline):
    with pytest.raises(ValueError, match="mode must be one of"):
        _build_with_outputs(local_config, isolated_cwd, new_pipeline, mode="dock")


def test_smoke_all_params(local_config, isolated_cwd, new_pipeline):
    content = _build(
        local_config, isolated_cwd, new_pipeline,
        center="(0, 0, 0)",
        size=20.0,
        autobox_add=5.0,
        exhaustiveness=12,
        num_modes=10,
        num_runs=2,
        seed=55,
        cnn_scoring="rescore",
        generate_conformers=True,
        num_conformers=10,
        energy_window=2.5,
        conformer_rmsd=1.0,
        cnn_score_threshold=0.4,
        rmsd_threshold=2.5,
        protonate=True,
        pH=7.0,
    )
    assert_substrings_in(content, [
        '"exhaustiveness": 12',
        '"num_modes": 10',
        '"num_runs": 2',
        '"seed": 55',
        '"cnn_scoring": "rescore"',
        '"generate_conformers": true',
        '"num_conformers": 10',
        '"energy_window": 2.5',
        '"cnn_score_threshold": 0.4',
        '"rmsd_threshold": 2.5',
        '"protonate": true',
        '"pH": 7.0',
    ])
