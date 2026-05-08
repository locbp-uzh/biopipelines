"""Parameter coverage for biopipelines.protein_mpnn.ProteinMPNN.

For each user-facing kwarg, drive a one-step Mock → ProteinMPNN pipeline and
assert that the kwarg's sentinel value appears in the emitted ``pipeline.sh``
next to the upstream CLI flag the wrapper is expected to emit. The test
runs at config time only — no GPU, no upstream binary execution.
"""

import pytest

from ._helpers import (
    assert_kwarg_emitted,
    assert_substrings_in,
    read_pipeline_sh,
)


pytestmark = pytest.mark.tool_parameters


def _build(local_config, isolated_cwd, new_pipeline, **pmpnn_kwargs):
    """Run a single Mock → ProteinMPNN pipeline and return pipeline.sh content."""
    from biopipelines.mock import Mock
    from biopipelines.protein_mpnn import ProteinMPNN

    pipeline = new_pipeline("pmpnn_params")
    with pipeline:
        m = Mock(
            ids=["d1", "d2"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            map_table_strategy="config",
        )
        ProteinMPNN(structures=m.streams.structures, **pmpnn_kwargs)
        script_path = pipeline.save()

    return read_pipeline_sh(script_path)


def test_num_sequences(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, num_sequences=7)
    assert_kwarg_emitted(content, "num_sequences", 7, flag="--num_seq_per_target 7")


def test_sampling_temp(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, sampling_temp=0.42)
    assert_kwarg_emitted(content, "sampling_temp", 0.42, flag="--sampling_temp 0.42")


def test_model_name(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, model_name="v_48_030")
    assert_kwarg_emitted(content, "model_name", "v_48_030", flag="--model_name v_48_030")


def test_soluble_model_off(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, soluble_model=False)
    assert "--use_soluble_model" not in content


def test_soluble_model_on(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, soluble_model=True)
    assert "--use_soluble_model" in content


def test_chain(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, chain="B")
    assert '"B"' in content or " B " in content, "chain value not in pipeline.sh"


# ── Track-1 additions ────────────────────────────────────────────────────────

def test_bias_AA_jsonl(local_config, isolated_cwd, new_pipeline):
    content = _build(
        local_config, isolated_cwd, new_pipeline, bias_AA_jsonl="/tmp/bias.jsonl"
    )
    assert_kwarg_emitted(
        content, "bias_AA_jsonl", "/tmp/bias.jsonl",
        flag="--bias_AA_jsonl /tmp/bias.jsonl",
    )


def test_omit_AA_jsonl(local_config, isolated_cwd, new_pipeline):
    content = _build(
        local_config, isolated_cwd, new_pipeline, omit_AA_jsonl="/tmp/omit.jsonl"
    )
    assert_kwarg_emitted(
        content, "omit_AA_jsonl", "/tmp/omit.jsonl",
        flag="--omit_AA_jsonl /tmp/omit.jsonl",
    )


def test_seed(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, seed=12345)
    assert_kwarg_emitted(content, "seed", 12345, flag="--seed 12345")


def test_ca_noise_std(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, ca_noise_std=0.25)
    assert_kwarg_emitted(content, "ca_noise_std", 0.25, flag="--backbone_noise 0.25")


def test_smoke_all_params(local_config, isolated_cwd, new_pipeline):
    """Exercise every kwarg simultaneously and assert no flag drops out."""
    content = _build(
        local_config, isolated_cwd, new_pipeline,
        num_sequences=3,
        sampling_temp=0.15,
        model_name="v_48_010",
        soluble_model=True,
        bias_AA_jsonl="/tmp/bias.jsonl",
        omit_AA_jsonl="/tmp/omit.jsonl",
        seed=99,
        ca_noise_std=0.05,
    )
    assert_substrings_in(content, [
        "--num_seq_per_target 3",
        "--sampling_temp 0.15",
        "--model_name v_48_010",
        "--use_soluble_model",
        "--bias_AA_jsonl /tmp/bias.jsonl",
        "--omit_AA_jsonl /tmp/omit.jsonl",
        "--seed 99",
        "--backbone_noise 0.05",
    ])
