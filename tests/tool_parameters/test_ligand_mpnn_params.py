"""Parameter coverage for biopipelines.ligand_mpnn.LigandMPNN."""

import pytest

from ._helpers import (
    assert_kwarg_emitted,
    assert_substrings_in,
    read_all_emitted_artifacts,
)


pytestmark = pytest.mark.tool_parameters


def _build(local_config, isolated_cwd, new_pipeline, **lmpnn_kwargs):
    from biopipelines.mock import Mock
    from biopipelines.ligand_mpnn import LigandMPNN
    from biopipelines.ligand import Ligand

    ligand_code = lmpnn_kwargs.pop("ligand", "ATP")

    pipeline = new_pipeline("lmpnn_params")
    with pipeline:
        m = Mock(
            ids=["s1", "s2"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            map_table_strategy="config",
        )
        LigandMPNN(structures=m.streams.structures, ligand=Ligand(code=ligand_code), **lmpnn_kwargs)
        script_path = pipeline.save()

    return read_all_emitted_artifacts(script_path)


def test_ligand(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, ligand="HEM")
    # The compounds-stream JSON is passed to the positions script, which reads
    # the residue `code` from it at runtime.
    assert "input_ligand.json" in content


def test_num_sequences_maps_to_batch_size(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, num_sequences=4)
    assert_kwarg_emitted(content, "num_sequences", 4, flag="--batch_size 4")


def test_num_batches(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, num_batches=3)
    assert_kwarg_emitted(content, "num_batches", 3, flag="--number_of_batches 3")


def test_design_within(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, design_within=8.0)
    assert_substrings_in(content, ["--ligand_mpnn_cutoff_for_score 8.0"])


def test_redesigned_table_reference_serializes(
    local_config, isolated_cwd, new_pipeline,
):
    from biopipelines.base_config import TableInfo

    reference = TableInfo(
        name="selections",
        path="/tmp/selections.csv",
        columns=["within"],
    ).within
    content = _build(
        local_config,
        isolated_cwd,
        new_pipeline,
        redesigned=reference,
    )
    assert_substrings_in(
        content,
        ["TABLE_REFERENCE:/tmp/selections.csv:within"],
        label="redesigned TableReference",
    )


def test_model(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, model="v_32_020")
    assert "v_32_020" in content


# ── Track-1 additions ────────────────────────────────────────────────────────

def test_temperature(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, temperature=0.3)
    assert_kwarg_emitted(content, "temperature", 0.3, flag="--temperature 0.3")


def test_bias_AA_per_residue(local_config, isolated_cwd, new_pipeline):
    content = _build(
        local_config, isolated_cwd, new_pipeline, bias_AA_per_residue="/tmp/bias.json"
    )
    assert "--bias_AA_per_residue" in content
    assert "/tmp/bias.json" in content


def test_seed(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, seed=4242)
    assert_kwarg_emitted(content, "seed", 4242, flag="--seed 4242")


def test_smoke_all_params(local_config, isolated_cwd, new_pipeline):
    content = _build(
        local_config, isolated_cwd, new_pipeline,
        ligand="ATP",
        num_sequences=2,
        num_batches=3,
        design_within=6.5,
        model="v_32_010",
        temperature=0.25,
        bias_AA_per_residue="/tmp/bias.json",
        seed=7,
    )
    assert_substrings_in(content, [
        "ATP",
        "--batch_size 2",
        "--number_of_batches 3",
        "--ligand_mpnn_cutoff_for_score",
        "v_32_010",
        "--temperature 0.25",
        "--bias_AA_per_residue",
        "--seed 7",
    ])


def test_pack_side_chains_off_by_default(local_config, isolated_cwd, new_pipeline):
    """No packing flags, and no structures stream, unless asked for."""
    content = _build(local_config, isolated_cwd, new_pipeline)
    assert "--pack_side_chains" not in content


def test_pack_side_chains_emits_flags(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, pack_side_chains=True)
    assert_substrings_in(content, [
        "--pack_side_chains 1",
        "--number_of_packs_per_design 1",
        "--pack_with_ligand_context 1",
    ])


def test_packs_per_design(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline,
                     pack_side_chains=True, packs_per_design=4)
    assert_kwarg_emitted(content, "packs_per_design", 4,
                         flag="--number_of_packs_per_design 4")


def test_pack_without_ligand_context(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline,
                     pack_side_chains=True, pack_with_ligand_context=False)
    assert_substrings_in(content, ["--pack_with_ligand_context 0"])


def test_packed_files_are_renamed_into_the_stream_folder(
    local_config, isolated_cwd, new_pipeline,
):
    """Upstream writes <struct>_packed_<seq>_<pack>.pdb; the stream ids carry no
    `_packed_` infix, so the script must rename as it moves them."""
    content = _build(local_config, isolated_cwd, new_pipeline, pack_side_chains=True)
    assert "/packed/*.pdb" in content
    assert "_packed_/_" in content
