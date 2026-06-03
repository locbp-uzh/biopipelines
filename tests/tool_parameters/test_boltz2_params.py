"""Parameter coverage for biopipelines.boltz2.Boltz2.

Boltz2 builds its YAML config at runtime via pipe_boltz_config_unified.py;
at config time the only artifact we can grep is pipeline.sh, where the
wrapper emits ``--<flag> '<json>'`` substrings on the runner invocation.
That's enough to verify each kwarg is wired through to the runner.
"""

import pytest

from ._helpers import (
    assert_kwarg_emitted,
    assert_substrings_in,
    read_pipeline_sh,
)


pytestmark = pytest.mark.tool_parameters


def _build(local_config, isolated_cwd, new_pipeline, **boltz_kwargs):
    from biopipelines.sequence import Sequence
    from biopipelines.boltz2 import Boltz2

    pipeline = new_pipeline("boltz_params")
    with pipeline:
        s = Sequence(seq="MKTAYIAKQRQISFVKSHFSRQLEERLGL", type="protein", ids="p1")
        Boltz2(proteins=s, **boltz_kwargs)
        script_path = pipeline.save()

    return read_pipeline_sh(script_path)


def test_recycling_steps(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, recycling_steps=5)
    assert "5" in content
    assert "recycling" in content.lower() or "--recycling" in content


def test_diffusion_samples(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, diffusion_samples=4)
    assert "4" in content
    assert "diffusion" in content.lower() or "--diffusion" in content


def test_use_potentials(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, use_potentials=True)
    assert "potential" in content.lower() or "use_potentials" in content


def test_output_format_mmcif(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline, output_format="mmcif")
    assert "mmcif" in content


def test_msa_server_used_when_no_msas(local_config, isolated_cwd, new_pipeline):
    content = _build(local_config, isolated_cwd, new_pipeline)
    assert "--use_msa_server" in content


def test_template(local_config, isolated_cwd, new_pipeline):
    content = _build(
        local_config, isolated_cwd, new_pipeline,
        template="/tmp/template.pdb",
        template_chain_ids=["A"],
        template_threshold=4.0,
    )
    assert "/tmp/template.pdb" in content
    assert "--template" in content


def test_pocket_constraints(local_config, isolated_cwd, new_pipeline):
    content = _build(
        local_config, isolated_cwd, new_pipeline,
        pocket_residues=[10, 20, 30],
        pocket_max_distance=8.0,
        pocket_force=True,
    )
    assert_substrings_in(content, [
        "--pocket-residues",
        "--pocket-max-distance 8.0",
        "--pocket-force",
    ])


def test_contacts(local_config, isolated_cwd, new_pipeline):
    content = _build(
        local_config, isolated_cwd, new_pipeline,
        contacts=[{"token1": ["A", 10], "token2": ["A", 20], "max_distance": 6.0}],
    )
    assert "--contacts" in content
    assert "token1" in content


def test_glycosylation(local_config, isolated_cwd, new_pipeline):
    content = _build(
        local_config, isolated_cwd, new_pipeline,
        glycosylation={"A": [25]},
    )
    assert "--glycosylation" in content


def test_covalent_linkage(local_config, isolated_cwd, new_pipeline):
    content = _build(
        local_config, isolated_cwd, new_pipeline,
        covalent_linkage={
            "chain": "A", "position": 30, "ligand_atom": "C1",
        },
    )
    assert "--covalent-linkage" in content


# ── Track-1 additions ────────────────────────────────────────────────────────

def test_disulfide_bonds(local_config, isolated_cwd, new_pipeline):
    content = _build(
        local_config, isolated_cwd, new_pipeline,
        disulfide_bonds=[{"token1": ["A", 12], "token2": ["A", 45]}],
    )
    assert "--disulfide-bonds" in content
    assert "12" in content and "45" in content


def test_metal_coord(local_config, isolated_cwd, new_pipeline):
    content = _build(
        local_config, isolated_cwd, new_pipeline,
        metal_coord=[{
            "atom1": ["A", 30, "ND1"],
            "atom2": ["B", 1, "ZN"],
        }],
    )
    assert "--metal-coord" in content
    assert "ZN" in content
    assert "ND1" in content


def test_smoke_all_params(local_config, isolated_cwd, new_pipeline):
    content = _build(
        local_config, isolated_cwd, new_pipeline,
        recycling_steps=3,
        diffusion_samples=2,
        use_potentials=True,
        output_format="pdb",
        template="/tmp/t.pdb",
        template_chain_ids=["A"],
        template_force=True,
        template_threshold=5.0,
        pocket_residues=[1, 2, 3],
        pocket_max_distance=7.0,
        pocket_force=True,
        contacts=[{"token1": ["A", 1], "token2": ["A", 2]}],
        glycosylation={"A": [10]},
        covalent_linkage={"chain": "A", "position": 5, "ligand_atom": "C1"},
        disulfide_bonds=[{"token1": ["A", 7], "token2": ["A", 9]}],
        metal_coord=[{"atom1": ["A", 11, "ND1"], "atom2": ["B", 1, "ZN"]}],
    )
    assert_substrings_in(content, [
        "--template",
        "--template-force",
        "--pocket-residues",
        "--pocket-force",
        "--contacts",
        "--glycosylation",
        "--covalent-linkage",
        "--disulfide-bonds",
        "--metal-coord",
    ])
