# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Parameter coverage for biopipelines.bioemu2.BioEmu2.

A single-step Pipeline is built and saved, and the emitted artifacts are read. BioEmu2 sends its parameters in `_configuration/config.yaml` and not on the command line, so that file is parsed and the generated script is only checked for what it must not contain.
The two runtime switches, BIOEMU2_SRUN_ENVIRONMENT and BIOEMU2_MSA_CACHE, are read where the helper runs. At configuration time the check is that neither is baked into the script: a hardcoded `srun --environment=` or cache path would make the variable a no-op.
"""

from pathlib import Path

import pytest
import yaml

from ._helpers import assert_substrings_in, read_all_emitted_artifacts


pytestmark = pytest.mark.tool_parameters

CONTACT = {"token1": {"asym_id": "A", "res_id": 10}, "token2": {"asym_id": "A", "res_id": 20}}
LINK = {"chain": "A", "position": 30, "protein_atom": "SG", "ligand_atom": "C1"}


def _build(new_pipeline, with_ligand=False, with_msas=False, **options):
    """(parsed config.yaml of the step, every emitted artifact as one string)."""
    from biopipelines.ligand import Ligand
    from biopipelines.bioemu2 import BioEmu2
    from biopipelines.mock import Mock
    from biopipelines.sequence import Sequence

    pipeline = new_pipeline("bioemu2_params")
    with pipeline:
        sequence = Sequence(seq="MKTAYIAKQRQISFVKSHFSRQLEERLGL", type="protein", ids="p1")
        if with_ligand:
            options["ligands"] = Ligand(smiles="CCO", ids="eth")
        if with_msas:
            options["msas"] = Mock(ids=["p1"], streams={"msas": {"format": "a3m", "file": "<id>.a3m"}}
                                   ).streams.msas
        BioEmu2(sequences=sequence, **options)
        script_path = pipeline.save()
    run_root = Path(script_path).resolve().parent.parent
    (config_path,) = [p for p in run_root.rglob("_configuration/config.yaml") if "BioEmu2" in str(p)]
    with open(config_path) as handle:
        return yaml.safe_load(handle), read_all_emitted_artifacts(script_path)


def test_num_samples(local_config, isolated_cwd, new_pipeline):
    config, _content = _build(new_pipeline, num_samples=25)
    assert config["num_samples"] == 25


def test_batch_size(local_config, isolated_cwd, new_pipeline):
    config, _content = _build(new_pipeline, batch_size=8)
    assert config["batch_size"] == 8


def test_seed(local_config, isolated_cwd, new_pipeline):
    config, _content = _build(new_pipeline, seed=4242)
    assert config["seed"] == 4242


def test_seed_absent_by_default(local_config, isolated_cwd, new_pipeline):
    config, _content = _build(new_pipeline)
    assert config["seed"] is None


def test_sampler(local_config, isolated_cwd, new_pipeline):
    config, _content = _build(new_pipeline, sampler="edm_gd")
    assert config["sampler"] == "edm_gd"


def test_physicality_filter(local_config, isolated_cwd, new_pipeline):
    config, _content = _build(new_pipeline, physicality_filter="no_filter")
    assert config["physicality_filter"] == "no_filter"


def test_contacts_are_written_unchanged(local_config, isolated_cwd, new_pipeline):
    config, _content = _build(new_pipeline, contacts=[CONTACT])
    assert config["contacts"] == [CONTACT]


def test_covalent_linkage_is_written_unchanged(local_config, isolated_cwd, new_pipeline):
    config, _content = _build(new_pipeline, with_ligand=True, covalent_linkage=LINK)
    assert config["covalent_linkage"] == LINK


def test_msas_flags_follow_the_parameter(local_config, isolated_cwd, new_pipeline):
    _config, plain = _build(new_pipeline)
    assert "--msas-json" not in plain and "--msa-cache-dir" not in plain
    _config, supplied = _build(new_pipeline, with_msas=True)
    assert_substrings_in(supplied, ["--msas-json", "--msa-cache-dir"])


def test_no_srun_wrapper_at_config_time(local_config, isolated_cwd, new_pipeline):
    """The EDF is the helper's runtime environment, never part of the generated script."""
    _config, content = _build(new_pipeline)
    assert "srun --environment=" not in content, (
        "the generated script hardcodes a Daint-specific srun wrapper; the image has to come from "
        "BIOEMU2_SRUN_ENVIRONMENT at runtime"
    )
    assert "BIOEMU2_SRUN_ENVIRONMENT" not in content


def test_no_hardcoded_msa_cache_path(local_config, isolated_cwd, new_pipeline):
    """bioemu2 reads BIOEMU2_MSA_CACHE itself; a path written into the script would make the override a no-op."""
    _config, content = _build(new_pipeline)
    assert ".cache/bioemu2/colabfold_msa" not in content


def test_smoke_all_params(local_config, isolated_cwd, new_pipeline):
    config, _content = _build(
        new_pipeline, with_ligand=True, num_samples=10, batch_size=4, seed=7, sampler="edm",
        physicality_filter="high_quality_per_residue", contacts=[CONTACT], covalent_linkage=LINK,
    )
    assert config == {"num_samples": 10, "batch_size": 4, "seed": 7, "sampler": "edm",
                      "physicality_filter": "high_quality_per_residue", "contacts": [CONTACT],
                      "covalent_linkage": LINK}
