"""BaseConfig.install() checks that a tool defines an install script by calling it with an empty folders dict.

A tool whose install script names a folder of its own (IPSAE, RFD4Proteina, ESM2) raises ValueError for that dict, so `Tool.install()` used to fail at pipeline build time with the message meant for a missing config key, even when the key was set. The probe now treats the ValueError as "defined"; the real folders reach `_install_script` when the step's script is generated.
"""

import pytest


def add_a_step(new_pipeline, job):
    """A pipeline that is already entered and holds one step, so leaving it can save."""
    from biopipelines.mock import Mock

    pipeline = new_pipeline(job)
    pipeline.__enter__()
    Mock(ids=["s1"], streams={"sequences": {"format": "csv", "values": ["MKT"]}})
    return pipeline


def test_a_tool_that_needs_its_own_folder_can_be_installed_in_a_pipeline(local_config, isolated_cwd, new_pipeline):
    from biopipelines.ipsae import IPSAE

    pipeline = add_a_step(new_pipeline, "install_ipsae")
    try:
        IPSAE.install()
    finally:
        pipeline.__exit__(None, None, None)


def test_a_folder_the_config_lacks_is_reported_when_the_install_step_is_generated(local_config, isolated_cwd, new_pipeline):
    from biopipelines.rfd4proteina import RFD4Proteina

    pipeline = add_a_step(new_pipeline, "install_rfd4")
    RFD4Proteina.install()
    with pytest.raises(ValueError, match="repositories.RFD4Proteina"):
        pipeline.__exit__(None, None, None)


def test_the_missing_folder_is_still_reported_by_the_install_script_itself():
    from biopipelines.ipsae import IPSAE

    with pytest.raises(ValueError, match="repositories.IPSAE"):
        IPSAE._install_script({}, "mamba")


def test_a_tool_without_an_install_script_still_reports_it(local_config, isolated_cwd, new_pipeline):
    from biopipelines.structure_cluster import StructureCluster

    pipeline = add_a_step(new_pipeline, "no_install")
    try:
        with pytest.raises(NotImplementedError, match="does not define installation steps"):
            StructureCluster.install()
    finally:
        pipeline.__exit__(None, None, None)
