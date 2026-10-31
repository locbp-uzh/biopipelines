"""config.daint.yaml routes only what can run on aarch64, and its docs say the same.

ESMFold2 stayed in the Daint config, in the README badges and in the list of tools verified by running after its environment became uninstallable: the wrapper needs esm >= 3.4 (torch 2.11), the EDF image ships torch 2.8, and the esm 3.3.0 commit that accepts that torch needs a transformers fork that no longer exists.
"""

import re
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parent.parent


def _daint():
    return yaml.safe_load((ROOT / "config.daint.yaml").read_text(encoding="utf-8"))


def test_esmfold2_has_no_route_on_daint():
    config = _daint()
    assert "ESMFold2" not in config["environments"]
    assert "ESMFold2" not in config["edf"]


def test_esmfold2_install_is_refused_by_name_under_daint(shipped_config):
    from biopipelines.esmfold2 import ESMFold2

    manager = shipped_config("daint")
    with pytest.raises(ValueError, match=r"ESMFold2: no environment is configured for it in config\.daint\.yaml"):
        ESMFold2._install_script({"biopipelines": str(ROOT)}, manager.get_env_manager())


def _readme_platforms(tool):
    row = re.search(rf'<td><sub><b>{re.escape(tool)}</b>.*?</tr>', (ROOT / "README.md").read_text(encoding="utf-8"), re.S).group(0)
    return re.findall(r'alt="(HPC aarch64|HPC x86-64|Colab) ok"', row)


@pytest.mark.parametrize("tool", ["ESMFold2", "BioEmu2"])
def test_no_aarch64_badge_for_a_tool_with_no_verified_daint_run(tool):
    assert "HPC aarch64" not in _readme_platforms(tool)


def test_the_list_of_verified_tools_does_not_carry_esmfold2_or_bioemu2():
    text = (ROOT / "skills" / "biopipelines" / "references" / "daint_backend.md").read_text(encoding="utf-8")
    verified = re.search(r"These \d+ are verified by running.*?\n", text).group(0)
    assert "ESMFold2" not in verified and "BioEmu2" not in verified
