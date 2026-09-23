"""BoltzGen's local 0.3.1 patches and per-code SMILES overrides.

The patches lived only in site-packages, so `bp-warm BoltzGen` or an env rebuild dropped them and BoltzGen silently behaved differently; `smiles_overrides` went into the job script as `export X="<smiles>"`, so a quote or `$` in a SMILES broke or injected into it.
"""

import os

import pytest

REPO = os.path.join(os.path.dirname(__file__), "..")


@pytest.fixture(autouse=True)
def _env(monkeypatch):
    from biopipelines.boltzgen import BoltzGen
    monkeypatch.setattr(BoltzGen, "_install_env", classmethod(lambda cls, mgr: "boltzgen"))


def _script(force_reinstall=False):
    from biopipelines.boltzgen import BoltzGen
    return BoltzGen._install_script({"biopipelines": "/repo"}, env_manager="mamba",
                                    force_reinstall=force_reinstall)


@pytest.mark.parametrize("force", [False, True])
def test_every_install_path_applies_or_confirms_the_patches(local_config, force):
    script = _script(force_reinstall=force)
    assert "/repo/environments/patches/boltzgen-0.3.1-locbp.patch" in script
    assert "-R --dry-run" in script and "--forward" in script


def test_the_already_installed_shortcut_still_checks_the_patches(local_config):
    script = _script()
    shortcut = script.split("already installed, skipping")[1].split("exit 0")[0]
    assert "BG_PATCH" in shortcut, "a rebuilt env reached through the shortcut kept no patches"


def test_boltzgen_is_pinned_to_the_version_the_patch_targets():
    pins = open(os.path.join(REPO, "environments", "boltzgen.pip.txt"), encoding="utf-8").read().split()
    assert "boltzgen==0.3.1" in pins


@pytest.mark.parametrize("overrides, message", [
    ({"lig": "CCO"}, "CCD-style code"),
    ({"LIG$": "CCO"}, "CCD-style code"),
    ({"LIG": "CC O"}, "single SMILES"),
])
def test_a_bad_override_is_refused_at_configuration(local_config, overrides, message):
    from biopipelines.boltzgen import BoltzGen
    tool = BoltzGen.__new__(BoltzGen)
    tool.smiles_overrides = overrides
    with pytest.raises(ValueError, match=message):
        tool.validate_params()


def test_the_export_quotes_the_smiles():
    source = open(os.path.join(REPO, "biopipelines", "boltzgen.py"), encoding="utf-8").read()
    assert "export BOLTZGEN_SMILES_{_code}={shlex.quote(_smi)}" in source


def test_a_missing_patch_binary_is_named_not_reported_as_a_version_mismatch(local_config):
    script = _script()
    assert "command -v patch" in script and "is not on PATH" in script
