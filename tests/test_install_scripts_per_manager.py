"""Every tool's ``_install_script`` must build under every SHIPPED config.

``install()`` runs at configuration time on the user's own machine, so a
manager-specific helper called without a guard raises there and never reaches
CI. ``OpenMM`` did exactly that: its antechamber check built a PATH from
``get_conda_env_root()``, which ``config_manager`` defines only for the
``venv`` manager, so a forced install raised ``KeyError: venv_root is
required`` on four of the five shipped variants -- and the message named a key
the variant does not use.

The invariant under test is not "every tool installs everywhere". A variant
that configures no environment for a tool is a real, documented condition and
the framework reports it as a ``ValueError`` naming the config file. What must
never happen is an internal error type -- ``KeyError``, ``AttributeError``,
``TypeError`` -- because that is a bug in the tool, not a statement about the
config.
"""
from __future__ import annotations

from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]

# Every config.<variant>.yaml the repo ships. daint is aarch64-only and
# deliberately configures a subset; it is included because the guard bug was
# invisible precisely there -- daint was the one variant that worked.
SHIPPED_VARIANTS = ["cluster", "colab", "container", "daint", "local"]

INTERNAL_ERRORS = (KeyError, AttributeError, TypeError, IndexError)


@pytest.fixture
def shipped_config(monkeypatch):
    """Load a repo-root ``config.<variant>.yaml`` and yield its ConfigManager."""
    from biopipelines.config_manager import ConfigManager

    def _reset():
        ConfigManager._instance = None
        ConfigManager._config = None
        ConfigManager._variant = None

    def _load(variant: str):
        path = REPO_ROOT / f"config.{variant}.yaml"
        assert path.exists(), f"Missing shipped config: {path}"
        _reset()
        monkeypatch.setattr(
            ConfigManager, "_get_config_path",
            classmethod(lambda cls, variant=None, _p=str(path): _p),
        )
        return ConfigManager(variant=variant)

    _reset()
    yield _load
    _reset()


def _installable_tools():
    """Every exported tool class that defines its own ``_install_script``."""
    import biopipelines
    from biopipelines.base_config import BaseConfig

    out = []
    for name in getattr(biopipelines, "__all__", []):
        cls = getattr(biopipelines, name, None)
        if not isinstance(cls, type) or not issubclass(cls, BaseConfig):
            continue
        if cls.__dict__.get("_install_script") is None:
            continue
        out.append((name, cls))
    return sorted(out, key=lambda p: p[0])


TOOLS = _installable_tools()


def test_at_least_one_tool_defines_an_install_script():
    """Guard the discovery -- an empty list would make every case below vacuous."""
    assert TOOLS, "no tool class defines its own _install_script"


@pytest.mark.parametrize("variant", SHIPPED_VARIANTS)
@pytest.mark.parametrize("name,cls", TOOLS, ids=[n for n, _ in TOOLS])
def test_install_script_never_raises_an_internal_error(
    name, cls, variant, shipped_config, tmp_path,
):
    """A config problem must be a ValueError naming the config, not a KeyError."""
    cm = shipped_config(variant)
    folders = {"install": str(tmp_path), "runtime": str(tmp_path), "logs": str(tmp_path)}

    try:
        script = cls._install_script(folders, cm.get_env_manager())
    except ValueError:
        return  # the framework's own diagnostic; the tool behaved correctly
    except INTERNAL_ERRORS as exc:
        pytest.fail(
            f"{name}._install_script raised {type(exc).__name__} under "
            f"config.{variant}.yaml (env_manager={cm.get_env_manager()!r}): {exc}"
        )

    assert script is None or isinstance(script, str)


@pytest.mark.parametrize("variant", SHIPPED_VARIANTS)
def test_openmm_install_script_builds_under_every_shipped_config(
    variant, shipped_config, tmp_path,
):
    """OpenMM specifically: the antechamber check is the guarded part.

    Named separately so the regression is legible in the failure list rather
    than buried in the per-tool sweep above.
    """
    from biopipelines.openmm import OpenMM

    cm = shipped_config(variant)
    try:
        script = OpenMM._install_script({"install": str(tmp_path)}, cm.get_env_manager())
    except ValueError:
        pytest.skip(f"config.{variant}.yaml configures no environment for OpenMM")

    assert isinstance(script, str) and script
    assert "antechamber" in script
    # Only the venv manager has a conda env root to build an explicit PATH from;
    # every other manager goes through the env run prefix instead.
    if cm.get_env_manager() == "venv":
        assert cm.get_conda_env_root() in script
    else:
        assert "/bin:$PATH" not in script
