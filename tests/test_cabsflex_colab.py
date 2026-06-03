"""CABSflex Colab-compatibility wiring tests.

CABSflex runs in an isolated micromamba py2.7 env (the framework's per-tool
sandbox). On Colab two things differ from the cluster, both verified live on a
Colab runtime before these tests were written:

  1. The py2.7 env must omit `modeller` — the only modeller builds on linux-64
     require Python >=3.6, unsatisfiable against cabs's python=2.7 pin, so the
     full solve fails. `environments/CABSflex.colab.yaml` drops it; the env then
     solves and `import CABS` passes.
  2. CABS's py2.7 matplotlib must use a headless backend. Colab exports
     MPLBACKEND=module://matplotlib_inline.backend_inline, which the py2.7 env
     cannot import and which crashes CABS at startup. The run script forces
     MPLBACKEND=Agg.

Because modeller is absent on Colab, aa_rebuild=True (which needs MODELLER) is
rejected at config time with a clear error rather than crashing mid-run.

These are config-time wiring checks (no CABS execution).
"""
from __future__ import annotations

import glob
import os
from unittest.mock import patch

import pytest


def _tool_script(pipeline_sh: str, tool_name: str) -> str:
    """Read the per-tool ``NNN_<Tool>.sh`` emitted next to pipeline.sh."""
    runtime_dir = os.path.dirname(pipeline_sh)
    matches = sorted(glob.glob(os.path.join(runtime_dir, f"*_{tool_name}.sh")))
    assert matches, f"no {tool_name} script in {runtime_dir}"
    with open(matches[0], "rb") as f:
        return f.read().decode("utf-8", errors="replace")


def test_colab_install_uses_modeller_free_env():
    """Under the colab scheduler, the install script targets
    CABSflex.colab.yaml (no modeller) via micromamba and drops the MODELLER
    license note."""
    from biopipelines.config_manager import ConfigManager
    from biopipelines.cabsflex import CABSflex

    with patch.object(ConfigManager, "get_scheduler", lambda self: "colab"), \
         patch.object(ConfigManager, "get_variant", lambda self: "colab"):
        script = CABSflex._install_script({"biopipelines": "/repo"},
                                          env_manager="micromamba")

    assert "CABSflex.colab.yaml" in script
    assert "micromamba env create" in script
    # Env-manager-agnostic skip check, not `conda list` (absent on Colab).
    assert "conda list" not in script
    assert 'micromamba run -n CABSflex python -c "import CABS"' in script
    # No MODELLER license echo on Colab (modeller isn't installed there).
    assert "KEY_MODELLER" not in script


def test_cluster_install_keeps_modeller_note():
    """Under a non-colab scheduler the MODELLER license note is retained and
    the variant-less / cluster yaml is used."""
    from biopipelines.config_manager import ConfigManager
    from biopipelines.cabsflex import CABSflex

    with patch.object(ConfigManager, "get_scheduler", lambda self: "slurm"), \
         patch.object(ConfigManager, "get_variant", lambda self: "cluster"):
        script = CABSflex._install_script({"biopipelines": "/repo"},
                                          env_manager="mamba")

    assert "KEY_MODELLER" in script
    assert "CABSflex.colab.yaml" not in script


def test_colab_env_yaml_omits_modeller():
    """The committed Colab env spec must not pin modeller as a dependency."""
    import yaml

    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    with open(os.path.join(repo_root, "environments", "CABSflex.colab.yaml")) as f:
        spec = yaml.safe_load(f)
    deps = [str(d) for d in spec["dependencies"]]
    assert not any("modeller" in d for d in deps), deps
    assert any(d.startswith("python=2.7") for d in deps), deps
    assert any(d == "cabs" or d.startswith("cabs") for d in deps), deps


def test_run_script_forces_headless_matplotlib_backend(
    local_config, isolated_cwd, new_pipeline,
):
    """The generated CABSflex run script exports MPLBACKEND=Agg before invoking
    CABS, so Colab's inline backend leak can't crash the py2.7 matplotlib."""
    from biopipelines.mock import Mock
    from biopipelines.cabsflex import CABSflex

    pipeline = new_pipeline("cabsflex_mpl")
    with pipeline:
        m = Mock(ids=["s1"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
        CABSflex(structures=m, num_models=2)
        script_path = pipeline.save()

    script = _tool_script(script_path, "CABSflex")
    assert "export MPLBACKEND=Agg" in script
    # The export precedes the CABSflex invocation.
    assert script.index("export MPLBACKEND=Agg") < script.index("CABSflex -i")


def test_aa_rebuild_rejected_on_colab(
    isolated_cwd, new_pipeline,
):
    """aa_rebuild=True must raise a clear error under the colab scheduler."""
    from biopipelines.config_manager import ConfigManager
    from biopipelines.mock import Mock
    from biopipelines.cabsflex import CABSflex

    with patch.object(ConfigManager, "get_scheduler", lambda self: "colab"):
        pipeline = new_pipeline("cabsflex_aa")
        with pytest.raises(ValueError, match="MODELLER"):
            with pipeline:
                m = Mock(ids=["s1"],
                         streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
                CABSflex(structures=m, num_models=2, aa_rebuild=True)


def test_aa_rebuild_allowed_off_colab(
    local_config, isolated_cwd, new_pipeline,
):
    """aa_rebuild=True is fine under a non-colab scheduler (cluster has MODELLER)."""
    from biopipelines.config_manager import ConfigManager
    from biopipelines.mock import Mock
    from biopipelines.cabsflex import CABSflex

    from biopipelines.pipeline import Resources

    with patch.object(ConfigManager, "get_scheduler", lambda self: "slurm"):
        pipeline = new_pipeline("cabsflex_aa_ok")
        with pipeline:
            Resources(time="1:00:00")
            m = Mock(ids=["s1"],
                     streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
            CABSflex(structures=m, num_models=2, aa_rebuild=True)
            script_path = pipeline.save()
        # Built without raising; the -A flag is emitted.
        assert "-A" in _tool_script(script_path, "CABSflex")
