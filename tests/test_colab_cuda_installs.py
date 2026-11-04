"""The Colab installs of ESMFold and Gnina against a Colab image that moved to CUDA 13 and gcc 13.

Found by running the tools on a Colab T4 on 2026-11-04. The image carries nvcc 13.0, gcc 13.3 and the CUDA 13 libraries:

* `ESMFold.install()` failed: the sokrypton openfold fork compiles a CUDA kernel whenever nvcc is visible, and torch
  2.1.2+cu121 refuses nvcc 13.0 ("detected CUDA version mismatches"); with nvcc 12.1 the host gcc 13 headers then fail
  ("_Float32 is undefined"). The Colab env now carries its own nvcc 12.1 and gcc 11, and the install points CUDA_HOME, CC
  and CXX at them. The recipe's comment, "no nvcc compile", had stopped being true.
* `Gnina.install()` reported COMPLETED although the prebuilt binary did not load (it links libcudart.so.12, libcufft.so.11,
  libcublas.so.12, libnvrtc.so.12 and libnvToolsExt.so.1, none of which the image has). The install now puts the CUDA 12
  wheels in the gnina env, the run script puts them on LD_LIBRARY_PATH, and the verification runs the binary.

Script-text checks, not executions: the live runs are in the changelog entries.
"""
from __future__ import annotations

import inspect
import os
from unittest.mock import patch

import yaml

from biopipelines.config_manager import ConfigManager
from biopipelines.esmfold import ESMFold
from biopipelines.gnina import Gnina

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FOLDERS = {"biopipelines": "/repo", "Gnina": "/opt/gnina"}


def _install(cls, variant):
    with patch.object(ConfigManager, "get_scheduler", lambda self: variant), \
         patch.object(ConfigManager, "get_variant", lambda self: variant):
        return cls._install_script(FOLDERS, env_manager="micromamba")


class TestESMFold:

    def test_the_colab_install_points_the_build_at_the_envs_own_cuda_12_toolchain(self):
        script = _install(ESMFold, "colab")
        assert 'export CUDA_HOME="$ESMFOLD_PREFIX"' in script
        assert 'export CC="$ESMFOLD_PREFIX/bin/x86_64-conda-linux-gnu-gcc"' in script
        assert 'export CXX="$ESMFOLD_PREFIX/bin/x86_64-conda-linux-gnu-g++"' in script
        assert script.index("export CUDA_HOME") < script.index("pip install"), "the exports must precede the pip phases"

    def test_the_cluster_install_does_not_touch_the_toolchain(self):
        script = _install(ESMFold, "cluster")
        assert "CUDA_HOME" not in script and "x86_64-conda-linux-gnu" not in script

    def test_the_install_check_imports_openfold_too(self):
        assert "import esm, openfold" in _install(ESMFold, "colab")

    def test_the_colab_env_carries_nvcc_12_1_and_gcc_11(self):
        with open(os.path.join(REPO_ROOT, "environments", "esmfold.colab.yaml"), encoding="utf-8") as handle:
            spec = yaml.safe_load(handle)
        deps = " ".join(d for d in spec["dependencies"] if isinstance(d, str))
        assert "cuda-nvcc=12.1" in deps and "gxx_linux-64=11" in deps and "gcc_linux-64=11" in deps
        assert any("cuda-12.1" in channel for channel in spec["channels"])


class TestGnina:

    def test_the_colab_install_adds_the_cuda_12_wheels_to_the_env(self):
        script = _install(Gnina, "colab")
        for wheel in ("nvidia-cuda-runtime-cu12", "nvidia-cufft-cu12", "nvidia-cublas-cu12",
                      "nvidia-cuda-nvrtc-cu12", "nvidia-nvtx-cu12"):
            assert wheel in script, wheel
        with patch.object(ConfigManager, "get_scheduler", lambda self: "colab"),              patch.object(ConfigManager, "get_variant", lambda self: "colab"):
            env = Gnina._install_env("micromamba")
        assert f"micromamba run -n {env} pip install" in script, "the env the config maps Gnina to, not a fixed name"

    def test_the_colab_install_runs_the_binary_before_it_calls_the_install_good(self):
        script = _install(Gnina, "colab")
        assert "--version" in script
        assert script.index("--version") < script.index('touch "$INSTALL_SUCCESS"', script.index("--version"))
        assert "binary does not load" in script

    def test_the_cluster_install_keeps_its_old_shape(self):
        script = _install(Gnina, "cluster")
        assert "nvidia-" not in script and "--version" not in script

    def test_the_run_script_puts_the_wheels_on_the_library_path_on_colab_only(self):
        block = Gnina._colab_ld_library_block()
        assert "site-packages/nvidia/*/lib" in block and "LD_LIBRARY_PATH" in block
        source = inspect.getsource(Gnina.generate_script)
        assert "_colab_ld_library_block" in source and '== "colab"' in source
