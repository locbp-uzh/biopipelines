"""The foundry env serves RFdiffusion3, RF3 and FoundryMPNN, and its install must prove they import.

An unpinned install resolved atomworks 3.x, where `rfd3.engine` and `rf3.inference_engines.rf3` no longer import, while the old checks (`import foundry`, `import rf3`, `import mpnn`) still passed. The pin only helps an env built after it, so the "already installed" shortcut has to run the engine import as well, and an env that fails it must have the spec applied again.
"""

import os
import shutil
import subprocess

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# tool name -> (engine import its install must run, files that count as its weights)
TOOLS = {
    "RFdiffusion3": ("import rfd3.engine", ["rfd3_latest.ckpt"]),
    "RF3": ("import rf3.inference_engines.rf3", ["rf3_latest.ckpt"]),
    "FoundryMPNN": ("import mpnn.inference_engines.mpnn", ["proteinmpnn_v1.pt", "ligandmpnn_v1.pt"]),
}

# Behaves as a `mamba` whose foundry env is stale while $FAKE_STALE exists, and that a `pip install` repairs.
FAKE_MAMBA = r'''
mamba() {
    echo "$*" >> "$FAKE_LOG"
    if [ "$1" = run ]; then
        shift 3
        tool="$1"; shift
        if [ "$tool" = pip ]; then rm -f "$FAKE_STALE"; return 0; fi
        if [ "$tool" = python ] && [ "$1" = -c ] && [ -n "$2" ] && [ -e "$FAKE_STALE" ]; then return 1; fi
    fi
    return 0
}
'''


def test_the_env_spec_keeps_atomworks_below_3():
    spec = open(os.path.join(REPO_ROOT, "environments", "foundry.pip.txt"), encoding="utf-8").read()
    lines = [ln.strip() for ln in spec.splitlines() if ln.strip() and not ln.startswith("#")]
    assert "atomworks[ml]<3" in lines


def _tool(name):
    import importlib

    return getattr(importlib.import_module("biopipelines"), name)


def _install_script(name, weights, force_reinstall=False):
    folders = {"biopipelines": "/repo", "RF3": f"{weights}/rf3", "FoundryMPNN": f"{weights}/mpnn", "RFdiffusion3": f"{weights}/rfd3"}
    return _tool(name)._install_script(folders, env_manager="mamba", force_reinstall=force_reinstall)


@pytest.mark.parametrize("name", TOOLS)
def test_each_install_proves_that_its_engine_imports(name, shipped_config):
    shipped_config("cluster")
    statement = TOOLS[name][0]
    script = _install_script(name, "/ck")
    assert f'python -c "{statement}"' in script, f"{name} does not verify with {statement}"


@pytest.mark.parametrize("name", TOOLS)
def test_the_already_installed_shortcut_runs_the_engine_import(name, shipped_config):
    """An env built before the pin has its weights and its env, and used to pass this check."""
    shipped_config("cluster")
    script = _install_script(name, "/ck")
    shortcut = script[script.index("# Check if already"):script.index("\nfi\n")]
    assert f'python -c "{TOOLS[name][0]}"' in shortcut


@pytest.mark.parametrize("name", TOOLS)
def test_an_env_that_fails_the_engine_import_gets_the_spec_applied_again(name, shipped_config):
    shipped_config("cluster")
    script = _install_script(name, "/ck")
    repair = script[script.index("elif ! "):script.index("else\n", script.index("elif ! "))]
    assert f'python -c "{TOOLS[name][0]}"' in repair
    assert 'pip install -r "/repo/environments/foundry.pip.txt"' in repair


@pytest.mark.skipif(shutil.which("bash") is None, reason="bash is not available")
@pytest.mark.parametrize("stale", [True, False])
@pytest.mark.parametrize("name", TOOLS)
def test_a_stale_env_is_repaired_and_a_good_one_is_left_alone(name, stale, shipped_config, tmp_path):
    shipped_config("cluster")
    folder = {"RFdiffusion3": "rfd3", "RF3": "rf3", "FoundryMPNN": "mpnn"}[name]
    weights = tmp_path / folder
    weights.mkdir()
    for file in TOOLS[name][1]:
        (weights / file).write_text("weights", encoding="utf-8")
    log, success, marker = tmp_path / "log.txt", tmp_path / "success", tmp_path / "stale"
    if stale:
        marker.write_text("", encoding="utf-8")
    script = "\n".join([
        f'export INSTALL_SUCCESS="{success.as_posix()}" FAKE_LOG="{log.as_posix()}" FAKE_STALE="{marker.as_posix()}"',
        FAKE_MAMBA,
        "(",
        _install_script(name, tmp_path.as_posix()),
        ")",
    ])
    result = subprocess.run([shutil.which("bash"), "-c", script], capture_output=True, text=True)
    commands = log.read_text(encoding="utf-8") if log.exists() else ""
    assert success.exists(), result.stdout + result.stderr
    assert ("pip install -r /repo/environments/foundry.pip.txt" in commands) is stale, commands
    assert "env create" not in commands
