"""The foundry env serves RFdiffusion3, RF3 and FoundryMPNN, and its install must prove they import.

An unpinned install resolved atomworks 3.x, where `rfd3.engine` and `rf3.inference_engines.rf3` no longer import, while the old checks (`import foundry`, `import rf3`, `import mpnn`) still passed.
"""

import os

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_the_env_spec_keeps_atomworks_below_3():
    spec = open(os.path.join(REPO_ROOT, "environments", "foundry.pip.txt"), encoding="utf-8").read()
    lines = [ln.strip() for ln in spec.splitlines() if ln.strip() and not ln.startswith("#")]
    assert "atomworks[ml]<3" in lines


def test_each_install_proves_that_its_engine_imports():
    from biopipelines.foundry_mpnn import FoundryMPNN
    from biopipelines.rf3 import RF3
    from biopipelines.rfdiffusion3 import RFdiffusion3

    folders = {"biopipelines": "/repo", "RF3": "/ck/rf3", "FoundryMPNN": "/ck/mpnn", "RFdiffusion3": "/ck/rfd3"}
    expected = {RFdiffusion3: "import rfd3.engine", RF3: "import rf3.inference_engines.rf3", FoundryMPNN: "import mpnn.inference_engines.mpnn"}
    for tool, statement in expected.items():
        script = tool._install_script(folders, env_manager="mamba")
        assert f'python -c "{statement}"' in script, f"{tool.__name__} does not verify with {statement}"
