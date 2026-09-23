# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""A handful of package modules are shared with the pipe scripts, which run on a compute node with no installed package: they put `biopipelines/` on `sys.path` and import the module by bare name. The module is then top-level, and `from .x import y` raises `ImportError: attempted relative import with no known parent package`.

The rest of the suite cannot see this, because a test imports `biopipelines` as a package, where the relative form works. That gap let a guardless `from .idset import ...` inside `combinatorics.predict_single_output_id` reach the cluster, where it produced zero Boltz2 configurations -- and Boltz2 then ran on zero inputs and raised nothing of its own.

Only the shared modules are checked. Wrappers and `pipeline.py` are configuration-time only, always imported as part of the package, and their relative imports are correct as they are.
"""

import ast
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PKG = ROOT / "biopipelines"
PIPE_SCRIPTS = ROOT / "pipe_scripts"

# Names a pipe script imports bare that are also modules here. Kept explicit rather than derived,
# because several bare imports in pipe_scripts (pymol, openmm, prolif, ...) are the external tool
# of that name and merely collide with a wrapper's filename.
SHARED_WITH_PIPE_SCRIPTS = ["combinatorics", "id_patterns", "ligand_utils", "nucleic_acids"]


def _intra_package_imports(path):
    """Module names this file imports from within the package, by either spelling."""
    names = set()
    tree = ast.parse(path.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module and node.level:
            names.add(node.module.split(".")[0])
        elif isinstance(node, ast.ImportFrom) and node.module and not node.level:
            head = node.module.split(".")[0]
            if head == "biopipelines" and len(node.module.split(".")) > 1:
                names.add(node.module.split(".")[1])
            elif (PKG / f"{head}.py").exists():
                names.add(head)
        elif isinstance(node, ast.Import):
            for alias in node.names:
                head = alias.name.split(".")[0]
                if (PKG / f"{head}.py").exists():
                    names.add(head)
    return {n for n in names if (PKG / f"{n}.py").exists()}


def _closure():
    """Every module reachable from the shared entry points, since a transitive import breaks the same way."""
    seen, stack = set(), list(SHARED_WITH_PIPE_SCRIPTS)
    while stack:
        name = stack.pop()
        if name in seen:
            continue
        seen.add(name)
        stack.extend(_intra_package_imports(PKG / f"{name}.py") - seen)
    return sorted(seen)


REACHABLE = _closure()


def test_the_entry_points_are_still_the_ones_pipe_scripts_import():
    """If a pipe script starts sharing another module, this list has to grow with it or the new module goes unchecked."""
    bare = set()
    for script in PIPE_SCRIPTS.glob("*.py"):
        tree = ast.parse(script.read_text(encoding="utf-8", errors="ignore"))
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.module and not node.level:
                head = node.module.split(".")[0]
                if (PKG / f"{head}.py").exists():
                    bare.add(head)
    missed = bare - set(REACHABLE)
    # Names that are really the external tool, not our wrapper of the same name.
    missed -= {"pymol", "openmm", "prolif", "plip", "posebusters", "admet_ai",
               "frame2seq", "bioemu", "openbabel", "utils"}
    assert not missed, f"pipe scripts import these package modules by bare name but they are unchecked: {sorted(missed)}"


def _run_standalone(code):
    """A fresh interpreter with only biopipelines/ on sys.path, exactly as a pipe script sets itself up."""
    script = f"import sys; sys.path.insert(0, {str(PKG)!r})\n" + textwrap.dedent(code)
    return subprocess.run([sys.executable, "-c", script], capture_output=True, text=True)


@pytest.mark.parametrize("module", REACHABLE)
def test_the_shared_modules_import_standalone(module):
    result = _run_standalone(f"import {module}")
    assert result.returncode == 0, f"{module} does not import with no package around it:\n{result.stderr}"


def test_combinatorics_predicts_an_id_when_loaded_standalone():
    """The call `pipe_boltz_config_unified.py` makes for every configuration it writes. It raised on the cluster, so the tool wrote none."""
    result = _run_standalone("""
        import combinatorics
        got = combinatorics.predict_single_output_id(
            proteins=("bundle", ["TrpR", "TrpR"], None, [], False),
            ligands=("each", ["lig1", "lig2", "lig3"], 1, [], False),
        )
        print("RESULT:" + str(got))
    """)
    assert result.returncode == 0, f"standalone call failed:\n{result.stderr}"
    assert "RESULT:" in result.stdout, result.stdout + result.stderr
