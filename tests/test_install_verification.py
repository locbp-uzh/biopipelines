"""An install may only report success on an env that imports every module the tool's pipe scripts import.

The expected modules come from an `ast` scan of the pipe scripts, so adding an import to a script without adding it to the install check fails here. A check that imported fewer modules passed on envs that could not run the tool: ESM2 verified torch and transformers while its pipe script also needs numpy, pandas and pyyaml.
"""

import ast
import functools
import importlib
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
FIRST_PARTY_FOLDERS = (REPO_ROOT / "biopipelines", REPO_ROOT / "pipe_scripts")
FOLDERS = {"biopipelines": "/repo", "pipe_scripts": "/repo/pipe_scripts"}

# Scripts a wrapper names that do not run in the tool's own env: the host env, or the import-free weights fetch.
RUNS_ELSEWHERE = {
    "Protenix": {"pipe_protenix_weights.py"},
    "RFD4Proteina": {"pipe_rfd4proteina_postprocess.py"},
    "ESMFold2": {"pipe_esmfold2_postprocessing.py"},
}

# tool -> (wrapper modules that name its pipe scripts, extra folders its install script needs)
TOOLS = {
    "ESM2": (("esm2",), {"ESM2": "/ck/esm2"}),
    "Peppr": (("peppr",), {}),
    "Foldseek": (("foldseek",), {}),
    "AiZynthFinder": (("aizynthfinder",), {"AiZynthFinder": "/ck/aizynthfinder"}),
    "REINVENT4": (("reinvent4",), {"REINVENT4": "/ck/reinvent4"}),
    "FLOWR": (("flowr",), {"FLOWR": "/ck/flowr"}),
    "PPIDdG": (("ppi_ddg",), {}),
    "LigandDdG": (("ligand_ddg",), {}),
    "Protenix": (("protenix",), {"ProtenixCache": "/cache/protenix"}),
    "RFD4Proteina": (("rfd4proteina",), {"RFD4Proteina": "/ck/rfd4"}),
    "BioEmu2": (("bioemu2",), {}),
    "ESMFold2": (("esmfold2",), {}),
    "IPSAE": (("ipsae",), {"IPSAE": "/ck/ipsae"}),
    "Embedding": (("embedding",), {}),
    "USalign": (("usalign", "pair_alignment"), {}),
    "TMalign": (("tmalign", "pair_alignment"), {}),
    "Foldcomp": (("foldcomp",), {}),
}

PYTHON_PAYLOAD = re.compile(r'python -c "([^"]*)"')


@functools.lru_cache(maxsize=None)
def parse_imports(path):
    """(module, level, names, inside_function) for every import statement of a file."""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    found = []

    def walk(node, inside):
        for child in ast.iter_child_nodes(node):
            deeper = inside or isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef))
            if isinstance(child, ast.Import):
                found.extend((alias.name, 0, (), deeper) for alias in child.names)
            elif isinstance(child, ast.ImportFrom):
                found.append((child.module or "", child.level, tuple(a.name for a in child.names), deeper))
            walk(child, deeper)

    walk(tree, False)
    return tuple(found)


def first_party_files(module, level, names, importer):
    """The files a first-party import loads, or None when the import is not first-party."""
    if level:
        base = importer.parent
        if not module:
            return [base / "__init__.py"] + [base / f"{n}.py" for n in names if (base / f"{n}.py").exists()]
        target = base / module.replace(".", "/")
        return [p for p in (target.with_suffix(".py"), target / "__init__.py") if p.exists()]
    top, *rest = module.split(".")
    if top == "biopipelines":
        files = [REPO_ROOT / "biopipelines" / "__init__.py"]
        if rest and (REPO_ROOT / "biopipelines" / f"{rest[0]}.py").exists():
            files.append(REPO_ROOT / "biopipelines" / f"{rest[0]}.py")
        return files
    for folder in FIRST_PARTY_FOLDERS:
        if (folder / f"{top}.py").exists():
            return [folder / f"{top}.py"]
    return None


def third_party_modules(path, whole_file=True, seen=None):
    """Top-level third-party modules a file needs: all of its own imports, then the module-level ones of the first-party files it loads."""
    seen = set() if seen is None else seen
    found = set()
    for module, level, names, inside in parse_imports(path):
        if inside and not whole_file:
            continue
        files = first_party_files(module, level, names, path)
        if files is None:
            top = module.split(".")[0]
            if top and top not in sys.stdlib_module_names:
                found.add(top)
            continue
        for file in files:
            if file not in seen:
                seen.add(file)
                found |= third_party_modules(file, False, seen)
    return found


def pipe_scripts_of(tool, wrapper_modules):
    """The pipe scripts a tool's wrapper modules name, minus those that run outside its env."""
    named = set()
    for module in wrapper_modules:
        text = (REPO_ROOT / "biopipelines" / f"{module}.py").read_text(encoding="utf-8")
        named |= set(re.findall(r"pipe_\w+\.py", text))
    existing = {n for n in named if (REPO_ROOT / "pipe_scripts" / n).exists()}
    return sorted(existing - RUNS_ELSEWHERE.get(tool, set()))


def expected_modules(tool, wrapper_modules):
    modules = set()
    for name in pipe_scripts_of(tool, wrapper_modules):
        modules |= third_party_modules(REPO_ROOT / "pipe_scripts" / name)
    return modules


def success_conditions(script):
    """The `if` condition that guards each `touch "$INSTALL_SUCCESS"`; a touch outside a `then` branch is an error."""
    lines = script.splitlines()
    stack, conditions, i = [], [], 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("if "):
            condition = line
            while not condition.endswith("; then") and i + 1 < len(lines):
                i += 1
                condition = condition.rstrip("\\").rstrip() + " " + lines[i].strip()
            stack.append([condition, "then"])
        elif line == "else" and stack:
            stack[-1][1] = "else"
        elif line == "fi" and stack:
            stack.pop()
        elif line == 'touch "$INSTALL_SUCCESS"':
            assert stack and stack[-1][1] == "then", "success is reported outside any verified branch"
            conditions.append(stack[-1][0])
        i += 1
    return conditions


def checked_modules(condition):
    """Top-level modules a condition imports in a `python -c` payload, or looks up with `find_spec` without importing."""
    names = set()
    for payload in PYTHON_PAYLOAD.findall(condition):
        for node in ast.walk(ast.parse(payload)):
            if isinstance(node, ast.Import):
                names |= {alias.name.split(".")[0] for alias in node.names}
            elif isinstance(node, ast.ImportFrom) and node.module:
                names.add(node.module.split(".")[0])
            elif isinstance(node, ast.Call) and node.args and isinstance(node.args[0], ast.Constant):
                function = node.func.attr if isinstance(node.func, ast.Attribute) else getattr(node.func, "id", "")
                if function == "find_spec":
                    names.add(str(node.args[0].value).split(".")[0])
    return names


def install_script(name, force_reinstall=False):
    _, extra_folders = TOOLS[name]
    tool = getattr(importlib.import_module("biopipelines"), name)
    return tool._install_script({**FOLDERS, **extra_folders}, env_manager="mamba", force_reinstall=force_reinstall)


@pytest.mark.parametrize("name", list(TOOLS))
def test_every_success_path_imports_what_the_pipe_scripts_import(name):
    wrapper_modules, _ = TOOLS[name]
    script = install_script(name)
    expected = expected_modules(name, wrapper_modules)
    assert expected, f"the scan found no module for {name}"
    conditions = success_conditions(script)
    assert conditions, f"{name}: the install script never reports success"
    for condition in conditions:
        missing = expected - checked_modules(condition)
        assert not missing, f"{name}: a success path does not import {sorted(missing)}:\n{condition}"


@pytest.mark.parametrize("force", [False, True])
@pytest.mark.parametrize("name", list(TOOLS))
def test_the_install_script_is_valid_bash(name, force, tmp_path):
    bash = shutil.which("bash")
    if bash is None:
        pytest.skip("bash is not available")
    path = tmp_path / "install.sh"
    path.write_text(install_script(name, force), encoding="utf-8", newline="\n")
    result = subprocess.run([bash, "-n", str(path)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("name", list(TOOLS))
def test_the_scripts_scanned_for_a_tool_are_the_ones_that_exist(name):
    wrapper_modules, _ = TOOLS[name]
    assert pipe_scripts_of(name, wrapper_modules), f"{name} names no pipe script that runs in its env"


def test_the_scan_reads_what_the_esm2_scripts_import():
    expected = expected_modules("ESM2", ("esm2",))
    assert {"numpy", "pandas", "yaml", "torch", "transformers", "huggingface_hub"} <= expected
    assert not expected & {"os", "sys", "argparse", "biopipelines"}


def test_a_module_loaded_only_through_a_first_party_helper_is_found():
    assert "pandas" in third_party_modules(REPO_ROOT / "pipe_scripts" / "pipe_ppi_ddg.py")


def test_a_check_that_lacks_a_module_is_reported():
    old_esm2_check = (
        'if mamba run -n esm2 python -c "import torch, transformers, huggingface_hub" >/dev/null 2>&1 '
        '&& [ -f "/w/REVISION" ]; then'
    )
    missing = expected_modules("ESM2", ("esm2",)) - checked_modules(old_esm2_check)
    assert missing == {"numpy", "pandas", "yaml"}


def test_a_lookup_without_importing_counts_as_checked():
    condition = 'if mamba run -n e python -c "import importlib.util as u; assert u.find_spec(\'atomworks\') and u.find_spec(\'rfproteina\')" >/dev/null 2>&1; then'
    assert checked_modules(condition) >= {"atomworks", "rfproteina", "importlib"}


def test_success_outside_a_verified_branch_is_refused():
    with pytest.raises(AssertionError, match="outside any verified branch"):
        success_conditions('if [ -f x ]; then\n    echo hi\nelse\n    touch "$INSTALL_SUCCESS"\nfi\n')
    with pytest.raises(AssertionError, match="outside any verified branch"):
        success_conditions('touch "$INSTALL_SUCCESS"\n')


def test_a_condition_on_continued_lines_is_read_whole():
    script = 'if [ -s a ] \\\n   && python -c "import numpy"; then\n    touch "$INSTALL_SUCCESS"\nfi\n'
    assert checked_modules(success_conditions(script)[0]) == {"numpy"}
