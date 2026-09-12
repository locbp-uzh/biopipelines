# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Guards the data_containers -> outputs -> base_config layering created when base_config.py was split.

The split is only worth having if it stays a layering: these tests fail if a future edit re-couples the
layers, drops a re-export that ~90 modules depend on, or lets __all__ promise a name the package lacks.
"""

import ast
import os

import pytest

import biopipelines
from biopipelines import base_config, data_containers, outputs


PACKAGE_DIR = os.path.dirname(os.path.abspath(biopipelines.__file__))

MOVED_TO_DATA_CONTAINERS = [
    "TableMetadata",
    "TableInfo",
    "IndexedTableContainer",
    "TableContainer",
    "StreamContainer",
    "resolve_table_reference",
]

MOVED_TO_OUTPUTS = ["StandardizedOutput", "ToolOutput"]


def _module_source(name):
    with open(os.path.join(PACKAGE_DIR, f"{name}.py"), encoding="utf-8") as handle:
        return handle.read()


def _imported_module_paths(name):
    """Every module path named by an import anywhere in the file, functions included.

    Both statement forms and every dotted segment, because `import biopipelines.base_config` and
    `from biopipelines.outputs import x` each hide the interesting name past the first segment.
    """
    found = set()
    for node in ast.walk(ast.parse(_module_source(name))):
        paths = []
        if isinstance(node, ast.Import):
            paths = [alias.name for alias in node.names]
        elif isinstance(node, ast.ImportFrom):
            # A relative `from . import x` has no module; the name itself is the module.
            paths = [node.module] if node.module else [alias.name for alias in node.names]
        for path in paths:
            if path:
                found.update(path.split("."))
    return found


@pytest.mark.parametrize("name", MOVED_TO_DATA_CONTAINERS)
def test_data_container_names_are_defined_in_data_containers(name):
    assert hasattr(data_containers, name)


@pytest.mark.parametrize("name", MOVED_TO_OUTPUTS)
def test_output_names_are_defined_in_outputs(name):
    assert hasattr(outputs, name)


@pytest.mark.parametrize("name", MOVED_TO_DATA_CONTAINERS + MOVED_TO_OUTPUTS)
def test_base_config_reexports_the_same_object(name):
    # ~90 modules import these from base_config; the re-export must be the identical object, not a copy.
    origin = data_containers if name in MOVED_TO_DATA_CONTAINERS else outputs
    assert getattr(base_config, name) is getattr(origin, name)


# outputs -> base_config is excluded here on purpose: it is legal under TYPE_CHECKING, and the
# next test governs that pair precisely.
@pytest.mark.parametrize("lower, higher", [
    ("data_containers", "outputs"),
    ("data_containers", "base_config"),
])
def test_lower_layer_never_imports_higher_layer(lower, higher):
    # Any import statement anywhere in the file counts, including one deferred inside a method.
    assert higher not in _imported_module_paths(lower)


def test_outputs_names_base_config_only_under_type_checking():
    # The `config: 'BaseConfig'` annotation must stay quoted, so the lower layer never imports the higher.
    tree = ast.parse(_module_source("outputs"))
    type_checking_nodes = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.If) and "TYPE_CHECKING" in ast.dump(node.test):
            type_checking_nodes.update(id(child) for child in ast.walk(node))
    for node in ast.walk(tree):
        names = []
        if isinstance(node, ast.Import):
            names = [alias.name for alias in node.names]
        elif isinstance(node, ast.ImportFrom) and node.module:
            names = [node.module]
        if any("base_config" in n.split(".") for n in names):
            assert id(node) in type_checking_nodes, (
                f"outputs.py line {node.lineno} imports base_config outside TYPE_CHECKING")


def test_package_all_has_no_phantom_names():
    missing = [name for name in biopipelines.__all__ if not hasattr(biopipelines, name)]
    assert missing == []
