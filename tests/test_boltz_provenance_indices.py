# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Boltz post-processing must not own an id-parsing rule.

`_prov_index` related a per-sample structure id to the '<base>' key combinatorics uses for provenance by doing `rpartition('_')` and requiring the suffix to be digits. That is the framework's parent tier restated inside one tool, and it diverged in two ways: a non-digit suffix silently lost its provenance columns, and a lookup that resolved through a deeper tier could not be reported, because only `id_map_utils` scores tiers.
"""

import importlib.util
import os

import pytest


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_boltz_postprocessing.py")


def _provenance_indices():
    """Import just the helper. The module body runs a whole post-processing pass, so it is read and the one function compiled out of it."""
    import ast

    source = open(SCRIPT, encoding="utf-8").read()
    tree = ast.parse(source)
    fn = next(n for n in tree.body
              if isinstance(n, ast.FunctionDef) and n.name == "_provenance_indices")
    namespace = {}
    from biopipelines.id_map_utils import get_mapped_ids
    namespace["get_mapped_ids"] = get_mapped_ids
    exec(compile(ast.Module(body=[fn], type_ignores=[]), SCRIPT, "exec"), namespace)
    return namespace["_provenance_indices"]


def test_the_local_id_parsing_helper_is_gone():
    source = open(SCRIPT, encoding="utf-8").read()
    assert "_prov_index" not in source
    assert "rpartition" not in source, "the parent-suffix rule belongs to id_map_utils"


def test_an_exact_id_finds_its_own_position():
    resolve = _provenance_indices()
    assert resolve(["p+l2", "p+l1"], ["p+l1", "p+l2"]) == {"p+l2": 1, "p+l1": 0}


def test_a_per_sample_suffix_resolves_to_its_parent():
    """With top_only=False boltz emits '<base>_<k>' while provenance is keyed by '<base>'."""
    resolve = _provenance_indices()
    assert resolve(["p+l1_0", "p+l2_3"], ["p+l1", "p+l2"]) == {"p+l1_0": 0, "p+l2_3": 1}


def test_a_non_digit_suffix_also_resolves():
    """The local copy required `isdigit()`, so this row silently got empty provenance columns."""
    resolve = _provenance_indices()
    assert resolve(["p+l1_top"], ["p+l1", "p+l2"]) == {"p+l1_top": 0}


def test_an_unrelated_id_resolves_to_nothing():
    resolve = _provenance_indices()
    assert resolve(["somethingelse"], ["p+l1", "p+l2"]) == {"somethingelse": None}
