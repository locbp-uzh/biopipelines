"""AiZynthFinder parameter validation and result-normalization regressions."""

import importlib.util
import os

import pytest

from biopipelines.aizynthfinder import _as_name_list

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HELPER = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_aizynthfinder.py")


def _load_helper():
    spec = importlib.util.spec_from_file_location("pipe_aizynthfinder", HELPER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("param", ["stocks", "expansion_policy", "filter_policy"])
def test_bare_string_rejected_instead_of_split_into_characters(param):
    with pytest.raises(ValueError, match="must be a list of names, not a string"):
        _as_name_list(param, "zinc")


def test_non_sequence_rejected():
    with pytest.raises(ValueError, match="must be a list or tuple"):
        _as_name_list("stocks", 5)


def test_lists_and_none_pass_through():
    assert _as_name_list("stocks", ["zinc", "emol"]) == ["zinc", "emol"]
    assert _as_name_list("stocks", ("zinc",)) == ["zinc"]
    assert _as_name_list("stocks", None) is None
    assert _as_name_list("stocks", []) is None


def test_missing_flag_is_not_in_stock():
    """bool(nan) is True, so an absent in_stock would mark a route solved."""
    as_bool = _load_helper()._as_bool
    assert as_bool(float("nan")) is False
    assert as_bool(None) is False


def test_present_flags_still_normalize():
    as_bool = _load_helper()._as_bool
    assert as_bool(True) is True
    assert as_bool("True") is True
    assert as_bool("no") is False
    assert as_bool(0) is False
