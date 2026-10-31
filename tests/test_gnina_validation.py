"""Gnina's configuration-time validation of `cnn_scoring` and of the docking-only parameters.

`cnn_scoring` was documented with `"rescore_only"`, which GNINA does not accept (upstream's values are none,
rescore, refinement, metrorescore, metrorefine and all), and nothing checked it: a user who followed the docs
queued a job whose every pair failed on the compute node with GNINA's argument error. `_validate_no_search_mode`
compared the docking-only parameters against copies of the signature defaults, the exact pattern that once stopped
Vina from constructing, so a changed default would have refused every `mode="score"` step.
"""

import inspect

import pytest

from biopipelines.gnina import Gnina

from tests.test_gnina_ids import _declared_stream

UPSTREAM_CNN_SCORING = ("none", "rescore", "refinement", "metrorescore", "metrorefine", "all")


def test_the_accepted_values_are_the_ones_gnina_documents():
    assert tuple(Gnina.CNN_SCORING) == UPSTREAM_CNN_SCORING


@pytest.mark.parametrize("value", UPSTREAM_CNN_SCORING)
def test_every_gnina_value_is_accepted(local_config, isolated_cwd, new_pipeline, value):
    # `rank_conformers_by` stays "auto", so only the value itself is under test.
    _declared_stream(new_pipeline, ["p1"], ["l1"], cnn_scoring=value)


@pytest.mark.parametrize("value", ["rescore_only", "Rescore", "", None, 1])
def test_a_value_gnina_does_not_accept_is_refused_with_the_valid_ones(
        local_config, isolated_cwd, new_pipeline, value):
    with pytest.raises(ValueError) as refused:
        _declared_stream(new_pipeline, ["p1"], ["l1"], cnn_scoring=value)
    message = str(refused.value)
    assert "cnn_scoring" in message
    for valid in UPSTREAM_CNN_SCORING:
        assert valid in message, f"the message must name {valid!r}: {message}"


def test_the_docstring_and_the_docs_name_the_real_values_only():
    docstring = inspect.getsource(Gnina.__init__)
    assert "rescore_only" not in docstring
    import pathlib
    docs = (pathlib.Path(__file__).resolve().parent.parent / "docs" / "tool" / "structure_prediction.md").read_text(
        encoding="utf-8")
    assert "rescore_only" not in docs
    for valid in UPSTREAM_CNN_SCORING:
        assert f'"{valid}"' in docs, f"docs/tool/structure_prediction.md does not list {valid!r}"


def _score_mode_gnina(**overrides):
    """A Gnina carrying only the attributes the no-search validator reads, all at their signature defaults."""
    parameters = inspect.signature(Gnina.__init__).parameters
    gnina = Gnina.__new__(Gnina)
    for name, parameter in parameters.items():
        if name not in ("self", "structures", "compounds", "kwargs"):
            setattr(gnina, name, parameter.default)
    gnina.mode = "score"
    gnina.autobox_ligand_stream = None
    gnina.autobox_ligand_path = None
    for key, value in overrides.items():
        setattr(gnina, key, value)
    return gnina


def test_defaults_are_accepted_in_score_mode():
    Gnina._validate_no_search_mode(_score_mode_gnina())


def test_a_docking_only_parameter_set_in_score_mode_is_still_refused():
    with pytest.raises(ValueError, match="exhaustiveness"):
        Gnina._validate_no_search_mode(_score_mode_gnina(exhaustiveness=32))


def test_a_changed_signature_default_does_not_refuse_score_mode(monkeypatch):
    """The regression: a default moved in the signature must not turn the inherited value into an offender."""
    signature = inspect.signature(Gnina.__init__)
    moved = signature.replace(parameters=[
        parameter.replace(default=16) if name == "exhaustiveness" else parameter
        for name, parameter in signature.parameters.items()])
    monkeypatch.setattr(Gnina.__init__, "__signature__", moved, raising=False)
    gnina = _score_mode_gnina()
    assert gnina.exhaustiveness == 16
    Gnina._validate_no_search_mode(gnina)
