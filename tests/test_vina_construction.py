"""Vina must construct with its inherited defaults.

_validate_vina_backend rejects Gnina's CNN parameters, which it detected by
comparing against hard-coded copies of Gnina's defaults. When 85d3134 changed
cnn_score_threshold's default from 0.5 to 0.0, the inherited value no longer
matched the copy and every Vina(...) raised at construction -- the tool was
unusable and nothing caught it, because no test instantiated it.
"""

import inspect

import pytest

from biopipelines.gnina import Gnina, Vina


def _gnina_default(name):
    return inspect.signature(Gnina.__init__).parameters[name].default


def _vina_with(**overrides):
    """A Vina carrying only the attributes the backend validator reads."""
    v = Vina.__new__(Vina)
    v.VINA_SCORING = Vina.VINA_SCORING
    v.scoring = "vina"
    v.cnn_scoring = _gnina_default("cnn_scoring")
    v.cnn_score_threshold = _gnina_default("cnn_score_threshold")
    for key, value in overrides.items():
        setattr(v, key, value)
    return v


def test_inherited_cnn_defaults_are_accepted():
    Vina._validate_vina_backend(_vina_with())


@pytest.mark.parametrize("name", ["cnn_scoring", "cnn_score_threshold"])
def test_cnn_parameters_are_still_rejected_when_set(name):
    sentinel = {"cnn_scoring": "none", "cnn_score_threshold": 0.7}[name]
    assert sentinel != _gnina_default(name), "sentinel must differ from the default"
    with pytest.raises(ValueError, match="no CNN scoring"):
        Vina._validate_vina_backend(_vina_with(**{name: sentinel}))


def test_validator_reads_defaults_from_gnina_not_from_a_copy():
    """The regression itself: a changed Gnina default must not break Vina.

    Asserting on the mechanism rather than on today's values, so this keeps
    holding whichever way the defaults move next.
    """
    src = inspect.getsource(Vina._validate_vina_backend)
    assert "signature(Gnina.__init__)" in src, (
        "the validator must derive the CNN defaults from Gnina's signature; "
        "hard-coded copies are what broke construction in the first place")


def test_scoring_enum_is_enforced():
    with pytest.raises(ValueError, match="scoring must be one of"):
        Vina._validate_vina_backend(_vina_with(scoring="rescore"))
