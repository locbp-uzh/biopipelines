# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""A forwarded key must be able to spell the upstream flag it stands for.

The key allowlist accepted an identifier or a dotted path, so a dash-spelled flag had no expressible form at all: `--msa-mode`, `--weighted-fit`, `--num-steps` could not be forwarded to any tool. That silently limited forwarding to whichever upstream parsers happen to use underscores -- NeuralPLexer's surface is almost entirely dashed, so its marker bought very little. None of these spellings is a Python identifier, so they all arrive the same way, as `Tool(**{"msa-mode": "single"})`.
"""

import pytest

from biopipelines.base_config import _validate_extra_arg, render_extra_args


@pytest.mark.parametrize("key", [
    "msa-mode", "weighted-fit", "num-steps",          # dash-spelled upstream flags
    "msa_mode", "num_recycle",                        # underscore-spelled ones still work
    "+potentials.substrate", "denoiser.noise_scale_ca",  # hydra overrides still work
])
def test_a_forwarded_key_may_spell_its_upstream_flag(key):
    _validate_extra_arg(key, "value")


@pytest.mark.parametrize("key", [
    "-x", "--msa-mode",   # a leading dash would let a forwarded key impersonate one of our own flags
    "-", "",              # nothing to name
    "a b", "a;b", "a|b",  # a key becomes a flag verbatim, so it may not carry shell syntax
])
def test_a_key_that_could_be_mistaken_for_something_else_is_refused(key):
    with pytest.raises(ValueError):
        _validate_extra_arg(key, "value")


def test_a_dashed_key_renders_as_the_flag_it_names():
    assert render_extra_args({"msa-mode": "single"}, "argparse") == ["--msa-mode", "single"]


def test_a_dashed_key_still_renders_as_a_bare_flag_when_true():
    assert render_extra_args({"weighted-fit": True}, "argparse") == ["--weighted-fit"]
