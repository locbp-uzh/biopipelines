# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Boltz2(single_sequence=True) reaches the config generator and sets `msa: empty` on every protein chain.

The wrapper emitted `--single-sequence` before the config script declared it, so argparse rejected the invocation, no YAML was written, and Boltz2 folded nothing while post-processing reported success.
These tests pin that the flag parses and that it, not the MSA backfill, is what writes `msa: empty`: the backfill acts only when some other chain has an alignment, so without the flag nothing marks a chain single-sequence.

Not covered: that Boltz accepts the resulting YAML, which needs the Boltz environment.
"""

import importlib.util
import os
import sys
import types

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_boltz_config_unified.py")

PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"
AXES = {"proteins": {"entity_type": "protein", "mode": "each", "static": [], "static_first": False,
                     "iterated": [{"id": "p1", "sequence": PROTEIN}, {"id": "p2", "sequence": "MKTAYIAKQRQ"}]}}


@pytest.fixture
def script(monkeypatch):
    monkeypatch.syspath_prepend(os.path.join(REPO_ROOT, "pipe_scripts"))
    spec = importlib.util.spec_from_file_location("pipe_boltz_config_unified", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _args(single_sequence):
    return types.SimpleNamespace(
        affinity=False, single_sequence=single_sequence, template=None, template_chains=None,
        glycosylation=None, covalent_linkage=None, contacts=None, disulfide_bonds=None, metal_coord=None,
        pocket_residues=None, pocket_binder=None, pocket_max_distance=None)


def test_the_parser_accepts_the_flag_the_wrapper_emits(script, monkeypatch):
    monkeypatch.setattr(sys, "argv", ["pipe_boltz_config_unified.py", "--combinatorics-config", "c.json",
                                      "--output-dir", "out", "--skipped-csv", "s.csv", "--single-sequence"])
    assert script.parse_arguments().single_sequence is True


def test_the_flag_is_off_unless_given(script, monkeypatch):
    monkeypatch.setattr(sys, "argv", ["pipe_boltz_config_unified.py", "--combinatorics-config", "c.json",
                                      "--output-dir", "out", "--skipped-csv", "s.csv"])
    assert script.parse_arguments().single_sequence is False


def test_every_protein_chain_is_marked_empty(script):
    configs, skipped = script.generate_configs(AXES, None, _args(True))
    assert not skipped and len(configs) == 2
    assert all(entry["protein"]["msa"] == "empty" for _, c in configs for entry in c["sequences"])


def test_without_the_flag_no_chain_is_marked(script):
    configs, _ = script.generate_configs(AXES, None, _args(False))
    assert all("msa" not in entry["protein"] for _, c in configs for entry in c["sequences"])
