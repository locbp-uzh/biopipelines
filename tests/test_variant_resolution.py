"""Variant auto-detection refuses to guess, and says which variant it resolved.

Every shipped ``config.*.yaml`` carries ``machine.username: ""``. The old detector treated a blank username as a wildcard and returned the alphabetically first one, so on a laptop with no scheduler it silently loaded ``config.cluster.yaml`` (``mamba`` + ``slurm``) and the first symptom was ``RuntimeError: Resources() must be called before adding tools`` -- an error that names nothing to do with the config. Detection now raises with the variants it found and the two ways to choose one, and every resolution route prints one line naming what it picked.

The explicit routes must keep working untouched: the env var, ``Pipeline(config=...)``, ``ConfigManager(variant=...)``, and the Colab branch (a real runtime signal, so it still wins silently).
"""

import getpass
import os
import sys

import pytest

from biopipelines import config_manager
from biopipelines.config_manager import (
    ConfigManager,
    _autodetect_variant,
    _autodetect_variant_with_source,
)


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@pytest.fixture
def unclaimed_host(monkeypatch):
    """A host no shipped config claims: no env var, not Colab, unknown user."""
    monkeypatch.delenv("BIOPIPELINES_CONFIG_VARIANT", raising=False)
    monkeypatch.setattr(getpass, "getuser", lambda: "no-such-biopipelines-user")


@pytest.fixture
def fresh_announcement(monkeypatch):
    """Re-arm the once-per-process announcement so a test can observe it."""
    monkeypatch.setattr(ConfigManager, "_variant_announced", False)
    monkeypatch.setattr(ConfigManager, "_variant", None)
    monkeypatch.setattr(ConfigManager, "_config", None)
    monkeypatch.setattr(ConfigManager, "_instance", None)


# ── refusal ──────────────────────────────────────────────────────────────────

def test_unmatched_autodetect_raises_instead_of_guessing(unclaimed_host):
    with pytest.raises(RuntimeError):
        _autodetect_variant()


def test_refusal_names_what_it_looked_for(unclaimed_host, record_case):
    with pytest.raises(RuntimeError) as excinfo:
        _autodetect_variant()
    message = str(excinfo.value)
    record_case(input="autodetect on an unclaimed host", expected="raises",
                actual="raises")
    assert "machine.username" in message
    assert REPO_ROOT in message


def test_refusal_lists_the_variants_it_found(unclaimed_host):
    with pytest.raises(RuntimeError) as excinfo:
        _autodetect_variant()
    message = str(excinfo.value)
    for variant in ("cluster", "container", "daint", "local"):
        assert variant in message, f"refusal does not name the {variant!r} variant"


def test_refusal_names_both_ways_to_choose(unclaimed_host):
    with pytest.raises(RuntimeError) as excinfo:
        _autodetect_variant()
    message = str(excinfo.value)
    assert "bp-config" in message
    assert "BIOPIPELINES_CONFIG_VARIANT" in message


def test_refusal_never_names_a_variant_as_chosen(unclaimed_host):
    """The point of the change: no config is loaded on a guess."""
    with pytest.raises(RuntimeError):
        ConfigManager()


def test_ambiguous_match_still_raises(monkeypatch, unclaimed_host):
    monkeypatch.setattr(
        config_manager, "_scan_variants",
        lambda root: (["alpha", "beta"], ["alpha", "beta"], "someone"),
    )
    with pytest.raises(RuntimeError, match="Ambiguous"):
        _autodetect_variant()


# ── the explicit routes still work ───────────────────────────────────────────

def test_env_var_wins_on_an_otherwise_unclaimed_host(monkeypatch, unclaimed_host):
    monkeypatch.setenv("BIOPIPELINES_CONFIG_VARIANT", "daint")
    variant, source = _autodetect_variant_with_source()
    assert variant == "daint"
    assert "BIOPIPELINES_CONFIG_VARIANT" in source


def test_explicit_config_manager_variant_needs_no_detection(unclaimed_host, monkeypatch):
    monkeypatch.setattr(ConfigManager, "_variant", None)
    monkeypatch.setattr(ConfigManager, "_config", None)
    assert ConfigManager(variant="local").get_variant() == "local"


def test_explicit_pipeline_config_needs_no_detection(unclaimed_host, monkeypatch,
                                                     isolated_cwd):
    monkeypatch.setattr(ConfigManager, "_variant", None)
    monkeypatch.setattr(ConfigManager, "_config", None)
    from biopipelines.pipeline import Pipeline

    pipeline = Pipeline(project="TestSuite", job="variant_route",
                        description="explicit variant route",
                        on_the_fly=False, local_output=True, config="local")
    assert ConfigManager.get_variant() == "local"
    assert pipeline.folders["output"]


def test_colab_still_wins_silently(unclaimed_host, monkeypatch):
    """A real runtime signal, so it must not be caught up in the refusal."""
    colab = type(sys)("google.colab")
    google = type(sys)("google")
    google.colab = colab
    monkeypatch.setitem(sys.modules, "google", google)
    monkeypatch.setitem(sys.modules, "google.colab", colab)
    variant, source = _autodetect_variant_with_source()
    assert variant == "colab"
    assert "Colab" in source


def test_matched_username_is_returned_with_its_reason(monkeypatch, unclaimed_host):
    monkeypatch.setattr(
        config_manager, "_scan_variants",
        lambda root: (["cluster", "local"], ["local"], "someone"),
    )
    variant, source = _autodetect_variant_with_source()
    assert variant == "local"
    assert "machine.username" in source and "someone" in source


# ── announcement ─────────────────────────────────────────────────────────────

def test_resolution_announces_variant_file_and_reason(fresh_announcement, capsys):
    ConfigManager(variant="local")
    line = capsys.readouterr().err.strip()
    assert line.count("\n") == 0, f"expected exactly one line, got: {line!r}"
    assert "local" in line
    assert "config.local.yaml" in line
    assert "explicit" in line


def test_announcement_reports_the_detection_reason(fresh_announcement, capsys,
                                                   monkeypatch):
    monkeypatch.setenv("BIOPIPELINES_CONFIG_VARIANT", "local")
    ConfigManager()
    assert "BIOPIPELINES_CONFIG_VARIANT" in capsys.readouterr().err


def test_announcement_prints_once_per_process(fresh_announcement, capsys):
    ConfigManager(variant="local")
    first = capsys.readouterr().err
    ConfigManager()
    ConfigManager.get_variant()
    ConfigManager(variant="cluster")
    ConfigManager.get_variant()
    assert first.strip()
    assert capsys.readouterr().err == ""


def test_announcement_goes_to_stderr_not_stdout(fresh_announcement, capsys):
    """``$(bp-config path)`` and friends must stay machine-readable."""
    ConfigManager(variant="local")
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "config.local.yaml" in captured.err
