import importlib.util
import logging
import sys
import types
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def load_pipe_script(name):
    sys.modules.setdefault("logomaker", types.ModuleType("logomaker"))
    spec = importlib.util.spec_from_file_location(name, ROOT / "pipe_scripts" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def one_position_frequency_row(original="A", mutant="C", mutant_prior=0.2):
    row = {"position": 1, "original": original}
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        row[aa] = 0.0
    row[original] = 1.0 - mutant_prior
    row[mutant] = mutant_prior
    return row


def one_position_correlation_row(original="A", mutant="C", correlation=1.0):
    row = {"position": 1, "original": original}
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        row[aa] = 0.0
    row[mutant] = correlation
    return row


def one_position_sample_count_row(original="A", mutant="C", count=2):
    row = {"position": 1, "original": original}
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        row[aa] = 0
    row[mutant] = count
    return row


def test_kappa_shrinks_bayesian_adjustment_when_counts_are_available():
    mod = load_pipe_script("pipe_bayesian_adjuster")
    logger = logging.getLogger("test")
    frequencies = pd.DataFrame([one_position_frequency_row()])
    correlations = pd.DataFrame([one_position_correlation_row()])
    sample_counts = pd.DataFrame([one_position_sample_count_row(count=2)])

    no_shrink_df, no_shrink_log = mod.apply_bayesian_adjustment(
        frequencies, correlations, "max", gamma=2.0, kappa=None, sample_counts_df=None, pseudocount=0.0, logger=logger
    )
    shrunk_df, shrunk_log = mod.apply_bayesian_adjustment(
        frequencies, correlations, "max", gamma=2.0, kappa=8.0, sample_counts_df=sample_counts, pseudocount=0.0, logger=logger
    )

    assert shrunk_df.loc[0, "C"] < no_shrink_df.loc[0, "C"]

    shrunk_c_log = next(row for row in shrunk_log if row["aa"] == "C")
    no_shrink_c_log = next(row for row in no_shrink_log if row["aa"] == "C")
    assert shrunk_c_log["raw_correlation"] == 1.0
    assert shrunk_c_log["correlation"] == 0.2
    assert shrunk_c_log["sample_count"] == 2
    assert shrunk_c_log["shrinkage_weight"] == 0.2
    assert no_shrink_c_log["correlation"] == 1.0


def test_kappa_requires_sample_counts_when_provided():
    mod = load_pipe_script("pipe_bayesian_adjuster")
    logger = logging.getLogger("test")
    frequencies = pd.DataFrame([one_position_frequency_row()])
    correlations = pd.DataFrame([one_position_correlation_row()])

    try:
        mod.apply_bayesian_adjustment(
            frequencies, correlations, "max", gamma=2.0, kappa=8.0,
            sample_counts_df=None, pseudocount=0.0, logger=logger
        )
    except ValueError as exc:
        assert "sample_counts_path is required when kappa is provided" in str(exc)
    else:
        raise AssertionError("Expected kappa > 0 without sample counts to raise")


def test_sequence_metric_correlation_2d_emits_separate_amino_acid_sample_counts():
    mod = load_pipe_script("pipe_sequence_metric_correlation")
    logger = logging.getLogger("test")
    merged = pd.DataFrame(
        [
            {"id": "s1", "sequence": "A", "score": 0.0},
            {"id": "s2", "sequence": "C", "score": 1.0},
            {"id": "s3", "sequence": "C", "score": 2.0},
            {"id": "s4", "sequence": "D", "score": 3.0},
        ]
    )

    corr_2d, sample_counts = mod.compute_correlation_2d(merged, "A", "score", logger)

    assert "n_A" not in corr_2d.columns
    assert corr_2d.loc[0, "original"] == "A"
    assert sample_counts.loc[0, "original"] == "A"
    assert sample_counts.loc[0, "A"] == 1
    assert sample_counts.loc[0, "C"] == 2
    assert sample_counts.loc[0, "D"] == 1


def test_bayesian_adjuster_accepts_legacy_wt_aa_correlation_schema():
    mod = load_pipe_script("pipe_bayesian_adjuster")
    logger = logging.getLogger("test")

    frequencies_path = ROOT / "tests" / "_tmp_frequencies.csv"
    correlations_path = ROOT / "tests" / "_tmp_correlations.csv"
    try:
        pd.DataFrame([one_position_frequency_row()]).to_csv(frequencies_path, index=False)
        legacy_corr = one_position_correlation_row()
        legacy_corr["wt_aa"] = legacy_corr.pop("original")
        pd.DataFrame([legacy_corr]).to_csv(correlations_path, index=False)

        _, correlations = mod.load_and_merge_tables(str(frequencies_path), str(correlations_path), logger)

        assert "original" in correlations.columns
        assert "wt_aa" not in correlations.columns
        assert correlations.loc[0, "original"] == "A"
    finally:
        frequencies_path.unlink(missing_ok=True)
        correlations_path.unlink(missing_ok=True)
