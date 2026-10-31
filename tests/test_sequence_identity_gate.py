# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""SequenceIdentityGate: the identity definition (matches over alignment length, gaps included), the nearest-reference and tie rules, NaN rows for sequences it cannot score, the exit codes, and the wrapper's validation and emitted bash."""

import glob
import importlib.util
import json
import os
import re
import sys

import pandas as pd
import pytest

pytest.importorskip("yaml")
pytest.importorskip("Bio")
import yaml  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "pipe_scripts", "pipe_sequence_identity_gate.py")


@pytest.fixture(scope="module")
def pipe():
    spec = importlib.util.spec_from_file_location("pipe_sequence_identity_gate", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ── identity ──────────────────────────────────────────────────────────────────

def test_an_identical_sequence_scores_one(pipe):
    assert pipe.global_identity(pipe.make_aligner(), "MKTAYIAKQR", "MKTAYIAKQR") == 1.0


def test_one_substitution_in_ten_scores_nine_tenths(pipe):
    assert pipe.global_identity(pipe.make_aligner(), "MKTAYIAKQR", "MKTWYIAKQR") == pytest.approx(0.9)


def test_gaps_count_in_the_denominator(pipe):
    """A design that is a prefix of a longer natural protein is not scored as that protein."""
    identity = pipe.global_identity(pipe.make_aligner(), "MKTAYIAKQR", "MKTAYIAKQRLLL")
    assert identity == pytest.approx(10 / 13)


def test_end_gaps_are_penalized_like_any_other_gap(pipe):
    """Needle's gap scores (open 10, extend 0.5) charge a 3-residue end gap 11; Needle itself would not charge it by default, and the docstring says so."""
    aligner = pipe.make_aligner()
    assert aligner.score("MKTAYIAKQR", "MKTAYIAKQR") - aligner.score("MKTAYIAKQR", "MKTAYIAKQRLLL") == pytest.approx(11.0)


def test_penalized_end_gaps_change_the_identity_of_a_design_cut_from_a_longer_protein(pipe):
    design = "QNIAFNHLCADFYHVELLCH"
    reference = "QNYQYVHIARQIRKFNHRCADIDSVEYNGHLPNLWGFYKT"
    free = pipe.make_aligner()
    free.end_gap_score = 0.0
    assert pipe.global_identity(pipe.make_aligner(), design, reference) == pytest.approx(0.325)
    assert pipe.global_identity(free, design, reference) == pytest.approx(0.25)


def test_the_docs_do_not_call_the_penalized_end_gaps_the_needle_defaults():
    from biopipelines import sequence_identity_gate

    for text in (sequence_identity_gate.SequenceIdentityGate.__doc__, open(SCRIPT, encoding="utf-8").read()):
        assert "the EMBOSS Needle defaults" not in text
        assert "end gaps" in text and "Needle" in text


def test_unrelated_sequences_score_low(pipe):
    assert pipe.global_identity(pipe.make_aligner(), "WWWWWWWWWW", "MKTAYIAKQR") < 0.2


def test_normalize_cleans_and_refuses(pipe):
    assert pipe.normalize(" mkt ayi*\n") == "MKTAYI"
    assert pipe.normalize("MK-TA") == "MKTA"
    for bad in ("", "   ", None, float("nan")):
        with pytest.raises(ValueError, match="empty"):
            pipe.normalize(bad)
    with pytest.raises(ValueError, match="BLOSUM62"):
        pipe.normalize("MKTUJ")


def test_the_nearest_reference_wins_and_the_first_wins_a_tie(pipe):
    aligner = pipe.make_aligner()
    references = [("far", "WWWWWWWWWW"), ("near", "MKTAYIAKQR"), ("twin", "MKTAYIAKQR")]
    identity, nearest = pipe.nearest_reference(aligner, "MKTAYIAKQR", references)
    assert (identity, nearest) == (1.0, "near")


# ── runtime ───────────────────────────────────────────────────────────────────

def _stream(tmp_path, name, sequences):
    """A value-based sequences stream: ids in the stream, sequences in its map table."""
    table = tmp_path / f"{name}.csv"
    pd.DataFrame({"id": list(sequences), "sequence": list(sequences.values())}).to_csv(table, index=False)
    path = tmp_path / f"{name}.json"
    path.write_text(json.dumps({"name": "sequences", "ids": list(sequences), "files": [],
                                "map_table": str(table), "format": "csv"}))
    return str(path)


def _run_main(pipe, monkeypatch, tmp_path, designs, references, threshold=0.4):
    config = tmp_path / "config.yaml"
    config.write_text(yaml.safe_dump({"threshold": threshold}))
    out = tmp_path / "tables" / "identity_gate.csv"
    local = tmp_path / "_execution" / "local_missing.csv"
    monkeypatch.setattr(sys, "argv", [
        "pipe_sequence_identity_gate.py", "--config-yaml", str(config),
        "--sequences-json", _stream(tmp_path, "designs", designs),
        "--references-json", _stream(tmp_path, "references", references),
        "--identity-gate-csv", str(out), "--local-missing-csv", str(local)])
    pipe.main()
    return pd.read_csv(out), pd.read_csv(local)


NATURALS = {"nat1": "MKTAYIAKQRQISFVKSHFSRQLEERLGL", "nat2": "GSHMLEDPVDAFKQLRELL"}


def test_the_table_reports_identity_nearest_reference_and_verdict(pipe, monkeypatch, tmp_path):
    designs = {"copy": NATURALS["nat1"], "novel": "WWWWWWWWWWWWWWWWWWWW"}
    table, missing = _run_main(pipe, monkeypatch, tmp_path, designs, NATURALS)
    assert list(table.columns) == ["id", "max_identity", "nearest_reference", "passed"]
    assert table["id"].tolist() == ["copy", "novel"]
    assert table.loc[0, "max_identity"] == 1.0 and table.loc[0, "nearest_reference"] == "nat1"
    assert table.loc[1, "max_identity"] < 0.4
    assert table["passed"].tolist() == [False, True]
    assert missing.empty


def test_the_threshold_is_a_strict_upper_bound(pipe, monkeypatch, tmp_path):
    designs = {"d": "MKTWYIAKQR"}
    refs = {"r": "MKTAYIAKQR"}
    table, _ = _run_main(pipe, monkeypatch, tmp_path, designs, refs, threshold=0.9)
    assert table.loc[0, "max_identity"] == pytest.approx(0.9)
    assert table.loc[0, "passed"] == False  # noqa: E712 - equal to the threshold does not pass


def test_a_sequence_that_cannot_be_scored_keeps_a_nan_row_and_a_failure(pipe, monkeypatch, tmp_path):
    designs = {"ok": NATURALS["nat2"], "bad": "MKTUJ"}
    table, missing = _run_main(pipe, monkeypatch, tmp_path, designs, NATURALS)
    assert table["id"].tolist() == ["ok", "bad"]
    assert table.loc[1, ["max_identity", "nearest_reference", "passed"]].isna().all()
    assert missing["id"].tolist() == ["bad"] and missing["kind"].tolist() == ["failure"]
    assert "BLOSUM62" in missing.loc[0, "cause"]


def test_every_design_failing_is_a_step_failure(pipe, monkeypatch, tmp_path):
    with pytest.raises(SystemExit) as raised:
        _run_main(pipe, monkeypatch, tmp_path, {"a": "MKTUJ", "b": "JJJ"}, NATURALS)
    assert raised.value.code == 1


def test_a_reference_that_cannot_be_scored_stops_the_step(pipe, monkeypatch, tmp_path, capsys):
    with pytest.raises(SystemExit) as raised:
        _run_main(pipe, monkeypatch, tmp_path, {"a": NATURALS["nat1"]}, {"sec": "MKTAUYIAKQR"})
    assert raised.value.code == 1
    assert "reference 'sec'" in capsys.readouterr().err
    assert not (tmp_path / "tables" / "identity_gate.csv").exists()


def test_a_repeated_reference_sequence_keeps_its_first_id(pipe, tmp_path):
    refs = pipe.load_references(_stream(tmp_path, "refs", {"first": "MKTAYIAKQR", "second": "MKTAYIAKQR",
                                                           "other": "GSHMLEDPVD"}))
    assert refs == [("first", "MKTAYIAKQR"), ("other", "GSHMLEDPVD")]


# ── wrapper ───────────────────────────────────────────────────────────────────

def _sequences(seqs):
    from biopipelines.sequence import Sequence
    return Sequence(seq=seqs, type="protein")


def test_the_gate_table_declares_its_columns(local_config, isolated_cwd, new_pipeline):
    from biopipelines.sequence_identity_gate import SequenceIdentityGate

    with new_pipeline("gate_columns"):
        out = SequenceIdentityGate(sequences=_sequences({"d1": "MKTAYIAKQR"}),
                                   references=_sequences({"r1": "GSHMLEDPVD"}))
    assert out.tables.identity_gate.info.columns == ["id", "max_identity", "nearest_reference", "passed"]
    assert out.tables.missing.info.columns == ["id", "removed_by", "kind", "cause"]


@pytest.mark.parametrize("kwargs, message", [
    ({"threshold": -0.1}, "threshold"),
    ({"threshold": 40}, "threshold"),
    ({"references": "naturals.fasta"}, "references must be"),
])
def test_bad_parameters_are_refused_at_configuration(local_config, isolated_cwd, new_pipeline, kwargs, message):
    from biopipelines.sequence_identity_gate import SequenceIdentityGate

    with new_pipeline("gate_bad"):
        arguments = {"sequences": _sequences({"d1": "MKTAYIAKQR"}), "references": _sequences({"r1": "GSHMLEDPVD"})}
        with pytest.raises(ValueError, match=message):
            SequenceIdentityGate(**dict(arguments, **kwargs))


def test_the_emitted_script_runs_the_pipe_with_every_declared_flag(local_config, isolated_cwd, new_pipeline):
    from biopipelines.sequence_identity_gate import SequenceIdentityGate

    pipeline = new_pipeline("gate_script")
    with pipeline:
        SequenceIdentityGate(sequences=_sequences({"d1": "MKTAYIAKQR"}),
                             references=_sequences({"r1": "GSHMLEDPVD"}), threshold=0.35)
        script_path = pipeline.save()
    step = sorted(glob.glob(os.path.join(os.path.dirname(script_path), "*_SequenceIdentityGate.sh")))[0]
    text = open(step, encoding="utf-8", errors="replace").read()
    call = text[text.index("pipe_sequence_identity_gate.py"):]
    call = call[:call.index("python ", 1)] if "python " in call[1:] else call
    emitted = set(re.findall(r"--[a-z][a-z-]+", call))
    declared = set(re.findall(r'add_argument\("(--[a-z-]+)"', open(SCRIPT, encoding="utf-8").read()))
    assert emitted == declared
    assert yaml.safe_load(open(pipeline.tools[-1].config_yaml))["threshold"] == 0.35
