# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""USalign and TMalign: the shared pair handling (pair table, id checks, parsing, the NaN rules), the pipe with the binary replaced by recorded output, and both wrappers.

The fixtures in tests/fixtures/pair_alignment are the real compact output (`-outfmt 1 -m -`) of US-align 20241108 and TM-align 20240303 as installed by bioconda, run on 101M against 1MBA (a sequential pair) and on 1AJK against 2AYH (a circular permutant). Nothing here runs the binaries.
"""

import glob
import importlib.util
import json
import math
import os
import re
import subprocess
import sys

import pandas as pd
import pytest

pytest.importorskip("yaml")
import yaml  # noqa: E402

from biopipelines import pair_alignment as pa  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "pipe_scripts", "pipe_pair_alignment.py")
FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures", "pair_alignment")

MODES = [{"label": pa.MODE_SEQUENTIAL, "args": ["-mm", "0"]},
         {"label": pa.MODE_CIRCULAR_PERMUTATION, "args": ["-mm", "3"]},
         {"label": pa.MODE_NON_SEQUENTIAL, "args": ["-mm", "5"]}]


def fixture_text(name):
    with open(os.path.join(FIXTURES, name), encoding="utf-8") as f:
        return f.read()


@pytest.fixture(scope="module")
def pipe():
    spec = importlib.util.spec_from_file_location("pipe_pair_alignment", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ── parsing real output ───────────────────────────────────────────────────────

def test_the_sequential_output_is_read_in_full():
    v = pa.parse_alignment(fixture_text("usalign_sequential.txt"))
    assert (v["length_1"], v["length_2"], v["aligned_length"], v["rmsd"]) == (154, 146, 143, 1.83)
    assert (v["tm_score_1"], v["tm_score_2"]) == (0.81453, 0.85377)
    assert v["aligned_1"].startswith("MVLSEGEWQLVLHVWAKVEADVAGHGQD") and v["aligned_2"].startswith("-SLSAAEADLAGKSW")
    assert v["translation_1"] == pytest.approx(-71.3397346333)
    assert [v["rotation_11"], v["rotation_12"], v["rotation_13"]] == pytest.approx(
        [0.0562537143, 0.6667196713, 0.7431826152])
    assert [v["rotation_21"], v["rotation_22"], v["rotation_23"]] == pytest.approx(
        [0.7560213042, -0.5146315959, 0.4044577952])
    assert [v["rotation_31"], v["rotation_32"], v["rotation_33"]] == pytest.approx(
        [0.6521252236, 0.5391096368, -0.5330042142])
    assert (v["translation_2"], v["translation_3"]) == pytest.approx((-66.1421187805, -53.2593104458))


def test_the_tm_score_of_each_structure_is_normalized_by_its_own_length():
    """US-align prints the structure_1 record first, with its length and the TM-score normalized by it."""
    v = pa.parse_alignment(fixture_text("usalign_sequential.txt"))
    assert v["length_1"] > v["length_2"] and v["tm_score_1"] < v["tm_score_2"]


@pytest.mark.parametrize("name", sorted(os.path.basename(p) for p in
                                        glob.glob(os.path.join(FIXTURES, "*_*.txt"))
                                        if "unparseable" not in p))
def test_every_recorded_output_is_internally_consistent(name):
    """The aligned strings agree with the reported lengths and aligned length, in every mode and for both binaries."""
    v = pa.parse_alignment(fixture_text(name))
    a, b = v["aligned_1"], v["aligned_2"]
    assert len(a) == len(b)
    assert len(a.replace("-", "").replace("*", "")) == v["length_1"]
    assert len(b.replace("-", "").replace("*", "")) == v["length_2"]
    assert sum(1 for x, y in zip(a, b) if x not in "-*" and y not in "-*") == v["aligned_length"]
    assert set(pa.ALIGNMENT_COLUMNS) - {"id", "structure_1", "structure_2", "mode"} == set(v)


def test_a_circular_permutation_is_marked_in_the_first_sequence_only():
    v = pa.parse_alignment(fixture_text("usalign_circular_permutation.txt"))
    assert v["aligned_1"].count("*") == 1 and "*" not in v["aligned_2"]
    assert (v["aligned_length"], v["tm_score_1"], v["tm_score_2"]) == (207, 0.94331, 0.93466)


def test_no_circular_permutation_returns_the_sequential_alignment():
    seq = pa.parse_alignment(fixture_text("usalign_sequential.txt"))
    none = pa.parse_alignment(fixture_text("usalign_no_circular_permutation.txt"))
    assert none == seq


def test_the_non_sequential_variants_are_read():
    fns = pa.parse_alignment(fixture_text("usalign_fns.txt"))
    sns = pa.parse_alignment(fixture_text("usalign_sns.txt"))
    assert (fns["aligned_length"], fns["tm_score_2"]) == (209, 0.94465)
    assert (sns["aligned_length"], sns["tm_score_2"]) == (209, 0.94389)
    assert "*" not in fns["aligned_1"]


def test_tmalign_and_usalign_agree_on_the_modes_they_share():
    for tm_name, us_name in (("tmalign_sequential.txt", "usalign_sequential.txt"),
                             ("tmalign_circular_permutation.txt", "usalign_circular_permutation.txt")):
        assert pa.parse_alignment(fixture_text(tm_name)) == pa.parse_alignment(fixture_text(us_name))


def test_windows_line_endings_are_read():
    text = fixture_text("usalign_circular_permutation.txt").replace("\n", "\r\n")
    assert pa.parse_alignment(text) == pa.parse_alignment(fixture_text("usalign_circular_permutation.txt"))


def test_output_without_an_alignment_is_refused():
    with pytest.raises(ValueError, match="no alignment"):
        pa.parse_alignment(fixture_text("usalign_unparseable.txt"))
    with pytest.raises(ValueError, match="no alignment"):
        pa.parse_alignment("")


@pytest.mark.parametrize("cut, message", [
    (lambda t: t.split("# Lali")[0], "Lali"),
    (lambda t: t.split("------ The rotation")[0], "rotation matrix"),
    (lambda t: re.sub(r"\n2 +-?\d+\.\d+ +-?\d+\.\d+ +-?\d+\.\d+ +-?\d+\.\d+", "", t), "rotation matrix"),
    (lambda t: t.replace("\n>1MBA", "\n#>1MBA"), "found 1 records"),
])
def test_a_truncated_output_is_refused_rather_than_read_as_zero(cut, message):
    with pytest.raises(ValueError, match=message):
        pa.parse_alignment(cut(fixture_text("usalign_sequential.txt")))


def test_two_aligned_sequences_of_different_length_are_refused():
    broken = fixture_text("usalign_sequential.txt").replace("-SLSAAEADLAGKSW", "SLSAAEADLAGKSW", 1)
    with pytest.raises(ValueError, match="differ in length"):
        pa.parse_alignment(broken)


# ── running the aligner ───────────────────────────────────────────────────────

def test_the_command_carries_the_files_the_mode_and_the_output_flags():
    cmd = pa.build_command("USalign", "a.pdb", "b.cif", ["-mm", "3"])
    assert cmd == ["USalign", "a.pdb", "b.cif", "-mm", "3", "-outfmt", "1", "-m", "-"]


def test_the_container_prefix_leads_the_command():
    cmd = pa.build_command("TMalign", "a.pdb", "b.pdb", ["-cp"], "apptainer exec -B /data img.sif ")
    assert cmd[:5] == ["apptainer", "exec", "-B", "/data", "img.sif"] and cmd[5] == "TMalign"


def test_an_unreadable_file_gives_the_aligners_own_words():
    """The binary exits 0 on a file it cannot parse and prints only a warning."""
    def run(cmd):
        return 0, fixture_text("usalign_unparseable.txt"), fixture_text("usalign_unparseable.stderr")
    with pytest.raises(RuntimeError, match=r"no alignment.*Cannot parse file: bad.pdb"):
        pa.align_one(["USalign"], run)


def test_a_crash_is_a_failure_even_when_it_wrote_output():
    def crashed(cmd):
        return -11, fixture_text("usalign_sequential.txt"), "Warning! Cannot parse file: bad.pdb. Chain number 0."
    with pytest.raises(RuntimeError, match=r"killed by signal 11: Warning!"):
        pa.align_one(["TMalign"], crashed)

    def failed(cmd):
        return 2, "", ""
    with pytest.raises(RuntimeError, match=r"^aligner exited 2$"):
        pa.align_one(["TMalign"], failed)


def test_run_command_returns_the_exit_code_and_both_streams(monkeypatch):
    def fake(cmd, **kwargs):
        assert kwargs["capture_output"] and kwargs["text"]
        return subprocess.CompletedProcess(cmd, 3, "out", "err")
    monkeypatch.setattr(pa.subprocess, "run", fake)
    assert pa.run_command(["x"]) == (3, "out", "err")


# ── pair table ────────────────────────────────────────────────────────────────

def ref(path, column):
    return f"TABLE_REFERENCE:{path}:{column}"


def write_pairs(tmp_path, rows, columns=("parent", "variant")):
    path = tmp_path / "pairs.csv"
    pd.DataFrame(rows, columns=list(columns)).to_csv(path, index=False)
    return str(path)


def test_pairs_are_distinct_and_in_first_seen_order(tmp_path):
    path = write_pairs(tmp_path, [("a", "b"), ("c", "d"), ("a", "b"), ("b", "a")])
    assert pa.read_pairs(ref(path, "parent"), ref(path, "variant")) == [("a", "b"), ("c", "d"), ("b", "a")]


def test_ids_are_read_as_text_and_trimmed(tmp_path):
    path = write_pairs(tmp_path, [("007", " NA "), ("1e5", "nan")])
    assert pa.read_pairs(ref(path, "parent"), ref(path, "variant")) == [("007", "NA"), ("1e5", "nan")]


def test_two_pairs_that_would_share_an_id_are_refused(tmp_path):
    path = write_pairs(tmp_path, [("x+y", "z"), ("x", "y+z")])
    with pytest.raises(ValueError, match=r"both have the id 'x\+y\+z'"):
        pa.read_pairs(ref(path, "parent"), ref(path, "variant"))


def test_a_column_the_table_lacks_is_named(tmp_path):
    path = write_pairs(tmp_path, [("a", "b")])
    with pytest.raises(ValueError, match=r"'child' is not in the pair table.*parent"):
        pa.read_pairs(ref(path, "parent"), ref(path, "child"))


def test_columns_of_two_tables_are_refused(tmp_path):
    first = write_pairs(tmp_path, [("a", "b")])
    other = tmp_path / "other.csv"
    other.write_text("parent,variant\na,b\n")
    with pytest.raises(ValueError, match="one table"):
        pa.read_pairs(ref(first, "parent"), ref(str(other), "variant"))


def test_a_missing_pair_table_is_an_error(tmp_path):
    with pytest.raises(FileNotFoundError, match="Pair table not found"):
        pa.read_pairs(ref(str(tmp_path / "none.csv"), "a"), ref(str(tmp_path / "none.csv"), "b"))


def test_an_id_absent_from_the_stream_is_told_apart_from_one_without_a_file():
    files, declared = {"a": "a.pdb"}, {"a", "b"}
    assert pa.resolve_structures("a", "a", files, declared) == ("a.pdb", "a.pdb")
    with pytest.raises(pa.PairError, match=r"structure_2 'z' is not in the structures stream"):
        pa.resolve_structures("a", "z", files, declared)
    with pytest.raises(pa.PairError, match=r"structure_1 'b' has no file .*removed upstream"):
        pa.resolve_structures("b", "a", files, declared)
    with pytest.raises(pa.PairError, match="structure_2 is empty"):
        pa.resolve_structures("a", "", files, declared)


# ── rows, NaN and failures ────────────────────────────────────────────────────

REPLAY = {"-mm 0": "usalign_sequential.txt", "-mm 3": "usalign_circular_permutation.txt",
          "-mm 5": "usalign_fns.txt", "-cp": "tmalign_circular_permutation.txt", "": "tmalign_sequential.txt"}


def replay(cmd):
    """Recorded output chosen by the mode flags of the command."""
    flags = " ".join(cmd[cmd.index("-outfmt") - 2:cmd.index("-outfmt")]) if "-mm" in cmd else (
        "-cp" if "-cp" in cmd else "")
    return 0, fixture_text(REPLAY[flags]), ""


def test_every_pair_and_mode_gives_one_row():
    pairs = [("a", "b"), ("c", "d")]
    files = {s: f"{s}.pdb" for s in "abcd"}
    rows, missing = pa.align_pairs(pairs, files, set(files), MODES, "USalign", "007_USalign",
                                   run=replay, log=lambda m: None)
    assert [r["id"] for r in rows] == [
        "a+b_sequential", "a+b_circular_permutation", "a+b_non_sequential",
        "c+d_sequential", "c+d_circular_permutation", "c+d_non_sequential"]
    assert missing == []
    assert [r["mode"] for r in rows[:3]] == ["sequential", "circular_permutation", "non_sequential"]
    assert rows[1]["tm_score_2"] == 0.93466 and rows[2]["aligned_length"] == 209
    assert (rows[0]["structure_1"], rows[0]["structure_2"]) == ("a", "b")


def test_a_pair_with_an_id_outside_the_stream_keeps_nan_rows_and_a_failure_per_row():
    files = {"a": "a.pdb", "b": "b.pdb"}
    rows, missing = pa.align_pairs([("a", "b"), ("a", "ghost")], files, set(files), MODES, "USalign",
                                   "007_USalign", run=replay, log=lambda m: None)
    assert len(rows) == 6
    ghost = [r for r in rows if r["structure_2"] == "ghost"]
    assert len(ghost) == 3
    for row in ghost:
        assert all(math.isnan(row[c]) for c in pa.NUMERIC_COLUMNS)
        assert row["aligned_1"] == "" and row["aligned_2"] == ""
        assert row["id"].startswith("a+ghost_")
    assert [m["id"] for m in missing] == [r["id"] for r in ghost]
    assert {m["kind"] for m in missing} == {"failure"} and {m["removed_by"] for m in missing} == {"007_USalign"}
    assert "'ghost' is not in the structures stream" in missing[0]["cause"]


def test_one_mode_failing_leaves_the_others_intact():
    def run(cmd):
        if "3" in cmd[cmd.index("-mm") + 1:cmd.index("-mm") + 2]:
            return 0, "", "Sequence is too short <3!: a.pdb"
        return replay(cmd)
    files = {"a": "a.pdb", "b": "b.pdb"}
    rows, missing = pa.align_pairs([("a", "b")], files, set(files), MODES, "USalign", "s",
                                   run=run, log=lambda m: None)
    assert [r["mode"] for r in rows] == ["sequential", "circular_permutation", "non_sequential"]
    assert math.isnan(rows[1]["tm_score_1"]) and not math.isnan(rows[0]["tm_score_1"])
    assert [m["id"] for m in missing] == ["a+b_circular_permutation"]
    assert "too short" in missing[0]["cause"]


def test_a_cause_is_cut_to_200_characters():
    files = {"a": "a.pdb", "b": "b.pdb"}
    rows, missing = pa.align_pairs([("a", "b")], files, set(files), MODES[:1], "USalign", "s",
                                   run=lambda cmd: (0, "", "x" * 500), log=lambda m: None)
    assert len(missing[0]["cause"]) == 200


# ── the TM-score difference ───────────────────────────────────────────────────

def tm_rows(pairs):
    """Alignment rows with chosen TM-scores: {(s1, s2): {mode: (tm_1, tm_2)}}."""
    rows = []
    for (s1, s2), modes in pairs.items():
        for mode, (tm_1, tm_2) in modes.items():
            row = pa.empty_row(s1, s2, mode)
            row.update(tm_score_1=tm_1, tm_score_2=tm_2)
            rows.append(row)
    return rows


def test_the_delta_is_permuted_minus_sequential_for_each_normalization():
    rows = tm_rows({("a", "b"): {"sequential": (0.30, 0.40), "circular_permutation": (0.90, 0.85)}})
    [delta] = pa.delta_rows(rows)
    assert delta["id"] == "a+b" and (delta["structure_1"], delta["structure_2"]) == ("a", "b")
    assert delta["delta_tm_1"] == pytest.approx(0.60) and delta["delta_tm_2"] == pytest.approx(0.45)


def test_the_delta_of_a_pair_without_a_value_is_nan_and_the_row_stays():
    rows = tm_rows({("a", "b"): {"sequential": (0.3, 0.4), "circular_permutation": (math.nan, math.nan)},
                    ("c", "d"): {"sequential": (0.5, 0.6), "circular_permutation": (0.5, 0.6)}})
    deltas = pa.delta_rows(rows)
    assert [d["id"] for d in deltas] == ["a+b", "c+d"]
    assert math.isnan(deltas[0]["delta_tm_1"]) and math.isnan(deltas[0]["delta_tm_2"])
    assert deltas[1]["delta_tm_1"] == 0.0


def test_the_delta_is_nan_when_a_mode_was_not_run():
    [delta] = pa.delta_rows(tm_rows({("a", "b"): {"sequential": (0.3, 0.4)}}))
    assert math.isnan(delta["delta_tm_1"])


def test_the_real_circular_permutant_gains_and_the_sequential_pair_gains_nothing():
    seq = pa.parse_alignment(fixture_text("tmalign_sequential.txt"))
    cp = pa.parse_alignment(fixture_text("tmalign_circular_permutation.txt"))
    none = pa.parse_alignment(fixture_text("tmalign_no_circular_permutation.txt"))
    assert cp["tm_score_2"] - seq["tm_score_2"] > 0.0
    assert none["tm_score_1"] - seq["tm_score_1"] == 0.0


# ── the pipe ──────────────────────────────────────────────────────────────────

def structures_json(tmp_path, ids, with_files=None):
    files = []
    for sid in ids:
        path = tmp_path / f"{sid}.pdb"
        if with_files is None or sid in with_files:
            path.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n")
        files.append(str(path))
    out = tmp_path / "structures.json"
    out.write_text(json.dumps({"name": "structures", "ids": list(ids), "files": files, "map_table": "",
                               "format": "pdb"}))
    return str(out)


def fake_subprocess(monkeypatch, outputs):
    """Replace subprocess.run in pair_alignment; `outputs(cmd)` returns (code, stdout, stderr)."""
    def fake(cmd, **kwargs):
        code, out, err = outputs(cmd)
        return subprocess.CompletedProcess(cmd, code, out, err)
    monkeypatch.setattr(pa.subprocess, "run", fake)


def run_pipe(pipe, monkeypatch, tmp_path, pairs, ids, binary="USalign", modes=MODES, delta=False, with_files=None):
    table = write_pairs(tmp_path, pairs)
    config = tmp_path / "config.yaml"
    config.write_text(yaml.safe_dump({"binary": binary, "modes": modes,
                                      "structure_1": ref(table, "parent"), "structure_2": ref(table, "variant")}))
    alignments = tmp_path / "tables" / "alignments.csv"
    local = tmp_path / "_execution" / "local_missing.csv"
    deltas = tmp_path / "tables" / "delta.csv"
    argv = ["pipe_pair_alignment.py", "--config-yaml", str(config),
            "--structures-json", structures_json(tmp_path, ids, with_files),
            "--alignments-csv", str(alignments), "--local-missing-csv", str(local)]
    if delta:
        argv += ["--delta-csv", str(deltas)]
    monkeypatch.setattr(sys, "argv", argv)
    pipe.main()
    out = {"alignments": pd.read_csv(alignments, keep_default_na=True), "missing": pd.read_csv(local)}
    if delta:
        out["delta"] = pd.read_csv(deltas)
    return out


def test_the_pipe_writes_one_row_per_pair_and_mode(pipe, monkeypatch, tmp_path):
    fake_subprocess(monkeypatch, replay)
    out = run_pipe(pipe, monkeypatch, tmp_path, [("a", "b"), ("b", "c")], ["a", "b", "c"])
    table = out["alignments"]
    assert list(table.columns) == pa.ALIGNMENT_COLUMNS
    assert table["id"].tolist() == [
        "a+b_sequential", "a+b_circular_permutation", "a+b_non_sequential",
        "b+c_sequential", "b+c_circular_permutation", "b+c_non_sequential"]
    assert table["tm_score_2"].tolist()[:3] == [0.85377, 0.93466, 0.94465]
    assert table["aligned_1"].str.contains(r"\*").tolist()[:3] == [False, True, False]
    assert out["missing"].empty and list(out["missing"].columns) == pa.MISSING_COLUMNS


def test_the_pipe_keeps_the_rows_of_a_pair_with_an_unknown_id_as_nan(pipe, monkeypatch, tmp_path):
    fake_subprocess(monkeypatch, replay)
    out = run_pipe(pipe, monkeypatch, tmp_path, [("a", "b"), ("a", "ghost")], ["a", "b"])
    table, missing = out["alignments"], out["missing"]
    assert len(table) == 6
    ghost = table[table["structure_2"] == "ghost"]
    assert len(ghost) == 3 and ghost["tm_score_1"].isna().all() and ghost["aligned_1"].isna().all()
    assert table[table["structure_2"] == "b"]["tm_score_1"].notna().all()
    assert sorted(missing["id"]) == sorted(ghost["id"]) and set(missing["kind"]) == {"failure"}


def test_length_columns_stay_whole_numbers_beside_a_failed_row(pipe, monkeypatch, tmp_path):
    fake_subprocess(monkeypatch, replay)
    run_pipe(pipe, monkeypatch, tmp_path, [("a", "b"), ("a", "ghost")], ["a", "b"])
    text = (tmp_path / "tables" / "alignments.csv").read_text()
    first = text.splitlines()[1].split(",")
    assert first[4:7] == ["154", "146", "143"]
    assert ".0," not in ",".join(first[4:7]) + ","


def test_a_structure_removed_upstream_is_named_as_such(pipe, monkeypatch, tmp_path):
    fake_subprocess(monkeypatch, replay)
    out = run_pipe(pipe, monkeypatch, tmp_path, [("a", "b"), ("c", "b")], ["a", "b", "c"], with_files={"a", "b"})
    assert out["missing"]["cause"].str.contains("removed upstream").all()
    assert sorted(out["missing"]["id"]) == sorted(i for i in out["alignments"]["id"] if i.startswith("c+b"))


def test_the_pipe_writes_the_delta_table_when_asked(pipe, monkeypatch, tmp_path):
    fake_subprocess(monkeypatch, replay)
    modes = [{"label": "sequential", "args": []}, {"label": "circular_permutation", "args": ["-cp"]}]
    out = run_pipe(pipe, monkeypatch, tmp_path, [("a", "b")], ["a", "b"], binary="TMalign", modes=modes, delta=True)
    assert list(out["alignments"].columns) == pa.ALIGNMENT_COLUMNS
    assert list(out["delta"].columns) == pa.DELTA_COLUMNS
    [delta] = out["delta"].to_dict("records")
    assert delta["id"] == "a+b"
    assert delta["delta_tm_1"] == pytest.approx(0.94331 - 0.81453) and delta["delta_tm_2"] == pytest.approx(0.93466 - 0.85377)


def test_no_delta_file_is_written_unless_asked(pipe, monkeypatch, tmp_path):
    fake_subprocess(monkeypatch, replay)
    run_pipe(pipe, monkeypatch, tmp_path, [("a", "b")], ["a", "b"])
    assert not (tmp_path / "tables" / "delta.csv").exists()


def test_every_alignment_failing_fails_the_step_but_keeps_the_table(pipe, monkeypatch, tmp_path):
    fake_subprocess(monkeypatch, lambda cmd: (0, "", "Warning! Cannot parse file: a.pdb. Chain number 0."))
    with pytest.raises(SystemExit) as raised:
        run_pipe(pipe, monkeypatch, tmp_path, [("a", "b")], ["a", "b"])
    assert raised.value.code == 1
    table = pd.read_csv(tmp_path / "tables" / "alignments.csv")
    assert len(table) == 3 and table["tm_score_1"].isna().all()
    assert len(pd.read_csv(tmp_path / "_execution" / "local_missing.csv")) == 3


def test_an_empty_pair_table_fails_the_step(pipe, monkeypatch, tmp_path):
    fake_subprocess(monkeypatch, replay)
    with pytest.raises(SystemExit) as raised:
        run_pipe(pipe, monkeypatch, tmp_path, [], ["a", "b"])
    assert raised.value.code == 1


def test_a_failed_alignment_among_good_ones_does_not_fail_the_step(pipe, monkeypatch, tmp_path):
    def outputs(cmd):
        if cmd[1].endswith("b.pdb"):
            return 0, "", "Sequence is too short <3!: b.pdb"
        return replay(cmd)
    fake_subprocess(monkeypatch, outputs)
    out = run_pipe(pipe, monkeypatch, tmp_path, [("a", "c"), ("b", "c")], ["a", "b", "c"])
    assert len(out["alignments"]) == 6 and len(out["missing"]) == 3


# ── wrappers ──────────────────────────────────────────────────────────────────

def _inputs():
    from biopipelines.mock import Mock
    structures = Mock(ids=["s1", "s2"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
    pairs = Mock(ids=["p"], tables={"pairs": {"columns": ["parent", "variant"]}})
    return structures, pairs


def test_both_wrappers_declare_the_columns_the_pipe_writes(local_config, isolated_cwd, new_pipeline):
    from biopipelines import TMalign, USalign

    with new_pipeline("pa_columns"):
        structures, pairs = _inputs()
        us = USalign(structures=structures, structure_1=pairs.tables.pairs.variant,
                     structure_2=pairs.tables.pairs.parent)
        tm = TMalign(structures=structures, structure_1=pairs.tables.pairs.variant,
                     structure_2=pairs.tables.pairs.parent)
    assert us.tables.alignments.info.columns == pa.ALIGNMENT_COLUMNS
    assert tm.tables.alignments.info.columns == pa.ALIGNMENT_COLUMNS
    assert tm.tables.delta.info.columns == pa.DELTA_COLUMNS
    assert "delta" not in us.tables._tables
    for out in (us, tm):
        assert out.tables.missing.info.columns == pa.MISSING_COLUMNS
        assert list(out.streams.items()) == []


def test_the_wrappers_declare_each_mode_the_aligner_runs(local_config, isolated_cwd, new_pipeline):
    from biopipelines import TMalign, USalign

    pipeline = new_pipeline("pa_modes")
    with pipeline:
        structures, pairs = _inputs()
        USalign(structures=structures, structure_1=pairs.tables.pairs.variant, structure_2=pairs.tables.pairs.parent)
        USalign(structures=structures, structure_1=pairs.tables.pairs.variant, structure_2=pairs.tables.pairs.parent,
                non_sequential_variant="sNS")
        TMalign(structures=structures, structure_1=pairs.tables.pairs.variant, structure_2=pairs.tables.pairs.parent)
    fns, sns, tm = [t for t in pipeline.tools if t.TOOL_NAME in ("USalign", "TMalign")]
    assert [m["label"] for m in fns.modes()] == ["sequential", "circular_permutation", "non_sequential"]
    assert [m["args"] for m in fns.modes()] == [["-mm", "0"], ["-mm", "3"], ["-mm", "5"]]
    assert sns.modes()[2]["args"] == ["-mm", "6"]
    assert [m["args"] for m in tm.modes()] == [[], ["-cp"]]
    assert fns.BINARY == "USalign" and tm.BINARY == "TMalign"


@pytest.mark.parametrize("kwargs, message", [
    ({"non_sequential_variant": "other"}, "non_sequential_variant"),
    ({"non_sequential_variant": ""}, "non_sequential_variant"),
])
def test_a_bad_variant_is_refused_at_configuration(local_config, isolated_cwd, new_pipeline, kwargs, message):
    from biopipelines import USalign

    with new_pipeline("pa_bad_variant"):
        structures, pairs = _inputs()
        with pytest.raises(ValueError, match=message):
            USalign(structures=structures, structure_1=pairs.tables.pairs.variant,
                    structure_2=pairs.tables.pairs.parent, **kwargs)


def test_pair_columns_must_be_references_into_one_table(local_config, isolated_cwd, new_pipeline):
    from biopipelines import TMalign
    from biopipelines.mock import Mock

    with new_pipeline("pa_bad_pairs"):
        structures, pairs = _inputs()
        other = Mock(ids=["q"], tables={"pairs": {"columns": ["parent", "variant"]}})
        with pytest.raises(ValueError, match="must reference a column of the pair table"):
            TMalign(structures=structures, structure_1="variant", structure_2=pairs.tables.pairs.parent)
        with pytest.raises(ValueError, match="same pair table"):
            TMalign(structures=structures, structure_1=other.tables.pairs.variant,
                    structure_2=pairs.tables.pairs.parent)
        TMalign(structures=structures, structure_1=(pairs.tables.pairs, "variant"),
                structure_2=(pairs.tables.pairs, "parent"))


def test_structures_must_be_a_stream(local_config, isolated_cwd, new_pipeline):
    from biopipelines import USalign

    with new_pipeline("pa_bad_structures"):
        _, pairs = _inputs()
        with pytest.raises(ValueError, match="structures must be DataStream or StandardizedOutput"):
            USalign(structures="a.pdb", structure_1=pairs.tables.pairs.variant, structure_2=pairs.tables.pairs.parent)


def _step_text(script_path, tool):
    step = sorted(glob.glob(os.path.join(os.path.dirname(script_path), f"*_{tool}.sh")))[0]
    return open(step, encoding="utf-8", errors="replace").read()


def test_the_emitted_script_runs_the_pipe_and_the_config_carries_the_references(local_config, isolated_cwd, new_pipeline):
    from biopipelines import USalign

    pipeline = new_pipeline("pa_script")
    with pipeline:
        structures, pairs = _inputs()
        USalign(structures=structures, structure_1=pairs.tables.pairs.variant, structure_2=pairs.tables.pairs.parent)
        script_path = pipeline.save()
    text = _step_text(script_path, "USalign")
    assert "pipe_pair_alignment.py" in text and "--container-prefix" in text
    assert "--alignments-csv" in text and "--local-missing-csv" in text and "--delta-csv" not in text
    tool = [t for t in pipeline.tools if t.TOOL_NAME == "USalign"][0]
    config = yaml.safe_load(open(tool.config_yaml))
    assert config["binary"] == "USalign"
    assert config["structure_1"].endswith(":variant") and config["structure_2"].endswith(":parent")
    assert config["structure_1"].startswith("TABLE_REFERENCE:")
    assert [m["label"] for m in config["modes"]] == ["sequential", "circular_permutation", "non_sequential"]


def test_every_pipe_flag_has_a_matching_argument(local_config, isolated_cwd, new_pipeline):
    from biopipelines import TMalign, USalign

    pipeline = new_pipeline("pa_flags")
    with pipeline:
        structures, pairs = _inputs()
        USalign(structures=structures, structure_1=pairs.tables.pairs.variant, structure_2=pairs.tables.pairs.parent)
        TMalign(structures=structures, structure_1=pairs.tables.pairs.variant, structure_2=pairs.tables.pairs.parent)
        script_path = pipeline.save()
    declared = set(re.findall(r'add_argument\("(--[a-z-]+)"', open(SCRIPT, encoding="utf-8").read()))
    emitted = {}
    for tool in ("USalign", "TMalign"):
        text = _step_text(script_path, tool)
        call = text[text.index("pipe_pair_alignment.py"):]
        call = call[:call.index("python ", 1)] if "python " in call[1:] else call
        emitted[tool] = set(re.findall(r"--[a-z][a-z-]+", call))
    assert emitted["TMalign"] == declared
    assert emitted["USalign"] == declared - {"--delta-csv"}


def test_tmalign_installs_through_the_usalign_environment(local_config, monkeypatch):
    from biopipelines import TMalign, USalign, base_config

    monkeypatch.setattr(base_config.ConfigManager, "get_environment", lambda self, name: "usalign")
    folders = {"biopipelines": "/repo"}
    assert TMalign._install_script(folders) == USalign._install_script(folders)
    text = USalign._install_script(folders)
    assert "USalign -h" in text and "TMalign -h" in text and "INSTALL_SUCCESS" in text


def test_the_two_wrappers_share_one_environment(local_config):
    from biopipelines import TMalign, USalign
    assert USalign.ENV_NAME == TMalign.ENV_NAME == "usalign"


# ── splitting a big pair table across sibling jobs ────────────────────────────

def test_a_pair_table_splits_with_panda_and_pool_across_parallel_jobs(local_config, isolated_cwd, new_pipeline):
    """The recipe in the USalign documentation: one Panda slice and one USalign per sibling job, merged by Pool."""
    from biopipelines import Pool, USalign
    from biopipelines.panda import Panda
    from biopipelines.pipeline import Parallel, Resources

    pipeline = new_pipeline("pa_split")
    with pipeline:
        structures, pairs = _inputs()
        Resources(time="00:10:00", memory="1GB")
        runs = []
        with Parallel():
            for k in range(3):
                Resources(time="00:30:00", memory="2GB")
                part = Panda(tables=pairs.tables.pairs, operations=[Panda.filter(f"index % 3 == {k}")])
                runs.append(USalign(structures=structures, structure_1=part.tables.result.variant,
                                    structure_2=part.tables.result.parent))
        Resources(time="00:10:00", memory="1GB")
        merged = Pool(runs=runs, keep_ids=True)
        script_path = pipeline.save()

    assert "alignments" in merged.tables._tables and "missing" in merged.tables._tables
    slices = [t for t in pipeline.tools if t.TOOL_NAME == "USalign"]
    assert len(slices) == 3
    references = {yaml.safe_load(open(t.config_yaml))["structure_1"] for t in slices}
    assert len(references) == 3
    assert all("_Panda" in ref_ or "Panda" in ref_ for ref_ in references)
    assert sorted(glob.glob(os.path.join(os.path.dirname(script_path), "*_USalign.sh")))
