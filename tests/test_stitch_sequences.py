"""Tests for StitchSequences marker-fill mode.

Marker keys (a substitutions key with no digits, e.g. "X") fill every template
position whose residue is one of the marker chars from a reference sequence at
the same index — used to repair gap-marked sequences (the 'X' padding from a
structure with missing residues) against a reference.

Two layers:
- Config-time: the wrapper classifies marker vs position keys, tags the marker
  in its config, predicts output ids, and rejects marker keys for indels.
- Runtime: run pipe_stitch_sequences.py end-to-end via a config JSON and inspect
  the produced sequences.csv / missing.csv — fill from a stream reference (matched
  by the template's own id), fill from a raw reference, and an unfilled marker
  position (source also a marker) left as-is and recorded in missing.
"""

import json
import os
import subprocess
import sys

import pandas as pd
import pytest


def _repo_root():
    import biopipelines.stitch_sequences as _s
    return os.path.dirname(os.path.dirname(os.path.abspath(_s.__file__)))


def _run_pipe(tmp_path, substitutions, template_seq, template_id="p1",
              indels=None, remove_duplicates=True):
    """Write a stitch config JSON and run the pipe script; return (rc, out_df, missing_df)."""
    tmpl = tmp_path / "tmpl.csv"
    pd.DataFrame([{"id": template_id, "sequence": template_seq}]).to_csv(tmpl, index=False)

    out_csv = tmp_path / "out.csv"
    missing_csv = tmp_path / "missing.csv"
    cfg = {
        "template": {"type": "tool_output", "sequences_file": str(tmpl),
                     "source_name": "t", "sequence_ids": [template_id]},
        "substitutions": substitutions,
        "indels": indels or {},
        "remove_duplicates": remove_duplicates,
        "output_csv": str(out_csv),
        "missing_csv": str(missing_csv),
        "step_tool_name": "StitchSequences",
        "upstream_missing_paths": [],
    }
    cfg_file = tmp_path / "cfg.json"
    cfg_file.write_text(json.dumps(cfg))

    script = os.path.join(_repo_root(), "pipe_scripts", "pipe_stitch_sequences.py")
    r = subprocess.run([sys.executable, script, "--config", str(cfg_file)],
                       capture_output=True, text=True)
    out_df = pd.read_csv(out_csv) if out_csv.exists() else None
    missing_df = pd.read_csv(missing_csv) if missing_csv.exists() else None
    return r.returncode, out_df, missing_df


def _marker_sub_stream(tmp_path, ref_seq, ref_id="p1", marker="X"):
    """A substitutions dict: marker key filled from a stream (CSV) reference."""
    ref = tmp_path / "ref.csv"
    pd.DataFrame([{"id": ref_id, "sequence": ref_seq}]).to_csv(ref, index=False)
    return {f"marker:{marker}": {
        "position_key": {"type": "marker", "markers": marker},
        "sequences": {"type": "tool_output", "sequences_file": str(ref),
                      "source_name": "r", "sequence_ids": [ref_id]},
    }}


def _marker_sub_raw(ref_seq, marker="X"):
    """A substitutions dict: marker key filled from a raw-string reference."""
    return {f"marker:{marker}": {
        "position_key": {"type": "marker", "markers": marker},
        "sequences": {"type": "raw", "sequences": [ref_seq]},
    }}


# ── config-time (wrapper) ────────────────────────────────────────────────────

def test_marker_key_classifier():
    from biopipelines.stitch_sequences import _is_marker_key
    assert _is_marker_key("X")
    assert _is_marker_key("XB")     # multi-char markers allowed
    assert _is_marker_key("-")
    assert not _is_marker_key("11-19")
    assert not _is_marker_key("145+147")
    assert not _is_marker_key("")


def test_marker_key_tagged_and_ids_predicted(local_config, isolated_cwd):
    from biopipelines.pipeline import Pipeline
    from biopipelines.stitch_sequences import StitchSequences
    from biopipelines.datastream import DataStream

    pipeline = Pipeline(project="TestSuite", job="stitch_marker",
                        on_the_fly=False, local_output=True, config="local")
    with pipeline:
        gappy = DataStream(name="sequences", ids=["p1"], files=[],
                           map_table="g.csv", format="csv")
        ref = DataStream(name="sequences", ids=["p1"], files=[],
                         map_table="r.csv", format="csv")
        out = StitchSequences(template=gappy, substitutions={"X": ref})

    # marker key tagged in config, not parsed as a position range
    tool = pipeline.tools[-1]
    assert tool.substitution_infos["marker:X"]["position_key"] == {
        "type": "marker", "markers": "X"}
    # one filled output per template id
    assert out.streams.sequences.ids == ["p1_1"]
    # missing table declared (carries unfilled-marker rows at runtime)
    assert hasattr(out.tables, "missing")


def test_marker_key_rejected_for_indels(local_config, isolated_cwd):
    from biopipelines.pipeline import Pipeline
    from biopipelines.stitch_sequences import StitchSequences

    Pipeline(project="TestSuite", job="stitch_indel_marker",
             on_the_fly=False, local_output=True, config="local")
    with pytest.raises(ValueError, match="only supported for substitutions"):
        StitchSequences(template="MKT", indels={"X": ["AA"]})


# ── runtime (pipe script) ────────────────────────────────────────────────────

def test_marker_fill_from_stream_reference(tmp_path):
    # template has X at positions 4 and 8; reference fills both
    rc, out, missing = _run_pipe(
        tmp_path,
        _marker_sub_stream(tmp_path, ref_seq="MKTAAGSWY"),
        template_seq="MKTXAGSXY",
    )
    assert rc == 0
    assert list(out["sequence"]) == ["MKTAAGSWY"]
    assert out["id"].iloc[0] == "p1_1"
    # nothing unfilled
    assert missing is None or missing.empty or "unfilled" not in str(missing.to_dict())


def test_marker_fill_from_raw_reference(tmp_path):
    rc, out, missing = _run_pipe(
        tmp_path,
        _marker_sub_raw("MKTLLGSAY"),
        template_seq="MKTXXGSAY",
    )
    assert rc == 0
    assert list(out["sequence"]) == ["MKTLLGSAY"]


def test_marker_unfilled_position_kept_not_removed(tmp_path):
    # template X at 4,5,8; reference is ALSO X at 5 -> position 5 unfillable
    rc, out, missing = _run_pipe(
        tmp_path,
        _marker_sub_stream(tmp_path, ref_seq="MKTAXGSWY"),
        template_seq="MKTXXGSXY",
    )
    assert rc == 0
    # 4 and 8 filled, 5 left as X — sequence is KEPT
    assert list(out["sequence"]) == ["MKTAXGSWY"]
    # crucially it must NOT be in missing.csv (that is a removal manifest; Load
    # would drop the very sequence we retained). Unfilled positions are logged only.
    assert missing is None or missing.empty or "p1_1" not in set(missing["id"])


def test_raw_string_marker_ref_predicted_id_matches_runtime(local_config, isolated_cwd, tmp_path):
    """Config-time predicted id must equal the runtime-emitted id for a raw-string
    template + raw-string marker reference (the seq vs seq_1_1 regression)."""
    from biopipelines.pipeline import Pipeline
    from biopipelines.stitch_sequences import StitchSequences

    pipeline = Pipeline(project="TestSuite", job="stitch_rawmarker",
                        on_the_fly=False, local_output=True, config="local")
    with pipeline:
        out = StitchSequences(template="ABXXEF", substitutions={"X": "ABQREF"})
    predicted = out.streams.sequences.ids

    # Run the same raw-string template + marker ref through the pipe script.
    import json, subprocess, sys
    o = tmp_path / "o2.csv"; m = tmp_path / "m2.csv"
    cfg = {"template": {"type": "raw", "sequences": ["ABXXEF"], "source_name": "t"},
           "substitutions": {"marker:X": {"position_key": {"type": "marker", "markers": "X"},
                             "sequences": {"type": "raw", "sequences": ["ABQREF"]}}},
           "indels": {}, "remove_duplicates": True, "output_csv": str(o),
           "missing_csv": str(m), "step_tool_name": "StitchSequences",
           "upstream_missing_paths": []}
    cf = tmp_path / "c2.json"; cf.write_text(json.dumps(cfg))
    script = os.path.join(_repo_root(), "pipe_scripts", "pipe_stitch_sequences.py")
    subprocess.run([sys.executable, script, "--config", str(cf)], capture_output=True, text=True)
    runtime = list(pd.read_csv(o)["id"])
    assert predicted == runtime == ["seq_1_1"]


def test_chained_stitch_predicted_id_matches_runtime(local_config, isolated_cwd, tmp_path):
    """Two-round stitching: the second round's predicted id must equal the
    runtime id. The predictor used to strip the first round's trailing _<n>
    (seq_1_1 -> predicted seq_1_1) while the runtime appends without stripping
    (seq_1_1 -> seq_1_1_1)."""
    import json, subprocess, sys
    from biopipelines.pipeline import Pipeline
    from biopipelines.stitch_sequences import StitchSequences
    from biopipelines.datastream import DataStream

    # round-1 output is a stream whose id already carries a _<n> suffix
    r1 = tmp_path / "r1.csv"
    pd.DataFrame([{"id": "seq_1_1", "sequence": "GGCDEF"}]).to_csv(r1, index=False)

    pipeline = Pipeline(project="TestSuite", job="stitch_chain",
                        on_the_fly=False, local_output=True, config="local")
    with pipeline:
        tmpl = DataStream(name="sequences", ids=["seq_1_1"], files=[],
                          map_table=str(r1), format="csv")
        out = StitchSequences(template=tmpl, substitutions={"3-4": ["TT"]})
    predicted = out.streams.sequences.ids

    # run the same round-2 config through the pipe script
    o = tmp_path / "r2.csv"; m = tmp_path / "r2m.csv"
    cfg = {"template": {"type": "tool_output", "sequences_file": str(r1),
                        "source_name": "t", "sequence_ids": ["seq_1_1"]},
           "substitutions": {"3-4": {"position_key": {"type": "fixed", "positions": "3-4"},
                             "sequences": {"type": "raw", "sequences": ["TT"]}}},
           "indels": {}, "remove_duplicates": True, "output_csv": str(o),
           "missing_csv": str(m), "step_tool_name": "StitchSequences",
           "upstream_missing_paths": []}
    cf = tmp_path / "c2.json"; cf.write_text(json.dumps(cfg))
    script = os.path.join(_repo_root(), "pipe_scripts", "pipe_stitch_sequences.py")
    subprocess.run([sys.executable, script, "--config", str(cf)], capture_output=True, text=True)
    runtime = list(pd.read_csv(o)["id"])
    assert predicted == runtime == ["seq_1_1_1"]


def test_matched_prediction_excludes_distant_sibling(local_config, isolated_cwd, tmp_path):
    """The matched-stream predictor must use the same closest_siblings_only=True
    grouping as the runtime: a distant id (Panda_2_1) is excluded for template
    Panda_29_1, so predicted ids match runtime (no phantom)."""
    import json, subprocess, sys
    from biopipelines.pipeline import Pipeline
    from biopipelines.stitch_sequences import StitchSequences
    from biopipelines.datastream import DataStream

    tmpl = tmp_path / "tmpl.csv"
    pd.DataFrame([{"id": "Panda_29_1", "sequence": "ABCDEF"}]).to_csv(tmpl, index=False)
    src = tmp_path / "src.csv"
    pd.DataFrame([{"id": "Panda_29_2", "sequence": "GGCDEF"},
                  {"id": "Panda_29_3", "sequence": "TTCDEF"},
                  {"id": "Panda_2_1", "sequence": "CCCDEF"}]).to_csv(src, index=False)

    pipeline = Pipeline(project="TestSuite", job="stitch_sib",
                        on_the_fly=False, local_output=True, config="local")
    with pipeline:
        t = DataStream(name="sequences", ids=["Panda_29_1"], files=[],
                       map_table=str(tmpl), format="csv")
        s = DataStream(name="sequences", ids=["Panda_29_2", "Panda_29_3", "Panda_2_1"],
                       files=[], map_table=str(src), format="csv")
        out = StitchSequences(template=t, substitutions={"1-2": s})
    predicted = sorted(out.streams.sequences.ids)

    o = tmp_path / "o.csv"; m = tmp_path / "m.csv"
    cfg = {"template": {"type": "tool_output", "sequences_file": str(tmpl),
                        "source_name": "t", "sequence_ids": ["Panda_29_1"]},
           "substitutions": {"1-2": {"position_key": {"type": "fixed", "positions": "1-2"},
                             "sequences": {"type": "tool_output", "sequences_file": str(src),
                             "source_name": "s",
                             "sequence_ids": ["Panda_29_2", "Panda_29_3", "Panda_2_1"]}}},
           "indels": {}, "remove_duplicates": False, "output_csv": str(o),
           "missing_csv": str(m), "step_tool_name": "StitchSequences",
           "upstream_missing_paths": []}
    cf = tmp_path / "c.json"; cf.write_text(json.dumps(cfg))
    script = os.path.join(_repo_root(), "pipe_scripts", "pipe_stitch_sequences.py")
    subprocess.run([sys.executable, script, "--config", str(cf)], capture_output=True, text=True)
    runtime = sorted(pd.read_csv(o)["id"])
    # Single tool-output axis -> output id is the matched source id itself
    # (provenance-carrying); distant sibling Panda_2_1 stays excluded.
    assert predicted == runtime == ["Panda_29_2", "Panda_29_3"]


def test_two_tool_axes_emit_plus_id_and_provenance(local_config, isolated_cwd, tmp_path):
    """Two tool-output substitution axes -> +-joined id carrying provenance, plus
    one sequences_k.id column per axis. Predictor and runtime agree, and each
    output id resolves back to the correct per-axis source id respectively."""
    import json, subprocess, sys
    from biopipelines.pipeline import Pipeline
    from biopipelines.stitch_sequences import StitchSequences
    from biopipelines.datastream import DataStream
    from biopipelines.id_map_utils import get_mapped_ids

    tmpl = tmp_path / "tmpl.csv"
    pd.DataFrame([{"id": "pdb_1_0", "sequence": "A" * 12}]).to_csv(tmpl, index=False)
    # Same id namespace across both axes — the case the +-id alone can't
    # disambiguate, so the provenance columns are what tells them apart.
    ax1 = tmp_path / "ax1.csv"
    pd.DataFrame([{"id": "pdb_1_1", "sequence": "C" * 12},
                  {"id": "pdb_1_2", "sequence": "D" * 12}]).to_csv(ax1, index=False)
    ax2 = tmp_path / "ax2.csv"
    pd.DataFrame([{"id": "pdb_1_1", "sequence": "E" * 12},
                  {"id": "pdb_1_2", "sequence": "F" * 12}]).to_csv(ax2, index=False)

    pipeline = Pipeline(project="TestSuite", job="stitch_two_axes",
                        on_the_fly=False, local_output=True, config="local")
    with pipeline:
        t = DataStream(name="sequences", ids=["pdb_1_0"], files=[],
                       map_table=str(tmpl), format="csv")
        s1 = DataStream(name="sequences", ids=["pdb_1_1", "pdb_1_2"], files=[],
                        map_table=str(ax1), format="csv")
        s2 = DataStream(name="sequences", ids=["pdb_1_1", "pdb_1_2"], files=[],
                        map_table=str(ax2), format="csv")
        out = StitchSequences(template=t, substitutions={"1-3": s1, "7-9": s2})
    predicted = sorted(out.streams.sequences.ids)
    assert out.tables.sequences.info.columns == [
        "id", "sequence", "sequences_1.id", "sequences_2.id"]

    o = tmp_path / "o.csv"; m = tmp_path / "m.csv"
    cfg = {"template": {"type": "tool_output", "sequences_file": str(tmpl),
                        "source_name": "t", "sequence_ids": ["pdb_1_0"]},
           "substitutions": {
               "1-3": {"position_key": {"type": "fixed", "positions": "1-3"},
                       "sequences": {"type": "tool_output", "sequences_file": str(ax1),
                                     "source_name": "s1",
                                     "sequence_ids": ["pdb_1_1", "pdb_1_2"]}},
               "7-9": {"position_key": {"type": "fixed", "positions": "7-9"},
                       "sequences": {"type": "tool_output", "sequences_file": str(ax2),
                                     "source_name": "s2",
                                     "sequence_ids": ["pdb_1_1", "pdb_1_2"]}}},
           "indels": {}, "remove_duplicates": False, "output_csv": str(o),
           "missing_csv": str(m), "step_tool_name": "StitchSequences",
           "upstream_missing_paths": []}
    cf = tmp_path / "c.json"; cf.write_text(json.dumps(cfg))
    script = os.path.join(_repo_root(), "pipe_scripts", "pipe_stitch_sequences.py")
    subprocess.run([sys.executable, script, "--config", str(cf)], capture_output=True, text=True)
    df = pd.read_csv(o)

    assert predicted == sorted(df["id"]) == [
        "pdb_1_1+pdb_1_1", "pdb_1_1+pdb_1_2", "pdb_1_2+pdb_1_1", "pdb_1_2+pdb_1_2"]
    assert list(df.columns) == ["id", "sequence", "sequences_1.id", "sequences_2.id"]
    # Provenance columns disambiguate the same-namespace axes by fixed position.
    row = df[df["id"] == "pdb_1_1+pdb_1_2"].iloc[0]
    assert row["sequences_1.id"] == "pdb_1_1"  # axis 1 source
    assert row["sequences_2.id"] == "pdb_1_2"  # axis 2 source
    # The +-id itself resolves to each axis's source id respectively.
    assert get_mapped_ids(["pdb_1_1+pdb_1_2"], ["pdb_1_1"])["pdb_1_1+pdb_1_2"] == "pdb_1_1"


def test_empty_option_list_rejected(local_config, isolated_cwd):
    from biopipelines.pipeline import Pipeline
    from biopipelines.stitch_sequences import StitchSequences

    Pipeline(project="TestSuite", job="stitch_empty",
             on_the_fly=False, local_output=True, config="local")
    with pytest.raises(ValueError, match="cannot be empty"):
        StitchSequences(template="ABCDEF", substitutions={"1-2": []})
    with pytest.raises(ValueError, match="cannot be empty"):
        StitchSequences(template="ABCDEF", indels={"1-2": []})


def test_marker_fill_leaves_non_marker_positions_untouched(tmp_path):
    # no X in template -> reference changes nothing
    rc, out, _ = _run_pipe(
        tmp_path,
        _marker_sub_stream(tmp_path, ref_seq="AAAAAAAAA"),
        template_seq="MKTAAGSWY",
    )
    assert rc == 0
    assert list(out["sequence"]) == ["MKTAAGSWY"]
