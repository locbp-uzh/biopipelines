"""Unit tests for biopipelines.stream_slicers — format-aware shared-file slicing."""

import pytest

from biopipelines.stream_slicers import (
    available_formats,
    get_slicer,
    register,
)


# ── registry ──────────────────────────────────────────────────────────────────

def test_available_formats_includes_fasta_and_csv(record_case):
    fmts = available_formats()
    record_case(input="available_formats()", expected=({"fasta", "fa", "csv"} <= set(fmts)),
                actual=fmts)
    assert "fasta" in fmts
    assert "fa" in fmts
    assert "csv" in fmts


def test_get_slicer_case_insensitive(record_case):
    fn1 = get_slicer("FASTA")
    fn2 = get_slicer("fasta")
    record_case(input="case-insensitive lookup", expected=True, actual=fn1 is fn2)
    assert fn1 is fn2


def test_get_slicer_unknown_raises(record_case):
    record_case(input="format=xyz", expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError, match="No shared-file slicer registered"):
        get_slicer("xyz")


# ── fasta ─────────────────────────────────────────────────────────────────────

_FASTA_SRC = """>a description text
SEQA1
SEQA2
>b
SEQB
>c another
SEQC
"""


def test_fasta_slice_keeps_subset_in_order(tmp_path, record_case):
    src = tmp_path / "src.fasta"
    dest = tmp_path / "out.fasta"
    src.write_text(_FASTA_SRC)

    get_slicer("fasta")(str(src), str(dest), ["c", "a"])

    out = dest.read_text()
    # Original order preserved: a then c
    expected_starts = [">a", ">c"]
    headers = [ln for ln in out.splitlines() if ln.startswith(">")]
    actual_starts = [h.split()[0] for h in headers]
    record_case(input="keep [c,a] from [a,b,c]",
                expected=expected_starts, actual=actual_starts)
    assert actual_starts == expected_starts
    assert "SEQA1" in out and "SEQA2" in out
    assert "SEQC" in out
    assert "SEQB" not in out
    assert ">b" not in out


def test_fasta_slice_empty_keep(tmp_path, record_case):
    src = tmp_path / "src.fasta"
    dest = tmp_path / "out.fasta"
    src.write_text(_FASTA_SRC)

    get_slicer("fa")(str(src), str(dest), [])

    out = dest.read_text()
    record_case(input="keep []", expected="", actual=out)
    assert out == ""


def test_fasta_slice_missing_ids_silently_skipped(tmp_path, record_case):
    """IDs not in the source file are silently ignored (caller decides policy)."""
    src = tmp_path / "src.fasta"
    dest = tmp_path / "out.fasta"
    src.write_text(_FASTA_SRC)

    get_slicer("fasta")(str(src), str(dest), ["a", "ghost"])

    out = dest.read_text()
    headers = [ln for ln in out.splitlines() if ln.startswith(">")]
    record_case(input="keep [a, ghost]", expected=1, actual=len(headers))
    assert len(headers) == 1
    assert headers[0].startswith(">a")


def test_fasta_slice_creates_parent_dir(tmp_path, record_case):
    src = tmp_path / "src.fasta"
    dest = tmp_path / "nested" / "deep" / "out.fasta"
    src.write_text(_FASTA_SRC)

    get_slicer("fasta")(str(src), str(dest), ["b"])

    record_case(input="nested dest", expected=True, actual=dest.exists())
    assert dest.exists()


# Header tokenizer hardening: comma/semicolon terminate the id, pipe is kept
# as part of the id (NCBI-style compound ids stay intact). Whitespace already
# tested above. These cases mirror real upstream headers from raw ProteinMPNN /
# LigandMPNN dumps (`>4LCD, score=...`) and NCBI-style fastas (`>sp|P12345|FOO`).

_FASTA_RICH = """>4LCD, score=1.5023, designed_chains=['E']
SEQA1
>T=0.1, sample=1, score=0.8533
SEQB
>sp|P12345|MYPROT_HUMAN My protein OS=Homo sapiens
SEQC
>id_semi;extra=meta
SEQD
"""


def test_fasta_slice_handles_comma_terminator(tmp_path, record_case):
    """`>4LCD, score=...` headers — id is `4LCD`, comma is not part of it."""
    src = tmp_path / "src.fasta"
    dest = tmp_path / "out.fasta"
    src.write_text(_FASTA_RICH)

    get_slicer("fasta")(str(src), str(dest), ["4LCD"])
    out = dest.read_text()
    headers = [ln for ln in out.splitlines() if ln.startswith(">")]
    record_case(input="keep [4LCD] from comma-separated header",
                expected=1, actual=len(headers))
    assert len(headers) == 1
    assert headers[0].startswith(">4LCD,")
    assert "SEQA1" in out


def test_fasta_slice_keeps_pipe_in_id(tmp_path, record_case):
    """NCBI-style `>sp|P12345|MYPROT_HUMAN ...` — full pipe-delimited id."""
    src = tmp_path / "src.fasta"
    dest = tmp_path / "out.fasta"
    src.write_text(_FASTA_RICH)

    get_slicer("fasta")(str(src), str(dest), ["sp|P12345|MYPROT_HUMAN"])
    out = dest.read_text()
    headers = [ln for ln in out.splitlines() if ln.startswith(">")]
    record_case(input="keep NCBI-style pipe id",
                expected=1, actual=len(headers))
    assert len(headers) == 1
    assert "SEQC" in out


def test_fasta_slice_handles_semicolon_terminator(tmp_path, record_case):
    src = tmp_path / "src.fasta"
    dest = tmp_path / "out.fasta"
    src.write_text(_FASTA_RICH)

    get_slicer("fasta")(str(src), str(dest), ["id_semi"])
    out = dest.read_text()
    headers = [ln for ln in out.splitlines() if ln.startswith(">")]
    record_case(input="keep [id_semi] from semicolon header",
                expected=1, actual=len(headers))
    assert len(headers) == 1
    assert "SEQD" in out


def test_fasta_slice_rename_preserves_suffix(tmp_path, record_case):
    """Renaming a kept record must keep the everything-after-id intact."""
    src = tmp_path / "src.fasta"
    dest = tmp_path / "out.fasta"
    src.write_text(_FASTA_RICH)

    get_slicer("fasta")(str(src), str(dest), ["4LCD"],
                       rename_map={"4LCD": "design_1"})
    out = dest.read_text()
    # Header should now be `>design_1, score=..., designed_chains=...`
    headers = [ln for ln in out.splitlines() if ln.startswith(">")]
    record_case(input="rename 4LCD -> design_1, comma-rich header",
                expected=(True, True, True),
                actual=(headers[0].startswith(">design_1,"),
                        "score=1.5023" in headers[0],
                        "designed_chains" in headers[0]))
    assert len(headers) == 1
    assert headers[0].startswith(">design_1,")
    assert "score=1.5023" in headers[0]
    assert "designed_chains" in headers[0]


def test_fasta_slice_rename_preserves_pipe_id_suffix(tmp_path, record_case):
    src = tmp_path / "src.fasta"
    dest = tmp_path / "out.fasta"
    src.write_text(_FASTA_RICH)

    get_slicer("fasta")(str(src), str(dest), ["sp|P12345|MYPROT_HUMAN"],
                       rename_map={"sp|P12345|MYPROT_HUMAN": "myprot"})
    out = dest.read_text()
    headers = [ln for ln in out.splitlines() if ln.startswith(">")]
    record_case(input="rename NCBI pipe id, preserve description",
                expected=(True, True),
                actual=(headers[0].startswith(">myprot "),
                        "OS=Homo sapiens" in headers[0]))
    assert headers[0].startswith(">myprot ")
    assert "OS=Homo sapiens" in headers[0]


# ── csv ───────────────────────────────────────────────────────────────────────

def test_csv_slice_keeps_subset(tmp_path, record_case):
    import pandas as pd

    src = tmp_path / "src.csv"
    dest = tmp_path / "out.csv"
    pd.DataFrame({"id": ["a", "b", "c"], "x": [1, 2, 3]}).to_csv(src, index=False)

    get_slicer("csv")(str(src), str(dest), ["a", "c"])

    out = pd.read_csv(dest)
    record_case(input="keep [a,c] from [a,b,c]",
                expected=(["a", "c"], [1, 3]),
                actual=(out["id"].tolist(), out["x"].tolist()))
    assert out["id"].tolist() == ["a", "c"]
    assert out["x"].tolist() == [1, 3]


def test_csv_slice_requires_id_column(tmp_path, record_case):
    import pandas as pd
    src = tmp_path / "src.csv"
    dest = tmp_path / "out.csv"
    pd.DataFrame({"name": ["a", "b"], "x": [1, 2]}).to_csv(src, index=False)

    record_case(input="csv without id column", expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError, match="requires an 'id' column"):
        get_slicer("csv")(str(src), str(dest), ["a"])


# ── registration ──────────────────────────────────────────────────────────────

def test_register_new_format(record_case):
    @register("test_xyz_format")
    def _impl(src, dest, kept_ids, rename_map=None):
        pass

    fn = get_slicer("TEST_XYZ_FORMAT")
    record_case(input="register('test_xyz_format')", expected=True, actual=fn is _impl)
    assert fn is _impl


def test_three_arg_slicer_still_callable(tmp_path, record_case):
    """A slicer registered without rename_map still works for callers that
    invoke with rename_map=None (the new default)."""
    @register("test_three_arg")
    def _three(src, dest, kept_ids):
        open(dest, "w").close()

    src = tmp_path / "a"; src.write_text("hi")
    dest = tmp_path / "b"
    fn = get_slicer("test_three_arg")
    fn(str(src), str(dest), ["x"])
    fn(str(src), str(dest), ["x"], rename_map=None)
    record_case(input="3-arg slicer with rename_map=None",
                expected=True, actual=dest.exists())
    assert dest.exists()


def test_three_arg_slicer_rejects_real_rename_map(tmp_path, record_case):
    """A slicer that doesn't know how to rename must raise when asked to."""
    @register("test_three_arg_b")
    def _three(src, dest, kept_ids):
        open(dest, "w").close()

    src = tmp_path / "a"; src.write_text("hi")
    dest = tmp_path / "b"
    fn = get_slicer("test_three_arg_b")
    record_case(input="rename_map={'x':'y'} on 3-arg slicer",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError, match="does not support rename_map"):
        fn(str(src), str(dest), ["x"], rename_map={"x": "y"})
