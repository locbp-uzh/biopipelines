"""Unit tests for biopipelines.datastream — __len__, __getitem__, iteration.

Also covers StandardizedOutput ID-based selection (c302877). Parametrized
so the XLSX report shows input IDs, expected slice, and actual slice.
"""

import pytest

from biopipelines.datastream import DataStream
from biopipelines.base_config import StandardizedOutput


# ── __len__ ───────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("ids, expected_len", [
    (["a", "b", "c"],                 3),
    (["5HG6_<0..49>"],                 50),
    (["prot_<0..4>[_<N><A V>]"],       5),   # lazy: prefix count only
    (["<0..1>_<A B>"],                 4),
])
def test_datastream_len(record_case, ids, expected_len):
    ds = DataStream(name="s", ids=list(ids))
    actual = len(ds)
    record_case(input=ids, expected=expected_len, actual=actual)
    assert actual == expected_len


# ── integer indexing ──────────────────────────────────────────────────────────

@pytest.mark.parametrize("ids, index, expected_id", [
    (["a", "b", "c"],      0, "a"),
    (["a", "b", "c"],      2, "c"),
    # Exercises the A3-fixed expand_at code path for single-pattern streams:
    (["base_<0..2>"],      0, "base_0"),
    (["base_<0..2>"],      1, "base_1"),
    (["base_<0..2>"],      2, "base_2"),
])
def test_datastream_integer_index(record_case, ids, index, expected_id):
    ds = DataStream(name="s", ids=list(ids))
    picked = ds[index]
    actual = picked.ids[0]
    record_case(input=(ids, index), expected=expected_id, actual=actual)
    assert actual == expected_id
    assert isinstance(picked, DataStream)


def test_datastream_integer_index_out_of_range(record_case):
    ds = DataStream(name="s", ids=["a", "b"])
    record_case(input=("ids=[a,b]", 5), expected="IndexError", actual="IndexError")
    with pytest.raises(IndexError):
        _ = ds[5]


def test_datastream_file_template_substituted(record_case):
    ds = DataStream(name="s", ids=["base_<0..2>"], files=["<id>.pdb"], format="pdb")
    first = ds[0]
    actual = (first.ids[0], first.files[0])
    expected = ("base_0", "base_0.pdb")
    record_case(input="ids=[base_<0..2>], files=[<id>.pdb]", expected=expected, actual=actual)
    assert actual == expected


# ── slice indexing ────────────────────────────────────────────────────────────

@pytest.mark.parametrize("ids, slc, expected_ids", [
    (["a", "b", "c", "d"],  slice(1, 3),    ["b", "c"]),
    (["base_<0..4>"],       slice(1, 4),    ["base_1", "base_2", "base_3"]),
    (["a", "b"],            slice(5, 10),   []),
])
def test_datastream_slice(record_case, ids, slc, expected_ids):
    ds = DataStream(name="s", ids=list(ids))
    sub = ds[slc]
    actual = list(sub.ids) if not sub.has_patterns() else list(sub.ids_expanded)
    # Slices always store literal IDs (DataStream stores ids_expanded[slc]).
    record_case(input=(ids, f"[{slc.start}:{slc.stop}]"),
                expected=expected_ids, actual=actual)
    assert actual == expected_ids


# ── iteration ─────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("ids, expected_iter_ids", [
    (["a", "b", "c"],        ["a", "b", "c"]),
    (["x_<0..2>"],           ["x_0", "x_1", "x_2"]),
])
def test_datastream_iteration(record_case, ids, expected_iter_ids):
    ds = DataStream(name="s", ids=list(ids))
    actual = [x.ids[0] for x in ds]
    record_case(input=ids, expected=expected_iter_ids, actual=actual)
    assert actual == expected_iter_ids


def test_iteration_count_matches_len(record_case):
    ds = DataStream(name="s", ids=["x_<0..9>"])
    actual = sum(1 for _ in ds)
    record_case(input="ids=[x_<0..9>]", expected=10, actual=actual)
    assert actual == len(ds) == 10


@pytest.mark.parametrize(
    "stream_format, allowed, expected",
    [
        ("pdb", ("pdb", "cif"), True),
        ("cif", ("pdb", "cif"), True),
        ("pdb|cif", ("pdb", "cif"), True),
        ("pdb|sdf", ("pdb", "cif"), False),
        ("sdf", ("pdb", "cif"), False),
    ],
)
def test_datastream_has_only_formats(record_case, stream_format, allowed, expected):
    ds = DataStream(name="s", ids=["x"], format=stream_format)
    actual = ds.has_only_formats(*allowed)
    record_case(input=(stream_format, allowed), expected=expected, actual=actual)
    assert actual is expected


def test_datastream_format_helpers_normalize_tokens(record_case):
    ds = DataStream(name="s", ids=["x"], format=" PDB | cif ")
    actual = (ds.formats, ds.has_format("pdb"), ds.has_only_formats("pdb", "cif"))
    expected = (("pdb", "cif"), True, True)
    record_case(input="format=' PDB | cif '", expected=expected, actual=actual)
    assert actual == expected


# ── StandardizedOutput ID-based selection (c302877) ───────────────────────────

def test_standardized_output_select_by_id(record_case):
    seqs = DataStream(name="sequences", ids=["CP1", "CP2", "CP3"])
    out = StandardizedOutput({"sequences": seqs})
    picked = out["CP2"]
    actual = list(picked.streams.sequences.ids)
    record_case(input="ids=[CP1,CP2,CP3], select=CP2", expected=["CP2"], actual=actual)
    assert actual == ["CP2"]
    assert isinstance(picked, StandardizedOutput)


def test_standardized_output_missing_id_raises(record_case):
    seqs = DataStream(name="sequences", ids=["CP1"])
    out = StandardizedOutput({"sequences": seqs})
    record_case(input="ids=[CP1], select=NOPE", expected="KeyError", actual="KeyError")
    with pytest.raises(KeyError):
        _ = out["NOPE"]


def test_standardized_output_multi_stream_selection(record_case):
    ids = ["a", "b", "c"]
    seqs = DataStream(name="sequences", ids=list(ids))
    structs = DataStream(
        name="structures",
        ids=list(ids),
        files=["a.pdb", "b.pdb", "c.pdb"],
        format="pdb",
    )
    out = StandardizedOutput({"sequences": seqs, "structures": structs})
    picked = out["b"]
    actual = (
        list(picked.streams.sequences.ids),
        list(picked.streams.structures.ids),
        picked.streams.structures.files,
    )
    expected = (["b"], ["b"], ["b.pdb"])
    record_case(input="2 streams, select=b", expected=expected, actual=actual)
    assert actual == expected


# ── Mock-driven: real files + map_table wiring ────────────────────────────────

def test_datastream_len_and_slice_via_mock(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """Mock-produced DataStream supports len, indexing, and slicing."""
    from biopipelines.mock import Mock

    pipeline = new_pipeline("ds_slice_mock")
    with pipeline:
        m = Mock(
            ids=["a", "b", "c", "d"],
            streams={"s": {"format": "pdb", "file": "<id>.pdb"}},
        )
        pipeline.save()

    s = m.streams.s
    sub = s[1:3]
    record_case(input="Mock ids=[a..d], [1:3]",
                expected=(4, ["b", "c"]),
                actual=(len(s), list(sub.ids)))
    assert len(s) == 4
    assert list(sub.ids) == ["b", "c"]


def test_datastream_from_mock_has_map_table_at_config_time(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """With strategy='config', Mock writes a map_table at config time."""
    import os
    from biopipelines.mock import Mock

    pipeline = new_pipeline("ds_maptable_mock")
    with pipeline:
        m = Mock(
            ids=["a", "b"],
            streams={"s": {"format": "pdb", "file": "<id>.pdb"}},
            map_table_strategy="config",
        )
        pipeline.save()

    map_path = m.streams.s.map_table
    record_case(input="map_table_strategy=config",
                expected=("exists", True),
                actual=("exists", os.path.exists(map_path)))
    assert os.path.exists(map_path)


# ── shared-file DataStream (files: str) ───────────────────────────────────────

def test_datastream_shared_construction(record_case):
    """files=str with N>1 ids is a valid shared-file stream."""
    ds = DataStream(name="fasta", ids=["a", "b", "c"], files="/tmp/shared.fasta",
                    map_table="", format="fasta")
    record_case(input="files='/tmp/shared.fasta', ids=[a,b,c]",
                expected=(True, "/tmp/shared.fasta"),
                actual=(ds.is_shared_file, ds.files))
    assert ds.is_shared_file is True
    assert ds.files == "/tmp/shared.fasta"
    assert len(ds) == 3


def test_datastream_shared_files_expanded_single(record_case):
    """files_expanded returns [shared_path] (length 1), not N copies."""
    ds = DataStream(name="fasta", ids=["a", "b", "c"], files="/tmp/s.fa",
                    map_table="", format="fasta")
    actual = ds.files_expanded
    record_case(input="3 ids, shared file",
                expected=["/tmp/s.fa"], actual=actual)
    assert actual == ["/tmp/s.fa"]


def test_datastream_shared_iter_preserves_str(record_case):
    """Iterating a shared stream yields sub-streams that keep files as str."""
    ds = DataStream(name="fasta", ids=["a", "b"], files="/tmp/s.fa",
                    map_table="", format="fasta")
    subs = list(ds)
    types = [type(s.files).__name__ for s in subs]
    paths = [s.files for s in subs]
    record_case(input="iter shared", expected=(["str", "str"], ["/tmp/s.fa", "/tmp/s.fa"]),
                actual=(types, paths))
    assert types == ["str", "str"]
    assert paths == ["/tmp/s.fa", "/tmp/s.fa"]


def test_datastream_shared_getitem_preserves_str(record_case):
    """Integer and slice indexing on a shared stream preserves str form."""
    ds = DataStream(name="fasta", ids=["a", "b", "c"], files="/tmp/s.fa",
                    map_table="", format="fasta")
    one = ds[1]
    two = ds[0:2]
    record_case(input="ds[1] and ds[0:2]",
                expected=("str", "str"),
                actual=(type(one.files).__name__, type(two.files).__name__))
    assert one.files == "/tmp/s.fa"
    assert two.files == "/tmp/s.fa"
    assert one.ids == ["b"]
    assert list(two.ids) == ["a", "b"]


def test_datastream_shared_to_dict_roundtrip(record_case):
    """to_dict / from_dict preserves the str form."""
    ds = DataStream(name="fasta", ids=["a", "b"], files="/tmp/s.fa",
                    map_table="", format="fasta")
    d = ds.to_dict()
    restored = DataStream.from_dict(d)
    record_case(input="roundtrip shared",
                expected=("str", "/tmp/s.fa"),
                actual=(type(restored.files).__name__, restored.files))
    assert restored.is_shared_file
    assert restored.files == "/tmp/s.fa"


def test_datastream_shared_filter_by_ids(record_case):
    """filter_by_ids on a shared stream keeps the shared path."""
    ds = DataStream(name="fasta", ids=["a", "b", "c"], files="/tmp/s.fa",
                    map_table="", format="fasta")
    sub = ds.filter_by_ids(["b", "c"])
    record_case(input="keep [b,c] from [a,b,c]",
                expected=(["b", "c"], "/tmp/s.fa"),
                actual=(list(sub.ids), sub.files))
    assert list(sub.ids) == ["b", "c"]
    assert sub.files == "/tmp/s.fa"


def test_create_map_table_shared_file(tmp_path, record_case):
    """create_map_table replicates a shared str path into every row."""
    from biopipelines.datastream import create_map_table
    import pandas as pd

    out = tmp_path / "shared_map.csv"
    create_map_table(str(out), ids=["a", "b", "c"], files="/tmp/s.fa")
    df = pd.read_csv(out)
    record_case(input="shared file 3 ids",
                expected=(3, ["/tmp/s.fa"] * 3),
                actual=(len(df), df["file"].tolist()))
    assert len(df) == 3
    assert df["file"].tolist() == ["/tmp/s.fa"] * 3
    assert df["id"].tolist() == ["a", "b", "c"]
