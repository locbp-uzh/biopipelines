"""Unit tests for the ReMap tool — several `onto` specifications.

These test the pure-logic path: given a source StandardizedOutput with per-ID
files, verify that ReMap computes the correct `id_mapping` for each onto mode
(basename, list, dict, DataStream). No pipeline / no bash-script emission.
"""

import os


def _make_source(tmp_path, ids, stream_name="structures"):
    """Build a minimal StandardizedOutput with per-ID .pdb files as the source."""
    from biopipelines.datastream import DataStream, create_map_table
    from biopipelines.base_config import StandardizedOutput

    tmp_path.mkdir(parents=True, exist_ok=True)
    files = [str(tmp_path / f"{sid}.pdb") for sid in ids]
    for f in files:
        open(f, "w").close()

    map_table = str(tmp_path / f"{stream_name}_map.csv")
    create_map_table(map_table, ids=ids, files=files)

    stream = DataStream(
        name=stream_name, ids=list(ids), files=files,
        format="pdb", map_table=map_table,
    )
    return StandardizedOutput({stream_name: stream})


def test_remap_onto_basename(tmp_path, record_case):
    """ReMap(source=..., onto='design') — basename + auto-numbering."""
    from biopipelines.remap import ReMap

    src = _make_source(tmp_path, ids=["a", "b", "c"])
    remapped = ReMap(source=src, onto="design")

    mapping_vals = list(remapped.id_mapping.values())
    record_case(input="ReMap onto='design' (3 IDs)",
                expected="3 entries mapped to design_* pattern",
                actual=f"{len(remapped.id_mapping)} entries, values={mapping_vals}")
    assert len(remapped.id_mapping) == 3
    assert all(v.startswith("design_") for v in mapping_vals)


def test_remap_onto_list(tmp_path, record_case):
    """ReMap(source=..., onto=[...]) — explicit new-ID list, length match."""
    from biopipelines.remap import ReMap

    src = _make_source(tmp_path, ids=["a", "b"])
    remapped = ReMap(source=src, onto=["alpha", "beta"])

    record_case(input="ReMap onto=['alpha','beta']",
                expected={"a": "alpha", "b": "beta"},
                actual=dict(remapped.id_mapping))
    assert remapped.id_mapping == {"a": "alpha", "b": "beta"}


def test_remap_onto_dict(tmp_path, record_case):
    """ReMap(source=..., onto={'a':'A', 'b':'B'}) — selective old->new mapping."""
    from biopipelines.remap import ReMap

    src = _make_source(tmp_path, ids=["a", "b"])
    remapped = ReMap(source=src, onto={"a": "A", "b": "B"})

    record_case(input="ReMap onto=dict",
                expected={"a": "A", "b": "B"},
                actual=dict(remapped.id_mapping))
    assert remapped.id_mapping == {"a": "A", "b": "B"}


def test_remap_onto_dict_partial(tmp_path, record_case):
    """ReMap(source=..., onto={'a':'A'}) — partial mapping only renames matching IDs."""
    from biopipelines.remap import ReMap

    src = _make_source(tmp_path, ids=["a", "b", "c"])
    remapped = ReMap(source=src, onto={"a": "A"})

    record_case(input="ReMap partial dict onto={'a':'A'}",
                expected={"a": "A"},
                actual=dict(remapped.id_mapping))
    assert remapped.id_mapping == {"a": "A"}


def test_remap_onto_datastream(tmp_path, record_case):
    """ReMap(source=seq_a, onto=stream_b) — align onto another DataStream."""
    from biopipelines.remap import ReMap

    src = _make_source(tmp_path / "a", ids=["a", "b"])
    dst = _make_source(tmp_path / "b", ids=["x", "y"])

    (tmp_path / "a").mkdir(exist_ok=True)
    (tmp_path / "b").mkdir(exist_ok=True)
    # _make_source already wrote into those dirs (mkdir no-op here is fine)

    remapped = ReMap(source=src, onto=dst.streams.structures)

    record_case(input="ReMap onto=DataStream(ids=[x,y])",
                expected={"a": "x", "b": "y"},
                actual=dict(remapped.id_mapping))
    assert remapped.id_mapping == {"a": "x", "b": "y"}


def test_remap_onto_standardized_output(tmp_path, record_case):
    """ReMap(source=a, onto=b) — align source onto another StandardizedOutput."""
    from biopipelines.remap import ReMap

    src = _make_source(tmp_path / "a", ids=["a", "b"])
    dst = _make_source(tmp_path / "b", ids=["x", "y"])

    remapped = ReMap(source=src, onto=dst)

    record_case(input="ReMap onto=StandardizedOutput(ids=[x,y])",
                expected={"a": "x", "b": "y"},
                actual=dict(remapped.id_mapping))
    assert remapped.id_mapping == {"a": "x", "b": "y"}


# ── ReMap driven by Mock — script-emission smoke ──────────────────────────────

def test_remap_via_mock_script_emitted(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """ReMap downstream of Mock emits a valid bash script with both tools."""
    from biopipelines.mock import Mock
    from biopipelines.remap import ReMap

    pipeline = new_pipeline("remap_mock_script")
    with pipeline:
        m = Mock(
            ids=["p1", "p2", "p3"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            map_table_strategy="config",
        )
        ReMap(source=m, onto="design")
        script_path = pipeline.save()

    record_case(input="Mock → ReMap onto='design'",
                expected="script with Mock+ReMap",
                actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Mock", "ReMap", "remap_mock_script")


# ── ReMap driven by basic input types — script-emission smoke ─────────────────
# These tests cover ReMap downstream of PDB / Ligand / Sequence (the non-Mock
# "basic input" tools). Config-time wiring only — no network fetches or tool
# binaries are exercised; execution happens at pipeline.sh run time.

def _wrap_stream(stream):
    """Wrap a single DataStream into a StandardizedOutput so ReMap can consume it.

    ReMap's ``source`` must be a StandardizedOutput. Some basic-input tools (PDB,
    Ligand) expose multiple streams of differing lengths (e.g. one sequences.csv
    vs. N structure files), which makes the raw ``tool.output`` unusable as a
    ReMap source. Isolating the per-ID stream avoids that mismatch.
    """
    from biopipelines.base_config import StandardizedOutput
    return StandardizedOutput({stream.name: stream})


def test_remap_source_pdb_script_emitted(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """ReMap downstream of PDB emits a valid script with PDB+ReMap."""
    from biopipelines.pdb import PDB
    from biopipelines.remap import ReMap

    pipeline = new_pipeline("remap_pdb_script")
    with pipeline:
        pdb = PDB(pdbs={"p1": "4ufc", "p2": "1ake"})
        ReMap(source=_wrap_stream(pdb.streams.structures), onto="design")
        script_path = pipeline.save()

    record_case(input="PDB.streams.structures → ReMap onto='design'",
                expected="script with PDB+ReMap",
                actual=os.path.basename(script_path))
    assert_valid_script(script_path, "PDB", "ReMap", "remap_pdb_script")


def test_remap_source_ligand_script_emitted(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """ReMap downstream of Ligand emits a valid script with Ligand+ReMap."""
    from biopipelines.ligand import Ligand
    from biopipelines.remap import ReMap

    pipeline = new_pipeline("remap_ligand_script")
    with pipeline:
        lig = Ligand(lookup={"l1": "ATP", "l2": "GDP"})
        ReMap(source=_wrap_stream(lig.streams.structures), onto="design")
        script_path = pipeline.save()

    record_case(input="Ligand.streams.structures → ReMap onto='design'",
                expected="script with Ligand+ReMap",
                actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Ligand", "ReMap", "remap_ligand_script")


def test_remap_source_sequence_script_emitted(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """ReMap downstream of Sequence emits a valid script with Sequence+ReMap."""
    from biopipelines.sequence import Sequence
    from biopipelines.remap import ReMap

    pipeline = new_pipeline("remap_sequence_script")
    with pipeline:
        seq = Sequence(seq={"s1": "MKVLWAALLV", "s2": "AETGFTSKLE"})
        ReMap(source=_wrap_stream(seq.streams.sequences), onto="design")
        script_path = pipeline.save()

    record_case(input="Sequence.streams.sequences → ReMap onto='design'",
                expected="script with Sequence+ReMap",
                actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Sequence", "ReMap", "remap_sequence_script")
