"""Unit tests for the Pool tool.

Pool gathers N StandardizedOutputs from parallel runs of the same upstream
tool into one combined StandardizedOutput. Tests exercise it on real
input-type tools (Sequence, Ligand) so the type-agnostic stream / table
iteration is covered against actual upstream payloads, not Mock stubs.

The runtime ``pipe_pool.py`` is exercised end-to-end where the local
pipeline can run; tests that only need to verify config-time output
shape stop at ``pipeline.save()`` and inspect the emitted DataStream
declarations directly.
"""
from __future__ import annotations

import os
from pathlib import Path

import pytest


# ── input-validation tests ────────────────────────────────────────────────────

def test_pool_requires_at_least_two_inputs(
    local_config, isolated_cwd, new_pipeline,
):
    from biopipelines.sequence import Sequence
    from biopipelines.pool import Pool

    pipeline = new_pipeline("pool_one_input")
    with pytest.raises(ValueError, match="at least 2"):
        with pipeline:
            s = Sequence(seq="MKTAY", ids="only")
            Pool(runs=[s])


def test_pool_rejects_non_standardized_output(
    local_config, isolated_cwd, new_pipeline,
):
    from biopipelines.sequence import Sequence
    from biopipelines.pool import Pool

    pipeline = new_pipeline("pool_bad_input")
    with pytest.raises(ValueError, match="StandardizedOutput"):
        with pipeline:
            s = Sequence(seq="MKTAY", ids="seq")
            Pool(runs=[s, "not_a_run"])


def test_pool_rejects_mismatched_streams(
    local_config, isolated_cwd, new_pipeline,
):
    """A Sequence (streams: sequences) and a Ligand (streams: compounds)
    do not share stream names — Pool should reject that pair."""
    from biopipelines.sequence import Sequence
    from biopipelines.ligand import Ligand
    from biopipelines.pool import Pool

    pipeline = new_pipeline("pool_mismatch")
    with pytest.raises(ValueError, match="stream names"):
        with pipeline:
            seq = Sequence(seq="MKTAY", ids="s")
            lig = Ligand("CCO", ids="l")
            Pool(runs=[seq, lig])


# ── happy path: pool two Sequence runs ────────────────────────────────────────

def test_pool_two_sequence_runs(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """Two Sequence runs with the SAME ids ('s1','s2') in each get
    composite ids ``s1_1, s2_1, s1_2, s2_2`` after pooling, exactly the
    case the user wants for parallel RFdiffusion runs."""
    from biopipelines.sequence import Sequence
    from biopipelines.pool import Pool

    pipeline = new_pipeline("pool_two_sequences")
    with pipeline:
        run_a = Sequence(seq=["MKTAY", "AETGF"], ids=["s1", "s2"], type="protein")
        run_b = Sequence(seq=["GGGGA", "ICCDE"], ids=["s1", "s2"], type="protein")
        pooled = Pool(runs=[run_a, run_b])
        pipeline.save()

    ids = list(pooled.streams.sequences.ids)
    record_case(
        input="Pool(Sequence(['s1','s2']), Sequence(['s1','s2']))",
        expected=["s1_1", "s2_1", "s1_2", "s2_2"],
        actual=ids,
    )
    assert ids == ["s1_1", "s2_1", "s1_2", "s2_2"]


def test_pool_three_sequence_runs_with_unique_ids(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """Pooling 3 runs each with unique ids: composite ids embed the pool
    index even when no original-id collision exists."""
    from biopipelines.sequence import Sequence
    from biopipelines.pool import Pool

    pipeline = new_pipeline("pool_three_sequences")
    with pipeline:
        a = Sequence(seq="MKTAY", ids="a")
        b = Sequence(seq="AETGF", ids="b")
        c = Sequence(seq="GGGGA", ids="c")
        pooled = Pool(runs=[a, b, c])
        pipeline.save()

    ids = list(pooled.streams.sequences.ids)
    record_case(
        input="Pool(Sequence('a'), Sequence('b'), Sequence('c'))",
        expected=["a_1", "b_2", "c_3"],
        actual=ids,
    )
    assert ids == ["a_1", "b_2", "c_3"]


def test_pool_no_config_time_map_and_content_table_unified(
    local_config, isolated_cwd, new_pipeline,
):
    """PDB's ``structures`` stream is content-bearing (it owns the
    ``structures`` table): its map_table is declared at the same file as the
    TableInfo, composite ids are predicted at config time, and — per the Map
    Table Contract — no CSV is written at config time."""
    from biopipelines.pdb import PDB
    from biopipelines.pool import Pool

    pipeline = new_pipeline("pool_provenance_column")
    with pipeline:
        run_a = PDB("1ubq", ids="p")
        run_b = PDB("1ubq", ids="p")
        pooled = Pool(runs=[run_a, run_b])
        pipeline.save()

    map_path = pooled.streams.structures.map_table
    assert map_path == pooled.tables.structures.info.path
    assert list(pooled.streams.structures.ids) == ["p_1", "p_2"]
    # No config-time map is written.
    assert not os.path.isfile(map_path)


def test_pool_lineage_only_stream_uses_map_csv_not_content_path(
    local_config, isolated_cwd, new_pipeline,
):
    """A stream that reuses a sibling's content CSV only as a metadata lookup
    (Sequence's ``fasta`` reusing ``sequences.csv``) is NOT content-bearing.
    Its pooled map_table must stay the lineage-only ``fasta/fasta_map.csv``,
    never a relocated/duplicated ``fasta/sequences.csv`` content table — and
    the real ``sequences`` content table must stay at ``sequences/sequences.csv``."""
    from biopipelines.sequence import Sequence
    from biopipelines.pool import Pool

    pipeline = new_pipeline("pool_lineage_only_fasta")
    with pipeline:
        a = Sequence(seq=["MK", "AY"], ids=["a", "b"], type="protein")
        b = Sequence(seq=["GF", "CD"], ids=["a", "b"], type="protein")
        pooled = Pool(runs=[a, b])
        pipeline.save()

    # The content-bearing `sequences` stream owns sequences/sequences.csv.
    seq_map = pooled.streams.sequences.map_table
    assert seq_map == pooled.tables.sequences.info.path
    assert os.path.basename(seq_map) == "sequences.csv"
    assert os.path.basename(os.path.dirname(seq_map)) == "sequences"

    # The lineage-only `fasta` stream stays at fasta/fasta_map.csv — it must
    # NOT be classified content-bearing nor relocate sequences.csv under fasta/.
    fasta_map = pooled.streams.fasta.map_table
    assert os.path.basename(fasta_map) == "fasta_map.csv"
    assert os.path.basename(os.path.dirname(fasta_map)) == "fasta"
    assert fasta_map != seq_map


# ── happy path: pool two Ligand runs ──────────────────────────────────────────

def test_pool_two_ligand_runs(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """Pool also works on Ligand outputs (compounds stream)."""
    from biopipelines.ligand import Ligand
    from biopipelines.pool import Pool

    pipeline = new_pipeline("pool_two_ligands")
    with pipeline:
        a = Ligand({"a": "CCO", "b": "C(=O)O"})
        b = Ligand({"a": "CC", "b": "CCC"})
        pooled = Pool(runs=[a, b])
        pipeline.save()

    # Ligand exposes its data via the `compounds` stream.
    assert "compounds" in pooled.streams._streams or hasattr(
        pooled.streams, "compounds"
    )
    ids = list(pooled.streams.compounds.ids)
    record_case(
        input="Pool(Ligand({'a','b'}), Ligand({'a','b'}))",
        expected=["a_1", "b_1", "a_2", "b_2"],
        actual=ids,
    )
    assert ids == ["a_1", "b_1", "a_2", "b_2"]


def test_pool_compounds_stream_map_is_the_content_table(
    local_config, isolated_cwd, new_pipeline,
):
    """The compounds stream is content-bearing: its map_table must be the
    same file as the pooled compounds TableInfo, and live at
    ``compounds/compounds.csv`` (the upstream convention), not a separate
    ``compounds_map.csv`` or ``tables/compounds.csv``."""
    from biopipelines.ligand import Ligand
    from biopipelines.pool import Pool

    pipeline = new_pipeline("pool_compounds_unify")
    with pipeline:
        a = Ligand({"a": "CCO", "b": "C(=O)O"})
        b = Ligand({"a": "CC", "b": "CCC"})
        pooled = Pool(runs=[a, b])
        pipeline.save()

    map_path = pooled.streams.compounds.map_table
    table_path = pooled.tables.compounds.info.path
    assert map_path == table_path
    assert os.path.basename(map_path) == "compounds.csv"
    assert os.path.basename(os.path.dirname(map_path)) == "compounds"


def test_pool_preserves_compounds_chemistry_at_runtime(
    local_config, isolated_cwd, new_pipeline, tmp_path,
):
    """End-to-end runtime check: after pipe_pool runs, the compounds stream
    map_table preserves the Ligand Contract domain columns (code, smiles)
    plus pool.path — so a downstream ligand-aware tool no longer hits a
    KeyError on streams.compounds."""
    import json
    import pandas as pd
    from biopipelines.ligand import Ligand
    from biopipelines.pool import Pool
    import pipe_scripts.pipe_pool as pipe_pool

    # Two upstream Ligand runs with real compounds.csv content tables.
    cols = ["id", "format", "code", "lookup", "source", "ccd", "cid",
            "cas", "smiles", "name", "formula", "file_path"]
    up = []
    for i, smis in enumerate([("CCO", "CC(=O)O"), ("CC", "CCC")], start=1):
        d = tmp_path / f"run{i}" / "compounds"
        d.mkdir(parents=True)
        csv = d / "compounds.csv"
        pd.DataFrame([
            {"id": "a", "code": "LIG", "smiles": smis[0], "format": "csv",
             "lookup": "", "source": "smiles", "ccd": "", "cid": "",
             "cas": "", "name": "", "formula": "", "file_path": ""},
            {"id": "b", "code": "LIG", "smiles": smis[1], "format": "csv",
             "lookup": "", "source": "smiles", "ccd": "", "cid": "",
             "cas": "", "name": "", "formula": "", "file_path": ""},
        ], columns=cols).to_csv(csv, index=False)
        up.append(str(csv))

    out_compounds = tmp_path / "out" / "compounds"
    out_compounds.mkdir(parents=True)
    combined = str(out_compounds / "compounds.csv")

    config = {
        "runs": [
            {"streams": [{"name": "compounds", "stream_key": "compounds",
                          "format": "csv", "ids": ["a", "b"], "files": [],
                          "map_table": up[i], "is_shared_file": False}],
             "tables": [{"name": "compounds", "path": up[i], "columns": cols}]}
            for i in range(2)
        ],
        "shared_streams": ["compounds"],
        "shared_tables": ["compounds"],
        "out_streams": {"compounds": {
            "stream_dir": str(out_compounds), "map_table": combined,
            "format": "csv"}},
        "out_tables": {"compounds": combined},
        "content_bearing_streams": ["compounds"],
        "content_file_cols": {"compounds": "file_path"},
        "output_folder": str(tmp_path / "out"),
        "recount_prefix": None,
    }
    config_json = tmp_path / "pool_config.json"
    config_json.write_text(json.dumps(config))

    import sys
    sys.argv = ["pipe_pool.py", str(config_json)]
    pipe_pool.main()

    df = pd.read_csv(combined)
    assert {"id", "code", "smiles", "pool.path"} <= set(df.columns)
    assert list(df["id"]) == ["a_1", "b_1", "a_2", "b_2"]
    assert list(df["pool.path"]) == [1, 1, 2, 2]
    assert list(df["smiles"]) == ["CCO", "CC(=O)O", "CC", "CCC"]
    # A content table must NOT carry the generic file/value scratch columns.
    assert "value" not in df.columns
    assert "file" not in df.columns


def test_pool_pdb_content_table_preserves_file_path_at_runtime(
    local_config, isolated_cwd, new_pipeline, tmp_path,
):
    """End-to-end runtime check for a content-bearing stream that ALSO carries
    per-id files (PDB's ``structures``): the pooled content table must keep the
    upstream ``file_path`` column (not a generic ``file``), and the structure
    files must be copied into the gather folder. Guards the Map Table Contract
    schema fidelity downstream readers (e.g. rcsb) rely on."""
    import json
    import pandas as pd
    import pipe_scripts.pipe_pool as pipe_pool

    cols = ["id", "pdb_id", "file_path", "format", "file_size", "source"]
    up = []
    for i in range(1, 3):
        d = tmp_path / f"run{i}" / "structures"
        d.mkdir(parents=True)
        pdb_file = d / "p.pdb"
        pdb_file.write_text("ATOM      1  N   MET A   1\n")
        csv = d / "structures.csv"
        pd.DataFrame([{
            "id": "p", "pdb_id": "1ubq", "file_path": str(pdb_file),
            "format": "pdb", "file_size": pdb_file.stat().st_size,
            "source": "rcsb",
        }], columns=cols).to_csv(csv, index=False)
        up.append((str(csv), str(pdb_file)))

    out_structs = tmp_path / "out" / "structures"
    out_structs.mkdir(parents=True)
    combined = str(out_structs / "structures.csv")

    config = {
        "runs": [
            {"streams": [{"name": "structures", "stream_key": "structures",
                          "format": "pdb", "ids": ["p"],
                          "files": [str(tmp_path / f"run{i+1}" / "structures" / "<id>.pdb")],
                          "map_table": up[i][0], "is_shared_file": False}],
             "tables": [{"name": "structures", "path": up[i][0], "columns": cols}]}
            for i in range(2)
        ],
        "shared_streams": ["structures"],
        "shared_tables": ["structures"],
        "out_streams": {"structures": {
            "stream_dir": str(out_structs), "map_table": combined,
            "format": "pdb"}},
        "out_tables": {"structures": combined},
        "content_bearing_streams": ["structures"],
        "content_file_cols": {"structures": "file_path"},
        "output_folder": str(tmp_path / "out"),
        "recount_prefix": None,
    }
    config_json = tmp_path / "pool_config.json"
    config_json.write_text(json.dumps(config))

    import sys
    sys.argv = ["pipe_pool.py", str(config_json)]
    pipe_pool.main()

    df = pd.read_csv(combined)
    assert "file_path" in df.columns
    assert "file" not in df.columns
    assert "value" not in df.columns
    assert list(df["id"]) == ["p_1", "p_2"]
    assert list(df["pool.path"]) == [1, 2]
    # The structure files were copied into the gather folder under composite ids.
    copied = [df["file_path"].iloc[0], df["file_path"].iloc[1]]
    assert all(os.path.exists(p) for p in copied)
    assert os.path.basename(copied[0]) == "p_1.pdb"
    assert os.path.basename(copied[1]) == "p_2.pdb"


# ── tables ────────────────────────────────────────────────────────────────────

def test_pool_propagates_table_schema(
    local_config, isolated_cwd, new_pipeline,
):
    """Every TableInfo present on input #1 should appear on the pooled
    output, with ``pool.path`` appended to the columns list."""
    from biopipelines.sequence import Sequence
    from biopipelines.pool import Pool

    pipeline = new_pipeline("pool_tables")
    with pipeline:
        a = Sequence(seq="MKTAY", ids="a")
        b = Sequence(seq="AETGF", ids="b")
        pooled = Pool(runs=[a, b])
        pipeline.save()

    # Sequence emits a `sequences` table (id | sequence | ...).
    assert hasattr(pooled.tables, "sequences")
    cols = list(pooled.tables.sequences.info.columns)
    assert "pool.path" in cols


# ── happy path: pool two PDB runs (file-based stream + multiple streams) ────

def test_pool_two_pdb_runs(
    local_config, isolated_cwd, new_pipeline,
):
    """PDB emits structures + sequences + compounds streams. Pool over
    two PDB runs must:

    * preserve all three streams,
    * track each via its native map_table column convention (PDB uses
      `file_path` for structures.csv, not `file`),
    * compute the *real* file extension at config time so the
      generated DataStream files don't carry a literal '.*' wildcard.
    """
    from biopipelines.pdb import PDB
    from biopipelines.pool import Pool

    pipeline = new_pipeline("pool_two_pdbs")
    with pipeline:
        run_a = PDB("1ubq", ids="p")
        run_b = PDB("1ubq", ids="p")
        pooled = Pool(runs=[run_a, run_b])
        pipeline.save()

    # Three streams should all be present.
    for name in ("structures", "sequences", "compounds"):
        assert hasattr(pooled.streams, name), f"missing stream {name}"

    # Composite ids on structures.
    s_ids = list(pooled.streams.structures.ids)
    assert s_ids == ["p_1", "p_2"]

    # The DataStream's files entry must NOT carry a literal '.*' that
    # would propagate into the per-id expansion. Either a concrete
    # extension is resolved at config time, or the framework's <id>
    # template form is used.
    s_files = list(pooled.streams.structures.files)
    for f in s_files:
        assert not f.endswith("p_1.*"), f"wildcard ext leaked into files: {f}"
        assert not f.endswith("p_2.*"), f"wildcard ext leaked into files: {f}"


# ── shape sanity: composite ids match the renumbering rule ───────────────────

def test_pool_id_format_is_orig_underscore_pool_idx(
    local_config, isolated_cwd, new_pipeline,
):
    """Verify the user-promised contract: composite ids are exactly
    ``<orig>_<pool_idx>`` with 1-based pool_idx, in input order."""
    from biopipelines.sequence import Sequence
    from biopipelines.pool import Pool

    pipeline = new_pipeline("pool_id_shape")
    with pipeline:
        runs = [
            Sequence(seq="MK", ids="design", type="protein"),
            Sequence(seq="AY", ids="design", type="protein"),
            Sequence(seq="GF", ids="design", type="protein"),
        ]
        pooled = Pool(runs=runs)
        pipeline.save()

    assert list(pooled.streams.sequences.ids) == ["design_1", "design_2", "design_3"]
