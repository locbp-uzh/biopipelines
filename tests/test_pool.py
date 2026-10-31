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


@pytest.mark.network
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


def test_pool_pdb_content_table_preserves_the_file_column_at_runtime(
    local_config, isolated_cwd, new_pipeline, tmp_path,
):
    """End-to-end runtime check for a content-bearing stream that ALSO carries
    per-id files (PDB's ``structures``): the pooled content table must keep the
    canonical ``file`` column, and the structure files must be copied into the
    gather folder. Guards the schema fidelity downstream readers rely on."""
    import json
    import pandas as pd
    import pipe_scripts.pipe_pool as pipe_pool

    cols = ["id", "pdb_id", "file", "format", "file_size", "source"]
    up = []
    for i in range(1, 3):
        d = tmp_path / f"run{i}" / "structures"
        d.mkdir(parents=True)
        pdb_file = d / "p.pdb"
        pdb_file.write_text("ATOM      1  N   MET A   1\n")
        csv = d / "structures.csv"
        pd.DataFrame([{
            "id": "p", "pdb_id": "1ubq", "file": str(pdb_file),
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
        "content_file_cols": {"structures": "file"},
        "output_folder": str(tmp_path / "out"),
        "recount_prefix": None,
    }
    config_json = tmp_path / "pool_config.json"
    config_json.write_text(json.dumps(config))

    import sys
    sys.argv = ["pipe_pool.py", str(config_json)]
    pipe_pool.main()

    df = pd.read_csv(combined)
    assert "file" in df.columns
    assert "value" not in df.columns
    assert list(df["id"]) == ["p_1", "p_2"]
    assert list(df["pool.path"]) == [1, 2]
    # The structure files were copied into the gather folder under composite ids.
    copied = [df["file"].iloc[0], df["file"].iloc[1]]
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

@pytest.mark.network
def test_pool_two_pdb_runs(
    local_config, isolated_cwd, new_pipeline,
):
    """PDB emits structures + sequences + compounds streams. Pool over
    two PDB runs must:

    * preserve all three streams,
    * track each through the canonical `file` column of its map_table,
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


# ── original.id and keep_ids ─────────────────────────────────────────────────
#
# The default suffix renames EVERY id. A step that maps a design to its parent
# with a raw ``sid.startswith(parent_id)`` then matches by accident: the pooled
# parent ``P_1`` is a prefix of the pooled design ``P_12_1`` but not of ``P_7_1``.

N_DESIGNS = 100
GROUPS = (1, 2, 3)


def _write_run_tables(tmp_path, name, ids_by_run):
    """One CSV per run, the shape a pooled upstream table has before Pool."""
    import pandas as pd

    runs = []
    for i, ids in enumerate(ids_by_run, start=1):
        d = tmp_path / f"{name}_run{i}"
        d.mkdir(parents=True)
        csv = d / f"{name}.csv"
        pd.DataFrame({"id": ids, "value": range(len(ids))}).to_csv(csv, index=False)
        runs.append([{"name": name, "path": str(csv), "columns": ["id", "value"]}])
    return runs


def _pool_table(tmp_path, name, ids_by_run, **options):
    """Pool the tables through pipe_pool's real table logic; return the result."""
    import pandas as pd
    import pipe_scripts.pipe_pool as pipe_pool

    out = tmp_path / f"{name}_pooled.csv"
    pipe_pool._gather_table(name, _write_run_tables(tmp_path, name, ids_by_run), str(out), **options)
    return pd.read_csv(out)


def test_prefix_matching_loses_pooled_designs_and_original_id_recovers_them(tmp_path):
    """Three parents pooled as P_1..P_3, 100 designs per group pooled as P_<n>_<group>.

    A prefix lookup finds only the designs whose number starts with the group's
    digit: 12, 11 and 11 of 100 (design 1, 10-19 and 100; 2 and 20-29; 3 and
    30-39). Looking up through ``original.id`` finds all of them.
    """
    parents = _pool_table(tmp_path, "parents", [["P"] for _ in GROUPS])
    designs = _pool_table(tmp_path, "designs",
                          [[f"P_{n}" for n in range(1, N_DESIGNS + 1)] for _ in GROUPS])

    assert list(parents["id"]) == ["P_1", "P_2", "P_3"]
    assert designs["id"].iloc[11] == "P_12_1"

    by_prefix, by_original = [], []
    for g in GROUPS:
        parent = parents[parents["pool.path"] == g].iloc[0]
        group = designs[designs["pool.path"] == g]
        by_prefix.append(int(group["id"].str.startswith(parent["id"]).sum()))
        # the design's own id minus its design number is the parent's pre-pool id
        by_original.append(int((group["original.id"].str.rsplit("_", n=1).str[0]
                                == parent["original.id"]).sum()))

    assert by_prefix == [12, 11, 11]
    assert by_original == [N_DESIGNS] * len(GROUPS)


def test_keep_ids_never_adds_a_group_suffix(tmp_path):
    """With distinct ids per run, keep_ids leaves every id exactly as written."""
    ids_by_run = [[f"P{g}_{n}" for n in range(1, N_DESIGNS + 1)] for g in GROUPS]
    designs = _pool_table(tmp_path, "designs", ids_by_run, keep_ids=True)

    assert list(designs["id"]) == [i for ids in ids_by_run for i in ids]
    assert list(designs["original.id"]) == list(designs["id"])
    assert list(designs["pool.path"][::N_DESIGNS]) == list(GROUPS)
    # the raw prefix lookup that failed for the suffixed ids now matches every design
    for g in GROUPS:
        group = designs[designs["pool.path"] == g]
        assert group["id"].str.startswith(f"P{g}").all()


def test_keep_ids_collision_in_a_table_names_the_ids(tmp_path):
    with pytest.raises(ValueError, match=r"table 'designs'.*dup_a \(runs 1, 2\).*dup_b \(runs 1, 2\)"):
        _pool_table(tmp_path, "designs", [["dup_a", "dup_b", "only_1"], ["dup_a", "dup_b"]],
                    keep_ids=True)


def test_keep_ids_collision_in_a_stream_map_is_raised_before_any_copy(tmp_path):
    import pipe_scripts.pipe_pool as pipe_pool

    runs = []
    for i in (1, 2):
        run_dir = tmp_path / f"run{i}"
        run_dir.mkdir()
        runs.append([{"name": "structures", "stream_key": "structures", "format": "pdb",
                      "ids": ["same"], "files": [], "map_table": "", "is_shared_file": False}])
    out_dir = tmp_path / "out"
    spec = {"stream_dir": str(out_dir), "map_table": str(out_dir / "structures_map.csv"),
            "format": "pdb"}
    with pytest.raises(ValueError, match=r"stream 'structures'.*same \(runs 1, 2\)"):
        pipe_pool._gather_stream("structures", runs, spec, keep_ids=True)
    assert not (out_dir / "structures_map.csv").exists()


def test_same_ids_are_fine_without_keep_ids(tmp_path):
    pooled = _pool_table(tmp_path, "designs", [["same"], ["same"]])
    assert list(pooled["id"]) == ["same_1", "same_2"]
    assert list(pooled["original.id"]) == ["same", "same"]


def test_original_id_is_emitted_on_every_stream_map_row_in_default_mode(tmp_path):
    import pandas as pd
    import pipe_scripts.pipe_pool as pipe_pool

    runs = []
    for i in (1, 2):
        run_dir = tmp_path / f"run{i}"
        run_dir.mkdir()
        pdb = run_dir / "p.pdb"
        pdb.write_text("ATOM\n")
        runs.append([{"name": "structures", "stream_key": "structures", "format": "pdb",
                      "ids": ["p"], "files": [str(pdb)], "map_table": "", "is_shared_file": False}])
    out_dir = tmp_path / "out"
    spec = {"stream_dir": str(out_dir), "map_table": str(out_dir / "structures_map.csv"),
            "format": "pdb"}
    pipe_pool._gather_stream("structures", runs, spec)

    df = pd.read_csv(out_dir / "structures_map.csv")
    assert list(df["id"]) == ["p_1", "p_2"]
    assert list(df["original.id"]) == ["p", "p"]
    assert list(df["pool.path"]) == [1, 2]


def test_renumber_ids_default_mode_keeps_an_upstream_original_id_and_suffixes_axis_ids():
    import pandas as pd
    import pipe_scripts.pipe_pool as pipe_pool

    plain = pd.DataFrame({"id": ["a"], "sequences.id": ["s"]})
    out = pipe_pool._renumber_ids(plain, 2)
    assert out.iloc[0].to_dict() == {"id": "a_2", "sequences.id": "s_2", "original.id": "a"}

    carried = pd.DataFrame({"id": ["m"], "original.id": ["parent"]})
    out = pipe_pool._renumber_ids(carried, 2)
    assert out.iloc[0].to_dict() == {"id": "m_2", "original.id": "parent"}


def test_a_stream_map_keeps_an_upstream_original_id_and_suffixes_axis_ids(tmp_path):
    """The stream path carries the column through the same rule as the tables."""
    import pandas as pd
    import pipe_scripts.pipe_pool as pipe_pool

    runs = []
    for i in (1, 2):
        run_dir = tmp_path / f"run{i}"
        run_dir.mkdir()
        pdb = run_dir / "m.pdb"
        pdb.write_text("ATOM\n")
        map_csv = run_dir / "structures_map.csv"
        pd.DataFrame({"id": ["m"], "file": [str(pdb)], "original.id": ["parent"],
                      "sequences.id": ["s"]}).to_csv(map_csv, index=False)
        runs.append([{"name": "structures", "stream_key": "structures", "format": "pdb",
                      "ids": ["m"], "files": [str(pdb)], "map_table": str(map_csv),
                      "is_shared_file": False}])
    out_dir = tmp_path / "out"
    spec = {"stream_dir": str(out_dir), "map_table": str(out_dir / "structures_map.csv"),
            "format": "pdb"}
    pipe_pool._gather_stream("structures", runs, spec)

    df = pd.read_csv(out_dir / "structures_map.csv")
    assert list(df["id"]) == ["m_1", "m_2"]
    assert list(df["original.id"]) == ["parent", "parent"]
    assert list(df["sequences.id"]) == ["s_1", "s_2"]


def test_renumber_ids_recount_mode_is_unchanged():
    import pandas as pd
    import pipe_scripts.pipe_pool as pipe_pool

    df = pd.DataFrame({"id": ["a", "b"], "sequences.id": ["s", "t"]})
    numbering = {(1, "x"): "r_1", (1, "y"): "r_2", (1, "z"): "r_3", (1, "w"): "r_4"}
    out = pipe_pool._renumber_ids(df, 2, recount_prefix="r", numbering=numbering)
    assert list(out["id"]) == ["r_5", "r_6"]
    assert list(out["original.id"]) == ["a", "b"]
    assert list(out["sequences.id"]) == ["s", "t"]


def test_pool_declares_the_original_id_column(local_config, isolated_cwd, new_pipeline):
    from biopipelines.sequence import Sequence
    from biopipelines.pool import Pool

    pipeline = new_pipeline("pool_original_id_column")
    with pipeline:
        pooled = Pool(runs=[Sequence(seq="MK", ids="a", type="protein"),
                            Sequence(seq="AY", ids="b", type="protein")])
        pipeline.save()

    cols = list(pooled.tables.sequences.info.columns)
    assert "original.id" in cols and "pool.path" in cols


def test_keep_ids_at_configuration_time_leaves_ids_untouched(local_config, isolated_cwd, new_pipeline):
    from biopipelines.sequence import Sequence
    from biopipelines.pool import Pool

    pipeline = new_pipeline("pool_keep_ids")
    with pipeline:
        pooled = Pool(runs=[Sequence(seq="MK", ids="a", type="protein"),
                            Sequence(seq="AY", ids="b", type="protein")],
                      keep_ids=True)
        pipeline.save()

    assert list(pooled.streams.sequences.ids) == ["a", "b"]


def test_keep_ids_collision_at_configuration_time_names_the_ids(local_config, isolated_cwd, new_pipeline):
    from biopipelines.sequence import Sequence
    from biopipelines.pool import Pool

    pipeline = new_pipeline("pool_keep_ids_collision")
    with pytest.raises(ValueError, match=r"keep_ids.*s1 \(runs 1, 2\).*s2 \(runs 1, 2\)"):
        with pipeline:
            Pool(runs=[Sequence(seq=["MK", "AY"], ids=["s1", "s2"], type="protein"),
                       Sequence(seq=["GF", "CD"], ids=["s1", "s2"], type="protein")],
                 keep_ids=True)


def test_keep_ids_and_recount_prefix_are_mutually_exclusive(local_config, isolated_cwd, new_pipeline):
    from biopipelines.sequence import Sequence
    from biopipelines.pool import Pool

    pipeline = new_pipeline("pool_keep_ids_recount")
    with pytest.raises(ValueError, match="mutually exclusive"):
        with pipeline:
            Pool(runs=[Sequence(seq="MK", ids="a", type="protein"),
                       Sequence(seq="AY", ids="b", type="protein")],
                 keep_ids=True, recount_prefix="r")


def _narrowed_batch(root, name, kept):
    """A finished ProteinMPNN-shaped step cut by hand to `kept` of its 10 designs, its manifest still naming all 10 and the source folder."""
    import csv
    import json

    folder = root / name / "003_ProteinMPNN"
    streams = {}
    for stream in ("sequences", "designs"):
        table = folder / stream / f"{stream}.csv"
        table.parent.mkdir(parents=True)
        with open(table, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["id", "sequence"])
            w.writerows([design, "MKTAY"] for design in kept)
        streams[stream] = {"name": stream, "ids": ["design_<1..10>"], "files": [],
                           "map_table": str(table), "format": "csv", "metadata": {}}
    (folder / ".expected_outputs.json").write_text(json.dumps({
        "tool_name": "ProteinMPNN", "tool_class": "ProteinMPNN",
        "output_structure": {**streams, "tables": {},
                             "output_folder": str(root / "source" / "003_ProteinMPNN")},
    }))
    return str(folder)


def test_pool_over_hand_narrowed_batches_declares_the_sum_of_the_batches(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """Two batches of one finished step pool to the designs they hold, not to runs x the producer's full list."""
    from biopipelines.load import Load
    from biopipelines.pool import Pool

    designs = [f"design_{i}" for i in range(1, 11)]
    batch_a = _narrowed_batch(isolated_cwd, "a", designs[:4])
    batch_b = _narrowed_batch(isolated_cwd, "b", designs[4:7])

    pipeline = new_pipeline("pool_narrowed_batches")
    with pipeline:
        pooled = Pool(runs=[Load(batch_a), Load(batch_b)])
        pipeline.save()

    expected = [f"design_{i}_1" for i in range(1, 5)] + [f"design_{i}_2" for i in range(5, 8)]
    actual = {name: list(getattr(pooled.streams, name).ids) for name in ("sequences", "designs")}
    record_case(input="Pool(Load(4 of 10), Load(3 of 10))",
                expected={"sequences": expected, "designs": expected}, actual=actual)
    assert actual == {"sequences": expected, "designs": expected}


def test_keep_ids_runs_end_to_end_through_pipe_pool_main(tmp_path):
    """The config Pool writes carries keep_ids; main() honors it for streams and tables."""
    import json
    import sys
    import pandas as pd
    import pipe_scripts.pipe_pool as pipe_pool

    cols = ["id", "code", "smiles"]
    up = []
    for i, ids in enumerate((["a", "b"], ["c"]), start=1):
        d = tmp_path / f"run{i}"
        d.mkdir()
        csv = d / "compounds.csv"
        pd.DataFrame({"id": ids, "code": "LIG", "smiles": "C"}, columns=cols).to_csv(csv, index=False)
        up.append((ids, str(csv)))
    out = tmp_path / "out" / "compounds"
    out.mkdir(parents=True)
    combined = str(out / "compounds.csv")
    config = {
        "runs": [{"streams": [{"name": "compounds", "stream_key": "compounds", "format": "csv",
                               "ids": ids, "files": [], "map_table": csv, "is_shared_file": False}],
                  "tables": [{"name": "compounds", "path": csv, "columns": cols}]}
                 for ids, csv in up],
        "shared_streams": ["compounds"], "shared_tables": ["compounds"],
        "out_streams": {"compounds": {"stream_dir": str(out), "map_table": combined, "format": "csv"}},
        "out_tables": {"compounds": combined},
        "content_bearing_streams": ["compounds"], "content_file_cols": {"compounds": "file_path"},
        "output_folder": str(tmp_path / "out"), "recount_prefix": None, "keep_ids": True,
    }
    config_json = tmp_path / "pool_config.json"
    config_json.write_text(json.dumps(config))

    sys.argv = ["pipe_pool.py", str(config_json)]
    pipe_pool.main()

    df = pd.read_csv(combined)
    assert list(df["id"]) == ["a", "b", "c"]
    assert list(df["original.id"]) == ["a", "b", "c"]
    assert list(df["pool.path"]) == [1, 1, 2]


def test_pool_over_rows_sharing_one_file_names_every_id_without_copying_per_row(tmp_path):
    """Several ids of one run point at the same file; each pooled id still gets its own name, and the bytes are not repeated.

    Panda used to write one copy per id for this shape and miss the file the stream declared; Pool reads each
    row's own `file` and links it, so it was never affected. This pins that.

    The extra `topology` column of the upstream rows is carried too, still naming the upstream file: nothing is lost,
    and Pool links its stream files to their upstream location the same way.
    """
    import json
    import sys
    import pandas as pd
    import pipe_scripts.pipe_pool as pipe_pool

    cols = ["id", "file", "topology"]
    maps = []
    for i in (1, 2):
        d = tmp_path / f"run{i}" / "trajectories"
        d.mkdir(parents=True)
        shared, own = d / "shared.xtc", d / "own.xtc"
        shared.write_bytes(b"shared")
        own.write_bytes(b"own")
        topology = d / "top.cif"
        topology.write_text(f"data_run{i}")
        csv = d / "trajectories_map.csv"
        pd.DataFrame([
            {"id": "a", "file": str(shared), "topology": str(topology)},
            {"id": "b", "file": str(shared), "topology": str(topology)},
            {"id": "c", "file": str(own), "topology": str(topology)},
        ], columns=cols).to_csv(csv, index=False)
        maps.append(str(csv))

    out = tmp_path / "out" / "trajectories"
    out.mkdir(parents=True)
    combined = str(out / "trajectories_map.csv")
    config = {
        "runs": [{"streams": [{"name": "trajectories", "stream_key": "trajectories", "format": "xtc",
                               "ids": ["a", "b", "c"],
                               "files": [str(tmp_path / f"run{i + 1}" / "trajectories" / "<id>.xtc")],
                               "map_table": maps[i], "is_shared_file": False}],
                  "tables": []} for i in range(2)],
        "shared_streams": ["trajectories"], "shared_tables": [],
        "out_streams": {"trajectories": {"stream_dir": str(out), "map_table": combined, "format": "xtc"}},
        "out_tables": {}, "content_bearing_streams": [], "content_file_cols": {},
        "output_folder": str(tmp_path / "out"), "recount_prefix": None,
    }
    config_json = tmp_path / "pool_config.json"
    config_json.write_text(json.dumps(config))

    sys.argv = ["pipe_pool.py", str(config_json)]
    pipe_pool.main()

    df = pd.read_csv(combined)
    assert list(df["id"]) == ["a_1", "b_1", "c_1", "a_2", "b_2", "c_2"]
    expected = [b"shared", b"shared", b"own"] * 2
    for path, content in zip(df["file"], expected):
        assert os.path.exists(path), path
        with open(path, "rb") as handle:
            assert handle.read() == content
    for path, run in zip(df["topology"], [1, 1, 1, 2, 2, 2]):
        assert os.path.exists(path), f"the carried topology column names a file that is not there: {path}"
        with open(path) as handle:
            assert handle.read() == f"data_run{run}", "each row keeps its own run's topology"

    probe = tmp_path / "probe"
    try:
        os.symlink(str(tmp_path / "run1"), str(probe))
    except (OSError, NotImplementedError):
        return
    assert len({(os.stat(p).st_dev, os.stat(p).st_ino) for p in df["file"]}) == 4, (
        "two runs of a shared and an own file are four files, not six")


def test_renumber_ids_leaves_a_blank_provenance_cell_blank():
    """A failed id keeps a blank row in a result table; its empty `<axis>.id` cell stays empty instead of becoming `nan_2`.

    Object dtype is pinned because that is what pandas 2 reads a column with a blank cell as, and what `astype(str)` turns into the text `nan`.
    """
    import pandas as pd
    import pipe_scripts.pipe_pool as pipe_pool

    with pd.option_context("future.infer_string", False):
        frame = pd.DataFrame({"id": ["a", "b"], "sequences.id": ["s", None]})
        out = pipe_pool._renumber_ids(frame, 2)

    assert list(out["id"]) == ["a_2", "b_2"]
    assert out["sequences.id"].iloc[0] == "s_2"
    assert pd.isna(out["sequences.id"].iloc[1])


# ── recount_prefix: one row has one number in every stream and table ─────────
#
# Under recount_prefix each stream and table used to number its own rows 1..N. An id the stream lacks at runtime
# (an upstream step filtered it) still has a blank row in the result table and a row in `missing`, so after that id
# the same row carried a different number in each file and a join on `id` paired unrelated rows.

def _write_recount_runs(tmp_path):
    """Two runs whose stream `kept` lacks the id `b` while `scores` and `missing` still name it."""
    import pandas as pd

    runs = []
    for i, (declared, present) in enumerate((
        (["a", "b", "c"], ["a", "c"]),
        (["a", "d"], ["a", "d"]),
    ), start=1):
        folder = tmp_path / f"run{i}"
        folder.mkdir()
        kept = folder / "kept_map.csv"
        pd.DataFrame({"id": present}).to_csv(kept, index=False)
        scores = folder / "scores.csv"
        pd.DataFrame({"id": declared, "score": range(len(declared))}).to_csv(scores, index=False)
        missing = folder / "missing.csv"
        pd.DataFrame({"id": [x for x in declared if x not in present], "removed_by": "T", "kind": "filter",
                      "cause": "c"}, columns=["id", "removed_by", "kind", "cause"]).to_csv(missing, index=False)
        runs.append({
            "streams": [{"name": "kept", "stream_key": "kept", "format": "csv", "ids": declared, "files": [],
                         "map_table": str(kept), "is_shared_file": False}],
            "tables": [{"name": name, "path": str(path), "columns": []}
                       for name, path in (("missing", missing), ("scores", scores))],
        })
    return runs


def test_recount_names_a_row_alike_in_the_stream_the_table_and_missing(tmp_path):
    import json
    import sys
    import pandas as pd
    import pipe_scripts.pipe_pool as pipe_pool

    out = tmp_path / "out"
    (out / "kept").mkdir(parents=True)
    config = {
        "runs": _write_recount_runs(tmp_path),
        "shared_streams": ["kept"], "shared_tables": ["missing", "scores"],
        "out_streams": {"kept": {"stream_dir": str(out / "kept"), "map_table": str(out / "kept_map.csv"), "format": "csv"}},
        "out_tables": {"missing": str(out / "missing.csv"), "scores": str(out / "scores.csv")},
        "content_bearing_streams": [], "content_file_cols": {}, "output_folder": str(out),
        "recount_prefix": "rc", "keep_ids": False,
        "recount_numbering": [[1, "a", "rc_1"], [1, "b", "rc_2"], [1, "c", "rc_3"], [2, "a", "rc_4"], [2, "d", "rc_5"]],
    }
    config_json = tmp_path / "pool_config.json"
    config_json.write_text(json.dumps(config))

    sys.argv = ["pipe_pool.py", str(config_json)]
    pipe_pool.main()

    kept, scores, missing = (pd.read_csv(out / name) for name in ("kept_map.csv", "scores.csv", "missing.csv"))
    names = {}
    for frame in (kept, scores, missing):
        for _, row in frame.iterrows():
            names.setdefault((row["pool.path"], row["original.id"]), set()).add(row["id"])
    assert all(len(ids) == 1 for ids in names.values()), f"a row is named differently in different files: {names}"
    assert list(kept["id"]) == ["rc_1", "rc_3", "rc_4", "rc_5"]
    assert list(scores["id"]) == ["rc_1", "rc_2", "rc_3", "rc_4", "rc_5"]
    assert list(missing["id"]) == ["rc_2"]


def test_recount_declares_the_ids_the_pipe_is_given(local_config, isolated_cwd, new_pipeline):
    """The numbering written to pool_config.json is the one the declared ids come from, so the two cannot drift."""
    import json
    from biopipelines.sequence import Sequence
    from biopipelines.pool import Pool

    pipeline = new_pipeline("pool_recount_numbering")
    with pipeline:
        pooled = Pool(runs=[Sequence(seq=["MK", "AY"], ids=["a", "b"], type="protein"),
                            Sequence(seq=["GF"], ids=["a"], type="protein")], recount_prefix="rc")
        pipeline.save()

    with open(pooled.producer.pool_config_json) as handle:
        numbering = json.load(handle)["recount_numbering"]
    assert numbering == [[1, "a", "rc_1"], [1, "b", "rc_2"], [2, "a", "rc_3"]]
    assert list(pooled.streams.sequences.ids) == ["rc_1", "rc_2", "rc_3"]
