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


def test_pool_emits_pool_path_column_in_map_table(
    local_config, isolated_cwd, new_pipeline,
):
    """The map_table CSV emitted at config time must carry a ``pool.path``
    column with the 1-based source-run index for every output id."""
    import pandas as pd
    from biopipelines.sequence import Sequence
    from biopipelines.pool import Pool

    pipeline = new_pipeline("pool_provenance_column")
    with pipeline:
        run_a = Sequence(seq=["MK", "AY"], ids=["a", "b"], type="protein")
        run_b = Sequence(seq=["GF", "CD"], ids=["a", "b"], type="protein")
        pooled = Pool(runs=[run_a, run_b])
        pipeline.save()

    map_path = pooled.streams.sequences.map_table
    assert os.path.isfile(map_path), f"map_table not written: {map_path}"
    df = pd.read_csv(map_path)
    assert "pool.path" in df.columns
    assert list(df["pool.path"]) == [1, 1, 2, 2]
    assert list(df["id"]) == ["a_1", "b_1", "a_2", "b_2"]


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
