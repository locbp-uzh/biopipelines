"""Smoke + runtime tests for the Panda tool in pool / no-pool modes.

Two layers:

- Config-time: ``Pipeline(...).save()`` emits a valid bash script when Panda
  is wired downstream of Sequence/Mock/etc. (tool registration, stream/table
  wiring, script emission).
- Runtime: actually execute ``pipe_sequence.py`` / ``pipe_ligand.py`` and
  ``pipe_panda.py`` end-to-end and inspect the resulting output folder —
  every pool-forwarded content table (``sequences.csv``, ``compounds.csv``)
  must retain its upstream content columns, and every stream map-table
  written alongside must have the expected ``id, file, <stream>.id``
  provenance schema. The runtime layer exists to catch regressions like the
  stream_map_target vs. pool_table_map overwrite bug where
  ``<panda>/sequences.csv`` was clobbered with a lineage-only schema and the
  ``sequence`` content column disappeared silently.
"""

import csv
import json
import os
import subprocess
import sys

import pytest


# ── runtime helpers ──────────────────────────────────────────────────────────

def _repo_root():
    import biopipelines.panda as _p
    return os.path.dirname(os.path.dirname(os.path.abspath(_p.__file__)))


def _run_pipe(name, config_path, config_flag=True):
    """Execute a pipe_*.py script with its JSON config.

    Most pipe scripts accept ``--config <path>``; a couple take the config
    path as a positional argument (``pipe_mock.py``). ``config_flag`` picks
    the calling convention.
    """
    script = os.path.join(_repo_root(), "pipe_scripts", f"pipe_{name}.py")
    args = [sys.executable, script]
    args += ["--config", str(config_path)] if config_flag else [str(config_path)]
    subprocess.run(args, check=True)


def _read_csv_rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def _read_csv_columns(path):
    with open(path, newline="") as f:
        return next(csv.reader(f))


# ── Panda: no-pool modes ──────────────────────────────────────────────────────

def test_panda_nopool_single_table_filter_sort(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda(tables=seq.tables.sequences, operations=[filter, sort]) — single table."""
    from biopipelines.sequence import Sequence
    from biopipelines.panda import Panda

    pipeline = new_pipeline("panda_nopool_single")
    with pipeline:
        seq = Sequence(seq=["MKTAY", "AETGF", "GGGGA"],
                       type="protein", ids=["a", "b", "c"])
        Panda(
            tables=seq.tables.sequences,
            operations=[
                Panda.filter("length >= 4"),
                Panda.sort("id", ascending=True),
            ],
        )
        script_path = pipeline.save()

    record_case(input="Panda nopool single-table",
                expected="script with Panda marker",
                actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Panda", "panda_nopool_single")


def test_panda_nopool_multi_table_merge(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda(tables=[a, b], operations=[merge]) — multi-table merge."""
    from biopipelines.sequence import Sequence
    from biopipelines.panda import Panda

    pipeline = new_pipeline("panda_nopool_merge")
    with pipeline:
        apo = Sequence(seq=["MKTAY", "AETGF"], type="protein", ids=["x", "y"])
        holo = Sequence(seq=["PPPPL", "QQQQL"], type="protein", ids=["x", "y"])
        Panda(
            tables=[apo.tables.sequences, holo.tables.sequences],
            operations=[Panda.merge(prefixes=["apo_", "holo_"])],
        )
        script_path = pipeline.save()

    record_case(input="Panda nopool merge 2 tables",
                expected="merge in script",
                actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Panda", "panda_nopool_merge")


def test_panda_nopool_multi_table_concat(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda(tables=[a, b, c], operations=[concat]) — multi-table concat."""
    from biopipelines.sequence import Sequence
    from biopipelines.panda import Panda

    pipeline = new_pipeline("panda_nopool_concat")
    with pipeline:
        s1 = Sequence(seq=["MKTAY"], type="protein", ids=["a"])
        s2 = Sequence(seq=["AETGF"], type="protein", ids=["b"])
        s3 = Sequence(seq=["GGGGA"], type="protein", ids=["c"])
        Panda(
            tables=[s1.tables.sequences, s2.tables.sequences, s3.tables.sequences],
            operations=[Panda.concat()],
        )
        script_path = pipeline.save()

    record_case(input="Panda nopool concat 3 tables",
                expected="concat in script",
                actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Panda", "panda_nopool_concat")


def test_panda_merge_requires_multi_table(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """Panda.merge on a single-table input raises ValueError."""
    from biopipelines.sequence import Sequence
    from biopipelines.panda import Panda

    pipeline = new_pipeline("panda_merge_invalid")
    record_case(input="Panda.merge with single table",
                expected="ValueError", actual="ValueError")
    with pipeline:
        seq = Sequence(seq="MKTAY", type="protein", ids="a")
        with pytest.raises(ValueError):
            Panda(tables=seq.tables.sequences, operations=[Panda.merge()])


# ── Panda: pool modes ─────────────────────────────────────────────────────────

def test_panda_single_pool(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda(pool=seq_out, ...) — single-pool filter keeps matching stream files."""
    from biopipelines.sequence import Sequence
    from biopipelines.panda import Panda

    pipeline = new_pipeline("panda_pool_single")
    with pipeline:
        seq = Sequence(seq=["MKTAY", "AETGF", "GGGGA"],
                       type="protein", ids=["a", "b", "c"])
        Panda(
            tables=seq.tables.sequences,
            operations=[Panda.filter("length >= 4")],
            pool=seq,
        )
        script_path = pipeline.save()

    record_case(input="Panda single-pool",
                expected="POOL MODE wired",
                actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Panda", "panda_pool_single")


def test_panda_multi_pool_concat(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda(pool=[p1, p2], operations=[concat]) — multi-pool concat + head."""
    from biopipelines.sequence import Sequence
    from biopipelines.panda import Panda

    pipeline = new_pipeline("panda_pool_multi")
    with pipeline:
        s1 = Sequence(seq=["MKTAY", "AETGF"], type="protein", ids=["a", "b"])
        s2 = Sequence(seq=["PPPPL", "QQQQL"], type="protein", ids=["c", "d"])
        Panda(
            tables=[s1.tables.sequences, s2.tables.sequences],
            operations=[
                Panda.concat(),
                Panda.sort("id"),
                Panda.head(2),
            ],
            pool=[s1, s2],
        )
        script_path = pipeline.save()

    record_case(input="Panda multi-pool concat+head",
                expected="multi-pool script",
                actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Panda", "panda_pool_multi")


# ── Panda driven by Mock (no Sequence shared-fasta quirk) ─────────────────────

def test_panda_nopool_via_mock(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda downstream of Mock's real per-ID table."""
    from biopipelines.mock import Mock
    from biopipelines.panda import Panda

    pipeline = new_pipeline("panda_mock_nopool")
    with pipeline:
        m = Mock(
            ids=["a", "b", "c"],
            streams={"s": {"format": "pdb", "file": "<id>.pdb"}},
            tables={"scores": {"columns": ["score"], "fill": {"score": 1.0}}},
        )
        Panda(
            tables=m.tables.scores,
            operations=[Panda.sort("id", ascending=True)],
        )
        script_path = pipeline.save()

    record_case(input="Panda nopool ← Mock",
                expected="script with Panda+Mock",
                actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Panda", "Mock", "panda_mock_nopool")


def test_panda_pool_via_mock(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda pool=[Mock1.stream, Mock2.stream] — two real per-ID streams."""
    from biopipelines.mock import Mock
    from biopipelines.panda import Panda

    pipeline = new_pipeline("panda_mock_pool")
    with pipeline:
        m1 = Mock(ids=["x1", "x2"],
                  streams={"s": {"format": "pdb", "file": "<id>.pdb"}},
                  tables={"scores": {"columns": ["score"], "fill": 1.0}})
        m2 = Mock(ids=["y1", "y2"],
                  streams={"s": {"format": "pdb", "file": "<id>.pdb"}},
                  tables={"scores": {"columns": ["score"], "fill": 2.0}})
        Panda(
            tables=[m1.tables.scores, m2.tables.scores],
            operations=[Panda.concat(), Panda.sort("id")],
            pool=[m1, m2],
        )
        script_path = pipeline.save()

    record_case(input="Panda pool ← 2 Mocks",
                expected="pool script with Mock",
                actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Panda", "Mock", "panda_mock_pool")


# ── Panda op-cycle tests: Mock → Panda(<op>) → Mock ──────────────────────────

def test_panda_cycle_filter_preserves_ids(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda(filter, pool) keeps the upstream IDs (no rename)."""
    from biopipelines.mock import Mock
    from biopipelines.panda import Panda

    pipeline = new_pipeline("panda_cycle_filter")
    with pipeline:
        m = Mock(
            ids=["p1", "p2", "p3"],
            streams={"s": {"format": "pdb", "file": "<id>.pdb"}},
            tables={"scores": {"columns": ["score"], "fill": {"score": 0.7}}},
        )
        p = Panda(
            tables=m.tables.scores,
            operations=[Panda.filter("score > 0.5")],
            pool=m,
        )
        Mock(source=p.streams.s, children="<1..2>",
             streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
        script_path = pipeline.save()

    ids = list(p.streams.s.ids)
    record_case(input="filter-only cycle", expected=["p1", "p2", "p3"], actual=ids)
    assert ids == ["p1", "p2", "p3"]
    assert_valid_script(script_path, "Mock", "Panda", "panda_cycle_filter")


def test_panda_cycle_sort_head_renames(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda(sort+head, pool) renames to fixed count; cycle is script-valid."""
    from biopipelines.mock import Mock
    from biopipelines.panda import Panda

    pipeline = new_pipeline("panda_cycle_sort_head")
    with pipeline:
        m = Mock(
            ids=["r1", "r2", "r3", "r4"],
            streams={"s": {"format": "pdb", "file": "<id>.pdb"}},
            tables={"scores": {"columns": ["score"], "fill": {"score": 0.8}}},
        )
        p = Panda(
            tables=m.tables.scores,
            operations=[Panda.sort("score", ascending=False), Panda.head(2)],
            pool=m,
        )
        Mock(source=p.streams.s, children="<1..2>",
             streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
        script_path = pipeline.save()

    ids = list(p.streams.s.ids)
    record_case(input="sort+head(2) cycle", expected=2, actual=len(ids))
    assert len(ids) == 2
    assert all("_" in i for i in ids)
    assert_valid_script(script_path, "Mock", "Panda", "panda_cycle_sort_head")


def test_panda_cycle_sort_tail_renames(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda(sort+tail, pool) also renames; bottom-N cycle."""
    from biopipelines.mock import Mock
    from biopipelines.panda import Panda

    pipeline = new_pipeline("panda_cycle_sort_tail")
    with pipeline:
        m = Mock(
            ids=["t1", "t2", "t3", "t4", "t5"],
            streams={"s": {"format": "pdb", "file": "<id>.pdb"}},
            tables={"scores": {"columns": ["score"], "fill": {"score": 0.4}}},
        )
        p = Panda(
            tables=m.tables.scores,
            operations=[Panda.sort("score"), Panda.tail(3)],
            pool=m,
        )
        Mock(source=p.streams.s, children="<A B>",
             streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
        script_path = pipeline.save()

    ids = list(p.streams.s.ids)
    record_case(input="sort+tail(3) cycle", expected=3, actual=len(ids))
    assert len(ids) == 3
    assert_valid_script(script_path, "Mock", "Panda", "panda_cycle_sort_tail")


def test_panda_cycle_filter_sort_head_combined(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda(filter + sort + head, pool) — the rename path still kicks in."""
    from biopipelines.mock import Mock
    from biopipelines.panda import Panda

    pipeline = new_pipeline("panda_cycle_filter_sort_head")
    with pipeline:
        m = Mock(
            ids=["u1", "u2", "u3", "u4", "u5", "u6"],
            streams={"s": {"format": "pdb", "file": "<id>.pdb"}},
            tables={"scores": {"columns": ["score"], "fill": {"score": 0.9}}},
        )
        p = Panda(
            tables=m.tables.scores,
            operations=[
                Panda.filter("score >= 0.5"),
                Panda.sort("score", ascending=False),
                Panda.head(3),
            ],
            pool=m,
        )
        Mock(source=p.streams.s, children="<1..2>",
             streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
        script_path = pipeline.save()

    ids = list(p.streams.s.ids)
    record_case(input="filter+sort+head(3) cycle", expected=3, actual=len(ids))
    assert len(ids) == 3
    assert_valid_script(
        script_path, "Mock", "Panda", "panda_cycle_filter_sort_head"
    )


def test_panda_cycle_sample_renames(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda(sample(n), pool) renames and fixes output count."""
    from biopipelines.mock import Mock
    from biopipelines.panda import Panda

    pipeline = new_pipeline("panda_cycle_sample")
    with pipeline:
        m = Mock(
            ids=["v1", "v2", "v3", "v4"],
            streams={"s": {"format": "pdb", "file": "<id>.pdb"}},
            tables={"scores": {"columns": ["score"], "fill": {"score": 0.6}}},
        )
        p = Panda(
            tables=m.tables.scores,
            operations=[Panda.sample(n=2, random_state=1)],
            pool=m,
        )
        Mock(source=p.streams.s, children="<A B>",
             streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
        script_path = pipeline.save()

    ids = list(p.streams.s.ids)
    record_case(input="sample(n=2) cycle", expected=2, actual=len(ids))
    assert len(ids) == 2
    assert_valid_script(script_path, "Mock", "Panda", "panda_cycle_sample")


# ── Runtime: pool-mode output content survival ───────────────────────────────
#
# These tests actually invoke pipe_sequence.py / pipe_ligand.py / pipe_panda.py
# and then read the Panda output folder to confirm:
#   (a) every pool_table_map'd content CSV (sequences.csv, compounds.csv)
#       still carries its upstream content columns after a round-trip through
#       Panda's pool mode, for both filter-only (no rename) and sort+head
#       (rename) cycles;
#   (b) every stream_map_target CSV (for streams with per-ID files) has the
#       canonical `id, file, <stream>.id` provenance schema and points at the
#       correct rows;
#   (c) the two are not at odds — when pool_table_map and stream_map_target
#       resolve to the same basename (e.g. `sequences.csv`), the content
#       must win, not the provenance.
#
# (c) is the regression lock for the 2026-04-23 bug where `sequences.csv`
# was silently overwritten and MutationProfiler failed with
# "Could not find id and sequence columns".


def test_panda_pool_sequence_preserves_content_columns(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """Sequence → Panda(pool=seq, filter) — post-execution `sequences.csv`
    must retain `id` and `sequence` (and `type`, `length`), not be
    overwritten with the stream-map schema `id, file, sequences.id`."""
    from biopipelines.sequence import Sequence
    from biopipelines.panda import Panda

    pipeline = new_pipeline("panda_pool_seq_content")
    with pipeline:
        seq = Sequence(seq=["MKTAYIAK", "AE", "GGGGALV"],
                       ids=["long1", "short", "long2"], type="protein")
        pan = Panda(
            tables=seq.tables.sequences,
            operations=[Panda.filter("length >= 5")],
            pool=seq,
        )
        pipeline.save()

    _run_pipe("sequence", os.path.join(seq.output_folder, "_configuration", "sequence_config.json"))
    _run_pipe("panda", os.path.join(pan.output_folder, "_configuration", "panda_config.json"))

    out_seq = os.path.join(pan.output_folder, "sequences", "sequences.csv")
    assert os.path.exists(out_seq), f"Panda did not emit {out_seq}"
    cols = _read_csv_columns(out_seq)
    rows = _read_csv_rows(out_seq)
    ids = sorted(r["id"] for r in rows)
    seqs = {r["id"]: r.get("sequence", "") for r in rows}

    record_case(
        input="Sequence[long1,short,long2] → Panda(pool, filter length>=5)",
        expected=({"id", "sequence"}, ["long1", "long2"],
                  {"long1": "MKTAYIAK", "long2": "GGGGALV"}),
        actual=(set(cols) & {"id", "sequence"}, ids, seqs),
    )
    assert {"id", "sequence"}.issubset(cols), (
        f"content columns lost from sequences.csv — got {cols}"
    )
    assert ids == ["long1", "long2"]
    assert seqs["long1"] == "MKTAYIAK"
    assert seqs["long2"] == "GGGGALV"


def test_panda_pool_sequence_rename_preserves_content_columns(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """Sequence → Panda(pool=seq, sort+head) triggers the rename path
    (`Panda_1`, `Panda_2`). Even then, `sequences.csv` must keep the
    `sequence` content column; `pool.id` provenance is expected as the
    only new column."""
    from biopipelines.sequence import Sequence
    from biopipelines.panda import Panda

    pipeline = new_pipeline("panda_pool_seq_rename")
    with pipeline:
        seq = Sequence(seq=["MKTAY", "AETGFLMK", "GG", "HHHHHKKK"],
                       ids=["a", "b", "c", "d"], type="protein")
        pan = Panda(
            tables=seq.tables.sequences,
            operations=[Panda.sort("length", ascending=False), Panda.head(2)],
            pool=seq,
        )
        pipeline.save()

    _run_pipe("sequence", os.path.join(seq.output_folder, "_configuration", "sequence_config.json"))
    _run_pipe("panda", os.path.join(pan.output_folder, "_configuration", "panda_config.json"))

    out_seq = os.path.join(pan.output_folder, "sequences", "sequences.csv")
    assert os.path.exists(out_seq)
    rows = _read_csv_rows(out_seq)
    cols = set(rows[0].keys()) if rows else set()
    ids = sorted(r["id"] for r in rows)
    parents = sorted(r.get("pool.id", "") for r in rows)
    longest_two_seqs = sorted(r["sequence"] for r in rows)

    record_case(
        input="Sequence[a,b,c,d lengths 5,8,2,8] → Panda(sort+head(2), pool)",
        expected=("sequence col kept", ["Panda_1", "Panda_2"],
                  ["b", "d"], ["AETGFLMK", "HHHHHKKK"]),
        actual=("sequence" in cols, ids, parents, longest_two_seqs),
    )
    assert "sequence" in cols, f"sequence column lost after rename — cols={cols}"
    assert ids == ["Panda_1", "Panda_2"]
    assert set(parents) == {"b", "d"}
    assert longest_two_seqs == ["AETGFLMK", "HHHHHKKK"]


def test_panda_pool_ligand_preserves_smiles_column(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """Ligand(smiles={...}) → Panda(pool=lig, filter) — `compounds.csv`
    must retain the `smiles` column (and the other content cols) after
    pool round-trip. Guards against the same stream_map_target
    overwrite for Ligand's `compounds` stream."""
    from biopipelines.ligand import Ligand
    from biopipelines.panda import Panda

    pipeline = new_pipeline("panda_pool_ligand_content")
    with pipeline:
        lig = Ligand(smiles={"ethanol": "CCO",
                             "acetate": "CC(=O)O",
                             "water": "O"})
        pan = Panda(
            tables=lig.tables.compounds,
            operations=[Panda.filter("id != 'water'")],
            pool=lig,
        )
        pipeline.save()

    _run_pipe("ligand", os.path.join(lig.output_folder, "_configuration", "fetch_config.json"))
    _run_pipe("panda", os.path.join(pan.output_folder, "_configuration", "panda_config.json"))

    out_cmp = os.path.join(pan.output_folder, "compounds", "compounds.csv")
    assert os.path.exists(out_cmp), f"Panda did not emit {out_cmp}"
    rows = _read_csv_rows(out_cmp)
    cols = set(rows[0].keys()) if rows else set()
    smiles_by_id = {r["id"]: r.get("smiles", "") for r in rows}

    record_case(
        input="Ligand(smiles={ethanol,acetate,water}) → Panda(pool, filter !=water)",
        expected=({"id", "smiles"}, {"ethanol": "CCO", "acetate": "CC(=O)O"}),
        actual=(cols & {"id", "smiles"}, smiles_by_id),
    )
    assert {"id", "smiles"}.issubset(cols), (
        f"content columns lost from compounds.csv — got {cols}"
    )
    assert smiles_by_id == {"ethanol": "CCO", "acetate": "CC(=O)O"}


def test_panda_pool_writes_stream_map_tables_with_provenance(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """Mock(structures stream with per-ID .pdb files) → Panda(pool, sort+head).

    The stream_map_target writer must still produce the structures
    map_table with `id, file, structures.id` provenance columns, pointing
    each renamed `Panda_N` back to its upstream parent. This locks in the
    behavior the bug-fix preserves for streams that genuinely need a map
    table (non-empty ``file_template``)."""
    from biopipelines.mock import Mock
    from biopipelines.panda import Panda

    pipeline = new_pipeline("panda_pool_stream_map")
    with pipeline:
        m = Mock(
            ids=["a", "b", "c", "d"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            tables={"scores": {"columns": ["score"],
                               "fill": {"score": 0.9}}},
        )
        pan = Panda(
            tables=m.tables.scores,
            operations=[Panda.sort("score", ascending=False), Panda.head(2)],
            pool=m,
        )
        pipeline.save()

    _run_pipe("mock", os.path.join(m.output_folder, "_configuration", "mock_config.json"),
              config_flag=False)
    _run_pipe("panda", os.path.join(pan.output_folder, "_configuration", "panda_config.json"))

    map_path = pan.streams.structures.map_table
    assert os.path.exists(map_path), f"stream map_table missing: {map_path}"
    rows = _read_csv_rows(map_path)
    cols = set(rows[0].keys()) if rows else set()
    rename_map = {r["id"]: r["structures.id"] for r in rows}

    record_case(
        input="Mock(4 ids) → Panda(sort+head(2), pool) — structures map_table",
        expected=({"id", "file", "structures.id"}, {"Panda_1", "Panda_2"}),
        actual=(cols, set(rename_map.keys())),
    )
    assert cols == {"id", "file", "structures.id"}, (
        f"unexpected stream map_table schema: {cols}"
    )
    assert set(rename_map.keys()) == {"Panda_1", "Panda_2"}
    # Each renamed ID must link to a concrete upstream parent.
    assert set(rename_map.values()).issubset({"a", "b", "c", "d"})


def test_panda_pool_content_and_stream_map_coexist(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """Sequence → Panda(pool, filter) emits BOTH a content `sequences.csv`
    (with `id, sequence, ...`) AND — if `sequences` is also emitted as a
    stream_map_target — something that carries provenance separately.

    The bug was: one path overwrites the other at the same filename.
    This test asserts that after the fix, `sequences.csv` is the content
    table (not a lineage stub), regardless of whether a sibling provenance
    file exists."""
    from biopipelines.sequence import Sequence
    from biopipelines.panda import Panda

    pipeline = new_pipeline("panda_pool_content_wins")
    with pipeline:
        seq = Sequence(seq=["MKTAY", "AETGFLMK"], ids=["a", "b"],
                       type="protein")
        pan = Panda(
            tables=seq.tables.sequences,
            operations=[Panda.filter("length >= 4")],
            pool=seq,
        )
        pipeline.save()

    _run_pipe("sequence", os.path.join(seq.output_folder, "_configuration", "sequence_config.json"))
    _run_pipe("panda", os.path.join(pan.output_folder, "_configuration", "panda_config.json"))

    out_seq = os.path.join(pan.output_folder, "sequences", "sequences.csv")
    rows = _read_csv_rows(out_seq)
    cols = set(rows[0].keys()) if rows else set()

    # The content table must win: `sequence` column is the strong signal.
    # The lineage-only schema would be {id, file, sequences.id} — explicitly
    # fail if we see that shape (no `sequence` key).
    is_content_schema = "sequence" in cols
    is_lineage_schema = cols == {"id", "file", "sequences.id"}

    record_case(
        input="Sequence → Panda(pool, filter) sequences.csv schema",
        expected=(True, False),
        actual=(is_content_schema, is_lineage_schema),
    )
    assert is_content_schema, (
        f"sequences.csv is a lineage stub, not a content table: cols={cols}"
    )
    assert not is_lineage_schema
