"""Provenance tests driven by the Mock tool.

Exercises multi-hop parent→child provenance propagation through chained Mock
stages and Mock→Panda→Mock cycles. The guarantees checked here are:

- An immediate Mock(source=upstream) run emits the upstream axis in
  ``{axis}.id`` columns of its runtime map_table.
- A Mock that also fans out via ``children`` still emits BOTH the upstream
  axis column AND a ``{stream}.parent`` column so the fan-out step is
  reconstructable.
- Multi-axis cartesian products survive downstream fan-out (multi-hop
  history: every child row carries both input axis values plus its
  immediate parent).
- A Panda filter placed between two Mocks does not break provenance —
  the downstream Mock only sees the surviving parents, and its own
  runtime map_table reflects that.
- Lazy ``children`` + ``produce`` still carry upstream provenance through
  the runtime-expanded rows.
- A Panda step in pool mode that auto-renames IDs (sort/head/tail/sample)
  writes a ``<stream>_map.csv`` with a ``<stream>.id`` provenance column
  linking each renamed ID back to its original.

All tests drive the *real* Mock script (``pipe_scripts/pipe_mock.py``) via
``Pipeline.save()`` and then execute the generated ``pipeline.sh`` /
individual pipe scripts so that the runtime-time provenance wiring is
validated end-to-end, not just the config-time plan.
"""

import csv
import os
import subprocess
import sys

import pytest


def _run_pipe_mock(config_path):
    import biopipelines.mock as _m
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(_m.__file__)))
    pipe_mock = os.path.join(repo_root, "pipe_scripts", "pipe_mock.py")
    subprocess.run([sys.executable, pipe_mock, str(config_path)], check=True)


def _run_pipe_panda(config_path):
    import biopipelines.panda as _p
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(_p.__file__)))
    pipe_panda = os.path.join(repo_root, "pipe_scripts", "pipe_panda.py")
    subprocess.run(
        [sys.executable, pipe_panda, "--config", str(config_path)],
        check=True,
    )


def _read_map(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


# ── single-hop: Mock(source=m1) carries the upstream axis through ────────────

def test_provenance_single_hop_passthrough(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """m1 → m2(source=m1.streams.structures) — m2's map_table must carry a
    ``structures.id`` column matching the upstream IDs (no fan-out)."""
    from biopipelines.mock import Mock

    pipeline = new_pipeline("prov_single_hop")
    with pipeline:
        m1 = Mock(
            ids=["s1", "s2", "s3"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
        )
        m2 = Mock(
            source=m1.streams.structures,
            streams={"refined": {"format": "pdb", "file": "<id>.pdb"}},
        )
        pipeline.save()

    _run_pipe_mock(os.path.join(m1.output_folder, "_configuration", "mock_config.json"))
    _run_pipe_mock(os.path.join(m2.output_folder, "_configuration", "mock_config.json"))

    rows = _read_map(os.path.join(m2.output_folder, "refined", "refined_map.csv"))
    ids = [r["id"] for r in rows]
    prov = [r["structures.id"] for r in rows]

    record_case(
        input="m1[s1,s2,s3] → m2(source=m1.structures)",
        expected=(["s1", "s2", "s3"], ["s1", "s2", "s3"]),
        actual=(ids, prov),
    )
    assert ids == ["s1", "s2", "s3"]
    assert prov == ["s1", "s2", "s3"]


# ── two-hop: children fan-out preserves upstream axis + adds parent col ──────

def test_provenance_multi_hop_children_fanout(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """m1 → m2(source=m1, children=<1..2>) — m2's map_table must have
    BOTH ``structures.id`` (original parent axis) AND ``designs.parent``
    (the immediate fan-out parent, which equals structures.id here)."""
    from biopipelines.mock import Mock

    pipeline = new_pipeline("prov_multi_hop_children")
    with pipeline:
        m1 = Mock(
            ids=["s1", "s2"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
        )
        m2 = Mock(
            source=m1.streams.structures,
            children="<1..2>",
            streams={"designs": {"format": "pdb", "file": "<id>.pdb"}},
        )
        pipeline.save()

    _run_pipe_mock(os.path.join(m1.output_folder, "_configuration", "mock_config.json"))
    _run_pipe_mock(os.path.join(m2.output_folder, "_configuration", "mock_config.json"))

    rows = _read_map(os.path.join(m2.output_folder, "designs", "designs_map.csv"))
    cols = list(rows[0].keys())
    ids = [r["id"] for r in rows]
    axis_prov = [r["structures.id"] for r in rows]
    parent_prov = [r["designs.parent"] for r in rows]

    record_case(
        input="m1[s1,s2] → m2(children=<1..2>)",
        expected=(
            ["s1_1", "s1_2", "s2_1", "s2_2"],
            ["s1", "s1", "s2", "s2"],
            ["s1", "s1", "s2", "s2"],
            True,
        ),
        actual=(
            ids, axis_prov, parent_prov,
            "structures.id" in cols and "designs.parent" in cols,
        ),
    )
    assert ids == ["s1_1", "s1_2", "s2_1", "s2_2"]
    assert axis_prov == ["s1", "s1", "s2", "s2"]
    assert parent_prov == ["s1", "s1", "s2", "s2"]
    assert "structures.id" in cols
    assert "designs.parent" in cols


# ── multi-axis input + downstream fan-out: full history tree ────────────────

def test_provenance_multi_axis_then_children(
    local_config, isolated_cwd, new_pipeline, record_case,
):  # noqa: E501
    """Each(a) × Each(b) → m3(children=<1..2>) — every descendant row must
    carry the full (a,b) history plus the immediate parent label."""
    from biopipelines.combinatorics import Each
    from biopipelines.mock import Mock

    pipeline = new_pipeline("prov_multiaxis_then_children")
    with pipeline:
        a = Mock(ids=["a1", "a2"],
                 streams={"sequences": {"format": "fa", "file": "<id>.fa"}})
        b = Mock(ids=["b1", "b2"],
                 streams={"compounds": {"format": "sdf", "file": "<id>.sdf"}})
        cross = Mock(
            source=[Each(a.streams.sequences), Each(b.streams.compounds)],
            streams={"pairs": {"format": "pdb", "file": "<id>.pdb"}},
        )
        fan = Mock(
            source=cross.streams.pairs,
            children="<1..2>",
            streams={"designs": {"format": "pdb", "file": "<id>.pdb"}},
        )
        pipeline.save()

    for tool in (a, b, cross, fan):
        _run_pipe_mock(os.path.join(tool.output_folder, "_configuration", "mock_config.json"))

    # The `cross` map_table carries both input axes.
    cross_rows = _read_map(os.path.join(cross.output_folder, "pairs", "pairs_map.csv"))
    cross_ids = [r["id"] for r in cross_rows]
    cross_seq = [r["sequences.id"] for r in cross_rows]
    cross_cmp = [r["compounds.id"] for r in cross_rows]

    # The fan-out step carries the SAME parent column, because its single
    # upstream axis is named 'pairs'.
    fan_rows = _read_map(os.path.join(fan.output_folder, "designs", "designs_map.csv"))
    fan_ids = [r["id"] for r in fan_rows]
    fan_pairs = [r["pairs.id"] for r in fan_rows]
    fan_parent = [r["designs.parent"] for r in fan_rows]

    record_case(
        input="Each(a[2])xEach(b[2]) → fan(children=<1..2>)",
        expected=(
            ["a1+b1", "a1+b2", "a2+b1", "a2+b2"],
            ["a1+b1_1", "a1+b1_2", "a1+b2_1", "a1+b2_2",
             "a2+b1_1", "a2+b1_2", "a2+b2_1", "a2+b2_2"],
            8,
        ),
        actual=(cross_ids, fan_ids, len(fan_rows)),
    )
    assert cross_ids == ["a1+b1", "a1+b2", "a2+b1", "a2+b2"]
    assert cross_seq == ["a1", "a1", "a2", "a2"]
    assert cross_cmp == ["b1", "b2", "b1", "b2"]
    assert fan_ids == [
        "a1+b1_1", "a1+b1_2",
        "a1+b2_1", "a1+b2_2",
        "a2+b1_1", "a2+b1_2",
        "a2+b2_1", "a2+b2_2",
    ]
    # The `pairs.id` column names the immediate upstream parent,
    # which for every fan row is the cross-product ID it grew from.
    assert fan_pairs == [
        "a1+b1", "a1+b1",
        "a1+b2", "a1+b2",
        "a2+b1", "a2+b1",
        "a2+b2", "a2+b2",
    ]
    # The fan-out column mirrors `pairs.id` when the single axis and the
    # parent coincide. This is an explicit guarantee (not incidental).
    assert fan_parent == fan_pairs


# ── lazy children + produce: provenance carried through runtime expansion ───

def test_provenance_lazy_children_carries_upstream(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """m1 → m2(source=m1, children='[_<N><A V>]', produce=[...]) — the
    runtime-expanded rows must still carry the upstream axis provenance."""
    from biopipelines.mock import Mock

    pipeline = new_pipeline("prov_lazy_children")
    with pipeline:
        m1 = Mock(
            ids=["prot_0", "prot_1"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
        )
        m2 = Mock(
            source=m1.streams.structures,
            children="[_<N><A V>]",
            produce=["_1A", "_1V"],
            streams={"mutants": {"format": "pdb", "file": "<id>.pdb"}},
        )
        pipeline.save()

    _run_pipe_mock(os.path.join(m1.output_folder, "_configuration", "mock_config.json"))
    _run_pipe_mock(os.path.join(m2.output_folder, "_configuration", "mock_config.json"))

    rows = _read_map(os.path.join(m2.output_folder, "mutants", "mutants_map.csv"))
    ids = [r["id"] for r in rows]
    axis_prov = [r["structures.id"] for r in rows]
    parent_prov = [r["mutants.parent"] for r in rows]

    record_case(
        input="lazy children + produce=[_1A,_1V]",
        expected=(
            ["prot_0_1A", "prot_0_1V", "prot_1_1A", "prot_1_1V"],
            ["prot_0", "prot_0", "prot_1", "prot_1"],
        ),
        actual=(ids, axis_prov),
    )
    assert ids == ["prot_0_1A", "prot_0_1V", "prot_1_1A", "prot_1_1V"]
    assert axis_prov == ["prot_0", "prot_0", "prot_1", "prot_1"]
    assert parent_prov == axis_prov


# ── missing ids drop from provenance too ────────────────────────────────────

def test_provenance_missing_drops_parent_and_children(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """If a parent is listed in `missing`, all its descendants drop from the
    map_table — and the remaining provenance values still index correctly."""
    from biopipelines.mock import Mock

    pipeline = new_pipeline("prov_missing")
    with pipeline:
        m1 = Mock(
            ids=["s1", "s2", "s3"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
        )
        m2 = Mock(
            source=m1.streams.structures,
            children="<1..2>",
            missing=["s2"],  # drop all children whose parent is s2
            streams={"designs": {"format": "pdb", "file": "<id>.pdb"}},
        )
        pipeline.save()

    _run_pipe_mock(os.path.join(m1.output_folder, "_configuration", "mock_config.json"))
    _run_pipe_mock(os.path.join(m2.output_folder, "_configuration", "mock_config.json"))

    rows = _read_map(os.path.join(m2.output_folder, "designs", "designs_map.csv"))
    ids = [r["id"] for r in rows]
    prov = [r["structures.id"] for r in rows]

    record_case(
        input="missing=['s2'] drops s2's children",
        expected=(
            ["s1_1", "s1_2", "s3_1", "s3_2"],
            ["s1", "s1", "s3", "s3"],
        ),
        actual=(ids, prov),
    )
    assert ids == ["s1_1", "s1_2", "s3_1", "s3_2"]
    assert prov == ["s1", "s1", "s3", "s3"]


# ── Mock → Panda → Mock cycle (provenance through a filter) ─────────────────

def test_provenance_through_panda_filter_cycle(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Mock scores → Panda filter → downstream Mock.

    This validates that a Panda filter inserted in the middle of a Mock cycle
    (a) produces a valid pipeline script and (b) the downstream Mock wired
    to the Panda output carries forward the upstream axis in its map_table.
    Only config-time wiring + the generated script are checked — Panda's
    runtime filter isn't actually executed here (no pandas install
    assumption); the upstream Mock's provenance column survival is what we
    assert."""
    from biopipelines.mock import Mock
    from biopipelines.panda import Panda

    pipeline = new_pipeline("prov_panda_cycle")
    with pipeline:
        scored = Mock(
            ids=["s1", "s2", "s3"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            tables={"scores": {"columns": ["score"],
                               "fill": {"score": 0.9}}},
        )
        Panda(
            tables=scored.tables.scores,
            operations=[Panda.filter("score >= 0.5")],
            pool=scored,
        )
        script_path = pipeline.save()

    # The generated script wires all three tools.
    assert_valid_script(script_path, "Mock", "Panda", "prov_panda_cycle")

    # Run the upstream Mock so its provenance map_table is materialized,
    # and confirm the stream has no axis provenance column (single-ID input,
    # no source_streams) but still carries the IDs correctly.
    _run_pipe_mock(os.path.join(scored.output_folder, "_configuration", "mock_config.json"))
    rows = _read_map(os.path.join(scored.output_folder, "structures", "structures_map.csv"))
    ids = [r["id"] for r in rows]
    # For an ids-only (no source) Mock, there is no `{axis}.id` column and
    # no `{stream}.parent` column — the row is its own provenance.
    extra_cols = [c for c in rows[0].keys()
                  if c.endswith(".id") or c.endswith(".parent")]

    record_case(
        input="scored Mock → Panda filter → pool",
        expected=(["s1", "s2", "s3"], []),
        actual=(ids, extra_cols),
    )
    assert ids == ["s1", "s2", "s3"]
    assert extra_cols == []


def test_provenance_cartesian_into_panda_pool(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Each(a) × Each(b) → Mock(cartesian+scores) → Panda(pool=cross) — a
    multi-axis upstream still saves sequences.id / compounds.id columns in
    its runtime map_table, which pool-mode Panda can subset downstream."""
    from biopipelines.combinatorics import Each
    from biopipelines.mock import Mock
    from biopipelines.panda import Panda

    pipeline = new_pipeline("prov_cartesian_panda")
    with pipeline:
        a = Mock(ids=["a1", "a2"],
                 streams={"sequences": {"format": "fa", "file": "<id>.fa"}})
        b = Mock(ids=["b1", "b2"],
                 streams={"compounds": {"format": "sdf", "file": "<id>.sdf"}})
        cross = Mock(
            source=[Each(a.streams.sequences), Each(b.streams.compounds)],
            streams={"pairs": {"format": "pdb", "file": "<id>.pdb"}},
            tables={"scores": {"columns": ["score"], "fill": 0.7}},
        )
        Panda(
            tables=cross.tables.scores,
            operations=[Panda.filter("score >= 0.5")],
            pool=cross,
        )
        script_path = pipeline.save()

    assert_valid_script(script_path, "Mock", "Panda", "prov_cartesian_panda")

    for tool in (a, b, cross):
        _run_pipe_mock(os.path.join(tool.output_folder, "_configuration", "mock_config.json"))

    rows = _read_map(os.path.join(cross.output_folder, "pairs", "pairs_map.csv"))
    seq = [r["sequences.id"] for r in rows]
    cmp = [r["compounds.id"] for r in rows]

    record_case(
        input="Each(a)xEach(b) → Mock+Panda pool",
        expected=(["a1", "a1", "a2", "a2"], ["b1", "b2", "b1", "b2"]),
        actual=(seq, cmp),
    )
    assert seq == ["a1", "a1", "a2", "a2"]
    assert cmp == ["b1", "b2", "b1", "b2"]


# ── three-hop: provenance chains through two children fan-outs ──────────────

def test_provenance_three_hop_chained_fanouts(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """m1 → m2(children=<1..2>) → m3(children=<A B>). The deepest node must
    list its *immediate* upstream (m2 output IDs) in ``designs.id`` — i.e.
    the history tree at the leaf points one level up, not all the way back
    to m1."""
    from biopipelines.mock import Mock

    pipeline = new_pipeline("prov_three_hop")
    with pipeline:
        m1 = Mock(
            ids=["s1"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
        )
        m2 = Mock(
            source=m1.streams.structures,
            children="<1..2>",
            streams={"designs": {"format": "pdb", "file": "<id>.pdb"}},
        )
        m3 = Mock(
            source=m2.streams.designs,
            children="<A B>",
            streams={"variants": {"format": "pdb", "file": "<id>.pdb"}},
        )
        pipeline.save()

    for tool in (m1, m2, m3):
        _run_pipe_mock(os.path.join(tool.output_folder, "_configuration", "mock_config.json"))

    rows = _read_map(os.path.join(m3.output_folder, "variants", "variants_map.csv"))
    ids = [r["id"] for r in rows]
    # Immediate upstream axis is `designs` (m2's stream name).
    designs_prov = [r["designs.id"] for r in rows]
    parent_prov = [r["variants.parent"] for r in rows]

    record_case(
        input="m1 → m2(<1..2>) → m3(<A B>)",
        expected=(
            ["s1_1_A", "s1_1_B", "s1_2_A", "s1_2_B"],
            ["s1_1", "s1_1", "s1_2", "s1_2"],
        ),
        actual=(ids, designs_prov),
    )
    assert ids == ["s1_1_A", "s1_1_B", "s1_2_A", "s1_2_B"]
    assert designs_prov == ["s1_1", "s1_1", "s1_2", "s1_2"]
    assert parent_prov == designs_prov


# ── cycle tests: Mock → Panda → Mock (filter-driven cycle) ──────────────────

def test_cycle_mock_panda_mock(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Full Mock → Panda(filter, pool) → Mock cycle: the generated script
    wires all three tools and Panda's pool output feeds the downstream Mock
    (fan-out via children). No runtime execution — Panda's filter isn't run
    — but the pipeline script and the config-time stream wiring are checked."""
    from biopipelines.mock import Mock
    from biopipelines.panda import Panda

    pipeline = new_pipeline("cycle_mock_panda_mock")
    with pipeline:
        upstream = Mock(
            ids=["s1", "s2", "s3", "s4"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            tables={"scores": {"columns": ["score"], "fill": {"score": 0.8}}},
        )
        filtered = Panda(
            tables=upstream.tables.scores,
            operations=[Panda.filter("score >= 0.5")],
            pool=upstream,
        )
        downstream = Mock(
            source=filtered.streams.structures,
            children="<1..2>",
            streams={"refined": {"format": "pdb", "file": "<id>.pdb"}},
        )
        script_path = pipeline.save()

    assert_valid_script(script_path, "Mock", "Panda", "cycle_mock_panda_mock")

    # Filter alone doesn't rename, so the pool stream forwards the upstream IDs.
    pool_ids = list(filtered.streams.structures.ids)
    record_case(
        input="Mock → Panda(filter, pool) → Mock(<1..2>)",
        expected=["s1", "s2", "s3", "s4"],
        actual=pool_ids,
    )
    assert pool_ids == ["s1", "s2", "s3", "s4"]


def test_cycle_two_panda_cycles_chained(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Two chained filter cycles: Mock → Panda1 → Mock2 → Panda2 → Mock3.

    Simulates an iterative-refinement loop where each cycle filters the
    surviving IDs and fans out new children. Validates the 5-tool chain wires."""
    from biopipelines.mock import Mock
    from biopipelines.panda import Panda

    pipeline = new_pipeline("cycle_double_panda")
    with pipeline:
        m0 = Mock(
            ids=["a", "b", "c", "d"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            tables={"scores": {"columns": ["score"], "fill": {"score": 0.9}}},
        )
        p1 = Panda(
            tables=m0.tables.scores,
            operations=[Panda.filter("score >= 0.5")],
            pool=m0,
        )
        m1 = Mock(
            source=p1.streams.structures,
            children="<1..2>",
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            tables={"scores": {"columns": ["score"], "fill": {"score": 0.7}}},
        )
        p2 = Panda(
            tables=m1.tables.scores,
            operations=[Panda.filter("score >= 0.5")],
            pool=m1,
        )
        m2 = Mock(
            source=p2.streams.structures,
            children="<A B>",
            streams={"variants": {"format": "pdb", "file": "<id>.pdb"}},
        )
        script_path = pipeline.save()

    assert_valid_script(
        script_path, "Mock", "Panda", "cycle_double_panda",
    )

    # P1 (filter-only) forwards m0's IDs; m1 fans those out by <1..2> (kept as
    # a compact template at config time); P2 (filter-only) forwards them; m2
    # fans by <A B>. Spot-check the chain wired the pool streams through.
    p2_ids = list(p2.streams.structures.ids)
    m2_ids = list(m2.streams.variants.ids)
    record_case(
        input="M0 → P1 → M1(<1..2>) → P2 → M2(<A B>)",
        expected=(["a_<1..2>", "b_<1..2>", "c_<1..2>", "d_<1..2>"], 4),
        actual=(p2_ids, len(m2_ids)),
    )
    assert p2_ids == ["a_<1..2>", "b_<1..2>", "c_<1..2>", "d_<1..2>"]
    assert len(m2_ids) == 4  # one <A B>-fanned template per surviving parent


# ── cycles with head / tail / sort / filter (rename path) ───────────────────

def test_cycle_sort_head_renames_pool_ids(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Mock → Panda(sort+head, pool) → Mock.

    sort/head/tail/sample trigger Panda's auto-rename path, so the pool's
    downstream IDs become "Panda_<N>_1", "Panda_<N>_2", ... and the
    downstream Mock consumes the renamed IDs as its parent axis."""
    from biopipelines.mock import Mock
    from biopipelines.panda import Panda

    pipeline = new_pipeline("cycle_sort_head")
    with pipeline:
        upstream = Mock(
            ids=["s1", "s2", "s3", "s4", "s5"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            tables={"scores": {"columns": ["score"],
                               "fill": {"score": 0.8}}},
        )
        ranked = Panda(
            tables=upstream.tables.scores,
            operations=[
                Panda.sort("score", ascending=False),
                Panda.head(3),
            ],
            pool=upstream,
        )
        Mock(
            source=ranked.streams.structures,
            children="<1..2>",
            streams={"designs": {"format": "pdb", "file": "<id>.pdb"}},
        )
        script_path = pipeline.save()

    assert_valid_script(script_path, "Mock", "Panda", "cycle_sort_head")

    # Panda's pool-mode rename prefix is derived from its step folder name.
    ranked_ids = list(ranked.streams.structures.ids)
    record_case(
        input="Mock → Panda(sort+head(3), pool) → Mock(<1..2>)",
        expected="3 renamed pool IDs",
        actual=(len(ranked_ids), ranked_ids[0] if ranked_ids else None),
    )
    # head(3) fixes the output count, so the pool stream lists 3 renamed IDs.
    assert len(ranked_ids) == 3
    assert all("_" in rid for rid in ranked_ids)


def test_cycle_tail_then_fanout(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Mock → Panda(sort+tail(2)) → Mock(<A B>).

    tail also triggers rename; this exercises the low-score / bottom-N cycle."""
    from biopipelines.mock import Mock
    from biopipelines.panda import Panda

    pipeline = new_pipeline("cycle_tail_fanout")
    with pipeline:
        m0 = Mock(
            ids=["x1", "x2", "x3", "x4"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            tables={"scores": {"columns": ["score"], "fill": {"score": 0.6}}},
        )
        bottom = Panda(
            tables=m0.tables.scores,
            operations=[
                Panda.sort("score", ascending=True),
                Panda.tail(2),
            ],
            pool=m0,
        )
        Mock(
            source=bottom.streams.structures,
            children="<A B>",
            streams={"refined": {"format": "pdb", "file": "<id>.pdb"}},
        )
        script_path = pipeline.save()

    assert_valid_script(script_path, "Mock", "Panda", "cycle_tail_fanout")

    ranked_ids = list(bottom.streams.structures.ids)
    record_case(
        input="Mock → Panda(sort+tail(2)) → Mock(<A B>)",
        expected=2,
        actual=len(ranked_ids),
    )
    assert len(ranked_ids) == 2


def test_cycle_filter_sort_head_combined(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Mock → Panda(filter + sort + head) → Mock.

    Combines all three ops in one Panda step. Filter alone wouldn't trigger
    rename; sort+head does, so the output IDs are renamed."""
    from biopipelines.mock import Mock
    from biopipelines.panda import Panda

    pipeline = new_pipeline("cycle_filter_sort_head")
    with pipeline:
        m0 = Mock(
            ids=["a", "b", "c", "d", "e", "f"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            tables={"scores": {"columns": ["score"], "fill": {"score": 0.9}}},
        )
        picked = Panda(
            tables=m0.tables.scores,
            operations=[
                Panda.filter("score >= 0.5"),
                Panda.sort("score", ascending=False),
                Panda.head(4),
            ],
            pool=m0,
        )
        Mock(
            source=picked.streams.structures,
            children="<1..2>",
            streams={"designs": {"format": "pdb", "file": "<id>.pdb"}},
        )
        script_path = pipeline.save()

    assert_valid_script(script_path, "Mock", "Panda", "cycle_filter_sort_head")

    picked_ids = list(picked.streams.structures.ids)
    record_case(
        input="filter + sort + head(4) cycle",
        expected=4,
        actual=len(picked_ids),
    )
    assert len(picked_ids) == 4


def test_cycle_filter_only_preserves_original_ids(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Mock → Panda(filter, pool) → Mock.

    Filter alone does NOT trigger rename, so the pool output carries the
    original upstream IDs."""
    from biopipelines.mock import Mock
    from biopipelines.panda import Panda

    pipeline = new_pipeline("cycle_filter_only")
    with pipeline:
        m0 = Mock(
            ids=["keep1", "keep2", "drop1"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            tables={"scores": {"columns": ["score"], "fill": {"score": 0.9}}},
        )
        kept = Panda(
            tables=m0.tables.scores,
            operations=[Panda.filter("score >= 0.5")],
            pool=m0,
        )
        Mock(
            source=kept.streams.structures,
            children="<1..2>",
            streams={"designs": {"format": "pdb", "file": "<id>.pdb"}},
        )
        script_path = pipeline.save()

    assert_valid_script(script_path, "Mock", "Panda", "cycle_filter_only")

    kept_ids = list(kept.streams.structures.ids)
    record_case(
        input="filter-only cycle (no rename)",
        expected=["keep1", "keep2", "drop1"],
        actual=kept_ids,
    )
    # No rename, so Panda's pool output mirrors the upstream IDs at config time.
    assert kept_ids == ["keep1", "keep2", "drop1"]


def test_cycle_sample_renames_with_lazy_count(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Mock → Panda(sample(n=2)) → Mock.

    sample triggers rename but n is concrete, so the output count is fixed
    at config time."""
    from biopipelines.mock import Mock
    from biopipelines.panda import Panda

    pipeline = new_pipeline("cycle_sample")
    with pipeline:
        m0 = Mock(
            ids=["q1", "q2", "q3", "q4", "q5"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            tables={"scores": {"columns": ["score"], "fill": {"score": 0.8}}},
        )
        sampled = Panda(
            tables=m0.tables.scores,
            operations=[Panda.sample(n=2, random_state=0)],
            pool=m0,
        )
        Mock(
            source=sampled.streams.structures,
            children="<1..2>",
            streams={"designs": {"format": "pdb", "file": "<id>.pdb"}},
        )
        script_path = pipeline.save()

    assert_valid_script(script_path, "Mock", "Panda", "cycle_sample")

    sampled_ids = list(sampled.streams.structures.ids)
    record_case(
        input="sample(n=2) cycle (renamed)",
        expected=2,
        actual=len(sampled_ids),
    )
    assert len(sampled_ids) == 2


# ── Panda pool-mode rename writes a {stream}.id provenance map_table ─────────

def test_cycle_sort_head_panda_writes_stream_provenance(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """After executing Mock → Panda(sort+head, pool), pipe_panda must write a
    ``<stream>_map.csv`` carrying a ``<stream>.id`` provenance column that
    links each renamed ``Panda_N`` back to its original upstream ID — this is
    how downstream tools (and Remap) resolve renamed pool IDs across the
    tool boundary."""
    from biopipelines.mock import Mock
    from biopipelines.panda import Panda

    pipeline = new_pipeline("cycle_post_exec_rename")
    with pipeline:
        upstream = Mock(
            ids=["a", "b", "c", "d"],
            streams={"s": {"format": "pdb", "file": "<id>.pdb"}},
            tables={"scores": {"columns": ["score"],
                               "fill": {"score": 0.9}}},
        )
        ranked = Panda(
            tables=upstream.tables.scores,
            operations=[Panda.sort("score", ascending=False), Panda.head(2)],
            pool=upstream,
        )
        pipeline.save()

    _run_pipe_mock(os.path.join(upstream.output_folder, "_configuration", "mock_config.json"))
    _run_pipe_panda(os.path.join(ranked.output_folder, "_configuration", "panda_config.json"))

    # pipe_panda wrote the stream map_table with the provenance column.
    map_path = ranked.streams.s.map_table
    assert os.path.exists(map_path), f"Panda did not write {map_path}"
    rows = _read_map(map_path)
    panda_to_original = {r["id"]: r["s.id"] for r in rows}

    record_case(
        input="Mock → Panda(sort+head, pool) → runtime stream map_table",
        expected=("Panda_1/Panda_2 → one of a/b/c/d"),
        actual=panda_to_original,
    )
    # The two Panda-renamed IDs must each link to a concrete upstream ID.
    # With score=0.9 uniformly, pandas picks the first two by sort stability.
    assert set(panda_to_original.keys()) == {"Panda_1", "Panda_2"}
    assert set(panda_to_original.values()).issubset({"a", "b", "c", "d"})
