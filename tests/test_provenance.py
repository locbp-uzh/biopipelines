"""Provenance / history-tree tests driven by the Mock tool.

Exercises multi-hop lineage propagation through chained Mock stages and
Mock→Panda→Mock cycles. The guarantees checked here are:

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

All tests drive the *real* Mock script (``pipe_scripts/pipe_mock.py``) via
``Pipeline.save()`` and then execute the generated ``pipeline.sh`` /
individual pipe scripts so that the runtime-time provenance wiring is
validated end-to-end, not just the config-time plan.
"""

import csv
import json
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


def _read_lineage(output_folder):
    """Read the pipeline's .lineage.csv from its output folder (one level
    above a tool's own output folder)."""
    lineage = os.path.join(output_folder, ".lineage.csv")
    with open(lineage, newline="") as f:
        reader = csv.reader(f)
        rows = list(reader)
    return rows[0], rows[1:]  # header, data_rows


# ── lineage plotting (columns = tools, nodes = IDs, arrows = lineage) ────────

_PLOTS_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "outputs", "lineage_plots"
)


def _plot_lineage(output_folder, title, png_name):
    """Render .lineage.csv as a PNG with columns per tool and arrows
    linking each row's consecutive non-empty cells.

    A '-' cell is drawn as a red "dropped" marker; every other cell is a
    node labeled by its ID. Duplicate IDs in the same tool column collapse
    to a single node. The image lands in tests/outputs/lineage_plots/.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return None

    header, rows = _read_lineage(output_folder)
    if not rows:
        return None

    # Column positions (x), per-column node ordering (y).
    col_x = {label: i for i, label in enumerate(header)}
    col_nodes = {label: [] for label in header}  # preserves insertion order
    col_node_y = {label: {} for label in header}

    def _register(label, node_id):
        if node_id not in col_node_y[label]:
            col_node_y[label][node_id] = len(col_nodes[label])
            col_nodes[label].append(node_id)

    # First pass: register nodes.
    for row in rows:
        for label, cell in zip(header, row):
            if cell:
                _register(label, cell)

    # Second pass: edges between consecutive non-empty cells in each row.
    edges = []  # list of ((label_a, id_a), (label_b, id_b), dropped_flag)
    for row in rows:
        prev = None
        for label, cell in zip(header, row):
            if not cell:
                continue
            if prev is None:
                prev = (label, cell)
                continue
            dropped = (cell == "-")
            edges.append((prev, (label, cell), dropped))
            if not dropped:
                prev = (label, cell)

    # Widen the per-column x-spacing so long labels don't overlap across
    # columns and arrows have room to breathe. Node x-positions are scaled
    # by `col_spacing`; arrow gaps are expressed in the same data units.
    max_label_len = max(
        (len(str(nid)) for nids in col_nodes.values() for nid in nids),
        default=1,
    )
    col_spacing = max(2.0, 0.18 * max_label_len + 1.2)

    max_rows = max((len(v) for v in col_nodes.values()), default=1)
    fig_w = max(5, col_spacing * len(header) + 1)
    fig_h = max(2.5, 0.6 * max_rows + 1.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_title(title, fontsize=11)
    ax.set_xlim(-col_spacing * 0.5, col_spacing * (len(header) - 0.5))
    ax.set_ylim(-max_rows - 0.5, 1.4)
    ax.axis("off")

    # Column headers.
    for label, x in col_x.items():
        ax.text(x * col_spacing, 0.8, label, ha="center", va="center",
                fontsize=10, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="#eef",
                          edgecolor="#557"))

    # Node geometry: half-width grows with label length so the arrow-end
    # offset clears even long IDs; half-height is fixed.
    def _node_xy(label, node_id):
        return col_x[label] * col_spacing, -col_node_y[label][node_id]

    def _node_halfwidth(node_id):
        # Roughly matches the rendered box width (fontsize 9, pad 0.25).
        return 0.05 * max(1, len(str(node_id))) + 0.18

    for label in header:
        for node_id in col_nodes[label]:
            x, y = _node_xy(label, node_id)
            is_drop = (node_id == "-")
            ax.text(
                x, y, node_id, ha="center", va="center", fontsize=9,
                bbox=dict(
                    boxstyle="round,pad=0.25",
                    facecolor=("#fdd" if is_drop else "#efe"),
                    edgecolor=("#c33" if is_drop else "#484"),
                ),
            )

    # Draw arrows: start/end offset by each node's half-width so the arrow
    # terminates at the box edge rather than piercing through it.
    for (la, ia), (lb, ib), dropped in edges:
        xa, ya = _node_xy(la, ia)
        xb, yb = _node_xy(lb, ib)
        gap_a = _node_halfwidth(ia) + 0.05
        gap_b = _node_halfwidth(ib) + 0.05
        color = "#c33" if dropped else "#666"
        ax.annotate(
            "", xy=(xb - gap_b, yb), xytext=(xa + gap_a, ya),
            arrowprops=dict(
                arrowstyle="->", color=color, lw=1.4,
                connectionstyle="arc3,rad=0.05",
            ),
        )

    os.makedirs(_PLOTS_DIR, exist_ok=True)
    png_path = os.path.join(_PLOTS_DIR, png_name)
    fig.tight_layout()
    fig.savefig(png_path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return png_path


@pytest.fixture
def plot_lineage(request):
    """Factory: render the current pipeline's lineage.csv to a PNG named
    after the test. Call inside a test AFTER the lineage CSV has been
    written (and, for post-exec variants, after regeneration)."""
    def _plot(output_folder, suffix=""):
        safe = request.node.name.replace("/", "_").replace("::", "_")
        if suffix:
            safe = f"{safe}__{suffix}"
        return _plot_lineage(output_folder, request.node.name + (
            f" ({suffix})" if suffix else ""
        ), f"{safe}.png")
    return _plot


# ── single-hop: Mock(source=m1) carries the upstream axis through ────────────

def test_provenance_single_hop_passthrough(
    local_config, isolated_cwd, new_pipeline, record_case, plot_lineage,
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
    plot_lineage(os.path.dirname(m1.output_folder))


# ── two-hop: children fan-out preserves upstream axis + adds parent col ──────

def test_provenance_multi_hop_children_fanout(
    local_config, isolated_cwd, new_pipeline, record_case, plot_lineage,
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
    plot_lineage(os.path.dirname(m1.output_folder))


# ── multi-axis input + downstream fan-out: full history tree ────────────────

def test_provenance_multi_axis_then_children(
    local_config, isolated_cwd, new_pipeline, record_case, plot_lineage,
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
    plot_lineage(os.path.dirname(a.output_folder))


# ── lazy children + produce: provenance carried through runtime expansion ───

def test_provenance_lazy_children_carries_upstream(
    local_config, isolated_cwd, new_pipeline, record_case, plot_lineage,
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
    plot_lineage(os.path.dirname(m1.output_folder))


# ── missing ids drop from provenance too ────────────────────────────────────

def test_provenance_missing_drops_parent_and_children(
    local_config, isolated_cwd, new_pipeline, record_case, plot_lineage,
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
    plot_lineage(os.path.dirname(m1.output_folder), "pre-regen")
    pipeline._generate_id_lineage_csv()
    plot_lineage(os.path.dirname(m1.output_folder), "post-exec")


# ── Mock → Panda → Mock cycle (history tree through a filter) ───────────────

def test_provenance_through_panda_filter_cycle(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
    plot_lineage,
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
    plot_lineage(os.path.dirname(scored.output_folder))


def test_provenance_cartesian_into_panda_pool(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
    plot_lineage,
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
    plot_lineage(os.path.dirname(a.output_folder))


# ── three-hop: provenance chains through two children fan-outs ──────────────

def test_provenance_three_hop_chained_fanouts(
    local_config, isolated_cwd, new_pipeline, record_case, plot_lineage,
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
    plot_lineage(os.path.dirname(m1.output_folder))


# ── post-execution lineage: dropped IDs must vanish ──────────────────────────

def _regenerate_lineage(pipeline):
    """Rerun the lineage writer against the already-materialized map_tables."""
    pipeline._generate_id_lineage_csv()


def test_lineage_before_and_after_execution_missing_drops(
    local_config, isolated_cwd, new_pipeline, record_case, plot_lineage,
):
    """m1 → m2(children=<1..2>, missing=['s2']).

    Before execution (config time), the lineage CSV shows s2 as a parent
    because `missing` is only applied at runtime. After executing the pipe
    scripts and regenerating the lineage, the runtime map_tables show only
    the surviving IDs, so s2's branch must disappear from the lineage."""
    from biopipelines.mock import Mock

    pipeline = new_pipeline("prov_lineage_runtime_missing")
    with pipeline:
        m1 = Mock(
            ids=["s1", "s2", "s3"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
        )
        m2 = Mock(
            source=m1.streams.structures,
            children="<1..2>",
            missing=["s2"],
            streams={"designs": {"format": "pdb", "file": "<id>.pdb"}},
        )
        pipeline.save()

    # Before: config-time lineage still lists s2 (missing not yet applied).
    output_root = os.path.dirname(m1.output_folder)
    header, pre_rows = _read_lineage(output_root)
    pre_ids_m1 = sorted({r[0] for r in pre_rows if r[0]})
    plot_lineage(output_root, "pre-exec")

    # Execute both Mock stages to materialize the runtime map_tables.
    _run_pipe_mock(os.path.join(m1.output_folder, "_configuration", "mock_config.json"))
    _run_pipe_mock(os.path.join(m2.output_folder, "_configuration", "mock_config.json"))

    # Regenerate lineage from the now-present runtime map_tables.
    _regenerate_lineage(pipeline)
    plot_lineage(output_root, "post-exec")

    _, post_rows = _read_lineage(output_root)
    post_pairs = [(r[0], r[1]) for r in post_rows]
    # IDs surviving to m2's output (filter out the dash-marker rows).
    post_m2_survivors = sorted({r[1] for r in post_rows
                                if r[1] and r[1] != "-"})
    # Filtered rows: ID present in m1 column with '-' in m2 column.
    filtered_rows = [(r[0], r[1]) for r in post_rows if r[1] == "-"]

    record_case(
        input="m1→m2(missing=['s2']); lineage before vs after execution",
        expected=(
            ["s1", "s2", "s3"],                        # pre: all three
            ["s1_1", "s1_2", "s3_1", "s3_2"],          # post: s2's kids gone
            [("s2", "-")],                             # s2 marked filtered
        ),
        actual=(
            sorted({r[0] for r in pre_rows if r[0]}),
            post_m2_survivors,
            filtered_rows,
        ),
    )
    # Pre-execution: all three parents present.
    assert pre_ids_m1 == ["s1", "s2", "s3"]
    # Post-execution: s2 appears only as a filtered row (m2 dropped its kids).
    assert ("s2", "-") in post_pairs
    assert ("s2", "s2_1") not in post_pairs
    assert ("s2", "s2_2") not in post_pairs
    # Surviving children are correctly paired back to their parents.
    assert post_m2_survivors == ["s1_1", "s1_2", "s3_1", "s3_2"]
    assert ("s1", "s1_1") in post_pairs and ("s1", "s1_2") in post_pairs
    assert ("s3", "s3_1") in post_pairs and ("s3", "s3_2") in post_pairs


def test_lineage_after_execution_expands_compact_patterns(
    local_config, isolated_cwd, new_pipeline, record_case, plot_lineage,
):
    """Before execution the lineage shows compact IDs like 's1_<1..2>'; after
    execution it must list the concrete runtime IDs ('s1_1', 's1_2')."""
    from biopipelines.mock import Mock

    pipeline = new_pipeline("prov_lineage_runtime_expand")
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

    output_root = os.path.dirname(m1.output_folder)
    _, pre_rows = _read_lineage(output_root)
    pre_m2 = sorted({r[1] for r in pre_rows if r[1]})
    plot_lineage(output_root, "pre-exec")

    _run_pipe_mock(os.path.join(m1.output_folder, "_configuration", "mock_config.json"))
    _run_pipe_mock(os.path.join(m2.output_folder, "_configuration", "mock_config.json"))
    _regenerate_lineage(pipeline)
    plot_lineage(output_root, "post-exec")

    _, post_rows = _read_lineage(output_root)
    post_m2 = sorted({r[1] for r in post_rows if r[1]})
    post_pairs = sorted((r[0], r[1]) for r in post_rows)

    record_case(
        input="compact vs runtime-expanded lineage IDs",
        expected=(
            ["s1_<1..2>", "s2_<1..2>"],
            ["s1_1", "s1_2", "s2_1", "s2_2"],
        ),
        actual=(pre_m2, post_m2),
    )
    assert pre_m2 == ["s1_<1..2>", "s2_<1..2>"]
    assert post_m2 == ["s1_1", "s1_2", "s2_1", "s2_2"]
    assert post_pairs == [
        ("s1", "s1_1"), ("s1", "s1_2"),
        ("s2", "s2_1"), ("s2", "s2_2"),
    ]


# ── cycle tests: Mock → Panda → Mock (filter-driven cycle) ──────────────────

def test_cycle_mock_panda_mock(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
    plot_lineage,
):
    """Full Mock → Panda(filter, pool) → Mock cycle.

    The downstream Mock consumes Panda's pool-mode output stream. The
    config-time lineage PNG should show all three columns (source Mock,
    Panda, downstream Mock) wired left-to-right. No runtime execution —
    Panda's filter isn't actually run — but the pipeline script and
    lineage wiring are checked."""
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

    header, rows = _read_lineage(os.path.dirname(upstream.output_folder))
    record_case(
        input="Mock → Panda(filter, pool) → Mock(<1..2>)",
        expected=3,
        actual=len(header),
    )
    # Three tool columns in lineage (source Mock, Panda, downstream Mock).
    assert len(header) == 3
    # Every source-Mock ID appears in column 0 of some row.
    col0 = {r[0] for r in rows if r[0]}
    assert {"s1", "s2", "s3", "s4"}.issubset(col0)
    plot_lineage(os.path.dirname(upstream.output_folder))


def test_cycle_two_panda_cycles_chained(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
    plot_lineage,
):
    """Two chained filter cycles: Mock → Panda1 → Mock2 → Panda2 → Mock3.

    Simulates an iterative-refinement loop where each cycle filters the
    surviving IDs and fans out new children. Validates the lineage PNG can
    render a 5-tool chain and that the script wires them all."""
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
        Mock(
            source=p2.streams.structures,
            children="<A B>",
            streams={"variants": {"format": "pdb", "file": "<id>.pdb"}},
        )
        script_path = pipeline.save()

    assert_valid_script(
        script_path, "Mock", "Panda", "cycle_double_panda",
    )

    header, _ = _read_lineage(os.path.dirname(m0.output_folder))
    record_case(
        input="M0 → P1 → M1(<1..2>) → P2 → M2(<A B>)",
        expected=5,
        actual=len(header),
    )
    assert len(header) == 5
    plot_lineage(os.path.dirname(m0.output_folder))


def test_cycle_panda_filter_then_fanout_executes(
    local_config, isolated_cwd, new_pipeline, record_case, plot_lineage,
):
    """Cycle where Mock upstream + Mock downstream are actually executed.

    Panda itself is not executed (no pandas assumption); we only need the
    upstream and downstream Mock runtime map_tables to exist so the
    post-execution lineage renders with concrete IDs on both ends."""
    from biopipelines.mock import Mock
    from biopipelines.panda import Panda

    pipeline = new_pipeline("cycle_exec_mocks")
    with pipeline:
        upstream = Mock(
            ids=["s1", "s2"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            tables={"scores": {"columns": ["score"], "fill": {"score": 0.7}}},
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
        pipeline.save()

    output_root = os.path.dirname(upstream.output_folder)
    plot_lineage(output_root, "pre-exec")

    # Run upstream and downstream Mocks (Panda's filter is not executed).
    _run_pipe_mock(os.path.join(upstream.output_folder, "_configuration", "mock_config.json"))
    _run_pipe_mock(os.path.join(downstream.output_folder, "_configuration", "mock_config.json"))
    _regenerate_lineage(pipeline)
    plot_lineage(output_root, "post-exec")

    _, post_rows = _read_lineage(output_root)
    # After execution, downstream Mock's runtime map_table gives concrete
    # fan-out IDs in its column — either the last or next-to-last
    # depending on whether Panda's column sits in between.
    all_ids = {cell for r in post_rows for cell in r if cell and cell != "-"}
    record_case(
        input="Mock→Panda→Mock cycle, exec upstream+downstream Mocks",
        expected="{s1_1,s1_2,s2_1,s2_2} ⊂ lineage",
        actual=sorted(all_ids),
    )
    assert {"s1_1", "s1_2", "s2_1", "s2_2"}.issubset(all_ids)


# ── cycles with head / tail / sort / filter (rename path) ───────────────────

def test_cycle_sort_head_renames_pool_ids(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
    plot_lineage,
):
    """Mock → Panda(sort+head, pool) → Mock.

    sort/head/tail/sample trigger Panda's auto-rename path, so the pool's
    downstream IDs become "Panda_<N>_1", "Panda_<N>_2", ... The lineage
    must still chain through, and the downstream Mock consumes the
    renamed IDs as its parent axis."""
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
    plot_lineage(os.path.dirname(upstream.output_folder))


def test_cycle_tail_then_fanout(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
    plot_lineage,
):
    """Mock → Panda(sort+tail(2)) → Mock(<A B>).

    tail also triggers rename; this exercises the low-score / bottom-N
    cycle path and checks the lineage renders the full chain."""
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
    plot_lineage(os.path.dirname(m0.output_folder))


def test_cycle_filter_sort_head_combined(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
    plot_lineage,
):
    """Mock → Panda(filter + sort + head) → Mock.

    Combines all three ops in one Panda step. Filter alone wouldn't trigger
    rename; sort+head does, so the output IDs are renamed. The downstream
    Mock's parent axis is the Panda-renamed IDs."""
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
    plot_lineage(os.path.dirname(m0.output_folder))


def test_cycle_filter_only_preserves_original_ids(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
    plot_lineage,
):
    """Mock → Panda(filter, pool) → Mock.

    Filter alone does NOT trigger rename, so the pool output carries the
    original upstream IDs. The downstream Mock's parent axis therefore
    shows the surviving upstream IDs directly in the lineage."""
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
    plot_lineage(os.path.dirname(m0.output_folder))


def test_cycle_sample_renames_with_lazy_count(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
    plot_lineage,
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
    plot_lineage(os.path.dirname(m0.output_folder))


# ── post-exec lineage: Panda rename links back to upstream IDs ──────────────

def test_cycle_sort_head_post_exec_lineage_links_rename(
    local_config, isolated_cwd, new_pipeline, record_case, plot_lineage,
):
    """After executing Mock → Panda(sort+head) → Mock, the regenerated
    lineage must link each renamed Panda_N back to its original upstream ID.

    This is what fails without a runtime-written `<stream>_map.csv` in
    Panda's pool mode: the result table carries `original_id`, but the
    pipeline's lineage regeneration reads only stream map_tables, so the
    rename map is invisible. pipe_panda now writes `<stream>_map.csv`
    with the `<stream>.id` provenance column at runtime."""
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

    output_root = os.path.dirname(upstream.output_folder)
    plot_lineage(output_root, "pre-exec")

    _run_pipe_mock(os.path.join(upstream.output_folder, "_configuration", "mock_config.json"))
    _run_pipe_panda(os.path.join(ranked.output_folder, "_configuration", "panda_config.json"))

    # Confirm pipe_panda wrote the stream map_table with the provenance col.
    map_path = ranked.streams.s.map_table
    assert os.path.exists(map_path), f"Panda did not write {map_path}"
    rows = _read_map(map_path)
    panda_to_original = {r["id"]: r["s.id"] for r in rows}

    # Regenerate lineage from runtime map_tables.
    pipeline._generate_id_lineage_csv()
    plot_lineage(output_root, "post-exec")

    _, post_rows = _read_lineage(output_root)
    post_pairs = {(r[0], r[1]) for r in post_rows}

    record_case(
        input="Mock → Panda(sort+head, pool) → post-exec lineage",
        expected=("Panda_1 ← original upstream ID", "lineage pair present"),
        actual=(panda_to_original, sorted(post_pairs)),
    )
    # The two Panda-renamed IDs must each link to a concrete upstream ID
    # (one of a/b/c/d). With score=0.9 uniformly, pandas picks the first
    # two by sort stability, so "a" and "b" survive.
    assert set(panda_to_original.keys()) == {"Panda_1", "Panda_2"}
    assert set(panda_to_original.values()).issubset({"a", "b", "c", "d"})
    # Lineage rows connect upstream Mock column to Panda column.
    for upstream_id, panda_id in panda_to_original.items():
        # post_pairs is (mock_id, panda_id) — ordered by column.
        assert (panda_id, upstream_id) in post_pairs, (
            f"Expected lineage pair ({panda_id}, {upstream_id}) in {post_pairs}"
        )
