"""Tests for bp-visualize — a self-contained HTML page for one tool's output.

Covers:
- locating a step's ToolOutputs manifest, including through a Folder() nesting
- ordering by <table>.<column>, ascending and descending
- a bare (unqualified) sort key being refused with the qualified spelling named
- --max-items overriding the structure viewer's own 5-item sampling
- --ids and --streams selection, and streams/tables agreeing on one order
- the written view not affecting the step's completion verdict
"""

import csv
import json
import os

import pytest


def _fixture_step(pipeline, folder=None):
    """A saved Mock step with 8 structures and a scores table, and its folder.

    Returns ``(step_folder, ids)``. The files and CSVs are written by hand rather
    than by running the step: bp-visualize reads the manifest and the on-disk
    artifacts, which is exactly what a half-finished cluster job presents.
    """
    from biopipelines.mock import Mock
    from biopipelines.pipeline import Folder

    ids = [f"s{i}" for i in range(1, 9)]
    with pipeline:
        if folder:
            with Folder(folder):
                m = Mock(ids=ids,
                         streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
                         tables={"scores": {"columns": ["plddt"], "fill": {"plddt": 0.0}}})
        else:
            m = Mock(ids=ids,
                     streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
                     tables={"scores": {"columns": ["plddt"], "fill": {"plddt": 0.0}}})
        pipeline.save()

    step = m.output_folder
    sdir = os.path.join(step, "structures")
    tdir = os.path.join(step, "tables")
    os.makedirs(sdir, exist_ok=True)
    os.makedirs(tdir, exist_ok=True)

    rows = []
    for i, sid in enumerate(ids):
        path = os.path.join(sdir, f"{sid}.pdb")
        with open(path, "w") as f:
            f.write(f"ATOM      1  CA  ALA A   1     {i:7.3f}   0.000   0.000"
                    f"  1.00  0.00           C\nEND\n")
        rows.append((sid, path, 50 + i * 5))

    with open(os.path.join(sdir, "structures_map.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "file"])
        for sid, path, _ in rows:
            w.writerow([sid, path])

    with open(os.path.join(tdir, "scores.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "plddt"])
        for sid, _, score in rows:
            w.writerow([sid, score])

    return step, ids


def _viewer_ids(html):
    """The ids the 3D viewer actually embedded, in order."""
    import re
    match = re.search(r'\["s\d+"(?:,\s*"s\d+")*\]', html)
    return json.loads(match.group(0)) if match else []


def _table_ids(html):
    """Every id cell rendered in the page's tables, in order."""
    import re
    return re.findall(r"<td>(s\d+)</td>", html)


# ── locating the step ────────────────────────────────────────────────────────

def test_locates_manifest_for_a_plain_step(
    renderers_config, isolated_cwd, new_pipeline, record_case,
):
    from biopipelines.visualize import locate_tool_outputs_json

    step, _ids = _fixture_step(new_pipeline("viz_locate"))
    _job_root, manifest = locate_tool_outputs_json(step)

    record_case(input="locate_tool_outputs_json(<step>)",
                expected=os.path.basename(step) + ".json",
                actual=os.path.basename(manifest))
    assert os.path.isfile(manifest)
    assert os.path.basename(manifest) == os.path.basename(step) + ".json"


def test_locates_manifest_through_a_folder_nesting(
    renderers_config, isolated_cwd, new_pipeline, record_case,
):
    """A Folder() puts the step one level deeper; the walk up must still find ToolOutputs/."""
    from biopipelines.visualize import locate_tool_outputs_json

    step, _ids = _fixture_step(new_pipeline("viz_folder"), folder="group")
    job_root, manifest = locate_tool_outputs_json(step)

    record_case(input="locate through Folder('group')",
                expected="manifest found beside the job root",
                actual=os.path.relpath(manifest, job_root))
    assert os.path.isfile(manifest)
    assert "group" in step.replace("\\", "/")


def test_rejects_a_folder_with_no_job_above_it(
    renderers_config, isolated_cwd, record_case,
):
    from biopipelines.visualize import VisualizeError, locate_tool_outputs_json

    stray = os.path.join(str(isolated_cwd), "not_a_job")
    os.makedirs(stray)
    with pytest.raises(VisualizeError) as excinfo:
        locate_tool_outputs_json(stray)

    record_case(input="locate on a folder outside any job",
                expected="VisualizeError naming ToolOutputs",
                actual=str(excinfo.value)[:60])
    assert "ToolOutputs" in str(excinfo.value)


# ── ordering ─────────────────────────────────────────────────────────────────

def test_descending_sort_orders_streams_and_tables_alike(
    renderers_config, isolated_cwd, new_pipeline, record_case,
):
    from biopipelines.visualize import write_view

    step, _ids = _fixture_step(new_pipeline("viz_desc"))
    page = write_view(step, sort_spec="scores.plddt", descending=True, max_items=3)
    html = open(page, encoding="utf-8").read()

    record_case(input="--descending scores.plddt --max-items 3",
                expected=["s8", "s7", "s6"], actual=_viewer_ids(html))
    assert _viewer_ids(html) == ["s8", "s7", "s6"]
    # Every table on the page shows the same three ids in the same order.
    assert set(_table_ids(html)) == {"s6", "s7", "s8"}
    assert _table_ids(html)[:3] == ["s8", "s7", "s6"]


def test_ascending_sort_reverses_the_order(
    renderers_config, isolated_cwd, new_pipeline, record_case,
):
    from biopipelines.visualize import write_view

    step, _ids = _fixture_step(new_pipeline("viz_asc"))
    page = write_view(step, sort_spec="scores.plddt", descending=False, max_items=3)
    html = open(page, encoding="utf-8").read()

    record_case(input="--ascending scores.plddt --max-items 3",
                expected=["s1", "s2", "s3"], actual=_viewer_ids(html))
    assert _viewer_ids(html) == ["s1", "s2", "s3"]


def test_bare_sort_key_is_refused_and_names_the_qualified_form(
    renderers_config, isolated_cwd, new_pipeline, record_case,
):
    """A bare column name must not be guessed at — two tables can share one column."""
    from biopipelines.visualize import VisualizeError, load_step, parse_sort_key

    step, _ids = _fixture_step(new_pipeline("viz_bare"))
    output, _manifest, _root = load_step(step)

    with pytest.raises(VisualizeError) as excinfo:
        parse_sort_key("plddt", output)

    message = str(excinfo.value)
    record_case(input="parse_sort_key('plddt')",
                expected="refused, suggesting scores.plddt", actual=message)
    assert "scores.plddt" in message
    assert "<table>.<column>" in message


def test_unknown_table_and_column_are_named(
    renderers_config, isolated_cwd, new_pipeline, record_case,
):
    from biopipelines.visualize import VisualizeError, load_step, ordered_ids, parse_sort_key

    step, _ids = _fixture_step(new_pipeline("viz_unknown"))
    output, _manifest, _root = load_step(step)

    with pytest.raises(VisualizeError) as bad_table:
        parse_sort_key("nosuch.plddt", output)
    with pytest.raises(VisualizeError) as bad_column:
        ordered_ids(output, "scores", "nope", True)

    record_case(input="unknown table / unknown column",
                expected="both refused",
                actual=(str(bad_table.value)[:40], str(bad_column.value)[:40]))
    assert "nosuch" in str(bad_table.value)
    assert "nope" in str(bad_column.value)


# ── selection ────────────────────────────────────────────────────────────────

def test_max_items_overrides_the_viewer_sampling(
    renderers_config, isolated_cwd, new_pipeline, record_case,
):
    """The structure viewer embeds 5 by default; an explicit cap must win."""
    from biopipelines.renderers import structures as structures_renderer
    from biopipelines.visualize import write_view

    step, _ids = _fixture_step(new_pipeline("viz_cap"))
    page = write_view(step, sort_spec="scores.plddt", descending=True, max_items=7)
    html = open(page, encoding="utf-8").read()

    record_case(input=f"--max-items 7 (default budget {structures_renderer.MAX_EMBEDDED})",
                expected=7, actual=len(_viewer_ids(html)))
    assert structures_renderer.MAX_EMBEDDED == 5
    assert len(_viewer_ids(html)) == 7


def test_no_selection_leaves_the_default_sampling_alone(
    renderers_config, isolated_cwd, new_pipeline, record_case,
):
    from biopipelines.renderers import structures as structures_renderer
    from biopipelines.visualize import write_view

    step, _ids = _fixture_step(new_pipeline("viz_default"))
    page = write_view(step)
    html = open(page, encoding="utf-8").read()

    record_case(input="no options (8 items, default budget)",
                expected=structures_renderer.MAX_EMBEDDED,
                actual=len(_viewer_ids(html)))
    assert len(_viewer_ids(html)) == structures_renderer.MAX_EMBEDDED


def test_cli_default_max_items_is_ten(record_case):
    """bp-visualize with no --max-items embeds 10.

    5 is the run page's per-node sampling budget and hides exactly one item of a
    6-item stream, which is the least useful place to stop. The library function keeps
    0 (no cap); the opinion lives in the command.
    """
    from biopipelines.visualize import DEFAULT_MAX_ITEMS, _build_parser

    parsed = _build_parser().parse_args(["/some/step"])
    record_case(input="bp-visualize <step> (no --max-items)",
                expected=10, actual=parsed.max_items)
    assert DEFAULT_MAX_ITEMS == 10
    assert parsed.max_items == 10
    assert _build_parser().parse_args(["/s", "--max-items", "0"]).max_items == 0


def test_ids_and_sort_together_use_the_sorted_order(
    renderers_config, isolated_cwd, new_pipeline, record_case,
):
    """With both, the order is the sort restricted to the ids — not the typed order."""
    from biopipelines.visualize import write_view

    step, _ids = _fixture_step(new_pipeline("viz_ids"))
    page = write_view(step, ids=["s3", "s5", "s1"],
                      sort_spec="scores.plddt", descending=True)
    html = open(page, encoding="utf-8").read()

    record_case(input="--ids s3,s5,s1 --descending scores.plddt",
                expected=["s5", "s3", "s1"], actual=_viewer_ids(html))
    assert _viewer_ids(html) == ["s5", "s3", "s1"]
    assert _table_ids(html)[:3] == ["s5", "s3", "s1"]


def test_cap_spares_a_stream_the_sort_does_not_rank(
    renderers_config, isolated_cwd, new_pipeline, record_case,
):
    """--max-items alongside a sort means the top N of what was ranked.

    Measured on a real Boltz2 step: ordered by confidence.confidence_score, it cut its
    6 sequences and 6 MSAs to 5 alongside the intended 30 structures to 5. The
    confidence table keys structure ids, so those six were dropped on no criterion.
    The streams are built directly here because Mock gives every stream the same ids,
    which is exactly the case where the distinction does not arise.
    """
    from biopipelines.datastream import DataStream
    from biopipelines.outputs import StandardizedOutput
    from biopipelines.visualize import apply_selection

    step, _ids = _fixture_step(new_pipeline("viz_unranked"))
    scores = os.path.join(step, "tables", "scores.csv")

    structures = DataStream(
        name="structures", format="pdb",
        ids=[f"s{i}" for i in range(1, 9)],
        files=[os.path.join(step, "structures", f"s{i}.pdb") for i in range(1, 9)],
    )
    # Protein-level ids, as a real Boltz2 sequences stream carries — none of them
    # appear in the structure-keyed scores table.
    sequences = DataStream(name="sequences", format="csv",
                           ids=[f"prot{i}" for i in range(1, 7)], files=[])

    from biopipelines.base_config import TableInfo
    output = StandardizedOutput({
        "structures": structures,
        "sequences": sequences,
        "tables": {"scores": TableInfo(name="scores", path=scores,
                                       columns=["id", "plddt"], description="")},
        "output_folder": step,
    })

    narrowed, notes = apply_selection(output, sort=("scores", "plddt"),
                                      descending=True, max_items=3)

    ranked = len(narrowed.streams.structures.ids_expanded)
    unranked = len(narrowed.streams.sequences.ids_expanded)
    record_case(input="--descending scores.plddt --max-items 3; sequences unranked",
                expected=(3, 6), actual=(ranked, unranked))
    assert ranked == 3, "the ranked stream should be capped"
    assert unranked == 6, "a stream the sort does not rank must be shown whole"
    assert any("does not rank it" in n for n in notes), notes


def test_cap_without_a_sort_applies_to_every_stream(
    renderers_config, isolated_cwd, new_pipeline, record_case,
):
    """With no ordering, --max-items is purely "keep the page small" and caps all."""
    from biopipelines.visualize import apply_selection, load_step

    step, _ids = _fixture_step(new_pipeline("viz_cap_nosort"))
    output, _manifest, _root = load_step(step)
    narrowed, _notes = apply_selection(output, max_items=2)

    count = len(narrowed.streams.structures.ids_expanded)
    record_case(input="--max-items 2 with no sort", expected=2, actual=count)
    assert count == 2


def test_streams_filter_drops_the_others(
    renderers_config, isolated_cwd, new_pipeline, record_case,
):
    from biopipelines.visualize import apply_selection, load_step

    step, _ids = _fixture_step(new_pipeline("viz_streams"))
    output, _manifest, _root = load_step(step)
    narrowed, _notes = apply_selection(output, streams=["structures"])

    names = sorted(narrowed.streams.keys())
    record_case(input="--streams structures", expected=["structures"], actual=names)
    assert names == ["structures"]


def test_page_is_self_contained(
    renderers_config, isolated_cwd, new_pipeline, record_case,
):
    """No external <script src> — the page must open from file:// off a cluster."""
    import re

    from biopipelines.visualize import write_view

    step, _ids = _fixture_step(new_pipeline("viz_selfcontained"))
    page = write_view(step, sort_spec="scores.plddt", descending=True, max_items=2)
    html = open(page, encoding="utf-8").read()

    external = re.findall(r'<(?:script|link)[^>]+(?:src|href)="https?://[^"]+"', html)
    record_case(input="write_view(...) external references",
                expected=[], actual=external)
    assert external == []


# ── the view must not disturb the step ───────────────────────────────────────

def test_page_styles_the_classes_its_renderers_emit(
    renderers_config, isolated_cwd, new_pipeline, record_case,
):
    """A page that emits bp-table but never styles it renders naked tables.

    pipeline_report._CSS covers the page chrome only; the renderer fragment classes
    are styled by StandardizedOutput._CSS, which render_page emits as a second style
    block. A shell carrying only the first looked nothing like the run page or the
    notebook view even though the fragments were byte-identical.
    """
    import re

    from biopipelines.visualize import write_view

    step, _ids = _fixture_step(new_pipeline("viz_css"))
    page = write_view(step, sort_spec="scores.plddt", descending=True, max_items=2)
    html = open(page, encoding="utf-8").read()

    emitted = set(re.findall(r'class="(bp-[a-z-]+)"', html))
    styled = {m for m in re.findall(r'\.(bp-[a-z-]+)\s*[{,: \[]', html)}
    unstyled = sorted(emitted - styled)

    record_case(input="fragment classes emitted vs styled",
                expected=[], actual=unstyled)
    assert emitted, "the page emitted no bp-* fragment classes at all"
    assert unstyled == [], f"emitted but never styled: {unstyled}"


def test_default_output_lands_in_extras(
    renderers_config, isolated_cwd, new_pipeline, record_case,
):
    from biopipelines.visualize import write_view

    step, _ids = _fixture_step(new_pipeline("viz_extras"))
    page = write_view(step)

    rel = os.path.relpath(page, step).replace("\\", "/")
    record_case(input="write_view(step) default path",
                expected="_extras/<step>_view.html", actual=rel)
    assert rel == f"_extras/{os.path.basename(step)}_view.html"


def test_view_in_extras_does_not_change_the_completion_verdict(
    renderers_config, isolated_cwd, new_pipeline, record_case,
):
    """The completion check is declaration-driven, so an extra file cannot fail a step.

    It reads the declared paths out of .expected_outputs.json and never enumerates a
    directory, which is what makes writing a view into _extras/ safe. Pinned here
    because "put the page beside the outputs" is only a good default while that holds.
    """
    import sys

    sys.path.insert(0, os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "pipe_scripts"))
    import pipe_check_completion

    from biopipelines.visualize import write_view

    step, _ids = _fixture_step(new_pipeline("viz_verdict"))
    expected_json = os.path.join(step, ".expected_outputs.json")
    assert os.path.isfile(expected_json), "fixture step wrote no expected-outputs manifest"
    with open(expected_json, encoding="utf-8") as f:
        expected = json.load(f)["output_structure"]

    before, _missing_before = pipe_check_completion.check_expected_outputs(
        expected, "Mock", step)
    page = write_view(step)
    after, missing_after = pipe_check_completion.check_expected_outputs(
        expected, "Mock", step)

    record_case(input="completion verdict before/after writing _extras/<step>_view.html",
                expected=(before, before), actual=(before, after))
    assert os.path.isfile(page)
    assert after == before, f"the view changed the verdict; missing={missing_after}"
