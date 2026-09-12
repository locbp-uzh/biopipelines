"""Tests for dataflow edge recovery and the per-run HTML pipeline page.

Covers:
- the producer back-reference `ToolOutput.output` stamps on outputs and streams
- which input shapes `recover_input_edges` recovers, and the one it cannot
- central population of the long-unused `input_sources`
- `pipeline.save()` writing RunTime/pipeline_graph.json + RunTime/pipeline.html
- the page being self-contained (no CDN / absolute URL fetches) and well formed
- the browser-library chain behind the 3D viewers: vendored copy inlined once, else the CDN tag under allow_external, else the metadata table
- regeneration from the saved graph, and honest degradation when an output
  inventory or a renderer is missing
"""

import json
import os
import re
from html.parser import HTMLParser

import pytest


_EXTERNAL_REF_RE = re.compile(r"""(?:src|href)\s*=\s*["']?\s*(?:https?:)?//""", re.I)

_VOID_TAGS = {
    "area", "base", "br", "col", "embed", "hr", "img", "input", "link",
    "meta", "param", "source", "track", "wbr",
}

# The fixture config declares no `renderers:` section, so tests needing real renderers pass this.
_REAL_RENDERERS = {
    "streams": {
        "_metadata": "renderers/streams.py",
        "pdb": "renderers/structures.py",
        "fasta": "renderers/fasta.py",
        "_default": "renderers/streams.py",
    },
    "tables": {"_default": "renderers/tables.py"},
}


class _TagBalance(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.stack = []
        self.errors = []

    def handle_starttag(self, tag, attrs):
        if tag not in _VOID_TAGS:
            self.stack.append(tag)

    def handle_endtag(self, tag):
        if tag in _VOID_TAGS:
            return
        if not self.stack:
            self.errors.append(f"stray </{tag}>")
            return
        if self.stack[-1] != tag:
            self.errors.append(f"</{tag}> closes <{self.stack[-1]}>")
            if tag in self.stack:
                del self.stack[self.stack.index(tag):]
            return
        self.stack.pop()


def _page_renderer():
    from biopipelines.pipeline import _load_page_renderer

    return _load_page_renderer()


def _mock_stream(name="structures", fmt="pdb"):
    return {name: {"format": fmt, "file": "<id>.pdb"}}


# ── the producer back-reference ───────────────────────────────────────────────

def test_output_and_streams_carry_their_producer(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    from biopipelines.mock import Mock

    pipeline = new_pipeline("page_producer")
    with pipeline:
        out = Mock(ids=["a"], streams=_mock_stream())
        producing_tool = pipeline.tools[0]
        record_case(
            input="Mock(ids=[a]).producer",
            expected=("Mock", "Mock"),
            actual=(out.producer.TOOL_NAME, out.streams.structures._producer.TOOL_NAME),
        )
        assert out.producer is producing_tool
        assert out.streams.structures._producer is producing_tool
        # Stamping the container is what would let the framework say `.streams.typo` matched nothing.
        assert out.streams._producer is producing_tool


def test_derived_outputs_inherit_the_producer(
    local_config, isolated_cwd, new_pipeline,
):
    from biopipelines.mock import Mock

    pipeline = new_pipeline("page_derived")
    with pipeline:
        out = Mock(ids=["a", "b"], streams=_mock_stream())
        producing_tool = pipeline.tools[0]

        assert out["a"].producer is producing_tool
        assert next(iter(out)).producer is producing_tool
        assert out.chunks(2)[0].producer is producing_tool


# ── which input shapes are recoverable ───────────────────────────────────────

def test_recovers_every_supported_input_shape(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    from biopipelines.combinatorics import Bundle, Each
    from biopipelines.mock import Mock
    from biopipelines.pool import Pool

    pipeline = new_pipeline("page_shapes")
    with pipeline:
        a = Mock(ids=["x1", "x2"], streams=_mock_stream())
        b = Mock(ids=["y1", "y2"], streams=_mock_stream())
        Mock(source=a.streams.structures, streams=_mock_stream("designs"))
        Mock(source=Bundle(a.streams.structures), streams=_mock_stream("designs"))
        Mock(source=[Each(a.streams.structures), Each(b.streams.structures)],
             streams=_mock_stream("designs"))
        Mock(source=next(iter(a)).streams.structures, streams=_mock_stream("designs"))
        Mock(source=a["x1"].streams.structures, streams=_mock_stream("designs"))
        Pool(runs=[a, b])

    edges = {(e["from_step"], e["to_step"]) for e in pipeline.dataflow_edges}
    expected = {(1, 3), (1, 4), (1, 5), (2, 5), (1, 6), (1, 7), (1, 8), (2, 8)}
    record_case(input="Mock/Bundle/Each/list/iter/getitem/Pool inputs",
                expected=sorted(expected), actual=sorted(edges))
    assert edges == expected

    pool_edges = [e for e in pipeline.dataflow_edges if e["to_step"] == 8]
    assert all(e["stream"] is None for e in pool_edges)
    assert all(e["shape"] == "StandardizedOutput" for e in pool_edges)
    assert all(e["arguments"] == ["runs"] for e in pool_edges)

    bare_stream_edge = next(e for e in pipeline.dataflow_edges if e["to_step"] == 3)
    assert bare_stream_edge["stream"] == "structures"
    assert bare_stream_edge["shape"] == "DataStream"


def test_datastream_derived_chunk_is_recovered(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """A stream derived inside DataStream keeps its producer, so the edge stands.

    This was pinned as a known gap: ``chunks`` / slicing / iterating a bare
    stream built a fresh DataStream without copying ``_producer``, so no edge
    was recovered and the page said in as many words that the two steps were
    "not a dependency" -- an affirmative false claim about a wired pair.
    """
    from biopipelines.mock import Mock

    pipeline = new_pipeline("page_chunk_gap")
    with pipeline:
        a = Mock(ids=["x1", "x2"], streams=_mock_stream())
        Mock(source=a.streams.structures.chunks(2)[0], streams=_mock_stream("designs"))

    edges = [e for e in pipeline.dataflow_edges if e["to_step"] == 2]
    record_case(input="Mock(source=stream.chunks(2)[0])",
                expected="one edge from step 1", actual=edges)
    assert len(edges) == 1, edges
    assert edges[0]["from_step"] == 1
    assert edges[0]["stream"] == "structures"


def test_a_filtered_stream_is_still_a_recovered_dependency(
    local_config, isolated_cwd, new_pipeline,
):
    """The commonest shape: a consumer taking a narrowed stream."""
    from biopipelines.mock import Mock

    pipeline = new_pipeline("page_filter_edge")
    with pipeline:
        a = Mock(ids=["x1", "x2", "x3"], streams=_mock_stream())
        Mock(source=a.streams.structures.filter_by_ids(["x1", "x3"]),
             streams=_mock_stream("designs"))

    edges = [e for e in pipeline.dataflow_edges if e["to_step"] == 2]
    assert len(edges) == 1, edges
    assert edges[0]["from_step"] == 1


def test_input_sources_is_populated_centrally(
    local_config, isolated_cwd, new_pipeline,
):
    """`input_sources` has existed on every tool since forever and 2 of 92 filled it."""
    from biopipelines.mock import Mock

    pipeline = new_pipeline("page_input_sources")
    with pipeline:
        a = Mock(ids=["x"], streams=_mock_stream())
        Mock(source=a.streams.structures, streams=_mock_stream("designs"))

    consumer = pipeline.tools[1]
    assert consumer.input_sources, "consumer recorded no input sources"
    record = next(iter(consumer.input_sources.values()))
    assert record["from_step"] == 1
    assert record["from_tool"] == "Mock"
    assert record["stream"] == "structures"


def test_a_tool_never_records_an_edge_to_itself(
    local_config, isolated_cwd, new_pipeline,
):
    from biopipelines.mock import Mock

    pipeline = new_pipeline("page_no_self_edge")
    with pipeline:
        Mock(ids=["x"], streams=_mock_stream())

    assert pipeline.dataflow_edges == []


# ── the graph record ─────────────────────────────────────────────────────────

def test_graph_records_batches_steps_and_edges(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    from biopipelines.mock import Mock
    from biopipelines.pipeline import Resources

    pipeline = new_pipeline("page_graph")
    with pipeline:
        Resources(cpus=1, memory="2GB")
        a = Mock(ids=["x"], streams=_mock_stream())
        Resources(cpus=4, memory="8GB")
        Mock(source=a.streams.structures, streams=_mock_stream("designs"))
        pipeline.save()

    graph_path = os.path.join(pipeline.folders["runtime"], "pipeline_graph.json")
    assert os.path.isfile(graph_path)
    graph = json.loads(open(graph_path, encoding="utf-8").read())

    record_case(input="two Resources() calls, two Mocks",
                expected=(2, 2, 1),
                actual=(len(graph["batches"]), len(graph["steps"]), len(graph["edges"])))
    assert [b["index"] for b in graph["batches"]] == [0, 1]
    assert graph["batches"][1]["parents"] == [0]
    assert [s["batch"] for s in graph["steps"]] == [0, 1]
    assert all(s["tool_version"] for s in graph["steps"])
    assert all(os.path.isfile(s["tool_outputs_json"]) for s in graph["steps"])
    for step in graph["steps"]:
        assert step["completed_marker"].endswith("_COMPLETED")
        assert step["script_file"].endswith(f"{step['script_basename']}.sh")


# ── the page ─────────────────────────────────────────────────────────────────

def test_save_writes_a_self_contained_page(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    from biopipelines.mock import Mock
    from biopipelines.pipeline import Resources

    pipeline = new_pipeline("page_selfcontained")
    with pipeline:
        Resources(cpus=1)
        a = Mock(ids=["x", "y"], streams=_mock_stream())
        Resources(cpus=2)
        Mock(source=a.streams.structures, streams=_mock_stream("designs"))
        pipeline.save()

    page_path = os.path.join(pipeline.folders["runtime"], "pipeline.html")
    assert os.path.isfile(page_path)
    page = open(page_path, encoding="utf-8").read()

    external = _EXTERNAL_REF_RE.findall(page)
    record_case(input="saved pipeline page", expected=[], actual=external)
    assert external == [], f"page fetches external resources: {external}"

    balance = _TagBalance()
    balance.feed(page)
    assert balance.errors == []
    assert balance.stack == []

    assert page.startswith("<!doctype html>")
    assert "no step has run" in page, "a page written before the run must say so"
    assert page.count('<details class="node"') == 2
    assert '<section class="lane"' in page
    assert "Batch 0" in page and "Batch 1" in page
    # Only recovered dataflow is drawn now: execution order is the order the steps are listed in, and
    # a dashed arrow for it was being read as a dependency.
    assert "recovered dataflow" in page
    assert "order only" not in page


def test_page_states_recovered_and_order_only_edges(
    local_config, isolated_cwd, new_pipeline,
):
    from biopipelines.mock import Mock

    pipeline = new_pipeline("page_edge_labels")
    with pipeline:
        a = Mock(ids=["x"], streams=_mock_stream())
        Mock(ids=["unrelated"], streams=_mock_stream("other"))
        Mock(source=a.streams.structures, streams=_mock_stream("designs"))
        pipeline.save()

    page = open(os.path.join(pipeline.folders["runtime"], "pipeline.html"),
                encoding="utf-8").read()
    wires = json.loads(
        re.search(r'<script id="edge-data" type="application/json">(.*?)</script>',
                  page, re.S).group(1)
    )
    flow = {(w["from"], w["to"]) for w in wires if w["kind"] == "flow"}
    order = {(w["from"], w["to"]) for w in wires if w["kind"] == "order"}

    assert flow == {(1, 3)}
    # Order edges skip pairs with a recovered edge, so a dependency is never overdrawn.
    assert order == {(1, 2), (2, 3)}


def test_regeneration_reproduces_the_page_from_the_saved_graph(
    local_config, isolated_cwd, new_pipeline,
):
    from biopipelines.mock import Mock
    from biopipelines.pipeline import regenerate_pipeline_page

    pipeline = new_pipeline("page_regenerate")
    with pipeline:
        a = Mock(ids=["x"], streams=_mock_stream())
        Mock(source=a.streams.structures, streams=_mock_stream("designs"))
        pipeline.save()

    runtime = pipeline.folders["runtime"]
    page_path = os.path.join(runtime, "pipeline.html")
    first = open(page_path, encoding="utf-8").read()
    assert regenerate_pipeline_page(runtime) == page_path
    second = open(page_path, encoding="utf-8").read()

    strip = lambda text: re.sub(r"Generated \S+ by", "Generated by", text)
    assert strip(first) == strip(second)


def test_regeneration_needs_the_graph(local_config, isolated_cwd, tmp_path):
    from biopipelines.pipeline import regenerate_pipeline_page

    with pytest.raises(FileNotFoundError, match="pipeline_graph.json"):
        regenerate_pipeline_page(str(tmp_path))


def test_missing_output_inventory_degrades_to_a_named_placeholder(
    local_config, isolated_cwd, new_pipeline,
):
    from biopipelines.mock import Mock
    from biopipelines.pipeline import regenerate_pipeline_page

    pipeline = new_pipeline("page_missing_inventory")
    with pipeline:
        Mock(ids=["x"], streams=_mock_stream())
        pipeline.save()

    runtime = pipeline.folders["runtime"]
    graph = json.loads(
        open(os.path.join(runtime, "pipeline_graph.json"), encoding="utf-8").read()
    )
    inventory = graph["steps"][0]["tool_outputs_json"]
    os.remove(inventory)

    page = open(regenerate_pipeline_page(runtime), encoding="utf-8").read()
    assert "no exported output inventory" in page
    assert '<details class="node"' in page, "the node must still render"


def test_the_on_the_fly_step_exports_the_inventory_before_refreshing():
    """Pins the call order in _execute_tool_on_the_fly: exporting after the page is written would leave the first refresh showing the placeholder."""
    import inspect

    from biopipelines.pipeline import Pipeline

    source = inspect.getsource(Pipeline._execute_tool_on_the_fly)
    assert "_export_tool_outputs" in source, "the on-the-fly step must export its inventory"
    assert source.index("_export_tool_outputs") < source.index("_write_pipeline_page")


def _pipeline_with_structures_on_disk(new_pipeline, monkeypatch, name, ids=("x",)):
    """A saved pipeline whose pdb stream has real files, so renderers/structures.py emits its viewer."""
    from biopipelines.config_manager import ConfigManager
    from biopipelines.mock import Mock

    monkeypatch.setattr(ConfigManager, "get_renderers_config",
                        lambda self: _REAL_RENDERERS)

    pipeline = new_pipeline(name)
    with pipeline:
        outs = [Mock(ids=[i], streams=_mock_stream(),
                     tables={"scores": {"columns": ["score"], "fill": {"score": 1.0}}},
                     map_table_strategy="config")
                for i in ids]
        pipeline.save()

    for out in outs:
        for path in out.streams.structures.files_expanded:
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, "w") as f:
                f.write("ATOM      1  CA  ALA A   1       0.000   0.000   0.000\nEND\n")

    return pipeline


def _render_without_vendored_libs(pipeline, monkeypatch, tmp_path, allow_external):
    """Render the saved graph with the vendored library directory pointed at an empty folder."""
    module = _page_renderer()
    empty = tmp_path / "no_vendor"
    empty.mkdir(exist_ok=True)
    monkeypatch.setattr(module, "_VENDOR_DIR", str(empty))

    graph_path = os.path.join(pipeline.folders["runtime"], "pipeline_graph.json")
    with open(graph_path, encoding="utf-8") as f:
        graph = json.load(f)
    out_path = str(tmp_path / f"page_external_{allow_external}.html")
    module.render_page(graph, out_path, allow_external=allow_external)
    return open(out_path, encoding="utf-8").read()


def test_the_vendored_library_is_inlined_once_and_keeps_the_page_offline(
    local_config, isolated_cwd, new_pipeline, monkeypatch, record_case,
):
    """Link 1 of the chain: a vendored copy makes the viewers work with no network, embedded once for the whole page."""
    from biopipelines.pipeline import regenerate_pipeline_page

    pipeline = _pipeline_with_structures_on_disk(
        new_pipeline, monkeypatch, "page_vendored_lib", ids=("x", "y", "z"))
    page = open(regenerate_pipeline_page(pipeline.folders["runtime"]), encoding="utf-8").read()

    embeds = page.count('<script data-bp-lib-embedded="3dmol">')
    viewers = len(re.findall(r'id="bp3d_[0-9a-f]+_viewer"', page))
    record_case(input="3 pdb streams with renderers/structures.py configured",
                expected=(1, 3, [], 0),
                actual=(embeds, viewers, _EXTERNAL_REF_RE.findall(page),
                        page.count("cdn.jsdelivr")))
    assert viewers == 3, "each structure node should still get its own viewer"
    assert embeds == 1, "the library must be embedded once per page, not once per node"
    assert _EXTERNAL_REF_RE.findall(page) == []
    assert "cdn.jsdelivr" not in page
    assert "needs a CDN" not in page
    assert "ATOM      1  CA  ALA" in page, "coordinates must still be embedded"
    assert "embedded from biopipelines/renderers/vendor/" in page


def test_without_a_vendored_copy_allow_external_keeps_the_cdn_tag(
    local_config, isolated_cwd, new_pipeline, monkeypatch, tmp_path, record_case,
):
    """Link 2 of the chain: no vendored copy but the opt-in flag, so the remote tag is shipped."""
    pipeline = _pipeline_with_structures_on_disk(
        new_pipeline, monkeypatch, "page_cdn_fallback")
    page = _render_without_vendored_libs(pipeline, monkeypatch, tmp_path, True)

    record_case(input="no vendored library, allow_external=True",
                expected=(1, 0, True),
                actual=(page.count("cdn.jsdelivr"),
                        page.count('data-bp-lib-embedded'),
                        "this page needs a network" in page))
    assert page.count("cdn.jsdelivr") == 1
    assert "data-bp-lib-embedded" not in page
    assert "this page needs a network" in page


def test_without_a_vendored_copy_or_the_flag_the_viewer_falls_back_to_metadata(
    local_config, isolated_cwd, new_pipeline, monkeypatch, tmp_path, record_case,
):
    """Link 3 of the chain: the pre-existing fallback, which still names what it dropped."""
    pipeline = _pipeline_with_structures_on_disk(
        new_pipeline, monkeypatch, "page_table_fallback")
    page = _render_without_vendored_libs(pipeline, monkeypatch, tmp_path, False)

    record_case(input="no vendored library, allow_external=False",
                expected=(True, [], True),
                actual=("needs a CDN" in page,
                        _EXTERNAL_REF_RE.findall(page),
                        "3Dmol" not in page))
    assert "needs a CDN" in page
    assert _EXTERNAL_REF_RE.findall(page) == []
    assert "3Dmol" not in page
    # The self-contained metadata renderer still ran, so the node is not empty.
    assert "structures_map.csv" in page


def test_the_vendored_library_matches_what_the_renderers_would_fetch(record_case):
    """The vendored file and the CDN URL the renderers fall back to must name one version."""
    import hashlib

    module = _page_renderer()
    path = module._vendored_lib("3dmol")
    record_case(input="renderers/vendor/3Dmol-min.js",
                expected="present, 2.5.2, sha256 7b26bfd8…",
                actual=path)
    assert path, "the vendored 3Dmol copy is missing; pages would need a network"

    blob = open(path, "rb").read()
    assert len(blob) == 524_222
    assert hashlib.sha256(blob).hexdigest() == (
        "7b26bfd8170372ea78ea06df4a45e4a55fa5f538c2dfd716166fedf31ebfce9a")
    # An inline <script> would be truncated by this substring, so page assembly escapes it; catch the need early.
    assert b"</script" not in blob.lower()

    for renderer in ("structures.py", "grids.py"):
        source = open(os.path.join(os.path.dirname(module.__file__), renderer),
                      encoding="utf-8").read()
        assert "3dmol@2.5.2/build/3Dmol-min.js" in source, (
            f"{renderer} falls back to a different 3Dmol version than the vendored copy")


def test_page_reuses_the_configured_renderers_for_node_bodies(
    local_config, isolated_cwd, new_pipeline, monkeypatch,
):
    from biopipelines.config_manager import ConfigManager
    from biopipelines.mock import Mock
    from biopipelines.pipeline import regenerate_pipeline_page

    monkeypatch.setattr(ConfigManager, "get_renderers_config",
                        lambda self: _REAL_RENDERERS)

    pipeline = new_pipeline("page_reuse_renderers")
    with pipeline:
        Mock(ids=["s1", "s2"], streams={"sequences": {"format": "csv",
                                                      "values": ["MKT", "AET"]}},
             tables={"scores": {"columns": ["score"], "fill": {"score": 2.5}}},
             map_table_strategy="config")
        pipeline.save()

    page = open(regenerate_pipeline_page(pipeline.folders["runtime"]),
                encoding="utf-8").read()
    # renderers/streams.py and renderers/tables.py markup, not a reimplementation.
    assert 'class="bp-table"' in page
    assert 'class="bp-section-title"' in page
    assert "s1" in page and "s2" in page


def test_page_names_the_gap_when_no_renderers_are_configured(
    local_config, isolated_cwd, new_pipeline,
):
    """The fixture config declares no `renderers:` section; the page must say so."""
    from biopipelines.mock import Mock

    pipeline = new_pipeline("page_no_renderers")
    with pipeline:
        Mock(ids=["x"], streams=_mock_stream())
        pipeline.save()

    page = open(os.path.join(pipeline.folders["runtime"], "pipeline.html"),
                encoding="utf-8").read()
    assert "no stream renderers are configured" in page


def test_completed_and_failed_markers_reach_the_page(
    local_config, isolated_cwd, new_pipeline,
):
    from biopipelines.mock import Mock
    from biopipelines.pipeline import regenerate_pipeline_page

    pipeline = new_pipeline("page_status")
    with pipeline:
        Mock(ids=["x"], streams=_mock_stream())
        Mock(ids=["y"], streams=_mock_stream("other"))
        pipeline.save()

    runtime = pipeline.folders["runtime"]
    graph = json.loads(
        open(os.path.join(runtime, "pipeline_graph.json"), encoding="utf-8").read()
    )
    open(graph["steps"][0]["completed_marker"], "w").close()
    open(graph["steps"][1]["failed_marker"], "w").close()

    page = open(regenerate_pipeline_page(runtime), encoding="utf-8").read()
    assert '<span class="pill completed">' in page
    assert '<span class="pill failed">' in page
    assert "left a FAILED marker" in page
    assert "no step has run" not in page


# ── file contents reach an inline <script>, so they cannot be trusted as text ──

def _load_renderer(rel_path):
    """Resolve the way production does: against the package, not the cwd."""
    import importlib.util
    import os
    import biopipelines
    from biopipelines.pipeline import _load_page_renderer  # noqa: F401  (proves the loader exists)
    abs_path = os.path.join(os.path.dirname(biopipelines.__file__), rel_path)
    spec = importlib.util.spec_from_file_location(rel_path.replace("/", "_"), abs_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_BREAKOUT = 'REMARK 1 </script><script>window.__PWNED=1</script>'


@pytest.mark.parametrize("rel_path", ["renderers/structures.py", "renderers/grids.py"])
def test_a_structure_file_cannot_close_the_viewer_script(rel_path):
    r"""Both viewers embed file contents (PDB/CIF/DX) into an inline `<script>` via `json.dumps`, which escapes quotes and backslashes but not `</script`. A REMARK carrying that sequence closed the script early and whatever followed ran as markup on open.

    `pipeline_report._library_scripts` already neutralizes the same sequence in the vendored library it inlines; this is the same rule applied to the data beside it. `<\/` is valid JSON for `/`, so the payload still round-trips as text in the viewer.
    """
    module = _load_renderer(rel_path)
    encoded = module._script_json([_BREAKOUT])
    assert "</script" not in encoded
    assert r"<\/script" in encoded
    assert json.loads(encoded) == [_BREAKOUT], "escaping must not change the data the viewer reads"


@pytest.mark.parametrize("rel_path", ["renderers/structures.py", "renderers/grids.py"])
def test_the_viewers_route_every_inline_payload_through_the_escaper(rel_path):
    """A future payload added with a bare `json.dumps` would reopen the hole silently, so no bare call survives outside the helper itself."""
    import biopipelines
    abs_path = os.path.join(os.path.dirname(biopipelines.__file__), rel_path)
    source = open(abs_path, encoding="utf-8").read()
    body = source.split("def _script_json", 1)[1].split("\n\n\n", 1)[1]
    assert "json.dumps" not in body, (
        f"{rel_path} embeds a payload with a bare json.dumps; use _script_json so `</script` is neutralized"
    )


# ── the page carries a few structures, so it is inspectable off the cluster ────

def test_a_structure_stream_embeds_the_head_plus_a_sample():
    """A page copied from the cluster to a laptop has no access to the compute node's filesystem, so a structure is only inspectable if its contents are inline -- which caps how many can go in. The head shows what the run produced; the sample is what shows whether quality holds across the run.

    Below the cap every structure is embedded, so a small run loses nothing.
    """
    module = _load_renderer("renderers/structures.py")
    assert module.MAX_EMBEDDED == module.SAMPLE_HEAD + module.SAMPLE_RANDOM

    for total in range(module.MAX_EMBEDDED + 1):
        assert module.sample_positions(total, seed="s") == list(range(total)), total

    picked = module.sample_positions(48, seed="structures:48")
    assert len(picked) == module.MAX_EMBEDDED
    assert picked[:module.SAMPLE_HEAD] == list(range(module.SAMPLE_HEAD))
    assert picked == sorted(picked), "positions stay in stream order so the viewer reads left to right"
    assert all(p >= module.SAMPLE_HEAD for p in picked[module.SAMPLE_HEAD:]), (
        "the sampled tail must come from beyond the head, or it is not a sample of the rest"
    )


def test_the_sample_is_stable_for_a_stream_and_differs_between_streams():
    """The page is regenerated after every completed step. A fresh sample each time would swap the structures under the reader and make two renders of an unchanged run differ, so the choice is seeded rather than random."""
    module = _load_renderer("renderers/structures.py")
    assert module.sample_positions(48, seed="structures:48") == module.sample_positions(48, seed="structures:48")
    assert module.sample_positions(48, seed="structures:48") != module.sample_positions(48, seed="designs:48"), (
        "two streams in one page should not sample identical positions"
    )


def test_a_twelve_structure_stream_embeds_five_and_says_so(
    local_config, isolated_cwd, new_pipeline, monkeypatch, record_case,
):
    """End to end: one stream of twelve structures, one viewer, five files inline, and the label saying which five.

    This is what makes the page worth copying off the cluster -- the reader can turn the structures themselves, not follow a link to a filesystem they cannot reach.
    """
    from biopipelines.config_manager import ConfigManager
    from biopipelines.mock import Mock
    from biopipelines.pipeline import regenerate_pipeline_page

    monkeypatch.setattr(ConfigManager, "get_renderers_config", lambda self: _REAL_RENDERERS)

    ids = [f"d{i:02d}" for i in range(12)]
    pipeline = new_pipeline("page_embed_sample")
    with pipeline:
        out = Mock(ids=ids, streams=_mock_stream(), map_table_strategy="config")
        pipeline.save()
    for n, path in enumerate(out.streams.structures.files_expanded):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            f.write(f"REMARK   1 STRUCTURE {n}\nATOM      1  CA  ALA A   1       0.000   0.000   0.000\nEND\n")

    page = open(regenerate_pipeline_page(pipeline.folders["runtime"]), encoding="utf-8").read()

    embedded = sorted(set(re.findall(r'"(d\d\d)"', page)))
    remarks = sorted(set(re.findall(r"STRUCTURE (\d+)", page)))
    record_case(input="12 structures on disk, one pdb stream",
                expected=("5 embedded", "label names head+sample"),
                actual=(f"{len(embedded)} embedded", embedded))
    assert len(remarks) == 5, f"expected 5 files inline, got {len(remarks)}: {remarks}"
    assert remarks[:3] == ["0", "1", "2"], f"the head must be the first three, got {remarks}"
    assert re.search(r"of 12 total — first 3 plus 2 sampled", page), "the label must say what was left out"
    assert len(re.findall(r'id="bp3d_[0-9a-f]+_viewer"', page)) == 1


def test_regenerating_a_page_with_viewers_is_byte_identical(
    local_config, isolated_cwd, new_pipeline, monkeypatch,
):
    """The cluster refreshes the page after every completed step, so two renders of an unchanged run must agree. `viewer_id` was `random.randint`, which made them differ on every element id -- invisible to the existing regeneration test because its fixture config declares no renderers and so emits no viewer at all.

    The page carries one deliberately varying field, the `Generated <timestamp>` footer (`pipeline_report.py:853`), so that one is normalized out rather than asserted on: comparing raw bytes made this test straddle a second boundary and fail roughly one run in ten.
    """
    from biopipelines.pipeline import regenerate_pipeline_page

    pipeline = _pipeline_with_structures_on_disk(
        new_pipeline, monkeypatch, "page_regen_viewers", ids=("x", "y", "z"))
    first = open(regenerate_pipeline_page(pipeline.folders["runtime"]), encoding="utf-8").read()
    second = open(regenerate_pipeline_page(pipeline.folders["runtime"]), encoding="utf-8").read()

    assert re.findall(r'id="bp3d_[0-9a-f]+_viewer"', first), "fixture emitted no viewer, so this proves nothing"
    stamp = re.compile(r"Generated \d{4}-\d\d-\d\dT\d\d:\d\d:\d\d")
    assert stamp.search(first), "the timestamp moved; this test would now be normalizing nothing"
    assert stamp.sub("Generated <stamp>", first) == stamp.sub("Generated <stamp>", second)
