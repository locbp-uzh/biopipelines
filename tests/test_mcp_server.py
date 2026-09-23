"""`bp_tools` is the contract three files now promise an agent: SKILL.md, llm/pipelines.md and
the index header all say "call bp_tools with a tool name and you get that tool's section".

These tests hold the server to that promise without a transport: they build the server object
and call the registered tool directly. The MCP SDK is an optional extra, so everything here
skips cleanly on an install that does not have it.
"""

import json

import pytest

from biopipelines import tool_docs

pytest.importorskip("mcp", reason="bp-mcp needs the MCP SDK: pip install -e '.[mcp]'")

from biopipelines import mcp_server  # noqa: E402


@pytest.fixture(scope="module")
def server():
    return mcp_server.build_server()


@pytest.fixture(autouse=True)
def no_saved_cluster(tmp_path_factory, monkeypatch):
    """Without this, a real saved connection sends every test's local tmp path to the cluster."""
    from biopipelines import remote
    monkeypatch.setattr(remote, "SETTINGS_FILE",
                        tmp_path_factory.mktemp("settings") / ".bp-mcp.json")


async def _call(server, **arguments):
    """The text an agent actually receives, out of the CallToolResult's content blocks."""
    result = await server.call_tool("bp_tools", arguments)
    return "\n".join(b.text for b in result.content if getattr(b, "type", None) == "text")


@pytest.mark.anyio
async def test_bp_tools_without_a_name_returns_the_catalog(server):
    text = await _call(server)
    assert "# BioPipelines tool index" in text
    for name in ("ProteinMPNN", "Boltz2", "RFdiffusion3", "PoseBusters"):
        assert f"**{name}**" in text, f"{name} missing from the catalog"


@pytest.mark.anyio
async def test_bp_tools_with_a_name_returns_that_tools_section(server):
    text = await _call(server, name="Boltz2")
    assert text.lstrip().startswith("### Boltz2")
    assert "ProteinMPNN" not in text.split("\n### ")[0][:200], "bled into a neighbouring section"


@pytest.mark.anyio
async def test_an_unknown_name_suggests_instead_of_inventing(server):
    text = await _call(server, name="ProteinMPN")
    assert "No BioPipelines tool named" in text
    assert "ProteinMPNN" in text, "a near-miss should be suggested, not left to the model"


@pytest.mark.anyio
async def test_tags_are_or_within_a_facet_and_and_across_them(server):
    """The query that motivated exclusion: protein binding tools that are not ligand tools."""
    text = await _call(server, tags=["protein", "binding"], exclude=["small-molecule"])
    assert "**Prodigy**" in text
    assert "**Boltz2**" not in text, "Boltz2 carries small-molecule and must be excluded"


@pytest.mark.anyio
async def test_a_retired_term_is_translated_rather_than_rejected(server):
    """`affinity` was the tag's name during design and is what an agent will type."""
    text = await _call(server, tags=["affinity"])
    assert "read 'affinity' as 'binding'" in text
    assert "**GEMS**" in text


@pytest.mark.anyio
async def test_an_unknown_tag_prints_the_whole_vocabulary(server):
    text = await _call(server, tags=["zzz"])
    assert "Not tags: zzz" in text
    for facet in ("action", "subject", "readout", "capability"):
        assert facet in text
    assert "covalent" in text


@pytest.mark.anyio
async def test_an_empty_result_explains_the_matching_rule(server):
    text = await _call(server, tags=["dock", "covalent"])
    assert "No tool matches" in text
    assert "OR-ed" in text and "AND-ed" in text


@pytest.mark.anyio
async def test_the_tool_is_advertised_with_a_usable_description(server):
    tools = await server.list_tools()
    entry = next(t for t in tools if t.name == "bp_tools")
    assert entry.description and "silently ignored" in entry.description, (
        "the description must warn that unknown kwargs are swallowed, not rejected")
    schema = entry.inputSchema if hasattr(entry, "inputSchema") else entry.input_schema
    assert "name" in json.dumps(schema)


def test_every_catalog_entry_can_be_fetched_by_name():
    """A pointer in the index that `section()` cannot resolve would strand the agent."""
    missing = [e["name"] for e in tool_docs.collect() if tool_docs.section(e["name"]) is None]
    assert not missing, f"listed in the index but not fetchable: {missing}"


def _build(job, steps):
    import io as _io
    (job / "Logs").mkdir(parents=True, exist_ok=True)
    for name, marker in steps:
        (job / name).mkdir(exist_ok=True)
        _io.open(job / "Logs" / f"{name}.log", "w", encoding="utf-8").write(f"line of {name}\n")
        if marker:
            _io.open(job / f"{name}_{marker}", "w", encoding="utf-8").write("")
    return job


@pytest.mark.anyio
async def test_bp_status_names_the_first_failure(server, tmp_path):
    job = _build(tmp_path / "run_001", [("001_A", "COMPLETED"), ("002_B", "FAILED")])
    result = await server.call_tool("bp_status", {"job": str(job)})
    text = "\n".join(b.text for b in result.content if getattr(b, "type", None) == "text")
    assert "1 completed, 1 failed" in text and "First failure: 002_B" in text


@pytest.mark.anyio
async def test_bp_logs_lists_what_is_available_when_the_step_is_wrong(server, tmp_path):
    job = _build(tmp_path / "run_001", [("001_A", "COMPLETED")])
    result = await server.call_tool("bp_logs", {"job": str(job), "step": "999_Nope"})
    text = "\n".join(b.text for b in result.content if getattr(b, "type", None) == "text")
    assert "No log for" in text and "001_A" in text, (
        "a wrong step name should point at the real ones, not just fail")


@pytest.mark.anyio
async def test_bp_project_scaffolds_then_surveys(server, tmp_path):
    project = tmp_path / "SNAP33"
    (project / "Binder_006").mkdir(parents=True)

    async def call(**args):
        result = await server.call_tool("bp_project", {"project_dir": str(project), **args})
        return "\n".join(b.text for b in result.content if getattr(b, "type", None) == "text")

    assert "project: created" in await call(action="scaffold")
    assert "Appended to HISTORY.md" in await call(action="history", title="Run 006")
    survey = await call(action="survey")
    assert "project=yes" in survey and "Binder_006" in survey and "Run 006" in survey


@pytest.mark.anyio
async def test_bp_project_history_without_a_title_says_so(server, tmp_path):
    result = await server.call_tool("bp_project",
                                    {"project_dir": str(tmp_path), "action": "history"})
    text = "\n".join(b.text for b in result.content if getattr(b, "type", None) == "text")
    assert "needs a `title`" in text


@pytest.mark.anyio
async def test_bp_project_rejects_an_unknown_action(server, tmp_path):
    result = await server.call_tool("bp_project",
                                    {"project_dir": str(tmp_path), "action": "rewrite"})
    text = "\n".join(b.text for b in result.content if getattr(b, "type", None) == "text")
    assert "Unknown action" in text, "an agent must not be able to invent a destructive action"


@pytest.mark.anyio
async def test_the_server_advertises_every_tool(server):
    names = {t.name for t in await server.list_tools()}
    assert names == {"bp_tools", "bp_setup", "bp_runs", "bp_status", "bp_logs",
                     "bp_lineage", "bp_table", "bp_submit", "bp_resubmit", "bp_cancel",
                     "bp_visualize", "bp_fetch", "bp_project", "bp_provenance",
                     "bp_reproduce", "bp_ancestry"}


@pytest.mark.anyio
async def test_run_reading_tools_take_a_host(server):
    """The cluster is the main platform: reading a run must not assume a local path.

    `bp_project` is in this list because it shipped without `host` and answered about a local
    path that did not exist — reporting a project full of runs as having no documents and no
    jobs, which is indistinguishable from an empty project.
    """
    for name in ("bp_runs", "bp_status", "bp_logs", "bp_setup", "bp_table", "bp_lineage",
                 "bp_project", "bp_provenance", "bp_reproduce", "bp_ancestry"):
        entry = next(t for t in await server.list_tools() if t.name == name)
        schema = entry.inputSchema if hasattr(entry, "inputSchema") else entry.input_schema
        assert "host" in json.dumps(schema), f"{name} cannot reach a cluster"


@pytest.mark.anyio
async def test_an_unreachable_host_is_reported_not_raised(server):
    result = await server.call_tool("bp_status",
                                    {"job": "run_001", "host": "no-such-host.invalid"})
    text = "\n".join(b.text for b in result.content if getattr(b, "type", None) == "text")
    assert "Could not read" in text and "no-such-host.invalid" in text


@pytest.mark.anyio
async def test_bp_setup_explains_what_to_do_when_nothing_is_configured(server):
    result = await server.call_tool("bp_setup", {})
    text = "\n".join(b.text for b in result.content if getattr(b, "type", None) == "text")
    assert "No cluster is configured yet" in text
    assert "~/.ssh/config" in text, "the alias lives in a file these tools cannot reach"


@pytest.mark.anyio
async def test_an_explicit_empty_host_forces_local_even_with_a_saved_cluster(server, tmp_path):
    """The escape hatch: `host=""` reads this machine, whatever is saved."""
    from biopipelines import remote
    remote.save_settings(host="cluster", repo="~/biopipelines-locbp")
    job = _build(tmp_path / "run_001", [("001_A", "COMPLETED")])
    result = await server.call_tool("bp_status", {"job": str(job), "host": ""})
    text = "\n".join(b.text for b in result.content if getattr(b, "type", None) == "text")
    assert "1 completed" in text


@pytest.mark.anyio
async def test_bp_submit_refuses_without_a_cluster(server):
    """Submitting is the one tool that spends money; it must not quietly fall back to local."""
    result = await server.call_tool("bp_submit", {"script": "my_pipelines/foo.py"})
    text = "\n".join(b.text for b in result.content if getattr(b, "type", None) == "text")
    assert "needs a cluster" in text and "bp_setup" in text


@pytest.mark.anyio
async def test_bp_submit_is_advertised_as_spending_compute(server):
    entry = next(t for t in await server.list_tools() if t.name == "bp_submit")
    assert "SPENDS COMPUTE" in entry.description
    assert "poll bp_status" in entry.description, (
        "an accepted job is not a finished job; the description has to say so")


@pytest.mark.anyio
async def test_bp_fetch_refuses_without_a_cluster(server):
    result = await server.call_tool("bp_fetch", {"remote_path": "/x/y.pdb", "local": "y.pdb"})
    text = "\n".join(b.text for b in result.content if getattr(b, "type", None) == "text")
    assert "needs a cluster" in text


@pytest.mark.anyio
async def test_bp_lineage_reports_attrition(server, tmp_path):
    import io as _io
    job = tmp_path / "run_001"
    (job / "004_Filter").mkdir(parents=True)
    _io.open(job / "004_Filter" / "missing.csv", "w", encoding="utf-8").write(
        'id,structure\nd1,a.pdb\nd2,b.pdb\n')
    result = await server.call_tool("bp_lineage", {"job": str(job)})
    text = "".join(b.text for b in result.content if getattr(b, "type", None) == "text")
    assert "2 ID(s) dropped" in text and "004_Filter" in text


@pytest.mark.anyio
async def test_bp_table_lists_then_reads(server, tmp_path):
    import io as _io
    job = tmp_path / "run_001"
    # The real layout: standalone tables under `tables/`, never beside the step. This fixture
    # used to put the CSV at the step root, which is why bp_table shipped unable to see one.
    (job / "001_A" / "tables").mkdir(parents=True)
    _io.open(job / "001_A" / "tables" / "analysis.csv", "w", encoding="utf-8").write(
        'id,plddt\nd1,0.9\nd2,0.8\n')

    async def call(**args):
        r = await server.call_tool("bp_table", {"job": str(job), "step": "001_A", **args})
        return "".join(b.text for b in r.content if getattr(b, "type", None) == "text")

    assert "analysis.csv" in await call()
    body = await call(table="analysis", limit=1)
    assert "1 of 2 rows" in body and "id,plddt" in body
    assert "Tables in that step" in await call(table="nope.csv")


@pytest.mark.anyio
async def test_the_tools_read_a_run_the_framework_actually_produced(server, produced_run):
    """End to end on a real tree: the case that would have caught the bp_table defect.

    Every other fixture here hand-builds the tree it then reads, so it tests the reader against
    its author's assumption. `produced_run` builds it with the framework, which is why this one
    could not have passed while `bp_table` listed `<step>/*.csv`.
    """
    job, steps = produced_run([
        {"streams": {"structures": 2}, "tables": {"metrics": [{"id": "d0", "n_res": 164}]}},
    ])
    step = next(iter(steps))

    async def call(tool, **arguments):
        result = await server.call_tool(tool, {"job": str(job), **arguments})
        return "\n".join(b.text for b in result.content if getattr(b, "type", None) == "text")

    listed = await call("bp_table", step=step)
    assert "tables/metrics.csv" in listed and "structures/structures_map.csv" in listed

    body = await call("bp_table", step=step, table="metrics")
    assert "id,n_res" in body and "164" in body

    assert "completed" in (await call("bp_status")).lower()
    assert "structures 2" in await call("bp_lineage")


@pytest.fixture
def settings(tmp_path, monkeypatch):
    """A settings file of this test's own, so a real saved cluster cannot leak in."""
    from biopipelines import remote
    monkeypatch.setattr(remote, "SETTINGS_FILE", tmp_path / ".bp-mcp.json")
    return remote


class TestTwoClusters:
    """A lab with access to both S3IT and Daint is the ordinary case.

    The flat settings file held one repo, one interpreter and no variant, and `defaults()`
    returned that repo for whichever host was named — so a call against Daint was sent S3IT's
    checkout path. That is the failure that answers confidently about the wrong machine.
    """

    async def call(self, server, **args):
        result = await server.call_tool("bp_setup", args)
        return "\n".join(b.text for b in result.content if getattr(b, "type", None) == "text")

    @pytest.mark.anyio
    async def test_adding_a_second_cluster_keeps_the_first(self, server, settings):
        settings.save_settings(host="cluster", repo="~/biopipelines-locbp")
        settings.save_settings(host="daint", repo="/scratch/me/bp", variant="daint")
        listing = await self.call(server)
        assert "cluster" in listing and "daint" in listing
        assert "repo=~/biopipelines-locbp" in listing and "repo=/scratch/me/bp" in listing
        assert "variant=daint" in listing

    @pytest.mark.anyio
    async def test_the_listing_marks_which_one_unnamed_calls_use(self, server, settings):
        settings.save_settings(host="cluster", repo="~/a")
        settings.save_settings(host="daint", repo="~/b", make_default=False)
        listing = await self.call(server)
        assert "cluster (default)" in listing and "daint (default)" not in listing

    @pytest.mark.anyio
    async def test_with_nothing_saved_it_says_how_to_start(self, server, settings):
        assert "No cluster is configured yet" in await self.call(server)

    def test_no_setting_is_inherited_across_hosts(self, settings):
        settings.save_settings(host="cluster", repo="~/biopipelines-locbp",
                               python="/opt/conda/bin/python")
        host, repo, python, variant, _prelude = settings.connection("daint")
        assert repo == settings.DEFAULT_REPO, "S3IT's checkout must not be sent to Daint"
        assert python == "python" and variant is None

    def test_a_pinned_variant_reaches_the_remote_command(self, settings):
        """Daint generates against `daint`; the framework would otherwise auto-detect."""
        ssh = settings.Ssh("daint", repo="/scratch/me/bp", variant="daint")
        assert ssh.env_prefix == "BIOPIPELINES_CONFIG_VARIANT=daint "
        assert settings.Ssh("cluster", repo="~/a").env_prefix == ""
