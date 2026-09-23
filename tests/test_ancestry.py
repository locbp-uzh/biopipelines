"""Per-ID ancestry: which design came from which input.

`lineage` counts IDs; this connects them. The rules under test are the ones in
`docs/design/id-ancestry.md`, and each of the four table shapes a real campaign produces has a
case here, because three of them cannot be answered by joining on a provenance column:

  * a concrete `<axis>.id` cell naming an upstream id          -> tier `recorded`
  * a combinatorial pattern (`3kzy_A_<1..100>_<1..2>`)         -> tier `matched`
  * no provenance column at all, pruned as recoverable         -> tier `matched`
  * neither resolves                                           -> reported as a root

The shapes and values are taken from `BinderDesign/tagged_binder_001`, an eight-step campaign
that ran on S3IT.
"""

import json

import pytest

from biopipelines import ancestry


class FakeFS:
    """The two reads ancestry performs, with nothing else implemented."""

    def __init__(self, tables, graph=None):
        self.tables = tables
        self.graph = graph
        self.id_column_calls = 0

    def join(self, *parts):
        return "/".join(str(p).rstrip("/") for p in parts if str(p) != "")

    def id_columns(self, root, names):
        self.id_column_calls += 1
        return self.tables

    def read_text(self, path):
        if self.graph is None:
            raise FileNotFoundError(path)
        return json.dumps(self.graph)


def graph_for(edges, steps):
    """A `pipeline_graph.json` in the shape `save()` writes."""
    return {
        "schema": 1,
        "steps": [{"execution_order": order, "tool": tool,
                   "output_folder": f"/out/job_001/{folder}"}
                  for order, tool, folder in steps],
        "edges": [{"from_step": a, "to_step": b, "stream": "structures"} for a, b in edges],
    }


@pytest.fixture
def campaign():
    """The real campaign's shapes: PDB -> RFdiffusion -> MPNN -> Boltz2 -> Panda, plus a tag."""
    tables = {
        "001_PDB/structures/structures_map.csv": [{"id": "3kzy_A"}],
        "002_Sequence/sequences/sequences.csv": [{"id": "tag"}],
        # concrete parent
        "003_RFdiffusion/structures/structures_map.csv": [
            {"id": "3kzy_A_50", "structures.id": "3kzy_A"},
            {"id": "3kzy_A_97", "structures.id": "3kzy_A"}],
        # concrete parent
        "004_ProteinMPNN/sequences/sequences.csv": [
            {"id": "3kzy_A_50_2", "structures.id": "3kzy_A_50"},
            {"id": "3kzy_A_97_1", "structures.id": "3kzy_A_97"}],
        # a PATTERN, not a parent -- the case a join cannot answer
        "005_Boltz2/structures/structures_map.csv": [
            {"id": "3kzy_A_50_2+tag", "proteins.id": "3kzy_A_<1..100>_<1..2>"},
            {"id": "3kzy_A_97_1+tag", "proteins.id": "3kzy_A_<1..100>_<1..2>"}],
        # concrete parent again
        "009_Panda/structures/structures_map.csv": [
            {"id": "9_Panda_1", "structures.id": "3kzy_A_50_2+tag"},
            {"id": "9_Panda_2", "structures.id": "3kzy_A_97_1+tag"}],
    }
    wiring = graph_for(
        edges=[(1, 3), (3, 4), (2, 5), (4, 5), (5, 9)],
        steps=[(1, "PDB", "001_PDB"), (2, "Sequence", "002_Sequence"),
               (3, "RFdiffusion", "003_RFdiffusion"), (4, "ProteinMPNN", "004_ProteinMPNN"),
               (5, "Boltz2", "005_Boltz2"), (9, "Panda", "009_Panda")])
    return FakeFS(tables, wiring)


def parents(graph, child):
    return sorted((e["parent_step"], e["parent"]) for e in graph["edges"] if e["child"] == child)


class TestTheTiers:

    def test_a_concrete_provenance_cell_is_taken_as_recorded(self, campaign):
        graph = ancestry.collect("/out/job_001", fs=campaign)
        edge = next(e for e in graph["edges"] if e["child"] == "3kzy_A_50")
        assert (edge["parent"], edge["tier"]) == ("3kzy_A", "recorded")

    def test_a_pattern_cell_never_becomes_a_parent(self, campaign):
        """`3kzy_A_<1..100>_<1..2>` is a combinatorial declaration, not an identity."""
        graph = ancestry.collect("/out/job_001", fs=campaign)
        assert not [e for e in graph["edges"] if "<" in e["parent"]]

    def test_the_pattern_case_still_finds_both_real_parents(self, campaign):
        """The case the whole module exists for: a join answers nothing here."""
        graph = ancestry.collect("/out/job_001", fs=campaign)
        assert parents(graph, "3kzy_A_50_2+tag") == [
            ("002_Sequence", "tag"), ("004_ProteinMPNN", "3kzy_A_50_2")]

    def test_a_matched_edge_carries_the_tier_that_answered(self, campaign):
        graph = ancestry.collect("/out/job_001", fs=campaign)
        edge = next(e for e in graph["edges"]
                    if e["child"] == "3kzy_A_50_2+tag" and e["parent"] == "tag")
        assert edge["tier"].startswith("matched:"), (
            "a reader who cannot see how a link was established cannot audit it")


class TestTheWiringScopesTheSearch:

    def test_a_source_step_is_given_no_parents(self, campaign):
        """The bug this pins: a step with no declared upstream once fell back to scanning
        every earlier step, and `002_Sequence` acquired `001_PDB` as a parent it never had."""
        graph = ancestry.collect("/out/job_001", fs=campaign)
        assert parents(graph, "tag") == []

    def test_a_grandparent_is_not_reported_as_a_parent(self, campaign):
        """`3kzy_A_50_2+tag` matches the RFdiffusion backbone as readily as the MPNN sequence."""
        assert ("003_RFdiffusion", "3kzy_A_50") not in parents(
            ancestry.collect("/out/job_001", fs=campaign), "3kzy_A_50_2+tag")

    def test_without_a_recorded_graph_the_grandparent_is_still_dropped(self, campaign):
        """An older run has no `pipeline_graph.json`; the transitive drop has to do the work."""
        blind = FakeFS(campaign.tables, graph=None)
        found = parents(ancestry.collect("/out/job_001", fs=blind), "3kzy_A_50_2+tag")
        assert ("003_RFdiffusion", "3kzy_A_50") not in found
        assert ("004_ProteinMPNN", "3kzy_A_50_2") in found


class TestTheWalkBack:

    def test_a_design_traces_to_the_structure_it_came_from(self, campaign):
        graph = ancestry.collect("/out/job_001", fs=campaign)
        found = ancestry.trace(graph, "9_Panda_1")
        seen = []

        def walk(node):
            seen.append(node["id"])
            for parent in node["parents"]:
                walk(parent)

        walk(found)
        assert seen[0] == "9_Panda_1"
        assert "3kzy_A" in seen, "the PDB entry is the question's whole point"
        assert "tag" in seen, "dropping the tag axis hides half a binder campaign"

    def test_an_unknown_id_is_reported_as_absent_not_as_rootless(self, campaign):
        assert ancestry.trace(ancestry.collect("/out/job_001", fs=campaign), "nope") is None

    def test_the_rendered_tree_names_the_step_for_every_generation(self, campaign):
        graph = ancestry.collect("/out/job_001", fs=campaign)
        text = "\n".join(ancestry.render_trace(ancestry.trace(graph, "9_Panda_1")))
        for step in ("005_Boltz2", "004_ProteinMPNN", "003_RFdiffusion", "001_PDB"):
            assert step in text


class TestWhatItRefusesToInvent:

    def test_an_id_with_no_resolvable_parent_is_reported_not_hidden(self):
        tables = {"001_A/structures/structures_map.csv": [{"id": "alpha"}],
                  "002_B/structures/structures_map.csv": [{"id": "unrelated_thing"}]}
        wiring = graph_for(edges=[(1, 2)],
                           steps=[(1, "A", "001_A"), (2, "B", "002_B")])
        graph = ancestry.collect("/out/job_001", fs=FakeFS(tables, wiring))
        assert graph["edges"] == []
        assert [u["id"] for u in graph["unresolved"]] == ["unrelated_thing"]

    def test_a_run_with_no_map_tables_says_so_rather_than_returning_empty(self):
        graph = ancestry.collect("/out/job_001", fs=FakeFS({}, None))
        assert "no ancestry" in ancestry.summarize(graph, job="job_001").lower()


class TestCost:

    def test_the_whole_run_costs_one_projection_call(self, campaign):
        """One ssh round trip, not one per step: a large campaign has hundreds."""
        ancestry.collect("/out/job_001", fs=campaign)
        assert campaign.id_column_calls == 1

    def test_rows_are_one_per_link(self, campaign):
        graph = ancestry.collect("/out/job_001", fs=campaign)
        rows = ancestry.to_rows(graph)
        assert len(rows) == len(graph["edges"])
        assert set(rows[0]) == {"child_step", "child", "parent_step", "parent", "tier"}


class TestFindingsFromThe151Review:

    def test_a_second_source_step_is_a_root_not_an_unresolved_id(self, campaign):
        graph = ancestry.collect("/out/job_001", fs=campaign)
        assert not any(u["step"] == "002_Sequence" for u in graph["unresolved"])
        assert "no traceable parent" not in ancestry.summarize(graph)

    def test_an_id_that_passes_through_unchanged_keeps_its_chain(self):
        """MPNN x_1_1 -> Boltz2 x_1_1 -> Filter x_1_1: keyed by id alone, the walk flattened
        the chain and listed the steps as siblings."""
        tables = {
            "002_RFdiffusion/structures/structures_map.csv": [{"id": "x_1"}],
            "003_ProteinMPNN/sequences/sequences.csv": [{"id": "x_1_1", "structures.id": "x_1"}],
            "004_Boltz2/structures/structures_map.csv": [{"id": "x_1_1", "proteins.id": "x_1_1"}],
            "005_Panda/structures/structures_map.csv": [{"id": "x_1_1", "structures.id": "x_1_1"}],
        }
        wiring = graph_for(edges=[(2, 3), (3, 4), (4, 5)],
                           steps=[(2, "RFdiffusion", "002_RFdiffusion"), (3, "ProteinMPNN", "003_ProteinMPNN"),
                                  (4, "Boltz2", "004_Boltz2"), (5, "Panda", "005_Panda")])
        found = ancestry.trace(ancestry.collect("/out/job_001", fs=FakeFS(tables, wiring)), "x_1_1")
        chain, node = [], found
        while node:
            chain.append(node["step"])
            node = node["parents"][0] if node["parents"] else None
        assert chain == ["005_Panda", "004_Boltz2", "003_ProteinMPNN", "002_RFdiffusion"]

    def test_an_echoed_stream_does_not_inflate_the_step_count(self):
        tables = {
            "001_PDB/structures/structures_map.csv": [{"id": "a"}],
            "002_Boltz2/structures/structures_map.csv": [{"id": "a+lig", "proteins.id": "a"}],
            "002_Boltz2/msas/msas_map.csv": [{"id": "a"}],
            "002_Boltz2/compounds/compounds_map.csv": [{"id": "lig"}],
        }
        text = ancestry.summarize(ancestry.collect("/out/job_001", fs=FakeFS(tables)))
        assert "002_Boltz2" in text and "    1 ID(s)" in text.split("002_Boltz2")[1].splitlines()[0]
