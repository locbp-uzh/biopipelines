"""Finding a tool by what it produces or consumes.

"Which tool gives me an RMSD" is the question a tag vocabulary answers worst: it would need one
tag per metric, hand-applied and drifting. The answer is already in the code — every tool
declares its streams, tables and columns in `get_output_files()` — so this search is derived
from verified source rather than from new metadata.

The matching rule is the load-bearing part: case-insensitive substring, so `rmsd` must find
`ligand_rmsd`, `RMSD_before` and `rmsd_to_ref`. A user who has to know the exact column name
already knows which tool they want.
"""

import pytest

from biopipelines import tool_schemas


class TestOutputsIndex:
    def test_the_schema_source_is_reachable(self):
        """If the AST extractor moves or breaks, the search silently returns nothing."""
        out = tool_schemas.outputs()
        assert len(out) > 50, f"only {len(out)} tools have schemas; the extractor is broken"

    def test_columns_streams_and_tables_are_all_indexed(self):
        out = tool_schemas.outputs()
        assert out["SASA"]["streams"], "SASA declares an accessibility stream"
        assert "sasa" in out["SASA"]["tables"]
        assert any("delta_sasa" == c for c in out["SASA"]["columns"])


class TestMatching:
    def test_a_fragment_finds_the_longer_column(self):
        """The whole point: `rmsd` must reach `ligand_rmsd` and `rmsd_to_ref`."""
        names = [n for n, _ in tool_schemas.search("rmsd")]
        assert {"PoseChange", "EnsembleAnalysis", "ConformationalChange"} <= set(names)

    def test_matching_ignores_case_in_both_directions(self):
        """ThermoMPNN's column is `ddG_pred`; a lowercase query must still find it."""
        assert "ThermoMPNN" in [n for n, _ in tool_schemas.search("ddg")]
        assert "ConformationalChange" in [n for n, _ in tool_schemas.search("RMSD")]

    def test_the_match_reports_which_name_hit(self):
        matched = dict(tool_schemas.search("rmsd"))["PoseChange"]
        assert any("ligand_rmsd" in v for v in matched.get("columns", []))

    def test_an_unmatched_query_returns_nothing_rather_than_everything(self):
        assert tool_schemas.search("zzznotacolumn") == []

    def test_an_empty_query_is_not_a_wildcard(self):
        assert tool_schemas.search("") == []
        assert tool_schemas.search("   ") == []


class TestInputs:
    def test_input_streams_come_from_the_parameter_types(self):
        """A parameter typed DataStream is a stream the tool consumes."""
        names = [n for n, _ in tool_schemas.search("structures", where="inputs")]
        assert len(names) > 20, "most tools take a structures stream"
        assert "AF2BIND" in names and "APBS" in names

    def test_a_renamed_input_is_still_found_by_fragment(self):
        matched = dict(tool_schemas.search("structures", where="inputs"))["ConformationalChange"]
        assert any("reference_structures" in v for v in matched["inputs"])

    def test_outputs_and_inputs_are_different_questions(self):
        """A tool that consumes compounds need not produce them."""
        produces = {n for n, _ in tool_schemas.search("compounds", where="outputs")}
        consumes = {n for n, _ in tool_schemas.search("compounds", where="inputs")}
        assert produces != consumes and (consumes - produces)


class TestSummary:
    def test_a_miss_explains_the_matching_rule(self):
        text = tool_schemas.summarize("zzznotacolumn")
        assert "case-insensitive substring" in text

    def test_a_hit_names_the_column_not_just_the_tool(self):
        text = tool_schemas.summarize("plddt")
        assert "AlphaFold" in text and "plddt" in text

    def test_the_summary_says_where_the_data_came_from(self):
        """An agent should know these are declared in source, not hand-maintained tags."""
        assert "not hand-maintained tags" in tool_schemas.summarize("rmsd")


class TestUnits:
    """A column name without its scale invites exactly the comparison that is wrong.

    Boltz2's `affinity_pred_value` is log10(IC50) where LOWER is stronger; GEMS's `pkd_pred` is a
    pKd where HIGHER is stronger. An agent handed both names and no units will rank one of them
    backwards, and the number it produces will look perfectly reasonable.
    """

    def test_the_documented_scale_is_returned(self):
        note = tool_schemas.units("Boltz2")
        assert "log10(IC50)" in note and "lower is stronger" in note.lower()

    def test_opposite_directions_are_both_captured(self):
        assert "higher is stronger" in tool_schemas.units("GEMS").lower()
        assert "lower is stronger" in tool_schemas.units("Boltz2").lower()

    def test_a_note_stops_at_the_next_block(self):
        """Without a trailing blank line the note used to swallow the following table row."""
        assert "removed_by" not in tool_schemas.units("GEMS")

    def test_a_tool_with_no_note_returns_empty_not_garbage(self):
        assert tool_schemas.units("Panda") == ""

    def test_units_reach_the_search_output(self):
        text = tool_schemas.summarize("affinity")
        assert "units:" in text
        assert "log10(IC50)" in text, "the scale must be visible without a second call"

    def test_undocumented_scales_say_so_rather_than_guessing(self):
        """DynamicBind and RTMScore are genuinely undocumented upstream; silence would be worse."""
        for tool in ("DynamicBind", "RTMScore"):
            note = tool_schemas.units(tool).lower()
            assert "does not" in note or "not state" in note or "not calibrated" in note, tool
