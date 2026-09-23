"""`head(n)` declares n outputs; a filter decides how many rows there are to fill them.

The two are set at different times. `Panda(rename="best", operations=[filter(...), head(10)])`
declares `best_1..best_10` at configuration time, because that is when the ids have to be known
for anything downstream to reference them. How many rows survive the filter is a runtime fact.

When fewer survive, the surplus slots are declared ids that nothing produced. They name no
upstream id, so no id-matching tier can excuse them from the missing manifest the dropped inputs
wrote — and the completion check then demands files that correctly have nothing to contain. This
is the same shape as the bug where a filter keeping nothing marked a step FAILED: the framework
reporting its own correct behaviour as a fault.

These tests pin the accounting: the unfilled slots are recorded under the ids the tool declared
them by, a full run records nothing, and the cause names the operation that set the count.
"""

import pytest

from pipe_scripts.pipe_panda import unfilled_slot_rows

DECLARED = ["best_1", "best_2", "best_3", "best_4", "best_5"]
HEAD_5 = [{"type": "filter", "params": {"expr": "score > 0"}},
          {"type": "head", "params": {"n": 5}}]


class TestWhatIsRecorded:

    def test_the_slots_no_row_reached_are_named(self):
        rows = unfilled_slot_rows(DECLARED, produced_count=2, step_tool_name="007_Panda",
                                  operations=HEAD_5)
        assert [r["id"] for r in rows] == ["best_3", "best_4", "best_5"], (
            "the unfilled slots must be recorded under the ids the tool declared, since that is "
            "the only spelling anything downstream knows them by")
        assert all(r["removed_by"] == "007_Panda" for r in rows)
        assert all(r["kind"] == "filter" for r in rows)

    def test_the_cause_names_the_operation_that_set_the_count(self):
        rows = unfilled_slot_rows(DECLARED, produced_count=2, step_tool_name="007_Panda",
                                  operations=HEAD_5)
        cause = rows[0]["cause"]
        assert "head(5)" in cause, f"the cause does not say what fixed the count: {cause!r}"
        assert "2 rows survived" in cause

    def test_the_tightest_limit_wins(self):
        """Two limiting operations in one chain: the smaller one decides the declared count."""
        operations = [{"type": "head", "params": {"n": 20}},
                      {"type": "sample", "params": {"n": 5}}]
        rows = unfilled_slot_rows(DECLARED, produced_count=1, step_tool_name="s",
                                  operations=operations)
        assert "sample(5)" in rows[0]["cause"]


class TestWhatIsNotRecorded:
    """Over-reporting here is the failure mode: a slot wrongly called missing is attrition
    invented by the framework, which is what a lineage layer must never do."""

    def test_a_full_run_records_nothing(self):
        assert unfilled_slot_rows(DECLARED, produced_count=5, step_tool_name="s",
                                  operations=HEAD_5) == []

    def test_more_rows_than_slots_records_nothing(self):
        assert unfilled_slot_rows(DECLARED, produced_count=9, step_tool_name="s",
                                  operations=HEAD_5) == []

    def test_a_step_declaring_nothing_records_nothing(self):
        assert unfilled_slot_rows([], produced_count=0, step_tool_name="s",
                                  operations=HEAD_5) == []

    def test_zero_survivors_accounts_for_every_slot(self):
        rows = unfilled_slot_rows(DECLARED, produced_count=0, step_tool_name="s",
                                  operations=HEAD_5)
        assert [r["id"] for r in rows] == DECLARED, (
            "a chain that kept nothing still declared five outputs, and all five are unfilled")

    def test_no_limiting_operation_still_produces_a_usable_cause(self):
        rows = unfilled_slot_rows(DECLARED, produced_count=1, step_tool_name="s", operations=[])
        assert rows and "declared" in rows[0]["cause"], (
            "a rename whose count came from somewhere other than head/tail/sample must still "
            "explain the slot rather than emit an empty cause")


class TestTheDeclarationItself:

    def test_only_pool_mode_declares_slots(self, local_config, isolated_cwd):
        """Outside pool mode the streams are not declared here, so there is nothing to account
        for; returning ids anyway would invent missing rows for a step that never promised them."""
        from biopipelines.panda import Panda

        tool = Panda.__new__(Panda)
        tool.use_pool_mode = False
        tool.pool_outputs = None
        tool.rename = "best"
        assert tool._renamed_output_ids() == []

    def test_an_unknown_count_declares_no_slots(self, local_config, isolated_cwd):
        from biopipelines.panda import Panda

        tool = Panda.__new__(Panda)
        tool.use_pool_mode = True
        tool.pool_outputs = ["a"]
        tool.rename = "best"
        tool._get_predicted_output_count = lambda: None
        assert tool._renamed_output_ids() == [], (
            "with no fixed count the ids are a lazy pattern, and a pattern has no slots to leave "
            "unfilled")
