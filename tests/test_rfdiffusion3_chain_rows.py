"""RFdiffusion3's `sequences` rows are per chain, and its declared ids have to say so.

`docs/design/multichain-sequences.md` sets the repo-wide rule: one `sequences` row is one polymer
chain. Both MPNN producers were given declared ids matching that. RFdiffusion3 was not, and the
mismatch stayed unreachable only because the wrapper rejected every contig containing `/` — so a
multi-chain design could not be expressed at all.

Once `/0` became the accepted chain-break token, a two-chain contig is the ordinary binder-design
spelling, and the mismatch became reachable through the documented API. Its failure is silent:
`select_ids(['bb_1'], ['bb_1_A','bb_1_B'])` returns nothing, so a Boltz2 chained off
`.streams.sequences` generates zero configs, boltz exits 0 over an empty directory, and the step
produces no structures while reporting no error of its own — the same signature as the incident
that cost four binder-design campaigns.
"""

import pytest

from biopipelines import id_patterns


@pytest.fixture
def rfd3(local_config, isolated_cwd):
    def _make(contig, **kwargs):
        from biopipelines.rfdiffusion3 import RFdiffusion3
        tool = RFdiffusion3.__new__(RFdiffusion3)
        tool.contig = contig
        tool.contig_reference = kwargs.get("contig_reference")
        return tool
    return _make


class TestTheDeclaredIdsMatchTheRowsWritten:
    """The runtime writer is `pipe_rfdiffusion3_postprocess.py`: `<structure>_<chain>` when a
    design has several chains, the bare id when it has one."""

    def test_a_chain_break_declares_a_lazy_chain_suffix(self, rfd3):
        tool = rfd3("A50-100,80-100,/0,A1-50")
        declared = tool._sequence_row_ids(["bb_1", "bb_2"])
        assert declared == ["bb_1[_<?>]", "bb_2[_<?>]"]

    def test_the_lazy_form_selects_the_rows_that_get_written(self, rfd3):
        declared = rfd3("A50-100,/0,A1-50")._sequence_row_ids(["bb_1"])
        assert id_patterns.select_ids(declared, ["bb_1_A", "bb_1_B"]) == ["bb_1_A", "bb_1_B"], (
            "this is the selection a downstream consumer performs; empty means zero configs")

    def test_it_also_covers_a_single_chain_design(self, rfd3):
        """The writer omits the suffix when there is one chain, and the same contig can do both."""
        declared = rfd3("A50-100,/0,A1-50")._sequence_row_ids(["bb_1"])
        assert id_patterns.select_ids(declared, ["bb_1"]) == ["bb_1"]

    def test_a_single_chain_contig_keeps_its_concrete_ids(self, rfd3):
        """No suffix is ever written for these, so a lazy id would trade precision for nothing."""
        assert rfd3("A50-100,80-100")._sequence_row_ids(["bb_1", "bb_2"]) == ["bb_1", "bb_2"]

    def test_a_runtime_contig_is_treated_as_multi_chain(self, rfd3):
        """A contig arriving per structure cannot be read at configuration time.

        Guessing single-chain there is the unsafe direction: it silently selects nothing, where
        guessing multi-chain costs only a lazy id on a stream that turns out to have one.
        """
        tool = rfd3("", contig_reference="tbl.contigs")
        assert tool._sequence_row_ids(["bb_1"]) == ["bb_1[_<?>]"]


class TestTheRegressionItself:

    def test_a_bare_structure_id_selects_nothing_against_chain_rows(self):
        """Why the fix is needed at all, pinned so nobody reintroduces the old declaration."""
        assert id_patterns.select_ids(["bb_1"], ["bb_1_A", "bb_1_B"]) == []


class TestJsonConfigContigs:
    """`json_config` is a documented way to give contigs, so its chain breaks count too."""

    def test_a_chain_break_inside_json_config_declares_chain_rows(self, rfd3):
        tool = rfd3("")
        tool.json_config = {"design_1": {"contig": "50-80,/0,A1-100"}}
        assert tool._sequence_row_ids(["bb_1"]) == ["bb_1[_<?>]"]

    def test_a_json_string_config_is_read_as_well(self, rfd3):
        tool = rfd3("")
        tool.json_config = '{"design_1": {"contig": "50-80,/0,A1-100"}}'
        assert tool._contig_has_chain_break()

    def test_a_single_chain_json_config_keeps_concrete_ids(self, rfd3):
        tool = rfd3("")
        tool.json_config = {"design_1": {"contig": "50-80"}}
        assert tool._sequence_row_ids(["bb_1"]) == ["bb_1"]
