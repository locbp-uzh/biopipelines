"""FoundryMPNN: ProteinMPNN / LigandMPNN through Foundry's `mpnn`.

What is worth pinning is where it differs from the two wrappers of the original repositories it overlaps: it takes
a JSON config rather than flags, writes one unseparated record for all chains of a design (so the chains are
recovered from the written structure), and has one id scheme, `<structure>_<n>`, that has to agree with the
declared ids under the same `chains` contract. The layouts parsed here are read from upstream's source, not from a
real run.
"""

import json
import os
import subprocess
import sys
from types import SimpleNamespace

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "pipe_scripts"))

from biopipelines.foundry_mpnn import FoundryMPNN
import pipe_foundry_mpnn as pipe  # noqa: E402


@pytest.fixture
def tool(local_config, isolated_cwd):
    def _make(**kwargs):
        made = FoundryMPNN.__new__(FoundryMPNN)
        made.structures_stream = kwargs.pop("structures_stream", [object()])
        defaults = dict(model_type="protein_mpnn", num_sequences=1, num_batches=1, temperature=0.1,
                        seed=None, structure_noise=None, designed_residues=None, fixed_residues=None,
                        designed_chains=None, fixed_chains=None, omit=None, bias=None, chains=None,
                        checkpoint_path=None, is_legacy_weights=True)
        for name, default in defaults.items():
            setattr(made, name, kwargs.pop(name, default))
        assert not kwargs, f"unexpected: {sorted(kwargs)}"
        return made
    return _make


class TestItRefusesWhatItCannotDo:

    def test_an_unknown_model_is_refused(self, tool):
        with pytest.raises(ValueError, match="model_type"):
            FoundryMPNN.validate_params(tool(model_type="solublempnn"))

    @pytest.mark.parametrize("name", ["num_sequences", "num_batches"])
    @pytest.mark.parametrize("bad", [0, -2, True, 1.5])
    def test_a_count_must_be_a_positive_integer(self, tool, name, bad):
        with pytest.raises(ValueError, match=name):
            FoundryMPNN.validate_params(tool(**{name: bad}))

    @pytest.mark.parametrize("bad", [0, -0.1, True])
    def test_the_temperature_must_be_positive(self, tool, bad):
        with pytest.raises(ValueError, match="temperature"):
            FoundryMPNN.validate_params(tool(temperature=bad))

    def test_designed_and_fixed_residues_are_mutually_exclusive(self, tool):
        with pytest.raises(ValueError, match="not both"):
            FoundryMPNN.validate_params(tool(designed_residues=["A1"], fixed_residues=["A2"]))

    def test_designed_and_fixed_chains_are_mutually_exclusive(self, tool):
        with pytest.raises(ValueError, match="not both"):
            FoundryMPNN.validate_params(tool(designed_chains=["A"], fixed_chains=["B"]))

    def test_a_residue_selection_and_a_chain_selection_cannot_be_mixed(self, tool):
        with pytest.raises(ValueError, match="cannot be mixed"):
            FoundryMPNN.validate_params(tool(designed_residues=["A1"], fixed_chains=["B"]))

    def test_a_chain_scope_alone_and_a_residue_scope_alone_are_accepted(self, tool):
        FoundryMPNN.validate_params(tool(designed_chains=["A", "B"]))
        FoundryMPNN.validate_params(tool(fixed_residues=["A35", "B40A"]))

    @pytest.mark.parametrize("bad", [["35"], ["A"], ["A-3"], []])
    def test_a_residue_id_needs_a_chain_and_a_number(self, tool, bad):
        with pytest.raises(ValueError, match="residue ids"):
            FoundryMPNN.validate_params(tool(designed_residues=bad))

    def test_a_chain_id_is_letters_only(self, tool):
        with pytest.raises(ValueError, match="chain ids"):
            FoundryMPNN.validate_params(tool(fixed_chains=["A1"]))

    def test_omit_and_bias_are_checked_for_shape(self, tool):
        with pytest.raises(ValueError, match="omit"):
            FoundryMPNN.validate_params(tool(omit="UNK"))
        with pytest.raises(ValueError, match="bias"):
            FoundryMPNN.validate_params(tool(bias={"ALA": "high"}))
        FoundryMPNN.validate_params(tool(omit=["UNK", "CYS"], bias={"ALA": -1.0}))

    def test_a_bad_chains_spec_is_refused_as_for_the_other_mpnns(self, tool):
        with pytest.raises(ValueError, match="chains"):
            FoundryMPNN.validate_params(tool(chains=["AB"]))

    def test_a_checkpoint_path_cannot_carry_shell_expansion(self, tool):
        with pytest.raises(ValueError, match="checkpoint_path"):
            FoundryMPNN.validate_params(tool(checkpoint_path="`x`.pt"))

    def test_a_selection_may_be_given_as_a_comma_string(self):
        from biopipelines.foundry_mpnn import as_id_list
        assert as_id_list("A35, B40,C52") == ["A35", "B40", "C52"]
        assert as_id_list(None) is None


class TestTheEngineConfig:

    BASE = {"model_type": "ligand_mpnn", "is_legacy_weights": True, "num_sequences": 4, "num_batches": 2,
            "temperature": 0.1, "chains": None}

    def test_the_global_keys_are_upstreams_own_and_both_outputs_are_written(self):
        built = pipe.build_config(self.BASE, "/s/a.pdb", "a", "/out")
        assert built["model_type"] == "ligand_mpnn" and built["is_legacy_weights"] is True
        assert built["out_directory"] == "/out"
        assert built["write_fasta"] is True and built["write_structures"] is True, \
            "the structure is how the chains are told apart, so it is always written"
        assert "checkpoint_path" not in built, "unset, upstream finds its registered weights itself"

    def test_the_input_names_the_structure_and_carries_the_batching(self):
        (entry,) = pipe.build_config(self.BASE, "/s/a.pdb", "a", "/out")["inputs"]
        assert entry["structure_path"] == "/s/a.pdb" and entry["name"] == "a"
        assert (entry["batch_size"], entry["number_of_batches"]) == (4, 2)

    def test_only_what_was_set_is_passed_per_input(self):
        (entry,) = pipe.build_config(self.BASE, "/s/a.pdb", "a", "/out")["inputs"]
        assert entry["temperature"] == 0.1
        for absent in ("seed", "structure_noise", "designed_residues", "omit", "bias"):
            assert absent not in entry

    def test_selections_omit_bias_and_the_checkpoint_are_passed_under_upstreams_spelling(self):
        config = dict(self.BASE, seed=5, structure_noise=0.2, designed_residues=["A1"],
                      omit=["UNK", "CYS"], bias={"ALA": -1.0}, checkpoint_path="/w/x.pt")
        built = pipe.build_config(config, "/s/a.pdb", "a", "/out")
        (entry,) = built["inputs"]
        assert built["checkpoint_path"] == "/w/x.pt"
        assert (entry["seed"], entry["structure_noise"], entry["designed_residues"]) == (5, 0.2, ["A1"])
        assert entry["omit"] == ["UNK", "CYS"] and entry["bias"] == {"ALA": -1.0}

    def test_the_wrapper_config_carries_the_chains_contract_and_leaves_unset_keys_out(self, tool):
        made = tool(chains=["A", "B"], seed=3, omit=["UNK"])
        config = FoundryMPNN._config_dict(made)
        assert config["chains"] == ["A", "B"] and config["seed"] == 3
        assert "bias" not in config and "checkpoint_path" not in config


class TestTheHeader:

    def test_a_current_header_carries_every_metric(self):
        ident, metrics = pipe.parse_header(
            ">p+l_b1_d2, confidence=0.5000, ligand_interface_confidence=0.4000, "
            "sequence_recovery=0.3000, ligand_interface_sequence_recovery=0.2000")
        assert ident == {"name": "p+l", "batch": 1, "design": 2}
        assert metrics == {"overall_confidence": 0.5, "ligand_confidence": 0.4, "seq_rec": 0.3,
                           "ligand_seq_rec": 0.2}

    def test_the_older_header_without_confidences_leaves_them_out(self):
        _ident, metrics = pipe.parse_header(">bb_1_b0_d0, sequence_recovery=0.1234")
        assert metrics == {"seq_rec": 0.1234}

    def test_the_name_may_itself_contain_underscores_and_batch_markers(self):
        ident, _ = pipe.parse_header(">bb_b3_x_1_b0_d4, sequence_recovery=0.1")
        assert ident == {"name": "bb_b3_x_1", "batch": 0, "design": 4}

    def test_a_header_that_is_not_ours_is_skipped(self):
        assert pipe.parse_header(">native, sequence_recovery=1.0")[0] is None

    def test_an_unparsable_value_is_left_out_rather_than_crashing(self):
        _ident, metrics = pipe.parse_header(">a_b0_d0, confidence=nan?, sequence_recovery=0.5")
        assert metrics == {"seq_rec": 0.5}

    def test_a_repeated_header_keeps_the_later_record(self, tmp_path):
        """Upstream appends to <name>.fa, so a retry leaves the earlier attempt's records in front."""
        fa = tmp_path / "a.fa"
        fa.write_text(">a_b0_d0, sequence_recovery=0.1\nAAA\n>a_b0_d0, sequence_recovery=0.1\nCCC\n")
        assert pipe.read_fasta(str(fa)) == [(">a_b0_d0, sequence_recovery=0.1", "CCC")]


def cif(residues):
    """A minimal mmCIF with the `_atom_site` loop upstream's writer produces; one CA per residue."""
    head = ("data_x\n#\n_mpnn_input.name  a\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n"
            "_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n"
            "_atom_site.label_seq_id\n_atom_site.pdbx_PDB_ins_code\n_atom_site.auth_seq_id\n"
            "_atom_site.auth_asym_id\n_atom_site.pdbx_PDB_model_num\n")
    rows = ""
    for index, (chain, number, comp) in enumerate(residues, start=1):
        group = "ATOM" if comp not in ("HEM", "MG") else "HETATM"
        atom = '"O5\'"' if comp == "HEM" else "CA"
        rows += f"{group} {index} {atom} {comp} {chain} {number} ? {number} {chain} 1\n"
    return head + rows + "#\n"


A, B = "ALA", "GLY"


class TestRecoveringTheChains:

    def test_only_standard_amino_acids_count_and_quoted_atom_names_do_not_break_the_read(self):
        text = cif([("A", 1, A), ("A", 2, B), ("A", 3, "HEM"), ("B", 1, A), ("C", 1, "MG")])
        assert pipe.chain_lengths(pipe.parse_atom_site(text)) == [("A", 2), ("B", 1)]

    def test_chains_come_out_in_order_of_appearance(self):
        text = cif([("B", 1, A), ("A", 1, A), ("A", 2, A)])
        assert pipe.chain_lengths(pipe.parse_atom_site(text)) == [("B", 1), ("A", 2)]

    def test_several_atoms_of_one_residue_count_once(self):
        text = cif([("A", 1, A), ("A", 1, A), ("A", 2, A)])
        assert pipe.chain_lengths(pipe.parse_atom_site(text)) == [("A", 2)]

    def test_a_second_model_is_not_counted_again(self):
        text = cif([("A", 1, A)]) .replace("#\n", "") + "ATOM 9 CA ALA A 1 ? 1 A 2\n#\n"
        assert pipe.chain_lengths(pipe.parse_atom_site(text)) == [("A", 1)]

    def test_a_structure_without_the_expected_columns_is_named_not_guessed(self):
        with pytest.raises(ValueError, match="_atom_site"):
            pipe.chain_lengths(pipe.parse_atom_site("data_x\n#\n"))

    def test_the_record_is_cut_at_the_chain_boundaries(self):
        assert pipe.split_chains("MKTAYG", [("A", 4), ("B", 2)]) == {"A": "MKTA", "B": "YG"}

    def test_a_record_that_disagrees_with_the_structure_is_refused_not_cut_anyway(self):
        """A silent wrong split would hand a downstream fold the wrong chains."""
        with pytest.raises(ValueError, match="6 residues but the written structure has 5"):
            pipe.split_chains("MKTAYG", [("A", 4), ("B", 1)])

    def test_none_expects_one_chain(self):
        assert pipe.select_chains({"A": "MK"}, None) == [(None, "MK")]
        with pytest.raises(ValueError, match="chains=\"all\""):
            pipe.select_chains({"A": "MK", "B": "TA"}, None)

    def test_all_keeps_every_chain_and_a_list_keeps_those(self):
        chains = {"A": "MK", "B": "TA"}
        assert pipe.select_chains(chains, "all") == [("A", "MK"), ("B", "TA")]
        assert pipe.select_chains(chains, ["B"]) == [("B", "TA")]

    def test_a_requested_chain_the_backbone_lacks_is_named(self):
        with pytest.raises(ValueError, match=r"\['Z'\]"):
            pipe.select_chains({"A": "MK"}, ["Z"])

    def test_only_a_setting_that_can_emit_several_chains_adds_a_suffix(self):
        assert [pipe.suffix_chains(x) for x in (None, ["A"], ["A", "B"], "all")] == [False, False, True, True]


def _write_run(out, name, designs, residues, batch_size=2, header_extra=""):
    """What `mpnn` leaves for one structure: <name>.fa plus a CIF per design."""
    out.mkdir(parents=True, exist_ok=True)
    fasta = ""
    for index, sequence in enumerate(designs):
        batch, design = divmod(index, batch_size)
        fasta += f">{name}_b{batch}_d{design}, {header_extra}sequence_recovery=0.{index + 1}\n{sequence}\n"
        (out / f"{name}_b{batch}_d{design}.cif").write_text(cif(residues))
    (out / f"{name}.fa").write_text(fasta)


CONFIG = {"num_sequences": 2, "num_batches": 1, "temperature": 0.1, "chains": None}


class TestCollectingAStructure:

    def test_a_single_chain_design_keeps_the_plain_ids(self, tmp_path):
        _write_run(tmp_path, "s", ["MKTA", "GGAA"], [("A", n, A) for n in range(1, 5)])
        sequences, designs, missing = pipe.collect_structure("s", "/in/s.pdb", str(tmp_path), CONFIG)
        assert [r["id"] for r in sequences] == ["s_1", "s_2"] == [d["id"] for d in designs]
        assert [r["sequence"] for r in sequences] == ["MKTA", "GGAA"]
        assert missing == []
        assert sequences[0]["design"] == "s_1" and sequences[0]["chain"] == ""
        assert (sequences[0]["sample"], sequences[0]["T"], sequences[0]["source_pdb"]) == (1, 0.1, "s")
        assert designs[1]["seq_rec"] == 0.2 and designs[0]["n_chains"] == 1

    def test_batches_are_numbered_batch_major(self, tmp_path):
        _write_run(tmp_path, "s", ["AA", "CC", "DD", "EE"], [("A", 1, A), ("A", 2, A)])
        config = dict(CONFIG, num_batches=2)
        sequences, _designs, _missing = pipe.collect_structure("s", "/in/s.pdb", str(tmp_path), config)
        assert [(r["id"], r["sequence"]) for r in sequences] == [
            ("s_1", "AA"), ("s_2", "CC"), ("s_3", "DD"), ("s_4", "EE")]

    def test_a_multi_chain_backbone_under_the_default_is_a_failure_not_a_fused_sequence(self, tmp_path):
        _write_run(tmp_path, "s", ["MKTA", "GGAA"], [("A", 1, A), ("A", 2, A), ("B", 1, A), ("B", 2, A)])
        sequences, designs, missing = pipe.collect_structure("s", "/in/s.pdb", str(tmp_path), CONFIG)
        assert sequences == [] and designs == []
        assert [m["id"] for m in missing] == ["s_1", "s_2"]
        assert all(m["kind"] == "failure" and "chains=" in m["cause"] for m in missing)

    def test_chains_all_gives_one_row_per_chain_under_a_chain_suffix(self, tmp_path):
        _write_run(tmp_path, "s", ["MKTA", "GGAA"], [("A", 1, A), ("A", 2, A), ("B", 1, A), ("B", 2, A)])
        sequences, designs, missing = pipe.collect_structure(
            "s", "/in/s.pdb", str(tmp_path), dict(CONFIG, chains="all"))
        assert [r["id"] for r in sequences] == ["s_1_A", "s_1_B", "s_2_A", "s_2_B"]
        assert [r["sequence"] for r in sequences] == ["MK", "TA", "GG", "AA"]
        assert [r["chain"] for r in sequences] == ["A", "B", "A", "B"]
        assert {r["design"] for r in sequences} == {"s_1", "s_2"}
        assert [d["n_chains"] for d in designs] == [2, 2] and missing == []

    def test_one_requested_chain_keeps_the_plain_ids(self, tmp_path):
        _write_run(tmp_path, "s", ["MKTA", "GGAA"], [("A", 1, A), ("A", 2, A), ("B", 1, A), ("B", 2, A)])
        sequences, _designs, _missing = pipe.collect_structure(
            "s", "/in/s.pdb", str(tmp_path), dict(CONFIG, chains=["B"]))
        assert [(r["id"], r["chain"], r["sequence"]) for r in sequences] == [
            ("s_1", "B", "TA"), ("s_2", "B", "AA")]

    def test_a_ligand_in_the_structure_is_context_not_a_chain(self, tmp_path):
        residues = [("A", 1, A), ("A", 2, A), ("B", 1, "HEM"), ("C", 1, "MG")]
        _write_run(tmp_path, "s", ["MK", "TA"], residues)
        sequences, _designs, missing = pipe.collect_structure("s", "/in/s.pdb", str(tmp_path), CONFIG)
        assert [r["id"] for r in sequences] == ["s_1", "s_2"] and missing == []

    def test_a_length_disagreement_fails_that_design_only(self, tmp_path):
        _write_run(tmp_path, "s", ["MKTA", "GGAAAA"], [("A", n, A) for n in range(1, 5)])
        sequences, designs, missing = pipe.collect_structure("s", "/in/s.pdb", str(tmp_path), CONFIG)
        assert [r["id"] for r in sequences] == ["s_1"] and [d["id"] for d in designs] == ["s_1"]
        assert [m["id"] for m in missing] == ["s_2"] and "6 residues" in missing[0]["cause"]

    def test_a_design_with_no_structure_cannot_have_its_chains_told_apart(self, tmp_path):
        _write_run(tmp_path, "s", ["MKTA", "GGAA"], [("A", n, A) for n in range(1, 5)])
        os.remove(tmp_path / "s_b0_d1.cif")
        sequences, _designs, missing = pipe.collect_structure("s", "/in/s.pdb", str(tmp_path), CONFIG)
        assert [r["id"] for r in sequences] == ["s_1"]
        assert [m["id"] for m in missing] == ["s_2"] and "no structure" in missing[0]["cause"]

    def test_a_design_upstream_never_wrote_is_missing_by_its_declared_id(self, tmp_path):
        _write_run(tmp_path, "s", ["MKTA"], [("A", n, A) for n in range(1, 5)])
        _sequences, _designs, missing = pipe.collect_structure("s", "/in/s.pdb", str(tmp_path), CONFIG)
        assert [m["id"] for m in missing] == ["s_2"]

    def test_no_fasta_at_all_fails_every_declared_design(self, tmp_path):
        sequences, designs, missing = pipe.collect_structure("s", "/in/s.pdb", str(tmp_path), CONFIG)
        assert sequences == designs == [] and [m["id"] for m in missing] == ["s_1", "s_2"]

    def test_a_record_for_another_structure_is_ignored(self, tmp_path):
        _write_run(tmp_path, "s", ["MKTA", "GGAA"], [("A", n, A) for n in range(1, 5)])
        with open(tmp_path / "s.fa", "a") as handle:
            handle.write(">other_b0_d0, sequence_recovery=0.9\nWWWW\n")
        sequences, _designs, _missing = pipe.collect_structure("s", "/in/s.pdb", str(tmp_path), CONFIG)
        assert [r["id"] for r in sequences] == ["s_1", "s_2"]

    def test_the_confidences_of_a_newer_release_reach_the_rows(self, tmp_path):
        _write_run(tmp_path, "s", ["MKTA", "GGAA"], [("A", n, A) for n in range(1, 5)],
                   header_extra="confidence=0.7000, ligand_interface_confidence=0.6000, ")
        sequences, designs, _missing = pipe.collect_structure("s", "/in/s.pdb", str(tmp_path), CONFIG)
        assert sequences[0]["overall_confidence"] == 0.7 and designs[0]["ligand_confidence"] == 0.6


class TestTheRunStep:
    """The whole script with `mpnn` replaced by a stand-in that writes what upstream writes."""

    @staticmethod
    def _stream(tmp_path, names):
        from biopipelines.datastream import DataStream
        files = []
        for name in names:
            path = tmp_path / f"{name}.pdb"
            path.write_text("ATOM\n")
            files.append(str(path))
        table = tmp_path / "structures_map.csv"
        table.write_text("id,file\n" + "".join(f"{n},{f}\n" for n, f in zip(names, files)))
        stream = DataStream(name="structures", ids=names, files=files, map_table=str(table), format="pdb")
        path = tmp_path / "structures.json"
        stream.save_json(str(path))
        return path

    def _run(self, tmp_path, monkeypatch, names, behaviour, config=None):
        import yaml
        monkeypatch.setattr(pipe, "run_mpnn", behaviour)
        config_yaml = tmp_path / "config.yaml"
        config_yaml.write_text(yaml.safe_dump(config or dict(
            CONFIG, model_type="protein_mpnn", is_legacy_weights=True)))
        out = tmp_path / "out"
        monkeypatch.setattr(sys, "argv", [
            "pipe", "--config-yaml", str(config_yaml),
            "--structures-json", str(self._stream(tmp_path, names)),
            "--out-dir", str(out / "001_FoundryMPNN" / "execution" / "mpnn"),
            "--sequences-csv", str(out / "sequences" / "sequences.csv"),
            "--designs-csv", str(out / "designs" / "designs.csv"),
            "--local-missing-csv", str(out / "001_FoundryMPNN" / "execution" / "local_missing.csv")])
        try:
            pipe.main()
            code = 0
        except SystemExit as stop:
            code = stop.code
        return out, code

    @staticmethod
    def _fake(failing=()):
        def behaviour(prefix, config_path):
            config = json.load(open(config_path))
            (entry,) = config["inputs"]
            if entry["name"] in failing:
                return "RuntimeError: could not parse the structure"
            from pathlib import Path
            _write_run(Path(config["out_directory"]), entry["name"], ["MKTA", "GGAA"],
                       [("A", n, A) for n in range(1, 5)])
            return None
        return behaviour

    def test_a_bad_structure_fails_its_own_designs_and_not_the_step(self, tmp_path, monkeypatch):
        import pandas as pd
        out, code = self._run(tmp_path, monkeypatch, ["good", "bad"], self._fake(failing=("bad",)))
        assert code == 0
        assert pd.read_csv(out / "sequences" / "sequences.csv")["id"].tolist() == ["good_1", "good_2"]
        assert pd.read_csv(out / "designs" / "designs.csv")["id"].tolist() == ["good_1", "good_2"]
        missing = pd.read_csv(out / "001_FoundryMPNN" / "execution" / "local_missing.csv")
        assert missing["id"].tolist() == ["bad_1", "bad_2"]
        assert set(missing["removed_by"]) == {"001_FoundryMPNN"} and set(missing["kind"]) == {"failure"}
        assert "could not parse" in missing["cause"].iloc[0]

    def test_every_structure_failing_exits_non_zero(self, tmp_path, monkeypatch):
        _out, code = self._run(tmp_path, monkeypatch, ["bad"], self._fake(failing=("bad",)))
        assert code == 1

    def test_a_stale_fasta_from_an_earlier_attempt_is_cleared_before_the_run(self, tmp_path, monkeypatch):
        stale = tmp_path / "out" / "001_FoundryMPNN" / "execution" / "mpnn"
        stale.mkdir(parents=True)
        (stale / "good.fa").write_text(">good_b0_d0, sequence_recovery=0.9\nWWWW\n")
        seen = {}

        def behaviour(prefix, config_path):
            seen["stale_present"] = (stale / "good.fa").exists()
            return self._fake()(prefix, config_path)
        self._run(tmp_path, monkeypatch, ["good"], behaviour)
        assert seen["stale_present"] is False

    def test_the_per_structure_config_is_written_where_the_run_can_be_reproduced(self, tmp_path, monkeypatch):
        self._run(tmp_path, monkeypatch, ["good"], self._fake())
        config = tmp_path / "out" / "001_FoundryMPNN" / "execution" / "mpnn" / "_configs" / "good.json"
        assert json.loads(config.read_text())["inputs"][0]["name"] == "good"


class TestTheDeclaredOutputs:

    @staticmethod
    def _built(local_config, isolated_cwd, **kwargs):
        from biopipelines.mock import Mock
        from biopipelines.pipeline import Pipeline
        with Pipeline(project="TestSuite", job="fmpnn", description="x", on_the_fly=False,
                      local_output=True, config="local") as pipeline:
            structures = Mock(ids=["s1", "s2"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
            FoundryMPNN(structures.streams.structures, **kwargs)
            tool = pipeline.tools[-1]
            script_path = pipeline.save()
        from pathlib import Path
        steps = sorted(Path(script_path).parent.glob("*_FoundryMPNN.sh"))
        return tool, "\n".join(step.read_text(encoding="utf-8", errors="replace") for step in steps)

    def test_designs_are_numbered_over_batches_and_sequences(self, local_config, isolated_cwd):
        tool, _script = self._built(local_config, isolated_cwd, num_sequences=3, num_batches=2)
        out = tool.get_output_files()
        assert out["designs"].ids == ["s1_<1..6>", "s2_<1..6>"]
        assert out["sequences"].ids == out["designs"].ids

    @pytest.mark.parametrize("chains,expected", [
        (None, ["s1_<1..2>", "s2_<1..2>"]),
        (["A"], ["s1_<1..2>", "s2_<1..2>"]),
        ("all", ["s1_<1..2>[_<?>]", "s2_<1..2>[_<?>]"]),
    ])
    def test_sequence_ids_follow_the_same_chains_contract_as_the_other_mpnns(
            self, local_config, isolated_cwd, chains, expected):
        from biopipelines.chain_rows import chain_row_ids
        tool, _script = self._built(local_config, isolated_cwd, num_sequences=2, chains=chains)
        out = tool.get_output_files()
        assert out["sequences"].ids == expected == chain_row_ids(out["designs"].ids, tool.chains)

    def test_several_chains_add_a_chain_suffix_to_the_declared_ids(self, local_config, isolated_cwd):
        tool, _script = self._built(local_config, isolated_cwd, chains=["A", "B"])
        assert all("A" in i and "B" in i for i in tool.get_output_files()["sequences"].ids)

    def test_the_streams_are_value_based_csvs_and_the_tables_are_declared(self, local_config, isolated_cwd):
        tool, _script = self._built(local_config, isolated_cwd)
        out = tool.get_output_files()
        for name in ("sequences", "designs"):
            assert out[name].format == "csv" and list(out[name].files) == []
        assert out["tables"]["sequences"].info.columns == [
            "id", "design", "structures.id", "source_pdb", "chain", "sequence", "sample", "T", "seed",
            "overall_confidence", "ligand_confidence", "seq_rec", "ligand_seq_rec"]
        assert out["tables"]["designs"].info.columns[:4] == ["id", "structures.id", "source_pdb", "n_chains"]
        assert out["tables"]["missing"].info.columns == ["id", "removed_by", "kind", "cause"]

    def test_the_sequence_table_columns_are_the_ones_the_pipe_writes(self):
        from biopipelines.foundry_mpnn import FoundryMPNN as tool_cls
        assert pipe.SEQUENCE_COLUMNS == [
            "id", "design", "structures.id", "source_pdb", "chain", "sequence", "sample", "T", "seed",
            "overall_confidence", "ligand_confidence", "seq_rec", "ligand_seq_rec"]
        assert pipe.DESIGN_COLUMNS == ["id", "structures.id", "source_pdb", "n_chains", "overall_confidence",
                                       "ligand_confidence", "seq_rec", "ligand_seq_rec"]
        assert tool_cls.TOOL_NAME == "FoundryMPNN"

    def test_every_flag_the_script_sends_is_known_to_the_pipe_script(self, local_config, isolated_cwd):
        import re
        _tool, script = self._built(local_config, isolated_cwd)
        block = script.split("pipe_foundry_mpnn.py")[1].split("\n\n")[0]
        sent = set(re.findall(r"--[a-z-]+", block))
        with open(os.path.join(REPO_ROOT, "pipe_scripts", "pipe_foundry_mpnn.py")) as handle:
            known = set(re.findall(r'add_argument\("(--[a-z-]+)"', handle.read()))
        assert sent <= known

    def test_the_pipe_runs_on_the_host_and_is_handed_the_container_prefix(self, local_config, isolated_cwd):
        _tool, script = self._built(local_config, isolated_cwd)
        assert '--container-prefix ""' in script
        assert "FOUNDRY_CHECKPOINT_DIRS" in script


class TestTheInstall:
    """The foundry environment is shared with RFdiffusion3; the weights and the shortcut are not."""

    @pytest.fixture(autouse=True)
    def _environment(self, local_config, monkeypatch):
        monkeypatch.setattr(FoundryMPNN, "_install_env", classmethod(lambda cls, manager: "foundry"))

    FOLDERS = {"biopipelines": "/repo", "FoundryMPNN": "/data/mpnn"}

    def test_it_downloads_both_mpnn_checkpoints_into_its_own_folder(self):
        script = FoundryMPNN._install_script(self.FOLDERS, env_manager="mamba")
        assert 'foundry install proteinmpnn ligandmpnn --checkpoint-dir "/data/mpnn"' in script
        assert "rfd3" not in script

    def test_the_shortcut_requires_both_checkpoints_not_just_the_environment(self):
        shortcut = FoundryMPNN._install_script(self.FOLDERS, env_manager="mamba").split(
            "already installed, skipping")[0]
        assert "proteinmpnn_*.pt" in shortcut and "ligandmpnn_*.pt" in shortcut

    def test_a_config_without_the_folder_is_refused_by_name(self):
        with pytest.raises(ValueError, match="repositories.FoundryMPNN"):
            FoundryMPNN._install_script({"biopipelines": "/repo"}, env_manager="mamba")
