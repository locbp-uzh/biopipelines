"""BioEmu2: what the wrapper writes for `bioemu2 sample-equilibrium`, reads back from it, and declares.

Nothing here runs the model or the CLI. The upstream flags, the system YAML, the MSA cache entry and the
physicality report are pinned from its source (github.com/microsoft/bioemu2: `cli.py`, `physicality.py`,
`msa/colabfold.py`, `checkpoints.py` and `docs/system.md`).

Not covered, because it needs the model or a cluster: that the sampler accepts the systems this builds, the
residue and atom naming of a ligand in a sampled structure, and that the MSA cache key below is the one the
installed release computes for a monomer. The trajectory split is run on a real single-chain `topology.cif`
when mdtraj is installed and skipped otherwise.
"""

import hashlib
import json
import os
import re
import shutil
import subprocess
import sys

import pandas as pd
import pytest
import yaml

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(REPO_ROOT, "pipe_scripts")
sys.path.insert(0, SCRIPTS)

from biopipelines.bioemu2 import (BioEmu2, COVALENT_LINKAGE_KEYS, MIN_PYTHON, PHYSICALITY_MODES,  # noqa: E402
                                  SAMPLERS)
from biopipelines.combinatorics import Bundle, CombinatoricsConfig, generate_combinatorics_config  # noqa: E402
from biopipelines.datastream import DataStream  # noqa: E402
import pipe_bioemu2 as pipe  # noqa: E402
import pipe_bioemu2_weights as weights  # noqa: E402

PIPE_SCRIPT = os.path.join(SCRIPTS, "pipe_bioemu2.py")
PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"


class Stream:
    """Stands in for a DataStream where only its truth value and length matter."""

    def __init__(self, n=1, fmt=None):
        self.n, self.format = n, fmt

    def __len__(self):
        return self.n

    def __bool__(self):
        return self.n > 0


@pytest.fixture
def tool(local_config, isolated_cwd):
    """A BioEmu2 with only the attributes validate_params and the section builders read."""
    defaults = dict(sequences=object(), ligands=None, msas=None, sequences_stream=Stream(),
                    ligands_stream=None, msas_stream=None, num_samples=100, batch_size=16, seed=None,
                    sampler="edm", contacts=None, covalent_linkage=None, physicality_filter="high_quality")

    def _make(**kwargs):
        made = BioEmu2.__new__(BioEmu2)
        for name, value in {**defaults, **kwargs}.items():
            setattr(made, name, value)
        return made
    return _make


CONTACT = {"token1": {"asym_id": "A", "res_id": 3}, "token2": {"asym_id": "B", "atom_name": "C1"},
           "max_distance": 6.0}
LINK = {"chain": "A", "position": 3, "protein_atom": "SG", "ligand_atom": "C1"}


class TestItRefusesWhatItCannotDo:
    """Each of these would otherwise run and return something quietly wrong."""

    def test_empty_sequences_are_refused(self, tool):
        with pytest.raises(ValueError, match="sequences"):
            BioEmu2.validate_params(tool(sequences_stream=Stream(0)))

    def test_ligands_given_but_empty_are_refused_instead_of_running_apo(self, tool):
        with pytest.raises(ValueError, match="apo"):
            BioEmu2.validate_params(tool(ligands=object(), ligands_stream=Stream(0)))

    @pytest.mark.parametrize("name", ["num_samples", "batch_size"])
    @pytest.mark.parametrize("value", [0, -1, 2.5, True, "10"])
    def test_bad_counts_are_refused(self, tool, name, value):
        with pytest.raises(ValueError, match=name):
            BioEmu2.validate_params(tool(**{name: value}))

    @pytest.mark.parametrize("seed", [1.5, True, "7"])
    def test_a_seed_that_is_not_an_integer_is_refused(self, tool, seed):
        with pytest.raises(ValueError, match="seed"):
            BioEmu2.validate_params(tool(seed=seed))

    def test_a_seed_and_none_are_accepted(self, tool):
        BioEmu2.validate_params(tool(seed=0))
        BioEmu2.validate_params(tool(seed=None))

    def test_an_unknown_sampler_is_refused_with_the_accepted_names(self, tool):
        with pytest.raises(ValueError, match="edm_fkc"):
            BioEmu2.validate_params(tool(sampler="euler"))

    def test_every_sampler_is_accepted(self, tool):
        assert SAMPLERS == ("edm", "edm_gd", "edm_fk", "edm_fkc")
        for name in SAMPLERS:
            BioEmu2.validate_params(tool(sampler=name))

    @pytest.mark.parametrize("name", ["medium_quality", "low_quality"])
    def test_the_physicality_modes_upstream_does_not_have_are_refused(self, tool, name):
        """Passed through, the CLI rejects them after the job has queued."""
        assert name not in PHYSICALITY_MODES
        with pytest.raises(ValueError, match="physicality_filter"):
            BioEmu2.validate_params(tool(physicality_filter=name))

    def test_the_modes_upstream_defines_are_the_ones_offered(self, tool):
        assert set(PHYSICALITY_MODES) == {"no_filter", "no_violations", "high_quality",
                                          "high_quality_per_residue"}
        for name in PHYSICALITY_MODES:
            BioEmu2.validate_params(tool(physicality_filter=name))

    def test_an_alignment_stream_that_is_not_a3m_or_csv_is_refused(self, tool):
        with pytest.raises(ValueError, match="a3m"):
            BioEmu2.validate_params(tool(msas=object(), msas_stream=Stream(1, "pdb")))

    @pytest.mark.parametrize("fmt", ["a3m", "csv", "A3M"])
    def test_a3m_and_boltz_csv_alignments_are_accepted(self, tool, fmt):
        BioEmu2.validate_params(tool(msas=object(), msas_stream=Stream(1, fmt)))


class TestContactValidation:

    def test_a_complete_contact_is_accepted(self, tool):
        BioEmu2.validate_params(tool(contacts=[CONTACT]))

    def test_max_distance_is_optional(self, tool):
        BioEmu2.validate_params(tool(contacts=[{k: v for k, v in CONTACT.items() if k != "max_distance"}]))

    @pytest.mark.parametrize("contacts", [[], {}, "x", [1]])
    def test_a_contacts_value_that_is_not_a_list_of_dicts_is_refused(self, tool, contacts):
        with pytest.raises(ValueError, match="contacts"):
            BioEmu2.validate_params(tool(contacts=contacts))

    @pytest.mark.parametrize("missing", ["token1", "token2"])
    def test_a_missing_token_is_named(self, tool, missing):
        with pytest.raises(ValueError, match=missing):
            BioEmu2.validate_params(tool(contacts=[{k: v for k, v in CONTACT.items() if k != missing}]))

    def test_a_token_needs_a_residue_or_an_atom(self, tool):
        bad = {**CONTACT, "token1": {"asym_id": "A"}}
        with pytest.raises(ValueError, match="res_id or an atom_name"):
            BioEmu2.validate_params(tool(contacts=[bad]))

    def test_a_token_key_upstream_does_not_know_is_refused(self, tool):
        bad = {**CONTACT, "token1": {"asym_id": "A", "res_id": 3, "chain": "A"}}
        with pytest.raises(ValueError, match="chain"):
            BioEmu2.validate_params(tool(contacts=[bad]))

    @pytest.mark.parametrize("distance", [0, -3.0, "6", True])
    def test_a_bad_max_distance_is_refused(self, tool, distance):
        with pytest.raises(ValueError, match="max_distance"):
            BioEmu2.validate_params(tool(contacts=[{**CONTACT, "max_distance": distance}]))


class TestCovalentLinkageValidation:

    @staticmethod
    def _with_ligand(tool, **kwargs):
        return tool(ligands=object(), ligands_stream=Stream(), **kwargs)

    def test_a_complete_linkage_is_accepted(self, tool):
        BioEmu2.validate_params(self._with_ligand(tool, covalent_linkage=LINK))

    @pytest.mark.parametrize("missing", COVALENT_LINKAGE_KEYS)
    def test_a_missing_key_is_named(self, tool, missing):
        link = {k: v for k, v in LINK.items() if k != missing}
        with pytest.raises(ValueError, match=missing):
            BioEmu2.validate_params(self._with_ligand(tool, covalent_linkage=link))

    def test_an_unknown_key_is_refused(self, tool):
        with pytest.raises(ValueError, match="ligand_chain"):
            BioEmu2.validate_params(self._with_ligand(tool, covalent_linkage={**LINK, "ligand_chain": "B"}))

    @pytest.mark.parametrize("position", [0, -1, "3", 2.0, True])
    def test_positions_count_from_one(self, tool, position):
        with pytest.raises(ValueError, match="position"):
            BioEmu2.validate_params(self._with_ligand(tool, covalent_linkage={**LINK, "position": position}))

    def test_a_linkage_without_ligands_is_refused(self, tool):
        with pytest.raises(ValueError, match="requires ligands"):
            BioEmu2.validate_params(tool(covalent_linkage=LINK))

    def test_a_linkage_that_is_not_a_dict_is_refused(self, tool):
        with pytest.raises(ValueError, match="dict"):
            BioEmu2.validate_params(self._with_ligand(tool, covalent_linkage=[LINK]))


class TestInputs:
    """Streams in, never bare values; the same normalization every other multi-axis tool applies."""

    @staticmethod
    def _inputs():
        from biopipelines.ligand import Ligand
        from biopipelines.sequence import Sequence
        return Sequence(seq=PROTEIN, type="protein", ids="p1"), Ligand(smiles="CCO", ids="eth")

    def test_a_bare_sequence_string_is_refused(self, local_config, isolated_cwd):
        with pytest.raises(ValueError, match="Cannot resolve str"):
            BioEmu2(sequences=PROTEIN)

    def test_a_bare_ligand_smiles_is_refused(self, local_config, isolated_cwd, new_pipeline):
        with new_pipeline("bioemu2_inputs"):
            seqs, _ = self._inputs()
            with pytest.raises(ValueError, match="Cannot resolve str"):
                BioEmu2(sequences=seqs, ligands="CCO")

    def test_a_tool_output_and_its_stream_are_both_accepted(self, local_config, isolated_cwd, new_pipeline):
        with new_pipeline("bioemu2_inputs"):
            seqs, lig = self._inputs()
            BioEmu2(sequences=seqs, ligands=lig)
            BioEmu2(sequences=seqs.streams.sequences, ligands=lig.streams.compounds)

    def test_each_and_bundle_wrappers_are_accepted(self, local_config, isolated_cwd, new_pipeline):
        with new_pipeline("bioemu2_inputs"):
            seqs, lig = self._inputs()
            BioEmu2(sequences=seqs, ligands=Bundle(lig))

    def test_the_old_precomputed_msa_spelling_still_binds_to_msas(self, local_config, isolated_cwd, new_pipeline):
        from biopipelines.mock import Mock
        with new_pipeline("bioemu2_inputs") as pipeline:
            seqs, _ = self._inputs()
            msas = Mock(ids=["p1"], streams={"msas": {"format": "a3m", "file": "<id>.a3m"}})
            BioEmu2(sequences=seqs, precomputed_msa=msas.streams.msas)
            made = pipeline.tools[-1]
        assert made.msas is not None and made.msas_stream is not None
        assert "precomputed_msa" in BioEmu2.DEPRECATED_ALIASES
        assert BioEmu2.PARAMETER_ALIASES["precomputed_msa"] == "msas"


class TestTheGeneratedScript:

    @staticmethod
    def _script(new_pipeline, ligands=True, with_msas=False, **kwargs):
        from pathlib import Path
        from biopipelines.ligand import Ligand
        from biopipelines.mock import Mock
        from biopipelines.sequence import Sequence
        pipeline = new_pipeline("bioemu2_script")
        with pipeline:
            seqs = Sequence(seq=PROTEIN, type="protein", ids="p1")
            options = dict(kwargs)
            if ligands:
                options["ligands"] = Ligand(smiles="CCO", ids="eth")
            if with_msas:
                options["msas"] = Mock(ids=["p1"], streams={"msas": {"format": "a3m", "file": "<id>.a3m"}}
                                       ).streams.msas
            BioEmu2(sequences=seqs, **options)
            made = pipeline.tools[-1]
            script_path = pipeline.save()
        steps = sorted(Path(script_path).parent.glob("*_BioEmu2.sh"))
        return made, "\n".join(step.read_text(encoding="utf-8", errors="replace") for step in steps)

    @staticmethod
    def _helper_flags(script):
        block = script.split("pipe_bioemu2.py")[1].split("\n\n")[0]
        return set(re.findall(r"(--[a-z0-9-]+)", block))

    @staticmethod
    def _declared_flags():
        with open(PIPE_SCRIPT) as handle:
            return set(re.findall(r'add_argument\("(--[a-z0-9-]+)"', handle.read()))

    def test_every_flag_the_wrapper_sends_is_known_to_the_pipe_script(self, local_config, isolated_cwd, new_pipeline):
        """A flag the wrapper sends and the pipe does not read would fail the step on the cluster."""
        _made, script = self._script(new_pipeline, with_msas=True)
        sent = self._helper_flags(script)
        assert sent and sent <= self._declared_flags(), sent - self._declared_flags()

    def test_every_required_flag_of_the_pipe_script_is_sent(self, local_config, isolated_cwd, new_pipeline):
        _made, script = self._script(new_pipeline)
        with open(PIPE_SCRIPT) as handle:
            required = set(re.findall(r'add_argument\("(--[a-z0-9-]+)"[^)]*required=True', handle.read()))
        assert required <= self._helper_flags(script), required - self._helper_flags(script)

    def test_the_msa_flags_are_sent_only_with_msas(self, local_config, isolated_cwd, new_pipeline):
        _made, plain = self._script(new_pipeline)
        assert "--msas-json" not in plain and "--msa-cache-dir" not in plain
        _made, with_msas = self._script(new_pipeline, with_msas=True)
        assert "--msas-json" in with_msas and "--msa-cache-dir" in with_msas

    def test_parameters_travel_in_the_yaml_not_on_the_command_line(self, local_config, isolated_cwd, new_pipeline):
        made, script = self._script(new_pipeline, seed=5, sampler="edm_gd", num_samples=7, batch_size=3,
                                    physicality_filter="no_filter", contacts=[CONTACT],
                                    covalent_linkage=LINK)
        for fragment in ("--num-samples", "--batch-size", "--sampler", "--seed", "--contacts", "--covalent"):
            assert fragment not in script
        with open(made.config_yaml) as handle:
            config = yaml.safe_load(handle)
        assert config == {"num_samples": 7, "batch_size": 3, "seed": 5, "sampler": "edm_gd",
                          "physicality_filter": "no_filter", "contacts": [CONTACT], "covalent_linkage": LINK}

    def test_the_defaults_are_written_with_empty_constraints(self, local_config, isolated_cwd, new_pipeline):
        made, _script = self._script(new_pipeline, ligands=False)
        with open(made.config_yaml) as handle:
            config = yaml.safe_load(handle)
        assert config["num_samples"] == 100 and config["batch_size"] == 16
        assert config["contacts"] is None and config["covalent_linkage"] is None and config["seed"] is None

    def test_the_combinatorics_config_names_both_axes(self, local_config, isolated_cwd, new_pipeline):
        made, _script = self._script(new_pipeline)
        config = CombinatoricsConfig.load(made.combinatorics_config_file)
        assert set(config.axes) == {"sequences", "ligands"}
        assert config.predicted_ids == ["p1+eth"]

    def test_an_apo_run_has_only_the_sequences_axis(self, local_config, isolated_cwd, new_pipeline):
        made, _script = self._script(new_pipeline, ligands=False)
        assert set(CombinatoricsConfig.load(made.combinatorics_config_file).axes) == {"sequences"}

    def test_the_local_failures_are_merged_into_the_missing_table_by_the_propagation_step(
            self, local_config, isolated_cwd, new_pipeline):
        made, script = self._script(new_pipeline)
        assert "pipe_propagate_missing.py" in script
        assert f'--local-missing "{made.local_missing_csv}"' in script
        assert f'--local-missing-csv "{made.local_missing_csv}"' in script
        assert "--upstream-missing" not in script
        assert made.local_missing_csv.endswith(os.path.join("_execution", "local_missing.csv"))

    def test_the_propagation_follows_the_helper_so_its_exit_status_is_captured(
            self, local_config, isolated_cwd, new_pipeline):
        _made, script = self._script(new_pipeline)
        assert script.index("pipe_bioemu2.py") < script.index("BP_MAIN_RC=$?") < script.index("pipe_propagate_missing.py")

    def test_the_helper_runs_on_the_host_and_only_the_container_prefix_is_passed_down(
            self, local_config, isolated_cwd, new_pipeline):
        _made, script = self._script(new_pipeline)
        assert '--container-prefix ""' in script
        assert "apptainer" not in script
        assert "WARNING: a container is configured" not in script

    def test_a_configured_container_reaches_the_helper_as_a_prefix(self, local_config, isolated_cwd, new_pipeline,
                                                                  monkeypatch):
        monkeypatch.setattr(BioEmu2, "container_prefix", lambda self: "apptainer exec --nv x.sif ")
        _made, script = self._script(new_pipeline)
        assert '--container-prefix "apptainer exec --nv x.sif "' in script
        assert "apptainer exec --nv x.sif python" not in script

    def test_the_work_directories_live_under_execution_not_extras(self, local_config, isolated_cwd, new_pipeline):
        made, script = self._script(new_pipeline)
        assert made.runs_folder.endswith(os.path.join("_execution", "runs"))
        assert "_extras" not in script


class TestTheDeclaredOutputs:

    @staticmethod
    def _tool(new_pipeline, sequences=("p1", "p2"), ligands=("l1",), num_samples=4, bundle=False):
        from biopipelines.ligand import Ligand
        from biopipelines.sequence import Sequence
        pipeline = new_pipeline("bioemu2_outputs")
        with pipeline:
            seqs = Sequence(seq=[PROTEIN] * len(sequences), type="protein", ids=list(sequences))
            options = {}
            if ligands:
                lig = Ligand(smiles=["CCO", "CCN"][:len(ligands)], ids=list(ligands))
                options["ligands"] = Bundle(lig) if bundle else lig
            BioEmu2(sequences=seqs, num_samples=num_samples, **options)
            return pipeline.tools[-1].get_output_files()

    def test_each_run_declares_num_samples_structures_under_the_framework_ids(self, local_config, isolated_cwd,
                                                                             new_pipeline):
        out = self._tool(new_pipeline)
        assert out["structures"].ids == ["p1+l1_<1..4>", "p2+l1_<1..4>"]
        assert out["trajectories"].ids == ["p1+l1", "p2+l1"]

    def test_an_apo_run_is_named_by_its_sequence(self, local_config, isolated_cwd, new_pipeline):
        out = self._tool(new_pipeline, ligands=())
        assert out["structures"].ids == ["p1_<1..4>", "p2_<1..4>"]

    def test_a_bundled_ligand_axis_is_one_system_per_sequence(self, local_config, isolated_cwd, new_pipeline):
        out = self._tool(new_pipeline, ligands=("l1", "l2"), bundle=True)
        assert out["trajectories"].ids == ["l1+l2+p1", "l1+l2+p2"]

    def test_the_stream_formats_and_files(self, local_config, isolated_cwd, new_pipeline):
        out = self._tool(new_pipeline)
        assert out["structures"].format == "pdb" and out["trajectories"].format == "xtc"
        assert out["structures"].files[0].endswith(os.path.join("structures", "<id>.pdb"))
        assert out["trajectories"].files[0].endswith(os.path.join("trajectories", "<id>.xtc"))
        assert out["structures"].map_table.endswith(os.path.join("structures", "structures_map.csv"))
        assert out["trajectories"].map_table.endswith(os.path.join("trajectories", "trajectories_map.csv"))

    def test_the_tables_carry_the_columns_the_pipe_script_writes(self, local_config, isolated_cwd, new_pipeline):
        tables = self._tool(new_pipeline)["tables"]
        assert tables["summary"].info.columns == pipe.SUMMARY_COLUMNS == [
            "id", "sequences.id", "ligands.id", "n_samples", "n_residues"]
        assert tables["physicality"].info.columns == pipe.PHYSICALITY_COLUMNS == [
            "id", "run.id", "frame", "status", "reason"]
        assert tables["missing"].info.columns == pipe.MISSING_COLUMNS == ["id", "removed_by", "kind", "cause"]

    def test_the_missing_table_is_the_standard_one(self, local_config, isolated_cwd, new_pipeline):
        tables = self._tool(new_pipeline)["tables"]
        assert tables["missing"].info.path.endswith(os.path.join("tables", "missing.csv"))

    def test_the_output_folder_is_returned(self, local_config, isolated_cwd, new_pipeline):
        assert "output_folder" in self._tool(new_pipeline)


def _axis(records, mode="each"):
    return {"mode": mode, "iterated": records if mode == "each" else [],
            "static": [] if mode == "each" else records, "static_first": False}


def _protein(rid, sequence=PROTEIN):
    return {"id": rid, "sequence": sequence}


class TestTheRunsAndTheirIds:

    def test_the_runtime_ids_are_the_ones_the_wrapper_predicts(self, tmp_path):
        """The completion check compares them, so a drift here fails a run that worked."""
        pd.DataFrame({"id": ["p1", "p2"], "sequence": [PROTEIN, PROTEIN]}).to_csv(tmp_path / "s.csv", index=False)
        pd.DataFrame({"id": ["l1", "l2"], "smiles": ["CCO", "CCN"]}).to_csv(tmp_path / "l.csv", index=False)
        seqs = DataStream(name="sequences", ids=["p1", "p2"], map_table=str(tmp_path / "s.csv"), format="csv")
        ligs = DataStream(name="compounds", ids=["l1", "l2"], map_table=str(tmp_path / "l.csv"), format="csv")
        for ligands in (ligs, Bundle(ligs)):
            path = str(tmp_path / "c.json")
            predicted = generate_combinatorics_config(
                path, sequences=(seqs, "sequences", "protein"), ligands=(ligands, "compounds", "ligand")
            ).predicted_ids
            runs = pipe.build_runs(pipe.axis_elements(CombinatoricsConfig.load(path)))
            assert [run["id"] for run in runs] == predicted

    def test_each_axis_gives_one_run_per_pair(self):
        runs = pipe.build_runs({"sequences": _axis([_protein("p1"), _protein("p2")]),
                                "ligands": _axis([{"id": "l1", "smiles": "CCO"}])})
        assert [run["id"] for run in runs] == ["p1+l1", "p2+l1"]
        assert runs[0]["provenance"] == {"sequences": "p1", "ligands": "l1"}

    def test_a_bundled_axis_is_in_every_run_and_in_its_provenance(self):
        runs = pipe.build_runs({"sequences": _axis([_protein("p1")]),
                                "ligands": _axis([{"id": "l1", "smiles": "CCO"}, {"id": "l2", "smiles": "CCN"}],
                                                 mode="bundle")})
        (run,) = runs
        assert run["id"] == "l1+l2+p1" and run["provenance"]["ligands"] == "l1+l2"

    def test_an_iterating_axis_with_no_records_gives_no_run_instead_of_crashing(self):
        assert pipe.build_runs({"sequences": _axis([])}) == []


class TestTheSystemYaml:

    @staticmethod
    def _system(records, contacts=None, covalent=None):
        return pipe.build_system(records, contacts, covalent)[0]

    def test_a_protein_and_a_smiles_ligand_take_chains_a_and_b(self):
        system = self._system({"sequences": [_protein("p1")], "ligands": [{"id": "l1", "smiles": "CCO"}]})
        assert system == {"entities": [
            {"entity_type": "protein", "sequence": PROTEIN, "asym_ids": ["A"]},
            {"entity_type": "ligand", "smiles": "CCO", "asym_ids": ["B"]}]}

    def test_an_apo_system_has_no_ligand_entity(self):
        assert len(self._system({"sequences": [_protein("p1")]})["entities"]) == 1

    def test_several_chains_continue_the_lettering_in_axis_order(self):
        system = self._system({"ligands": [{"id": "l1", "smiles": "CCO"}, {"id": "l2", "ccd": "ATP"}],
                               "sequences": [_protein("p1"), _protein("p2")]})
        assert [(e["entity_type"], e["asym_ids"]) for e in system["entities"]] == [
            ("protein", ["A"]), ("protein", ["B"]), ("ligand", ["C"]), ("ligand", ["D"])]

    def test_a_ccd_ligand_is_passed_by_code_as_boltz2_and_protenix_pass_it(self):
        system = self._system({"sequences": [_protein("p1")],
                               "ligands": [{"id": "atp", "ccd": "ATP", "smiles": "C", "source": "rcsb"}]})
        assert system["entities"][1] == {"entity_type": "ligand", "ccd_code": "ATP", "asym_ids": ["B"]}

    def test_a_smiles_with_a_backslash_survives_the_yaml_round_trip(self, tmp_path):
        """The old writer hand-built `smiles: "..."`, where `\\C` is an invalid escape and upstream could not parse it."""
        smiles = r"C/C=C\C(=O)[O-]"
        system = self._system({"sequences": [_protein("p1")], "ligands": [{"id": "l1", "smiles": smiles}]})
        path = tmp_path / "system.yaml"
        pipe.write_system_yaml(str(path), system)
        assert yaml.safe_load(path.read_text())["entities"][1]["smiles"] == smiles

    def test_a_ligand_without_chemistry_fails_its_run_instead_of_running_apo(self):
        """A code-only Ligand has an empty smiles; the old runner dropped the entity and sampled the protein alone."""
        with pytest.raises(pipe.RunError, match="no chemistry"):
            self._system({"sequences": [_protein("p1")], "ligands": [{"id": "hetatm", "code": "LIG"}]})

    @pytest.mark.parametrize("sequence", [float("nan"), "", "   "])
    def test_an_empty_sequence_fails_its_run_by_id(self, sequence):
        with pytest.raises(pipe.RunError, match="p_empty"):
            self._system({"sequences": [_protein("p_empty", sequence)]})

    def test_contacts_are_forwarded_unchanged_and_not_aliased(self):
        contacts = [dict(CONTACT)]
        system = self._system({"sequences": [_protein("p1")]}, contacts=contacts)
        assert system["contacts"] == contacts and system["contacts"] is not contacts

    def test_more_chains_than_letters_fails_the_run(self):
        records = {"sequences": [_protein(f"p{i}") for i in range(27)]}
        with pytest.raises(pipe.RunError, match="at most 26"):
            self._system(records)


class TestCovalentBond:

    LIGAND = {"id": "atp", "ccd": "ATP"}

    def _bond(self, link=LINK, ligand=LIGAND, sequence="MKCAYIAK"):
        return pipe.build_system({"sequences": [_protein("p1", sequence)], "ligands": [ligand]}, None, link)[0]

    def test_it_becomes_an_upstream_bonds_entry_and_not_a_covalent_linkage_key(self):
        """`covalent_linkage` is not in upstream's system schema; `bonds` with atom1 and atom2 is."""
        system = self._bond()
        assert "covalent_linkage" not in system
        assert system["bonds"] == [{
            "atom1": {"name": "SG", "res_name": "CYS", "res_id": 3, "asym_id": "A"},
            "atom2": {"name": "C1", "res_name": "ATP", "res_id": 1, "asym_id": "B"}}]

    def test_the_residue_name_comes_from_the_sequence(self):
        system = self._bond(link={**LINK, "position": 2, "protein_atom": "NZ"})
        assert system["bonds"][0]["atom1"]["res_name"] == "LYS"

    def test_the_ccd_code_is_used_even_when_the_row_would_otherwise_be_sent_as_smiles(self):
        system = self._bond(ligand={"id": "atp", "ccd": "ATP", "smiles": "C", "source": "user"})
        assert system["entities"][1] == {"entity_type": "ligand", "ccd_code": "ATP", "asym_ids": ["B"]}

    def test_a_smiles_only_ligand_fails_the_run_with_the_reason(self):
        with pytest.raises(pipe.RunError, match="CCD code"):
            self._bond(ligand={"id": "eth", "smiles": "CCO"})

    def test_a_position_past_the_end_of_the_chain_fails_the_run(self):
        with pytest.raises(pipe.RunError, match="outside chain A"):
            self._bond(link={**LINK, "position": 99})

    def test_a_chain_that_is_not_a_protein_chain_fails_the_run(self):
        with pytest.raises(pipe.RunError, match="not a protein chain"):
            self._bond(link={**LINK, "chain": "B"})

    def test_an_unidentified_residue_is_named_unk(self):
        assert self._bond(sequence="MKXAYIAK")["bonds"][0]["atom1"]["res_name"] == "UNK"


class TestTheCommand:

    CONFIG = {"num_samples": 8, "batch_size": 4, "sampler": "edm_gd", "physicality_filter": "no_filter",
              "seed": None}

    def test_the_flags_upstream_spells(self):
        command = pipe.build_command("", "/s.yaml", "/out", self.CONFIG)
        assert command == ["bioemu2", "sample-equilibrium", "--system", "/s.yaml", "--out", "/out",
                           "--num-samples", "8", "--batch-size", "4", "--sampler", "edm_gd",
                           "--physicality-filter", "no_filter"]

    def test_the_seed_is_passed_only_when_set(self):
        assert "--seed" not in pipe.build_command("", "s", "o", self.CONFIG)
        command = pipe.build_command("", "s", "o", {**self.CONFIG, "seed": 0})
        assert command[command.index("--seed") + 1] == "0"

    def test_a_container_prefix_goes_in_front_of_the_binary_only(self):
        command = pipe.build_command("apptainer exec --nv x.sif ", "s", "o", self.CONFIG)
        assert command[:5] == ["apptainer", "exec", "--nv", "x.sif", "bioemu2"]


def _upstream_cache_name(sequence):
    """Upstream's key, restated from msa/colabfold.py so the test does not share the code under test."""
    request = {"sequences": [sequence], "use_pairing": False, "pairing_strategy": "greedy", "use_env": True}
    return hashlib.sha256(json.dumps(request, sort_keys=True).encode("utf-8")).hexdigest() + ".json"


class TestAlignments:

    A3M = f">101\n{PROTEIN}\n>hit_1\n{PROTEIN[:10]}ab{PROTEIN[10:]}\n>hit_2\n{'-' * 5}{PROTEIN[5:]}\n"

    def test_the_cache_entry_has_the_name_and_checksum_upstream_validates(self):
        """Upstream refuses an entry whose sha256 is not that of the JSON list of alignments, and does not search again."""
        name, entry = pipe.msa_cache_entry(PROTEIN, self.A3M)
        assert name == _upstream_cache_name(PROTEIN)
        assert entry["a3ms"] == [self.A3M]
        assert entry["sha256"] == hashlib.sha256(json.dumps([self.A3M]).encode("utf-8")).hexdigest()
        assert entry["sha256"] != hashlib.sha256(self.A3M.encode("utf-8")).hexdigest()
        assert entry["request"] == {"sequences": [PROTEIN], "use_pairing": False,
                                    "pairing_strategy": "greedy", "use_env": True}

    def test_the_query_is_the_first_record_without_insertions_or_gaps(self):
        assert pipe.a3m_query(self.A3M) == PROTEIN
        assert pipe.a3m_query(">q\nMK-tT\nAY\n>h\nGG\n") == "MKTAY"

    def test_a_boltz_csv_becomes_an_a3m_with_the_query_first(self, tmp_path):
        pd.DataFrame({"key": [0, 1], "sequence": [PROTEIN, PROTEIN[:5] + "-" * 3]}).to_csv(tmp_path / "m.csv",
                                                                                           index=False)
        a3m = pipe.load_alignment(str(tmp_path / "m.csv"))
        assert a3m.splitlines()[0] == ">101" and pipe.a3m_query(a3m) == PROTEIN

    def test_a_missing_or_unknown_alignment_file_fails_the_run(self, tmp_path):
        with pytest.raises(pipe.RunError, match="does not exist"):
            pipe.load_alignment(str(tmp_path / "gone.a3m"))
        (tmp_path / "m.sto").write_text("x")
        with pytest.raises(pipe.RunError, match="neither"):
            pipe.load_alignment(str(tmp_path / "m.sto"))

    @staticmethod
    def _msas_json(tmp_path, rows, columns):
        table = tmp_path / "msas_map.csv"
        pd.DataFrame(rows, columns=columns).to_csv(table, index=False)
        stream = tmp_path / "msas.json"
        stream.write_text(json.dumps({"map_table": str(table)}))
        return str(stream)

    def test_an_alignment_is_found_by_the_sequence_id_it_was_built_for(self, tmp_path):
        a3m = tmp_path / "x.a3m"
        stream = self._msas_json(tmp_path, [["m1", "p1", str(a3m)]], ["id", "sequences.id", "file"])
        lookup = pipe.load_msa_lookup(stream)
        assert lookup["p1"] == lookup["m1"] == os.path.abspath(a3m)

    def test_an_alignment_without_a_sequence_column_is_found_by_its_own_id(self, tmp_path):
        stream = self._msas_json(tmp_path, [["p1", str(tmp_path / "x.a3m")]], ["id", "file"])
        assert set(pipe.load_msa_lookup(stream)) == {"p1"}

    def test_a_supplied_alignment_is_written_to_the_step_cache(self, tmp_path):
        (tmp_path / "p1.a3m").write_text(self.A3M)
        system, sources = pipe.build_system({"sequences": [_protein("p1")]}, None, None)
        staged = pipe.stage_msa(sources, {"p1": str(tmp_path / "p1.a3m")}, str(tmp_path / "cache"))
        assert staged == ["p1"]
        entry = json.loads((tmp_path / "cache" / _upstream_cache_name(PROTEIN)).read_text())
        assert entry["a3ms"] == [self.A3M]

    def test_an_alignment_for_another_sequence_fails_the_run_instead_of_crashing_upstream(self, tmp_path):
        (tmp_path / "p1.a3m").write_text(">101\nAAAAAAAA\n")
        _system, sources = pipe.build_system({"sequences": [_protein("p1")]}, None, None)
        with pytest.raises(pipe.RunError, match="not its sequence"):
            pipe.stage_msa(sources, {"p1": str(tmp_path / "p1.a3m")}, str(tmp_path / "cache"))
        assert not (tmp_path / "cache").exists()

    def test_a_sequence_without_an_alignment_is_left_to_the_colabfold_search(self, tmp_path):
        _system, sources = pipe.build_system({"sequences": [_protein("p1")]}, None, None)
        assert pipe.stage_msa(sources, {"other": "x.a3m"}, str(tmp_path / "cache")) == []
        assert not (tmp_path / "cache").exists()

    def test_a_system_with_several_protein_chains_refuses_a_supplied_alignment(self, tmp_path):
        (tmp_path / "p1.a3m").write_text(self.A3M)
        _system, sources = pipe.build_system({"sequences": [_protein("p1"), _protein("p2")]}, None, None)
        with pytest.raises(pipe.RunError, match="one protein chain"):
            pipe.stage_msa(sources, {"p1": str(tmp_path / "p1.a3m")}, str(tmp_path / "cache"))


class TestThePhysicalityReport:

    def test_the_report_is_read_by_frame(self, tmp_path):
        path = tmp_path / "physicality.json"
        path.write_text(json.dumps([{"frame": 0, "status": "passed", "reason": None},
                                    {"frame": 1, "status": "rejected", "reason": "bad_bonds_backbone"}]))
        report = pipe.read_physicality(str(path))
        assert report[1] == {"frame": 1, "status": "rejected", "reason": "bad_bonds_backbone"}

    def test_an_absent_or_unreadable_report_is_empty_rather_than_an_error(self, tmp_path):
        assert pipe.read_physicality(str(tmp_path / "none.json")) == {}
        (tmp_path / "bad.json").write_text("{not json")
        assert pipe.read_physicality(str(tmp_path / "bad.json")) == {}
        (tmp_path / "old.json").write_text(json.dumps({"frames": [{"clash": True}]}))
        assert pipe.read_physicality(str(tmp_path / "old.json")) == {}


SINGLE_CHAIN_CIF = """data_topology
#
loop_
_chem_comp_bond.pdbx_ordinal
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
1   MET CA  C   ? ? ?
2   MET C   O   ? ? ?
3   MET CA  CB  ? ? ?
4   MET N   CA  ? ? ?
5   MET CB  CG  ? ? ?
6   MET SD  CE  ? ? ?
7   MET CG  SD  ? ? ?
#
loop_
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.label_entity_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
_atom_site.id
ATOM N N   . MET A 1  . 1  MET A N   1 2.7527437    -5.8161097   15.860241   1 1
ATOM C CA  . MET A 1  . 1  MET A CA  1 2.800047     -5.19038     14.56017    1 2
ATOM C C   . MET A 1  . 1  MET A C   1 2.9121616    -6.137147    13.36379    1 3
ATOM O O   . MET A 1  . 1  MET A O   1 3.7524195    -6.99292     13.227924   1 4
ATOM C CB  . MET A 1  . 1  MET A CB  1 3.9507031    -4.144405    14.456848   1 5
ATOM C CG  . MET A 1  . 1  MET A CG  1 3.7267213    -3.102971    13.402943   1 6
ATOM S SD  . MET A 1  . 1  MET A SD  1 2.4838057    -1.802271    13.745588   1 7
ATOM C CE  . MET A 1  . 1  MET A CE  1 3.107978     -0.5251236   12.631086   1 8
"""


class TestTheTrajectorySplit:
    """A single-chain apo system is what failed on the cluster, in every run: mdtraj gives a CIF atom a text serial and its PDB writer only takes integers when the system has one chain."""

    def test_every_atom_gets_an_integer_serial_from_its_index(self):
        class Atom:
            def __init__(self, index):
                self.index, self.serial = index, str(index + 1)

        class Topology:
            atoms = [Atom(0), Atom(1), Atom(2)]

        pipe.number_atoms(Topology)
        assert [atom.serial for atom in Topology.atoms] == [1, 2, 3]

    def test_a_single_chain_cif_is_written_as_one_pdb_per_frame(self, tmp_path):
        md = pytest.importorskip("mdtraj")
        cif = tmp_path / "topology.cif"
        cif.write_text(SINGLE_CHAIN_CIF)
        frames = md.load(str(cif))
        frames = md.join([frames, frames])
        frames.save_xtc(str(tmp_path / "samples.xtc"))

        written, n_residues = pipe.split_trajectory(str(cif), str(tmp_path / "samples.xtc"), "run",
                                                    str(tmp_path / "structures"))

        assert [name for name, _ in written] == ["run_1", "run_2"]
        assert n_residues == 1
        atoms = [line for line in open(written[0][1]).read().splitlines() if line.startswith("ATOM")]
        assert [int(line[6:11]) for line in atoms] == list(range(1, 9))


class TestTheRunner:
    """`main` end to end with the sampler and the trajectory split replaced: ids, tables, missing rows and exit status."""

    @staticmethod
    def _world(tmp_path, sequences=("ok", "bad"), ligands=None, msas=None):
        pd.DataFrame({"id": list(sequences), "sequence": [PROTEIN] * len(sequences)}).to_csv(
            tmp_path / "s.csv", index=False)
        axes = {"sequences": (DataStream(name="sequences", ids=list(sequences), map_table=str(tmp_path / "s.csv"),
                                         format="csv"), "sequences", "protein")}
        if ligands:
            pd.DataFrame({"id": list(ligands), "smiles": ["CCO"] * len(ligands)}).to_csv(tmp_path / "l.csv",
                                                                                         index=False)
            axes["ligands"] = (DataStream(name="compounds", ids=list(ligands), map_table=str(tmp_path / "l.csv"),
                                          format="csv"), "compounds", "ligand")
        generate_combinatorics_config(str(tmp_path / "c.json"), **axes)
        (tmp_path / "config.yaml").write_text(yaml.safe_dump({
            "num_samples": 3, "batch_size": 2, "seed": None, "sampler": "edm",
            "physicality_filter": "high_quality", "contacts": None, "covalent_linkage": None}))
        out = tmp_path / "out"
        argv = ["pipe_bioemu2.py", "--config-yaml", str(tmp_path / "config.yaml"),
                "--combinatorics-config", str(tmp_path / "c.json"), "--runs-dir", str(out / "runs"),
                "--structures-dir", str(out / "structures"),
                "--structures-map-csv", str(out / "structures" / "m.csv"),
                "--trajectories-dir", str(out / "trajectories"),
                "--trajectories-map-csv", str(out / "trajectories" / "m.csv"),
                "--summary-csv", str(out / "tables" / "summary.csv"),
                "--physicality-csv", str(out / "tables" / "physicality.csv"),
                "--local-missing-csv", str(out / "005_BioEmu2" / "_execution" / "local_missing.csv")]
        if msas:
            argv += ["--msas-json", str(msas), "--msa-cache-dir", str(out / "msa_cache")]
        return out, argv

    @staticmethod
    def _fake(monkeypatch, fail=(), seen=None):
        def sample(command, env):
            out_dir = command[command.index("--out") + 1]
            run_id = os.path.basename(out_dir)
            if seen is not None:
                seen.append({"id": run_id, "command": command, "env": env,
                             "stale": os.path.exists(os.path.join(out_dir, "stale.txt"))})
            if run_id in fail:
                raise pipe.RunError("bioemu2 sample-equilibrium exited with status 1; see the step log")
            for name in ("samples.xtc", "topology.cif"):
                open(os.path.join(out_dir, name), "w").write(name)
            with open(os.path.join(out_dir, "physicality.json"), "w") as handle:
                json.dump([{"frame": 0, "status": "passed", "reason": None},
                           {"frame": 1, "status": "rejected", "reason": "bad_bonds_backbone"}], handle)

        def split(topology, xtc, run_id, structures_dir):
            os.makedirs(structures_dir, exist_ok=True)
            frames = []
            for k in (1, 2, 3):
                path = os.path.join(structures_dir, f"{run_id}_{k}.pdb")
                open(path, "w").write("ATOM")
                frames.append((f"{run_id}_{k}", path))
            return frames, 12

        monkeypatch.setattr(pipe, "sample", sample)
        monkeypatch.setattr(pipe, "split_trajectory", split)

    @staticmethod
    def _main(monkeypatch, argv):
        monkeypatch.setattr(sys, "argv", argv)
        try:
            pipe.main()
        except SystemExit as stop:
            return stop.code
        return 0

    def test_a_failed_run_is_in_the_table_the_missing_file_and_absent_from_the_maps(self, tmp_path, monkeypatch):
        out, argv = self._world(tmp_path)
        self._fake(monkeypatch, fail=("bad",))
        assert self._main(monkeypatch, argv) == 0

        summary = pd.read_csv(out / "tables" / "summary.csv")
        assert list(summary["id"]) == ["ok", "bad"]
        assert list(summary.columns) == pipe.SUMMARY_COLUMNS
        assert summary.loc[0, "n_samples"] == 3 and summary.loc[0, "n_residues"] == 12
        assert summary.loc[1, ["n_samples", "n_residues"]].isna().all()
        assert summary.loc[1, "sequences.id"] == "bad"

        assert list(pd.read_csv(out / "trajectories" / "m.csv")["id"]) == ["ok"]
        assert list(pd.read_csv(out / "structures" / "m.csv")["id"]) == ["ok_1", "ok_2", "ok_3"]

        missing = pd.read_csv(out / "005_BioEmu2" / "_execution" / "local_missing.csv")
        assert list(missing["id"]) == ["bad"] and list(missing["kind"]) == ["failure"]
        assert list(missing["removed_by"]) == ["005_BioEmu2"]
        assert "step log" in missing.loc[0, "cause"]
        assert list(missing.columns) == pipe.MISSING_COLUMNS

    def test_every_run_failing_exits_non_zero_and_still_writes_every_artifact(self, tmp_path, monkeypatch):
        out, argv = self._world(tmp_path)
        self._fake(monkeypatch, fail=("ok", "bad"))
        assert self._main(monkeypatch, argv) == 1
        assert len(pd.read_csv(out / "tables" / "summary.csv")) == 2
        assert len(pd.read_csv(out / "005_BioEmu2" / "_execution" / "local_missing.csv")) == 2
        assert pd.read_csv(out / "structures" / "m.csv").empty

    def test_no_run_at_all_exits_non_zero(self, tmp_path, monkeypatch):
        """An empty upstream stream must not be rubber-stamped as a step that worked."""
        out, argv = self._world(tmp_path, sequences=())
        self._fake(monkeypatch)
        assert self._main(monkeypatch, argv) == 1

    def test_the_physicality_table_has_a_row_for_every_structure_with_nan_where_unreported(
            self, tmp_path, monkeypatch):
        out, argv = self._world(tmp_path, sequences=("ok",))
        self._fake(monkeypatch)
        assert self._main(monkeypatch, argv) == 0
        table = pd.read_csv(out / "tables" / "physicality.csv")
        assert list(table.columns) == pipe.PHYSICALITY_COLUMNS
        assert list(table["id"]) == ["ok_1", "ok_2", "ok_3"] and list(table["frame"]) == [0, 1, 2]
        assert list(table["status"][:2]) == ["passed", "rejected"] and pd.isna(table.loc[2, "status"])
        assert table.loc[1, "reason"] == "bad_bonds_backbone" and pd.isna(table.loc[0, "reason"])
        assert set(table["run.id"]) == {"ok"}

    def test_provenance_columns_name_the_axes_that_made_each_run(self, tmp_path, monkeypatch):
        out, argv = self._world(tmp_path, sequences=("p1",), ligands=("l1", "l2"))
        self._fake(monkeypatch)
        assert self._main(monkeypatch, argv) == 0
        summary = pd.read_csv(out / "tables" / "summary.csv")
        assert list(summary["id"]) == ["p1+l1", "p1+l2"]
        assert list(summary["ligands.id"]) == ["l1", "l2"] and set(summary["sequences.id"]) == {"p1"}
        structures = pd.read_csv(out / "structures" / "m.csv")
        assert set(structures["run.id"]) == {"p1+l1", "p1+l2"}
        assert {"sequences.id", "ligands.id"} <= set(structures.columns)

    def test_a_stale_run_directory_cannot_pass_for_new_output(self, tmp_path, monkeypatch):
        out, argv = self._world(tmp_path, sequences=("ok",))
        stale = out / "runs" / "ok"
        stale.mkdir(parents=True)
        (stale / "stale.txt").write_text("from an earlier submission")
        seen = []
        self._fake(monkeypatch, seen=seen)
        assert self._main(monkeypatch, argv) == 0
        assert seen[0]["stale"] is False

    def test_a_container_prefix_prefixes_the_binary(self, tmp_path, monkeypatch):
        out, argv = self._world(tmp_path, sequences=("ok",))
        seen = []
        self._fake(monkeypatch, seen=seen)
        assert self._main(monkeypatch, argv + ["--container-prefix", "apptainer exec x.sif "]) == 0
        assert seen[0]["command"][:3] == ["apptainer", "exec", "x.sif"]

    def test_an_unexpected_error_is_recorded_with_its_type(self, tmp_path, monkeypatch):
        out, argv = self._world(tmp_path, sequences=("ok", "bad"))
        self._fake(monkeypatch)
        real = pipe.build_system

        def explode(records, contacts, covalent):
            if records["sequences"][0]["id"] == "bad":
                raise KeyError("boom")
            return real(records, contacts, covalent)
        monkeypatch.setattr(pipe, "build_system", explode)
        assert self._main(monkeypatch, argv) == 0
        missing = pd.read_csv(out / "005_BioEmu2" / "_execution" / "local_missing.csv")
        assert missing.loc[0, "cause"].startswith("KeyError")

    def test_a_supplied_alignment_redirects_only_that_runs_msa_cache(self, tmp_path, monkeypatch):
        (tmp_path / "ok.a3m").write_text(f">101\n{PROTEIN}\n")
        table = tmp_path / "msas_map.csv"
        pd.DataFrame([["ok", "ok", str(tmp_path / "ok.a3m")]], columns=["id", "sequences.id", "file"]).to_csv(
            table, index=False)
        (tmp_path / "msas.json").write_text(json.dumps({"map_table": str(table)}))
        out, argv = self._world(tmp_path, sequences=("ok", "bad"), msas=tmp_path / "msas.json")
        seen = []
        self._fake(monkeypatch, seen=seen)
        monkeypatch.delenv("BIOEMU2_MSA_CACHE", raising=False)
        assert self._main(monkeypatch, argv) == 0
        by_id = {entry["id"]: entry for entry in seen}
        assert by_id["ok"]["env"]["BIOEMU2_MSA_CACHE"] == str(out / "msa_cache")
        assert "BIOEMU2_MSA_CACHE" not in by_id["bad"]["env"]
        assert (out / "msa_cache" / _upstream_cache_name(PROTEIN)).is_file()


class TestTheInstall:

    ENV_YAML = os.path.join(REPO_ROOT, "environments", "bioemu2.yaml")
    PIP_TXT = os.path.join(REPO_ROOT, "environments", "bioemu2.pip.txt")

    @staticmethod
    def _script(**kwargs):
        return BioEmu2._install_script({"biopipelines": REPO_ROOT}, env_manager="mamba", **kwargs)

    def test_the_env_spec_and_the_install_script_name_the_same_python(self):
        """The first cluster install built Python 3.11 from a stale spec, and pip then found no bioemu2 for it."""
        with open(self.ENV_YAML) as handle:
            spec = yaml.safe_load(handle)
        pins = [d for d in spec["dependencies"] if isinstance(d, str) and d.split("=")[0] == "python"]
        assert pins == [f"python={MIN_PYTHON[0]}.{MIN_PYTHON[1]}"]
        script = self._script()
        assert f"sys.version_info >= {MIN_PYTHON}" in script
        assert f"{MIN_PYTHON[0]}.{MIN_PYTHON[1]}" in script
        assert "environments/bioemu2.yaml" in script

    def test_the_python_floor_is_upstreams_requires_python(self):
        """bioemu2 0.1.0rc0 declares requires-python >=3.12 on PyPI; a lower pin finds no distribution."""
        assert MIN_PYTHON >= (3, 12)

    def test_the_requirement_names_a_pre_release_so_pip_picks_the_only_published_version(self):
        """Without a pre-release specifier or --pre, pip skips a package whose only release is a release candidate."""
        from packaging.requirements import Requirement
        lines = [l.strip() for l in open(self.PIP_TXT) if l.strip() and not l.lstrip().startswith("#")]
        (requirement,) = [Requirement(l) for l in lines if l.lower().startswith("bioemu2")]
        assert requirement.specifier.contains("0.1.0rc0")
        assert requirement.specifier.contains("0.1.0")
        assert "cuda" in requirement.extras
        assert not any(l.startswith("--pre") for l in lines), "a global --pre would admit every dependency's pre-releases"

    def test_the_pipe_scripts_own_requirements_are_installed(self):
        lines = [l.strip() for l in open(self.PIP_TXT) if l.strip() and not l.lstrip().startswith("#")]
        assert "pandas" in lines

    def test_an_env_left_by_an_earlier_spec_is_rebuilt_instead_of_failing_on_create(self):
        script = self._script()
        guard = "if mamba run -n bioemu2 python -c \"\" >/dev/null 2>&1 && ! mamba run -n bioemu2 python -c"
        assert guard in script
        assert script.index(guard) < script.index("mamba env remove -n bioemu2 -y") < script.index("mamba env create")

    def test_an_env_with_a_good_python_gets_its_requirements_installed_into_it(self):
        script = self._script()
        assert f'pip install -r "{REPO_ROOT}/environments/bioemu2.pip.txt"' in script

    def test_the_checkpoint_is_fetched_and_a_failed_fetch_warns_instead_of_failing_the_install(self):
        script = self._script()
        assert "pipe_bioemu2_weights.py" in script
        assert "WARNING: the BioEmu2 checkpoint could not be fetched" in script

    def test_the_verification_checks_python_imports_and_the_cli_not_just_one_import(self):
        verify = self._script().split("# Verify installation")[1]
        for fragment in (f"sys.version_info >= {MIN_PYTHON}", "import bioemu2, mdtraj, torch",
                         "bioemu2 sample-equilibrium --help", 'touch "$INSTALL_SUCCESS"', "exit 1"):
            assert fragment in verify, fragment

    def test_a_forced_reinstall_removes_the_env_and_skips_the_already_installed_exit(self):
        forced = self._script(force_reinstall=True)
        assert "already installed, skipping" not in forced
        assert "Removing existing bioemu2 env for reinstall" in forced
        assert "already installed, skipping" in self._script()

    @pytest.mark.skipif(shutil.which("bash") is None, reason="needs bash")
    @pytest.mark.parametrize("force", [False, True])
    def test_the_script_is_valid_bash(self, force):
        done = subprocess.run(["bash", "-n"], input=self._script(force_reinstall=force), text=True,
                              capture_output=True)
        assert done.returncode == 0, done.stderr


class TestTheCheckpointFetch:
    """Pinned from checkpoints.py: two files per named checkpoint, from microsoft/bioemu2, through hf_hub_download."""

    @staticmethod
    def _fake_hub(monkeypatch, fail=False):
        calls = []
        module = type(sys)("huggingface_hub")

        def hf_hub_download(**kwargs):
            calls.append(kwargs)
            if fail:
                raise OSError("offline")
            return "/cache/" + kwargs["filename"].replace("/", "_")
        module.hf_hub_download = hf_hub_download
        monkeypatch.setitem(sys.modules, "huggingface_hub", module)
        return calls

    def test_both_files_of_the_named_checkpoint_are_requested_from_the_upstream_repository(self, monkeypatch):
        calls = self._fake_hub(monkeypatch)
        weights.fetch("bioemu2-preview", None)
        assert [c["filename"] for c in calls] == ["checkpoints/bioemu2-preview/checkpoint.ckpt",
                                                  "checkpoints/bioemu2-preview/config.yaml"]
        assert {c["repo_id"] for c in calls} == {"microsoft/bioemu2"}
        assert all(c["local_files_only"] is False for c in calls)

    def test_check_only_never_reaches_the_network(self, monkeypatch):
        calls = self._fake_hub(monkeypatch)
        monkeypatch.setattr(sys, "argv", ["pipe_bioemu2_weights.py", "--check-only", "--cache-dir", "/c"])
        weights.main()
        assert all(c["local_files_only"] is True and c["cache_dir"] == "/c" for c in calls)

    def test_a_failed_download_exits_non_zero(self, monkeypatch):
        self._fake_hub(monkeypatch, fail=True)
        monkeypatch.setattr(sys, "argv", ["pipe_bioemu2_weights.py"])
        with pytest.raises(SystemExit) as stop:
            weights.main()
        assert stop.value.code == 1

    def test_the_default_checkpoint_is_the_one_upstream_loads(self):
        assert weights.DEFAULT_CHECKPOINT == "bioemu2-preview"
