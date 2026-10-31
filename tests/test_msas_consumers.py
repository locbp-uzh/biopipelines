# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Every consumer of an `msas` stream finds its alignments through `id_map_utils.MsaLookup`, and a protein without one fails alone.

The predictors used to carry one copy of the lookup each (`load_msa_lookup`, `load_msa_mappings`, `msa_paths_for`), and three of them matched on `id` alone, so an alignment from MMseqs2 matched no chain and the fold ran single-sequence in silence. These tests pin the shared function's use, the per-id failure (a blank row of the result table and a `failure` row of `missing`, the rest of the batch running), and the exact earlier bug.

Not covered, because it needs a cluster: the predictors accepting the alignments, and a real batch.
"""

import ast
import json
import os
import subprocess
import sys
import types

import pandas as pd
import pytest

from biopipelines.biopipelines_io import read_skipped_queries, write_skipped_queries
from biopipelines.id_map_utils import MsaLookup
from tests._helpers import load_pipe_script

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(REPO_ROOT, "pipe_scripts")

PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"


def _module(name):
    return load_pipe_script(name, f"{name}.py")


def _write_table(tmp_path, rows, name="msas_map.csv"):
    path = tmp_path / name
    pd.DataFrame(rows).to_csv(path, index=False)
    return str(path)


def _a3m(tmp_path, name="p1"):
    path = tmp_path / f"{name}.a3m"
    path.write_text(f">{name}\n{PROTEIN}\n")
    return str(path)


def _axis(entity_type, records):
    return {"entity_type": entity_type, "mode": "each", "iterated": records, "static": [],
            "static_first": False}


TWO_PROTEINS = {"proteins": _axis("protein", [{"id": "p1", "sequence": PROTEIN},
                                              {"id": "p2", "sequence": "MKTAYIAKQRQ"}])}


def _boltz_args():
    return types.SimpleNamespace(
        affinity=False, single_sequence=False, template=None, template_chains=None, glycosylation=None,
        covalent_linkage=None, contacts=None, disulfide_bonds=None, metal_coord=None,
        pocket_residues=None, pocket_binder=None, pocket_max_distance=None)


class TestTheSkippedQueryFile:

    def test_what_a_config_step_writes_the_postprocess_reads(self, tmp_path):
        path = str(tmp_path / "x" / "skipped.csv")
        write_skipped_queries(path, {"a": "no alignment in msas for a (ids present: b)", "c": "boom"})
        assert read_skipped_queries(path) == {"a": "no alignment in msas for a (ids present: b)", "c": "boom"}

    def test_no_file_is_no_skips(self, tmp_path):
        assert read_skipped_queries(None) == {} and read_skipped_queries(str(tmp_path / "none.csv")) == {}


class TestAlphaFoldWritesTheSequenceId:

    def test_its_msas_table_has_no_repeated_column(self, tmp_path):
        folder = tmp_path / "msas"
        folder.mkdir()
        (folder / "p1.a3m").write_text(f">101\n{PROTEIN}\n")
        pd.DataFrame([{"id": "p1", "sequence": PROTEIN}]).to_csv(tmp_path / "queries.csv", index=False)
        done = subprocess.run([sys.executable, os.path.join(SCRIPTS, "pipe_alphafold_msas.py"), str(folder),
                               str(tmp_path / "queries.csv"), str(tmp_path / "msas_map.csv")],
                              capture_output=True, text=True)
        assert done.returncode == 0, done.stdout + done.stderr
        table = pd.read_csv(tmp_path / "msas_map.csv")
        assert table.columns.tolist() == ["id", "sequence", "file"] and table["id"].tolist() == ["p1"]


class TestAMissingAlignmentFailsOnlyItsOwnId:

    def test_openfold3_skips_the_query_and_the_postprocess_leaves_a_blank_row_and_a_failure(self, tmp_path):
        config = _module("pipe_openfold3_config")
        msas = MsaLookup([{"id": "p1", "file": _a3m(tmp_path)}])
        alignments = config.Alignments(msas, str(tmp_path / "staged"))
        payload, skipped = config.build_queries(TWO_PROTEINS, alignments, None, None)
        assert list(payload["queries"]) == ["p1"]
        assert skipped == {"p2": "no alignment in msas for p2 (ids present: p1)"}
        assert payload["queries"]["p1"]["chains"][0]["main_msa_file_paths"][0].endswith("colabfold_main.a3m")

        write_skipped_queries(str(tmp_path / "skipped.csv"), skipped)
        pred = tmp_path / "pred" / "p1" / "seed_1"
        pred.mkdir(parents=True)
        (pred / "p1_seed_1_sample_1_model.cif").write_text("data_p1\n")
        (pred / "p1_seed_1_sample_1_confidences_aggregated.json").write_text(
            json.dumps({"sample_ranking_score": 0.5}))
        (tmp_path / "queries.json").write_text(json.dumps(payload))
        out = tmp_path / "out"
        done = subprocess.run(
            [sys.executable, os.path.join(SCRIPTS, "pipe_openfold3_postprocess.py"),
             "--prediction-folder", str(tmp_path / "pred"), "--combinatorics-config", str(tmp_path / "none.json"),
             "--queries-json", str(tmp_path / "queries.json"), "--skipped-csv", str(tmp_path / "skipped.csv"),
             "--structures-dir", str(out / "structures"), "--structures-map-csv", str(out / "structures" / "m.csv"),
             "--confidence-csv", str(out / "c.csv"), "--local-missing-csv", str(out / "x.csv")],
            capture_output=True, text=True)
        assert done.returncode == 0, done.stdout + done.stderr
        confidence = pd.read_csv(out / "c.csv")
        assert confidence["id"].tolist() == ["p1", "p2"]
        assert confidence.loc[confidence["id"] == "p2", "plddt"].isna().all()
        missing = pd.read_csv(out / "x.csv")
        assert missing["id"].tolist() == ["p2"] and missing["kind"].tolist() == ["failure"]
        assert missing["cause"][0].startswith("no alignment in msas for p2")

    def test_openfold3_stages_only_the_alignments_that_are_asked_for(self, tmp_path):
        config = _module("pipe_openfold3_config")
        msas = MsaLookup([{"id": "p1", "file": _a3m(tmp_path)}, {"id": "other", "file": _a3m(tmp_path, "other")}])
        alignments = config.Alignments(msas, str(tmp_path / "staged"))
        config.build_queries({"proteins": _axis("protein", [{"id": "p1", "sequence": PROTEIN}])},
                             alignments, None, None)
        assert os.listdir(tmp_path / "staged") == ["p1"]

    def test_rf3_skips_the_example_and_the_postprocess_leaves_a_blank_row_and_a_failure(self, tmp_path):
        config = _module("pipe_rf3_config")
        msas = MsaLookup([{"id": "p1_msa", "sequences.id": "p1", "file": _a3m(tmp_path)}])
        examples, skipped = config.build_examples(TWO_PROTEINS, msas)
        assert [e["name"] for e in examples] == ["p1"] and examples[0]["components"][0]["msa_path"].endswith("p1.a3m")
        assert skipped == {"p2": "no alignment in msas for p2 (ids present: p1_msa)"}

        write_skipped_queries(str(tmp_path / "skipped.csv"), skipped)
        (tmp_path / "queries.json").write_text(json.dumps(examples))
        out = tmp_path / "out"
        done = subprocess.run(
            [sys.executable, os.path.join(SCRIPTS, "pipe_rf3_postprocess.py"),
             "--prediction-folder", str(tmp_path / "pred"), "--queries-json", str(tmp_path / "queries.json"),
             "--skipped-csv", str(tmp_path / "skipped.csv"), "--structures-dir", str(out / "structures"),
             "--structures-map-csv", str(out / "structures" / "m.csv"), "--confidence-csv", str(out / "c.csv"),
             "--pae-dir", str(out / "pae"), "--pae-map-csv", str(out / "pae" / "pae_map.csv"),
             "--local-missing-csv", str(out / "001_RF3" / "execution" / "local_missing.csv"),
             "--diffusion-batch-size", "2"],
            capture_output=True, text=True)
        # p1 produced no structure here, so the step reports total failure; the rows are what is under test.
        assert done.returncode == 1 and "no RF3 prediction produced a structure" in done.stdout
        missing = pd.read_csv(out / "001_RF3" / "execution" / "local_missing.csv")
        assert sorted(missing["id"]) == ["p1", "p2"]
        assert missing.loc[missing["id"] == "p2", "cause"].iloc[0].startswith("no alignment in msas for p2")
        assert sorted(pd.read_csv(out / "c.csv")["id"]) == ["p1", "p2"]

    def test_rf3_without_msas_builds_every_example_as_before(self):
        config = _module("pipe_rf3_config")
        examples, skipped = config.build_examples(TWO_PROTEINS, None)
        assert [e["name"] for e in examples] == ["p1", "p2"] and not skipped
        assert all("msa_path" not in c for e in examples for c in e["components"])

    def test_protenix_names_the_protein_and_the_ids_that_are_present(self, tmp_path):
        config = _module("pipe_protenix_config")
        msas = MsaLookup([{"id": "p1", "file": _a3m(tmp_path)}])
        jobs, skipped = config.build_queries(TWO_PROTEINS, msas, None)
        assert list(jobs) == ["p1"]
        assert skipped == {"p2": "no alignment in msas for p2 (ids present: p1)"}

    def test_boltz2_skips_the_config_instead_of_folding_the_chain_single_sequence(self, tmp_path):
        config = _module("pipe_boltz_config_unified")
        msas = MsaLookup([{"id": "p1", "file": _a3m(tmp_path)}])
        configs, skipped = config.generate_configs(TWO_PROTEINS, msas, _boltz_args())
        assert [cid for cid, _ in configs] == ["p1"]
        assert configs[0][1]["sequences"][0]["protein"]["msa"].endswith("p1.a3m")
        assert skipped == {"p2": "no alignment in msas for p2 (ids present: p1)"}

    def test_boltz2_without_msas_leaves_the_server_to_search(self):
        config = _module("pipe_boltz_config_unified")
        configs, skipped = config.generate_configs(TWO_PROTEINS, None, _boltz_args())
        assert len(configs) == 2 and not skipped
        assert all("msa" not in entry["protein"] for _, c in configs for entry in c["sequences"])

    def test_boltz2_postprocess_leaves_a_blank_row_and_a_failure_for_a_skipped_config(self, tmp_path):
        prediction = tmp_path / "pred" / "boltz_results_x" / "predictions" / "p1"
        prediction.mkdir(parents=True)
        (prediction / "p1_model_0.pdb").write_text("ATOM\n")
        (prediction / "confidence_p1_model_0.json").write_text(json.dumps({"confidence_score": 0.9, "ptm": 0.8}))
        pd.DataFrame({"id": ["p1"]}).to_csv(tmp_path / "ids.csv", index=False)
        write_skipped_queries(str(tmp_path / "skipped.csv"), {"p2": "no alignment in msas for p2 (ids present: p1)"})
        out = tmp_path / "out"
        out.mkdir()
        (out / "structures").mkdir()
        missing = out / "001_Boltz2" / "execution" / "local_missing.csv"
        done = subprocess.run(
            [sys.executable, os.path.join(SCRIPTS, "pipe_boltz_postprocessing.py"), str(tmp_path / "pred"),
             str(out), str(tmp_path / "ids.csv"), "--structures-folder", str(out / "structures"),
             "--confidence-csv", str(out / "confidence.csv"), "--msas-folder", str(out / "msas"),
             "--msas-csv", str(out / "msas_map.csv"), "--sequences-csv", str(out / "sequences.csv"),
             "--pae-folder", str(out / "pae"), "--pae-map", str(out / "pae_map.csv"),
             "--skipped-csv", str(tmp_path / "skipped.csv"), "--local-missing-csv", str(missing)],
            capture_output=True, text=True)
        assert done.returncode == 0, done.stdout + done.stderr
        confidence = pd.read_csv(out / "confidence.csv")
        assert confidence["id"].tolist() == ["p1", "p2"]
        assert confidence.loc[confidence["id"] == "p2", "confidence_score"].isna().all()
        rows = pd.read_csv(missing)
        assert rows["id"].tolist() == ["p2"] and rows["kind"].tolist() == ["failure"]
        assert rows["removed_by"].tolist() == ["001_Boltz2"]

    def test_alphafold_copies_each_query_its_alignment_and_records_the_ones_without(self, tmp_path):
        copy = _module("pipe_boltz_msa_copy")
        table = _write_table(tmp_path, [{"id": "p1_msa", "sequences.id": "p1", "file": _a3m(tmp_path)}])
        pd.DataFrame(TWO_PROTEINS["proteins"]["iterated"]).to_csv(tmp_path / "queries.csv", index=False)
        out = tmp_path / "staged"
        missing = tmp_path / "001_AlphaFold" / "execution" / "local_missing.csv"
        copied = copy.copy_msa_files(table, str(out), str(tmp_path / "queries.csv"), str(missing))
        assert copied == 1 and os.listdir(out) == ["p1.a3m"]
        rows = pd.read_csv(missing)
        assert rows["id"].tolist() == ["p2"] and rows["kind"].tolist() == ["failure"]
        assert rows["removed_by"].tolist() == ["001_AlphaFold"]
        assert rows["cause"][0] == "no alignment in msas for p2 (ids present: p1_msa)"

    def test_boltz2_recycling_names_each_copy_after_the_sequence_it_belongs_to(self, tmp_path):
        copy = _module("pipe_boltz_msa_copy")
        table = _write_table(tmp_path, [{"id": "p1_msa", "sequences.id": "p1", "file": _a3m(tmp_path)},
                                        {"id": "p2", "file": _a3m(tmp_path, "p2")}])
        out = tmp_path / "recycled"
        assert copy.copy_msa_files(table, str(out)) == 2
        assert sorted(os.listdir(out)) == ["p1.a3m", "p2.a3m"]

    def test_alphafold_keeps_a_blank_confidence_row_for_the_id_without_an_alignment(self, tmp_path):
        folding = tmp_path / "folding"
        folding.mkdir()
        missing = tmp_path / "missing.csv"
        pd.DataFrame([{"id": "p2", "removed_by": "001_AlphaFold", "kind": "failure", "cause": "x"}]).to_csv(
            missing, index=False)
        done = subprocess.run(
            [sys.executable, os.path.join(SCRIPTS, "pipe_alphafold_confidence.py"), str(folding),
             str(tmp_path / "structures"), str(tmp_path / "confidence.csv"), "--missing-csv", str(missing)],
            capture_output=True, text=True)
        assert done.returncode == 0, done.stdout + done.stderr
        confidence = pd.read_csv(tmp_path / "confidence.csv")
        assert confidence["id"].tolist() == ["p2"] and confidence["plddt"].isna().all()

    def test_esmfold2_fails_the_complex_not_the_batch(self, tmp_path):
        inference = _module("pipe_esmfold2_inference")
        table = _write_table(tmp_path, TWO_PROTEINS["proteins"]["iterated"], "proteins.csv")
        config = {"axes": {"proteins": {"entity_type": "protein", "mode": "each",
                                        "sources": [{"path": table, "iterate": True}]}}}
        msas = MsaLookup([{"id": "p1", "file": _a3m(tmp_path)}])
        built = dict(inference.build_complexes(config, msas))
        assert built["p1"][0]["msa_path"].endswith("p1.a3m") and built["p1"][0]["msa_error"] is None
        assert built["p2"][0]["msa_error"] == "no alignment in msas for p2 (ids present: p1)"
        with pytest.raises(ValueError, match="no alignment in msas for p2"):
            inference.to_structure_input(built["p2"], None, None, 1024, {})

    def test_esmfold2_still_matches_by_sequence_when_the_ids_differ(self, tmp_path):
        inference = _module("pipe_esmfold2_inference")
        table = _write_table(tmp_path, [{"id": "p1", "sequence": PROTEIN}], "proteins.csv")
        config = {"axes": {"proteins": {"entity_type": "protein", "mode": "each",
                                        "sources": [{"path": table, "iterate": True}]}}}
        msas = MsaLookup([{"id": "renamed", "sequence": PROTEIN, "file": _a3m(tmp_path)}])
        (_cid, chains), = list(inference.build_complexes(config, msas))
        assert chains[0]["msa_path"].endswith("p1.a3m") and chains[0]["msa_error"] is None


class TestTheEarlierBugWithTheNamesItWasReportedFor:
    """A legacy-named map (`<id>_msa` rows, the id in `sequences.id`) and the protein id, for the three predictors that matched on `id` alone."""

    ONE = {"proteins": _axis("protein", [{"id": "p1", "sequence": PROTEIN}])}

    @pytest.fixture
    def legacy(self, tmp_path):
        return MsaLookup.load(_write_table(
            tmp_path, [{"id": "p1_msa", "sequences.id": "p1", "sequence": "", "file": _a3m(tmp_path)}]))

    def test_openfold3(self, tmp_path, legacy):
        config = _module("pipe_openfold3_config")
        payload, skipped = config.build_queries(self.ONE, config.Alignments(legacy, str(tmp_path / "s")), None, None)
        assert not skipped and "main_msa_file_paths" in payload["queries"]["p1"]["chains"][0]

    def test_rf3(self, legacy):
        config = _module("pipe_rf3_config")
        examples, skipped = config.build_examples(self.ONE, legacy)
        assert not skipped and examples[0]["components"][0]["msa_path"].endswith("p1.a3m")

    def test_protenix(self, legacy):
        config = _module("pipe_protenix_config")
        jobs, skipped = config.build_queries(self.ONE, legacy, None)
        assert not skipped and jobs["p1"]["sequences"][0]["proteinChain"]["unpairedMsaPath"].endswith("p1.a3m")


CONSUMERS = [
    "pipe_openfold3_config.py", "pipe_rf3_config.py", "pipe_protenix_config.py", "pipe_boltz_config_unified.py",
    "pipe_boltz_msa_copy.py", "pipe_esmfold2_inference.py", "pipe_bioemu2.py", "pipe_mutagenesis_msa.py",
]
OLD_PER_TOOL_FUNCTIONS = {"load_msa_lookup", "load_msa_mappings", "load_msa_by_sequence", "msa_paths_for",
                          "stage_alignments"}


def _tree(path):
    with open(path, encoding="utf-8") as handle:
        return ast.parse(handle.read())


class TestOneLookupForEveryConsumer:

    @pytest.mark.parametrize("script", CONSUMERS)
    def test_the_consumer_builds_its_lookup_through_the_shared_function(self, script):
        tree = _tree(os.path.join(SCRIPTS, script))
        imported = {alias.name for node in ast.walk(tree)
                    if isinstance(node, ast.ImportFrom) and node.module == "biopipelines.id_map_utils"
                    for alias in node.names}
        assert "MsaLookup" in imported, f"{script} does not import MsaLookup from biopipelines.id_map_utils"

    @pytest.mark.parametrize("script", CONSUMERS)
    def test_the_per_tool_copies_are_gone(self, script):
        tree = _tree(os.path.join(SCRIPTS, script))
        defined = {node.name for node in ast.walk(tree) if isinstance(node, ast.FunctionDef)}
        assert not defined & OLD_PER_TOOL_FUNCTIONS, f"{script} still defines {sorted(defined & OLD_PER_TOOL_FUNCTIONS)}"

    # pipe_bioemu2 writes `sequences.id` as a provenance column of its own outputs, which is not matching.
    @pytest.mark.parametrize("script", [name for name in CONSUMERS if name != "pipe_bioemu2.py"])
    def test_the_consumer_does_not_read_the_provenance_column_itself(self, script):
        source = open(os.path.join(SCRIPTS, script), encoding="utf-8").read()
        assert "sequences.id" not in source, f"{script} matches on sequences.id itself instead of asking MsaLookup"


# Each wrapper, the pipe scripts it calls with the skipped-query plumbing, and the format of the msas it takes.
WIRING = {
    "AlphaFold": ("alphafold", ["pipe_boltz_msa_copy.py", "pipe_alphafold_confidence.py"], "a3m"),
    "Boltz2": ("boltz2", ["pipe_boltz_config_unified.py", "pipe_boltz_postprocessing.py"], "csv"),
    "OpenFold3": ("openfold3", ["pipe_openfold3_config.py", "pipe_openfold3_postprocess.py"], "a3m"),
    "RF3": ("rf3", ["pipe_rf3_config.py", "pipe_rf3_postprocess.py"], "a3m"),
}


class TestTheWrapperWiring:

    @staticmethod
    def _built(new_pipeline, name, with_msas):
        import importlib
        from pathlib import Path
        from biopipelines.mock import Mock
        module, _scripts, fmt = WIRING[name]
        tool_class = getattr(importlib.import_module(f"biopipelines.{module}"), name)
        pipeline = new_pipeline(f"msas_wiring_{name}_{with_msas}")
        with pipeline:
            seqs = Mock(ids=["p1"], streams={"sequences": {"format": "csv", "file": "<id>.csv"}})
            kwargs = {"proteins": seqs.streams.sequences}
            if with_msas:
                msas = Mock(ids=["p1"], streams={"msas": {"format": fmt, "file": f"<id>.{fmt}"}})
                kwargs["msas"] = msas if name == "AlphaFold" else msas.streams.msas
            tool = tool_class(**kwargs)
            script_path = pipeline.save()
        steps = sorted(Path(script_path).parent.glob(f"*_{name}.sh"))
        return tool, "\n".join(step.read_text(encoding="utf-8", errors="replace") for step in steps)

    @pytest.mark.parametrize("name", list(WIRING))
    def test_every_flag_the_wrapper_sends_is_known_to_the_pipe_script(self, local_config, isolated_cwd,
                                                                      new_pipeline, name):
        import re
        _tool, script = self._built(new_pipeline, name, with_msas=True)
        for pipe in WIRING[name][1]:
            block = script.split(pipe)[1].split("\n\n")[0]
            sent = set(re.findall(r"(?<![\w-])--[a-z][a-z-]*", block))
            source = open(os.path.join(SCRIPTS, pipe), encoding="utf-8").read()
            known = set(re.findall(r"""add_argument\(\s*['"](--[a-z-]+)['"]""", source))
            assert sent <= known, f"{pipe}: {sorted(sent - known)} are sent but not read"

    @pytest.mark.parametrize("name", list(WIRING))
    def test_supplying_msas_declares_the_missing_table_that_excuses_the_ids_left_out(
            self, local_config, isolated_cwd, new_pipeline, name):
        tool, script = self._built(new_pipeline, name, with_msas=True)
        assert "missing" in list(tool.tables.keys())
        assert "skipped_" in script or "--missing-csv" in script

    @pytest.mark.parametrize("name", ["AlphaFold", "Boltz2"])
    def test_without_msas_nothing_new_is_declared_or_sent(self, local_config, isolated_cwd, new_pipeline, name):
        tool, script = self._built(new_pipeline, name, with_msas=False)
        assert "missing" not in list(tool.tables.keys())
        assert "local_missing.csv" not in script and "--missing-csv" not in script
