# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Foldseek: parameter validation, the emitted bash, the output columns, and the runtime half (staging, command lines, cluster TSV parsing, the table) with the binary replaced by a fake that writes a cluster TSV."""

import glob
import importlib.util
import json
import os
import shutil
import sys

import pandas as pd
import pytest

pytest.importorskip("yaml")
import yaml  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "pipe_scripts", "pipe_foldseek.py")


@pytest.fixture(scope="module")
def pipe():
    spec = importlib.util.spec_from_file_location("pipe_foldseek", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CFG = {"tm_thresholds": [0.5, 0.8], "multimer": False, "alignment_type": 1, "min_seq_id": 0.0,
       "evalue": None, "coverage": None, "cov_mode": 0,
       "multimer_tm_threshold": None, "interface_lddt_threshold": None, "threads": 1}


# ── runtime half ──────────────────────────────────────────────────────────────

def test_threshold_tag_is_the_percent_form(pipe):
    assert [pipe.threshold_tag(t) for t in (0.5, 0.65, 0.8, 0.05)] == ["tm50", "tm65", "tm80", "tm05"]


def test_entry_names_resolve_to_the_staged_index(pipe):
    assert pipe.entry_index("12") == 12
    assert pipe.entry_index("12_A") == 12
    assert pipe.entry_index("3.pdb_B") == 3
    assert pipe.entry_index("junk") is None


def test_the_single_chain_command_carries_the_threshold_and_no_default_min_seq_id(pipe):
    cmd = pipe.build_command("in", "out/res", "out/tmp", 0.7, CFG)
    assert cmd[:5] == ["foldseek", "easy-cluster", "in", "out/res", "out/tmp"]
    assert cmd[cmd.index("--tmscore-threshold") + 1] == "0.7"
    assert cmd[cmd.index("--alignment-type") + 1] == "1"
    assert "--min-seq-id" not in cmd and "--multimer-tm-threshold" not in cmd


def test_min_seq_id_is_passed_when_set(pipe):
    cmd = pipe.build_command("in", "o/res", "o/tmp", 0.5, dict(CFG, min_seq_id=0.3))
    assert cmd[cmd.index("--min-seq-id") + 1] == "0.3"


def test_unset_evalue_and_coverage_emit_no_flag_and_cov_mode_zero_stays(pipe):
    for multimer in (False, True):
        cfg = dict(CFG, multimer=multimer, multimer_tm_threshold=0.65, interface_lddt_threshold=0.65)
        cmd = pipe.build_command("in", "o/res", "o/tmp", 0.5, cfg)
        assert "-e" not in cmd and "-c" not in cmd
        assert cmd[cmd.index("--cov-mode") + 1] == "0"


def test_evalue_and_coverage_are_passed_when_set(pipe):
    cmd = pipe.build_command("in", "o/res", "o/tmp", 0.5, dict(CFG, evalue=1e-5, coverage=0.9))
    assert cmd[cmd.index("-e") + 1] == "1e-05"
    assert cmd[cmd.index("-c") + 1] == "0.9"


def test_a_small_coverage_is_written_without_an_exponent(pipe):
    cmd = pipe.build_command("in", "o/res", "o/tmp", 0.5, dict(CFG, coverage=1e-5))
    assert cmd[cmd.index("-c") + 1] == "0.00001"


def test_each_of_evalue_and_coverage_is_passed_alone(pipe):
    only_e = pipe.build_command("in", "o/res", "o/tmp", 0.5, dict(CFG, evalue=0.01))
    assert "-e" in only_e and "-c" not in only_e
    only_c = pipe.build_command("in", "o/res", "o/tmp", 0.5, dict(CFG, coverage=0.0))
    assert only_c[only_c.index("-c") + 1] == "0.0" and "-e" not in only_c


def test_cov_mode_is_passed_as_set(pipe):
    for mode in range(6):
        cmd = pipe.build_command("in", "o/res", "o/tmp", 0.5, dict(CFG, cov_mode=mode))
        assert cmd[cmd.index("--cov-mode") + 1] == str(mode)


def test_the_new_flags_reach_the_multimer_command_beside_its_own_thresholds(pipe):
    cfg = dict(CFG, multimer=True, multimer_tm_threshold=0.65, interface_lddt_threshold=0.65,
               evalue=0.1, coverage=0.8, cov_mode=2)
    cmd = pipe.build_command("in", "o/res", "o/tmp", 0.6, cfg)
    assert cmd[1] == "easy-multimercluster"
    assert cmd[cmd.index("-e") + 1] == "0.1" and cmd[cmd.index("-c") + 1] == "0.8"
    assert cmd[cmd.index("--cov-mode") + 1] == "2"
    assert cmd[cmd.index("--chain-tm-threshold") + 1] == "0.6"


def test_a_config_written_before_the_new_keys_still_builds_the_old_command(pipe):
    old = {k: v for k, v in CFG.items() if k not in ("evalue", "coverage", "cov_mode")}
    cmd = pipe.build_command("in", "o/res", "o/tmp", 0.5, old)
    assert cmd == pipe.build_command("in", "o/res", "o/tmp", 0.5, CFG)
    assert cmd[cmd.index("--cov-mode") + 1] == "0" and "-e" not in cmd and "-c" not in cmd


def test_the_multimer_command_sweeps_the_chain_threshold_and_pins_the_rest(pipe):
    cfg = dict(CFG, multimer=True, multimer_tm_threshold=0.65, interface_lddt_threshold=0.6)
    cmd = pipe.build_command("in", "o/res", "o/tmp", 0.6, cfg)
    assert cmd[1] == "easy-multimercluster"
    assert cmd[cmd.index("--chain-tm-threshold") + 1] == "0.6"
    assert cmd[cmd.index("--multimer-tm-threshold") + 1] == "0.65"
    assert cmd[cmd.index("--interface-lddt-threshold") + 1] == "0.6"
    assert "--tmscore-threshold" not in cmd


def test_the_container_prefix_leads_the_command(pipe):
    cmd = pipe.build_command("in", "o/res", "o/tmp", 0.5, CFG, "apptainer exec -B /data img.sif ")
    assert cmd[:5] == ["apptainer", "exec", "-B", "/data", "img.sif"]
    assert cmd[5] == "foldseek"


def test_members_map_to_their_representatives(pipe):
    index_to_id = {0: "a", 1: "b", 2: "c"}
    pairs = [("0", "0"), ("0", "1"), ("2", "2")]
    assigned, problems = pipe.assign_representatives(pairs, index_to_id, multimer=False)
    assert assigned == {"a": "a", "b": "a", "c": "c"}
    assert problems == {}


def test_a_structure_foldseek_did_not_report_is_a_problem_not_a_guess(pipe):
    assigned, problems = pipe.assign_representatives([("0", "0")], {0: "a", 1: "b"}, multimer=False)
    assert assigned == {"a": "a"}
    assert "not reported" in problems["b"]


def test_a_multi_chain_file_is_refused_in_single_chain_mode_only(pipe):
    pairs = [("0_A", "0_A"), ("0_A", "0_B")]
    assigned, problems = pipe.assign_representatives(pairs, {0: "a"}, multimer=False)
    assert assigned == {} and "multimer=True" in problems["a"]
    assigned, problems = pipe.assign_representatives([("0", "0")], {0: "a"}, multimer=True)
    assert assigned == {"a": "a"} and problems == {}


def test_clusters_are_numbered_by_size_then_input_order(pipe):
    assigned = {"a": "a", "b": "c", "c": "c", "d": "d"}
    assert pipe.number_clusters(assigned, ["a", "b", "c", "d"]) == {"c": 1, "a": 2, "d": 3}


def _structures_json(tmp_path, ids):
    files = []
    for sid in ids:
        path = tmp_path / f"{sid}.pdb"
        path.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n")
        files.append(str(path))
    stream = {"name": "structures", "ids": list(ids), "files": files, "map_table": "", "format": "pdb"}
    out = tmp_path / "structures.json"
    out.write_text(json.dumps(stream))
    return str(out)


def _run_main(pipe, monkeypatch, tmp_path, ids, fake, cfg=None):
    monkeypatch.setattr(pipe.os, "symlink", lambda src, dst: shutil.copyfile(src, dst))
    monkeypatch.setattr(pipe, "run_foldseek", fake)
    config = tmp_path / "config.yaml"
    config.write_text(yaml.safe_dump(cfg or CFG))
    out = tmp_path / "tables" / "tm_clusters.csv"
    local = tmp_path / "_execution" / "local_missing.csv"
    argv = ["pipe_foldseek.py", "--config-yaml", str(config),
            "--structures-json", _structures_json(tmp_path, ids),
            "--scratch-dir", str(tmp_path / "_execution"),
            "--tm-clusters-csv", str(out), "--local-missing-csv", str(local)]
    monkeypatch.setattr(sys, "argv", argv)
    pipe.main()
    return pd.read_csv(out), pd.read_csv(local)


def _fake_cluster(groups_by_tm):
    """A stand-in for the binary: writes <prefix>_cluster.tsv from {tm: [[member indices]]}, first member leading."""
    def fake(cmd):
        tm = float(cmd[cmd.index("--tmscore-threshold") + 1])
        prefix = cmd[3]
        with open(f"{prefix}_cluster.tsv", "w") as f:
            for group in groups_by_tm[tm]:
                for member in group:
                    f.write(f"{group[0]}\t{member}\n")
    return fake


def test_the_table_holds_one_row_per_structure_with_each_threshold(pipe, monkeypatch, tmp_path):
    fake = _fake_cluster({0.5: [[0, 1], [2]], 0.8: [[0], [1], [2]]})
    table, missing = _run_main(pipe, monkeypatch, tmp_path, ["a", "b", "c"], fake)
    assert list(table.columns) == ["id", "cluster_tm50", "representative_tm50", "n_clusters_tm50",
                                   "cluster_tm80", "representative_tm80", "n_clusters_tm80"]
    assert table["id"].tolist() == ["a", "b", "c"]
    assert table["cluster_tm50"].tolist() == [1, 1, 2]
    assert table["representative_tm50"].tolist() == ["a", "a", "c"]
    assert table["n_clusters_tm50"].tolist() == [2, 2, 2]
    assert table["cluster_tm80"].tolist() == [1, 2, 3]
    assert table["n_clusters_tm80"].tolist() == [3, 3, 3]
    assert missing.empty


def test_an_unreported_structure_keeps_its_row_with_nan_and_a_failure(pipe, monkeypatch, tmp_path):
    fake = _fake_cluster({0.5: [[0, 1]], 0.8: [[0], [1]]})
    table, missing = _run_main(pipe, monkeypatch, tmp_path, ["a", "b", "c"], fake)
    assert table["id"].tolist() == ["a", "b", "c"]
    assert table.loc[2, ["cluster_tm50", "representative_tm50", "n_clusters_tm50"]].isna().all()
    assert missing["id"].tolist() == ["c"] and missing["kind"].tolist() == ["failure"]
    assert "TM 0.5" in missing.loc[0, "cause"] and "TM 0.8" in missing.loc[0, "cause"]


def test_a_failing_threshold_fails_the_step_but_keeps_the_other_columns(pipe, monkeypatch, tmp_path):
    good = _fake_cluster({0.5: [[0, 1]]})

    def fake(cmd):
        if "0.8" in cmd:
            raise RuntimeError("foldseek exited 1: boom")
        good(cmd)

    with pytest.raises(SystemExit) as raised:
        _run_main(pipe, monkeypatch, tmp_path, ["a", "b"], fake)
    assert raised.value.code == 1
    table = pd.read_csv(tmp_path / "tables" / "tm_clusters.csv")
    missing = pd.read_csv(tmp_path / "_execution" / "local_missing.csv")
    assert table["cluster_tm50"].tolist() == [1, 1]
    assert table["cluster_tm80"].isna().all()
    assert sorted(missing["id"]) == ["a", "b"]


def test_every_threshold_failing_is_a_step_failure(pipe, monkeypatch, tmp_path):
    def fake(cmd):
        raise RuntimeError("foldseek exited 1: boom")

    with pytest.raises(SystemExit) as raised:
        _run_main(pipe, monkeypatch, tmp_path, ["a", "b"], fake)
    assert raised.value.code == 1
    assert pd.read_csv(tmp_path / "tables" / "tm_clusters.csv")["id"].tolist() == ["a", "b"]


def test_every_threshold_run_receives_the_alignment_filters(pipe, monkeypatch, tmp_path):
    seen = []
    write = _fake_cluster({0.5: [[0]], 0.8: [[0]]})

    def fake(cmd):
        seen.append(cmd)
        write(cmd)

    cfg = dict(CFG, evalue=0.5, coverage=0.7, cov_mode=2)
    _run_main(pipe, monkeypatch, tmp_path, ["a"], fake, cfg)
    assert len(seen) == 2
    for cmd in seen:
        assert cmd[cmd.index("-e") + 1] == "0.5" and cmd[cmd.index("-c") + 1] == "0.7"
        assert cmd[cmd.index("--cov-mode") + 1] == "2"


def test_the_raw_cluster_tsv_stays_in_the_execution_folder(pipe, monkeypatch, tmp_path):
    fake = _fake_cluster({0.5: [[0]], 0.8: [[0]]})
    _run_main(pipe, monkeypatch, tmp_path, ["a"], fake)
    kept = glob.glob(str(tmp_path / "_execution" / "foldseek_tm*" / "res_cluster.tsv"))
    assert len(kept) == 2


# ── wrapper ───────────────────────────────────────────────────────────────────

def _structures():
    from biopipelines.mock import Mock
    return Mock(ids=["s1", "s2"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})


def test_default_parameters_declare_the_chaperone_thresholds(local_config, isolated_cwd, new_pipeline):
    from biopipelines.foldseek import Foldseek

    pipeline = new_pipeline("fs_defaults")
    with pipeline:
        out = Foldseek(structures=_structures())
    assert pipeline.tools[-1].tm_thresholds == [0.5, 0.6, 0.7, 0.8]
    columns = out.tables.tm_clusters.info.columns
    assert columns[0] == "id" and len(columns) == 1 + 3 * 4
    assert "cluster_tm50" in columns and "n_clusters_tm80" in columns and "representative_tm60" in columns


def test_the_output_columns_follow_the_thresholds(local_config, isolated_cwd, new_pipeline):
    from biopipelines.foldseek import Foldseek

    with new_pipeline("fs_columns"):
        out = Foldseek(structures=_structures(), tm_thresholds=[0.75, 0.55])
    assert out.tables.tm_clusters.info.columns == [
        "id", "cluster_tm55", "representative_tm55", "n_clusters_tm55",
        "cluster_tm75", "representative_tm75", "n_clusters_tm75"]


@pytest.mark.parametrize("kwargs, message", [
    ({"tm_thresholds": []}, "at least one"),
    ({"tm_thresholds": [0.0]}, r"in \(0, 1\)"),
    ({"tm_thresholds": [1.0]}, r"in \(0, 1\)"),
    ({"tm_thresholds": [0.555]}, "two decimals"),
    ({"tm_thresholds": [0.5, 0.5]}, "distinct"),
    ({"alignment_type": 0}, "alignment_type"),
    ({"alignment_type": 3}, "alignment_type"),
    ({"min_seq_id": 1.5}, "min_seq_id"),
    ({"evalue": 0}, "evalue"),
    ({"evalue": -1e-3}, "evalue"),
    ({"evalue": float("inf")}, "evalue"),
    ({"evalue": float("nan")}, "evalue"),
    ({"coverage": -0.1}, "coverage"),
    ({"coverage": 1.1}, "coverage"),
    ({"coverage": float("nan")}, "coverage"),
    ({"cov_mode": -1}, "cov_mode"),
    ({"cov_mode": 6}, "cov_mode"),
    ({"multimer": True, "cov_mode": 3}, "cov_mode"),
    ({"multimer": True, "cov_mode": 5}, "cov_mode"),
    ({"threads": 0}, "threads"),
    ({"multimer": True, "min_seq_id": 0.3}, "single-chain"),
    ({"multimer_tm_threshold": 0.7}, "multimer=True"),
    ({"interface_lddt_threshold": 0.5}, "multimer=True"),
    ({"multimer": True, "multimer_tm_threshold": 1.2}, "multimer_tm_threshold"),
])
def test_bad_parameters_are_refused_at_configuration(local_config, isolated_cwd, new_pipeline, kwargs, message):
    from biopipelines.foldseek import Foldseek

    with new_pipeline("fs_bad"):
        with pytest.raises(ValueError, match=message):
            Foldseek(structures=_structures(), **kwargs)


@pytest.mark.parametrize("kwargs", [
    {},
    {"evalue": 1e-3, "coverage": 0.9, "cov_mode": 1},
    {"evalue": 10, "coverage": 0.0},
    {"coverage": 1.0, "cov_mode": 5},
    {"multimer": True, "evalue": 0.1, "coverage": 0.8, "cov_mode": 2},
    {"multimer": True, "cov_mode": 0},
])
def test_good_alignment_filters_are_accepted(local_config, isolated_cwd, new_pipeline, kwargs):
    from biopipelines.foldseek import Foldseek

    with new_pipeline("fs_good"):
        Foldseek(structures=_structures(), **kwargs)


def test_the_alignment_filters_default_to_nothing_new(local_config, isolated_cwd, new_pipeline):
    from biopipelines.foldseek import Foldseek

    pipeline = new_pipeline("fs_filter_defaults")
    with pipeline:
        Foldseek(structures=_structures())
    fs = pipeline.tools[-1]
    assert fs.evalue is None and fs.coverage is None and fs.cov_mode == 0


def test_the_filters_are_recorded_in_the_config_and_the_display(local_config, isolated_cwd, new_pipeline):
    from biopipelines.foldseek import Foldseek

    pipeline = new_pipeline("fs_filter_record")
    with pipeline:
        Foldseek(structures=_structures(), evalue=1e-3, coverage=0.9, cov_mode=1)
        Foldseek(structures=_structures())
        pipeline.save()
    set_step, default_step = [t for t in pipeline.tools if t.TOOL_NAME == "Foldseek"]
    config = yaml.safe_load(open(set_step.config_yaml))
    assert (config["evalue"], config["coverage"], config["cov_mode"]) == (1e-3, 0.9, 1)
    config = yaml.safe_load(open(default_step.config_yaml))
    assert (config["evalue"], config["coverage"], config["cov_mode"]) == (None, None, 0)
    shown = "\n".join(set_step.get_config_display())
    assert "E-VALUE: 0.001" in shown and "COVERAGE: 0.9 (cov-mode 1)" in shown
    assert "E-VALUE: Foldseek default" in "\n".join(default_step.get_config_display())


def test_multimer_defaults_resolve_to_the_documented_values(local_config, isolated_cwd, new_pipeline):
    from biopipelines.foldseek import Foldseek

    pipeline = new_pipeline("fs_multimer")
    with pipeline:
        Foldseek(structures=_structures(), multimer=True)
    fs = pipeline.tools[-1]
    assert fs.multimer_tm_threshold == 0.65 and fs.interface_lddt_threshold == 0.65


def _step_text(script_path):
    step = sorted(glob.glob(os.path.join(os.path.dirname(script_path), "*_Foldseek.sh")))[0]
    return open(step, encoding="utf-8", errors="replace").read()


def test_the_emitted_script_runs_the_pipe_with_the_container_prefix(local_config, isolated_cwd, new_pipeline):
    from biopipelines.foldseek import Foldseek

    pipeline = new_pipeline("fs_script")
    with pipeline:
        out = Foldseek(structures=_structures(), tm_thresholds=[0.5], threads=4)
        script_path = pipeline.save()
    text = _step_text(script_path)
    assert "pipe_foldseek.py" in text
    assert "--container-prefix" in text and "--tm-clusters-csv" in text and "--local-missing-csv" in text
    config = yaml.safe_load(open(pipeline.tools[-1].config_yaml))
    assert config["tm_thresholds"] == [0.5] and config["threads"] == 4 and config["multimer"] is False
    assert out.tables.missing.info.columns == ["id", "removed_by", "kind", "cause"]


def test_every_pipe_flag_has_a_matching_argument(local_config, isolated_cwd, new_pipeline):
    import re
    from biopipelines.foldseek import Foldseek

    pipeline = new_pipeline("fs_flags")
    with pipeline:
        Foldseek(structures=_structures())
        script_path = pipeline.save()
    text = _step_text(script_path)
    call = text[text.index("pipe_foldseek.py"):]
    call = call[:call.index("python ", 1)] if "python " in call[1:] else call
    emitted = set(re.findall(r"--[a-z][a-z-]+", call))
    declared = set(re.findall(r'add_argument\("(--[a-z-]+)"', open(SCRIPT, encoding="utf-8").read()))
    assert emitted == declared
