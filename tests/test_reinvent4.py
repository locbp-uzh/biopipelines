"""REINVENT4 parameter validation, predicted outputs, install script and runtime helper."""

import importlib.util
import json
import os
import shutil
import subprocess
import sys

import pandas as pd
import pytest

from biopipelines.reinvent4 import (
    CONDITIONAL_GENERATORS,
    GENERATORS,
    GENERATOR_PRIOR,
    PRIOR_MD5,
    SIMILARITY_PRIOR,
    UPSTREAM_COLUMNS,
    REINVENT4,
    prior_file,
)

compound_columns = REINVENT4.compound_columns
sample_columns = REINVENT4.sample_columns

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HELPER = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_reinvent4.py")

# What `reinvent` writes per generator, header copied from upstream's run_sampling.py HEADERS.
UPSTREAM_HEADERS = {
    "reinvent": "SMILES,SMILES_state,NLL",
    "libinvent": "SMILES,SMILES_state,Scaffold,R-groups,NLL",
    "linkinvent": "SMILES,SMILES_state,Warheads,Linker,NLL",
    "mol2mol": "SMILES,SMILES_state,Input_SMILES,Tanimoto,NLL",
}


def _load_helper():
    spec = importlib.util.spec_from_file_location("pipe_reinvent4", HELPER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _compounds(ids=("s1", "s2")):
    from biopipelines.datastream import DataStream
    return DataStream(name="compounds", ids=list(ids), files=[],
                      map_table="compounds.csv", format="csv")


def _tool(generator="reinvent", with_inputs=None, **kwargs):
    from biopipelines.reinvent4 import REINVENT4
    if with_inputs is None:
        with_inputs = generator in CONDITIONAL_GENERATORS
    compounds = _compounds() if with_inputs else None
    tool = REINVENT4(generator=generator, compounds=compounds, **kwargs)
    tool.folders = {"REINVENT4": "/data/r4", "pipe_scripts": "/repo/pipe_scripts"}
    tool.output_folder = os.path.join(os.sep, "runs", "Job_001", "003_REINVENT4")
    return tool


# ── parameter validation ─────────────────────────────────────────────────────

@pytest.mark.parametrize("generator", GENERATORS)
def test_every_generator_validates_with_its_defaults(generator):
    _tool(generator).validate_params()


def test_unknown_generator_is_rejected():
    with pytest.raises(ValueError, match="generator must be one of"):
        _tool("pepinvent", with_inputs=False).validate_params()


def test_de_novo_refuses_compounds():
    with pytest.raises(ValueError, match="takes no compounds"):
        _tool("reinvent", with_inputs=True).validate_params()


@pytest.mark.parametrize("generator", sorted(CONDITIONAL_GENERATORS))
def test_conditional_generators_need_compounds(generator):
    with pytest.raises(ValueError, match="needs compounds"):
        _tool(generator, with_inputs=False).validate_params()


def test_compounds_must_be_a_stream():
    from biopipelines.reinvent4 import REINVENT4
    with pytest.raises(ValueError, match="DataStream or StandardizedOutput"):
        REINVENT4(generator="libinvent", compounds="[*:0]CC")


@pytest.mark.parametrize("value", [0, -3, True, 2.5, "10"])
def test_num_compounds_must_be_a_positive_integer(value):
    with pytest.raises(ValueError, match="num_compounds must be an integer >= 1"):
        _tool("reinvent", num_compounds=value).validate_params()


@pytest.mark.parametrize("value", [0, -1, True, 1.5])
def test_seed_must_be_a_positive_integer(value):
    """Upstream reads `args.seed or config seed`, so 0 would silently mean no seed."""
    with pytest.raises(ValueError, match="seed must be an integer >= 1"):
        _tool("reinvent", seed=value).validate_params()


@pytest.mark.parametrize("generator", ["reinvent", "libinvent", "linkinvent"])
@pytest.mark.parametrize("kwargs", [{"similarity": "high_similarity"},
                                    {"sample_strategy": "beamsearch"},
                                    {"temperature": 0.8}])
def test_mol2mol_settings_are_refused_for_the_other_generators(generator, kwargs):
    with pytest.raises(ValueError, match="only applies to generator='mol2mol'"):
        _tool(generator, **kwargs).validate_params()


def test_unknown_similarity_and_strategy_are_rejected():
    with pytest.raises(ValueError, match="similarity must be one of"):
        _tool("mol2mol", similarity="low").validate_params()
    with pytest.raises(ValueError, match="sample_strategy must be one of"):
        _tool("mol2mol", sample_strategy="greedy").validate_params()


@pytest.mark.parametrize("value", [0, -1.0, True, "hot"])
def test_temperature_must_be_positive(value):
    with pytest.raises(ValueError, match="temperature must be a number > 0"):
        _tool("mol2mol", temperature=value).validate_params()


def test_temperature_is_refused_with_beam_search():
    with pytest.raises(ValueError, match="multinomial sampling only"):
        _tool("mol2mol", sample_strategy="beamsearch", temperature=0.5).validate_params()


def test_mol2mol_accepts_its_settings():
    tool = _tool("mol2mol", similarity="scaffold", sample_strategy="multinomial", temperature=1.2, seed=3)
    tool.validate_params()


# ── priors ───────────────────────────────────────────────────────────────────

def test_every_selectable_prior_has_a_checksum_and_no_file_is_unused():
    selectable = set(GENERATOR_PRIOR.values()) | set(SIMILARITY_PRIOR.values())
    assert selectable == set(PRIOR_MD5)
    assert all(len(md5) == 32 and set(md5) <= set("0123456789abcdef") for md5 in PRIOR_MD5.values())


@pytest.mark.parametrize("generator, expected", [
    ("reinvent", "reinvent.prior"),
    ("libinvent", "libinvent.prior"),
    ("linkinvent", "linkinvent.prior"),
    ("mol2mol", "mol2mol_medium_similarity.prior"),
])
def test_default_prior_per_generator(generator, expected):
    assert prior_file(generator) == expected


@pytest.mark.parametrize("similarity, expected", [
    ("similarity", "mol2mol_similarity.prior"),
    ("high_similarity", "mol2mol_high_similarity.prior"),
    ("scaffold_generic", "mol2mol_scaffold_generic.prior"),
    ("mmp", "mol2mol_mmp.prior"),
    ("pubchem", "pubchem_ecfp4_with_count_with_rank_reinvent4_dict_voc.prior"),
])
def test_similarity_selects_the_mol2mol_prior(similarity, expected):
    assert prior_file("mol2mol", similarity) == expected


# ── generated config and predicted outputs ───────────────────────────────────

def test_pipe_config_points_at_the_installed_prior():
    cfg = _tool("mol2mol", similarity="mmp", seed=5)._pipe_config()
    assert cfg["model_file"] == os.path.join("/data/r4", "priors", "mol2mol_mmp.prior")
    assert cfg["seed"] == 5
    assert cfg["id_prefix"] is None
    assert cfg["compounds_json"].endswith("compounds_ds.json")
    assert cfg["inputs_csv"].endswith(os.path.join("tables", "inputs.csv"))


def test_pipe_config_names_the_missing_repository_folder():
    tool = _tool("reinvent")
    tool.folders = {"pipe_scripts": "/repo/pipe_scripts"}
    with pytest.raises(ValueError, match="folders.repositories.REINVENT4"):
        tool._pipe_config()


def test_de_novo_config_names_ids_after_the_pipeline_and_reads_no_inputs():
    cfg = _tool("reinvent")._pipe_config()
    assert cfg["id_prefix"] == "Job_001"
    assert cfg["compounds_json"] is None
    assert cfg["inputs_csv"] is None
    assert cfg["sample_strategy"] is None and cfg["temperature"] is None


@pytest.mark.parametrize("generator", GENERATORS)
def test_compounds_stream_is_a_value_based_csv_with_lazy_ids(generator):
    out = _tool(generator).get_output_files()
    stream = out["compounds"]
    assert stream.format == "csv"
    assert list(stream.files) == []
    assert stream.map_table.endswith(os.path.join("compounds", "compounds.csv"))
    if generator == "reinvent":
        assert list(stream.ids) == ["Job_001_[<#>]"]
    else:
        assert list(stream.ids) == ["s1_[<#>]", "s2_[<#>]"]


@pytest.mark.parametrize("generator", GENERATORS)
def test_declared_tables_match_the_helper_columns(generator):
    tables = _tool(generator).get_output_files()["tables"]
    assert tables["compounds"].info.columns == compound_columns(generator)
    assert tables["samples"].info.columns == sample_columns(generator)
    assert tables["compounds"].info.path == _tool(generator).get_output_files()["compounds"].map_table
    assert ("inputs" in tables) == (generator in CONDITIONAL_GENERATORS)
    assert tables["missing"].info.columns == ["id", "removed_by", "kind", "cause"]


def test_canonical_compounds_schema_comes_first():
    for generator in GENERATORS:
        assert compound_columns(generator)[:5] == ["id", "format", "code", "smiles", "ccd"]


@pytest.mark.parametrize("generator, extra", [
    ("reinvent", []),
    ("libinvent", ["scaffold", "r_groups"]),
    ("linkinvent", ["warheads", "linker"]),
    ("mol2mol", ["input_smiles", "tanimoto"]),
])
def test_samples_columns_follow_what_upstream_reports(generator, extra):
    columns = sample_columns(generator)
    assert columns[-len(extra):] == extra if extra else columns[-1] == "nll"
    assert ("compounds.id" in columns) == (generator in CONDITIONAL_GENERATORS)
    assert set(UPSTREAM_COLUMNS[generator].values()) <= set(columns)


@pytest.mark.parametrize("generator", GENERATORS)
def test_upstream_columns_cover_the_header_reinvent_writes(generator):
    written = UPSTREAM_HEADERS[generator].split(",")
    assert set(UPSTREAM_COLUMNS[generator]) == set(written) - {"SMILES_state"}


# ── generated script ─────────────────────────────────────────────────────────

def _script(tmp_path, generator):
    tool = _tool(generator)
    tool.output_folder = str(tmp_path / "Job_001" / "003_REINVENT4")
    os.makedirs(tool.configuration_folder, exist_ok=True)
    script = tool.generate_script(str(tmp_path / "step.sh"))
    with open(tool.config_json) as f:
        return script, json.load(f), tool


@pytest.mark.parametrize("generator", GENERATORS)
def test_script_runs_the_helper_with_its_config(tmp_path, local_config, generator):
    script, cfg, tool = _script(tmp_path, generator)
    assert f'python "{tool.helper_py}"' in script
    assert f'--config "{tool.config_json}"' in script
    assert "--container-prefix" in script
    assert cfg["generator"] == generator
    assert "pipe_propagate_missing.py" in script
    assert os.path.exists(tool.compounds_json) == (generator in CONDITIONAL_GENERATORS)


# ── install ──────────────────────────────────────────────────────────────────

def _install(monkeypatch, force=False, **kwargs):
    from biopipelines.reinvent4 import REINVENT4
    monkeypatch.setattr(REINVENT4, "_install_env", classmethod(lambda cls, mgr: "REINVENT4"))
    return REINVENT4._install_script(
        {"biopipelines": "/repo", "REINVENT4": "/data/r4"},
        env_manager="mamba", force_reinstall=force, **kwargs)


def test_install_uses_upstreams_installer_without_optional_dependencies(local_config, monkeypatch):
    script = _install(monkeypatch)
    assert "python install.py cu126 -d none" in script
    assert 'cd "/data/r4"' in script


def test_install_processor_is_a_keyword(local_config, monkeypatch):
    assert "python install.py cpu -d none" in _install(monkeypatch, processor="cpu")


@pytest.mark.parametrize("processor", ["cu126; rm -rf ~", "", "cu 126", None, 126])
def test_install_processor_cannot_carry_shell(local_config, monkeypatch, processor):
    with pytest.raises(ValueError, match="processor must be a PyTorch build tag"):
        _install(monkeypatch, processor=processor)


def test_install_verifies_the_cli_starts_and_every_prior_is_present(local_config, monkeypatch):
    script = _install(monkeypatch)
    assert script.count("reinvent --version") == 2
    for name in PRIOR_MD5:
        assert script.count(f'[ -s "/data/r4/priors/{name}" ]') == 2


def test_install_pins_the_zenodo_record_and_checks_every_md5(local_config, monkeypatch):
    script = _install(monkeypatch)
    assert "https://zenodo.org/api/records/15641297/files" in script
    for name, md5 in PRIOR_MD5.items():
        assert f'fetch_prior "{name}" "{md5}" || exit 1' in script
    assert "checksum mismatch" in script


def test_install_keeps_an_existing_env_so_a_failed_download_can_be_retried(local_config, monkeypatch):
    script = _install(monkeypatch)
    guard = script.split("not recreating it")[1].split("cd \"/data/r4\"")[0]
    assert "env create" in guard
    assert script.index("not recreating it") < script.index("env create")


def test_forced_reinstall_removes_the_env_before_the_guard(local_config, monkeypatch):
    script = _install(monkeypatch, force=True)
    assert "already installed, skipping" not in script
    assert script.index("env remove") < script.index("not recreating it")


@pytest.mark.parametrize("force", [False, True])
def test_install_script_is_valid_bash(local_config, monkeypatch, tmp_path, force):
    bash = shutil.which("bash")
    if bash is None:
        pytest.skip("bash is not available")
    path = tmp_path / "install.sh"
    path.write_text(_install(monkeypatch, force=force), encoding="utf-8", newline="\n")
    result = subprocess.run([bash, "-n", str(path)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


def test_install_probe_without_folders_does_not_raise(local_config, monkeypatch):
    """BaseConfig.install() calls _install_script({}, manager) just to see that it is defined."""
    from biopipelines.reinvent4 import REINVENT4
    monkeypatch.setattr(REINVENT4, "_install_env", classmethod(lambda cls, mgr: "REINVENT4"))
    assert REINVENT4._install_script({}, "mamba")


# ── runtime helper: configuration, parsing ───────────────────────────────────

def _cfg(**overrides):
    cfg = {"generator": "reinvent", "model_file": "/data/r4/priors/reinvent.prior",
           "num_compounds": 20, "sample_strategy": None, "temperature": None, "seed": None}
    cfg.update(overrides)
    return cfg


def test_toml_literals():
    helper = _load_helper()
    assert helper.toml_value(True) == "true"
    assert helper.toml_value(False) == "false"
    assert helper.toml_value(20) == "20"
    assert helper.toml_value(0.5) == "0.5"
    assert helper.toml_value('a"b\\c') == '"a\\"b\\\\c"'


def test_sampling_toml_round_trips_for_de_novo():
    tomllib = pytest.importorskip("tomllib")
    text = _load_helper().build_sampling_toml(_cfg(), None, "/run/denovo.csv", "cpu")
    parsed = tomllib.loads(text)
    assert parsed == {
        "run_type": "sampling",
        "device": "cpu",
        "parameters": {"model_file": "/data/r4/priors/reinvent.prior", "output_file": "/run/denovo.csv",
                       "num_smiles": 20, "unique_molecules": True},
    }


def test_sampling_toml_carries_the_input_and_mol2mol_settings():
    tomllib = pytest.importorskip("tomllib")
    cfg = _cfg(generator="mol2mol", sample_strategy="multinomial", temperature=2)
    parsed = tomllib.loads(_load_helper().build_sampling_toml(cfg, "/run/s1.smi", "/run/s1.csv", "cuda:0"))
    params = parsed["parameters"]
    assert parsed["device"] == "cuda:0"
    assert params["smiles_file"] == "/run/s1.smi"
    assert params["sample_strategy"] == "multinomial"
    assert params["temperature"] == 2.0 and isinstance(params["temperature"], float)


def test_unset_mol2mol_settings_are_left_to_upstream():
    tomllib = pytest.importorskip("tomllib")
    parsed = tomllib.loads(_load_helper().build_sampling_toml(_cfg(generator="mol2mol"), "/a.smi", "/a.csv", "cpu"))
    assert "sample_strategy" not in parsed["parameters"]
    assert "temperature" not in parsed["parameters"]


def _write_csv(path, generator, rows):
    with open(path, "w") as f:
        f.write(UPSTREAM_HEADERS[generator] + "\n")
        for row in rows:
            f.write(",".join(row) + "\n")


def test_read_samples_maps_each_generator_to_the_emitted_columns(tmp_path):
    helper = _load_helper()
    cases = {
        "reinvent": (["CCO", "1", "9.5"], {"smiles": "CCO", "nll": 9.5}),
        "libinvent": (["c1ccccc1C", "1", "[*]c1ccccc1", "[*]C|[*]N", "12.25"],
                      {"smiles": "c1ccccc1C", "nll": 12.25, "scaffold": "[*]c1ccccc1", "r_groups": "[*]C|[*]N"}),
        "linkinvent": (["CCOCC", "1", "*C|*C", "*O*", "7.0"],
                       {"smiles": "CCOCC", "nll": 7.0, "warheads": "*C|*C", "linker": "*O*"}),
        "mol2mol": (["CCN", "1", "CCO", "0.5", "3.75"],
                    {"smiles": "CCN", "nll": 3.75, "input_smiles": "CCO", "tanimoto": 0.5}),
    }
    for generator, (record, expected) in cases.items():
        path = tmp_path / f"{generator}.csv"
        _write_csv(path, generator, [record])
        rows = helper.read_samples(str(path), generator)
        assert rows == [expected]
        assert set(rows[0]) == set(UPSTREAM_COLUMNS[generator].values())


def test_read_samples_keeps_smiles_that_pandas_would_read_as_missing(tmp_path):
    helper = _load_helper()
    path = tmp_path / "r.csv"
    _write_csv(path, "reinvent", [["N", "1", "4.0"], ["NA", "1", "5.0"]])
    assert [r["smiles"] for r in helper.read_samples(str(path), "reinvent")] == ["N", "NA"]


def test_read_samples_header_only_is_no_rows(tmp_path):
    path = tmp_path / "r.csv"
    _write_csv(path, "libinvent", [])
    assert _load_helper().read_samples(str(path), "libinvent") == []


def test_read_samples_rejects_a_header_it_does_not_know(tmp_path):
    path = tmp_path / "r.csv"
    path.write_text("SMILES,NLL\nCCO,1.0\n")
    with pytest.raises(ValueError, match="lacks the column"):
        _load_helper().read_samples(str(path), "mol2mol")


def test_read_samples_rejects_a_non_numeric_nll(tmp_path):
    path = tmp_path / "r.csv"
    _write_csv(path, "reinvent", [["CCO", "1", "high"]])
    with pytest.raises(ValueError):
        _load_helper().read_samples(str(path), "reinvent")


def test_input_text():
    helper = _load_helper()
    assert helper.input_text("  [*:0]CC[*:1] ") == "[*:0]CC[*:1]"
    assert helper.input_text("") is None
    assert helper.input_text(None) is None
    assert helper.input_text(float("nan")) is None


def test_last_line_of_a_log(tmp_path):
    helper = _load_helper()
    log = tmp_path / "x.log"
    log.write_text("start\nValueError: Tokens {'Q'} are not supported\n\n")
    assert helper.last_line(str(log)).startswith("ValueError: Tokens")
    assert helper.last_line(str(tmp_path / "absent.log")) == "no log"


def test_cuda_probe_reads_the_answer_and_refuses_noise(monkeypatch):
    helper = _load_helper()

    class Done:
        def __init__(self, stdout):
            self.stdout = stdout

    monkeypatch.setattr(helper.subprocess, "run", lambda *a, **k: Done("warning\nTrue\n"))
    assert helper.cuda_available("") is True
    monkeypatch.setattr(helper.subprocess, "run", lambda *a, **k: Done("False\n"))
    assert helper.cuda_available("") is False
    monkeypatch.setattr(helper.subprocess, "run", lambda *a, **k: Done("maybe\n"))
    with pytest.raises(RuntimeError, match="could not tell"):
        helper.cuda_available("")


# ── runtime helper: a whole run against a stand-in for `reinvent` ────────────

def _fake_reinvent(helper, monkeypatch, behavior):
    """Replace the subprocess calls: the CUDA probe, and `reinvent` writing its CSV from the TOML."""
    tomllib = pytest.importorskip("tomllib")
    calls = []

    class Done:
        def __init__(self, returncode=0, stdout=""):
            self.returncode = returncode
            self.stdout = stdout

    def run(cmd, **kwargs):
        if "-c" in cmd:
            return Done(stdout="False\n")
        calls.append(list(cmd))
        with open(cmd[-1], "rb") as f:
            params = tomllib.load(f)["parameters"]
        smiles_file = params.get("smiles_file")
        seed_smiles = open(smiles_file).read().strip() if smiles_file else None
        kwargs["stdout"].write(f"log for {seed_smiles}\n")
        result = behavior(seed_smiles, params)
        if isinstance(result, int):
            kwargs["stdout"].write("RuntimeError: boom\n")
            return Done(returncode=result)
        _write_csv(params["output_file"], helper_generator[0], result)
        return Done()

    helper_generator = ["reinvent"]
    monkeypatch.setattr(helper.subprocess, "run", run)
    return calls, helper_generator


def _run_helper(tmp_path, monkeypatch, generator, behavior, inputs=None, extra=None):
    helper = _load_helper()
    calls, holder = _fake_reinvent(helper, monkeypatch, behavior)
    holder[0] = generator

    out = tmp_path / "Job_001" / "003_REINVENT4"
    paths = {
        "compounds_csv": str(out / "compounds" / "compounds.csv"),
        "samples_csv": str(out / "tables" / "samples.csv"),
        "inputs_csv": str(out / "tables" / "inputs.csv"),
        "local_missing_csv": str(out / "execution" / "local_missing.csv"),
        "sampling_dir": str(out / "execution" / "sampling"),
    }
    cfg = _cfg(generator=generator, **paths)
    cfg["compounds_json"] = None
    cfg["id_prefix"] = "Job_001"
    if generator not in CONDITIONAL_GENERATORS:
        cfg["inputs_csv"] = None
    if inputs is not None:
        from biopipelines.datastream import DataStream
        map_csv = tmp_path / "in_compounds.csv"
        pd.DataFrame({"id": list(inputs), "smiles": list(inputs.values())}).to_csv(map_csv, index=False)
        DataStream(name="compounds", ids=list(inputs), files=[], map_table=str(map_csv),
                   format="csv").save_json(str(tmp_path / "in.json"))
        cfg["compounds_json"] = str(tmp_path / "in.json")
    cfg.update(extra or {})
    cfg_path = tmp_path / "cfg.json"
    cfg_path.write_text(json.dumps(cfg))

    monkeypatch.setattr(sys, "argv", ["pipe_reinvent4.py", "--config", str(cfg_path)])
    code = 0
    try:
        helper.main()
    except SystemExit as e:
        code = e.code
    return code, paths, calls


def test_de_novo_run_writes_the_stream_and_samples(tmp_path, monkeypatch):
    rows = [["CCO", "1", "8.5"], ["CCN", "1", "9.0"], ["CCC", "1", "9.25"]]
    code, paths, calls = _run_helper(tmp_path, monkeypatch, "reinvent", lambda s, p: rows)
    assert code == 0
    compounds = pd.read_csv(paths["compounds_csv"], keep_default_na=False)
    assert list(compounds.columns) == compound_columns("reinvent")
    assert list(compounds["id"]) == ["Job_001_1", "Job_001_2", "Job_001_3"]
    assert set(compounds["format"]) == {"smiles"} and set(compounds["code"]) == {"LIG"}
    samples = pd.read_csv(paths["samples_csv"])
    assert list(samples.columns) == sample_columns("reinvent")
    assert list(samples["nll"]) == [8.5, 9.0, 9.25]
    assert not os.path.exists(paths["inputs_csv"])
    assert len(calls) == 1 and calls[0][0] == "reinvent"
    assert pd.read_csv(paths["local_missing_csv"]).empty


def test_seed_flag_is_passed_through(tmp_path, monkeypatch):
    code, _, calls = _run_helper(tmp_path, monkeypatch, "reinvent",
                                 lambda s, p: [["CCO", "1", "1.0"]], extra={"seed": 11})
    assert code == 0
    assert calls[0][:3] == ["reinvent", "--seed", "11"]


def test_de_novo_with_no_valid_molecule_fails_the_step(tmp_path, monkeypatch):
    code, paths, _ = _run_helper(tmp_path, monkeypatch, "reinvent", lambda s, p: [])
    assert code == 1
    assert pd.read_csv(paths["compounds_csv"]).empty


def test_de_novo_reinvent_failure_fails_the_step(tmp_path, monkeypatch):
    code, _, _ = _run_helper(tmp_path, monkeypatch, "reinvent", lambda s, p: 1)
    assert code == 1


def test_each_input_keeps_its_own_ids_and_provenance(tmp_path, monkeypatch):
    def behavior(seed, params):
        assert params["smiles_file"].endswith(".smi")
        return {"CCO": [["CCN", "1", "CCO", "0.7", "2.5"], ["CCC", "1", "CCO", "0.6", "3.5"]],
                "c1ccccc1": [["c1ccncc1", "1", "c1ccccc1", "0.55", "4.25"]]}[seed]

    code, paths, calls = _run_helper(tmp_path, monkeypatch, "mol2mol", behavior,
                                     inputs={"a": "CCO", "b": "c1ccccc1"})
    assert code == 0
    assert len(calls) == 2
    compounds = pd.read_csv(paths["compounds_csv"], keep_default_na=False)
    assert list(compounds["id"]) == ["a_1", "a_2", "b_1"]
    assert list(compounds["compounds.id"]) == ["a", "a", "b"]
    samples = pd.read_csv(paths["samples_csv"])
    assert list(samples.columns) == sample_columns("mol2mol")
    assert list(samples["tanimoto"]) == [0.7, 0.6, 0.55]
    inputs = pd.read_csv(paths["inputs_csv"])
    assert list(inputs["id"]) == ["a", "b"] and list(inputs["num_compounds"]) == [2, 1]
    assert pd.read_csv(paths["local_missing_csv"]).empty


def test_one_failing_input_fails_alone(tmp_path, monkeypatch):
    def behavior(seed, params):
        return 1 if seed == "bad" else [["CCN", "1", "CCO", "0.7", "2.5"]]

    code, paths, _ = _run_helper(tmp_path, monkeypatch, "mol2mol", behavior,
                                 inputs={"a": "CCO", "b": "bad"})
    assert code == 0
    assert list(pd.read_csv(paths["compounds_csv"])["id"]) == ["a_1"]
    inputs = pd.read_csv(paths["inputs_csv"])
    assert list(inputs["id"]) == ["a", "b"]
    assert inputs.loc[1, "num_compounds"] != inputs.loc[1, "num_compounds"]
    missing = pd.read_csv(paths["local_missing_csv"])
    assert list(missing["id"]) == ["b"] and list(missing["kind"]) == ["failure"]
    assert "boom" in missing.loc[0, "cause"]
    assert missing.loc[0, "removed_by"] == "003_REINVENT4"


def test_an_input_with_no_valid_molecule_is_a_filter_not_a_failure(tmp_path, monkeypatch):
    def behavior(seed, params):
        return [] if seed == "CCO" else [["CCN", "1", "x", "0.7", "2.5"]]

    code, paths, _ = _run_helper(tmp_path, monkeypatch, "mol2mol", behavior,
                                 inputs={"a": "CCO", "b": "CCC"})
    assert code == 0
    inputs = pd.read_csv(paths["inputs_csv"])
    assert list(inputs["num_compounds"]) == [0, 1]
    missing = pd.read_csv(paths["local_missing_csv"])
    assert list(missing["id"]) == ["a"] and list(missing["kind"]) == ["filter"]


def test_all_inputs_failing_fails_the_step_but_keeps_the_tables(tmp_path, monkeypatch):
    code, paths, _ = _run_helper(tmp_path, monkeypatch, "libinvent", lambda s, p: 2,
                                 inputs={"a": "[*:0]CC[*:1]", "b": "[*:0]CCC[*:1]"})
    assert code == 1
    assert len(pd.read_csv(paths["inputs_csv"])) == 2
    assert set(pd.read_csv(paths["local_missing_csv"])["kind"]) == {"failure"}


def test_every_input_yielding_nothing_still_exits_zero(tmp_path, monkeypatch):
    code, paths, _ = _run_helper(tmp_path, monkeypatch, "linkinvent", lambda s, p: [],
                                 inputs={"a": "*C|*N"})
    assert code == 0
    assert set(pd.read_csv(paths["local_missing_csv"])["kind"]) == {"filter"}


def test_input_without_smiles_is_dropped_and_none_usable_fails_the_step(tmp_path, monkeypatch):
    code, paths, calls = _run_helper(tmp_path, monkeypatch, "mol2mol", lambda s, p: [],
                                     inputs={"a": ""})
    assert code == 1
    assert calls == []
    missing = pd.read_csv(paths["local_missing_csv"])
    assert list(missing["kind"]) == ["filter"] and "no smiles" in missing.loc[0, "cause"]
    assert list(pd.read_csv(paths["inputs_csv"])["id"]) == ["a"]


def test_input_with_whitespace_is_a_failure_that_never_reaches_reinvent(tmp_path, monkeypatch):
    code, paths, calls = _run_helper(tmp_path, monkeypatch, "mol2mol",
                                     lambda s, p: [["CCN", "1", "CCO", "0.7", "2.5"]],
                                     inputs={"a": "CCO", "b": "CC O"})
    assert code == 0
    assert len(calls) == 1
    missing = pd.read_csv(paths["local_missing_csv"])
    assert list(missing["id"]) == ["b"] and list(missing["kind"]) == ["failure"]
