"""AiZynthFinder parameter validation, generated config and result-parsing regressions."""

import copy
import gzip
import importlib.util
import json
import os
import sys

import pandas as pd
import pytest
import yaml

from biopipelines.aizynthfinder import _as_name_list, _scorer_column

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HELPER = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_aizynthfinder.py")


def _load_helper():
    spec = importlib.util.spec_from_file_location("pipe_aizynthfinder", HELPER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _compounds(ids=("c1", "c2")):
    from biopipelines.datastream import DataStream
    return DataStream(name="compounds", ids=list(ids), files=[],
                      map_table="compounds.csv", format="csv")


def _tool(**kwargs):
    from biopipelines.aizynthfinder import AiZynthFinder
    tool = AiZynthFinder(compounds=_compounds(), **kwargs)
    tool.folders = {"AiZynthFinder": "/data/azf"}
    tool.output_folder = "/runs/001_AiZynthFinder"
    return tool


@pytest.mark.parametrize("param", ["stocks", "expansion_policy", "filter_policy"])
def test_bare_string_rejected_instead_of_split_into_characters(param):
    with pytest.raises(ValueError, match="must be a list of names, not a string"):
        _as_name_list(param, "zinc")


def test_non_sequence_rejected():
    with pytest.raises(ValueError, match="must be a list or tuple"):
        _as_name_list("stocks", 5)


def test_lists_and_none_pass_through():
    assert _as_name_list("stocks", ["zinc", "emol"]) == ["zinc", "emol"]
    assert _as_name_list("stocks", ("zinc",)) == ["zinc"]
    assert _as_name_list("stocks", None) is None
    assert _as_name_list("stocks", []) is None


def test_missing_flag_is_not_in_stock():
    """bool(nan) is True, so an absent in_stock would mark a route solved."""
    as_bool = _load_helper()._as_bool
    assert as_bool(float("nan")) is False
    assert as_bool(None) is False


def test_present_flags_still_normalize():
    as_bool = _load_helper()._as_bool
    assert as_bool(True) is True
    assert as_bool("True") is True
    assert as_bool("no") is False
    assert as_bool(0) is False


# ── install ──────────────────────────────────────────────────────────────────

def _install(monkeypatch, force=False):
    from biopipelines.aizynthfinder import AiZynthFinder
    monkeypatch.setattr(AiZynthFinder, "_install_env", classmethod(lambda cls, mgr: "aizynthfinder"))
    return AiZynthFinder._install_script(
        {"biopipelines": "/repo", "AiZynthFinder": "/data/azf"},
        env_manager="mamba", force_reinstall=force)


def test_install_verifies_the_cli_starts_not_just_the_package(local_config, monkeypatch):
    """`import aizynthfinder` passes with a broken dependency chain; aizynthcli imports it all."""
    script = _install(monkeypatch)
    assert script.count("aizynthcli --help") == 2
    assert 'import aizynthfinder"' not in script
    assert "download_public_data ." in script


def test_install_keeps_an_existing_env_so_a_failed_download_can_be_retried(local_config, monkeypatch):
    script = _install(monkeypatch)
    guard = script.split("not recreating it")[1].split("mkdir -p")[0]
    assert "env create" in guard
    assert script.index("not recreating it") < script.index("env create")


@pytest.mark.parametrize("force", [False, True])
def test_install_script_is_valid_bash(local_config, monkeypatch, tmp_path, force):
    import shutil
    import subprocess
    bash = shutil.which("bash")
    if bash is None:
        pytest.skip("bash is not available")
    path = tmp_path / "install.sh"
    path.write_text(_install(monkeypatch, force=force), encoding="utf-8", newline="\n")
    result = subprocess.run([bash, "-n", str(path)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


def test_forced_reinstall_removes_the_env_before_the_guard(local_config, monkeypatch):
    script = _install(monkeypatch, force=True)
    assert "already installed, skipping" not in script
    assert script.index("env remove") < script.index("not recreating it")


# ── parameter validation ─────────────────────────────────────────────────────

def test_defaults_add_nothing_to_the_generated_config():
    cfg = _tool()._pipe_config()
    assert cfg["base_config"] == os.path.join("/data/azf", "config.yml")
    for key in ("search_algorithm", "max_transforms", "route_scorer", "filter_cutoff",
                "frozen_substructures", "stock_stop_criteria"):
        assert cfg[key] is None
    assert cfg["extra_scorers"] == []


@pytest.mark.parametrize("algorithm, path", [
    ("mcts", "mcts"),
    ("retrostar", "aizynthfinder.search.retrostar.search_tree.SearchTree"),
    ("breadth_first", "aizynthfinder.search.breadth_first.search_tree.SearchTree"),
    ("dfpn", "aizynthfinder.search.dfpn.search_tree.SearchTree"),
])
def test_algorithm_resolves_to_the_class_upstream_loads(algorithm, path):
    tool = _tool(search_algorithm=algorithm)
    tool.validate_params()
    assert tool._pipe_config()["search_algorithm"] == path


def test_unknown_algorithm_is_refused():
    with pytest.raises(ValueError, match="search_algorithm"):
        _tool(search_algorithm="astar").validate_params()


@pytest.mark.parametrize("kwargs", [
    {"filter_policy": ["uspto"]},
    {"filter_cutoff": 0.1},
    {"frozen_substructures": ["c1ccccc1"]},
])
def test_dfpn_refuses_every_filter_option(kwargs):
    with pytest.raises(ValueError, match="dfpn"):
        _tool(search_algorithm="dfpn", **kwargs).validate_params()


@pytest.mark.parametrize("value", [0, -1, 2.5, True, "6"])
def test_max_transforms_must_be_a_positive_integer(value):
    with pytest.raises(ValueError, match="max_transforms"):
        _tool(max_transforms=value).validate_params()


@pytest.mark.parametrize("value", [0, -0.1, 1.5, True, "0.05"])
def test_filter_cutoff_must_lie_in_zero_one(value):
    with pytest.raises(ValueError, match="filter_cutoff"):
        _tool(filter_cutoff=value).validate_params()


def test_frozen_substructures_reject_a_bare_string_and_blanks():
    with pytest.raises(ValueError, match="not a string"):
        _tool(frozen_substructures="c1ccccc1")
    with pytest.raises(ValueError, match="frozen_substructures"):
        _tool(frozen_substructures=["c1ccccc1", " "]).validate_params()


def test_scorers_come_from_the_catalog_unless_a_config_defines_more():
    _tool(scorers=["route cost", "sum of prices"]).validate_params()
    with pytest.raises(ValueError, match="own config file"):
        _tool(scorers=["my scorer"]).validate_params()
    _tool(scorers=["my scorer"], config="my.yml").validate_params()


def test_always_computed_scorers_are_refused_as_extra_columns():
    with pytest.raises(ValueError, match="always computed"):
        _tool(scorers=["number of reactions"]).validate_params()


def test_scorer_bare_string_and_duplicates_are_refused():
    with pytest.raises(ValueError, match="not a string"):
        _tool(scorers="route cost")
    with pytest.raises(ValueError, match="more than once"):
        _tool(scorers=["route cost", "route cost"]).validate_params()


def test_scorer_column_may_not_shadow_an_existing_routes_column():
    with pytest.raises(ValueError, match="already exists"):
        _tool(scorers=["Score"], config="my.yml").validate_params()


def test_route_scorer_accepts_the_always_computed_scorers_without_a_new_column():
    tool = _tool(route_scorer="number of reactions")
    tool.validate_params()
    assert tool.extra_scorer_names == []


def test_route_scorer_outside_the_catalog_needs_a_config():
    with pytest.raises(ValueError, match="route_scorer"):
        _tool(route_scorer="my scorer").validate_params()
    _tool(route_scorer="my scorer", config="my.yml").validate_params()


def test_route_scorer_gets_its_own_column_once():
    tool = _tool(route_scorer="route cost", scorers=["route cost", "sum of prices"])
    assert tool.extra_scorer_names == ["route cost", "sum of prices"]
    tool = _tool(route_scorer="route cost", scorers=["sum of prices"])
    assert tool.extra_scorer_names == ["sum of prices", "route cost"]


def test_scorer_column_spelling():
    assert _scorer_column("route cost") == "route_cost"
    assert _scorer_column("average template occurrence") == "average_template_occurrence"
    assert _scorer_column("number of pre-cursors") == "number_of_pre_cursors"


@pytest.mark.parametrize("criteria", [
    {"mass": 250},
    {"weight": 0},
    {"weight": "250"},
    {"counts": {}},
    {"counts": {"C": -1}},
    {"counts": {"6": 3}},
    {"price": True},
    {},
    "weight",
])
def test_bad_stock_stop_criteria_are_refused(criteria):
    with pytest.raises(ValueError, match="stock_stop_criteria"):
        _tool(stock_stop_criteria=criteria).validate_params()


def test_good_stock_stop_criteria_pass():
    _tool(stock_stop_criteria={"weight": 250.5, "price": 10, "counts": {"C": 12, "Cl": 0}}).validate_params()


def test_existing_range_checks_still_apply():
    with pytest.raises(ValueError, match="min_routes"):
        _tool(min_routes=30, max_routes=25).validate_params()
    with pytest.raises(ValueError, match="nproc"):
        _tool(nproc=0).validate_params()


def test_routes_table_declares_one_column_per_extra_scorer():
    columns = _tool(scorers=["route cost"], route_scorer="average template occurrence"
                    ).get_output_files()["tables"]["routes"].info.columns
    assert columns[:8] == ["id", "compounds.id", "route", "score", "is_solved", "number_of_steps",
                           "number_of_precursors", "number_of_precursors_in_stock"]
    assert columns[8:] == ["route_cost", "average_template_occurrence"]


def test_existing_tables_keep_their_columns_by_default():
    from biopipelines.aizynthfinder import ROUTE_COLUMNS, TARGET_COLUMNS, PRECURSOR_COLUMNS
    tables = _tool().get_output_files()["tables"]
    assert tables["routes"].info.columns == ROUTE_COLUMNS
    assert tables["targets"].info.columns == TARGET_COLUMNS
    assert tables["precursors"].info.columns == PRECURSOR_COLUMNS


def test_pipe_config_carries_the_resolved_options():
    cfg = _tool(search_algorithm="retrostar", max_transforms=4, route_scorer="route cost",
                scorers=["fraction in stock"], filter_cutoff=0.2,
                frozen_substructures=["c1ccccc1"], stock_stop_criteria={"weight": 300}
                )._pipe_config()
    assert cfg["search_algorithm"].endswith("retrostar.search_tree.SearchTree")
    assert cfg["max_transforms"] == 4
    assert cfg["route_scorer"] == "route cost"
    assert cfg["extra_scorers"] == [
        {"name": "fraction in stock", "column": "fraction_in_stock", "class": "FractionInStockScorer"},
        {"name": "route cost", "column": "route_cost", "class": "RouteCostScorer"},
    ]
    assert cfg["filter_cutoff"] == 0.2
    assert cfg["frozen_substructures"] == ["c1ccccc1"]
    assert cfg["stock_stop_criteria"] == {"weight": 300}


def test_custom_scorer_gets_a_column_but_no_class_to_load():
    cfg = _tool(scorers=["my scorer"], config="my.yml")._pipe_config()
    assert cfg["base_config"] == "my.yml"
    assert cfg["extra_scorers"] == [{"name": "my scorer", "column": "my_scorer", "class": None}]


def test_generate_script_writes_the_config_the_pipe_reads(local_config, tmp_path):
    from biopipelines.aizynthfinder import AiZynthFinder
    from biopipelines.datastream import DataStream

    csv = tmp_path / "compounds.csv"
    csv.write_text("id,smiles\nc1,CCO\nc2,CCN\n")
    stream = DataStream(name="compounds", ids=["c1", "c2"], files=[], map_table=str(csv), format="csv")
    tool = AiZynthFinder(compounds=stream, route_scorer="route cost", search_algorithm="retrostar", nproc=2)
    tool.validate_params()
    tool.output_folder = str(tmp_path / "001_AiZynthFinder")
    tool.folders = {"AiZynthFinder": "/data/azf", "pipe_scripts": "/repo/pipe_scripts"}

    script = tool.generate_script("unused.sh")

    with open(tool.config_json, encoding="utf-8") as f:
        written = json.load(f)
    assert written == json.loads(json.dumps(tool._pipe_config()))
    assert written["search_algorithm"].endswith("retrostar.search_tree.SearchTree")
    assert written["extra_scorers"][0]["column"] == "route_cost"
    assert "pipe_aizynthfinder.py" in script and "aizynthfinder_config.json" in script
    assert "--local-missing" in script


def test_display_lists_only_what_was_set():
    assert not any("ALGORITHM" in line for line in _tool().get_config_display())
    lines = _tool(search_algorithm="dfpn", max_transforms=3).get_config_display()
    assert "SEARCH ALGORITHM: dfpn" in lines and "MAX TRANSFORMS: 3" in lines


# ── overlay onto the installed config ────────────────────────────────────────

PUBLIC_CONFIG = {
    "expansion": {"uspto": ["/d/uspto_model.onnx", "/d/uspto_templates.csv.gz"],
                  "ringbreaker": ["/d/rb_model.onnx", "/d/rb_templates.csv.gz"]},
    "filter": {"uspto": "/d/uspto_filter_model.onnx"},
    "stock": {"zinc": "/d/zinc_stock.hdf5"},
}


def _pipe_cfg(**overrides):
    cfg = {"min_routes": 5, "max_routes": 25, "time_limit": None, "iteration_limit": None,
           "search_algorithm": None, "max_transforms": None, "route_scorer": None,
           "extra_scorers": [], "filter_cutoff": None, "frozen_substructures": None,
           "stock_stop_criteria": None}
    cfg.update(overrides)
    return cfg


def test_overlay_without_options_only_adds_route_counts():
    """Backward compatibility: the 2.1 overlay wrote search limits and route counts, nothing else."""
    build = _load_helper().build_search_config
    config = build(copy.deepcopy(PUBLIC_CONFIG), _pipe_cfg())
    assert config["expansion"] == PUBLIC_CONFIG["expansion"]
    assert config["filter"] == PUBLIC_CONFIG["filter"]
    assert config["stock"] == PUBLIC_CONFIG["stock"]
    assert config["search"] == {}
    assert config["post_processing"] == {"min_routes": 5, "max_routes": 25}
    assert "scorer" not in config


def test_overlay_does_not_touch_the_input_config():
    build = _load_helper().build_search_config
    source = copy.deepcopy(PUBLIC_CONFIG)
    build(source, _pipe_cfg(filter_cutoff=0.3, frozen_substructures=["C=O"],
                            stock_stop_criteria={"weight": 200}))
    assert source == PUBLIC_CONFIG


def test_overlay_limits_algorithm_and_depth():
    build = _load_helper().build_search_config
    config = build(copy.deepcopy(PUBLIC_CONFIG), _pipe_cfg(
        time_limit=60, iteration_limit=200, max_transforms=4,
        search_algorithm="aizynthfinder.search.dfpn.search_tree.SearchTree", min_routes=2, max_routes=9))
    assert config["search"] == {
        "time_limit": 60, "iteration_limit": 200, "max_transforms": 4,
        "algorithm": "aizynthfinder.search.dfpn.search_tree.SearchTree"}
    assert config["post_processing"] == {"min_routes": 2, "max_routes": 9}


def test_overlay_keeps_a_search_setting_the_user_config_already_has():
    build = _load_helper().build_search_config
    base = {"search": {"algorithm": "mcts", "algorithm_config": {"C": 1.0}, "max_transforms": 9}}
    config = build(base, _pipe_cfg(time_limit=30))
    assert config["search"] == {"algorithm": "mcts", "algorithm_config": {"C": 1.0},
                                "max_transforms": 9, "time_limit": 30}


def test_overlay_ranks_routes_by_one_scorer_and_loads_the_extras():
    build = _load_helper().build_search_config
    config = build(copy.deepcopy(PUBLIC_CONFIG), _pipe_cfg(
        route_scorer="route cost",
        extra_scorers=[{"name": "route cost", "column": "route_cost", "class": "RouteCostScorer"},
                       {"name": "my scorer", "column": "my_scorer", "class": None}]))
    assert config["post_processing"]["route_scorers"] == ["route cost"]
    assert config["scorer"] == {"RouteCostScorer": None}


def test_overlay_keeps_the_arguments_of_a_scorer_the_user_configured():
    build = _load_helper().build_search_config
    base = dict(PUBLIC_CONFIG, scorer={"RouteCostScorer": {"reaction_cost": 5}})
    config = build(base, _pipe_cfg(extra_scorers=[
        {"name": "route cost", "column": "route_cost", "class": "RouteCostScorer"}]))
    assert config["scorer"] == {"RouteCostScorer": {"reaction_cost": 5}}


def test_overlay_sets_the_cutoff_on_the_quick_filter_shortcut():
    build = _load_helper().build_search_config
    config = build(copy.deepcopy(PUBLIC_CONFIG), _pipe_cfg(filter_cutoff=0.3))
    assert config["filter"] == {"uspto": {"model": "/d/uspto_filter_model.onnx", "filter_cutoff": 0.3}}


def test_overlay_sets_the_cutoff_only_on_quick_filters():
    build = _load_helper().build_search_config
    base = dict(PUBLIC_CONFIG, filter={
        "a": {"type": "quick-filter", "model": "/m1"},
        "b": {"model": "/m2"},
        "c": {"type": "reactants_count"},
        "d": {"type": "frozen_substructure", "smarts_list": ["C"]},
    })
    filters = build(base, _pipe_cfg(filter_cutoff=0.4))["filter"]
    assert filters["a"]["filter_cutoff"] == 0.4 and filters["b"]["filter_cutoff"] == 0.4
    assert "filter_cutoff" not in filters["c"] and "filter_cutoff" not in filters["d"]


def test_overlay_refuses_a_cutoff_when_no_quick_filter_exists():
    build = _load_helper().build_search_config
    with pytest.raises(ValueError, match="no quick-filter"):
        build({"expansion": {}}, _pipe_cfg(filter_cutoff=0.1))


def test_overlay_adds_the_frozen_substructure_filter():
    build = _load_helper().build_search_config
    config = build(copy.deepcopy(PUBLIC_CONFIG), _pipe_cfg(frozen_substructures=["c1ccccc1", "C(=O)N"]))
    assert config["filter"]["uspto"] == PUBLIC_CONFIG["filter"]["uspto"]
    assert config["filter"]["frozen_substructure"] == {
        "type": "frozen_substructure", "smarts_list": ["c1ccccc1", "C(=O)N"]}


def test_overlay_refuses_to_overwrite_a_filter_of_the_same_name():
    build = _load_helper().build_search_config
    base = dict(PUBLIC_CONFIG, filter={"frozen_substructure": {"type": "frozen_substructure", "smarts_list": ["C"]}})
    with pytest.raises(ValueError, match="already defines"):
        build(base, _pipe_cfg(frozen_substructures=["N"]))


def test_overlay_adds_stock_stop_criteria_beside_the_stocks():
    build = _load_helper().build_search_config
    config = build(copy.deepcopy(PUBLIC_CONFIG), _pipe_cfg(stock_stop_criteria={"weight": 250, "counts": {"C": 10}}))
    assert config["stock"] == {"zinc": "/d/zinc_stock.hdf5",
                               "stop_criteria": {"weight": 250, "counts": {"C": 10}}}


def test_overlay_merges_stop_criteria_over_the_ones_in_the_config():
    build = _load_helper().build_search_config
    base = {"stock": {"zinc": "/z", "stop_criteria": {"price": 10, "weight": 500}}}
    config = build(base, _pipe_cfg(stock_stop_criteria={"weight": 250}))
    assert config["stock"]["stop_criteria"] == {"price": 10, "weight": 250}


def test_overlay_survives_a_yaml_round_trip():
    build = _load_helper().build_search_config
    config = build(copy.deepcopy(PUBLIC_CONFIG), _pipe_cfg(
        route_scorer="route cost", filter_cutoff=0.1, frozen_substructures=["C=O"],
        stock_stop_criteria={"weight": 250},
        extra_scorers=[{"name": "route cost", "column": "route_cost", "class": "RouteCostScorer"}]))
    assert yaml.safe_load(yaml.safe_dump(config, sort_keys=False)) == config


def test_filter_selection_adds_the_frozen_filter_to_an_explicit_choice_only():
    selection = _load_helper()._filter_selection
    assert selection(_pipe_cfg(filter_policy=["uspto"], frozen_substructures=["C"])) == ["uspto", "frozen_substructure"]
    assert selection(_pipe_cfg(filter_policy=["uspto", "frozen_substructure"], frozen_substructures=["C"])) == [
        "uspto", "frozen_substructure"]
    assert selection(_pipe_cfg(filter_policy=None, frozen_substructures=["C"])) == []
    assert selection(_pipe_cfg(filter_policy=["uspto"])) == ["uspto"]


def test_invalid_smarts_stops_the_step_before_the_search():
    pytest.importorskip("rdkit")
    validate = _load_helper()._validate_smarts
    validate(["c1ccccc1", "[OH]C=O"])
    with pytest.raises(ValueError, match="not a valid SMARTS"):
        validate(["C(("])


# ── aizynthcli command line ──────────────────────────────────────────────────

def test_cli_command_passes_selections_and_caps_the_processes(monkeypatch, tmp_path):
    helper = _load_helper()
    seen = {}
    monkeypatch.setattr(helper.subprocess, "run", lambda cmd, check: seen.update(cmd=cmd))
    cfg = _pipe_cfg(targets_smi="t.smi", cli_output=str(tmp_path / "o.json.gz"), nproc=16,
                    stocks=["zinc"], expansion_policy=["uspto", "ringbreaker"],
                    filter_policy=["uspto"], frozen_substructures=["C"])
    helper._run_cli(cfg, "search.yml", 3, "")
    cmd = seen["cmd"]
    assert cmd[:7] == ["aizynthcli", "--config", "search.yml", "--smiles", "t.smi", "--output", cfg["cli_output"]]
    assert cmd[cmd.index("--stocks") + 1] == "zinc"
    assert cmd[cmd.index("--policy") + 1:cmd.index("--policy") + 3] == ["uspto", "ringbreaker"]
    assert cmd[cmd.index("--filter") + 1:cmd.index("--filter") + 3] == ["uspto", "frozen_substructure"]
    assert cmd[cmd.index("--nproc") + 1] == "3"


def test_cli_command_leaves_defaults_to_aizynthcli_when_nothing_is_set(monkeypatch, tmp_path):
    helper = _load_helper()
    seen = {}
    monkeypatch.setattr(helper.subprocess, "run", lambda cmd, check: seen.update(cmd=cmd))
    helper._run_cli(_pipe_cfg(targets_smi="t.smi", cli_output=str(tmp_path / "o.json.gz")), "search.yml", 3, "")
    assert not {"--stocks", "--policy", "--filter", "--nproc"} & set(seen["cmd"])


# ── aizynthcli output (4.x format) ───────────────────────────────────────────

def _leaf(smiles, in_stock):
    return {"type": "mol", "hide": False, "smiles": smiles, "is_chemical": True, "in_stock": in_stock}


def _tree(target, leaves, scores):
    """A reaction tree the way ReactionTree._build_dict plus dict_with_extra write it."""
    return {
        "type": "mol", "hide": False, "smiles": target, "is_chemical": True, "in_stock": False,
        "children": [{
            "type": "reaction", "hide": False, "smiles": "A.B>>C", "is_reaction": True,
            "metadata": {"policy_name": "uspto", "template_code": 7},
            "children": [_leaf(smiles, flag) for smiles, flag in leaves],
        }],
        "scores": scores,
        "metadata": {"created_at_iteration": 3},
    }


def _scores(**extra):
    scores = {"state score": 0.97, "number of reactions": 1, "number of pre-cursors": 2,
              "number of pre-cursors in stock": 2}
    scores.update(extra)
    return scores


def _row(target, trees, is_solved=True):
    return {
        "target": target, "search_time": 1.5, "first_solution_time": 0.2,
        "first_solution_iteration": 3, "number_of_nodes": 12, "max_transforms": 2,
        "max_children": 5, "number_of_routes": 7, "number_of_solved_routes": 2,
        "top_score": 0.97, "is_solved": is_solved, "number_of_steps": 1,
        "number_of_precursors": 2, "number_of_precursors_in_stock": 2,
        "precursors_in_stock": "CC, O", "precursors_not_in_stock": "",
        "precursors_availability": "zinc;zinc", "policy_used_counts": {"uspto": 3},
        "profiling": {"expansion_calls": 4, "reactants_generations": 9},
        "stock_info": {"CC": ["zinc"], "O": ["zinc"]}, "trees": trees,
    }


def _write_cli_output(path, rows):
    """What aizynthcli's save_datafile does: DataFrame.to_json(orient='table')."""
    payload = pd.DataFrame(rows).to_json(orient="table")
    with gzip.open(path, "wt") as f:
        f.write(payload)


REAL_OUTPUT = os.path.join(REPO_ROOT, "tests", "fixtures", "aizynthfinder_4_4_1_cli_output.json.gz")


def test_real_aizynthcli_4_4_1_output_parses_end_to_end(tmp_path):
    """Written by aizynthcli 4.4.1 on toy targets (an ester, an amide, a nitrile; a fourth, unsanitizable SMILES was skipped).

    The search used a toy expansion and stock, `route_scorers: [route cost]` and the four extra scorers, so this pins the real row keys, the tree layout and the scores dict.
    """
    pytest.importorskip("rdkit")
    helper = _load_helper()
    rows = helper._read_cli_output(REAL_OUTPUT)
    assert len(rows) == 3
    assert {"index", "target", "is_solved", "top_score", "number_of_routes", "number_of_solved_routes",
            "number_of_nodes", "search_time", "first_solution_time", "first_solution_iteration",
            "trees", "stock_info", "policy_used_counts", "profiling"} <= set(rows[0])

    smiles = ["CCOC(C)=O", "C(C)(C)(C)(C)C", "CCNC(C)=O", "c1ccccc1C#N"]
    ids = ["t1", "t2", "t3", "t4"]
    targets = [{"id": i, "smiles": s, "canonical": helper._canonical_smiles(s)} for i, s in zip(ids, smiles)]
    valid = [t for t in targets if t["canonical"]]
    assert [t["id"] for t in valid] == ["t1", "t3", "t4"]
    assert [row["target"] for row in rows] == [t["canonical"] for t in valid]

    extra = [{"name": "route cost", "column": "route_cost", "class": "RouteCostScorer"},
             {"name": "fraction in stock", "column": "fraction_in_stock", "class": "FractionInStockScorer"}]
    pairs, skipped = helper._align_results(valid, rows)
    assert skipped == []
    converted = {t["id"]: helper._convert_result(t, row, extra, str(tmp_path), "001_AiZynthFinder")
                 for t, row in pairs}

    ester = converted["t1"]
    assert ester["target"]["is_solved"] is True and ester["target"]["top_score"] == 3.5
    best = ester["routes"][0]
    assert (best["number_of_steps"], best["number_of_precursors"], best["number_of_precursors_in_stock"]) == (1, 2, 2)
    assert best["route_cost"] == 3.5 and best["fraction_in_stock"] == 1.0
    assert [p["smiles"] for p in ester["precursors"]] == ["CCO", "CC(=O)O"]

    amide = converted["t3"]
    assert [r["route_cost"] for r in amide["routes"]] == [3.5, 14.75]
    assert [r["is_solved"] for r in amide["routes"]] == [True, False]

    nitrile = converted["t4"]
    assert nitrile["target"]["is_solved"] is False
    assert nitrile["target"]["first_solution_time"] == 0 and nitrile["target"]["first_solution_iteration"] == 0
    assert nitrile["routes"][0]["number_of_steps"] == 0 and nitrile["routes"][0]["is_solved"] is False


def test_cli_output_rows_come_from_the_table_envelope(tmp_path):
    rows = [_row("CCO", [_tree("CCO", [("CC", True), ("O", True)], _scores())])]
    out = str(tmp_path / "output.json.gz")
    _write_cli_output(out, rows)
    parsed = _load_helper()._read_cli_output(out)
    assert len(parsed) == 1
    assert parsed[0]["target"] == "CCO"
    assert parsed[0]["trees"][0]["scores"]["state score"] == 0.97
    assert parsed[0]["trees"][0]["children"][0]["children"][1]["in_stock"] is True


def test_cli_output_with_no_rows_parses_to_an_empty_list(tmp_path):
    out = str(tmp_path / "output.json.gz")
    with gzip.open(out, "wt") as f:
        f.write(pd.DataFrame({"target": []}).to_json(orient="table"))
    assert _load_helper()._read_cli_output(out) == []


def test_convert_result_reads_scores_by_name_and_counts_leaves(tmp_path):
    helper = _load_helper()
    extra = [{"name": "route cost", "column": "route_cost", "class": "RouteCostScorer"},
             {"name": "sum of prices", "column": "sum_of_prices", "class": "PriceSumScorer"}]
    trees = [_tree("CCO", [("CC", True), ("O", False)],
                   _scores(**{"route cost": 3.4, "number of pre-cursors in stock": 1})),
             _tree("CCO", [("CC", True), ("O", True)], _scores(**{"route cost": 2.9}))]
    out = helper._convert_result({"id": "c1", "smiles": "CCO"}, _row("CCO", trees), extra,
                                 str(tmp_path), "001_AiZynthFinder")

    assert out["target"]["number_of_routes"] == 7
    assert out["target"]["number_of_returned_routes"] == 2
    assert out["target"]["is_solved"] is True
    first, second = out["routes"]
    assert first["id"] == "c1_1" and first["route"] == 1
    assert first["score"] == 0.97 and first["number_of_steps"] == 1
    assert first["is_solved"] is False and second["is_solved"] is True
    assert first["number_of_precursors_in_stock"] == 1
    assert first["route_cost"] == 3.4 and second["route_cost"] == 2.9
    assert first["sum_of_prices"] is None
    assert [p["id"] for p in out["precursors"]] == ["c1_1_1", "c1_1_2", "c1_2_1", "c1_2_2"]
    assert [p["in_stock"] for p in out["precursors"]] == [True, False, True, True]
    assert (tmp_path / "c1_1.json").is_file() and (tmp_path / "c1_2.json").is_file()
    assert out["missing"] == []


def test_convert_result_without_trees_reports_a_filter(tmp_path):
    out = _load_helper()._convert_result(
        {"id": "c1", "smiles": "CCO"}, _row("CCO", [], is_solved=False), [], str(tmp_path), "001_AiZynthFinder")
    assert out["target"]["number_of_returned_routes"] == 0
    assert out["routes"] == [] and out["precursors"] == []
    assert out["missing"] == [{"id": "c1", "removed_by": "001_AiZynthFinder", "kind": "filter",
                               "cause": "no retrosynthetic route found"}]


def test_convert_result_takes_counts_from_the_tree_when_scores_lack_them(tmp_path):
    tree = _tree("CCO", [("CC", True), ("O", True)], {"state score": 0.5})
    out = _load_helper()._convert_result({"id": "c1", "smiles": "CCO"}, _row("CCO", [tree]), [],
                                         str(tmp_path), "s")
    route = out["routes"][0]
    assert route["number_of_steps"] == 1 and route["number_of_precursors"] == 2
    assert route["number_of_precursors_in_stock"] == 2


# ── aligning rows to targets ─────────────────────────────────────────────────

def _targets(*canonical):
    return [{"id": f"c{i}", "smiles": s, "canonical": s} for i, s in enumerate(canonical, start=1)]


def test_align_pairs_rows_to_targets_by_order_when_counts_agree():
    align = _load_helper()._align_results
    targets = _targets("CC", "CCO")
    pairs, skipped = align(targets, [{"target": "CC"}, {"target": "CCO"}])
    assert [t["id"] for t, _ in pairs] == ["c1", "c2"] and skipped == []


def test_align_names_the_target_aizynthcli_skipped():
    align = _load_helper()._align_results
    targets = _targets("CC", "CCO", "CCN")
    pairs, skipped = align(targets, [{"target": "CC"}, {"target": "CCN"}])
    assert [t["id"] for t, _ in pairs] == ["c1", "c3"]
    assert [t["id"] for t in skipped] == ["c2"]


def test_align_handles_a_skipped_first_and_last_target():
    align = _load_helper()._align_results
    targets = _targets("CC", "CCO", "CCN")
    pairs, skipped = align(targets, [{"target": "CCO"}])
    assert [t["id"] for t, _ in pairs] == ["c2"]
    assert [t["id"] for t in skipped] == ["c1", "c3"]


def test_align_with_duplicate_targets_keeps_order():
    align = _load_helper()._align_results
    targets = _targets("CC", "CC", "CCO")
    pairs, skipped = align(targets, [{"target": "CC", "n": 1}, {"target": "CCO"}])
    assert [t["id"] for t, _ in pairs] == ["c1", "c3"] and [t["id"] for t in skipped] == ["c2"]


def test_align_refuses_rows_that_match_no_target():
    align = _load_helper()._align_results
    with pytest.raises(RuntimeError, match="cannot align"):
        align(_targets("CC", "CCO", "CCN"), [{"target": "CC"}, {"target": "OTHER"}])
    with pytest.raises(RuntimeError, match="cannot align"):
        align(_targets("CC"), [{"target": "CC"}, {"target": "CC"}])


def test_canonical_smiles_matches_what_upstream_reports_as_target():
    pytest.importorskip("rdkit")
    canonical = _load_helper()._canonical_smiles
    assert canonical("OCC") == "CCO"
    assert canonical("C1=CC=CC=C1") == "c1ccccc1"
    # Explicit hydrogens survive: upstream parses unsanitized, so none are removed.
    assert canonical("[H]OC") == "[H]OC"


def test_canonical_smiles_is_none_for_what_upstream_cannot_sanitize():
    pytest.importorskip("rdkit")
    canonical = _load_helper()._canonical_smiles
    assert canonical("not a smiles") is None
    assert canonical("C(C)(C)(C)(C)C") is None
    assert canonical("n1(c(=O)[nH]c2c(c1=O)c(c(s2)C(=O)N(C)C)C)c1c(N(=O)O)cccc1") is None


# ── the whole helper on a synthetic aizynthcli ───────────────────────────────

def _run_helper(tmp_path, monkeypatch, smiles, cli_rows, **cfg_extra):
    """Run main() with aizynthcli replaced by a function that writes `cli_rows(targets)`."""
    from biopipelines.datastream import DataStream
    helper = _load_helper()
    step = tmp_path / "001_AiZynthFinder"
    (step / "_configuration").mkdir(parents=True)
    ids = [f"c{i}" for i in range(1, len(smiles) + 1)]
    map_table = step / "compounds.csv"
    pd.DataFrame({"id": ids, "smiles": smiles}).to_csv(map_table, index=False)
    DataStream(name="compounds", ids=ids, files=[], map_table=str(map_table), format="csv"
               ).save_json(str(step / "_configuration" / "compounds_ds.json"))
    base = step / "config.yml"
    base.write_text(yaml.safe_dump(PUBLIC_CONFIG))

    cfg = _pipe_cfg(
        compounds_json=str(step / "_configuration" / "compounds_ds.json"),
        base_config=str(base),
        search_config_yaml=str(step / "_configuration" / "search_config.yml"),
        targets_smi=str(step / "_configuration" / "targets.smi"),
        cli_output=str(step / "_execution" / "aizynth_output.json.gz"),
        routes_dir=str(step / "routes"),
        routes_map_csv=str(step / "routes" / "routes_map.csv"),
        precursors_csv=str(step / "precursors" / "precursors.csv"),
        targets_csv=str(step / "tables" / "targets.csv"),
        routes_csv=str(step / "tables" / "routes.csv"),
        local_missing_csv=str(step / "_execution" / "local_missing.csv"),
        stocks=None, expansion_policy=None, filter_policy=None, nproc=None, **cfg_extra)
    cfg_path = step / "_configuration" / "aizynthfinder_config.json"
    cfg_path.write_text(json.dumps(cfg))

    def fake_cli(cfg_, search_config, n_targets, prefix):
        with open(cfg_["targets_smi"]) as f:
            sent = [line.strip() for line in f]
        assert len(sent) == n_targets
        os.makedirs(os.path.dirname(cfg_["cli_output"]), exist_ok=True)
        _write_cli_output(cfg_["cli_output"], cli_rows(sent))

    monkeypatch.setattr(helper, "_run_cli", fake_cli)
    monkeypatch.setattr(sys, "argv", ["pipe_aizynthfinder.py", "--config", str(cfg_path)])
    helper.main()
    read = lambda name: pd.read_csv(step / name)
    return step, read


def _good(target):
    return _row(target, [_tree(target, [("CC", True), ("O", True)], _scores(**{"route cost": 2.5}))])


def test_main_unsanitizable_target_no_longer_breaks_the_row_alignment(tmp_path, monkeypatch):
    """aizynthcli drops a target it cannot set up, which used to abort the step on the row count."""
    pytest.importorskip("rdkit")
    extra = [{"name": "route cost", "column": "route_cost", "class": "RouteCostScorer"}]
    step, read = _run_helper(
        tmp_path, monkeypatch, ["OCC", "C(C)(C)(C)(C)C", "CCN"],
        lambda sent: [_good("CCO"), _good("CCN")], extra_scorers=extra, route_scorer="route cost")

    assert (step / "_configuration" / "targets.smi").read_text().split() == ["OCC", "CCN"]
    targets = read("tables/targets.csv")
    assert list(targets["id"]) == ["c1", "c2", "c3"]
    assert targets.loc[targets["id"] == "c2", "smiles"].item() == "C(C)(C)(C)(C)C"
    assert pd.isna(targets.loc[targets["id"] == "c2", "number_of_routes"].item())
    assert targets.loc[targets["id"] == "c3", "number_of_routes"].item() == 7

    routes = read("tables/routes.csv")
    assert list(routes["compounds.id"]) == ["c1", "c3"]
    assert list(routes["route_cost"]) == [2.5, 2.5]
    assert set(read("precursors/precursors.csv")["compounds.id"]) == {"c1", "c3"}

    missing = read("_execution/local_missing.csv")
    assert list(missing["id"]) == ["c2"] and list(missing["kind"]) == ["failure"]

    written = yaml.safe_load((step / "_configuration" / "search_config.yml").read_text())
    assert written["post_processing"]["route_scorers"] == ["route cost"]
    assert written["scorer"] == {"RouteCostScorer": None}


def test_main_reports_a_target_the_cli_skipped_without_guessing(tmp_path, monkeypatch):
    pytest.importorskip("rdkit")
    step, read = _run_helper(
        tmp_path, monkeypatch, ["CCO", "CCN", "CCC"],
        lambda sent: [_good("CCO"), _good("CCC")])
    targets = read("tables/targets.csv")
    assert list(targets["id"]) == ["c1", "c2", "c3"]
    assert pd.isna(targets.loc[targets["id"] == "c2", "is_solved"].item())
    missing = read("_execution/local_missing.csv")
    assert list(missing["id"]) == ["c2"] and list(missing["kind"]) == ["failure"]
    assert list(read("tables/routes.csv")["compounds.id"]) == ["c1", "c3"]


def test_main_default_run_keeps_the_2_1_tables(tmp_path, monkeypatch):
    pytest.importorskip("rdkit")
    step, read = _run_helper(tmp_path, monkeypatch, ["CCO", "CCN"],
                             lambda sent: [_good("CCO"), _row("CCN", [], is_solved=False)])
    from biopipelines.aizynthfinder import ROUTE_COLUMNS, TARGET_COLUMNS, PRECURSOR_COLUMNS
    assert list(read("tables/targets.csv").columns) == TARGET_COLUMNS
    assert list(read("tables/routes.csv").columns) == ROUTE_COLUMNS
    assert list(read("precursors/precursors.csv").columns) == PRECURSOR_COLUMNS
    targets = read("tables/targets.csv")
    assert list(targets["is_solved"]) == [True, False]
    missing = read("_execution/local_missing.csv")
    assert list(missing["cause"]) == ["no retrosynthetic route found"]
    assert list(missing["kind"]) == ["filter"]


def test_main_exits_non_zero_when_no_target_can_be_searched(tmp_path, monkeypatch):
    pytest.importorskip("rdkit")
    with pytest.raises(SystemExit) as exit_info:
        _run_helper(tmp_path, monkeypatch, ["C(C)(C)(C)(C)C", "xyz"], lambda sent: [])
    assert exit_info.value.code == 1


def test_main_refuses_an_invalid_frozen_pattern_before_searching(tmp_path, monkeypatch):
    pytest.importorskip("rdkit")
    with pytest.raises(ValueError, match="not a valid SMARTS"):
        _run_helper(tmp_path, monkeypatch, ["CCO"], lambda sent: [_good("CCO")],
                    frozen_substructures=["C(("])
