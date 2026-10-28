"""RFD4Proteina: contig and condition handling, the spec it writes, ids and streams, the generated script.

The model itself needs a GPU and a CUDA 13 environment, so nothing here runs it. What can be checked without it is everything the wrapper decides: that the spec matches what upstream documents, that a contig upstream would reject is refused at construction, that the ids the wrapper declares are the ids the post-processing writes, and that every value reaching `rfd4-generate` is spelled the way its CLI documents.
"""

import importlib.util
import json
import pathlib
import subprocess
import sys

import pytest

yaml = pytest.importorskip("yaml")

from biopipelines.rfd4proteina import (
    RFD4Proteina, build_conditions, hotspot_selects, merge_conditions, normalize_conditions,
    parse_contig, recycle_strategy,
)

ROOT = pathlib.Path(__file__).resolve().parent.parent
PIPE_SCRIPTS = ROOT / "pipe_scripts"

# docs/how-to/input_prep/examples/condition_spec_ppi.yaml of RFD4-Proteina (BSD-3-Clause, Copyright (c) 2026 Institute for Protein Design, University of Washington).
UPSTREAM_PPI_EXAMPLE = """
pdl1_binder:
  input: ../../../examples/inputs/5o45_cropped.cif.zst
  contig: "50-120, /0, A17-131"
  conditions:
    C_CRDGRP:
      1:
        - select: "sel('A/*/*/[N,CA,C,O]') and not has_nan_coord()"
    C_SEQ:
      true:
        - select: "sel('A')"
    C_HOT:
      true:
        - select: "sel('A/*/56/[CG,OH]')"
        - select: "sel('A/*/115/[CG,SD]')"
        - select: "sel('A/*/123/[CD2,OH]')"
"""

PD_L1_HOTSPOTS = "A56:CG,OH;A115:CG,SD;A123:CD2,OH"


def _load_script(name):
    spec = importlib.util.spec_from_file_location(name, PIPE_SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _structures(ids=("s1",)):
    from biopipelines.mock import Mock
    return Mock(ids=list(ids), streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
                map_table_strategy="config")


def _build(local_config, isolated_cwd, new_pipeline, name="rfd4", ids=("s1",), **kwargs):
    """The tool's output and its step script, built inside a real pipeline."""
    if "contigs" not in kwargs:
        kwargs.setdefault("contig", "50-120, /0, A17-131")
    kwargs.setdefault("checkpoint", "/w/rfd4")
    pipeline = new_pipeline(name)
    with pipeline:
        built = RFD4Proteina(structures=_structures(ids), **kwargs)
        script_path = pipeline.save()
    runtime = pathlib.Path(script_path).parent
    matches = sorted(runtime.glob("*_RFD4Proteina.sh"))
    assert matches, f"no RFD4Proteina step script under {runtime}"
    return built, matches[0].read_text(encoding="utf-8", errors="replace")


# ── the contig ────────────────────────────────────────────────────────────────

def test_the_upstream_binder_contig_parses():
    tokens = parse_contig("50-120, /0, A17-131")
    assert [t.kind for t in tokens] == ["expandable", "chainbreak", "preserved"]
    assert (tokens[0].seg_min, tokens[0].seg_max) == (50, 120)
    assert (tokens[2].chain, tokens[2].start, tokens[2].end) == ("A", 17, 131)


def test_named_and_typed_segments_and_multiletter_chains_parse():
    tokens = parse_contig("A42-48, E1:2-20, A52, E2:10R, /0, 20P, AB7-9")
    kinds = [t.kind for t in tokens]
    assert kinds.count("preserved") == 3 and kinds.count("expandable") == 3
    assert [t.chain for t in tokens if t.kind == "preserved"] == ["A", "A", "AB"]
    assert [t.name for t in tokens if t.kind == "expandable"] == ["E1", "E2", ""]


@pytest.mark.parametrize("contig, message", [
    ("", "contig is required"),
    (None, "contig is required"),
    ("50-120, /0, A17-131,", "empty contig token"),
    ("50-120, /0, A17-", "malformed preserved segment"),
    ("50-120, /0, 17A", "malformed expandable segment"),
    ("50-120, /0, A131-17", "end < start"),
    ("E1:20-10, A1-5", "end < start"),
    ("50-120", "no preserved segment"),
    ("A10-20, A15-30", "overlap"),
    ("A10-20, B1-5", "mixes preserved segments"),
    ("A1-5, /0, 10, /0, A9-12", "both named 'A'"),
    ("A1-5, E1:3, E1:4", "duplicate expandable-segment"),
    ("MASKED:3, A1-5", "MASKED"),
    ("E1:, A1-5", "malformed expandable segment"),
])
def test_a_contig_upstream_would_reject_is_refused(contig, message):
    with pytest.raises(ValueError, match=message):
        parse_contig(contig)


def test_the_wrapper_refuses_a_bad_contig_at_construction(local_config, isolated_cwd, new_pipeline):
    with new_pipeline("badcontig"):
        with pytest.raises(ValueError, match="no preserved segment"):
            RFD4Proteina(structures=_structures(), contig="50-120", checkpoint="/w/rfd4")


def test_contigs_is_a_synonym_for_contig(local_config, isolated_cwd, new_pipeline):
    """RFdiffusion, RFdiffusion2 and RFdiffusionAllAtom spell it `contigs`; forwarded to hydra it would be silently dropped."""
    built, script = _build(local_config, isolated_cwd, new_pipeline, name="alias",
                           contigs="50-120, /0, A17-131")
    assert "contigs=" not in script
    assert yaml.safe_load(open(built.producer.config_yaml))["spec"]["contig"] == "50-120, /0, A17-131"


# ── hotspots ──────────────────────────────────────────────────────────────────

def test_hotspots_become_atom_selections():
    assert hotspot_selects(PD_L1_HOTSPOTS) == [
        "sel('A/*/56/[CG,OH]')", "sel('A/*/115/[CG,SD]')", "sel('A/*/123/[CD2,OH]')"]


def test_whole_residue_and_single_atom_hotspots():
    assert hotspot_selects("A67,A89") == ["sel('A/*/67')", "sel('A/*/89')"]
    assert hotspot_selects({"A67": "NE", "B5": ""}) == ["sel('A/*/67/NE')", "sel('B/*/5')"]


def test_no_hotspots_means_no_selections():
    assert hotspot_selects(None) == [] and hotspot_selects("") == [] and hotspot_selects({}) == []


@pytest.mark.parametrize("bad", ["67", "A:CA", "A67:C A", "A67:O5'", 5, ["A67"]])
def test_a_malformed_hotspot_is_refused(bad):
    with pytest.raises(ValueError):
        hotspot_selects(bad)


# ── conditions, compared with what upstream documents ────────────────────────

def test_the_binder_recipe_equals_upstreams_ppi_example():
    expected = yaml.safe_load(UPSTREAM_PPI_EXAMPLE)["pdl1_binder"]
    tokens = parse_contig(expected["contig"])
    assert build_conditions("protein_binder", tokens, PD_L1_HOTSPOTS) == expected["conditions"]


def test_a_multichain_target_is_one_internal_coordinate_group():
    tokens = parse_contig("40-80, /0, A1-100, /0, B5-60")
    conditions = build_conditions("protein_binder", tokens)
    assert list(conditions["C_CRDGRP"]) == [1]
    assert [e["select"] for e in conditions["C_CRDGRP"][1]] == [
        "sel('A/*/*/[N,CA,C,O]') and not has_nan_coord()", "sel('B/*/*/[N,CA,C,O]') and not has_nan_coord()"]
    assert [e["select"] for e in conditions["C_SEQ"][True]] == ["sel('A')", "sel('B')"]


def test_a_small_molecule_binder_adds_no_conditions_of_its_own():
    """Upstream sequence-conditions atomized tokens itself and the SAM benchmark leaves the pose free."""
    tokens = parse_contig("150, /0, B386")
    assert build_conditions("small_molecule_binder", tokens) == {}


def test_motif_scaffolding_fixes_every_preserved_segment():
    tokens = parse_contig("A10-25, 20-40, A60, 20-40")
    conditions = build_conditions("motif_scaffold", tokens)
    expected = [{"select": "sel('A/*/[10-25]')"}, {"select": "sel('A/*/60')"}]
    assert conditions == {"C_CRD": {True: expected}, "C_SEQ": {True: expected}}


def test_user_conditions_come_after_the_recipe():
    """Upstream resolves overlaps last write wins, in declaration order."""
    tokens = parse_contig("50-120, /0, A17-131")
    extra = {"C_SEQ": {False: [{"select": "sel('A/*/[20-30]')"}]},
             "C_RAS": {"BURIED": [{"select": "sel('B')"}]}}
    conditions = build_conditions("protein_binder", tokens, None, extra)
    assert list(conditions["C_SEQ"]) == [True, False]
    assert conditions["C_RAS"] == {"BURIED": [{"select": "sel('B')"}]}


def test_hotspots_sit_between_the_recipe_and_the_user_conditions():
    tokens = parse_contig("50-120, /0, A17-131")
    extra = {"C_HOT": {True: [{"select": "sel('A/*/99')"}]}}
    conditions = build_conditions("protein_binder", tokens, "A56", extra)
    assert [e["select"] for e in conditions["C_HOT"][True]] == ["sel('A/*/56')", "sel('A/*/99')"]


@pytest.mark.parametrize("bad", [
    [],
    {"C_HOT": []},
    {"C_HOT": {True: {"select": "sel('A')"}}},
    {"C_HOT": {True: []}},
    {"C_HOT": {True: [{"select": ""}]}},
    {"C_HOT": {True: [{"selection": "sel('A')"}]}},
    {"C_HOT": {True: [{"select": "sel('A')", "color": "red"}]}},
    {"C_HOT": {(1, 2): [{"select": "sel('A')"}]}},
    {5: {True: [{"select": "sel('A')"}]}},
])
def test_a_malformed_conditions_mapping_is_refused(bad):
    with pytest.raises(ValueError):
        normalize_conditions(bad)


def test_a_group_label_is_kept_for_index_conditions():
    given = {"C_IDX": {False: [{"select": "sel('A/*/[10-12]')", "group": "linker_a"}]}}
    assert normalize_conditions(given) == given


def test_merge_does_not_touch_its_inputs():
    base = {"C_SEQ": {True: [{"select": "sel('A')"}]}}
    extra = {"C_SEQ": {True: [{"select": "sel('B')"}]}}
    merged = merge_conditions(base, extra)
    assert len(merged["C_SEQ"][True]) == 2 and len(base["C_SEQ"][True]) == 1 and len(extra["C_SEQ"][True]) == 1


# ── thinking budget ───────────────────────────────────────────────────────────

@pytest.mark.parametrize("task, thinking, expected", [
    ("protein_binder", "low", "low/front_loaded"),
    ("protein_binder", "medium", "medium/front_loaded"),
    ("protein_binder", "high", "high/front_loaded"),
    ("small_molecule_binder", "medium", "medium/midpeak"),
    ("motif_scaffold", "high", "high/midpeak"),
])
def test_the_recycle_strategy_follows_task_and_budget(task, thinking, expected):
    assert recycle_strategy(task, thinking) == expected


# ── construction ──────────────────────────────────────────────────────────────

@pytest.mark.parametrize("kwargs, message", [
    ({"task": "dna_binder"}, "task must be one of"),
    ({"thinking": "extreme"}, "thinking must be one of"),
    ({"num_designs": 0}, "num_designs"),
    ({"num_designs": True}, "num_designs"),
    ({"batch_size": 0}, "batch_size"),
    ({"seed": -1}, "seed"),
    ({"seed": "7"}, "seed"),
    ({"metrics_tags": "design"}, "list of tag strings"),
    ({"metrics_tags": ["design;rm"]}, "must be a string of letters"),
    ({"metrics_tags": ["ppi_refolding"]}, "refolding oracle"),
    ({"metrics_tags": ["antibody"]}, "refolding oracle"),
    ({"checkpoint": '/w/"x'}, "checkpoint"),
    ({"inputs": "other.yaml"}, "set by RFD4Proteina itself"),
    ({"compress_outputs": True}, "set by RFD4Proteina itself"),
    ({"+hydra.run.dir": "x"}, "set by RFD4Proteina itself"),
    ({"conditions": {"C_HOT": []}}, "conditions"),
    ({"hotspots": "67"}, "hotspot"),
])
def test_bad_parameters_are_refused_at_construction(local_config, isolated_cwd, new_pipeline, kwargs, message):
    with new_pipeline("bad"):
        with pytest.raises(ValueError, match=message):
            RFD4Proteina(**{"structures": _structures(), "contig": "50-120, /0, A17-131",
                            "checkpoint": "/w/rfd4", **kwargs})


def test_an_empty_structures_stream_is_refused(local_config, isolated_cwd, new_pipeline):
    from biopipelines.datastream import DataStream
    empty = DataStream(name="structures", ids=[], files=[], map_table="", format="pdb")
    with new_pipeline("empty"):
        _structures()
        with pytest.raises(ValueError, match="structures"):
            RFD4Proteina(structures=empty, contig="50-120, /0, A17-131", checkpoint="/w/rfd4")


def test_a_non_stream_input_is_refused(local_config, isolated_cwd, new_pipeline):
    with new_pipeline("raw"):
        _structures()
        with pytest.raises(ValueError, match="DataStream or StandardizedOutput"):
            RFD4Proteina(structures="7KDL.pdb", contig="50-120, /0, A17-131")


# ── ids and streams ───────────────────────────────────────────────────────────

def test_the_declared_ids_count_designs_from_one(local_config, isolated_cwd, new_pipeline):
    built, _ = _build(local_config, isolated_cwd, new_pipeline, name="ids", num_designs=3, ids=("t1", "t2"))
    structures = built.streams.structures
    assert list(structures.ids) == ["t1_<1..3>", "t2_<1..3>"]
    assert structures.format == "pdb"
    assert set(structures.ids_expanded) == {"t1_1", "t1_2", "t1_3", "t2_1", "t2_2", "t2_3"}
    assert pathlib.Path(structures.map_table).name == "structures_map.csv"


def test_a_multichain_contig_declares_the_lazy_chain_suffix(local_config, isolated_cwd, new_pipeline):
    """The sequences stream is per chain: a bare `s1_1` selector would match nothing against rows `s1_1_A`, `s1_1_B`."""
    built, _ = _build(local_config, isolated_cwd, new_pipeline, name="chains")
    sequences = built.streams.sequences
    assert list(sequences.ids) == ["s1_<1..1>[_<?>]"]
    assert sequences.format == "csv" and list(sequences.files) == []


def test_a_single_chain_contig_declares_plain_ids(local_config, isolated_cwd, new_pipeline):
    built, _ = _build(local_config, isolated_cwd, new_pipeline, name="onechain",
                      contig="A10-25, 20-40, A60-70", task="motif_scaffold", num_designs=2)
    assert list(built.streams.sequences.ids) == ["s1_<1..2>"]


def test_the_tables_are_declared(local_config, isolated_cwd, new_pipeline):
    built, _ = _build(local_config, isolated_cwd, new_pipeline, name="tables")
    tables = built.tables
    assert tables.structures.info.columns == ["id", "source_id", "design", "cif", "pdb", "n_chains", "n_residues"]
    assert tables.sequences.info.columns == ["id", "source_id", "source_pdb", "chain", "sequence", "length"]
    assert tables.metrics.info.columns == ["id", "source_id"]
    assert tables.missing.info.columns == ["id", "removed_by", "kind", "cause"]


# ── the spec the wrapper writes ──────────────────────────────────────────────

def test_the_config_holds_the_binder_entry_upstream_documents(local_config, isolated_cwd, new_pipeline):
    built, _ = _build(local_config, isolated_cwd, new_pipeline, name="cfg", hotspots=PD_L1_HOTSPOTS, num_designs=4)
    config = yaml.safe_load(open(built.producer.config_yaml))
    expected = yaml.safe_load(UPSTREAM_PPI_EXAMPLE)["pdl1_binder"]
    assert config["num_designs"] == 4
    assert config["spec"] == {"contig": expected["contig"], "conditions": expected["conditions"]}


def test_a_spec_without_conditions_carries_only_the_contig(local_config, isolated_cwd, new_pipeline):
    built, _ = _build(local_config, isolated_cwd, new_pipeline, name="nocond",
                      contig="150, /0, B386", task="small_molecule_binder")
    assert yaml.safe_load(open(built.producer.config_yaml))["spec"] == {"contig": "150, /0, B386"}


def test_the_config_round_trips_boolean_and_integer_keys(local_config, isolated_cwd, new_pipeline):
    """`true:` and `1:` must come back as a bool and an int; quoting them would change what upstream reads."""
    built, _ = _build(local_config, isolated_cwd, new_pipeline, name="keys")
    text = open(built.producer.config_yaml).read()
    assert "\n      true:" in text and "\n      1:" in text
    spec = yaml.safe_load(text)["spec"]
    assert list(spec["conditions"]["C_SEQ"]) == [True] and list(spec["conditions"]["C_CRDGRP"]) == [1]


def test_the_builder_writes_one_entry_per_input_in_upstreams_key_order():
    inputs = _load_script("pipe_rfd4proteina_inputs")
    template = {"contig": "50-120, /0, A17-131", "conditions": {"C_SEQ": {True: [{"select": "sel('A')"}]}}}
    spec = inputs.build_spec(template, {"t1": "/x/t1.cif", "t2": "/x/t2.cif"})
    assert list(spec) == ["t1", "t2"]
    assert list(spec["t1"]) == ["input", "contig", "conditions"]
    assert spec["t2"]["input"] == "/x/t2.cif" and spec["t2"]["contig"] == template["contig"]


def test_numeric_looking_ids_stay_strings_in_the_written_spec(tmp_path):
    """An id like 1234 or 1E58 would be read back as a number if it were written unquoted."""
    inputs = _load_script("pipe_rfd4proteina_inputs")
    path = tmp_path / "spec.yaml"
    with open(path, "w") as f:
        yaml.safe_dump(inputs.build_spec({"contig": "A1-5"}, {"1234": "/x/a.cif", "1E58": "/x/b.cif"}),
                       f, sort_keys=False, width=1000)
    assert list(yaml.safe_load(path.read_text())) == ["1234", "1E58"]


@pytest.mark.parametrize("path, expected", [
    ("/x/a.cif", True), ("/x/a.CIF", True), ("/x/a.cif.gz", True), ("/x/a.bcif", True),
    ("/x/a.cif.zst", True), ("/x/a.pdb", False), ("/x/a.ent", False), ("/x/cif", False),
])
def test_only_cif_like_files_skip_the_conversion(path, expected):
    assert _load_script("pipe_rfd4proteina_inputs").is_cif_like(path) is expected


@pytest.mark.parametrize("design_id, valid", [
    ("7KDL", True), ("t1+lig2", True), ("a_b", True), ("", False), (".", False), ("..", False),
    ("a/b", False), ("a\\b", False), ("a\0b", False),
])
def test_ids_that_cannot_name_a_directory_are_refused(design_id, valid):
    assert _load_script("pipe_rfd4proteina_inputs").valid_design_id(design_id) is valid


def _datastream_json(tmp_path, files, ids):
    path = tmp_path / "structures.json"
    path.write_text(json.dumps({"name": "structures", "format": "pdb", "ids": ids, "files": files}))
    return path


def _run_inputs(tmp_path, ids, files, template=None):
    config = tmp_path / "config.yaml"
    config.write_text(yaml.safe_dump({"num_designs": 2, "spec": template or {"contig": "50-120, /0, A17-131"}}))
    return subprocess.run(
        [sys.executable, str(PIPE_SCRIPTS / "pipe_rfd4proteina_inputs.py"),
         "--config-yaml", str(config), "--structures-json", str(_datastream_json(tmp_path, files, ids)),
         "--inputs-dir", str(tmp_path / "inputs"), "--spec-yaml", str(tmp_path / "spec.yaml"),
         "--manifest-json", str(tmp_path / "manifest.json")],
        capture_output=True, text=True)


def test_cif_inputs_are_referenced_in_place(tmp_path):
    for sid in ("t1", "t2"):
        (tmp_path / f"{sid}.cif").write_text("data_x\n")
    result = _run_inputs(tmp_path, ["t1", "t2"], [str(tmp_path / "t1.cif"), str(tmp_path / "t2.cif")])
    assert result.returncode == 0, result.stderr
    spec = yaml.safe_load((tmp_path / "spec.yaml").read_text())
    assert spec["t1"]["input"] == str(tmp_path / "t1.cif") and spec["t2"]["contig"] == "50-120, /0, A17-131"
    manifest = json.loads((tmp_path / "manifest.json").read_text())
    assert manifest == {"attempted": ["t1", "t2"], "built": ["t1", "t2"], "failed": {}}


@pytest.mark.skipif(importlib.util.find_spec("atomworks") is not None,
                    reason="AtomWorks is installed here, so a PDB would be converted rather than fail")
def test_a_failed_conversion_is_recorded_and_the_others_still_run(tmp_path):
    (tmp_path / "ok.cif").write_text("data_x\n")
    (tmp_path / "bad.pdb").write_text("ATOM      1  CA  ALA A   1\n")
    result = _run_inputs(tmp_path, ["ok", "bad"], [str(tmp_path / "ok.cif"), str(tmp_path / "bad.pdb")])
    assert result.returncode == 0
    manifest = json.loads((tmp_path / "manifest.json").read_text())
    assert manifest["built"] == ["ok"] and list(manifest["failed"]) == ["bad"]
    assert list(yaml.safe_load((tmp_path / "spec.yaml").read_text())) == ["ok"]


@pytest.mark.skipif(importlib.util.find_spec("atomworks") is not None,
                    reason="AtomWorks is installed here, so a PDB would be converted rather than fail")
def test_every_input_failing_is_a_step_failure(tmp_path):
    (tmp_path / "bad.pdb").write_text("ATOM      1  CA  ALA A   1\n")
    result = _run_inputs(tmp_path, ["bad"], [str(tmp_path / "bad.pdb")])
    assert result.returncode == 1 and "all 1 input(s) failed" in result.stderr
    assert not (tmp_path / "spec.yaml").exists()
    assert json.loads((tmp_path / "manifest.json").read_text())["built"] == []


def test_a_pdb_is_parsed_without_the_annotations_atomworks_refuses_for_it(monkeypatch, tmp_path):
    """AtomWorks raises on a PDB parsed with standard annotations, which upstream's DISTILLATION preset requests."""
    import types

    class Preset(dict):
        def replace(self, **changes):
            return Preset({**self, **changes})

    def parse(path, config):
        if config["load_standard_annotations"]:
            raise ValueError("load_standard_annotations=True is not supported for PDB files. Use CIF format instead.")
        return {"asym_unit": ["atoms"]}

    written = {}
    modules = {
        "atomworks": types.ModuleType("atomworks"), "atomworks.io": types.ModuleType("atomworks.io"),
        "atomworks.io.utils": types.ModuleType("atomworks.io.utils"),
        "atomworks.io.utils.io_utils": types.ModuleType("atomworks.io.utils.io_utils"),
        "rfproteina": types.ModuleType("rfproteina"), "rfproteina.utils": types.ModuleType("rfproteina.utils"),
        "rfproteina.utils.parse_presets": types.ModuleType("rfproteina.utils.parse_presets"),
    }
    modules["atomworks.io"].parse = parse
    modules["atomworks.io.utils.io_utils"].to_cif_file = lambda atoms, path: written.update(path=path, atoms=atoms)
    modules["rfproteina.utils.parse_presets"].DISTILLATION = Preset(build_assembly="all", load_standard_annotations=True)
    for name, module in modules.items():
        monkeypatch.setitem(sys.modules, name, module)

    _load_script("pipe_rfd4proteina_inputs").convert_to_cif("in.pdb", str(tmp_path / "out.cif"))
    assert written == {"path": str(tmp_path / "out.cif"), "atoms": "atoms"}


def test_an_input_the_contig_does_not_fit_fails_alone_with_upstreams_message(monkeypatch, tmp_path):
    """rfd4-generate stops at the first example that raises, so a contig one input cannot hold must not reach it."""
    import types

    message = "Preserved segment 'A1-40' runs outside chain A, whose residues are 1-20."

    def build(spec, *, spec_dir, example_id):
        if example_id == "short":
            raise ValueError(message)

    for name in ("rfproteina", "rfproteina.inference_utils", "rfproteina.inference_utils.condition_spec"):
        monkeypatch.setitem(sys.modules, name, types.ModuleType(name))
    sys.modules["rfproteina.inference_utils.condition_spec"].build_annotated_atom_array_from_condition_spec = build

    for sid in ("short", "ubq"):
        (tmp_path / f"{sid}.cif").write_text("data_x\n")
    (tmp_path / "config.yaml").write_text(yaml.safe_dump({"num_designs": 1, "spec": {"contig": "50-70, /0, A1-40"}}))
    structures = _datastream_json(tmp_path, [str(tmp_path / "short.cif"), str(tmp_path / "ubq.cif")], ["short", "ubq"])
    monkeypatch.setattr(sys, "argv", [
        "inputs", "--config-yaml", str(tmp_path / "config.yaml"), "--structures-json", str(structures),
        "--inputs-dir", str(tmp_path / "inputs"), "--spec-yaml", str(tmp_path / "spec.yaml"),
        "--manifest-json", str(tmp_path / "manifest.json")])
    _load_script("pipe_rfd4proteina_inputs").main()

    manifest = json.loads((tmp_path / "manifest.json").read_text())
    assert manifest["built"] == ["ubq"] and manifest["failed"] == {"short": message}
    assert list(yaml.safe_load((tmp_path / "spec.yaml").read_text())) == ["ubq"]


def test_an_empty_input_stream_is_a_step_failure(tmp_path):
    result = _run_inputs(tmp_path, ["gone"], [str(tmp_path / "gone.cif")])
    assert result.returncode == 1 and "no input structures" in result.stderr


# ── the generated step script ────────────────────────────────────────────────

def test_the_script_runs_inputs_then_the_model_then_post_processing(local_config, isolated_cwd, new_pipeline):
    _, script = _build(local_config, isolated_cwd, new_pipeline, name="order", seed=7, num_designs=3, batch_size=2)
    inputs = script.index("pipe_rfd4proteina_inputs.py")
    model = script.index("rfd4-generate")
    post = script.index("pipe_rfd4proteina_postprocess.py")
    propagate = script.index("pipe_propagate_missing.py")
    assert inputs < model < post < propagate
    assert 'hf_ckpt_path="/w/rfd4"' in script
    assert "num_replicates=3" in script and "batch_size=2" in script and "seed=7" in script
    assert "recycle_strategy=medium/front_loaded" in script
    assert "condition_spec.yaml" in script and "hydra.run.dir=" in script


def test_scoring_is_off_unless_tags_are_given(local_config, isolated_cwd, new_pipeline):
    _, script = _build(local_config, isolated_cwd, new_pipeline, name="nometrics")
    assert "'~metrics'" in script and "metrics_tags" not in script


def test_inline_metric_tags_are_quoted_as_a_hydra_list(local_config, isolated_cwd, new_pipeline):
    _, script = _build(local_config, isolated_cwd, new_pipeline, name="metrics",
                       metrics_tags=["design", "hotspot_conditioning"])
    assert "metrics_tags='[design,hotspot_conditioning]'" in script and "'~metrics'" not in script


def test_the_seed_is_omitted_when_none_is_given(local_config, isolated_cwd, new_pipeline):
    _, script = _build(local_config, isolated_cwd, new_pipeline, name="noseed")
    assert "seed=" not in script


def test_a_forwarded_override_reaches_the_command_and_is_announced(local_config, isolated_cwd, new_pipeline):
    _, script = _build(local_config, isolated_cwd, new_pipeline, name="fwd",
                       **{"generation.args.nsteps": 80, "compile_lengths": "null"})
    assert '"generation.args.nsteps=80"' in script and '"compile_lengths=null"' in script
    assert "Forwarding unrecognized argument(s) to RFD4Proteina" in script


def test_the_thinking_budget_selects_the_task_schedule(local_config, isolated_cwd, new_pipeline):
    _, script = _build(local_config, isolated_cwd, new_pipeline, name="think", thinking="high",
                       contig="150, /0, B386", task="small_molecule_binder")
    assert "recycle_strategy=high/midpeak" in script


def test_the_weights_revision_is_printed_to_the_log(local_config, isolated_cwd, new_pipeline):
    _, script = _build(local_config, isolated_cwd, new_pipeline, name="rev")
    assert '"/w/rfd4/REVISION"' in script and "Weights revision" in script


def test_a_missing_bundle_is_reported_instead_of_run(local_config, isolated_cwd, new_pipeline):
    _, script = _build(local_config, isolated_cwd, new_pipeline, name="nobundle")
    assert '"/w/rfd4/model.safetensors"' in script and "weights not found" in script


def test_no_weights_location_is_refused_with_the_fix(local_config, isolated_cwd, new_pipeline):
    pipeline = new_pipeline("noweights")
    with pytest.raises(ValueError, match="no weights bundle configured"):
        with pipeline:
            RFD4Proteina(structures=_structures(), contig="50-120, /0, A17-131")
            pipeline.save()


def test_a_configured_weights_folder_is_the_default(local_config, isolated_cwd, new_pipeline):
    pipeline = new_pipeline("folderdefault")
    with pipeline:
        built = RFD4Proteina(structures=_structures(), contig="50-120, /0, A17-131")
        built.producer.folders = {**built.producer.folders, "RFD4Proteina": "/data/rfd4_proteina"}
        script = built.producer.generate_script("x.sh")
    assert 'hf_ckpt_path="/data/rfd4_proteina"' in script


def test_the_container_prefix_leads_only_the_model_command(local_config, isolated_cwd, new_pipeline, monkeypatch):
    from biopipelines.config_manager import ConfigManager
    monkeypatch.setattr(ConfigManager, "get_container_executor", lambda self: "apptainer")
    pipeline = new_pipeline("container")
    with pipeline:
        built = RFD4Proteina(structures=_structures(), contig="50-120, /0, A17-131", checkpoint="/w/rfd4")
        built.producer.folders = {**built.producer.folders, "container:RFD4Proteina": "/img/rfd4.sif",
                                  "__container_binds__": ["/w"]}
        script = built.producer.generate_script("x.sh")
    line = next(l for l in script.splitlines() if "rfd4-generate" in l and "echo" not in l)
    assert line.strip().startswith("apptainer exec --nv -B /w /img/rfd4.sif rfd4-generate")
    assert "apptainer" not in script.split("pipe_rfd4proteina_inputs.py")[1].split("RFD4_INPUTS_RC")[0]


# ── post-processing ───────────────────────────────────────────────────────────

def test_replicates_are_renumbered_from_one():
    post = _load_script("pipe_rfd4proteina_postprocess")
    assert post.upstream_id("binder", 1) == "binder_0" and post.design_id("binder", 1) == "binder_1"
    assert post.upstream_id("binder", 3) == "binder_2" and post.design_id("binder", 3) == "binder_3"


def test_an_empty_generated_file_is_not_a_design(tmp_path):
    post = _load_script("pipe_rfd4proteina_postprocess")
    folder = tmp_path / "inference_samples" / "t1_0"
    folder.mkdir(parents=True)
    path = folder / "t1_0_generated.cif"
    path.write_text("")
    assert post.locate_generated_cif(str(tmp_path), "t1_0") is None
    path.write_text("data_x\n")
    assert post.locate_generated_cif(str(tmp_path), "t1_0") == str(path)


def _metrics_csv(run_dir, name, rows):
    import pandas as pd
    folder = pathlib.Path(run_dir) / "inference_metrics"
    folder.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(folder / name, index=False)


def test_metrics_are_rekeyed_to_design_ids_and_keep_the_run_record(tmp_path):
    post = _load_script("pipe_rfd4proteina_postprocess")
    _metrics_csv(tmp_path, "per_example_metrics.csv", [
        {"example_id": "t1_0", "designed_structure_path": "/x", "input_json_path": "/spec.yaml",
         "generation.args.nsteps": 60, "seed": 7, "clash_metrics.frac_atoms_clashing_1p5A": 0.01, "ckpt_path": "/w"},
        {"example_id": "t1_1", "designed_structure_path": "/y", "input_json_path": "/spec.yaml",
         "generation.args.nsteps": 60, "seed": 7, "clash_metrics.frac_atoms_clashing_1p5A": 0.02, "ckpt_path": "/w"},
    ])
    table = post.build_metrics_table(["t1_1", "t1_2", "t1_3"], {"t1_1": "t1_0", "t1_2": "t1_1", "t1_3": "t1_2"},
                                     {"t1_1": "t1", "t1_2": "t1", "t1_3": "t1"},
                                     post.read_inference_metrics(str(tmp_path)))
    assert list(table.columns) == ["id", "source_id", "seed", "clash_metrics.frac_atoms_clashing_1p5A", "ckpt_path"]
    assert table["clash_metrics.frac_atoms_clashing_1p5A"].tolist()[:2] == [0.01, 0.02]
    assert table["id"].tolist() == ["t1_1", "t1_2", "t1_3"]
    assert table["seed"].isna().tolist() == [False, False, True]


def test_without_scoring_the_metrics_table_still_lists_every_design(tmp_path):
    post = _load_script("pipe_rfd4proteina_postprocess")
    assert post.read_inference_metrics(str(tmp_path)) is None
    table = post.build_metrics_table(["t1_1", "t1_2"], {}, {"t1_1": "t1", "t1_2": "t1"}, None)
    assert table["id"].tolist() == ["t1_1", "t1_2"] and list(table.columns) == ["id", "source_id"]


def test_re_evaluation_files_do_not_duplicate_a_design(tmp_path):
    post = _load_script("pipe_rfd4proteina_postprocess")
    _metrics_csv(tmp_path, "per_example_metrics.csv", [{"example_id": "t1_0", "m.x": 1.0}])
    _metrics_csv(tmp_path, "per_example_metrics_reeval.csv", [{"example_id": "t1_0", "m.x": 2.0}])
    frame = post.read_inference_metrics(str(tmp_path))
    assert frame["example_id"].tolist() == ["t1_0"] and frame["m.x"].tolist() == [2.0]


def _write_design_cif(path, chains):
    """A designed CIF as Biopython writes it: `chains` maps a chain id to its residue names, one C-alpha each."""
    from Bio.PDB import MMCIFIO
    from Bio.PDB.StructureBuilder import StructureBuilder
    builder = StructureBuilder()
    builder.init_structure("d")
    builder.init_model(0)
    builder.init_seg(" ")
    serial = 1
    for chain_id, residues in chains.items():
        builder.init_chain(chain_id)
        for number, name in enumerate(residues, start=1):
            builder.init_residue(name, " ", number, " ")
            builder.init_atom("CA", [float(serial), 0.0, 0.0], 0.0, 1.0, " ", " CA ", serial, "C")
            serial += 1
    io = MMCIFIO()
    io.set_structure(builder.get_structure())
    path.parent.mkdir(parents=True, exist_ok=True)
    io.save(str(path))


def _run_postprocess(tmp_path, manifest, num_designs):
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    return subprocess.run(
        [sys.executable, str(PIPE_SCRIPTS / "pipe_rfd4proteina_postprocess.py"),
         "--run-dir", str(tmp_path / "run"), "--manifest-json", str(tmp_path / "manifest.json"),
         "--num-designs", str(num_designs), "--structures-dir", str(tmp_path / "structures"),
         "--structures-map-csv", str(tmp_path / "structures" / "structures_map.csv"),
         "--sequences-csv", str(tmp_path / "sequences" / "sequences.csv"),
         "--structures-csv", str(tmp_path / "tables" / "structures.csv"),
         "--metrics-csv", str(tmp_path / "tables" / "metrics.csv"),
         "--local-missing-csv", str(tmp_path / "007_RFD4Proteina" / "execution" / "local_missing.csv")],
        capture_output=True, text=True)


def test_designs_become_pdbs_sequences_and_table_rows(tmp_path):
    pytest.importorskip("Bio")
    import pandas as pd
    _write_design_cif(tmp_path / "run" / "inference_samples" / "t1_0" / "t1_0_generated.cif",
                      {"A": ["ALA", "GLY", "SER"], "B": ["LEU", "LEU"]})
    _write_design_cif(tmp_path / "run" / "inference_samples" / "t1_1" / "t1_1_generated.cif",
                      {"A": ["VAL", "VAL", "VAL"], "B": ["LEU", "LEU"]})
    result = _run_postprocess(tmp_path, {"attempted": ["t1"], "built": ["t1"], "failed": {}}, 2)
    assert result.returncode == 0, result.stderr

    assert (tmp_path / "structures" / "t1_1.pdb").is_file() and (tmp_path / "structures" / "t1_2.pdb").is_file()
    stream = pd.read_csv(tmp_path / "structures" / "structures_map.csv")
    assert stream["id"].tolist() == ["t1_1", "t1_2"] and stream["structures.id"].tolist() == ["t1", "t1"]

    sequences = pd.read_csv(tmp_path / "sequences" / "sequences.csv")
    assert sequences["id"].tolist() == ["t1_1_A", "t1_1_B", "t1_2_A", "t1_2_B"]
    assert sequences["sequence"].tolist() == ["AGS", "LL", "VVV", "LL"]

    table = pd.read_csv(tmp_path / "tables" / "structures.csv")
    assert table["id"].tolist() == ["t1_1", "t1_2"] and table["n_residues"].tolist() == [5, 5]
    assert table["n_chains"].tolist() == [2, 2] and table["design"].tolist() == [1, 2]
    metrics = pd.read_csv(tmp_path / "tables" / "metrics.csv")
    assert metrics["id"].tolist() == ["t1_1", "t1_2"] and metrics["source_id"].tolist() == ["t1", "t1"]
    assert pd.read_csv(tmp_path / "007_RFD4Proteina" / "execution" / "local_missing.csv").empty


def test_a_design_that_was_never_written_is_a_failure_row_not_a_gap(tmp_path):
    pytest.importorskip("Bio")
    import pandas as pd
    _write_design_cif(tmp_path / "run" / "inference_samples" / "t1_0" / "t1_0_generated.cif", {"A": ["ALA"]})
    result = _run_postprocess(tmp_path, {"attempted": ["t1"], "built": ["t1"], "failed": {}}, 2)
    assert result.returncode == 0, result.stderr

    table = pd.read_csv(tmp_path / "tables" / "structures.csv")
    assert table["id"].tolist() == ["t1_1", "t1_2"] and table["n_residues"].isna().tolist() == [False, True]
    assert pd.read_csv(tmp_path / "structures" / "structures_map.csv")["id"].tolist() == ["t1_1"]
    missing = pd.read_csv(tmp_path / "007_RFD4Proteina" / "execution" / "local_missing.csv")
    assert missing["id"].tolist() == ["t1_2"] and missing["kind"].tolist() == ["failure"]
    assert missing["removed_by"].tolist() == ["007_RFD4Proteina"]
    assert "no generated CIF" in missing["cause"][0]


def test_a_failed_input_gets_failure_rows_for_each_of_its_designs(tmp_path):
    pytest.importorskip("Bio")
    import pandas as pd
    _write_design_cif(tmp_path / "run" / "inference_samples" / "ok_0" / "ok_0_generated.cif", {"A": ["ALA"]})
    manifest = {"attempted": ["ok", "bad"], "built": ["ok"], "failed": {"bad": "could not convert"}}
    result = _run_postprocess(tmp_path, manifest, 2)
    assert result.returncode == 0, result.stderr
    missing = pd.read_csv(tmp_path / "007_RFD4Proteina" / "execution" / "local_missing.csv")
    assert set(missing["id"]) == {"ok_2", "bad_1", "bad_2"}
    assert all("input preparation failed: could not convert" in c for c in missing[missing["id"].str.startswith("bad")]["cause"])
    table = pd.read_csv(tmp_path / "tables" / "structures.csv")
    assert table["id"].tolist() == ["ok_1", "ok_2", "bad_1", "bad_2"]


def test_a_step_that_produced_nothing_fails(tmp_path):
    pytest.importorskip("Bio")
    result = _run_postprocess(tmp_path, {"attempted": ["t1"], "built": ["t1"], "failed": {}}, 2)
    assert result.returncode == 1 and "all 2 expected design(s) failed" in result.stderr


def test_a_step_with_no_prepared_input_fails(tmp_path):
    pytest.importorskip("Bio")
    result = _run_postprocess(tmp_path, {"attempted": [], "built": [], "failed": {}}, 2)
    assert result.returncode == 1 and "no designs expected" in result.stderr


def test_a_missing_manifest_fails_with_the_cause(tmp_path):
    pytest.importorskip("Bio")
    result = subprocess.run(
        [sys.executable, str(PIPE_SCRIPTS / "pipe_rfd4proteina_postprocess.py"),
         "--run-dir", str(tmp_path), "--manifest-json", str(tmp_path / "none.json"), "--num-designs", "1",
         "--structures-dir", str(tmp_path), "--structures-map-csv", str(tmp_path / "a.csv"),
         "--sequences-csv", str(tmp_path / "b.csv"), "--structures-csv", str(tmp_path / "c.csv"),
         "--metrics-csv", str(tmp_path / "d.csv"), "--local-missing-csv", str(tmp_path / "e.csv")],
        capture_output=True, text=True)
    assert result.returncode == 1 and "inputs manifest not found" in result.stderr


# ── the weights download ──────────────────────────────────────────────────────

def test_the_pinned_release_is_read_from_the_packages_own_config():
    weights = _load_script("pipe_rfd4proteina_weights")
    config = (
        "defaults:\n  - _self_\n  - recycle_strategy: thinking_medium\n"
        "hf_ckpt_path: nvidia/RFD4-Proteina  # HF repo\n"
        "hf_ckpt_revision: a06e2509c7033b010fc43614e9fd840303588502\n"
        "ckpt_path: ${oc.env:DATA_PATH}/flow_checkpoints\n")
    assert weights.pinned_release(config) == ("nvidia/RFD4-Proteina", "a06e2509c7033b010fc43614e9fd840303588502")


def test_an_unpinned_release_resolves_later_not_to_a_guess():
    weights = _load_script("pipe_rfd4proteina_weights")
    assert weights.pinned_release("hf_ckpt_path: nvidia/RFD4-Proteina\nhf_ckpt_revision: null\n") == (
        "nvidia/RFD4-Proteina", None)


def test_a_config_without_a_repository_is_refused():
    with pytest.raises(ValueError, match="hf_ckpt_path"):
        _load_script("pipe_rfd4proteina_weights").pinned_release("seed: 1\n")


def test_the_bundle_is_exactly_the_files_upstreams_loader_reads():
    weights = _load_script("pipe_rfd4proteina_weights")
    assert {"config.json", "model.safetensors"} <= set(weights.BUNDLE_FILES)
    assert not any("runtime" in name or ".sif" in name for name in weights.BUNDLE_FILES)


# ── install ───────────────────────────────────────────────────────────────────

@pytest.fixture
def install_env(monkeypatch):
    monkeypatch.setattr(RFD4Proteina, "_install_env", classmethod(lambda cls, manager, env_name=None: "RFD4Proteina"))


def _install(force=False, folders=None):
    folders = {"biopipelines": str(ROOT), "pipe_scripts": str(PIPE_SCRIPTS),
               "RFD4Proteina": "/data/rfd4_proteina", **(folders or {})}
    return RFD4Proteina._install_script(folders, env_manager="mamba", force_reinstall=force)


def test_the_install_creates_the_env_then_downloads_the_weights_then_verifies(local_config, install_env):
    script = _install()
    assert "mamba env create" in script
    assert script.index("env create") < script.index("pipe_rfd4proteina_weights.py") < script.index("# Verify installation")
    assert '--dest "/data/rfd4_proteina"' in script
    assert "rfd4-generate" in script and "model.safetensors" in script


def test_the_pip_layers_install_pip_first_then_the_cuda_extra(local_config, install_env):
    script = _install()
    first = script.index("RFD4Proteina.pip.1.txt")
    second = script.index("RFD4Proteina.pip.2.txt")
    assert first < second
    assert "rfd4-proteina[all-cuda13]" in (ROOT / "environments" / "RFD4Proteina.pip.2.txt").read_text()
    assert "pip-26.2.2" in (ROOT / "environments" / "RFD4Proteina.pip.1.txt").read_text()


def test_the_env_spec_is_python_312():
    spec = yaml.safe_load((ROOT / "environments" / "RFD4Proteina.yaml").read_text())
    assert spec["name"] == "RFD4Proteina" and "python=3.12" in spec["dependencies"]


def test_an_installed_env_and_bundle_short_circuit_the_install(local_config, install_env):
    skip = _install().split("already installed, skipping")[0]
    assert "model.safetensors" in skip and "importlib.metadata" in skip
    assert "already installed, skipping" not in _install(force=True)


def test_force_reinstall_removes_the_env_first(local_config, install_env):
    assert "env remove -n RFD4Proteina" in _install(force=True)


def test_an_install_without_a_weights_folder_names_the_config_key(local_config, install_env):
    with pytest.raises(ValueError, match="folders.repositories.RFD4Proteina"):
        RFD4Proteina._install_script({"biopipelines": str(ROOT)}, env_manager="mamba")


def test_every_variant_but_daint_maps_the_env_and_the_weights_folder():
    for variant in ("cluster", "colab", "container"):
        config = yaml.safe_load((ROOT / f"config.{variant}.yaml").read_text())
        assert config["environments"]["RFD4Proteina"] == "RFD4Proteina", variant
        assert config["folders"]["repositories"]["RFD4Proteina"], variant
    daint = yaml.safe_load((ROOT / "config.daint.yaml").read_text())
    assert "RFD4Proteina" not in daint["environments"], "RFD4-Proteina ships x86-64 wheels only"
