# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""DiffDock, DynamicBind and NeuralPLexer: one spelling of a pose id, and the parents come from a column.

A pose id is `<protein>+<ligand>_rank<N>`. Its two axes used to be recovered at the post-processing step by splitting the pair id on the first `+`, which is wrong as soon as an id already carries one (a Boltz2 complex `prot+lig0` docked against `LIG1` gives `prot+lig0+LIG1`, split into `prot` and `lig0+LIG1`). The producer knows both axes when it writes the pair, so it records them and the consumer reads them back.

The pair is composed by `combinatorics.predict_single_output_id` at runtime and `predict_output_ids_with_provenance` at configuration time, and the rank suffix is spelled once, by `combinatorics.ranked_id`, so the declared ids and the written ids cannot drift.
"""

import csv
import importlib.util
import os
import re
import subprocess
import sys

import pytest

from biopipelines import combinatorics
from biopipelines.datastream import DataStream, create_map_table

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PIPE_SCRIPTS = os.path.join(REPO_ROOT, "pipe_scripts")

# Both a plain id and an id that already carries the axis separator, as a Boltz2 complex does.
PROTEINS = ["prot1", "prot2+lig0"]
LIGANDS = ["LIG1", "LIG2"]
RANKS = 2


def _pipe(name):
    spec = importlib.util.spec_from_file_location(f"{name}_ids_under_test", os.path.join(PIPE_SCRIPTS, f"{name}.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _structures(tmp_path):
    folder = tmp_path / "proteins"
    folder.mkdir(exist_ok=True)
    for pid in PROTEINS:
        (folder / f"{pid}.pdb").write_text("ATOM      1  N   ALA A   1       0.0   0.0   0.0\n")
    template = str(folder / "<id>.pdb")
    table = str(tmp_path / "structures_map.csv")
    create_map_table(table, ids=list(PROTEINS), files=[template])
    return DataStream(name="structures", ids=list(PROTEINS), files=[template], map_table=table, format="pdb")


def _compounds(tmp_path):
    table = str(tmp_path / "compounds_map.csv")
    create_map_table(table, ids=list(LIGANDS), additional_columns={"smiles": ["C" * (n + 1) for n in range(len(LIGANDS))]})
    return DataStream(name="compounds", ids=list(LIGANDS), map_table=table, format="csv")


def _compound_files(tmp_path):
    """A compounds stream with one SDF per id, which staging copies without needing RDKit."""
    folder = tmp_path / "ligands"
    folder.mkdir(exist_ok=True)
    for lid in LIGANDS:
        (folder / f"{lid}.sdf").write_text("fake ligand\n")
    template = str(folder / "<id>.sdf")
    table = str(tmp_path / "compound_files_map.csv")
    create_map_table(table, ids=list(LIGANDS), files=[template])
    return DataStream(name="compounds", ids=list(LIGANDS), files=[template], map_table=table, format="sdf")


def _expected_poses():
    """(pose id, protein id, ligand id) for every pair and rank, in declared order."""
    return [(f"{p}+{l}_rank{n}", p, l)
            for p in PROTEINS for l in LIGANDS for n in range(1, RANKS + 1)]


def _rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def _provenance(rows):
    return {row["id"]: (row["structures.id"], row["compounds.id"]) for row in rows}


# ── the shared suffix spelling ────────────────────────────────────────────────

def test_ranked_id_appends_the_rank_suffix_to_the_parent():
    assert combinatorics.ranked_id("prot1+LIG1", 3) == "prot1+LIG1_rank3"


def test_ranked_id_keeps_a_parent_that_carries_the_axis_separator_whole():
    assert combinatorics.ranked_id("prot2+lig0+LIG1", 1) == "prot2+lig0+LIG1_rank1"


@pytest.mark.parametrize("rank", [0, -1])
def test_ranked_id_refuses_a_rank_below_one(rank):
    with pytest.raises(ValueError, match="rank"):
        combinatorics.ranked_id("p+l", rank)


# ── DiffDock ──────────────────────────────────────────────────────────────────

def _diffdock(new_pipeline, tmp_path):
    from biopipelines.diffdock import DiffDock

    pipeline = new_pipeline("docking_ids_diffdock")
    with pipeline:
        dock = DiffDock(structures=_structures(tmp_path), compounds=_compounds(tmp_path),
                        samples_per_complex=RANKS)
        pipeline.save()
    return dock.producer, list(dock.streams.structures.ids)


def _diffdock_run(config, monkeypatch):
    monkeypatch.setattr("sys.argv", [
        "pipe_diffdock_build_csv.py",
        "--structures-json", config.structures_json, "--compounds-json", config.compounds_json,
        "--output-csv", config.input_csv])
    _pipe("pipe_diffdock_build_csv").main()

    # inference writes one index<row>_<name> folder per CSV row, in row order
    for row_index in range(len(_rows(config.input_csv))):
        folder = os.path.join(config.raw_out_folder, f"index{row_index}_pair")
        os.makedirs(folder, exist_ok=True)
        open(os.path.join(folder, "rank1.sdf"), "w").write("pose\n")
        for rank in range(1, RANKS + 1):
            open(os.path.join(folder, f"rank{rank}_confidence-0.{rank}.sdf"), "w").write("pose\n")

    monkeypatch.setattr("sys.argv", [
        "pipe_diffdock_postprocess.py",
        "--raw-out", config.raw_out_folder, "--input-csv", config.input_csv,
        "--structures-folder", config.stream_folder("structures"),
        "--structures-map", config.structures_map,
        "--confidence-csv", config.confidence_csv, "--missing-csv", config.missing_csv,
        "--samples-per-complex", str(RANKS)])
    _pipe("pipe_diffdock_postprocess").main()
    return _rows(config.structures_map), _rows(config.confidence_csv)


def test_diffdock_writes_exactly_the_declared_pose_ids(local_config, isolated_cwd, new_pipeline, monkeypatch, tmp_path):
    config, declared = _diffdock(new_pipeline, tmp_path)
    structures, _ = _diffdock_run(config, monkeypatch)

    assert [r["id"] for r in structures] == declared
    assert declared == [pose for pose, _, _ in _expected_poses()]


def test_diffdock_provenance_names_each_axis_even_when_an_id_carries_a_plus(
    local_config, isolated_cwd, new_pipeline, monkeypatch, tmp_path,
):
    config, _ = _diffdock(new_pipeline, tmp_path)
    structures, confidence = _diffdock_run(config, monkeypatch)

    expected = {pose: (p, l) for pose, p, l in _expected_poses()}
    assert _provenance(structures) == expected
    assert _provenance(confidence) == expected


# ── NeuralPLexer ──────────────────────────────────────────────────────────────

def _neuralplexer(new_pipeline, tmp_path):
    from biopipelines.neuralplexer import NeuralPLexer

    pipeline = new_pipeline("docking_ids_neuralplexer")
    with pipeline:
        dock = NeuralPLexer(structures=_structures(tmp_path), compounds=_compound_files(tmp_path),
                            n_samples=RANKS, chunk_size=1)
        pipeline.save()
    return dock.producer, list(dock.streams.structures.ids)


def _neuralplexer_run(config, monkeypatch):
    done = subprocess.run(
        [sys.executable, os.path.join(PIPE_SCRIPTS, "pipe_neuralplexer_stage.py"),
         "--structures-json", config.structures_json, "--compounds-json", config.compounds_json,
         "--staging-folder", config.staging_folder, "--pairs-csv", config.pairs_csv],
        capture_output=True, text=True)
    assert done.returncode == 0, done.stdout + done.stderr

    # upstream writes one folder per staged pair with the ranked prot/lig files inside
    for pair_id in sorted(os.listdir(config.staging_folder)):
        folder = os.path.join(config.raw_out_folder, pair_id)
        os.makedirs(folder, exist_ok=True)
        for rank in range(1, RANKS + 1):
            open(os.path.join(folder, f"prot_rank{rank}_plddt0.{rank}.pdb"), "w").write("ATOM\n")
            open(os.path.join(folder, f"lig_rank{rank}_plddt0.{rank}.sdf"), "w").write("lig\n")

    postprocess = _pipe("pipe_neuralplexer_postprocess")
    # merging reads the ligand through RDKit; the ids are what this test is about
    postprocess.merge_protein_ligand_pdb = lambda prot, lig, out: open(out, "w").write("END\n")
    monkeypatch.setattr("sys.argv", [
        "pipe_neuralplexer_postprocess.py",
        "--raw-out", config.raw_out_folder, "--pairs-csv", config.pairs_csv,
        "--structures-folder", config.stream_folder("structures"),
        "--structures-map", config.structures_map,
        "--confidence-csv", config.confidence_csv, "--missing-csv", config.missing_csv,
        "--n-samples", str(RANKS)])
    postprocess.main()
    return _rows(config.structures_map), _rows(config.confidence_csv)


def test_neuralplexer_writes_exactly_the_declared_pose_ids(local_config, isolated_cwd, new_pipeline, monkeypatch, tmp_path):
    config, declared = _neuralplexer(new_pipeline, tmp_path)
    structures, _ = _neuralplexer_run(config, monkeypatch)

    assert sorted(r["id"] for r in structures) == sorted(declared)
    assert declared == [pose for pose, _, _ in _expected_poses()]


def test_neuralplexer_provenance_names_each_axis_even_when_an_id_carries_a_plus(
    local_config, isolated_cwd, new_pipeline, monkeypatch, tmp_path,
):
    config, _ = _neuralplexer(new_pipeline, tmp_path)
    structures, confidence = _neuralplexer_run(config, monkeypatch)

    expected = {pose: (p, l) for pose, p, l in _expected_poses()}
    assert _provenance(structures) == expected
    assert _provenance(confidence) == expected


# ── DynamicBind ───────────────────────────────────────────────────────────────

def _dynamicbind(new_pipeline, tmp_path):
    from biopipelines.dynamicbind import DynamicBind

    pipeline = new_pipeline("docking_ids_dynamicbind")
    with pipeline:
        dock = DynamicBind(structures=_structures(tmp_path), compounds=_compounds(tmp_path),
                           num_samples=RANKS, num_saved=RANKS)
        pipeline.save()
    return dock.producer, list(dock.streams.structures.ids)


def _dynamicbind_run(config, monkeypatch):
    monkeypatch.setattr("sys.argv", [
        "pipe_dynamicbind_build_csv.py",
        "--structures-json", config.structures_json, "--compounds-json", config.compounds_json,
        "--output-dir", config.ligand_csvs_dir])
    _pipe("pipe_dynamicbind_build_csv").main()

    for csv_name in sorted(os.listdir(config.ligand_csvs_dir)):
        prot_id = csv_name[:-4]
        for lig_idx in range(len(LIGANDS)):
            folder = os.path.join(config.results_root, prot_id, f"index0_idx_{lig_idx}")
            os.makedirs(folder, exist_ok=True)
            for rank in range(1, RANKS + 1):
                open(os.path.join(folder, f"rank{rank}_ligand_lddt0.9{rank}_affinity5.0{rank}_relaxed.sdf"), "w").write("pose\n")

    monkeypatch.setattr("sys.argv", [
        "pipe_dynamicbind_postprocess.py",
        "--results-root", config.results_root, "--ligand-csvs-dir", config.ligand_csvs_dir,
        "--structures-folder", config.stream_folder("structures"),
        "--structures-map", config.structures_map,
        "--affinity-csv", config.affinity_csv, "--missing-csv", config.missing_csv,
        "--num-saved", str(RANKS)])
    _pipe("pipe_dynamicbind_postprocess").main()
    return _rows(config.structures_map), _rows(config.affinity_csv)


def test_dynamicbind_writes_exactly_the_declared_pose_ids(local_config, isolated_cwd, new_pipeline, monkeypatch, tmp_path):
    config, declared = _dynamicbind(new_pipeline, tmp_path)
    structures, _ = _dynamicbind_run(config, monkeypatch)

    assert sorted(r["id"] for r in structures) == sorted(declared)
    assert declared == [pose for pose, _, _ in _expected_poses()]


def test_dynamicbind_provenance_names_each_axis_even_when_an_id_carries_a_plus(
    local_config, isolated_cwd, new_pipeline, monkeypatch, tmp_path,
):
    config, _ = _dynamicbind(new_pipeline, tmp_path)
    structures, affinity = _dynamicbind_run(config, monkeypatch)

    expected = {pose: (p, l) for pose, p, l in _expected_poses()}
    assert _provenance(structures) == expected
    assert _provenance(affinity) == expected


# ── nobody spells the suffix or splits the pair any more ──────────────────────

SOURCES = [
    "biopipelines/diffdock.py", "biopipelines/dynamicbind.py", "biopipelines/neuralplexer.py",
    "pipe_scripts/pipe_diffdock_build_csv.py", "pipe_scripts/pipe_diffdock_postprocess.py",
    "pipe_scripts/pipe_dynamicbind_build_csv.py", "pipe_scripts/pipe_dynamicbind_postprocess.py",
    "pipe_scripts/pipe_neuralplexer_stage.py", "pipe_scripts/pipe_neuralplexer_postprocess.py",
]


@pytest.mark.parametrize("relative", SOURCES)
def test_no_source_spells_a_rank_suffix_or_splits_a_pair_id(relative):
    with open(os.path.join(REPO_ROOT, relative), encoding="utf-8") as f:
        code = f.read()
    assert not re.search(r"\}_rank\{", code), "the rank suffix is spelled by combinatorics.ranked_id"
    assert not re.search(r"\{(prot|protein)_id\}\+\{(lig|ligand)_id\}", code), "the pair id is composed by combinatorics"
    assert not re.search(r"\.split\(\s*[\"']\+[\"']", code), "a pair id is not split back into its axes"
