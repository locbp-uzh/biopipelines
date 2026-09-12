# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""What DiffDock declares it will produce is what DiffDock produces.

DiffDock predicts its pair ids at configuration time through `predict_output_ids_with_provenance`, which honors Bundle and Each, while `pipe_diffdock_build_csv.py` composes `complex_name` itself as a flat product of `f"{prot_id}+{lig_id}"`. A hand-rolled two-axis join has no notion of bundle versus each, so on paper the two sides can drift: a bundled axis should contribute one collapsed prefix (`LIG1+LIG2+prot1`) where a flat product iterates it.

Measured, they cannot drift here. DiffDock's constructor accepts only a DataStream or a StandardizedOutput and raises on anything else, so a Bundle never reaches `get_output_files`, both axes are always iterated, and a flat product is exactly the pair-ID rule for that case. That is a contract rather than an accident: DiffDock's `inference.py` runs one prediction per CSV row, one row per (protein, ligand) pair, so a bundled axis has no meaning to express.

These cases pin the whole invariant, not just today's behavior. `test_predict_equals_produce_for_every_axis_mode` compares for real whenever the tool accepts the input, so the day someone lets a Bundle through without routing the runtime via `biopipelines.combinatorics.predict_single_output_id`, the bundled cases start comparing and fail.
"""

import csv
import importlib.util
import os

import pytest

from biopipelines.combinatorics import Bundle, Each, predict_output_ids_with_provenance
from biopipelines.datastream import DataStream, create_map_table


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

PROTEIN_IDS = ["prot1", "prot2"]
LIGAND_IDS = ["LIG1", "LIG2"]

# A bare input is the documented spelling of an iterated axis, so "each" needs no wrapper.
AXIS_MODES = {"each": lambda stream: stream, "bundle": Bundle}


def _pipe_build_csv():
    spec = importlib.util.spec_from_file_location(
        "pipe_diffdock_build_csv_under_test",
        os.path.join(REPO_ROOT, "pipe_scripts", "pipe_diffdock_build_csv.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _protein_stream(tmp_path, ids=PROTEIN_IDS):
    folder = tmp_path / "proteins"
    folder.mkdir(exist_ok=True)
    for pid in ids:
        (folder / f"{pid}.pdb").write_text("ATOM      1  N   ALA A   1       0.0   0.0   0.0\n")
    template = str(folder / "<id>.pdb")
    table = str(tmp_path / "structures_map.csv")
    create_map_table(table, ids=list(ids), files=[template])
    return DataStream(name="structures", ids=list(ids), files=[template],
                      map_table=table, format="pdb")


def _ligand_stream(tmp_path, ids=LIGAND_IDS):
    table = str(tmp_path / "compounds_map.csv")
    smiles = ["C" * (n + 1) for n in range(len(ids))]
    create_map_table(table, ids=list(ids), additional_columns={"smiles": smiles})
    return DataStream(name="compounds", ids=list(ids), map_table=table, format="csv")


def _unique_pair_ids(ranked_ids):
    """Strip the `_rank<N>` tail the wrapper multiplies each pair id by."""
    pairs = []
    for rid in ranked_ids:
        pair = rid.rsplit("_rank", 1)[0]
        if pair not in pairs:
            pairs.append(pair)
    return pairs


def _declare(new_pipeline, job, structures, compounds, samples=2):
    """Build the tool and return (config, declared pair ids), or (None, error) if it refuses the input."""
    from biopipelines.diffdock import DiffDock

    pipeline = new_pipeline(job)
    try:
        with pipeline:
            dock = DiffDock(structures=structures, compounds=compounds,
                            samples_per_complex=samples)
            pipeline.save()
    except ValueError as exc:
        if "must be DataStream or StandardizedOutput" not in str(exc):
            raise
        return None, str(exc)
    return dock.producer, _unique_pair_ids(dock.streams.structures.ids)


def _produce(config, monkeypatch):
    """Run the real build-csv step on the JSONs the wrapper itself wrote, and return its complex names."""
    monkeypatch.setattr("sys.argv", [
        "pipe_diffdock_build_csv.py",
        "--structures-json", config.structures_json,
        "--compounds-json", config.compounds_json,
        "--output-csv", config.input_csv,
    ])
    _pipe_build_csv().main()
    with open(config.input_csv, newline="") as f:
        return [row["complex_name"] for row in csv.DictReader(f)]


# ── every axis-mode combination ───────────────────────────────────────────────

@pytest.mark.parametrize("structures_mode, compounds_mode", [
    ("each", "each"),
    ("bundle", "each"),
    ("each", "bundle"),
    ("bundle", "bundle"),
])
def test_predict_equals_produce_for_every_axis_mode(
    local_config, isolated_cwd, new_pipeline, monkeypatch, tmp_path,
    structures_mode, compounds_mode,
):
    """Either DiffDock refuses the axis mode outright, or what it declares is what its runtime produces."""
    structures = AXIS_MODES[structures_mode](_protein_stream(tmp_path))
    compounds = AXIS_MODES[compounds_mode](_ligand_stream(tmp_path))

    config, declared = _declare(
        new_pipeline, f"dd_{structures_mode}_{compounds_mode}", structures, compounds)

    if config is None:
        # A mode the wrapper rejects never reaches the runtime, so there is nothing to disagree about.
        assert "must be DataStream or StandardizedOutput" in declared
        return

    produced = _produce(config, monkeypatch)
    assert produced == declared, (
        "DiffDock declared pair ids its runtime does not compose; route "
        "pipe_diffdock_build_csv through combinatorics.predict_single_output_id")


# ── why a flat product suffices today ─────────────────────────────────────────

@pytest.mark.parametrize("wrapper", [Bundle, Each])
@pytest.mark.parametrize("axis", ["structures", "compounds"])
def test_the_wrapper_refuses_a_combinatorics_wrapper(local_config, tmp_path, axis, wrapper):
    """DiffDock docks one pair per CSV row, so both axes must iterate; the refusal is what keeps the two composers in step."""
    from biopipelines.diffdock import DiffDock

    inputs = {"structures": _protein_stream(tmp_path), "compounds": _ligand_stream(tmp_path)}
    inputs[axis] = wrapper(inputs[axis])

    with pytest.raises(ValueError, match=f"{axis} must be DataStream or StandardizedOutput"):
        DiffDock(**inputs)


def test_a_flat_product_could_not_express_a_bundled_axis(tmp_path):
    """The counterfactual, kept as measured evidence: were a Bundle to get through, the hand-rolled join would be wrong."""
    proteins = _protein_stream(tmp_path)
    ligands = _ligand_stream(tmp_path)

    bundled, _ = predict_output_ids_with_provenance(
        structures=(Bundle(proteins), "structures"),
        compounds=(ligands, "compounds"))
    flat_product = [f"{pid}+{lid}" for pid in PROTEIN_IDS for lid in LIGAND_IDS]

    assert bundled == ["prot1+prot2+LIG1", "prot1+prot2+LIG2"]
    assert bundled != flat_product


# ── the ordinary case, end to end ─────────────────────────────────────────────

def test_the_iterated_pair_ids_are_the_complex_names(
    local_config, isolated_cwd, new_pipeline, monkeypatch, tmp_path,
):
    config, declared = _declare(
        new_pipeline, "dd_pairs", _protein_stream(tmp_path), _ligand_stream(tmp_path))

    assert declared == ["prot1+LIG1", "prot1+LIG2", "prot2+LIG1", "prot2+LIG2"]
    assert _produce(config, monkeypatch) == declared


def test_a_single_pair_needs_no_join(
    local_config, isolated_cwd, new_pipeline, monkeypatch, tmp_path,
):
    """One protein and one ligand still carry both axes in the id, as the multi-axis rule requires."""
    config, declared = _declare(
        new_pipeline, "dd_single",
        _protein_stream(tmp_path, ["prot1"]), _ligand_stream(tmp_path, ["LIG1"]))

    assert declared == ["prot1+LIG1"]
    assert _produce(config, monkeypatch) == declared
