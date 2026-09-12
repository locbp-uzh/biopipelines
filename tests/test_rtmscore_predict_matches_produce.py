# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Nothing predicts RTMScore's output ids, and nothing can.

`pipe_rtmscore.py` composes a `<protein>+<ligand>` pair id at runtime, which is the same hand-rolled two-axis join that broke Mock. It cannot break the same way here, for two reasons this file pins.

`get_output_files` declares tables only -- no streams -- and `pipe_check_completion._collect_declared_output_ids` reads declared ids from stream categories, skipping `tables`. So RTMScore contributes zero declared ids, the pair id is compared against nothing, and there is no predicted set for it to disagree with. The pair id is the scores table's join key; the "+" convention is what keeps `id_map_utils` able to match it back to either input axis.

Both axes are also always iterated: the constructor takes only a `DataStream` or a `StandardizedOutput`, so a `Bundle`/`Each` wrapper is refused before anything runs. Those assertions are the alarm. If the gate ever widens, or if RTMScore ever grows a predicted stream, the flat product needs the treatment Mock got -- axis mode into the config, `predict_single_output_id` at runtime -- and these tests are what say so.

One thing to know before wiring a predictor: RTMScore's id set is not knowable at configuration time. When the two input streams carry identical id lists the runner takes the diagonal instead of the product, because a ligand carved out of each complex shares its structure's id and only the paired rows are meaningful. That branch reads the id lists, so any config-time predictor would have to reproduce the same data-dependent decision rather than assume a product.
"""

import importlib.util
import os

import pytest

from biopipelines.combinatorics import Bundle, Each, predict_single_output_id
from biopipelines.datastream import DataStream, create_map_table
from biopipelines.rtmscore import RTMScore


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

AXIS_MODES = [
    ("each", "each", Each, Each),
    ("bundle", "each", Bundle, Each),
    ("each", "bundle", Each, Bundle),
    ("bundle", "bundle", Bundle, Bundle),
]


def _load(name, relative):
    spec = importlib.util.spec_from_file_location(name, os.path.join(REPO_ROOT, *relative))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _runner():
    return _load("pipe_rtmscore_under_test", ("pipe_scripts", "pipe_rtmscore.py"))


def _stream(tmp_path, name, ids, fmt, ext):
    table = tmp_path / f"{name}_{fmt}_map.csv"
    create_map_table(str(table), ids=list(ids))
    for one in ids:
        (tmp_path / f"{one}.{ext}").write_text("ATOM\n")
    return DataStream(name=name, ids=list(ids), map_table=str(table), format=fmt,
                      files=[str(tmp_path / f"<id>.{ext}")])


@pytest.fixture
def proteins(tmp_path):
    return _stream(tmp_path, "structures", ["prot1", "prot2"], "pdb", "pdb")


@pytest.fixture
def ligands(tmp_path):
    # A ligand coordinate stream is a structures stream: that is what `ligands.streams.structures` hands over.
    return _stream(tmp_path, "structures", ["LIG1", "LIG2", "LIG3"], "sdf", "sdf")


def _produce(protein_ids, ligand_ids):
    """Replay the runner's composition over the cross-product branch."""
    compose = _runner().compose_pair_id
    return [compose(pid, lid) for pid in protein_ids for lid in ligand_ids]


# ── nothing predicts these ids ────────────────────────────────────────────────

def test_rtmscore_declares_tables_but_no_streams(proteins, ligands, record_case):
    tool = RTMScore(structures=proteins, ligands_3d=ligands)
    tool.validate_params()
    out = tool.get_output_files()

    stream_categories = [k for k in out if k not in ("tables", "output_folder")]
    record_case(input="RTMScore(structures[2], ligands_3d[3]).get_output_files()",
                expected="tables only, no stream categories", actual=sorted(out.keys()))
    assert sorted(out.keys()) == ["output_folder", "tables"]
    assert stream_categories == []
    assert sorted(out["tables"]) == ["missing", "scores"]


def test_the_completion_check_derives_no_declared_ids_from_rtmscore(proteins, ligands, record_case):
    """With no streams there is nothing for the pair id to be checked against, so it cannot be marked FAILED for a count mismatch."""
    tool = RTMScore(structures=proteins, ligands_3d=ligands)
    tool.validate_params()
    expected_outputs = {"tables": {name: info.to_dict() for name, info in
                                   tool.get_output_files()["tables"].items()},
                        "output_folder": tool.output_folder}

    checker = _load("pipe_check_completion_under_test",
                    ("pipe_scripts", "pipe_check_completion.py"))
    declared_ids, _map_tables = checker._collect_declared_output_ids(expected_outputs)

    record_case(input="_collect_declared_output_ids(RTMScore expected outputs)",
                expected=[], actual=declared_ids)
    assert declared_ids == []


# ── why no axis mode can reach the runner ─────────────────────────────────────

@pytest.mark.parametrize("structures_mode, ligand_mode, wrap_structures, wrap_ligand", AXIS_MODES)
def test_a_wrapped_axis_is_refused_at_construction(
    proteins, ligands, structures_mode, ligand_mode, wrap_structures, wrap_ligand, record_case,
):
    with pytest.raises(ValueError) as excinfo:
        RTMScore(structures=wrap_structures(proteins), ligands_3d=wrap_ligand(ligands))

    record_case(input=f"structures={structures_mode}, ligands_3d={ligand_mode}",
                expected="ValueError naming the accepted input types",
                actual=str(excinfo.value))
    assert "DataStream or StandardizedOutput" in str(excinfo.value)


@pytest.mark.parametrize("wrapper", [Bundle, Each])
def test_a_wrapped_ligand_axis_alone_is_refused_too(proteins, ligands, wrapper):
    """The structures gate fires first, so check the ligand gate on its own with a bare structures axis."""
    with pytest.raises(ValueError, match="ligands_3d must be a DataStream"):
        RTMScore(structures=proteins, ligands_3d=wrapper(ligands)).validate_params()


# ── the composed pair id is the framework's product id ────────────────────────

def test_the_pair_id_matches_the_shared_composer_for_two_iterated_axes(
    proteins, ligands, record_case,
):
    through_composer = [
        predict_single_output_id(
            structures=("each", list(proteins.ids), i, [], False),
            ligands=("each", list(ligands.ids), j, [], False))
        for i in range(len(proteins.ids))
        for j in range(len(ligands.ids))
    ]
    produced = _produce(proteins.ids, ligands.ids)

    record_case(input="structures[2] x ligands[3], both iterated",
                expected=through_composer, actual=produced)
    assert produced == ["prot1+LIG1", "prot1+LIG2", "prot1+LIG3",
                        "prot2+LIG1", "prot2+LIG2", "prot2+LIG3"]
    assert produced == through_composer


def test_matching_id_lists_take_the_diagonal_not_the_product(tmp_path, record_case):
    """The runner's own branch: identical id lists mean a carved ligand per complex, where only the paired rows score against the right protein."""
    ids = ["cplx1", "cplx2", "cplx3"]
    prot_pairs = [(one, f"{one}.pdb") for one in ids]
    lig_pairs = [(one, f"{one}.sdf") for one in ids]

    paired = ([pid for pid, _ in prot_pairs] == [lid for lid, _ in lig_pairs]
              and len(prot_pairs) > 1)
    compose = _runner().compose_pair_id
    diagonal = [compose(pp[0], lp[0]) for pp, lp in zip(prot_pairs, lig_pairs)]

    record_case(input="structures ids == ligands ids, 3 entries",
                expected="3 diagonal ids, not 9 product ids", actual=diagonal)
    assert paired
    assert diagonal == ["cplx1+cplx1", "cplx2+cplx2", "cplx3+cplx3"]
