# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""What NeuralPLexer declares it will produce is what its staging produces.

`get_output_files` predicts pair ids through `predict_output_ids_with_provenance`, which honors `Bundle` and `Each`; `pipe_neuralplexer_stage.py` re-composes them by hand as `f"{prot_id}+{lig_id}"` in a nested loop, and `pipe_neuralplexer_postprocess.py` splits them back apart on the first `+`. A hand-rolled two-axis join has no notion of bundle versus each, so the two sides can only be trusted to agree while every axis is iterated.

They do agree, because the constructor refuses a `Bundle` (and a bare `Each`) on either axis: `structures` and `compounds` accept a `DataStream` or a `StandardizedOutput` and nothing else. That refusal is the whole reason the flat product at runtime is correct, so it is pinned here alongside the predict-equals-produce comparison. Widen the constructor to accept a wrapper and these tests fail together, which is the signal to move the runtime onto `predict_single_output_id` the way `pipe_mock.py` does.
"""

import importlib.util
import os
import re
import subprocess
import sys

import pytest

pytest.importorskip("rdkit")

from biopipelines.combinatorics import Bundle, Each
from biopipelines.datastream import DataStream, create_map_table
from biopipelines.neuralplexer import NeuralPLexer


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STAGE_SCRIPT = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_neuralplexer_stage.py")

PROTEINS = ["prot1", "prot2"]
LIGANDS = ["LIG1", "LIG2", "LIG3"]

_RANK_SUFFIX = re.compile(r"_rank\d+$")

# One glycine plus a one-heavy-atom HETATM: enough for staging to copy and for RDKit to embed.
SCAFFOLD_PDB = """ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  GLY A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   GLY A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   GLY A   1       1.251   2.390   0.000  1.00  0.00           O
HETATM    5  C1  LIG A   2       5.000   5.000   5.000  1.00  0.00           C
END
"""


def _structures_stream(tmp_path):
    folder = tmp_path / "scaffolds"
    folder.mkdir()
    for pid in PROTEINS:
        (folder / f"{pid}.pdb").write_text(SCAFFOLD_PDB)
    template = str(folder / "<id>.pdb")
    table = tmp_path / "structures_map.csv"
    create_map_table(str(table), ids=list(PROTEINS), files=[template])
    return DataStream(name="structures", ids=list(PROTEINS), files=[template],
                      map_table=str(table), format="pdb")


def _compounds_stream(tmp_path):
    table = tmp_path / "compounds_map.csv"
    create_map_table(str(table), ids=list(LIGANDS),
                     additional_columns={"code": ["LIG"] * len(LIGANDS),
                                         "smiles": ["C"] * len(LIGANDS)})
    return DataStream(name="compounds", ids=list(LIGANDS), files=[],
                      map_table=str(table), format="csv")


def _configured_tool(tmp_path, structures, compounds):
    tool = NeuralPLexer(structures=structures, compounds=compounds, n_samples=2, chunk_size=1)
    tool.configure_inputs({"NeuralPLexer": str(tmp_path / "repo")})
    return tool


def _declared_pair_ids(tool):
    """The pair ids behind the declared structures stream, with the per-rank suffix removed."""
    seen, pairs = set(), []
    for ranked_id in tool.get_output_files()["structures"].ids:
        pair = _RANK_SUFFIX.sub("", ranked_id)
        assert pair != ranked_id, f"declared id {ranked_id!r} carries no _rank<N> suffix"
        if pair not in seen:
            seen.add(pair)
            pairs.append(pair)
    return pairs


def _staged_pair_ids(tmp_path, structures, compounds):
    """Run the real staging script and return the pair folders it created."""
    structures_json = tmp_path / "structures.json"
    compounds_json = tmp_path / "compounds.json"
    structures.save_json(str(structures_json))
    compounds.save_json(str(compounds_json))
    staging = tmp_path / "staging"

    result = subprocess.run(
        [sys.executable, STAGE_SCRIPT,
         "--structures-json", str(structures_json),
         "--compounds-json", str(compounds_json),
         "--staging-folder", str(staging)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, f"staging failed: {result.stderr}"
    return sorted(os.listdir(staging))


def _postprocess():
    spec = importlib.util.spec_from_file_location(
        "pipe_neuralplexer_postprocess_under_test",
        os.path.join(REPO_ROOT, "pipe_scripts", "pipe_neuralplexer_postprocess.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ── predict vs produce ────────────────────────────────────────────────────────

def test_predict_equals_produce_for_iterated_axes(tmp_path, record_case):
    structures = _structures_stream(tmp_path)
    compounds = _compounds_stream(tmp_path)

    declared = _declared_pair_ids(_configured_tool(tmp_path, structures, compounds))
    staged = _staged_pair_ids(tmp_path, structures, compounds)

    record_case(input=f"structures={PROTEINS}, compounds={LIGANDS}",
                expected=sorted(declared), actual=staged)
    assert declared == ["prot1+LIG1", "prot1+LIG2", "prot1+LIG3",
                        "prot2+LIG1", "prot2+LIG2", "prot2+LIG3"]
    assert staged == sorted(declared)


def test_postprocess_recovers_both_axis_ids_from_every_pair_id(tmp_path):
    """The pair id must survive the round trip: staging joins on '+', post-processing splits on it."""
    structures = _structures_stream(tmp_path)
    compounds = _compounds_stream(tmp_path)
    parse_pair_id = _postprocess().parse_pair_id

    for pair_id in _staged_pair_ids(tmp_path, structures, compounds):
        prot_id, lig_id = parse_pair_id(pair_id)
        assert prot_id in PROTEINS
        assert lig_id in LIGANDS


# ── every axis-mode combination ───────────────────────────────────────────────

@pytest.mark.parametrize("structures_mode", [Each, Bundle], ids=["Each", "Bundle"])
@pytest.mark.parametrize("compounds_mode", [Each, Bundle], ids=["Each", "Bundle"])
def test_a_wrapped_axis_is_refused_at_construction(tmp_path, structures_mode, compounds_mode):
    structures = structures_mode(_structures_stream(tmp_path))
    compounds = compounds_mode(_compounds_stream(tmp_path))

    with pytest.raises(ValueError) as excinfo:
        NeuralPLexer(structures=structures, compounds=compounds)
    assert "structures must be DataStream or StandardizedOutput" in str(excinfo.value)


@pytest.mark.parametrize("mode", [Each, Bundle], ids=["Each", "Bundle"])
def test_a_wrapped_compounds_axis_alone_is_refused(tmp_path, mode):
    with pytest.raises(ValueError) as excinfo:
        NeuralPLexer(structures=_structures_stream(tmp_path),
                     compounds=mode(_compounds_stream(tmp_path)))
    assert "compounds must be DataStream or StandardizedOutput" in str(excinfo.value)
