# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""The ESMFold2 `pae` stream: the matrix written per kept structure, its map, and its declaration.

The inference script keeps one sample per complex (best by pLDDT or ipTM), so the property worth
pinning is that the matrix on disk belongs to the structure on disk, not to some other sample. The
model is replaced by fakes (torch and esm are stubbed in sys.modules), so none of
this needs a GPU or the weights.

Not covered: a real `res.pae` tensor. esm 3.4.2 builds it as `output["pae"][i]`, a (tokens x tokens) tensor
in Angstrom (read from source, not run). The script treats anything but a square 2-D matrix as unwritable
and says so rather than guessing.
"""

import contextlib
import csv
import importlib.util
import json
import os
import subprocess
import sys
import types
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INFERENCE = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_esmfold2_inference.py")
POSTPROCESS = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_esmfold2_postprocessing.py")

SEQ = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def inference():
    return _load(INFERENCE, "pipe_esmfold2_inference")


# ── fake model stack ──────────────────────────────────────────────────────────

class _FakeResult:
    """What ESMFold2InputBuilder.fold returns per diffusion sample."""

    def __init__(self, plddt, pae, iptm=0.5):
        self.plddt = np.full(4, plddt / 100.0)
        self.ptm = 0.7
        self.iptm = iptm
        self.pae = pae
        self.complex = types.SimpleNamespace(to_mmcif=lambda: f"data_plddt_{plddt}\n")


class _Stub:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


@pytest.fixture
def fold_with(monkeypatch):
    """Install the fakes; returns a function that sets what each fold() call yields."""
    queued = {"results": []}

    class FakeModel:
        @classmethod
        def from_pretrained(cls, name, **kwargs):
            return cls()

        def eval(self):
            return self

    class FakeBuilder:
        def fold(self, model, spi, **kwargs):
            if isinstance(queued["results"], Exception):
                raise queued["results"]
            return queued["results"]

    builder = types.SimpleNamespace(ProteinInput=_Stub, DNAInput=_Stub, RNAInput=_Stub,
                                    LigandInput=_Stub, Modification=_Stub, CovalentBond=_Stub,
                                    StructurePredictionInput=_Stub)
    fake_torch = types.ModuleType("torch")
    fake_torch.no_grad = contextlib.nullcontext
    fake_torch.cuda = types.SimpleNamespace(is_available=lambda: False)
    modules = {
        "torch": fake_torch,
        "esm": types.ModuleType("esm"),
        "esm.utils": types.ModuleType("esm.utils"),
        "esm.utils.structure": types.ModuleType("esm.utils.structure"),
        "esm.utils.msa": types.ModuleType("esm.utils.msa"),
        "esm.models": types.ModuleType("esm.models"),
        "esm.models.esmfold2": types.ModuleType("esm.models.esmfold2"),
    }
    modules["esm.utils.structure"].input_builder = builder
    modules["esm.utils.msa"].MSA = _Stub
    modules["esm.models.esmfold2"].ESMFold2InputBuilder = FakeBuilder
    modules["esm.models.esmfold2"].EsmFold2Model = FakeModel
    for name, module in modules.items():
        monkeypatch.setitem(sys.modules, name, module)

    def configure(results):
        queued["results"] = results
    return configure


def _run_inference(module, tmp_path, monkeypatch, *flags):
    """Fold one monomer, "p1", and return the output directory."""
    sequences = tmp_path / "sequences.csv"
    sequences.write_text(f"id,sequence\np1,{SEQ}\n")
    config = tmp_path / "comb.json"
    config.write_text(json.dumps({"axes": {"proteins": {
        "mode": "each", "entity_type": "protein", "sources": [str(sequences)]}}}))
    out = tmp_path / "structures"
    monkeypatch.setattr(sys, "argv", [
        "pipe_esmfold2_inference.py", "--combinatorics-config", str(config),
        "--output-dir", str(out), "--num-diffusion-samples", "2", *flags])
    module.main()
    return out


def _matrix(seed, n=6):
    return np.random.default_rng(seed).uniform(0, 32, size=(n, n)).astype(np.float32)


# ── the matrix ────────────────────────────────────────────────────────────────

class TestTheMatrix:

    def test_a_numpy_matrix_passes_through_as_float32(self, inference):
        out = inference.pae_matrix(_matrix(0).astype(np.float64))
        assert out.dtype == np.float32 and out.shape == (6, 6)

    def test_a_torch_tensor_is_converted_even_when_half_precision(self, inference):
        """The result's tensor is bfloat16 on some builds, which numpy cannot take directly."""
        torch = pytest.importorskip("torch")
        tensor = torch.tensor(_matrix(1), dtype=torch.bfloat16)
        out = inference.pae_matrix(tensor)
        assert out.dtype == np.float32 and out.shape == (6, 6)

    def test_a_duck_typed_tensor_goes_through_detach_float_cpu_numpy(self, inference):
        class Tensor:
            def __init__(self, array):
                self.array = array

            def detach(self):
                return self

            def float(self):
                return self

            def cpu(self):
                return self

            def numpy(self):
                return self.array

        assert np.array_equal(inference.pae_matrix(Tensor(_matrix(2))), _matrix(2))

    @pytest.mark.parametrize("shape", [(3, 6, 6), (6,), (5, 6)])
    def test_anything_but_a_square_matrix_is_refused(self, inference, shape):
        with pytest.raises(ValueError, match="square token x token"):
            inference.pae_matrix(np.zeros(shape))

    def test_the_file_holds_the_matrix_under_the_key_pae(self, inference, tmp_path):
        """IPSAE reads Boltz npz files as np.load(path)['pae']; the layout has to be the same."""
        path = inference.write_pae(str(tmp_path / "pae"), "p1", _matrix(3))
        assert path == str(tmp_path / "pae" / "p1.npz")
        with np.load(path) as archive:
            assert list(archive.files) == ["pae"]
            assert np.array_equal(archive["pae"], _matrix(3))

    def test_a_refused_matrix_leaves_no_folder_behind(self, inference, tmp_path):
        with pytest.raises(ValueError):
            inference.write_pae(str(tmp_path / "pae"), "p1", np.zeros((2, 3)))
        assert not (tmp_path / "pae").exists()


# ── the inference loop ────────────────────────────────────────────────────────

class TestTheInferenceLoop:

    def test_the_matrix_on_disk_belongs_to_the_structure_that_was_kept(
            self, inference, fold_with, tmp_path, monkeypatch):
        weak, strong = _matrix(10), _matrix(11)
        fold_with([_FakeResult(60.0, weak), _FakeResult(80.0, strong)])
        out = _run_inference(inference, tmp_path, monkeypatch,
                             "--include-pae", "--pae-dir", str(tmp_path / "pae"))
        assert (out / "p1.cif").read_text() == "data_plddt_80.0\n"
        with np.load(tmp_path / "pae" / "p1.npz") as archive:
            assert np.array_equal(archive["pae"], strong)
        scores = json.loads((out / "ESMFold2_scores.json").read_text())
        assert scores[0]["max_pae"] == pytest.approx(float(strong.max()))

    def test_every_sample_gets_its_own_matrix_when_all_are_kept(
            self, inference, fold_with, tmp_path, monkeypatch):
        weak, strong = _matrix(12), _matrix(13)
        fold_with([_FakeResult(60.0, weak), _FakeResult(80.0, strong)])
        out = _run_inference(inference, tmp_path, monkeypatch, "--all-samples",
                             "--include-pae", "--pae-dir", str(tmp_path / "pae"))
        assert (out / "p1_1.cif").read_text() == "data_plddt_80.0\n"
        assert (out / "p1_2.cif").read_text() == "data_plddt_60.0\n"
        with np.load(tmp_path / "pae" / "p1_1.npz") as first, \
                np.load(tmp_path / "pae" / "p1_2.npz") as second:
            assert np.array_equal(first["pae"], strong)
            assert np.array_equal(second["pae"], weak)

    def test_without_a_pae_dir_nothing_extra_is_written(
            self, inference, fold_with, tmp_path, monkeypatch):
        """The default run must be byte-for-byte what it was before the stream existed."""
        fold_with([_FakeResult(60.0, _matrix(14)), _FakeResult(80.0, _matrix(15))])
        out = _run_inference(inference, tmp_path, monkeypatch, "--include-pae")
        assert sorted(p.name for p in out.iterdir()) == ["ESMFold2_failures.json", "ESMFold2_scores.json", "p1.cif"]
        assert "max_pae" in json.loads((out / "ESMFold2_scores.json").read_text())[0]

    def test_a_result_without_pae_still_yields_its_structure(
            self, inference, fold_with, tmp_path, monkeypatch):
        fold_with([_FakeResult(60.0, None), _FakeResult(80.0, None)])
        out = _run_inference(inference, tmp_path, monkeypatch,
                             "--include-pae", "--pae-dir", str(tmp_path / "pae"))
        assert (out / "p1.cif").exists()
        assert not (tmp_path / "pae").exists()
        assert "max_pae" not in json.loads((out / "ESMFold2_scores.json").read_text())[0]

    def test_an_unwritable_matrix_costs_the_pae_file_not_the_structure(
            self, inference, fold_with, tmp_path, monkeypatch, capsys):
        fold_with([_FakeResult(60.0, np.zeros((2, 3))), _FakeResult(80.0, np.zeros((2, 3)))])
        out = _run_inference(inference, tmp_path, monkeypatch,
                             "--include-pae", "--pae-dir", str(tmp_path / "pae"))
        assert (out / "p1.cif").exists()
        assert not (tmp_path / "pae" / "p1.npz").exists()
        assert "no PAE matrix written for p1" in capsys.readouterr().err
        assert json.loads((out / "ESMFold2_scores.json").read_text())[0]["id"] == "p1"


# ── the post-processing map ───────────────────────────────────────────────────

def _postprocess(tmp_path, *extra):
    """Two scored ids, p1 and p2, both with a structure; returns the CompletedProcess."""
    predictions = tmp_path / "structures"
    predictions.mkdir(exist_ok=True)
    for cid in ("p1", "p2"):
        (predictions / f"{cid}.cif").write_text("data_x\n")
    (predictions / "ESMFold2_scores.json").write_text(json.dumps([
        {"id": "p1", "plddt": 80.0, "ptm": 0.7, "iptm": 0.5, "max_pae": 12.0},
        {"id": "p2", "plddt": 70.0, "ptm": 0.6, "iptm": 0.4, "max_pae": 14.0}]))
    config = tmp_path / "comb.json"
    config.write_text(json.dumps({"predicted_ids": ["p1", "p2"],
                                  "provenance": {"proteins": ["a", "b"]}}))
    return subprocess.run(
        [sys.executable, POSTPROCESS, "--combinatorics-config", str(config),
         "--predictions-dir", str(predictions),
         "--structures-map", str(tmp_path / "out" / "structures_map.csv"),
         "--confidence-csv", str(tmp_path / "out" / "confidence.csv"), *extra],
        capture_output=True, text=True)


class TestThePaeMap:

    def _pae_dir(self, tmp_path, ids):
        pae_dir = tmp_path / "pae"
        pae_dir.mkdir()
        for cid in ids:
            np.savez_compressed(pae_dir / f"{cid}.npz", pae=_matrix(0))
        return pae_dir

    def test_it_lists_only_the_ids_whose_matrix_exists(self, tmp_path):
        pae_dir = self._pae_dir(tmp_path, ["p1"])
        done = _postprocess(tmp_path, "--pae-dir", str(pae_dir),
                            "--pae-map", str(pae_dir / "pae_map.csv"))
        assert done.returncode == 0, done.stderr
        rows = pd.read_csv(pae_dir / "pae_map.csv", dtype=str)
        assert rows["id"].tolist() == ["p1"]
        assert rows["file"].tolist() == [str(pae_dir / "p1.npz")]
        assert "p2" in done.stderr

    def test_it_carries_the_provenance_columns_the_structures_map_carries(self, tmp_path):
        pae_dir = self._pae_dir(tmp_path, ["p1", "p2"])
        _postprocess(tmp_path, "--pae-dir", str(pae_dir), "--pae-map", str(pae_dir / "pae_map.csv"))
        pae_rows = pd.read_csv(pae_dir / "pae_map.csv", dtype=str)
        structure_rows = pd.read_csv(tmp_path / "out" / "structures_map.csv", dtype=str)
        assert pae_rows.columns.tolist() == structure_rows.columns.tolist()
        assert pae_rows["proteins.id"].tolist() == ["a", "b"]

    def test_the_map_is_written_even_when_no_matrix_exists(self, tmp_path):
        pae_dir = self._pae_dir(tmp_path, [])
        done = _postprocess(tmp_path, "--pae-dir", str(pae_dir),
                            "--pae-map", str(pae_dir / "pae_map.csv"))
        assert done.returncode == 0, done.stderr
        assert len(pd.read_csv(pae_dir / "pae_map.csv")) == 0

    def test_without_the_flags_the_other_outputs_are_unchanged(self, tmp_path):
        plain = tmp_path / "plain"
        plain.mkdir()
        with_pae = tmp_path / "with_pae"
        with_pae.mkdir()
        _postprocess(plain)
        pae_dir = self._pae_dir(with_pae, ["p1", "p2"])
        _postprocess(with_pae, "--pae-dir", str(pae_dir), "--pae-map", str(pae_dir / "pae_map.csv"))
        for name in ("structures_map.csv", "confidence.csv"):
            plain_text = (plain / "out" / name).read_text().replace(str(plain), "")
            pae_text = (with_pae / "out" / name).read_text().replace(str(with_pae), "")
            assert plain_text == pae_text
        assert not (plain / "pae").exists()


# ── the declaration ───────────────────────────────────────────────────────────

def _fold(new_pipeline, **kwargs):
    """Declare one ESMFold2 step in a local pipeline; returns (output, per-step script text)."""
    from biopipelines.esmfold2 import ESMFold2
    from biopipelines.sequence import Sequence
    pipeline = new_pipeline("esmfold2_pae")
    with pipeline:
        tool = ESMFold2(proteins=Sequence(seq=SEQ, type="protein", ids="p1"), **kwargs)
        script_path = pipeline.save()
    steps = sorted(f for f in Path(script_path).parent.glob("*_ESMFold2.sh"))
    return tool, "\n".join(f.read_text(encoding="utf-8", errors="replace") for f in steps)


class TestTheDeclaration:

    def test_the_stream_is_absent_unless_pae_is_requested(self, local_config, isolated_cwd, new_pipeline):
        tool, script = _fold(new_pipeline)
        assert len(tool.streams.pae) == 0
        assert "--pae-dir" not in script and "--pae-map" not in script

    def test_the_stream_mirrors_the_structures_ids(self, local_config, isolated_cwd, new_pipeline):
        tool, _script = _fold(new_pipeline, include_pae=True)
        assert tool.streams.pae.ids == tool.streams.structures.ids == ["p1"]
        assert tool.streams.pae.format == "npz"
        assert tool.streams.pae.files[0].endswith(os.path.join("pae", "<id>.npz"))
        assert tool.streams.pae.map_table.endswith(os.path.join("pae", "pae_map.csv"))

    def test_every_diffusion_sample_is_declared_when_all_are_kept(
            self, local_config, isolated_cwd, new_pipeline):
        tool, _script = _fold(new_pipeline, include_pae=True, top_only=False, num_diffusion_samples=3)
        assert tool.streams.pae.ids == tool.streams.structures.ids == ["p1_<1..3>"]

    def test_the_scripts_hand_the_same_folder_to_inference_and_post_processing(
            self, local_config, isolated_cwd, new_pipeline):
        tool, script = _fold(new_pipeline, include_pae=True)
        folder = os.path.dirname(tool.streams.pae.map_table)
        assert script.count(f'--pae-dir "{folder}"') == 2
        assert f'--pae-map "{tool.streams.pae.map_table}"' in script
        assert "--include-pae" in script

    def test_the_confidence_table_is_unchanged(self, local_config, isolated_cwd, new_pipeline):
        without, _ = _fold(new_pipeline)
        with_pae, _ = _fold(new_pipeline, include_pae=True)
        assert "max_pae" not in without.tables.confidence.info.columns
        assert with_pae.tables.confidence.info.columns == [
            "id", "file", "plddt", "ptm", "iptm", "max_pae"]


# ── a complex that fails to fold ──────────────────────────────────────────────

class TestAComplexThatFailsToFold:
    """One id that raises must leave its trace: a NaN row in the result table and a `failure` row in `missing`."""

    def test_the_inference_script_records_why_a_complex_failed(
            self, inference, fold_with, tmp_path, monkeypatch):
        fold_with(RuntimeError("CUDA out of memory"))
        with pytest.raises(SystemExit):
            _run_inference(inference, tmp_path, monkeypatch)
        failures = json.loads((tmp_path / "structures" / "ESMFold2_failures.json").read_text())
        assert failures == {"p1": "CUDA out of memory"}

    def _postprocess_with_a_failure(self, tmp_path):
        predictions = tmp_path / "structures"
        predictions.mkdir()
        (predictions / "p1.cif").write_text("data_x\n")
        (predictions / "ESMFold2_scores.json").write_text(json.dumps(
            [{"id": "p1", "plddt": 80.0, "ptm": 0.7, "iptm": 0.5}]))
        (predictions / "ESMFold2_failures.json").write_text(
            json.dumps({"p2": "position 60 is past the end of chain A"}))
        config = tmp_path / "comb.json"
        config.write_text(json.dumps({"predicted_ids": ["p1", "p2"], "provenance": {"proteins": ["a", "b"]}}))
        done = subprocess.run(
            [sys.executable, POSTPROCESS, "--combinatorics-config", str(config),
             "--predictions-dir", str(predictions),
             "--structures-map", str(tmp_path / "out" / "structures_map.csv"),
             "--confidence-csv", str(tmp_path / "out" / "confidence.csv"),
             "--local-missing-csv", str(tmp_path / "exec" / "local_missing.csv")],
            capture_output=True, text=True)
        assert done.returncode == 0, done.stderr

    def test_the_result_table_keeps_a_row_for_it_with_nothing_but_the_id(self, tmp_path):
        self._postprocess_with_a_failure(tmp_path)
        table = pd.read_csv(tmp_path / "out" / "confidence.csv")
        assert table["id"].tolist() == ["p1", "p2"]
        assert table.loc[1, ["plddt", "ptm", "iptm"]].isna().all()

    def test_the_stream_map_lists_only_what_was_produced(self, tmp_path):
        self._postprocess_with_a_failure(tmp_path)
        assert pd.read_csv(tmp_path / "out" / "structures_map.csv")["id"].tolist() == ["p1"]

    def test_the_local_missing_file_names_it_as_a_failure_with_the_reason(self, tmp_path):
        self._postprocess_with_a_failure(tmp_path)
        rows = pd.read_csv(tmp_path / "exec" / "local_missing.csv")
        assert rows.to_dict("records") == [{"id": "p2", "removed_by": "ESMFold2", "kind": "failure",
                                            "cause": "position 60 is past the end of chain A"}]

    def test_the_local_missing_file_is_written_even_when_nothing_failed(self, tmp_path):
        done = _postprocess(tmp_path, "--local-missing-csv", str(tmp_path / "exec" / "local_missing.csv"))
        assert done.returncode == 0, done.stderr
        assert list(pd.read_csv(tmp_path / "exec" / "local_missing.csv").columns) == [
            "id", "removed_by", "kind", "cause"]

    def test_the_missing_table_is_always_declared(self, local_config, isolated_cwd, new_pipeline):
        tool, script = _fold(new_pipeline)
        assert tool.tables.missing.info.columns == ["id", "removed_by", "kind", "cause"]
        assert "--local-missing-csv" in script
