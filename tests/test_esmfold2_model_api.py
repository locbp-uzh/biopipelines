# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""The esm >= 3.4 model API that ESMFold2's install check and inference script rely on.

esm 3.3.0 kept the model in a transformers fork (`transformers.models.esmfold2`); esm 3.4 ships it as `esm.models.esmfold2.EsmFold2Model`, loaded with `from_pretrained(name, device=...)` and driven through `ESMFold2InputBuilder.fold`.
These tests pin that the script touches only the esm API, loads the model where the user asked, hands every documented knob to `fold`, and reads the result fields of `MolecularComplexResult` as esm 3.4.2 declares them.
torch and esm are stubbed in `sys.modules` and `transformers` is made unimportable, so none of this needs a GPU, the weights or either package.

Not covered: that the real `EsmFold2Model.from_pretrained` accepts the current `biohub/ESMFold2` checkpoint, and the numbers a real fold returns.
Both need a GPU run in an esm >= 3.4 environment.
"""

import contextlib
import importlib.util
import json
import os
import sys
import types
from dataclasses import dataclass
from typing import Any, Optional

import numpy as np
import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INFERENCE = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_esmfold2_inference.py")

SEQ = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"


@dataclass
class _Result:
    """The fields of MolecularComplexResult in esm 3.4.2 that the script can touch."""

    complex: Any
    plddt: Any = None
    ptm: Optional[float] = None
    iptm: Optional[float] = None
    pae: Any = None
    pde: Any = None
    distogram: Any = None
    pair_chains_iptm: Any = None
    num_tokens: Optional[int] = None


class _Stub:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def _result(plddt=0.8, iptm=0.5, pae=None, pair_chains=None):
    return _Result(
        complex=types.SimpleNamespace(to_mmcif=lambda: "data_x\n"),
        plddt=np.full(5, plddt, dtype=np.float32),
        ptm=0.7, iptm=iptm,
        pae=np.full((5, 5), 4.0, dtype=np.float32) if pae is None else pae,
        pair_chains_iptm=pair_chains, num_tokens=5)


@pytest.fixture
def esm_api(monkeypatch):
    """Install the fakes; returns the record of load and fold calls and a setter for what fold yields."""
    record = types.SimpleNamespace(loads=[], folds=[], evals=[], results=[_result()])

    class FakeModel:
        @classmethod
        def from_pretrained(cls, name, device=None, **kwargs):
            record.loads.append((name, device))
            return cls()

        def eval(self):
            record.evals.append(True)
            return self

    class FakeBuilder:
        def fold(self, model, spi, **kwargs):
            record.folds.append(kwargs)
            return list(record.results)

    builder = types.SimpleNamespace(
        ProteinInput=_Stub, DNAInput=_Stub, RNAInput=_Stub, LigandInput=_Stub,
        Modification=_Stub, CovalentBond=_Stub, StructurePredictionInput=_Stub)
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
    # None in sys.modules makes `import transformers` raise ImportError, whether or not it is installed.
    for name in ("transformers", "transformers.models", "transformers.models.esmfold2",
                 "transformers.models.esmfold2.modeling_esmfold2"):
        monkeypatch.setitem(sys.modules, name, None)
    return record


@pytest.fixture
def inference():
    spec = importlib.util.spec_from_file_location("pipe_esmfold2_inference", INFERENCE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _axis(tmp_path, name, entity_type, rows):
    table = tmp_path / f"{name}.csv"
    table.write_text("id,sequence\n" + "".join(f"{i},{s}\n" for i, s in rows))
    return {"mode": "each", "entity_type": entity_type, "sources": [str(table)]}


def _run(module, tmp_path, monkeypatch, *flags, with_dna=False):
    """Fold protein p1 (plus a DNA chain when with_dna, which makes it a complex); returns the output directory."""
    axes = {"proteins": _axis(tmp_path, "proteins", "protein", [("p1", SEQ)])}
    if with_dna:
        axes["ssDNA"] = _axis(tmp_path, "ssdna", "ssdna", [("d1", "GATAGCGCTATC")])
    config = tmp_path / "comb.json"
    config.write_text(json.dumps({"axes": axes}))
    out = tmp_path / "structures"
    monkeypatch.setattr(sys, "argv", [
        "pipe_esmfold2_inference.py", "--combinatorics-config", str(config),
        "--output-dir", str(out), *flags])
    module.main()
    return out


class TestLoadingTheModel:

    def test_it_comes_from_esm_and_goes_straight_to_the_gpu(self, inference, esm_api, tmp_path, monkeypatch):
        _run(inference, tmp_path, monkeypatch, "--model-name", "lab/ESMFold2-custom")
        assert esm_api.loads == [("lab/ESMFold2-custom", "cuda")]
        assert esm_api.evals == [True]

    def test_the_default_checkpoint_is_the_hub_one(self, inference, esm_api, tmp_path, monkeypatch):
        _run(inference, tmp_path, monkeypatch)
        assert esm_api.loads == [("biohub/ESMFold2", "cuda")]

    def test_nothing_is_imported_from_transformers(self, inference, esm_api, tmp_path, monkeypatch):
        """The fixture makes `import transformers` raise, so a run that completes never asked for it."""
        out = _run(inference, tmp_path, monkeypatch)
        assert (out / "p1.cif").exists()

    def test_an_esm_without_the_model_stops_with_the_reinstall_advice(
            self, inference, esm_api, tmp_path, monkeypatch):
        """esm 3.3.0 has ESMFold2InputBuilder but not EsmFold2Model, and an existing env skips install."""
        del sys.modules["esm.models.esmfold2"].EsmFold2Model
        with pytest.raises(SystemExit) as stopped:
            _run(inference, tmp_path, monkeypatch)
        message = str(stopped.value)
        assert "esm >= 3.4" in message and "force_reinstall=True" in message
        assert esm_api.loads == [] and esm_api.folds == []


class TestFolding:

    def test_every_seed_is_folded_with_the_knobs_the_wrapper_sets(
            self, inference, esm_api, tmp_path, monkeypatch):
        _run(inference, tmp_path, monkeypatch, "--num-seeds", "2", "--seed-offset", "5",
             "--num-loops", "3", "--num-sampling-steps", "7",
             "--num-diffusion-samples", "1", "--msa-max-depth", "256")
        assert [call["seed"] for call in esm_api.folds] == [5, 6]
        for call in esm_api.folds:
            assert call["num_loops"] == 3 and call["num_sampling_steps"] == 7
            assert call["num_diffusion_samples"] == 1
            assert call["msa_max_depth"] == 256

    def test_the_default_msa_depth_is_the_one_upstream_defaults_to(
            self, inference, esm_api, tmp_path, monkeypatch):
        _run(inference, tmp_path, monkeypatch)
        assert esm_api.folds[0]["msa_max_depth"] == 1024

    def test_a_single_result_is_accepted_as_well_as_a_list(
            self, inference, esm_api, tmp_path, monkeypatch):
        """fold() returns one result for one diffusion sample and a list for several."""
        class Single:
            def fold(self, model, spi, **kwargs):
                return _result()
        sys.modules["esm.models.esmfold2"].ESMFold2InputBuilder = Single
        out = _run(inference, tmp_path, monkeypatch)
        assert (out / "p1.cif").exists()


class TestReadingTheResult:

    def _scores(self, inference, esm_api, tmp_path, monkeypatch, *flags, **result):
        esm_api.results = [_result(**result)]
        out = _run(inference, tmp_path, monkeypatch, *flags)
        return json.loads((out / "ESMFold2_scores.json").read_text())

    def test_plddt_is_rescaled_to_percent_and_ptm_and_iptm_are_read_as_declared(
            self, inference, esm_api, tmp_path, monkeypatch):
        row, = self._scores(inference, esm_api, tmp_path, monkeypatch, plddt=0.8, iptm=0.5)
        assert row["id"] == "p1"
        assert row["plddt"] == pytest.approx(80.0)
        assert row["ptm"] == pytest.approx(0.7) and row["iptm"] == pytest.approx(0.5)
        assert "max_pae" not in row

    def test_max_pae_is_the_largest_entry_of_the_token_by_token_matrix(
            self, inference, esm_api, tmp_path, monkeypatch):
        pae = np.full((5, 5), 4.0, dtype=np.float32)
        pae[1, 3] = 21.5
        row, = self._scores(inference, esm_api, tmp_path, monkeypatch, "--include-pae", pae=pae)
        assert row["max_pae"] == pytest.approx(21.5)

    def test_the_chain_pair_matrix_becomes_one_column_per_unordered_pair(
            self, inference, esm_api, tmp_path, monkeypatch):
        pairs = np.array([[0.0, 0.6, 0.2], [0.6, 0.0, 0.3], [0.2, 0.3, 0.0]], dtype=np.float32)
        row, = self._scores(inference, esm_api, tmp_path, monkeypatch, pair_chains=pairs)
        assert row["iptm_chain_0_1"] == pytest.approx(0.6)
        assert row["iptm_chain_0_2"] == pytest.approx(0.2)
        assert row["iptm_chain_1_2"] == pytest.approx(0.3)
        assert "iptm_chain_1_0" not in row

    def test_a_monomer_is_ranked_by_plddt(self, inference, esm_api, tmp_path, monkeypatch):
        esm_api.results = [_result(plddt=0.5, iptm=0.9), _result(plddt=0.9, iptm=0.1)]
        out = _run(inference, tmp_path, monkeypatch, "--num-diffusion-samples", "2")
        row, = json.loads((out / "ESMFold2_scores.json").read_text())
        assert row["plddt"] == pytest.approx(90.0)

    def test_a_complex_is_ranked_by_iptm(self, inference, esm_api, tmp_path, monkeypatch):
        esm_api.results = [_result(plddt=0.9, iptm=0.1), _result(plddt=0.5, iptm=0.9)]
        out = _run(inference, tmp_path, monkeypatch, "--num-diffusion-samples", "2", with_dna=True)
        row, = json.loads((out / "ESMFold2_scores.json").read_text())
        assert row["iptm"] == pytest.approx(0.9)


class TestTheInstallCheck:

    @pytest.fixture(autouse=True)
    def _env(self, monkeypatch):
        from biopipelines.esmfold2 import ESMFold2
        monkeypatch.setattr(ESMFold2, "_install_env", classmethod(lambda cls, mgr: "esmfold2"))

    def _script(self, **kwargs):
        from biopipelines.esmfold2 import ESMFold2
        return ESMFold2._install_script({"biopipelines": "/repo"}, env_manager="mamba", **kwargs)

    @pytest.mark.parametrize("force", [False, True])
    def test_it_imports_the_esm_model_class_not_the_transformers_one(self, local_config, force):
        script = self._script(force_reinstall=force)
        assert "from esm.models.esmfold2 import EsmFold2Model, ESMFold2InputBuilder" in script
        assert "transformers.models" not in script

    def test_a_failed_import_names_the_minimum_version(self, local_config):
        script = self._script()
        assert "needs esm >= 3.4" in script.split("verification failed")[1].split("exit 1")[0]
