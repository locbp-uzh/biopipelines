"""Vina score and minimize mode must hand forwarded kwargs to vina without a NameError.

`_run_vina_score` read `config.get("extra_args")` although `config` is not one of its
parameters, so every Vina `mode="score"` or `mode="minimize"` complex raised NameError, which the
per-complex handler turned into a `failure` row: the step ran and scored nothing.
"""

import importlib.util
import pathlib
import subprocess

import pytest

REPO = pathlib.Path(__file__).resolve().parent.parent


@pytest.fixture(scope="module")
def pipe_gnina():
    """Load the pipe script directly; pipe_scripts is not an importable package."""
    spec = importlib.util.spec_from_file_location("pipe_gnina", REPO / "pipe_scripts" / "pipe_gnina.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_vina_score_passes_the_forwarded_tokens_to_vina(pipe_gnina, monkeypatch, tmp_path):
    seen = {}

    def fake_run(binary, receptor, ligand, out, box, **kwargs):
        seen.update(kwargs)
        return subprocess.CompletedProcess([binary], 0, "Estimated Free Energy of Binding   : -6.5 (kcal/mol)\n", "")

    backend = pipe_gnina.vina_backend
    monkeypatch.setattr(backend, "receptor_pdbqt", lambda src, dst: dst)
    monkeypatch.setattr(backend, "ligand_pdbqt", lambda src, dst: dst)
    monkeypatch.setattr(backend, "box_args", lambda *a, **k: ["--center_x", "0"])
    monkeypatch.setattr(backend, "run", fake_run)
    monkeypatch.setattr(backend, "parse_score_stdout", lambda text: -6.5)

    scores, pose = pipe_gnina._run_vina_score(
        "vina", "protein.pdb", "lig.sdf", str(tmp_path), "complex_1", "score", "vina", 1, 60,
        str(tmp_path / "out.sdf"), extra_args=["--energy_range", "5"])

    assert seen["extra_args"] == ["--energy_range", "5"]
    assert scores == {"vina_affinity": -6.5}
    assert pose == "lig.sdf"
