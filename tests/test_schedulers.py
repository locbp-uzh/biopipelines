"""Tests for the LSF and PBS scheduler backends.

Two layers:

* Backend unit tests (``biopipelines.schedulers``): directive text for
  headers, GPU specs (incl. best-effort translation + warnings),
  dependencies, options, email — no Pipeline involved.
* End-to-end script generation: build a multi-batch DAG (with a Parallel
  block) under an ``lsf`` / ``pbs`` fixture config and assert the emitted
  ``<scheduler>_batch*.sh`` carry native directives and never leak
  ``#SBATCH``.

The SLURM backend's byte-for-byte output is covered by test_parallel.py /
test_service.py, so it is not re-checked here.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from biopipelines.schedulers import get_backend, BATCH_SCHEDULERS


# ── backend unit tests ────────────────────────────────────────────────────────

def test_known_backends():
    assert BATCH_SCHEDULERS == ("slurm", "lsf", "pbs")


def test_get_backend_rejects_non_batch():
    for name in ("colab", "none", "bogus"):
        with pytest.raises(ValueError):
            get_backend(name)


def test_env_var_overrides_variant_detection(monkeypatch):
    # BIOPIPELINES_CONFIG_VARIANT (the same var the submit wrapper reads) must
    # win, so the bash and Python sides agree on the active config variant.
    from biopipelines.config_manager import _autodetect_variant
    monkeypatch.setenv("BIOPIPELINES_CONFIG_VARIANT", "lsf")
    assert _autodetect_variant() == "lsf"
    monkeypatch.delenv("BIOPIPELINES_CONFIG_VARIANT", raising=False)
    # Without the env var it falls back to runtime detection (not "lsf").
    assert _autodetect_variant() != "lsf"


def test_memory_to_mb():
    from biopipelines.schedulers import _memory_to_mb
    assert _memory_to_mb("16GB") == 16000
    assert _memory_to_mb("512MB") == 512
    assert _memory_to_mb("8000") == 8000
    assert _memory_to_mb("1.5GB") == 1500
    assert _memory_to_mb("weird") is None


def test_hms_to_hhmm():
    from biopipelines.schedulers import _hms_to_hhmm
    assert _hms_to_hhmm("24:00:00") == "24:00"
    assert _hms_to_hhmm("4:30:00") == "4:30"
    assert _hms_to_hhmm("90:00") == "90:00"
    assert _hms_to_hhmm("1-00:00:00") == "24:00"  # day prefix folded into hours
    # Unparseable walltime is passed through unchanged (best-effort contract).
    assert _hms_to_hhmm("12h") == "12h"


@pytest.mark.parametrize("name,opts_key", [
    ("slurm", "slurm_options"),
    ("lsf", "lsf_options"),
    ("pbs", "pbs_options"),
])
def test_options_key(name, opts_key):
    assert get_backend(name).options_key == opts_key


def test_lsf_headers_and_gpu():
    b = get_backend("lsf")
    # Memory -> MB on -M plus a matching rusage; walltime -> [HH:]MM (no seconds).
    assert b.header_directives("16GB", "12:00:00", "job.out") == (
        '#BSUB -M 16000\n#BSUB -R "rusage[mem=16000]"\n#BSUB -W 12:00\n#BSUB -o job.out'
    )
    # named model translates cleanly, no warning
    directive, warns = b.gpu_directive("A100", 2)
    assert directive == '#BSUB -gpu "num=2:gmodel=A100"'
    assert warns == []
    # generic count
    assert b.gpu_directive("gpu", 1)[0] == '#BSUB -gpu "num=1"'
    # constraint-style spec falls back to count-only + warning
    directive, warns = b.gpu_directive("high-memory", 1)
    assert directive == '#BSUB -gpu "num=1"'
    assert warns and "high-memory" in warns[0]


def test_pbs_headers_and_gpu():
    b = get_backend("pbs")
    assert b.header_directives("16GB", "12:00:00", "job.out") == (
        "#PBS -l mem=16GB\n#PBS -l walltime=12:00:00\n#PBS -o job.out"
    )
    directive, warns = b.gpu_directive("A100", 2)
    assert directive == "#PBS -l select=1:ngpus=2:gpu_model=A100"
    assert warns == []
    directive, warns = b.gpu_directive("80GB|96GB", 1)
    assert directive == "#PBS -l select=1:ngpus=1"
    assert warns and "80GB|96GB" in warns[0]


def test_lsf_dependency_and_options_and_email():
    b = get_backend("lsf")
    dep = b.dependency_directive(["<JOBID_BATCH_001>", "<JOBID_BATCH_002>"], ["<JOBID_BATCH_003>"])
    assert dep == (
        "\n#BSUB -w 'done(<JOBID_BATCH_001>) && done(<JOBID_BATCH_002>) "
        "&& started(<JOBID_BATCH_003>)'"
    )
    assert b.dependency_directive([], []) == ""
    assert b.extra_options({"q": "normal", "cpus": 4}) == "\n#BSUB -q normal\n#BSUB -n 4"
    assert b.email_directive("x@y.com") == "\n#BSUB -N\n#BSUB -u x@y.com"
    assert b.email_directive("") == ""


def test_pbs_dependency_and_options_and_email():
    b = get_backend("pbs")
    dep = b.dependency_directive(["<JOBID_BATCH_001>", "<JOBID_BATCH_002>"], ["<JOBID_BATCH_003>"])
    # Types separated by ',', job ids within a type by ':'.
    assert dep == "\n#PBS -W depend=afterok:<JOBID_BATCH_001>:<JOBID_BATCH_002>,after:<JOBID_BATCH_003>"
    assert b.dependency_directive([], []) == ""
    assert b.extra_options({"q": "normal", "cpus": 4}) == "\n#PBS -q normal\n#PBS -l ncpus=4"
    assert b.email_directive("x@y.com") == "\n#PBS -m ae\n#PBS -M x@y.com"


# ── end-to-end script generation ──────────────────────────────────────────────

FIXTURES_DIR = Path(__file__).parent / "fixtures"


@pytest.fixture
def scheduler_config(monkeypatch):
    """Point ConfigManager at an lsf/pbs fixture for the duration of a test."""
    from biopipelines.config_manager import ConfigManager

    def _activate(name: str):
        config_path = FIXTURES_DIR / f"config.{name}_local.yaml"
        assert config_path.exists(), f"Missing fixture: {config_path}"
        ConfigManager._instance = None
        ConfigManager._config = None
        ConfigManager._variant = None
        monkeypatch.setattr(
            ConfigManager, "_get_config_path",
            classmethod(lambda cls, variant=None: str(config_path)),
        )
        return config_path

    yield _activate

    ConfigManager._instance = None
    ConfigManager._config = None
    ConfigManager._variant = None


def _make_pipeline(name: str, job: str):
    from biopipelines.pipeline import Pipeline
    return Pipeline(
        project="TestSuite",
        job=job,
        description=f"Scheduler test: {job}",
        on_the_fly=False,
        local_output=True,
        config=f"{name}_local",
    )


@pytest.mark.parametrize("name", ["lsf", "pbs"])
def test_multi_batch_native_directives(name, scheduler_config, isolated_cwd):
    """A Parallel-block DAG emits native dependency directives, no #SBATCH."""
    scheduler_config(name)
    from biopipelines.pipeline import Resources, Parallel
    from biopipelines.mock import Mock

    pipeline = _make_pipeline(name, "multi")
    with pipeline:
        Resources()
        Mock(ids=["pre"], streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
        with Parallel():
            Resources()
            Mock(ids=["s1"], streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
            Resources()
            Mock(ids=["s2"], streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
        Resources()
        Mock(ids=["post"], streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
        pipeline.save()
        pipeline.generate_job_scripts()

    assert pipeline.batch_parents == [[], [0], [0], [1, 2]]

    runtime = Path(pipeline.folders["runtime"])
    batch2 = (runtime / f"{name}_batch2.sh").read_text(encoding="utf-8")
    batch4 = (runtime / f"{name}_batch4.sh").read_text(encoding="utf-8")
    # No SLURM leakage anywhere.
    for i in range(1, 5):
        assert "#SBATCH" not in (runtime / f"{name}_batch{i}.sh").read_text(encoding="utf-8")

    if name == "lsf":
        assert "#BSUB -w 'done(<JOBID_BATCH_001>)'" in batch2
        assert "done(<JOBID_BATCH_002>) && done(<JOBID_BATCH_003>)" in batch4
    else:
        assert "#PBS -W depend=afterok:<JOBID_BATCH_001>" in batch2
        assert "afterok:<JOBID_BATCH_002>:<JOBID_BATCH_003>" in batch4


@pytest.mark.parametrize("name", ["lsf", "pbs"])
def test_options_passthrough(name, scheduler_config, isolated_cwd):
    """Native knobs given to Resources() land under the right options key."""
    scheduler_config(name)
    from biopipelines.pipeline import Resources
    from biopipelines.mock import Mock

    pipeline = _make_pipeline(name, "opts")
    with pipeline:
        Resources(q="normal")
        Mock(ids=["a"], streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
        pipeline.save()
        pipeline.generate_job_scripts()

    assert pipeline.batch_resources[0][f"{name}_options"] == {"q": "normal"}
    script = (Path(pipeline.folders["runtime"]) / f"{name}.sh").read_text(encoding="utf-8")
    if name == "lsf":
        assert "#BSUB -q normal" in script
    else:
        assert "#PBS -q normal" in script
