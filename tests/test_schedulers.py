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
test_service.py, so it is not re-checked here beyond the integer ``gpu``
count, which those files never pass.
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
    import getpass

    from biopipelines.config_manager import _autodetect_variant
    monkeypatch.setenv("BIOPIPELINES_CONFIG_VARIANT", "lsf")
    assert _autodetect_variant() == "lsf"

    # Without the env var, and with a user no shipped config claims, detection
    # refuses rather than guessing a variant that merely might fit.
    monkeypatch.delenv("BIOPIPELINES_CONFIG_VARIANT", raising=False)
    monkeypatch.setattr(getpass, "getuser", lambda: "no-such-biopipelines-user")
    with pytest.raises(RuntimeError, match="BIOPIPELINES_CONFIG_VARIANT"):
        _autodetect_variant()


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


@pytest.mark.parametrize("spec", [1, 2, "2", " 4 "])
def test_slurm_integer_gpu_spec_is_a_count(spec):
    # An integer or all-digit string is a count, not a model: it used to raise AttributeError in `.startswith`.
    from biopipelines.schedulers import SlurmBackend
    directive, warns = SlurmBackend().gpu_directive(spec, 1, single_gpu_type=False)
    assert directive == f"#SBATCH --gpus={int(spec)}"
    assert warns == []


@pytest.mark.parametrize("spec, expected", [(2, 2), ("3", 3), (1, 1)])
def test_slurm_integer_gpu_count_survives_a_single_gpu_type_site(spec, expected):
    # The single_gpu_type branch ran before the count branch and asked for one GPU however many were requested.
    from biopipelines.schedulers import SlurmBackend
    directive, warns = SlurmBackend().gpu_directive(spec, 1, single_gpu_type=True)
    assert directive == f"#SBATCH --gres=gpu:{expected}"
    assert warns == []


def test_slurm_single_gpu_type_still_uses_the_gpus_argument_for_a_non_count():
    from biopipelines.schedulers import SlurmBackend
    assert SlurmBackend().gpu_directive("A100", 3, single_gpu_type=True)[0] == "#SBATCH --gres=gpu:3"


def test_16gb_is_a_gpu_memory_constraint_not_a_model():
    # "16GB" is in the Resources() docstring but was absent from MEMORY_SPECS, so it emitted `--gpus=16GB:1`.
    from biopipelines.schedulers import SlurmBackend
    directive, warns = SlurmBackend().gpu_directive("16GB", 1)
    assert directive == '#SBATCH --gpus=1\n#SBATCH --constraint="GPUMEM16GB"'
    assert warns == []


def test_every_documented_gpu_memory_is_a_memory_spec():
    """Resources() advertises these, so each must reach the memory branch rather than the typed one."""
    import re
    from biopipelines import pipeline as pipeline_module
    from biopipelines.schedulers import SlurmBackend
    # Only the gpu: line — the memory: line next to it is node RAM, a different quantity.
    gpu_line = next(line for line in (pipeline_module.Resources.__doc__ or "").splitlines()
                    if line.strip().startswith("gpu:"))
    documented = set(re.findall(r'"(\d+GB)"', gpu_line))
    assert documented, "the Resources() docstring no longer lists GPU memories"
    assert documented <= set(SlurmBackend.MEMORY_SPECS)


def test_the_memory_tokens_are_the_scheduler_memory_specs():
    # Two copies of one list drift; resource_requirements now reads the scheduler's.
    from biopipelines.resource_requirements import MEMORY_TOKENS
    from biopipelines.schedulers import SchedulerBackend
    assert MEMORY_TOKENS is SchedulerBackend.MEMORY_SPECS


def test_slurm_gpu_spec_none_still_declines_and_a_model_is_still_typed():
    from biopipelines.schedulers import SlurmBackend
    b = SlurmBackend()
    assert b.gpu_directive(None, 1)[0] == ""
    assert b.gpu_directive("none", 1)[0] == ""
    assert b.gpu_directive("A100", 2)[0] == "#SBATCH --gpus=A100:2"


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


@pytest.mark.parametrize("spec", [2, "2"])
@pytest.mark.parametrize("name, expected", [
    ("slurm", "#SBATCH --gpus=2"),
    ("lsf", '#BSUB -gpu "num=2"'),
    ("pbs", "#PBS -l select=1:ngpus=2"),
])
def test_an_integer_gpu_request_is_a_count_on_every_backend(name, expected, spec):
    """`gpu=2` used to reach `.startswith` on LSF and PBS and raise inside pipeline generation."""
    assert get_backend(name).gpu_directive(spec, 1) == (expected, [])


@pytest.mark.parametrize("name", ["lsf", "pbs"])
def test_140gb_is_a_memory_not_a_model(name):
    directive, warns = get_backend(name).gpu_directive("140GB", 1)
    assert "140GB" not in directive
    assert warns and "140GB" in warns[0]


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


def test_slurm_pipeline_with_an_integer_gpu_generates_its_batch_script(
    slurm_local_config, isolated_cwd, new_slurm_pipeline,
):
    """`Resources(gpu=1)` once crashed inside generation: the run directory was written and no job submitted."""
    from biopipelines.pipeline import Resources
    from biopipelines.mock import Mock

    pipeline = new_slurm_pipeline("int_gpu")
    with pipeline:
        Resources(gpu=1)
        Mock(ids=["a"], streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
        pipeline.save()
        pipeline.generate_job_scripts()

    script = (Path(pipeline.folders["runtime"]) / "slurm.sh").read_text(encoding="utf-8")
    assert "#SBATCH --gpus=1\n" in script
    assert "--gpus=None" not in script
