"""Tests for ``Parallel(pack=N)`` / ``Run()`` — node packing via job steps.

On a machine that allocates whole nodes, N sibling jobs waste the rest of
each node. A packed block emits ONE job whose tasks are concurrent job
steps, N per node. Covers the header geometry, the per-tool step prefix,
per-task CPU shares, and the misuse errors.

Also covers the EDF prefix in the multi-batch path, which previously
applied only to single-batch pipelines.
"""
from __future__ import annotations

from pathlib import Path

import pytest

STREAM = {"out": {"format": "pdb", "file": "<id>.pdb"}}


def _batch_files(pipeline):
    runtime = Path(pipeline.folders["runtime"])
    return runtime, (runtime / "slurm.sh").read_text(encoding="utf-8")


# ── header geometry ───────────────────────────────────────────────────────────

def test_pack_derives_node_count_from_tasks_per_node(
    slurm_packed_config, isolated_cwd, new_packed_pipeline,
):
    """10 tasks at pack=4 spread 4/4/2 over 3 nodes."""
    from biopipelines.pipeline import Resources, Parallel, Run
    from biopipelines.mock import Mock

    pipeline = new_packed_pipeline("pack_geometry")
    with pipeline:
        with Parallel(pack=4):
            Resources(gpu="gh", time="06:00:00")
            for i in range(10):
                with Run():
                    Mock(ids=[f"s{i}"], streams=STREAM)
        pipeline.save()
        pipeline.generate_job_scripts()

    _, batch1 = _batch_files(pipeline)
    assert "#SBATCH --nodes=3" in batch1
    assert "#SBATCH --ntasks=10" in batch1
    assert "#SBATCH --ntasks-per-node=4" in batch1
    # node_exclusive fixture: the node's GPUs are requested once, not per task,
    # so a CPU-only task alongside GPU ones cannot overcommit them.
    assert "#SBATCH --gres=gpu:4" in batch1
    assert "--gpus-per-task" not in batch1.split("--distribution")[0]
    assert "#SBATCH --distribution=block" in batch1


def test_pack_cpu_task_alongside_gpu_tasks_does_not_overcommit(
    slurm_packed_config, isolated_cwd, new_packed_pipeline,
):
    """A CPU-only task plus 4 GPU tasks fits a 4-GPU node.

    --gpus-per-task multiplies by task count, so 5 tasks would ask for 5 GPUs
    and SLURM rejects the allocation outright ("Requested node configuration is
    not available"). The whole node is allocated anyway when it is exclusive.
    """
    from biopipelines.pipeline import Resources, Parallel, Run
    from biopipelines.mock import Mock

    pipeline = new_packed_pipeline("pack_mixed_gpus")
    with pipeline:
        with Parallel(pack=5):
            Resources(gpu="gh", time="06:00:00")
            with Run(cpus=160, gpus=0, memory="720GB"):
                Mock(ids=["server"], streams=STREAM)
            for i in range(4):
                with Run(cpus=32, gpus=1, memory="20GB"):
                    Mock(ids=[f"g{i}"], streams=STREAM)
        pipeline.save()
        pipeline.generate_job_scripts()

    runtime, batch1 = _batch_files(pipeline)
    assert "#SBATCH --nodes=1" in batch1
    assert "#SBATCH --ntasks=5" in batch1
    assert "#SBATCH --gres=gpu:4" in batch1
    assert "--gpus-per-task=5" not in batch1
    # The per-step assignment still distinguishes the CPU task from the GPU ones.
    body = (runtime / "pipeline.sh").read_text(encoding="utf-8")
    assert "--gres=none" in body
    assert "--gpus-per-task=0" not in body
    assert "--mem=720000M" in body
    assert body.count("--mem=20000M") == 4
    assert "--mem-per-cpu=" not in body
    assert "--gpus-per-task=1" in body
    # A batch-level --gpus total would let SLURM pick a different spread.
    assert "#SBATCH --gpus=" not in batch1


def test_pack_fewer_tasks_than_pack_uses_one_node(
    slurm_packed_config, isolated_cwd, new_packed_pipeline,
):
    from biopipelines.pipeline import Resources, Parallel, Run
    from biopipelines.mock import Mock

    pipeline = new_packed_pipeline("pack_small")
    with pipeline:
        with Parallel(pack=4):
            Resources()
            for i in range(2):
                with Run():
                    Mock(ids=[f"s{i}"], streams=STREAM)
        pipeline.save()
        pipeline.generate_job_scripts()

    _, batch1 = _batch_files(pipeline)
    assert "#SBATCH --nodes=1" in batch1
    assert "#SBATCH --ntasks=2" in batch1
    assert "#SBATCH --ntasks-per-node=2" in batch1


# ── body: one step per tool ───────────────────────────────────────────────────

def test_each_tool_is_its_own_step_with_explicit_cpus(
    slurm_packed_config, isolated_cwd, new_packed_pipeline,
):
    """Every tool gets an srun step; the CPU share divides the node by pack.

    A step with no --cpus-per-task silently gets ONE core, which runs
    correctly and very slowly, so the flag must always be present.
    """
    from biopipelines.pipeline import Resources, Parallel, Run
    from biopipelines.mock import Mock

    pipeline = new_packed_pipeline("pack_steps")
    with pipeline:
        with Parallel(pack=4):
            Resources()
            for i in range(2):
                with Run():
                    Mock(ids=[f"s{i}"], streams=STREAM)
        pipeline.save()
        pipeline.generate_job_scripts()

    body = (Path(pipeline.folders["runtime"]) / "pipeline.sh").read_text(encoding="utf-8")
    # 288 cores / pack 4 = 72
    assert body.count("srun --exclusive --exact --nodes=1 --ntasks=1 --gpus-per-task=1 --cpus-per-task=72") == 2
    assert body.count(") &") == 2
    assert "wait" in body


def test_failed_packed_tasks_make_the_batch_fail(
    slurm_packed_config, isolated_cwd, new_packed_pipeline,
):
    """After waiting for every sibling, any failed task must make the packed
    script exit non-zero so SLURM and afterok dependencies see the failure."""
    from biopipelines.pipeline import Resources, Parallel, Run
    from biopipelines.mock import Mock

    pipeline = new_packed_pipeline("pack_failure_status")
    with pipeline:
        with Parallel(pack=4):
            Resources()
            with Run():
                Mock(ids=["a"], streams=STREAM)
        pipeline.save()
        pipeline.generate_job_scripts()

    body = (Path(pipeline.folders["runtime"]) / "pipeline.sh").read_text(encoding="utf-8")
    failure_block = body.split('if [ "$_pack_failed" -gt 0 ]; then', 1)[1].split("fi", 1)[0]
    assert 'echo "$_pack_failed of ${#_pack_pids[@]} packed tasks failed"' in failure_block
    assert "exit 1" in failure_block


def test_run_cpus_overrides_derived_share(
    slurm_packed_config, isolated_cwd, new_packed_pipeline,
):
    from biopipelines.pipeline import Resources, Parallel, Run
    from biopipelines.mock import Mock

    pipeline = new_packed_pipeline("pack_cpus")
    with pipeline:
        with Parallel(pack=4):
            Resources()
            with Run(cpus=16):
                Mock(ids=["a"], streams=STREAM)
        pipeline.save()
        pipeline.generate_job_scripts()

    body = (Path(pipeline.folders["runtime"]) / "pipeline.sh").read_text(encoding="utf-8")
    assert "--cpus-per-task=16" in body
    assert "--cpus-per-task=72" not in body


def test_run_may_mix_containerized_and_plain_tools(
    slurm_packed_config, isolated_cwd, new_packed_pipeline,
):
    """Per-tool EDF selection survives packing: the container flag rides on
    each tool's own step, since nesting a container step inside a task step
    fails to bind CPUs."""
    from biopipelines.pipeline import Resources, Parallel, Run
    from biopipelines.mock import Mock
    from biopipelines.panda import Panda

    pipeline = new_packed_pipeline("pack_mixed")
    with pipeline:
        with Parallel(pack=4):
            Resources()
            with Run():
                m = Mock(ids=["a"], streams=STREAM, tables={"info": {"columns": ["id"]}})
                Panda(                                # no edf: entry
                    tables=m.tables.info,
                    operations=[Panda.sort("id", ascending=True)],
                )
        pipeline.save()
        pipeline.generate_job_scripts()

    body = (Path(pipeline.folders["runtime"]) / "pipeline.sh").read_text(encoding="utf-8")
    assert "--environment=/fixtures/edf/ngc.toml" in body
    # Both tools are steps; only one carries the container flag.
    assert body.count("srun --exclusive") == 2
    assert body.count("--environment=") == 1


# ── misuse ────────────────────────────────────────────────────────────────────

def test_run_outside_packed_block_raises(
    slurm_packed_config, isolated_cwd, new_packed_pipeline,
):
    from biopipelines.pipeline import Resources, Run
    from biopipelines.mock import Mock

    pipeline = new_packed_pipeline("run_alone")
    with pipeline:
        Resources()
        Mock(ids=["a"], streams=STREAM)
        with pytest.raises(RuntimeError, match="only valid inside"):
            with Run():
                pass
        pipeline.save()


def test_second_resources_in_packed_block_raises(
    slurm_packed_config, isolated_cwd, new_packed_pipeline,
):
    """Resources() describes the shared allocation, so it is called once."""
    from biopipelines.pipeline import Resources, Parallel, Run
    from biopipelines.mock import Mock

    pipeline = new_packed_pipeline("pack_two_resources")
    with pipeline:
        with pytest.raises(RuntimeError, match="single Resources"):
            with Parallel(pack=4):
                Resources()
                with Run():
                    Mock(ids=["a"], streams=STREAM)
                Resources()


def test_pack_requires_positive_int():
    from biopipelines.pipeline import Parallel

    for bad in (0, -1, 2.5, True, "4"):
        with pytest.raises(ValueError, match="positive integer"):
            Parallel(pack=bad)


def test_run_rejects_invalid_cpus_and_gpus():
    """Both reach an srun directive verbatim, so a bad value would otherwise
    surface only after the allocation is granted. bool is an int subclass, so
    ``gpus=True`` would format as ``--gpus-per-task=True``."""
    from biopipelines.pipeline import Run

    for bad in (0, -1, 2.5, True, "8"):
        with pytest.raises(ValueError, match="positive integer"):
            Run(cpus=bad)
    for bad in (-1, 2.5, True, "1"):
        with pytest.raises(ValueError, match="non-negative integer"):
            Run(gpus=bad)
    # gpus=0 is how a CPU-only task declines GRES — it must stay legal.
    Run(cpus=32, gpus=0)
    Run()


# ── the node-geometry contract ────────────────────────────────────────────────

def test_pack_refuses_a_non_exclusive_machine(
    slurm_local_config, isolated_cwd, new_slurm_pipeline,
):
    """``get_node_exclusive`` documents that pack refuses to engage without it.

    The slurm_local fixture declares no node geometry, so packing there would
    emit steps with no --cpus-per-task at all — one core each, silently.
    """
    from biopipelines.pipeline import Resources, Parallel, Run
    from biopipelines.mock import Mock

    pipeline = new_slurm_pipeline("pack_non_exclusive")
    with pipeline:
        Resources()
        Mock(ids=["pre"], streams=STREAM)
        with pytest.raises(RuntimeError, match="node_exclusive"):
            with Parallel(pack=4):
                Resources()
                with Run():
                    Mock(ids=["a"], streams=STREAM)


def test_pack_requires_cores_per_node(
    slurm_packed_config, isolated_cwd, new_packed_pipeline, monkeypatch,
):
    """Without it there is nothing to divide between the tasks."""
    from biopipelines.config_manager import ConfigManager
    from biopipelines.pipeline import Resources, Parallel, Run
    from biopipelines.mock import Mock

    monkeypatch.setattr(ConfigManager, "get_cores_per_node", lambda self: None)

    pipeline = new_packed_pipeline("pack_no_cores")
    with pipeline:
        Resources()
        Mock(ids=["pre"], streams=STREAM)
        with pytest.raises(RuntimeError, match="cores_per_node"):
            with Parallel(pack=4):
                Resources()
                with Run():
                    Mock(ids=["a"], streams=STREAM)


def test_pack_with_gpu_tasks_requires_gpus_per_node(
    slurm_packed_config, isolated_cwd, new_packed_pipeline, monkeypatch,
):
    """On an exclusive node the GPUs are requested once, so the count must be known."""
    from biopipelines.config_manager import ConfigManager
    from biopipelines.pipeline import Resources, Parallel, Run
    from biopipelines.mock import Mock

    monkeypatch.setattr(ConfigManager, "get_gpus_per_node", lambda self: None)

    pipeline = new_packed_pipeline("pack_no_gpus_per_node")
    with pipeline:
        with pytest.raises(RuntimeError, match="gpus_per_node"):
            with Parallel(pack=4):
                Resources(gpu="gh")
                with Run(gpus=1):
                    Mock(ids=["a"], streams=STREAM)


def test_pack_without_gpu_tasks_ignores_gpus_per_node(
    slurm_packed_config, isolated_cwd, new_packed_pipeline, monkeypatch,
):
    """A CPU-only packed block never requests GRES, so the count is irrelevant."""
    from biopipelines.config_manager import ConfigManager
    from biopipelines.pipeline import Resources, Parallel, Run
    from biopipelines.mock import Mock

    monkeypatch.setattr(ConfigManager, "get_gpus_per_node", lambda self: None)

    pipeline = new_packed_pipeline("pack_cpu_only")
    with pipeline:
        with Parallel(pack=4):
            Resources()
            with Run(gpus=0):
                Mock(ids=["a"], streams=STREAM)
        pipeline.save()
        pipeline.generate_job_scripts()

    body = (Path(pipeline.folders["runtime"]) / "pipeline.sh").read_text(encoding="utf-8")
    assert "--gres=none" in body


def test_pack_refuses_a_non_slurm_scheduler(
    pbs_local_config, isolated_cwd,
):
    """Packing emits srun job steps, which only SLURM has."""
    from biopipelines.pipeline import Pipeline, Resources, Parallel, Run
    from biopipelines.mock import Mock

    pipeline = Pipeline(
        project="TestSuite", job="pack_pbs", on_the_fly=False,
        local_output=True, config="pbs_local",
    )
    with pipeline:
        Resources()
        Mock(ids=["pre"], streams=STREAM)
        with pytest.raises(RuntimeError, match="SLURM-specific"):
            with Parallel(pack=4):
                Resources()
                with Run():
                    Mock(ids=["a"], streams=STREAM)


def test_pack_refuses_scheduler_none_even_with_exclusive_geometry(
    local_config, isolated_cwd, monkeypatch,
):
    """The inert-script backend falls back to SLURM, but packing must inspect
    the configured scheduler rather than that script-generation fallback."""
    from biopipelines.config_manager import ConfigManager
    from biopipelines.pipeline import Pipeline, Resources, Parallel, Run
    from biopipelines.mock import Mock

    monkeypatch.setattr(ConfigManager, "get_node_exclusive", lambda self: True)
    monkeypatch.setattr(ConfigManager, "get_cores_per_node", lambda self: 64)
    monkeypatch.setattr(ConfigManager, "get_gpus_per_node", lambda self: 4)

    pipeline = Pipeline(
        project="TestSuite", job="pack_none", on_the_fly=False,
        local_output=True, config="local",
    )
    with pipeline:
        Resources()
        Mock(ids=["pre"], streams=STREAM)
        with pytest.raises(RuntimeError, match="scheduler is 'none'"):
            with Parallel(pack=4):
                Resources()
                with Run():
                    Mock(ids=["a"], streams=STREAM)


def test_unpacked_parallel_is_unchanged(
    slurm_packed_config, isolated_cwd, new_packed_pipeline,
):
    """Parallel() with no pack still opens one batch per sibling."""
    from biopipelines.pipeline import Resources, Parallel
    from biopipelines.mock import Mock

    pipeline = new_packed_pipeline("pack_absent")
    with pipeline:
        Resources()
        Mock(ids=["pre"], streams=STREAM)
        with Parallel():
            Resources()
            Mock(ids=["s1"], streams=STREAM)
            Resources()
            Mock(ids=["s2"], streams=STREAM)
        pipeline.save()
        pipeline.generate_job_scripts()

    assert pipeline.batch_parents == [[], [0], [0]]
    body = (Path(pipeline.folders["runtime"]) / "pipeline_batch2.sh").read_text(encoding="utf-8")
    assert "srun --exclusive" not in body


# ── EDF in the multi-batch path (regression) ──────────────────────────────────

def test_edf_prefix_applies_in_multi_batch_pipelines(
    slurm_packed_config, isolated_cwd, new_packed_pipeline,
):
    """The container prefix used to be emitted only for single-batch
    pipelines, so every EDF tool in a chained pipeline ran uncontainerized."""
    from biopipelines.pipeline import Resources
    from biopipelines.mock import Mock

    pipeline = new_packed_pipeline("edf_multibatch")
    with pipeline:
        Resources()
        Mock(ids=["a"], streams=STREAM)
        Resources()
        Mock(ids=["b"], streams=STREAM)
        pipeline.save()
        pipeline.generate_job_scripts()

    runtime = Path(pipeline.folders["runtime"])
    for n in (1, 2):
        body = (runtime / f"pipeline_batch{n}.sh").read_text(encoding="utf-8")
        assert "srun --environment=/fixtures/edf/ngc.toml" in body, f"batch {n} missing EDF"
