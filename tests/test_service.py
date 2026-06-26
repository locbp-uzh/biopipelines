"""Tests for ``Service()`` and the ``after:`` (started, not afterok) edge.

A Service block launches a single long-running daemon batch (e.g. an MSA
server). The batch keeps its normal upstream ``afterok`` dependency, but the
first batch AFTER the block gets a SLURM ``after:`` dependency on it — it
starts once the daemon is RUNNING rather than after it finishes. Nothing
``afterok``-waits on the daemon, so the pipeline never blocks on its exit.

Covers:

* The after: edge lands on the batch following the block; the daemon keeps its
  afterok chain parent.
* The daemon is absent from every downstream afterok dependency.
* A Service block with zero or >1 batch raises.
* Nested Service(), Service inside Parallel, and Dependencies() inside Service
  all raise.

Like test_parallel.py, these use the SLURM-flavoured config so the framework
emits slurm_batch*.sh whose dependency lines we inspect directly.
"""
from __future__ import annotations

from pathlib import Path

import pytest


def test_service_block_emits_after_edge_on_next_batch(
    slurm_local_config, isolated_cwd, new_slurm_pipeline,
):
    """A -> Service:B -> C  =>  B is afterok:A, C is after:B."""
    from biopipelines.pipeline import Resources, Service
    from biopipelines.mock import Mock

    pipeline = new_slurm_pipeline("service_basic")
    with pipeline:
        Resources()
        Mock(ids=["a"], streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
        with Service():
            Resources()
            Mock(ids=["srv"], streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
        Resources()
        Mock(ids=["c"], streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
        pipeline.save()
        pipeline.generate_job_scripts()

    # Batches: 0 = A, 1 = daemon B, 2 = C.
    # afterok parents: A=[], B=[0], C=[]  (C's only predecessor is via after:)
    assert pipeline.batch_parents == [[], [0], []]
    assert pipeline.batch_after_parents == [[], [], [1]]

    runtime = Path(pipeline.folders["runtime"])
    batch2 = (runtime / "slurm_batch2.sh").read_text(encoding="utf-8")  # daemon
    batch3 = (runtime / "slurm_batch3.sh").read_text(encoding="utf-8")  # consumer
    assert "afterok:<JOBID_BATCH_001>" in batch2          # daemon afterok:A
    assert "after:<JOBID_BATCH_002>" in batch3            # consumer after:daemon
    assert "afterok:" not in batch3                       # never afterok the daemon


def test_service_daemon_absent_from_downstream_afterok(
    slurm_local_config, isolated_cwd, new_slurm_pipeline,
):
    """A step after the consumer afterok-chains off the consumer, never the daemon."""
    from biopipelines.pipeline import Resources, Service
    from biopipelines.mock import Mock

    pipeline = new_slurm_pipeline("service_chain")
    with pipeline:
        Resources()
        Mock(ids=["a"], streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
        with Service():
            Resources()
            Mock(ids=["srv"], streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
        Resources()
        Mock(ids=["c"], streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
        Resources()
        Mock(ids=["d"], streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
        pipeline.save()
        pipeline.generate_job_scripts()

    # D (batch 3) chains off C (batch 2) via afterok; daemon (1) is never a parent.
    assert pipeline.batch_parents == [[], [0], [], [2]]
    assert pipeline.batch_after_parents == [[], [], [1], []]

    runtime = Path(pipeline.folders["runtime"])
    batch4 = (runtime / "slurm_batch4.sh").read_text(encoding="utf-8")
    assert "afterok:<JOBID_BATCH_003>" in batch4
    assert "<JOBID_BATCH_002>" not in batch4   # daemon never referenced downstream


def test_service_empty_block_raises(
    slurm_local_config, isolated_cwd, new_slurm_pipeline,
):
    from biopipelines.pipeline import Resources, Service

    pipeline = new_slurm_pipeline("service_empty")
    with pytest.raises(RuntimeError, match="exactly one batch"):
        with pipeline:
            Resources()
            with Service():
                pass


def test_service_two_batches_raises(
    slurm_local_config, isolated_cwd, new_slurm_pipeline,
):
    from biopipelines.pipeline import Resources, Service
    from biopipelines.mock import Mock

    pipeline = new_slurm_pipeline("service_two")
    with pytest.raises(RuntimeError, match="exactly one batch"):
        with pipeline:
            Resources()
            with Service():
                Resources()
                Mock(ids=["s1"], streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
                Resources()
                Mock(ids=["s2"], streams={"out": {"format": "pdb", "file": "<id>.pdb"}})


def test_nested_service_raises(
    slurm_local_config, isolated_cwd, new_slurm_pipeline,
):
    from biopipelines.pipeline import Resources, Service

    pipeline = new_slurm_pipeline("service_nested")
    with pytest.raises(RuntimeError, match="Nested"):
        with pipeline:
            Resources()
            with Service():
                with Service():
                    pass


def test_service_inside_parallel_raises(
    slurm_local_config, isolated_cwd, new_slurm_pipeline,
):
    from biopipelines.pipeline import Resources, Parallel, Service

    pipeline = new_slurm_pipeline("service_in_parallel")
    with pytest.raises(RuntimeError, match="Parallel"):
        with pipeline:
            Resources()
            with Parallel():
                with Service():
                    pass


def test_dependencies_inside_service_raises(
    slurm_local_config, isolated_cwd, new_slurm_pipeline,
):
    from biopipelines.pipeline import Resources, Service, Dependencies

    pipeline = new_slurm_pipeline("deps_in_service")
    with pytest.raises(RuntimeError, match="Service"):
        with pipeline:
            Resources()
            with Service():
                Dependencies("12345")
