"""Tests for ``Parallel()`` and the underlying batch-DAG machinery.

Covers:

* The chain default (no Parallel block): batch N depends on batch N-1.
* Inside a Parallel block: each Resources() opens a sibling sharing the
  pre-block "anchor" parent.
* After a Parallel block: the next Resources() (or the next tool added
  without a Resources call) opens a fan-in batch depending on every
  sibling.
* Tools added as kwargs (e.g. Mock(structures=Mock(...))) land in the
  same iteration batch as the outer tool, not their own sibling.
* Empty Parallel block: no siblings opened; chain default resumes.
* Dependencies() inside Parallel: raises.
* Nested Parallel(): raises.

Tests use the SLURM-flavoured config so the framework emits
``slurm_batch*.sh`` artifacts whose dependency-line placeholder format we
inspect directly.
"""
from __future__ import annotations

import os
from pathlib import Path

import pytest


# ── chain-default regression (no Parallel) ────────────────────────────────────

def test_chain_default_two_batches_depends_on_one(
    slurm_local_config, isolated_cwd, new_slurm_pipeline,
):
    """Without Parallel, batch 2 should depend on batch 1 only."""
    from biopipelines.pipeline import Resources
    from biopipelines.mock import Mock

    pipeline = new_slurm_pipeline("chain_two")
    with pipeline:
        Resources()
        Mock(ids=["a"], streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
        Resources()
        Mock(ids=["b"], streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
        pipeline.save()
        # slurm_batch*.sh are emitted by .generate_job_scripts(), not .save().
        pipeline.generate_job_scripts()

    assert pipeline.batch_parents == [[], [0]]

    runtime = Path(pipeline.folders["runtime"])
    batch1 = (runtime / "slurm_batch1.sh").read_text(encoding="utf-8")
    batch2 = (runtime / "slurm_batch2.sh").read_text(encoding="utf-8")
    assert "<JOBID_BATCH_" not in batch1
    assert "afterok:<JOBID_BATCH_001>" in batch2
    assert "<JOBID_BATCH_002>" not in batch2  # no spurious self-dep


# ── single Parallel block ─────────────────────────────────────────────────────

def test_parallel_block_opens_sibling_batches(
    slurm_local_config, isolated_cwd, new_slurm_pipeline,
):
    """Two Resources() calls inside a Parallel block should each depend on
    the same anchor (batch 0), not on each other."""
    from biopipelines.pipeline import Resources, Parallel
    from biopipelines.mock import Mock

    pipeline = new_slurm_pipeline("parallel_block")
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
        # slurm_batch*.sh are emitted by .generate_job_scripts(), not .save().
        pipeline.generate_job_scripts()

    # Batches: 0 = pre, 1 = sibling 1, 2 = sibling 2, 3 = post.
    # Parents: [], [0], [0], [1, 2].
    assert pipeline.batch_parents == [[], [0], [0], [1, 2]]

    runtime = Path(pipeline.folders["runtime"])
    batch2 = (runtime / "slurm_batch2.sh").read_text(encoding="utf-8")
    batch3 = (runtime / "slurm_batch3.sh").read_text(encoding="utf-8")
    batch4 = (runtime / "slurm_batch4.sh").read_text(encoding="utf-8")
    assert "afterok:<JOBID_BATCH_001>" in batch2
    assert "afterok:<JOBID_BATCH_001>" in batch3
    assert "afterok:<JOBID_BATCH_002>:<JOBID_BATCH_003>" in batch4


def test_parallel_kwargs_land_in_iteration_batch(
    slurm_local_config, isolated_cwd, new_slurm_pipeline,
):
    """Entity tools used as kwargs (registered before the outer tool)
    must land in the same sibling batch as the outer tool, not open
    their own. The Resources() call at the top of the iteration is the
    one and only batch opener."""
    from biopipelines.pipeline import Resources, Parallel
    from biopipelines.mock import Mock

    pipeline = new_slurm_pipeline("parallel_kwargs")
    with pipeline:
        with Parallel():
            for label in ("a", "b"):
                Resources(memory="42GB", time="1:00:00")
                # Mock-as-kwarg registers first; the outer Mock should
                # join the same batch instead of opening a new sibling.
                inner = Mock(ids=[f"{label}_in"], streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
                Mock(ids=[label], structures=inner,
                     streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
        pipeline.save()

    # Two siblings, no anchor (Parallel sat at top of pipeline → empty parents).
    assert pipeline.batch_parents == [[], []]
    # Each sibling holds two Mocks (inner + outer).
    starts = pipeline.batch_start_indices
    assert starts == [0, 2]
    assert len(pipeline.tools) == 4


# ── empty Parallel block ──────────────────────────────────────────────────────

def test_empty_parallel_block_falls_back_to_chain(
    slurm_local_config, isolated_cwd, new_slurm_pipeline,
):
    """A Parallel block with no Resources() inside should not change the
    dependency graph: the next Resources() resumes the chain default."""
    from biopipelines.pipeline import Resources, Parallel
    from biopipelines.mock import Mock

    pipeline = new_slurm_pipeline("empty_parallel")
    with pipeline:
        Resources()
        Mock(ids=["a"], streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
        with Parallel():
            pass
        Resources()
        Mock(ids=["b"], streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
        pipeline.save()

    assert pipeline.batch_parents == [[], [0]]


# ── Dependencies inside Parallel raises ───────────────────────────────────────

def test_dependencies_inside_parallel_raises(
    slurm_local_config, isolated_cwd, new_slurm_pipeline,
):
    from biopipelines.pipeline import Resources, Parallel, Dependencies

    pipeline = new_slurm_pipeline("deps_in_parallel")
    with pytest.raises(RuntimeError, match="Parallel"):
        with pipeline:
            Resources()
            with Parallel():
                Dependencies("12345")


# ── nested Parallel raises ────────────────────────────────────────────────────

def test_nested_parallel_raises(
    slurm_local_config, isolated_cwd, new_slurm_pipeline,
):
    from biopipelines.pipeline import Resources, Parallel

    pipeline = new_slurm_pipeline("nested_parallel")
    with pytest.raises(RuntimeError, match="Nested"):
        with pipeline:
            Resources()
            with Parallel():
                with Parallel():
                    pass


# ── auto fan-in: tool after Parallel without Resources() ──────────────────────

def test_post_parallel_tool_without_resources_auto_opens_fan_in(
    slurm_local_config, isolated_cwd, new_slurm_pipeline,
):
    """If the user adds a tool right after a Parallel block without an
    explicit Resources() call, the fan-in batch is opened automatically
    (inheriting resources from the most recent sibling)."""
    from biopipelines.pipeline import Resources, Parallel
    from biopipelines.mock import Mock

    pipeline = new_slurm_pipeline("post_parallel_auto")
    with pipeline:
        Resources()
        Mock(ids=["pre"], streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
        with Parallel():
            Resources()
            Mock(ids=["s1"], streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
            Resources()
            Mock(ids=["s2"], streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
        # No explicit Resources() here — Mock should still produce a fan-in batch.
        Mock(ids=["post"], streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
        pipeline.save()

    assert pipeline.batch_parents == [[], [0], [0], [1, 2]]
