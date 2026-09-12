"""Tests for ``Parallel(pack=N)`` in a MULTI-batch pipeline.

Every case in ``test_parallel_pack.py`` declares one ``Resources()`` and so
takes the single-batch path. A second ``Resources()`` routes generation
through ``_generate_multi_batch`` instead, which reads ``_packed_batches``
through different code — and got the container shape wrong, so every
multi-batch packed pipeline raised before a script was written.

Also covers the marker a packed task's later steps get when an earlier one
in the same task fails: exiting the subshell skips them, and a step holding
no marker at all reads as "still queued" rather than "never ran".
"""
from __future__ import annotations

from pathlib import Path

import pytest

STREAM = {"out": {"format": "pdb", "file": "<id>.pdb"}}


def _batch_script(pipeline, n):
    """The SBATCH wrapper; the body it invokes lives in ``pipeline_batch<N>.sh``."""
    return (Path(pipeline.folders["runtime"]) / f"slurm_batch{n}.sh").read_text(encoding="utf-8")


def _batch_body(pipeline, n):
    return (Path(pipeline.folders["runtime"]) / f"pipeline_batch{n}.sh").read_text(encoding="utf-8")


def _build(pipeline, packed_first: bool = False):
    """Two ``Resources()`` blocks, one of them a packed ``Parallel``."""
    from biopipelines.pipeline import Resources, Parallel, Run
    from biopipelines.mock import Mock

    with pipeline:
        if packed_first:
            with Parallel(pack=2):
                Resources(gpu="gh", time="06:00:00")
                for i in range(3):
                    with Run():
                        Mock(ids=[f"p{i}"], streams=STREAM)
            Resources(time="01:00:00")
            Mock(ids=["tail"], streams=STREAM)
        else:
            Resources(time="01:00:00")
            Mock(ids=["head"], streams=STREAM)
            with Parallel(pack=2):
                Resources(gpu="gh", time="06:00:00")
                for i in range(3):
                    with Run():
                        Mock(ids=[f"p{i}"], streams=STREAM)
        pipeline.save()
        pipeline.generate_job_scripts()


def test_multi_batch_packed_pipeline_generates(
    slurm_packed_config, isolated_cwd, new_packed_pipeline,
):
    """A packed block beside a second Resources() writes its scripts.

    ``_packed_batches`` maps a batch index to ``{"pack": N, "runs": [...]}``;
    iterating that dict yields its keys, so building the batch's tool list
    raised ``TypeError: string indices must be integers``.
    """
    pipeline = new_packed_pipeline("pack_multi_batch")
    _build(pipeline)

    runtime = Path(pipeline.folders["runtime"])
    batches = sorted(runtime.glob("slurm_batch*.sh"))
    assert len(batches) == 2, [b.name for b in batches]

    packed = _batch_script(pipeline, 2)
    assert "#SBATCH --ntasks=3" in packed
    assert "#SBATCH --ntasks-per-node=2" in packed
    assert "_pack_pids=()" in _batch_body(pipeline, 2)


def test_multi_batch_packed_pipeline_generates_when_packed_block_is_first(
    slurm_packed_config, isolated_cwd, new_packed_pipeline,
):
    """Batch order must not matter: the packed block can be batch 1."""
    pipeline = new_packed_pipeline("pack_multi_batch_first")
    _build(pipeline, packed_first=True)

    assert "#SBATCH --ntasks=3" in _batch_script(pipeline, 1)
    assert "_pack_pids=()" in _batch_body(pipeline, 1)
    # The unpacked tail batch carries no pack machinery.
    assert "_pack_pids=()" not in _batch_body(pipeline, 2)


def test_multi_batch_packed_batch_traps_kill_for_its_packed_tools(
    slurm_packed_config, isolated_cwd, new_packed_pipeline,
):
    """The kill trap names the packed tasks' steps, not the batch slice.

    The trap's marker list is what the broken comprehension fed, so a
    regression there is silent unless the emitted names are asserted.
    """
    pipeline = new_packed_pipeline("pack_multi_batch_trap")
    _build(pipeline)

    packed = _batch_body(pipeline, 2)
    assert "trap _bp_on_kill TERM XCPU" in packed
    # The trap's marker list is built from the same comprehension the bug broke.
    assert packed.count("_COMPLETED") >= 3


def test_packed_task_marks_its_later_steps_failed(
    slurm_packed_config, isolated_cwd, new_packed_pipeline,
):
    """A step that fails inside a packed Run() marks the rest of its task.

    The guard exits the subshell so ``wait`` sees the failure, which means
    every later step in that task never runs and never writes a marker --
    and an unmarked step reads as "still queued", not as "never ran".
    """
    from biopipelines.pipeline import Resources, Parallel, Run
    from biopipelines.mock import Mock

    pipeline = new_packed_pipeline("pack_successor_markers")
    with pipeline:
        with Parallel(pack=2):
            Resources(gpu="gh", time="06:00:00")
            with Run():
                Mock(ids=["a"], streams=STREAM)
                Mock(ids=["b"], streams=STREAM)
                Mock(ids=["c"], streams=STREAM)
        pipeline.save()
        pipeline.generate_job_scripts()

    script = (Path(pipeline.folders["runtime"]) / "pipeline.sh").read_text(encoding="utf-8")
    first, second, third = pipeline.tools
    markers = {t: pipeline._completion_marker_paths(t) for t in (first, second, third)}

    # The guard is everything between this step's invocation and the next one's.
    guard = script.split(first.script_basename + ".sh", 1)[1].split(second.script_basename + ".sh", 1)[0]

    for tool in (first, second, third):
        assert markers[tool]["failed"] in guard, f"{tool.script_basename} unmarked by the first step's guard"
    assert "exit 1" in guard

    # Marking is guarded on the step not having already completed.
    for tool in (second, third):
        assert f'[ -f "{markers[tool]["completed"]}" ] || touch' in guard

    # The last step has no successors, so its guard marks only itself. Bound the
    # slice at the subshell close, or the kill trap's MARKERS array is caught too.
    tail = script.split(third.script_basename + ".sh", 1)[1].split(") &", 1)[0]
    assert markers[third]["failed"] in tail
    assert markers[first]["failed"] not in tail
    assert markers[second]["failed"] not in tail
