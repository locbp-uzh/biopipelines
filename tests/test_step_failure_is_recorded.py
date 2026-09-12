# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""A failing step used to be invisible three times over: it wrote no marker, the job kept going, and the scheduler reported success.

Every step runs as `step.sh 2>&1 | tee log`, and a pipeline's exit status in bash is the *last* command's -- tee, which always succeeds. So a tool that exited non-zero left `sacct` saying COMPLETED. Worse, a tool that exits from inside its own emitted bash (several do, on an upstream tool's failure) never reaches the block that writes COMPLETED or FAILED, so the step stayed unmarked -- and unmarked reads as "not run yet" everywhere, including the run page.

Measured on a real cluster run: CABSflex aborted on its first structure, wrote neither marker, never attempted the second structure, and SLURM reported COMPLETED. These tests hold the generated bash to reading `${PIPESTATUS[0]}` rather than tee's status, to writing the FAILED marker itself when the tool did not, and to failing the job at the end.
"""

import os
import re
import shutil

import pytest


def _script(new_pipeline, name, n_tools=2, **pipeline_kwargs):
    from biopipelines.mock import Mock

    pipeline = new_pipeline(name, **pipeline_kwargs)
    with pipeline:
        prev = None
        for i in range(n_tools):
            kwargs = {"source": prev.streams.structures} if prev else {"ids": ["a"]}
            prev = Mock(streams={"structures": {"format": "pdb", "file": "<id>.pdb"}}, **kwargs)
        path = pipeline.save()
    return open(path, encoding="utf-8").read()


def test_the_guard_reads_the_step_status_not_tees(local_config, isolated_cwd, new_pipeline):
    script = _script(new_pipeline, "guard_pipestatus")
    # A step script is `<order>_<Tool>.sh`; config.sh is also piped through tee and is not a step.
    tees = [ln for ln in script.splitlines() if "| tee " in ln and re.search(r"\d+_\w+\.sh ", ln)]
    assert tees, "no step invocations found"
    assert script.count("${PIPESTATUS[0]}") >= len(tees), (
        "every step invocation needs its own PIPESTATUS check, or tee's zero exit status hides the failure"
    )


def test_a_failing_step_gets_a_failed_marker_written_for_it(local_config, isolated_cwd, new_pipeline):
    """The tool may exit before its own completion check runs, so the pipeline writes the marker when the tool did not -- and never overwrites a COMPLETED one it did write."""
    script = _script(new_pipeline, "guard_marker")
    assert re.search(r'\[ -f "[^"]*_COMPLETED" \] \|\| touch "[^"]*_FAILED"', script), script[-2000:]


def test_the_job_fails_when_any_step_failed(local_config, isolated_cwd, new_pipeline):
    script = _script(new_pipeline, "guard_exit")
    assert "BP_FAILED_STEPS=0" in script
    assert "BP_FAILED_STEPS=$((BP_FAILED_STEPS + 1))" in script
    assert '"${BP_FAILED_STEPS:-0}" -ne 0' in script
    assert "exit 1" in script


def test_a_failed_marker_alone_also_fails_the_job(local_config, isolated_cwd, new_pipeline):
    """A tool can run to the end, find its outputs missing, write its own FAILED marker and still exit 0 -- measured on AlphaFold. The PIPESTATUS guard cannot see that one, so the markers are counted too."""
    script = _script(new_pipeline, "guard_marker_count")
    assert "BP_MY_MARKERS=(" in script
    assert re.search(r'\[ -f "\$_m" \] && BP_FAILED_MARKERS=\$\(\(BP_FAILED_MARKERS \+ 1\)\)', script)
    assert '"${BP_FAILED_MARKERS:-0}" -ne 0' in script


def test_the_marker_count_reaches_a_step_inside_a_folder(local_config, isolated_cwd, new_pipeline):
    """The count was a flat glob of the output root, but the marker goes beside the step folder: a step inside Folder(...) writes `<root>/<folder>/001_Mock_FAILED`, one level down, and an internal step writes under `.internal/`. So the count was structurally zero for both, and a grouped step that wrote FAILED and exited 0 still let afterok release every dependent batch.

    Asserting on the marker path the pipeline itself computes, rather than on a hardcoded depth, is what makes this hold if the layout moves again.
    """
    from biopipelines.mock import Mock
    from biopipelines.pipeline import Folder

    pipeline = new_pipeline("guard_grouped")
    with pipeline:
        with Folder("grouped"):
            tool = Mock(ids=["a"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
        path = pipeline.save()
    script = open(path, encoding="utf-8").read()

    marker = pipeline._completion_marker_paths(pipeline.tools[-1])["failed"]
    root = pipeline.folders["output"]
    assert marker.startswith(root) and "grouped" in marker[len(root):], (
        f"expected the marker below the output root, got {marker!r} under {root!r}"
    )
    # The script names the marker path itself, so no glob depth can miss it.
    assert marker in script, f"the grouped step's marker is not in the counted set: {marker}"


def test_the_page_is_refreshed_before_the_job_gives_up(local_config, isolated_cwd, new_pipeline):
    """A failed run is exactly when its page is worth having, so the exit comes after the refresh."""
    script = _script(new_pipeline, "guard_order")
    assert script.index("regenerate_pipeline_page") < script.index('"${BP_FAILED_STEPS:-0}"')


def test_every_step_is_guarded_not_just_the_last(local_config, isolated_cwd, new_pipeline):
    script = _script(new_pipeline, "guard_all", n_tools=4)
    assert script.count("${PIPESTATUS[0]}") >= 4
    assert script.count("touch ") >= 4


def test_sibling_batches_do_not_count_each_others_markers(local_config, isolated_cwd, new_pipeline):
    """Batches under `Parallel()` share one afterok parent, so they run at the same time. Counting `*_FAILED` across the whole output root made every sibling see every other sibling's failure: one bad batch failed its healthy siblings too, and each of those cancelled its own afterok successors.

    Each batch script now names its own steps' marker paths, so the count is scoped by construction rather than by where the files happen to sit.
    """
    from biopipelines.mock import Mock
    from biopipelines.pipeline import Parallel, Resources

    pipeline = new_pipeline("guard_sibling_scope")
    with pipeline:
        with Parallel():
            Resources(time="1:00:00")
            first = Mock(ids=["a"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
            Resources(time="1:00:00")
            second = Mock(source=first.streams.structures,
                          streams={"designs": {"format": "pdb", "file": "<id>.pdb"}})
        pipeline.save()
        pipeline.generate_job_scripts()

    runtime = pipeline.folders["runtime"]
    batch_scripts = sorted(f for f in os.listdir(runtime)
                           if f.startswith("slurm_batch") and f.endswith(".sh"))
    assert len(batch_scripts) >= 2, batch_scripts

    # Keyed by index, not TOOL_NAME: both steps here are Mock.
    markers = [pipeline._completion_marker_paths(t)["failed"] for t in pipeline.tools]
    assert len(set(markers)) == len(pipeline.tools), f"steps must have distinct markers: {markers}"

    for name in batch_scripts:
        body = open(os.path.join(runtime, name), encoding="utf-8").read()
        if "BP_MY_MARKERS=(" not in body:
            continue
        counted = [m for m in markers if m in body]
        assert len(counted) == 1, (
            f"{name} counts {len(counted)} markers; each batch must count only its own step"
        )


def _run_trap_probe(pipeline, name):
    """Run the kill trap in real bash and deliver SIGTERM from inside the script.

    Windows has no POSIX signals: `Popen.send_signal(SIGTERM)` becomes TerminateProcess, which kills bash without the trap ever running. Signalling from inside is what the script sees under SLURM anyway -- the trap body and the marker-pair expansions are the part that can be wrong, and this exercises both on any host with bash.
    """
    import subprocess

    trap = "\n".join(pipeline._kill_trap_lines(pipeline.tools))
    assert "trap _bp_on_kill TERM XCPU" in trap
    # Drop the page refresh: this is about the marker, and the refresh needs a full run on disk.
    body = "\n".join(ln for ln in trap.splitlines() if "python -c" not in ln)

    script = os.path.join(pipeline.folders["runtime"], f"_{name}.sh")
    with open(script, "w", newline="\n") as f:
        f.write("#!/bin/bash\n" + body + "\nkill -TERM $$\nsleep 5\n")

    return subprocess.run(["bash", script], capture_output=True, timeout=60)


@pytest.mark.skipif(not shutil.which("bash"), reason="needs bash to execute the trap")
def test_a_scheduler_kill_leaves_a_marker_behind(local_config, isolated_cwd, new_pipeline):
    """SLURM TIMEOUT and `scancel` kill the job mid-step, so the PIPESTATUS guard never runs and the step writes nothing. No marker reads as "not run yet" everywhere, including the page -- a dead step and a queued one looked identical, which is the defect the marker work exists to close.

    SIGTERM precedes SIGKILL on both, so trapping it is enough to leave a record. An OOM kill is SIGKILL and cannot be trapped; that one is still only visible in `sacct`.

    This runs the trap rather than grepping for it: a trap that is syntactically present but never fires records nothing, and the parameter expansions that split the marker pairs are exactly the kind of bash that looks right and is not.
    """
    from biopipelines.mock import Mock

    pipeline = new_pipeline("guard_kill_trap")
    with pipeline:
        Mock(ids=["a"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
        pipeline.save()

    markers = pipeline._completion_marker_paths(pipeline.tools[-1])
    os.makedirs(os.path.dirname(markers["failed"]), exist_ok=True)
    assert not os.path.exists(markers["failed"])

    proc = _run_trap_probe(pipeline, "trap_probe")
    out = (proc.stdout + proc.stderr).decode("utf-8", "replace")

    assert os.path.isfile(markers["failed"]), f"SIGTERM left no FAILED marker; script said: {out}"
    assert proc.returncode != 0, "a killed job must not report success"
    assert "timeout or scancel" in out


@pytest.mark.skipif(not shutil.which("bash"), reason="needs bash to execute the trap")
def test_the_kill_trap_does_not_overwrite_a_finished_step(local_config, isolated_cwd, new_pipeline):
    """A step that already finished keeps its verdict: the trap marks only what has neither marker, so a late kill cannot turn a completed run into a failed one."""
    from biopipelines.mock import Mock

    pipeline = new_pipeline("guard_kill_keeps_completed")
    with pipeline:
        Mock(ids=["a"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
        pipeline.save()

    markers = pipeline._completion_marker_paths(pipeline.tools[-1])
    os.makedirs(os.path.dirname(markers["completed"]), exist_ok=True)
    open(markers["completed"], "w").close()

    _run_trap_probe(pipeline, "trap_probe_completed")

    assert not os.path.exists(markers["failed"]), "a completed step must keep its verdict"
