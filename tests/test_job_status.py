"""Reading a run's state off its completion markers.

The case that matters is the half-failed campaign: an agent that can only see stdout has to
grep `Logs/` to learn that step 4 of 9 failed. These tests pin the two claims the summary
makes about a failure — that the run stopped, or that it did not — because getting either
backwards sends the agent hunting for a cascade that never happened.
"""

import io

from biopipelines import job_status, run_log
from biopipelines.remote import LocalFS


def build(job, steps):
    """steps: [(name, marker or None)] — marker None leaves the step folder with no marker."""
    (job / "Logs").mkdir(parents=True, exist_ok=True)
    for name, marker in steps:
        (job / name).mkdir(exist_ok=True)
        io.open(job / "Logs" / f"{name}.log", "w", encoding="utf-8").write(f"log of {name}\n")
        if marker:
            io.open(job / f"{name}_{marker}", "w", encoding="utf-8").write("")
    return job


def test_steps_are_listed_in_execution_order(tmp_path):
    job = build(tmp_path / "run_001", [("002_Boltz2", "COMPLETED"),
                                       ("001_RFdiffusion", "COMPLETED"),
                                       ("010_PoseBusters", "COMPLETED")])
    assert [s["step"] for s in job_status.steps(job)] == [
        "001_RFdiffusion", "002_Boltz2", "010_PoseBusters"]


def test_a_tool_name_containing_underscores_is_parsed_whole(tmp_path):
    job = build(tmp_path / "run_001", [("003_PLM_Sol", "COMPLETED")])
    step = job_status.steps(job)[0]
    assert step["tool"] == "PLM_Sol" and step["index"] == "003"


def test_each_marker_maps_to_its_status(tmp_path):
    job = build(tmp_path / "run_001", [("001_A", "COMPLETED"), ("002_B", "FAILED"),
                                       ("003_C", "WARNING"), ("004_D", None)])
    assert {s["step"]: s["status"] for s in job_status.steps(job)} == {
        "001_A": "completed", "002_B": "failed", "003_C": "warning", "004_D": "pending"}


def test_a_missing_job_directory_is_reported_not_raised(tmp_path):
    state = job_status.status(tmp_path / "never_ran")
    assert state["exists"] is False and state["steps"] == []
    assert "No run directory" in job_status.summarize(state)


def test_summary_says_the_run_stopped_when_everything_after_is_pending(tmp_path):
    job = build(tmp_path / "run_001", [("001_A", "COMPLETED"), ("002_B", "FAILED"),
                                       ("003_C", None), ("004_D", None)])
    text = job_status.summarize(job_status.status(job))
    assert "First failure: 002_B" in text
    assert "the run stopped there" in text


def test_summary_says_the_pipeline_continued_when_later_steps_ran(tmp_path):
    """A failed step does not always halt the rest — the real mock_run_001 behaves this way."""
    job = build(tmp_path / "run_001", [("001_A", "COMPLETED"), ("002_B", "FAILED"),
                                       ("003_C", "COMPLETED")])
    text = job_status.summarize(job_status.status(job))
    assert "did not halt the pipeline" in text
    assert "stopped there" not in text


def test_counts_group_by_status(tmp_path):
    job = build(tmp_path / "run_001", [("001_A", "COMPLETED"), ("002_B", "COMPLETED"),
                                       ("003_C", "FAILED")])
    assert job_status.status(job)["counts"] == {"completed": 2, "failed": 1}


def test_the_run_header_is_surfaced_in_the_summary(tmp_path):
    job = build(tmp_path / "run_001", [("001_A", "COMPLETED")])
    run_log.header(job, pipeline="x.py")
    state = job_status.status(job)
    assert state["header"]["action"] == "run"
    assert "biopipelines=" in job_status.summarize(state)


def test_log_returns_the_tail_with_counts(tmp_path):
    job = build(tmp_path / "run_001", [("001_A", "COMPLETED")])
    io.open(job / "Logs" / "001_A.log", "w", encoding="utf-8").write(
        "\n".join(str(i) for i in range(500)))
    found = job_status.log(job, "001_A", tail=10)
    assert found["lines"] == 500 and found["shown"] == 10
    assert found["text"].splitlines()[-1] == "499"


def test_a_log_can_be_addressed_by_tool_name_alone(tmp_path):
    job = build(tmp_path / "run_001", [("007_Boltz2", "COMPLETED")])
    assert job_status.log(job, "Boltz2")["step"] == "007_Boltz2"


def test_an_ambiguous_tool_name_resolves_to_nothing(tmp_path):
    """Two Boltz2 steps: guessing which one the user meant would be worse than saying no."""
    job = build(tmp_path / "run_001", [("001_Boltz2", "COMPLETED"), ("005_Boltz2", "COMPLETED")])
    assert job_status.log(job, "Boltz2") is None
    assert job_status.log(job, "005_Boltz2")["step"] == "005_Boltz2"


def test_a_log_with_undecodable_bytes_still_reads(tmp_path):
    job = build(tmp_path / "run_001", [("001_A", "COMPLETED")])
    io.open(job / "Logs" / "001_A.log", "wb").write(b"fine\n\xff\xfe bad bytes\n")
    assert "fine" in job_status.log(job, "001_A")["text"]


class TestSuffixedSteps:
    """A `Suffix` lands on the output folder and the log, but not on the completion marker.

    Found live on a cluster run: one job has 36 markers like `001_PDB_COMPLETED` beside
    36 folders like `001_PDB_stitched_001/`. Treating both as steps reported 72 — 36 real and
    36 invented as "pending" — which is exactly the false signal that sends someone hunting for
    work that never existed.
    """

    def suffixed(self, tmp_path):
        job = tmp_path / "run_001"
        (job / "Logs").mkdir(parents=True)
        for index, tool, marker in (("001", "PDB", "COMPLETED"),
                                    ("002", "AlphaFold", "FAILED")):
            folder = f"{index}_{tool}_stitched_001"
            (job / folder).mkdir()
            io.open(job / f"{index}_{tool}_{marker}", "w", encoding="utf-8").write("")
            io.open(job / "Logs" / f"{folder}.log", "w", encoding="utf-8").write("out\n")
        return job

    def test_a_suffixed_folder_does_not_become_a_second_step(self, tmp_path):
        steps = job_status.steps(self.suffixed(tmp_path))
        assert [s["step"] for s in steps] == ["001_PDB", "002_AlphaFold"]
        assert [s["status"] for s in steps] == ["completed", "failed"]

    def test_the_suffixed_folder_is_recorded_against_its_step(self, tmp_path):
        assert job_status.steps(self.suffixed(tmp_path))[0]["folder"] == "001_PDB_stitched_001"

    def test_the_suffixed_log_is_found_by_the_unsuffixed_step_name(self, tmp_path):
        found = job_status.log(self.suffixed(tmp_path), "002_AlphaFold")
        assert found is not None and found["step"] == "002_AlphaFold_stitched_001"

    def test_counts_match_the_number_of_markers(self, tmp_path):
        assert job_status.status(self.suffixed(tmp_path))["counts"] == {
            "completed": 1, "failed": 1}


def test_a_folder_with_no_marker_is_still_a_pending_step(tmp_path):
    """The suffix fix must not swallow a step that genuinely has not finished."""
    job = tmp_path / "run_001"
    (job / "Logs").mkdir(parents=True)
    (job / "001_PDB").mkdir()
    io.open(job / "001_PDB_COMPLETED", "w", encoding="utf-8").write("")
    (job / "002_Boltz2").mkdir()
    assert {s["step"]: s["status"] for s in job_status.steps(job)} == {
        "001_PDB": "completed", "002_Boltz2": "pending"}


class TestAgainstAFrameworkBuiltRun:
    """One case built by the framework rather than by this file's memory of it.

    The fixtures here happen to be right — markers at the job root, logs in `Logs/` — but they
    are right by having been checked against a cluster once, not by construction. Its sibling
    reader `lineage.tables` was wrong the same way and nothing noticed. This case notices.
    """

    def test_markers_and_logs_are_read_where_the_framework_puts_them(self, produced_run):
        job, steps = produced_run([
            {"streams": {"structures": 2}},
            {"streams": {"sequences": 4}, "status": "FAILED"},
        ])
        found = job_status.steps(job, fs=LocalFS())
        assert [s["step"] for s in found] == sorted(steps)
        assert [s["status"] for s in found] == ["completed", "failed"]
        assert all(s["has_log"] for s in found), "every step wrote a log the reader must find"

    def test_the_failure_is_the_one_reported(self, produced_run):
        job, _steps = produced_run([{"streams": {"structures": 1}, "status": "FAILED"}])
        state = job_status.status(job, fs=LocalFS())
        assert state["failed"] == ["001_Mock"] and state["first_failure"] == "001_Mock"
