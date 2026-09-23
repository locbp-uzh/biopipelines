"""The operational log replaces `llm/log.sh`, so it has to hold up where that wrapper did not:
it must exist on every backend, live with the run rather than in a dated file, and never take
a run down with it when the filesystem says no.
"""

import io
import json

import pytest

from biopipelines import run_log


def test_records_land_next_to_the_run(tmp_path):
    job = tmp_path / "Examples" / "Binders_001"
    run_log.record(job, "submitted", job_id="4711")
    assert (job / "_operations.jsonl").exists(), "the log must be inside the job folder"
    assert run_log.read(job) == [
        {"ts": run_log.read(job)[0]["ts"], "action": "submitted", "job_id": "4711"}]


def test_appends_rather_than_overwrites(tmp_path):
    job = tmp_path / "job"
    for state in ("submitted", "running", "completed"):
        run_log.record(job, state)
    assert [e["action"] for e in run_log.read(job)] == ["submitted", "running", "completed"]


def test_timestamps_are_utc_and_unambiguous(tmp_path):
    job = tmp_path / "job"
    entry = run_log.record(job, "submitted")
    assert entry["ts"].endswith("+00:00"), (
        "a laptop and a cluster in different zones must produce comparable timestamps")


def test_command_records_argv_and_exit_code(tmp_path):
    job = tmp_path / "job"
    run_log.command(job, ["ssh", "cluster", "squeue -u me"], exit_code=0, host="cluster")
    entry = run_log.read(job)[0]
    assert entry["action"] == "command"
    assert entry["argv"] == ["ssh", "cluster", "squeue -u me"]
    assert entry["exit_code"] == 0 and entry["host"] == "cluster"


def test_command_truncates_output(tmp_path):
    job = tmp_path / "job"
    run_log.command(job, ["ssh", "cluster", "cat huge.log"], output="x" * 5000)
    recorded = run_log.read(job)[0]["output"]
    assert len(recorded) < 5000 and "more chars" in recorded, (
        "the log is a breadcrumb trail; full output belongs in <Job>/Logs/")


def test_none_fields_are_dropped(tmp_path):
    job = tmp_path / "job"
    run_log.command(job, ["ls"], exit_code=None, output=None)
    assert set(run_log.read(job)[0]) == {"ts", "action", "argv"}


def test_a_malformed_line_does_not_break_reading(tmp_path):
    """A log truncated by a killed job must still be readable."""
    job = tmp_path / "job"
    run_log.record(job, "submitted")
    with io.open(run_log.path_for(job), "a", encoding="utf-8") as handle:
        handle.write('{"action": "truncated mid-writ\n')
    run_log.record(job, "completed")
    assert [e["action"] for e in run_log.read(job)] == ["submitted", "completed"]


def test_reading_a_run_with_no_log_is_empty_not_an_error(tmp_path):
    assert run_log.read(tmp_path / "never-ran") == []


def test_a_write_failure_never_takes_the_run_down(tmp_path, monkeypatch):
    """Logging is not worth failing a six-hour campaign over."""
    monkeypatch.setattr(run_log.pathlib.Path, "mkdir",
                        lambda *a, **k: (_ for _ in ()).throw(OSError("read-only")))
    assert run_log.record(tmp_path / "job", "submitted") is None


def test_find_runs_orders_newest_first(tmp_path):
    import os
    import time
    for name, age in (("old", 100), ("new", 0)):
        job = tmp_path / name
        run_log.record(job, "submitted")
        stamp = time.time() - age
        os.utime(run_log.path_for(job), (stamp, stamp))
    assert [p.name for p in run_log.find_runs(tmp_path)] == ["new", "old"]


def test_records_are_valid_json_lines(tmp_path):
    """The file has to be machine-readable by something that is not us."""
    job = tmp_path / "job"
    run_log.record(job, "submitted", parameters={"gpu": "A100", "designs": 100})
    run_log.command(job, ["scp", "a", "b"], exit_code=1)
    for line in io.open(run_log.path_for(job), encoding="utf-8"):
        assert json.loads(line)


class TestBufferThenFlush:
    """`bp_submit` knows the submission parameters before the job folder exists."""

    def test_records_made_before_bind_are_held_then_written_in_order(self, tmp_path):
        log = run_log.RunLog()
        log.record("submitting", designs=500)
        log.command(["ssh", "cluster", "sbatch run.sh"], exit_code=0)
        assert run_log.read(tmp_path / "job") == [], "nothing should be on disk yet"
        assert len(log.pending) == 2

        log.bind(tmp_path / "job")
        assert [e["action"] for e in run_log.read(tmp_path / "job")] == ["submitting", "command"]
        assert log.pending == []

    def test_a_buffered_record_keeps_the_time_it_happened(self, tmp_path):
        log = run_log.RunLog()
        log.record("submitting")
        stamped = log.pending[0][1]["ts"]
        log.bind(tmp_path / "job")
        assert run_log.read(tmp_path / "job")[0]["ts"] == stamped, (
            "the useful timestamp is when it happened, not when the folder appeared")

    def test_records_after_bind_go_straight_to_disk(self, tmp_path):
        log = run_log.RunLog(tmp_path / "job")
        log.record("running")
        assert [e["action"] for e in run_log.read(tmp_path / "job")] == ["running"]
        assert log.pending == []

    def test_constructing_with_a_directory_binds_immediately(self, tmp_path):
        # Kept as a string, not a Path: a run on a cluster has a POSIX path that WindowsPath
        # would mangle into backslashes.
        assert run_log.RunLog(tmp_path / "job").job_dir == str(tmp_path / "job")


class TestRunHeader:
    """The first question when a rerun behaves differently is what changed around it."""

    def test_header_records_the_environment_that_produced_the_run(self, tmp_path):
        job = tmp_path / "job"
        run_log.header(job, pipeline="snap33_solubilize.py")
        entry = run_log.read(job)[0]
        assert entry["action"] == "run"
        assert entry["pipeline"] == "snap33_solubilize.py"
        for field in ("biopipelines", "host", "python"):
            assert entry.get(field), f"{field} missing from the run header"

    def test_header_reports_the_config_variant(self, tmp_path, monkeypatch):
        monkeypatch.setenv("BIOPIPELINES_CONFIG_VARIANT", "daint")
        run_log.header(tmp_path / "job")
        assert run_log.read(tmp_path / "job")[0]["config_variant"] == "daint"

    def test_a_buffered_header_captures_state_at_submission(self, tmp_path, monkeypatch):
        monkeypatch.setenv("BIOPIPELINES_CONFIG_VARIANT", "cluster")
        log = run_log.RunLog()
        log.header(pipeline="x.py")
        monkeypatch.setenv("BIOPIPELINES_CONFIG_VARIANT", "colab")
        log.bind(tmp_path / "job")
        assert run_log.read(tmp_path / "job")[0]["config_variant"] == "cluster", (
            "the header describes the checkout as it was at submission, not at flush")

    def test_header_survives_not_being_in_a_git_repo(self, tmp_path, monkeypatch):
        monkeypatch.setattr(run_log, "_git_commit", lambda: None)
        run_log.header(tmp_path / "job")
        assert "commit" not in run_log.read(tmp_path / "job")[0], "None fields are dropped"


class TestRemoteRuns:
    """The record must live with the run, and the run is usually on a cluster.

    The first real submission wrote nothing: `record` used local pathlib against a cluster path,
    and its never-raise contract turned that into silence. Both halves are pinned here.
    """

    class FakeRemoteFS:
        def __init__(self):
            self.files = {}

        def join(self, *parts):
            return "/".join(str(p).rstrip("/") for p in parts if str(p) != "")

        def is_file(self, path):
            return str(path) in self.files

        def append_text(self, path, text):
            self.files[str(path)] = self.files.get(str(path), "") + text
            return True

        def read_text(self, path):
            return self.files[str(path)]

    def test_a_record_is_written_through_the_given_filesystem(self):
        fs = self.FakeRemoteFS()
        run_log.record("/shares/x/Bio/Proj/run_001", "submitted", fs=fs, job_id="4711")
        assert "/shares/x/Bio/Proj/run_001/_operations.jsonl" in fs.files
        assert run_log.read("/shares/x/Bio/Proj/run_001", fs=fs)[0]["job_id"] == "4711"

    def test_a_posix_path_is_not_mangled_into_backslashes(self):
        fs = self.FakeRemoteFS()
        run_log.record("/shares/x/run_001", "submitted", fs=fs)
        assert list(fs.files) == ["/shares/x/run_001/_operations.jsonl"]

    def test_the_buffer_flushes_through_the_remote_filesystem(self):
        fs = self.FakeRemoteFS()
        log = run_log.RunLog(fs=fs)
        log.header(pipeline="foo.py")
        log.record("submitted", job_id="4711")
        assert fs.files == {}
        log.bind("/shares/x/run_001")
        assert [e["action"] for e in run_log.read("/shares/x/run_001", fs=fs)] == [
            "run", "submitted"]

    def test_a_filesystem_that_cannot_write_returns_none_rather_than_raising(self):
        class Refuses(TestRemoteRuns.FakeRemoteFS):
            def append_text(self, path, text):
                raise OSError("read-only")

        assert run_log.record("/x/run_001", "submitted", fs=Refuses()) is None
