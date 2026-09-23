"""Submitting a pipeline, and recording what was submitted against the run it creates.

`bp_submit` is the first tool that changes something, so it is the first that writes to a run's
operations log — and it is where `RunLog`'s buffer earns its place: the submission parameters
are known before `./submit` runs, the job directory it creates is not known until the output is
parsed.

The output shapes asserted here are the ones `submit` actually prints (`Runtime directory:`,
`Batch N submitted: Job ID`, and its explicit failure lines). If that script's wording changes,
these tests are what notices.
"""

import io

import pytest

from biopipelines import remote, run_log, submission

SUCCESS = """Generating pipeline: my_pipelines/foo.py
Found 1 pipeline(s) to submit
==========================================
Processing pipeline: foo
Runtime directory: /shares/x/Bio/Proj/foo_003/RunTime
Detected multi-batch pipeline
Submitting batch 001...
  Job name: foo_b1
✓ Batch 001 submitted: Job ID 4711
Submitting batch 002...
✓ Batch 002 submitted: Job ID 4712
Pipeline foo: 2 batch(es) submitted
"""

GEN_FAILED = """Generating pipeline: my_pipelines/foo.py
Pipeline generation failed with exit code 1
Re-run with -v to see full output
"""

SUB_FAILED = """Runtime directory: /shares/x/Bio/Proj/foo_003/RunTime
✗ Submission failed (slurm): sbatch: error: Invalid partition name specified
"""


class FakeSsh:
    """Stands in for a cluster: runs commands, and is also the filesystem the log writes to.

    That double duty is the point — the record has to land on the machine the run is on, which
    is what the first real submission got wrong.
    """

    def __init__(self, output="", code=0, upload_code=0, listing=None, dirs=()):
        self.host, self.repo, self.timeout = "cluster", "~/biopipelines", 30
        self.output, self.code, self.upload_code = output, code, upload_code
        self.commands, self.uploads, self.files = [], [], {}
        # `listing` is {path: [(name, is_dir)]}, for the tools that look before they act.
        self.listing, self.dirs = listing or {}, set(dirs)

    def run(self, command, timeout=None, stdin=None):
        self.commands.append(command)
        return self.code, self.output, ""

    def upload(self, local, remote_path, timeout=None):
        self.uploads.append((str(local), remote_path))
        return self.upload_code, "", "" if self.upload_code == 0 else "Permission denied"

    # --- filesystem half -------------------------------------------------
    def join(self, *parts):
        return "/".join(str(p).rstrip("/") for p in parts if str(p) != "")

    def is_file(self, path):
        return str(path) in self.files

    def is_dir(self, path):
        return str(path) in self.dirs or str(path) in self.listing

    def listdir(self, path):
        return list(self.listing.get(str(path), []))

    def append_text(self, path, text):
        self.files[str(path)] = self.files.get(str(path), "") + text
        return True

    def read_text(self, path):
        return self.files[str(path)]


class TestParse:
    def test_the_run_directory_is_the_parent_of_runtime(self):
        assert submission.parse(SUCCESS)["jobs"] == ["/shares/x/Bio/Proj/foo_003"]

    def test_every_batch_job_id_is_captured(self):
        assert submission.parse(SUCCESS)["job_ids"] == ["4711", "4712"]

    def test_a_clean_submission_is_ok(self):
        assert submission.parse(SUCCESS)["ok"] is True

    def test_generation_failure_is_reported_not_silently_empty(self):
        result = submission.parse(GEN_FAILED)
        assert result["ok"] is False and result["jobs"] == []
        assert "pipeline generation failed (exit 1)" in result["failures"]

    def test_a_scheduler_rejection_is_not_success_despite_a_run_directory(self):
        """The directory exists but nothing was queued — the dangerous shape to get wrong."""
        result = submission.parse(SUB_FAILED)
        assert result["jobs"] == ["/shares/x/Bio/Proj/foo_003"]
        assert result["ok"] is False
        assert "Invalid partition name" in result["failures"][0]

    def test_no_runtime_directories_is_a_named_failure(self):
        result = submission.parse("No runtime directories found in pipeline output\n")
        assert result["ok"] is False and "nothing was generated" in result["failures"][0]


class TestSubmit:
    def test_the_header_and_submission_land_in_the_new_run_directory(self, tmp_path):
        job = tmp_path / "foo_003"
        ssh = FakeSsh(SUCCESS.replace("/shares/x/Bio/Proj/foo_003", str(job)))
        result = submission.submit(ssh, "my_pipelines/foo.py")

        assert result["ok"] is True
        entries = run_log.read(job, fs=ssh)
        assert [e["action"] for e in entries] == ["run", "submitted"]
        assert entries[0]["pipeline"] == "my_pipelines/foo.py"
        assert entries[1]["job_ids"] == ["4711", "4712"]

    def test_records_carry_the_time_they_happened_not_the_flush(self, tmp_path):
        job = tmp_path / "foo_003"
        ssh = FakeSsh(SUCCESS.replace("/shares/x/Bio/Proj/foo_003", str(job)))
        submission.submit(ssh, "my_pipelines/foo.py")
        first, second = run_log.read(job, fs=ssh)
        assert first["ts"] <= second["ts"]

    def test_an_upload_is_copied_before_submitting_and_recorded(self, tmp_path):
        script = tmp_path / "foo.py"
        script.write_text("print('hi')", encoding="utf-8")
        job = tmp_path / "foo_003"
        ssh = FakeSsh(SUCCESS.replace("/shares/x/Bio/Proj/foo_003", str(job)))

        submission.submit(ssh, "my_pipelines/foo.py", upload=script)
        assert ssh.uploads == [(str(script), "~/biopipelines/my_pipelines/foo.py")]
        assert [e["action"] for e in run_log.read(job, fs=ssh)] == ["run", "command", "submitted"]

    def test_a_failed_upload_does_not_submit(self, tmp_path):
        ssh = FakeSsh(SUCCESS, upload_code=1)
        result = submission.submit(ssh, "my_pipelines/foo.py", upload=tmp_path / "foo.py")
        assert result["ok"] is False
        assert ssh.commands == [], "nothing should run on the cluster after a failed copy"
        assert "could not copy" in result["failures"][0]

    def test_a_nonzero_exit_with_no_parsed_failure_is_still_a_failure(self, tmp_path):
        ssh = FakeSsh("", code=127)
        result = submission.submit(ssh, "my_pipelines/foo.py")
        assert result["ok"] is False and "./submit exited 127" in result["failures"]

    def test_the_command_runs_in_the_repo(self, tmp_path):
        job = tmp_path / "foo_003"
        ssh = FakeSsh(SUCCESS.replace("/shares/x/Bio/Proj/foo_003", str(job)))
        submission.submit(ssh, "my_pipelines/foo.py")
        assert ssh.commands[0].startswith("cd ~/biopipelines && ./submit"), (
            "a quoted tilde would not expand on the remote shell")


class TestSummarize:
    def test_success_says_nothing_has_run_yet(self):
        text = submission.summarize({**submission.parse(SUCCESS), "output": SUCCESS},
                                    host="cluster")
        assert "Submitted on cluster" in text
        assert "Nothing has run yet" in text, (
            "the scheduler accepting a job is not the job succeeding")

    def test_failure_shows_the_reason_and_the_tail(self):
        text = submission.summarize({**submission.parse(SUB_FAILED), "output": SUB_FAILED},
                                    host="cluster")
        assert "Submission failed on cluster" in text
        assert "Invalid partition name" in text
        assert "Last lines of ./submit" in text


def test_read_only_tools_write_nothing_to_the_log(tmp_path):
    """Only actions that change something earn a record — that is why log.sh was retired."""
    from biopipelines import job_status
    job = tmp_path / "run_001"
    (job / "Logs").mkdir(parents=True)
    io.open(job / "001_A_COMPLETED", "w", encoding="utf-8").write("")
    (job / "001_A").mkdir()

    job_status.status(job)
    job_status.log(job, "001_A")
    assert run_log.read(job) == [], "a status poll must not add to the run's record"


def test_the_record_is_written_to_the_cluster_not_this_machine(tmp_path):
    """The first real submission produced no record: the log used local pathlib against a
    cluster path, and never-raise turned the failure into silence."""
    ssh = FakeSsh(SUCCESS)
    submission.submit(ssh, "my_pipelines/foo.py")
    assert ssh.files, "nothing was written to the cluster filesystem"
    assert list(ssh.files) == ["/shares/x/Bio/Proj/foo_003/_operations.jsonl"]
    assert run_log.read("/shares/x/Bio/Proj/foo_003") == [], (
        "and nothing should have been written locally")



RUNTIME = "/runs/job_001/RunTime"


def _ssh(scripts, output="", code=0, log=None):
    """A cluster whose run has these batch scripts, and optionally an operations log."""
    ssh = FakeSsh(output=output, code=code,
                  listing={RUNTIME: [(name, False) for name in scripts]})
    if log is not None:
        ssh.files["/runs/job_001/_operations.jsonl"] = log
    return ssh


class TestResubmit:
    """A campaign that fails at step 4 of 9 is the ordinary case.

    Until there was a tool for it, recovering meant handing the user a shell command — the
    moment the premise of driving BioPipelines through tools visibly breaks.
    """

    def test_a_single_batch_needs_no_script_named(self):
        ssh = _ssh(["slurm.sh", "config.sh"], output="Batch 1 submitted: Job ID 991\n")
        result = submission.resubmit(ssh, "/runs/job_001")
        assert result["ok"] and result["script"] == "slurm.sh" and result["job_ids"] == ["991"]

    def test_config_sh_is_not_mistaken_for_a_batch_script(self):
        assert submission.job_scripts(_ssh(["slurm.sh", "config.sh"]), "/runs/job_001") == \
            ["slurm.sh"]

    def test_several_batches_must_be_disambiguated(self):
        """Resubmitting the wrong batch spends compute on work that already succeeded."""
        result = submission.resubmit(_ssh(["slurm_batch1.sh", "slurm_batch2.sh"]), "/runs/job_001")
        assert not result["ok"]
        assert result["available"] == ["slurm_batch1.sh", "slurm_batch2.sh"]
        assert "name the one" in result["failures"][0]

    def test_a_named_script_that_does_not_exist_lists_the_real_ones(self):
        result = submission.resubmit(_ssh(["slurm_batch1.sh"]), "/runs/job_001",
                                     script="slurm_batch9.sh")
        assert not result["ok"] and result["available"] == ["slurm_batch1.sh"]

    def test_a_run_with_no_batch_script_says_so(self):
        result = submission.resubmit(_ssh([]), "/runs/job_001")
        assert not result["ok"] and "nothing to resubmit" in result["failures"][0]

    def test_dependencies_are_stripped_unless_asked(self):
        """The script on disk names the original run's ids, long since finished or aged out."""
        ssh = _ssh(["slurm.sh"], output="submitted: Job ID 7\n")
        submission.resubmit(ssh, "/runs/job_001")
        assert not any("--keep-dependencies" in c for c in ssh.commands)
        submission.resubmit(ssh, "/runs/job_001", keep_dependencies=True)
        assert any("--keep-dependencies" in c for c in ssh.commands)

    def test_it_is_recorded_against_the_run(self):
        ssh = _ssh(["slurm.sh"], output="submitted: Job ID 7\n")
        submission.resubmit(ssh, "/runs/job_001")
        written = ssh.files["/runs/job_001/_operations.jsonl"]
        assert '"action": "resubmitted"' in written and '"7"' in written

    def test_a_failing_resubmit_is_reported_not_recorded_as_success(self):
        ssh = _ssh(["slurm.sh"], output="sbatch: error: Invalid partition\n", code=1)
        result = submission.resubmit(ssh, "/runs/job_001")
        assert not result["ok"]
        assert "Invalid partition" in submission.summarize_action(result, host="cluster")


def _queued(ssh, *ids):
    """Make `squeue` report these ids as still held, the way a live cluster would."""
    plain = ssh.run

    def run(command, timeout=None, stdin=None):
        if command.startswith("squeue"):
            ssh.commands.append(command)
            return 0, "".join(f"{i}\n" for i in ids), ""
        return plain(command, timeout=timeout, stdin=stdin)

    ssh.run = run
    return ssh


class TestCancel:
    def test_ids_come_from_the_runs_own_log(self):
        """Making the caller find the ids first is what pushes an agent back to a shell."""
        ssh = _queued(_ssh(["slurm.sh"], log='{"action": "submitted", "job_ids": ["6085300", "6085301"]}\n'
                                             '{"action": "status"}\n'), "6085300", "6085301")
        result = submission.cancel(ssh, job_dir="/runs/job_001")
        assert result["ok"] and result["job_ids"] == ["6085300", "6085301"]
        assert any("scancel" in c and "6085300" in c for c in ssh.commands)

    def test_explicit_ids_win_over_the_log(self):
        ssh = _queued(_ssh(["slurm.sh"], log='{"action": "submitted", "job_ids": ["1"]}\n'), "1", "42")
        result = submission.cancel(ssh, job_dir="/runs/job_001", job_ids=["42"])
        assert result["job_ids"] == ["42"]

    def test_no_ids_anywhere_is_an_error_not_a_bare_scancel(self):
        """A `scancel` with no ids would be a very bad thing to send by accident."""
        ssh = _ssh(["slurm.sh"])
        result = submission.cancel(ssh, job_dir="/runs/job_001")
        assert not result["ok"] and result["job_ids"] == []
        assert not any(c.strip() == "scancel" for c in ssh.commands)

    def test_the_cancellation_is_recorded(self):
        ssh = _queued(_ssh(["slurm.sh"], log='{"action": "submitted", "job_ids": ["9"]}\n'), "9")
        submission.cancel(ssh, job_dir="/runs/job_001")
        assert '"action": "cancelled"' in ssh.files["/runs/job_001/_operations.jsonl"]

    def test_the_summary_warns_that_a_cancelled_step_reads_as_pending(self):
        ssh = _queued(_ssh(["slurm.sh"]), "42")
        result = submission.cancel(ssh, job_ids=["42"])
        text = submission.summarize_action(result, host="cluster")
        assert "'pending'" in text and "rather than becoming 'failed'" in text, (
            "a cancelled step writes no FAILED marker, and an agent that does not know that "
            "will report the cancel as a crash")

    def test_a_dry_run_names_the_live_ids_and_cancels_nothing(self):
        ssh = _queued(_ssh(["slurm.sh"], log='{"action": "submitted", "job_ids": ["7", "8"]}\n'), "8")
        result = submission.cancel(ssh, job_dir="/runs/job_001", dry_run=True)
        assert result["job_ids"] == ["8"] and result["finished"] == ["7"]
        assert not any(c.startswith("scancel") for c in ssh.commands)
        assert "confirm=True" in submission.summarize_action(result)

    def test_finished_ids_are_left_alone(self):
        """scancel on a dead id exits non-zero and made a successful cancel read as a failure."""
        ssh = _queued(_ssh(["slurm.sh"], log='{"action": "submitted", "job_ids": ["7", "8"]}\n'), "8")
        submission.cancel(ssh, job_dir="/runs/job_001")
        assert [c for c in ssh.commands if c.startswith("scancel")] == ["scancel 8"]


class TestSubmitRefusesRawShell:
    def test_extra_is_limited_to_the_flags_submit_takes(self):
        ssh = FakeSsh()
        result = submission.submit(ssh, "my_pipelines/x.py", extra="; rm -rf ~")
        assert not result["ok"] and not ssh.commands

    def test_verbose_becomes_dash_v(self):
        ssh = FakeSsh(output="Runtime directory: /runs/job_001/RunTime\nsubmitted: Job ID 5\n")
        submission.submit(ssh, "my_pipelines/x.py", verbose=True)
        assert any("./submit -v " in c for c in ssh.commands)

    def test_an_upload_cannot_escape_the_repo(self):
        ssh = FakeSsh()
        result = submission.submit(ssh, "../.bashrc", upload="x.py")
        assert not result["ok"] and not ssh.uploads

    def test_a_timeout_is_reported_as_an_unknown_outcome(self):
        ssh = FakeSsh()
        ssh.run = lambda *a, **k: (_ for _ in ()).throw(submission.RemoteError("ssh cluster timed out after 600s"))
        result = submission.submit(ssh, "my_pipelines/x.py")
        assert result.get("outcome_unknown") and "before submitting again" in result["failures"][0]


def test_the_submit_header_files_the_laptop_under_submitted_from():
    """bp_submit runs on the laptop; its version and commit are not the run's identity."""
    import json
    ssh = FakeSsh(output="Runtime directory: /runs/job_001/RunTime\nsubmitted: Job ID 5\n")
    submission.submit(ssh, "my_pipelines/x.py")
    header = json.loads(ssh.files["/runs/job_001/_operations.jsonl"].splitlines()[0])
    assert header["action"] == "run" and "submitted_from" in header
    assert "commit" not in header and header["host"] == "cluster"


def test_an_array_task_id_is_recognised_as_live():
    ssh = _queued(FakeSsh(), "123_4", "456+0")
    assert submission.live_job_ids(ssh, ["123_4", "456", "789"]) == ["123_4", "456"]


def test_an_unreachable_controller_is_not_read_as_nothing_to_cancel():
    ssh = FakeSsh()
    ssh.run = lambda *a, **k: (1, "", "slurm_load_jobs error: Unable to contact slurm controller")
    with pytest.raises(remote.RemoteError, match="squeue failed"):
        submission.live_job_ids(ssh, ["123"])
