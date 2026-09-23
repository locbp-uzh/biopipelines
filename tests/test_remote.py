"""The ssh layer: reading a cluster's runs from the machine the agent runs on.

The cluster is the main platform, so these paths matter more than the local ones. No cluster is
available in CI, so `subprocess.run` is replaced with a recorder: what matters here is the
command that would be sent, how its output is parsed, and that a failure becomes a message a
tool can return rather than a traceback.
"""

import subprocess

import pytest

from biopipelines import job_status, remote


class FakeSsh:
    """Answers canned output per command substring, and records everything asked."""

    def __init__(self, answers, default=(0, "", "")):
        self.answers = answers
        self.default = default
        self.calls = []

    def __call__(self, argv, **kwargs):
        self.calls.append(argv)
        command = argv[-1]
        for needle, result in self.answers.items():
            if needle in command:
                code, out, err = result
                return subprocess.CompletedProcess(argv, code, out, err)
        code, out, err = self.default
        return subprocess.CompletedProcess(argv, code, out, err)


@pytest.fixture
def fake(monkeypatch):
    def install(answers, default=(0, "", "")):
        recorder = FakeSsh(answers, default)
        monkeypatch.setattr(remote.subprocess, "run", recorder)
        return recorder
    return install


def test_filesystem_picks_local_without_a_host_and_ssh_with_one():
    assert isinstance(remote.filesystem(None), remote.LocalFS)
    assert isinstance(remote.filesystem("cluster"), remote.Ssh)


def test_the_remote_resolves_its_own_output_root(fake):
    """The local config expands <username> for the wrong user, so the cluster must answer."""
    recorder = fake({"ConfigManager": (0, "BP_OUTPUT_ROOT=/shares/group/user/Bio\n", "")})
    assert remote.Ssh("cluster").output_root() == "/shares/group/user/Bio"
    sent = recorder.calls[0][-1]
    assert sent.startswith("cd ~/biopipelines &&"), (
        "a quoted tilde is not expanded by the remote shell")


def test_output_root_failure_is_a_remote_error_not_a_traceback(fake):
    fake({"ConfigManager": (1, "", "ModuleNotFoundError: No module named 'biopipelines'")})
    with pytest.raises(remote.RemoteError, match="ModuleNotFoundError"):
        remote.Ssh("cluster").output_root()


def test_listdir_reads_trailing_slashes_as_directories(fake):
    fake({"ls -Ap": (0, "001_Mock/\n001_Mock_COMPLETED\nLogs/\n_operations.jsonl\n", "")})
    assert remote.Ssh("cluster").listdir("/runs/job_001") == [
        ("001_Mock", True), ("001_Mock_COMPLETED", False),
        ("Logs", True), ("_operations.jsonl", False)]


def test_listdir_of_a_missing_directory_is_empty(fake):
    fake({}, default=(2, "", "No such file or directory"))
    assert remote.Ssh("cluster").listdir("/nope") == []


def test_read_text_failure_names_the_host_and_path(fake):
    fake({"cat": (1, "", "cat: /x: No such file")})
    with pytest.raises(remote.RemoteError, match="cannot read /x on cluster"):
        remote.Ssh("cluster").read_text("/x")


def test_check_reports_an_unreachable_host_rather_than_raising(fake):
    fake({"echo ok": (255, "", "ssh: Could not resolve hostname cluster")})
    ok, detail = remote.Ssh("cluster").check()
    assert ok is False and "Could not resolve hostname" in detail


def test_a_missing_ssh_binary_is_a_remote_error(monkeypatch):
    monkeypatch.setattr(remote.subprocess, "run",
                        lambda *a, **k: (_ for _ in ()).throw(FileNotFoundError()))
    with pytest.raises(remote.RemoteError, match="ssh is not on PATH"):
        remote.Ssh("cluster").run("echo ok")


def test_a_timeout_names_the_host_and_the_limit(monkeypatch):
    def boom(*a, **k):
        raise subprocess.TimeoutExpired(cmd="ssh", timeout=30)
    monkeypatch.setattr(remote.subprocess, "run", boom)
    with pytest.raises(remote.RemoteError, match="ssh cluster timed out after 30s"):
        remote.Ssh("cluster").run("echo ok")


def test_remote_paths_are_posix_regardless_of_the_local_platform():
    """Building a cluster path with os.path.join on Windows would produce backslashes."""
    assert remote.Ssh("cluster").join("/shares/user/Bio", "Proj", "run_001") == \
        "/shares/user/Bio/Proj/run_001"


class TestResolveJobDir:
    def test_an_absolute_path_is_used_as_given(self):
        fs = remote.LocalFS()
        assert remote.resolve_job_dir(fs, "/shares/x/Bio/Proj/run_001") == "/shares/x/Bio/Proj/run_001"

    def test_a_project_slash_job_resolves_under_the_root(self, fake):
        """`bp_runs` prints "2 runs in Template", so `Template/example_002` is what comes next.

        Returned unchanged it was a relative path resolving to nothing, and bp_status reported
        "No run directory" for a run bp_runs had just listed.
        """
        fake({"ConfigManager": (0, "BP_OUTPUT_ROOT=/shares/x/Bio\n", "")})
        fs = remote.Ssh("cluster")
        assert remote.resolve_job_dir(fs, "Proj/run_001") == "/shares/x/Bio/Proj/run_001"

    def test_a_named_project_is_joined_under_the_root(self, fake):
        fake({"ConfigManager": (0, "BP_OUTPUT_ROOT=/shares/x/Bio\n", "")})
        fs = remote.Ssh("cluster")
        assert remote.resolve_job_dir(fs, "run_001", project="Proj") == "/shares/x/Bio/Proj/run_001"

    def test_without_a_project_the_projects_are_searched(self, fake):
        fake({"ConfigManager": (0, "BP_OUTPUT_ROOT=/shares/x/Bio\n", ""),
              "for d in */": (0, "Beta/run_001\n", "")},
             default=(1, "", ""))
        fs = remote.Ssh("cluster")
        assert remote.resolve_job_dir(fs, "run_001") == "/shares/x/Bio/Beta/run_001"

    def test_a_job_name_in_two_projects_is_refused_naming_both(self, fake):
        """Picking the first silently is how a cancel lands on another project's run."""
        fake({"ConfigManager": (0, "BP_OUTPUT_ROOT=/shares/x/Bio\n", ""),
              "for d in */": (0, "Alpha/run_001\nBeta/run_001\n", "")},
             default=(1, "", ""))
        with pytest.raises(remote.RemoteError, match="Alpha/run_001.*Beta/run_001"):
            remote.resolve_job_dir(remote.Ssh("cluster"), "run_001")

    def test_an_unfound_job_still_returns_a_path_to_report(self, fake):
        fake({"ConfigManager": (0, "BP_OUTPUT_ROOT=/shares/x/Bio\n", ""), "ls -Ap": (0, "", "")},
             default=(1, "", ""))
        assert remote.resolve_job_dir(remote.Ssh("cluster"), "run_001") == "/shares/x/Bio/run_001"


def test_job_status_reads_a_run_over_the_remote_filesystem(fake):
    """The whole point: the same status logic answering for a run on a login node."""
    fake({"ls -Ap /runs/job_001/Logs": (0, "001_A.log\n002_B.log\n", ""),
          "ls -Ap /runs/job_001 ": (0, "001_A/\n001_A_COMPLETED\n002_B/\n002_B_FAILED\nLogs/\n", ""),
          "test -d /runs/job_001": (0, "", ""),
          "test -f": (1, "", ""),
          "cat /runs/job_001/Logs/002_B.log": (0, "boom\ntraceback\n", "")},
         default=(1, "", ""))
    fs = remote.Ssh("cluster")

    state = job_status.status("/runs/job_001", fs=fs)
    assert state["counts"] == {"completed": 1, "failed": 1}
    assert state["first_failure"] == "002_B"
    assert "First failure: 002_B" in job_status.summarize(state)

    assert job_status.log("/runs/job_001", "002_B", fs=fs)["text"] == "boom\ntraceback"


def test_ssh_output_is_decoded_as_utf8_not_the_local_codepage(monkeypatch):
    """Found live: `text=True` alone uses cp1252 on Windows and dies on the first non-ASCII
    byte in a cluster log, taking down subprocess's reader thread rather than returning."""
    monkeypatch.setattr(remote.subprocess, "run",
                        lambda argv, **k: subprocess.CompletedProcess(argv, 0, "Å ok \x81".encode("latin-1"), b""))
    code, out, _ = remote.Ssh("cluster").run("echo ok")
    assert code == 0 and "ok" in out


def test_stdin_is_sent_as_bytes_so_windows_cannot_add_crlf(monkeypatch):
    seen = {}

    def record(argv, **kwargs):
        seen.update(kwargs)
        return subprocess.CompletedProcess(argv, 0, b"", b"")

    monkeypatch.setattr(remote.subprocess, "run", record)
    remote.Ssh("cluster").run("cat >> x", stdin="a\nb\n")
    assert seen["input"] == b"a\nb\n" and "text" not in seen


def test_the_host_cannot_become_an_ssh_option():
    with pytest.raises(remote.RemoteError, match="not a valid ssh host"):
        remote.Ssh("-oProxyCommand=touch /tmp/x")


def test_ssh_passes_batchmode_and_ends_options_before_the_host(monkeypatch):
    seen = {}
    monkeypatch.setattr(remote.subprocess, "run",
                        lambda argv, **k: seen.setdefault("argv", argv) and subprocess.CompletedProcess(argv, 0, b"", b""))
    remote.Ssh("cluster").run("true")
    argv = seen["argv"]
    assert "BatchMode=yes" in argv and argv[argv.index("--") + 1] == "cluster"


def test_an_ssh_failure_is_not_reported_as_an_absent_path(fake):
    fake({}, default=(255, "", "Permission denied (publickey)"))
    with pytest.raises(remote.RemoteError, match="Permission denied"):
        remote.Ssh("cluster").is_dir("/runs/job_001")


@pytest.mark.parametrize("path, quoted", [
    ("$SCRATCH/biopipelines", '"${SCRATCH:?}"/biopipelines'),
    ("${SCRATCH}", '"${SCRATCH:?}"'),
    ("~/x y", "~/'x y'"),
    ("/a;b", "'/a;b'"),
])
def test_a_leading_variable_stays_expandable_and_the_rest_is_quoted(path, quoted):
    assert remote.quote(path) == quoted


@pytest.mark.parametrize("path", ["x;touch pwned", "../.bashrc", "a/../../b", "$(id)", "a b"])
def test_an_unsafe_scp_path_is_refused(path):
    with pytest.raises(remote.RemoteError, match="refusing remote path"):
        remote.validate_scp_path(path)


def test_none_streams_become_empty_strings(monkeypatch):
    monkeypatch.setattr(remote.subprocess, "run",
                        lambda argv, **k: subprocess.CompletedProcess(argv, 0, None, None))
    assert remote.Ssh("cluster").run("echo ok") == (0, "", "")


def test_the_root_snippet_does_not_import_output_root(monkeypatch):
    """It must run against an older cluster checkout: the first live S3IT attempt failed with
    `ImportError: cannot import name 'output_root'` because the remote predated that helper."""
    assert "output_root" not in remote._ROOT_SNIPPET
    assert "ConfigManager" in remote._ROOT_SNIPPET


def test_the_root_answer_is_taken_from_the_sentinel_line(fake):
    fake({"ConfigManager": (0, "[biopipelines] config variant 'cluster': ...\n"
                               "BP_OUTPUT_ROOT=/shares/x/Bio\n", "")})
    assert remote.Ssh("cluster").output_root() == "/shares/x/Bio"


class TestSavedConnection:
    """A first-time user with only MCP tools cannot discover their cluster's alias or repo path,
    and repeating both on every call is noise. `bp_setup` verifies a pair once and saves it."""

    @pytest.fixture(autouse=True)
    def isolated(self, tmp_path, monkeypatch):
        monkeypatch.setattr(remote, "SETTINGS_FILE", tmp_path / ".bp-mcp.json")

    def test_nothing_saved_yields_no_host_and_the_default_repo(self):
        assert remote.defaults() == (None, remote.DEFAULT_REPO)

    def test_a_saved_pair_fills_in_omitted_arguments(self):
        remote.save_settings(host="cluster", repo="~/biopipelines-locbp")
        assert remote.defaults() == ("cluster", "~/biopipelines-locbp")

    def test_an_explicit_argument_wins_over_the_saved_one(self):
        remote.save_settings(host="cluster", repo="~/a")
        assert remote.defaults(host="daint", repo="~/b") == ("daint", "~/b")

    def test_saving_merges_rather_than_replacing(self):
        remote.save_settings(host="cluster", repo="~/a")
        remote.save_settings(host="cluster", python="python3")
        entry = remote.settings_for("cluster")
        assert entry == {"repo": "~/a", "python": "python3"}

    def test_a_second_cluster_does_not_replace_the_first(self):
        """A lab with two clusters is the ordinary case, not an edge case."""
        remote.save_settings(host="cluster", repo="~/biopipelines-locbp")
        remote.save_settings(host="daint", repo="/scratch/me/biopipelines", variant="daint")
        assert remote.settings_for("cluster")["repo"] == "~/biopipelines-locbp"
        assert remote.settings_for("daint")["variant"] == "daint"

    def test_no_repo_is_ever_inherited_across_hosts(self):
        """The defect this shape exists to remove.

        The flat file carried one repo for whichever host was named, so a call against Daint
        got S3IT's `~/biopipelines-locbp` — a path that resolves to nothing there, or worse to
        something else's tree. A wrong path that happens to exist is the failure that produces
        a confident answer about the wrong machine.
        """
        remote.save_settings(host="cluster", repo="~/biopipelines-locbp")
        assert remote.defaults("daint") == ("daint", remote.DEFAULT_REPO)

    def test_each_host_keeps_its_own_interpreter_and_variant(self):
        """Daint's login nodes have no `python` at all, and it needs its own config variant."""
        remote.save_settings(host="daint", repo="/scratch/me/bp",
                             python="/scratch/me/venv/bin/python", variant="daint")
        remote.save_settings(host="cluster", repo="~/bp")
        assert remote.connection("daint")[2:4] == ("/scratch/me/venv/bin/python", "daint")
        assert remote.connection("cluster")[2:4] == ("python", None)

    def test_the_default_host_is_what_an_unnamed_call_uses(self):
        remote.save_settings(host="cluster", repo="~/a")
        remote.save_settings(host="daint", repo="~/b")       # newest becomes default
        assert remote.defaults()[0] == "daint"
        remote.save_settings(host="cluster", repo="~/a")     # switching back is a save
        assert remote.defaults()[0] == "cluster"

    def test_saving_without_making_it_default_leaves_the_default_alone(self):
        remote.save_settings(host="cluster", repo="~/a")
        remote.save_settings(host="daint", repo="~/b", make_default=False)
        assert remote.defaults()[0] == "cluster"
        assert remote.settings_for("daint")["repo"] == "~/b"

    def test_a_flat_settings_file_is_read_as_one_host(self):
        """Nobody should have to run bp_setup again because the file grew a level."""
        remote.SETTINGS_FILE.write_text(
            '{"host": "cluster", "repo": "~/biopipelines-locbp", "python": "python"}',
            encoding="utf-8")
        assert remote.defaults() == ("cluster", "~/biopipelines-locbp")
        assert remote.settings_for("cluster")["python"] == "python"

    def test_a_corrupt_settings_file_is_ignored_not_fatal(self):
        remote.SETTINGS_FILE.write_text("{not json", encoding="utf-8")
        assert remote.load_settings() == {"default": None, "hosts": {}}
        assert remote.defaults() == (None, remote.DEFAULT_REPO)


class TestDiagnose:
    """The stages are ordered so the first failure names the thing to fix, rather than surfacing
    as a config error three steps later."""

    def test_an_unreachable_alias_stops_at_the_first_stage(self, fake):
        fake({"echo ok": (255, "", "ssh: Could not resolve hostname cluster")})
        ok, stages = remote.diagnose("cluster")
        assert ok is False and len(stages) == 1
        assert stages[0][0] == "ssh alias reachable"

    def test_a_missing_repo_is_named_as_a_missing_repo(self, fake):
        fake({"echo ok": (0, "ok", "")}, default=(1, "", ""))
        ok, stages = remote.diagnose("cluster", repo="~/wrong-name")
        assert ok is False
        assert stages[-1][0] == "repository at ~/wrong-name"
        assert "clone it there" in stages[-1][2]

    def test_an_unimportable_framework_stops_before_the_config(self, fake):
        fake({"echo ok": (0, "ok", ""),
              "test -d": (0, "present", ""),
              "__version__": (1, "", "ModuleNotFoundError: No module named 'biopipelines'")},
             default=(1, "", ""))
        ok, stages = remote.diagnose("cluster")
        assert ok is False
        assert stages[-1][0] == "biopipelines importable"
        assert "ModuleNotFoundError" in stages[-1][2]

    def test_a_healthy_cluster_reports_every_stage(self, fake):
        fake({"echo ok": (0, "ok", ""),
              "test -d": (0, "present", ""),
              "__version__": (0, "1.4.0", ""),
              "ConfigManager": (0, "BP_OUTPUT_ROOT=/shares/x/Bio\n", ""),
              "ls -Ap": (0, "Alpha/\nBeta/\n", "")})
        ok, stages = remote.diagnose("cluster")
        assert ok is True
        assert [name for name, _, _ in stages] == [
            "ssh alias reachable", "repository at ~/biopipelines", "biopipelines importable",
            "console scripts on the interpreter's path",
            "config variant resolves", "output root readable"]
        assert stages[-1][2] == "2 projects"


class TestDownload:
    """Pulling results back: a run's structures often exist only on the cluster."""

    def test_a_directory_is_copied_recursively_without_being_told(self, monkeypatch):
        argv_seen = []

        def record(argv, **kwargs):
            argv_seen.append(argv)
            # `is_dir` runs first as an ssh call, then scp.
            return subprocess.CompletedProcess(argv, 0, "", "")

        monkeypatch.setattr(remote.subprocess, "run", record)
        remote.Ssh("cluster").download("/shares/x/run_001/005_ProteinMPNN", "./out")
        assert "-r" in argv_seen[-1], "a step folder must not be copied as a single file"

    def test_a_file_is_copied_without_r(self, monkeypatch):
        calls = []

        def record(argv, **kwargs):
            calls.append(argv)
            # test -d fails -> not a directory
            return subprocess.CompletedProcess(argv, 1 if "test -d" in argv[-1] else 0, "", "")

        monkeypatch.setattr(remote.subprocess, "run", record)
        remote.Ssh("cluster").download("/shares/x/run_001/a.pdb", "a.pdb")
        assert "-r" not in calls[-1]

    def test_an_explicit_recursive_flag_skips_the_probe(self, monkeypatch):
        calls = []
        monkeypatch.setattr(remote.subprocess, "run",
                            lambda argv, **k: (calls.append(argv),
                                               subprocess.CompletedProcess(argv, 0, "", ""))[1])
        remote.Ssh("cluster").download("/x", "y", recursive=False)
        assert len(calls) == 1 and calls[0][0] == "scp"

    def test_a_missing_scp_is_a_remote_error(self, monkeypatch):
        monkeypatch.setattr(remote.subprocess, "run",
                            lambda *a, **k: (_ for _ in ()).throw(FileNotFoundError()))
        with pytest.raises(remote.RemoteError, match="scp is not on PATH"):
            remote.Ssh("cluster").download("/x", "y", recursive=False)


def test_walk_files_is_one_round_trip(fake):
    """views._walk listed every folder separately: one ssh connection per directory."""
    recorder = fake({"find ": (0, "/r/step/top.png\n/r/step/_extras/i.png\n", "")})
    found = remote.Ssh("cluster").walk_files("/r/step", depth=2)
    assert found == [("/r/step/top.png", "top.png"), ("/r/step/_extras/i.png", "i.png")]
    assert len(recorder.calls) == 1
    assert "-maxdepth 3" in recorder.calls[0][-1] and "! -name _extras -prune" in recorder.calls[0][-1]
