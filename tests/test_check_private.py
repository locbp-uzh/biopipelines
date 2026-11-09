"""`versions/check_private.py` names every `#PRIVATE` file a publish would carry.

The checker runs against a real repository built in a tmp_path, because what it
asserts is a property of git's object graph: a merge publishes everything
reachable from HEAD and not from the target, deleted files included.
"""

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
CHECKER = ROOT / "versions" / "check_private.py"


def git(repo, *args):
    done = subprocess.run(["git", *args], cwd=repo, capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    return done.stdout


def commit(repo, message):
    git(repo, "add", "-A")
    git(repo, "-c", "user.email=t@t", "-c", "user.name=t", "commit", "-q",
        "--no-verify", "-m", message)


@pytest.fixture
def repo(tmp_path):
    """A repo whose `base` branch is the publish target and whose HEAD is ahead of it."""
    repo = tmp_path / "r"
    repo.mkdir()
    git(repo, "init", "-q", "-b", "base")
    (repo / "public.py").write_text("print('hello')\n")
    commit(repo, "base")
    git(repo, "checkout", "-q", "-b", "work")
    return repo


def check(repo, *args):
    return subprocess.run([sys.executable, str(CHECKER), *args],
                          cwd=repo, capture_output=True, text=True, timeout=120)


class TestWhatIsFlagged:

    def test_a_clean_range_passes(self, repo):
        (repo / "more.py").write_text("x = 1\n")
        commit(repo, "ordinary work")
        done = check(repo, "base")
        assert done.returncode == 0, done.stdout + done.stderr
        assert "no #PRIVATE files" in done.stdout

    def test_a_marked_file_at_head_is_named(self, repo):
        (repo / "secret.py").write_text("#PRIVATE\nx = 1\n")
        commit(repo, "add a private file")
        done = check(repo, "base")
        assert done.returncode == 1
        assert "secret.py" in done.stderr

    def test_a_marked_file_deleted_before_head_is_still_named(self, repo):
        """The reason the whole range is checked: a merge publishes the ancestry too,
        so a file deleted in a later commit stays readable in the published history."""
        (repo / "secret.py").write_text("#PRIVATE\nx = 1\n")
        commit(repo, "add a private file")
        (repo / "secret.py").unlink()
        commit(repo, "delete it again")
        assert not (repo / "secret.py").exists()
        done = check(repo, "base")
        assert done.returncode == 1, done.stdout
        assert "secret.py" in done.stderr

    def test_a_marked_version_of_a_surviving_file_is_still_named(self, repo):
        """The marker was removed in a later commit; the marked blob is still published."""
        path = repo / "notes.py"
        path.write_text("#PRIVATE\nx = 1\n")
        commit(repo, "private for a while")
        path.write_text("x = 1\n")
        commit(repo, "unmarked now")
        done = check(repo, "base")
        assert done.returncode == 1
        assert "notes.py" in done.stderr

    def test_a_file_that_only_mentions_the_marker_is_not_flagged(self, repo):
        """Line-start, so prose about the convention is not itself private."""
        (repo / "doc.md").write_text("A file whose first line is #PRIVATE never ships.\n")
        commit(repo, "documentation")
        done = check(repo, "base")
        assert done.returncode == 0, done.stderr

    def test_an_indented_marker_still_counts(self, repo):
        (repo / "conf.yaml").write_text("block:\n  #PRIVATE\n  key: value\n")
        commit(repo, "indented marker")
        assert check(repo, "base").returncode == 1

    def test_a_binary_file_does_not_crash_the_check(self, repo):
        (repo / "blob.bin").write_bytes(bytes(range(256)) * 8)
        commit(repo, "a binary file")
        done = check(repo, "base")
        assert done.returncode == 0, done.stdout + done.stderr

    def test_a_file_already_on_the_target_is_not_reported(self, repo):
        """It is published already; the check is about what this range would add."""
        git(repo, "checkout", "-q", "base")
        (repo / "old_secret.py").write_text("#PRIVATE\nx = 1\n")
        commit(repo, "marked, but already on the target")
        git(repo, "checkout", "-q", "work")
        git(repo, "merge", "-q", "--no-edit", "base")
        (repo / "more.py").write_text("x = 2\n")
        commit(repo, "ordinary work")
        done = check(repo, "base")
        assert done.returncode == 0, done.stdout + done.stderr


class TestTheListing:

    def test_list_names_the_marked_files_of_the_working_tree(self, repo):
        (repo / "secret.py").write_text("#PRIVATE\nx = 1\n")
        (repo / "plain.py").write_text("x = 1\n")
        commit(repo, "one of each")
        done = check(repo, "--list")
        assert done.returncode == 0, done.stderr
        assert done.stdout.split() == ["secret.py"]


class TestItScalesWithTheBlobsNotTheCommits:

    def test_a_range_of_many_files_over_many_commits_finishes_quickly(self, repo):
        """Reading each path once per commit cost 211 x 385 git invocations on the real
        repository, twenty minutes for work that takes seconds. The timeout is the
        assertion: a reader that re-reads unchanged blobs cannot meet it.

        It also pins the deadlock this rewrite first walked into: the request list is
        handed to `git cat-file --batch` through a file, because writing it down a pipe
        while not yet reading git's answers makes writer and reader wait on each other.
        """
        for n in range(40):
            (repo / f"f{n:03d}.py").write_text(f"x = {n}\n")
        commit(repo, "forty files")
        for c in range(25):
            (repo / f"f{c:03d}.py").write_text(f"x = {c}  # touched in commit {c}\n")
            commit(repo, f"commit {c}")
        done = check(repo, "base")
        assert done.returncode == 0, done.stdout + done.stderr

    def test_a_marked_file_is_found_among_many(self, repo):
        for n in range(40):
            (repo / f"f{n:03d}.py").write_text(f"x = {n}\n")
        (repo / "f017.py").write_text("#PRIVATE\nx = 17\n")
        commit(repo, "forty files, one of them private")
        done = check(repo, "base")
        assert done.returncode == 1
        assert "f017.py" in done.stderr
        assert "f016.py" not in done.stderr
