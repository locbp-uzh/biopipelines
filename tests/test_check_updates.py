"""`check-updates` had no tests, and it runs on every first submit of the day.

Its whole contract is negative: it must never fetch, never write to the
worktree, never need python, and never exit non-zero — a notice that blocks a
job is worse than no notice. Those are exactly the properties that rot
silently, so they are asserted here rather than described in the header.

`git` is stubbed on PATH where the test needs to control the answer, so these
neither reach GitHub nor depend on what the real release branch happens to be.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
CHECK_UPDATES = REPO_ROOT / "check-updates"

BASH = shutil.which("bash")

pytestmark = pytest.mark.skipif(
    BASH is None,
    reason="check-updates is a POSIX shell script and needs a real bash to drive",
)


def _fake_bin(tmp_path: Path, **scripts: str) -> Path:
    """A PATH directory holding stub executables."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir(exist_ok=True)
    for name, body in scripts.items():
        path = bin_dir / name
        path.write_text("#!/usr/bin/env bash\n" + body, encoding="utf-8")
        path.chmod(0o755)
    return bin_dir


def _run(tmp_path, git_body: str, *, drop_network_tools: bool = True):
    """Run check-updates with a stubbed git, and no downloader by default."""
    bin_dir = _fake_bin(tmp_path, git=git_body)
    if drop_network_tools:
        # An empty curl/wget keeps the optional feed call from reaching GitHub.
        _fake_bin(tmp_path, curl="exit 1\n", wget="exit 1\n")

    # A minimal PATH: the stubs, plus what bash itself needs.
    env = dict(os.environ, PATH=f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    return subprocess.run(
        [BASH, str(CHECK_UPDATES)], capture_output=True, text=True, env=env, cwd=str(tmp_path),
    )


def _git_dir(tmp_path) -> str:
    """A per-test stand-in for the repo's .git directory."""
    d = tmp_path / "fakegit"
    d.mkdir(exist_ok=True)
    return str(d).replace("\\", "/")


# A stub git that answers ls-remote with a sha and ancestry with the given code.
#
# `rev-parse` has to answer too. check-updates derives its once-per-day sentinel
# from `rev-parse --absolute-git-dir`; a stub that stays silent yields an empty
# GIT_DIR, so the sentinel becomes "/.bp-check-updates-day" -- which a ROOT
# shell (CI) can really create, after which every later case in the same job
# hits the throttle and exits silently. Pointing it at the test's own tmp_path
# keeps the cases isolated from each other.
def _git_stub(ancestor_exit: int, git_dir: str, sha: str = "a" * 40) -> str:
    return (
        'for a in "$@"; do\n'
        '  case "$a" in\n'
        f'    ls-remote) echo "{sha}\trefs/heads/main"; exit 0 ;;\n'
        f'    merge-base) exit {ancestor_exit} ;;\n'
        f'    rev-parse) echo "{git_dir}"; exit 0 ;;\n'
        '    fetch|pull) echo "FETCHED" >&2; exit 1 ;;\n'
        "  esac\n"
        "done\n"
        "exit 0\n"
    )


def test_it_never_fetches(tmp_path):
    """A fetch would write to the worktree of a job that is about to start."""
    proc = _run(tmp_path, _git_stub(ancestor_exit=0, git_dir=_git_dir(tmp_path)))

    assert "FETCHED" not in proc.stderr, "check-updates invoked git fetch"
    assert "FETCHED" not in proc.stdout


def test_it_exits_zero_even_when_everything_fails(tmp_path):
    """Any failure must be silent: the notice must never block a job."""
    proc = _run(tmp_path, "exit 127\n")

    assert proc.returncode == 0, proc.stderr


def test_it_exits_zero_when_git_is_absent(tmp_path):
    """No git at all is the ordinary case on a bare login shell."""
    bin_dir = _fake_bin(tmp_path, curl="exit 1\n", wget="exit 1\n")
    env = dict(os.environ, PATH=str(bin_dir))
    proc = subprocess.run(
        [BASH, str(CHECK_UPDATES)], capture_output=True, text=True, env=env, cwd=str(tmp_path),
    )

    assert proc.returncode == 0


def test_up_to_date_checkout_is_silent(tmp_path):
    """merge-base --is-ancestor succeeding means the release line is behind us."""
    proc = _run(tmp_path, _git_stub(ancestor_exit=0, git_dir=_git_dir(tmp_path)))

    assert proc.returncode == 0
    assert "behind" not in proc.stdout.lower()


def test_a_behind_checkout_is_reported(tmp_path):
    """Ancestry failing means main's HEAD is not in our history."""
    proc = _run(tmp_path, _git_stub(ancestor_exit=1, git_dir=_git_dir(tmp_path)))

    assert proc.returncode == 0
    assert proc.stdout.strip(), "a checkout behind the release line printed nothing"


def test_an_empty_ls_remote_is_silent(tmp_path):
    """Network failure gives an empty sha; that is 'cannot tell', not 'behind'."""
    # rev-parse must answer here too, or this passes because the run was
    # throttled rather than because the empty sha was handled.
    stub = (
        'for a in "$@"; do\n'
        '  case "$a" in\n'
        '    ls-remote) exit 0 ;;\n'
        f'    rev-parse) echo "{_git_dir(tmp_path)}"; exit 0 ;;\n'
        "  esac\n"
        "done\n"
        "exit 0\n"
    )
    proc = _run(tmp_path, stub)

    assert proc.returncode == 0
    assert "behind" not in proc.stdout.lower()
    # It reached the network step; it just could not tell.
    assert "Checking for BioPipelines updates" in proc.stdout


def test_it_does_not_require_python(tmp_path):
    """It runs from ./submit before any environment is activated.

    Comments are stripped first: the header explains at length that it must not
    need python or jq, and matching that would be the test failing on its own
    rationale.
    """
    code = "\n".join(
        line for line in CHECK_UPDATES.read_text(encoding="utf-8").splitlines()
        if not line.lstrip().startswith("#")
    )

    for token in ("python", "python3", "jq"):
        assert token not in code, (
            f"check-updates invokes {token}, which is not available where it runs"
        )


def test_the_timeout_call_is_guarded(tmp_path):
    """`timeout` is GNU coreutils and absent from a stock macOS."""
    source = CHECK_UPDATES.read_text(encoding="utf-8")

    assert "command -v timeout" in source, "timeout is used without checking it exists"
    assert "$TIMEOUT_CMD" in source


def test_it_works_when_timeout_is_missing(tmp_path):
    """With no `timeout` on PATH the ls-remote still runs, just unbounded."""
    bin_dir = _fake_bin(
        tmp_path,
        git=_git_stub(ancestor_exit=1, git_dir=_git_dir(tmp_path)),
        curl="exit 1\n",
        wget="exit 1\n",
        # `command -v timeout` must fail, so PATH below is only our stub dir
        # plus the interpreter's own; a stub named `timeout` would defeat that.
    )
    env = dict(os.environ, PATH=str(bin_dir))
    proc = subprocess.run(
        [BASH, str(CHECK_UPDATES)], capture_output=True, text=True, env=env, cwd=str(tmp_path),
    )

    assert proc.returncode == 0
    assert "timeout: command not found" not in proc.stderr


def _git_stub_with_origin(origin: str, git_dir: str, asked_file: str,
                          ancestor_exit: int = 1, sha: str = "b" * 40) -> str:
    """Like `_git_stub`, but answers `remote get-url origin` and records the URL that

    `ls-remote` was actually asked for. It records to a file rather than stderr: the
    script sends ls-remote's stderr to /dev/null, so a stderr probe sees nothing.
    """
    return (
        'for a in "$@"; do\n'
        '  case "$a" in\n'
        f'    remote) echo "{origin}"; exit 0 ;;\n'
        f'    ls-remote) echo "$4" > "{asked_file}"; echo "{sha}\trefs/heads/main"; exit 0 ;;\n'
        f'    merge-base) exit {ancestor_exit} ;;\n'
        f'    rev-parse) echo "{git_dir}"; exit 0 ;;\n'
        "  esac\n"
        "done\n"
        "exit 0\n"
    )


def _asked(tmp_path) -> str:
    """Where the stub records the URL, in a form bash can write to."""
    return str(tmp_path / "asked.txt").replace("\\", "/")


class TestWhichRemoteIsAsked:
    """The lab's GitLab `main` is the line an internal clone tracks; GitHub's is a mirror.

    Asking GitHub for every checkout answered the wrong question internally: GitHub's HEAD is
    an ancestor of an internal clone's, so the check stayed silent while the line that had
    actually moved went unmentioned.
    """

    def test_the_checkouts_own_origin_is_queried(self, tmp_path):
        origin = "https://gitlab.uzh.ch/locbp/public/biopipelines-locbp.git"
        asked = Path(_asked(tmp_path))
        proc = _run(tmp_path, _git_stub_with_origin(origin, _git_dir(tmp_path), str(asked)))
        assert asked.read_text().strip() == origin, "ls-remote did not use the configured origin"
        assert "github.com" not in proc.stdout

    def test_the_notice_names_that_remote(self, tmp_path):
        origin = "https://gitlab.uzh.ch/locbp/public/biopipelines-locbp.git"
        proc = _run(tmp_path, _git_stub_with_origin(origin, _git_dir(tmp_path), _asked(tmp_path)))
        assert "gitlab.uzh.ch/locbp/public/biopipelines-locbp" in proc.stdout
        assert ".git" not in proc.stdout.split("branch")[0].split("available")[-1]

    def test_a_checkout_with_no_origin_falls_back_to_github(self, tmp_path):
        """An export or a tarball still has a release line worth checking."""
        proc = _run(tmp_path, _git_stub(ancestor_exit=1, git_dir=_git_dir(tmp_path)))
        assert "github.com/locbp-uzh/biopipelines" in proc.stdout

    def test_credentials_in_a_remote_url_never_reach_the_output(self, tmp_path):
        """A token pasted into a remote would otherwise land in every job's log."""
        origin = "https://oauth2:glpat-SECRETTOKEN@gitlab.uzh.ch/locbp/public/bp.git"
        proc = _run(tmp_path, _git_stub_with_origin(origin, _git_dir(tmp_path), _asked(tmp_path)))
        assert "glpat-SECRETTOKEN" not in proc.stdout
        assert "oauth2" not in proc.stdout
        assert "gitlab.uzh.ch/locbp/public/bp" in proc.stdout

    def test_a_password_containing_an_at_sign_is_stripped_whole(self, tmp_path):
        """Stripping to the first @ printed the rest of the password in front of the host."""
        origin = "https://user:pa@ssw0rd@gitlab.uzh.ch/locbp/public/bp.git"
        proc = _run(tmp_path, _git_stub_with_origin(origin, _git_dir(tmp_path), _asked(tmp_path)))
        assert "ssw0rd" not in proc.stdout
        assert "gitlab.uzh.ch/locbp/public/bp" in proc.stdout

    def test_an_ssh_remote_is_shown_as_host_and_path(self, tmp_path):
        origin = "git@gitlab.uzh.ch:locbp/public/biopipelines-locbp.git"
        proc = _run(tmp_path, _git_stub_with_origin(origin, _git_dir(tmp_path), _asked(tmp_path)))
        assert "gitlab.uzh.ch/locbp/public/biopipelines-locbp" in proc.stdout


def test_it_never_waits_for_a_password(tmp_path):
    """An internal remote can prompt, and a prompt inside ./submit hangs the job.

    `timeout` bounds this only where coreutils exists, so refusing to prompt is the real guard.
    """
    source = CHECK_UPDATES.read_text(encoding="utf-8")
    assert "GIT_TERMINAL_PROMPT=0" in source
    assert "BatchMode=yes" in source
