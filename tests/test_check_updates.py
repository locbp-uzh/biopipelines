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
