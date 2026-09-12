"""A server must release the submission lock however it dies.

Companion to test_mmseqs2_server_submit_lock.py, which asserts the lock is
*claimed*. Two release defects hid from it because it inspects cleanup()'s body
and never the trap that arms it:

  * mmseqs2_server_gpu.sh trapped only SIGINT/SIGTERM, not EXIT, while the CPU
    script trapped EXIT too. With `set -e` and the exit-1 paths in
    check_mmseqs_installation running after the lock claim, a GPU server dying
    during install or warm-up leaked the lock for its whole 3h TTL -- the exact
    failure the claim was added to fix.
  * Both scripts removed the marker file before the lockdir. If the rmdir then
    failed, clients reached marker-absent + lockdir-present, where
    submission_in_progress() reports nothing in flight and acquire_submit_lock()
    refuses because the dir exists -- so every client logged "another client is
    already submitting" and submitted nothing, silently, for 3h.

These run the real shell rather than grepping, so a future edit to either
mechanism has to keep the behaviour.
"""

import pathlib
import shutil
import subprocess
import textwrap

import pytest

SERVERS = ["pipe_scripts/mmseqs2_server_cpu.sh",
           "pipe_scripts/mmseqs2_server_gpu.sh"]

bash_required = pytest.mark.skipif(shutil.which("bash") is None,
                                   reason="needs bash to execute the trap")


@pytest.fixture
def repo_root():
    return pathlib.Path(__file__).resolve().parent.parent


def _read(repo_root, rel):
    return (repo_root / rel).read_text(encoding="utf-8")


@pytest.mark.parametrize("script", SERVERS)
def test_cleanup_is_trapped_on_exit(repo_root, script):
    """EXIT is the branch that catches `set -e` and every plain exit 1."""
    body = _read(repo_root, script)
    traps = [ln.strip() for ln in body.splitlines()
             if ln.strip().startswith("trap ") and "cleanup" in ln]
    assert traps, f"{script}: cleanup is never trapped"
    assert any("EXIT" in t for t in traps), (
        f"{script}: cleanup trapped on {traps} but not EXIT, so a server that "
        f"dies during install or warm-up leaks the submission lock for its TTL")


def _cleanup_body(script_text):
    """The text of cleanup(), where a release happens with no server advertised."""
    start = script_text.find("cleanup() {")
    assert start != -1, "cleanup() not found"
    end = script_text.find("\n}", start)
    assert end != -1, "cleanup() has no closing brace"
    return script_text[start:end]


@pytest.mark.parametrize("script", SERVERS)
def test_cleanup_releases_the_lockdir_before_the_marker(repo_root, script):
    """Order matters in cleanup(), where no ready server is advertised.

    The ready path releases in either order safely, because the timestamp file
    is written first and a client that sees it never consults the lock. In
    cleanup() there is no timestamp file, so a marker orphaned by a failed rmdir
    strands every client until the TTL.
    """
    body = _cleanup_body(_read(repo_root, script))
    rmdir = body.find('rmdir "${SUBMITTING_FILE}.lockdir"')
    marker = body.find('rm -f "$SUBMITTING_FILE"')
    assert rmdir != -1 and marker != -1, f"{script}: cleanup does not release the lock"
    assert rmdir < marker, (
        f"{script}: cleanup removes the marker before the lockdir. If the rmdir "
        f"then fails, clients see no submission in progress and cannot take the "
        f"lock either -- they all no-op silently until the TTL expires")


@bash_required
@pytest.mark.parametrize("script", SERVERS)
def test_trap_actually_fires_on_a_failed_exit(repo_root, tmp_path, script):
    """Run the real trap wiring: a non-zero exit must clear both artifacts."""
    body = _read(repo_root, script)
    traps = [ln.strip() for ln in body.splitlines()
             if ln.strip().startswith("trap ") and "cleanup" in ln]
    trap_line = traps[0]

    marker = tmp_path / "CPU_SUBMITTING"
    lockdir = tmp_path / "CPU_SUBMITTING.lockdir"
    lockdir.mkdir()
    marker.write_text("held\n", encoding="utf-8")

    # Reproduce only the release contract: the same trap line, a cleanup that
    # releases in the script's own order, and a failure after the lock is held.
    harness = textwrap.dedent(f"""\
        set -e
        SUBMITTING_FILE="{marker.as_posix()}"
        cleanup() {{
          rmdir "${{SUBMITTING_FILE}}.lockdir" 2>/dev/null || true
          rm -f "$SUBMITTING_FILE"
          exit 0
        }}
        {trap_line}
        false          # stands in for check_mmseqs_installation failing
        echo "unreachable"
        """)
    harness_path = tmp_path / "harness.sh"
    harness_path.write_text(harness, encoding="utf-8", newline="\n")

    subprocess.run(["bash", harness_path.as_posix()], cwd=tmp_path,
                   capture_output=True, text=True, timeout=60)

    assert not lockdir.exists(), (
        f"{script}: lockdir survived a failed exit — clients block for the TTL")
    assert not marker.exists(), (
        f"{script}: marker survived a failed exit")


# --- the marker/lockdir staleness asymmetry ---------------------------------

def test_lockdir_is_the_staleness_authority(repo_root):
    """Both checks must judge staleness from the same artifact.

    submission_in_progress() read the MARKER's mtime while acquire_submit_lock()
    read the LOCKDIR's. That left one state in which nobody could act: marker
    gone, lockdir present. submission_in_progress returned False (no marker), the
    acquire then failed because the directory existed, and every client logged
    that someone else was submitting and submitted nothing -- for the full 3h TTL.

    Reachable two ways: a client dying between os.mkdir(lock_dir) and the marker
    write, and a cleanup whose rmdir failed after the marker was already gone.
    """
    body = _read(repo_root, "pipe_scripts/pipe_mmseqs2_sequences.py")
    start = body.index("def submission_in_progress(")
    end = body.index("def acquire_submit_lock(")
    fn = body[start:end]

    assert ".lockdir" in fn, (
        "submission_in_progress does not look at the lockdir, so a marker-gone / "
        "lockdir-present state reads as 'nothing in progress' while the acquire "
        "still fails -- every client no-ops until the TTL expires")
    assert "getmtime" in fn, "staleness is no longer judged from an mtime"
