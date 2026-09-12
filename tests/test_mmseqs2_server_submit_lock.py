"""A server must hold the submission lock while its index loads.

The lock exists so a client does not submit a second server while one is on its
way. A client-started server inherits the lock the client took before
submitting, but a server started directly -- `Service()`, or by hand -- took
nothing, so for the several minutes its index loads a client sees neither
CPU_SERVER (written only once resident) nor CPU_SUBMITTING, and submits a
duplicate.

That is not hypothetical: in msa_short_001 the client probed at 10:54:25, eight
seconds into a 271-second load, and submitted a second 800 GB server. Both fed
the same queue so no work was lost, but one whole allocation was wasted.
"""

import os
import re
import shutil
import subprocess

import pytest

SERVERS = ["pipe_scripts/mmseqs2_server_cpu.sh",
           "pipe_scripts/mmseqs2_server_gpu.sh"]


def _read(repo_root, rel):
    return (repo_root / rel).read_text(encoding="utf-8")


def _trap_registration(body):
    """The top-level ``trap cleanup ...`` line, or None.

    Anchored at column 0 so the ``trap - EXIT ...`` disarm *inside* cleanup() (which is indented) is never mistaken for the registration.
    """
    return re.search(r"^trap\s+cleanup\s+(?P<signals>.+?)\s*$", body, re.M)


def _line_of(body, index):
    return body[:index].count("\n") + 1


@pytest.fixture
def repo_root():
    import pathlib
    return pathlib.Path(__file__).resolve().parent.parent


@pytest.mark.parametrize("script", SERVERS)
def test_server_claims_the_lock_before_loading(repo_root, script):
    """The claim must happen before the index load, not after."""
    body = _read(repo_root, script)

    claim = body.find('mkdir "${SUBMITTING_FILE}.lockdir"')
    assert claim != -1, "server never claims the submission lock"

    # The claim has to precede whatever makes the index resident, otherwise the
    # window it is meant to cover is already over.
    # The load is whatever call precedes the ready advertisement; anchoring on
    # that is robust to the two servers naming their load steps differently.
    load = body.find("> \"$SERVER_TIMESTAMP_FILE\"")
    assert load != -1, "could not locate the ready advertisement"
    assert claim < load, "lock is claimed after the server already advertised ready"


@pytest.mark.parametrize("script", SERVERS)
def test_claim_is_mkdir_not_o_excl(repo_root, script):
    """mkdir is the atomic primitive over NFS; O_EXCL file creation is not."""
    body = _read(repo_root, script)
    claim_line = [l for l in body.splitlines()
                  if 'mkdir "${SUBMITTING_FILE}.lockdir"' in l]
    assert claim_line, "no lockdir claim found"
    # Must be the condition of an if, so a failed claim is handled rather than
    # aborting the server under `set -e`.
    assert any(l.strip().startswith("if mkdir") for l in claim_line), \
        "claim must be guarded by `if` so losing the race is not fatal"


@pytest.mark.parametrize("script", SERVERS)
def test_ready_path_releases_the_lock(repo_root, script):
    """Going ready must clear both the marker and the lockdir."""
    body = _read(repo_root, script)
    # Anchor on the line that WRITES the timestamp, not the comment describing it.
    ready = body.find("> \"$SERVER_TIMESTAMP_FILE\"")
    assert ready != -1, "no ready-advertisement found"

    # cleanup() holds the same two lines, and is defined before the ready path so the trap can be armed before the claim, so searching from the ready advertisement onwards cannot match its copies by accident.
    cleanup_at = body.find("cleanup() {")
    assert cleanup_at != -1, "no cleanup() found"
    assert cleanup_at < ready, (
        "cleanup() must be defined before the ready path, so the trap can be "
        "armed before the lock is claimed")
    after = body[ready:]
    assert 'rmdir "${SUBMITTING_FILE}.lockdir"' in after, "ready path leaves the lockdir"
    assert 'rm -f "$SUBMITTING_FILE"' in after, "ready path leaves the marker"

    # A client probing between the two writes sees neither marker and submits a duplicate.
    assert ready < body.index('rmdir "${SUBMITTING_FILE}.lockdir"', ready), \
        "the lock is released before the server advertises ready"


@pytest.mark.parametrize("script", SERVERS)
def test_cleanup_releases_the_lock(repo_root, script):
    """A server that dies before ready must not block clients for the TTL.

    Without this the lock outlives the process by its 3 h staleness TTL, so one
    crashed server silently suppresses every client's server for an afternoon --
    strictly worse than the duplicate submission the lock prevents.
    """
    body = _read(repo_root, script)
    m = re.search(r"^cleanup\(\)\s*\{(.*?)^\}", body, re.S | re.M)
    assert m, "no cleanup() function found"
    fn = m.group(1)
    assert 'rm -f "$SUBMITTING_FILE"' in fn, "cleanup leaves the submission marker"
    assert 'rmdir "${SUBMITTING_FILE}.lockdir"' in fn, "cleanup leaves the lockdir"

    # A cleanup() nobody registers as a trap releases nothing; asserting only on its body is what let a lock leak through CI.
    trap = _trap_registration(body)
    assert trap, "cleanup() is never registered with `trap`, so it never runs"
    signals = set(trap.group("signals").split())
    for signal in ("EXIT", "SIGINT", "SIGTERM"):
        assert signal in signals, (
            f"trap must cover {signal} or that exit path leaks the lock; got {sorted(signals)}"
        )


@pytest.mark.parametrize("script", SERVERS)
def test_cleanup_trap_is_armed_before_the_lock_is_claimed(repo_root, script):
    """The trap must cover the window in which the lock is actually held.

    The lock is claimed before the index load and released once the server is ready. If the trap is only armed after that release, every failure during the load -- the entire window the lock exists to cover -- leaks it.
    """
    body = _read(repo_root, script)
    claim = body.find('if mkdir "${SUBMITTING_FILE}.lockdir"')
    assert claim != -1, "server never claims the submission lock"
    trap = _trap_registration(body)
    assert trap, "cleanup() is never registered with `trap`"

    assert trap.start() < claim, (
        f"trap armed at line {_line_of(body, trap.start())} but the lock is claimed at "
        f"line {_line_of(body, claim)}: a death during the load leaks the lockdir"
    )


def test_lock_ttl_outlives_a_realistic_load(repo_root):
    """The TTL must exceed how long an index load actually takes.

    Measured: 271 s to lock uniref30 (220 GB) + envdb (527 GB), plus SLURM queue
    time ahead of it. A TTL near that would let a client steal the lock from a
    server that is still legitimately loading.
    """
    src = _read(repo_root, "pipe_scripts/pipe_mmseqs2_sequences.py")
    m = re.search(r"SUBMIT_LOCK_TTL_SECONDS\s*=\s*([^\n#]+)", src)
    assert m, "SUBMIT_LOCK_TTL_SECONDS not found"
    ttl = eval(m.group(1).strip(), {"__builtins__": {}})
    assert ttl >= 1800, f"TTL {ttl}s is too short for a ~271s load plus queue time"


# ── the trap, executed rather than grepped ───────────────────────────────────
#
# The greps above pin the textual ordering; these run the real scripts, driving a server to the point where it holds the lock and then making it die inside the window the lock covers. This is the check that survives any future re-layout.

bash_required = pytest.mark.skipif(shutil.which("bash") is None,
                                   reason="needs bash to execute the trap")


def _stub(path, body="exit 0\n"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/bash\n" + body, encoding="utf-8", newline="\n")
    path.chmod(0o755)


def _server_env(tmp_path, stub_dir):
    """Folders and PATH that let a server script run with no cluster present."""
    shared = tmp_path / "shared"
    (shared / "bin").mkdir(parents=True)
    for name in ("db", "data", "colabfold/colabfold-conda/bin"):
        (tmp_path / name).mkdir(parents=True, exist_ok=True)
    for tool in ("colabfold_search", "mmseqs"):
        _stub(tmp_path / "colabfold/colabfold-conda/bin" / tool, 'echo stub\n')
    _stub(stub_dir / "biopipelines-config", "exit 1\n")
    env = dict(os.environ)
    env.update({
        "PATH": stub_dir.as_posix() + os.pathsep + env.get("PATH", ""),
        "MMSEQS2_SHARED_FOLDER": shared.as_posix(),
        "MMSEQS2_DB_DIR": (tmp_path / "db").as_posix(),
        "BIOPIPELINES_DATA_DIR": (tmp_path / "data").as_posix(),
        "COLABFOLD_DIR": (tmp_path / "colabfold").as_posix(),
        "MMSEQS2_USE_SHM": "0",
        # Set so the script does not need /proc/meminfo or nproc.
        "SLURM_MEM_PER_NODE": "1024",
        "SLURM_CPUS_PER_TASK": "2",
    })
    return shared, env


@bash_required
def test_cpu_server_releases_the_lock_when_it_dies_during_the_load(repo_root, tmp_path):
    """A death inside the index load must not strand the lockdir for the TTL.

    The lock is claimed before the load and released only once the server is ready, so the load is the whole window the lock exists to cover. Here the .idx files are absent, which is exactly the ``exit 1`` inside ``lock_one_file`` -- a death after the claim and before the release.
    """
    stub_dir = tmp_path / "stub"
    shared, env = _server_env(tmp_path, stub_dir)
    # A pre-placed vmtouch means ensure_vmtouch does not try to build one.
    _stub(shared / "bin" / "vmtouch")

    proc = subprocess.run(
        ["bash", (repo_root / "pipe_scripts/mmseqs2_server_cpu.sh").as_posix()],
        env=env, cwd=tmp_path.as_posix(), capture_output=True, text=True, timeout=180)

    log = proc.stdout + proc.stderr
    assert "Claimed submission lock" in log, f"never reached the claim:\n{log}"
    assert "index not found" in log, f"never reached the load:\n{log}"
    assert not (shared / "CPU_SUBMITTING.lockdir").exists(), (
        f"lockdir survived a death during the load -- every client is "
        f"suppressed until the TTL expires:\n{log}")
    assert not (shared / "CPU_SUBMITTING").exists(), f"marker survived:\n{log}"
    # The failure must still reach SLURM as a failure, not as COMPLETED.
    assert proc.returncode != 0, "cleanup masked the failure with exit 0"


@bash_required
def test_gpu_server_releases_the_lock_when_it_dies_during_warm_up(repo_root, tmp_path):
    """Same window on the GPU side: install + DB warm-up, before the ready path."""
    stub_dir = tmp_path / "stub"
    shared, env = _server_env(tmp_path, stub_dir)
    _stub(stub_dir / "nvidia-smi", 'echo "StubGPU, 1, 2, 3"\n')
    # check_mmseqs_installation downloads when mmseqs is absent, so a failing wget is the `set -e` abort that test_mmseqs2_server_lock_release only simulates.
    _stub(stub_dir / "wget", 'echo "wget: stub failure" >&2\nexit 1\n')
    env["MMSEQS2_DIR"] = (tmp_path / "mmseqs2").as_posix()
    (tmp_path / "mmseqs2").mkdir()

    proc = subprocess.run(
        ["bash", (repo_root / "pipe_scripts/mmseqs2_server_gpu.sh").as_posix()],
        env=env, cwd=tmp_path.as_posix(), capture_output=True, text=True, timeout=180)

    log = proc.stdout + proc.stderr
    assert "Claimed submission lock" in log, f"never reached the claim:\n{log}"
    assert "downloading" in log, f"never reached the install step:\n{log}"
    assert not (shared / "GPU_SUBMITTING.lockdir").exists(), (
        f"lockdir survived a death during warm-up:\n{log}")
    assert not (shared / "GPU_SUBMITTING").exists(), f"marker survived:\n{log}"
    assert proc.returncode != 0, "cleanup masked the failure with exit 0"


@bash_required
@pytest.mark.parametrize("script", SERVERS)
def test_a_server_that_never_went_ready_leaves_a_peers_advertisement_alone(
        repo_root, tmp_path, script):
    """cleanup() must not delete a CPU_SERVER/GPU_SERVER it did not write.

    Arming the trap before the load means cleanup() now runs on deaths that happen while a *different* server is already serving, and removing that peer's advertisement would recreate the very "client sees no server" state the lock exists to prevent -- so the removal is gated on having advertised.
    """
    body = _read(repo_root, script)
    m = re.search(r"^cleanup\(\)\s*\{(.*?)^\}", body, re.S | re.M)
    assert m, "no cleanup() found"
    fn = m.group(1)
    assert 'rm -f "$SERVER_TIMESTAMP_FILE"' in fn, "cleanup no longer clears its own advertisement"
    guard = re.search(r'if \[\[ "\$SERVER_READY" == "1" \]\]; then\s*\n\s*rm -f "\$SERVER_TIMESTAMP_FILE"',
                      fn)
    assert guard, "the advertisement removal is not gated on SERVER_READY"
    # And the flag has to actually be raised where the advertisement is written.
    ready = body.find('> "$SERVER_TIMESTAMP_FILE"')
    assert body.find("SERVER_READY=1", ready) - ready < 120, \
        "SERVER_READY is not set on the ready path"
