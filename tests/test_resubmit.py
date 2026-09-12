"""`resubmit` had no tests at all, and it is the script people reach for when a
run has already failed.

`submit` resolves the `<JOBID_BATCH_NNN>` placeholders with `sed -i`, so the
script on disk permanently names the job ids of the run that produced it.
Handing that straight back to the scheduler is wrong in every case resubmit
exists for: a parent that failed or was cancelled leaves the job pending on a
dependency that can never be satisfied, and a parent whose id has aged out is
either rejected or treated as already satisfied depending on site config.

The scheduler binaries are stubbed by a shell script on PATH that records its
argv and its stdin, so these assert what resubmit actually hands over.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
RESUBMIT = REPO_ROOT / "resubmit"

# Verified under git-bash on Windows as well as Linux CI, so this skips only
# where there is genuinely no bash, rather than on platform.
pytestmark = pytest.mark.skipif(
    shutil.which("bash") is None,
    reason="resubmit is a POSIX shell script and needs a real bash to drive",
)

SLURM_SCRIPT = """#!/usr/bin/bash
#SBATCH --job-name=old
#SBATCH --dependency=afterok:12345
#SBATCH --mem=4GB
#SBATCH -d afterany:99999
echo body
"""

LSF_SCRIPT = """#!/usr/bin/bash
#BSUB -J old
#BSUB -w done(4242)
#BSUB -M 4096
echo body
"""

# -W carries other PBS attributes too; only depend= may be removed.
PBS_SCRIPT = """#!/usr/bin/bash
#PBS -N old
#PBS -W depend=afterok:777.server
#PBS -W group_list=chem
#PBS -l walltime=01:00:00
echo body
"""


def _stub_scheduler(bin_dir: Path, name: str) -> Path:
    """A fake sbatch/bsub/qsub that records argv and stdin, then succeeds."""
    record = bin_dir / f"{name}.record"
    stub = bin_dir / name
    stub.write_text(
        "#!/usr/bin/env bash\n"
        f'echo "ARGV: $*" > "{record}"\n'
        f'echo "--- STDIN ---" >> "{record}"\n'
        f'cat >> "{record}"\n'
        "exit 0\n",
        encoding="utf-8",
    )
    stub.chmod(0o755)
    return record


def _run(tmp_path, script_text, scheduler, *args):
    run_dir = tmp_path / "JobFolder" / "RunTime"
    run_dir.mkdir(parents=True)
    job_script = run_dir / "slurm_batch1.sh"
    job_script.write_text(script_text, encoding="utf-8")

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    record = _stub_scheduler(bin_dir, scheduler)

    env = dict(os.environ, PATH=f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    proc = subprocess.run(
        ["bash", str(RESUBMIT), *args, str(job_script)],
        capture_output=True, text=True, env=env,
    )
    submitted = record.read_text(encoding="utf-8") if record.exists() else ""
    return proc, submitted, job_script


def _stdin_of(submitted: str) -> str:
    return submitted.split("--- STDIN ---\n", 1)[1] if "--- STDIN ---" in submitted else ""


def test_slurm_dependency_directives_are_stripped(tmp_path):
    """Both spellings: `--dependency=` and the short `-d`."""
    proc, submitted, _ = _run(tmp_path, SLURM_SCRIPT, "sbatch")

    assert proc.returncode == 0, proc.stderr
    body = _stdin_of(submitted)
    assert "--dependency" not in body
    assert "-d afterany" not in body
    # Everything else survives.
    assert "#SBATCH --mem=4GB" in body
    assert "echo body" in body


def test_removed_directives_are_echoed(tmp_path):
    """A dependency is never dropped silently."""
    proc, _, _ = _run(tmp_path, SLURM_SCRIPT, "sbatch")

    assert "Dependencies: STRIPPED" in proc.stdout
    assert "afterok:12345" in proc.stdout
    assert "afterany:99999" in proc.stdout


def test_keep_dependencies_retains_them(tmp_path):
    """The one case that wants them: a parent batch still queued."""
    proc, submitted, _ = _run(tmp_path, SLURM_SCRIPT, "sbatch", "--keep-dependencies")

    body = _stdin_of(submitted)
    assert "#SBATCH --dependency=afterok:12345" in body
    assert "#SBATCH -d afterany:99999" in body
    assert "Dependencies: KEPT" in proc.stdout


def test_the_on_disk_script_is_never_modified(tmp_path):
    """Fed through stdin, so the record of what the original run submitted stands."""
    _, _, job_script = _run(tmp_path, SLURM_SCRIPT, "sbatch")

    assert job_script.read_text(encoding="utf-8") == SLURM_SCRIPT


def test_lsf_w_directive_is_stripped(tmp_path):
    proc, submitted, _ = _run(tmp_path, LSF_SCRIPT, "bsub")

    assert proc.returncode == 0, proc.stderr
    body = _stdin_of(submitted)
    assert "-w done(4242)" not in body
    assert "#BSUB -M 4096" in body


def test_pbs_strips_depend_but_keeps_other_w_attributes(tmp_path):
    """`-W` also carries group_list and umask; only depend= is a dependency."""
    proc, submitted, _ = _run(tmp_path, PBS_SCRIPT, "qsub")

    assert proc.returncode == 0, proc.stderr
    body = _stdin_of(submitted)
    assert "depend=afterok" not in body
    assert "#PBS -W group_list=chem" in body, "an unrelated -W attribute was removed"
    assert "#PBS -l walltime=01:00:00" in body


def test_a_script_with_no_dependencies_says_so(tmp_path):
    clean = "#!/usr/bin/bash\n#SBATCH --job-name=old\n#SBATCH --mem=4GB\necho body\n"
    proc, submitted, _ = _run(tmp_path, clean, "sbatch")

    assert "Dependencies: none in this script." in proc.stdout
    assert "echo body" in _stdin_of(submitted)


def test_a_script_with_no_directives_is_refused(tmp_path):
    """Without a directive there is no way to tell which scheduler this is."""
    proc, _, _ = _run(tmp_path, "#!/usr/bin/bash\necho body\n", "sbatch")

    assert proc.returncode != 0
    assert "cannot tell the scheduler" in proc.stdout + proc.stderr
