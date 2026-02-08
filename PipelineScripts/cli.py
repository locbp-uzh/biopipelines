"""
Command-line entry points for BioPipelines.

Provides biopipelines-submit and biopipelines-run commands that work
from any directory by locating the repository root automatically.
"""

import os
import sys
import subprocess


def _get_repo_root():
    """Get the biopipelines repository root (parent of PipelineScripts/)."""
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _run_script(script_name):
    """Run a bash script from the repo root, passing through all CLI arguments."""
    repo_root = _get_repo_root()
    script_path = os.path.join(repo_root, script_name)

    if not os.path.exists(script_path):
        print(f"ERROR: {script_name} not found at {script_path}")
        sys.exit(1)

    # Run the script from the repo root so that os.getcwd() resolves correctly
    result = subprocess.run(
        ["bash", script_path] + sys.argv[1:],
        cwd=repo_root
    )
    sys.exit(result.returncode)


def submit():
    """Entry point for biopipelines-submit."""
    _run_script("submit")


def run():
    """Entry point for biopipelines-run."""
    _run_script("run")
