# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""A step's own exit status decides its marker.

The completion check judged declared paths by existence alone and wrote its marker before the footer looked at the command's status; `pipeline.py`'s guard then found a COMPLETED marker and declined to write FAILED. A tool that died after creating its declared files therefore left the run green.
"""

import json
import os
import subprocess
import sys

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CHECK = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_check_completion.py")


def _step(tmp_path, table_rows, ids=("a",)):
    """A step folder whose one declared table holds `table_rows` data rows."""
    folder = os.path.join(str(tmp_path), "007_Prodigy")
    os.makedirs(os.path.join(folder, "tables"), exist_ok=True)
    table = os.path.join(folder, "tables", "affinity.csv")
    with open(table, "w") as handle:
        handle.write("id,affinity\n")
        for i in range(table_rows):
            handle.write(f"row{i},-9.1\n")
    expected = {
        "structures": {"ids": list(ids), "files": [], "map_table": "", "format": "csv"},
        "tables": {"affinity": {"name": "affinity", "path": table,
                                "columns": ["id", "affinity"], "description": ""}},
        "output_folder": folder,
    }
    manifest = os.path.join(folder, ".expected_outputs.json")
    json.dump(expected, open(manifest, "w"))
    return folder, manifest


def _run(folder, manifest, main_rc):
    return subprocess.run(
        [sys.executable, CHECK, folder, "Prodigy", manifest, "--main-rc", str(main_rc)],
        capture_output=True, text=True)


def _markers(tmp_path):
    return sorted(f for f in os.listdir(str(tmp_path)) if f.startswith("007_Prodigy_"))


def test_a_nonzero_exit_is_failed_however_complete_the_outputs_look(record_case, tmp_path):
    folder, manifest = _step(tmp_path, table_rows=3)
    result = _run(folder, manifest, main_rc=1)
    markers = _markers(tmp_path)
    record_case(input="rc=1, all declared outputs present",
                expected=["007_Prodigy_FAILED"], actual=markers)
    assert markers == ["007_Prodigy_FAILED"]
    assert result.returncode == 1


def test_a_zero_exit_with_rows_is_completed(record_case, tmp_path):
    folder, manifest = _step(tmp_path, table_rows=3)
    result = _run(folder, manifest, main_rc=0)
    markers = _markers(tmp_path)
    record_case(input="rc=0, 3 rows", expected=["007_Prodigy_COMPLETED"], actual=markers)
    assert markers == ["007_Prodigy_COMPLETED"]
    assert result.returncode == 0


def test_the_generated_footer_passes_the_step_status_through(record_case):
    """Without this argument the check cannot know the command died."""
    import inspect
    from biopipelines.base_config import BaseConfig
    source = inspect.getsource(BaseConfig.generate_completion_check_footer)
    record_case(input="generate_completion_check_footer source",
                expected="--main-rc \"${BP_MAIN_RC:-0}\"",
                actual="--main-rc" in source)
    assert '--main-rc "${{BP_MAIN_RC:-0}}"' in source
