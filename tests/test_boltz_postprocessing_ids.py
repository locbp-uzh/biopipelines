# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Boltz post-processing maps a prediction folder to its id by exact match, or refuses it.

The config step writes one `<id>.yaml` per expected id and Boltz names its prediction folder after that file, so a folder is either exactly one of the expected ids or it is not one of this step's predictions (a leftover from an earlier run in the same execution folder). The script also accepted a folder whose name merely started or ended with `<expected id>+`, which attributed a foreign prediction to an expected id and let it overwrite that id's structure, scores and alignment without a word.
"""

import os
import subprocess
import sys

import pandas as pd
import pytest

from tests._helpers import BOLTZ_POST, boltz_prediction


def _postprocess(root, ids):
    for folder in ("_configuration", "structures", "tables", "sequences", "msas", "_extras", "pae"):
        (root / folder).mkdir(exist_ok=True)
    sequence_ids = root / "_configuration" / "sequence_ids.csv"
    pd.DataFrame({"id": ids}).to_csv(sequence_ids, index=False)
    args = [BOLTZ_POST, root / "execution", root, sequence_ids,
            "--structures-map", root / "structures" / "structures_map.csv",
            "--structures-folder", root / "structures",
            "--msas-folder", root / "msas",
            "--confidence-csv", root / "tables" / "confidence.csv",
            "--affinity-csv", root / "tables" / "affinity.csv",
            "--sequences-csv", root / "sequences" / "sequences.csv",
            "--msas-csv", root / "msas" / "msas_map.csv",
            "--scores-info", root / "_extras" / "scores_info.txt"]
    return subprocess.run([sys.executable, *[str(a) for a in args]], capture_output=True, text=True)


def test_every_exact_folder_maps_to_its_own_id(tmp_path):
    for config_id in ("p1+l1", "p1+l2"):
        boltz_prediction(tmp_path, config_id)

    done = _postprocess(tmp_path, ["p1+l1", "p1+l2"])

    assert done.returncode == 0, done.stdout + done.stderr
    assert sorted(os.listdir(tmp_path / "structures")) == ["p1+l1.pdb", "p1+l2.pdb", "structures_map.csv"]


@pytest.mark.parametrize("expected, foreign", [
    (["p1+l1", "p1+l2"], "p1+l1+x"),
    (["A", "B"], "X+A"),
])
def test_a_folder_that_only_contains_an_expected_id_is_refused(tmp_path, expected, foreign):
    for config_id in (*expected, foreign):
        boltz_prediction(tmp_path, config_id)

    done = _postprocess(tmp_path, expected)

    assert done.returncode != 0
    assert f"Cannot map config_id '{foreign}'" in done.stdout + done.stderr
