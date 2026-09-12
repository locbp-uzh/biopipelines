# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""fpocket writes every pocket metric on a tab-indented line, so a parser that splits the raw line on its first tab reads an empty key and drops the lot. The fixture below is copied verbatim from an fpocket run on 1ubq, tabs and stray colons included, because the whitespace is the thing under test."""

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

INFO_TXT = "\n".join([
    "Pocket 1 :",
    "\tScore : \t0.433",
    "\tDruggability Score : \t0.235",
    "\tNumber of Alpha Spheres : \t50",
    "\tVolume : \t550.097",
    "\tMean alpha sphere radius :\t3.685",
    "\tHydrophobicity score:\t-2.111",
    "\tVolume score: \t 4.000",
    "",
    "Pocket 2 :",
    "\tDruggability Score : \t0.010",
    "\tVolume : \t120.500",
    "",
])


def load_pipe_fpocket():
    spec = importlib.util.spec_from_file_location("pipe_fpocket", ROOT / "pipe_scripts" / "pipe_fpocket.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["pipe_fpocket"] = module
    spec.loader.exec_module(module)
    return module


def write_info(tmp_path):
    path = tmp_path / "1ubq_info.txt"
    path.write_text(INFO_TXT, encoding="utf-8")
    return str(path)


def test_pocket_metrics_survive_the_leading_tab(tmp_path):
    pockets = load_pipe_fpocket().parse_info(write_info(tmp_path))
    assert len(pockets) == 2
    assert pockets[0]["pocket_idx"] == 1
    assert pockets[0]["druggability_score"] == 0.235
    assert pockets[0]["volume"] == 550.097
    assert pockets[0]["number_of_alpha_spheres"] == 50
    assert pockets[1]["pocket_idx"] == 2
    assert pockets[1]["druggability_score"] == 0.010


def test_keys_carry_no_colon_or_trailing_underscore(tmp_path):
    pockets = load_pipe_fpocket().parse_info(write_info(tmp_path))
    for key in pockets[0]:
        assert ":" not in key
        assert not key.endswith("_")
    assert pockets[0]["hydrophobicity_score"] == -2.111
    assert pockets[0]["volume_score"] == 4.0


def test_the_columns_the_tables_read_are_the_ones_produced(tmp_path):
    """The pockets table reads `druggability_score` and falls back to `volume`; the summary maxes over the same keys. A rename in the parser that missed these two would put empty cells in the tables again without failing anything."""
    pockets = load_pipe_fpocket().parse_info(write_info(tmp_path))
    top_drug = max((p.get("druggability_score", 0) for p in pockets), default=0)
    top_vol = max((p.get("real_volume_(approximation)", p.get("volume", 0)) for p in pockets), default=0)
    assert top_drug == 0.235
    assert top_vol == 550.097
