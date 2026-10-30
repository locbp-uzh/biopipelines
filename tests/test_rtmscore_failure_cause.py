# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""RTMScore: the `cause` of a failure row names the error, not the head of a traceback.

A failed `rtmscore.py` run raised `RuntimeError("rtmscore.py failed (exit 1): <last 1500 characters of stderr>")` and the row stored `str(e)[:200]`, the first 200 characters of that tail. The tail of a Python traceback starts in the middle of a chained exception, so the stored cause was a fragment of stack and the actual reason ("The graph of pocket cannot be generated") sat at the end, cut off.
"""

import importlib.util
import os
import sys

import pytest

pd = pytest.importorskip("pandas")

from biopipelines.datastream import DataStream, create_map_table

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

TRACEBACK_TAIL = (
    "etCoords'\n\nDuring handling of the above exception, another exception occurred:\n\n"
    "Traceback (most recent call last):\n"
    "  File \"RTMScore/data/data.py\", line 120, in __init__\n    extract_pocket(prot, reflig, cutoff,\n"
    "TypeError: atoms must be an object with coordinate sets\n\n"
    "During handling of the above exception, another exception occurred:\n\n"
    "Traceback (most recent call last):\n"
    "  File \"example/rtmscore.py\", line 228, in <module>\n    main()\n"
    "ValueError: The graph of pocket cannot be generated"
)


def _stream(tmp_path, label, ids, ext):
    name = "structures"
    table = tmp_path / f"{label}_map.csv"
    create_map_table(str(table), ids=list(ids))
    for one in ids:
        (tmp_path / f"{one}.{ext}").write_text("ATOM\n")
    return DataStream(name=name, ids=list(ids), map_table=str(table), format=ext,
                      files=[str(tmp_path / f"<id>.{ext}")])


@pytest.fixture
def rtmscore(monkeypatch):
    spec = importlib.util.spec_from_file_location("pipe_rtmscore_cause_test", os.path.join(REPO_ROOT, "pipe_scripts", "pipe_rtmscore.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, "write_ligand_sdf", lambda src, dst, smiles: open(dst, "w").write("sdf\n"))
    return module


def _run(rtmscore, monkeypatch, tmp_path, error):
    proteins = _stream(tmp_path, "proteins", ["prot+one"], "pdb")
    ligands = _stream(tmp_path, "ligands", ["lig_a"], "sdf")
    proteins.save_json(str(tmp_path / "proteins.json"))
    ligands.save_json(str(tmp_path / "ligands.json"))

    def fail(*args, **kwargs):
        raise error

    monkeypatch.setattr(rtmscore, "run_rtmscore", fail)
    monkeypatch.setattr(sys, "argv", [
        "pipe_rtmscore.py",
        "--structures-json", str(tmp_path / "proteins.json"),
        "--compounds-json", str(tmp_path / "ligands.json"),
        "--rtmscore-script", "rtmscore.py",
        "--model-path", "model.pth",
        "--scratch-dir", str(tmp_path / "scratch"),
        "--scores-csv", str(tmp_path / "scores.csv"),
        "--missing-csv", str(tmp_path / "001_RTMScore" / "missing.csv"),
    ])
    with pytest.raises(SystemExit) as exit_info:
        rtmscore.main()
    return exit_info.value.code, pd.read_csv(tmp_path / "001_RTMScore" / "missing.csv")


def test_a_failed_run_is_a_failure_row_whose_cause_is_the_final_error_line(rtmscore, monkeypatch, tmp_path, record_case):
    error = RuntimeError(f"rtmscore.py failed (exit 1): {TRACEBACK_TAIL}")

    code, missing = _run(rtmscore, monkeypatch, tmp_path, error)

    record_case(input="rtmscore.py exits 1 with a chained traceback ending in a ValueError",
                expected="ValueError: The graph of pocket cannot be generated", actual=missing.loc[0, "cause"])
    assert code == 1
    assert list(missing["id"]) == ["prot+one+lig_a"]
    assert list(missing["kind"]) == ["failure"]
    assert missing.loc[0, "cause"] == "ValueError: The graph of pocket cannot be generated"


def test_a_one_line_error_is_kept_whole(rtmscore, monkeypatch, tmp_path):
    _, missing = _run(rtmscore, monkeypatch, tmp_path, RuntimeError("rtmscore.py produced no output CSV at /x/scores.csv"))

    assert missing.loc[0, "cause"] == "rtmscore.py produced no output CSV at /x/scores.csv"


def test_a_long_final_line_is_cut_to_200_characters(rtmscore, monkeypatch, tmp_path):
    _, missing = _run(rtmscore, monkeypatch, tmp_path, RuntimeError("first\n" + "x" * 500))

    assert missing.loc[0, "cause"] == "x" * 200
