"""The committed templates have to run for whoever copies them.

`my_pipelines/_template.py` calls `Scripting("_template.py")`, which resolves against
`my_scripts/`. For a while it named `template.py`, which existed in nobody's checkout: the
folder is gitignored except for a short whitelist, so the file a fresh clone needs was simply
absent and the template failed at its only real step. Nothing caught it, because a template is
documentation that happens to be executable and no test ran it.

These tests close that: every script a committed template names must itself be committed, and
the script template must produce the shapes it claims.
"""

import importlib.util
import pathlib
import re
import subprocess

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
TEMPLATES = [ROOT / "my_pipelines" / "_template.py", ROOT / "my_pipelines" / "_template.ipynb"]


def tracked(rel):
    """Is this path committed? `my_scripts/` and `my_pipelines/` are ignored but for a whitelist.

    A missing `git` is raised as itself rather than as a FileNotFoundError from inside a helper.
    `python:*-slim` has no git binary, so these tests died with a traceback that named neither
    the test nor the cause, and main stayed red for a day before anyone read it closely enough
    to see that the failure was about the runner rather than about the templates.
    """
    try:
        done = subprocess.run(["git", "ls-files", "--error-unmatch", rel],
                              cwd=ROOT, capture_output=True, text=True)
    except FileNotFoundError:
        raise AssertionError(
            "`git` is not on PATH, so this test cannot tell whether the templates are "
            "committed. Install git in the CI image — skipping would report success for a "
            "check that never ran.") from None
    return done.returncode == 0


@pytest.mark.parametrize("template", TEMPLATES, ids=lambda p: p.name)
def test_a_template_only_names_scripts_that_ship_with_it(template):
    named = set(re.findall(r'Scripting\(\\?"([^"\\]+)\\?"', template.read_text(encoding="utf-8")))
    assert named, f"{template.name} no longer calls Scripting; update this test"
    for script in named:
        path = ROOT / "my_scripts" / script
        assert path.exists(), f"{template.name} calls Scripting({script!r}), absent from my_scripts/"
        assert tracked(f"my_scripts/{script}"), (
            f"my_scripts/{script} exists here but is not committed, so a fresh clone lacks it "
            f"— add it to the whitelist in .gitignore")


@pytest.mark.parametrize("rel", ["my_pipelines/_template.py", "my_pipelines/_template.ipynb",
                                 "my_scripts/_template.py"])
def test_the_templates_themselves_are_committed(rel):
    assert tracked(rel), f"{rel} is not committed; check the whitelist in .gitignore"


class TestTheScriptTemplate:
    """It is the file a user copies, so its two functions must work as written."""

    @staticmethod
    def load():
        spec = importlib.util.spec_from_file_location(
            "script_template", ROOT / "my_scripts" / "_template.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    @staticmethod
    def structures(tmp_path):
        pdb = tmp_path / "s1.pdb"
        pdb.write_text(
            "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
            "ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C\n"
            "ATOM      3  CA  GLY A   2       2.000   0.000   0.000  1.00  0.00           C\n"
            "HETATM    4  O   HOH A 101      9.000   0.000   0.000  1.00  0.00           O\n")

        class Stream:
            ids = ["s1"]
            format = "pdb"

            def iterate(self):
                return [("s1", str(pdb))]

        return Stream(), pdb

    def test_configuration_declares_a_stream_and_a_table(self, tmp_path):
        stream, _pdb = self.structures(tmp_path)
        declared = self.load().configuration({"structures": stream})
        assert set(declared) == {"structures", "metrics"}
        assert declared["structures"].format == "pdb"
        assert declared["structures"].ids == ["s1"], "ids must come from the input, not be invented"
        assert declared["metrics"].columns == ["id", "n_residues"]

    def test_execution_fills_both_and_passes_the_file_through(self, tmp_path):
        stream, pdb = self.structures(tmp_path)
        out_dir = tmp_path / "out"

        class Out:
            def __init__(self):
                self.rows = []

            def file(self, _id, name):
                out_dir.mkdir(exist_ok=True)
                return str(out_dir / name)

            def row(self, r):
                self.rows.append(r)

        class Outs(dict):
            def drop(self, _id, cause=""):
                raise AssertionError("the template drops nothing")

        outs = Outs(structures=Out(), metrics=Out())
        self.load().execution({"structures": stream}, outs)

        assert outs["metrics"].rows == [{"id": "s1", "n_residues": 2}], "CA count, waters excluded"
        assert (out_dir / "s1.pdb").read_text() == pdb.read_text(), "pass-through must be verbatim"
