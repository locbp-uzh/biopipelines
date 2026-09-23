"""DSSP recovers from mkdssp's duplicate-`refine` defect instead of failing the run.

mkdssp 4.6.1 converts PDB input to mmCIF internally, and for some entries that conversion emits
two `refine` rows with the same key, which its validator rejects. 4UFC is one: it failed twice
on S3IT before the runner learned to retry without `REMARK 3`.

The fix is invisible until it is needed, which is exactly why it needs a test — nothing else
would notice it rotting.
"""

import importlib.util
import io
import pathlib

import pytest

pytest.importorskip("pandas")

ROOT = pathlib.Path(__file__).resolve().parent.parent
RUNNER = ROOT / "pipe_scripts" / "pipe_dssp.py"

DUPLICATE_KEY = ("Error while setting validator in datablock 4UFC\n"
                 " >> Duplicate Key violation, cat: refine values: aniso_B[1][1]: \"-0.55\"")


@pytest.fixture(scope="module")
def dssp():
    spec = importlib.util.spec_from_file_location("pipe_dssp", RUNNER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_pdb(path, remark3=True):
    lines = ["HEADER    HYDROLASE                               16-MAR-15   4UFC\n"]
    if remark3:
        lines += [f"REMARK   3   LINE {i}\n" for i in range(5)]
    lines += ["REMARK 465 MISSING RESIDUES\n",
              "ATOM      1  N   MET A   1      10.000  10.000  10.000  1.00 20.00           N\n",
              "END\n"]
    io.open(path, "w", encoding="utf-8").write("".join(lines))
    return str(path)


class TestStripping:
    def test_only_remark_3_is_removed(self, dssp, tmp_path):
        cleaned = dssp.strip_refinement_remarks(write_pdb(tmp_path / "in.pdb"))
        body = io.open(cleaned, encoding="utf-8").read()
        assert "REMARK   3" not in body
        assert "REMARK 465" in body, "other REMARK categories must survive"
        assert "ATOM      1" in body and "HEADER" in body

    def test_a_file_without_remark_3_is_passed_through_untouched(self, dssp, tmp_path):
        original = write_pdb(tmp_path / "in.pdb", remark3=False)
        assert dssp.strip_refinement_remarks(original) == original

    def test_the_copy_is_not_written_beside_the_results(self, dssp, tmp_path):
        """A sanitization artifact in a declared output stream is what run 005 produced."""
        original = write_pdb(tmp_path / "in.pdb")
        cleaned = dssp.strip_refinement_remarks(original)
        assert pathlib.Path(cleaned).parent != pathlib.Path(original).parent

    def test_an_unreadable_file_is_returned_unchanged(self, dssp, tmp_path):
        missing = str(tmp_path / "nope.pdb")
        assert dssp.strip_refinement_remarks(missing) == missing


class TestRetry:
    """`run_dssp` drives a binary, so the binary is faked: it fails on input carrying REMARK 3."""

    def fake_mkdssp(self, dssp, monkeypatch, calls):
        def run(argv, capture_output=None, text=None):
            path, out = argv[-2], argv[-1]
            calls.append(path)
            body = io.open(path, encoding="utf-8").read()
            if "REMARK   3" in body:
                return type("R", (), {"returncode": 1, "stdout": "", "stderr": DUPLICATE_KEY})()
            io.open(out, "w", encoding="utf-8").write("  #  RESIDUE\n    1 A M  H\n")
            return type("R", (), {"returncode": 0, "stdout": "", "stderr": ""})()

        monkeypatch.setattr(dssp.subprocess, "run", run)

    def test_a_duplicate_key_failure_is_retried_without_remark_3(self, dssp, tmp_path,
                                                                monkeypatch, capsys):
        calls = []
        self.fake_mkdssp(dssp, monkeypatch, calls)
        out = tmp_path / "out.dssp"
        dssp.run_dssp("mkdssp", write_pdb(tmp_path / "in.pdb"), str(out))

        assert out.exists() and out.stat().st_size > 0
        assert len(calls) > 1, "the first attempt should have failed and been retried"
        assert "REMARK 3 stripped" in capsys.readouterr().out, (
            "a silent input rewrite should still be visible in the log")

    def test_a_file_that_needs_no_retry_is_run_once(self, dssp, tmp_path, monkeypatch):
        calls = []
        self.fake_mkdssp(dssp, monkeypatch, calls)
        dssp.run_dssp("mkdssp", write_pdb(tmp_path / "in.pdb", remark3=False),
                      str(tmp_path / "out.dssp"))
        assert len(calls) == 1, "the common path must not pay for the workaround"

    def test_a_failure_the_retry_cannot_fix_still_raises(self, dssp, tmp_path, monkeypatch):
        def always_fails(argv, capture_output=None, text=None):
            return type("R", (), {"returncode": 1, "stdout": "", "stderr": "broken"})()

        monkeypatch.setattr(dssp.subprocess, "run", always_fails)
        with pytest.raises(RuntimeError, match="broken"):
            dssp.run_dssp("mkdssp", write_pdb(tmp_path / "in.pdb"), str(tmp_path / "out.dssp"))

    def test_the_retry_reports_its_own_error_and_removes_its_copy(self, dssp, tmp_path, monkeypatch):
        """The retry's failure was dropped for the first attempt's, and each retry leaked a temp dir."""
        seen = []

        def fails_differently(argv, capture_output=None, text=None):
            seen.append(argv[-2])
            body = io.open(argv[-2], encoding="utf-8").read()
            reason = "Duplicate Key violation, cat: refine" if "REMARK   3" in body else "second reason"
            return type("R", (), {"returncode": 1, "stdout": "", "stderr": reason})()

        monkeypatch.setattr(dssp.subprocess, "run", fails_differently)
        with pytest.raises(RuntimeError, match="second reason"):
            dssp.run_dssp("mkdssp", write_pdb(tmp_path / "in.pdb"), str(tmp_path / "out.dssp"))
        copies = [p for p in seen if "dssp_norefine_" in p]
        assert copies and not any(pathlib.Path(p).parent.exists() for p in copies)
