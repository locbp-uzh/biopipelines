"""Scripting resolves a bare script filename against the configured folder.

Scripting("foo.py") should find <scripts>/foo.py without an absolute path, while
an existing path given as-is always wins. A missing script raises listing every
place searched.
"""

import os

import pytest

from biopipelines.config_manager import ConfigManager
from biopipelines.scripting import Scripting


@pytest.fixture
def scripts_folder(monkeypatch, tmp_path):
    """Point ConfigManager.get_scripts_folder() at a temp dir."""
    folder = tmp_path / "scripts"
    folder.mkdir()
    monkeypatch.setattr(
        ConfigManager, "get_scripts_folder",
        lambda self: str(folder),
    )
    return folder


def _resolver():
    """A Scripting instance whose __init__ is bypassed (we only test the
    pure-resolution method, which needs no pipeline context)."""
    return Scripting.__new__(Scripting)


def test_bare_filename_found_in_scripts_folder(scripts_folder):
    script = scripts_folder / "step.py"
    script.write_text("def configuration(i): return {}\n", encoding="utf-8")
    resolved = _resolver()._resolve_script_path("step.py")
    assert resolved == os.path.abspath(str(script))


def test_existing_relative_path_wins_over_scripts_folder(scripts_folder, tmp_path, monkeypatch):
    # A file that exists as given (here via cwd) is used directly, even when a
    # same-named file also sits in the scripts folder.
    (scripts_folder / "dup.py").write_text("# in scripts folder\n", encoding="utf-8")
    workdir = tmp_path / "work"
    workdir.mkdir()
    (workdir / "dup.py").write_text("# in cwd\n", encoding="utf-8")
    monkeypatch.chdir(workdir)
    resolved = _resolver()._resolve_script_path("dup.py")
    assert resolved == os.path.abspath(str(workdir / "dup.py"))


def test_missing_script_raises_listing_candidates(scripts_folder):
    with pytest.raises(ValueError) as exc:
        _resolver()._resolve_script_path("nope.py")
    msg = str(exc.value)
    assert "nope.py" in msg
    # The scripts-folder candidate is listed (repr escapes backslashes on
    # Windows, so match the folder's leaf name rather than the full path).
    assert os.path.basename(str(scripts_folder)) in msg


def _standardized_output(streams):
    from biopipelines.base_config import StandardizedOutput, StreamContainer
    o = StandardizedOutput.__new__(StandardizedOutput)
    o.streams = StreamContainer(streams)
    return o


def _ds(name, n_ids):
    from biopipelines.datastream import DataStream
    ids = [f"{name}_{i}" for i in range(n_ids)]
    return DataStream(name=name, ids=ids, files=[f"{i}.pdb" for i in ids], format="pdb")


def test_pick_stream_exact_key_match():
    s = Scripting.__new__(Scripting)
    out = _standardized_output({"structures": _ds("structures", 2), "sequences": None})
    assert s._pick_stream("structures", out).name == "structures"


def test_pick_stream_single_stream_fallback():
    # Key does not match a stream name, but there's exactly one non-empty stream.
    s = Scripting.__new__(Scripting)
    out = _standardized_output({"sequences": _ds("sequences", 3), "structures": None})
    assert s._pick_stream("seqs", out).name == "sequences"


def test_pick_stream_ambiguous_raises():
    s = Scripting.__new__(Scripting)
    out = _standardized_output({"structures": _ds("structures", 2), "sequences": _ds("sequences", 2)})
    with pytest.raises(ValueError, match="multiple non-empty streams"):
        s._pick_stream("nope", out)


def test_absolute_path_used_directly(tmp_path, monkeypatch):
    # No scripts folder configured; an absolute path to an existing file works.
    monkeypatch.setattr(ConfigManager, "get_scripts_folder", lambda self: None)
    script = tmp_path / "abs.py"
    script.write_text("x = 1\n", encoding="utf-8")
    resolved = _resolver()._resolve_script_path(str(script))
    assert resolved == os.path.abspath(str(script))
