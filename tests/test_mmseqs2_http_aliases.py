"""Remote MSA-server harvest: members must map to queries by alias, not position."""

import importlib.util
import io
import json
import os
import sys
import tarfile
import urllib.request

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HELPER = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_mmseqs2_sequences.py")


def _load_helper():
    spec = importlib.util.spec_from_file_location("pipe_mmseqs2_sequences", HELPER)
    module = importlib.util.module_from_spec(spec)
    sys.modules["pipe_mmseqs2_sequences"] = module
    spec.loader.exec_module(module)
    return module


class _Response:
    def __init__(self, payload):
        self._payload = payload

    def read(self):
        return self._payload

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _serve(monkeypatch, blob):
    def fake_urlopen(request, timeout=None):
        url = request.full_url
        if url.endswith("/ticket/msa"):
            return _Response(json.dumps({"id": "t1"}).encode())
        if "/ticket/" in url:
            return _Response(json.dumps({"status": "COMPLETE"}).encode())
        return _Response(blob)

    monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)


def _tarball(names):
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        for name in names:
            payload = f"MSA_{os.path.splitext(name)[0]}".encode()
            info = tarfile.TarInfo(name)
            info.size = len(payload)
            tar.addfile(info, io.BytesIO(payload))
    return buf.getvalue()


def test_double_digit_batch_maps_each_member_to_its_own_query(tmp_path, monkeypatch):
    """Members sort lexically ("10" before "2"), so a positional rename would hand
    most queries another sequence's alignment once the batch passes ten."""
    module = _load_helper()
    n = 12
    _serve(monkeypatch, _tarball([f"{a}.a3m" for a in range(n)]))
    seqs = [(f"seq_{i}", "AAAA") for i in range(n)]

    assert module.submit_batch_http("https://msa.example", seqs, str(tmp_path)) is True

    for alias in range(n):
        got = (tmp_path / f"{alias}.a3m").read_text()
        assert got == f"MSA_{alias}", f"alias {alias} received {got}"


def test_unmappable_member_name_fails_loudly(tmp_path, monkeypatch):
    """A server naming results by query header rather than alias must error, not
    silently misalign."""
    module = _load_helper()
    _serve(monkeypatch, _tarball(["seq_a.a3m", "seq_b.a3m"]))
    seqs = [("seq_a", "AAAA"), ("seq_b", "CCCC")]

    assert module.submit_batch_http("https://msa.example", seqs, str(tmp_path)) is False


def test_out_of_range_alias_fails_loudly(tmp_path, monkeypatch):
    module = _load_helper()
    _serve(monkeypatch, _tarball(["0.a3m", "7.a3m"]))
    seqs = [("seq_a", "AAAA"), ("seq_b", "CCCC")]

    assert module.submit_batch_http("https://msa.example", seqs, str(tmp_path)) is False


def test_server_url_requires_a3m_output_format():
    """The ColabFold protocol returns a3m; a csv request would have the harvester
    look for files the remote path never writes."""
    from biopipelines.mmseqs2 import MMseqs2

    with pytest.raises(ValueError, match="requires output_format='a3m'"):
        MMseqs2(sequences=["MKV"], server_url="https://msa.example",
                output_format="csv").validate_params()
