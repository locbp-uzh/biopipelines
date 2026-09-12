"""Remote MSA-server harvest: the ColabFold result format, as the server sends it.

The server returns one member per DATABASE — uniref.a3m always, and
bfd.mgnify30.metaeuk30.smag30.a3m under mode=env — each concatenating every
query's alignment separated by NUL bytes, in submission order. It does NOT
return one member per query. An earlier version of this file asserted the
opposite and passed while the remote path could never work against the real
endpoint.
"""

import importlib.util
import io
import os
import sys
import tarfile
import urllib.request

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HELPER = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_mmseqs2_sequences.py")
NUL = bytes([0])
UNIREF = "uniref.a3m"
ENVDB = "bfd.mgnify30.metaeuk30.smag30.a3m"


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
        url = request.full_url if hasattr(request, "full_url") else str(request)
        if "/ticket/msa" in url:
            return _Response(b'{"id": "TICKET"}')
        if "/ticket/" in url:
            return _Response(b'{"status": "COMPLETE"}')
        return _Response(blob)

    monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)


def _tarball(members):
    """members: {filename: bytes}."""
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        for name, payload in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(payload)
            tar.addfile(info, io.BytesIO(payload))
    return buf.getvalue()


def _query(alias, seq):
    """One query's alignment: its own >101+i header, then the query sequence."""
    return f">{101 + alias}\n{seq}\n>hit_{alias}\n{seq}".encode()


def _concat(prefix, n):
    """One database member: n query alignments separated by NUL."""
    return NUL.join(_query(i, "AAAA") for i in range(n))


def _concat_permuted(n, order):
    """The same alignments, delivered in the server's own internal order."""
    return NUL.join(_query(i, "AAAA") for i in order)


def test_queries_are_split_out_of_the_concatenated_databases(tmp_path, monkeypatch):
    """Each query gets its own slice, and both databases are joined into it."""
    module = _load_helper()
    n = 12
    _serve(monkeypatch, _tarball({UNIREF: _concat("UNIREF", n),
                                  ENVDB: _concat("ENV", n)}))
    seqs = [(f"seq_{i}", "AAAA") for i in range(n)]

    assert module.submit_batch_http("https://msa.example", seqs, str(tmp_path)) is True

    for alias in range(n):
        got = (tmp_path / f"{alias}.a3m").read_text()
        assert got.startswith(f">{101 + alias}\n"), f"alias {alias} got: {got[:40]!r}"
        assert got.count(f">{101 + alias}\n") == 2, "both databases must be joined in"


def test_chunks_are_paired_by_header_not_position(tmp_path, monkeypatch):
    """The server returns alignments in its OWN order, not submission order.

    Pairing by position attached 98 of 100 alignments to the wrong design in
    fullmsa_001 -- every fold then fell back to a dummy MSA. The query header
    is the only trustworthy link back to the submitted sequence.
    """
    module = _load_helper()
    n = 6
    order = [3, 0, 5, 1, 4, 2]
    _serve(monkeypatch, _tarball({UNIREF: _concat_permuted(n, order)}))
    seqs = [(f"seq_{i}", "AAAA") for i in range(n)]

    assert module.submit_batch_http("https://msa.example", seqs, str(tmp_path)) is True

    for alias in range(n):
        got = (tmp_path / f"{alias}.a3m").read_text()
        assert got.startswith(f">{101 + alias}\n"), (
            f"alias {alias} received a neighbour's alignment: {got[:40]!r}")


def test_alignment_disagreeing_with_the_submission_is_refused(tmp_path, monkeypatch):
    """A header that survives but a sequence that does not match is not usable."""
    module = _load_helper()
    payload = NUL.join([_query(0, "AAAA"), _query(1, "WWWW")])
    _serve(monkeypatch, _tarball({UNIREF: payload}))
    seqs = [("seq_a", "AAAA"), ("seq_b", "CCCC")]

    assert module.submit_batch_http("https://msa.example", seqs, str(tmp_path)) is False


def test_uniref_alone_still_succeeds(tmp_path, monkeypatch):
    """A server without the environmental database is usable, just shallower."""
    module = _load_helper()
    _serve(monkeypatch, _tarball({UNIREF: _concat("UNIREF", 2)}))
    seqs = [("seq_a", "AAAA"), ("seq_b", "AAAA")]

    assert module.submit_batch_http("https://msa.example", seqs, str(tmp_path)) is True
    assert (tmp_path / "0.a3m").read_text().startswith(">101\n")


def test_missing_uniref_fails_loudly(tmp_path, monkeypatch):
    """No recognised database member means no alignment — do not invent one."""
    module = _load_helper()
    _serve(monkeypatch, _tarball({"something_else.a3m": b"X"}))
    seqs = [("seq_a", "AAAA")]

    assert module.submit_batch_http("https://msa.example", seqs, str(tmp_path)) is False


def test_wrong_query_count_fails_loudly(tmp_path, monkeypatch):
    """If the split does not match the submission, the pairing is unknowable."""
    module = _load_helper()
    _serve(monkeypatch, _tarball({UNIREF: _concat("UNIREF", 3)}))
    seqs = [("seq_a", "AAAA"), ("seq_b", "CCCC")]

    assert module.submit_batch_http("https://msa.example", seqs, str(tmp_path)) is False


def test_server_url_accepts_csv_output_format():
    """csv is legal with server_url: a3m comes off the wire and is converted."""
    from biopipelines.mmseqs2 import MMseqs2

    for fmt in ("csv", "a3m"):
        MMseqs2(sequences=["MKV"], server_url="https://msa.example",
                output_format=fmt).validate_params()
