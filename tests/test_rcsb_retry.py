"""A network blip is not evidence that a PDB id is wrong.

`PDB(...)` checks RCSB while the pipeline is being constructed, and every exception from that
request — timeout, connection reset, 5xx, 429 — used to become a `ValueError` that aborted the
build. So an unrelated outage presented as a bad PDB id, and the suite inherited the same
fragility: `conftest`'s reachability probe retries three times with exponential backoff, while the
code it guards retried not at all. A single blip minutes after a successful probe failed the run.

The distinction these tests protect is between the two answers RCSB can give. A 404 is a fact about
the id and must still stop the pipeline. Anything else is a fact about the network and must not.
"""

import pytest

from biopipelines import pdb as pdb_module


class Response:
    def __init__(self, status=200, payload=None):
        self.status_code = status
        self._payload = payload if payload is not None else {}

    def raise_for_status(self):
        if self.status_code >= 400:
            import requests
            raise requests.HTTPError(f"{self.status_code}", response=self)

    def json(self):
        return self._payload


def patch_requests(monkeypatch, responses):
    """Serve `responses` in order; an exception instance is raised instead of returned."""
    calls = []

    def fake_get(url, timeout=None):
        calls.append(url)
        item = responses[min(len(calls) - 1, len(responses) - 1)]
        if isinstance(item, Exception):
            raise item
        return item

    import requests
    monkeypatch.setattr(requests, "get", fake_get)
    monkeypatch.setattr(pdb_module.time if hasattr(pdb_module, "time") else __import__("time"),
                        "sleep", lambda _s: None, raising=False)
    return calls


class TestATransportFailureIsRetried:

    def test_it_recovers_when_a_later_attempt_succeeds(self, monkeypatch):
        import requests
        monkeypatch.setattr("time.sleep", lambda _s: None)
        calls = patch_requests(monkeypatch, [
            requests.ConnectionError("reset"),
            Response(200, {"rcsb_entry_info": {}}),
        ])
        assert pdb_module._rcsb_entry("https://data.rcsb.org/x") == {"rcsb_entry_info": {}}
        assert len(calls) == 2, "the first failure must be retried, not surfaced"

    def test_it_gives_up_as_unreachable_not_as_absent(self, monkeypatch):
        import requests
        monkeypatch.setattr("time.sleep", lambda _s: None)
        calls = patch_requests(monkeypatch, [requests.ReadTimeout("slow")])
        with pytest.raises(pdb_module._RcsbUnreachable):
            pdb_module._rcsb_entry("https://data.rcsb.org/x")
        assert len(calls) == 3, "three attempts before giving up, matching the conftest probe"

    def test_a_server_error_is_transport_not_absence(self, monkeypatch):
        monkeypatch.setattr("time.sleep", lambda _s: None)
        patch_requests(monkeypatch, [Response(503)])
        with pytest.raises(pdb_module._RcsbUnreachable):
            pdb_module._rcsb_entry("https://data.rcsb.org/x")

    @pytest.mark.parametrize("status", [400, 403])
    def test_a_refusal_is_not_retried(self, monkeypatch, status):
        monkeypatch.setattr("time.sleep", lambda _s: None)
        calls = patch_requests(monkeypatch, [Response(status)])
        with pytest.raises(pdb_module._RcsbUnreachable, match=f"HTTP {status}"):
            pdb_module._rcsb_entry("https://data.rcsb.org/x")
        assert len(calls) == 1

    def test_a_rate_limit_is_retried(self, monkeypatch):
        monkeypatch.setattr("time.sleep", lambda _s: None)
        calls = patch_requests(monkeypatch, [Response(429), Response(200, {"ok": 1})])
        assert pdb_module._rcsb_entry("https://data.rcsb.org/x") == {"ok": 1} and len(calls) == 2


class TestAMissingEntryStillStopsTheBuild:

    def test_a_404_raises_immediately_and_is_not_retried(self, monkeypatch):
        monkeypatch.setattr("time.sleep", lambda _s: None)
        calls = patch_requests(monkeypatch, [Response(404)])
        with pytest.raises(pdb_module._RcsbNotFound):
            pdb_module._rcsb_entry("https://data.rcsb.org/x")
        assert len(calls) == 1, (
            "retrying a 404 spends three round trips to be told the same true thing")

    def test_the_wrapper_turns_only_a_404_into_a_hard_error(self, monkeypatch, local_config,
                                                            isolated_cwd):
        """The whole point: one of these aborts pipeline construction, the other must not."""
        from biopipelines.pdb import PDB

        tool = PDB.__new__(PDB)
        tool.predicted_compound_ids = []

        monkeypatch.setattr(pdb_module, "_rcsb_entry",
                            lambda *_a, **_k: (_ for _ in ()).throw(pdb_module._RcsbNotFound("x")))
        with pytest.raises(ValueError, match="not found on RCSB"):
            tool._check_rcsb_exists("9XYZ")

        monkeypatch.setattr(pdb_module, "_rcsb_entry",
                            lambda *_a, **_k: (_ for _ in ()).throw(pdb_module._RcsbUnreachable("x")))
        assert tool._check_rcsb_exists("1UBQ") is False, (
            "an unreachable RCSB must not abort the build; it means the ligand check could not "
            "run, not that the entry is absent")


def test_the_entry_is_fetched_once_per_check(monkeypatch, local_config, isolated_cwd):
    """`_check_rcsb_exists` used to GET the entry, then call a helper that GET the same URL again."""
    from biopipelines.pdb import PDB

    tool = PDB.__new__(PDB)
    tool.predicted_compound_ids = []
    seen = []
    monkeypatch.setattr(pdb_module, "_rcsb_entry",
                        lambda url, *a, **k: (seen.append(url), {"rcsb_entry_info": {}})[1])
    tool._check_rcsb_exists("1UBQ")
    assert len(seen) == 1, f"the entry was requested {len(seen)} times for one check"
