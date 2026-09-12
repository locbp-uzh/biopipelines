"""BindingData must send the SMILES it was given.

ChEMBL's similarity and substructure endpoints take the query pattern in the URL
*path*. It was interpolated raw, and requests keeps '#' and '/' in its safe set,
so:

  * 'C#N' produced .../similarity/C?limit=N#N/85.json -- everything from the '#'
    is a client-side fragment and never reaches the server. The request that goes
    out is a similarity search for plain carbon, with the threshold gone.
  * 'C/C=C/C' produced .../similarity/C/C=C/C/85.json, splitting the pattern into
    extra path segments.

Either way the wrong molecule is queried, and because the failure surfaces as an
HTTP error the compound is filed under `missing` as "no matching affinity record"
-- a service problem reported as a scientific negative.

These tests build the real PreparedRequest offline, so they assert on the bytes
that would leave the machine without touching the network.
"""

import importlib.util
import pathlib

import pytest

requests = pytest.importorskip("requests")

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"


def _load_module():
    """Import the pipe script directly; pipe_scripts is not a package."""
    path = pathlib.Path(__file__).resolve().parent.parent / "pipe_scripts" / "pipe_binding_data.py"
    spec = importlib.util.spec_from_file_location("pipe_binding_data", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def mod():
    return _load_module()


def _sent_url(endpoint, params=None):
    """The URL requests would actually transmit, fragment stripped as on the wire."""
    prepared = requests.Request("GET", f"{CHEMBL}/{endpoint}", params=params or {}).prepare()
    return prepared.url


# --- the encoder ------------------------------------------------------------

@pytest.mark.parametrize("smiles", [
    "C#N",                          # nitrile: the '#' case
    "C/C=C/C",                      # stereo slashes
    "CC(=O)Oc1ccccc1C(=O)O",        # aspirin, no special characters
    "C[Si](C)(C)c1ccccc1",          # silicon, brackets
    "O=C(N)c1ccc(cc1)S(=O)(=O)N",   # sulfonamide
    "c1ccc2c(c1)[nH]c1ccccc12",     # carbazole
    "[Na+].[Cl-]",                  # salt: '.' and charges
    "C%10CCCCC%10",                 # ring closure above 9: the '%' case
])
def test_encoded_smiles_survives_the_path(mod, smiles):
    """Whatever the SMILES, the encoded segment must round-trip intact."""
    from urllib.parse import unquote

    segment = mod.smiles_path_segment(smiles)
    url = _sent_url(f"similarity/{segment}/85.json", {"limit": 10})

    assert "#" not in url, f"{smiles!r}: a '#' survived and truncates the request"
    assert url.startswith(f"{CHEMBL}/similarity/"), url
    # exactly three path segments after /data: similarity / pattern / 85.json
    tail = url.split("/data/", 1)[1].split("?", 1)[0]
    parts = tail.split("/")
    assert len(parts) == 3, f"{smiles!r}: pattern split the path into {parts}"
    assert unquote(parts[1]) == smiles, f"{smiles!r}: round-trip gave {unquote(parts[1])!r}"
    assert parts[2] == "85.json", f"{smiles!r}: threshold lost, tail is {parts[2]!r}"


def test_raw_nitrile_demonstrates_the_bug_being_fixed(mod):
    """Guard the guard: unencoded, 'C#N' really does lose everything after '#'."""
    raw = _sent_url("similarity/C#N/85.json", {"limit": 10})
    assert "#" in raw and "85.json" not in raw.split("#")[0], (
        "the unencoded form no longer loses the threshold, so this test no longer "
        "demonstrates anything -- check whether requests changed its safe set")

    encoded = _sent_url(f"similarity/{mod.smiles_path_segment('C#N')}/85.json",
                        {"limit": 10})
    assert "85.json" in encoded and "#" not in encoded


def test_substructure_pattern_is_encoded_too(mod):
    segment = mod.smiles_path_segment("C#N")
    url = _sent_url(f"substructure/{segment}.json", {"limit": 10})
    assert "#" not in url
    assert "substructure/C%23N.json" in url


# --- politeness -------------------------------------------------------------

def test_a_user_agent_is_declared(mod):
    """Free academic endpoints should be able to attribute their traffic."""
    assert "biopipelines" in mod.USER_AGENT.lower()
    assert "python-urllib" not in mod.USER_AGENT.lower()
    assert "http" in mod.USER_AGENT, "no contact URL in the User-Agent"


def test_requests_are_rate_limited_and_retried(mod):
    assert mod.REQUEST_INTERVAL_S > 0, "no gap between calls to a shared service"
    assert mod.MAX_ATTEMPTS > 1, "a 429 would be recorded as a permanent failure"


def test_throttling_is_retried_not_reported_as_no_data(mod, monkeypatch):
    """429 must not become 'this compound has no known binder'."""
    calls = []

    class _Resp:
        def __init__(self, status):
            self.status_code = status
            self.headers = {"Retry-After": "0"}

        def raise_for_status(self):
            if self.status_code >= 400:
                raise requests.HTTPError(f"{self.status_code}")

        def json(self):
            return {"molecules": []}

    def fake_get(url, params=None, timeout=None, headers=None):
        calls.append(headers)
        return _Resp(429 if len(calls) == 1 else 200)

    monkeypatch.setattr(requests, "get", fake_get)
    monkeypatch.setattr(mod.time, "sleep", lambda s: None)

    response = mod._get(f"{CHEMBL}/molecule.json", {"limit": 1})

    assert len(calls) == 2, f"a 429 was not retried (calls={len(calls)})"
    assert response.status_code == 200
    assert all(h and "biopipelines" in h.get("User-Agent", "").lower() for h in calls)


def test_a_real_error_still_raises(mod, monkeypatch):
    """Retrying must not swallow a genuine 404."""
    class _Resp:
        status_code = 404
        headers = {}

        def raise_for_status(self):
            raise requests.HTTPError("404")

    monkeypatch.setattr(requests, "get",
                        lambda url, **kw: _Resp())
    monkeypatch.setattr(mod.time, "sleep", lambda s: None)

    with pytest.raises(requests.HTTPError):
        mod._get(f"{CHEMBL}/molecule.json", {"limit": 1})


# --- failure vs. absence ----------------------------------------------------

@pytest.mark.parametrize("raw,expected", [
    (">30000", (">", 30000.0)),
    ("<=5", ("<=", 5.0)),
    (" 40", ("=", 40.0)),
    ("40", ("=", 40.0)),
    ("~12.5", ("~", 12.5)),
    ("n/a", ("", None)),
    ("", ("", None)),
    (">abc", ("", None)),
])
def test_split_relation_does_not_assert_an_unread_equality(mod, raw, expected):
    """An unparseable value must not come back as '=' over a missing number."""
    assert mod.split_relation(raw) == expected


@pytest.mark.parametrize("payload", [
    {"getTargetByCompound": ["not", "an", "object"]},
    {"getTargetByCompound": "error: rate limited"},
    {"getTargetByCompound": {"message": "service unavailable"}},
])
def test_unrecognised_bindingdb_payload_raises(mod, monkeypatch, payload):
    """A schema or service problem must not read as 'no known binder'.

    Returning [] here would file the compound under `missing` with
    kind="filter" / "no matching affinity record", which is a scientific claim
    about the molecule rather than a report about the service.
    """
    import json as _json

    class _Resp:
        status_code = 200
        text = _json.dumps(payload)

        def raise_for_status(self):
            pass

    monkeypatch.setattr(requests, "get", lambda url, **kw: _Resp())
    monkeypatch.setattr(mod.time, "sleep", lambda s: None)

    config = {"match": "exact", "similarity": 1.0, "affinity_types": ["Ki"]}
    with pytest.raises(ValueError, match="unexpected BindingDB payload"):
        mod.bindingdb_rows("cmp_1", "CCO", config)


def test_genuinely_empty_bindingdb_answer_is_still_no_rows(mod, monkeypatch):
    """The real 'nothing matched' answer must stay distinguishable from a fault."""
    class _Resp:
        status_code = 200
        text = ""

        def raise_for_status(self):
            pass

    monkeypatch.setattr(requests, "get", lambda url, **kw: _Resp())
    monkeypatch.setattr(mod.time, "sleep", lambda s: None)

    config = {"match": "exact", "similarity": 1.0, "affinity_types": ["Ki"]}
    assert mod.bindingdb_rows("cmp_1", "CCO", config) == []
