#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Runtime helper for BindingData.

Reads SMILES from the compounds DataStream's map_table, queries ChEMBL and/or
BindingDB for measured affinity records against protein targets, and writes the
per-record affinities table plus a per-target aggregate.
"""

import os
import sys
import json
import argparse
from datetime import datetime, timezone
import time
from urllib.parse import quote

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from biopipelines.biopipelines_io import load_datastream, iterate_values
from biopipelines.binding_data import CHEMBL_API_URL, BINDINGDB_API_URL

TOOL_NAME = "BindingData"

# Columns of the affinities table, in declared order.
AFFINITY_COLUMNS = [
    "id", "smiles", "source", "source_molecule_id", "target_id", "target_name",
    "organism", "affinity_type", "relation", "affinity_nm", "pchembl_value",
    "assay_id", "assay_description", "assay_type", "document_id", "match_similarity",
]

# BindingDB reports the relation glued onto the value ('>30000', ' 40').
RELATIONS = (">=", "<=", ">", "<", "~", "=")


# Both endpoints are free academic services. Identify ourselves so their operators
# can attribute and debug the traffic, and leave a gap between calls.
USER_AGENT = ("biopipelines/1.0 (LOCBP, University of Zurich; "
              "https://github.com/locbp-uzh/biopipelines-locbp)")
REQUEST_INTERVAL_S = 0.34      # ~3 requests/second
MAX_ATTEMPTS = 4
_last_request = [0.0]


def smiles_path_segment(smiles):
    """Percent-encode a SMILES for use as a URL *path* segment.

    ChEMBL's similarity/substructure endpoints take the pattern in the path, and a
    raw SMILES breaks it: requests keeps '#' and '/' in its safe set, so 'C#N'
    truncates at the '#' (everything after it becomes a client-side fragment and is
    never sent, threshold included) and 'C/C=C/C' splits into extra path segments.
    Both cases silently query for the wrong molecule.
    """
    return quote(str(smiles), safe="")


def _get(url, params):
    """One GET, rate-limited and retried, with a generous timeout."""
    import requests

    for attempt in range(MAX_ATTEMPTS):
        wait = REQUEST_INTERVAL_S - (time.monotonic() - _last_request[0])
        if wait > 0:
            time.sleep(wait)
        _last_request[0] = time.monotonic()

        response = requests.get(url, params=params, timeout=120,
                                headers={"User-Agent": USER_AGENT})
        # 429/503 mean "come back later", not "no such compound". Retrying keeps a
        # throttled batch from being recorded as a batch of scientific negatives.
        if response.status_code in (429, 503) and attempt < MAX_ATTEMPTS - 1:
            delay = float(response.headers.get("Retry-After") or 2 ** attempt)
            print(f"  {response.status_code} from {url.split('/')[2]}, retrying in "
                  f"{delay:.0f}s", file=sys.stderr)
            time.sleep(delay)
            continue
        response.raise_for_status()
        return response
    response.raise_for_status()
    return response


def _to_float(value):
    """Parse a numeric field, yielding None for absent or malformed values."""
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def chembl_release():
    """The ChEMBL release that answered, or "" when the status endpoint declines.

    These tables are a snapshot of a database that grows with every release, so
    two runs months apart give different numbers with nothing in the output to
    tell them apart. Recording the release is what makes a result citable.
    """
    try:
        response = _get(f'{CHEMBL_API_URL}/status.json', {})
        return str(response.json().get('chembl_db_version') or '')
    except Exception as e:
        # Provenance is worth recording but not worth failing a whole run over.
        print(f'  could not read the ChEMBL release: {e}', file=sys.stderr)
        return ''


def resolve_chembl_molecules(smiles, match, similarity, max_records):
    """Return [(molecule_chembl_id, match_similarity)] for one query SMILES."""
    if match == "exact":
        endpoint = "molecule.json"
        params = {"molecule_structures__canonical_smiles__flexmatch": smiles,
                  "limit": max_records}
    else:
        # The similarity/substructure endpoints take the pattern in the path, so it
        # has to be percent-encoded -- see smiles_path_segment.
        percent = int(round(similarity * 100))
        pattern = smiles_path_segment(smiles)
        stem = (f"similarity/{pattern}/{percent}" if match == "similarity"
                else f"substructure/{pattern}")
        endpoint = f"{stem}.json"
        params = {"limit": max_records}

    response = _get(f"{CHEMBL_API_URL}/{endpoint}", params)
    molecules = response.json().get("molecules", [])
    resolved = []
    for molecule in molecules:
        chembl_id = molecule.get("molecule_chembl_id")
        if chembl_id:
            resolved.append((chembl_id, _to_float(molecule.get("similarity"))))
    return resolved


def fetch_chembl_activities(chembl_id, affinity_types, max_records):
    """Return the standardized activity records for one ChEMBL molecule."""
    params = {
        "molecule_chembl_id": chembl_id,
        "standard_type__in": ",".join(affinity_types),
        "limit": max_records,
    }
    response = _get(f"{CHEMBL_API_URL}/activity.json", params)
    return response.json().get("activities", [])


def chembl_rows(compound_id, smiles, config):
    """Collect the affinities table rows a ChEMBL lookup yields for one compound."""
    molecules = resolve_chembl_molecules(
        smiles, config["match"], config["similarity"], config["max_records"]
    )

    rows = []
    for chembl_id, match_similarity in molecules:
        for activity in fetch_chembl_activities(
            chembl_id, config["affinity_types"], config["max_records"]
        ):
            # ChEMBL harmonises units per standard_type; anything not already in
            # nM is a different quantity (e.g. a percentage) rather than a
            # convertible one, so it is not an affinity record we can report.
            if (activity.get("standard_units") or "") != "nM":
                continue
            rows.append({
                "id": compound_id,
                "smiles": smiles,
                "source": "chembl",
                "source_molecule_id": chembl_id,
                "target_id": activity.get("target_chembl_id"),
                "target_name": activity.get("target_pref_name"),
                "organism": activity.get("target_organism"),
                "affinity_type": activity.get("standard_type"),
                "relation": activity.get("standard_relation"),
                "affinity_nm": _to_float(activity.get("standard_value")),
                "pchembl_value": _to_float(activity.get("pchembl_value")),
                "assay_id": activity.get("assay_chembl_id"),
                "assay_description": activity.get("assay_description"),
                "assay_type": activity.get("assay_type"),
                "document_id": activity.get("document_chembl_id"),
                "match_similarity": match_similarity,
            })
    return rows


def split_relation(affinity):
    """Split BindingDB's glued affinity string into (relation, value_nM).

    An unparseable value yields ("", None) rather than ("=", None): reporting '='
    over a number that was never read asserts an equality nobody measured.
    """
    text = str(affinity).strip()
    for relation in RELATIONS:
        if text.startswith(relation):
            value = _to_float(text[len(relation):].strip())
            return (relation, value) if value is not None else ("", None)
    value = _to_float(text)
    return ("=", value) if value is not None else ("", None)


def bindingdb_rows(compound_id, smiles, config):
    """Collect the affinities table rows a BindingDB lookup yields for one compound."""
    cutoff = 1.0 if config["match"] == "exact" else config["similarity"]
    response = _get(f"{BINDINGDB_API_URL}/getTargetByCompound",
                    {"smiles": smiles, "cutoff": cutoff, "response": "application/json"})

    # The service answers with an empty body when nothing matches, and names its
    # payload after a different endpoint; read the single top-level value rather
    # than that misleading key.
    text = response.text.strip()
    if not text:
        return []
    payload = json.loads(text)
    if not isinstance(payload, dict) or not payload:
        return []
    body = next(iter(payload.values()))
    if not isinstance(body, dict):
        # A shape we do not recognise is a service or schema problem, not an answer.
        # Returning [] here would file the compound under `missing` as "no matching
        # affinity record" -- reporting a failure as a scientific negative.
        raise ValueError(
            f"unexpected BindingDB payload for {compound_id}: top-level value is "
            f"{type(body).__name__}, expected an object with 'bdb.affinities'")
    if "bdb.affinities" not in body:
        raise ValueError(
            f"unexpected BindingDB payload for {compound_id}: no 'bdb.affinities' key "
            f"(got {sorted(body)[:6]})")
    records = body.get("bdb.affinities") or []

    rows = []
    for record in records:
        affinity_type = record.get("bdb.affinity_type")
        if affinity_type not in config["affinity_types"]:
            continue
        relation, value = split_relation(record.get("bdb.affinity"))
        rows.append({
            "id": compound_id,
            "smiles": smiles,
            "source": "bindingdb",
            "source_molecule_id": record.get("bdb.monomerid"),
            "target_id": None,
            "target_name": record.get("bdb.target"),
            "organism": record.get("bdb.species"),
            "affinity_type": affinity_type,
            "relation": relation,
            "affinity_nm": value,
            "pchembl_value": None,
            "assay_id": None,
            "assay_description": None,
            "assay_type": None,
            "document_id": None,
            "match_similarity": None,
        })
    return rows


def apply_filters(rows, config):
    """Drop records outside the requested potency / organism window."""
    max_affinity = config["max_affinity"]
    organism = (config["organism"] or "").strip().lower()

    kept = []
    for row in rows:
        if max_affinity is not None:
            value = row["affinity_nm"]
            # A record with no parsable value cannot be shown to satisfy a
            # potency ceiling, so an explicit ceiling excludes it.
            if value is None or value > max_affinity:
                continue
        if organism and organism not in (row["organism"] or "").lower():
            continue
        kept.append(row)
    return kept


TARGET_COLUMNS = ["target_id", "source", "target_name", "organism",
                  "n_compounds", "n_records", "best_affinity_nm"]


def build_targets(affinities):
    """Aggregate the record-level rows into one row per target and source."""
    resolved = affinities[affinities["target_id"].notna()
                          | affinities["target_name"].notna()]
    if resolved.empty:
        return pd.DataFrame(columns=TARGET_COLUMNS)

    # Group on the target's identity where the source supplies one: distinct
    # ChEMBL targets share a pref_name (human and sheep COX-1 both read
    # "Prostaglandin G/H synthase 1"), so keying on the name merges them.
    # BindingDB returns no target id, leaving the name as its only key.
    resolved = resolved.assign(
        _key=resolved["target_id"].where(resolved["target_id"].notna(),
                                         resolved["target_name"])
    )

    targets = resolved.groupby(["_key", "source"], dropna=False).agg(
        target_id=("target_id", "first"),
        target_name=("target_name", "first"),
        organism=("organism", "first"),
        n_compounds=("id", "nunique"),
        n_records=("id", "size"),
        best_affinity_nm=("affinity_nm", "min"),
    ).reset_index()

    return targets[TARGET_COLUMNS]


def main():
    parser = argparse.ArgumentParser(description="Fetch measured binding affinities")
    parser.add_argument("--config", required=True, help="Path to binding_data_config.json")
    args = parser.parse_args()

    with open(args.config) as f:
        config = json.load(f)

    compounds_ds = load_datastream(config["compounds_json"])

    fetchers = {"chembl": chembl_rows, "bindingdb": bindingdb_rows}

    all_rows = []
    missing_rows = []
    attempted = 0
    failed = 0

    for compound_id, values in iterate_values(compounds_ds, columns=["smiles"]):
        smiles = values["smiles"]
        if smiles is None or (isinstance(smiles, float) and pd.isna(smiles)) or not str(smiles).strip():
            print(f"WARNING: {compound_id} has no SMILES, skipping", file=sys.stderr)
            missing_rows.append({"id": compound_id, "removed_by": TOOL_NAME,
                                 "kind": "filter", "cause": "no smiles"})
            continue

        smiles = str(smiles).strip()
        attempted += 1

        compound_rows = []
        errors = []
        for source in config["sources"]:
            try:
                compound_rows.extend(fetchers[source](compound_id, smiles, config))
            except Exception as exc:
                errors.append(f"{source}: {exc}")
                print(f"WARNING: {compound_id} failed against {source}: {exc}",
                      file=sys.stderr)

        if errors and not compound_rows:
            failed += 1
            missing_rows.append({"id": compound_id, "removed_by": TOOL_NAME,
                                 "kind": "failure", "cause": "; ".join(errors)})
            continue

        kept = apply_filters(compound_rows, config)
        if kept:
            all_rows.extend(kept)
        else:
            # The id was queried successfully and simply has no matching record.
            # It stays in the table as an all-NaN row so the table remains a
            # complete matrix over the input ids.
            all_rows.append({column: None for column in AFFINITY_COLUMNS}
                            | {"id": compound_id, "smiles": smiles})
            missing_rows.append({"id": compound_id, "removed_by": TOOL_NAME,
                                 "kind": "filter", "cause": "no matching affinity record"})

    if attempted == 0:
        print("ERROR: no compounds with a SMILES were attempted", file=sys.stderr)
        sys.exit(1)

    if failed == attempted:
        print(f"ERROR: all {attempted} compounds failed to query", file=sys.stderr)
        sys.exit(1)

    # Provenance columns: which services answered, which ChEMBL release, and when.
    queried_utc = datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
    release = chembl_release() if 'chembl' in config['sources'] else ''

    affinities = pd.DataFrame(all_rows, columns=AFFINITY_COLUMNS)
    os.makedirs(os.path.dirname(config["affinities_csv"]), exist_ok=True)
    affinities["queried_utc"] = queried_utc
    affinities["chembl_release"] = release
    affinities.to_csv(config["affinities_csv"], index=False)

    targets = build_targets(affinities)
    targets["queried_utc"] = queried_utc
    targets["chembl_release"] = release
    targets.to_csv(config["targets_csv"], index=False)

    local_missing = config["local_missing_csv"]
    os.makedirs(os.path.dirname(local_missing), exist_ok=True)
    pd.DataFrame(missing_rows, columns=["id", "removed_by", "kind", "cause"]).to_csv(
        local_missing, index=False)

    measured = int(affinities["affinity_type"].notna().sum())
    print(f"Wrote {measured} affinity records over {len(targets)} targets "
          f"for {attempted} compounds to {config['affinities_csv']}")


if __name__ == "__main__":
    main()
