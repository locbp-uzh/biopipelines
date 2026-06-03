#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""UniProt fetch helper. Reads a config JSON with a list of accessions, hits the
UniProt REST API, and writes a content-bearing sequences CSV + annotations CSV.

Endpoints used:
  - GET https://rest.uniprot.org/uniprotkb/{acc}.fasta  (sequence)
  - GET https://rest.uniprot.org/uniprotkb/{acc}.json   (annotations)
"""

import argparse
import json
import os
import sys
import time
from typing import List

import pandas as pd
import requests

UNIPROT = "https://rest.uniprot.org/uniprotkb"
SEQ_COLUMNS = ["id", "sequence", "type", "length"]
ANN_COLUMNS = ["id", "name", "organism", "taxonomy_id", "length",
               "go_terms", "pfam", "ec_number", "reviewed"]


def fetch_fasta(acc: str):
    """Return (id, sequence) for the canonical record of a single accession,
    or None if the accession is not found. The returned id is `acc` itself
    so downstream rows are keyed to the declared accessions."""
    url = f"{UNIPROT}/{acc}.fasta"
    r = requests.get(url, timeout=30)
    if r.status_code == 404:
        return None
    r.raise_for_status()
    seq_parts: List[str] = []
    seen_header = False
    for line in r.text.splitlines():
        if line.startswith(">"):
            if seen_header:
                break
            seen_header = True
            continue
        seq_parts.append(line.strip())
    if not seen_header:
        return None
    return acc, "".join(seq_parts)


def fetch_annotations(acc: str) -> dict:
    url = f"{UNIPROT}/{acc}.json"
    r = requests.get(url, timeout=30)
    if r.status_code == 404:
        return {}
    r.raise_for_status()
    data = r.json()

    name = ""
    pname = data.get("proteinDescription", {}).get("recommendedName", {})
    if pname:
        name = pname.get("fullName", {}).get("value", "")

    organism = data.get("organism", {}).get("scientificName", "")
    tax_id = data.get("organism", {}).get("taxonId", "")
    length = data.get("sequence", {}).get("length", "")
    reviewed = data.get("entryType", "").lower().startswith("uniprotkb reviewed")

    go = []
    pfam = []
    for ref in data.get("uniProtKBCrossReferences", []):
        db = ref.get("database", "")
        if db == "GO":
            go.append(ref.get("id", ""))
        elif db == "Pfam":
            pfam.append(ref.get("id", ""))

    ec_numbers = []
    for ec in pname.get("ecNumbers", []) if pname else []:
        ec_numbers.append(ec.get("value", ""))

    return {
        "name": name,
        "organism": organism,
        "taxonomy_id": tax_id,
        "length": length,
        "go_terms": ";".join(go),
        "pfam": ";".join(pfam),
        "ec_number": ";".join(ec_numbers),
        "reviewed": reviewed,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", required=True)
    p.add_argument("--sequences-csv", required=True)
    p.add_argument("--annotations-csv", required=True)
    args = p.parse_args()

    with open(args.config, encoding="utf-8") as f:
        cfg = json.load(f)
    accessions = cfg["accessions"]

    seq_rows, ann_rows, failed = [], [], []
    for acc in accessions:
        try:
            pair = fetch_fasta(acc)
            if pair is None:
                print(f"WARNING: {acc} returned no FASTA record", file=sys.stderr)
                failed.append(acc)
                continue
            sid, seq = pair
            seq_rows.append({"id": sid, "sequence": seq, "type": "protein", "length": len(seq)})

            ann = fetch_annotations(acc)
            ann_rows.append({"id": acc, **{k: ann.get(k, "") for k in
                                           ["name", "organism", "taxonomy_id", "length",
                                            "go_terms", "pfam", "ec_number", "reviewed"]}})
            print(f"  {acc}: {ann.get('name', '<no name>')}")
        except requests.RequestException as e:
            print(f"WARNING: {acc} fetch failed: {e}", file=sys.stderr)
            failed.append(acc)
        time.sleep(0.05)

    os.makedirs(os.path.dirname(args.sequences_csv), exist_ok=True)
    os.makedirs(os.path.dirname(args.annotations_csv), exist_ok=True)
    pd.DataFrame(seq_rows, columns=SEQ_COLUMNS).to_csv(args.sequences_csv, index=False)
    pd.DataFrame(ann_rows, columns=ANN_COLUMNS).to_csv(args.annotations_csv, index=False)
    print(f"Sequences: {args.sequences_csv} ({len(seq_rows)} rows)")
    print(f"Annotations: {args.annotations_csv} ({len(ann_rows)} rows)")

    if failed:
        print(f"Failed: {len(failed)}/{len(accessions)}: {failed}", file=sys.stderr)
    if not seq_rows:
        sys.exit(1)


if __name__ == "__main__":
    main()
