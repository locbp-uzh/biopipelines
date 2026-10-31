#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""The provenance of each query a config step builds, carried to the postprocess that writes the result rows.

The config step knows which input of every axis a query was built from: the same selections `predict_single_output_id` composes the id from.
The postprocess reads those cells back by the query's exact id and stamps them as `<axis>.id` columns, so a `+`-joined id is never split to recover its parts.
"""

import csv
import os
from typing import Dict, List, Optional

COLUMN_SUFFIX = ".id"


def write_query_provenance(path: str, provenance: Dict[str, Dict[str, str]]) -> None:
    """One row per query: its id, then one `<axis>.id` cell per provenance key (`predict_single_provenance`)."""
    keys = list(dict.fromkeys(key for cells in provenance.values() for key in cells))
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["id", *[key + COLUMN_SUFFIX for key in keys]])
        for query_id, cells in provenance.items():
            writer.writerow([query_id, *[cells.get(key, "") for key in keys]])


def read_query_provenance(path: Optional[str]) -> Dict[str, Dict[str, str]]:
    """{query id: {"<axis>.id": input id}}; empty when there is no file or it has no provenance column."""
    if not path or not os.path.exists(path):
        return {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        columns = [name for name in (reader.fieldnames or []) if name != "id"]
        if not columns:
            return {}
        return {row["id"]: {name: row[name] for name in columns} for row in reader}


def provenance_columns(provenance: Dict[str, Dict[str, str]]) -> List[str]:
    """The `<axis>.id` column names, in the order the queries first carried them."""
    return list(dict.fromkeys(name for cells in provenance.values() for name in cells))
