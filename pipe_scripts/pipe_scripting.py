#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Execution-phase runner for the Scripting tool.

Loads the user's script, builds runtime input proxies (resolved DataStreams and table references) and an outputs writer, calls ``execution(inputs, outputs)``, then materializes the declared stream map_tables and standalone tables.

Runs in one process under the step's env with the repo on ``sys.path``, so the user's ``execution`` may import biopipelines freely when that env carries its deps. Every declared output — file stream, value stream, or standalone table — is a handle reached by ``outputs[name]``; all three accumulate rows and are written once at the end. ``execution`` returns nothing.
"""

import argparse
import importlib.util
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pandas as pd

from biopipelines.biopipelines_io import (
    load_datastream, iterate_files, iterate_values,
    load_table, lookup_table_value,
)


class _StreamInput:
    """Runtime view of a stream input: iterate files or values, resolve one."""

    def __init__(self, ds):
        self._ds = ds
        self.ids = ds.ids_expanded
        self.name = ds.name
        self.format = ds.format
        self.map_table = ds.map_table
        self.metadata = ds.metadata

    def iterate(self):
        """(id, file_path) for file-based streams; (id, value_dict) otherwise."""
        if self._ds.files:
            yield from iterate_files(self._ds)
        else:
            yield from iterate_values(self._ds)


class _TableInput:
    """Runtime view of a table-column input: per-id lookup over one column."""

    def __init__(self, reference):
        self._table, self._column = load_table(reference)
        self.ids = list(self._table["id"]) if "id" in self._table.columns else []

    def value(self, item_id):
        return lookup_table_value(self._table, item_id, self._column)


class _TableFullInput:
    """Runtime view of a whole-table input: all rows and columns.

    ``ids`` are the table's id column; ``columns`` its columns. ``value(id, col)`` looks one cell up (framework id-matching), ``row(id)`` the whole row as a dict, and ``rows()`` iterates ``(id, row_dict)``.
    """

    def __init__(self, path):
        self._table, _ = load_table(path)
        self.columns = list(self._table.columns)
        self.ids = list(self._table["id"]) if "id" in self._table.columns else []

    def value(self, item_id, column):
        return lookup_table_value(self._table, item_id, column)

    def row(self, item_id):
        return {c: lookup_table_value(self._table, item_id, c) for c in self.columns}

    def rows(self):
        for _, r in self._table.iterrows():
            yield (r["id"], dict(r))


class _OutputHandle:
    """One declared output (file stream, value stream, or table).

    All three accumulate rows and write one CSV at ``target`` on close. A file stream additionally allocates per-id paths under ``folder`` via ``file()``.
    """

    def __init__(self, name, kind, target, folder=None):
        self.name = name
        self.kind = kind  # "file_stream" | "value_stream" | "table"
        self._target = target
        self._folder = folder
        self._rows = []
        self._adopted = None  # a CSV handed over wholesale via `= path`

    def file(self, id, filename):
        """Allocate ``filename`` under the stream folder for id ``id``.

        Returns the absolute path to write into and records the ``{id, file}`` map row. File-based streams only.
        """
        if self.kind != "file_stream":
            raise ValueError(
                f"output '{self.name}' is not a file-based stream; use .row()/.add()"
            )
        os.makedirs(self._folder, exist_ok=True)
        path = os.path.join(self._folder, filename)
        self._rows.append({"id": id, "file": path})
        return path

    def add(self, id, file=None, value=None, **kwargs):
        """Append one row keyed by ``id`` with any of file/value/extra columns."""
        row = {"id": id}
        if file is not None:
            row["file"] = file
        if value is not None:
            row["value"] = value
        row.update(kwargs)
        self._rows.append(row)

    def row(self, mapping):
        """Append one row from a mapping (must include ``id``)."""
        m = dict(mapping)
        if "id" not in m:
            raise ValueError(f"output '{self.name}': row() mapping must include 'id'")
        self._rows.append(m)

    def dataframe(self, df):
        """Hand over a whole DataFrame as this output's rows."""
        self._rows.extend(df.to_dict("records"))

    def _write(self, columns=None):
        if self._adopted is not None:
            df = pd.read_csv(self._adopted)
        elif self._rows:
            df = pd.DataFrame(self._rows)
        else:
            df = pd.DataFrame(columns=columns or ["id"])
        os.makedirs(os.path.dirname(self._target), exist_ok=True)
        df.to_csv(self._target, index=False)


class _Outputs:
    """Mapping of output name -> ``_OutputHandle``, reached via ``outputs[name]``.

    Assigning ``outputs[name] = "path.csv"`` adopts that CSV wholesale (copied to the managed target on close) instead of accumulating rows. ``drop(id)`` records an intentionally-filtered id so the completion check excuses it.
    """

    def __init__(self, manifest):
        self._handles = {}
        for name, info in manifest["streams"].items():
            kind = "value_stream" if info["format"] == "csv" else "file_stream"
            self._handles[name] = _OutputHandle(
                name, kind, info["map_table"], folder=info["folder"]
            )
        for name, info in manifest["tables"].items():
            self._handles[name] = _OutputHandle(name, "table", info["path"])
        self._table_columns = {n: i["columns"] for n, i in manifest["tables"].items()}
        self._removed_by = manifest.get("tool_name", "Scripting")
        self._missing_csv = manifest.get("missing_csv")
        self._dropped = []

    def __getitem__(self, name):
        if name not in self._handles:
            raise ValueError(f"'{name}' is not a declared output")
        return self._handles[name]

    def __setitem__(self, name, csv_path):
        if name not in self._handles:
            raise ValueError(f"'{name}' is not a declared output")
        if not isinstance(csv_path, str):
            raise ValueError(
                f"outputs['{name}'] = ... expects a CSV path string, got "
                f"{type(csv_path).__name__}"
            )
        self._handles[name]._adopted = csv_path

    def drop(self, id, cause=""):
        """Record ``id`` as intentionally filtered out.

        Its absent outputs are then excused by the completion check instead of being flagged as a failure. ``cause`` is a free-text reason.
        """
        self._dropped.append(
            {"id": id, "removed_by": self._removed_by, "kind": "filter", "cause": cause}
        )

    def _write_missing(self):
        if not self._missing_csv:
            return
        os.makedirs(os.path.dirname(self._missing_csv), exist_ok=True)
        columns = ["id", "removed_by", "kind", "cause"]
        df = pd.DataFrame(self._dropped) if self._dropped else pd.DataFrame(columns=columns)
        df.to_csv(self._missing_csv, index=False)

    def _write_all(self):
        for name, handle in self._handles.items():
            handle._write(columns=self._table_columns.get(name))


def _load_user_module(script_path):
    spec = importlib.util.spec_from_file_location("_bp_scripting_user", script_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _build_inputs(inputs_manifest):
    inputs = {}
    for name, spec in inputs_manifest.items():
        if spec["kind"] == "stream":
            inputs[name] = _StreamInput(load_datastream(spec["payload"]))
        elif spec["kind"] == "table_full":
            inputs[name] = _TableFullInput(spec["payload"])
        else:
            inputs[name] = _TableInput(spec["payload"])
    return inputs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--script", required=True)
    parser.add_argument("--inputs", required=True)
    parser.add_argument("--outputs", required=True)
    args = parser.parse_args()

    with open(args.inputs) as f:
        inputs_manifest = json.load(f)
    with open(args.outputs) as f:
        out_manifest = json.load(f)

    module = _load_user_module(args.script)
    inputs = _build_inputs(inputs_manifest)
    outputs = _Outputs(out_manifest)

    module.execution(inputs, outputs)
    outputs._write_all()
    outputs._write_missing()


if __name__ == "__main__":
    main()
