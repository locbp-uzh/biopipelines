#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Execution-phase runner for the Scripting tool.

Loads the user's script, builds runtime input proxies (resolved DataStreams and
table references) and an outputs writer, calls ``execution(inputs, outputs)``,
then materializes the declared stream map_tables and standalone tables.
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
        self.format = ds.format

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
        self.ids = []

    def value(self, item_id):
        return lookup_table_value(self._table, item_id, self._column)


class _Outputs:
    """Writer the user's execution() uses to place files and accumulate rows.

    - ``file(stream, name)`` returns a path under the stream folder (created on
      demand) for the user to write a per-id file into.
    - ``row(table, mapping)`` accumulates a row for a declared standalone table.
    The user's execution() returns the stream map rows; this object owns the
    standalone-table rows.
    """

    def __init__(self, manifest):
        self._streams = manifest["streams"]
        self._tables = manifest["tables"]
        self._table_rows = {name: [] for name in self._tables}

    def file(self, stream, name):
        if stream not in self._streams:
            raise ValueError(f"'{stream}' is not a declared output stream")
        folder = self._streams[stream]["folder"]
        os.makedirs(folder, exist_ok=True)
        return os.path.join(folder, name)

    def row(self, table, mapping):
        if table not in self._tables:
            raise ValueError(f"'{table}' is not a declared output table")
        self._table_rows[table].append(dict(mapping))


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
        else:
            inputs[name] = _TableInput(spec["payload"])
    return inputs


def _write_stream_maps(declared_rows, streams_manifest):
    for name, info in streams_manifest.items():
        rows = declared_rows.get(name, [])
        map_path = info["map_table"]
        os.makedirs(os.path.dirname(map_path), exist_ok=True)
        if rows:
            pd.DataFrame(rows).to_csv(map_path, index=False)
        else:
            pd.DataFrame(columns=["id"]).to_csv(map_path, index=False)


def _write_tables(outputs, tables_manifest):
    for name, info in tables_manifest.items():
        rows = outputs._table_rows[name]
        path = info["path"]
        os.makedirs(os.path.dirname(path), exist_ok=True)
        columns = info["columns"]
        if rows:
            pd.DataFrame(rows).to_csv(path, index=False)
        else:
            pd.DataFrame(columns=columns).to_csv(path, index=False)


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

    declared_rows = module.execution(inputs, outputs)
    if declared_rows is None:
        declared_rows = {}
    if not isinstance(declared_rows, dict):
        raise ValueError(
            f"execution() must return a dict of stream rows or None, "
            f"got {type(declared_rows).__name__}"
        )

    _write_stream_maps(declared_rows, out_manifest["streams"])
    _write_tables(outputs, out_manifest["tables"])


if __name__ == "__main__":
    main()
