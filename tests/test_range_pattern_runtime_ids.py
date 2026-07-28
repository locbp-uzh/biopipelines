"""A range pattern states the design space; at runtime the map_table states what exists.

ids_expanded already read the map for LAZY ids (bracket segments, `design[_<chain>]`), but
a RANGE pattern (`12_Panda_1_<1..200>_<1..2>`) took the arithmetic branch and re-declared
every id the design space allows -- including the ones an upstream filter had dropped.

MMseqs2 died on that: its sequences stream declared 2400 ids from six range patterns while
the map held 159 rows, so materializing the filtered CSV raised

    KeyError: ids absent from map_table for 'sequences': ['12_Panda_1_1_1', ...]

and Boltz2 then failed for want of MSAs. The strictness there is correct -- the stream was
lying about its contents.
"""

import json

import pandas as pd

from biopipelines.biopipelines_io import load_datastream, write_filtered_map_table
from biopipelines.datastream import DataStream


def _runtime_stream(tmp_path, pattern, rows):
    mp = tmp_path / "map.csv"
    pd.DataFrame({"id": rows, "sequence": ["A" * (i + 1) for i in range(len(rows))]}).to_csv(mp, index=False)
    p = tmp_path / "ds.json"
    p.write_text(json.dumps({"name": "sequences", "ids": [pattern], "files": [],
                             "map_table": str(mp), "format": "csv"}))
    return load_datastream(str(p))


def test_range_pattern_resolves_to_map_rows(tmp_path):
    survivors = ["12_Panda_1_3_1", "12_Panda_1_3_2", "12_Panda_1_7_1"]
    ds = _runtime_stream(tmp_path, "12_Panda_1_<1..200>_<1..2>", survivors)
    assert ds.ids_expanded == survivors, "must not re-declare filtered-out ids"


def test_filtered_map_table_materializes(tmp_path):
    """The consumer path that MMseqs2 uses must now succeed."""
    survivors = ["12_Panda_1_3_1", "12_Panda_1_7_1"]
    ds = _runtime_stream(tmp_path, "12_Panda_1_<1..200>_<1..2>", survivors)
    out = str(tmp_path / "filtered.csv")
    write_filtered_map_table(ds, out)
    assert pd.read_csv(out)["id"].tolist() == survivors


def test_config_time_expansion_is_unchanged():
    """Without a map there is nothing to consult; the pattern still expands arithmetically."""
    ds = DataStream(name="sequences", ids=["x_<1..3>"], files=[], map_table="", format="csv")
    assert ds.ids_expanded == ["x_1", "x_2", "x_3"]
