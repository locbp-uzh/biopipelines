"""RFdiffusion3 must accept a per-structure contig, not only one shared string.

`contig` was typed `str` and copied verbatim into every entry of the inputs
JSON. A contig names a chain range (`1,A1-135,1`), so one literal is only ever
correct for structures that share an end residue. Over a set whose chains ended
anywhere from 118 to 146, a single string would have kept the wrong range on
almost every design — silently, since RFD3 accepts any well-formed contig.

The other three RFdiffusion wrappers already took a table-column reference here;
this is the fourth.
"""

import csv
import json
import os
import subprocess
import sys

import pytest

from biopipelines.biopipelines_io import TableReference
from biopipelines.base_config import TableInfo, resolve_table_reference


def _norm(value):
    """What RFdiffusion3 stores: ("literal", str) or ("table", token).

    The wrapper used to carry its own copy of this normalization, which read
    ``value[0].info.path`` and so AttributeError'd on the documented
    ``(path, "column")`` pair. It now calls the shared helper.
    """
    resolved = resolve_table_reference(value, "contig")
    return ("literal", resolved) if isinstance(resolved, str) else ("table", str(resolved))


def test_literal_stays_literal():
    assert _norm("A1-135") == ("literal", "A1-135")
    assert _norm("") == ("literal", "")


def test_table_reference_forms():
    tr = TableReference("/x/c.csv", "contig_1")
    assert _norm(tr) == ("table", "TABLE_REFERENCE:/x/c.csv:contig_1")
    # (TableInfo, "column") tuple — the LigandMPNN/ThermoMPNN convention.
    info = TableInfo(name="c", path="/x/c.csv", columns=["id", "contig_2"], description="")
    assert _norm((info, "contig_2")) == ("table", "TABLE_REFERENCE:/x/c.csv:contig_2")


def test_plain_path_column_pair_is_accepted():
    """The documented (path, "column") form. The wrapper's own copy raised
    AttributeError on it, because it assumed a TableInfo and read .info.path."""
    assert _norm(("/x/c.csv", "contig_3")) == ("table", "TABLE_REFERENCE:/x/c.csv:contig_3")


def test_non_string_non_reference_rejected():
    with pytest.raises(ValueError):
        _norm(123)


def test_builder_resolves_one_contig_per_structure(tmp_path):
    """The end each design keeps must come from its own row."""
    for sid in ("dA", "dB"):
        (tmp_path / f"{sid}.pdb").write_text("ATOM      1  CA  ALA A   1\n")

    table = tmp_path / "contigs.csv"
    with open(table, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "contig_1"])
        w.writerow(["dA", "1,A1-135,1"])
        w.writerow(["dB", "1,A1-142,1"])

    ds = tmp_path / "ds.json"
    ds.write_text(json.dumps({
        "name": "structures", "format": "pdb", "ids": ["dA", "dB"],
        "files": [str(tmp_path / "dA.pdb"), str(tmp_path / "dB.pdb")],
    }))
    tmpl = tmp_path / "tmpl.json"
    tmpl.write_text(json.dumps({"num_designs": 2}))
    out = tmp_path / "out.json"

    script = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                          "pipe_scripts", "pipe_rfdiffusion3_build_inputs.py")
    subprocess.run(
        [sys.executable, script, "--template", str(tmpl), "--output", str(out),
         "--structures-json", str(ds),
         "--contig-reference", f"TABLE_REFERENCE:{table}:contig_1"],
        check=True, capture_output=True)

    cfg = json.loads(out.read_text())
    assert cfg["dA"]["contig"] == "1,A1-135,1"
    assert cfg["dB"]["contig"] == "1,A1-142,1"


def test_builder_rejects_an_id_with_no_contig(tmp_path):
    """An unmatched id must fail loudly, not fall back to a shared default."""
    (tmp_path / "dA.pdb").write_text("ATOM      1  CA  ALA A   1\n")
    table = tmp_path / "contigs.csv"
    with open(table, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "contig_1"])
        w.writerow(["someone_else", "1,A1-135,1"])

    ds = tmp_path / "ds.json"
    ds.write_text(json.dumps({
        "name": "structures", "format": "pdb", "ids": ["dA"],
        "files": [str(tmp_path / "dA.pdb")],
    }))
    tmpl = tmp_path / "tmpl.json"
    tmpl.write_text(json.dumps({}))

    script = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                          "pipe_scripts", "pipe_rfdiffusion3_build_inputs.py")
    r = subprocess.run(
        [sys.executable, script, "--template", str(tmpl),
         "--output", str(tmp_path / "out.json"), "--structures-json", str(ds),
         "--contig-reference", f"TABLE_REFERENCE:{table}:contig_1"],
        capture_output=True, text=True)
    assert r.returncode != 0
    # lookup_table_value raises on the unmatched id before the empty-cell check.
    assert "dA" in r.stderr


def test_reference_without_pdb_is_rejected():
    """A per-PDB contig is keyed by pdb id, so it needs the structures."""
    import inspect
    src = inspect.getsource(
        __import__("biopipelines.rfdiffusion3", fromlist=["x"]).RFdiffusion3.validate_params)
    assert "per-PDB contig reference requires an input structure" in src


def test_reference_is_not_baked_into_the_shared_template(local_config, isolated_cwd):
    """The template is copied to every entry; a reference there would broadcast
    one design's contig onto all of them.

    Asserted on the constructed object rather than by grepping __init__, so a
    refactor that keeps the behavior passes and one that loses it does not.
    """
    from biopipelines.datastream import DataStream
    from biopipelines.rfdiffusion3 import RFdiffusion3

    pdb = DataStream(name="structures", ids=["dA"], files=["/tmp/dA.pdb"], format="pdb")

    literal = RFdiffusion3(pdb=pdb, contig="A1-135")
    assert literal.contig == "A1-135"
    assert literal.contig_reference is None

    ref = RFdiffusion3(pdb=pdb, contig=TableReference("/x/c.csv", "contig_1"))
    assert ref.contig == "", "a table token reached the shared template"
    assert ref.contig_reference == "TABLE_REFERENCE:/x/c.csv:contig_1"
