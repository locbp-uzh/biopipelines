"""Tests for mmCIF parsing in the shared structure reader (pdb_parser).

parse_pdb_file() dispatches .cif/_atom_site input to parse_cif_file(); these pin
the two failure modes that line-based naive parsing reintroduces: loop rows that
wrap across physical lines, and CIF null markers (./?) that must not block the
auth_*->label_* fallback. parse_models_file() shares the same loop reader, so the
multi-model path is checked too."""

import os

from biopipelines.pdb_parser import (
    parse_pdb_file, parse_cif_file, parse_models_file, read_cif_loop, _cif_value,
    _tokenize_cif_line,
)


_HEADER = """data_test
loop_
_atom_site.group_PDB
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_atom_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.type_symbol
_atom_site.pdbx_PDB_model_num
"""


def _write(tmp_path, name, body):
    p = os.path.join(str(tmp_path), name)
    with open(p, "w") as f:
        f.write(_HEADER + body)
    return p


def test_single_line_rows(tmp_path):
    body = (
        "ATOM CA ALA A 1 CA ALA A 50 1.0 2.0 3.0 C 1\n"
        "ATOM CB ALA A 1 CB ALA A 50 1.5 2.5 3.5 C 1\n"
    )
    atoms = parse_pdb_file(_write(tmp_path, "x.cif", body))
    assert len(atoms) == 2
    # auth_seq_id (50) is used, not label_seq_id (1).
    assert atoms[0].res_num == 50
    assert atoms[0].res_name == "ALA"
    assert atoms[0].chain == "A"


def test_wrapped_row_across_lines(tmp_path):
    """A single loop row split across physical lines must still parse (Finding 1)."""
    body = (
        "ATOM CA ALA A 1 CA ALA A 50\n"
        "1.0 2.0 3.0 C 1\n"
        "ATOM CB ALA A 1 CB ALA A 50 1.5 2.5 3.5 C 1\n"
    )
    atoms = parse_cif_file(_write(tmp_path, "wrap.cif", body))
    assert len(atoms) == 2
    assert {a.atom_name for a in atoms} == {"CA", "CB"}
    assert atoms[0].res_num == 50


def test_null_marker_falls_back_to_label(tmp_path):
    """auth_seq_id '?' must fall back to label_seq_id, not drop the atom (Finding 2)."""
    body = "ATOM CA ALA A 7 CA ALA A ? 1.0 2.0 3.0 C 1\n"
    atoms = parse_cif_file(_write(tmp_path, "null.cif", body))
    assert len(atoms) == 1
    assert atoms[0].res_num == 7  # from label_seq_id


def test_null_chain_falls_back(tmp_path):
    body = "ATOM CA ALA A 7 CA ALA . 7 1.0 2.0 3.0 C 1\n"
    atoms = parse_cif_file(_write(tmp_path, "nullchain.cif", body))
    assert len(atoms) == 1
    assert atoms[0].chain == "A"  # auth '.' -> label_asym_id


def test_multi_model_cif(tmp_path):
    body = (
        "ATOM CA ALA A 1 CA ALA A 1 1.0 2.0 3.0 C 1\n"
        "ATOM CA GLY A 2 CA GLY A 2 4.0 5.0 6.0 C 1\n"
        "ATOM CA ALA A 1 CA ALA A 1 1.1 2.1 3.1 C 2\n"
        "ATOM CA GLY A 2 CA GLY A 2 4.1 5.1 6.1 C 2\n"
    )
    models = parse_models_file(_write(tmp_path, "multi.cif", body), records=("ATOM",))
    assert len(models) == 2
    assert [len(m) for m in models] == [2, 2]
    assert abs(models[0][0].x - 1.0) < 1e-9
    assert abs(models[1][0].x - 1.1) < 1e-9


def test_record_filter_excludes_hetatm(tmp_path):
    body = (
        "ATOM CA ALA A 1 CA ALA A 1 1.0 2.0 3.0 C 1\n"
        "HETATM C1 LIG B 1 C1 LIG B 901 9.0 9.0 9.0 C 1\n"
    )
    models = parse_models_file(_write(tmp_path, "het.cif", body), records=("ATOM",))
    assert len(models) == 1
    assert all(a.res_name == "ALA" for a in models[0])


def test_atom_name_prefers_auth(tmp_path):
    """atom_name uses auth_atom_id (matches fixed-column PDB) when it differs."""
    # label_atom_id=CA1, auth_atom_id=CA -> expect CA.
    body = "ATOM CA1 ALA A 1 CA ALA A 1 1.0 2.0 3.0 C 1\n"
    atoms = parse_cif_file(_write(tmp_path, "name.cif", body))
    assert len(atoms) == 1
    assert atoms[0].atom_name == "CA"


def test_unquoted_prime_atom_id(tmp_path):
    """An unquoted nucleic-acid atom id like O5' must parse, not abort the loop."""
    body = (
        "ATOM O5' G A 1 O5' G A 1 1.0 2.0 3.0 O 1\n"
        "ATOM C1' G A 1 C1' G A 1 1.5 2.5 3.5 C 1\n"
    )
    atoms = parse_cif_file(_write(tmp_path, "rna.cif", body))
    assert {a.atom_name for a in atoms} == {"O5'", "C1'"}


def test_tokenize_prime_and_quotes():
    # Unquoted prime stays one token (shlex would raise here).
    assert _tokenize_cif_line("ATOM O5' G A 1") == ["ATOM", "O5'", "G", "A", "1"]
    # A double-quoted token with an internal prime keeps the prime.
    assert _tokenize_cif_line('a "O5\'" b') == ["a", "O5'", "b"]
    # A quoted token with an internal space is one token.
    assert _tokenize_cif_line("a 'two words' b") == ["a", "two words", "b"]


def test_cif_value_skips_nulls():
    assert _cif_value({"a": "?", "b": "5"}, "a", "b") == "5"
    assert _cif_value({"a": ".", "b": "X"}, "a", "b") == "X"
    assert _cif_value({"a": "."}, "a") is None


def test_read_cif_loop_token_accumulation():
    text = _HEADER + "ATOM CA ALA A 1 CA ALA A 50\n1.0 2.0 3.0 C 1\n"
    recs = read_cif_loop(text, "_atom_site.")
    assert len(recs) == 1
    assert recs[0]["auth_seq_id"] == "50"
    assert recs[0]["Cartn_x"] == "1.0"
