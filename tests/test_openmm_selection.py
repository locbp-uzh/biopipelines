"""OpenMM's frozen/mobile selections: a bare chain is the whole chain, and selecting nothing is an error."""

import pytest

from biopipelines.sele_utils import sele_to_chains_and_list, sele_to_list


@pytest.mark.parametrize("sele, chains, residues", [
    ("B", {"B"}, []),
    ("A2-3+B", {"B"}, [("A", 2), ("A", 3)]),
    ("A1+C", {"C"}, [("A", 1)]),
    ("b10-11+c", {"c"}, [("b", 10), ("b", 11)]),
])
def test_bare_chain_letter_selects_the_whole_chain(sele, chains, residues, record_case):
    got = sele_to_chains_and_list(sele)
    record_case(input=sele, expected=(chains, residues), actual=got)
    assert got == (chains, residues)


@pytest.mark.parametrize("sele", ["", "nan", "-", "foo", "A1-2-3", "A1+xyz", "A1-x", "BB"])
def test_strict_parse_refuses_what_would_select_nothing(sele, record_case):
    with pytest.raises(ValueError) as exc:
        sele_to_chains_and_list(sele)
    record_case(input=sele, expected="ValueError", actual=str(exc.value))


def test_lenient_parse_is_unchanged(capsys):
    assert sele_to_list("B+A2") == [("A", 2)]
    assert "Could not parse 'B'" in capsys.readouterr().out
    assert sele_to_list("b10") == []
    with pytest.raises(ValueError):
        sele_to_list("A1-x")


def test_strict_errors_name_the_bad_token():
    for sele, token in [("b1+A1-x", "A1-x"), ("A1-2-3", "A1-2-3")]:
        with pytest.raises(ValueError, match=f"Malformed range '{token}'"):
            sele_to_chains_and_list(sele)


def _pdb(isolated_cwd):
    path = isolated_cwd / "s1.pdb"
    path.write_text(
        "HEADER test\n"
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "END\n")
    return str(path)


@pytest.mark.parametrize("param", ["frozen_selection", "mobile_selection"])
def test_openmm_accepts_a_bare_chain(param, local_config, isolated_cwd, new_pipeline):
    import biopipelines as bp
    from biopipelines.entities import PDB

    pipeline = new_pipeline(f"openmm_chain_{param}")
    with pipeline:
        bp.OpenMM(structures=PDB(_pdb(isolated_cwd)), **{param: "B"})
    tool = [t for t in pipeline.tools if t.TOOL_NAME == "OpenMM"][-1]
    assert getattr(tool, param) == "B"


@pytest.mark.parametrize("param", ["frozen_selection", "mobile_selection"])
def test_openmm_refuses_an_unparseable_selection(param, local_config, isolated_cwd, new_pipeline,
                                                 record_case):
    import biopipelines as bp
    from biopipelines.entities import PDB

    with pytest.raises(ValueError) as exc:
        with new_pipeline(f"openmm_bad_{param}"):
            bp.OpenMM(structures=PDB(_pdb(isolated_cwd)), **{param: "A1+chainB"})
    record_case(input=f"{param}='A1+chainB'", expected="ValueError", actual=str(exc.value))
    assert param in str(exc.value)
