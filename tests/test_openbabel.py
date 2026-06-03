"""Config-time stream-prediction tests for the OpenBabel converter.

OpenBabel routes a 3-D conversion (convert_3d) to a structures stream and 1-D
notations (convert_1d) to columns on the compounds stream. The obabel binary is
not exercised here (no local install); these assert the predicted stream shapes
and the validation rules, which is the framework contract.
"""

import pytest


def _compounds(ids=("aspirin",)):
    from biopipelines.datastream import DataStream
    return DataStream(name="compounds", ids=list(ids), files=[],
                      map_table="compounds.csv", format="csv")


def test_convert_3d_emits_structures_and_compounds_passthrough(record_case):
    from biopipelines.openbabel import OpenBabel

    ob = OpenBabel(compounds=_compounds(), convert_3d="sdf")
    ob.validate_params()
    out = ob.get_output_files()
    structures, compounds = out["structures"], out["compounds"]
    record_case(input="OpenBabel(compounds, convert_3d='sdf')",
                expected=("sdf structures with files", "csv compounds files=[]"),
                actual=(structures.format, compounds.format))
    assert structures.format == "sdf"
    assert len(structures.files) == 1 and structures.files[0].endswith("<id>.sdf")
    assert structures.map_table
    # chemistry passthrough: value-based compounds
    assert compounds.format == "csv"
    assert list(compounds.files) == []


def test_convert_1d_adds_compounds_columns_no_structures(record_case):
    from biopipelines.openbabel import OpenBabel

    ob = OpenBabel(compounds=_compounds(), convert_1d=["smi", "inchi"])
    ob.validate_params()
    out = ob.get_output_files()
    record_case(input="OpenBabel(compounds, convert_1d=['smi','inchi'])",
                expected=("no structures", "csv compounds"),
                actual=(len(out["structures"].files), out["compounds"].format))
    assert len(out["structures"].files) == 0
    assert out["compounds"].format == "csv"
    assert list(out["compounds"].files) == []


def test_requires_at_least_one_conversion(record_case):
    from biopipelines.openbabel import OpenBabel

    record_case(input="OpenBabel(compounds) — no convert_3d/convert_1d",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError, match="convert_3d.*convert_1d"):
        OpenBabel(compounds=_compounds()).validate_params()


def test_convert_3d_rejects_line_format(record_case):
    from biopipelines.openbabel import OpenBabel

    record_case(input="OpenBabel(convert_3d='smi')",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError, match="convert_3d"):
        OpenBabel(compounds=_compounds(), convert_3d="smi").validate_params()


def test_convert_1d_rejects_3d_format(record_case):
    from biopipelines.openbabel import OpenBabel

    record_case(input="OpenBabel(convert_1d=['sdf'])",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError, match="convert_1d"):
        OpenBabel(compounds=_compounds(), convert_1d=["sdf"]).validate_params()
