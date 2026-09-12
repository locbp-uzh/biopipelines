"""The documented `(TableInfo, "column")` form must bind on every tool that advertises it."""

import pytest


def _csv(tmp_path):
    import pandas as pd

    path = tmp_path / "positions.csv"
    pd.DataFrame({"id": ["s1"], "designed": ["A1-5"]}).to_csv(path, index=False)
    return str(path)


def _pdb(isolated_cwd):
    path = isolated_cwd / "s1.pdb"
    path.write_text(
        "HEADER test\n"
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "END\n")
    return str(path)


def _built(pipeline, tool_name):
    return [t for t in pipeline.tools if t.TOOL_NAME == tool_name][-1]


@pytest.mark.parametrize("tool_name", ["ProteinMPNN", "LigandMPNN"])
def test_tuple_from_a_table_binds_the_same_reference_as_attribute_access(
        tool_name, local_config, isolated_cwd, new_pipeline, tmp_path, record_case):
    import biopipelines as bp
    from biopipelines.entities import PDB, Table

    kwargs = {"ligand": "LIG"} if tool_name == "LigandMPNN" else {}
    tool_cls = getattr(bp, tool_name)

    pipeline = new_pipeline(f"tuple_{tool_name.lower()}")
    with pipeline:
        table = Table(_csv(tmp_path), table_name="positions").tables.positions
        tool_cls(structures=PDB(_pdb(isolated_cwd)),
                 redesigned=(table, "designed"), **kwargs)
        pipeline.save()
    from_tuple = str(_built(pipeline, tool_name).redesigned)

    pipeline2 = new_pipeline(f"attr_{tool_name.lower()}")
    with pipeline2:
        table2 = Table(_csv(tmp_path), table_name="positions").tables.positions
        tool_cls(structures=PDB(_pdb(isolated_cwd)),
                 redesigned=table2.designed, **kwargs)
        pipeline2.save()
    from_attr = str(_built(pipeline2, tool_name).redesigned)

    record_case(input=f"{tool_name}(redesigned=(table, 'designed'))",
                expected=from_attr, actual=from_tuple)
    assert from_tuple.startswith("TABLE_REFERENCE:")
    assert from_tuple == from_attr


@pytest.mark.parametrize("tool_name", ["ProteinMPNN", "LigandMPNN"])
def test_tuple_from_an_upstream_tools_table_binds(
        tool_name, local_config, isolated_cwd, new_pipeline, record_case):
    import biopipelines as bp

    kwargs = {"ligand": "LIG"} if tool_name == "LigandMPNN" else {}
    tool_cls = getattr(bp, tool_name)

    pipeline = new_pipeline(f"upstream_{tool_name.lower()}")
    with pipeline:
        rfd = bp.RFdiffusion(contigs="50-100", num_designs=2)
        tool_cls(structures=rfd,
                 redesigned=(rfd.tables.structures, "designed"), **kwargs)
        pipeline.save()

    bound = str(_built(pipeline, tool_name).redesigned)
    record_case(input=f"{tool_name}(redesigned=(rfd.tables.structures, 'designed'))",
                expected="TABLE_REFERENCE:...:designed", actual=bound)
    assert bound.startswith("TABLE_REFERENCE:")
    assert bound.endswith(":designed")


@pytest.mark.parametrize("tool_name", ["ProteinMPNN", "LigandMPNN"])
def test_a_literal_selection_and_the_default_are_untouched(
        tool_name, local_config, isolated_cwd, new_pipeline):
    import biopipelines as bp
    from biopipelines.entities import PDB

    kwargs = {"ligand": "LIG"} if tool_name == "LigandMPNN" else {}
    tool_cls = getattr(bp, tool_name)

    pipeline = new_pipeline(f"literal_{tool_name.lower()}")
    with pipeline:
        pdbs = PDB(_pdb(isolated_cwd))
        tool_cls(structures=pdbs, redesigned="A1-5", **kwargs)
        tool_cls(structures=pdbs, **kwargs)
        pipeline.save()

    built = [t for t in pipeline.tools if t.TOOL_NAME == tool_name]
    assert built[0].redesigned == "A1-5"
    assert built[1].redesigned == ""


def test_a_malformed_reference_is_refused_at_construction(local_config, tmp_path):
    """A bad column reference must name the parameter, not surface as a json TypeError from generate_script."""
    from biopipelines.base_config import resolve_table_reference

    with pytest.raises(ValueError, match="column of type int"):
        resolve_table_reference(("some/path.csv", 5), "redesigned")
    with pytest.raises(ValueError, match="first element must be a TableInfo or path"):
        resolve_table_reference((object(), "designed"), "redesigned")
    with pytest.raises(ValueError, match="redesigned must be a string"):
        resolve_table_reference(object(), "redesigned")


def test_the_passthrough_forms_are_returned_unchanged(local_config):
    from biopipelines.base_config import resolve_table_reference
    from biopipelines.biopipelines_io import TableReference

    ref = TableReference("m.csv", "designed")
    assert resolve_table_reference(None) is None
    assert resolve_table_reference("A1-5") == "A1-5"
    assert resolve_table_reference(ref) is ref
    assert str(resolve_table_reference(("m.csv", "designed"))) == str(ref)


def test_openmm_selections_accept_the_documented_tuple(
        local_config, isolated_cwd, new_pipeline, tmp_path, record_case):
    """`mobile_selection` / `frozen_selection` advertise the tuple; a validator used to refuse it."""
    import biopipelines as bp
    from biopipelines.entities import PDB, Table

    pipeline = new_pipeline("openmm_tuple")
    with pipeline:
        table = Table(_csv(tmp_path), table_name="positions").tables.positions
        bp.OpenMM(structures=PDB(_pdb(isolated_cwd)), mobile_selection=(table, "designed"))
        bp.OpenMM(structures=PDB(_pdb(isolated_cwd)), frozen_selection=(table, "designed"))
        pipeline.save()

    built = [t for t in pipeline.tools if t.TOOL_NAME == "OpenMM"]
    record_case(input="OpenMM(mobile_selection=(table, 'designed'))",
                expected="TABLE_REFERENCE:...:designed",
                actual=str(built[0].mobile_selection))
    assert str(built[0].mobile_selection).endswith(":designed")
    assert str(built[1].frozen_selection).endswith(":designed")


def test_ligand_atom_selector_restrict_to_accepts_the_documented_tuple(
        local_config, isolated_cwd, new_pipeline, tmp_path):
    import biopipelines as bp
    from biopipelines.entities import PDB, Table

    pipeline = new_pipeline("las_tuple")
    with pipeline:
        table = Table(_csv(tmp_path), table_name="positions").tables.positions
        bp.LigandAtomSelector(structures=PDB(_pdb(isolated_cwd)), ligand="LIG",
                              atoms="C61+S57", restrict_to=(table, "designed"))
        bp.LigandAtomSelector(structures=PDB(_pdb(isolated_cwd)), ligand="LIG",
                              atoms="C61+S57")
        pipeline.save()

    built = [t for t in pipeline.tools if t.TOOL_NAME == "LigandAtomSelector"]
    assert str(built[0].restrict_to_selection).endswith(":designed")
    assert built[1].restrict_to_selection is None
