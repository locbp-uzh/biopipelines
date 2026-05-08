"""Unit tests for the Ligand input-type tool.

Covers the four lookup paths (CCD/RCSB, PubChem CID/CAS/name detection), SMILES
input, .txt SMILES file loading, dict-keyed IDs, residue-code defaulting, and
integration with Mock and Panda downstream.

Network calls during configure_inputs are allowed — Ligand checks RCSB/PubChem
for remote lookups at config time.
"""

import os

import pytest


# ── lookup construction ──────────────────────────────────────────────────────

def test_ligand_single_ccd(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Ligand(lookup='ATP') — smoke script, default id and residue code = 'ATP'."""
    from biopipelines.ligand import Ligand

    pipeline = new_pipeline("lig_single_ccd")
    with pipeline:
        lig = Ligand(lookup="ATP")
        script_path = pipeline.save()

    ids = list(lig.streams.structures.ids)
    record_case(input="Ligand(lookup='ATP')",
                expected=["ATP"], actual=ids)
    assert_valid_script(script_path, "Ligand")
    assert ids == ["ATP"]


def test_ligand_list_of_codes(record_case):
    from biopipelines.ligand import Ligand

    lig = Ligand(lookup=["ATP", "GDP", "HEM"])
    record_case(input="Ligand(lookup=['ATP','GDP','HEM'])",
                expected=["ATP", "GDP", "HEM"],
                actual=list(lig.custom_ids))
    assert lig.custom_ids == ["ATP", "GDP", "HEM"]


def test_ligand_dict_uses_keys_as_ids(record_case):
    from biopipelines.ligand import Ligand

    lig = Ligand(lookup={"cofactor": "ATP", "gnp": "GDP"})
    record_case(input="Ligand(lookup={'cofactor':'ATP','gnp':'GDP'})",
                expected=(["cofactor", "gnp"], ["ATP", "GDP"]),
                actual=(list(lig.custom_ids), list(lig.lookup_values)))
    assert lig.custom_ids == ["cofactor", "gnp"]
    assert lig.lookup_values == ["ATP", "GDP"]


def test_ligand_custom_ids_and_codes(record_case):
    """Explicit ids and codes override defaults."""
    from biopipelines.ligand import Ligand

    lig = Ligand(lookup=["ATP", "GDP"],
                 ids=["lig1", "lig2"],
                 codes=["A", "G"])
    record_case(input="Ligand(lookup=[...], ids=[...], codes=[...])",
                expected=(["lig1", "lig2"], ["A", "G"]),
                actual=(list(lig.custom_ids), list(lig.residue_codes)))
    assert lig.custom_ids == ["lig1", "lig2"]
    assert lig.residue_codes == ["A", "G"]


# ── SMILES input ─────────────────────────────────────────────────────────────

def test_ligand_single_smiles(record_case):
    """smiles='CCO' — default id is smiles1, default code is LIG."""
    from biopipelines.ligand import Ligand

    lig = Ligand(smiles="CCO")
    record_case(input="Ligand(smiles='CCO')",
                expected=(["smiles1"], ["LIG"], ["CCO"]),
                actual=(list(lig.custom_ids), list(lig.residue_codes),
                        list(lig.smiles_values)))
    assert lig.custom_ids == ["smiles1"]
    assert lig.residue_codes == ["LIG"]
    assert lig.smiles_values == ["CCO"]


def test_ligand_smiles_dict_uses_keys_as_ids(record_case):
    from biopipelines.ligand import Ligand

    lig = Ligand(smiles={"ethanol": "CCO", "methane": "C"})
    record_case(input="Ligand(smiles={'ethanol':'CCO','methane':'C'})",
                expected=(["ethanol", "methane"], ["CCO", "C"]),
                actual=(list(lig.custom_ids), list(lig.smiles_values)))
    assert lig.custom_ids == ["ethanol", "methane"]
    assert lig.smiles_values == ["CCO", "C"]


def test_ligand_smiles_from_txt_file(isolated_cwd, record_case):
    """SMILES per line in a .txt file; IDs derived from file basename + index."""
    from biopipelines.ligand import Ligand

    txt_path = isolated_cwd / "mylibs.txt"
    txt_path.write_text("CCO\nCC(=O)O\nc1ccccc1\n")

    lig = Ligand(lookup=str(txt_path))
    record_case(input="Ligand(lookup='mylibs.txt')",
                expected=(["mylibs1", "mylibs2", "mylibs3"],
                          ["CCO", "CC(=O)O", "c1ccccc1"]),
                actual=(list(lig.custom_ids), list(lig.smiles_values)))
    assert lig.custom_ids == ["mylibs1", "mylibs2", "mylibs3"]
    assert lig.smiles_values == ["CCO", "CC(=O)O", "c1ccccc1"]


def test_ligand_combines_lookup_and_smiles(record_case):
    """Both lookup and smiles can be supplied together."""
    from biopipelines.ligand import Ligand

    lig = Ligand(lookup=["ATP"], smiles=["CCO"])
    record_case(input="Ligand(lookup=['ATP'], smiles=['CCO'])",
                expected=(["ATP", "smiles1"], ["ATP", "LIG"]),
                actual=(list(lig.custom_ids), list(lig.residue_codes)))
    assert lig.custom_ids == ["ATP", "smiles1"]
    assert lig.residue_codes == ["ATP", "LIG"]


# ── lookup-type detection ────────────────────────────────────────────────────

def test_ligand_detect_lookup_types():
    """_detect_lookup_type classifies CCD / CID / CAS / name correctly."""
    from biopipelines.ligand import Ligand

    lig = Ligand(lookup=["ATP"])  # any construction works; we use the helper
    assert lig._detect_lookup_type("ATP") == "ccd"
    assert lig._detect_lookup_type("HEM") == "ccd"
    assert lig._detect_lookup_type("2244") == "cid"
    assert lig._detect_lookup_type("50-78-2") == "cas"
    assert lig._detect_lookup_type("aspirin") == "name"
    assert lig._detect_lookup_type("caffeine") == "name"


# ── option validation ────────────────────────────────────────────────────────

def test_ligand_invalid_source_raises(record_case):
    from biopipelines.ligand import Ligand

    record_case(input="Ligand(source='nowhere')",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError, match="source"):
        Ligand(lookup="ATP", source="nowhere")


def test_ligand_invalid_output_format_raises(record_case):
    from biopipelines.ligand import Ligand

    record_case(input="Ligand(output_format='sdf')",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError, match="output_format"):
        Ligand(lookup="ATP", output_format="sdf")


def test_ligand_output_format_cif(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """output_format='cif' propagates through the streams.structures format."""
    from biopipelines.ligand import Ligand

    pipeline = new_pipeline("lig_cif")
    with pipeline:
        lig = Ligand(lookup="ATP", output_format="cif")
        pipeline.save()

    record_case(input="Ligand(output_format='cif')",
                expected="cif", actual=lig.streams.structures.format)
    assert lig.streams.structures.format == "cif"


def test_ligand_no_lookup_or_smiles_raises(record_case):
    from biopipelines.ligand import Ligand

    record_case(input="Ligand() — neither lookup nor smiles",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError, match="lookup.*smiles"):
        Ligand()


def test_ligand_ids_length_mismatch_raises(record_case):
    from biopipelines.ligand import Ligand

    record_case(input="len(ids)=1 vs 2 lookups",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError, match="Length mismatch"):
        Ligand(lookup=["ATP", "GDP"], ids=["only_one"])


def test_ligand_codes_length_mismatch_raises(record_case):
    from biopipelines.ligand import Ligand

    record_case(input="len(codes)=1 vs 2 lookups",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError, match="Length mismatch"):
        Ligand(lookup=["ATP", "GDP"], codes=["A"])


def test_ligand_empty_dict_raises(record_case):
    from biopipelines.ligand import Ligand

    record_case(input="Ligand(lookup={})",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError, match="cannot be empty"):
        Ligand(lookup={})


# ── integration: Ligand + Mock ───────────────────────────────────────────────

def test_ligand_feeds_mock(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Mock consumes Ligand.streams.compounds (value-stream with a map_table)
    and fans out with children=<A B>.

    Uses .streams.compounds rather than .streams.structures because Mock's
    source resolution requires a map_table, which Ligand only attaches to the
    compounds stream.
    """
    from biopipelines.ligand import Ligand
    from biopipelines.mock import Mock

    pipeline = new_pipeline("lig_feeds_mock")
    with pipeline:
        lig = Ligand(lookup={"cofactor": "ATP", "nucleotide": "GDP"})
        m = Mock(
            source=lig.streams.compounds,
            children="<A B>",
            streams={"poses": {"format": "pdb", "file": "<id>.pdb"}},
        )
        script_path = pipeline.save()

    ids = list(m.streams.poses.ids)
    record_case(input="Ligand.compounds → Mock(children=<A B>)",
                expected=["cofactor_<A B>", "nucleotide_<A B>"], actual=ids)
    assert ids == ["cofactor_<A B>", "nucleotide_<A B>"]
    assert_valid_script(script_path, "Ligand", "Mock")


# ── integration: Ligand + Panda ──────────────────────────────────────────────

def test_ligand_feeds_panda_on_compounds_table(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda sorts Ligand's compounds table."""
    from biopipelines.ligand import Ligand
    from biopipelines.panda import Panda

    pipeline = new_pipeline("lig_feeds_panda")
    with pipeline:
        lig = Ligand(lookup=["ATP", "GDP"])
        Panda(
            tables=lig.tables.compounds,
            operations=[Panda.sort("id", ascending=True)],
        )
        script_path = pipeline.save()

    record_case(input="Ligand → Panda(sort compounds)",
                expected="script with Ligand+Panda",
                actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Ligand", "Panda", "lig_feeds_panda")


def test_ligand_smiles_feeds_panda(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """SMILES-based Ligand + Panda filter on the compounds table."""
    from biopipelines.ligand import Ligand
    from biopipelines.panda import Panda

    pipeline = new_pipeline("lig_smiles_feeds_panda")
    with pipeline:
        lig = Ligand(smiles={"ethanol": "CCO", "acetate": "CC(=O)O"})
        Panda(
            tables=lig.tables.compounds,
            operations=[Panda.filter("id != 'ethanol'")],
        )
        script_path = pipeline.save()

    record_case(input="Ligand(smiles=dict) → Panda(filter)",
                expected="script",
                actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Ligand", "Panda")


def test_ligand_panda_head(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda.head(n) on Ligand's compounds table."""
    from biopipelines.ligand import Ligand
    from biopipelines.panda import Panda

    pipeline = new_pipeline("lig_panda_head")
    with pipeline:
        lig = Ligand(lookup=["ATP", "GDP", "HEM"])
        Panda(
            tables=lig.tables.compounds,
            operations=[Panda.head(2)],
        )
        script_path = pipeline.save()

    record_case(input="Ligand → Panda(head(2))",
                expected="script", actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Ligand", "Panda")


def test_ligand_panda_tail(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda.tail(n) on Ligand's compounds table."""
    from biopipelines.ligand import Ligand
    from biopipelines.panda import Panda

    pipeline = new_pipeline("lig_panda_tail")
    with pipeline:
        lig = Ligand(lookup=["ATP", "GDP", "HEM"])
        Panda(
            tables=lig.tables.compounds,
            operations=[Panda.tail(1)],
        )
        script_path = pipeline.save()

    record_case(input="Ligand → Panda(tail(1))",
                expected="script", actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Ligand", "Panda")


def test_ligand_panda_filter_sort_head_chain(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """filter + sort + head chain on Ligand.tables.compounds."""
    from biopipelines.ligand import Ligand
    from biopipelines.panda import Panda

    pipeline = new_pipeline("lig_panda_fsh")
    with pipeline:
        lig = Ligand(lookup=["ATP", "GDP", "HEM", "NAD"])
        Panda(
            tables=lig.tables.compounds,
            operations=[
                Panda.filter("id != 'NAD'"),
                Panda.sort("id", ascending=True),
                Panda.head(2),
            ],
        )
        script_path = pipeline.save()

    record_case(input="Ligand → Panda(filter+sort+head)",
                expected="chained script",
                actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Ligand", "Panda")


def test_ligand_panda_select_columns(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """select_columns trims the compounds table to id + smiles."""
    from biopipelines.ligand import Ligand
    from biopipelines.panda import Panda

    pipeline = new_pipeline("lig_panda_select")
    with pipeline:
        lig = Ligand(lookup=["ATP", "GDP"])
        Panda(
            tables=lig.tables.compounds,
            operations=[Panda.select_columns(["id", "smiles"])],
        )
        script_path = pipeline.save()

    record_case(input="Ligand → Panda(select_columns id,smiles)",
                expected="script", actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Ligand", "Panda")


def test_ligand_panda_rename(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """rename remaps code → residue_code on the compounds table."""
    from biopipelines.ligand import Ligand
    from biopipelines.panda import Panda

    pipeline = new_pipeline("lig_panda_rename")
    with pipeline:
        lig = Ligand(lookup=["ATP"])
        Panda(
            tables=lig.tables.compounds,
            operations=[Panda.rename({"code": "residue_code"})],
        )
        script_path = pipeline.save()

    record_case(input="Ligand → Panda(rename code→residue_code)",
                expected="script", actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Ligand", "Panda")


def test_ligand_panda_sample(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda.sample on compounds."""
    from biopipelines.ligand import Ligand
    from biopipelines.panda import Panda

    pipeline = new_pipeline("lig_panda_sample")
    with pipeline:
        lig = Ligand(lookup=["ATP", "GDP", "HEM", "NAD"])
        Panda(
            tables=lig.tables.compounds,
            operations=[Panda.sample(n=2, random_state=11)],
        )
        script_path = pipeline.save()

    record_case(input="Ligand → Panda(sample n=2)",
                expected="script", actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Ligand", "Panda")
