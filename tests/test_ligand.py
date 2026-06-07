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
    assert lig._detect_lookup_type("A1LU6") == "ccd"   # extended 5-char CCD
    assert lig._detect_lookup_type("2244") == "cid"
    assert lig._detect_lookup_type("50-78-2") == "cas"
    assert lig._detect_lookup_type("aspirin") == "name"
    assert lig._detect_lookup_type("caffeine") == "name"
    # Short lowercase compound names are PubChem names, not CCD codes.
    assert lig._detect_lookup_type("water") == "name"
    assert lig._detect_lookup_type("urea") == "name"


# ── option validation ────────────────────────────────────────────────────────

def test_auth_ligand_field_prefers_ccd_for_authoritative_rows():
    from biopipelines.ligand_utils import auth_ligand_field

    assert auth_ligand_field({
        "format": "ccd",
        "source": "",
        "ccd": "ATP",
        "smiles": "C1=NC=NC2=C1N=CN2",
    }) == ("ccd", "ATP")
    assert auth_ligand_field({
        "format": "smiles",
        "source": "rcsb",
        "ccd": "ATP",
        "smiles": "C1=NC=NC2=C1N=CN2",
    }) == ("ccd", "ATP")
    assert auth_ligand_field({
        "format": "smiles",
        "source": "pubchem",
        "ccd": "ATP",
        "smiles": "CCO",
    }) == ("smiles", "CCO")


def test_ligand_invalid_source_raises(record_case):
    from biopipelines.ligand import Ligand

    record_case(input="Ligand(source='nowhere')",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError, match="source"):
        Ligand(lookup="ATP", source="nowhere")


def test_ligand_output_format_removed_raises(record_case):
    """output_format is no longer user-selectable; passing it points to OpenBabel."""
    from biopipelines.ligand import Ligand

    record_case(input="Ligand(output_format='cif')",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError, match="output_format was removed"):
        Ligand(lookup="ATP", output_format="cif")


def test_ligand_structures_default_sdf(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """A fetched/generated Ligand emits its structures stream as SDF."""
    from biopipelines.ligand import Ligand

    pipeline = new_pipeline("lig_sdf")
    with pipeline:
        lig = Ligand(lookup="ATP")
        pipeline.save()

    record_case(input="Ligand(lookup='ATP') structures format",
                expected="sdf", actual=lig.streams.structures.format)
    assert lig.streams.structures.format == "sdf"


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


# ── code-only construction ───────────────────────────────────────────────────

def test_ligand_code_only_basic(record_case):
    """Ligand(code='ZIT') — code-only mode: single compounds row, no structures."""
    from biopipelines.ligand import Ligand

    lig = Ligand(code="ZIT")
    record_case(input="Ligand(code='ZIT')",
                expected=(["ZIT"], ["ZIT"], True),
                actual=(list(lig.custom_ids), list(lig.residue_codes), lig.code_only))
    assert lig.code_only is True
    assert lig.custom_ids == ["ZIT"]
    assert lig.residue_codes == ["ZIT"]
    assert lig.lookup_values == [] and lig.smiles_values == []


def test_ligand_code_only_streams(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """Code-only mode: compounds is value-based csv; no structures stream."""
    from biopipelines.ligand import Ligand

    pipeline = new_pipeline("lig_code_only")
    with pipeline:
        lig = Ligand(code="ZIT")
        pipeline.save()

    compounds = lig.streams.compounds
    structures = lig.streams.structures
    record_case(input="Ligand(code='ZIT') streams",
                expected=("csv", [], 0),
                actual=(compounds.format, list(compounds.files), len(structures)))
    assert compounds.format == "csv"
    assert list(compounds.files) == []
    assert compounds.map_table
    assert len(structures) == 0  # no coordinate files in code-only mode


def test_ligand_code_mutually_exclusive(record_case):
    from biopipelines.ligand import Ligand

    record_case(input="Ligand(code=, smiles=) — mutually exclusive",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError, match="mutually exclusive"):
        Ligand(code="ZIT", smiles="CCO")
    with pytest.raises(ValueError, match="mutually exclusive"):
        Ligand(code="ZIT", lookup="ATP")
    with pytest.raises(ValueError, match="mutually exclusive"):
        Ligand(code="ZIT", codes="LIG")


# ── structures stream now carries a map_table ────────────────────────────────

def test_ligand_structures_stream_has_map_and_template(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """A normal Ligand emits a structures stream with a map_table and an
    <id>-templated file under the structures/ folder (not compounds/)."""
    from biopipelines.ligand import Ligand

    pipeline = new_pipeline("lig_structures_map")
    with pipeline:
        lig = Ligand(lookup="ATP")
        pipeline.save()

    structures = lig.streams.structures
    files = list(structures.files)
    record_case(input="Ligand(lookup='ATP') structures stream",
                expected=("map_table set", "<id> template", "under structures/"),
                actual=(bool(structures.map_table), files, files))
    assert structures.map_table  # now has its own lineage map
    assert len(files) == 1 and files[0].endswith("<id>.sdf")
    # coordinate files live in the structures/ stream folder, not compounds/
    assert "structures" in str(structures.map_table)


# ── HETATM extraction (structures= path) ─────────────────────────────────────

def _load_pipe_ligand():
    import importlib.util
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    path = os.path.join(repo, "pipe_scripts", "pipe_ligand.py")
    spec = importlib.util.spec_from_file_location("pipe_ligand", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_extract_hetatm_block_pdb(tmp_path, record_case):
    """The PDB carver keeps only the named ligand's HETATM lines (first copy)."""
    pl = _load_pipe_ligand()
    pdb = tmp_path / "complex.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A   1      11.000  11.000  11.000  1.00  0.00           C\n"
        "HETATM    2  C1  LIG A 200      20.000  20.000  20.000  1.00  0.00           C\n"
        "HETATM    3  O1  LIG A 200      21.000  21.000  21.000  1.00  0.00           O\n"
        "HETATM    4  C1  LIG B 200      30.000  30.000  30.000  1.00  0.00           C\n"
    )
    block = pl.extract_hetatm_block(str(pdb), "LIG")
    record_case(input="extract_hetatm_block(pdb,'LIG')",
                expected="2 HETATM lines (first copy)",
                actual=block)
    assert block is not None
    lines = [l for l in block.splitlines() if l.startswith("HETATM")]
    assert len(lines) == 2  # only chain A copy
    assert pl.extract_hetatm_block(str(pdb), "NOPE") is None


def test_extract_hetatm_block_cif(tmp_path, record_case):
    """The CIF carver selects _atom_site rows by comp_id, keeping the first copy,
    and supports an extended (5-char) code."""
    pl = _load_pipe_ligand()
    cif = tmp_path / "complex.cif"
    cif.write_text(
        "data_complex\n"
        "loop_\n"
        "_atom_site.group_PDB\n"
        "_atom_site.id\n"
        "_atom_site.label_atom_id\n"
        "_atom_site.label_comp_id\n"
        "_atom_site.auth_asym_id\n"
        "_atom_site.auth_seq_id\n"
        "_atom_site.Cartn_x\n"
        "_atom_site.Cartn_y\n"
        "_atom_site.Cartn_z\n"
        "ATOM   1 CA ALA A 1 11.0 11.0 11.0\n"
        "HETATM 2 C1 A7ZK B 200 20.0 20.0 20.0\n"
        "HETATM 3 O1 A7ZK B 200 21.0 21.0 21.0\n"
        "HETATM 4 C1 A7ZK C 200 30.0 30.0 30.0\n"
        "#\n"
    )
    block = pl.extract_hetatm_block_cif(str(cif), "A7ZK")
    record_case(input="extract_hetatm_block_cif(cif,'A7ZK')",
                expected="2 HETATM rows (first copy), valid mini-cif",
                actual=block)
    assert block is not None
    rows = [l for l in block.splitlines() if l.startswith("HETATM")]
    assert len(rows) == 2  # only chain B copy
    assert block.startswith("data_A7ZK")
    assert "_atom_site.label_comp_id" in block
    assert pl.extract_hetatm_block_cif(str(cif), "NOPE") is None


def test_extract_hetatm_block_cif_quoted_values(tmp_path, record_case):
    """The CIF carver matches comp_id correctly even when a row carries a
    space-containing quoted value (shlex-aware split, not a naive .split())."""
    pl = _load_pipe_ligand()
    cif = tmp_path / "quoted.cif"
    # label_atom_id values are single-quoted and contain a space.
    cif.write_text(
        "data_x\n"
        "loop_\n"
        "_atom_site.group_PDB\n"
        "_atom_site.id\n"
        "_atom_site.label_atom_id\n"
        "_atom_site.label_comp_id\n"
        "_atom_site.auth_asym_id\n"
        "_atom_site.auth_seq_id\n"
        "_atom_site.Cartn_x\n"
        "_atom_site.Cartn_y\n"
        "_atom_site.Cartn_z\n"
        "HETATM 1 'C 1' LIG A 200 20.0 20.0 20.0\n"
        "HETATM 2 'O 1' LIG A 200 21.0 21.0 21.0\n"
    )
    block = pl.extract_hetatm_block_cif(str(cif), "LIG")
    record_case(input="extract_hetatm_block_cif quoted atom ids",
                expected="2 rows matched despite spaces in quoted field",
                actual=block)
    assert block is not None
    rows = [l for l in block.splitlines() if l.startswith("HETATM")]
    assert len(rows) == 2
