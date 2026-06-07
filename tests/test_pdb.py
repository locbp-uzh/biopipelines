"""Unit tests for the PDB input-type tool.

Covers construction paths (single string, list, dict, folder, upstream stream),
option handling (convert, ids, operations), validation errors, and integration
with Mock and Panda downstream.

RCSB network calls are intentionally allowed — PDB hits the REST API during
configure_inputs to look up ligand metadata and to verify downloads.
"""

import os

import pytest


# ── single + list + dict construction ────────────────────────────────────────

def test_pdb_single_rcsb_code(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    from biopipelines.pdb import PDB

    pipeline = new_pipeline("pdb_single")
    with pipeline:
        p = PDB(pdbs="4ufc")
        script_path = pipeline.save()

    ids = list(p.streams.structures.ids)
    record_case(input="PDB(pdbs='4ufc')",
                expected=["4ufc"], actual=ids)
    assert_valid_script(script_path, "PDB")
    assert ids == ["4ufc"]


def test_pdb_list_of_codes(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    from biopipelines.pdb import PDB

    pipeline = new_pipeline("pdb_list")
    with pipeline:
        p = PDB(pdbs=["4ufc", "1ake"])
        pipeline.save()

    ids = list(p.streams.structures.ids)
    record_case(input="PDB(pdbs=['4ufc','1ake'])",
                expected=["4ufc", "1ake"], actual=ids)
    assert ids == ["4ufc", "1ake"]


def test_pdb_dict_uses_keys_as_ids(record_case):
    """Dict {id: pdb_code} — keys become IDs. Constructed outside a Pipeline to
    inspect raw attributes (pipeline context returns a StandardizedOutput)."""
    from biopipelines.pdb import PDB

    p = PDB(pdbs={"POI1": "4ufc", "POI2": "1ake"})
    record_case(input="PDB(pdbs={'POI1':'4ufc','POI2':'1ake'})",
                expected=(["POI1", "POI2"], ["4ufc", "1ake"]),
                actual=(list(p.custom_ids), list(p.pdb_ids)))
    assert p.custom_ids == ["POI1", "POI2"]
    assert p.pdb_ids == ["4ufc", "1ake"]


def test_pdb_custom_ids_list(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """PDB(pdbs=[...], ids=[...]) renames independently."""
    from biopipelines.pdb import PDB

    pipeline = new_pipeline("pdb_custom_ids")
    with pipeline:
        p = PDB(pdbs=["4ufc", "1ake"], ids=["target1", "target2"])
        pipeline.save()

    ids = list(p.streams.structures.ids)
    record_case(input="PDB(pdbs=[...], ids=['target1','target2'])",
                expected=["target1", "target2"], actual=ids)
    assert ids == ["target1", "target2"]


# ── folder loading ───────────────────────────────────────────────────────────

def test_pdb_loads_from_folder(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """PDB(pdbs=<folder>) scans the folder and uses basenames as IDs."""
    from biopipelines.pdb import PDB

    folder = isolated_cwd / "my_pdbs"
    folder.mkdir()
    (folder / "struct1.pdb").write_text("HEADER test1\nEND\n")
    (folder / "struct2.pdb").write_text("HEADER test2\nEND\n")

    pipeline = new_pipeline("pdb_folder")
    with pipeline:
        p = PDB(pdbs=str(folder))
        pipeline.save()

    ids = sorted(list(p.streams.structures.ids))
    record_case(input="PDB(pdbs=<folder with 2 .pdbs>)",
                expected=["struct1", "struct2"], actual=ids)
    assert ids == ["struct1", "struct2"]


# ── convert option ───────────────────────────────────────────────────────────

def test_pdb_convert_pdb(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    from biopipelines.pdb import PDB

    pipeline = new_pipeline("pdb_convert_pdb")
    with pipeline:
        p = PDB(pdbs="4ufc", convert="pdb")
        pipeline.save()

    record_case(input="PDB(convert='pdb')",
                expected="pdb", actual=p.streams.structures.format)
    assert p.streams.structures.format == "pdb"


def test_pdb_convert_cif(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    from biopipelines.pdb import PDB

    pipeline = new_pipeline("pdb_convert_cif")
    with pipeline:
        p = PDB(pdbs="4ufc", convert="cif")
        pipeline.save()

    record_case(input="PDB(convert='cif')",
                expected="cif", actual=p.streams.structures.format)
    assert p.streams.structures.format == "cif"


def test_pdb_convert_none_is_mixed_format(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """convert=None leaves format as 'pdb|cif' when a download is needed."""
    from biopipelines.pdb import PDB

    pipeline = new_pipeline("pdb_convert_none")
    with pipeline:
        p = PDB(pdbs="4ufc")  # convert defaults to None
        pipeline.save()

    record_case(input="PDB(convert=None)",
                expected="pdb|cif", actual=p.streams.structures.format)
    assert p.streams.structures.format == "pdb|cif"


def test_pdb_invalid_convert_raises(record_case):
    from biopipelines.pdb import PDB

    record_case(input="PDB(convert='xyz')",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError, match="convert"):
        PDB(pdbs="4ufc", convert="xyz")


# ── operations (PDB.rename) ──────────────────────────────────────────────────

def test_pdb_rename_operation(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """PDB.rename creates a PDBOperation that rides along in to_dict."""
    from biopipelines.pdb import PDB, PDBOperation

    op = PDB.rename("LIG", ":L:")
    record_case(input="PDB.rename('LIG', ':L:')",
                expected=("rename", "LIG", ":L:"),
                actual=(op.op_type, op.params["old"], op.params["new"]))
    assert isinstance(op, PDBOperation)
    assert op.op_type == "rename"
    assert op.params == {"old": "LIG", "new": ":L:"}

    # Outside a pipeline to keep raw-attribute access to `operations`.
    p = PDB("4ufc", PDB.rename("LIG", ":L:"))
    assert len(p.operations) == 1
    assert p.operations[0].op_type == "rename"


def test_pdb_rotate_bond_operation(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """PDB.rotate_bond creates a PDBOperation carrying atoms + angle."""
    from biopipelines.pdb import PDB, PDBOperation

    op = PDB.rotate_bond("LIG.C60", "LIG.C70", 180)
    record_case(input="PDB.rotate_bond('LIG.C60', 'LIG.C70', 180)",
                expected=("rotate_bond", "LIG.C60", "LIG.C70", 180.0),
                actual=(op.op_type, op.params["atom1"], op.params["atom2"],
                        op.params["angle"]))
    assert isinstance(op, PDBOperation)
    assert op.op_type == "rotate_bond"
    assert op.params == {"atom1": "LIG.C60", "atom2": "LIG.C70", "angle": 180.0}
    # angle is coerced to float for clean JSON serialization.
    assert isinstance(op.params["angle"], float)


def test_pdb_rejects_non_operation_positional_args(record_case):
    from biopipelines.pdb import PDB

    record_case(input="PDB('4ufc', 'not_an_op')",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError, match="PDBOperation"):
        PDB("4ufc", "not_an_op")


# ── validation errors ───────────────────────────────────────────────────────

def test_pdb_ids_length_mismatch_raises(record_case):
    from biopipelines.pdb import PDB

    record_case(input="len(ids)=1 vs len(pdbs)=2",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError, match="Length mismatch"):
        PDB(pdbs=["4ufc", "1ake"], ids=["only_one"])


def test_pdb_unsupported_type_raises(record_case):
    from biopipelines.pdb import PDB

    record_case(input="PDB(pdbs=42) — int",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError, match="string, list"):
        PDB(pdbs=42)


def test_pdb_empty_dict_raises(record_case):
    from biopipelines.pdb import PDB

    record_case(input="PDB(pdbs={})",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError, match="at least one"):
        PDB(pdbs={})


# ── integration: PDB + Mock ──────────────────────────────────────────────────

def test_pdb_feeds_mock(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Mock consumes PDB.streams.structures and fans out with children=<1..3>."""
    from biopipelines.pdb import PDB
    from biopipelines.mock import Mock

    pipeline = new_pipeline("pdb_feeds_mock")
    with pipeline:
        p = PDB(pdbs={"A": "4ufc", "B": "1ake"})
        m = Mock(
            source=p.streams.structures,
            children="<1..3>",
            streams={"designs": {"format": "pdb", "file": "<id>.pdb"}},
        )
        script_path = pipeline.save()

    ids = list(m.streams.designs.ids)
    record_case(input="PDB → Mock(children=<1..3>)",
                expected=["A_<1..3>", "B_<1..3>"], actual=ids)
    assert ids == ["A_<1..3>", "B_<1..3>"]
    assert_valid_script(script_path, "PDB", "Mock")


# ── integration: PDB + Panda ─────────────────────────────────────────────────

def test_pdb_feeds_panda_on_sequences_table(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda sorts PDB's sequences table."""
    from biopipelines.pdb import PDB
    from biopipelines.panda import Panda

    pipeline = new_pipeline("pdb_feeds_panda")
    with pipeline:
        p = PDB(pdbs=["4ufc", "1ake"])
        Panda(
            tables=p.tables.sequences,
            operations=[Panda.sort("id", ascending=True)],
        )
        script_path = pipeline.save()

    record_case(input="PDB → Panda(sort on sequences table)",
                expected="script with PDB+Panda",
                actual=os.path.basename(script_path))
    assert_valid_script(script_path, "PDB", "Panda", "pdb_feeds_panda")


def test_pdb_feeds_panda_on_structures_table(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda can also sort the structures table (file_size, source, etc.)."""
    from biopipelines.pdb import PDB
    from biopipelines.panda import Panda

    pipeline = new_pipeline("pdb_feeds_panda_struct")
    with pipeline:
        p = PDB(pdbs=["4ufc", "1ake"])
        Panda(
            tables=p.tables.structures,
            operations=[
                Panda.sort("file_size", ascending=False),
                Panda.head(1),
            ],
        )
        script_path = pipeline.save()

    record_case(input="PDB → Panda(sort structures by file_size + head)",
                expected="script",
                actual=os.path.basename(script_path))
    assert_valid_script(script_path, "PDB", "Panda")


def test_pdb_panda_filter(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda.filter on PDB.tables.structures — expression over 'source' column."""
    from biopipelines.pdb import PDB
    from biopipelines.panda import Panda

    pipeline = new_pipeline("pdb_panda_filter")
    with pipeline:
        p = PDB(pdbs=["4ufc", "1ake"])
        Panda(
            tables=p.tables.structures,
            operations=[Panda.filter("source == 'rcsb_download'")],
        )
        script_path = pipeline.save()

    record_case(input="PDB → Panda(filter by source)",
                expected="script", actual=os.path.basename(script_path))
    assert_valid_script(script_path, "PDB", "Panda")


def test_pdb_panda_tail(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda.tail on PDB's sequences table."""
    from biopipelines.pdb import PDB
    from biopipelines.panda import Panda

    pipeline = new_pipeline("pdb_panda_tail")
    with pipeline:
        p = PDB(pdbs=["4ufc", "1ake"])
        Panda(
            tables=p.tables.sequences,
            operations=[Panda.tail(1)],
        )
        script_path = pipeline.save()

    record_case(input="PDB → Panda(tail(1))",
                expected="script", actual=os.path.basename(script_path))
    assert_valid_script(script_path, "PDB", "Panda")


def test_pdb_panda_filter_sort_head_chain(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Full filter+sort+head chain over PDB.tables.structures."""
    from biopipelines.pdb import PDB
    from biopipelines.panda import Panda

    pipeline = new_pipeline("pdb_panda_fsh")
    with pipeline:
        p = PDB(pdbs=["4ufc", "1ake"])
        Panda(
            tables=p.tables.structures,
            operations=[
                Panda.filter("file_size > 0"),
                Panda.sort("file_size", ascending=False),
                Panda.head(1),
            ],
        )
        script_path = pipeline.save()

    record_case(input="PDB → Panda(filter+sort+head)",
                expected="chained script",
                actual=os.path.basename(script_path))
    assert_valid_script(script_path, "PDB", "Panda")


def test_pdb_panda_select_columns(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """select_columns keeps only id + format from the structures table."""
    from biopipelines.pdb import PDB
    from biopipelines.panda import Panda

    pipeline = new_pipeline("pdb_panda_select")
    with pipeline:
        p = PDB(pdbs=["4ufc"])
        Panda(
            tables=p.tables.structures,
            operations=[Panda.select_columns(["id", "format"])],
        )
        script_path = pipeline.save()

    record_case(input="PDB → Panda(select_columns id,format)",
                expected="script", actual=os.path.basename(script_path))
    assert_valid_script(script_path, "PDB", "Panda")


def test_pdb_panda_rename_column(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """rename remaps pdb_id → rcsb_id."""
    from biopipelines.pdb import PDB
    from biopipelines.panda import Panda

    pipeline = new_pipeline("pdb_panda_rename")
    with pipeline:
        p = PDB(pdbs=["4ufc"])
        Panda(
            tables=p.tables.structures,
            operations=[Panda.rename({"pdb_id": "rcsb_id"})],
        )
        script_path = pipeline.save()

    record_case(input="PDB → Panda(rename pdb_id→rcsb_id)",
                expected="script", actual=os.path.basename(script_path))
    assert_valid_script(script_path, "PDB", "Panda")


# ── chain semantics: auto / all / list / explicit letter ─────────────────────

def test_pdb_chain_auto_emits_literal_sequence_ids(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """chain='auto' (default) yields one literal id per input on both streams."""
    from biopipelines.pdb import PDB

    pipeline = new_pipeline("pdb_chain_auto")
    with pipeline:
        p = PDB(pdbs=["4ufc", "1ake"])  # default chain='auto'
        pipeline.save()

    seq_ids = list(p.streams.sequences.ids)
    struct_ids = list(p.streams.structures.ids)
    record_case(input="PDB(['4ufc','1ake']) chain='auto'",
                expected={"sequences": ["4ufc", "1ake"], "structures": ["4ufc", "1ake"]},
                actual={"sequences": seq_ids, "structures": struct_ids})
    assert seq_ids == ["4ufc", "1ake"]
    assert struct_ids == ["4ufc", "1ake"]


def test_pdb_chain_all_emits_lazy_sequence_ids(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """chain='all' declares lazy <id>[_<chain>] sequences (cardinality at runtime)."""
    from biopipelines.pdb import PDB

    pipeline = new_pipeline("pdb_chain_all")
    with pipeline:
        p = PDB(pdbs=["1a3n"], chain="all")
        pipeline.save()

    seq_ids = list(p.streams.sequences.ids)
    struct_ids = list(p.streams.structures.ids)
    record_case(input="PDB('1a3n', chain='all')",
                expected={"sequences": ["1a3n[_<chain>]"], "structures": ["1a3n"]},
                actual={"sequences": seq_ids, "structures": struct_ids})
    assert seq_ids == ["1a3n[_<chain>]"]
    # split_chains defaults to False → structures stream stays single-id literal
    assert struct_ids == ["1a3n"]


def test_pdb_chain_list_emits_literal_per_chain_ids(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """chain=['A','C'] yields literal <id>_A, <id>_C sequences (cardinality known)."""
    from biopipelines.pdb import PDB

    pipeline = new_pipeline("pdb_chain_list")
    with pipeline:
        p = PDB(pdbs=["1a3n"], chain=["A", "C"])
        pipeline.save()

    seq_ids = list(p.streams.sequences.ids)
    struct_ids = list(p.streams.structures.ids)
    record_case(input="PDB('1a3n', chain=['A','C'])",
                expected={"sequences": ["1a3n_A", "1a3n_C"], "structures": ["1a3n"]},
                actual={"sequences": seq_ids, "structures": struct_ids})
    assert seq_ids == ["1a3n_A", "1a3n_C"]
    assert struct_ids == ["1a3n"]  # split_chains=False → single structure file


def test_pdb_chain_list_split_chains_emits_per_chain_structures(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """chain=['A','B'] + split_chains=True splits structures stream too."""
    from biopipelines.pdb import PDB

    pipeline = new_pipeline("pdb_chain_list_split")
    with pipeline:
        p = PDB(pdbs=["4ake"], chain=["A", "B"], split_chains=True)
        pipeline.save()

    seq_ids = list(p.streams.sequences.ids)
    struct_ids = list(p.streams.structures.ids)
    record_case(input="PDB('4ake', chain=['A','B'], split_chains=True)",
                expected={"sequences": ["4ake_A", "4ake_B"], "structures": ["4ake_A", "4ake_B"]},
                actual={"sequences": seq_ids, "structures": struct_ids})
    assert seq_ids == ["4ake_A", "4ake_B"]
    assert struct_ids == ["4ake_A", "4ake_B"]


def test_pdb_chain_all_split_chains_emits_lazy_structures(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """chain='all' + split_chains=True declares both streams lazy."""
    from biopipelines.pdb import PDB

    pipeline = new_pipeline("pdb_chain_all_split")
    with pipeline:
        p = PDB(pdbs=["4lcd"], chain="all", split_chains=True)
        pipeline.save()

    seq_ids = list(p.streams.sequences.ids)
    struct_ids = list(p.streams.structures.ids)
    record_case(input="PDB('4lcd', chain='all', split_chains=True)",
                expected={"sequences": ["4lcd[_<chain>]"], "structures": ["4lcd[_<chain>]"]},
                actual={"sequences": seq_ids, "structures": struct_ids})
    assert seq_ids == ["4lcd[_<chain>]"]
    assert struct_ids == ["4lcd[_<chain>]"]


def test_pdb_chain_explicit_letter_emits_single_literal_id(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """chain='A' yields a single literal id and leaves the structure file filtered."""
    from biopipelines.pdb import PDB

    pipeline = new_pipeline("pdb_chain_letter")
    with pipeline:
        p = PDB(pdbs=["4lcd"], chain="E")
        pipeline.save()

    seq_ids = list(p.streams.sequences.ids)
    struct_ids = list(p.streams.structures.ids)
    record_case(input="PDB('4lcd', chain='E')",
                expected={"sequences": ["4lcd"], "structures": ["4lcd"]},
                actual={"sequences": seq_ids, "structures": struct_ids})
    assert seq_ids == ["4lcd"]
    assert struct_ids == ["4lcd"]


def test_pdb_split_chains_requires_multi_chain_selector(record_case):
    """split_chains=True is meaningless with single-chain selectors."""
    from biopipelines.pdb import PDB

    record_case(input="PDB('4ufc', chain='auto', split_chains=True)",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError, match="split_chains"):
        PDB(pdbs="4ufc", chain="auto", split_chains=True)
    with pytest.raises(ValueError, match="split_chains"):
        PDB(pdbs="4ufc", chain="A", split_chains=True)


def test_pdb_chain_empty_list_raises(record_case):
    """chain=[] is not a valid selector."""
    from biopipelines.pdb import PDB

    record_case(input="PDB(pdbs='4ufc', chain=[])",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError, match="empty"):
        PDB(pdbs="4ufc", chain=[])


def test_pdb_chain_list_with_non_string_raises(record_case):
    """chain list elements must be non-empty strings."""
    from biopipelines.pdb import PDB

    record_case(input="PDB(pdbs='4ufc', chain=['A', 1])",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError, match="non-empty"):
        PDB(pdbs="4ufc", chain=["A", 1])
    with pytest.raises(ValueError, match="non-empty"):
        PDB(pdbs="4ufc", chain=["A", ""])


def test_pdb_chain_list_split_propagates_to_dict(record_case):
    """to_dict round-trips chain list and split_chains."""
    from biopipelines.pdb import PDB

    p = PDB(pdbs="4ake", chain=["A", "B"], split_chains=True)
    d = p.to_dict()
    record_case(input="PDB('4ake', chain=['A','B'], split_chains=True).to_dict()",
                expected=(["A", "B"], True),
                actual=(d["tool_params"]["chain"], d["tool_params"]["split_chains"]))
    assert d["tool_params"]["chain"] == ["A", "B"]
    assert d["tool_params"]["split_chains"] is True


# ── remove-operation selection parsing ───────────────────────────────────────

def _load_pipe_pdb():
    import importlib.util
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    path = os.path.join(repo, "pipe_scripts", "pipe_pdb.py")
    spec = importlib.util.spec_from_file_location("pipe_pdb", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_parse_residue_selection_ranges_and_names(record_case):
    """Range-shaped tokens parse as ranges; the rest as residue names; they mix."""
    pp = _load_pipe_pdb()
    ranges, names = pp._parse_residue_selection("NAP+A1-83+EDO+90-95")
    record_case(input="_parse_residue_selection('NAP+A1-83+EDO+90-95')",
                expected=([("A", 1, 83), (None, 90, 95)], {"NAP", "EDO"}),
                actual=(ranges, names))
    assert ranges == [("A", 1, 83), (None, 90, 95)]
    assert names == {"NAP", "EDO"}


def test_parse_residue_selection_digit_bearing_names(record_case):
    """A CCD code with digits ('9DP', 'A7ZK') is a residue name, not a range."""
    pp = _load_pipe_pdb()
    ranges, names = pp._parse_residue_selection("9DP+A7ZK+A1-83")
    record_case(input="_parse_residue_selection('9DP+A7ZK+A1-83')",
                expected=([("A", 1, 83)], {"9DP", "A7ZK"}),
                actual=(ranges, names))
    assert ranges == [("A", 1, 83)]
    assert names == {"9DP", "A7ZK"}


def test_apply_remove_operation_digit_bearing_name(record_case):
    """PDB.remove('9DP') drops the digit-bearing CCD ligand by name."""
    pp = _load_pipe_pdb()
    pdb = (
        "ATOM      1  CA  ALA A   1      11.0  11.0  11.0  1.00  0.00           C\n"
        "HETATM    2  PA  9DP A 500      20.0  20.0  20.0  1.00  0.00           P\n"
        "HETATM    3  C1  STI A 700      40.0  40.0  40.0  1.00  0.00           C\n"
    )
    out = pp.apply_remove_operation(pdb, "pdb", "9DP")
    record_case(input="apply_remove_operation(pdb,'9DP')",
                expected="9DP gone; ALA, STI kept", actual=out)
    assert "9DP" not in out
    assert "ALA" in out and "STI" in out


def test_parse_residue_selection_b12_collision(record_case):
    """'B12'-shaped token: a present residue name wins over the chain reading;
    absent, it falls back to a chain-qualified range (B, 12, 12)."""
    pp = _load_pipe_pdb()
    as_name = pp._parse_residue_selection("B12", present_names={"B12"})
    as_range = pp._parse_residue_selection("B12", present_names={"ALA"})
    record_case(input="_parse_residue_selection('B12', present/absent)",
                expected="present -> name {B12}; absent -> range (B,12,12)",
                actual=(as_name, as_range))
    assert as_name == ([], {"B12"})
    assert as_range == ([("B", 12, 12)], set())


def test_apply_remove_operation_b12_ligand_vs_chain(record_case):
    """End to end: when B12 is a bound ligand it is removed as a name; the same
    token against a structure with no B12 ligand removes chain-B residue 12."""
    pp = _load_pipe_pdb()
    # B12 present as a HETATM ligand -> name removal (the chain-A copy stays).
    pdb_lig = (
        "ATOM      1  CA  ALA A  12      11.0  11.0  11.0  1.00  0.00           C\n"
        "HETATM    2  CO  B12 B 500      20.0  20.0  20.0  1.00  0.00          CO\n"
    )
    out_lig = pp.apply_remove_operation(pdb_lig, "pdb", "B12")
    assert "B12" not in out_lig and "ALA" in out_lig

    # No B12 ligand -> 'B12' reads as chain B residue 12 (the chain-A res 12 stays).
    pdb_nolig = (
        "ATOM      1  CA  ALA A  12      11.0  11.0  11.0  1.00  0.00           C\n"
        "ATOM      2  CA  GLY B  12      20.0  20.0  20.0  1.00  0.00           C\n"
    )
    out_nolig = pp.apply_remove_operation(pdb_nolig, "pdb", "B12")
    record_case(input="apply_remove_operation B12 ligand-vs-chain",
                expected="ligand: B12 gone/ALA kept; no-ligand: chain-B res12 gone/chain-A kept",
                actual=(out_lig, out_nolig))
    kept = [l for l in out_nolig.split("\n") if l.startswith("ATOM")]
    assert len(kept) == 1 and kept[0][21] == "A"


def test_parse_residue_selection_pure_name(record_case):
    """A bare residue name (the 3QWI case: PDB.remove('NAP')) parses, no crash."""
    pp = _load_pipe_pdb()
    ranges, names = pp._parse_residue_selection("NAP")
    record_case(input="_parse_residue_selection('NAP')",
                expected=([], {"NAP"}), actual=(ranges, names))
    assert ranges == []
    assert names == {"NAP"}


def test_parse_residue_selection_unparseable_token_raises(record_case):
    """A digit-bearing but malformed token gives a clear error (not int(''))."""
    pp = _load_pipe_pdb()
    with pytest.raises(ValueError, match="unparseable token"):
        pp._parse_residue_selection("1-x")
    record_case(input="_parse_residue_selection('1-x')",
                expected="ValueError(unparseable token)",
                actual="ValueError(unparseable token)")


def test_apply_remove_operation_by_name(record_case):
    """PDB.remove('NAP') drops the named HETATM group regardless of remove_hetatm,
    while leaving protein and other ligands intact."""
    pp = _load_pipe_pdb()
    pdb = (
        "ATOM      1  CA  ALA A   1      11.0  11.0  11.0  1.00  0.00           C\n"
        "HETATM    2  PA  NAP A 500      20.0  20.0  20.0  1.00  0.00           P\n"
        "HETATM    3  C1  EDO A 600      30.0  30.0  30.0  1.00  0.00           C\n"
        "HETATM    4  C1  STI A 700      40.0  40.0  40.0  1.00  0.00           C\n"
    )
    # remove_hetatm=False would normally keep all HETATM, but a NAMED removal
    # still drops NAP.
    out = pp.apply_remove_operation(pdb, "pdb", "NAP", remove_hetatm=False)
    record_case(input="apply_remove_operation(pdb,'NAP',remove_hetatm=False)",
                expected="NAP gone; ALA, EDO, STI kept",
                actual=out)
    assert "NAP" not in out
    assert "ALA" in out and "EDO" in out and "STI" in out

    # Mixed name+name removal ('NAP+EDO') drops both.
    out2 = pp.apply_remove_operation(pdb, "pdb", "NAP+EDO")
    assert "NAP" not in out2 and "EDO" not in out2
    assert "STI" in out2 and "ALA" in out2


def test_apply_remove_operation_drops_dependent_records(record_case):
    """Removing a residue also drops ANISOU/CONECT/LINK records that reference its
    atoms: a CONECT keeping >=1 partner is rewritten; one left with only its base
    serial (no bonded partner) is dropped, not emitted as a lone-base record."""
    pp = _load_pipe_pdb()
    pdb = (
        "ATOM      1  CA  ALA A   1      11.0  11.0  11.0  1.00  0.00           C\n"
        "ANISOU    1  CA  ALA A   1      100  100  100    0    0    0       C\n"
        "HETATM    2  PA  NAP A 500      20.0  20.0  20.0  1.00  0.00           P\n"
        "ANISOU    2  PA  NAP A 500      200  200  200    0    0    0       P\n"
        "HETATM    3  C1  STI A 700      40.0  40.0  40.0  1.00  0.00           C\n"
        "HETATM    4  C2  STI A 700      41.0  40.0  40.0  1.00  0.00           C\n"
        "LINK         PA  NAP A 500                 C1  STI A 700     1555 1555  1.5\n"
        "CONECT    2    3\n"      # base removed -> dropped
        "CONECT    3    2    4\n"  # partner 2 removed, 4 survives -> rewritten
        "CONECT    4    2\n"      # only partner removed, lone base left -> dropped
    )
    out = pp.apply_remove_operation(pdb, "pdb", "NAP")
    lines = out.split("\n")
    record_case(input="apply_remove_operation(pdb,'NAP') with deps",
                expected="NAP atom+ANISOU+LINK gone; CONECT 2 & 4 dropped; CONECT 3 -> '3 4'",
                actual=out)
    # NAP's atom and its ANISOU are gone; serial-2 ANISOU gone.
    assert not any(l.startswith(("ATOM", "HETATM", "ANISOU")) and pp.field_res_name(l).upper() == "NAP"
                   for l in lines if len(l) > 20)
    assert not any(l.startswith("ANISOU") and l[6:11].strip() == "2" for l in lines)
    # LINK referencing NAP is dropped entirely.
    assert not any(l.startswith("LINK") for l in lines)
    # Only CONECT 3 survives, rewritten to drop serial 2 and keep serial 4.
    conects = [l for l in lines if l.startswith("CONECT")]
    assert len(conects) == 1
    refs = [conects[0][i:i + 5].strip() for i in range(6, len(conects[0]), 5)]
    refs = [r for r in refs if r]
    assert refs == ["3", "4"]  # base + one surviving partner, no lone base, no dead serial 2
    # The surviving STI atoms are untouched.
    assert sum(1 for l in lines if l.startswith("HETATM") and pp.field_res_name(l) == "STI") == 2
