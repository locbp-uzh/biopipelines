"""SelectionTransfer carries a residue selection across trims, renumbering, chain renames and circular permutations.

A wrong map is silent downstream: ProteinMPNN redesigns whatever residues it is handed. These tests pin the mapping on synthetic structures whose correspondence is known by construction.
"""

import importlib.util
import os
import subprocess
import sys

import pandas as pd
import pytest

pytest.importorskip("Bio.Align")

from biopipelines.datastream import DataStream  # noqa: E402

REPO = os.path.join(os.path.dirname(__file__), "..")
SCRIPT = os.path.join(REPO, "pipe_scripts", "pipe_selection_transfer.py")

REF = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQ"

ONE_TO_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU",
    "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE",
    "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}


@pytest.fixture(scope="module")
def st():
    spec = importlib.util.spec_from_file_location("pipe_selection_transfer", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def chains(seq, chain="A", start=1):
    return {chain: [(start + i, aa) for i, aa in enumerate(seq)]}


def write_pdb(path, seq, chain="A", start=1):
    lines = []
    for i, aa in enumerate(seq):
        lines.append(f"ATOM  {i + 1:5d}  CA  {ONE_TO_THREE[aa]} {chain}{start + i:4d}    "
                     f"{3.8 * i:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\nEND\n")


def run(st, residues, ref, tgt):
    return st.transfer(residues, ref, tgt, st.make_aligner())


def test_a_trimmed_renumbered_renamed_target_maps_by_sequence(st):
    mapped, unmapped, identity = run(st, [("A", 10), ("A", 11), ("A", 12)],
                                     chains(REF), chains(REF[5:60], chain="B", start=1))
    assert mapped == [("B", 5), ("B", 6), ("B", 7)]
    assert unmapped == []
    assert identity == 1.0


def test_a_circular_permutant_with_a_linker_maps_both_halves(st):
    permuted = REF[30:] + "GSGG" + REF[:30]
    mapped, unmapped, _ = run(st, [("A", 5), ("A", 40)], chains(REF), chains(permuted))
    tail = len(REF) - 30
    assert mapped == sorted([("A", tail + 4 + 5), ("A", 40 - 30)])
    assert unmapped == []


def test_point_mutations_still_map_and_lower_identity(st):
    mutant = REF[:20] + "W" + REF[21:50] + "A" + REF[51:]
    mapped, _, identity = run(st, [("A", 21), ("A", 51)], chains(REF), chains(mutant))
    assert mapped == [("A", 21), ("A", 51)]
    assert identity < 1.0


def test_a_deleted_residue_is_reported_not_guessed(st):
    deleted = REF[:30] + REF[35:]
    mapped, unmapped, _ = run(st, [("A", 29), ("A", 33), ("A", 40)], chains(REF), chains(deleted))
    assert unmapped == [("A", 33)]
    assert mapped == [("A", 29), ("A", 35)]


def test_a_chainless_selection_needs_a_single_chain_reference(st):
    mapped, _, _ = run(st, [("", 3)], chains(REF), chains(REF))
    assert mapped == [("A", 3)]
    two = {**chains(REF), **chains(REF, chain="B")}
    with pytest.raises(ValueError, match="chainless"):
        run(st, [("", 3)], two, chains(REF))


def test_a_residue_absent_from_the_reference_is_an_error(st):
    with pytest.raises(ValueError, match="not in the reference"):
        run(st, [("A", 500)], chains(REF), chains(REF))


def test_the_pipe_pairs_targets_with_their_reference_and_reads_the_table(st, tmp_path):
    folder = tmp_path / "003_SelectionTransfer"
    for sub in ("_configuration", "_execution", "tables", "inputs"):
        (folder / sub).mkdir(parents=True)
    inputs = folder / "inputs"
    write_pdb(inputs / "d1.pdb", REF)
    write_pdb(inputs / "d2.pdb", REF[10:])
    write_pdb(inputs / "d1_1.pdb", REF[3:], chain="B")
    write_pdb(inputs / "d2_1.pdb", REF[20:] + "GSG" + REF[:20])

    pd.DataFrame({"id": ["d1", "d2"], "designed": ["A8-9", "A12"]}).to_csv(inputs / "rfd.csv", index=False)
    DataStream(name="structures", ids=["d1", "d2"], files=[str(inputs / "<id>.pdb")],
               format="pdb").save_json(str(folder / "_configuration" / "reference.json"))
    DataStream(name="structures", ids=["d1_1", "d2_1"], files=[str(inputs / "<id>.pdb")],
               format="pdb").save_json(str(folder / "_configuration" / "structures.json"))
    config = folder / "_configuration" / "config.json"
    config.write_text(f'{{"selection": {{"table": "{(inputs / "rfd.csv").as_posix()}", "column": "designed"}}}}')

    out = folder / "tables" / "selections.csv"
    subprocess.run([sys.executable, SCRIPT,
                    "--config-json", str(config),
                    "--structures-json", str(folder / "_configuration" / "structures.json"),
                    "--reference-json", str(folder / "_configuration" / "reference.json"),
                    "--selections-csv", str(out),
                    "--local-missing-csv", str(folder / "_execution" / "local_missing.csv")],
                   check=True, cwd=REPO)

    df = pd.read_csv(out).set_index("id")
    assert df.loc["d1_1", "selection"] == "B5-6"
    assert df.loc["d1_1", "reference.id"] == "d1"
    # d2 starts at REF[10] numbered 1, so A12 is REF[21]; in d2_1 that is position 2 of REF[20:].
    assert df.loc["d2_1", "selection"] == "A2"
    assert df.loc["d2_1", "reference.id"] == "d2"


# ── pairing a target with its reference ──────────────────────────────────────

def test_a_sibling_reference_is_not_a_match(st):
    """`des_3_1` shares only the ancestor `des` with `des_1`, whose selection sits on another design."""
    pairing = st.pair_references(["des_3_1"], ["des_1", "des_2"], 2)
    assert pairing == {"des_3_1": (None, "des_1")}


def test_the_exact_parent_and_child_tiers_still_match(st):
    references = ["des_1", "des_2"]
    assert st.pair_references(["des_2"], references, 2)["des_2"] == ("des_2", None)
    assert st.pair_references(["des_2_1"], references, 2)["des_2_1"] == ("des_2", None)
    assert st.pair_references(["des_2"], ["des_1_1", "des_2_1"], 2)["des_2"] == ("des_2_1", None)


def test_the_provenance_of_the_target_map_names_the_reference(st, tmp_path):
    table = tmp_path / "structures_map.csv"
    table.write_text("id,file,structures.id\nPanda_1,a.pdb,des_2\n")
    assert st.pair_references(["Panda_1"], ["des_1", "des_2"], 2, [str(table)])["Panda_1"] == ("des_2", None)


def test_one_declared_reference_is_the_reference_of_every_target(st):
    assert st.pair_references(["a_1", "b_1"], ["ref"], 1) == {"a_1": ("ref", None), "b_1": ("ref", None)}


def test_a_reference_filtered_upstream_is_not_broadcast_to_the_other_targets(st):
    """The stream declares two ids but only `des_1` has a file; counting the files would hand its selection to `des_2_1`."""
    pairing = st.pair_references(["des_1_1", "des_2_1"], ["des_1"], 2)
    assert pairing["des_1_1"] == ("des_1", None)
    assert pairing["des_2_1"][0] is None


def _stage(tmp_path, references, targets, absent=()):
    """A step folder with a reference and a target stream; `references` and `targets` map id to sequence, `absent` ids have no file."""
    folder = tmp_path / "003_SelectionTransfer"
    for sub in ("_configuration", "_execution", "tables", "inputs"):
        (folder / sub).mkdir(parents=True)
    inputs = folder / "inputs"
    for sid, seq in {**references, **targets}.items():
        if sid not in absent:
            write_pdb(inputs / f"{sid}.pdb", seq)
    DataStream(name="structures", ids=list(references), files=[str(inputs / "<id>.pdb")],
               format="pdb").save_json(str(folder / "_configuration" / "reference.json"))
    DataStream(name="structures", ids=list(targets), files=[str(inputs / "<id>.pdb")],
               format="pdb").save_json(str(folder / "_configuration" / "structures.json"))
    config = folder / "_configuration" / "config.json"
    config.write_text('{"selection": {"literal": "A8-9"}}')
    return folder, config


def _run_pipe(folder, config, check=True):
    out = folder / "tables" / "selections.csv"
    result = subprocess.run([sys.executable, SCRIPT,
                             "--config-json", str(config),
                             "--structures-json", str(folder / "_configuration" / "structures.json"),
                             "--reference-json", str(folder / "_configuration" / "reference.json"),
                             "--selections-csv", str(out),
                             "--local-missing-csv", str(folder / "_execution" / "local_missing.csv")],
                            cwd=REPO, capture_output=True, text=True)
    if check:
        assert result.returncode == 0, result.stderr
    return result, out


def test_a_target_whose_own_reference_is_absent_is_unmatched_not_given_a_siblings(tmp_path):
    folder, config = _stage(tmp_path,
                            references={"des_1": REF, "des_2": REF[10:]},
                            targets={"des_1_1": REF[3:], "des_3_1": REF[3:]})
    _, out = _run_pipe(folder, config)

    df = pd.read_csv(out).set_index("id")
    assert df.loc["des_1_1", "reference.id"] == "des_1" and df.loc["des_1_1", "selection"] == "A5-6"
    assert pd.isna(df.loc["des_3_1", "reference.id"]) and pd.isna(df.loc["des_3_1", "selection"])
    missing = pd.read_csv(folder / "_execution" / "local_missing.csv")
    assert list(missing["id"]) == ["des_3_1"] and list(missing["kind"]) == ["unmatched"]
    assert "'des_1' is only a sibling" in missing.loc[0, "cause"]
    assert list(missing["removed_by"]) == ["003_SelectionTransfer"]


def test_a_reference_with_a_missing_file_is_not_broadcast_by_the_pipe(tmp_path):
    folder, config = _stage(tmp_path,
                            references={"des_1": REF, "des_2": REF},
                            targets={"des_1_1": REF[3:], "des_2_1": REF[3:]}, absent={"des_2"})
    _, out = _run_pipe(folder, config)
    df = pd.read_csv(out).set_index("id")
    assert df.loc["des_1_1", "reference.id"] == "des_1" and pd.isna(df.loc["des_2_1", "reference.id"])
    assert list(pd.read_csv(folder / "_execution" / "local_missing.csv")["id"]) == ["des_2_1"]


def test_every_target_unmatched_is_a_step_failure_with_the_tables_written(tmp_path):
    folder, config = _stage(tmp_path, references={"des_1": REF, "des_2": REF}, targets={"des_3_1": REF})
    result, out = _run_pipe(folder, config, check=False)
    assert result.returncode == 1 and "all 1 id(s) failed" in result.stderr
    assert out.exists() and len(pd.read_csv(folder / "_execution" / "local_missing.csv")) == 1


def test_a_target_carrying_the_reference_twice_fails_and_says_why(st):
    """A tandem fusion aligns each reference residue twice across the doubled reference; the id fails instead of mapping one copy at random."""
    with pytest.raises(ValueError, match="aligned twice"):
        run(st, [("A", 10)], chains(REF), chains(REF + REF))


def test_a_tandem_target_is_a_failure_row_of_the_step(tmp_path):
    folder, config = _stage(tmp_path, references={"des_1": REF}, targets={"des_1_1": REF[3:], "des_1_2": REF + REF})
    _, out = _run_pipe(folder, config)
    missing = pd.read_csv(folder / "_execution" / "local_missing.csv")
    assert list(missing["id"]) == ["des_1_2"] and list(missing["kind"]) == ["failure"]
    assert "aligned twice" in missing.loc[0, "cause"]


def _saved_config(root) -> dict:
    import glob
    import json
    matches = glob.glob(os.path.join(str(root), "**", "*_SelectionTransfer", "_configuration",
                                     "selection_transfer_config.json"), recursive=True)
    assert len(matches) == 1, matches
    with open(matches[0]) as f:
        return json.load(f)["selection"]


def test_a_column_reference_is_serialized_for_the_pipe(local_config, isolated_cwd, new_pipeline):
    from biopipelines.mock import Mock
    from biopipelines.selection_transfer import SelectionTransfer

    pipeline = new_pipeline("st_ref")
    with pipeline:
        ref = Mock(ids=["d1"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
                   tables={"structures": {"columns": ["designed"], "fill": {"designed": "A1-5"}}})
        tgt = Mock(ids=["d1_1"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
        st = SelectionTransfer(selection=ref.tables.structures.designed, reference=ref, structures=tgt)
        pipeline.save()

    saved = _saved_config(isolated_cwd)
    assert saved["column"] == "designed" and saved["table"].endswith("structures.csv")
    assert "selection" in st.tables.selections.info.columns


@pytest.mark.parametrize("bad", ["", "A10-", "loop1"])
def test_a_malformed_literal_is_refused_at_configuration_time(local_config, isolated_cwd, new_pipeline, bad):
    from biopipelines.mock import Mock
    from biopipelines.selection_transfer import SelectionTransfer

    pipeline = new_pipeline("st_bad")
    with pipeline:
        m = Mock(ids=["d1"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
        with pytest.raises(ValueError, match="not a residue selection"):
            SelectionTransfer(selection=bad, reference=m, structures=m)
