# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""The chain contract: one `sequences` row is one polymer chain.

Covers the producer split in `pipe_fa_to_csv_fasta.py`, the id schemes `chains=` selects, and the `Grouped` axis that folds a design's chain rows back into one prediction.
"""

import json
import os
import subprocess
import sys

import pandas as pd
import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONVERTER = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_fa_to_csv_fasta.py")

AA3 = {'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE', 'G': 'GLY',
       'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN',
       'P': 'PRO', 'Q': 'GLN'}

CHAIN_A = "ACDEFGHIKL"   # 10 residues
CHAIN_B = "MNPQ"         # 4 residues


def _pdb(path, chains):
    """A CA-only PDB carrying the given {chain: sequence}."""
    lines = []
    serial = 1
    for chain, seq in chains.items():
        for i, aa in enumerate(seq, 1):
            lines.append(
                f"ATOM  {serial:5d}  CA  {AA3[aa]} {chain}{i:4d}    "
                f"{i * 3.8:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C")
            serial += 1
    open(path, "w").write("\n".join(lines) + "\nEND\n")


def _run_converter(tmp_path, chains_spec, records, pdb_chains, separator="/"):
    """Convert one .fa holding `records` designs, and return (sequences, designs, missing)."""
    seqs_dir = os.path.join(tmp_path, "seqs")
    out_dir = os.path.join(tmp_path, "out")
    os.makedirs(seqs_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)

    pdb_path = os.path.join(tmp_path, "bb.pdb")
    _pdb(pdb_path, pdb_chains)

    native = separator.join(pdb_chains.values())
    lines = [">bb, score=1.0, seq_recovery=0.5", native]
    for n, segments in enumerate(records, start=1):
        lines.append(f">T=0.1, sample={n}, score=0.9, seq_recovery=0.42")
        lines.append(separator.join(segments))
    open(os.path.join(seqs_dir, "bb.fa"), "w").write("\n".join(lines) + "\n")

    ds = os.path.join(tmp_path, "ds.json")
    json.dump({"name": "structures", "ids": ["bb"], "files": [pdb_path],
               "map_table": "", "format": "pdb"}, open(ds, "w"))

    seq_csv = os.path.join(out_dir, "sequences.csv")
    designs_csv = os.path.join(out_dir, "designs.csv")
    missing_csv = os.path.join(out_dir, "missing.csv")
    cmd = [sys.executable, CONVERTER, seqs_dir, seq_csv,
           os.path.join(out_dir, "sequences.fasta"),
           "--ds-json", ds, "--designs-csv", designs_csv,
           "--missing-csv", missing_csv, "--step-tool-name", "004_ProteinMPNN"]
    if chains_spec:
        cmd += ["--chains", chains_spec]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr

    read = lambda p: pd.read_csv(p) if os.path.exists(p) else pd.DataFrame()
    return read(seq_csv), read(designs_csv), read(missing_csv)


# ── the producer split ────────────────────────────────────────────────────────

def test_single_chain_backbone_keeps_the_ids_it_always_had(record_case, tmp_path):
    seqs, designs, missing = _run_converter(
        str(tmp_path), None, [["AAAAAAAAAA"], ["DDDDDDDDDD"]], {"A": CHAIN_A})
    record_case(input="1 chain, chains=None", expected=["bb_1", "bb_2"],
                actual=list(seqs["id"]))
    assert list(seqs["id"]) == ["bb_1", "bb_2"]
    assert list(designs["id"]) == ["bb_1", "bb_2"]
    assert missing.empty


def test_multi_chain_without_chains_is_a_named_failure(record_case, tmp_path):
    seqs, designs, missing = _run_converter(
        str(tmp_path), None, [["AAAAAAAAAA", "CCCC"]], {"A": CHAIN_A, "B": CHAIN_B})
    record_case(input="2 chains, chains=None", expected="failure row naming the chains",
                actual=list(missing["cause"]))
    assert seqs.empty and designs.empty
    assert list(missing["kind"]) == ["failure"]
    assert "A+B" in missing["cause"].iloc[0]
    # Not excused by the completion check, so the step is reported FAILED.
    assert missing["removed_by"].iloc[0] == "004_ProteinMPNN"


def test_one_named_chain_emits_that_chain_under_the_unsuffixed_id(record_case, tmp_path):
    seqs, _, _ = _run_converter(
        str(tmp_path), "A", [["AAAAAAAAAA", "CCCC"]], {"A": CHAIN_A, "B": CHAIN_B})
    record_case(input='chains="A"', expected=[("bb_1", "AAAAAAAAAA")],
                actual=list(zip(seqs["id"], seqs["sequence"])))
    assert list(seqs["id"]) == ["bb_1"]
    assert list(seqs["sequence"]) == ["AAAAAAAAAA"]
    assert list(seqs["chain"]) == ["A"]


def test_named_chains_emit_one_suffixed_row_each(record_case, tmp_path):
    seqs, designs, _ = _run_converter(
        str(tmp_path), "A,B", [["AAAAAAAAAA", "CCCC"]], {"A": CHAIN_A, "B": CHAIN_B})
    record_case(input='chains="A,B"', expected=["bb_1_A", "bb_1_B"], actual=list(seqs["id"]))
    assert list(seqs["id"]) == ["bb_1_A", "bb_1_B"]
    assert list(seqs["sequence"]) == ["AAAAAAAAAA", "CCCC"]
    assert list(seqs["design"]) == ["bb_1", "bb_1"]
    # One design however many chains it split into.
    assert list(designs["id"]) == ["bb_1"]
    assert list(designs["n_chains"]) == [2]


def test_all_resolves_chain_letters_from_the_backbone(record_case, tmp_path):
    seqs, _, _ = _run_converter(
        str(tmp_path), "all", [["AAAAAAAAAA", "CCCC"]], {"A": CHAIN_A, "B": CHAIN_B})
    record_case(input='chains="all"', expected=["A", "B"], actual=list(seqs["chain"]))
    assert list(seqs["chain"]) == ["A", "B"]
    assert list(seqs["id"]) == ["bb_1_A", "bb_1_B"]


def test_ligandmpnn_colon_separator_splits_the_same_way(record_case, tmp_path):
    seqs, _, _ = _run_converter(
        str(tmp_path), "A,B", [["AAAAAAAAAA", "CCCC"]], {"A": CHAIN_A, "B": CHAIN_B},
        separator=":")
    record_case(input="':' separator", expected=["bb_1_A", "bb_1_B"], actual=list(seqs["id"]))
    assert list(seqs["id"]) == ["bb_1_A", "bb_1_B"]


def test_segments_that_do_not_match_the_backbone_are_dropped_not_mislabelled(record_case, tmp_path):
    # Three segments against a two-chain backbone: no assignment is defensible.
    seqs, _, missing = _run_converter(
        str(tmp_path), "all", [["AAAAAAAAAA", "CCCC", "GG"]], {"A": CHAIN_A, "B": CHAIN_B})
    record_case(input="3 segments, 2 chains", expected="dropped", actual=len(seqs))
    assert seqs.empty
    assert list(missing["kind"]) == ["failure"]
    assert "cannot match" in missing["cause"].iloc[0]


def test_a_homodimers_identical_chains_both_survive(record_case, tmp_path):
    # Dedup is per design, not per chain, or the second copy would vanish.
    seqs, _, _ = _run_converter(
        str(tmp_path), "A,B", [["AAAAAAAAAA", "AAAAAAAAAA"]], {"A": CHAIN_A, "B": CHAIN_A})
    record_case(input="homodimer", expected=2, actual=len(seqs))
    assert len(seqs) == 2
    assert list(seqs["chain"]) == ["A", "B"]


# ── the id schemes the wrapper declares ──────────────────────────────────────

@pytest.mark.parametrize("chains,expected", [
    (None, ["bb_<1..2>"]),
    ("A", ["bb_<1..2>"]),
    (["A", "B"], ["bb_<1..2>_<A B>"]),
    ("all", ["bb_<1..2>[_<?>]"]),
])
def test_declared_sequence_ids_follow_the_chains_setting(record_case, chains, expected):
    from biopipelines.chain_rows import chain_row_ids, normalize_chains
    actual = chain_row_ids(["bb_<1..2>"], normalize_chains(chains))
    record_case(input=f"chains={chains!r}", expected=expected, actual=actual)
    assert actual == expected


def test_only_all_is_lazy(record_case):
    from biopipelines.chain_rows import chain_row_ids, normalize_chains
    from biopipelines.id_patterns import is_lazy
    lazy = {repr(c): is_lazy(chain_row_ids(["bb_<1..2>"], normalize_chains(c))[0])
            for c in (None, "A", ["A", "B"], "all")}
    record_case(input="is_lazy per chains setting", expected="only 'all'", actual=lazy)
    assert lazy == {"None": False, "'A'": False, "['A', 'B']": False, "'all'": True}


@pytest.mark.parametrize("bad", ["", [], ["AB"], ["A", "A"], 3])
def test_invalid_chains_are_refused_at_construction(record_case, bad):
    from biopipelines.chain_rows import validate_chains, normalize_chains
    with pytest.raises(ValueError):
        validate_chains(normalize_chains(bad))
    record_case(input=f"chains={bad!r}", expected="ValueError", actual="ValueError")


# ── the Grouped axis ─────────────────────────────────────────────────────────

def test_grouped_iterates_group_ids_not_member_ids(record_case):
    from biopipelines.combinatorics import Grouped, predict_output_ids
    from biopipelines.datastream import DataStream

    members = DataStream(name="sequences", ids=["d_1_A", "d_1_B", "d_2_A", "d_2_B"],
                         files=[], map_table="/tmp/seq.csv", format="csv")
    groups = DataStream(name="designs", ids=["d_1", "d_2"],
                        files=[], map_table="/tmp/designs.csv", format="csv")
    actual = predict_output_ids(proteins=(Grouped(members, groups=groups), "sequences"))
    record_case(input="Grouped over 4 chain rows in 2 designs",
                expected=["d_1", "d_2"], actual=actual)
    assert actual == ["d_1", "d_2"]


def test_grouped_inside_bundle_keeps_a_static_partner(record_case):
    from biopipelines.combinatorics import Bundle, Grouped, predict_output_ids
    from biopipelines.datastream import DataStream

    members = DataStream(name="sequences", ids=["d_1_A", "d_1_B", "d_2_A", "d_2_B"],
                         files=[], map_table="/tmp/seq.csv", format="csv")
    groups = DataStream(name="designs", ids=["d_1", "d_2"],
                        files=[], map_table="/tmp/designs.csv", format="csv")
    tag = DataStream(name="sequences", ids=["tag"], files=[],
                     map_table="/tmp/tag.csv", format="csv")
    actual = predict_output_ids(
        proteins=(Bundle(Grouped(members, groups=groups), tag), "sequences"))
    record_case(input="Bundle(Grouped(...), tag)",
                expected=["d_1+tag", "d_2+tag"], actual=actual)
    assert actual == ["d_1+tag", "d_2+tag"]


def test_a_bare_stream_cannot_be_grouped_without_groups(record_case):
    from biopipelines.combinatorics import Grouped, predict_output_ids
    from biopipelines.datastream import DataStream

    members = DataStream(name="sequences", ids=["d_1_A"], files=[],
                         map_table="/tmp/seq.csv", format="csv")
    with pytest.raises(ValueError, match="groups="):
        predict_output_ids(proteins=(Grouped(members), "sequences"))
    record_case(input="Grouped(bare DataStream)", expected="ValueError", actual="ValueError")


def test_group_membership_survives_composed_and_renamed_ids(record_case):
    """The ids need not look alike: matching walks the parent chain, not the spelling."""
    from biopipelines.id_map_utils import get_mapped_ids

    groups = ["chihuahua+labradudor_1", "chihuahua+labradudor_2"]
    members = ["chihuahua+labradudor_1_A", "chihuahua+labradudor_1_B",
               "chihuahua+labradudor_2_A"]
    actual = get_mapped_ids(groups, members, unique=False)
    record_case(input="'+'-composed bundle ids", expected="2 + 1 members",
                actual={k: len(v) for k, v in actual.items()})
    assert actual[groups[0]] == members[:2]
    assert actual[groups[1]] == members[2:]


# ── the consumers, end to end ────────────────────────────────────────────────

def _grouped_combinatorics(tmp_path):
    """Two two-chain designs plus a fixed partner, as a combinatorics config."""
    import csv
    base = str(tmp_path)
    with open(os.path.join(base, "sequences.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "design", "chain", "sequence"])
        w.writerow(["d_1_A", "d_1", "A", "AAAAAAAAAA"])
        w.writerow(["d_1_B", "d_1", "B", "CCCC"])
        w.writerow(["d_2_A", "d_2", "A", "DDDDDDDDDD"])
        w.writerow(["d_2_B", "d_2", "B", "EEEE"])
    with open(os.path.join(base, "designs.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "n_chains"]); w.writerow(["d_1", 2]); w.writerow(["d_2", 2])
    with open(os.path.join(base, "tag.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "sequence"]); w.writerow(["tag", "CHRVVQGDLDVG"])

    cfg = {"axes": {"proteins": {
        "name": "proteins", "mode": "bundle", "entity_type": "protein", "sources": [
            {"path": os.path.join(base, "sequences.csv"), "iterate": True, "order": 0,
             "ids": ["d_<1 2>_<A B>"],
             "group_by": os.path.join(base, "designs.csv"), "group_ids": ["d_<1 2>"]},
            {"path": os.path.join(base, "tag.csv"), "iterate": False, "order": 1},
        ]}}, "predicted_ids": ["d_1+tag", "d_2+tag"], "provenance": {}}
    path = os.path.join(base, "comb.json")
    json.dump(cfg, open(path, "w"))
    return path


def test_boltz_writes_one_chain_entry_per_group_member(record_case, tmp_path):
    """The exact failure the contract exists for: a fused sequence became a UNK residue."""
    import yaml
    cfg = _grouped_combinatorics(tmp_path)
    out = os.path.join(str(tmp_path), "out")
    pipe = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_boltz_config_unified.py")
    r = subprocess.run([sys.executable, pipe, "--combinatorics-config", cfg,
                        "--output-dir", out], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr

    written = yaml.safe_load(open(os.path.join(out, "d_1+tag.yaml")))
    entries = [(e["protein"]["id"], e["protein"]["sequence"]) for e in written["sequences"]]
    record_case(input="Bundle(Grouped(2-chain design), tag)",
                expected=[("A", "AAAAAAAAAA"), ("B", "CCCC"), ("C", "CHRVVQGDLDVG")],
                actual=entries)
    assert entries == [("A", "AAAAAAAAAA"), ("B", "CCCC"), ("C", "CHRVVQGDLDVG")]
    # No entry may carry a chain separator: that is what Boltz turns into UNK.
    assert not any("/" in seq or ":" in seq for _, seq in entries)
    # The partner appears once, not once per chain of the design.
    assert [seq for _, seq in entries].count("CHRVVQGDLDVG") == 1


def test_alphafold_colon_joins_a_group_then_its_partner(record_case, tmp_path):
    cfg = _grouped_combinatorics(tmp_path)
    out_csv = os.path.join(str(tmp_path), "queries.csv")
    pipe = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_alphafold_queries.py")
    r = subprocess.run([sys.executable, pipe, "--combinatorics-config", cfg,
                        "--output", out_csv], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr

    rows = pd.read_csv(out_csv)
    actual = dict(zip(rows["id"], rows["sequence"]))
    record_case(input="Bundle(Grouped(...), tag) into ColabFold",
                expected="chains colon-joined, partner last", actual=actual)
    assert actual["d_1+tag"] == "AAAAAAAAAA:CCCC:CHRVVQGDLDVG"
    assert actual["d_2+tag"] == "DDDDDDDDDD:EEEE:CHRVVQGDLDVG"


# ── grouping at any ancestor level ───────────────────────────────────────────

def _streams():
    from biopipelines.datastream import DataStream
    from biopipelines.outputs import StandardizedOutput
    members = DataStream(name="sequences", ids=["A_1_1", "A_1_2", "A_2_1", "A_2_2"],
                         files=[], map_table="/x/s.csv", format="csv")
    designs = DataStream(name="designs", ids=["A_1", "A_2"],
                         files=[], map_table="/x/d.csv", format="csv")
    backbone = StandardizedOutput({
        "structures": DataStream(name="structures", ids=["A"], files=["<id>.pdb"],
                                 map_table="/x/st.csv", format="pdb"),
        "tables": {}, "output_folder": "/x"})
    return members, designs, backbone


def test_groups_can_be_the_parent_or_any_ancestor(record_case):
    """`A_1_1` groups under its design `A_1` or its backbone `A` — same call, different groups."""
    from biopipelines.combinatorics import Grouped, predict_output_ids
    members, designs, backbone = _streams()
    by_design = predict_output_ids(proteins=(Grouped(members, groups=designs), "sequences"))
    by_backbone = predict_output_ids(proteins=(Grouped(members, groups=backbone), "sequences"))
    record_case(input="4 rows, grouped two ways",
                expected=(["A_1", "A_2"], ["A"]), actual=(by_design, by_backbone))
    assert by_design == ["A_1", "A_2"]
    assert by_backbone == ["A"]


def test_a_grouping_output_need_not_carry_the_consumed_stream_name(record_case):
    """An RFdiffusion output groups `sequences` by its `structures`; that must not be refused."""
    from biopipelines.combinatorics import Grouped, predict_output_ids
    members, _, backbone = _streams()
    whole = predict_output_ids(proteins=(Grouped(members, groups=backbone), "sequences"))
    named = predict_output_ids(
        proteins=(Grouped(members, groups=backbone.streams.structures), "sequences"))
    record_case(input="groups=<output with only 'structures'>", expected=["A"],
                actual=(whole, named))
    assert whole == named == ["A"]


def test_an_ambiguous_grouping_output_is_refused_not_guessed(record_case):
    from biopipelines.combinatorics import Grouped, predict_output_ids
    from biopipelines.datastream import DataStream
    from biopipelines.outputs import StandardizedOutput
    members, _, _ = _streams()
    ambiguous = StandardizedOutput({
        "structures": DataStream(name="structures", ids=["A"], files=["<id>.pdb"],
                                 map_table="/x/st.csv", format="pdb"),
        "sequences": DataStream(name="sequences", ids=["A_1", "A_2"], files=[],
                                map_table="/x/sq.csv", format="csv"),
        "tables": {}, "output_folder": "/x"})
    with pytest.raises(ValueError, match="disagree on ids"):
        predict_output_ids(proteins=(Grouped(members, groups=ambiguous), "sequences"))
    record_case(input="groups=<output whose streams disagree>", expected="ValueError",
                actual="ValueError")


def test_runtime_membership_reaches_a_grandparent(record_case):
    from biopipelines.id_map_utils import get_mapped_ids
    actual = get_mapped_ids(["A", "B"], ["A_1_1", "A_2_1", "B_1_1"], unique=False)
    record_case(input="grandparent keys", expected={"A": 2, "B": 1},
                actual={k: len(v) for k, v in actual.items()})
    assert actual == {"A": ["A_1_1", "A_2_1"], "B": ["B_1_1"]}


# ── one `chains` parameter, and per-chain selections ─────────────────────────

FIXED_POS = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_pmpnn_fixed_positions.py")


def _run_fixed_positions(tmp_path, fixed, designed, fixed_chain, pdb_chains):
    base = str(tmp_path)
    os.makedirs(base, exist_ok=True)
    pdb = os.path.join(base, "bb.pdb")
    _pdb(pdb, pdb_chains)
    ds = os.path.join(base, "ds.json")
    json.dump({"name": "structures", "ids": ["bb"], "files": [pdb],
               "map_table": "", "format": "pdb"}, open(ds, "w"))
    cfg = os.path.join(base, "cfg.json")
    json.dump({"structures_json": ds, "FIXED": fixed, "DESIGNED": designed,
               "FIXED_CHAIN": fixed_chain,
               "fixed_jsonl_file": os.path.join(base, "fixed.jsonl"),
               "sele_csv_file": os.path.join(base, "sele.csv")}, open(cfg, "w"))
    result = subprocess.run([sys.executable, FIXED_POS, cfg], capture_output=True, text=True)
    return result, base


def test_a_chainless_selection_on_a_multichain_structure_raises(record_case, tmp_path):
    """It used to pick the first chain and say so in a print nobody reads."""
    result, _ = _run_fixed_positions(tmp_path, "-", "3-5", "auto",
                                     {"A": CHAIN_A, "B": CHAIN_B})
    record_case(input="chainless redesigned, 2 chains", expected="ValueError",
                actual=result.returncode)
    assert result.returncode != 0
    assert "which chain do the residues belong to" in result.stderr


def test_a_per_chain_dict_selects_on_each_chain(record_case, tmp_path):
    result, base = _run_fixed_positions(tmp_path, "-", {"A": "3-5", "B": "2"}, "auto",
                                        {"A": CHAIN_A, "B": CHAIN_B})
    assert result.returncode == 0, result.stderr
    fixed = json.loads(open(os.path.join(base, "fixed.jsonl")).read())
    record_case(input='redesigned={"A": "3-5", "B": "2"}',
                expected={"A": [1, 2, 6, 7, 8, 9, 10], "B": [1, 3, 4]}, actual=fixed["bb"])
    # Positions are 1-based ranks within the chain, and each chain's complement is its own.
    assert fixed["bb"] == {"A": [1, 2, 6, 7, 8, 9, 10], "B": [1, 3, 4]}


def test_a_multichain_summary_is_chain_aware_and_a_single_chain_one_is_not(record_case, tmp_path):
    """Inter-tool selection columns must carry chains; the single-chain form is unchanged."""
    _, multi = _run_fixed_positions(tmp_path / "m", "-", {"A": "3-5", "B": "2"}, "auto",
                                    {"A": CHAIN_A, "B": CHAIN_B})
    _, single = _run_fixed_positions(tmp_path / "s", "-", "3-5", "A", {"A": CHAIN_A})
    multi_row = pd.read_csv(os.path.join(multi, "sele.csv")).iloc[0]
    single_row = pd.read_csv(os.path.join(single, "sele.csv")).iloc[0]
    record_case(input="sele.csv summaries",
                expected=("A1-2+A6-10+B1+B3-4", "1-2+6-10"),
                actual=(multi_row["fixed"], single_row["fixed"]))
    assert multi_row["fixed"] == "A1-2+A6-10+B1+B3-4"
    assert multi_row["mobile"] == "A3-5+B2"
    assert single_row["fixed"] == "1-2+6-10"


@pytest.mark.parametrize("tool_name", ["ProteinMPNN", "LigandMPNN", "SolubleMPNN"])
def test_chain_is_gone_and_only_bound_as_a_deprecated_alias(record_case, tool_name):
    import inspect
    import biopipelines
    tool = getattr(biopipelines, tool_name)
    params = list(inspect.signature(tool.__init__).parameters)
    record_case(input=f"{tool_name} signature", expected="chains present, chain absent",
                actual=[p for p in params if p == "chain" or p == "chains"])
    assert "chains" in params
    assert "chain" not in params
    assert tool.PARAMETER_ALIASES.get("chain") == "chains"
    assert "chain" in tool.DEPRECATED_ALIASES


def test_the_retired_spelling_still_binds(record_case):
    from biopipelines.protein_mpnn import ProteinMPNN
    from biopipelines.datastream import DataStream
    ds = DataStream(name="structures", ids=["bb"], files=["/x/bb.pdb"],
                    map_table="/x/m.csv", format="pdb")
    bound = {spelling: ProteinMPNN(structures=ds, chain=spelling).chains
             for spelling in ("A", "auto")}
    record_case(input="chain= on a new wrapper", expected={"A": ["A"], "auto": None},
                actual=bound)
    # "auto" was the old "work it out from the structure", which is what None means now.
    assert bound == {"A": ["A"], "auto": None}


def test_a_per_chain_dict_is_validated_at_construction(record_case):
    from biopipelines.data_containers import resolve_table_reference
    with pytest.raises(ValueError, match="single-character chain ids"):
        resolve_table_reference({"AB": "1-5"}, "redesigned")
    record_case(input='redesigned={"AB": ...}', expected="ValueError", actual="ValueError")


# ── LigandMPNN's chainless-position chain ────────────────────────────────────

LMPNN_POS = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_lmpnn_runtime_positions.py")


def _run_lmpnn_positions(tmp_path, pdb_chains, designed="3-5", default_chain="auto"):
    base = str(tmp_path)
    os.makedirs(base, exist_ok=True)
    pdb = os.path.join(base, "bb.pdb")
    _pdb(pdb, pdb_chains)
    ds = os.path.join(base, "ds.json")
    json.dump({"name": "structures", "ids": ["bb"], "files": [pdb],
               "map_table": "", "format": "pdb"}, open(ds, "w"))
    cfg = os.path.join(base, "cfg.json")
    out = os.path.join(base, "out.json")
    json.dump({"structures_json": ds, "input_source": "selection", "input_table": "-",
               "fixed_positions": "-", "designed_positions": designed, "ligand": "-",
               "design_within": "5.0", "output_file": out,
               "default_chain": default_chain}, open(cfg, "w"))
    result = subprocess.run([sys.executable, LMPNN_POS, cfg], capture_output=True, text=True)
    written = json.load(open(out)) if os.path.exists(out) else None
    return result, written


def test_a_chainless_position_takes_the_structures_own_chain(record_case, tmp_path):
    """The old default was the literal "A", which was wrong whenever the chain was not A."""
    result, written = _run_lmpnn_positions(tmp_path, {"B": CHAIN_A})
    assert result.returncode == 0, result.stderr
    actual = written["bb"]["redesigned_option"]
    record_case(input="single chain named B, redesigned='3-5'",
                expected='--redesigned_residues "B3 B4 B5"', actual=actual)
    assert actual == '--redesigned_residues "B3 B4 B5"'


def test_lmpnn_chainless_positions_on_a_multichain_structure_raise(record_case, tmp_path):
    result, written = _run_lmpnn_positions(tmp_path, {"A": CHAIN_A, "B": CHAIN_B})
    record_case(input="2 chains, chainless redesigned", expected="ValueError",
                actual=result.returncode)
    assert result.returncode != 0 and written is None
    assert "which chain do the residues belong to" in result.stderr


def test_a_chain_qualified_position_is_untouched_on_a_multichain_structure(record_case, tmp_path):
    """Qualified positions are unambiguous, so several chains is not an error."""
    result, written = _run_lmpnn_positions(tmp_path, {"A": CHAIN_A, "B": CHAIN_B},
                                           designed="A3-5+B2")
    assert result.returncode == 0, result.stderr
    actual = written["bb"]["redesigned_option"]
    record_case(input="2 chains, redesigned='A3-5+B2'",
                expected='--redesigned_residues "A3 A4 A5 B2"', actual=actual)
    assert actual == '--redesigned_residues "A3 A4 A5 B2"'


def test_a_chain_qualified_string_needs_no_dict(record_case, tmp_path):
    """The dict is for attaching a different source per chain, not a workaround.

    ProteinMPNN's script aliased the shared chain-aware parser away and defined a local
    `sele_to_list` of the same name that stripped the chain, so a qualified selection
    silently collapsed onto one chain while LigandMPNN handled it correctly.
    """
    _, qualified = _run_fixed_positions(tmp_path / "q", "-", "A3-5+B2", "auto",
                                        {"A": CHAIN_A, "B": CHAIN_B})
    _, via_dict = _run_fixed_positions(tmp_path / "d", "-", {"A": "3-5", "B": "2"}, "auto",
                                       {"A": CHAIN_A, "B": CHAIN_B})
    a = json.loads(open(os.path.join(qualified, "fixed.jsonl")).read())
    b = json.loads(open(os.path.join(via_dict, "fixed.jsonl")).read())
    record_case(input='"A3-5+B2" vs {"A": "3-5", "B": "2"}', expected="identical",
                actual=(a, b))
    assert a == b == {"bb": {"A": [1, 2, 6, 7, 8, 9, 10], "B": [1, 3, 4]}}


def test_both_mpnn_scripts_agree_on_a_qualified_selection(record_case, tmp_path):
    """The two used to disagree; that is the whole point of the shared parser."""
    _, pmpnn_dir = _run_fixed_positions(tmp_path / "p", "-", "A3-5", "auto",
                                        {"A": CHAIN_A, "B": CHAIN_B})
    pmpnn = json.loads(open(os.path.join(pmpnn_dir, "fixed.jsonl")).read())["bb"]
    _, lmpnn = _run_lmpnn_positions(tmp_path / "l", {"A": CHAIN_A, "B": CHAIN_B},
                                    designed="A3-5")
    record_case(input='redesigned="A3-5" on a 2-chain backbone',
                expected="both place it on chain A",
                actual=(pmpnn, lmpnn["bb"]["redesigned_option"]))
    # ProteinMPNN fixes A's complement; LigandMPNN names the redesigned residues directly.
    assert pmpnn["A"] == [1, 2, 6, 7, 8, 9, 10]
    assert lmpnn["bb"]["redesigned_option"] == '--redesigned_residues "A3 A4 A5"'


def test_grouped_accepts_a_bare_stream_from_a_tool(record_case, isolated_cwd, new_pipeline):
    """`Grouped(pmpnn.streams.sequences)` is the spelling that reads most naturally.

    A DataStream carries `_producer`, the tool that made it, so the sibling `designs`
    stream is reachable without the caller naming it.
    """
    from biopipelines.pipeline import Resources, Save
    from biopipelines import ProteinMPNN, Mock, Grouped
    from biopipelines.combinatorics import predict_output_ids

    with new_pipeline("grouped_bare"):
        Resources()
        bb = Mock(streams={"structures": {"format": "pdb"}}, ids=["bb1", "bb2"])
        pm = ProteinMPNN(structures=bb.streams.structures, num_sequences=2, chains=["A", "B"])
        bare = predict_output_ids(proteins=(Grouped(pm.streams.sequences), "sequences"))
        whole = predict_output_ids(proteins=(Grouped(pm), "sequences"))
        Save()

    record_case(input="Grouped(stream) vs Grouped(tool)",
                expected=["bb1_<1..2>", "bb2_<1..2>"], actual=(bare, whole))
    assert bare == whole == ["bb1_<1..2>", "bb2_<1..2>"]


def test_a_stream_with_no_producer_still_says_what_to_pass(record_case):
    from biopipelines.combinatorics import Grouped, predict_output_ids
    from biopipelines.datastream import DataStream
    orphan = DataStream(name="sequences", ids=["d_1_A"], files=[],
                        map_table="/x/s.csv", format="csv")
    with pytest.raises(ValueError, match="groups="):
        predict_output_ids(proteins=(Grouped(orphan), "sequences"))
    record_case(input="Grouped(stream with no producer)", expected="ValueError naming groups=",
                actual="ValueError")


def test_a_grouped_item_in_a_bare_list_is_unwrapped():
    from biopipelines.combinatorics import Grouped, _unwrap_sources
    from biopipelines.datastream import DataStream
    members = DataStream(name="sequences", ids=["A_1_A", "A_1_B"], files=[],
                         map_table="/x/s.csv", format="csv")
    designs = DataStream(name="designs", ids=["A_1"], files=[], map_table="/x/d.csv", format="csv")
    other = DataStream(name="sequences", ids=["tag"], files=[], map_table="/x/t.csv", format="csv")
    mode, entries = _unwrap_sources([Grouped(members, groups=designs), other], "sequences")
    assert mode == "each"
    assert [bool(e.get("group_by")) for e in entries] == [True, False]


def test_esmfold2_folds_a_group_as_one_complex(record_case, tmp_path):
    """ESMFold2 accepted Grouped at configuration time; its runtime used to fold each chain row alone."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "pipe_esmfold2_inference", os.path.join(REPO_ROOT, "pipe_scripts", "pipe_esmfold2_inference.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    config = json.load(open(_grouped_combinatorics(tmp_path)))
    complexes = {cid: [c["sequence"] for c in chains]
                 for cid, chains in module.build_complexes(config, {})}
    record_case(input="Bundle(Grouped(2-chain design), tag) into ESMFold2",
                expected="one complex per design, partner last", actual=complexes)
    assert complexes == {"d_1+tag": ["AAAAAAAAAA", "CCCC", "CHRVVQGDLDVG"],
                         "d_2+tag": ["DDDDDDDDDD", "EEEE", "CHRVVQGDLDVG"]}


def test_equal_length_chains_are_assigned_by_sequence_not_by_alphabet(tmp_path):
    """Chains of equal length match both candidate orders by length; alphabetical always won,
    so a producer writing appearance order (B before A) had its two sequences swapped."""
    seqs, _designs, _missing = _run_converter(
        str(tmp_path), "A,B", [["GGGGA", "KKKKA"]], {"B": "GGGGM", "A": "KKKKM"})
    by_id = dict(zip(seqs["id"], seqs["sequence"]))
    assert by_id["bb_1_B"] == "GGGGA" and by_id["bb_1_A"] == "KKKKA"


def test_naming_several_chains_does_not_attach_a_chainless_selection_to_the_first():
    """chains=["A", "B"] used to set the default chain to "A", so an unqualified fixed= went to
    chain A silently; several named chains are exactly the ambiguous case."""
    from biopipelines.chain_rows import positions_chain
    assert positions_chain(["A", "B"]) == "auto"
    assert positions_chain(["B"]) == "B"
    assert positions_chain("all") == "auto"
