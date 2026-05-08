"""Unit tests for the Sequence input-type tool.

Covers the config-time construction paths (raw string, list, dict, CSV path,
FASTA path, RCSB PDB code — network; type auto-detection, default IDs,
validation errors) and integration with Mock and Panda downstream.
"""

import os

import pytest


# ── raw-sequence construction ────────────────────────────────────────────────

def test_sequence_single_string_default_id(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    from biopipelines.sequence import Sequence

    pipeline = new_pipeline("seq_single_string")
    with pipeline:
        s = Sequence(seq="MKTAYIAKQRQISFVKSHFSRQLE")
        script_path = pipeline.save()

    ids = list(s.streams.sequences.ids)
    record_case(input="Sequence(seq='MKTAY...')",
                expected=(["seq1"], "protein"),
                actual=(ids, s.tables.sequences.info.name))
    assert_valid_script(script_path, "Sequence")
    assert ids == ["seq1"]


def test_sequence_list_with_custom_ids(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    from biopipelines.sequence import Sequence

    pipeline = new_pipeline("seq_list_ids")
    with pipeline:
        s = Sequence(seq=["MKTAY", "AETGF", "GGGGA"],
                     ids=["a", "b", "c"], type="protein")
        pipeline.save()

    ids = list(s.streams.sequences.ids)
    record_case(input="Sequence(list, ids=[a,b,c])",
                expected=["a", "b", "c"], actual=ids)
    assert ids == ["a", "b", "c"]


def test_sequence_list_default_ids(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """Default IDs when none provided: seq1, seq2, seq3."""
    from biopipelines.sequence import Sequence

    pipeline = new_pipeline("seq_default_ids")
    with pipeline:
        s = Sequence(seq=["MKTAY", "AETGF", "GGGGA"], type="protein")
        pipeline.save()

    ids = list(s.streams.sequences.ids)
    record_case(input="Sequence(list, ids=None)",
                expected=["seq1", "seq2", "seq3"], actual=ids)
    assert ids == ["seq1", "seq2", "seq3"]


def test_sequence_dict_uses_keys_as_ids(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """{id: seq} dict — keys become IDs, values become sequences."""
    from biopipelines.sequence import Sequence

    pipeline = new_pipeline("seq_dict")
    with pipeline:
        s = Sequence(seq={"proteinA": "MKTAY", "proteinB": "AETGF"},
                     type="protein")
        pipeline.save()

    ids = list(s.streams.sequences.ids)
    record_case(input="Sequence(seq={proteinA:..., proteinB:...})",
                expected=["proteinA", "proteinB"], actual=ids)
    assert ids == ["proteinA", "proteinB"]


def test_sequence_dict_ids_kwarg_is_ignored(
    local_config, isolated_cwd, new_pipeline, record_case, capsys,
):
    """When seq is a dict, ids= is ignored (dict keys win)."""
    from biopipelines.sequence import Sequence

    pipeline = new_pipeline("seq_dict_ids_ignored")
    with pipeline:
        s = Sequence(seq={"A": "MKTAY", "B": "AETGF"},
                     ids=["X", "Y"], type="protein")
        pipeline.save()

    ids = list(s.streams.sequences.ids)
    record_case(input="dict + ids= (ids ignored)",
                expected=["A", "B"], actual=ids)
    assert ids == ["A", "B"]


# ── type detection ───────────────────────────────────────────────────────────

def test_sequence_auto_detects_dna(record_case):
    """Construct outside a Pipeline to inspect raw Sequence attributes."""
    from biopipelines.sequence import Sequence

    s = Sequence(seq="ACGTACGTACGT")
    record_case(input="Sequence(seq='ACGTACGT', type='auto')",
                expected=["dna"], actual=s.detected_types)
    assert s.detected_types == ["dna"]


def test_sequence_auto_detects_rna(record_case):
    from biopipelines.sequence import Sequence

    s = Sequence(seq="ACGUACGUACGU")
    record_case(input="Sequence(seq='ACGUACGU', type='auto')",
                expected=["rna"], actual=s.detected_types)
    assert s.detected_types == ["rna"]


def test_sequence_auto_detects_protein(record_case):
    """Anything outside ACGT/ACGU falls through to protein."""
    from biopipelines.sequence import Sequence

    s = Sequence(seq="MKTAYIAKQRQ")
    record_case(input="Sequence(seq='MKTAY...', type='auto')",
                expected=["protein"], actual=s.detected_types)
    assert s.detected_types == ["protein"]


def test_sequence_auto_mixed_types_in_list(record_case):
    """A mixed list auto-detects each entry independently."""
    from biopipelines.sequence import Sequence

    s = Sequence(seq=["ACGT", "ACGU", "MKTAY"])
    record_case(input="list with DNA/RNA/protein",
                expected=["dna", "rna", "protein"],
                actual=s.detected_types)
    assert s.detected_types == ["dna", "rna", "protein"]


# ── CSV loading ──────────────────────────────────────────────────────────────

def test_sequence_loads_from_csv(isolated_cwd, record_case):
    """CSV with 'id' and 'sequence' columns — custom_ids come from the file."""
    from biopipelines.sequence import Sequence

    csv_path = isolated_cwd / "seqs.csv"
    csv_path.write_text("id,sequence\np1,MKTAY\np2,AETGF\np3,GGGGA\n")

    s = Sequence(seq=str(csv_path), type="protein")
    record_case(input="Sequence(seq='seqs.csv')",
                expected=(["p1", "p2", "p3"], ["MKTAY", "AETGF", "GGGGA"], True),
                actual=(list(s.custom_ids), list(s.sequences), s.from_csv))
    assert s.custom_ids == ["p1", "p2", "p3"]
    assert s.sequences == ["MKTAY", "AETGF", "GGGGA"]
    assert s.from_csv is True


def test_sequence_csv_preserves_extra_columns(isolated_cwd, record_case):
    """CSV columns other than id/sequence ride through as extra_columns."""
    from biopipelines.sequence import Sequence

    csv_path = isolated_cwd / "seqs_extra.csv"
    csv_path.write_text(
        "id,sequence,source,chain\n"
        "p1,MKTAY,rcsb,A\n"
        "p2,AETGF,user,B\n"
    )

    s = Sequence(seq=str(csv_path), type="protein")
    extra_keys = sorted(s.extra_columns.keys())
    record_case(input="CSV with extra columns source,chain",
                expected=["chain", "source"], actual=extra_keys)
    assert extra_keys == ["chain", "source"]
    assert s.extra_columns["source"] == ["rcsb", "user"]


def test_sequence_csv_missing_id_column_raises(isolated_cwd, record_case):
    from biopipelines.sequence import Sequence

    csv_path = isolated_cwd / "bad.csv"
    csv_path.write_text("seqid,sequence\np1,MKTAY\n")

    record_case(input="CSV without 'id' column",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError, match="id"):
        Sequence(seq=str(csv_path), type="protein")


def test_sequence_csv_missing_sequence_column_raises(isolated_cwd, record_case):
    from biopipelines.sequence import Sequence

    csv_path = isolated_cwd / "bad.csv"
    csv_path.write_text("id,seq\np1,MKTAY\n")

    record_case(input="CSV without 'sequence' column",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError, match="sequence"):
        Sequence(seq=str(csv_path), type="protein")


# ── FASTA loading ────────────────────────────────────────────────────────────

def test_sequence_loads_from_fasta(isolated_cwd, record_case):
    """FASTA file — headers become IDs, sequences are the blocks that follow."""
    from biopipelines.sequence import Sequence

    fasta_path = isolated_cwd / "seqs.fasta"
    fasta_path.write_text(
        ">p1|Chain A|test\nMKTAY\n"
        ">p2|Chain B|test\nAETGF\n"
    )

    s = Sequence(seq=str(fasta_path), type="protein")
    record_case(input="Sequence(seq='seqs.fasta')",
                expected=(True, 2, ["MKTAY", "AETGF"]),
                actual=(s.from_fasta, len(s.sequences), list(s.sequences)))
    assert s.from_fasta is True
    assert len(s.sequences) == 2
    assert s.sequences == ["MKTAY", "AETGF"]


def test_sequence_fasta_missing_file_raises(isolated_cwd, record_case):
    from biopipelines.sequence import Sequence

    record_case(input="Sequence(seq='absent.fasta') — absolute path",
                expected="FileNotFoundError", actual="FileNotFoundError")
    # Use an absolute path so the tool doesn't defer to configure_inputs.
    missing = str(isolated_cwd / "absent.fasta")
    # Make it absolute and non-existent but tagged as .fasta — tool defers to
    # configure_inputs, so we call the loader directly.
    s = Sequence.__new__(Sequence)
    with pytest.raises(FileNotFoundError):
        s._load_from_fasta(missing)


# ── RCSB PDB-code path (hits network) ────────────────────────────────────────

def test_sequence_fetches_from_rcsb_pdb_code(record_case):
    """Passing a 4-char alphanumeric code fetches chains from the RCSB FASTA API."""
    from biopipelines.sequence import Sequence

    s = Sequence(seq="4EQ7")
    record_case(input="Sequence(seq='4EQ7') — RCSB fetch",
                expected=(True, ">=1"),
                actual=(s.from_pdb_code, len(s.sequences)))
    assert s.from_pdb_code is True
    assert len(s.sequences) >= 1
    assert all(seq for seq in s.sequences)


# ── validation errors ───────────────────────────────────────────────────────

def test_sequence_empty_list_raises(record_case):
    from biopipelines.sequence import Sequence

    record_case(input="Sequence(seq=[])",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError):
        Sequence(seq=[])


def test_sequence_ids_length_mismatch_raises(record_case):
    from biopipelines.sequence import Sequence

    record_case(input="len(ids)=2 vs len(seq)=3",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError, match="Length mismatch"):
        Sequence(seq=["MKTAY", "AETGF", "GGGGA"],
                 ids=["a", "b"], type="protein")


def test_sequence_invalid_type_raises(record_case):
    from biopipelines.sequence import Sequence

    record_case(input="Sequence(type='peptide')",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError, match="Invalid type"):
        Sequence(seq="MKTAY", type="peptide")


def test_sequence_invalid_protein_chars_caught_by_validate(record_case):
    """Invalid characters outside the 20 canonical AAs are rejected at construction.

    validate_params() runs from BaseConfig.__init__, so Sequence(seq='MK*TAY',
    type='protein') itself raises — no separate .validate_params() call needed.
    """
    from biopipelines.sequence import Sequence

    record_case(input="Sequence(seq='MK*TAY', type='protein')",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError, match="invalid characters"):
        Sequence(seq="MK*TAY", type="protein")


# ── integration: Sequence + Mock ─────────────────────────────────────────────

def test_sequence_feeds_mock(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Mock consumes Sequence.streams.sequences, fans out with children=<1..2>."""
    from biopipelines.sequence import Sequence
    from biopipelines.mock import Mock

    pipeline = new_pipeline("seq_feeds_mock")
    with pipeline:
        s = Sequence(seq={"protA": "MKTAY", "protB": "AETGF"}, type="protein")
        m = Mock(
            source=s.streams.sequences,
            children="<1..2>",
            streams={"designs": {"format": "pdb", "file": "<id>.pdb"}},
        )
        script_path = pipeline.save()

    ids = list(m.streams.designs.ids)
    record_case(input="Sequence → Mock(children=<1..2>)",
                expected=["protA_<1..2>", "protB_<1..2>"], actual=ids)
    assert ids == ["protA_<1..2>", "protB_<1..2>"]
    assert_valid_script(script_path, "Sequence", "Mock")


# ── integration: Sequence + Panda ────────────────────────────────────────────

def test_sequence_feeds_panda_filter_sort(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda filters on the sequences table produced by Sequence."""
    from biopipelines.sequence import Sequence
    from biopipelines.panda import Panda

    pipeline = new_pipeline("seq_feeds_panda")
    with pipeline:
        s = Sequence(seq=["MKTAY", "AETGFLMK", "GG"],
                     ids=["short1", "long1", "short2"], type="protein")
        Panda(
            tables=s.tables.sequences,
            operations=[
                Panda.filter("length >= 5"),
                Panda.sort("id", ascending=True),
            ],
        )
        script_path = pipeline.save()

    record_case(input="Sequence → Panda(filter+sort)",
                expected="script with Sequence+Panda",
                actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Sequence", "Panda", "seq_feeds_panda")


def test_sequence_panda_head(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda.head(n) keeps the first n rows of the sequences table."""
    from biopipelines.sequence import Sequence
    from biopipelines.panda import Panda

    pipeline = new_pipeline("seq_panda_head")
    with pipeline:
        s = Sequence(seq=["MKTAY", "AETGFL", "GGGGA", "HHHHK"],
                     ids=["a", "b", "c", "d"], type="protein")
        Panda(
            tables=s.tables.sequences,
            operations=[Panda.head(2)],
        )
        script_path = pipeline.save()

    record_case(input="Sequence → Panda(head(2))",
                expected="script",
                actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Sequence", "Panda")


def test_sequence_panda_tail(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda.tail(n) keeps the last n rows."""
    from biopipelines.sequence import Sequence
    from biopipelines.panda import Panda

    pipeline = new_pipeline("seq_panda_tail")
    with pipeline:
        s = Sequence(seq=["MKTAY", "AETGFL", "GGGGA", "HHHHK"],
                     ids=["a", "b", "c", "d"], type="protein")
        Panda(
            tables=s.tables.sequences,
            operations=[Panda.tail(2)],
        )
        script_path = pipeline.save()

    record_case(input="Sequence → Panda(tail(2))",
                expected="script", actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Sequence", "Panda")


def test_sequence_panda_sort_descending(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda.sort(by='length', ascending=False) — sort by a numeric column."""
    from biopipelines.sequence import Sequence
    from biopipelines.panda import Panda

    pipeline = new_pipeline("seq_panda_sort_desc")
    with pipeline:
        s = Sequence(seq=["MK", "MKTAY", "MKT"],
                     ids=["a", "b", "c"], type="protein")
        Panda(
            tables=s.tables.sequences,
            operations=[Panda.sort("length", ascending=False)],
        )
        script_path = pipeline.save()

    record_case(input="Sequence → Panda(sort by length desc)",
                expected="script", actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Sequence", "Panda")


def test_sequence_panda_filter_sort_head_chain(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """filter + sort + head in a single Panda call."""
    from biopipelines.sequence import Sequence
    from biopipelines.panda import Panda

    pipeline = new_pipeline("seq_panda_fsh")
    with pipeline:
        s = Sequence(seq=["MK", "MKTAY", "MKTAYLMK", "MKTAYL"],
                     ids=["a", "b", "c", "d"], type="protein")
        Panda(
            tables=s.tables.sequences,
            operations=[
                Panda.filter("length >= 5"),
                Panda.sort("length", ascending=False),
                Panda.head(2),
            ],
        )
        script_path = pipeline.save()

    record_case(input="Sequence → Panda(filter+sort+head)",
                expected="chained ops script",
                actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Sequence", "Panda")


def test_sequence_panda_select_drop_columns(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """select_columns then drop_columns compose cleanly."""
    from biopipelines.sequence import Sequence
    from biopipelines.panda import Panda

    pipeline = new_pipeline("seq_panda_cols")
    with pipeline:
        s = Sequence(seq=["MKTAY", "AETGFL"],
                     ids=["a", "b"], type="protein")
        Panda(
            tables=s.tables.sequences,
            operations=[
                Panda.select_columns(["id", "sequence", "length"]),
                Panda.drop_columns(["length"]),
            ],
        )
        script_path = pipeline.save()

    record_case(input="Sequence → Panda(select+drop columns)",
                expected="script", actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Sequence", "Panda")


def test_sequence_panda_rename(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda.rename remaps a column name."""
    from biopipelines.sequence import Sequence
    from biopipelines.panda import Panda

    pipeline = new_pipeline("seq_panda_rename")
    with pipeline:
        s = Sequence(seq=["MKTAY", "AETGFL"],
                     ids=["a", "b"], type="protein")
        Panda(
            tables=s.tables.sequences,
            operations=[Panda.rename({"length": "n_residues"})],
        )
        script_path = pipeline.save()

    record_case(input="Sequence → Panda(rename length→n_residues)",
                expected="script", actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Sequence", "Panda")


def test_sequence_panda_sample(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda.sample(n, random_state) fixes the subsample size."""
    from biopipelines.sequence import Sequence
    from biopipelines.panda import Panda

    pipeline = new_pipeline("seq_panda_sample")
    with pipeline:
        s = Sequence(seq=["MK", "MKTAY", "MKTAYL", "MKTAYLMK"],
                     ids=["a", "b", "c", "d"], type="protein")
        Panda(
            tables=s.tables.sequences,
            operations=[Panda.sample(n=2, random_state=7)],
        )
        script_path = pipeline.save()

    record_case(input="Sequence → Panda(sample n=2)",
                expected="script", actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Sequence", "Panda")


def test_sequence_panda_drop_duplicates(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda.drop_duplicates collapses repeated rows."""
    from biopipelines.sequence import Sequence
    from biopipelines.panda import Panda

    pipeline = new_pipeline("seq_panda_dedup")
    with pipeline:
        s = Sequence(seq=["MKTAY", "MKTAY", "AETGF"],
                     ids=["a", "b", "c"], type="protein")
        Panda(
            tables=s.tables.sequences,
            operations=[Panda.drop_duplicates(subset="sequence")],
        )
        script_path = pipeline.save()

    record_case(input="Sequence → Panda(drop_duplicates on sequence)",
                expected="script", actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Sequence", "Panda")
