"""Unit tests for the Table input-type tool.

Table is a passive loader: it reads a CSV or Excel file at config time, records
columns+row count, and exposes a TableInfo that downstream tools reference.
Tests cover the CSV path, the Excel-to-CSV conversion, column-reference access,
and integration with Panda (primary consumer) and a Mock (co-resident in a
pipeline so the combined script still emits correctly).
"""

import os

import pandas as pd
import pytest


# ── CSV loading ──────────────────────────────────────────────────────────────

def test_table_loads_csv(isolated_cwd, record_case):
    """Raw-attribute inspection — construct outside Pipeline context."""
    from biopipelines.table import Table

    csv_path = isolated_cwd / "metrics.csv"
    csv_path.write_text(
        "id,score,rmsd\n"
        "a,0.9,1.1\n"
        "b,0.7,1.5\n"
        "c,0.8,0.9\n"
    )

    t = Table(str(csv_path))
    record_case(input="Table('metrics.csv')",
                expected=(["id", "score", "rmsd"], 3, "data"),
                actual=(list(t.table_columns), t.table_count, t.table_name))
    assert t.table_columns == ["id", "score", "rmsd"]
    assert t.table_count == 3
    assert t.table_name == "data"


def test_table_script_emitted(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Smoke: Table inside a Pipeline produces a valid script."""
    from biopipelines.table import Table

    csv_path = isolated_cwd / "metrics.csv"
    csv_path.write_text("id,score\na,0.9\nb,0.7\n")

    pipeline = new_pipeline("table_script")
    with pipeline:
        Table(str(csv_path))
        script_path = pipeline.save()

    record_case(input="Table in pipeline",
                expected="script with Table",
                actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Table")


def test_table_custom_name_and_description(isolated_cwd, record_case):
    from biopipelines.table import Table

    csv_path = isolated_cwd / "scores.csv"
    csv_path.write_text("id,plddt\na,0.92\nb,0.85\n")

    t = Table(str(csv_path), name="metrics",
              description="pLDDT from AlphaFold")
    record_case(input="Table(..., name='metrics', description=...)",
                expected=("metrics", "pLDDT from AlphaFold"),
                actual=(t.table_name, t.table_description))
    assert t.table_name == "metrics"
    assert t.table_description == "pLDDT from AlphaFold"


def test_table_exposes_tableinfo_via_standardized_output(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """Inside a Pipeline, Table(...) returns a StandardizedOutput exposing
    `.tables.<name>` with the expected TableInfo."""
    from biopipelines.table import Table
    from biopipelines.base_config import StandardizedOutput

    csv_path = isolated_cwd / "scores.csv"
    csv_path.write_text("id,score\na,0.9\nb,0.8\n")

    pipeline = new_pipeline("table_stdout")
    with pipeline:
        t = Table(str(csv_path), name="scores")
        pipeline.save()

    record_case(input="Table(...).tables.scores (StandardizedOutput access)",
                expected=("scores", ["id", "score"]),
                actual=(t.tables.scores.info.name,
                        list(t.tables.scores.info.columns)))
    assert isinstance(t, StandardizedOutput)
    assert t.tables.scores.info.name == "scores"
    assert list(t.tables.scores.info.columns) == ["id", "score"]


def test_table_column_reference_string(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    """tool.tables.x.col yields a TABLE_REFERENCE string (same contract as Mock)."""
    from biopipelines.table import Table

    csv_path = isolated_cwd / "scores.csv"
    csv_path.write_text("id,score\na,0.9\nb,0.8\n")

    pipeline = new_pipeline("table_colref")
    with pipeline:
        t = Table(str(csv_path), name="scores")
        pipeline.save()

    ref = t.tables.scores.score
    ref_str = str(ref)
    record_case(input="Table(...).tables.scores.score",
                expected=("TABLE_REFERENCE:", ":score"),
                actual=(ref_str[:16], ref_str[-6:]))
    assert ref_str.startswith("TABLE_REFERENCE:")
    assert ref_str.endswith(":score")
    assert ref.column == "score"


# ── Excel loading + conversion ───────────────────────────────────────────────

def test_table_loads_excel_and_converts_to_csv(isolated_cwd, record_case):
    """Excel (.xlsx) is converted to a sibling .csv; table_path points to the CSV.

    openpyxl is a core dependency — this test runs unconditionally.
    """
    from biopipelines.table import Table

    xlsx_path = isolated_cwd / "scores.xlsx"
    df = pd.DataFrame({"id": ["a", "b"], "score": [0.9, 0.8]})
    df.to_excel(xlsx_path, index=False)

    t = Table(str(xlsx_path))
    csv_sibling = str(xlsx_path).replace(".xlsx", ".csv")
    record_case(input="Table('scores.xlsx')",
                expected=(True, ["id", "score"], 2),
                actual=(os.path.exists(csv_sibling),
                        list(t.table_columns), t.table_count))
    assert os.path.exists(csv_sibling)
    assert t.table_path == csv_sibling
    assert t.table_columns == ["id", "score"]


# ── validation errors ───────────────────────────────────────────────────────

def test_table_missing_file_raises(local_config, isolated_cwd, new_pipeline, record_case):
    from biopipelines.table import Table

    # Table defers existence checks to configure_inputs (so it can resolve
    # bare filenames against the pipeline's tables/ folder). The error
    # surfaces when the Table is registered into a Pipeline, not at
    # construction time — drive it through a Pipeline to trigger.
    record_case(input="Table('absent.csv') inside Pipeline",
                expected="FileNotFoundError", actual="FileNotFoundError")
    with pytest.raises(FileNotFoundError):
        with new_pipeline("table_missing"):
            Table("absent.csv")


def test_table_rejects_unsupported_extension(isolated_cwd, record_case):
    from biopipelines.table import Table

    bad_path = isolated_cwd / "notes.txt"
    bad_path.write_text("id,score\na,0.9\n")

    record_case(input="Table('notes.txt') — unsupported",
                expected="ValueError", actual="ValueError")
    with pytest.raises(ValueError, match="CSV or Excel"):
        Table(str(bad_path))


# ── integration: Table + Panda ───────────────────────────────────────────────

def test_table_feeds_panda_filter_sort(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda filters and sorts a Table loaded from a CSV."""
    from biopipelines.table import Table
    from biopipelines.panda import Panda

    csv_path = isolated_cwd / "metrics.csv"
    csv_path.write_text(
        "id,score,rmsd\n"
        "a,0.9,1.1\n"
        "b,0.7,1.5\n"
        "c,0.8,0.9\n"
        "d,0.4,2.0\n"
    )

    pipeline = new_pipeline("table_feeds_panda")
    with pipeline:
        t = Table(str(csv_path), name="metrics")
        Panda(
            tables=t.tables.metrics,
            operations=[
                Panda.filter("score > 0.6"),
                Panda.sort("score", ascending=False),
            ],
        )
        script_path = pipeline.save()

    record_case(input="Table → Panda(filter+sort)",
                expected="script with Table+Panda",
                actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Table", "Panda", "table_feeds_panda")


def test_table_multi_table_panda_merge(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Two Table inputs merged by Panda.merge — exercises multi-table wiring."""
    from biopipelines.table import Table
    from biopipelines.panda import Panda

    a_path = isolated_cwd / "scoresA.csv"
    b_path = isolated_cwd / "scoresB.csv"
    a_path.write_text("id,plddt\nx,0.9\ny,0.8\n")
    b_path.write_text("id,rmsd\nx,1.1\ny,1.4\n")

    pipeline = new_pipeline("table_multi_panda")
    with pipeline:
        ta = Table(str(a_path), name="a")
        tb = Table(str(b_path), name="b")
        Panda(
            tables=[ta.tables.a, tb.tables.b],
            operations=[Panda.merge(prefixes=["A_", "B_"])],
        )
        script_path = pipeline.save()

    record_case(input="Table+Table → Panda.merge",
                expected="merge script",
                actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Table", "Panda")


# ── integration: Table + Mock (co-resident) ─────────────────────────────────

def test_table_coresident_with_mock(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Table and Mock in the same pipeline both land in the generated script.

    Table has no DataStream, so it can't source a Mock directly — but both tools
    should register cleanly side by side and the combined script should emit.
    """
    from biopipelines.table import Table
    from biopipelines.mock import Mock

    csv_path = isolated_cwd / "metrics.csv"
    csv_path.write_text("id,score\na,0.9\nb,0.8\n")

    pipeline = new_pipeline("table_and_mock")
    with pipeline:
        Table(str(csv_path), name="metrics")
        Mock(
            ids=["a", "b"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
        )
        script_path = pipeline.save()

    record_case(input="Table + Mock (co-resident)",
                expected="script with Table+Mock",
                actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Table", "Mock", "table_and_mock")


def test_table_plus_mock_feeds_panda(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda merges a Table-loaded CSV with Mock's per-ID scores table."""
    from biopipelines.table import Table
    from biopipelines.mock import Mock
    from biopipelines.panda import Panda

    csv_path = isolated_cwd / "ext_metrics.csv"
    csv_path.write_text("id,external\na,0.1\nb,0.2\n")

    pipeline = new_pipeline("table_mock_panda")
    with pipeline:
        t = Table(str(csv_path), name="ext")
        m = Mock(
            ids=["a", "b"],
            streams={"s": {"format": "pdb", "file": "<id>.pdb"}},
            tables={"scores": {"columns": ["score"], "fill": {"score": 0.5}}},
        )
        Panda(
            tables=[m.tables.scores, t.tables.ext],
            operations=[Panda.merge(prefixes=["mock_", "ext_"])],
        )
        script_path = pipeline.save()

    record_case(input="Mock.scores + Table.ext → Panda.merge",
                expected="script",
                actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Table", "Mock", "Panda",
                        "table_mock_panda")


def test_table_panda_head(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda.head on a Table."""
    from biopipelines.table import Table
    from biopipelines.panda import Panda

    csv_path = isolated_cwd / "m.csv"
    csv_path.write_text(
        "id,score\na,0.9\nb,0.7\nc,0.8\nd,0.4\n"
    )

    pipeline = new_pipeline("table_panda_head")
    with pipeline:
        t = Table(str(csv_path))
        Panda(tables=t.tables.data, operations=[Panda.head(2)])
        script_path = pipeline.save()

    record_case(input="Table → Panda(head(2))",
                expected="script", actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Table", "Panda")


def test_table_panda_tail(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda.tail on a Table."""
    from biopipelines.table import Table
    from biopipelines.panda import Panda

    csv_path = isolated_cwd / "m.csv"
    csv_path.write_text(
        "id,score\na,0.9\nb,0.7\nc,0.8\nd,0.4\n"
    )

    pipeline = new_pipeline("table_panda_tail")
    with pipeline:
        t = Table(str(csv_path))
        Panda(tables=t.tables.data, operations=[Panda.tail(1)])
        script_path = pipeline.save()

    record_case(input="Table → Panda(tail(1))",
                expected="script", actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Table", "Panda")


def test_table_panda_filter_sort_head_chain(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """filter + sort + head chain on a Table."""
    from biopipelines.table import Table
    from biopipelines.panda import Panda

    csv_path = isolated_cwd / "m.csv"
    csv_path.write_text(
        "id,score\na,0.9\nb,0.3\nc,0.8\nd,0.5\ne,0.7\n"
    )

    pipeline = new_pipeline("table_panda_fsh")
    with pipeline:
        t = Table(str(csv_path))
        Panda(
            tables=t.tables.data,
            operations=[
                Panda.filter("score >= 0.5"),
                Panda.sort("score", ascending=False),
                Panda.head(2),
            ],
        )
        script_path = pipeline.save()

    record_case(input="Table → Panda(filter+sort+head)",
                expected="chained script",
                actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Table", "Panda")


def test_table_panda_select_and_drop_columns(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """select_columns + drop_columns compose over a Table."""
    from biopipelines.table import Table
    from biopipelines.panda import Panda

    csv_path = isolated_cwd / "m.csv"
    csv_path.write_text(
        "id,score,rmsd,note\na,0.9,1.1,first\nb,0.7,1.5,second\n"
    )

    pipeline = new_pipeline("table_panda_cols")
    with pipeline:
        t = Table(str(csv_path))
        Panda(
            tables=t.tables.data,
            operations=[
                Panda.select_columns(["id", "score", "rmsd"]),
                Panda.drop_columns(["rmsd"]),
            ],
        )
        script_path = pipeline.save()

    record_case(input="Table → Panda(select_columns + drop_columns)",
                expected="script", actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Table", "Panda")


def test_table_panda_rename(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """rename remaps score → plddt on a Table."""
    from biopipelines.table import Table
    from biopipelines.panda import Panda

    csv_path = isolated_cwd / "m.csv"
    csv_path.write_text("id,score\na,0.9\nb,0.7\n")

    pipeline = new_pipeline("table_panda_rename")
    with pipeline:
        t = Table(str(csv_path))
        Panda(
            tables=t.tables.data,
            operations=[Panda.rename({"score": "plddt"})],
        )
        script_path = pipeline.save()

    record_case(input="Table → Panda(rename score→plddt)",
                expected="script", actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Table", "Panda")


def test_table_panda_sample(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda.sample on a Table."""
    from biopipelines.table import Table
    from biopipelines.panda import Panda

    csv_path = isolated_cwd / "m.csv"
    csv_path.write_text(
        "id,score\na,0.9\nb,0.7\nc,0.8\nd,0.4\n"
    )

    pipeline = new_pipeline("table_panda_sample")
    with pipeline:
        t = Table(str(csv_path))
        Panda(
            tables=t.tables.data,
            operations=[Panda.sample(n=2, random_state=3)],
        )
        script_path = pipeline.save()

    record_case(input="Table → Panda(sample n=2)",
                expected="script", actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Table", "Panda")


def test_table_panda_drop_duplicates(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """Panda.drop_duplicates on a subset column."""
    from biopipelines.table import Table
    from biopipelines.panda import Panda

    csv_path = isolated_cwd / "m.csv"
    csv_path.write_text(
        "id,score\na,0.9\nb,0.9\nc,0.7\n"
    )

    pipeline = new_pipeline("table_panda_dedup")
    with pipeline:
        t = Table(str(csv_path))
        Panda(
            tables=t.tables.data,
            operations=[Panda.drop_duplicates(subset=["score"])],
        )
        script_path = pipeline.save()

    record_case(input="Table → Panda(drop_duplicates on score)",
                expected="script", actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Table", "Panda")


def test_table_panda_concat_three(
    local_config, isolated_cwd, new_pipeline, assert_valid_script, record_case,
):
    """concat three Table tables into one."""
    from biopipelines.table import Table
    from biopipelines.panda import Panda

    a = isolated_cwd / "a.csv"; a.write_text("id,v\na,1\n")
    b = isolated_cwd / "b.csv"; b.write_text("id,v\nb,2\n")
    c = isolated_cwd / "c.csv"; c.write_text("id,v\nc,3\n")

    pipeline = new_pipeline("table_panda_concat")
    with pipeline:
        ta = Table(str(a), name="a")
        tb = Table(str(b), name="b")
        tc = Table(str(c), name="c")
        Panda(
            tables=[ta.tables.a, tb.tables.b, tc.tables.c],
            operations=[Panda.concat()],
        )
        script_path = pipeline.save()

    record_case(input="3×Table → Panda.concat",
                expected="script", actual=os.path.basename(script_path))
    assert_valid_script(script_path, "Table", "Panda")
