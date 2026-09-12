"""Every tool's documented output schema must match the schema its code declares.

`docs/tool/*.md` tells the reader which streams a tool emits, which tables it writes and
what columns each table has. None of that was checked by anything, and an audit found only
16 of 77 tool sections correct: `missing` tables documented with three of their four
canonical columns, whole tables and streams omitted, columns elided behind `...`, a `PyMOL`
`images` stream that is really called `renders`, and a `ThermoMPNN` `profile` stream with no
implementation at all.

The ground truth is `get_output_files()`. `versions/extract_output_schemas.py` reads it with
an abstract interpreter over the AST; this module reads the docs back and compares the two.

The contract each doc section is held to:

  * A ``**Streams**`` block lists exactly the stream names the code returns as top-level
    `DataStream` keys. `DataStream.empty(...)` placeholders are *not* streams (`Reduce` and
    `UniProt` return empty `sequences`/`compounds` so downstream attribute access works, and
    neither documents them), so listing one is an error.
  * A ``**Tables**`` block lists exactly the table names the code declares, and each table's
    column row matches its `TableInfo(columns=...)` in order.
  * A column whose *name* comes from a parameter (`Distance(metric_name=...)`) is written
    `{parameter_name}` in the docs and matched positionally.
  * Where the code's column list is genuinely open-ended (extended by a comprehension or
    grown in a loop, e.g. `RDKit`'s per-descriptor columns) the docs must contain every
    column that *is* statically known, and may list more.
  * Where a tool computes its stream or table *names* from its inputs (`Panda`, `Pool`,
    `ReMap`, `Load`, `Table`, `Scripting`, `ExtractMetrics`, `Consensus`, `Plot`), only the
    statically-declared names are enforced; the rest can only be described in prose.

Exemptions are explicit and named, the way `test_registry_consistency.py` does it, so adding
one is a decision rather than an oversight.
"""

import io
import pathlib
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "versions"))

import extract_output_schemas as schemas  # noqa: E402  (needs the sys.path entry above)


@pytest.fixture(scope="module")
def code_schemas():
    """{TOOL_NAME: schema} read out of biopipelines/*.py."""
    found = schemas.extract_all()
    assert len(found) > 50, f"only found {len(found)} tools; the extractor is broken"
    return found


@pytest.fixture(scope="module")
def doc_sections():
    """{section title: parsed **Streams** / **Tables** blocks} read out of docs/tool/*.md."""
    found = schemas.parse_doc_sections()
    assert len(found) > 50, f"only found {len(found)} doc sections; the parser is broken"
    return found


def _ids(mapping):
    return sorted(mapping)


@pytest.mark.parametrize("tool", _ids(schemas.extract_all()))
def test_documented_schema_matches_code(tool, code_schemas, doc_sections):
    """The tool's **Streams** and **Tables** blocks agree with its get_output_files()."""
    problems = schemas.compare(tool, code_schemas[tool], doc_sections.get(tool))
    if problems:
        pytest.fail(
            f"{tool}: docs/tool/ disagrees with {code_schemas[tool]['source']}.\n"
            "Ground truth is the code; run `python versions/extract_output_schemas.py "
            f"{tool}` to see what it declares.\n\n" + "\n".join("  - " + p for p in problems))


def test_every_tool_has_a_doc_section(code_schemas, doc_sections):
    """A public tool with outputs must have a docs/tool/ section to document them in."""
    undocumented = sorted(
        tool for tool, schema in code_schemas.items()
        if tool not in doc_sections and (schema["streams"] or schema["tables"]))
    assert not undocumented, (
        "no docs/tool/*.md section documents: " + ", ".join(undocumented)
        + " (add one, or add the tool to extract_output_schemas.INTERNAL if it is not "
          "part of the public API)")


def test_no_doc_section_invents_a_tool(code_schemas, doc_sections):
    """A section carrying a schema block must correspond to a real tool."""
    orphans = sorted(
        title for title, doc in doc_sections.items()
        if title not in code_schemas and title not in schemas.DOC_NO_SCHEMA
        and (doc["has_streams_block"] or doc["has_tables_block"]))
    assert not orphans, (
        "docs/tool/ sections declare streams/tables but match no TOOL_NAME: "
        + ", ".join(orphans)
        + " (rename the heading, add an alias to extract_output_schemas.DOC_ALIASES, or "
          "list it in DOC_NO_SCHEMA if it documents a variant of another tool)")


def test_missing_table_schema_is_canonical(code_schemas):
    """`missing` has one schema, fixed in BaseConfig.missing_table_info().

    Doc drift here was the single most common defect (45 tools declare the table; the docs
    got it right 7 times), so pin the code side too: a tool that spells its own `missing`
    columns out by hand must spell out the same four.
    """
    wrong = {
        tool: entry["columns"]
        for tool, schema in code_schemas.items()
        for name, entry in schema["tables"].items()
        if name == "missing" and entry["variants"]
        and tuple(entry["columns"]) != schemas.MISSING_COLUMNS
    }
    # RCSB reuses its `failed` schema for `missing` rather than the canonical one. That is a code-side inconsistency, not a docs one; the docs match what it declares today.
    assert wrong == {"RCSB": ["pdb_id", "error_message", "source", "attempted_path"]}, (
        "unexpected `missing` schemas (canonical is "
        + " | ".join(schemas.MISSING_COLUMNS) + "): "
        + "; ".join(f"{t}: {' | '.join(c)}" for t, c in sorted(wrong.items())))


def test_doc_sections_declare_a_block_when_the_code_has_output(code_schemas, doc_sections):
    """A tool with streams or tables must actually carry the corresponding block."""
    gaps = []
    for tool, schema in code_schemas.items():
        doc = doc_sections.get(tool)
        if doc is None:
            continue
        if schema["streams"] and not doc["has_streams_block"]:
            gaps.append(f"{tool}: no **Streams** block")
        if schema["tables"] and not doc["has_tables_block"]:
            gaps.append(f"{tool}: no **Tables** block")
    assert not gaps, "docs/tool/ sections missing a schema block:\n  " + "\n  ".join(gaps)


def test_extractor_resolves_the_schemas_it_claims_to(code_schemas):
    """Guard the extractor itself: only the known-dynamic tools may be unresolved.

    Without this, a refactor that made `get_output_files()` unreadable to the interpreter
    would silently turn every check above into a no-op.
    """
    # Tools that build stream or table names from their inputs, so the names cannot be known from the source alone. Each is documented in prose instead.
    DYNAMIC = {
        "Consensus": "the output stream takes the input stream's name",
        "ExtractMetrics": "one table per requested metric",
        "Load": "inherits the loaded tool's whole shape",
        "Panda": "pool mode mirrors the pool's streams; `result` columns follow the op chain",
        "Plot": "one table per plot operation",
        "Pool": "mirrors the shared streams/tables of its runs",
        "ReMap": "mirrors the source streams/tables under new ids",
        "Scripting": "whatever the user's configuration() declares",
        "Table": "the user supplies the columns",
    }
    unexpected = sorted(
        tool for tool, schema in code_schemas.items()
        if (schema["dynamic_streams"] or schema["dynamic_tables"]) and tool not in DYNAMIC)
    assert not unexpected, (
        "these tools' stream/table names no longer resolve statically: "
        + ", ".join(unexpected)
        + " — either the code moved the declaration out of reach of "
          "versions/extract_output_schemas.py (teach the interpreter), or the tool really is "
          "input-driven now (add it to DYNAMIC with the reason)")

    still_static = sorted(t for t in DYNAMIC if t in code_schemas
                          and not code_schemas[t]["dynamic_streams"]
                          and not code_schemas[t]["dynamic_tables"])
    assert not still_static, (
        "these tools now declare their names statically and should be checked in full: "
        + ", ".join(still_static) + " (drop them from DYNAMIC)")


def test_repo_wide_check_is_clean(code_schemas, doc_sections):
    """The same comparison the CLI runs, as one summary assertion."""
    problems, orphans = schemas.compare_all()
    total = sum(len(v) for v in problems.values()) + len(orphans)
    assert total == 0, (
        f"{total} doc/code schema mismatch(es) across {len(problems)} tool(s). "
        "Run `python versions/extract_output_schemas.py --check` for the full list.\n"
        + "\n".join("  - " + p for tool in sorted(problems) for p in problems[tool]))


def test_docs_and_code_are_both_actually_being_read():
    """A smoke test for the two halves, so a silent parse failure cannot pass everything."""
    aggrescan = schemas.extract_tool_schema("Aggrescan3D")
    assert set(aggrescan["streams"]) == {"structures", "images", "aggregation"}
    assert aggrescan["tables"]["aggregation_all"]["columns"] == ["id", "chain", "resi", "score"]

    docs = schemas.parse_doc_sections()
    assert docs["Aggrescan3D"]["streams"] == ["aggregation", "structures", "images"]
    assert docs["Aggrescan3D"]["tables"]["aggregation_all"]["rows"] == [
        ["id", "chain", "resi", "score"]]


def test_docs_tool_files_are_all_scanned():
    """Every docs/tool/*.md file contributes sections, so none is silently skipped."""
    seen = {doc["file"] for doc in schemas.parse_doc_sections().values()}
    on_disk = {f"docs/tool/{p.name}" for p in sorted((ROOT / "docs" / "tool").glob("*.md"))}
    assert on_disk and on_disk <= seen, f"never parsed: {sorted(on_disk - seen)}"


def test_extractor_columns_are_not_silently_empty(code_schemas):
    """A declared table must resolve to some column list, or be a known-dynamic one."""
    UNRESOLVED = {
        ("Panda", "result"): "columns are whatever the operation chain leaves",
    }
    empty = sorted(
        (tool, name) for tool, schema in code_schemas.items()
        for name, entry in schema["tables"].items()
        if not [c for c in entry["columns"] if c is not None]
        and (tool, name) not in UNRESOLVED)
    assert not empty, (
        "these tables resolved to no column names at all: "
        + ", ".join(f"{t}.{n}" for t, n in empty)
        + " (teach versions/extract_output_schemas.py to follow the declaration, or record "
          "it in UNRESOLVED with the reason)")


def test_thermompnn_has_no_profile_stream():
    """ThermoMPNN's docs described a `sequences=`/`profile` feature that does not exist.

    It was removed from the docs rather than implemented; this pins the removal so the
    prose cannot creep back without the code. Delete this test when the feature lands.
    """
    source = io.open(ROOT / "biopipelines" / "thermompnn.py", encoding="utf-8").read()
    assert "sequences" not in source and "profile" not in source, (
        "thermompnn.py now mentions sequences/profile: if the resi-csv profile output was "
        "implemented, document it and drop this test")
    section = schemas.parse_doc_sections()["ThermoMPNN"]
    assert "profile" not in section["streams"], (
        "docs/tool/analysis.md documents a ThermoMPNN `profile` stream again, but "
        "biopipelines/thermompnn.py still has no implementation")
