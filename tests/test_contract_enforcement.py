# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Tests for the contract-enforcement layer.

Covers, in order:

* the severity ladder — per-check env var, global env var, ``SEVERITY`` table, ``DEFAULT_SEVERITY`` — including that an invalid environment value falls through instead of being honoured, and that the global switch outranks the table;
* ``report()`` behaviour per severity, and the dedup cache that keeps a per-id loop from printing hundreds of identical lines;
* each individual check's true and false cases, including the empty-placeholder exemption that keeps "this tool emits none of these" from reading as "value-based";
* ``check_no_unknown_kwargs`` as reached through a real ``BaseConfig.__init__``, i.e. after the subclass has bound its own named parameters;
* the ``tool_context`` plumbing that lets a warning name the offending tool while a ``DataStream`` stays constructible standalone;
* the ``USER_STREAM_NAMES`` opt-out for tools whose stream names come from their user;
* that a ``ContractViolation`` raised inside ``get_output_files()`` reaches the caller instead of being swallowed by ``pipeline.py``'s output-population handlers, while every other exception is still swallowed as before.
"""

import pytest

from biopipelines import contract_enforcement as ce
from biopipelines.contract_enforcement import ContractViolation, Violation


@pytest.fixture(autouse=True)
def _clean_enforcement_state(monkeypatch):
    """The dedup cache and the env-var ladder are process-global."""
    for var in list(__import__("os").environ):
        if var.startswith("BIOPIPELINES_ENFORCE"):
            monkeypatch.delenv(var, raising=False)
    ce.reset_reported()
    yield
    ce.reset_reported()


def _violation(check="stream_name", message="boom"):
    return Violation(check=check, message=message)


# ── severity ladder ───────────────────────────────────────────────────────────

def test_severity_falls_back_to_the_table(record_case):
    """With no environment set, the SEVERITY table decides."""
    record_case(input="severity_for('stream_name') with no env",
                expected=ce.SEVERITY["stream_name"], actual=ce.severity_for("stream_name"))
    assert ce.severity_for("stream_name") == ce.SEVERITY["stream_name"]


def test_severity_of_an_unlisted_check_is_the_default(record_case):
    actual = ce.severity_for("not_a_registered_check")
    record_case(input="severity_for('not_a_registered_check')",
                expected=ce.DEFAULT_SEVERITY, actual=actual)
    assert actual == ce.DEFAULT_SEVERITY


def test_check_specific_env_var_wins(monkeypatch, record_case):
    monkeypatch.setenv("BIOPIPELINES_ENFORCE_STREAM_NAME", "raise")
    actual = ce.severity_for("stream_name")
    record_case(input="ENFORCE_STREAM_NAME=raise", expected="raise", actual=actual)
    assert actual == ce.RAISE


@pytest.mark.parametrize("name", sorted(["forwarded_kwargs", "deprecated_alias"]))
def test_a_row_reporting_a_working_construction_can_still_be_promoted(monkeypatch, name):
    """These report a supported construction rather than a fault, so they ship at warn. Naming one is a deliberate choice about that check and still raises."""
    monkeypatch.setenv(f"BIOPIPELINES_ENFORCE_{name.upper()}", "raise")
    assert ce.severity_for(name) == ce.RAISE


def test_the_global_switch_is_gone(monkeypatch):
    """`BIOPIPELINES_ENFORCE` promised to silence everything and to make everything fatal, and did neither. Setting it now changes nothing."""
    for value in ("raise", "off"):
        monkeypatch.setenv("BIOPIPELINES_ENFORCE", value)
        assert {ce.severity_for(name) for name in ce.SEVERITY} == {ce.WARN}


def test_a_check_specific_truthy_value_is_accepted(monkeypatch):
    """`=1` must not be read as "unrecognized, ignore": the caller believes enforcement is on."""
    monkeypatch.setenv("BIOPIPELINES_ENFORCE_STREAM_NAME", "1")
    assert ce.severity_for("stream_name") == ce.RAISE


def test_an_unrecognized_value_falls_through_and_says_so(monkeypatch, capsys):
    """Falling through is right for a value that means nothing, but doing it in silence is the same failure the aliases exist to prevent, so it announces itself once."""
    ce._warned_unrecognized.clear()
    monkeypatch.setenv("BIOPIPELINES_ENFORCE_COMPOUNDS_FORMAT", "sometimes")
    assert ce.severity_for("compounds_format") == ce.SEVERITY["compounds_format"]
    assert "is not a severity and was ignored" in capsys.readouterr().err


def test_the_unrecognized_notice_is_printed_once_per_variable(monkeypatch, capsys):
    ce._warned_unrecognized.clear()
    monkeypatch.setenv("BIOPIPELINES_ENFORCE_COMPOUNDS_FORMAT", "sometimes")
    ce.severity_for("compounds_format")
    capsys.readouterr()
    ce.severity_for("compounds_format")
    assert capsys.readouterr().err == ""


# ── report() ──────────────────────────────────────────────────────────────────

def test_report_of_none_is_a_pass():
    ce.report(None)  # a check that found nothing must not be a special case


def test_report_off_is_silent(monkeypatch, capsys, record_case):
    monkeypatch.setenv("BIOPIPELINES_ENFORCE_STREAM_NAME", "off")
    ce.report(_violation())
    err = capsys.readouterr().err
    record_case(input="report() at off", expected="", actual=err)
    assert err == ""


def test_report_warn_prints_to_stderr(monkeypatch, capsys, record_case):
    monkeypatch.setenv("BIOPIPELINES_ENFORCE_STREAM_NAME", "warn")
    ce.report(_violation(message="unregistered stream name 'x'"))
    err = capsys.readouterr().err
    record_case(input="report() at warn",
                expected="[contract:stream_name] on stderr", actual=err.strip())
    assert "[contract:stream_name] unregistered stream name 'x'" in err


def test_report_raise_raises_contract_violation(monkeypatch, record_case):
    monkeypatch.setenv("BIOPIPELINES_ENFORCE_STREAM_NAME", "raise")
    with pytest.raises(ContractViolation, match="unregistered"):
        ce.report(_violation(message="unregistered stream name 'x'"))
    record_case(input="report() at raise", expected="ContractViolation", actual="ContractViolation")


def test_hint_is_rendered_under_the_message(record_case):
    rendered = Violation(check="c", message="m", hint="h").render()
    record_case(input="Violation(hint='h').render()",
                expected="[contract:c] m\n    h", actual=rendered)
    assert rendered == "[contract:c] m\n    h"


def test_warnings_are_deduped_by_rendered_message(monkeypatch, capsys, record_case):
    monkeypatch.setenv("BIOPIPELINES_ENFORCE_STREAM_NAME", "warn")
    for _ in range(5):
        ce.report(_violation(message="same"))
    lines = [ln for ln in capsys.readouterr().err.splitlines() if ln.strip()]
    record_case(input="the same violation reported 5x", expected=1, actual=len(lines))
    assert len(lines) == 1


def test_dedup_does_not_hide_a_different_offender(monkeypatch, capsys, record_case):
    monkeypatch.setenv("BIOPIPELINES_ENFORCE_STREAM_NAME", "warn")
    ce.report(_violation(message="stream 'a' (in ToolA)"))
    ce.report(_violation(message="stream 'a' (in ToolB)"))
    lines = [ln for ln in capsys.readouterr().err.splitlines() if ln.strip()]
    record_case(input="same stream name from two tools", expected=2, actual=len(lines))
    assert len(lines) == 2


def test_reset_reported_clears_the_cache(monkeypatch, capsys):
    monkeypatch.setenv("BIOPIPELINES_ENFORCE_STREAM_NAME", "warn")
    ce.report(_violation(message="same"))
    ce.reset_reported()
    ce.report(_violation(message="same"))
    lines = [ln for ln in capsys.readouterr().err.splitlines() if ln.strip()]
    assert len(lines) == 2


def test_check_all_reports_every_violation(monkeypatch, capsys, record_case):
    monkeypatch.setenv("BIOPIPELINES_ENFORCE", "warn")
    ce.check_all(_violation("stream_name", "one"), None, _violation("compounds_format", "two"))
    err = capsys.readouterr().err
    record_case(input="check_all(v1, None, v2)", expected=("one" in err, "two" in err),
                actual=("one" in err, "two" in err))
    assert "one" in err and "two" in err


def test_check_all_raises_on_the_first_raise_severity_violation(monkeypatch):
    monkeypatch.setenv("BIOPIPELINES_ENFORCE_STREAM_NAME", "raise")
    monkeypatch.setenv("BIOPIPELINES_ENFORCE_COMPOUNDS_FORMAT", "raise")
    with pytest.raises(ContractViolation, match="one"):
        ce.check_all(_violation("stream_name", "one"), _violation("compounds_format", "two"))


# ── check_stream_name ─────────────────────────────────────────────────────────

def test_known_stream_name_passes():
    assert ce.check_stream_name("structures") is None


def test_empty_stream_name_passes():
    """An unnamed placeholder stream carries no name to check."""
    assert ce.check_stream_name("") is None


def test_unregistered_stream_name_is_a_violation(record_case):
    v = ce.check_stream_name("strucutres")
    record_case(input="check_stream_name('strucutres')",
                expected="stream_name violation", actual=v.check)
    assert v is not None and v.check == "stream_name"
    assert "strucutres" in v.message


def test_where_is_rendered_into_the_message(record_case):
    v = ce.check_stream_name("strucutres", where="MyTool")
    record_case(input="check_stream_name('strucutres', where='MyTool')",
                expected="(in MyTool)", actual=v.message)
    assert "(in MyTool)" in v.message


# ── check_compounds_format ────────────────────────────────────────────────────

def test_compounds_as_csv_passes():
    assert ce.check_compounds_format("compounds", [], "csv") is None


def test_compounds_format_is_case_insensitive():
    assert ce.check_compounds_format("compounds", [], "CSV") is None


def test_compounds_as_sdf_is_a_violation(record_case):
    v = ce.check_compounds_format("compounds", ["a.sdf"], "sdf")
    record_case(input="compounds stream with format='sdf'",
                expected="compounds_format violation", actual=v.check)
    assert v.check == "compounds_format"


def test_a_non_compounds_stream_may_be_sdf():
    """The contract is about the compounds stream, not about sdf."""
    assert ce.check_compounds_format("structures", ["a.sdf"], "sdf") is None


# ── check_value_based_format ──────────────────────────────────────────────────

def test_value_based_stream_as_csv_passes():
    assert ce.check_value_based_format("binding", [], "csv", ids=["a"], map_table="m.csv") is None


def test_value_based_stream_with_a_non_csv_format_is_a_violation(record_case):
    v = ce.check_value_based_format("binding", [], "pdb", ids=["a"], map_table="m.csv")
    record_case(input="files=[] with format='pdb'",
                expected="value_based_format violation", actual=v.check)
    assert v.check == "value_based_format"


def test_a_stream_with_files_is_not_value_based():
    assert ce.check_value_based_format("structures", ["a.pdb"], "pdb", ids=["a"]) is None


def test_a_shared_file_stream_is_not_value_based():
    """files as a non-empty str is the shared-file form, not the value-based one."""
    assert ce.check_value_based_format("sequences", "all.fasta", "fasta", ids=["a"]) is None


def test_an_empty_placeholder_is_exempt(record_case):
    """No files, no ids, no map_table means 'this tool emits none of these'."""
    actual = ce.check_value_based_format("structures", [], "pdb", ids=[], map_table="")
    record_case(input="DataStream.empty-shaped stream",
                expected=None, actual=actual)
    assert actual is None


def test_an_empty_placeholder_with_ids_is_not_exempt():
    """Ids with no files means the content lives in the map_table, i.e. csv."""
    assert ce.check_value_based_format("structures", [], "pdb", ids=["a"], map_table="") is not None


def test_an_empty_placeholder_with_a_map_table_is_not_exempt():
    assert ce.check_value_based_format("structures", [], "pdb", ids=[], map_table="m.csv") is not None


# ── check_stream: the two format checks never double-report ───────────────────

def test_a_non_csv_compounds_stream_reports_once(monkeypatch, capsys, record_case):
    """A compounds stream is value-based by definition; two lines for one mistake would train readers to ignore them."""
    monkeypatch.setenv("BIOPIPELINES_ENFORCE", "warn")
    ce.check_stream("compounds", [], "sdf", ids=["a"], map_table="m.csv")
    lines = [ln for ln in capsys.readouterr().err.splitlines() if ln.startswith("[contract:")]
    record_case(input="check_stream('compounds', files=[], format='sdf')",
                expected=["compounds_format"], actual=[ln.split("]")[0][10:] for ln in lines])
    assert len(lines) == 1 and "compounds_format" in lines[0]


# ── tool_context ──────────────────────────────────────────────────────────────

def test_no_context_by_default(record_case):
    ctx = ce.current_context()
    record_case(input="current_context() outside any tool",
                expected=("", False), actual=(ctx.where, ctx.user_stream_names))
    assert ctx.where == "" and ctx.user_stream_names is False


def test_tool_context_supplies_where(monkeypatch, capsys, record_case):
    monkeypatch.setenv("BIOPIPELINES_ENFORCE", "warn")
    with ce.tool_context("MyTool"):
        ce.check_stream("strucutres", ["a.pdb"], "pdb", ids=["a"])
    err = capsys.readouterr().err
    record_case(input="check_stream inside tool_context('MyTool')",
                expected="(in MyTool)", actual=err.strip())
    assert "(in MyTool)" in err


def test_an_explicit_where_beats_the_context(monkeypatch, capsys):
    monkeypatch.setenv("BIOPIPELINES_ENFORCE", "warn")
    with ce.tool_context("Outer"):
        ce.check_stream("strucutres", ["a.pdb"], "pdb", ids=["a"], where="Explicit")
    err = capsys.readouterr().err
    assert "(in Explicit)" in err and "Outer" not in err


def test_tool_context_is_restored_on_exit_and_on_error():
    with ce.tool_context("A"):
        assert ce.current_context().where == "A"
        with pytest.raises(RuntimeError):
            with ce.tool_context("B"):
                assert ce.current_context().where == "B"
                raise RuntimeError("boom")
        assert ce.current_context().where == "A"
    assert ce.current_context().where == ""


def test_user_stream_names_skips_only_the_name_check(monkeypatch, capsys, record_case):
    monkeypatch.setenv("BIOPIPELINES_ENFORCE", "warn")
    with ce.tool_context("MyTool", user_stream_names=True):
        ce.check_stream("whatever_the_user_wants", [], "pdb", ids=["a"], map_table="m.csv")
    lines = [ln for ln in capsys.readouterr().err.splitlines() if ln.startswith("[contract:")]
    record_case(input="unregistered name + non-csv value-based stream, user_stream_names=True",
                expected=["value_based_format"], actual=[ln.split("]")[0][10:] for ln in lines])
    assert len(lines) == 1 and "value_based_format" in lines[0]


# ── check_no_unknown_kwargs ───────────────────────────────────────────────────

def test_reserved_kwargs_are_accepted(record_case):
    actual = ce.check_no_unknown_kwargs("MyTool", dict.fromkeys(ce.RESERVED_KWARGS))
    record_case(input="every reserved framework key", expected=None, actual=actual)
    assert actual is None


def test_no_kwargs_is_a_pass():
    assert ce.check_no_unknown_kwargs("MyTool", {}) is None


def test_an_unknown_kwarg_is_a_violation(record_case):
    v = ce.check_no_unknown_kwargs("MyTool", {"num_designs": 64, "name": "x"})
    record_case(input="MyTool(num_designs=64, name='x')",
                expected="unknown_kwargs violation naming num_designs", actual=v.message)
    assert v.check == "unknown_kwargs"
    assert "'num_designs'" in v.message and "'name'" not in v.message
    assert "MyTool" in v.message


def test_several_unknown_kwargs_are_listed_in_one_violation(record_case):
    v = ce.check_no_unknown_kwargs("MyTool", {"b": 1, "a": 2})
    record_case(input="MyTool(b=1, a=2)", expected="'a', 'b'", actual=v.message)
    assert "'a', 'b'" in v.message  # sorted, so the message is stable for dedup


def test_an_exempt_tool_is_not_flagged(monkeypatch, record_case):
    monkeypatch.setattr(ce, "UNKNOWN_KWARGS_EXEMPT_TOOLS", {"Forwarder"})
    actual = ce.check_no_unknown_kwargs("Forwarder", {"anything": 1})
    record_case(input="an exempt tool with an unknown kwarg", expected=None, actual=actual)
    assert actual is None


# ── the check as reached through a real tool constructor ──────────────────────

def _mock_kwargs(**extra):
    return dict(ids=["a"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}}, **extra)


def test_a_tools_own_named_parameters_are_not_flagged(local_config, isolated_cwd, capsys):
    """BaseConfig.__init__ runs after the subclass bound its named parameters, so a correctly spelled call must be silent."""
    from biopipelines.mock import Mock

    Mock(**_mock_kwargs(map_table_strategy="config", missing=["a"]))
    err = capsys.readouterr().err
    assert "unknown_kwargs" not in err


def test_a_typo_in_a_tool_constructor_is_flagged(local_config, isolated_cwd, capsys, record_case):
    from biopipelines.mock import Mock

    Mock(**_mock_kwargs(map_table_strategie="config"))
    err = capsys.readouterr().err
    record_case(input="Mock(map_table_strategie='config')",
                expected="unknown_kwargs naming Mock", actual=err.strip())
    assert "[contract:unknown_kwargs] Mock got unknown constructor parameter(s): 'map_table_strategie'" in err


def test_framework_keys_reach_a_tool_without_a_warning(local_config, isolated_cwd, capsys):
    from biopipelines.mock import Mock

    Mock(**_mock_kwargs(name="job", resources={"memory": "1GB"}, dependencies=[], _internal=True))
    assert "unknown_kwargs" not in capsys.readouterr().err


def test_a_promoted_unknown_kwarg_raises_from_the_constructor(
    monkeypatch, local_config, isolated_cwd,
):
    from biopipelines.mock import Mock

    monkeypatch.setenv("BIOPIPELINES_ENFORCE_UNKNOWN_KWARGS", "raise")
    with pytest.raises(ContractViolation, match="typoo_param"):
        Mock(**_mock_kwargs(typoo_param=1))


def test_get_output_files_names_the_tool_in_a_stream_warning(
    local_config, isolated_cwd, capsys, record_case,
):
    """The whole point of the context: the warning says which tool to fix."""
    from biopipelines.mock import Mock

    tool = Mock(ids=["a"], streams={"compounds": {"format": "sdf", "file": "<id>.sdf"}})
    tool.get_output_files()
    err = capsys.readouterr().err
    record_case(input="Mock(streams={'compounds': {'format': 'sdf'}})",
                expected="(in Mock)", actual=err.strip())
    assert "compounds_format" in err and "(in Mock)" in err


def test_a_user_named_stream_is_not_flagged_by_a_user_stream_names_tool(
    local_config, isolated_cwd, capsys, record_case,
):
    """Mock(streams={...}) and Scripting let the caller name the streams, which no registry can know in advance."""
    from biopipelines.mock import Mock

    assert Mock.USER_STREAM_NAMES is True
    tool = Mock(ids=["a"], streams={"my_own_name": {"format": "pdb", "file": "<id>.pdb"}})
    tool.get_output_files()
    err = capsys.readouterr().err
    record_case(input="Mock(streams={'my_own_name': ...}).get_output_files()",
                expected="no stream_name warning", actual=err.strip() or "no stream_name warning")
    assert "stream_name" not in err


def test_the_same_stream_name_outside_a_user_stream_names_tool_is_flagged(capsys):
    """The opt-out is scoped to the declaring tool, not to the name."""
    from biopipelines.datastream import DataStream

    DataStream(name="my_own_name", ids=["a"], files=["a.pdb"], format="pdb")
    assert "stream_name" in capsys.readouterr().err


def test_a_datastream_is_usable_standalone(capsys, record_case):
    """Pipe scripts import DataStream directly at runtime; with no tool context the check still runs, just without a name to blame."""
    from biopipelines.datastream import DataStream

    stream = DataStream(name="structures", ids=["a"], files=["a.pdb"], format="pdb")
    err = capsys.readouterr().err
    record_case(input="DataStream(name='structures', ...) with no tool context",
                expected=(["a"], ""), actual=(stream.ids, err))
    assert stream.ids == ["a"] and err == ""


# ── pipeline.py must not swallow a ContractViolation ──────────────────────────

def test_a_contract_violation_surfaces_from_tool_construction(
    monkeypatch, local_config, isolated_cwd, record_case,
):
    """Nearly every DataStream is built under one of pipeline.py's bare `except Exception` handlers around get_output_files(). Swallowing a promoted check there would turn it into a warning plus empty outputs — a worse failure than warn, reported from the wrong place."""
    from biopipelines.mock import Mock

    monkeypatch.setenv("BIOPIPELINES_ENFORCE_COMPOUNDS_FORMAT", "raise")
    pipeline = _pipeline("contract_raise")
    with pytest.raises(ContractViolation, match="compounds stream declared"):
        with pipeline:
            Mock(ids=["a"], streams={"compounds": {"format": "sdf", "file": "<id>.sdf"}})
    record_case(input="Mock(compounds as sdf) in a Pipeline, compounds_format=raise",
                expected="ContractViolation", actual="ContractViolation")


def test_every_other_exception_is_still_swallowed(
    monkeypatch, local_config, isolated_cwd, capsys, record_case,
):
    """The handlers exist because a tool's outputs may not be computable until its dependencies resolve. Re-raising ContractViolation must not change that."""
    from biopipelines.mock import Mock

    def _boom(self):
        raise RuntimeError("outputs not ready yet")

    tool = Mock(ids=["a"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
    monkeypatch.setattr(Mock, "get_output_files", _boom)
    _pipeline("contract_swallow")._auto_register(tool)
    out = capsys.readouterr().out
    record_case(input="get_output_files() raising RuntimeError",
                expected="a warning, no exception", actual=out.strip().splitlines()[-1:])
    assert "Warning: Could not immediately populate outputs for Mock" in out


def test_a_contract_violation_is_not_swallowed_at_the_same_handler(
    monkeypatch, local_config, isolated_cwd,
):
    """Same handler, same call, different exception type."""
    from biopipelines.mock import Mock

    def _boom(self):
        raise ContractViolation("[contract:compounds_format] promoted")

    tool = Mock(ids=["a"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
    monkeypatch.setattr(Mock, "get_output_files", _boom)
    with pytest.raises(ContractViolation, match="promoted"):
        _pipeline("contract_reraise")._auto_register(tool)


def test_every_get_output_files_handler_reraises_contract_violation(record_case):
    """Structural guard: a handler added later without the re-raise would silently re-break promotion, and only an end-to-end run would notice."""
    import ast
    import pathlib

    source = pathlib.Path(__import__("biopipelines.pipeline", fromlist=["x"]).__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"), str(source))

    unguarded = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Try):
            continue
        calls = [n for n in ast.walk(ast.Module(body=node.body, type_ignores=[]))
                 if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)
                 and n.func.attr == "get_output_files"]
        if not calls:
            continue
        names = [h.type.id for h in node.handlers
                 if isinstance(h.type, ast.Name)]
        if "ContractViolation" not in names:
            unguarded.append(node.lineno)

    record_case(input="every try/except around get_output_files() in pipeline.py",
                expected=[], actual=unguarded)
    assert unguarded == [], (
        f"pipeline.py lines {unguarded} swallow a ContractViolation raised while "
        "building a DataStream")


# ── DataStream indexing must not fabricate ids ────────────────────────────────
# Home is tests/test_datastream.py, which belongs to another wave while datastream.py is this one's.

def test_indexing_a_lazy_stream_keeps_the_id_symbolic(record_case):
    """`DataStream.__getitem__` used to guard on has_patterns(), which is True for a lazy id, so expand_at() substituted the in-bracket slot and produced 'prot[_N]' — neither a pattern nor a real id, matching no map_table row and interpolating into a path with a literal bracket."""
    from biopipelines.datastream import DataStream

    stream = DataStream(name="structures", ids=["prot[_<?>]"],
                        files=["prot[_<?>].pdb"], format="pdb")
    actual = stream[0].ids
    record_case(input="DataStream(ids=['prot[_<?>]'])[0].ids",
                expected=["prot[_<?>]"], actual=actual)
    assert actual == ["prot[_<?>]"]


def test_indexing_a_deterministic_pattern_still_expands(record_case):
    """The O(1) expand_at path is the point of the guard; keep it working."""
    from biopipelines.datastream import DataStream

    stream = DataStream(name="structures", ids=["p_<0..3>"],
                        files=["p_<0..3>.pdb"], format="pdb")
    actual = (len(stream), stream[2].ids)
    record_case(input="DataStream(ids=['p_<0..3>'])[2].ids", expected=(4, ["p_2"]), actual=actual)
    assert actual == (4, ["p_2"])


def _pipeline(job):
    from biopipelines.pipeline import Pipeline

    return Pipeline(project="TestSuite", job=job, description=f"Contract: {job}",
                    on_the_fly=False, local_output=True, config="local")
