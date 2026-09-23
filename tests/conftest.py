"""Shared pytest fixtures + reporting hooks for the BioPipelines test suite."""

import csv
import os
import re
from pathlib import Path

import pytest


FIXTURES_DIR = Path(__file__).parent / "fixtures"
REPO_ROOT = Path(__file__).resolve().parent.parent
# Write results to tests/_results/ (gitignored). Previously written into
# tests/ directly, where the .csv/.xlsx ended up tracked and every run
# dirtied the working tree.
RESULTS_DIR = Path(__file__).resolve().parent / "_results"
RESULTS_CSV = RESULTS_DIR / "test_results.csv"
RESULTS_XLSX = RESULTS_DIR / "test_results.xlsx"

_CASE_STASH_KEY = pytest.StashKey[dict]()


# ── network tests: an RCSB outage must not read as a regression ──────────────

# One probe per session against a structure that has existed since 1987. Only a
# TRANSPORT failure skips: a 4xx/5xx from a reachable RCSB is a real answer and
# must still fail the suite, or this fixture would hide the thing it guards.
_RCSB_PROBE_URL = "https://data.rcsb.org/rest/v1/core/entry/1UBQ"


@pytest.fixture(scope="session")
def rcsb_reachable():
    """(ok, reason). Retries, so one blip does not skip 25 tests."""
    import time

    try:
        import requests
    except ImportError as exc:  # pragma: no cover - requests is a hard dependency
        return False, f"requests is not importable: {exc}"

    last = ""
    for attempt in range(3):
        try:
            response = requests.get(_RCSB_PROBE_URL, timeout=15)
        except requests.exceptions.RequestException as exc:
            last = f"{type(exc).__name__}: {exc}"
            time.sleep(2 ** attempt)
            continue
        if response.status_code < 500:
            # Reachable. A 4xx here is RCSB answering, which is all this proves.
            return True, ""
        last = f"HTTP {response.status_code}"
        time.sleep(2 ** attempt)
    return False, last


@pytest.fixture(autouse=True)
def _skip_when_rcsb_is_unreachable(request):
    """Skip a `network` test when RCSB cannot be reached at all.

    These tests gate the mirror: a red test stage on GitLab blocks the push to
    the public GitHub repo, so an RCSB outage would stop a release on a tree
    with nothing wrong with it. Skipping is the lesser harm, and `-ra` prints
    the reason so a skipped run is never mistaken for a passing one.
    """
    if request.node.get_closest_marker("network") is None:
        return
    ok, reason = request.getfixturevalue("rcsb_reachable")
    if not ok:
        pytest.skip(f"RCSB unreachable ({reason}); cannot tell a regression from an outage")


# ── default config variant ───────────────────────────────────────────────────

@pytest.fixture(autouse=True)
def _default_config_variant(monkeypatch):
    """Pin the variant the suite has always implicitly run under.

    Auto-detection used to fall back to the alphabetically first config whose ``machine.username`` was blank, which on a clean checkout is ``config.cluster.yaml``. Every test that constructs a tool without one of the ``*_config`` fixtures therefore loaded the shipped cluster config — and would have loaded a different one on a machine that had claimed a variant in its ``.config.<variant>.yaml`` overlay. Detection now refuses rather than guessing, so the dependency is named here instead: same file as before, on every machine. Tests that exercise detection itself override or delete this variable.
    """
    monkeypatch.setenv("BIOPIPELINES_CONFIG_VARIANT", "cluster")


# ── active-pipeline isolation ─────────────────────────────────────────────────

@pytest.fixture(autouse=True)
def _reset_active_pipeline():
    """Ensure the module-level ``_active_pipeline`` context is empty between
    tests. ``on_the_fly=True`` pipelines constructed in a test without a
    ``with`` block (notably shell-safety tests) leak the active context and
    make any subsequent test that instantiates a tool return a Standardized-
    Output or raise 'Cannot nest Pipeline contexts'."""
    try:
        from biopipelines.pipeline import _active_pipeline
    except ImportError:
        yield
        return
    _active_pipeline.set(None)
    yield
    _active_pipeline.set(None)


# ── local_config / isolated_cwd fixtures ──────────────────────────────────────

@pytest.fixture
def local_config(monkeypatch, tmp_path):
    """Point ConfigManager at tests/fixtures/config.local.yaml and reset its
    singleton state so each test gets a fresh load.

    Yields the absolute path of the fixture config file.
    """
    from biopipelines.config_manager import ConfigManager

    config_path = FIXTURES_DIR / "config.local.yaml"
    assert config_path.exists(), f"Missing fixture: {config_path}"

    ConfigManager._instance = None
    ConfigManager._config = None
    ConfigManager._variant = None

    monkeypatch.setattr(
        ConfigManager, "_get_config_path",
        classmethod(lambda cls, variant=None: str(config_path)),
    )

    yield config_path

    ConfigManager._instance = None
    ConfigManager._config = None
    ConfigManager._variant = None


@pytest.fixture
def slurm_local_config(monkeypatch, tmp_path):
    """Like ``local_config`` but pointing at ``config.slurm_local.yaml`` so the
    framework's SLURM script-generation path is exercised. Tests that need to
    inspect ``slurm_batch*.sh`` artifacts use this instead of ``local_config``."""
    from biopipelines.config_manager import ConfigManager

    config_path = FIXTURES_DIR / "config.slurm_local.yaml"
    assert config_path.exists(), f"Missing fixture: {config_path}"

    ConfigManager._instance = None
    ConfigManager._config = None
    ConfigManager._variant = None

    monkeypatch.setattr(
        ConfigManager, "_get_config_path",
        classmethod(lambda cls, variant=None: str(config_path)),
    )

    yield config_path

    ConfigManager._instance = None
    ConfigManager._config = None
    ConfigManager._variant = None


@pytest.fixture
def renderers_config(monkeypatch, tmp_path):
    """Points at ``config.renderers_local.yaml``, the local fixture plus a
    ``renderers:`` block. ``config.local.yaml`` declares none, so a page built under
    it has no viewers at all -- which is what silenced an earlier viewer regression."""
    from biopipelines.config_manager import ConfigManager

    config_path = FIXTURES_DIR / "config.renderers_local.yaml"
    assert config_path.exists(), f"Missing fixture: {config_path}"

    ConfigManager._instance = None
    ConfigManager._config = None
    ConfigManager._variant = None

    monkeypatch.setattr(
        ConfigManager, "_get_config_path",
        classmethod(lambda cls, variant=None: str(config_path)),
    )

    yield config_path

    ConfigManager._instance = None
    ConfigManager._config = None
    ConfigManager._variant = None


@pytest.fixture
def pbs_local_config(monkeypatch, tmp_path):
    """Points at ``config.pbs_local.yaml``, for tests asserting that a
    SLURM-only feature declines to engage on another scheduler."""
    from biopipelines.config_manager import ConfigManager

    config_path = FIXTURES_DIR / "config.pbs_local.yaml"
    assert config_path.exists(), f"Missing fixture: {config_path}"

    ConfigManager._instance = None
    ConfigManager._config = None
    ConfigManager._variant = None

    monkeypatch.setattr(
        ConfigManager, "_get_config_path",
        classmethod(lambda cls, variant=None: str(config_path)),
    )

    yield config_path

    ConfigManager._instance = None
    ConfigManager._config = None
    ConfigManager._variant = None


@pytest.fixture
def isolated_cwd(tmp_path, monkeypatch):
    """Run a test with cwd set to an isolated temp directory."""
    monkeypatch.chdir(tmp_path)
    return tmp_path


@pytest.fixture
def new_pipeline():
    """Factory for a minimal local-config Pipeline used by smoke tests."""
    from biopipelines.pipeline import Pipeline

    def _make(job: str):
        return Pipeline(
            project="TestSuite",
            job=job,
            description=f"Smoke test: {job}",
            on_the_fly=False,
            local_output=True,
            config="local",
        )
    return _make


@pytest.fixture
def slurm_packed_config(monkeypatch, tmp_path):
    """SLURM config declaring exclusive whole-node allocation, so
    ``Parallel(pack=N)`` engages and emits job steps."""
    from biopipelines.config_manager import ConfigManager

    config_path = FIXTURES_DIR / "config.slurm_packed.yaml"
    assert config_path.exists(), f"Missing fixture: {config_path}"

    ConfigManager._instance = None
    ConfigManager._config = None
    ConfigManager._variant = None

    monkeypatch.setattr(
        ConfigManager, "_get_config_path",
        classmethod(lambda cls, variant=None: str(config_path)),
    )

    yield config_path

    ConfigManager._instance = None
    ConfigManager._config = None
    ConfigManager._variant = None


@pytest.fixture
def new_packed_pipeline():
    """Factory for a Pipeline on the packed SLURM fixture."""
    from biopipelines.pipeline import Pipeline

    def _make(job: str):
        return Pipeline(
            project="TestSuite",
            job=job,
            description=f"Smoke test: {job}",
            on_the_fly=False,
            local_output=True,
            config="slurm_packed",
        )
    return _make


@pytest.fixture
def new_slurm_pipeline():
    """Factory for a minimal SLURM-config Pipeline. Needs ``slurm_local_config``
    to be active so ConfigManager loads the SLURM-flavoured fixture."""
    from biopipelines.pipeline import Pipeline

    def _make(job: str):
        return Pipeline(
            project="TestSuite",
            job=job,
            description=f"Smoke test: {job}",
            on_the_fly=False,
            local_output=True,
            config="slurm_local",
        )
    return _make


@pytest.fixture
def config_variant(monkeypatch):
    """Load any ``tests/fixtures/config.<variant>.yaml`` and yield its ConfigManager.

    The per-manager fixtures (mamba/conda/micromamba/venv/container) exist only to exercise environment resolution and container prefixing, which ``config.local.yaml`` cannot reach: it is ``env_manager: pip`` with no ``environments:`` block, so every tool resolves to no environment at all.
    """
    from biopipelines.config_manager import ConfigManager

    def _reset():
        ConfigManager._instance = None
        ConfigManager._config = None
        ConfigManager._variant = None

    def _load(variant: str):
        config_path = FIXTURES_DIR / f"config.{variant}.yaml"
        assert config_path.exists(), f"Missing fixture: {config_path}"
        _reset()
        monkeypatch.setattr(
            ConfigManager, "_get_config_path",
            classmethod(lambda cls, variant=None, _p=str(config_path): _p),
        )
        # Named explicitly: auto-detection now refuses when no config claims the
        # current user, and a fixture variant never does.
        return ConfigManager(variant=variant)

    _reset()
    yield _load
    _reset()


# ── generated-script structure helpers ────────────────────────────────────────

_TOOLS_HEADER_RE = re.compile(r"^#\s*Tools:\s*(.+)$", re.MULTILINE)

# Steps are emitted as <runtime>/NNN_<Tool>[_<?>].sh and invoked by path.
_STEP_PATH_RE = re.compile(r"(?:^|[/\\])(\d{3})_([A-Za-z0-9]+)(?:_[^\s/\\]*)?\.sh")


def declared_tools(content: str) -> list:
    """Tool names listed in the generated script's ``# Tools:`` comment header."""
    match = _TOOLS_HEADER_RE.search(content)
    if not match:
        return []
    return [name.strip() for name in match.group(1).split(",") if name.strip()]


def script_body(content: str) -> str:
    """The script with comment-only lines removed.

    Every tool name appears in the ``# Tools:`` header, so a bare ``name in content`` test is satisfied by that comment even when the script invokes nothing; structural assertions must run against the body.
    """
    return "\n".join(
        line for line in content.splitlines() if not line.lstrip().startswith("#")
    )


def step_invocation_re(tool_name: str):
    """Regex matching an emitted step-script invocation for ``tool_name``."""
    return re.compile(
        r"(?:^|[/\\])\d{3}_" + re.escape(tool_name) + r"(?:_[^\s/\\]*)?\.sh"
    )


def emitted_tool_order(content: str) -> list:
    """Tool names in the order their step scripts are invoked."""
    return [m.group(2) for m in _STEP_PATH_RE.finditer(script_body(content))]


@pytest.fixture
def assert_valid_script():
    """Assert that a saved pipeline.sh exists and actually invokes its tools.

    A marker naming a declared tool is checked against the emitted step invocation rather than a bare substring, and every tool the script declares must be invoked -- otherwise the ``# Tools:`` comment alone satisfies the assertion and a script with zero invocations passes.
    """
    def _check(script_path: str, *markers: str):
        assert os.path.isfile(script_path), f"pipeline.sh missing: {script_path}"
        content = open(script_path, encoding="utf-8").read()
        assert content.startswith("#!/bin/bash"), "missing shebang"
        assert len(content) > 200, "script suspiciously short"

        declared = declared_tools(content)
        assert declared, "generated script declares no tools in its '# Tools:' header"

        body = script_body(content)

        for tool in declared:
            assert step_invocation_re(tool).search(body), (
                f"tool {tool!r} is declared in the '# Tools:' header but is never "
                f"invoked in the script body"
            )

        # The step scripts the body invokes must have been written to disk.
        for path in re.findall(r"(\S*[/\\]\d{3}_[A-Za-z0-9_]+\.sh)", body):
            assert os.path.isfile(path), f"invoked step script does not exist: {path}"

        for marker in markers:
            if marker in declared:
                assert step_invocation_re(marker).search(body), (
                    f"expected tool {marker!r} to be invoked, but the script body "
                    f"has no step invocation for it"
                )
            else:
                assert marker in body, (
                    f"expected marker {marker!r} missing from script body "
                    f"(comment lines excluded)"
                )
    return _check


@pytest.fixture
def tool_order_in_script():
    """Return the real step-invocation order of tools in a generated script."""
    def _order(script_path: str) -> list:
        return emitted_tool_order(open(script_path, encoding="utf-8").read())
    return _order


# ── record_case fixture (adds input/expected/actual to the report) ────────────

@pytest.fixture
def record_case(request):
    """Let a test record its (input, expected, actual) so the XLSX report can
    show what was tested, not just the pass/fail flag.

    Usage inside a test:
        def test_foo(record_case):
            record_case(input="a_<0..1>", expected=["a_0", "a_1"], actual=actual)
            assert actual == ["a_0", "a_1"]

    For parametrized tests, call this once inside the test body.
    """
    def _record(*, input, expected, actual):
        node = request.node
        node.stash[_CASE_STASH_KEY] = {
            "input": _short_repr(input),
            "expected": _short_repr(expected),
            "actual": _short_repr(actual),
            "matched": expected == actual,
        }
    return _record


def _short_repr(value, limit: int = 200) -> str:
    r = repr(value)
    if len(r) > limit:
        return r[: limit - 3] + "..."
    return r


# ── per-test result capture + post-session CSV/XLSX report ────────────────────

_RESULT_ROWS: list[dict] = []


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    """Capture each test's call-phase outcome, duration, and recorded case."""
    outcome = yield
    report = outcome.get_result()

    # Only log one row per test; prefer the call phase, but fall back to a
    # non-passing setup phase so collection/setup errors still show up.
    if report.when != "call" and not (report.when == "setup" and report.outcome != "passed"):
        return

    nodeid = report.nodeid
    parts = nodeid.split("::")
    module = parts[0] if parts else nodeid
    if len(parts) == 3:
        classname, testname = parts[1], parts[2]
    elif len(parts) == 2:
        classname, testname = "", parts[1]
    else:
        classname, testname = "", nodeid

    case = item.stash.get(_CASE_STASH_KEY, None) if hasattr(item, "stash") else None

    _RESULT_ROWS.append({
        "module": module,
        "class": classname,
        "test": testname,
        "outcome": report.outcome,
        "duration_s": round(report.duration, 4),
        "input": case["input"] if case else "",
        "expected": case["expected"] if case else "",
        "actual": case["actual"] if case else "",
        "matched": ("yes" if case["matched"] else "no") if case else "",
        "error_detail": _short_repr(report.longreprtext, 500) if report.outcome == "failed" else "",
        "phase": report.when,
        "nodeid": nodeid,
    })


def pytest_sessionfinish(session, exitstatus):
    """Write collected test outcomes to outputs/test_results.{csv,xlsx}."""
    if not _RESULT_ROWS:
        return

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    columns = [
        "module", "class", "test", "outcome", "duration_s",
        "input", "expected", "actual", "matched",
        "error_detail", "phase", "nodeid",
    ]

    with RESULTS_CSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(_RESULT_ROWS)

    try:
        from openpyxl import Workbook
        from openpyxl.styles import Font, PatternFill, Alignment
    except ImportError:
        return

    wb = Workbook()
    ws = wb.active
    ws.title = "results"
    ws.append(columns)
    header_font = Font(bold=True)
    for cell in ws[1]:
        cell.font = header_font

    outcome_fills = {
        "passed": PatternFill("solid", fgColor="C6EFCE"),
        "failed": PatternFill("solid", fgColor="FFC7CE"),
        "skipped": PatternFill("solid", fgColor="FFEB9C"),
        "error": PatternFill("solid", fgColor="FFC7CE"),
    }
    outcome_col_idx = columns.index("outcome") + 1

    wrap = Alignment(wrap_text=True, vertical="top")
    wrap_cols = {"input", "expected", "actual", "error_detail"}

    for row in _RESULT_ROWS:
        ws.append([row[c] for c in columns])
        r_idx = ws.max_row
        fill = outcome_fills.get(row["outcome"])
        if fill is not None:
            ws.cell(row=r_idx, column=outcome_col_idx).fill = fill
        for i, col in enumerate(columns, start=1):
            if col in wrap_cols:
                ws.cell(row=r_idx, column=i).alignment = wrap

    for i, col in enumerate(columns, start=1):
        max_len = max([len(str(r[col])) for r in _RESULT_ROWS] + [len(col)])
        ws.column_dimensions[ws.cell(row=1, column=i).column_letter].width = min(max_len + 2, 60)

    ws.freeze_panes = "A2"

    # Summary sheet with counts per outcome.
    summary = wb.create_sheet("summary")
    summary.append(["outcome", "count"])
    for cell in summary[1]:
        cell.font = header_font
    counts: dict[str, int] = {}
    for r in _RESULT_ROWS:
        counts[r["outcome"]] = counts.get(r["outcome"], 0) + 1
    for outcome, count in sorted(counts.items()):
        summary.append([outcome, count])
    summary.append([])
    summary.append(["recorded cases", sum(1 for r in _RESULT_ROWS if r["input"])])
    summary.append(["total tests", len(_RESULT_ROWS)])

    wb.save(RESULTS_XLSX)


# --- a run tree built by the framework, not by a fixture's memory ----------------------------
#
# The readers (`job_status`, `lineage`, the MCP tools) are the only tests that consume an output
# tree instead of producing one, so they have no producer to disagree with -- and three of them
# invented a layout that did not exist. `bp_table` shipped listing `<step>/*.csv`, which is
# empty for every real run: standalone tables go to `<step>/tables/` and a stream's map table
# into the stream's own folder. Both the code and its fixture were written from the same wrong
# picture in the same hour, so they agreed and CI stayed green.
#
# This factory takes every path from the framework -- `stream_map_path`, `tables_folder`,
# `_compute_log_file_path` -- so a fixture cannot hold an opinion about the layout. If the
# framework moves a folder, the readers' expectations fail here, loudly, instead of quietly
# agreeing with a stale assumption.

@pytest.fixture
def produced_run(local_config, isolated_cwd):
    """Factory: a real run tree on disk, laid out by the framework. Returns (job_dir, steps).

    `steps` is a list of dicts:
      {"streams": {name: n_ids}, "tables": {name: [row dicts]}, "status": "COMPLETED"|"FAILED"|None}

    Every path written here comes from what the tool *declared* -- `stream.map_table`,
    `TableInfo.path`, `_compute_log_file_path()` -- so the fixture cannot hold an opinion about
    where things go. That is the whole point: the reader tests are the only ones that consume a
    tree rather than produce one, and when they invent the tree they are testing their own
    assumption. `bp_table` shipped listing `<step>/*.csv` -- empty for every real run -- because
    its fixture agreed with it.
    """
    from biopipelines.mock import Mock
    from biopipelines.pipeline import Pipeline

    def _make(steps, job="produced"):
        pipeline = Pipeline(project="TestSuite", job=job, description="fixture",
                            on_the_fly=False, local_output=True, config="local")
        made = []
        with pipeline:
            for spec in steps:
                widest = max(spec.get("streams", {}).values(), default=1)
                ids = [f"d{i}" for i in range(widest)]
                streams = {name: {"format": "pdb", "file": "<id>.pdb"}
                           for name in spec.get("streams", {})}
                tables = {name: {"columns": list(rows[0]) if rows else ["id"]}
                          for name, rows in spec.get("tables", {}).items()}
                made.append((Mock(ids=ids, streams=streams or None, tables=tables or None), spec))
            pipeline.save()

        job_dir = os.path.dirname(made[0][0].output_folder) if made else ""
        produced = {}
        for out, spec in made:
            step = os.path.basename(out.output_folder)
            produced[step] = out

            for name, count in spec.get("streams", {}).items():
                stream = getattr(out.streams, name)
                # These come back as reference proxies; str() is what makes them a path.
                map_table = str(stream.map_table)
                folder = os.path.dirname(map_table)
                os.makedirs(folder, exist_ok=True)
                with open(map_table, "w", encoding="utf-8", newline="") as handle:
                    handle.write("id,file,value\n")
                    for i in range(count):
                        path = os.path.join(folder, f"d{i}.pdb")
                        open(path, "w", encoding="utf-8").write("ATOM\n")
                        handle.write(f"d{i},{path},\n")

            for name, rows in spec.get("tables", {}).items():
                # `TableInfo.path` comes back as a pipeline TableReference, whose str() is a
                # reference token, not a path -- the folder itself is the authority.
                target = os.path.join(out._producer.tables_folder, f"{name}.csv")
                os.makedirs(os.path.dirname(target), exist_ok=True)
                columns = list(rows[0]) if rows else ["id"]
                with open(target, "w", encoding="utf-8", newline="") as handle:
                    handle.write(",".join(columns) + "\n")
                    for row in rows:
                        handle.write(",".join(str(row[c]) for c in columns) + "\n")

            status = spec.get("status", "COMPLETED")
            if status:
                open(os.path.join(job_dir, f"{step}_{status}"), "w", encoding="utf-8").write("")
            log = out._producer._compute_log_file_path()
            os.makedirs(os.path.dirname(log), exist_ok=True)
            open(log, "w", encoding="utf-8").write(f"=== {step} ===\n")
        return job_dir, produced

    return _make
