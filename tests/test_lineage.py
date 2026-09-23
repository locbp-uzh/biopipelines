"""Where the designs went — the counts behind "I started with 500 and have 12".

Derived from tables the framework already writes: `<stream>_map.csv` for the IDs a step
produced, `missing.csv` for the ones that entered and did not come out. The arithmetic is the
whole value, so these tests pin it rather than the formatting.
"""

import io

import pytest

from biopipelines import lineage
from biopipelines.remote import LocalFS


def build(job, steps):
    """steps: {folder: {"streams": {name: rows}, "missing": rows}}

    Writes the layout the framework really produces — a stream's map table inside the stream's
    own folder, standalone tables under `tables/`. An earlier version of this helper put both
    beside the step, and `lineage.tables` was written against that shape and so reported a real
    run as having written nothing.
    """
    for folder, spec in steps.items():
        (job / folder).mkdir(parents=True, exist_ok=True)
        for name, rows in spec.get("streams", {}).items():
            (job / folder / name).mkdir(parents=True, exist_ok=True)
            body = "id,file,value\n" + "".join(f"d{i},f{i}.pdb,\n" for i in range(rows))
            io.open(job / folder / name / f"{name}_map.csv", "w", encoding="utf-8").write(body)
        if spec.get("missing"):
            (job / folder / "tables").mkdir(parents=True, exist_ok=True)
            # Ids distinct per step: the same id in two steps' tables is a propagated drop, counted once.
            body = "id,structure\n" + "".join(
                f"{folder}_d{i},f{i}.pdb\n" for i in range(spec["missing"]))
            io.open(job / folder / "tables" / "missing.csv", "w", encoding="utf-8").write(body)
    return job


def test_counts_exclude_the_csv_header(tmp_path):
    job = build(tmp_path / "run_001", {"001_RFdiffusion": {"streams": {"structures": 500}}})
    assert lineage.collect(job, fs=LocalFS())[0]["streams"]["structures"] == 500


def test_steps_are_ordered_by_index_not_alphabetically(tmp_path):
    job = build(tmp_path / "run_001", {
        "010_Boltz2": {"streams": {"structures": 2}},
        "002_ProteinMPNN": {"streams": {"sequences": 3}}})
    assert [s["folder"] for s in lineage.collect(job, fs=LocalFS())] == [
        "002_ProteinMPNN", "010_Boltz2"]


def test_a_step_with_neither_table_is_omitted(tmp_path):
    """Reporting zero for a tool that writes no stream would read as total attrition."""
    job = build(tmp_path / "run_001", {"001_RFdiffusion": {"streams": {"structures": 5}}})
    (job / "002_PyMOL").mkdir()
    assert [s["folder"] for s in lineage.collect(job, fs=LocalFS())] == ["001_RFdiffusion"]


def test_dropped_ids_are_counted_per_step(tmp_path):
    job = build(tmp_path / "run_001", {
        "004_Filter": {"missing": 26},
        "008_Filter": {"missing": 2}})
    assert {s["folder"]: s["missing"] for s in lineage.collect(job, fs=LocalFS())} == {
        "004_Filter": 26, "008_Filter": 2}


def test_a_step_reports_both_what_it_made_and_what_it_lost(tmp_path):
    job = build(tmp_path / "run_001",
                {"004_Filter": {"streams": {"structures": 74}, "missing": 26}})
    step = lineage.collect(job, fs=LocalFS())[0]
    assert step["produced"] == 74 and step["missing"] == 26


def test_produced_is_the_largest_stream_when_a_step_writes_several(tmp_path):
    job = build(tmp_path / "run_001",
                {"001_A": {"streams": {"structures": 10, "sequences": 40}}})
    assert lineage.collect(job, fs=LocalFS())[0]["produced"] == 40


class TestSummary:
    def test_it_names_the_largest_single_loss(self, tmp_path):
        job = build(tmp_path / "run_001", {
            "004_Filter": {"missing": 26},
            "022_Filter": {"missing": 40},
            "031_Filter": {"missing": 34}})
        text = lineage.summarize(lineage.collect(job, fs=LocalFS()), job="run_001")
        assert "100 ID(s) dropped across 3 step(s)" in text
        assert "largest single loss is 40 at 022_Filter" in text
        assert 'bp_table(step="022_Filter", table="missing.csv")' in text

    def test_a_clean_run_says_so_rather_than_staying_silent(self, tmp_path):
        job = build(tmp_path / "run_001", {"001_A": {"streams": {"structures": 5}}})
        assert "No step reported dropped IDs." in lineage.summarize(
            lineage.collect(job, fs=LocalFS()))

    def test_a_large_campaign_keeps_only_the_steps_that_lost_ids(self, tmp_path):
        """253 steps buries the attrition, which is the only part anyone reads."""
        spec = {f"{i:03d}_Boltz2": {"streams": {"structures": 9}} for i in range(1, 60)}
        spec["030_Filter"] = {"missing": 7}
        text = lineage.summarize(lineage.collect(build(tmp_path / "r", spec), fs=LocalFS()))
        assert "030_Filter" in text
        assert "001_Boltz2" not in text
        assert "step(s) with no dropped IDs not listed" in text

    def test_a_run_with_no_lineage_tables_explains_itself(self):
        assert "nothing to trace" in lineage.summarize([], job="run_001")


class TestTables:
    def test_tables_lists_only_csvs(self, tmp_path):
        job = build(tmp_path / "run_001", {"001_A": {"streams": {"structures": 2}}})
        io.open(job / "001_A" / "structures" / "out.pdb", "w", encoding="utf-8").write("x")
        assert lineage.tables(job, "001_A", fs=LocalFS()) == ["structures/structures_map.csv"]

    def test_tables_finds_what_a_real_step_wrote(self, tmp_path):
        """The live shape: a stream's map beside its files, standalone tables under `tables/`.

        `Template/example_002/002_Scripting` held exactly this and `bp_table` answered "wrote no
        CSV tables" — a confident denial about a step whose outputs the completion check had
        just counted.
        """
        job = build(tmp_path / "run_001",
                    {"002_Scripting": {"streams": {"structures": 1}, "missing": 2}})
        io.open(job / "002_Scripting" / "tables" / "metrics.csv", "w", encoding="utf-8").write(
            "id,n_residues\n168L,164\n")
        assert lineage.tables(job, "002_Scripting", fs=LocalFS()) == [
            "structures/structures_map.csv", "tables/metrics.csv", "tables/missing.csv"]

    def test_the_frameworks_own_folders_are_not_results(self, tmp_path):
        """`_configuration/scripting_dropped.csv` is plumbing; listing it invites reading it."""
        job = build(tmp_path / "run_001", {"002_Scripting": {"streams": {"structures": 1}}})
        (job / "002_Scripting" / "_configuration").mkdir(parents=True, exist_ok=True)
        io.open(job / "002_Scripting" / "_configuration" / "scripting_dropped.csv", "w",
                encoding="utf-8").write("id\n")
        assert lineage.tables(job, "002_Scripting", fs=LocalFS()) == [
            "structures/structures_map.csv"]

    @pytest.mark.parametrize("asked", ["missing", "missing.csv", "tables/missing.csv", "MISSING"])
    def test_a_table_is_found_however_it_is_named(self, tmp_path, asked):
        """Nobody should have to know which folder the framework chose for a table."""
        job = build(tmp_path / "run_001", {"001_A": {"streams": {"structures": 1}, "missing": 3}})
        found = lineage.read_table(job, "001_A", asked, fs=LocalFS())
        assert "error" not in found, found
        assert found["table"] == "tables/missing.csv" and found["rows"] == 3

    def test_read_table_truncates_and_reports_the_true_row_count(self, tmp_path):
        job = build(tmp_path / "run_001", {"001_A": {"streams": {"structures": 500}}})
        found = lineage.read_table(job, "001_A", "structures_map.csv", limit=5, fs=LocalFS())
        assert found["rows"] == 500 and found["shown"] == 5 and len(found["lines"]) == 5
        assert found["header"] == "id,file,value"

    def test_a_wrong_table_name_lists_the_real_ones(self, tmp_path):
        job = build(tmp_path / "run_001", {"001_A": {"streams": {"structures": 2}}})
        found = lineage.read_table(job, "001_A", "nope.csv", fs=LocalFS())
        assert "error" in found and found["available"] == ["structures/structures_map.csv"]


class TestNestedTables:
    """Map tables sit in the stream's own subfolder, not beside it.

    Found live: `pdb_dssp_003` wrote `002_DSSP/dssp/dssp_map.csv`, `002_DSSP/ss/ss_map.csv` and
    `002_DSSP/tables/missing.csv`, and a one-level search reported the run as having no lineage
    at all — a confident "nothing to trace" about a run that traced fine.
    """

    def nested(self, tmp_path):
        job = tmp_path / "run_001"
        for sub, name, rows in (("dssp", "dssp_map.csv", 1),
                                ("ss", "ss_map.csv", 1),
                                ("tables", "missing.csv", 3)):
            (job / "002_DSSP" / sub).mkdir(parents=True, exist_ok=True)
            header = "id,file,value\n" if name.endswith("_map.csv") else "id,structure\n"
            body = header + "".join(f"d{i},f{i}\n" for i in range(rows))
            io.open(job / "002_DSSP" / sub / name, "w", encoding="utf-8").write(body)
        return job

    def test_streams_in_subfolders_are_found(self, tmp_path):
        step = lineage.collect(self.nested(tmp_path), fs=LocalFS())[0]
        assert step["streams"] == {"dssp": 1, "ss": 1}

    def test_a_nested_missing_table_is_attributed_to_the_step(self, tmp_path):
        step = lineage.collect(self.nested(tmp_path), fs=LocalFS())[0]
        assert step["folder"] == "002_DSSP" and step["missing"] == 3

    def test_the_run_is_not_reported_as_traceless(self, tmp_path):
        text = lineage.summarize(lineage.collect(self.nested(tmp_path), fs=LocalFS()))
        assert "nothing to trace" not in text
        assert "dssp 1, ss 1" in text


class TestAgainstAFrameworkBuiltRun:
    """The anchor case: a tree the framework laid out, not one this file imagined.

    Every other case here builds its own fixture, and for a while that fixture and the code
    under test shared the same wrong picture of where a step keeps its tables — so they agreed,
    CI was green, and `bp_table` answered "wrote no CSV tables" about every real run. This case
    cannot drift: `produced_run` takes each path from what the tool declared.
    """

    def test_tables_and_lineage_agree_with_what_was_written(self, produced_run):
        job, steps = produced_run([
            {"streams": {"structures": 3}, "tables": {"metrics": [{"id": "d0", "n_res": 164}]}},
        ])
        step = next(iter(steps))
        assert lineage.tables(job, step, fs=LocalFS()) == [
            "structures/structures_map.csv", "tables/metrics.csv"]

        found = lineage.read_table(job, step, "metrics", fs=LocalFS())
        assert found["rows"] == 1 and found["header"] == "id,n_res"

        counted = lineage.collect(job, fs=LocalFS())
        assert counted[0]["streams"] == {"structures": 3}


def _missing(job, folder, rows):
    (job / folder / "tables").mkdir(parents=True, exist_ok=True)
    body = "id,removed_by,kind,cause\n" + "".join(f"{i},{by},filter,x\n" for i, by in rows)
    io.open(job / folder / "tables" / "missing.csv", "w", encoding="utf-8").write(body)


def test_a_drop_propagated_downstream_is_counted_once_where_it_happened(tmp_path):
    """Every step merges its upstream missing.csv, so a raw row count charged one filter's
    drops again to each step after it and named the last step as the largest loss."""
    job = build(tmp_path / "run_001", {
        "002_Panda": {"streams": {"sequences": 6}},
        "003_Boltz2": {"streams": {"structures": 6}},
        "004_Panda": {"streams": {"structures": 5}}})
    _missing(job, "002_Panda", [("d1", "002_Panda"), ("d2", "002_Panda")])
    _missing(job, "003_Boltz2", [("d1", "002_Panda"), ("d2", "002_Panda")])
    _missing(job, "004_Panda", [("d1", "002_Panda"), ("d2", "002_Panda"), ("d7", "004_Panda")])
    steps = {s["folder"]: s for s in lineage.collect(job, fs=LocalFS())}
    assert [steps[f]["missing"] for f in ("002_Panda", "003_Boltz2", "004_Panda")] == [2, 0, 1]
    text = lineage.summarize(list(steps.values()))
    assert "3 ID(s) dropped across 2 step(s)" in text and "largest single loss is 2 at 002_Panda" in text


def test_tables_reads_only_the_step_it_lists(tmp_path):
    """It used to count rows in every CSV of the run to list one step, once per step on the provenance page."""
    job = tmp_path / "run_001"
    (job / "003_ExtractMetrics").mkdir(parents=True)
    (job / "003_ExtractMetrics" / "plddt.csv").write_text("id,v\na,1\n")
    (job / "003_ExtractMetrics" / "tables").mkdir()
    (job / "003_ExtractMetrics" / "tables" / "missing.csv").write_text("id\n")
    (job / "004_Other" / "tables").mkdir(parents=True)
    (job / "004_Other" / "tables" / "big.csv").write_text("id\n")

    seen = []

    class Counting(LocalFS):
        def tally(self, root, names, include_root=False):
            seen.append(str(root))
            return super().tally(root, names, include_root=include_root)

    assert lineage.tables(str(job), "003_ExtractMetrics", fs=Counting()) == ["plddt.csv", "tables/missing.csv"]
    assert seen == [str(job / "003_ExtractMetrics")]
