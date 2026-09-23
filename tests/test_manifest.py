"""A run has to say what produced it without being asked.

`Pipeline(debug=True)` already captured the environment, and captured it well. But it is opt-in, and the moment you want it is the moment after a result surprised you — by which time the run is over and the flag was off. These tests hold the manifest to three claims: it is written whether or not anyone asked, it records the defaults that actually determined the output rather than only the kwargs someone typed, and two identical runs hash the same so a real difference is visible as a difference.

The stability claim is the fragile one. A parameter value can be a DataStream or an upstream tool, and a default `repr()` of either embeds a memory address; a manifest that hashed those would give every run a new hash and quietly answer "everything changed" forever.
"""

import json
import os

import pytest

from biopipelines import manifest


def build_pipeline(job, sampling_temp=0.1):
    """A real two-step pipeline: a Mock source feeding a real ProteinMPNN."""
    from biopipelines.mock import Mock
    from biopipelines.pipeline import Pipeline
    from biopipelines.protein_mpnn import ProteinMPNN

    pipeline = Pipeline(project="TestSuite", job=job, description="manifest fixture",
                        on_the_fly=False, local_output=True, config="local")
    with pipeline:
        source = Mock(ids=["d0", "d1"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
        ProteinMPNN(structures=source.streams.structures, num_sequences=4,
                    sampling_temp=sampling_temp)
        pipeline.save()
    return pipeline


@pytest.fixture
def pipelines(local_config, isolated_cwd):
    return build_pipeline


class TestItIsWrittenWithoutBeingAsked:

    def test_save_writes_a_manifest_with_debug_off(self, pipelines):
        pipeline = pipelines("unasked")
        assert pipeline.debug is False, "the point of this test is that debug is off"
        written = manifest.path_for(pipeline.folders["output"])
        assert os.path.exists(written), (
            "save() produced no manifest.json; provenance you have to enable is provenance you "
            "do not have when the run surprises you")

    def test_the_manifest_is_valid_json_naming_every_step(self, pipelines):
        pipeline = pipelines("shape")
        with open(manifest.path_for(pipeline.folders["output"]), encoding="utf-8") as handle:
            record = json.load(handle)
        assert record["schema"] == manifest.SCHEMA
        assert [tool["tool"] for tool in record["tools"]] == ["Mock", "ProteinMPNN"]
        assert record["hash"].startswith("sha256:")


class TestResolvedParameters:
    """`passed` is what someone typed. `resolved` is what decided the output."""

    def test_defaults_nobody_passed_are_recorded(self, pipelines):
        pipeline = pipelines("defaults")
        record = manifest.read(pipeline.folders["output"])
        mpnn = [tool for tool in record["tools"] if tool["tool"] == "ProteinMPNN"][0]

        passed, resolved = mpnn["parameters"]["passed"], mpnn["parameters"]["resolved"]
        assert "num_sequences" in passed and passed["num_sequences"] == 4, (
            "`passed` must be the call the user wrote — a subclass binds its named parameters "
            "before BaseConfig sees them, so reading `params` reports nothing was passed")
        assert "structures" in passed, "a stream argument is part of the call too"
        assert "model_name" not in passed, "the fixture does not pass model_name"
        assert resolved["model_name"] == "v_48_020", (
            "the default that chose the weights is absent, so the manifest cannot answer which "
            "model produced these sequences")
        assert resolved["num_sequences"] == 4, "a value that was passed must survive resolution"

    def test_a_passed_value_beats_the_default_even_when_unbound(self, pipelines):
        """ProteinMPNN consumes `structures` without binding it, so the attribute read misses.

        Falling back to the signature default there records the run as having used `None` for the
        input it was actually given -- a resolved value contradicting the passed one beside it.
        """
        record = manifest.read(pipelines("unbound").folders["output"])
        mpnn = [tool for tool in record["tools"] if tool["tool"] == "ProteinMPNN"][0]
        parameters = mpnn["parameters"]
        assert "structures" in parameters["passed"]
        assert parameters["resolved"]["structures"] == parameters["passed"]["structures"]
        assert parameters["resolved"]["structures"] is not None

    def test_a_stream_argument_is_named_not_reprd(self, pipelines):
        pipeline = pipelines("stream")
        text = open(manifest.path_for(pipeline.folders["output"]), encoding="utf-8").read()
        assert "0x" not in text, (
            "a memory address reached the manifest; every run would then hash differently and "
            "the comparison would report a change that did not happen")


class TestRenderingIsStable:
    """Guards the branch a real pipeline happens not to reach, which is where the trap lives.

    Streams and upstream tools are caught by named branches, so a fixture built from them proves
    nothing about the fallback. An unrecognized object is the one that would arrive carrying an
    address.
    """

    def test_an_unrecognized_object_collapses_to_its_type(self):
        class Opaque:
            pass

        assert manifest._jsonable(Opaque()) == "<Opaque>"

    def test_two_instances_of_the_same_class_render_identically(self):
        class Opaque:
            pass

        assert manifest._jsonable(Opaque()) == manifest._jsonable(Opaque()), (
            "two equivalent runs would hash differently on a value neither of them chose")

    def test_a_nested_value_is_rendered_too(self):
        class Opaque:
            pass

        rendered = manifest._jsonable({"b": [Opaque()], "a": 1})
        assert rendered == {"a": 1, "b": ["<Opaque>"]}


class TestTheHashMeansSomething:

    def test_two_identical_runs_hash_the_same(self, pipelines):
        first = manifest.read(pipelines("twin_a").folders["output"])
        second = manifest.read(pipelines("twin_b").folders["output"])
        assert first["created"] != second["created"] or first["job"] != second["job"], (
            "the two runs must differ in something volatile or this proves nothing")
        assert first["hash"] == second["hash"], (
            "the same pipeline hashed differently; a hash that changes on its own cannot flag a "
            "change that matters")

    def test_a_changed_default_changes_the_hash(self, pipelines):
        cool = manifest.read(pipelines("cool", sampling_temp=0.1).folders["output"])
        warm = manifest.read(pipelines("warm", sampling_temp=0.3).folders["output"])
        assert cool["hash"] != warm["hash"], "a changed sampling temperature left the hash alone"

    def test_compare_names_the_parameter_that_moved(self, pipelines):
        cool = manifest.read(pipelines("cmp_cool", sampling_temp=0.1).folders["output"])
        warm = manifest.read(pipelines("cmp_warm", sampling_temp=0.3).folders["output"])
        differences = manifest.compare(cool, warm)
        assert any("sampling_temp" in line for line in differences), (
            f"compare() did not name sampling_temp; it said {differences}")
        assert all("created" not in line for line in differences), (
            "a timestamp is not a reason two runs differ")

    def test_identical_manifests_compare_clean(self, pipelines):
        record = manifest.read(pipelines("clean").folders["output"])
        assert manifest.compare(record, dict(record)) == []


def build_table_pipeline(job, threshold=0.5):
    """Mock -> ProteinMPNN -> Panda over a table: the wiring 1.5.0 could not record."""
    from biopipelines.mock import Mock
    from biopipelines.pipeline import Pipeline
    from biopipelines.protein_mpnn import ProteinMPNN
    from biopipelines.panda import Panda

    pipeline = Pipeline(project="TestSuite", job=job, description="table fixture",
                        on_the_fly=False, local_output=True, config="local")
    with pipeline:
        source = Mock(ids=["d0", "d1"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
        mpnn = ProteinMPNN(structures=source.streams.structures, num_sequences=4)
        Panda(tables=mpnn.tables.sequences, operations=[Panda.filter(f"score > {threshold}")])
        pipeline.save()
    return pipeline


class TestTablesAndOperations:
    """A table argument answers every attribute probe, which crashed the whole manifest in 1.5.0."""

    def test_a_panda_over_a_table_still_writes_a_manifest(self, local_config, isolated_cwd):
        pipeline = build_table_pipeline("panda_table")
        record = manifest.read(pipeline.folders["output"])
        assert record is not None, f"no manifest: {manifest.write.last_error}"
        panda = [tool for tool in record["tools"] if tool["tool"] == "Panda"][0]
        assert panda["parameters"]["passed"]["tables"].startswith("<table sequences ")

    def test_the_filter_expression_is_part_of_the_hash(self, local_config, isolated_cwd):
        loose = manifest.read(build_table_pipeline("loose", threshold=0.5).folders["output"])
        strict = manifest.read(build_table_pipeline("strict", threshold=0.9).folders["output"])
        assert loose["hash"] != strict["hash"]
        assert any("0.9" in line for line in manifest.compare(loose, strict))

    def test_a_table_reference_renders_without_an_absolute_path(self):
        from biopipelines.data_containers import TableInfo
        table = TableInfo("sequences", "/abs/root/job_001/002_ProteinMPNN/sequences/sequences.csv")
        assert manifest._jsonable(table) == "<table sequences 002_ProteinMPNN/sequences/sequences.csv>"
        assert manifest._jsonable(table.score) == "<column score of 002_ProteinMPNN/sequences/sequences.csv>"

    def test_an_unrenderable_value_degrades_to_its_type(self):
        class Hostile:
            TOOL_NAME = "Hostile"
            output_folder = object()
        assert manifest._jsonable(Hostile()) == "<Hostile>"


class TestGrouped:

    def test_regrouping_changes_the_rendering(self):
        from biopipelines.combinatorics import Grouped
        from biopipelines.datastream import DataStream
        rows = DataStream(name="sequences", ids=["a_A", "a_B"], files=[], format="csv")
        by_design = DataStream(name="designs", ids=["a"], files=[], format="csv")
        by_backbone = DataStream(name="structures", ids=["bb"], files=[], format="pdb")
        assert (manifest._jsonable(Grouped(rows, groups=by_design))
                != manifest._jsonable(Grouped(rows, groups=by_backbone)))


class TestCompareNeverHidesAHashDifference:

    def test_forwarded_arguments_are_named(self):
        base = {"hash": "sha256:a", "tools": [{"order": 1, "tool": "X", "parameters": {
            "passed": {}, "resolved": {}, "forwarded": {"--seed": 1}}}]}
        other = {"hash": "sha256:b", "tools": [{"order": 1, "tool": "X", "parameters": {
            "passed": {}, "resolved": {}, "forwarded": {"--seed": 2}}}]}
        lines = manifest.compare(base, other)
        assert any("forwarded --seed" in line for line in lines)

    def test_an_unlisted_difference_is_still_reported(self):
        lines = manifest.compare({"hash": "sha256:a", "tools": []}, {"hash": "sha256:b", "tools": []})
        assert lines and lines[0].startswith("hash:")

    def test_an_environment_without_a_digest_is_not_called_compared(self):
        left = {"tools": [{"environments": ["A", "B"]}], "environments_resolved": {"A": "x"}}
        right = {"tools": [{"environments": ["A", "B"]}], "environments_resolved": {"A": "x"}}
        assert manifest.environments_were_compared(left, right) is False


class TestInputFileContent:

    def test_the_same_path_with_other_content_hashes_differently(self, tmp_path):
        target = tmp_path / "in.pdb"
        target.write_text("ATOM 1\n")
        first = manifest._input_files({"structures": str(target)})
        target.write_text("ATOM 2\n")
        second = manifest._input_files({"structures": str(target)})
        assert first["structures"].startswith("in.pdb:") and first != second


class TestReadSaysWhatWentWrong:

    def test_an_absent_manifest_is_none(self, tmp_path):
        assert manifest.read(str(tmp_path)) is None

    def test_a_corrupt_manifest_raises_instead_of_reading_as_absent(self, tmp_path):
        """None is reported as "this run predates the manifest", which sends the reader looking for a version problem."""
        (tmp_path / manifest.FILENAME).write_text("{not json", encoding="utf-8")
        with pytest.raises(ValueError, match="not valid JSON"):
            manifest.read(str(tmp_path))


class TestReviewOfTheFixes:

    def test_a_dataclass_field_named_type_keeps_its_value(self):
        from biopipelines.panda import Panda
        rendered = manifest._jsonable(Panda.filter("score > 0.5"))
        assert rendered["type"] == "filter" and rendered["__class__"] == "Operation"

    def test_a_short_sha_and_a_missing_dirty_flag_are_not_differences(self):
        old = {"hash": "sha256:a", "commit": "abc1234", "tools": []}
        new = {"hash": "sha256:b", "commit": "abc1234def5678", "dirty": False, "tools": []}
        assert all(not line.startswith(("commit", "dirty")) for line in manifest.compare(old, new))

    def test_hashing_stops_at_its_file_budget(self, tmp_path):
        paths = []
        for i in range(205):
            p = tmp_path / f"f{i}.pdb"
            p.write_text("x")
            paths.append(str(p))
        digests = manifest._input_files({"structures": paths})["structures"]
        assert sum(":size=" in d for d in digests) == 5
