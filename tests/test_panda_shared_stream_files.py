"""Panda pool mode over a stream whose ids are coarser than the table it filters by.

BioEmu2 declares `structures` per conformer (`<run>_<k>.pdb`) and `trajectories` per run
(`<run>.xtc`, one file for all `num_samples` conformers of that run). Filtering the per-conformer
`physicality` table therefore selects conformer ids that no `trajectories` row carries; the id
matcher resolves each of them to its run (the parent tier), and the extraction used to copy the
run's whole trajectory once per surviving conformer under the conformer's name.

That left the declared `<run>.xtc` unwritten, so the step was marked FAILED, and it wrote one full
copy of a multi-megabyte trajectory per conformer. These tests pin what a many-ids-to-one-file
stream must do instead: one copy per distinct source file, under the id the stream declares it by,
with every surviving conformer still resolving to it through the map table.
"""

import csv
import json
import os
import subprocess
import sys

import pytest

PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"
XTC = b"XTC-" + b"x" * 4096


def _repo_root():
    import biopipelines.panda as _p
    return os.path.dirname(os.path.dirname(os.path.abspath(_p.__file__)))


def _run_panda(tool):
    script = os.path.join(_repo_root(), "pipe_scripts", "pipe_panda.py")
    config = os.path.join(tool.output_folder, "_configuration", "panda_config.json")
    subprocess.run([sys.executable, script, "--config", config], check=True)


def _check_completion(tool):
    sys.path.insert(0, os.path.join(_repo_root(), "pipe_scripts"))
    import pipe_check_completion

    with open(os.path.join(tool.output_folder, ".expected_outputs.json")) as handle:
        manifest = json.load(handle)
    return pipe_check_completion.check_expected_outputs(
        manifest["output_structure"], manifest["tool_name"], tool.output_folder)


def _rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path, columns, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def _stage_bioemu2_outputs(be, runs, num_samples, rejected=()):
    """Write what pipe_bioemu2.py leaves behind, for `runs` runs of `num_samples` conformers each."""
    structures_dir = os.path.join(be.output_folder, "structures")
    trajectories_dir = os.path.join(be.output_folder, "trajectories")
    os.makedirs(structures_dir, exist_ok=True)
    os.makedirs(trajectories_dir, exist_ok=True)

    structure_rows, trajectory_rows, physicality_rows = [], [], []
    for run in runs:
        trajectory = os.path.join(trajectories_dir, f"{run}.xtc")
        topology = os.path.join(trajectories_dir, f"{run}_topology.cif")
        with open(trajectory, "wb") as handle:
            handle.write(XTC + run.encode())
        with open(topology, "w") as handle:
            handle.write(f"data_{run}\n")
        trajectory_rows.append({"id": run, "file": trajectory, "topology": topology})
        for k in range(1, num_samples + 1):
            frame_id = f"{run}_{k}"
            path = os.path.join(structures_dir, f"{frame_id}.pdb")
            with open(path, "w") as handle:
                handle.write(f"ATOM  {frame_id}\n")
            structure_rows.append({"id": frame_id, "file": path, "run.id": run,
                                   "sequences.id": run, "ligands.id": ""})
            physicality_rows.append({"id": frame_id, "run.id": run, "frame": k - 1,
                                     "status": "rejected" if frame_id in rejected else "passed",
                                     "reason": ""})

    _write_csv(be.streams.structures.map_table,
               ["id", "file", "run.id", "sequences.id", "ligands.id"], structure_rows)
    _write_csv(be.streams.trajectories.map_table, ["id", "file", "topology"], trajectory_rows)
    _write_csv(be.tables.physicality.info.path,
               ["id", "run.id", "frame", "status", "reason"], physicality_rows)
    _write_csv(be.tables.summary.info.path,
               ["id", "sequences.id", "ligands.id", "n_samples", "n_residues"],
               [{"id": run, "sequences.id": run, "ligands.id": "", "n_samples": num_samples,
                 "n_residues": 33} for run in runs])
    _write_csv(be.tables.missing.info.path, ["id", "removed_by", "kind", "cause"], [])


def _bioemu2_then_panda(new_pipeline, operations, table="physicality", runs=("p1", "p2"),
                        num_samples=3, rejected=(), **panda_kwargs):
    from biopipelines.bioemu2 import BioEmu2
    from biopipelines.panda import Panda
    from biopipelines.sequence import Sequence

    pipeline = new_pipeline("panda_shared_files")
    with pipeline:
        seqs = Sequence(seq=[PROTEIN] * len(runs), type="protein", ids=list(runs))
        be = BioEmu2(sequences=seqs, num_samples=num_samples)
        pan = Panda(tables=getattr(be.tables, table), operations=operations, pool=be, **panda_kwargs)
        pipeline.save()

    _stage_bioemu2_outputs(be, runs, num_samples, rejected)
    _run_panda(pan)
    return be, pan


class TestAPerRunTrajectoryFilteredByConformer:

    def test_the_step_is_not_failed_by_a_trajectory_it_never_renamed(
            self, local_config, isolated_cwd, new_pipeline):
        from biopipelines.panda import Panda
        _be, pan = _bioemu2_then_panda(
            new_pipeline, [Panda.filter("status == 'passed'")], rejected={"p1_2"})

        ok, missing = _check_completion(pan)
        assert ok, f"a filtered per-conformer table failed the step over the trajectories: {missing}"

    def test_each_trajectory_is_copied_once_not_once_per_conformer(
            self, local_config, isolated_cwd, new_pipeline):
        from biopipelines.panda import Panda
        _be, pan = _bioemu2_then_panda(new_pipeline, [Panda.filter("status == 'passed'")])

        written = sorted(os.listdir(os.path.join(pan.output_folder, "trajectories")))
        copies = [name for name in written if name.endswith(".xtc")]
        assert copies == ["p1.xtc", "p2.xtc"], (
            f"expected one trajectory per run, found {copies}; a copy per conformer is what the "
            "3 conformers of each run used to write")

    def test_every_surviving_conformer_still_reaches_its_trajectory(
            self, local_config, isolated_cwd, new_pipeline):
        from biopipelines.panda import Panda
        _be, pan = _bioemu2_then_panda(
            new_pipeline, [Panda.filter("status == 'passed'")], rejected={"p1_2"})

        rows = {row["id"]: row for row in _rows(pan.streams.trajectories.map_table)}
        assert set(rows) == {"p1", "p2"}
        for run, row in rows.items():
            assert os.path.exists(row["file"]), row
            assert os.path.dirname(row["file"]) == os.path.join(pan.output_folder, "trajectories")
            with open(row["file"], "rb") as handle:
                assert handle.read() == XTC + run.encode(), "the copy is not the run's own trajectory"


class TestARunWhoseEveryConformerWasFiltered:

    def test_its_trajectory_is_accounted_for_instead_of_failing_the_step(
            self, local_config, isolated_cwd, new_pipeline):
        from biopipelines.panda import Panda
        _be, pan = _bioemu2_then_panda(
            new_pipeline, [Panda.filter("status == 'passed'")], rejected={"p2_1", "p2_2", "p2_3"})

        ok, missing = _check_completion(pan)
        assert ok, f"a run with nothing left failed the step over its own trajectory: {missing}"
        assert sorted(os.listdir(os.path.join(pan.output_folder, "trajectories"))) == [
            "p1.xtc", "p1_topology.cif", "trajectories_map.csv"]

        lost = [row for row in _rows(os.path.join(pan.output_folder, "tables", "missing.csv"))
                if row["id"] == "p2"]
        assert [row["kind"] for row in lost] == ["filter"], (
            "the run is gone because the filter removed all of it, which is a filter row")

    def test_a_run_that_kept_one_conformer_still_keeps_its_trajectory(
            self, local_config, isolated_cwd, new_pipeline):
        from biopipelines.panda import Panda
        _be, pan = _bioemu2_then_panda(
            new_pipeline, [Panda.filter("status == 'passed'")], rejected={"p1_1", "p1_3"})

        ok, missing = _check_completion(pan)
        assert ok, missing
        assert os.path.exists(os.path.join(pan.output_folder, "trajectories", "p1.xtc"))


def _can_link(directory):
    source = os.path.join(directory, "probe_source")
    with open(source, "w") as handle:
        handle.write("x")
    for make in (os.symlink, os.link):
        target = os.path.join(directory, f"probe_{make.__name__}")
        try:
            make(source, target)
            return True
        except (OSError, NotImplementedError):
            continue
    return False


def _independent_files(paths):
    """How many distinct files these names resolve to, following links."""
    return len({(os.stat(path).st_dev, os.stat(path).st_ino) for path in paths})


class TestARenamedSelection:

    def test_every_new_id_has_a_trajectory_and_each_run_is_written_once(
            self, local_config, isolated_cwd, new_pipeline, tmp_path):
        from biopipelines.panda import Panda
        _be, pan = _bioemu2_then_panda(
            new_pipeline, [Panda.sort("frame"), Panda.head(4)], rename="sel")

        ok, missing = _check_completion(pan)
        assert ok, missing

        rows = _rows(pan.streams.trajectories.map_table)
        assert [row["id"] for row in rows] == ["sel_1", "sel_2", "sel_3", "sel_4"]
        for row in rows:
            assert os.path.exists(row["file"]), row
            with open(row["file"], "rb") as handle:
                assert handle.read() in (XTC + b"p1", XTC + b"p2")
        if _can_link(str(tmp_path)):
            assert _independent_files(row["file"] for row in rows) == 2, (
                "four conformers of two runs should be backed by two trajectories on disk")


class TestSeveralIdsSharingOneFile:
    """The shape without BioEmu2: a file map in which some ids carry their own file and some share one."""

    @staticmethod
    def _stream(tmp_path):
        source = tmp_path / "upstream"
        source.mkdir()
        shared, own = source / "shared.xtc", source / "own.xtc"
        shared.write_bytes(b"shared-content")
        own.write_bytes(b"own-content")
        return {"traj": {"a": str(shared), "b": str(shared), "c": str(shared), "d": str(own)}}

    @staticmethod
    def _extract(file_map, out, rename_map=None, ids=("a", "b", "c", "d")):
        from pipe_scripts.pipe_panda import extract_pool_data_for_filtered_ids
        os.makedirs(out, exist_ok=True)
        return extract_pool_data_for_filtered_ids(
            list(ids), "pool", str(out), rename_map=rename_map, file_map=file_map)

    def test_every_id_has_its_file_and_the_shared_one_is_stored_once(self, tmp_path):
        out = tmp_path / "out"
        extracted = self._extract(self._stream(tmp_path), out)

        names = sorted(os.path.basename(path) for path in extracted["traj"])
        assert names == ["a.xtc", "b.xtc", "c.xtc", "d.xtc"]
        for name, content in (("a", b"shared-content"), ("b", b"shared-content"),
                              ("c", b"shared-content"), ("d", b"own-content")):
            assert (out / f"{name}.xtc").read_bytes() == content
        if _can_link(str(tmp_path)):
            assert _independent_files(extracted["traj"]) == 2

    def test_a_rename_gives_each_new_id_the_shared_file_as_well(self, tmp_path):
        out = tmp_path / "out"
        rename = {"a": "best_1", "b": "best_2", "c": "best_3", "d": "best_4"}
        extracted = self._extract(self._stream(tmp_path), out, rename_map=rename)

        assert sorted(os.path.basename(path) for path in extracted["traj"]) == [
            "best_1.xtc", "best_2.xtc", "best_3.xtc", "best_4.xtc"]
        assert (out / "best_3.xtc").read_bytes() == b"shared-content"
        assert (out / "best_4.xtc").read_bytes() == b"own-content"

    def test_a_second_run_over_the_same_folder_does_not_write_through_a_link(self, tmp_path):
        out = tmp_path / "out"
        file_map = self._stream(tmp_path)
        self._extract(file_map, out)
        self._extract(file_map, out)

        assert (out / "a.xtc").read_bytes() == b"shared-content"
        assert (out / "b.xtc").read_bytes() == b"shared-content"
        assert (out / "d.xtc").read_bytes() == b"own-content"

    def test_a_link_left_by_an_earlier_run_is_replaced_not_followed(self, tmp_path):
        """`b.xtc` was a link to `a.xtc`; when it is now another file's own name, it must not overwrite `a.xtc`."""
        if not _can_link(str(tmp_path)):
            pytest.skip("the filesystem supports neither symlinks nor hard links")
        out = tmp_path / "out"
        file_map = self._stream(tmp_path)
        self._extract(file_map, out)

        other = tmp_path / "upstream" / "other.xtc"
        other.write_bytes(b"other-content")
        file_map["traj"]["b"] = str(other)
        self._extract(file_map, out)

        assert (out / "a.xtc").read_bytes() == b"shared-content"
        assert (out / "b.xtc").read_bytes() == b"other-content"


class TestIdsOnACoarserStream:

    @staticmethod
    def _per_run(tmp_path, runs):
        source = tmp_path / "upstream"
        source.mkdir(exist_ok=True)
        file_map = {"traj": {}}
        for run in runs:
            path = source / f"{run}.xtc"
            path.write_bytes(run.encode())
            file_map["traj"][run] = str(path)
        return file_map

    def test_conformer_ids_resolve_to_one_file_per_run_under_the_runs_own_id(self, tmp_path):
        from pipe_scripts.pipe_panda import extract_pool_data_for_filtered_ids
        out = tmp_path / "out"
        out.mkdir()
        extracted = extract_pool_data_for_filtered_ids(
            ["p1_1", "p1_3", "p2_2"], "pool", str(out),
            file_map=self._per_run(tmp_path, ["p1", "p2"]))

        assert sorted(os.path.basename(path) for path in extracted["traj"]) == ["p1.xtc", "p2.xtc"]
        assert sorted(os.listdir(out)) == ["p1.xtc", "p2.xtc"]

    def test_the_multi_pool_path_agrees(self, tmp_path):
        import pandas as pd
        from pipe_scripts.pipe_panda import PANDA_SOURCE, extract_from_multiple_pools
        out = tmp_path / "out"
        out.mkdir()
        frame = pd.DataFrame({"id": ["p1_1", "p1_2", "q1_1"], PANDA_SOURCE: [0, 0, 1]})
        extracted = extract_from_multiple_pools(
            frame, ["pool0", "pool1"], str(out),
            file_maps=[self._per_run(tmp_path, ["p1"]), self._per_run(tmp_path, ["q1"])])

        assert sorted(os.path.basename(path) for path in extracted["traj"]) == ["p1.xtc", "q1.xtc"]
        assert sorted(os.listdir(out)) == ["p1.xtc", "q1.xtc"]

    def test_lost_stream_ids_names_only_runs_the_filter_emptied(self):
        from pipe_scripts.pipe_panda import lost_stream_ids
        stream = {"p1": "p1.xtc", "p2": "p2.xtc", "p3": "p3.xtc", "z": "z.xtc"}
        original = ["p1_1", "p1_2", "p2_1", "p2_2"]
        lost = lost_stream_ids(original, ["p1_1", "p2_2"], stream)
        assert lost == [], "both runs kept a conformer"
        assert lost_stream_ids(original, ["p1_1", "p1_2"], stream) == ["p2"], (
            "p3 and z were never reached by a table id, so the filter did not remove them")

    def test_a_stream_id_that_is_itself_a_table_id_is_left_to_the_ordinary_route(self):
        from pipe_scripts.pipe_panda import lost_stream_ids
        assert lost_stream_ids(["a", "b"], ["a"], {"a": "a.pdb", "b": "b.pdb"}) == []


def _topology_files(pan):
    """{id: topology path} from the pooled trajectories map, failing if the column is gone."""
    rows = _rows(pan.streams.trajectories.map_table)
    assert rows and "topology" in rows[0], (
        f"the pooled trajectories map lost its topology column: {list(rows[0]) if rows else rows}")
    return {row["id"]: row["topology"] for row in rows}


class TestAnExtraFileColumnOfTheStreamMap:
    """A trajectory without its topology cannot be read, so the second file a map row names travels with it."""

    def test_the_topology_of_each_run_is_copied_and_the_column_points_at_the_copy(
            self, local_config, isolated_cwd, new_pipeline):
        from biopipelines.panda import Panda
        be, pan = _bioemu2_then_panda(new_pipeline, [Panda.filter("status == 'passed'")])

        topology = _topology_files(pan)
        assert sorted(topology) == ["p1", "p2"]
        folder = os.path.join(pan.output_folder, "trajectories")
        for run, path in topology.items():
            assert os.path.dirname(path) == folder, f"{run}: the column still names the upstream file {path}"
            with open(path) as handle:
                assert handle.read() == f"data_{run}\n"

    def test_a_rename_gives_each_new_id_a_readable_topology_stored_once_per_run(
            self, local_config, isolated_cwd, new_pipeline, tmp_path):
        from biopipelines.panda import Panda
        _be, pan = _bioemu2_then_panda(
            new_pipeline, [Panda.sort("frame"), Panda.head(4)], rename="sel")

        topology = _topology_files(pan)
        assert sorted(topology) == ["sel_1", "sel_2", "sel_3", "sel_4"]
        for path in topology.values():
            assert os.path.dirname(path) == os.path.join(pan.output_folder, "trajectories")
            with open(path) as handle:
                assert handle.read() in ("data_p1\n", "data_p2\n")
        if _can_link(str(tmp_path)):
            assert _independent_files(topology.values()) == 2

    def test_a_run_that_kept_one_conformer_keeps_its_topology(
            self, local_config, isolated_cwd, new_pipeline):
        from biopipelines.panda import Panda
        _be, pan = _bioemu2_then_panda(
            new_pipeline, [Panda.filter("status == 'passed'")], rejected={"p1_1", "p1_3"})

        with open(_topology_files(pan)["p1"]) as handle:
            assert handle.read() == "data_p1\n"

    def test_every_pool_of_a_multi_pool_step_brings_its_own_topologies(
            self, local_config, isolated_cwd, new_pipeline):
        from biopipelines.bioemu2 import BioEmu2
        from biopipelines.panda import Panda
        from biopipelines.sequence import Sequence

        pipeline = new_pipeline("panda_shared_files_multi")
        with pipeline:
            first = BioEmu2(sequences=Sequence(seq=PROTEIN, type="protein", ids="p1"), num_samples=3)
            second = BioEmu2(sequences=Sequence(seq=PROTEIN, type="protein", ids="q1"), num_samples=3)
            pan = Panda(tables=[first.tables.physicality, second.tables.physicality],
                        operations=[Panda.concat(), Panda.filter("status == 'passed'")],
                        pool=[first, second])
            pipeline.save()
        _stage_bioemu2_outputs(first, ("p1",), 3)
        _stage_bioemu2_outputs(second, ("q1",), 3)
        _run_panda(pan)

        topology = _topology_files(pan)
        assert sorted(topology) == ["p1", "q1"]
        for run, path in topology.items():
            with open(path) as handle:
                assert handle.read() == f"data_{run}\n"


class TestCarryingTheExtraColumnsOfAMapRow:

    @staticmethod
    def _carry(tmp_path, columns, row_id="r", staged=None, name="traj"):
        from pipe_scripts.pipe_panda import stage_row_columns
        out = tmp_path / "out"
        out.mkdir(exist_ok=True)
        return stage_row_columns(name, row_id, columns, str(out), {} if staged is None else staged), out

    def test_a_value_naming_a_file_is_copied_and_rewritten_others_are_kept(self, tmp_path):
        top = tmp_path / "top.cif"
        top.write_text("data_x")
        carried, out = self._carry(tmp_path, {"topology": str(top), "n_atoms": 120,
                                              "label": "closed", "empty": "", "run.id": "p1"})

        assert carried["topology"] == str(out / "r_topology.cif")
        assert (out / "r_topology.cif").read_text() == "data_x"
        assert (carried["n_atoms"], carried["label"], carried["empty"], carried["run.id"]) == (
            120, "closed", "", "p1")

    def test_a_value_that_only_looks_like_a_path_is_left_alone(self, tmp_path):
        carried, out = self._carry(tmp_path, {"gone": str(tmp_path / "deleted.cif"),
                                              "relative": "top.cif", "folder": str(tmp_path)})

        assert carried == {"gone": str(tmp_path / "deleted.cif"), "relative": "top.cif",
                           "folder": str(tmp_path)}
        assert list(out.iterdir()) == []

    def test_rows_sharing_a_source_store_it_once(self, tmp_path):
        top = tmp_path / "top.cif"
        top.write_text("data_x")
        staged = {}
        first, out = self._carry(tmp_path, {"topology": str(top)}, row_id="a", staged=staged)
        second, _ = self._carry(tmp_path, {"topology": str(top)}, row_id="b", staged=staged)

        assert sorted(path.name for path in out.iterdir()) == ["a_topology.cif", "b_topology.cif"]
        assert (out / "b_topology.cif").read_text() == "data_x"
        if _can_link(str(tmp_path)):
            assert _independent_files([first["topology"], second["topology"]]) == 1
