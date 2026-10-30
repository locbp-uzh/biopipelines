# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""PLACER: a run that gives no model is a `failure` row of `missing`, and a run that cannot work is not started.

On the cluster a structure without the ligand's residue code made `predict_ligand` segfault inside PLACER. The driver (one process for every run) died with it, the wrapper stopped before the post-processing, and the step wrote no stream, no scores and no `missing.csv`, although the runs before it had finished. And a run that raised a Python error was only printed: the post-processing found its units from the model files, so an id with no file was in no table.

Nothing here runs PLACER; the driver is run against a stand-in `PLACER` package and the post-processing on synthetic runs.
"""

import json
import os
import subprocess
import sys
import textwrap

import pytest

pd = pytest.importorskip("pandas")

from biopipelines.datastream import DataStream, create_map_table
from biopipelines.placer import PLACER
from tests.test_id_construction_pair_sites import SCRIPTS, _run

WITH_LIGAND = (
    "ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N\n"
    "HETATM    2  C1  LDP A 500       5.000   5.000   5.000  1.00  0.00           C\n"
    "END\n")
WITHOUT_LIGAND = (
    "ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N\n"
    "END\n")

STAND_IN = textwrap.dedent('''
    import json, os

    class PLACERinput:
        def __init__(self):
            self.ligands = []
        def pdb(self, path): self.path = path
        def cif(self, path): self.path = path
        def name(self, name): self.id = name
        def predict_ligand(self, codes): self.ligands = codes
        def predict_multi(self, flag): pass
        def exclude_sm(self, flag): pass

    class PLACER:
        started = []
        def run(self, placer_input, nsamples):
            PLACER.started.append(placer_input.id)
            with open(os.path.join(os.environ["STAND_IN_LOG"]), "a") as log:
                log.write(placer_input.id + "\\n")
            if "raises" in placer_input.id:
                raise RuntimeError("stand-in PLACER refuses this run")
            return [object()] * nsamples

    class protocol:
        @staticmethod
        def dump_output(outputs, prefix, rerank=None):
            with open(prefix + "_model.pdb", "w") as handle:
                handle.write("".join(f"MODEL {n}\\nATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N\\nENDMDL\\n"
                                     for n in range(1, len(outputs) + 1)))
            with open(prefix + ".csv", "w") as handle:
                handle.write("model_idx,plddt\\n" + "".join(f"{n},0.5\\n" for n in range(1, len(outputs) + 1)))
''')


def _streams(tmp_path, structures):
    """Write the structures stream and a one-ligand compounds stream the driver reads."""
    ids = []
    for name, text in structures.items():
        (tmp_path / f"{name}.pdb").write_text(text)
        ids.append(name)
    table = tmp_path / "structures_map.csv"
    create_map_table(str(table), ids=ids)
    DataStream(name="structures", ids=ids, map_table=str(table), format="pdb",
               files=[str(tmp_path / "<id>.pdb")]).save_json(str(tmp_path / "structures.json"))
    compounds = tmp_path / "compounds.csv"
    create_map_table(str(compounds), ids=["l_1"], additional_columns={"code": ["LDP"]})
    DataStream(name="compounds", ids=["l_1"], map_table=str(compounds), format="csv",
               files=[]).save_json(str(tmp_path / "ligand.json"))


def _drive(tmp_path, structures):
    _streams(tmp_path, structures)
    repo = tmp_path / "placer_repo"
    repo.mkdir()
    (repo / "PLACER.py").write_text(STAND_IN)
    log = tmp_path / "started.log"
    log.write_text("")
    env = dict(os.environ, STAND_IN_LOG=str(log), PYTHONPATH=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    runs = tmp_path / "runs"
    done = subprocess.run(
        [sys.executable, os.path.join(SCRIPTS, "pipe_placer_driver.py"), "--placer-repo", str(repo),
         "--structures-json", str(tmp_path / "structures.json"), "--runs-folder", str(runs), "--nsamples", "2",
         "--mode", "ligand", "--ligand-json", str(tmp_path / "ligand.json")],
        capture_output=True, text=True, env=env)
    return done, runs, log.read_text().split()


class TestTheDriverDoesNotStartARunThatCannotWork:

    def test_a_structure_without_the_ligand_code_is_not_given_to_placer(self, tmp_path):
        done, runs, started = _drive(tmp_path, {"holo": WITH_LIGAND, "apo+1": WITHOUT_LIGAND})

        assert done.returncode == 0, done.stderr
        assert started == ["holo+l_1"], "PLACER segfaults on a structure that has no residue of the ligand's code"
        assert (runs / "holo+l_1_model.pdb").is_file()
        failures = json.loads((runs / "failures.json").read_text())
        assert list(failures) == ["apo+1+l_1"]
        assert "LDP" in failures["apo+1+l_1"]

    def test_a_run_that_raises_is_recorded_with_its_error(self, tmp_path):
        done, runs, started = _drive(tmp_path, {"holo": WITH_LIGAND, "raises": WITH_LIGAND})

        assert done.returncode == 0, done.stderr
        failures = json.loads((runs / "failures.json").read_text())
        assert list(failures) == ["raises+l_1"]
        assert "stand-in PLACER refuses this run" in failures["raises+l_1"]

    def test_every_run_failing_exits_non_zero_and_still_names_each(self, tmp_path):
        done, runs, started = _drive(tmp_path, {"apo": WITHOUT_LIGAND})

        assert done.returncode != 0
        assert list(json.loads((runs / "failures.json").read_text())) == ["apo+l_1"]


class TestThePostProcessingNamesEveryRunThatGaveNoModel:

    def _collect(self, tmp_path, finished, manifest, failures=None):
        runs = tmp_path / "runs"
        runs.mkdir()
        for unit in finished:
            (runs / f"{unit}_model.pdb").write_text(
                "MODEL 1\nATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N\nENDMDL\n")
            (runs / f"{unit}.csv").write_text("model_idx,plddt\n1,0.5\n")
        (runs / "units.json").write_text(json.dumps(manifest))
        if failures is not None:
            (runs / "failures.json").write_text(json.dumps(failures))
        (tmp_path / "ligand.json").write_text(json.dumps(
            {"name": "compounds", "ids": ["l_1"], "files": [], "map_table": str(tmp_path / "compounds.csv"), "format": "csv"}))
        (tmp_path / "compounds.csv").write_text("id,code\nl_1,LDP\n")
        step = tmp_path / "003_PLACER"
        done = _run("pipe_placer_postprocess.py", "--mode", "ligand", "--runs-folder", str(runs),
                    "--structures-folder", str(step / "structures"),
                    "--structures-map", str(step / "structures" / "structures_map.csv"),
                    "--scores-csv", str(step / "tables" / "scores.csv"),
                    "--missing-csv", str(step / "tables" / "missing.csv"),
                    "--ligand-json", str(tmp_path / "ligand.json"),
                    "--compounds-map", str(step / "compounds" / "compounds_map.csv"))
        return done, step

    MANIFEST = {
        "holo+l_1": {"structures.id": "holo", "compounds.id": "l_1"},
        "apo+l_1": {"structures.id": "apo", "compounds.id": "l_1"},
        "dead+l_1": {"structures.id": "dead", "compounds.id": "l_1"},
    }

    def test_a_run_with_a_recorded_error_is_a_failure_row_with_that_error(self, tmp_path):
        done, step = self._collect(tmp_path, ["holo+l_1"], self.MANIFEST,
                                   failures={"apo+l_1": "ValueError: no residue LDP in the structure",
                                             "dead+l_1": "RuntimeError: out of memory"})

        assert done.returncode == 0, done.stderr
        missing = pd.read_csv(step / "tables" / "missing.csv")
        assert sorted(missing["id"]) == ["apo+l_1", "dead+l_1"]
        assert set(missing["kind"]) == {"failure"}
        assert set(missing["removed_by"]) == {"003_PLACER"}
        causes = dict(zip(missing["id"], missing["cause"]))
        assert causes["apo+l_1"] == "ValueError: no residue LDP in the structure"
        table = pd.read_csv(step / "structures" / "structures_map.csv")
        assert list(table["id"]) == ["holo+l_1_1"]

    def test_a_run_the_driver_never_reached_is_a_failure_row_too(self, tmp_path):
        """The driver died (a segmentation fault) before it recorded anything."""
        done, step = self._collect(tmp_path, ["holo+l_1"], self.MANIFEST)

        assert done.returncode == 0, done.stderr
        missing = pd.read_csv(step / "tables" / "missing.csv")
        assert sorted(missing["id"]) == ["apo+l_1", "dead+l_1"]
        assert all("no model" in cause for cause in missing["cause"])


class TestTheWrapperCollectsWhatFinished:

    def test_the_post_processing_runs_when_the_driver_exits_non_zero(self, tmp_path):
        """A driver that was killed after some runs must not stop the step before the post-processing."""
        structures = DataStream(name="structures", ids=["s"], map_table=str(tmp_path / "s.csv"), format="pdb",
                                files=[str(tmp_path / "<id>.pdb")])
        create_map_table(str(tmp_path / "s.csv"), ids=["s"])
        ligands = DataStream(name="compounds", ids=["l"], map_table=str(tmp_path / "l.csv"), format="csv", files=[])
        create_map_table(str(tmp_path / "l.csv"), ids=["l"], additional_columns={"code": ["LDP"]})
        tool = PLACER(structures=structures, ligand=ligands)
        tool.output_folder = str(tmp_path / "003_PLACER")
        tool.folders = {"PLACER": "/placer", "pipe_scripts": "/repo/pipe_scripts"}
        script = tool.generate_script("unused.sh")

        driver_at = script.index("pipe_placer_driver.py")
        post_at = script.index("pipe_placer_postprocess.py")
        between = script[driver_at:post_at]
        assert "exit 1" not in between, "the step must not stop between the driver and the post-processing"
        assert "$?" in between, "the driver's status is read"
