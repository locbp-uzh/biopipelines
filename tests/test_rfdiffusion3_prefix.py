"""`prefix=` has to reach the output ids, or sibling steps collide silently.

With a PDB input the design key was the input pdb id alone and `prefix` was consulted only in the de-novo branch, so a packed run with one step per design target -- same parent PDB, different contigs -- had every step emit `<pdb>_1 ... <pdb>_N`. Nothing failed: each step writes its own folder, so the collision only appears once the steps are pooled, where N steps land on N ids and the rest are lost. One real run put 3,700 designs onto 100 ids.

The key is now `<prefix>_<pdb id>`, which the declared ids, the runtime inputs JSON and the provenance regex all have to agree on.
"""

import csv
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]


@pytest.fixture
def rfd3(local_config, isolated_cwd):
    def _make(**kwargs):
        from biopipelines import PDB, RFdiffusion3
        from biopipelines.pipeline import Pipeline, Resources
        with Pipeline(project="TestSuite", job="rfd3_prefix"):
            Resources(gpu=None)
            parent = PDB("4AKE")
            out = RFdiffusion3(pdb=parent, contig="A1-50,10-20,A51-100",
                               num_designs=3, **kwargs)
        return out.producer
    return _make


class TestTheDeclaredIdsCarryThePrefix:

    def test_a_prefix_prepends_to_every_declared_id(self, rfd3):
        ids = rfd3(prefix="op7").get_output_files()["structures"].ids
        assert all("op7_" in str(i) for i in ids), ids

    def test_without_a_prefix_the_ids_are_unchanged(self, rfd3):
        ids = rfd3().get_output_files()["structures"].ids
        assert all("4AKE" in str(i) for i in ids), ids

    def test_two_steps_on_one_parent_no_longer_collide(self, rfd3):
        first = rfd3(prefix="op7").get_output_files()["structures"].ids
        second = rfd3(prefix="op8").get_output_files()["structures"].ids
        assert set(map(str, first)).isdisjoint(map(str, second))


class TestTheRuntimeKeysMatchTheDeclaredIds:
    """The builder writes the foundry inputs JSON; its keys become the output filenames."""

    def _build(self, tmp_path, pdb_id, key_prefix):
        template = tmp_path / "template.json"
        template.write_text(json.dumps({"contig": "A1-50,10-20,A51-100"}))
        pdb = tmp_path / f"{pdb_id}.pdb"
        pdb.write_text("ATOM\n")
        ds = tmp_path / "ds.json"
        ds.write_text(json.dumps({"name": "structures", "ids": [pdb_id],
                                  "files": [str(pdb)], "format": "pdb"}))
        out = tmp_path / "inputs.json"
        cmd = [sys.executable, str(REPO / "pipe_scripts" / "pipe_rfdiffusion3_build_inputs.py"),
               "--template", str(template), "--output", str(out),
               "--structures-json", str(ds)]
        if key_prefix:
            cmd += ["--key-prefix", key_prefix]
        subprocess.run(cmd, check=True, capture_output=True)
        return json.loads(out.read_text())

    def test_the_key_is_prefixed(self, tmp_path):
        assert list(self._build(tmp_path, "4AKE", "op7")) == ["op7_4AKE"]

    def test_the_bare_id_is_kept_without_a_prefix(self, tmp_path):
        assert list(self._build(tmp_path, "4AKE", None)) == ["4AKE"]


class TestProvenanceStillRecoversTheParent:
    """structures_map.csv joins each design back to the PDB it came from; the prefix sits in front of that id, so the regex has to take it off as well as the design suffix."""

    def _regex(self, rfd3, prefix):
        script = rfd3(prefix=prefix)._generate_script_update_structures_map()
        return re.search(r'--suffix-regex "([^"]+)"', script).group(1)

    def test_the_parent_id_survives_a_prefix(self, rfd3):
        pattern = self._regex(rfd3, "op7")
        assert re.sub(pattern, "", "op7_4AKE_3") == "4AKE"

    def test_a_prefix_with_a_regex_metacharacter_is_taken_literally(self, rfd3):
        pattern = self._regex(rfd3, "38-40.55-59")
        assert re.sub(pattern, "", "38-40.55-59_4AKE_3") == "4AKE"
        assert re.sub(pattern, "", "38-40X55-59_4AKE_3") == "38-40X55-59_4AKE"

    def test_the_regex_reaches_the_helper_intact_through_bash(self, rfd3, tmp_path):
        """The pattern is interpolated into a double-quoted shell argument."""
        pattern = self._regex(rfd3, "op7")
        folder = tmp_path / "structures"
        folder.mkdir()
        (folder / "op7_4AKE_1.pdb").write_text("ATOM\n")
        out = tmp_path / "structures_map.csv"
        subprocess.run(
            f'{sys.executable} "{REPO / "pipe_scripts" / "pipe_update_structures_map.py"}"'
            f' --structures-map "{out}" --output-folder "{folder}"'
            f' --provenance-from-suffix "structures.id" --suffix-regex "{pattern}"',
            shell=True, check=True, capture_output=True)
        rows = list(csv.DictReader(out.open()))
        assert rows[0]["structures.id"] == "4AKE", rows
