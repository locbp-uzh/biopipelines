"""The parent of a design is recorded by exact id, never recovered from the id text.

The RFdiffusion family and HBDesigner write `structures_map.csv` at runtime with a `structures.id` column naming the input each design was made from. `pipe_update_structures_map.py` used to derive that column by stripping a regex the wrapper handed it from every file id, so any file whose name did not follow the wrapper's shape got a parent that was merely the name with its last `_<digits>` cut off: a stray `design_2_7.pdb` claimed the parent `design_2`, silently. The parent is now looked up by exact id in the map `combinatorics.multiplied_id_parents` builds from the inputs that really exist at runtime (the input stream's own ids, which a lazy bracket resolves against the upstream map_table), and a file that is no input's output carries no parent and is named on stderr.

What must hold, and what each test pins:
  - the shared function maps each output id to its own parent when parents and outputs both carry digits and underscores (`5HG6_1_3` is an output of `5HG6_1` and a parent in its own right);
  - the ids it builds are the ids the wrapper declares at configuration time;
  - the wrapper-generated command, run for real over a folder of the declared files, writes the same ids and the right parent for each, and a stray file keeps a blank parent;
  - a folder in which nothing matches fails loudly instead of writing a column of blanks;
  - no wrapper or pipe script spells a suffix regex any more.
"""

import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "pipe_scripts" / "pipe_update_structures_map.py"

# Digits and underscores in the parent: a strip of the last `_<digits>` gives the wrong parent for a name that is not an output.
PARENTS = ["design_2_7", "5HG6_1_3", "5HG6_1"]
DESIGNS = 3


class TestTheSharedFunction:

    def test_each_output_maps_to_its_own_parent(self):
        from biopipelines.combinatorics import multiplied_id_parents

        mapped = multiplied_id_parents(PARENTS, "<1..3>")
        assert mapped == {
            "design_2_7_1": "design_2_7", "design_2_7_2": "design_2_7", "design_2_7_3": "design_2_7",
            "5HG6_1_3_1": "5HG6_1_3", "5HG6_1_3_2": "5HG6_1_3", "5HG6_1_3_3": "5HG6_1_3",
            "5HG6_1_1": "5HG6_1", "5HG6_1_2": "5HG6_1", "5HG6_1_3": "5HG6_1",
        }

    def test_an_output_that_is_also_a_parent_keeps_the_parent_it_was_made_from(self):
        from biopipelines.combinatorics import multiplied_id_parents

        mapped = multiplied_id_parents(PARENTS, "<1..3>")
        assert mapped["5HG6_1_3"] == "5HG6_1"
        assert mapped["5HG6_1_3_2"] == "5HG6_1_3"

    def test_the_ids_are_the_declared_ids_expanded_in_order(self):
        from biopipelines import id_patterns
        from biopipelines.combinatorics import generate_multiplied_ids_pattern, multiplied_id_parents

        declared = id_patterns.expand_ids(generate_multiplied_ids_pattern(PARENTS, "<1..3>"))
        assert list(multiplied_id_parents(PARENTS, "<1..3>")) == declared

    def test_a_key_prefix_names_the_output_and_leaves_the_parent_bare(self):
        from biopipelines.combinatorics import multiplied_id_parents

        mapped = multiplied_id_parents(["5HG6_1"], "d<1..2>_m<1..2>", key_prefix="op7")
        assert mapped == {
            "op7_5HG6_1_d1_m1": "5HG6_1", "op7_5HG6_1_d1_m2": "5HG6_1",
            "op7_5HG6_1_d2_m1": "5HG6_1", "op7_5HG6_1_d2_m2": "5HG6_1",
        }


def _stream_json(tmp_path, ids, map_table=""):
    path = tmp_path / "parents.json"
    path.write_text(json.dumps({"name": "structures", "ids": ids, "files": ["<id>.pdb"],
                                "map_table": str(map_table), "format": "pdb"}))
    return path


def _run_update(tmp_path, parents_json, files, suffix="<1..3>", key_prefix=None):
    folder = tmp_path / "structures"
    folder.mkdir(exist_ok=True)
    for name in files:
        (folder / f"{name}.pdb").write_text("ATOM\n")
    out = tmp_path / "structures_map.csv"
    cmd = [sys.executable, str(SCRIPT), "--structures-map", str(out), "--output-folder", str(folder),
           "--parents-json", str(parents_json), "--suffix-pattern", suffix,
           "--parent-column", "structures.id"]
    if key_prefix:
        cmd += ["--key-prefix", key_prefix]
    done = subprocess.run(cmd, capture_output=True, text=True)
    rows = dict(pd.read_csv(out, dtype=str, keep_default_na=False)[["id", "structures.id"]].values) \
        if out.exists() else None
    return done, rows


class TestThePipeScript:

    def test_a_stray_file_is_not_given_a_parent_by_string_surgery(self, tmp_path):
        """`design_2_7.pdb` is no input's output; stripping its `_7` would call its parent `design_2`."""
        parents = _stream_json(tmp_path, PARENTS)
        done, rows = _run_update(tmp_path, parents, ["design_2_7_1", "5HG6_1_3_2", "design_2_7"])
        assert done.returncode == 0, done.stderr
        assert rows == {"design_2_7_1": "design_2_7", "5HG6_1_3_2": "5HG6_1_3", "design_2_7": ""}
        assert "design_2_7" in done.stderr and "WARNING" in done.stderr

    def test_a_lazy_parent_is_resolved_from_the_upstream_map_at_runtime(self, tmp_path):
        """The configuration cannot name these parents; the runtime map of the upstream step does."""
        upstream = tmp_path / "upstream_map.csv"
        upstream.write_text("id,file\nbase_1,a.pdb\nbase_12,b.pdb\n")
        parents = _stream_json(tmp_path, ["base[_<#>]"], map_table=upstream)
        done, rows = _run_update(tmp_path, parents, ["base_1_3", "base_12_1"])
        assert done.returncode == 0, done.stderr
        assert rows == {"base_1_3": "base_1", "base_12_1": "base_12"}

    def test_a_parent_an_upstream_filter_dropped_has_no_outputs(self, tmp_path):
        """A range pattern states the design space; the upstream map says which ids exist."""
        upstream = tmp_path / "upstream_map.csv"
        upstream.write_text("id,file\np_1,a.pdb\np_3,b.pdb\n")
        parents = _stream_json(tmp_path, ["p_<1..3>"], map_table=upstream)
        done, rows = _run_update(tmp_path, parents, ["p_1_1", "p_3_2", "p_2_1"])
        assert done.returncode == 0, done.stderr
        assert rows == {"p_1_1": "p_1", "p_3_2": "p_3", "p_2_1": ""}

    def test_a_prefix_is_taken_off_the_output_and_not_off_the_parent(self, tmp_path):
        parents = _stream_json(tmp_path, ["op7_4AKE", "4AKE"])
        done, rows = _run_update(tmp_path, parents, ["op7_op7_4AKE_1", "op7_4AKE_2"],
                                 suffix="<1..3>", key_prefix="op7")
        assert done.returncode == 0, done.stderr
        assert rows == {"op7_op7_4AKE_1": "op7_4AKE", "op7_4AKE_2": "4AKE"}

    def test_a_folder_in_which_nothing_matches_fails_instead_of_writing_blanks(self, tmp_path):
        """Config and runtime disagreeing on the id shape would otherwise blank the whole column."""
        parents = _stream_json(tmp_path, PARENTS)
        done, rows = _run_update(tmp_path, parents, ["unrelated_1", "unrelated_2"])
        assert done.returncode != 0
        assert rows is None, "the map must not be written with a column of blanks"
        assert "unrelated_1" in done.stderr

    def test_a_step_without_parents_writes_the_map_as_before(self, tmp_path):
        folder = tmp_path / "structures"
        folder.mkdir()
        (folder / "a.pdb").write_text("ATOM\n")
        out = tmp_path / "structures_map.csv"
        subprocess.run([sys.executable, str(SCRIPT), "--structures-map", str(out),
                        "--output-folder", str(folder)], check=True, capture_output=True)
        assert list(pd.read_csv(out)["id"]) == ["a"]


def _plain(parent, designs=DESIGNS):
    return {f"{parent}_{n}": parent for n in range(1, designs + 1)}


def _expected_plain():
    return {i: p for parent in PARENTS for i, p in _plain(parent).items()}


def _expected_rfd3_prefixed():
    return {f"op7_{p}_d{d}_m{m}": p for p in PARENTS for d in (1, 2, 3) for m in (1, 2)}


def _cases():
    def rfdiffusion(ds):
        from biopipelines.rfdiffusion import RFdiffusion
        return RFdiffusion(pdb=ds, contigs="A1-50/10-20", num_designs=DESIGNS)

    def rfdiffusion2(ds):
        from biopipelines.ligand import Ligand
        from biopipelines.rfdiffusion2 import RFdiffusion2
        return RFdiffusion2(pdb=ds, ligand=Ligand(codes="LIG"), contigs="A1-50/10-20", num_designs=DESIGNS)

    def allatom(ds):
        from biopipelines.ligand import Ligand
        from biopipelines.rfdiffusion_allatom import RFdiffusionAllAtom
        return RFdiffusionAllAtom(pdb=ds, ligand=Ligand(codes="LIG"), contigs="A1-50/10-20", num_designs=DESIGNS)

    def rfd3(ds):
        from biopipelines.rfdiffusion3 import RFdiffusion3
        return RFdiffusion3(pdb=ds, contig="A1-50,10-20,A51-100", num_designs=DESIGNS)

    def rfd3_prefixed(ds):
        from biopipelines.rfdiffusion3 import RFdiffusion3
        return RFdiffusion3(pdb=ds, contig="A1-50,10-20,A51-100", num_designs=DESIGNS,
                            num_models=2, prefix="op7")

    def hbdesigner(ds):
        from biopipelines.hbdesigner import HBDesigner
        return HBDesigner(structures=ds, top_k=DESIGNS)

    return {
        "RFdiffusion": (rfdiffusion, _expected_plain),
        "RFdiffusion2": (rfdiffusion2, _expected_plain),
        "RFdiffusionAllAtom": (allatom, _expected_plain),
        "RFdiffusion3": (rfd3, _expected_plain),
        "RFdiffusion3 prefix and models": (rfd3_prefixed, _expected_rfd3_prefixed),
        "HBDesigner": (hbdesigner, _expected_plain),
    }


CASES = _cases()


@pytest.fixture
def producer(local_config, isolated_cwd):
    def _make(case):
        from biopipelines.datastream import DataStream
        from biopipelines.pipeline import Pipeline, Resources

        build, _expected = CASES[case]
        with Pipeline(project="TestSuite", job="map_parents", description="x",
                      on_the_fly=False, local_output=True, config="local"):
            Resources(gpu=None)
            ds = DataStream(name="structures", ids=list(PARENTS), files=["<id>.pdb"],
                            map_table="", format="pdb")
            built = build(ds)
        return built.tool if hasattr(built, "tool") else built.producer
    return _make


@pytest.mark.parametrize("case", list(CASES))
class TestTheWrappersAgreeWithTheirDeclaredIds:

    def test_the_command_carries_no_suffix_regex(self, producer, case):
        script = producer(case)._generate_script_update_structures_map()
        assert "--suffix-regex" not in script and "--provenance-from-suffix" not in script, script
        assert "--parents-json" in script and "--suffix-pattern" in script, script

    def test_runtime_provenance_equals_the_declared_ids_and_their_parents(self, producer, case, tmp_path):
        """Run the generated command over a folder holding exactly the declared files, plus one stray."""
        from biopipelines import id_patterns

        cfg = producer(case)
        cfg.generate_script("unused.sh")
        expected = CASES[case][1]()
        declared = id_patterns.expand_ids(cfg.get_output_files()["structures"].ids)
        assert sorted(declared) == sorted(expected)

        folder = Path(cfg.stream_folder("structures"))
        folder.mkdir(parents=True, exist_ok=True)
        for design_id in declared + ["design_2_7"]:
            (folder / f"{design_id}.pdb").write_text("ATOM\n")

        line = next(l for l in cfg._generate_script_update_structures_map().splitlines()
                    if l.startswith("python "))
        cmd = f'"{sys.executable}" ' + line[len("python "):]
        done = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        assert done.returncode == 0, done.stderr

        written = pd.read_csv(cfg.stream_map_path("structures"), dtype=str, keep_default_na=False)
        by_id = dict(zip(written["id"], written["structures.id"]))
        assert {i: p for i, p in by_id.items() if i != "design_2_7"} == expected
        assert by_id["design_2_7"] == ""


class TestNoSuffixRegexSurvives:
    """The option is gone from the helper and from every module that drove it."""

    @pytest.mark.parametrize("folder", ["biopipelines", "pipe_scripts"])
    def test_no_module_spells_a_suffix_regex(self, folder):
        pattern = re.compile(r"suffix[_-]regex|provenance[_-]from[_-]suffix")
        offenders = [str(p.relative_to(REPO)) for p in (REPO / folder).rglob("*.py")
                     if pattern.search(p.read_text(encoding="utf-8", errors="ignore"))]
        assert not offenders, offenders
