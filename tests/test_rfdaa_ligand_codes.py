# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""RFdiffusionAllAtom: each ligand's own residue code reaches the run made for it.

The wrapper resolved the code of the FIRST ligand of the compounds stream once and passed it to `inference.ligand=` for every input structure, so a second ligand with another code (prepared ligands `lig_a` = UNL and `lig_b` = ASP) was designed against `UNL` and RFdiffusion-AllAtom looked for a residue the file does not hold. The code is now resolved inside the loop, for the structure's own id when the compounds stream has that id, and for every structure when the stream holds one ligand (broadcast over several proteins). A structure that matches none of several ligands was still given the first ligand's code in silence; it is now a `failure` row of `missing` naming the ligands it could have been, and the other structures go on.

Nothing here runs RFdiffusion-AllAtom; the resolver is driven on synthetic streams and the script is inspected.
"""

import os
import subprocess
import sys

import pytest

from biopipelines.datastream import DataStream, create_map_table

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESOLVER = os.path.join(REPO_ROOT, "pipe_scripts", "resolve_rfdaa_ligand_code.py")

LIGANDS = {"lig_a": "UNL", "lig_b": "ASP"}


def _compounds(tmp_path, codes=LIGANDS):
    table = tmp_path / "compounds_map.csv"
    create_map_table(str(table), ids=list(codes), additional_columns={"code": list(codes.values())})
    stream = DataStream(name="compounds", ids=list(codes), files=[], map_table=str(table), format="csv")
    stream.save_json(str(tmp_path / "ligand.json"))
    return stream


def _structures(tmp_path, ids):
    folder = tmp_path / "prepared"
    folder.mkdir(exist_ok=True)
    for sid in ids:
        (folder / f"{sid}.pdb").write_text("END\n")
    template = str(folder / "<id>.pdb")
    table = tmp_path / "structures_map.csv"
    create_map_table(str(table), ids=list(ids), files=[template])
    return DataStream(name="structures", ids=list(ids), files=[template], map_table=str(table), format="pdb")


def _resolve(tmp_path, struct_id):
    return subprocess.run([sys.executable, RESOLVER, str(tmp_path / "ligand.json"), struct_id],
                          capture_output=True, text=True)


def _code_for(tmp_path, struct_id):
    done = _resolve(tmp_path, struct_id)
    assert done.returncode == 0, done.stderr
    return done.stdout.strip()


class TestResolver:

    def test_a_structure_named_after_a_ligand_gets_that_ligands_code(self, tmp_path):
        _compounds(tmp_path)
        assert _code_for(tmp_path, "lig_a") == "UNL"
        assert _code_for(tmp_path, "lig_b") == "ASP"

    def test_a_structure_that_matches_none_of_several_ligands_is_refused_not_given_the_first(self, tmp_path):
        _compounds(tmp_path)
        done = _resolve(tmp_path, "protein_1")
        assert done.returncode != 0
        assert done.stdout.strip() == ""
        assert "protein_1" in done.stderr and "lig_a" in done.stderr and "lig_b" in done.stderr

    def test_one_ligand_is_broadcast_over_every_structure(self, tmp_path):
        _compounds(tmp_path, codes={"9DP": "9DP"})
        assert [_code_for(tmp_path, sid) for sid in ("a_1_3", "b_2")] == ["9DP", "9DP"]


class TestScript:

    def _body(self, new_pipeline, tmp_path):
        import pathlib
        from biopipelines.rfdiffusion_allatom import RFdiffusionAllAtom

        pipeline = new_pipeline("rfdaa_ligand_codes")
        with pipeline:
            RFdiffusionAllAtom(pdb=_structures(tmp_path, list(LIGANDS)), ligand=_compounds(tmp_path),
                               contigs="60-60", num_designs=2)
            script_path = pipeline.save()
        return next(pathlib.Path(script_path).parent.glob("*_RFdiffusionAllAtom.sh")).read_text(encoding="utf-8", errors="replace")

    def test_the_code_is_resolved_for_each_structure_inside_the_loop(self, local_config, isolated_cwd, new_pipeline,
                                                                       tmp_path):
        body = self._body(new_pipeline, tmp_path)
        loop = body[body.index("for STRUCT_ID in"):]
        assert 'resolve_rfdaa_ligand_code.py' in loop
        resolved = loop.index("resolve_rfdaa_ligand_code.py")
        assert loop.index("LIGAND_CODE=") < loop.index("AA_OPTIONS=(")
        assert resolved < loop.index("AA_OPTIONS=(")
        assert '"$STRUCT_ID"' in loop[resolved:resolved + 400]

    def test_the_first_ligand_is_no_longer_resolved_once_before_the_loop(self, local_config, isolated_cwd,
                                                                         new_pipeline, tmp_path):
        body = self._body(new_pipeline, tmp_path)
        assert "LIGAND_ID=" not in body[:body.index("for STRUCT_ID in")]

    def test_the_loop_skips_a_structure_whose_code_cannot_be_resolved(self, local_config, isolated_cwd, new_pipeline,
                                                                       tmp_path):
        body = self._body(new_pipeline, tmp_path)
        loop = body[body.index("for STRUCT_ID in"):]
        head = loop[:loop.index("OPTS=$(")]
        assert 'LIGAND_ERR="' in head and '/$STRUCT_ID.stderr"' in head
        assert "if ! LIGAND_CODE=$(python" in head and '2> "$LIGAND_ERR"' in head
        assert 'cat "$LIGAND_ERR" >&2' in head and "continue" in head

    def test_the_unconditional_run_still_reads_the_first_ligand(self, local_config, isolated_cwd, new_pipeline,
                                                                  tmp_path):
        import pathlib
        from biopipelines.rfdiffusion_allatom import RFdiffusionAllAtom

        pipeline = new_pipeline("rfdaa_unconditional")
        with pipeline:
            RFdiffusionAllAtom(ligand=_compounds(tmp_path), contigs="60-60", num_designs=2)
            script_path = pipeline.save()
        body = next(pathlib.Path(script_path).parent.glob("*_RFdiffusionAllAtom.sh")).read_text(
            encoding="utf-8", errors="replace")
        assert "LIGAND_ID=" in body and "LIGAND_CODE=" in body
        assert "resolve_rfdaa_ligand_code.py" not in body


class TestAStructureMatchingNoLigandFailsAlone:

    def test_the_resolver_message_is_the_cause_of_its_failure_row(self, tmp_path):
        import importlib.util

        _compounds(tmp_path)
        spec = importlib.util.spec_from_file_location(
            "update_map", os.path.join(REPO_ROOT, "pipe_scripts", "pipe_update_structures_map.py"))
        update_map = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(update_map)
        stderr_dir = tmp_path / "007_RFdiffusionAllAtom" / "_execution" / "stderr"
        stderr_dir.mkdir(parents=True)
        (stderr_dir / "protein_1.stderr").write_text(_resolve(tmp_path, "protein_1").stderr)
        rows = update_map.record_failures(
            str(tmp_path / "007_RFdiffusionAllAtom" / "_execution" / "local_missing.csv"),
            {"protein_1_1": "protein_1"}, {"protein_1"}, set(), stderr_dir=str(stderr_dir))
        assert [row["kind"] for row in rows] == ["failure"]
        assert "protein_1" in rows[0]["cause"] and "lig_a" in rows[0]["cause"]
