"""One structure LigandMPNN cannot design must not take the batch down.

A shell loop calls upstream `run.py` once per structure. An empty or unparseable PDB made that call exit non-zero, and nothing recorded it: the `missing` table had no row for the structure, and its raw `<id>.fa` was declared but never written, so the completion check marked the whole step FAILED while the good structures' sequences sat on disk. ProteinMPNN records the same case under the STRUCTURE id (kind="failure"), carries on, and its step completes.

These run the generated script end to end under bash, with only upstream's `run.py` replaced by a stand-in that fails on a structure without atoms.
"""

import json
import os
import pathlib
import re
import shutil
import subprocess
import sys

import pandas as pd
import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
POSITIONS = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_lmpnn_runtime_positions.py")
RESOLVER = os.path.join(REPO_ROOT, "pipe_scripts", "resolve_lmpnn_positions.py")

bash_required = pytest.mark.skipif(shutil.which("bash") is None, reason="runs the generated script under bash")

GOOD_PDB = "".join(
    f"ATOM  {i:5d}  CA  ALA A{i:4d}    {i * 3.8:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C\n"
    for i in range(1, 7)) + "END\n"
BINARY_JUNK = b"\x81\x9d\xff\xfe\x00junk\x80"

# Stands in for upstream: exits 1 on a structure without atoms, writes nothing for one marked SILENT, otherwise writes the native record plus the designs.
FAKE_RUN = '''
import argparse, os, sys
p = argparse.ArgumentParser()
for flag in ("--model_type", "--checkpoint_ligand_mpnn", "--batch_size", "--number_of_batches",
             "--ligand_mpnn_cutoff_for_score", "--out_folder", "--pdb_path"):
    p.add_argument(flag)
args, _ = p.parse_known_args()
text = open(args.pdb_path, errors="replace").read()
if "SILENT" in text:
    sys.exit(0)
if "ATOM" not in text:
    sys.stderr.write("Traceback (most recent call last):\\nIndexError: no residues in " + os.path.basename(args.pdb_path) + "\\n")
    sys.exit(1)
name = os.path.basename(args.pdb_path)[:-4]
os.makedirs(os.path.join(args.out_folder, "seqs"), exist_ok=True)
lines = [">" + name + ", T=0.1, seed=1, num_res=6", "AAAAAA"]
for i in range(1, int(args.batch_size) * int(args.number_of_batches) + 1):
    seq = "".join(chr(65 + (ord(c) + i) % 20) for c in name * 2)[:6]
    lines += [">" + name + ", id=%d, T=0.1, seed=1, overall_confidence=0.5, ligand_confidence=0.4, seq_rec=0.3" % i, seq]
open(os.path.join(args.out_folder, "seqs", name + ".fa"), "w").write("\\n".join(lines) + "\\n")
'''


def _posix_paths(script):
    """Git-bash reads `C:\\a\\b` with the backslashes stripped; spell every drive path with slashes."""
    return re.sub(r"(?:[A-Za-z]:|\.)\\[^\s\"']*", lambda m: m.group(0).replace("\\", "/"), script)


def _build(tmp_path, new_pipeline, files, **kwargs):
    """A LigandMPNN step over `files` ({id: bytes-or-str}); the caller's cwd is `tmp_path`."""
    from biopipelines.datastream import DataStream
    from biopipelines.ligand import Ligand
    from biopipelines.ligand_mpnn import LigandMPNN

    inputs = tmp_path / "inputs"
    inputs.mkdir(exist_ok=True)
    for sid, content in files.items():
        data = content if isinstance(content, bytes) else content.encode()
        (inputs / f"{sid}.pdb").write_bytes(data)
    map_csv = tmp_path / "structures_map.csv"
    pd.DataFrame({"id": list(files), "file": [str(inputs / f"{i}.pdb") for i in files]}).to_csv(map_csv, index=False)
    stream = DataStream(name="structures", ids=list(files), files=[str(inputs / "<id>.pdb")],
                        map_table=str(map_csv), format="pdb")

    pipeline = new_pipeline("lmpnn_fail")
    with pipeline:
        LigandMPNN(structures=stream, ligand=Ligand(codes="STI").streams.compounds,
                   redesigned="3-5", **kwargs)
    return next(t for t in pipeline.tools if t.TOOL_NAME == "LigandMPNN"), pipeline


def _run_step(tool, pipeline, tmp_path):
    """Execute the generated step script; returns (exit code, combined output)."""
    import glob

    script = glob.glob(pipeline.folders["runtime"] + "/*LigandMPNN*.sh")[0]
    text = _posix_paths(open(script, encoding="utf-8", errors="replace").read())
    os.makedirs(tool.lmpnn_folder, exist_ok=True)
    with open(os.path.join(tool.lmpnn_folder, "run.py"), "w", newline="\n") as f:
        f.write(FAKE_RUN)

    shim = tmp_path / "bin"
    shim.mkdir(exist_ok=True)
    wrapper = shim / "python"
    # Windows python writes CRLF, which would stick to every id the script captures from it.
    interpreter = pathlib.Path(sys.executable).as_posix()
    wrapper.write_text(f'#!/usr/bin/env bash\n"{interpreter}" "$@" | tr -d "\\r"\nexit ${{PIPESTATUS[0]}}\n', newline="\n")
    wrapper.chmod(0o755)
    runnable = tmp_path / "step.sh"
    runnable.write_text(text, newline="\n")

    env = dict(os.environ, PATH=str(shim).replace("\\", "/") + os.pathsep + os.environ["PATH"], PYTHONIOENCODING="utf-8")
    proc = subprocess.run(["bash", runnable.as_posix()], cwd=tmp_path, capture_output=True, text=True,
                          errors="replace", env=env, timeout=900)
    return proc.returncode, proc.stdout + proc.stderr


def _missing(tool):
    return pd.read_csv(tool.missing_csv) if os.path.exists(tool.missing_csv) else None


def _fresh_pipeline(job):
    from biopipelines.pipeline import Pipeline
    return Pipeline(project="TestSuite", job=job, description=job, on_the_fly=False,
                    local_output=True, config="local")


@pytest.fixture(scope="module")
def module_config():
    """`local_config` for module-scoped fixtures, which the function-scoped one cannot serve."""
    from _pytest.monkeypatch import MonkeyPatch
    from biopipelines.config_manager import ConfigManager

    config = os.path.join(REPO_ROOT, "tests", "fixtures", "config.local.yaml")
    mp = MonkeyPatch()

    def reset():
        ConfigManager._instance = None
        ConfigManager._config = None
        ConfigManager._variant = None

    reset()
    mp.setattr(ConfigManager, "_get_config_path", classmethod(lambda cls, variant=None: config))
    yield
    mp.undo()
    reset()


def _step_in(tmp_path, files, **kwargs):
    """Build and run a step in its own directory, restoring the working directory afterwards."""
    previous = os.getcwd()
    os.chdir(tmp_path)
    try:
        tool, pipeline = _build(tmp_path, _fresh_pipeline, files, **kwargs)
        rc, output = _run_step(tool, pipeline, tmp_path)
    finally:
        os.chdir(previous)
    return tool, rc, output


@pytest.fixture(scope="module")
def partial_run(module_config, tmp_path_factory):
    """Two good structures with an empty, a binary and a silent one between them."""
    files = {"good_a": GOOD_PDB, "empty": "", "binary": BINARY_JUNK, "silent": "ATOM SILENT\n", "good_b": GOOD_PDB}
    return _step_in(tmp_path_factory.mktemp("lmpnn_partial"), files, num_sequences=2)


@pytest.fixture(scope="module")
def all_failed_run(module_config, tmp_path_factory):
    return _step_in(tmp_path_factory.mktemp("lmpnn_all_failed"), {"empty_a": "", "empty_b": ""})


# ── the generated step, run end to end ────────────────────────────────────────

@bash_required
def test_unparsable_structures_are_listed_in_missing_under_their_own_id(partial_run):
    tool, _, output = partial_run
    missing = _missing(tool)
    assert missing is not None, output
    failures = missing[missing["kind"] == "failure"].set_index("id")
    assert sorted(failures.index) == ["binary", "empty", "silent"], (list(missing["id"]), output)
    assert "no designed sequences" in failures.loc["silent", "cause"]
    step = os.path.basename(tool.output_folder)
    assert set(failures["removed_by"]) == {step}
    assert all(isinstance(c, str) and len(c) > 10 for c in failures["cause"]), failures["cause"].tolist()


@bash_required
def test_the_loop_goes_on_after_a_failure(partial_run):
    """`good_b` comes after both junk files in the stream; its design must exist."""
    tool, _, output = partial_run
    sequences = pd.read_csv(tool.queries_csv)
    assert sorted(sequences["structures.id"].unique()) == ["good_a", "good_b"], output
    assert sorted(sequences["id"]) == ["good_a_1", "good_a_2", "good_b_1", "good_b_2"]


@bash_required
def test_a_partial_failure_exits_zero_and_the_step_completes(partial_run):
    tool, rc, output = partial_run
    parent = os.path.dirname(tool.output_folder)
    step_number = os.path.basename(tool.output_folder).split("_")[0]
    assert rc == 0, output
    assert os.path.exists(os.path.join(parent, f"{step_number}_LigandMPNN_COMPLETED")), output
    assert not os.path.exists(os.path.join(parent, f"{step_number}_LigandMPNN_FAILED")), output


@bash_required
def test_every_structure_failing_exits_non_zero_and_lists_each(all_failed_run):
    tool, rc, output = all_failed_run
    missing = _missing(tool)
    assert rc != 0, output
    assert sorted(missing[missing["kind"] == "failure"]["id"]) == ["empty_a", "empty_b"], output
    parent = os.path.dirname(tool.output_folder)
    step_number = os.path.basename(tool.output_folder).split("_")[0]
    assert os.path.exists(os.path.join(parent, f"{step_number}_LigandMPNN_FAILED"))


# ── the pieces ────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("attempted,failed,expected", [
    (0, 0, 1),   # nothing was attempted
    (3, 3, 1),   # every attempted structure failed
    (3, 2, 0),   # some structures were designed
    (3, 0, 0),
])
@bash_required
def test_exit_status_follows_the_tally(local_config, isolated_cwd, new_pipeline, tmp_path,
                                       attempted, failed, expected):
    """Checklist section 10: non-zero only when nothing was attempted or every attempt failed."""
    tool, _ = _build(tmp_path, new_pipeline, {"s": GOOD_PDB})
    block = tool._generate_script_exit_status()
    harness = f"LMPNN_ATTEMPTED={attempted}\nLMPNN_FAILED={failed}\nBP_MAIN_RC=0\n{block}\nexit $BP_MAIN_RC\n"
    path = tmp_path / "rule.sh"
    path.write_text(harness, newline="\n")
    proc = subprocess.run(["bash", path.as_posix()], capture_output=True, text=True, timeout=60)
    assert proc.returncode == expected, proc.stdout + proc.stderr


def test_a_failed_convert_step_is_not_hidden_by_the_tally(local_config, isolated_cwd, new_pipeline, tmp_path):
    tool, _ = _build(tmp_path, new_pipeline, {"s": GOOD_PDB})
    block = tool._generate_script_exit_status()
    harness = f"LMPNN_ATTEMPTED=3\nLMPNN_FAILED=0\nBP_MAIN_RC=7\n{block}\nexit $BP_MAIN_RC\n"
    path = tmp_path / "rule.sh"
    path.write_text(harness, newline="\n")
    proc = subprocess.run(["bash", path.as_posix()], capture_output=True, text=True, timeout=60)
    assert proc.returncode == 7


def test_fasta_is_one_shared_file_so_a_failed_structure_leaves_no_declared_gap(
        local_config, isolated_cwd, new_pipeline, tmp_path):
    """Per-structure `<id>.fa` files were declared, and a local failure row does not excuse an absent file."""
    tool, _ = _build(tmp_path, new_pipeline, {"a": GOOD_PDB, "b": GOOD_PDB})
    fasta = tool.get_output_files()["fasta"]
    assert fasta.files == tool.queries_fasta, fasta.files
    assert "<id>" not in str(fasta.files)


@pytest.mark.parametrize("pack_side_chains", [False, True])
def test_the_completion_check_on_a_step_with_one_failed_structure(
        local_config, isolated_cwd, new_pipeline, tmp_path, pack_side_chains):
    """Every structure but `b` has its outputs; `b` has a local failure row.

    The packed structures are per-id files, which a local failure row does not excuse, so that stream is the one a failed structure still leaves unmet.
    """
    sys.path.insert(0, os.path.join(REPO_ROOT, "pipe_scripts"))
    import pipe_check_completion as check

    tool, _ = _build(tmp_path, new_pipeline, {"a": GOOD_PDB, "b": GOOD_PDB}, pack_side_chains=pack_side_chains)
    declared = json.load(open(os.path.join(tool.output_folder, ".expected_outputs.json")))["output_structure"]

    os.makedirs(os.path.dirname(tool.queries_csv), exist_ok=True)
    os.makedirs(os.path.dirname(tool.designs_csv), exist_ok=True)
    os.makedirs(tool.seqs_folder, exist_ok=True)
    os.makedirs(os.path.dirname(tool.missing_csv), exist_ok=True)
    pd.DataFrame({"id": ["a_1"], "sequence": ["AAAA"]}).to_csv(tool.queries_csv, index=False)
    pd.DataFrame({"id": ["a_1"]}).to_csv(tool.designs_csv, index=False)
    open(tool.queries_fasta, "w").write(">a_1\nAAAA\n")
    open(os.path.join(tool.seqs_folder, "a.fa"), "w").write(">a\nAAAA\n")
    step = os.path.basename(tool.output_folder)
    pd.DataFrame([{"id": "b", "removed_by": step, "kind": "failure", "cause": "run.py exited 1"}]
                 ).to_csv(tool.missing_csv, index=False)

    if pack_side_chains:
        os.makedirs(tool.packed_folder, exist_ok=True)
        open(os.path.join(tool.packed_folder, "a_1_1.pdb"), "w").write("ATOM\n")

    ok, missing = check.check_expected_outputs(declared, "LigandMPNN", tool.output_folder)
    assert ok is (not pack_side_chains), missing
    assert list(missing) == (["structures"] if pack_side_chains else [])


# ── the position helper, which parses every structure before the loop ─────────

def _positions(tmp_path, files, designed="3-5", default_chain="auto"):
    base = str(tmp_path)
    paths = {}
    for sid, content in files.items():
        paths[sid] = os.path.join(base, f"{sid}.pdb")
        data = content if isinstance(content, bytes) else content.encode()
        with open(paths[sid], "wb") as f:
            f.write(data)
    ds = os.path.join(base, "ds.json")
    json.dump({"name": "structures", "ids": list(files), "files": list(paths.values()),
               "map_table": "", "format": "pdb"}, open(ds, "w"))
    cfg = os.path.join(base, "cfg.json")
    out = os.path.join(base, "out.json")
    json.dump({"structures_json": ds, "input_source": "selection", "input_table": "-",
               "fixed_positions": "-", "designed_positions": designed, "ligand": "-",
               "design_within": "5.0", "output_file": out, "default_chain": default_chain}, open(cfg, "w"))
    result = subprocess.run([sys.executable, POSITIONS, cfg], capture_output=True, text=True)
    return result, (json.load(open(out)) if os.path.exists(out) else None), out


@pytest.mark.parametrize("designed,default_chain", [("3-5", "auto"), ("-", "A")])
def test_an_unreadable_structure_does_not_stop_the_position_helper(tmp_path, designed, default_chain):
    """Both parse sites: the chain lookup for a chainless selection, and the fix-everything fallback."""
    result, written, _ = _positions(tmp_path, {"good": GOOD_PDB, "bad": BINARY_JUNK},
                                    designed=designed, default_chain=default_chain)
    assert result.returncode == 0, result.stderr
    assert "error" in written["bad"] and "bad.pdb" in written["bad"]["error"], written
    assert "error" not in written["good"]


def test_the_resolver_fails_for_a_structure_the_helper_could_not_read(tmp_path):
    _, _, out = _positions(tmp_path, {"good": GOOD_PDB, "bad": BINARY_JUNK})
    bad = subprocess.run([sys.executable, RESOLVER, out, "bad"], capture_output=True, text=True)
    good = subprocess.run([sys.executable, RESOLVER, out, "good"], capture_output=True, text=True)
    assert bad.returncode != 0 and "bad.pdb" in bad.stderr
    assert good.returncode == 0 and good.stdout.split("\0")[:2] == ["--redesigned_residues", "A3 A4 A5"]
