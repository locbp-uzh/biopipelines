# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Foldcomp: entry and range parsing, the cut and renumbering of a decompressed model, the runtime half end to end with the binary replaced by a stand-in that follows upstream's `decompress --id-list` behavior, and the wrapper (validation, outputs, emitted script, install bash).

The decompressed model in tests/fixtures/foldcomp (a .txt, since *.pdb is git-ignored) is real output of Foldcomp 1.0.0: 1CRN compressed and decompressed through the Python package, whose binding calls the same writeAtomCoordinatesToPDB as the command line. The real `foldcomp` binary was not available to the suite, so the command line itself is exercised only through the stand-in.
"""

import glob
import importlib.util
import json
import os
import pathlib
import re
import sys

import pandas as pd
import pytest

pytest.importorskip("yaml")
import yaml  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "pipe_scripts", "pipe_foldcomp.py")
REAL_MODEL = os.path.join(ROOT, "tests", "fixtures", "foldcomp", "crambin_decompressed.txt")


@pytest.fixture(scope="module")
def pipe():
    spec = importlib.util.spec_from_file_location("pipe_foldcomp", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def real_lines():
    with open(REAL_MODEL) as f:
        return f.read().splitlines()


def real_atoms():
    return [line for line in real_lines() if line.startswith("ATOM")]


def format_atom(serial, atom, resname, chain, number, xyz, bfactor, element):
    """One ATOM line laid out as Foldcomp's writeAtomCoordinatesToPDB does."""
    name = f"{atom:<4}" if len(atom) == 4 else f" {atom:<3}"
    x, y, z = (f"{c:8.3f}" for c in xyz)
    return f"ATOM  {serial:5d} {name} {resname:>3} {chain}{number:4d}    {x}{y}{z}  1.00{bfactor:6.2f}          {element:>2}  "


def synthetic_model(first, last, chain="A", b_of=lambda n: float(n % 100), title="synthetic"):
    """A one-chain Foldcomp-format model of residues first..last, B-factor b_of(residue number)."""
    lines, serial = [f"TITLE     {title}"], 0
    for number in range(first, last + 1):
        atoms = [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")] + ([("OXT", "O")] if number == last else [])
        for atom, element in atoms:
            serial += 1
            lines.append(format_atom(serial, atom, "ALA", chain, number, (number * 1.5, serial * 0.1, 0.0), b_of(number), element))
    lines.append(f"TER   {serial + 1:5d}      ALA {chain}{last:4d}")
    return "\n".join(lines) + "\n"


def atoms_of(text):
    return [line for line in text.splitlines() if line.startswith("ATOM")]


# ── real output format ────────────────────────────────────────────────────────

def test_the_atom_formatter_reproduces_real_foldcomp_lines():
    for line in real_atoms():
        rebuilt = format_atom(int(line[6:11]), line[12:16].strip(), line[17:20].strip(), line[21], int(line[22:26]),
                              (float(line[30:38]), float(line[38:46]), float(line[46:54])), float(line[60:66]), line[76:78].strip())
        assert rebuilt == line


# ── entries and positions ─────────────────────────────────────────────────────

def test_a_bare_accession_expands_to_the_foldcomp_name(pipe):
    assert pipe.parse_entry("A0A009F5K6", 4) == "AF-A0A009F5K6-F1-model_v4"
    assert pipe.parse_entry("P69905", 6) == "AF-P69905-F1-model_v6"
    assert pipe.parse_entry("B1XK97", 4) == "AF-B1XK97-F1-model_v4"


def test_a_full_afdb_name_is_used_as_given_which_is_how_a_fragment_is_asked_for(pipe):
    assert pipe.parse_entry("AF-Q8W3K0-F2-model_v4", 4) == "AF-Q8W3K0-F2-model_v4"


@pytest.mark.parametrize("entry", ["p69905", "P1234", "AF-P69905", "AF-P69905-F1-model_v4.pdb", "../P69905", "P69905 ", "bad entry!", ""])
def test_anything_else_is_refused_before_it_reaches_a_file_name(pipe, entry):
    with pytest.raises(ValueError, match="neither a UniProt accession"):
        pipe.parse_entry(entry, 4)


@pytest.mark.parametrize("value, expected", [(12, 12), ("12", 12), (12.0, 12), ("12.0", 12), (" 7 ", 7)])
def test_positions_arrive_as_int_text_or_float(pipe, value, expected):
    assert pipe.parse_position(value, "start") == expected


@pytest.mark.parametrize("value, message", [(float("nan"), "empty"), (None, "empty"), ("", "empty"), ("abc", "not an integer"), (1.5, "not an integer")])
def test_a_position_that_is_not_an_integer_is_refused(pipe, value, message):
    with pytest.raises(ValueError, match=message):
        pipe.parse_position(value, "end")


def test_ids_derive_from_the_entry_as_written_and_the_range(pipe):
    requests = pipe.build_requests(["A0A009F5K6", "AF-Q8W3K0-F2-model_v4"], [12, "201"], [150, 300.0], 4)
    assert [r.id for r in requests] == ["A0A009F5K6_12-150", "AF-Q8W3K0-F2-model_v4_201-300"]
    assert [r.name for r in requests] == ["AF-A0A009F5K6-F1-model_v4", "AF-Q8W3K0-F2-model_v4"]
    assert all(not r.error for r in requests)


def test_a_bad_row_keeps_its_place_with_its_error_and_the_others_are_untouched(pipe):
    requests = pipe.build_requests(["A0A009F5K6", "A0A009F5K6", "A0A009F5K6", "nope", "", "A0A009F5K6"],
                                   [1, 0, 30, 1, 1, float("nan")], [10, 5, 20, 5, 5, 9], 4)
    assert [r.id for r in requests] == ["A0A009F5K6_1-10", "A0A009F5K6_0-5", "A0A009F5K6_30-20",
                                        "nope_1-5", "row4", "A0A009F5K6_nan-9"]
    errors = [r.error for r in requests]
    assert errors[0] == ""
    assert "at least 1" in errors[1] and "after end" in errors[2] and "neither a UniProt" in errors[3]
    assert "entry is empty" in errors[4] and "start is empty" in errors[5]


def test_a_repeated_domain_collapses_to_its_first_row(pipe, capsys):
    requests = pipe.build_requests(["A0A009F5K6", "A0A009F5K6", "A0A009F5K6"], [1, 1, 2], [10, 10, 10], 4)
    assert [r.id for r in requests] == ["A0A009F5K6_1-10", "A0A009F5K6_2-10"]
    assert "more than once" in capsys.readouterr().err


# ── the cut ───────────────────────────────────────────────────────────────────

def test_a_range_is_cut_renumbered_from_one_and_the_offset_recovers_the_original(pipe):
    atoms = real_atoms()
    cut = pipe.cut_domain(atoms, 5, 20, "title")
    assert (cut.length, cut.offset, cut.entry_length, cut.clipped) == (16, 4, 46, False)
    kept = [line for line in atoms if 5 <= int(line[22:26]) <= 20]
    written = atoms_of(cut.text)
    assert len(written) == len(kept)
    assert [int(line[6:11]) for line in written] == list(range(1, len(kept) + 1))
    for before, after in zip(kept, written):
        assert int(after[22:26]) + cut.offset == int(before[22:26])
        assert after[11:22] == before[11:22] and after[26:] == before[26:]
    assert int(written[0][22:26]) == 1 and int(written[-1][22:26]) == 16


def test_the_whole_entry_comes_back_unchanged_apart_from_the_title(pipe):
    lines = real_lines()
    cut = pipe.cut_domain(real_atoms(), 1, 46, "whole")
    assert cut.offset == 0 and cut.length == 46 and cut.clipped is False
    assert cut.text.splitlines() == ["TITLE     whole"] + lines[1:]


def test_the_terminator_follows_the_last_kept_residue(pipe):
    text = pipe.cut_domain(real_atoms(), 5, 20, "t").text.splitlines()
    n_atoms = len(atoms_of("\n".join(text)))
    assert text[0] == "TITLE     t"
    assert text[-1] == f"TER   {n_atoms + 1:5d}      {text[-2][17:20]} A  16"


def test_the_plddt_is_the_mean_over_residues_not_over_atoms(pipe):
    cut = pipe.cut_domain(real_atoms(), 5, 20, "t")
    per_residue = {int(line[22:26]): float(line[60:66]) for line in real_atoms() if 5 <= int(line[22:26]) <= 20}
    assert cut.plddt_mean == round(sum(per_residue.values()) / 16, 2)
    # Residues carry different numbers of atoms, so an atom mean would differ.
    atom_mean = sum(float(line[60:66]) for line in real_atoms() if 5 <= int(line[22:26]) <= 20) / len(
        [line for line in real_atoms() if 5 <= int(line[22:26]) <= 20])
    assert cut.plddt_mean != round(atom_mean, 2)


def test_a_range_past_the_end_is_clipped_and_flagged(pipe):
    cut = pipe.cut_domain(real_atoms(), 40, 100, "t")
    assert (cut.length, cut.offset, cut.clipped) == (7, 39, True)
    written = atoms_of(cut.text)
    assert int(written[-1][22:26]) == 7 and "OXT" in written[-1]


def test_a_range_outside_the_model_fails_and_names_its_length(pipe):
    with pytest.raises(ValueError, match=r"range 50-60 lies outside the entry \(residues 1-46, length 46\)"):
        pipe.cut_domain(real_atoms(), 50, 60, "t")


def test_a_fragment_keeps_its_own_numbering_in_the_offset_and_the_clip(pipe):
    atoms = atoms_of(synthetic_model(201, 300))
    inside = pipe.cut_domain(atoms, 250, 260, "t")
    assert (inside.length, inside.offset, inside.entry_length, inside.clipped) == (11, 249, 100, False)
    clipped = pipe.cut_domain(atoms, 150, 210, "t")
    assert (clipped.length, clipped.offset, clipped.clipped) == (10, 200, True)
    with pytest.raises(ValueError, match=r"residues 201-300, length 100"):
        pipe.cut_domain(atoms, 1, 200, "t")


def test_a_model_with_two_chains_is_refused(pipe):
    atoms = atoms_of(synthetic_model(1, 5, chain="A")) + atoms_of(synthetic_model(1, 5, chain="B"))
    with pytest.raises(ValueError, match="2 chains"):
        pipe.cut_domain(atoms, 1, 3, "t")


# ── command line and diagnostics ──────────────────────────────────────────────

def test_the_command_decompresses_every_listed_name_in_one_call(pipe):
    cmd = pipe.build_command("ids.txt", "/db/afdb_swissprot_v4", "out", 4)
    assert cmd == ["foldcomp", "decompress", "-t", "4", "--id-list", "ids.txt", "/db/afdb_swissprot_v4", "out"]


def test_the_container_prefix_leads_the_command(pipe):
    cmd = pipe.build_command("ids.txt", "db", "out", 1, "apptainer exec -B /data img.sif ")
    assert cmd[:5] == ["apptainer", "exec", "-B", "/data", "img.sif"] and cmd[5:7] == ["foldcomp", "decompress"]


def test_not_found_warnings_are_read_from_stderr(pipe):
    stderr = "[Warning] AF-A-F1-model_v4 not found in database.\nnoise\n[Warning] AF-B-F1-model_v4 not found in database.\n"
    assert pipe.parse_not_found(stderr) == {"AF-A-F1-model_v4", "AF-B-F1-model_v4"}


def test_an_absent_name_is_explained_by_the_warning_or_by_the_line_that_mentions_it(pipe):
    stderr = "[Error] processing db entry AF-C-F1-model_v4 failed.\n[Error] again AF-C-F1-model_v4\n"
    mentions = pipe.parse_mentions(stderr)
    assert list(mentions) == ["AF-C-F1-model_v4"] and mentions["AF-C-F1-model_v4"].startswith("[Error] processing")
    assert "is not in the Foldcomp database db" in pipe.explain_absence("AF-A-F1-model_v4", "/x/db", {"AF-A-F1-model_v4"}, {})
    assert "processing db entry" in pipe.explain_absence("AF-C-F1-model_v4", "/x/db", set(), mentions)
    assert pipe.explain_absence("AF-D-F1-model_v4", "/x/db", set(), {}) == "foldcomp wrote no structure for AF-D-F1-model_v4"


def test_a_prefix_without_its_parts_is_not_a_database(pipe, tmp_path):
    prefix = tmp_path / "db"
    prefix.write_text("data")
    with pytest.raises(RuntimeError, match=r"db\.dbtype.*db\.index.*db\.lookup|not a Foldcomp database"):
        pipe.check_database(str(prefix))
    for suffix in (".dbtype", ".index", ".lookup"):
        (tmp_path / f"db{suffix}").write_text("x")
    pipe.check_database(str(prefix))


# ── runtime half: stand-in binary ─────────────────────────────────────────────

# Follows src/main.cpp + src/input_processor.h of Foldcomp 1.0.0: names are looked up in <db>.lookup, a missing one only warns on stderr and the exit status stays 0, output is <out>/<name>.pdb.
STAND_IN = '''
import json, os, sys
args = sys.argv[1:]
assert args[0] == "foldcomp" and args[1] == "decompress", args
ids_file = args[args.index("--id-list") + 1]
positional = [a for i, a in enumerate(args[2:], 2) if not a.startswith("-") and args[i - 1] not in ("-t", "--id-list")]
database, out_dir = positional
with open(database + ".calls", "a") as log:
    log.write(json.dumps(args) + "\\n")
if os.environ.get("FOLDCOMP_STAND_IN_EXIT"):
    sys.stderr.write("[Error] stand-in failure\\n")
    sys.exit(int(os.environ["FOLDCOMP_STAND_IN_EXIT"]))
with open(ids_file) as f:
    names = [line.rstrip("\\n") for line in f]
with open(database + ".lookup") as f:
    known = {line.split("\\t")[1] for line in f}
with open(database) as f:
    stored = json.load(f)
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)
for name in names:
    if name not in known:
        sys.stderr.write("[Warning] " + name + " not found in database.\\n")
        continue
    with open(os.path.join(out_dir, name + ".pdb"), "w") as f:
        f.write(stored[name])
'''


def make_database(tmp_path, entries):
    """A tiny database folder in Foldcomp's layout: data, .index, .dbtype and .lookup."""
    folder = tmp_path / "database"
    folder.mkdir(exist_ok=True)
    prefix = folder / "tiny_v4"
    prefix.write_text(json.dumps(entries))
    (folder / "tiny_v4.index").write_text("".join(f"{k}\t{k * 10}\t100\n" for k, _ in enumerate(entries)))
    (folder / "tiny_v4.dbtype").write_bytes(bytes([12, 0, 0, 0]))
    (folder / "tiny_v4.lookup").write_text("".join(f"{k}\t{name}\t0\n" for k, name in enumerate(entries)))
    return str(prefix)


def stand_in_prefix(tmp_path):
    script = tmp_path / "stand_in.py"
    script.write_text(STAND_IN)
    return f'"{pathlib.Path(sys.executable).as_posix()}" "{script.as_posix()}" '


CRAMBIN = "".join(line + "\n" for line in real_lines())
DATABASE_ENTRIES = {
    "AF-A0A009F5K6-F1-model_v4": CRAMBIN,
    "AF-A0A009E921-F1-model_v4": synthetic_model(1, 120, b_of=lambda n: 50.0 + n % 40),
    "AF-A0A009ECR5-F2-model_v4": synthetic_model(201, 300, b_of=lambda n: 90.0),
}
MEMBERS = {
    "uniprot": ["A0A009F5K6", "A0A009F5K6", "A0A009F5K6", "A0A009E921", "A0A009E9H3", "A0A009F5K6",
                "AF-A0A009ECR5-F2-model_v4", "bad entry!", "A0A009E921"],
    "structure_start": [5, 40, 60, 1, 1, 5, 210, 1, 0],
    "structure_end": [20, 100, 80, 120, 10, 20, 230, 5, 5],
}


def run_main(pipe, monkeypatch, tmp_path, members, database, container_prefix=None, extra_cfg=None):
    members_csv = tmp_path / "members.csv"
    pd.DataFrame(members).to_csv(members_csv, index=False)
    cfg = {"entries": f"TABLE_REFERENCE:{members_csv}:uniprot",
           "start": f"TABLE_REFERENCE:{members_csv}:structure_start",
           "end": f"TABLE_REFERENCE:{members_csv}:structure_end",
           "database": database, "model_version": 4, "threads": 2}
    cfg.update(extra_cfg or {})
    config = tmp_path / "config.yaml"
    config.write_text(yaml.safe_dump(cfg))
    step = tmp_path / "job" / "004_Foldcomp"
    argv = ["pipe_foldcomp.py", "--config-yaml", str(config),
            "--structures-dir", str(step / "structures"),
            "--structures-map-csv", str(step / "structures" / "structures.csv"),
            "--domains-csv", str(step / "tables" / "domains.csv"),
            "--scratch-dir", str(step / "_execution"),
            "--local-missing-csv", str(step / "_execution" / "local_missing.csv"),
            "--container-prefix", stand_in_prefix(tmp_path) if container_prefix is None else container_prefix]
    monkeypatch.setattr(sys, "argv", argv)
    pipe.main()


def read_outputs(tmp_path):
    step = tmp_path / "job" / "004_Foldcomp"
    return (pd.read_csv(step / "tables" / "domains.csv"),
            pd.read_csv(step / "structures" / "structures.csv"),
            pd.read_csv(step / "_execution" / "local_missing.csv"),
            step)


def test_a_mixed_table_yields_one_row_per_domain_and_every_failure_kind(pipe, monkeypatch, tmp_path):
    run_main(pipe, monkeypatch, tmp_path, MEMBERS, make_database(tmp_path, DATABASE_ENTRIES))
    domains, structures, missing, step = read_outputs(tmp_path)

    assert list(domains.columns) == ["id", "entry", "database_name", "start", "end", "length", "offset",
                                     "entry_length", "plddt_mean", "found", "clipped"]
    assert domains["id"].tolist() == ["A0A009F5K6_5-20", "A0A009F5K6_40-100", "A0A009F5K6_60-80", "A0A009E921_1-120",
                                      "A0A009E9H3_1-10", "AF-A0A009ECR5-F2-model_v4_210-230", "bad entry!_1-5",
                                      "A0A009E921_0-5"]
    by_id = domains.set_index("id")
    assert by_id.loc["A0A009F5K6_5-20", ["length", "offset", "entry_length", "found", "clipped"]].tolist() == [16, 4, 46, True, False]
    assert by_id.loc["A0A009F5K6_40-100", ["length", "offset", "clipped"]].tolist() == [7, 39, True]
    assert by_id.loc["A0A009E921_1-120", ["length", "offset", "entry_length"]].tolist() == [120, 0, 120]
    fragment = by_id.loc["AF-A0A009ECR5-F2-model_v4_210-230"]
    assert (fragment["length"], fragment["offset"], fragment["entry_length"], fragment["plddt_mean"]) == (21, 209, 100, 90.0)
    assert by_id.loc["A0A009F5K6_5-20", "database_name"] == "AF-A0A009F5K6-F1-model_v4"

    failed = by_id[~by_id["found"]]
    assert sorted(failed.index) == ["A0A009E921_0-5", "A0A009E9H3_1-10", "A0A009F5K6_60-80", "bad entry!_1-5"]
    assert failed[["length", "offset", "plddt_mean", "clipped"]].isna().all().all()
    assert failed.loc["A0A009F5K6_60-80", "entry_length"] != failed.loc["A0A009F5K6_60-80", "entry_length"]

    assert structures["id"].tolist() == ["A0A009F5K6_5-20", "A0A009F5K6_40-100", "A0A009E921_1-120",
                                         "AF-A0A009ECR5-F2-model_v4_210-230"]
    assert all(pathlib.Path(p).exists() for p in structures["file"])
    assert sorted(p.name for p in (step / "structures").glob("*.pdb")) == sorted(f"{i}.pdb" for i in structures["id"])

    assert missing["kind"].unique().tolist() == ["failure"] and (missing["removed_by"] == "004_Foldcomp").all()
    cause = missing.set_index("id")["cause"]
    assert "is not in the Foldcomp database tiny_v4" in cause["A0A009E9H3_1-10"]
    assert "residues 1-46, length 46" in cause["A0A009F5K6_60-80"]
    assert "at least 1" in cause["A0A009E921_0-5"] and "neither a UniProt" in cause["bad entry!_1-5"]


def test_the_database_is_read_in_one_pass_over_each_distinct_valid_entry(pipe, monkeypatch, tmp_path):
    database = make_database(tmp_path, DATABASE_ENTRIES)
    run_main(pipe, monkeypatch, tmp_path, MEMBERS, database)
    calls = [json.loads(line) for line in open(database + ".calls")]
    assert len(calls) == 1 and calls[0][2:4] == ["-t", "2"]
    names = open(tmp_path / "job" / "004_Foldcomp" / "_execution" / "foldcomp_ids.txt").read().split()
    assert names == ["AF-A0A009F5K6-F1-model_v4", "AF-A0A009E921-F1-model_v4", "AF-A0A009E9H3-F1-model_v4",
                     "AF-A0A009ECR5-F2-model_v4"]


def test_a_written_domain_is_the_renumbered_cut_of_the_stored_model(pipe, monkeypatch, tmp_path):
    run_main(pipe, monkeypatch, tmp_path, MEMBERS, make_database(tmp_path, DATABASE_ENTRIES))
    _, _, _, step = read_outputs(tmp_path)
    text = (step / "structures" / "A0A009F5K6_5-20.pdb").read_text()
    assert text.splitlines()[0] == "TITLE     AF-A0A009F5K6-F1-model_v4 residues 5-20"
    original = next(line for line in real_atoms() if int(line[22:26]) == 5)
    assert atoms_of(text)[0] == f"{original[:6]}    1{original[11:22]}   1{original[26:]}"
    assert text.splitlines()[-1].startswith("TER")


def test_the_decompressed_entries_are_not_kept_but_the_ids_and_the_log_are(pipe, monkeypatch, tmp_path):
    run_main(pipe, monkeypatch, tmp_path, MEMBERS, make_database(tmp_path, DATABASE_ENTRIES))
    scratch = tmp_path / "job" / "004_Foldcomp" / "_execution"
    assert (scratch / "foldcomp_ids.txt").exists() and (scratch / "foldcomp.log").exists()
    assert not (scratch / "foldcomp_out").exists()
    assert "AF-A0A009E9H3-F1-model_v4 not found in database" in (scratch / "foldcomp.log").read_text()


def test_every_domain_failing_keeps_the_complete_table_and_fails_the_step(pipe, monkeypatch, tmp_path):
    members = {"uniprot": ["A0A009E9H3", "A0A009F5K6"], "structure_start": [1, 90], "structure_end": [5, 99]}
    with pytest.raises(SystemExit) as raised:
        run_main(pipe, monkeypatch, tmp_path, members, make_database(tmp_path, DATABASE_ENTRIES))
    assert raised.value.code == 1
    domains, structures, missing, _ = read_outputs(tmp_path)
    assert domains["id"].tolist() == ["A0A009E9H3_1-5", "A0A009F5K6_90-99"] and not domains["found"].any()
    assert structures.empty and sorted(missing["id"]) == sorted(domains["id"])


def test_a_prefix_that_is_not_a_database_fails_every_domain_with_the_reason(pipe, monkeypatch, tmp_path):
    with pytest.raises(SystemExit) as raised:
        run_main(pipe, monkeypatch, tmp_path, MEMBERS, str(tmp_path / "nowhere" / "db"))
    assert raised.value.code == 1
    domains, structures, missing, _ = read_outputs(tmp_path)
    assert len(domains) == 8 and structures.empty
    assert "not a Foldcomp database" in missing.set_index("id").loc["A0A009F5K6_5-20", "cause"]


def test_a_foldcomp_crash_fails_every_domain_with_its_exit_status(pipe, monkeypatch, tmp_path):
    monkeypatch.setenv("FOLDCOMP_STAND_IN_EXIT", "3")
    with pytest.raises(SystemExit):
        run_main(pipe, monkeypatch, tmp_path, MEMBERS, make_database(tmp_path, DATABASE_ENTRIES))
    _, structures, missing, step = read_outputs(tmp_path)
    assert structures.empty
    assert "foldcomp exited 3" in missing.set_index("id").loc["A0A009F5K6_5-20", "cause"]
    assert "foldcomp exited 3" in (step / "_execution" / "foldcomp.log").read_text()


def test_a_table_with_no_rows_or_a_missing_column_stops_the_step(pipe, monkeypatch, tmp_path):
    database = make_database(tmp_path, DATABASE_ENTRIES)
    with pytest.raises(SystemExit) as raised:
        run_main(pipe, monkeypatch, tmp_path, {"uniprot": [], "structure_start": [], "structure_end": []}, database)
    assert raised.value.code == 1
    with pytest.raises(SystemExit) as raised:
        run_main(pipe, monkeypatch, tmp_path, {"uniprot": ["A0A009F5K6"], "structure_start": [1]}, database)
    assert raised.value.code == 1


def test_entries_that_are_all_invalid_never_call_foldcomp(pipe, monkeypatch, tmp_path):
    database = make_database(tmp_path, DATABASE_ENTRIES)
    with pytest.raises(SystemExit):
        run_main(pipe, monkeypatch, tmp_path, {"uniprot": ["nope"], "structure_start": [1], "structure_end": [5]}, database)
    assert not os.path.exists(database + ".calls")


# ── wrapper ───────────────────────────────────────────────────────────────────

def _members(columns=("uniprot", "structure_start", "structure_end"), name="members"):
    from biopipelines.mock import Mock
    return Mock(ids=["t"], tables={name: {"columns": list(columns)}})


def _build(new_pipeline, job, **overrides):
    from biopipelines.foldcomp import Foldcomp

    pipeline = new_pipeline(job)
    with pipeline:
        members = _members()
        kwargs = dict(entries=members.tables.members.uniprot, start=members.tables.members.structure_start,
                      end=members.tables.members.structure_end, database="/data/foldcomp/afdb_rep_v4")
        kwargs.update(overrides)
        out = Foldcomp(**kwargs)
    return pipeline, out


def test_defaults_and_declared_outputs(local_config, isolated_cwd, new_pipeline):
    pipeline, out = _build(new_pipeline, "fc_defaults")
    tool = pipeline.tools[-1]
    assert (tool.model_version, tool.threads) == (4, 1)
    assert out.tables.domains.info.columns == ["id", "entry", "database_name", "start", "end", "length", "offset",
                                               "entry_length", "plddt_mean", "found", "clipped"]
    assert out.tables.missing.info.columns == ["id", "removed_by", "kind", "cause"]
    stream = out.streams.structures
    assert stream.format == "pdb"


def test_the_structures_stream_is_lazy_and_covers_ids_with_underscores():
    from biopipelines import id_patterns
    from biopipelines.foldcomp import DOMAIN_ID_PATTERN

    ids = ["A0A009F5K6_12-150", "AF-Q8W3K0-F2-model_v4_201-300"]
    assert id_patterns.is_lazy(DOMAIN_ID_PATTERN)
    assert id_patterns.select_ids([DOMAIN_ID_PATTERN], ids) == ids


def test_a_path_and_column_tuple_names_the_entries(local_config, isolated_cwd, new_pipeline):
    from biopipelines.foldcomp import Foldcomp

    with new_pipeline("fc_tuple"):
        out = Foldcomp(entries=("members.csv", "uniprot"), start=("members.csv", "s"), end=("members.csv", "e"), database="db")
    assert out is not None


@pytest.mark.parametrize("overrides, message", [
    ({"entries": "A0A009F5K6"}, "must be a column of a table"),
    ({"database": ""}, "path prefix"),
    ({"model_version": 0}, "model_version"),
    ({"threads": 0}, "threads"),
])
def test_bad_parameters_are_refused_at_configuration(local_config, isolated_cwd, new_pipeline, overrides, message):
    with pytest.raises(ValueError, match=message):
        _build(new_pipeline, "fc_bad", **overrides)


@pytest.mark.parametrize("database", ["/data/$(touch pwned)/afdb", "/data/`id`/afdb", '/data/"x/afdb', "C:\\data\\afdb"])
def test_a_database_path_that_would_be_expanded_or_would_break_the_script_is_refused(
        local_config, isolated_cwd, new_pipeline, database):
    """The path is echoed inside double quotes in the generated bash, where $(...), a backtick and a quote act."""
    with pytest.raises(ValueError, match="'database'.*would break the generated shell script"):
        _build(new_pipeline, "fc_unsafe_database", database=database)


def test_an_ordinary_database_path_still_builds_and_its_script_echoes_it_plainly(local_config, isolated_cwd, new_pipeline):
    pipeline, _ = _build(new_pipeline, "fc_safe_database", database="/data/foldcomp-v4/afdb_rep_v4")
    assert 'echo "Extracting domains from afdb_rep_v4"' in _step_text(pipeline.save())


def test_the_domain_id_is_the_multiplied_id_of_the_entry_and_its_range_and_the_unnamed_row_keeps_its_name(pipe):
    from biopipelines.combinatorics import generate_multiplied_ids

    requests = pipe.build_requests(["A0A009F5K6", ""], [12, 3], [150, 9], 4)
    assert requests[0].id == generate_multiplied_ids(["A0A009F5K6"], ["12-150"])[0][0] == "A0A009F5K6_12-150"
    assert requests[1].id == "row1"
    assert "generate_multiplied_ids" in open(pipe.__file__, encoding="utf-8").read()


def test_columns_of_different_tables_are_refused(local_config, isolated_cwd, new_pipeline):
    from biopipelines.foldcomp import Foldcomp

    with new_pipeline("fc_two_tables"):
        a, b = _members(name="first"), _members(name="second")
        with pytest.raises(ValueError, match="same table"):
            Foldcomp(entries=a.tables.first.uniprot, start=b.tables.second.structure_start,
                     end=a.tables.first.structure_end, database="db")


def test_one_column_serving_two_roles_is_refused(local_config, isolated_cwd, new_pipeline):
    from biopipelines.foldcomp import Foldcomp

    with new_pipeline("fc_same_column"):
        m = _members()
        with pytest.raises(ValueError, match="three different columns"):
            Foldcomp(entries=m.tables.members.uniprot, start=m.tables.members.structure_start,
                     end=m.tables.members.structure_start, database="db")


def _step_text(script_path):
    step = sorted(glob.glob(os.path.join(os.path.dirname(script_path), "*_Foldcomp.sh")))[0]
    return open(step, encoding="utf-8", errors="replace").read()


def test_the_emitted_script_runs_the_pipe_and_writes_the_parameters(local_config, isolated_cwd, new_pipeline):
    pipeline, out = _build(new_pipeline, "fc_script", threads=4, model_version=6)
    script_path = pipeline.save()
    text = _step_text(script_path)
    assert "pipe_foldcomp.py" in text and "--container-prefix" in text and "--local-missing-csv" in text
    assert "pipe_propagate_missing.py" in text
    config = yaml.safe_load(open(pipeline.tools[-1].config_yaml))
    assert config["threads"] == 4 and config["model_version"] == 6 and config["database"] == "/data/foldcomp/afdb_rep_v4"
    for key, column in (("entries", "uniprot"), ("start", "structure_start"), ("end", "structure_end")):
        assert config[key].startswith("TABLE_REFERENCE:") and config[key].endswith(f":{column}")
    assert len({config[k].rsplit(":", 1)[0] for k in ("entries", "start", "end")}) == 1


def test_every_pipe_flag_has_a_matching_argument(local_config, isolated_cwd, new_pipeline):
    pipeline, _ = _build(new_pipeline, "fc_flags")
    text = _step_text(pipeline.save())
    call = text[text.index("pipe_foldcomp.py"):]
    call = call[:call.index("pipe_propagate_missing.py")]
    emitted = set(re.findall(r"--[a-z][a-z-]+", call))
    declared = set(re.findall(r'add_argument\("(--[a-z-]+)"', open(SCRIPT, encoding="utf-8").read()))
    assert emitted == declared


# ── install ───────────────────────────────────────────────────────────────────

@pytest.fixture
def cluster_config(monkeypatch):
    """The shipped config.cluster.yaml, whose `environments:` block names the env the install builds."""
    from biopipelines.config_manager import ConfigManager

    def _reset():
        ConfigManager._instance = None
        ConfigManager._config = None
        ConfigManager._variant = None

    _reset()
    path = os.path.join(ROOT, "config.cluster.yaml")
    monkeypatch.setattr(ConfigManager, "_get_config_path", classmethod(lambda cls, variant=None: path))
    yield ConfigManager(variant="cluster")
    _reset()


def test_the_install_builds_the_env_and_verifies_the_binary(cluster_config):
    from biopipelines.foldcomp import Foldcomp

    script = Foldcomp._install_script({"biopipelines": "/bp"}, "mamba")
    assert "foldcomp --version" in script and "INSTALL_SUCCESS" in script
    assert "afdb_swissprot_v4" not in script


def test_the_swissprot_download_is_opt_in_and_fetches_the_four_database_files(cluster_config):
    from biopipelines.foldcomp import Foldcomp

    script = Foldcomp._install_script({"biopipelines": "/bp"}, "mamba", swissprot_dir="/data/foldcomp")
    for suffix in ("", ".index", ".dbtype", ".lookup"):
        assert f'"{suffix}"' in script
    assert "https://opendata.mmseqs.org/foldcomp/afdb_swissprot_v4$SUFFIX" in script
    assert 'FILE="/data/foldcomp/afdb_swissprot_v4$SUFFIX"' in script
    assert "curl -fL -C -" in script


def test_a_download_folder_that_would_break_the_shell_is_refused(cluster_config):
    from biopipelines.foldcomp import Foldcomp

    with pytest.raises(ValueError, match="swissprot_dir"):
        Foldcomp._install_script({"biopipelines": "/bp"}, "mamba", swissprot_dir='/data/$(rm -rf x)')


def test_the_environment_entry_exists_in_the_variants_it_runs_under():
    for variant in ("cluster", "colab", "container"):
        config = yaml.safe_load(open(os.path.join(ROOT, f"config.{variant}.yaml"), encoding="utf-8"))
        assert config["environments"]["Foldcomp"] == "foldcomp"
    assert os.path.exists(os.path.join(ROOT, "environments", "foldcomp.yaml"))
