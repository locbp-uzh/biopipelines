"""IPSAE wraps the upstream ipsae.py, so what the wrapper owns is staging, parsing, the NaN rules and the per-id summary.

Upstream prints a score table and exits 0 even when an input is unusable, and it prints pDockQ and ipTM computed from zeros when pLDDT or the predictor's ipTM is absent; the tests pin how each of those is reported.
"""

import importlib.util
import json
import math
import os
import sys
import urllib.request

import numpy as np
import pytest

pytest.importorskip("yaml")

ROOT = os.path.join(os.path.dirname(__file__), "..")
SCRIPT = os.path.join(ROOT, "pipe_scripts", "pipe_ipsae.py")
UPSTREAM_OUTPUT = os.path.join(os.path.dirname(__file__), "fixtures", "ipsae_upstream_output.txt")
UPSTREAM_URL = "https://raw.githubusercontent.com/DunbrackLab/IPSAE/main/ipsae.py"

N_A, N_B = 12, 10


@pytest.fixture(scope="module")
def pipe():
    spec = importlib.util.spec_from_file_location("pipe_ipsae", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_pdb(path, chains=(("A", N_A, 0.0), ("B", N_B, 6.0))):
    """Two ALA chains, CA and CB per residue, chain B laid 3 A beside chain A."""
    lines, serial = [], 1
    for chain, length, y in chains:
        for i in range(length):
            for name, dy in (("CA", 0.0), ("CB", -1.5 if y else 1.5)):
                lines.append(f"ATOM  {serial:5d} {name:^4s} ALA {chain}{i + 1:4d}    "
                             f"{3.8 * i:8.3f}{y + dy:8.3f}{0.0:8.3f}  1.00 80.00           C")
                serial += 1
    path.write_text("\n".join(lines) + "\nEND\n")
    return str(path)


def write_cif(path, chains=(("A", N_A, 0.0), ("B", N_B, 6.0))):
    """The same two chains as mmCIF, with the atom_site fields upstream reads."""
    fields = ["group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id", "label_seq_id",
              "auth_seq_id", "pdbx_PDB_ins_code", "label_asym_id", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy",
              "label_entity_id", "auth_asym_id", "auth_comp_id", "B_iso_or_equiv", "pdbx_PDB_model_num"]
    lines = ["data_model", "loop_"] + [f"_atom_site.{f}" for f in fields]
    serial = 1
    for chain, length, y in chains:
        for i in range(length):
            for name, dy in (("CA", 0.0), ("CB", -1.5 if y else 1.5)):
                lines.append(f"ATOM {serial} C {name} . ALA {i + 1} {i + 1} ? {chain} "
                             f"{3.8 * i:.3f} {y + dy:.3f} 0.000 1 1 {chain} ALA 80.00 1")
                serial += 1
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def pae_matrix():
    """Intra-chain 1 A; A-aligned error on B is 2 A, B-aligned error on A is 6 A, so the two directions differ."""
    pae = np.full((N_A + N_B, N_A + N_B), 1.0)
    pae[:N_A, N_A:] = 2.0
    pae[N_A:, :N_A] = 6.0
    return pae


def write_npz(path):
    np.savez(path, pae=pae_matrix())
    return str(path)


def write_json(path, **extra):
    path.write_text(json.dumps({"pae": pae_matrix().tolist(), **extra}))
    return str(path)


def write_stream(path, name, ids, files, fmt):
    path.write_text(json.dumps({"name": name, "ids": ids, "files": files, "map_table": "", "format": fmt}))
    return str(path)


def row(chain1, chain2, type_="max", ipsae=0.5, **over):
    base = {"chain1": chain1, "chain2": chain2, "type": type_, "ipsae": ipsae, "ipsae_d0chn": 0.4,
            "ipsae_d0dom": 0.45, "iptm_af": 0.7, "iptm_d0chn": 0.3, "pdockq": 0.2, "pdockq2": 0.1, "lis": 0.25}
    base.update(over)
    return base


def trio(chain1, chain2, ipsae):
    """Both directions and the max, as upstream prints a chain pair."""
    return [row(chain1, chain2, "asym", ipsae - 0.1), row(chain2, chain1, "asym", ipsae),
            row(chain1, chain2, "max", ipsae)]


# ── parsing what upstream wrote ──────────────────────────────────────────────

def test_the_upstream_score_table_parses_into_typed_columns(pipe):
    rows = pipe.parse_ipsae_output(open(UPSTREAM_OUTPUT).read())
    assert [(r["chain1"], r["chain2"], r["type"]) for r in rows] == [("A", "B", "asym"), ("B", "A", "asym"), ("A", "B", "max")]
    first = rows[0]
    assert set(first) == set(pipe.PAIR_COLUMNS) - {"id", "pae_cutoff", "dist_cutoff"}
    assert isinstance(first["ipsae"], float) and isinstance(first["n0res"], int) and isinstance(first["nres1"], int)
    assert rows[2]["ipsae"] == max(rows[0]["ipsae"], rows[1]["ipsae"])


def test_a_directional_pair_scores_differently_in_each_direction(pipe):
    """The fixture PAE is asymmetric on purpose; a parser that merged the directions would lose this."""
    asym_ab, asym_ba, _ = pipe.parse_ipsae_output(open(UPSTREAM_OUTPUT).read())
    assert asym_ab["ipsae"] != asym_ba["ipsae"]


def test_output_without_a_score_table_is_an_error(pipe):
    with pytest.raises(ValueError, match="no score table"):
        pipe.parse_ipsae_output("")


def test_a_renamed_upstream_column_is_refused_not_guessed(pipe):
    text = open(UPSTREAM_OUTPUT).read().replace("pDockQ2", "pDockQ_2")
    with pytest.raises(ValueError, match="pDockQ2"):
        pipe.parse_ipsae_output(text)


def test_a_row_with_the_wrong_field_count_is_an_error(pipe):
    text = open(UPSTREAM_OUTPUT).read().rstrip() + "\nA B 10\n"
    with pytest.raises(ValueError, match="fields"):
        pipe.parse_ipsae_output(text)


def test_a_single_chain_structure_parses_to_no_rows(pipe):
    """Upstream still writes the header when there is no chain pair."""
    header = [l for l in open(UPSTREAM_OUTPUT).read().splitlines() if l.startswith("Chn1")][0]
    assert pipe.parse_ipsae_output("\n" + header + "\n") == []


# ── what the PAE file does not carry ─────────────────────────────────────────

def test_an_npz_carries_neither_plddt_nor_iptm(pipe, tmp_path):
    assert pipe.confidence_available(write_npz(tmp_path / "p.npz"), ".pdb") == (False, False)


def test_a_colabfold_json_is_checked_key_by_key(pipe, tmp_path):
    assert pipe.confidence_available(write_json(tmp_path / "a.json", plddt=[90.0], iptm=0.8), ".pdb") == (True, True)
    assert pipe.confidence_available(write_json(tmp_path / "b.json"), ".pdb") == (False, False)


def test_an_alphafold3_json_has_plddt_but_no_staged_summary(pipe, tmp_path):
    assert pipe.confidence_available(write_json(tmp_path / "c.json", atom_plddts=[90.0]), ".cif") == (True, False)


def test_unavailable_columns_become_nan_and_the_rest_are_untouched(pipe):
    masked = pipe.mask_unavailable([row("A", "B")], has_plddt=False, has_iptm=False)[0]
    assert math.isnan(masked["pdockq"]) and math.isnan(masked["pdockq2"]) and math.isnan(masked["iptm_af"])
    assert masked["ipsae"] == 0.5 and masked["lis"] == 0.25
    kept = pipe.mask_unavailable([row("A", "B")], has_plddt=True, has_iptm=True)[0]
    assert kept["pdockq"] == 0.2 and kept["iptm_af"] == 0.7


# ── the per-id summary ───────────────────────────────────────────────────────

def test_the_summary_is_the_max_row_of_the_best_target(pipe):
    rows = trio("A", "B", 0.3) + trio("A", "C", 0.8) + trio("B", "C", 0.9)
    summary = pipe.summarize(rows, "A", None)
    assert summary["target_chain"] == "C" and summary["ipsae"] == 0.8 and summary["binder_chain"] == "A"
    assert set(summary) == set(pipe.SUMMARY_COLUMNS) - {"id"}


def test_the_binder_is_found_whichever_side_upstream_printed_it_on(pipe):
    rows = trio("A", "B", 0.4) + trio("B", "C", 0.9)
    assert pipe.summarize(rows, "C", ["B"])["ipsae"] == 0.9


def test_target_chains_restrict_the_choice(pipe):
    rows = trio("A", "B", 0.3) + trio("A", "C", 0.8)
    assert pipe.summarize(rows, "A", ["B"])["target_chain"] == "B"


def test_a_tie_goes_to_the_first_target_listed(pipe):
    rows = trio("A", "B", 0.5) + trio("A", "C", 0.5)
    assert pipe.summarize(rows, "A", ["C", "B"])["target_chain"] == "C"


def test_a_missing_binder_chain_names_the_chains_present(pipe):
    with pytest.raises(ValueError, match="binder_chain 'Z'.*A, B"):
        pipe.summarize(trio("A", "B", 0.5), "Z", None)


def test_a_missing_target_chain_is_an_error(pipe):
    with pytest.raises(ValueError, match="target chain"):
        pipe.summarize(trio("A", "B", 0.5), "A", ["Q"])


def test_a_structure_with_only_the_binder_chain_is_an_error(pipe):
    with pytest.raises(ValueError, match="binder_chain 'A'"):
        pipe.summarize([], "A", None)


# ── staging and the step loop, against a stand-in for upstream ───────────────

STAND_IN = '''
import json, sys
pae_file, structure_file, pae_cutoff, dist_cutoff = sys.argv[1:5]
mode = json.load(open(pae_file)).get("mode") if pae_file.endswith(".json") else "ok"
stem = structure_file.rsplit(".", 1)[0]
if mode == "crash":
    sys.exit(3)
open(f"{stem}_10_10.txt", "w").close()
if mode == "silent":
    print("AF2 PAE file does not exist: nowhere")
    sys.exit()
open(f"{stem}_10_10.txt", "w").write(open(%r).read())
'''


@pytest.fixture
def stand_in(tmp_path):
    path = tmp_path / "ipsae.py"
    path.write_text(STAND_IN % UPSTREAM_OUTPUT)
    return str(path)


def run_main(pipe, monkeypatch, tmp_path, stand_in, ids, pae_files, binder="A", targets=None, structure_ext="pdb"):
    structures = [write_pdb(tmp_path / f"{i}.{structure_ext}") for i in ids]
    cfg = tmp_path / "config.yaml"
    cfg.write_text(json.dumps({"binder_chain": binder, "target_chains": targets, "pae_cutoff": 10.0, "dist_cutoff": 10.0}))
    out = tmp_path / "out"
    argv = ["pipe_ipsae.py", "--config-yaml", str(cfg),
            "--structures-json", write_stream(tmp_path / "s.json", "structures", ids, structures, structure_ext),
            "--pae-json", write_stream(tmp_path / "p.json", "pae", ids, pae_files, "json"),
            "--ipsae-script", stand_in, "--runs-dir", str(out / "runs"),
            "--ipsae-csv", str(out / "tables" / "ipsae.csv"), "--pairs-csv", str(out / "tables" / "ipsae_pairs.csv"),
            "--local-missing-csv", str(out / "_execution" / "local_missing.csv")]
    monkeypatch.setattr(sys, "argv", argv)
    return out


def test_one_row_per_id_with_a_failed_id_as_nan_and_a_missing_entry(pipe, monkeypatch, tmp_path, stand_in):
    import pandas as pd

    pae = [write_json(tmp_path / "a.json", mode="ok"), write_json(tmp_path / "b.json", mode="crash"),
           write_json(tmp_path / "c.json", mode="silent")]
    out = run_main(pipe, monkeypatch, tmp_path, stand_in, ["a", "b", "c"], pae)
    pipe.main()

    summary = pd.read_csv(out / "tables" / "ipsae.csv")
    assert list(summary.columns) == pipe.SUMMARY_COLUMNS
    assert list(summary["id"]) == ["a", "b", "c"]
    assert summary.loc[0, "target_chain"] == "B" and not math.isnan(summary.loc[0, "ipsae"])
    assert summary.loc[1:, "ipsae"].isna().all() and summary.loc[1:, "target_chain"].isna().all()
    assert list(summary["binder_chain"]) == ["A", "A", "A"]

    pairs = pd.read_csv(out / "tables" / "ipsae_pairs.csv")
    assert list(pairs.columns) == pipe.PAIR_COLUMNS
    assert set(pairs["id"]) == {"a"} and len(pairs) == 3
    assert (pairs["pae_cutoff"] == 10.0).all() and (pairs["dist_cutoff"] == 10.0).all()

    missing = pd.read_csv(out / "_execution" / "local_missing.csv")
    assert list(missing["id"]) == ["b", "c"] and set(missing["kind"]) == {"failure"}
    assert "exited 3" in missing.loc[0, "cause"]
    assert "AF2 PAE file does not exist" in missing.loc[1, "cause"]


def test_pdockq_and_iptm_are_nan_when_the_json_carries_no_plddt_or_iptm(pipe, monkeypatch, tmp_path, stand_in):
    import pandas as pd

    out = run_main(pipe, monkeypatch, tmp_path, stand_in, ["a"], [write_json(tmp_path / "a.json", mode="ok")])
    pipe.main()
    summary = pd.read_csv(out / "tables" / "ipsae.csv")
    assert summary[["pdockq", "pdockq2", "iptm_af"]].isna().all().all()
    assert not summary[["ipsae", "ipsae_d0chn", "lis"]].isna().any().any()


def test_pdockq_and_iptm_survive_when_the_json_carries_them(pipe, monkeypatch, tmp_path, stand_in):
    import pandas as pd

    pae = [write_json(tmp_path / "a.json", mode="ok", plddt=[90.0] * (N_A + N_B), iptm=0.8)]
    out = run_main(pipe, monkeypatch, tmp_path, stand_in, ["a"], pae)
    pipe.main()
    summary = pd.read_csv(out / "tables" / "ipsae.csv")
    assert not summary[["pdockq", "pdockq2", "iptm_af"]].isna().any().any()


def test_the_upstream_outputs_are_kept_and_the_staged_inputs_removed(pipe, monkeypatch, tmp_path, stand_in):
    out = run_main(pipe, monkeypatch, tmp_path, stand_in, ["a"], [write_json(tmp_path / "a.json", mode="ok")])
    pipe.main()
    assert sorted(os.listdir(out / "runs" / "a")) == ["model_10_10.txt"]


def test_every_id_failing_is_a_step_failure_but_the_tables_are_still_written(pipe, monkeypatch, tmp_path, stand_in):
    out = run_main(pipe, monkeypatch, tmp_path, stand_in, ["a"], [write_json(tmp_path / "a.json", mode="crash")])
    with pytest.raises(SystemExit) as raised:
        pipe.main()
    assert raised.value.code == 1
    assert os.path.exists(out / "tables" / "ipsae.csv") and os.path.exists(out / "_execution" / "local_missing.csv")


def test_no_declared_id_having_a_file_is_a_step_failure(pipe, monkeypatch, tmp_path, stand_in):
    out = run_main(pipe, monkeypatch, tmp_path, stand_in, ["a"], [write_json(tmp_path / "a.json", mode="ok")])
    os.remove(tmp_path / "a.pdb")
    with pytest.raises(SystemExit) as raised:
        pipe.main()
    assert raised.value.code == 1


def test_an_id_absent_from_the_pae_stream_is_a_failure_naming_it(pipe, monkeypatch, tmp_path, stand_in):
    import pandas as pd

    out = run_main(pipe, monkeypatch, tmp_path, stand_in, ["a", "b"], [write_json(tmp_path / "a.json", mode="ok")] * 2)
    with open(tmp_path / "p.json", "w") as fh:
        json.dump({"name": "pae", "ids": ["a"], "files": [str(tmp_path / "a.json")], "map_table": "", "format": "json"}, fh)
    pipe.main()
    missing = pd.read_csv(out / "_execution" / "local_missing.csv")
    assert list(missing["id"]) == ["b"] and "'b'" in missing.loc[0, "cause"]


def test_a_structure_that_is_not_pdb_or_cif_is_a_failure(pipe, monkeypatch, tmp_path, stand_in):
    import pandas as pd

    pae = [write_json(tmp_path / "a.json", mode="ok"), write_json(tmp_path / "b.json", mode="ok")]
    out = run_main(pipe, monkeypatch, tmp_path, stand_in, ["a", "b"], pae)
    structures = [str(tmp_path / "a.pdb"), str(tmp_path / "b.xyz")]
    os.rename(tmp_path / "b.pdb", tmp_path / "b.xyz")
    write_stream(tmp_path / "s.json", "structures", ["a", "b"], structures, "pdb")
    pipe.main()
    missing = pd.read_csv(out / "_execution" / "local_missing.csv")
    assert list(missing["id"]) == ["b"] and "neither PDB nor mmCIF" in missing.loc[0, "cause"]


# ── the wrapper ──────────────────────────────────────────────────────────────

def predictor(Mock, pae_format="npz", with_pae=True):
    streams = {"structures": {"format": "pdb", "file": "<id>.pdb"}}
    if with_pae:
        streams["pae"] = {"format": pae_format, "file": "<id>." + pae_format}
    return Mock(ids=["s1", "s2"], streams=streams)


def build(new_pipeline, **overrides):
    from biopipelines.ipsae import IPSAE
    from biopipelines.mock import Mock

    pipeline = new_pipeline("ipsae")
    with pipeline:
        source = predictor(Mock)
        IPSAE(structures=source, pae=source, **overrides)
        script_path = pipeline.save()
    return pipeline.tools[-1], script_path


def test_a_predictor_output_is_accepted_for_both_inputs(local_config, isolated_cwd, new_pipeline):
    tool, _ = build(new_pipeline)
    assert tool.pae_stream.name == "pae" and tool.structures_stream.name == "structures"
    assert tool.binder_chain == "A" and tool.target_chains is None
    assert (tool.pae_cutoff, tool.dist_cutoff) == (10.0, 10.0)


def test_a_predictor_without_a_pae_stream_is_refused_naming_it(local_config, isolated_cwd, new_pipeline):
    from biopipelines.ipsae import IPSAE
    from biopipelines.mock import Mock

    with new_pipeline("ipsae_nopae"):
        source = predictor(Mock, with_pae=False)
        with pytest.raises(ValueError, match="no 'pae' stream.*structures"):
            IPSAE(structures=source, pae=source)


def test_a_bare_value_for_pae_is_refused(local_config, isolated_cwd, new_pipeline):
    from biopipelines.ipsae import IPSAE
    from biopipelines.mock import Mock

    with new_pipeline("ipsae_bare"):
        with pytest.raises(ValueError, match="pae must be"):
            IPSAE(structures=predictor(Mock), pae="scores.json")


@pytest.mark.parametrize("overrides,message", [
    ({"binder_chain": ""}, "binder_chain"),
    ({"binder_chain": "A B"}, "whitespace"),
    ({"target_chains": []}, "non-empty list"),
    ({"target_chains": ["A"]}, "must not contain binder_chain"),
    ({"target_chains": ["B", "B"]}, "more than once"),
    ({"target_chains": [""]}, "target_chains"),
    ({"pae_cutoff": 0}, "pae_cutoff"),
    ({"pae_cutoff": -3.0}, "pae_cutoff"),
    ({"dist_cutoff": "10"}, "dist_cutoff"),
    ({"dist_cutoff": True}, "dist_cutoff"),
])
def test_bad_parameters_fail_at_configuration(local_config, isolated_cwd, new_pipeline, overrides, message):
    with pytest.raises(ValueError, match=message):
        build(new_pipeline, **overrides)


def test_a_pae_stream_in_an_unsupported_format_is_refused(local_config, isolated_cwd, new_pipeline):
    from biopipelines.ipsae import IPSAE
    from biopipelines.mock import Mock

    with new_pipeline("ipsae_fmt"):
        source = predictor(Mock, pae_format="csv")
        with pytest.raises(ValueError, match="npz.*json"):
            IPSAE(structures=source, pae=source)


@pytest.mark.parametrize("fmt", ["pdb", "cif", "mmcif"])
def test_every_structure_format_a_predictor_emits_is_accepted(local_config, isolated_cwd, new_pipeline, fmt):
    """Boltz2(output_format="mmcif") advertises its structures as `mmcif`, OpenFold3 as `cif`."""
    from biopipelines.ipsae import IPSAE
    from biopipelines.mock import Mock

    pipeline = new_pipeline("ipsae_structure_format")
    with pipeline:
        source = Mock(ids=["s1"], streams={"structures": {"format": fmt, "file": "<id>." + ("pdb" if fmt == "pdb" else "cif")},
                                           "pae": {"format": "npz", "file": "<id>.npz"}})
        IPSAE(structures=source, pae=source)
        pipeline.save()
    assert pipeline.tools[-1].structures_stream.format == fmt


def test_a_single_target_chain_string_is_one_chain(local_config, isolated_cwd, new_pipeline):
    tool, _ = build(new_pipeline, target_chains="B")
    assert tool.target_chains == ["B"]


def test_the_generated_script_runs_the_helper_with_every_flag_it_accepts(local_config, isolated_cwd, new_pipeline):
    tool, _ = build(new_pipeline, binder_chain="H", target_chains=["L"], pae_cutoff=12.0, dist_cutoff=8.0)
    script = tool.generate_script("unused")
    assert "pipe_ipsae.py" in script
    for flag in ("--config-yaml", "--structures-json", "--pae-json", "--ipsae-script", "--runs-dir",
                 "--ipsae-csv", "--pairs-csv", "--local-missing-csv"):
        assert flag in script
    import yaml
    assert yaml.safe_load(open(tool.config_yaml)) == {
        "binder_chain": "H", "target_chains": ["L"], "pae_cutoff": 12.0, "dist_cutoff": 8.0}
    assert os.path.exists(tool.structures_json) and os.path.exists(tool.pae_json)
    assert "pipe_propagate_missing.py" in script and tool.local_missing_csv in script


def test_every_flag_the_wrapper_emits_is_declared_by_the_pipe(local_config, isolated_cwd, new_pipeline):
    import re

    tool, _ = build(new_pipeline)
    emitted = set(re.findall(r"^\s+(--[a-z-]+)", tool.generate_script("unused").split("pipe_ipsae.py")[1].split("pipe_propagate")[0], re.M))
    declared = set(re.findall(r'add_argument\("(--[a-z-]+)"', open(SCRIPT).read()))
    assert emitted == declared


def test_the_declared_outputs_are_two_tables_and_missing_with_the_pipes_columns(local_config, isolated_cwd, new_pipeline, pipe):
    tool, _ = build(new_pipeline)
    outputs = tool.get_output_files()
    tables = outputs["tables"]
    assert set(tables) == {"ipsae", "ipsae_pairs", "missing"}
    assert tables["ipsae"].info.columns == pipe.SUMMARY_COLUMNS and tables["ipsae_pairs"].info.columns == pipe.PAIR_COLUMNS
    assert tables["missing"].info.columns == ["id", "removed_by", "kind", "cause"]
    assert "structures" not in outputs and "pae" not in outputs


# ── installation ─────────────────────────────────────────────────────────────

def test_the_install_downloads_the_upstream_script_unpinned_and_verifies_it():
    from biopipelines.ipsae import IPSAE

    script = IPSAE._install_script({"IPSAE": "/data/ipsae"}, "mamba")
    assert "raw.githubusercontent.com/DunbrackLab/IPSAE/main/ipsae.py" in script
    assert 'touch "$INSTALL_SUCCESS"' in script and "numpy" in script
    assert "if [ -s \"/data/ipsae/ipsae.py\" ]" in script


def test_force_reinstall_drops_the_already_installed_shortcut():
    from biopipelines.ipsae import IPSAE

    assert "already installed" not in IPSAE._install_script({"IPSAE": "/data/ipsae"}, "mamba", force_reinstall=True)


def test_the_install_names_the_missing_config_key():
    from biopipelines.ipsae import IPSAE

    with pytest.raises(ValueError, match="repositories.IPSAE"):
        IPSAE._install_script({}, "mamba")


# ── against the real upstream script ─────────────────────────────────────────

def _upstream_script(tmp_path_factory):
    folder = tmp_path_factory.mktemp("upstream")
    path = folder / "ipsae.py"
    try:
        with urllib.request.urlopen(UPSTREAM_URL, timeout=20) as response:
            path.write_bytes(response.read())
    except OSError as exc:
        pytest.skip(f"{UPSTREAM_URL} unreachable ({exc}); cannot run upstream ipsae.py")
    return str(path)


@pytest.fixture(scope="module")
def real_upstream(tmp_path_factory):
    return _upstream_script(tmp_path_factory)


@pytest.mark.network
@pytest.mark.parametrize("flavor", ["boltz_pdb_npz", "boltz_cif_npz", "colabfold_pdb_json", "alphafold3_cif_json"])
def test_the_real_upstream_script_is_staged_run_and_parsed(pipe, real_upstream, tmp_path, flavor):
    n_atoms = 2 * (N_A + N_B)
    structure = (write_cif(tmp_path / "x.cif") if "cif" in flavor else write_pdb(tmp_path / "x.pdb"))
    pae = {"boltz_pdb_npz": lambda: write_npz(tmp_path / "x.npz"),
           "boltz_cif_npz": lambda: write_npz(tmp_path / "x.npz"),
           "colabfold_pdb_json": lambda: write_json(tmp_path / "x.json", plddt=[90.0] * (N_A + N_B), iptm=0.8),
           "alphafold3_cif_json": lambda: write_json(tmp_path / "x.json", atom_plddts=[90.0] * n_atoms)}[flavor]()
    rows = pipe.run_ipsae(real_upstream, str(tmp_path / "run"), structure, pae, 10.0, 10.0)

    assert [(r["chain1"], r["chain2"], r["type"]) for r in rows] == [("A", "B", "asym"), ("B", "A", "asym"), ("A", "B", "max")]
    assert all(0.0 <= r["ipsae"] <= 1.0 for r in rows)
    assert rows[0]["ipsae"] > rows[1]["ipsae"]
    summary = pipe.summarize(rows, "A", None)
    assert summary["target_chain"] == "B" and summary["ipsae"] == rows[2]["ipsae"]
    assert all(math.isnan(r["iptm_af"]) for r in rows) == (flavor != "colabfold_pdb_json")
    assert all(math.isnan(r["pdockq"]) for r in rows) == flavor.endswith("npz")
    if flavor == "colabfold_pdb_json":
        assert all(r["iptm_af"] == 0.8 for r in rows)
