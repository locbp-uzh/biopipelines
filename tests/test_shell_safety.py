"""Shell-safety validation tests.

Covers the three layers added in response to the reviewer's observation
that user-controlled strings reach generated bash:

  1. Pipeline(project=..., job=...)   -> _validate_identifier
  2. Free-form tool strings           -> _validate_freeform_string
  3. config.<variant>.yaml values     -> ConfigManager._validate_shell_safety

The bar is: invalid input raises ValueError at Python-construction / load
time, naming the offending field, before any shell script is written.
"""

import os
import pytest
import yaml


# ── 1. Pipeline identifier validation ────────────────────────────────────────

@pytest.mark.parametrize("bad_value", [
    "has space",
    'has"quote',
    "has$dollar",
    "has`backtick",
    "has;semicolon",
    "..",
    ".",
    "-rf",
    "",
])
def test_pipeline_rejects_unsafe_project(local_config, isolated_cwd, bad_value):
    from biopipelines.pipeline import Pipeline
    with pytest.raises(ValueError, match="project"):
        Pipeline(
            project=bad_value, job="ok", description="d",
            on_the_fly=False, local_output=True, config="local",
        )


@pytest.mark.parametrize("bad_value", [
    'has"quote', "has$dollar", "has`backtick", "..", ".", "-rf", "",
])
def test_pipeline_rejects_unsafe_job(local_config, isolated_cwd, bad_value):
    from biopipelines.pipeline import Pipeline
    with pytest.raises(ValueError, match="job"):
        Pipeline(
            project="ok", job=bad_value, description="d",
            on_the_fly=False, local_output=True, config="local",
        )


def test_pipeline_accepts_safe_identifiers(local_config, isolated_cwd):
    from biopipelines.pipeline import Pipeline
    Pipeline(
        project="Proj.1-A_2", job="run_03",
        description="fine", on_the_fly=False, local_output=True, config="local",
    )


# ── 2. Description escaping at emission time ─────────────────────────────────

def test_pipeline_description_is_escaped_in_generated_script(
    local_config, isolated_cwd,
):
    """A description containing bash-special characters must not terminate
    the surrounding double-quoted echo literal."""
    from biopipelines.pipeline import Pipeline
    from biopipelines.sequence import Sequence

    tricky = 'breaks " and $VAR and `cmd` and \\slash'
    pipeline = Pipeline(
        project="TestSuite", job="escape",
        description=tricky,
        on_the_fly=False, local_output=True, config="local",
    )
    with pipeline:
        Sequence(seq="MKTAYIA", type="protein", ids="demo")
        pipeline.save()

    # The echo "Description: ..." line is emitted into config.sh, which is
    # piped into pipeline.sh at runtime.
    runtime_dir = os.path.join(pipeline.folders["runtime"])
    config_sh = os.path.join(runtime_dir, "config.sh")
    content = open(config_sh, encoding="utf-8").read()

    # Raw unescaped chars would break the surrounding "..." literal; assert
    # the escaped forms are what actually reached the file.
    assert 'breaks \\"' in content, f"\" not escaped in:\n{content}"
    assert "\\$VAR" in content, f"$ not escaped in:\n{content}"
    assert "\\`cmd\\`" in content, f"` not escaped in:\n{content}"


# ── 3. _validate_freeform_string helper ──────────────────────────────────────

@pytest.mark.parametrize("bad_value", [
    'has"quote', "has$dollar", "has`backtick", "has\\backslash",
])
def test_freeform_string_rejects_unsafe_chars(bad_value):
    from biopipelines.base_config import _validate_freeform_string
    with pytest.raises(ValueError, match="foo"):
        _validate_freeform_string("foo", bad_value)


@pytest.mark.parametrize("ok_value", [
    "plain", "with space", "with/slash", "a-b_c.d", "", "<placeholder>",
])
def test_freeform_string_accepts_safe_values(ok_value):
    from biopipelines.base_config import _validate_freeform_string
    _validate_freeform_string("foo", ok_value)  # must not raise


def test_freeform_string_accepts_none():
    from biopipelines.base_config import _validate_freeform_string
    _validate_freeform_string("foo", None)  # must not raise


def test_freeform_string_rejects_non_string():
    from biopipelines.base_config import _validate_freeform_string
    with pytest.raises(ValueError, match="must be a string"):
        _validate_freeform_string("foo", 42)


# ── 3b. Per-tool wiring (spot checks) ────────────────────────────────────────

def test_contacts_rejects_unsafe_ligand_name(local_config, isolated_cwd):
    """Contacts.ligand must be a compounds stream, not a string. A bare string
    (which could otherwise carry shell-injection chars) is rejected outright at
    construction; the residue code is resolved from the stream's `code` column."""
    from biopipelines.pipeline import Pipeline
    from biopipelines.mock import Mock
    from biopipelines.contacts import Contacts

    pipeline = Pipeline(
        project="TestSuite", job="contacts",
        on_the_fly=True, local_output=True, config="local",
    )
    structures = Mock(
        ids=["x"],
        streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
    )
    with pytest.raises(ValueError, match="ligand"):
        Contacts(
            structures=structures.streams.structures,
            ligand='LIG"; rm -rf /',
        )


def test_af2bind_rejects_unsafe_chain(local_config, isolated_cwd):
    """AF2BIND.chain is interpolated into an echo line and a CLI arg; unsafe
    chars must be rejected by the tool's validate_params()."""
    from biopipelines.pipeline import Pipeline
    from biopipelines.mock import Mock
    from biopipelines.af2bind import AF2BIND

    Pipeline(
        project="TestSuite", job="af2bind",
        on_the_fly=True, local_output=True, config="local",
    )
    structures = Mock(
        ids=["x"],
        streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
    )
    with pytest.raises(ValueError, match="chain"):
        AF2BIND(
            structures=structures.streams.structures,
            chain='A`rm -rf /`',
        )


def test_rfdiffusion_allatom_rejects_unsafe_contigs(local_config, isolated_cwd):
    """RFdiffusionAllAtom.contigs is interpolated into an unquoted CLI arg.
    Shell metas must be rejected at config time."""
    from biopipelines.pipeline import Pipeline
    from biopipelines.rfdiffusion_allatom import RFdiffusionAllAtom

    Pipeline(
        project="TestSuite", job="rfdaa",
        on_the_fly=True, local_output=True, config="local",
    )
    from biopipelines.ligand import Ligand
    with pytest.raises(ValueError, match="contigs"):
        RFdiffusionAllAtom(ligand=Ligand(code="LIG"), contigs="A1-50 `rm -rf /`")


def test_rfdiffusion2_rejects_unsafe_contigs(local_config, isolated_cwd):
    """RFdiffusion2.contigs is interpolated into an unquoted Hydra CLI arg.
    Shell metas must be rejected at config time."""
    from biopipelines.pipeline import Pipeline
    from biopipelines.rfdiffusion2 import RFdiffusion2
    from biopipelines.datastream import DataStream
    from biopipelines.ligand import Ligand

    Pipeline(
        project="TestSuite", job="rfd2",
        on_the_fly=True, local_output=True, config="local",
    )
    pdb = DataStream(name="structures", ids=["x"], files=["<id>.pdb"],
                     map_table="", format="pdb")
    with pytest.raises(ValueError, match="contigs"):
        RFdiffusion2(pdb=pdb, ligand=Ligand(code="LIG"),
                     contigs="A1-50 `rm -rf /`")


def test_hbdesigner_rejects_unsafe_guide_res(local_config, isolated_cwd):
    """HBDesigner.guide_res literal reaches the generated bash through the
    constraints resolver invocation. Shell metas must be rejected at config
    time."""
    from biopipelines.pipeline import Pipeline
    from biopipelines.hbdesigner import HBDesigner
    from biopipelines.datastream import DataStream

    Pipeline(
        project="TestSuite", job="hbdes",
        on_the_fly=True, local_output=True, config="local",
    )
    pdb = DataStream(name="structures", ids=["x"], files=["<id>.pdb"],
                     map_table="", format="pdb")
    with pytest.raises(ValueError, match="guide_res"):
        HBDesigner(structures=pdb, guide_res="A12 `rm -rf /`")


def test_ligand_rejects_unsafe_code(local_config, isolated_cwd):
    """Ligand(code=...) reaches the generated bash; the residue-code contract
    (1-5 alphanumeric, extended CCD) is enforced at construction, rejecting shell
    metacharacters. This is where the former per-tool ligand_code injection
    check now lives — Ligand is the sole validator of the code."""
    from biopipelines.pipeline import Pipeline
    from biopipelines.ligand import Ligand

    Pipeline(
        project="TestSuite", job="ligand",
        on_the_fly=True, local_output=True, config="local",
    )
    with pytest.raises(ValueError, match="alphanumeric"):
        Ligand(code='LIG"; rm -rf /')


def test_table_rejects_unsafe_name(local_config, isolated_cwd, tmp_path):
    """Table.name is echoed directly into bash; unsafe chars must be
    rejected by validate_params()."""
    import pandas as pd
    from biopipelines.pipeline import Pipeline
    from biopipelines.table import Table

    csv = tmp_path / "data.csv"
    pd.DataFrame({"id": ["a"], "v": [1]}).to_csv(csv, index=False)

    Pipeline(
        project="TestSuite", job="tbl",
        on_the_fly=True, local_output=True, config="local",
    )
    with pytest.raises(ValueError, match="name"):
        Table(path=str(csv), name='bad"name')


# ── 4. ConfigManager shell-safety at load time ───────────────────────────────

def _write_config(tmp_path, mutator):
    """Copy the fixture config, apply `mutator(dict) -> None`, write back,
    and point ConfigManager at the copy. Returns the config path."""
    from biopipelines.config_manager import ConfigManager

    src = os.path.join(
        os.path.dirname(__file__), "fixtures", "config.local.yaml"
    )
    with open(src) as f:
        data = yaml.safe_load(f)
    mutator(data)

    dst = tmp_path / "config.local.yaml"
    with open(dst, "w") as f:
        yaml.safe_dump(data, f)

    ConfigManager._instance = None
    ConfigManager._config = None
    ConfigManager._variant = None
    return str(dst)


def _load_with_config_path(monkeypatch, config_path):
    from biopipelines.config_manager import ConfigManager
    monkeypatch.setattr(
        ConfigManager, "_get_config_path",
        classmethod(lambda cls, variant=None: config_path),
    )
    ConfigManager._instance = None
    ConfigManager._config = None
    ConfigManager._variant = None
    return ConfigManager(variant="local")


@pytest.mark.parametrize("field,value,expected_in_msg", [
    ("folders.base.home",        '/home/"oops',       "folders.base.home"),
    ("folders.base.home",        "/home/$(whoami)",   "folders.base.home"),
    ("folders.base.home",        "/home/`whoami`",    "folders.base.home"),
    ("folders.infrastructure.pipe_scripts",
                                 "/path/with\\back",  "folders.infrastructure.pipe_scripts"),
])
def test_config_rejects_unsafe_folder_value(
    tmp_path, monkeypatch, field, value, expected_in_msg,
):
    """Shell-unsafe characters anywhere in folders.* must fail at load time."""
    section, subsection, key = field.split(".")

    def mutate(cfg):
        cfg[section][subsection][key] = value

    path = _write_config(tmp_path, mutate)

    with pytest.raises(ValueError, match=expected_in_msg):
        _load_with_config_path(monkeypatch, path)


def test_config_rejects_unsafe_slurm_module(tmp_path, monkeypatch):
    def mutate(cfg):
        cfg["machine"]["scheduler"]["modules"] = ["ok_mod", "bad$mod"]
    path = _write_config(tmp_path, mutate)
    with pytest.raises(ValueError, match=r"scheduler\.modules\[1\]"):
        _load_with_config_path(monkeypatch, path)


def test_config_rejects_unsafe_env_manager(tmp_path, monkeypatch):
    """env_manager.name is interpolated *unquoted* into `<mgr> activate <env>`,
    so any value outside the allowlist — including command-injection attempts
    like 'mamba; rm -rf /' — must be rejected at load time."""
    def mutate(cfg):
        cfg["machine"]["env_manager"]["name"] = "mamba; rm -rf /"
    path = _write_config(tmp_path, mutate)
    with pytest.raises(ValueError, match="machine.env_manager.name"):
        _load_with_config_path(monkeypatch, path)


def test_config_rejects_unknown_scheduler(tmp_path, monkeypatch):
    def mutate(cfg):
        cfg["machine"]["scheduler"]["name"] = "cobalt"
    path = _write_config(tmp_path, mutate)
    with pytest.raises(ValueError, match="machine.scheduler.name"):
        _load_with_config_path(monkeypatch, path)


@pytest.mark.parametrize("name", ["slurm", "lsf", "pbs", "colab", "none"])
def test_config_accepts_known_schedulers(tmp_path, monkeypatch, name):
    def mutate(cfg):
        cfg["machine"]["scheduler"]["name"] = name
    path = _write_config(tmp_path, mutate)
    cm = _load_with_config_path(monkeypatch, path)
    assert cm.get_scheduler() == name


def test_config_accepts_safe_paths_with_spaces_and_placeholders(
    tmp_path, monkeypatch,
):
    """Spaces, forward slashes, and <placeholder> syntax are all legitimate
    in folder paths and must NOT be rejected."""
    def mutate(cfg):
        cfg["folders"]["base"]["home"] = "/home/user with space"
        cfg["folders"]["infrastructure"]["cache"] = "<biopipelines>/cache"
    path = _write_config(tmp_path, mutate)
    # Must not raise:
    _load_with_config_path(monkeypatch, path)


# ── 4. Generated scripts quote paths so a spaced workspace doesn't break bash ──

def test_generated_scripts_quote_paths_under_spaced_workspace(
    local_config, tmp_path, monkeypatch,
):
    """A workspace whose path contains spaces (e.g. OneDrive paths) must not
    produce broken bash: every `python <script>` invocation the framework emits
    has to carry a quoted script path. Regression for the completion-check /
    missing-propagation snippets that interpolated paths unquoted."""
    import re

    spaced = tmp_path / "dir with space"
    spaced.mkdir()
    monkeypatch.chdir(spaced)

    from biopipelines.pipeline import Pipeline
    from biopipelines.sequence import Sequence

    pipeline = Pipeline(
        project="TestSuite", job="spaced", description="spaced-path regression",
        on_the_fly=False, local_output=True, config="local",
    )
    with pipeline:
        Sequence(seq="MKTAYIAKQRQISFVKSHFSRQLEERLGL", type="protein", ids="demo")
        script_path = pipeline.save()

    runtime_dir = os.path.join(os.path.dirname(script_path))
    scripts = [
        os.path.join(runtime_dir, f)
        for f in os.listdir(runtime_dir)
        if f.endswith(".sh")
    ]
    assert scripts, "no generated scripts found"

    # Any script path passed to `python` must be quoted; an unquoted spaced
    # path word-splits in bash. A correctly-quoted call reads `python "...`;
    # an unquoted one reads `python /abs/path...`. Flag the latter.
    bad = re.compile(r'python\s+[^"\'$\s-]')
    for sh in scripts:
        content = open(sh, encoding="utf-8").read()
        for line in content.splitlines():
            stripped = line.strip()
            assert not bad.search(stripped), (
                f"unquoted script path in {os.path.basename(sh)}: {stripped!r}"
            )


# ── 5. HBDesigner emits its command without eval ─────────────────────────────

def _hbdesigner_run_block(ids):
    """Emit HBDesigner's per-input run loop for a structures stream of `ids`."""
    import pathlib
    from biopipelines.pipeline import Pipeline
    from biopipelines.hbdesigner import HBDesigner
    from biopipelines.datastream import DataStream

    pipeline = Pipeline(project="TestSuite", job="hbdes", description="d",
                        on_the_fly=False, local_output=True, config="local")
    with pipeline:
        pdb = DataStream(name="structures", ids=ids, files=["<id>.pdb"],
                         map_table="", format="pdb")
        HBDesigner(structures=pdb)
        script_path = pipeline.save()
    step = sorted(pathlib.Path(script_path).parent.glob("*_HBDesigner.sh"))[0]
    body = step.read_text(encoding="utf-8")
    start = body.index("for STRUCT_ID in")
    return body[start:body.index("\ndone", start) + len("\ndone")]


def test_hbdesigner_does_not_eval_its_command(local_config, isolated_cwd):
    """`eval` re-parses the assembled line, so any runtime value reaching it is
    a second parse away from executing. The command must go out as a bash array."""
    import re

    block = _hbdesigner_run_block(["x"])
    # Match `eval` as a command word only; a tmp_path name can contain the substring.
    evals = [l for l in block.splitlines() if re.match(r"\s*(if\s+!\s+)?eval\b", l)]
    assert not evals, f"HBDesigner still emits through eval: {evals}"
    assert 'HBDES_OPTIONS=(' in block, f"expected a bash array:\n{block}"
    assert 'run_hbdesigner "${HBDES_OPTIONS[@]}"' in block, (
        f"array must be expanded one-element-per-argument:\n{block}"
    )


@pytest.mark.skipif(not __import__("shutil").which("bash"),
                    reason="needs bash to replay the emitted loop")
def test_hbdesigner_hostile_id_executes_nothing(local_config, isolated_cwd, tmp_path):
    """A stream id carrying `";cmd;"` reaches $OUT_DIR. Under the former
    `eval run_hbdesigner $HBDES_OPTIONS` it broke out of the quoting and ran
    `cmd`; replayed against the emitted loop it must stay one argument."""
    import subprocess

    hostile = 'x";>PWNED_ID;"'  # no whitespace, so the loop's word splitting keeps it whole
    block = _hbdesigner_run_block([hostile])

    lines = ['run_hbdesigner() { printf "ARGV:"; for a in "$@"; do printf " <%s>" "$a"; done; echo; }']
    for line in block.splitlines():
        s = line.strip()
        if s.startswith("for STRUCT_ID in"):
            lines.append("for STRUCT_ID in $(printf '%s\n' \"$HOSTILE\"); do")
        elif s.startswith("INPUT_PDB="):
            lines.append('    INPUT_PDB=/tmp/in.pdb')
        elif s.startswith("OUT_DIR="):
            lines.append('    OUT_DIR="/tmp/exec/$STRUCT_ID"')
        elif s.startswith("OPTS="):
            lines.append('    OPTS=$(printf "\n\n\n")')
        elif s.startswith("mkdir -p") or s.startswith("rm -f") or s.startswith("touch "):
            continue
        else:
            lines.append(line)
    harness = tmp_path / "replay.sh"
    harness.write_text("\n".join(lines) + "\n", encoding="utf-8")

    result = subprocess.run(
        ["bash", str(harness)], capture_output=True, text=True,
        cwd=str(tmp_path), env=dict(os.environ, HOSTILE=hostile),
    )
    assert not os.path.exists(tmp_path / "PWNED_ID"), (
        f"injected command ran; stdout={result.stdout!r} stderr={result.stderr!r}"
    )
    assert f'<--out_dir> </tmp/exec/{hostile}>' in result.stdout, (
        f"the id must survive as one argument; got {result.stdout!r}"
    )


# ── 6. LigandMPNN: the last eval in emitted bash ──────────────────────────────

def test_ligandmpnn_does_not_eval_its_command(local_config, isolated_cwd, new_pipeline):
    """LigandMPNN assembled `eval python run.py … --pdb_path '"$PDB_FILE"' $FIXED_OPTION` -- the same escaped-quote-under-eval construction HBDesigner shed. `$PDB_FILE` is a runtime-resolved path and the position options carried a flag plus a quoted value, which is why the eval was there: it split them. A bash array does that without a second parse."""
    import glob

    from biopipelines.ligand import Ligand
    from biopipelines.ligand_mpnn import LigandMPNN
    from biopipelines.mock import Mock

    pipeline = new_pipeline("lmpnn_noeval")
    with pipeline:
        seeds = Mock(ids=["s1"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
        ligand = Ligand(codes="STI")
        LigandMPNN(structures=seeds.streams.structures, ligand=ligand.streams.compounds)

    script = glob.glob(pipeline.folders["runtime"] + "/*LigandMPNN*.sh")[0]
    body = open(script, encoding="utf-8").read()
    assert "eval " not in body, "the command is evaluated a second time"
    assert 'LMPNN_ARGS=(' in body
    assert 'python run.py "${LMPNN_ARGS[@]}"' in body


def test_a_hostile_path_stays_one_argument_for_ligandmpnn():
    """Replays the emitted shape: a path that closes its own quoting used to reach a second parse."""
    import subprocess
    import tempfile

    workdir = tempfile.mkdtemp()
    script = (
        f'cd {workdir}\n'
        'PDB_FILE=\'x";>PWNED;"\'\n'
        "POSITION_OPTIONS=(--fixed_residues 'A10 A11')\n"
        'LMPNN_ARGS=("--model_type" "ligand_mpnn" --pdb_path "$PDB_FILE" "${POSITION_OPTIONS[@]}")\n'
        'printf "ARG <%s>\n" "${LMPNN_ARGS[@]}"\n'
    )
    result = subprocess.run(["bash", "-c", script], capture_output=True, text=True)
    import os

    assert "PWNED" not in os.listdir(workdir), "the injected command ran"
    assert 'ARG <x";>PWNED;">' in result.stdout, "the path was split instead of staying one argument"
    assert "ARG <A10 A11>" in result.stdout, "a residue list with a space was word-split"
