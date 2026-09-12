"""Environment resolution and container prefixing, per env manager.

Every other fixture in this suite is ``env_manager: pip`` with no ``environments:`` block. Under that combination ``_load_environments`` sets ``self.environments = []`` and ``get_activate_command`` returns "", so the generated scripts contain no activation at all -- the whole environment layer was passing its tests by never being reached. These tests drive the four real managers (mamba, conda, micromamba, venv) and the container path against fixtures that do configure them.
"""

import pytest

NEWLINE = chr(10)


# name -> (activate command, python-resolution command) for env "MockEnv".
MANAGERS = {
    "mamba_local": (
        "mamba activate MockEnv",
        '$(mamba run -n MockEnv python -c "import sys; print(sys.executable)")',
    ),
    "conda_local": (
        "conda activate MockEnv",
        '$(conda run -n MockEnv python -c "import sys; print(sys.executable)")',
    ),
    "micromamba_local": (
        "micromamba activate MockEnv",
        '$(micromamba run -n MockEnv python -c "import sys; print(sys.executable)")',
    ),
}


@pytest.mark.parametrize("variant", sorted(MANAGERS))
def test_activation_command_per_manager(variant, config_variant, record_case):
    """Each conda-family manager activates by name with its own binary."""
    expected_activate, expected_python = MANAGERS[variant]
    cm = config_variant(variant)

    actual = (cm.get_activate_command("MockEnv"), cm.get_env_python_command("MockEnv"))
    record_case(input=f"{variant}: env 'MockEnv'",
                expected=[expected_activate, expected_python], actual=list(actual))
    assert cm.get_env_manager() == variant.replace("_local", "")
    assert actual[0] == expected_activate
    assert actual[1] == expected_python


def test_venv_activation_is_a_source_of_the_env_bin(config_variant):
    """venv has no name registry: the name only becomes a path via venv_root."""
    cm = config_variant("venv_local")
    assert cm.get_env_manager() == "venv"

    venv_path = cm.get_venv_path("MockEnv")
    assert venv_path.endswith("/_stub_venvs/MockEnv"), venv_path
    assert cm.get_activate_command("MockEnv") == f'source "{venv_path}/bin/activate"'
    assert cm.get_env_python_command("MockEnv") == f'"{venv_path}/bin/python"'


def test_venv_absolute_env_name_passes_through(config_variant):
    cm = config_variant("venv_local")
    assert cm.get_venv_path("/opt/envs/prebuilt") == "/opt/envs/prebuilt"


def test_pip_manager_activates_nothing(config_variant):
    """Why the environment layer was untested: pip resolves to no activation."""
    cm = config_variant("local")
    assert cm.get_env_manager() == "pip"
    assert cm.get_activate_command("MockEnv") == ""
    assert cm.get_environment("Mock") is None


@pytest.mark.parametrize("variant", sorted(MANAGERS) + ["venv_local"])
def test_environments_block_names_the_tool_env(variant, config_variant):
    """A configured name must win over the "biopipelines" default."""
    cm = config_variant(variant)
    assert cm.get_environment("Mock") == "MockEnv"
    assert cm.get_environment("Sequence") == "SeqEnv"
    # A tool absent from the block gets no name, and so falls back downstream.
    assert cm.get_environment("NotATool") is None


@pytest.mark.parametrize("variant", sorted(MANAGERS) + ["venv_local"])
def test_tool_resolves_its_configured_environment(
    variant, config_variant, isolated_cwd,
):
    """Under a real manager a tool must carry the configured env name."""
    config_variant(variant)
    from biopipelines.pipeline import Pipeline
    from biopipelines.mock import Mock

    pipeline = Pipeline(project="TestSuite", job="env_res", description="env resolution",
                        on_the_fly=False, local_output=True, config=variant)
    with pipeline:
        Mock(ids=["a"],
             streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
             map_table_strategy="config")

    # The tool object behind the StandardizedOutput carries the resolved envs.
    assert pipeline.tools[0].environments == ["MockEnv"], pipeline.tools[0].environments


@pytest.mark.parametrize("variant", sorted(MANAGERS) + ["venv_local"])
def test_emitted_step_script_activates_the_configured_env(
    variant, config_variant, isolated_cwd,
):
    """End-to-end: the activation actually lands in the generated bash."""
    import pathlib

    cm = config_variant(variant)
    from biopipelines.pipeline import Pipeline
    from biopipelines.mock import Mock

    pipeline = Pipeline(project="TestSuite", job="env_emit", description="env emission",
                        on_the_fly=False, local_output=True, config=variant)
    with pipeline:
        Mock(ids=["a"],
             streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
             map_table_strategy="config")
        script_path = pipeline.save()

    runtime_dir = pathlib.Path(script_path).parent
    step_scripts = sorted(runtime_dir.glob("*_Mock.sh"))
    assert step_scripts, f"no Mock step script under {runtime_dir}"
    body = step_scripts[0].read_text(encoding="utf-8")

    assert 'echo "Requested: MockEnv"' in body, "the step never announces its env"
    assert cm.get_activate_command("MockEnv") in body, (
        f"{variant}: step script is missing the activation line "
        f"{cm.get_activate_command('MockEnv')!r}"
    )


def test_pip_step_script_has_no_activation(local_config, isolated_cwd, new_pipeline):
    """The baseline this suite used everywhere: nothing to activate, nothing emitted."""
    import pathlib

    pipeline = new_pipeline("env_pip")
    with pipeline:
        from biopipelines.mock import Mock
        Mock(ids=["a"],
             streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
             map_table_strategy="config")
        script_path = pipeline.save()

    runtime_dir = pathlib.Path(script_path).parent
    body = sorted(runtime_dir.glob("*_Mock.sh"))[0].read_text(encoding="utf-8")
    for manager in ("mamba activate", "conda activate", "micromamba activate",
                    "bin/activate"):
        assert manager not in body, f"pip mode must not emit {manager!r}"


# ── container prefixing ───────────────────────────────────────────────────────

def test_container_prefix_wraps_the_command(config_variant, isolated_cwd):
    """A configured image must produce an `<executor> exec` prefix with binds."""
    cm = config_variant("container_local")
    assert cm.get_container_executor() == "apptainer"

    from biopipelines.pipeline import Pipeline
    from biopipelines.mock import Mock

    pipeline = Pipeline(project="TestSuite", job="container", description="container",
                        on_the_fly=False, local_output=True, config="container_local")
    with pipeline:
        Mock(ids=["a"],
             streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
             map_table_strategy="config")
        pipeline.save()

    tool = pipeline.tools[0]
    assert tool.uses_container(), \
        "a containers: entry for Mock must reach folders as container:Mock"

    prefix = tool.container_prefix()
    assert prefix.startswith("apptainer exec --nv "), prefix
    assert prefix.endswith("mock.sif "), prefix
    assert "-B " in prefix, f"bind mounts missing, outputs would be invisible: {prefix}"


def test_no_container_means_no_prefix(local_config, isolated_cwd, new_pipeline):
    """container_executor: none must leave commands unwrapped."""
    pipeline = new_pipeline("no_container")
    with pipeline:
        from biopipelines.mock import Mock
        Mock(ids=["a"],
             streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
             map_table_strategy="config")
        pipeline.save()

    tool = pipeline.tools[0]
    assert not tool.uses_container()
    assert tool.container_prefix() == ""


# ── the environments: block redirects the install too ─────────────────────────

def test_install_script_honours_the_configured_environment(config_variant, record_case, tmp_path):
    """An `environments:` entry must redirect the install, not just the activation.

    The redirected env gets a real spec file here, because that is what makes the redirect work end to end: a spec is named after the environment it defines, so pointing a tool at `MyDSSPEnv` without writing `environments/MyDSSPEnv.yaml` leaves nothing to create. Asserting the `env create` line alone proved only that the command mentioned the new name.
    """
    cm = config_variant("mamba_local")
    assert cm.get_environment("DSSP") == "MyDSSPEnv", "fixture must redirect DSSP"

    repo = tmp_path / "bp"
    (repo / "environments").mkdir(parents=True)
    spec = NEWLINE.join(["name: MyDSSPEnv", "channels: [conda-forge]", "dependencies: [python=3.11]", ""])
    (repo / "environments" / "MyDSSPEnv.yaml").write_text(spec, encoding="utf-8")

    from biopipelines.dssp import DSSP
    script = DSSP._install_script({"biopipelines": str(repo)}, env_manager="mamba")

    record_case(input="environments: {DSSP: MyDSSPEnv}, env_manager: mamba",
                expected="install addresses MyDSSPEnv",
                actual=[l for l in script.splitlines() if "mamba " in l])
    assert "MyDSSPEnv" in script, (
        "install script builds a different env than the one the pipeline activates"
    )
    # Not just mentioned: every env-addressing command must name it, and none may fall back to the wrapper's own literal.
    assert "mamba env create" in script and "-n MyDSSPEnv" in script, script
    assert "mamba run -n MyDSSPEnv" in script, script
    assert f"run -n {DSSP.ENV_NAME} " not in script, (
        f"install still addresses the hardcoded {DSSP.ENV_NAME!r} env"
    )


def test_install_and_runtime_env_names_agree(config_variant):
    """The env the install builds is the env the pipeline activates.

    Replaces a test that pinned the old split (install built "dssp" while the runtime activated the configured name). Both sides are now resolved from the same `environments:` entry, so the check is that they match rather than that they differ.
    """
    cm = config_variant("mamba_local")
    from biopipelines.dssp import DSSP

    configured = cm.get_environment("DSSP")
    assert configured == "MyDSSPEnv"
    assert cm.get_activate_command(configured) == "mamba activate MyDSSPEnv"

    install_script = DSSP._install_script({"biopipelines": "/tmp/bp"},
                                          env_manager="mamba")
    assert DSSP._install_env("mamba") == configured
    assert f"mamba run -n {configured}" in install_script
    assert "mamba run -n dssp " not in install_script


def test_install_env_is_resolved_the_same_way_for_every_env_building_tool(
    config_variant,
):
    """The resolution lives in one place, so a tool cannot drift from it.

    Guards the shared helper rather than 42 wrappers: every tool that declares an ENV_NAME must return the configured entry from ``_install_env``, and the two tools the fixture redirects prove the lookup is not just echoing ENV_NAME back.
    """
    cm = config_variant("mamba_local")
    from biopipelines import dssp, mock, panda

    assert dssp.DSSP._install_env("mamba") == "MyDSSPEnv" != dssp.DSSP.ENV_NAME
    assert panda.Panda._install_env("mamba") == cm.get_environment("Panda")
    # BaseConfig's own default is empty: a tool must declare its env, not inherit one.
    from biopipelines.base_config import BaseConfig
    assert BaseConfig.ENV_NAME == ""
    assert mock.Mock.ENV_NAME == "", "a base-env tool declares no env of its own"


def test_missing_entry_raises_under_a_real_manager(config_variant):
    """No hardcoded fallback: an unconfigured tool fails at install time, by name."""
    cm = config_variant("mamba_local")
    assert cm.get_environment("XTB") is None, "fixture must not configure XTB"

    from biopipelines.xtb import XTB
    with pytest.raises(ValueError) as excinfo:
        XTB._install_script({"biopipelines": "/tmp/bp"}, env_manager="mamba")

    message = str(excinfo.value)
    # The message has to be actionable on its own: which tool, which config, what to add.
    assert "XTB" in message
    assert "config.mamba_local.yaml" in message
    assert "mamba_local" in message
    assert "environments:" in message
    assert XTB.ENV_NAME == "xtb", "the literal still belongs on the class, for pip"


def test_missing_entry_does_not_raise_under_pip(config_variant):
    """pip mode is exempt and keeps emitting the tool's own literal."""
    cm = config_variant("local")
    assert cm.get_env_manager() == "pip"
    assert cm.get_environment("XTB") is None

    from biopipelines.xtb import XTB
    script = XTB._install_script({"biopipelines": "/tmp/bp"}, env_manager="pip")
    assert XTB._install_env("pip") == XTB.ENV_NAME == "xtb"
    assert "xtb" in script


def test_explicit_env_name_argument_wins(config_variant):
    """A caller-supplied name outranks both the config and the raise.

    ``venv_local`` configures no OpenMM entry, so resolution would raise; the explicit argument is what keeps a caller that pins a name (an installer, a warm-up) working.
    """
    cm = config_variant("venv_local")
    from biopipelines.openmm import OpenMM

    assert cm.get_environment("OpenMM") is None
    assert OpenMM._install_env("venv", "CallerEnv") == "CallerEnv"
    with pytest.raises(ValueError):
        OpenMM._install_env("venv")

    script = OpenMM._install_script({"biopipelines": "/tmp/bp"},
                                    env_manager="venv", env_name="CallerEnv")
    assert cm.get_venv_path("CallerEnv") in script
    assert cm.get_venv_path(OpenMM.ENV_NAME) not in script


# ── the naming standard ───────────────────────────────────────────────────────

def test_pip_mode_reads_the_config_before_falling_back():
    """`env_manager: pip` returned ENV_NAME before consulting the config at all, so an environments: entry was ignored there -- the one manager where the wrapper's own literal still won. It is exempt from the raise, not from the config."""
    import inspect

    from biopipelines.base_config import BaseConfig

    source = inspect.getsource(BaseConfig._install_env)
    config_read = source.index("get_environment")
    pip_fallback = source.index('env_manager == "pip"', source.index("if env_name"))
    assert config_read < pip_fallback, "pip mode must consult the config before its fallback"


def test_pip_mode_refuses_to_build_an_env_instead_of_emitting_fake_commands():
    """It used to emit `pip env create -f ...` and `pip run -n ...`; neither is a pip command, so the script could only fail later and less clearly."""
    from biopipelines.base_config import BaseConfig

    block = BaseConfig._env_install_block("SomeEnv", "pip", "/repo")
    assert "pip env create" not in block
    assert "pip run -n" not in block
    assert "env_manager" in block and "exit 1" in block


def test_an_env_with_no_spec_falls_through_to_assuming_it_exists():
    """A spec file is named after the environment it defines, but an env may legitimately be built by an upstream installer instead: RF3's `modelforge` has no spec on cluster, colab or container. The check has to be in the emitted bash, not in Python -- the script is generated locally and runs on the cluster, so `folders["biopipelines"]` is a remote path no config-time `os.path.isfile` can see."""
    from biopipelines.base_config import BaseConfig

    block = BaseConfig._env_install_block("modelforge", "mamba", "/remote/repo")
    assert 'elif [ -f "/remote/repo/environments/modelforge.yaml" ]' in block
    assert "already exists" in block
    # Both spec spellings still get their create; only the third branch is new.
    assert block.count("env create -f") >= 2
