"""Every fixture AND shipped config must supply every folder key live code hard-subscripts.

Tools read repository and cache paths as ``self.folders["<Tool>"]``, not ``.get(...)``, so a key the fixture omits is a ``KeyError`` at config time rather than a soft fallback. That is a fixture gap, not a product bug, and it made whole tools unconstructible under the test configs: before this test, ``config.local.yaml`` was missing ``LASErMPNN``, ``RFdiffusion2`` and ``ColabFoldDatabases``, and the four scheduler fixtures were additionally missing ``DiffDock``, ``DynamicBind``, ``GEMS``, ``NeuralPLexer``, ``PLACER``, ``PocketGen``, ``RTMScore``, ``AF2BIND``, ``AlphaFoldParams`` and ``MMseqs2Server``.

Wave 2 completed the fixtures, which made the suite green while leaving the *shipped* repo-root configs untouched -- so ``ProteinMPNN(...)`` still raised ``KeyError: 'ProteinMPNN'`` on ``config.local.yaml``, and eleven more tools did the same on ``config.daint.yaml``. The shipped configs are therefore checked here too, against the same derived set: a fixture being complete says nothing about what a new user actually loads.

The required set is derived by scanning the package rather than listed here, so a new hard subscript fails this test instead of surfacing as a KeyError in whichever test happens to construct that tool next.
"""

import glob
import os
import re

import pytest


PACKAGE_GLOB = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "biopipelines", "*.py"
)
FIXTURE_GLOB = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "fixtures", "config.*.yaml")
SHIPPED_GLOB = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "config.*.yaml"
)

_SUBSCRIPT_RE = re.compile(r'folders\[\s*"([^"]+)"\s*\]')


def _required_folder_keys():
    """Folder keys the package reads as a hard subscript."""
    keys = set()
    for path in glob.glob(PACKAGE_GLOB):
        with open(path, encoding="utf-8") as handle:
            keys.update(_SUBSCRIPT_RE.findall(handle.read()))
    assert keys, "found no folders[...] subscripts -- the scan is broken, not the fixtures"
    return keys


def _runtime_supplied_keys(config_variant, isolated_cwd):
    """Keys a live Pipeline injects itself, which no config needs to declare."""
    config_variant("local")
    from biopipelines.pipeline import Pipeline

    supplied = _config_folder_keys(config_variant, "local")
    pipeline = Pipeline(project="TestSuite", job="folder_keys", description="folder keys",
                        on_the_fly=False, local_output=True, config="local")
    return set(pipeline.folders) - supplied


def _config_folder_keys(config_variant, variant):
    cm = config_variant(variant)
    keys = set()
    for section in cm.get_folder_config().values():
        if isinstance(section, dict):
            keys.update(section)
    return keys


def _fixture_variants():
    names = []
    for path in sorted(glob.glob(FIXTURE_GLOB)):
        base = os.path.basename(path)
        names.append(base[len("config."):-len(".yaml")])
    return names


FIXTURE_VARIANTS = _fixture_variants()


def test_fixture_variants_were_discovered():
    assert FIXTURE_VARIANTS, "no config.*.yaml fixtures found"
    assert "local" in FIXTURE_VARIANTS


@pytest.mark.parametrize("variant", FIXTURE_VARIANTS)
def test_fixture_supplies_every_required_folder_key(
    variant, config_variant, isolated_cwd, record_case,
):
    """A missing key is a config-time KeyError, so every fixture must be complete."""
    required = _required_folder_keys()
    runtime_keys = _runtime_supplied_keys(config_variant, isolated_cwd)
    supplied = _config_folder_keys(config_variant, variant)

    missing = sorted(required - supplied - runtime_keys)
    record_case(input=f"config.{variant}.yaml vs {len(required)} hard subscripts",
                expected=[], actual=missing)
    assert not missing, (
        f"config.{variant}.yaml omits folder keys that live code subscripts "
        f"directly, so constructing those tools raises KeyError: {missing}"
    )


# ── the shipped repo-root configs, held to the same requirement ─────────────

def _shipped_variants():
    names = []
    for path in sorted(glob.glob(SHIPPED_GLOB)):
        base = os.path.basename(path)
        names.append(base[len("config."):-len(".yaml")])
    return names


SHIPPED_VARIANTS = _shipped_variants()


def _committed_folder_keys(variant):
    """Folder keys a shipped config declares, read straight off disk.

    Not via ConfigManager on purpose: that deep-merges the gitignored ``.config.<variant>.yaml``, so a developer's own overlay could satisfy a requirement a fresh clone does not.
    """
    import yaml

    path = os.path.join(os.path.dirname(SHIPPED_GLOB), f"config.{variant}.yaml")
    with open(path, encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    keys = set()
    for section in (data.get("folders") or {}).values():
        if isinstance(section, dict):
            keys.update(section)
    return keys


def test_shipped_variants_were_discovered():
    assert set(SHIPPED_VARIANTS) >= {"cluster", "colab", "container", "daint", "local"}


@pytest.mark.parametrize("variant", SHIPPED_VARIANTS)
def test_shipped_config_supplies_every_required_folder_key(
    variant, config_variant, isolated_cwd, record_case,
):
    """A new user loads one of these, not a fixture, so completeness is theirs too."""
    required = _required_folder_keys()
    runtime_keys = _runtime_supplied_keys(config_variant, isolated_cwd)
    supplied = _committed_folder_keys(variant)

    missing = sorted(required - supplied - runtime_keys)
    record_case(input=f"config.{variant}.yaml vs {len(required)} hard subscripts",
                expected=[], actual=missing)
    assert not missing, (
        f"config.{variant}.yaml omits folder keys that live code subscripts "
        f"directly, so constructing those tools raises KeyError: {missing}"
    )


# ── the tools whose keys were missing, actually constructed ───────────────────

def test_lasermpnn_constructs_under_the_local_fixture(local_config, isolated_cwd, new_pipeline):
    """Used to raise KeyError: 'LASErMPNN' from folders["LASErMPNN"]."""
    from biopipelines.lasermpnn import LASErMPNN
    from biopipelines.mock import Mock

    pipeline = new_pipeline("lasermpnn_keys")
    with pipeline:
        m = Mock(ids=["s1"],
                 streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
                 map_table_strategy="config")
        LASErMPNN(structures=m.streams.structures)
        script_path = pipeline.save()
    assert os.path.isfile(script_path)


@pytest.mark.parametrize("tool_name", ["DiffDock", "RTMScore", "GEMS", "PLACER"])
def test_scheduler_fixture_constructs_tools_it_used_to_lack(
    tool_name, slurm_local_config, isolated_cwd, new_slurm_pipeline,
):
    """Each of these raised KeyError under config.slurm_local.yaml."""
    from biopipelines.pipeline import Resources
    from biopipelines.mock import Mock
    from biopipelines.ligand import Ligand

    pipeline = new_slurm_pipeline(f"{tool_name.lower()}_keys")
    with pipeline:
        Resources()
        m = Mock(
            ids=["s1"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"},
                     "compounds": {"format": "csv", "smiles": "CCO"}},
            map_table_strategy="config",
        )
        if tool_name == "DiffDock":
            from biopipelines.diffdock import DiffDock
            DiffDock(structures=m, compounds=Ligand(smiles="CCO", ids="c1", codes="ETH"))
        elif tool_name == "RTMScore":
            from biopipelines.rtmscore import RTMScore
            RTMScore(structures=m, ligands=m)
        elif tool_name == "GEMS":
            from biopipelines.gems import GEMS
            GEMS(structures=m, ligands=m)
        else:
            from biopipelines.placer import PLACER
            PLACER(structures=m, target_res="A50", exclude_sm=True)
        script_path = pipeline.save()
    assert os.path.isfile(script_path)
