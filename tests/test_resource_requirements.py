# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""A tool declares which GPUs it does not support, and generating a batch warns when its request can land on one.

DynamicBind and PocketGen run on an A100 and stop on an H100 (their torch builds have no sm_90 kernels), but
`Resources(gpu="any")` can land on an H100 and nothing said so until the CUDA error. The declaration is the class
attribute `RESOURCE_REQUIREMENTS`; the judgement of a request against it is pure and is tested on its own, then through
the batch-script generation.
"""

from pathlib import Path

import pytest

from biopipelines.base_config import BaseConfig
from biopipelines.config_manager import ConfigManager
from biopipelines.resource_requirements import gpu_warnings, resolve_request, validate_requirements

S3IT = {"L4": "", "A100": "", "H100": "", "H200": "140GB"}
DYNAMICBIND = {"gpu": {"incompatible": ["H100", "H200"], "use": "A100"}}
FIXTURES_DIR = Path(__file__).parent / "fixtures"


class TestTheDefault:

    def test_a_tool_declares_nothing_unless_it_says_so(self):
        assert BaseConfig.RESOURCE_REQUIREMENTS == {}

    def test_nothing_declared_never_warns(self):
        assert gpu_warnings("Tool", "H100", {}, S3IT) == []
        assert gpu_warnings("Tool", "any", {"gpu": {}}, S3IT) == []

    def test_the_two_tools_with_a_verified_limit_declare_it(self):
        from biopipelines.dynamicbind import DynamicBind
        from biopipelines.pocketgen import PocketGen

        for tool in (DynamicBind, PocketGen):
            assert tool.RESOURCE_REQUIREMENTS == {
                "gpu": {"incompatible": ["H100", "H200", "GH200"], "use": "A100"}}, tool.TOOL_NAME


class TestAnExplicitModel:

    def test_a_request_naming_an_incompatible_model_warns_without_any_inventory(self):
        messages = gpu_warnings("DynamicBind", "H100", DYNAMICBIND, {})

        assert len(messages) == 1
        assert "DynamicBind" in messages[0] and "H100" in messages[0]
        assert 'Resources(gpu="H100")' in messages[0]
        assert 'Resources(gpu="A100")' in messages[0], "the fix is named"

    def test_the_model_name_is_compared_without_regard_to_case(self):
        assert len(gpu_warnings("DynamicBind", "h100", DYNAMICBIND, {})) == 1

    def test_a_supported_model_is_silent(self):
        assert gpu_warnings("DynamicBind", "A100", DYNAMICBIND, S3IT) == []

    def test_none_is_silent(self):
        assert gpu_warnings("DynamicBind", "none", DYNAMICBIND, S3IT) == []
        assert gpu_warnings("DynamicBind", None, DYNAMICBIND, S3IT) == []
        assert gpu_warnings("DynamicBind", "", DYNAMICBIND, S3IT) == []

    def test_a_token_that_is_neither_a_model_nor_a_memory_is_not_judged(self):
        assert gpu_warnings("DynamicBind", "banana", DYNAMICBIND, S3IT) == []


class TestARequestThatPinsNoModel:

    def test_any_may_land_on_every_incompatible_model_of_the_inventory(self):
        messages = gpu_warnings("DynamicBind", "any", DYNAMICBIND, S3IT)

        assert len(messages) == 1
        assert "may land on H100, H200" in messages[0]
        assert 'Resources(gpu="any")' in messages[0]
        assert 'Resources(gpu="A100")' in messages[0]

    @pytest.mark.parametrize("request_", ["gpu", "any", 1, "2"])
    def test_every_spelling_of_any_gpu_is_judged_alike(self, request_):
        assert "may land on H100, H200" in gpu_warnings("DynamicBind", request_, DYNAMICBIND, S3IT)[0]

    def test_without_an_inventory_only_a_named_model_is_judged(self):
        assert gpu_warnings("DynamicBind", "any", DYNAMICBIND, {}) == []
        assert gpu_warnings("DynamicBind", "80GB", DYNAMICBIND, {}) == []

    def test_an_exclusion_leaves_what_is_not_excluded(self):
        assert "may land on H100, H200" in gpu_warnings("DynamicBind", "!L4", DYNAMICBIND, S3IT)[0]
        assert gpu_warnings("DynamicBind", "!H100", DYNAMICBIND, S3IT)[0].count("H200") == 1
        assert "H100" not in gpu_warnings("DynamicBind", "!H100", DYNAMICBIND, S3IT)[0].split("may land on")[1].split(",")[0]

    def test_a_memory_request_may_land_on_a_model_whose_memory_is_not_stated(self):
        """H100 has no memory in the inventory, so any memory request can be answered by it; H200 is 140 GB."""
        assert "may land on H100" in gpu_warnings("DynamicBind", "80GB", DYNAMICBIND, S3IT)[0]
        assert "H200" not in gpu_warnings("DynamicBind", "80GB", DYNAMICBIND, S3IT)[0]
        assert "H200" in gpu_warnings("DynamicBind", "140GB", DYNAMICBIND, S3IT)[0]

    def test_alternatives_are_judged_as_a_union(self):
        assert "H200" in gpu_warnings("DynamicBind", "80GB|140GB", DYNAMICBIND, S3IT)[0]
        assert gpu_warnings("DynamicBind", "A100|L4", DYNAMICBIND, S3IT) == []

    def test_high_memory_leaves_out_the_models_it_cannot_mean(self):
        named, possible = resolve_request("high-memory", {"L4": "24GB", "A100": "80GB", "H200": "140GB"})
        assert named == [] and possible == ["A100"]

    def test_a_request_that_cannot_reach_an_incompatible_model_is_silent(self):
        assert gpu_warnings("DynamicBind", "A100", DYNAMICBIND, {"A100": "80GB", "H100": "80GB"}) == []
        assert gpu_warnings("DynamicBind", "24GB", DYNAMICBIND, {"L4": "24GB", "H100": "80GB"}) == []


class TestACompatibleList:

    ONLY_A100 = {"gpu": {"compatible": ["A100"]}}

    def test_every_model_outside_the_list_warns(self):
        assert len(gpu_warnings("Tool", "H100", self.ONLY_A100, S3IT)) == 1
        assert len(gpu_warnings("Tool", "L4", self.ONLY_A100, S3IT)) == 1
        assert "may land on L4, H100, H200" in gpu_warnings("Tool", "any", self.ONLY_A100, S3IT)[0]

    def test_a_model_inside_the_list_is_silent(self):
        assert gpu_warnings("Tool", "A100", self.ONLY_A100, S3IT) == []

    def test_the_fix_is_the_first_compatible_model_when_no_use_is_given(self):
        assert 'Resources(gpu="A100")' in gpu_warnings("Tool", "H100", self.ONLY_A100, S3IT)[0]

    def test_the_fix_is_the_first_listed_model_not_the_first_alphabetically(self):
        declared = {"gpu": {"compatible": ["V100", "A100"]}}
        assert 'Resources(gpu="V100")' in gpu_warnings("Tool", "H100", declared, S3IT)[0]


class TestASingleGpuTypeSite:
    """On a node-exclusive site the scheduler emits an untyped `gpu:N`, so every request lands on the inventory."""

    DAINT = {"GH200": ""}
    POCKETGEN = {"gpu": {"incompatible": ["GH200"], "use": "A100"}}

    @pytest.mark.parametrize("request_", ["gh", "any", "A100", "4"])
    def test_every_request_is_judged_against_the_inventory(self, request_):
        messages = gpu_warnings("PocketGen", request_, self.POCKETGEN, self.DAINT, single_gpu_type=True)
        assert len(messages) == 1 and "may land on GH200" in messages[0]

    def test_none_is_still_not_judged(self):
        assert gpu_warnings("PocketGen", "none", self.POCKETGEN, self.DAINT, single_gpu_type=True) == []

    def test_a_typed_site_still_ignores_an_unknown_token(self):
        assert gpu_warnings("PocketGen", "gh", self.POCKETGEN, self.DAINT) == []


class TestTheDeclarationIsChecked:

    def test_a_misspelled_key_is_refused_instead_of_ignored(self):
        with pytest.raises(ValueError, match="incompatable"):
            validate_requirements("Tool", {"gpu": {"incompatable": ["H100"]}})

    def test_a_model_list_must_be_a_list_of_names(self):
        with pytest.raises(ValueError, match="list of GPU model names"):
            validate_requirements("Tool", {"gpu": {"incompatible": "H100"}})

    def test_a_tool_class_with_a_bad_declaration_cannot_be_defined(self):
        with pytest.raises(ValueError, match="RESOURCE_REQUIREMENTS"):
            class Broken(BaseConfig):
                TOOL_NAME = "Broken"
                RESOURCE_REQUIREMENTS = ["H100"]

    def test_other_resource_keys_are_left_alone(self):
        validate_requirements("Tool", {"gpu": {"incompatible": ["H100"]}, "memory": "64GB"})


class TestTheInventoryOfTheConfig:

    def _inventory(self, monkeypatch, machine):
        monkeypatch.setattr(ConfigManager, "_get_machine_config", lambda self: machine)
        return ConfigManager().get_gpu_models()

    def test_an_absent_key_is_an_empty_inventory(self, monkeypatch):
        assert self._inventory(monkeypatch, {}) == {}

    def test_a_model_without_memory_maps_to_an_empty_string(self, monkeypatch):
        assert self._inventory(monkeypatch, {"gpu_models": {"A100": None, "H200": "140GB"}}) == {"A100": "", "H200": "140GB"}

    def test_a_list_is_refused_with_the_shape_it_wants(self, monkeypatch):
        with pytest.raises(ValueError, match="machine.gpu_models must map"):
            self._inventory(monkeypatch, {"gpu_models": ["A100"]})

    def test_the_shipped_configs_declare_what_the_repo_documents(self):
        import yaml

        root = Path(__file__).parent.parent
        cluster = yaml.safe_load((root / "config.cluster.yaml").read_text(encoding="utf-8"))["machine"]["gpu_models"]
        daint = yaml.safe_load((root / "config.daint.yaml").read_text(encoding="utf-8"))["machine"]["gpu_models"]
        assert set(cluster) == {"L4", "A100", "H100", "H200"}
        assert cluster["H200"] == "140GB", "the only size the repo states, in the SLURM backend"
        assert set(daint) == {"GH200"}


# ── through the batch-script generation ──────────────────────────────────────

@pytest.fixture
def slurm(monkeypatch):
    config_path = FIXTURES_DIR / "config.slurm_local.yaml"
    ConfigManager._instance = None
    ConfigManager._config = None
    ConfigManager._variant = None
    monkeypatch.setattr(ConfigManager, "_get_config_path", classmethod(lambda cls, variant=None: str(config_path)))
    monkeypatch.setattr(ConfigManager, "get_gpu_models", lambda self: dict(S3IT))
    from biopipelines.mock import Mock
    monkeypatch.setattr(Mock, "RESOURCE_REQUIREMENTS", DYNAMICBIND)
    yield
    ConfigManager._instance = None
    ConfigManager._config = None
    ConfigManager._variant = None


def _generate(job, requests):
    """Generate the scripts of a pipeline with one Mock step per GPU request, each in its own batch."""
    from biopipelines.mock import Mock
    from biopipelines.pipeline import Pipeline, Resources

    pipeline = Pipeline(project="TestSuite", job=job, description="requirements", on_the_fly=False,
                        local_output=True, config="slurm_local")
    with pipeline:
        for index, request in enumerate(requests):
            Resources(gpu=request)
            Mock(ids=[f"s{index}"], streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
        pipeline.save()
        pipeline.generate_job_scripts()
    return pipeline


class TestTheWarningWhenTheScriptsAreGenerated:

    def test_a_request_that_may_land_on_an_unsupported_gpu_warns(self, slurm, isolated_cwd, capsys):
        _generate("any", ["any"])

        out = capsys.readouterr().out
        assert "Warning: Mock: Resources(gpu=\"any\") may land on H100, H200" in out
        assert 'Resources(gpu="A100")' in out

    def test_a_supported_request_adds_no_warning(self, slurm, isolated_cwd, capsys):
        _generate("a100", ["A100"])

        assert "does not support" not in capsys.readouterr().out

    def test_each_batch_is_judged_by_its_own_request_once(self, slurm, isolated_cwd, capsys):
        _generate("batches", ["A100", "H100", "any"])

        lines = [line for line in capsys.readouterr().out.splitlines() if "does not support" in line]
        assert len(lines) == 2
        assert 'Resources(gpu="H100")' in lines[0] and 'Resources(gpu="any")' in lines[1]

    def test_a_tool_that_appears_twice_in_a_batch_is_named_once(self, slurm, isolated_cwd, capsys):
        from biopipelines.mock import Mock
        from biopipelines.pipeline import Pipeline, Resources

        pipeline = Pipeline(project="TestSuite", job="twice", description="requirements", on_the_fly=False,
                            local_output=True, config="slurm_local")
        with pipeline:
            Resources(gpu="H100")
            Mock(ids=["a"], streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
            Mock(ids=["b"], streams={"out": {"format": "pdb", "file": "<id>.pdb"}})
            pipeline.save()
            pipeline.generate_job_scripts()

        assert capsys.readouterr().out.count("does not support H100") == 1

    def test_a_node_exclusive_site_warns_for_its_own_request(self, slurm, isolated_cwd, capsys, monkeypatch):
        """Daint pipelines ask for `gpu="gh"`, which names no model, and every node is a GH200."""
        monkeypatch.setattr(ConfigManager, "get_node_exclusive", lambda self: True)
        monkeypatch.setattr(ConfigManager, "get_gpu_models", lambda self: {"GH200": ""})
        from biopipelines.mock import Mock
        monkeypatch.setattr(Mock, "RESOURCE_REQUIREMENTS", {"gpu": {"incompatible": ["GH200"], "use": "A100"}})
        _generate("daint", ["gh"])

        assert 'Warning: Mock: Resources(gpu="gh") may land on GH200' in capsys.readouterr().out

    def test_the_generated_scripts_do_not_depend_on_the_warning(self, slurm, isolated_cwd, monkeypatch):
        from biopipelines.mock import Mock

        def scripts(pipeline):
            runtime = Path(pipeline.folders["runtime"])
            return {path.name: path.read_text(encoding="utf-8").replace(str(runtime), "<RUNTIME>")
                    for path in sorted(runtime.glob("slurm*.sh"))}

        warned = scripts(_generate("same_a", ["any"]))
        monkeypatch.setattr(Mock, "RESOURCE_REQUIREMENTS", {})
        silent = scripts(_generate("same_b", ["any"]))

        normalize = lambda texts: {name: text.replace("same_a", "same").replace("same_b", "same") for name, text in texts.items()}
        assert normalize(warned) == normalize(silent)
