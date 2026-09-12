"""A parameter's alternative spelling must bind the parameter, not fall through to **kwargs.

Four near-misses between our own sibling tools, all of which failed silently:

* `RFdiffusion3(contigs=...)` — its three sibling wrappers take `contigs`, RFdiffusion3 types `contig`. Since RFdiffusion3 forwards untyped kwargs to its hydra CLI, the plural was reported as a probable misspelling and then passed through as a bogus override: `self.contig` stayed `''`, the motif was discarded and the job ran to completion producing unconditioned de-novo backbones.
* `Ligand(code=...)` vs `Ligand(codes=...)` — one letter apart, mutually exclusive, different objects. `code` is now `codes`.
* `Fuse(name=...)` and `Table(name=...)` — `name` is framework-reserved (`BaseConfig` reads it as the job name), so declaring it as their own parameter meant the job name was silently lost. Their own parameters are now `prefix` and `table_name`, and `name` reaches the framework again.

The mechanism is `BaseConfig.PARAMETER_ALIASES`, resolved before the tool's own `__init__` binds anything — which is what keeps an alias out of the leftover-kwargs logic that would otherwise report it as a typo or forward it. `DEPRECATED_ALIASES` marks which of those spellings are retired rather than first-class, so a synonym stays silent while a retired spelling reports one line through `contract_enforcement`.

`Structure` / `Compound` (the PDB / Ligand aliases) are covered at the bottom: those are whole-tool aliases, not parameter ones, and deliberately carry no deprecation.
"""

import os
import pathlib

import pytest

from biopipelines import contract_enforcement as ce


def _step_script(script_path, tool_name):
    runtime = pathlib.Path(script_path).parent
    matches = sorted(runtime.glob(f"*_{tool_name}.sh"))
    assert matches, f"no {tool_name} step script under {runtime}"
    return matches[0].read_text(encoding="utf-8")


# ── the mechanism ─────────────────────────────────────────────────────────────

def test_an_alias_target_that_no_parameter_names_is_refused_at_class_creation():
    """An alias landing back in **kwargs would be the silence this whole mechanism removes, so it fails when the class is defined rather than when a user calls it."""
    from biopipelines.base_config import BaseConfig

    with pytest.raises(ValueError, match="PARAMETER_ALIASES"):
        class BrokenAliasTool(BaseConfig):
            TOOL_NAME = "BrokenAliasTool"
            PARAMETER_ALIASES = {"olde": "no_such_parameter"}

            def __init__(self, real=None, **kwargs):
                super().__init__(**kwargs)


def test_a_deprecation_for_an_undeclared_alias_is_refused_at_class_creation():
    from biopipelines.base_config import BaseConfig

    with pytest.raises(ValueError, match="DEPRECATED_ALIASES"):
        class MislabelledAliasTool(BaseConfig):
            TOOL_NAME = "MislabelledAliasTool"
            PARAMETER_ALIASES = {"olde": "real"}
            DEPRECATED_ALIASES = ("olde", "never_mapped")

            def __init__(self, real=None, **kwargs):
                super().__init__(**kwargs)


# ── item 1: contigs on RFdiffusion3 ───────────────────────────────────────────

def test_contigs_binds_contig_on_rfdiffusion3(local_config, record_case):
    from biopipelines.rfdiffusion3 import RFdiffusion3

    tool = RFdiffusion3(contigs="A50-100,80-100")
    record_case(input="RFdiffusion3(contigs='A50-100,80-100')",
                expected="A50-100,80-100", actual=tool.contig)
    assert tool.contig == "A50-100,80-100"


def test_contigs_is_not_forwarded_to_hydra(local_config, record_case):
    """RFdiffusion3 renders untyped kwargs as `key=value` overrides; a bound alias must not be among them."""
    from biopipelines.rfdiffusion3 import RFdiffusion3

    tool = RFdiffusion3(contigs="A50-100,80-100")
    record_case(input="RFdiffusion3(contigs=...).extra_args",
                expected={}, actual=tool.extra_args)
    assert tool.extra_args == {}
    assert tool.extra_args_tokens() == []


def test_contigs_reaches_the_emitted_script_as_the_contig(
        local_config, isolated_cwd, new_pipeline):
    """End to end: the motif must be in the artifacts, and no `contigs=` override on the command line."""
    from biopipelines.rfdiffusion3 import RFdiffusion3

    pipeline = new_pipeline("rfd3_contigs_alias")
    with pipeline:
        RFdiffusion3(contigs="A50-100,80-100", num_designs=1)
        script_path = pipeline.save()

    template = pathlib.Path(pipeline.tools[0].json_template).read_text(encoding="utf-8")
    assert '"contig": "A50-100,80-100"' in template
    assert "contigs=" not in _step_script(script_path, "RFdiffusion3")


def test_contig_still_works_on_rfdiffusion3(local_config):
    from biopipelines.rfdiffusion3 import RFdiffusion3

    assert RFdiffusion3(contig="A50-100,80-100").contig == "A50-100,80-100"


def test_contigs_is_a_synonym_and_says_nothing(local_config, capsys):
    """`contigs` is a first-class synonym, not a retired spelling: warning on every use would be noise."""
    from biopipelines.rfdiffusion3 import RFdiffusion3

    ce.reset_reported()
    RFdiffusion3(contigs="A50-100,80-100")
    err = capsys.readouterr().err
    assert "deprecated" not in err
    assert "did you mean" not in err


def test_both_spellings_at_once_is_refused(local_config):
    from biopipelines.rfdiffusion3 import RFdiffusion3

    with pytest.raises(ValueError, match="another spelling"):
        RFdiffusion3(contig="A1-50", contigs="A1-50")


# ── item 2: Ligand codes ────────────────────────────────────────────────

def test_codes_builds_a_code_only_ligand(local_config, record_case):
    from biopipelines.ligand import Ligand

    lig = Ligand(codes="ZIT")
    record_case(input="Ligand(codes='ZIT')",
                expected=(True, ["ZIT"]), actual=(lig.code_only, lig.residue_codes))
    assert lig.code_only is True
    assert lig.residue_codes == ["ZIT"]


def test_the_retired_code_spelling_still_builds_the_same_ligand(local_config, record_case):
    from biopipelines.ligand import Ligand

    lig = Ligand(code="ZIT")
    record_case(input="Ligand(code='ZIT') — retired spelling",
                expected=(True, ["ZIT"]), actual=(lig.code_only, lig.residue_codes))
    assert lig.code_only is True
    assert lig.residue_codes == ["ZIT"]


def test_the_retired_code_spelling_reports_a_deprecation(local_config, capsys, record_case):
    from biopipelines.ligand import Ligand

    ce.reset_reported()
    Ligand(code="ZIT")
    err = capsys.readouterr().err
    record_case(input="Ligand(code='ZIT') stderr",
                expected="deprecated; use codes=", actual=err.strip()[:80])
    assert "deprecated_alias" in err
    assert "code= is a synonym for codes= and will soon be deprecated" in err


def test_the_retired_code_spelling_is_not_reported_as_a_typo(local_config, capsys):
    """Before the alias existed the leftover-kwargs check called `code` a misspelling of `codes` — the wrong diagnosis, since the two are different objects."""
    from biopipelines.ligand import Ligand

    ce.reset_reported()
    Ligand(code="ZIT")
    err = capsys.readouterr().err
    assert "did you mean" not in err
    assert "unknown constructor parameter" not in err


def test_codes_still_carries_residue_codes_on_a_real_ligand(local_config, record_case):
    """The parameter that stayed: `codes` names the residue on a ligand that has chemistry."""
    from biopipelines.ligand import Ligand

    lig = Ligand(smiles="CCO", ids="etoh", codes="LIG")
    record_case(input="Ligand(smiles='CCO', codes='LIG')",
                expected=(False, ["LIG"]), actual=(lig.code_only, lig.residue_codes))
    assert lig.code_only is False
    assert lig.residue_codes == ["LIG"]


def test_both_ligand_spellings_at_once_is_refused(local_config):
    from biopipelines.ligand import Ligand

    with pytest.raises(ValueError, match="another spelling"):
        Ligand(code="ZIT", codes="ZIT")


# ── item 3: name is the framework's, on Fuse and Table ────────────────────────

def _two_sequence_streams():
    from biopipelines.datastream import DataStream

    return [DataStream(name="sequences", ids=[i], files=[], map_table=f"{i}.csv",
                       format="csv") for i in ("a", "b")]


def test_fuse_prefix_labels_the_construct(local_config, record_case):
    from biopipelines.fuse import Fuse

    fuse = Fuse(_two_sequence_streams(), prefix="mybinder")
    record_case(input="Fuse(..., prefix='mybinder')",
                expected="mybinder", actual=fuse.prefix)
    assert fuse.prefix == "mybinder"
    assert fuse._get_job_base() == "mybinder"


def test_fuse_name_is_only_the_job_name(local_config, record_case):
    """Fuse has no alias from `name`: the framework keeps that key, so `name` means the job name here exactly as it does on every other tool, and the construct label has its own spelling."""
    from biopipelines.fuse import Fuse

    fuse = Fuse(_two_sequence_streams(), name="jobname")
    record_case(input="Fuse(..., name='jobname')",
                expected=("", "jobname"), actual=(fuse.prefix, fuse.job_name))
    assert fuse.prefix == ""
    assert fuse.job_name == "jobname"


def test_fuse_prefix_and_name_are_independent(local_config):
    from biopipelines.fuse import Fuse

    fuse = Fuse(_two_sequence_streams(), prefix="label", name="jobname")
    assert (fuse.prefix, fuse.job_name) == ("label", "jobname")


def test_fuse_declares_no_parameter_aliases(local_config):
    """Pins the decision: Fuse deliberately has no alias mechanism, so `name` cannot be captured away from the framework again."""
    from biopipelines.fuse import Fuse

    assert Fuse.PARAMETER_ALIASES == {}
    assert Fuse.DEPRECATED_ALIASES == ()


def _csv(tmp_path):
    import pandas as pd

    path = tmp_path / "metrics.csv"
    pd.DataFrame({"id": ["a"], "plddt": [90.0]}).to_csv(path, index=False)
    return str(path)


def test_table_table_name_identifies_the_table(local_config, tmp_path, record_case):
    from biopipelines.table import Table

    table = Table(_csv(tmp_path), table_name="metrics")
    record_case(input="Table(..., table_name='metrics')",
                expected="metrics", actual=table.table_name)
    assert table.table_name == "metrics"


def test_table_name_sets_the_job_name_and_still_identifies_the_table(
        local_config, tmp_path, record_case):
    from biopipelines.table import Table

    table = Table(_csv(tmp_path), name="metrics")
    record_case(input="Table(..., name='metrics')",
                expected=("metrics", "metrics"),
                actual=(table.table_name, table.job_name))
    assert table.table_name == "metrics"
    assert table.job_name == "metrics"


def test_table_name_alongside_table_name_is_only_the_job_name(local_config, tmp_path):
    from biopipelines.table import Table

    table = Table(_csv(tmp_path), table_name="metrics", name="jobname")
    assert (table.table_name, table.job_name) == ("metrics", "jobname")


def test_table_name_reports_a_deprecation(local_config, tmp_path, capsys):
    from biopipelines.table import Table

    ce.reset_reported()
    Table(_csv(tmp_path), name="metrics")
    err = capsys.readouterr().err
    assert "name= is a synonym for table_name= and will soon be deprecated" in err


def test_the_deprecated_table_name_is_still_shell_safety_checked(local_config, tmp_path):
    """The alias must not become a way past a validator: the value is checked wherever it was bound from."""
    from biopipelines.table import Table

    with pytest.raises(ValueError):
        Table(_csv(tmp_path), name='bad"name')


# ── the deprecation line is severity-gated like every other contract ──────────

def test_a_deprecation_can_be_silenced(local_config, monkeypatch, capsys):
    from biopipelines.ligand import Ligand

    monkeypatch.setenv("BIOPIPELINES_ENFORCE_DEPRECATED_ALIAS", "off")
    ce.reset_reported()
    Ligand(code="ZIT")
    err = capsys.readouterr().err
    # Only this check is silenced; an unrelated one on the same call still reports.
    assert "deprecated_alias" not in err


def test_a_deprecation_can_be_promoted_to_an_error(local_config, monkeypatch):
    from biopipelines.ligand import Ligand

    monkeypatch.setenv("BIOPIPELINES_ENFORCE_DEPRECATED_ALIAS", "raise")
    ce.reset_reported()
    with pytest.raises(ce.ContractViolation, match="deprecated"):
        Ligand(code="ZIT")


# ── item 5: Structure and Compound are the same tools under stream names ──────

def test_structure_is_pdb_and_compound_is_ligand():
    import biopipelines as bp

    assert bp.Structure is bp.PDB
    assert bp.Compound is bp.Ligand


def test_both_spellings_are_exported():
    import biopipelines as bp

    for name in ("PDB", "Structure", "Ligand", "Compound"):
        assert name in bp.__all__, f"{name} missing from biopipelines.__all__"
        assert getattr(bp, name, None) is not None


def test_an_alias_defines_no_second_tool_name():
    """A bare alias — not a subclass — is what keeps the registry, the docs index and the config's environments/folders keys single-entry for one tool."""
    import biopipelines as bp

    assert bp.Structure.TOOL_NAME == "PDB"
    assert bp.Compound.TOOL_NAME == "Ligand"


def test_a_pipeline_built_through_structure_keeps_the_pdb_folder_name(
        local_config, isolated_cwd, new_pipeline, record_case):
    """The documented cost of keeping TOOL_NAME: outputs are named for the original tool whichever spelling was used."""
    import biopipelines as bp

    pdb_file = isolated_cwd / "M0584_1ldm.pdb"
    pdb_file.write_text("HEADER test\nEND\n")

    pipeline = new_pipeline("structure_alias")
    with pipeline:
        bp.Structure(pdbs=str(pdb_file))
        pipeline.save()

    folder = os.path.basename(pipeline.tools[0].output_folder)
    record_case(input="Structure(pdbs=...) step folder", expected="001_PDB", actual=folder)
    assert folder == "001_PDB"


def test_a_pipeline_built_through_compound_keeps_the_ligand_folder_name(
        local_config, isolated_cwd, new_pipeline):
    import biopipelines as bp

    pipeline = new_pipeline("compound_alias")
    with pipeline:
        bp.Compound(codes="ZIT")
        pipeline.save()

    assert os.path.basename(pipeline.tools[0].output_folder) == "001_Ligand"
