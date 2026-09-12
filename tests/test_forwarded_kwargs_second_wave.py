"""The second wave of forwarding tools: an untyped kwarg must reach the emitted step script.

`tests/test_forwarded_kwargs.py` pins the renderer, the shell-safety layer and the checks. What it cannot pin for a newly converted tool is the only thing that makes the marker worth setting: that the rendered tokens are actually joined onto the upstream command line the wrapper emits, rather than merely stored on `self.extra_args`. Every test here therefore greps the generated `NNN_<Tool>.sh`, not `pipeline.sh` and not the config folder.

The seven tools covered are the ones whose wrapper's own bash carries the upstream command: AlphaFold (colabfold_batch), LASErMPNN (run_batch_inference, through its free-form `--run-options` channel), Boltz2 (`boltz predict`), NeuralPLexer (neuralplexer-inference), Aggrescan3D (`aggrescan`), CABSflex (`CABSflex`) and DynamicBind (run_single_protein_inference.py).

Known limit, deliberate: `render_extra_args` preserves underscores, so a dash-spelled upstream flag (`--msa-mode`, `--weighted-fit`) has no kwarg spelling — `_SAFE_EXTRA_KEY_RE` refuses a key with a dash in it. On the dash-style CLIs here (colabfold_batch, neuralplexer-inference, part of CABSflex) only single-word flags are reachable.
"""

import pathlib

import pytest


def _step_script(script_path, tool_name):
    runtime = pathlib.Path(script_path).parent
    matches = sorted(runtime.glob(f"*_{tool_name}.sh"))
    assert matches, f"no {tool_name} step script under {runtime}"
    return matches[0].read_text(encoding="utf-8", errors="replace")


def _structures_mock():
    from biopipelines.mock import Mock
    return Mock(
        ids=["s1"],
        streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
        map_table_strategy="config",
    )


# ── the six second-wave tools ─────────────────────────────────────────────────

def test_alphafold_forwards_onto_colabfold_batch(local_config, isolated_cwd, new_pipeline, record_case):
    from biopipelines.alphafold import AlphaFold
    from biopipelines.sequence import Sequence

    pipeline = new_pipeline("fwd_alphafold")
    with pipeline:
        s = Sequence(seq="MKTAYIAKQRQISFVKSHFSRQLEERLGL", type="protein", ids="p1")
        AlphaFold(proteins=s, templates=True)
        script = _step_script(pipeline.save(), "AlphaFold")

    line = next(l for l in script.splitlines() if "colabfold_batch" in l and "--templates" in l)
    record_case(input="AlphaFold(templates=True)", expected='"--templates" on the colabfold_batch line',
                actual=line.strip()[:120])
    assert '"--templates"' in line


def test_alphafold_forwards_a_flag_with_a_value(local_config, isolated_cwd, new_pipeline):
    from biopipelines.alphafold import AlphaFold
    from biopipelines.sequence import Sequence

    pipeline = new_pipeline("fwd_alphafold_value")
    with pipeline:
        s = Sequence(seq="MKTAYIAKQRQISFVKSHFSRQLEERLGL", type="protein", ids="p1")
        AlphaFold(proteins=s, host_url="https://example.org/api")
        script = _step_script(pipeline.save(), "AlphaFold")

    assert '"--host_url" "https://example.org/api"' in script


def test_alphafold_options_echo_stays_readable(local_config, isolated_cwd, new_pipeline):
    """The forwarded tokens are announced by their own echo; the typed "Options:" echo keeps the bash quoting out."""
    from biopipelines.alphafold import AlphaFold
    from biopipelines.sequence import Sequence

    pipeline = new_pipeline("fwd_alphafold_echo")
    with pipeline:
        s = Sequence(seq="MKTAYIAKQRQISFVKSHFSRQLEERLGL", type="protein", ids="p1")
        AlphaFold(proteins=s, templates=True)
        script = _step_script(pipeline.save(), "AlphaFold")

    assert "Forwarding unrecognized argument(s) to AlphaFold: --templates" in script
    options_echo = next(l for l in script.splitlines() if l.startswith('echo "Options:'))
    assert "--templates" not in options_echo


def test_lasermpnn_forwards_through_run_options(local_config, isolated_cwd, new_pipeline, record_case):
    """LASErMPNN's channel is the free-form `--run-options` string the pipe script shlex-splits into argv."""
    from biopipelines.lasermpnn import LASErMPNN

    pipeline = new_pipeline("fwd_lasermpnn")
    with pipeline:
        m = _structures_mock()
        LASErMPNN(structures=m.streams.structures, bb_noise=0.1)
        script = _step_script(pipeline.save(), "LASErMPNN")

    line = next(l for l in script.splitlines() if "--run-options" in l)
    record_case(input="LASErMPNN(bb_noise=0.1)", expected="--bb_noise 0.1 inside --run-options",
                actual=line.strip()[:160])
    assert "--bb_noise 0.1" in line


def test_lasermpnn_quotes_a_forwarded_value_that_contains_a_space(local_config, isolated_cwd, new_pipeline):
    """shlex.split would otherwise turn one value into two argv tokens."""
    import shlex

    from biopipelines.lasermpnn import LASErMPNN

    pipeline = new_pipeline("fwd_lasermpnn_space")
    with pipeline:
        m = _structures_mock()
        LASErMPNN(structures=m.streams.structures, note="two words")
        script = _step_script(pipeline.save(), "LASErMPNN")

    line = next(l for l in script.splitlines() if "--run-options" in l)
    run_options = line.split("--run-options", 1)[1].strip().rstrip("\\").strip()
    argv = shlex.split(shlex.split(run_options)[0])
    assert argv[-2:] == ["--note", "two words"]


def test_boltz2_forwards_onto_boltz_predict(local_config, isolated_cwd, new_pipeline, record_case):
    from biopipelines.boltz2 import Boltz2
    from biopipelines.sequence import Sequence

    pipeline = new_pipeline("fwd_boltz2")
    with pipeline:
        s = Sequence(seq="MKTAYIAKQRQISFVKSHFSRQLEERLGL", type="protein", ids="p1")
        Boltz2(proteins=s, step_scale=1.5)
        script = _step_script(pipeline.save(), "Boltz2")

    line = next(l for l in script.splitlines() if "boltz predict" in l)
    record_case(input="Boltz2(step_scale=1.5)", expected='"--step_scale" "1.5"',
                actual=line.strip()[:160])
    assert '"--step_scale" "1.5"' in line


def test_neuralplexer_forwards_onto_the_inference_command(local_config, isolated_cwd, new_pipeline, record_case):
    from biopipelines.ligand import Ligand
    from biopipelines.neuralplexer import NeuralPLexer

    pipeline = new_pipeline("fwd_neuralplexer")
    with pipeline:
        s = _structures_mock()
        NeuralPLexer(structures=s, compounds=Ligand(smiles="CCO", ids="c1", codes="ETH"),
                     **{"seed": 7})
        script = _step_script(pipeline.save(), "NeuralPLexer")

    block = script.split("neuralplexer-inference", 1)[1].split("rc=$?", 1)[0]
    # One token per continuation line, the shape RFdiffusion3 already emits.
    tokens = [l.strip().rstrip("\\").strip() for l in block.splitlines()]
    record_case(input="NeuralPLexer(seed=7)", expected='"--seed" and "7" as argv tokens',
                actual=tokens[-3:])
    assert '"--seed"' in tokens and '"7"' in tokens
    assert tokens.index('"--seed"') + 1 == tokens.index('"7"')


def test_aggrescan3d_forwards_onto_the_aggrescan_binary(local_config, isolated_cwd, new_pipeline, record_case):
    from biopipelines.aggrescan3d import Aggrescan3D

    pipeline = new_pipeline("fwd_aggrescan3d")
    with pipeline:
        s = _structures_mock()
        Aggrescan3D(structures=s.streams.structures, dynamic=True)
        script = _step_script(pipeline.save(), "Aggrescan3D")

    line = next(l for l in script.splitlines() if "aggrescan -i" in l)
    record_case(input="Aggrescan3D(dynamic=True)", expected='"--dynamic" on the aggrescan line',
                actual=line.strip()[:140])
    assert '"--dynamic"' in line


def test_aggrescan3d_forwards_in_the_parallel_branch_too(local_config, isolated_cwd, new_pipeline):
    """Two structures with max_parallel>1 emit the other of the wrapper's two run blocks; both read the same flags string."""
    from biopipelines.mock import Mock
    from biopipelines.aggrescan3d import Aggrescan3D

    pipeline = new_pipeline("fwd_aggrescan3d_par")
    with pipeline:
        m = Mock(ids=["s1", "s2"],
                 streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
                 map_table_strategy="config")
        Aggrescan3D(structures=m.streams.structures, max_parallel=2, dynamic=True)
        script = _step_script(pipeline.save(), "Aggrescan3D")

    assert "PIDS=()" in script  # the parallel branch
    assert '"--dynamic"' in next(l for l in script.splitlines() if "aggrescan -i" in l)


def test_cabsflex_forwards_onto_the_cabsflex_binary(local_config, isolated_cwd, new_pipeline, record_case):
    from biopipelines.cabsflex import CABSflex

    pipeline = new_pipeline("fwd_cabsflex")
    with pipeline:
        m = _structures_mock()
        CABSflex(structures=m.streams.structures, verbose=3)
        script = _step_script(pipeline.save(), "CABSflex")

    line = next(l for l in script.splitlines() if "CABSflex -i" in l)
    record_case(input="CABSflex(verbose=3)", expected='"--verbose" "3"', actual=line.strip()[:140])
    assert '"--verbose" "3"' in line


def test_dynamicbind_forwards_onto_the_inference_script(local_config, isolated_cwd, new_pipeline, record_case):
    from biopipelines.dynamicbind import DynamicBind
    from biopipelines.ligand import Ligand

    pipeline = new_pipeline("fwd_dynamicbind")
    with pipeline:
        s = _structures_mock()
        DynamicBind(structures=s, compounds=Ligand(smiles="CCO", ids="c1", codes="ETH"),
                    no_final_step_noise=True)
        script = _step_script(pipeline.save(), "DynamicBind")

    block = script.split("run_single_protein_inference", 1)[1].split("rc=$?", 1)[0]
    tokens = [l.strip().rstrip("\\").strip() for l in block.splitlines()]
    record_case(input="DynamicBind(no_final_step_noise=True)",
                expected='"--no_final_step_noise" as an argv token', actual=tokens[-3:])
    assert '"--no_final_step_noise"' in tokens


# ── the shared guarantees, once per tool ──────────────────────────────────────

_TOOLS = ["AlphaFold", "LASErMPNN", "Boltz2", "NeuralPLexer", "Aggrescan3D",
          "CABSflex", "DynamicBind"]


@pytest.mark.parametrize("tool_name", _TOOLS)
def test_each_second_wave_tool_declares_the_argparse_dialect(tool_name):
    """Every upstream command here is `--flag value`, so none of them may claim the hydra dialect."""
    import biopipelines  # noqa: F401  (imports every wrapper)
    from biopipelines.base_config import BaseConfig

    def find(cls):
        for sub in cls.__subclasses__():
            if sub.TOOL_NAME == tool_name and sub.__module__.startswith("biopipelines."):
                return sub
            found = find(sub)
            if found is not None:
                return found
        return None

    tool = find(BaseConfig)
    assert tool is not None, f"{tool_name} not found among BaseConfig subclasses"
    assert tool.FORWARD_UNKNOWN_KWARGS == "argparse"


def test_a_dash_spelled_upstream_flag_is_reachable():
    """This used to be a documented limit: the key allowlist took only an identifier, so `--msa-mode`, `--weighted-fit` and `--num-steps` could not be forwarded to any tool at all. That silently confined forwarding to whichever upstream parsers happen to use underscores, and NeuralPLexer's surface is almost entirely dashed. An inner dash is now allowed; a leading one is still refused, so a forwarded key cannot impersonate one of the framework's own flags."""
    from biopipelines.base_config import _validate_extra_arg, render_extra_args

    for key in ("msa-mode", "weighted-fit", "num-steps"):
        _validate_extra_arg(key, "value")
    assert render_extra_args({"msa-mode": "single_sequence"}, "argparse") == [
        "--msa-mode", "single_sequence"]
    with pytest.raises(ValueError):
        _validate_extra_arg("--msa-mode", "value")
