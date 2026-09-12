"""An untyped constructor kwarg must either reach the upstream command line or be named as a typo.

Every tool constructor ends in ``**kwargs`` and ``BaseConfig.__init__`` parked whatever was left in ``self.params`` — which nothing ever read. An advanced user reaching for an upstream flag the wrapper does not type got exactly the same silence as someone who misspelled a real parameter: the value vanished between Python and bash.

Two behaviours are pinned here:

* a tool that sets ``FORWARD_UNKNOWN_KWARGS`` renders the leftover keys onto its upstream command line through one shared renderer, after the same shell-safety validation every other user string goes through;
* a leftover key that is lexically close to a real parameter name is reported as a probable misspelling naming the candidate, whether or not the tool forwards — because a near-miss on a forwarding tool is still almost certainly a typo, and forwarding it only moves the failure onto a compute node.
"""

import pathlib

import os
import pytest

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

from biopipelines import contract_enforcement as ce
from biopipelines.base_config import (
    EXTRA_ARGS_DIALECTS, _validate_extra_arg, render_extra_args,
)


# ── the renderer ──────────────────────────────────────────────────────────────

def test_argparse_dialect_renders_flag_and_value(record_case):
    actual = render_extra_args({"omit_AAs": "CX"}, "argparse")
    record_case(input="omit_AAs='CX' (argparse)", expected=["--omit_AAs", "CX"], actual=actual)
    assert actual == ["--omit_AAs", "CX"]


def test_hydra_dialect_renders_a_dotted_override():
    assert render_extra_args({"inference": "x"}, "hydra") == ["inference=x"]


def test_underscores_are_preserved():
    """Both upstream families take the underscore verbatim; rewriting it to a dash would invent a flag that does not exist."""
    assert render_extra_args({"num_seq_per_target": 8}, "argparse") == ["--num_seq_per_target", "8"]
    assert render_extra_args({"seed_offset": 3}, "hydra") == ["seed_offset=3"]


def test_true_is_a_bare_flag_under_argparse():
    assert render_extra_args({"score_only": True}, "argparse") == ["--score_only"]


def test_true_is_explicit_under_hydra(record_case):
    """A bare token is not a valid hydra override, so the value has to be spelled out."""
    actual = render_extra_args({"inference.cyclic": True}, "hydra")
    record_case(input="hydra True", expected=["inference.cyclic=True"], actual=actual)
    assert actual == ["inference.cyclic=True"]


@pytest.mark.parametrize("dialect", EXTRA_ARGS_DIALECTS)
def test_false_and_none_are_omitted(dialect):
    """"Off" is the upstream default the wrapper never touched; emitting a flag would change it."""
    assert render_extra_args({"a": False, "b": None}, dialect) == []


def test_a_list_repeats_an_argparse_flag():
    assert render_extra_args({"chain": ["A", "B"]}, "argparse") == [
        "--chain", "A", "--chain", "B"]


def test_a_list_is_one_hydra_override(record_case):
    """A repeated hydra override keeps only the last value, so repeating one would drop elements silently."""
    actual = render_extra_args({"ppi.hotspot_res": ["A30", "A33"]}, "hydra")
    record_case(input="hydra list", expected=["ppi.hotspot_res=[A30,A33]"], actual=actual)
    assert actual == ["ppi.hotspot_res=[A30,A33]"]


def test_an_unknown_dialect_is_refused():
    with pytest.raises(ValueError, match="dialect"):
        render_extra_args({"a": 1}, "getopt")


def _step_script(script_path, tool_name):
    runtime = pathlib.Path(script_path).parent
    matches = sorted(runtime.glob(f"*_{tool_name}.sh"))
    assert matches, f"no {tool_name} step script under {runtime}"
    return matches[0].read_text(encoding="utf-8", errors="replace")


def _forwarding_tools():
    """{TOOL_NAME: class} for every wrapper that declares the marker, found rather than listed.

    A hand-kept list is what let the previous gate drift: it named ten tools, three of which do build a command line, and none of the six that were actually forwarding. Everything downstream of this function is therefore parametrized over what the tree contains right now.
    """
    import biopipelines  # noqa: F401  (imports every wrapper)
    from biopipelines.base_config import BaseConfig

    found = {}

    def walk(cls):
        for sub in cls.__subclasses__():
            # A class whose __init_subclass__ raised still lingers in __subclasses__, so the throwaway ones defined in this file would otherwise be counted.
            if sub.__module__.startswith("biopipelines.") and sub.FORWARD_UNKNOWN_KWARGS:
                # SolubleMPNN shares ProteinMPNN's TOOL_NAME and its script, so the parent answers for both.
                found.setdefault(sub.TOOL_NAME, sub)
            walk(sub)

    walk(BaseConfig)
    return found


# ── shell safety (developer_manual.md "Shell Safety", layer 4) ────────────────

@pytest.mark.parametrize("value", ['$(whoami)', 'a`id`b', 'say "hi"', "back\\slash"])
def test_an_unsafe_forwarded_value_is_refused(value):
    """The denylist the manual's fourth layer names, applied to a value that reaches bash."""
    with pytest.raises(ValueError):
        _validate_extra_arg("flag", value)


def test_an_unsafe_element_inside_a_list_is_refused():
    with pytest.raises(ValueError):
        _validate_extra_arg("flag", ["ok", "$(id)"])


def test_a_safe_value_passes():
    for value in ["A30,A33", 42, 1.5, True, None, ["A", "B"]]:
        _validate_extra_arg("flag", value)


def test_a_value_that_cannot_become_argv_is_refused():
    """A dict has no command-line spelling; parking it would be the silence this whole change removes."""
    with pytest.raises(ValueError, match="cannot be forwarded"):
        _validate_extra_arg("flag", {"a": 1})


def test_a_non_identifier_key_is_refused():
    with pytest.raises(ValueError, match="identifier"):
        _validate_extra_arg("a;rm -rf /", "x")


@pytest.mark.parametrize("key", ["inference.str_self_cond", "+potentials.substrate"])
def test_a_dotted_hydra_key_is_allowed(key):
    """A hydra override is dotted, and `Tool(**{"denoiser.noise_scale_ca": 1})` is how one reaches a constructor; refusing it would close the escape hatch for the whole RFdiffusion family."""
    _validate_extra_arg(key, True)


@pytest.mark.parametrize("key", ["a b", ".x", "x.", "1x", "a$b"])
def test_a_malformed_key_is_refused(key):
    with pytest.raises(ValueError, match="identifier"):
        _validate_extra_arg(key, True)


def test_a_dotted_key_reaches_a_hydra_tool(local_config, isolated_cwd, new_pipeline):
    from biopipelines.rfdiffusion import RFdiffusion

    pipeline = new_pipeline("fwd_dotted")
    with pipeline:
        RFdiffusion(contigs="100-100", num_designs=1,
                    **{"denoiser.noise_scale_seq": 0.5})
        script_path = pipeline.save()
    assert '"denoiser.noise_scale_seq=0.5"' in _step_script(script_path, "RFdiffusion")


# ── the similarity cutoff ─────────────────────────────────────────────────────

# (unknown key, real parameter names, expected candidate or None) — the pairs the cutoff was tuned against: `omit_AAs` is a genuine ProteinMPNN flag scoring 0.76 against a typed name, so it must forward, while `code`/`codes` at 0.89 must not.
_SIMILARITY_CASES = [
    ("contig", ["contigs"], "contigs"),
    ("code", ["codes"], "codes"),
    ("num_desgins", ["num_designs"], "num_designs"),
    ("symetry", ["symmetry"], "symmetry"),
    ("omit_AAs", ["omit_AA_jsonl", "bias_AA_jsonl"], None),
    ("sampling_temp", ["temperature"], None),
    ("bias_by_res_jsonl", ["bias_AA_jsonl"], None),
    ("model_type", ["model"], None),
]


@pytest.mark.parametrize("key,known,expected", _SIMILARITY_CASES)
def test_similarity_cutoff(key, known, expected, record_case):
    actual = ce.probable_typos([key], known).get(key)
    record_case(input=f"{key!r} against {known}", expected=expected, actual=actual)
    assert actual == expected


def test_sampling_temp_and_temperature_are_not_lexically_close(record_case):
    """They mean the same thing but share almost no letters; no lexical cutoff can pair them, and pretending otherwise would drag real flags in with it."""
    import difflib
    ratio = difflib.SequenceMatcher(a="sampling_temp", b="temperature").ratio()
    record_case(input="sampling_temp vs temperature",
                expected=f"< {ce.TYPO_SIMILARITY_CUTOFF}", actual=round(ratio, 3))
    assert ratio < ce.TYPO_SIMILARITY_CUTOFF


def test_a_framework_key_is_a_typo_candidate_too():
    assert ce.probable_typos(["dependancies"], []) == {"dependancies": "dependencies"}


# ── the checks ────────────────────────────────────────────────────────────────

def test_a_near_miss_names_the_candidate_on_a_non_forwarding_tool(record_case):
    v = ce.check_no_unknown_kwargs("MyTool", {"num_desgins": 8}, known=["num_designs"])
    record_case(input="MyTool(num_desgins=8)", expected="names 'num_designs'", actual=v.message)
    assert "did you mean 'num_designs'?" in v.message


def test_a_near_miss_is_reported_on_a_forwarding_tool_too(record_case):
    v = ce.check_probable_typo("RFdiffusion", {"num_desgins": 8}, known=["num_designs"])
    record_case(input="RFdiffusion(num_desgins=8)",
                expected="probable misspelling", actual=v.message)
    assert v.check == "unknown_kwargs"
    assert "unknown parameter 'num_desgins'; did you mean 'num_designs'?" in v.message


def test_a_non_similar_key_on_a_forwarding_tool_is_informational(record_case):
    v = ce.check_forwarded_kwargs("RFdiffusion", {"denoiser": "x"}, known=["num_designs"])
    record_case(input="RFdiffusion(denoiser='x')",
                expected="forwarding to RFdiffusion", actual=v.message)
    assert v.check == "forwarded_kwargs"
    assert "is forwarding to RFdiffusion: 'denoiser'" in v.message


def test_a_near_miss_is_not_also_announced_as_a_forward():
    """One mistake, one line."""
    assert ce.check_forwarded_kwargs(
        "RFdiffusion", {"num_desgins": 8}, known=["num_designs"]) is None


def test_a_forwarding_tool_never_emits_the_no_effect_warning(capsys):
    ce.reset_reported()
    ce.check_kwargs("RFdiffusion", {"denoiser": "x"}, known=["num_designs"], forwards=True)
    err = capsys.readouterr().err
    assert "never read" not in err and "forwarding to RFdiffusion" in err


# ── refusing a marker that cannot work ────────────────────────────────────────

def test_a_wrapper_that_never_renders_the_tokens_cannot_declare_forwarding(record_case):
    """OpenMM's step script is `python pipe_openmm.py --config <json>`, so the pipe script builds argv and a token the wrapper does not place has nothing to join."""
    from biopipelines.openmm import OpenMM

    assert ce.renders_forwarded_tokens(OpenMM) is False
    with pytest.raises(ce.ContractViolation, match="rendered and then dropped"):
        ce.assert_can_forward(OpenMM)
    record_case(input="assert_can_forward(OpenMM)", expected="raises", actual="raises")


@pytest.mark.parametrize("tool", sorted(_forwarding_tools()))
def test_a_forwarding_tool_passes_the_gate(tool):
    assert ce.assert_can_forward(_forwarding_tools()[tool]) is None


def test_declaring_the_marker_on_an_in_process_tool_fails_at_class_definition():
    """Loud at import, not silently accepted and then dropped at runtime."""
    from biopipelines.openmm import OpenMM

    with pytest.raises(ce.ContractViolation, match="rendered and then dropped"):
        class _Bad(OpenMM):
            FORWARD_UNKNOWN_KWARGS = "argparse"


def test_echoing_the_tokens_is_not_forwarding_them():
    """`extra_args_echo` names the tokens in the step log and nowhere else, so it must not be read as evidence that they reach a command line."""
    from biopipelines.base_config import BaseConfig

    with pytest.raises(ce.ContractViolation, match="rendered and then dropped"):
        class _EchoOnly(BaseConfig):
            TOOL_NAME = "EchoOnly"
            FORWARD_UNKNOWN_KWARGS = "argparse"

            def generate_script(self):
                return self.extra_args_echo()


def test_an_inherited_renderer_counts_for_the_subclass():
    """SolubleMPNN adds no script of its own; answering on its own source alone would refuse a marker that works."""
    from biopipelines.protein_mpnn import SolubleMPNN

    assert ce.renders_forwarded_tokens(SolubleMPNN) is True


def test_an_unknown_dialect_fails_at_class_definition():
    from biopipelines.base_config import BaseConfig

    with pytest.raises(ValueError, match="FORWARD_UNKNOWN_KWARGS"):
        class _Bad(BaseConfig):
            TOOL_NAME = "Whatever"
            FORWARD_UNKNOWN_KWARGS = "getopt"


# ── the tools that forward today keep forwarding ──────────────────────────────

def test_the_converted_tools_still_forward(record_case):
    """A floor, not a census. Losing the marker is a regression and is caught here; gaining one is admitted by `test_every_forwarding_tool_puts_its_token_on_a_command_line`, which makes the new tool prove itself instead of asking this list to be edited. An equality assertion here is exactly the shape that went stale before."""
    found = set(_forwarding_tools())
    converted = {"RFdiffusion", "RFdiffusion2", "RFdiffusion3", "RFdiffusionAllAtom",
                 "ProteinMPNN", "LigandMPNN",
                 "AlphaFold", "LASErMPNN", "Boltz2", "NeuralPLexer",
                 "Aggrescan3D", "CABSflex", "DynamicBind"}
    record_case(input="every BaseConfig subclass", expected=sorted(converted),
                actual=sorted(found))
    # SolubleMPNN subclasses ProteinMPNN and shares its TOOL_NAME, so it forwards by inheritance rather than through a line of its own.
    assert converted <= found


# ── end to end: every marker, discovered, must reach a command line ───────────

# A key no wrapper types and nothing is lexically close to, so it is forwarded rather than reported as a typo, and a value that cannot be confused with anything a step script says on its own.
PROBE_KEY = "zz_forward_probe"
PROBE_VALUE = "BPFWDPROBE"


def _mock(streams, ids=("s1",)):
    from biopipelines.mock import Mock
    return Mock(ids=list(ids), streams=streams, map_table_strategy="config")


def _structures():
    return _mock({"structures": {"format": "pdb", "file": "<id>.pdb"}})


def _protein():
    from biopipelines.sequence import Sequence
    return Sequence(seq="MKTAYIAKQRQISFVKSHFSRQLEERLGL", type="protein", ids="p1")


def _ligand(**kwargs):
    from biopipelines.ligand import Ligand
    return Ligand(**(kwargs or {"code": "HEM"}))


# The minimum each tool needs in order to emit a step script at all. Nothing here is about forwarding; it is the cost of driving a real wrapper end to end, and a tool that gains the marker needs a row so the test below can build it.
FORWARDING_INPUTS = {
    "Aggrescan3D": lambda: dict(structures=_structures().streams.structures),
    "AlphaFold": lambda: dict(proteins=_protein()),
    "Boltz2": lambda: dict(proteins=_protein()),
    "CABSflex": lambda: dict(structures=_structures().streams.structures),
    "Gnina": lambda: dict(structures=_structures(),
                          compounds=_ligand(smiles="CCO", ids="c1", codes="ETH")),
    "DynamicBind": lambda: dict(structures=_structures(),
                                compounds=_ligand(smiles="CCO", ids="c1", codes="ETH")),
    "LASErMPNN": lambda: dict(structures=_structures().streams.structures),
    "Vina": lambda: dict(structures=_structures(),
                         compounds=_ligand(smiles="CCO", ids="c1", codes="ETH")),
    "LigandMPNN": lambda: dict(structures=_structures(), ligand=_ligand(), num_sequences=1),
    "NeuralPLexer": lambda: dict(structures=_structures(),
                                 compounds=_ligand(smiles="CCO", ids="c1", codes="ETH")),
    "ProteinMPNN": lambda: dict(structures=_structures(), num_sequences=1),
    "RFdiffusion": lambda: dict(contigs="100-100", num_designs=1),
    "RFdiffusion2": lambda: dict(pdb=_structures(), ligand=_ligand(), contigs="A1-50",
                                 num_designs=1),
    "RFdiffusion3": lambda: dict(contig="100-100", num_designs=1),
    "RFdiffusionAllAtom": lambda: dict(ligand=_ligand(), contigs="100-200", num_designs=1),
}


def _command_lines(body):
    """The step script with its `echo` lines removed.

    `extra_args_echo` prints the very tokens this test looks for, so a script that announces the forwarded arguments and then drops them would pass a plain substring search over the whole file. That is the failure the marker exists to prevent, so the announcement is not allowed to answer for the command.
    """
    return "\n".join(l for l in body.splitlines() if not l.strip().startswith("echo "))


def test_every_forwarding_tool_has_a_way_to_be_built():
    """The row a new marker needs. Fails the moment someone declares one without giving the test below something to construct, rather than letting the tool go unchecked."""
    missing = sorted(set(_forwarding_tools()) - set(FORWARDING_INPUTS))
    assert not missing, (
        f"{missing} declare FORWARD_UNKNOWN_KWARGS but have no row in "
        f"FORWARDING_INPUTS, so nothing checks that their forwarded tokens reach a "
        f"command line. Add the constructor kwargs each one needs to emit a step script."
    )


def _consumes_extra_args(script_body):
    """Whether the pipe scripts this step script drives read the forwarded tokens back out.

    Follows one level of `import pipe_x` because a tool can be driven by another tool's pipe script -- Vina runs through `pipe_gnina.py`, which delegates to `pipe_vina_backend`. Guessing the script from the tool's name gets that wrong.
    """
    import glob
    import os
    import re

    names = set(re.findall(r"(pipe_[A-Za-z0-9_]+)[.]py", script_body))
    seen, pending = set(), list(names)
    while pending:
        name = pending.pop()
        if name in seen:
            continue
        seen.add(name)
        path = os.path.join(_REPO_ROOT, "pipe_scripts", name + ".py")
        if not os.path.isfile(path):
            continue
        body = open(path, encoding="utf-8", errors="ignore").read()
        if 'config.get("extra_args")' in body or "config.get('extra_args')" in body:
            return True
        pending.extend(re.findall(r"import (pipe_[A-Za-z0-9_]+)", body))
    return False


def _missing_from_config_route(config, tokens, script_body):
    """Tokens that reach neither a command line nor a config JSON some pipe script reads back.

    A wrapper whose bash is `python pipe_<tool>.py --config <json>` has no command line of its own for a token to join, so it carries the tokens in that JSON and the pipe script appends them to the argv list it hands to subprocess.run. No shell is involved there, so nothing re-parses them. Both halves have to hold: a JSON nobody reads drops the tokens exactly as silently as no JSON at all.
    """
    import glob
    import json
    import os

    folder = getattr(config, "output_folder", "") or ""
    carried = []
    for path in glob.glob(os.path.join(folder, "**", "*.json"), recursive=True):
        try:
            with open(path, encoding="utf-8") as handle:
                written = json.load(handle)
        except (OSError, ValueError):
            continue
        if isinstance(written, dict) and written.get("extra_args"):
            carried.append(written["extra_args"])

    if not any(all(token in entry for token in tokens) for entry in carried):
        return list(tokens)
    return [] if _consumes_extra_args(script_body) else list(tokens)


@pytest.mark.parametrize("tool", sorted(_forwarding_tools()))
def test_every_forwarding_tool_puts_its_token_on_a_command_line(
    tool, local_config, isolated_cwd, new_pipeline, record_case,
):
    """The marker's whole claim, checked against what the wrapper actually emits.

    Storing the value on `self.extra_args` is free and proves nothing: the previous gate was a list of tool names that turned out to name none of the tools that forward, so a marker on a wrapper whose bash is `python pipe_<tool>.py --config <json>` was accepted and its tokens dropped. Discovery is the point — a tool that gains the marker is tested by having gained it.
    """
    cls = _forwarding_tools()[tool]
    pipeline = new_pipeline(f"fwd_{tool}")
    with pipeline:
        built = cls(**FORWARDING_INPUTS[tool](), **{PROBE_KEY: PROBE_VALUE})
        script_path = pipeline.save()

    # Inside a Pipeline a wrapper hands back its output, not itself.
    config = getattr(built, "producer", built)
    tokens = config.extra_args_tokens()
    assert tokens, f"{tool} rendered no tokens for {PROBE_KEY}={PROBE_VALUE!r}"

    commands = _command_lines(_step_script(script_path, tool))
    missing = [t for t in tokens if t not in commands]
    if missing:
        # The second legitimate route: a wrapper whose bash is `python pipe_<tool>.py --config <json>`
        # puts the tokens in that JSON, and the pipe script appends them to the argv list it hands to
        # subprocess.run. No shell is involved there, so nothing re-parses them.
        missing = _missing_from_config_route(config, tokens, _step_script(script_path, tool))
    record_case(input=f"{tool}({PROBE_KEY}={PROBE_VALUE!r})",
                expected=f"{tokens} on a command line", actual=f"missing {missing}")
    assert not missing, (
        f"{tool} declares FORWARD_UNKNOWN_KWARGS ({cls.FORWARD_UNKNOWN_KWARGS}) but "
        f"{missing} reach no command line in its step script — only `echo` lines, or "
        f"nothing at all. Interpolate self.extra_args_bash() into the upstream command "
        f"the wrapper emits, or drop the marker."
    )


def test_a_forwarded_kwarg_reaches_the_rfdiffusion_command(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    from biopipelines.rfdiffusion import RFdiffusion

    pipeline = new_pipeline("fwd_rfd")
    with pipeline:
        RFdiffusion(contigs="100-100", num_designs=1, denoiser="x")
        script_path = pipeline.save()
    body = _step_script(script_path, "RFdiffusion")
    record_case(input="RFdiffusion(denoiser='x')", expected='"denoiser=x" in RFD_OPTIONS',
                actual='"denoiser=x"' in body)
    assert '"denoiser=x"' in body
    assert "Forwarding unrecognized argument(s) to RFdiffusion: denoiser=x" in body


@pytest.mark.network
def test_a_forwarded_kwarg_reaches_the_proteinmpnn_command(
    local_config, isolated_cwd, new_pipeline, record_case,
):
    from biopipelines.pdb import PDB
    from biopipelines.protein_mpnn import ProteinMPNN

    pipeline = new_pipeline("fwd_pmpnn")
    with pipeline:
        structures = PDB(pdbs="1CRN")
        ProteinMPNN(structures=structures, num_sequences=1, omit_AAs="CX")
        script_path = pipeline.save()
    body = _step_script(script_path, "ProteinMPNN")
    record_case(input="ProteinMPNN(omit_AAs='CX')", expected='--omit_AAs "CX" on the run line',
                actual='"--omit_AAs" "CX"' in body)
    assert '"--omit_AAs" "CX"' in body
    assert "Forwarding unrecognized argument(s) to ProteinMPNN: --omit_AAs CX" in body


def test_an_unsafe_forwarded_value_fails_at_construction(
    local_config, isolated_cwd, new_pipeline,
):
    """At the user's desk, not as shell injection on a compute node."""
    from biopipelines.rfdiffusion import RFdiffusion

    pipeline = new_pipeline("fwd_unsafe")
    with pytest.raises(ValueError, match="break the generated shell script"):
        with pipeline:
            RFdiffusion(contigs="100-100", num_designs=1, denoiser="$(rm -rf /)")


SHELL_METACHARACTER_VALUES = [
    "CX;id",
    "CX|id",
    "CX && id",
    "CX > /tmp/out",
    "CX < /etc/passwd",
    "CX (id)",
    "CX" + chr(10) + "id",
    "CX" + chr(13) + "id",
]


@pytest.mark.parametrize("value", SHELL_METACHARACTER_VALUES)
def test_a_shell_metacharacter_in_a_forwarded_value_is_refused(value):
    """LigandMPNN emits its command through `eval`, which re-parses the line and strips one layer of quoting, so the double quotes the renderer puts around a forwarded value do not survive to contain a `;`. The freeform denylist covers only the four characters that break double-quoted interpolation, and every character here passes it."""
    with pytest.raises(ValueError, match="shell metacharacter"):
        _validate_extra_arg("omit_AAs", value)


@pytest.mark.parametrize("value", [
    "CX", "A10-A20", "0.75", "chain A", "/path/to/file.pdb", "a,b,c", "key=value", "-1.5e3",
])
def test_a_legitimate_forwarded_value_still_passes(value):
    _validate_extra_arg("omit_AAs", value)


def test_a_non_forwarding_tool_still_parks_and_warns(
    local_config, isolated_cwd, capsys,
):
    from biopipelines.mock import Mock

    Mock(ids=["a"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
         some_upstream_flag=1)
    err = capsys.readouterr().err
    assert "[contract:unknown_kwargs] Mock got unknown constructor parameter(s)" in err
    assert "forwarding" not in err
