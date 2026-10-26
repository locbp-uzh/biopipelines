"""Upstream `missing`-manifest propagation is the footer's job, not each wrapper's.

A tool that consumes ids declares, at configuration time, the output files it expects. When an
upstream step legitimately filtered some ids away, those ids produce no output here — correctly —
and `pipe_check_completion` excuses the absence only if this tool owns a `tables/missing.csv`
naming them. While that was opt-in per wrapper, 38 of 89 were silent and reported legitimate
absences as lost outputs, wrote a `_FAILED` marker and stranded every queued job behind them on
`DependencyNeverSatisfied`. Six separate commits fixed one instance each; the last was
Aggrescan3D, where 118 of 2,799 designs rejected by a geometric gate failed a step whose 2,682
real outputs were all on disk.

What must hold, and what each test here pins:
  - a tool handed a filtered input gets the block and the declaration with no code of its own,
    whether it kept the input source or only the stream;
  - an id generator gets neither, so no tool promises a manifest nothing writes;
  - a tool that owns its manifest keeps it, because the propagation helper rewrites its
    destination and would erase rows the tool recorded itself;
  - the step's real exit status survives, which is why the block the footer emits must not
    re-capture `$?` after the footer already saved it.
"""

import json
import os

import pytest

pytest.importorskip("yaml")


def _config(built):
    """The BaseConfig behind whatever a wrapper returned."""
    return built.tool if hasattr(built, "tool") else built


def _generate(built):
    """(script, tables declared in .expected_outputs.json) for one constructed tool."""
    cfg = _config(built)
    script = cfg.generate_script("unused.sh")
    manifest = os.path.join(cfg.output_folder, ".expected_outputs.json")
    with open(manifest, encoding="utf-8") as fh:
        declared = json.load(fh)["output_structure"].get("tables", {})
    return script, declared


def _blocks(script):
    return script.count("pipe_propagate_missing.py")


@pytest.fixture
def sources(local_config, isolated_cwd):
    """A filtering step and the unfiltered generator feeding it, inside a live pipeline."""
    from biopipelines.mock import Mock
    from biopipelines.panda import Panda
    from biopipelines.pipeline import Pipeline

    with Pipeline(project="TestSuite", job="prop", description="x", on_the_fly=False,
                  local_output=True, config="local"):
        generator = Mock(
            ids=["a", "b", "c"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"},
                     "sequences": {"format": "csv", "values": ["AAA", "BB", "CCC"]}},
            tables={"scores": {"columns": ["score"], "fill": {"score": 1.0}}},
        )
        yield generator, Panda(pool=generator, operations=[Panda.filter("score > 0.5")])


# ── a silent wrapper is covered with no code of its own ──────────────────────

def test_a_wrapper_that_wires_nothing_still_propagates(sources, record_case):
    """Reduce declares no `missing` table and emits no block of its own."""
    from biopipelines.reduce import Reduce

    _generator, filtered = sources
    script, declared = _generate(Reduce(structures=filtered))

    actual = (_blocks(script), "missing" in declared)
    record_case(input="Reduce(structures=<filtered>)", expected=(1, True), actual=actual)
    assert actual == (1, True)


def test_the_declared_manifest_is_this_tools_own_file(sources):
    """`pipe_check_completion` reads tables.missing.path, and the block writes that path."""
    from biopipelines.reduce import Reduce

    _generator, filtered = sources
    built = Reduce(structures=filtered)
    script, declared = _generate(built)

    path = declared["missing"]["path"] if isinstance(declared["missing"], dict) \
        else declared["missing"].path
    assert path == _config(built).table_path("missing")
    assert path.replace("\\", "/") in script.replace("\\", "/")


def test_a_wrapper_keeping_only_the_stream_is_covered(sources):
    """A `DataStream` carries no `.tables`; it resolves through the producer that stamped it.

    This is the shape Aggrescan3D had to be edited for — it kept `structures.streams.structures`
    and dropped the source, so nothing downstream could find the manifest.
    """
    from biopipelines.reduce import Reduce

    _generator, filtered = sources
    script, declared = _generate(Reduce(structures=filtered.streams.structures))

    assert (_blocks(script), "missing" in declared) == (1, True)


def test_discovery_looks_inside_bundle_and_each(sources):
    from biopipelines.combinatorics import Bundle
    from biopipelines.esmfold2 import ESMFold2

    _generator, filtered = sources
    script, declared = _generate(ESMFold2(proteins=Bundle(filtered.streams.sequences)))

    assert (_blocks(script), "missing" in declared) == (1, True)


# ── an id generator must be left alone ───────────────────────────────────────

@pytest.mark.parametrize("by_stream", [False, True], ids=["source", "stream"])
def test_an_unfiltered_input_gets_no_manifest(sources, by_stream, record_case):
    """Every stream carries a canonical `missing_table` stamp whether or not the file is written.

    Trusting the stamp alone would make an id generator look like a filter, and the tool would
    then promise a table nothing ever writes.
    """
    from biopipelines.reduce import Reduce

    generator, _filtered = sources
    structures = generator.streams.structures if by_stream else generator
    script, declared = _generate(Reduce(structures=structures))

    actual = (_blocks(script), "missing" in declared)
    record_case(input=f"Reduce(structures=<generator {'stream' if by_stream else 'source'}>)",
                expected=(0, False), actual=actual)
    assert actual == (0, False)


# ── a tool that owns its manifest keeps it ───────────────────────────────────

@pytest.mark.parametrize("tool_name", [
    # Merge upstream rows with their own failures inside their pipe script.
    "DSSP", "Prodigy", "FPocket", "ProteinMPNN",
    # Re-keys each upstream row into its own id space, so a rename via `ids=` still excuses.
    "PDB",
    # Record only their OWN failures and read no upstream manifest, so they signal nothing
    # about which helper they use -- the declaration is the only thing that marks them.
    "ThermoMPNN", "PLM_Sol", "Frame2Seq", "VespaG", "LASErMPNN", "HBDesigner",
])
def test_a_tool_that_writes_its_own_manifest_gets_no_central_block(sources, tool_name):
    """`pipe_propagate_missing.py` rewrites its destination, so a central block would erase these.

    The rows at stake carry `kind="failure"`, which `pipe_check_completion` must NOT excuse. An
    upstream-only merge dropping them would turn a real failure into an expected absence — worse
    than the bug this central path fixes.
    """
    from biopipelines.dssp import DSSP
    from biopipelines.fpocket import FPocket
    from biopipelines.frame2seq import Frame2Seq
    from biopipelines.hbdesigner import HBDesigner
    from biopipelines.lasermpnn import LASErMPNN
    from biopipelines.pdb import PDB
    from biopipelines.plm_sol import PLM_Sol
    from biopipelines.prodigy import Prodigy
    from biopipelines.protein_mpnn import ProteinMPNN
    from biopipelines.thermompnn import ThermoMPNN
    from biopipelines.vespag import VespaG

    _generator, filtered = sources
    build = {
        "DSSP": lambda: DSSP(structures=filtered),
        "Prodigy": lambda: Prodigy(structures=filtered),
        "FPocket": lambda: FPocket(structures=filtered),
        "ProteinMPNN": lambda: ProteinMPNN(structures=filtered),
        "PDB": lambda: PDB(pdbs=filtered, keep_chains="A"),
        "ThermoMPNN": lambda: ThermoMPNN(structures=filtered),
        "PLM_Sol": lambda: PLM_Sol(sequences=filtered),
        "Frame2Seq": lambda: Frame2Seq(structures=filtered),
        "VespaG": lambda: VespaG(sequences=filtered),
        "LASErMPNN": lambda: LASErMPNN(structures=filtered),
        "HBDesigner": lambda: HBDesigner(structures=filtered),
    }[tool_name]
    script, declared = _generate(build())

    assert _blocks(script) == 0, f"{tool_name} owns its manifest; the central block would erase it"
    assert "missing" in declared, f"{tool_name} must still declare the table it writes"


def test_a_tool_emitting_its_own_block_does_not_get_a_second(sources):
    from biopipelines.aggrescan3d import Aggrescan3D

    _generator, filtered = sources
    script, declared = _generate(Aggrescan3D(structures=filtered))

    assert _blocks(script) == 1
    assert "missing" in declared


def test_generating_twice_gives_the_same_script(sources):
    """The ownership flag is per-generation state; a stale one would drop or double the block."""
    from biopipelines.aggrescan3d import Aggrescan3D
    from biopipelines.reduce import Reduce

    _generator, filtered = sources
    for built in (Reduce(structures=filtered), Aggrescan3D(structures=filtered)):
        cfg = _config(built)
        first = cfg.generate_script("unused.sh")
        assert cfg.generate_script("unused.sh") == first
        assert cfg.generate_script("unused.sh") == first
        assert _blocks(first) == 1


# ── the step's own exit status survives the block ────────────────────────────

def test_the_footers_block_runs_after_the_status_is_captured(sources):
    """The footer saves `$?` first, so the block cannot turn a failed step into a clean one."""
    from biopipelines.reduce import Reduce

    _generator, filtered = sources
    script, _declared = _generate(Reduce(structures=filtered))

    capture = script.index("BP_MAIN_RC=${BP_MAIN_RC:-$?}")
    assert capture < script.index("pipe_propagate_missing.py")
    # And it captures once: a second would report the propagation's own success.
    assert script.count("BP_MAIN_RC=$?") == 0
