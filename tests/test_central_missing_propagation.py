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


# ── a manifest the footer writes reaches every later step ───────────────────

@pytest.mark.parametrize("by_stream", [False, True], ids=["source", "stream"])
@pytest.mark.parametrize("silent_steps", [1, 2])
def test_a_footer_written_manifest_reaches_the_next_step(sources, by_stream, silent_steps,
                                                          record_case):
    """Panda.filter -> Reduce [-> Reduce] -> Reduce: the last step still excuses the filtered ids.

    The footer writes `missing.csv` for a silent tool but never adds it to that tool's
    `get_output_files()`, so discovery used to stop one step after the filter, and every step
    behind a silent wrapper counted the filtered ids as lost outputs and wrote `_FAILED`.
    """
    from biopipelines.reduce import Reduce

    _generator, upstream = sources
    for _ in range(silent_steps):
        upstream = Reduce(structures=upstream)
    previous = _config(upstream)
    script, declared = _generate(Reduce(structures=upstream.streams.structures if by_stream
                                        else upstream))

    actual = (_blocks(script), "missing" in declared,
              os.path.dirname(previous.table_path("missing")).replace("\\", "/")
              in script.replace("\\", "/"))
    record_case(input=f"Panda.filter -> {silent_steps} x Reduce -> Reduce "
                      f"({'stream' if by_stream else 'source'})",
                expected=(1, True, True), actual=actual)
    assert actual == (1, True, True)


def test_a_tool_owning_its_manifest_reads_a_footer_written_one(sources):
    """Panda.filter -> Reduce -> DSSP: DSSP merges the upstream rows in its own pipe script."""
    from biopipelines.dssp import DSSP
    from biopipelines.reduce import Reduce

    _generator, filtered = sources
    reduced = Reduce(structures=filtered)
    script, declared = _generate(DSSP(structures=reduced))

    expected = _config(reduced).table_path("missing").replace("\\", "/")
    assert f'--upstream-missing "{expected}"' in script.replace("\\", "/")
    assert _blocks(script) == 0 and "missing" in declared


def test_an_upstream_that_cannot_describe_itself_is_named(sources, capsys):
    """Discovery propagates nothing for such a producer, so it must say which one and why."""
    from biopipelines.base_config import BaseConfig
    from biopipelines.reduce import Reduce

    _generator, filtered = sources
    reduced = Reduce(structures=filtered)
    stream = reduced.streams.structures
    upstream = _config(reduced)
    upstream.__dict__.pop("_written_missing_cache", None)

    def broken():
        raise RuntimeError("cannot list outputs")
    # Set by hand, not by monkeypatch: the pipeline's exit still needs the real method
    upstream.get_output_files = broken
    try:
        assert BaseConfig._missing_path_of(stream) is None
    finally:
        del upstream.get_output_files
    out = capsys.readouterr().out
    assert "Warning: Reduce could not describe its outputs" in out and "cannot list outputs" in out


def test_an_unfiltered_chain_still_gets_no_manifest(sources):
    """generator -> Reduce -> Reduce: nothing upstream writes a manifest, so nothing is promised."""
    from biopipelines.reduce import Reduce

    generator, _filtered = sources
    script, declared = _generate(Reduce(structures=Reduce(structures=generator)))

    assert (_blocks(script), "missing" in declared) == (0, False)


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
    # Merges the upstream manifests inside its collector (see test_upstream_missing_owners.py for the tools that merge through `local_missing`).
    "HBDesigner",
])
def test_a_tool_that_writes_its_own_manifest_gets_no_central_block(sources, tool_name):
    """`pipe_propagate_missing.py` rewrites its destination, so a central block would erase these.

    The rows at stake carry `kind="failure"`, which `pipe_check_completion` must NOT excuse. An
    upstream-only merge dropping them would turn a real failure into an expected absence — worse
    than the bug this central path fixes.
    """
    from biopipelines.dssp import DSSP
    from biopipelines.fpocket import FPocket
    from biopipelines.hbdesigner import HBDesigner
    from biopipelines.pdb import PDB
    from biopipelines.prodigy import Prodigy
    from biopipelines.protein_mpnn import ProteinMPNN

    _generator, filtered = sources
    build = {
        "DSSP": lambda: DSSP(structures=filtered),
        "Prodigy": lambda: Prodigy(structures=filtered),
        "FPocket": lambda: FPocket(structures=filtered),
        "ProteinMPNN": lambda: ProteinMPNN(structures=filtered),
        "PDB": lambda: PDB(pdbs=filtered, keep_chains="A"),
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


def test_the_footers_block_runs_in_the_biopipelines_env(sources):
    """The propagator imports pandas, which a tool's own env need not provide."""
    from biopipelines.reduce import Reduce

    _generator, filtered = sources
    built = Reduce(structures=filtered)
    script, _declared = _generate(built)

    footer = script[script.index("BP_MAIN_RC=${BP_MAIN_RC:-$?}"):]
    activation = _config(built).activate_environment(name="biopipelines").strip()
    assert activation and activation in footer
    assert footer.index(activation) < footer.index("pipe_propagate_missing.py")


# ── the filtered map_table excuses what upstream dropped ─────────────────────

def test_the_filtered_map_table_block_names_the_upstream_manifest(sources):
    """`input_sources` was accepted, documented and then dropped on the floor.

    Without `--excused` the materializer raises `ids absent from map_table` for every id the
    upstream step filtered away, whenever the fallback inside `_excused_ids_for` cannot find the
    manifest by itself — a `Load`ed stream, or a tool that writes it off the conventional path.
    """
    from biopipelines.dna_encoder import DNAEncoder

    _generator, filtered = sources
    built = DNAEncoder(sequences=filtered.streams.sequences)
    cfg = _config(built)
    script = cfg.generate_script("unused.sh")

    manifests = cfg._collect_upstream_missing_paths(filtered.streams.sequences)
    assert manifests, "the filtered input carries no upstream manifest to excuse against"

    block = script[script.index("materialize_filtered_map_table.py"):]
    block = block[:block.index("\n")]
    assert "--excused" in block
    for path in manifests:
        assert path.replace("\\", "/") in block.replace("\\", "/")


def test_an_unfiltered_input_adds_no_excused_flag(sources):
    """No manifest upstream, no flag: the materializer must still raise on a genuinely absent id."""
    from biopipelines.dna_encoder import DNAEncoder

    generator, _filtered = sources
    cfg = _config(DNAEncoder(sequences=generator.streams.sequences))
    script = cfg.generate_script("unused.sh")

    block = script[script.index("materialize_filtered_map_table.py"):]
    assert "--excused" not in block[:block.index("\n")]
