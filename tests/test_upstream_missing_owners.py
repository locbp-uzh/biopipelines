# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""A tool that writes its own `missing.csv` still carries the rows an upstream filter wrote.

Declaring a `missing` table makes the footer's central propagation stand down, because `pipe_propagate_missing.py` rewrites its destination and would erase the tool's own failure rows. The tools that then wrote only their own failures read no upstream manifest, so a parent dropped by `Panda.filter` left its declared outputs unexcused and the step was marked FAILED with every surviving output on disk (HBDesigner showed it on the cluster).

The framework's standard is the `local_missing` idiom: the pipe writes its own rows to `_execution/local_missing.csv` and `generate_missing_propagation(..., local_missing=...)` merges them with every upstream manifest into `tables/missing.csv`, so neither side is erased. Each wrapper here must either do that or hand its pipe the upstream manifests (`--upstream-missing`).
"""

import json
import os
import re

import pytest

pytest.importorskip("yaml")

from biopipelines.datastream import DataStream, create_map_table


@pytest.fixture
def filtered(local_config, isolated_cwd):
    """A structures source and a sequences source, each with one id dropped by an upstream filter."""
    from biopipelines.mock import Mock
    from biopipelines.panda import Panda
    from biopipelines.pipeline import Pipeline

    with Pipeline(project="TestSuite", job="owners", description="x", on_the_fly=False,
                  local_output=True, config="local"):
        structures = Mock(
            ids=["a", "b", "c"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            tables={"scores": {"columns": ["score"], "fill": {"score": 1.0}}},
        )
        sequences = Mock(
            ids=["a", "b", "c"],
            streams={"sequences": {"format": "csv", "values": ["AAA", "BB", "CCC"]}},
            tables={"scores": {"columns": ["score"], "fill": {"score": 1.0}}},
        )
        yield {
            "structures": Panda(pool=structures, operations=[Panda.filter("score > 0.5")]),
            "sequences": Panda(pool=sequences, operations=[Panda.filter("score > 0.5")]),
        }


def _compounds(tmp_path):
    table = tmp_path / "compounds_map.csv"
    create_map_table(str(table), ids=["l1"], additional_columns={"code": ["LIG"], "smiles": ["CCO"]})
    return DataStream(name="compounds", ids=["l1"], files=[], map_table=str(table), format="csv")


def _build(name, filtered, tmp_path):
    from biopipelines.diffdock import DiffDock
    from biopipelines.dynamicbind import DynamicBind
    from biopipelines.frame2seq import Frame2Seq
    from biopipelines.lasermpnn import LASErMPNN
    from biopipelines.neuralplexer import NeuralPLexer
    from biopipelines.placer import PLACER
    from biopipelines.plm_sol import PLM_Sol
    from biopipelines.pocketgen import PocketGen
    from biopipelines.thermompnn import ThermoMPNN
    from biopipelines.vespag import VespaG

    structures, sequences = filtered["structures"], filtered["sequences"]
    ligands = _compounds(tmp_path)
    return {
        "DiffDock": lambda: DiffDock(structures=structures, compounds=ligands),
        "DynamicBind": lambda: DynamicBind(structures=structures, compounds=ligands),
        "NeuralPLexer": lambda: NeuralPLexer(structures=structures, compounds=ligands),
        "PocketGen": lambda: PocketGen(structures=structures, ligand=ligands),
        "PLACER": lambda: PLACER(structures=structures, ligand=ligands),
        "Frame2Seq": lambda: Frame2Seq(structures=structures),
        "LASErMPNN": lambda: LASErMPNN(structures=structures),
        "ThermoMPNN": lambda: ThermoMPNN(structures=structures),
        "PLM_Sol": lambda: PLM_Sol(sequences=sequences),
        "VespaG": lambda: VespaG(sequences=sequences),
    }[name]()


OWNERS = ["DiffDock", "DynamicBind", "NeuralPLexer", "PocketGen", "PLACER",
          "Frame2Seq", "LASErMPNN", "ThermoMPNN", "PLM_Sol", "VespaG"]


def _script_and_tables(built):
    tool = built.tool if hasattr(built, "tool") else built
    script = tool.generate_script("unused.sh")
    with open(os.path.join(tool.output_folder, ".expected_outputs.json"), encoding="utf-8") as handle:
        declared = json.load(handle)["output_structure"].get("tables", {})
    return tool, script, declared


@pytest.mark.parametrize("name", OWNERS)
def test_the_upstream_manifest_reaches_the_tools_own_missing_table(name, filtered, tmp_path):
    tool, script, declared = _script_and_tables(_build(name, filtered, tmp_path))
    canonical = tool.table_path("missing").replace("\\", "/")
    flat = script.replace("\\", "/")

    assert "missing" in declared, f"{name} must still declare the table it writes"
    handed_to_pipe = "--upstream-missing" in flat
    merged = "pipe_propagate_missing.py" in flat and "--local-missing" in flat
    assert handed_to_pipe or merged, f"{name} reads no upstream `missing` manifest"
    if merged:
        # The pipe must not write the canonical file itself: the merge rewrites it and would erase its rows.
        commands = re.split(r"\n\s*\n", flat)
        writers = [c for c in commands if f'"{canonical}"' in c]
        assert writers and all("pipe_propagate_missing.py" in c for c in writers), (
            f"{name}: a command other than the merge writes {canonical}")
        assert flat.count("BP_MAIN_RC=$?") == 1


@pytest.mark.parametrize("name", OWNERS)
def test_the_merge_runs_after_the_pipe_and_before_the_footer(name, filtered, tmp_path):
    _tool, script, _declared = _script_and_tables(_build(name, filtered, tmp_path))
    flat = script.replace("\\", "/")
    if "pipe_propagate_missing.py" not in flat:
        pytest.skip("hands the manifest to its pipe instead")
    assert flat.index("BP_MAIN_RC=$?") < flat.index("pipe_propagate_missing.py") < flat.index("Checking outputs")
