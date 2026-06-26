"""Tests for internal tools, the Folder() context manager, and the shared
input resolver (resolve_basic_input) that powers the ligand="LIG" shorthand.

Covers:

* Internal tools (``_internal=True``) are routed under .internal/ in the
  output tree and under RunTime/.internal, Logs/.internal, ToolOutputs/.internal.
* Public step numbering skips internal tools; execution order does not.
* The generated pipeline.sh runs internal and public tools in execution order.
* resolve_basic_input: StandardizedOutput / DataStream / string / error.
* ligand="LIG" creates one internal Ligand + one public consumer; an
  explicitly-constructed Ligand stays public; a malformed code is rejected.
* Folder() nests public outputs, leaves numbering monotonic, reserves
  .internal, validates the name, and never affects internal placement.
"""
from __future__ import annotations

import os
import sys

import pytest

from biopipelines._layout import INTERNAL_FOLDER

PDB_STREAM = {"structures": {"format": "pdb", "file": "<id>.pdb"}}


def _rel(tool, pipeline):
    return os.path.relpath(tool.output_folder, pipeline.folders["output"])


def _by_name(pipeline, name):
    return [t for t in pipeline.tools if t.TOOL_NAME == name]


# ── internal tools: numbering + layout ────────────────────────────────────────

def test_internal_tool_layout_and_numbering(local_config, isolated_cwd, new_pipeline):
    from biopipelines.pipeline import Resources
    from biopipelines.mock import Mock

    p = new_pipeline("internal_layout")
    with p:
        Resources()
        a = Mock(ids=["x"], name="a")
        b = Mock(ids=["y"], name="b", _internal=True)
        c = Mock(ids=["z"], name="c")

    a_t, b_t, c_t = p.tools
    assert (a_t.internal, b_t.internal, c_t.internal) == (False, True, False)
    # Execution order counts everything.
    assert [t.execution_order for t in p.tools] == [1, 2, 3]
    # Public numbering skips the internal tool: a=1, c=2.
    assert a_t.public_step == 1 and c_t.public_step == 2
    assert b_t.public_step is None and b_t.internal_order == 1
    # Output folders.
    assert _rel(a_t, p) == "001_Mock"
    assert _rel(b_t, p) == os.path.join(INTERNAL_FOLDER, "001_Mock")
    assert _rel(c_t, p) == "002_Mock"
    # Script basenames (internal nests under .internal/).
    assert a_t.script_basename == "001_Mock"
    assert b_t.script_basename == os.path.join(INTERNAL_FOLDER, "001_Mock")
    assert c_t.script_basename == "002_Mock"


def test_internal_scripts_logs_nest_under_dot_internal(local_config, isolated_cwd, new_pipeline):
    from biopipelines.pipeline import Resources
    from biopipelines.mock import Mock

    p = new_pipeline("internal_runtime")
    with p:
        Resources()
        Mock(ids=["x"], name="a")
        Mock(ids=["y"], name="b", _internal=True)

    runtime, logs = p.folders["runtime"], p.folders["logs"]
    assert os.path.exists(os.path.join(runtime, "001_Mock.sh"))
    assert os.path.exists(os.path.join(runtime, INTERNAL_FOLDER, "001_Mock.sh"))
    # Logs/.internal dir is created at registration even before execution.
    assert os.path.isdir(os.path.join(logs, INTERNAL_FOLDER))
    tool_outputs = os.path.join(p.folders["output"], "ToolOutputs")
    assert os.path.exists(os.path.join(tool_outputs, INTERNAL_FOLDER, "001_Mock.json"))


def test_pipeline_sh_runs_internal_before_consumer_in_exec_order(
    local_config, isolated_cwd, new_pipeline
):
    from biopipelines.pipeline import Resources
    from biopipelines.mock import Mock

    p = new_pipeline("internal_exec_order")
    with p:
        Resources()
        Mock(ids=["x"], name="a")
        Mock(ids=["y"], name="b", _internal=True)
        Mock(ids=["z"], name="c")

    with open(os.path.join(p.folders["runtime"], "pipeline.sh")) as f:
        text = f.read()
    pos_a = text.index(os.path.join("RunTime", "001_Mock.sh"))
    pos_b = text.index(os.path.join("RunTime", INTERNAL_FOLDER, "001_Mock.sh"))
    pos_c = text.index(os.path.join("RunTime", "002_Mock.sh"))
    assert pos_a < pos_b < pos_c


# ── resolve_basic_input ───────────────────────────────────────────────────────

def test_resolve_basic_input_branches(local_config, isolated_cwd, new_pipeline):
    from biopipelines.input_standardization import resolve_basic_input
    from biopipelines.datastream import DataStream
    from biopipelines.ligand import Ligand

    p = new_pipeline("resolver")
    with p:
        from biopipelines.pipeline import Resources
        Resources()
        # String -> internal Ligand -> compounds stream.
        ds = resolve_basic_input("LIG", Ligand, "compounds", "code")
        assert isinstance(ds, DataStream)
        assert ds.name == "compounds"
        # DataStream passes through unchanged.
        assert resolve_basic_input(ds, Ligand, "compounds", "code") is ds
        # None passes through when allowed; raises otherwise.
        assert resolve_basic_input(None, Ligand, "compounds", "code") is None
        with pytest.raises(ValueError):
            resolve_basic_input(None, Ligand, "compounds", "code", allow_none=False)
        # A non-string, non-stream object raises.
        with pytest.raises(ValueError):
            resolve_basic_input(123, Ligand, "compounds", "code")


# ── ligand="LIG" shorthand on a real consumer ─────────────────────────────────

def test_ligand_string_creates_internal_ligand(local_config, isolated_cwd, new_pipeline):
    from biopipelines.pipeline import Resources
    from biopipelines.mock import Mock
    from biopipelines.ligand_mpnn import LigandMPNN

    p = new_pipeline("ligand_sugar")
    with p:
        Resources()
        s = Mock(ids=["prot"], streams=PDB_STREAM, name="s")
        LigandMPNN(structures=s, ligand="STI", num_sequences=1)

    ligands = _by_name(p, "Ligand")
    assert len(ligands) == 1
    lig = ligands[0]
    assert lig.internal is True
    assert lig.residue_codes == ["STI"]
    # The internal Ligand runs before the public consumer.
    lmpnn = _by_name(p, "LigandMPNN")[0]
    assert lig.execution_order < lmpnn.execution_order
    assert lmpnn.public_step == 2 and lig.public_step is None
    assert lmpnn.ligand_stream is not None and lmpnn.ligand_stream.name == "compounds"


def test_explicit_ligand_stays_public(local_config, isolated_cwd, new_pipeline):
    from biopipelines.pipeline import Resources
    from biopipelines.mock import Mock
    from biopipelines.ligand_mpnn import LigandMPNN
    from biopipelines.ligand import Ligand

    p = new_pipeline("ligand_explicit")
    with p:
        Resources()
        lig = Ligand(code="STI")
        s = Mock(ids=["prot"], streams=PDB_STREAM, name="s")
        LigandMPNN(structures=s, ligand=lig, num_sequences=1)

    assert _by_name(p, "Ligand")[0].internal is False


def test_ligand_string_rejects_unsafe_code(local_config, isolated_cwd, new_pipeline):
    from biopipelines.pipeline import Resources
    from biopipelines.mock import Mock
    from biopipelines.ligand_mpnn import LigandMPNN

    p = new_pipeline("ligand_bad")
    with p:
        Resources()
        s = Mock(ids=["prot"], streams=PDB_STREAM, name="s")
        with pytest.raises(ValueError):
            LigandMPNN(structures=s, ligand="ST I", num_sequences=1)


# ── Folder() ──────────────────────────────────────────────────────────────────

def test_folder_nests_and_keeps_monotonic_numbering(local_config, isolated_cwd, new_pipeline):
    from biopipelines.pipeline import Resources, Folder
    from biopipelines.mock import Mock

    p = new_pipeline("folder_basic")
    with p:
        Resources()
        a = Mock(ids=["x"], name="a")
        with Folder("group"):
            b = Mock(ids=["y"], name="b")
        c = Mock(ids=["z"], name="c")

    a_t, b_t, c_t = p.tools
    assert _rel(a_t, p) == "001_Mock"
    assert _rel(b_t, p) == os.path.join("group", "002_Mock")
    assert _rel(c_t, p) == "003_Mock"


def test_folder_nesting(local_config, isolated_cwd, new_pipeline):
    from biopipelines.pipeline import Resources, Folder
    from biopipelines.mock import Mock

    p = new_pipeline("folder_nested")
    with p:
        Resources()
        with Folder("a"):
            with Folder("b"):
                t = Mock(ids=["x"], name="t")

    assert _rel(p.tools[0], p) == os.path.join("a", "b", "001_Mock")


def test_folder_reserved_and_validated(local_config, isolated_cwd, new_pipeline):
    from biopipelines.pipeline import Resources, Folder
    from biopipelines.mock import Mock

    p = new_pipeline("folder_reserved")
    with p:
        Resources()
        Mock(ids=["x"], name="a")  # keep the pipeline non-empty for save()
        with pytest.raises(ValueError):
            with Folder(INTERNAL_FOLDER):
                pass
        with pytest.raises(ValueError):
            with Folder(".."):
                pass
        with pytest.raises(ValueError):
            with Folder("a/b"):
                pass


def test_folder_does_not_affect_internal_placement(local_config, isolated_cwd, new_pipeline):
    from biopipelines.pipeline import Resources, Folder
    from biopipelines.mock import Mock

    p = new_pipeline("folder_internal")
    with p:
        Resources()
        with Folder("group"):
            Mock(ids=["y"], name="b", _internal=True)

    assert _rel(p.tools[0], p) == os.path.join(INTERNAL_FOLDER, "001_Mock")


@pytest.mark.skipif(
    sys.platform.startswith("win"),
    reason="on-the-fly runs bash on the generated .sh; MSYS/Git-bash mangles "
    "the embedded Windows paths. Runs on POSIX CI. See test_pipeline_generation.",
)
def test_folder_download_zips_on_the_fly(local_config, isolated_cwd):
    import zipfile
    from biopipelines.pipeline import Pipeline, Folder
    from biopipelines.mock import Mock

    p = Pipeline("TestSuite", "folder_dl", on_the_fly=True, local_output=True, config="local")
    with Folder("Results") as results:
        Mock(ids=["x"], streams=PDB_STREAM, name="a")
    archive = results.download()  # not in Colab -> writes the zip, prints path

    assert archive.endswith("Results.zip") and os.path.exists(archive)
    assert any(n.startswith("001_Mock/") for n in zipfile.ZipFile(archive).namelist())


def test_folder_download_rejected_in_submit_mode(local_config, isolated_cwd, new_pipeline):
    from biopipelines.pipeline import Resources, Folder
    from biopipelines.mock import Mock

    p = new_pipeline("folder_dl_submit")
    with p:
        Resources()
        with Folder("Results") as results:
            Mock(ids=["y"], streams=PDB_STREAM, name="b")
        with pytest.raises(RuntimeError):
            results.download()


# ── multi-batch config display labels (internal vs public) ────────────────────

def test_multi_batch_config_display_labels_internal_and_public(
    slurm_local_config, isolated_cwd, new_slurm_pipeline
):
    from pathlib import Path
    from biopipelines.pipeline import Resources
    from biopipelines.mock import Mock
    from biopipelines.ligand_mpnn import LigandMPNN

    p = new_slurm_pipeline("mbatch_labels")
    with p:
        Resources()
        s = Mock(ids=["prot"], streams=PDB_STREAM, name="s")
        Resources()
        LigandMPNN(structures=s, ligand="STI", num_sequences=1)
        p.save()
        p.generate_job_scripts()

    runtime = Path(p.folders["runtime"])
    text = "".join(c.read_text(encoding="utf-8") for c in runtime.glob("config_batch*.sh"))
    # Public tools labeled by public step; the auto-created Ligand is internal.
    assert "Step 001: Mock" in text
    assert "Step 002: LigandMPNN" in text
    assert "Internal 001: Ligand" in text
    # No public tool is ever labeled Step 003 (the internal tool doesn't bump it).
    assert "Step 003" not in text
