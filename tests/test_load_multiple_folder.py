"""LoadMultiple nests its created Load steps under Folder(folder) by default.

Builds two fake tool-output dirs (valid .expected_outputs.json), then loads
them with LoadMultiple inside a pipeline and checks where the Load steps land:
under <output>/LoadMultiple/ by default, at the job root when folder=None.
"""
from __future__ import annotations

import json
import os


def _make_tool_output(parent, folder_name, tool="Mock"):
    """Write a minimal loadable tool-output dir and return its path."""
    p = os.path.join(parent, folder_name)
    os.makedirs(os.path.join(p, "structures"), exist_ok=True)
    open(os.path.join(p, "structures", "x.pdb"), "w").close()
    json.dump(
        {
            "tool_name": tool,
            "tool_class": tool,
            "job_name": "src",
            "output_structure": {
                "output_folder": p,
                "structures": {"x": os.path.join(p, "structures", "x.pdb")},
            },
        },
        open(os.path.join(p, ".expected_outputs.json"), "w"),
    )
    return p


def _step_rels(pipeline):
    """Relative output_folder of each Load *step* registered in the pipeline.

    NB: the StandardizedOutputs returned by LoadMultiple keep pointing at the
    original source data; it is the Load tool steps (pipeline.tools) whose
    output folders reflect the Folder() nesting in the new pipeline.
    """
    return [
        os.path.relpath(t.output_folder, pipeline.folders["output"]).replace(os.sep, "/")
        for t in pipeline.tools
        if t.TOOL_NAME == "Load"
    ]


def test_loadmultiple_nests_load_steps(local_config, isolated_cwd, new_pipeline):
    from biopipelines.pipeline import Resources
    from biopipelines.load import LoadMultiple

    src = os.path.join(isolated_cwd, "src_job")
    _make_tool_output(src, "001_Mock_a")
    _make_tool_output(src, "002_Mock_b")

    p = new_pipeline("lm_nest")
    with p:
        Resources()
        loaded = LoadMultiple(src, tool="Mock", validate_files=False)

    assert set(loaded.keys()) == {"001_Mock_a", "002_Mock_b"}
    rels = _step_rels(p)
    assert len(rels) == 2
    for r in rels:
        assert r.startswith("LoadMultiple/"), r


def test_loadmultiple_folder_none_keeps_root(local_config, isolated_cwd, new_pipeline):
    from biopipelines.pipeline import Resources
    from biopipelines.load import LoadMultiple

    src = os.path.join(isolated_cwd, "src_job")
    _make_tool_output(src, "001_Mock_a")

    p = new_pipeline("lm_root")
    with p:
        Resources()
        LoadMultiple(src, tool="Mock", folder=None, validate_files=False)

    rels = _step_rels(p)
    assert len(rels) == 1
    assert not rels[0].startswith("LoadMultiple/"), rels[0]


def test_loadmultiple_accepts_four_digit_step_prefixes(local_config, isolated_cwd, new_pipeline):
    from biopipelines.pipeline import Resources
    from biopipelines.load import LoadMultiple

    src = os.path.join(isolated_cwd, "src_job")
    _make_tool_output(src, "001_Mock_a")
    _make_tool_output(src, "0999_Mock_b")
    _make_tool_output(src, "1000_Mock_c")
    _make_tool_output(src, "10000_Mock_d")

    p = new_pipeline("lm_4digit")
    with p:
        Resources()
        loaded = LoadMultiple(src, tool="Mock", validate_files=False)

    assert set(loaded.keys()) == {"001_Mock_a", "0999_Mock_b", "1000_Mock_c", "10000_Mock_d"}
