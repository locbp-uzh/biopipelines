"""Parameter coverage for biopipelines.pymol.PyMOL.

PyMOL is a declarative op-API, not a CLI mirror. The wrapper serialises
the list of operations into a JSON config consumed at runtime by
pipe_pymol.py. We assert each op-type and its salient parameters land in
that JSON.
"""

import pytest

from ._helpers import assert_substrings_in, read_all_emitted_artifacts


pytestmark = pytest.mark.tool_parameters


def _run_pipeline(local_config, isolated_cwd, new_pipeline, *ops, session="my_session"):
    from biopipelines.mock import Mock
    from biopipelines.pymol import PyMOL

    pipeline = new_pipeline("pymol_params")
    with pipeline:
        Mock(
            ids=["s1", "s2"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            map_table_strategy="config",
        )
        # Resolve `<structures>` placeholders inside ops by building real ones
        # against the Mock above. Caller passes ops that reference the Mock's
        # streams.structures by closure if needed.
        PyMOL(*ops, session=session)
        script_path = pipeline.save()
    return read_all_emitted_artifacts(script_path)


def test_load_op(local_config, isolated_cwd, new_pipeline):
    from biopipelines.mock import Mock
    from biopipelines.pymol import PyMOL

    pipeline = new_pipeline("pymol_params")
    with pipeline:
        m = Mock(
            ids=["s1"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            map_table_strategy="config",
        )
        PyMOL(PyMOL.Load(structures=m.streams.structures))
        script_path = pipeline.save()
    content = read_all_emitted_artifacts(script_path)
    assert '"op": "load"' in content


def test_color_op(local_config, isolated_cwd, new_pipeline):
    from biopipelines.mock import Mock
    from biopipelines.pymol import PyMOL

    pipeline = new_pipeline("pymol_params")
    with pipeline:
        m = Mock(
            ids=["s1"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            map_table_strategy="config",
        )
        PyMOL(
            PyMOL.Load(structures=m.streams.structures),
            PyMOL.Color(color="cyan", structures=m.streams.structures),
        )
        script_path = pipeline.save()
    content = read_all_emitted_artifacts(script_path)
    assert '"op": "color"' in content
    assert '"color": "cyan"' in content


def test_show_op(local_config, isolated_cwd, new_pipeline):
    from biopipelines.mock import Mock
    from biopipelines.pymol import PyMOL

    pipeline = new_pipeline("pymol_params")
    with pipeline:
        m = Mock(
            ids=["s1"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            map_table_strategy="config",
        )
        PyMOL(
            PyMOL.Load(structures=m.streams.structures),
            PyMOL.Show(representation="cartoon"),
        )
        script_path = pipeline.save()
    content = read_all_emitted_artifacts(script_path)
    assert '"op": "show"' in content


def test_align_op(local_config, isolated_cwd, new_pipeline):
    from biopipelines.mock import Mock
    from biopipelines.pymol import PyMOL

    pipeline = new_pipeline("pymol_params")
    with pipeline:
        m = Mock(
            ids=["s1", "s2"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            map_table_strategy="config",
        )
        PyMOL(
            PyMOL.Load(structures=m.streams.structures),
            PyMOL.Align(method="cealign"),
        )
        script_path = pipeline.save()
    content = read_all_emitted_artifacts(script_path)
    assert '"op": "align"' in content
    assert '"method": "cealign"' in content


def test_session_name(local_config, isolated_cwd, new_pipeline):
    from biopipelines.mock import Mock
    from biopipelines.pymol import PyMOL

    pipeline = new_pipeline("pymol_params")
    with pipeline:
        m = Mock(
            ids=["s1"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            map_table_strategy="config",
        )
        PyMOL(PyMOL.Load(structures=m.streams.structures), session="custom_view")
        script_path = pipeline.save()
    content = read_all_emitted_artifacts(script_path)
    assert '"session_name": "custom_view"' in content


def test_smoke_multi_op(local_config, isolated_cwd, new_pipeline):
    from biopipelines.mock import Mock
    from biopipelines.pymol import PyMOL

    pipeline = new_pipeline("pymol_params")
    with pipeline:
        m = Mock(
            ids=["s1", "s2"],
            streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
            map_table_strategy="config",
        )
        PyMOL(
            PyMOL.Load(structures=m.streams.structures),
            PyMOL.Color(color="magenta", structures=m.streams.structures),
            PyMOL.Show(representation="cartoon"),
            PyMOL.Align(method="align"),
            session="multi",
        )
        script_path = pipeline.save()
    content = read_all_emitted_artifacts(script_path)
    assert_substrings_in(content, [
        '"op": "load"',
        '"op": "color"',
        '"op": "show"',
        '"op": "align"',
        '"session_name": "multi"',
        '"color": "magenta"',
    ])
