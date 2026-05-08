"""Parameter coverage for biopipelines.pool.Pool.

Pool has no user-facing kwargs of its own (the constructor accepts
``*runs: StandardizedOutput``; everything else is internal bookkeeping).
The kwarg-coverage test for Pool therefore reduces to a smoke check that
a 2-run pool emits a pipe_pool.py invocation in the per-step script and
that the gather output folder reaches pipeline.sh.
"""

import pytest

from ._helpers import (
    assert_substrings_in,
    read_pipeline_sh,
)


pytestmark = pytest.mark.tool_parameters


def test_pool_emits_pipe_pool_call(local_config, isolated_cwd, new_pipeline):
    from biopipelines.sequence import Sequence
    from biopipelines.pool import Pool

    pipeline = new_pipeline("pool_params")
    with pipeline:
        a = Sequence(seq="MKTAY", ids="a")
        b = Sequence(seq="AETGF", ids="b")
        Pool(a, b)
        script_path = pipeline.save()

    content = read_pipeline_sh(script_path)
    assert_substrings_in(
        content,
        ["pipe_pool.py", "pool_config.json"],
        label="Pool per-step bash",
    )
