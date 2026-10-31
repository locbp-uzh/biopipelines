# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""MutationComposer: the two ways the step fails are the two the documentation names.

The docstring and `docs/tool/sequence_design.md` said the step "fails only when it wrote no sequence at all". The pipe also raises when it wrote sequences but every one of them is the reference (a mode that finds no mutation at or above `min_frequency` writes the reference under the description `none`), which points at the frequency data or the threshold, not at a short table. Both are documented now; the behavior is unchanged.
"""

import os
import re

import pandas as pd

from tests.test_mutation_composer_missing import NUM_SEQUENCES, _configured_step, _run

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_a_run_whose_every_sequence_is_the_reference_fails_the_step_but_writes_them(local_config, isolated_cwd,
                                                                                    new_pipeline):
    config_path, config, script, step_id = _configured_step(
        new_pipeline, "top_mutations", "average", 1, False, min_frequency=0.9)

    done = _run(config_path)

    assert done.returncode != 0
    assert "Failed to generate any mutations" in done.stdout + done.stderr
    written = pd.read_csv(config["sequences_output"], dtype=str, keep_default_na=False)
    assert len(written) == NUM_SEQUENCES and set(written["mutations"]) == {"none"}
    assert len(pd.read_csv(config["local_missing_csv"])) == 0


def _statements():
    from biopipelines.mutation_composer import MutationComposer

    docs = open(os.path.join(REPO_ROOT, "docs", "tool", "sequence_design.md"), encoding="utf-8").read()
    docs_line = next(line for line in docs.splitlines() if "The step declares `num_sequences` ids" in line)
    doc = MutationComposer.__init__.__doc__
    docstring = next(line for line in doc.splitlines() if line.strip().startswith("missing: id | removed_by"))
    return {"docs": docs_line, "docstring": docstring}


def test_both_ways_the_step_fails_are_documented():
    for where, text in _statements().items():
        assert "no sequence at all" in text, where
        assert re.search(r"every sequence (it wrote )?is (identical to )?the reference", text), (where, text)
        assert "fails only when" not in text, where
