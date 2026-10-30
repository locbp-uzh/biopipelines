# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""MutationComposer: the sequence ids it declares are the ids it writes, for every mode and strategy.

With no `prefix`, the wrapper declared `mut_001`, `rand_001`, `hotspot_001` or `top_001` from the mode alone, while the script chose its own default per mode and per combination strategy (`stack_001`, `rr_001`, `wr_stack_001`, `hs_rr_001`, `top_stack_001`, ...) wherever two or more tables were combined with `stack` or `round_robin`. A pipeline wired against the declared ids then found none of the sequences. The wrapper now resolves the prefix once and hands it to the script, and both build every id through `combinatorics.numbered_id`.
"""

import os
import re
import subprocess
import sys

import pytest

pd = pytest.importorskip("pandas")

from biopipelines import combinatorics

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PIPE = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_mutation_composer.py")

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
REFERENCE = "MKTA"
NUM_SEQUENCES = 3


def _frequency_table(path, bias=0):
    """Four positions, every substitution at 0.05 (the wild type at 0.8), a little more for one residue per table."""
    rows = []
    for position, original in enumerate(REFERENCE, start=1):
        row = {"position": position, "original": original}
        for aa in AMINO_ACIDS:
            row[aa] = 0.8 if aa == original else 0.05
        row[AMINO_ACIDS[(position + bias) % 20]] += 0.01 * (bias + 1)
        rows.append(row)
    pd.DataFrame(rows).to_csv(path, index=False)
    return str(path)


def _compose(new_pipeline, tmp_path, mode, strategy, tables, prefix=""):
    from biopipelines.mutation_composer import MutationComposer

    paths = [_frequency_table(tmp_path / f"freq{i}.csv", bias=i) for i in range(tables)]
    pipeline = new_pipeline(f"composer_{mode}_{strategy}_{tables}")
    with pipeline:
        composer = MutationComposer(
            frequencies=paths if tables > 1 else paths[0],
            num_sequences=NUM_SEQUENCES, mode=mode, combination_strategy=strategy,
            hotspot_count=3, random_seed=0, prefix=prefix)
        pipeline.save()
    return composer.producer, list(composer.streams.sequences.ids)


def _written_ids(config):
    done = subprocess.run([sys.executable, PIPE, "--config", config.config_file], capture_output=True, text=True)
    assert done.returncode == 0, done.stdout + done.stderr
    return list(pd.read_csv(config.sequences_csv)["id"])


@pytest.mark.parametrize("mode", ["single_point", "weighted_random", "hotspot_focused", "top_mutations"])
@pytest.mark.parametrize("strategy, tables", [
    ("average", 1), ("average", 2), ("maximum", 2), ("stack", 1), ("stack", 2), ("round_robin", 2),
])
def test_the_written_ids_are_the_declared_ones(local_config, isolated_cwd, new_pipeline, tmp_path, mode, strategy, tables):
    config, declared = _compose(new_pipeline, tmp_path, mode, strategy, tables)

    written = _written_ids(config)

    assert len(written) == NUM_SEQUENCES
    assert written == declared[:len(written)]


@pytest.mark.parametrize("strategy, tables", [("average", 1), ("stack", 2), ("round_robin", 2)])
def test_a_given_prefix_names_every_sequence(local_config, isolated_cwd, new_pipeline, tmp_path, strategy, tables):
    config, declared = _compose(new_pipeline, tmp_path, "single_point", strategy, tables, prefix="scan")

    assert declared == ["scan_001", "scan_002", "scan_003"]
    assert _written_ids(config) == declared


def test_the_default_prefix_follows_the_mode(local_config, isolated_cwd, new_pipeline, tmp_path):
    _, declared = _compose(new_pipeline, tmp_path, "weighted_random", "stack", 2)

    assert declared == ["rand_001", "rand_002", "rand_003"]


# ── numbered_id, the one place the id is spelled ──────────────────────────────

def test_numbered_id_pads_the_number_to_three_digits():
    assert combinatorics.numbered_id("mut", 7) == "mut_007"
    assert combinatorics.numbered_id("mut", 1234) == "mut_1234"


def test_numbered_id_refuses_an_empty_prefix():
    with pytest.raises(ValueError, match="prefix"):
        combinatorics.numbered_id("", 1)


def test_the_script_no_longer_chooses_a_prefix_of_its_own():
    with open(PIPE, encoding="utf-8") as f:
        code = f.read()
    assert not re.search(r"if prefix else", code), "the wrapper resolves the prefix once"
    assert not re.search(r"f[\"']\{prefix\}_", code), "ids come from combinatorics.numbered_id"
