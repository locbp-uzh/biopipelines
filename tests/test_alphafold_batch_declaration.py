"""A fold step fed a hand-narrowed batch is judged on the ids it was given, not on its producer's full list.

A finished ProteinMPNN step can be cut into batches by copying its folder and narrowing the
tables (`sequences.csv`, `designs.csv`) while the manifest stays whole. `sequences` and
`designs` are value-based streams (`files == []`), which `Load` used to propagate with the
producer's full id list. `AlphaFold` declares `structures`, `pae` and `msas` from those ids, so
it declared every design of the producer, folded the few rows of the batch, and the completion
check marked a step that did exactly what it was asked as FAILED.
"""

import csv
import json
import os
import sys

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "pipe_scripts"))

import pipe_check_completion  # noqa: E402
from biopipelines import id_patterns  # noqa: E402

FULL = [f"design_{i}" for i in range(1, 11)]
BATCH = FULL[:4]


def _write_csv(path, columns, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(columns)
        w.writerows(rows)


@pytest.fixture
def batch_folder(tmp_path):
    """A ProteinMPNN-shaped folder narrowed to BATCH, its manifest still naming the 10 designs and the source folder."""
    folder = str(tmp_path / "batch" / "003_ProteinMPNN")
    _write_csv(os.path.join(folder, "sequences", "sequences.csv"), ["id", "sequence"],
               [[i, "MKTAY"] for i in BATCH])
    _write_csv(os.path.join(folder, "designs", "designs.csv"), ["id", "source"],
               [[i, "backbone"] for i in BATCH])
    fasta = os.path.join(folder, "sequences", "sequences.fasta")
    with open(fasta, "w") as f:
        f.writelines(f">{i}\nMKTAY\n" for i in BATCH)
    missing = os.path.join(folder, "tables", "missing.csv")
    _write_csv(missing, ["id", "removed_by", "kind", "cause"], [])

    def stream(name, files, map_table, fmt):
        return {"name": name, "ids": ["design_<1..10>"], "files": files,
                "map_table": map_table, "format": fmt, "metadata": {}}

    sequences_csv = os.path.join(folder, "sequences", "sequences.csv")
    manifest = {
        "tool_name": "ProteinMPNN",
        "tool_class": "ProteinMPNN",
        "output_structure": {
            "sequences": stream("sequences", [], sequences_csv, "csv"),
            "designs": stream("designs", [], os.path.join(folder, "designs", "designs.csv"), "csv"),
            "fasta": stream("fasta", fasta, sequences_csv, "fasta"),
            "tables": {"missing": {"name": "missing", "path": missing,
                                   "columns": ["id", "removed_by", "kind", "cause"],
                                   "description": ""}},
            "output_folder": str(tmp_path / "source" / "003_ProteinMPNN"),
        },
    }
    with open(os.path.join(folder, ".expected_outputs.json"), "w") as f:
        json.dump(manifest, f)
    return folder


@pytest.fixture
def fold(local_config, isolated_cwd, batch_folder):
    """The AlphaFold step over the batch, with the manifest its bash would pass to the completion check."""
    from biopipelines.alphafold import AlphaFold
    from biopipelines.load import Load
    from biopipelines.pipeline import Pipeline

    with Pipeline(project="TestSuite", job="af_batch", description="x",
                  on_the_fly=False, local_output=True, config="local"):
        batch = Load(batch_folder)
        built = AlphaFold(proteins=batch)
        tool = built.tool if hasattr(built, "tool") else built
        tool.generate_script("unused.sh")
        with open(os.path.join(tool.output_folder, ".expected_outputs.json")) as f:
            declared = json.load(f)["output_structure"]
        yield tool, declared


def _produce(declared, ids):
    """Write what a finished fold step leaves: one file per id and stream, and every declared table."""
    for name in ("structures", "pae", "msas"):
        template = declared[name]["files"][0]
        for item_id in ids:
            path = template.replace("<id>", item_id)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            open(path, "w").write("x\n")
    for info in declared["tables"].values():
        _write_csv(info["path"], ["id"], [[i] for i in ids])


def _check(tool, declared):
    return pipe_check_completion.check_expected_outputs(declared, "AlphaFold", tool.output_folder)


def test_a_narrowed_batch_declares_its_own_ids_per_stream(fold, record_case):
    """Four designs in the batch declare four structures, four pae and four msas, not ten."""
    _tool, declared = fold
    actual = {name: len(id_patterns.expand_ids(declared[name]["ids"]))
              for name in ("structures", "pae", "msas")}
    record_case(input="AlphaFold(Load(batch of 4 of 10 designs))",
                expected={"structures": 4, "pae": 4, "msas": 4}, actual=actual)
    for name in ("structures", "pae", "msas"):
        assert id_patterns.expand_ids(declared[name]["ids"]) == BATCH


def test_a_step_that_folded_its_batch_completes(fold, record_case):
    tool, declared = fold
    _produce(declared, BATCH)

    found, total = pipe_check_completion.tally_declared_outputs(declared)
    ok, missing = _check(tool, declared)

    record_case(input="4 pdb + 4 pae + 4 a3m + tables written", expected=(True, {}, total),
                actual=(ok, missing, found))
    assert (ok, missing) == (True, {})
    assert found == total


def test_a_batch_structure_that_was_not_produced_still_fails_the_step(fold):
    """Narrowing the declaration must not excuse a fold that really went missing."""
    tool, declared = fold
    _produce(declared, BATCH)
    os.remove(declared["structures"]["files"][0].replace("<id>", BATCH[2]))

    ok, missing = _check(tool, declared)

    assert not ok
    assert list(missing) == ["structures"]
    assert [os.path.basename(path) for path in missing["structures"]] == [f"{BATCH[2]}.pdb"]
