# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""MSA declares every column its pipe can write.

The `msas` table declared three columns while `pipe_msa.py` builds the frame with four and lets
`prune_redundant_provenance_columns` drop `sequences.id` only when the row id already spells it. An input
whose id does not spell its source therefore wrote a column absent from the declaration, so `.tables`
under-reported the schema and a consumer projecting onto the declared columns lost the provenance.
"""

import json
import os
import subprocess
import sys

import pandas as pd

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PIPE = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_msa.py")
PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"


def _a3m(folder, name):
    path = folder / f"{name}.a3m"
    path.write_text(f"#{len(PROTEIN)}\t1\n>101\n{PROTEIN}\n>UniRef100_X\n{PROTEIN}\n")
    return str(path)


def _convert(tmp_path, rows):
    """Run the converter over `rows` and return the frame it wrote."""
    step = tmp_path / "007_MSA"
    table = tmp_path / "in.csv"
    pd.DataFrame(rows).to_csv(table, index=False)
    config = tmp_path / "config.json"
    config.write_text(json.dumps({
        "input_msa_table": str(table), "convert": "csv", "output_folder": str(step / "msas"),
        "output_msas_csv": str(step / "msas" / "msas_map.csv"),
        "local_missing_csv": str(step / "_execution" / "local_missing.csv")}))
    done = subprocess.run([sys.executable, PIPE, "--config", str(config)], capture_output=True, text=True)
    assert done.returncode == 0, done.stdout + done.stderr
    return pd.read_csv(step / "msas" / "msas_map.csv")


def _declared_columns(isolated_cwd, new_pipeline):
    from biopipelines.mock import Mock
    from biopipelines.msa import MSA

    pipeline = new_pipeline("msa_columns")
    with pipeline:
        upstream = Mock(ids=["aln_1"], streams={"msas": {"format": "a3m", "file": "<id>.a3m"}},
                        tables={"msas": {"columns": ["sequence", "file"]}})
        return MSA(upstream, convert="csv").tables.msas.info.columns


def test_the_carried_provenance_column_is_declared(tmp_path, local_config, isolated_cwd, new_pipeline):
    """An id that does not spell its source keeps `sequences.id`, so the declaration has to name it."""
    written = _convert(tmp_path, [{"id": "aln_1", "sequences.id": "p1", "file": _a3m(tmp_path, "p1")}])
    assert "sequences.id" in written.columns, "the pruning kept the provenance, as this case intends"
    assert set(written.columns) <= set(_declared_columns(isolated_cwd, new_pipeline))


def test_a_pruned_table_is_still_within_the_declaration(tmp_path, local_config, isolated_cwd, new_pipeline):
    """When the row id already spells the source the column is pruned; the declaration is a superset, not an equality."""
    written = _convert(tmp_path, [{"id": "p1", "sequences.id": "p1", "file": _a3m(tmp_path, "p1")}])
    assert "sequences.id" not in written.columns
    assert set(written.columns) <= set(_declared_columns(isolated_cwd, new_pipeline))
