# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""The id of an `msas` row is the id of the sequence it belongs to, and one function finds it.

MMseqs2 and MSA wrote `<id>_msa` for a stream that is 1:1 with its sequences, so every consumer had to special-case the suffix, and three predictors that matched on `id` alone folded single-sequence without a word. The rule (docs/developer_manual.md, "ID Generation and Provenance"): a 1:1 derived stream keeps the id, and a consumer asks `id_map_utils.MsaLookup`, which also reads the maps written before the rule.

Not covered, because it needs a cluster: a real MMseqs2 search, and the predictors accepting the alignments this hands them.
"""

import ast
import json
import os
import subprocess
import sys

import pandas as pd
import pytest

from biopipelines import id_map_utils
from biopipelines.id_map_utils import MissingAlignment, MsaLookup, prune_redundant_provenance_columns
from tests._helpers import load_pipe_script

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(REPO_ROOT, "pipe_scripts")
PACKAGE = os.path.join(REPO_ROOT, "biopipelines")

PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"

# The three layouts a map can have: written now, written with the provenance column, written with the suffix only.
LAYOUTS = {
    "named_after_the_sequence": [{"id": "p1", "file": "ALN"}],
    "provenance_column": [{"id": "p1_msa", "sequences.id": "p1", "file": "ALN"}],
    "legacy_suffix_only": [{"id": "p1_msa", "file": "ALN"}],
}


def _module(name):
    return load_pipe_script(name, f"{name}.py")


def _write_table(tmp_path, rows, name="msas_map.csv"):
    path = tmp_path / name
    pd.DataFrame(rows).to_csv(path, index=False)
    return str(path)


def _stream_json(tmp_path, table):
    path = tmp_path / "msas.json"
    path.write_text(json.dumps({"map_table": table}))
    return str(path)


def _a3m(tmp_path, name="p1"):
    path = tmp_path / f"{name}.a3m"
    path.write_text(f">{name}\n{PROTEIN}\n")
    return str(path)


# -- the shared lookup -------------------------------------------------------------

class TestTheSharedLookup:

    @pytest.mark.parametrize("layout", list(LAYOUTS))
    def test_exact_provenance_and_legacy_tables_resolve_identically(self, tmp_path, layout):
        rows = [{**row, "file": "/x/p1.a3m"} for row in LAYOUTS[layout]]
        lookup = MsaLookup.load(_write_table(tmp_path, rows))
        assert lookup.path("p1") == "/x/p1.a3m"
        assert lookup.path("p2") is None

    @pytest.mark.parametrize("how", ["stream_json", "stream_dict", "table_csv"])
    def test_the_stream_or_its_table_is_accepted(self, tmp_path, how):
        table = _write_table(tmp_path, [{"id": "p1", "file": "/x/p1.a3m"}])
        source = {"stream_json": _stream_json(tmp_path, table),
                  "stream_dict": {"map_table": table},
                  "table_csv": table}[how]
        lookup = MsaLookup.load(source)
        assert (lookup.path("p1"), lookup.path("p2")) == ("/x/p1.a3m", None)

    def test_the_exact_id_wins_over_a_provenance_or_legacy_row(self, tmp_path):
        rows = [{"id": "p1_msa", "sequences.id": "p1", "file": "/x/legacy.a3m"},
                {"id": "p1", "file": "/x/exact.a3m"}]
        assert MsaLookup.load(_write_table(tmp_path, rows)).path("p1") == "/x/exact.a3m"

    def test_provenance_wins_over_the_legacy_suffix(self, tmp_path):
        rows = [{"id": "p1_msa", "file": "/x/legacy.a3m"},
                {"id": "row7", "sequences.id": "p1", "file": "/x/prov.a3m"}]
        assert MsaLookup.load(_write_table(tmp_path, rows)).path("p1") == "/x/prov.a3m"

    def test_nothing_looser_than_these_tiers_is_tried(self, tmp_path):
        """The child, parent and sibling tiers of get_mapped_ids would hand `prot` the alignment of `prot_1`."""
        rows = [{"id": "prot_1", "file": "/x/one.a3m"}, {"id": "prot", "file": "/x/base.a3m"}]
        lookup = MsaLookup.load(_write_table(tmp_path, rows[:1]))
        assert lookup.path("prot") is None and lookup.path("prot_2") is None
        assert MsaLookup.load(_write_table(tmp_path, rows, "two.csv")).path("prot") == "/x/base.a3m"

    def test_the_sequence_is_the_last_resort_and_only_when_given(self, tmp_path):
        table = _write_table(tmp_path, [{"id": "other", "sequence": PROTEIN, "file": "/x/o.a3m"}])
        lookup = MsaLookup.load(table)
        assert lookup.path("p1") is None
        assert lookup.path("p1", PROTEIN) == "/x/o.a3m"
        assert lookup.path("p1", "MKT") is None

    def test_a_row_without_a_file_answers_nobody(self, tmp_path):
        lookup = MsaLookup.load(_write_table(tmp_path, [{"id": "p1", "file": ""}, {"id": "p2", "file": "/x/p2.a3m"}]))
        assert lookup.path("p1") is None and lookup.ids == ["p2"]

    def test_a_missing_protein_is_named_with_the_ids_that_are_present(self, tmp_path):
        lookup = MsaLookup.load(_write_table(tmp_path, [{"id": "p1", "file": "/x/p1.a3m"}]))
        with pytest.raises(MissingAlignment, match=r"^no alignment in msas for p9 \(ids present: p1\)$"):
            lookup.require("p9")

    def test_the_cause_lists_a_few_ids_and_counts_the_rest(self, tmp_path):
        rows = [{"id": f"p{i}", "file": f"/x/{i}.a3m"} for i in range(8)]
        cause = MsaLookup.load(_write_table(tmp_path, rows)).cause("zz")
        assert cause == "no alignment in msas for zz (ids present: p0, p1, p2, p3, p4, and 3 more)"

    def test_entries_name_the_sequence_each_row_belongs_to(self, tmp_path):
        rows = [{"id": "p1_msa", "sequences.id": "p1", "file": "/x/1.a3m"}, {"id": "p2", "file": "/x/2.a3m"}]
        assert MsaLookup.load(_write_table(tmp_path, rows)).entries == [("p1", "/x/1.a3m"), ("p2", "/x/2.a3m")]

    def test_a_table_that_is_not_there_or_has_no_file_column_is_refused(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            MsaLookup.load(str(tmp_path / "gone.csv"))
        with pytest.raises(ValueError, match="no 'file' column"):
            MsaLookup.load(_write_table(tmp_path, [{"id": "p1", "sequence": PROTEIN}]))
        with pytest.raises(ValueError, match="names no map_table"):
            MsaLookup.load({"ids": ["p1"]})

    def test_the_lookup_imports_with_no_package_around_it(self):
        """The pipe scripts that share it may put biopipelines/ on sys.path and import it by bare name."""
        code = ("import sys; sys.path.insert(0, %r)\n"
                "import id_map_utils\n"
                "rows = [{'id': 'p1_msa', 'file': '/x/p1.a3m'}]\n"
                "print('RESULT:' + str(id_map_utils.MsaLookup(rows).path('p1')))\n" % PACKAGE)
        done = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
        assert done.returncode == 0, done.stderr
        assert "RESULT:/x/p1.a3m" in done.stdout

    def test_the_framework_prunes_a_provenance_column_that_only_repeats_the_id(self):
        """Why the producers do not write `sequences.id` beside an id that already is the sequence id."""
        repeated = pd.DataFrame({"id": ["p1", "p2"], "sequences.id": ["p1", "p2"], "file": ["a", "b"]})
        assert "sequences.id" not in prune_redundant_provenance_columns(repeated).columns
        renamed = pd.DataFrame({"id": ["Panda_1"], "sequences.id": ["LID_1"], "file": ["a"]})
        assert "sequences.id" in prune_redundant_provenance_columns(renamed).columns


# -- producers ---------------------------------------------------------------------

class TestMMseqs2WritesTheSequenceId:

    def test_every_branch_that_writes_a_row_keeps_the_sequence_id(self, tmp_path):
        mm = _module("pipe_mmseqs2_sequences")
        a3m = tmp_path / "raw.a3m"
        a3m.write_text(f"#33\t1\n>101\n{PROTEIN}\n>hit\n{PROTEIN}\n")
        csv_in = tmp_path / "raw.csv"
        pd.DataFrame({"key": [-1, -1], "sequence": [PROTEIN, PROTEIN]}).to_csv(csv_in, index=False)

        (row_csv,) = mm.convert_a3m_to_csv_format(str(a3m), "p1", str(tmp_path / "p1.csv"))
        (row_a3m,) = mm.write_a3m_output(str(a3m), "p1", str(tmp_path / "p1.a3m"))
        (row_plain,) = mm.process_csv_output(str(csv_in), "p1", str(tmp_path / "out.csv"))
        for row in (row_csv, row_a3m, row_plain):
            assert row["id"] == "p1"
            assert set(row) == {"id", "sequence", "file"}

    def test_a_reused_alignment_and_the_empty_table_follow_the_same_rule(self, tmp_path):
        script = os.path.join(SCRIPTS, "pipe_mmseqs2_sequences.py")
        pd.DataFrame([{"id": "p1", "sequence": PROTEIN}]).to_csv(tmp_path / "in.csv", index=False)
        out_dir = tmp_path / "msas"
        out_dir.mkdir()
        (out_dir / "p1.a3m").write_text(f">101\n{PROTEIN}\n")
        done = subprocess.run(
            [sys.executable, script, str(tmp_path / "in.csv"), str(out_dir / "msas_map.csv"), "unused.sh",
             "--output_format", "a3m", "--server_dir", str(tmp_path / "server"),
             "--server_url", "https://example.invalid", "--missing_csv", str(tmp_path / "missing.csv")],
            capture_output=True, text=True)
        assert done.returncode == 0, done.stdout + done.stderr
        table = pd.read_csv(out_dir / "msas_map.csv")
        assert table.columns.tolist() == ["id", "sequence", "file"]
        assert table["id"].tolist() == ["p1"]
        assert MsaLookup.load(str(out_dir / "msas_map.csv")).path("p1") == str(out_dir / "p1.a3m")

    def test_the_declared_ids_are_the_sequence_ids_and_the_table_has_no_repeated_column(self):
        from biopipelines.mmseqs2 import MMseqs2

        class Ids:
            ids = ["p1", "p2"]

            def __len__(self):
                return 2

        tool = MMseqs2.__new__(MMseqs2)
        tool.sequences_stream, tool.raw_sequences, tool.output_format = Ids(), None, "a3m"
        tool.output_msa_csv, tool.missing_csv, tool.output_folder = "/o/msas.csv", "/o/missing.csv", "/o"
        tool.stream_path = lambda stream, name: f"/o/{stream}/{name}"
        declared = tool.get_output_files()
        assert list(declared["msas"].ids) == ["p1", "p2"]
        assert declared["tables"]["msas"].info.columns == ["id", "sequence", "file"]


class TestTheMsaToolKeepsTheInputId:

    @staticmethod
    def _run(tmp_path, rows, convert):
        table = _write_table(tmp_path, rows, "in.csv")
        config = tmp_path / "config.json"
        config.write_text(json.dumps({
            "input_msa_table": table, "convert": convert,
            "output_folder": str(tmp_path / "msas"), "output_msas_csv": str(tmp_path / "out.csv"),
            "local_missing_csv": str(tmp_path / "007_MSA" / "_execution" / "local_missing.csv")}))
        done = subprocess.run([sys.executable, os.path.join(SCRIPTS, "pipe_msa.py"), "--config", str(config)],
                              capture_output=True, text=True)
        assert done.returncode == 0, done.stdout + done.stderr
        return pd.read_csv(tmp_path / "out.csv")

    def test_ids_pass_through_and_no_suffix_is_added(self, tmp_path):
        a3m = _a3m(tmp_path)
        out = self._run(tmp_path, [{"id": "p1", "file": a3m, "sequence": PROTEIN}], "csv")
        assert out["id"].tolist() == ["p1"]
        assert out.columns.tolist() == ["id", "sequence", "file"]
        assert os.path.basename(out["file"][0]) == "p1.csv"

    def test_a_legacy_row_keeps_its_id_and_a_chain_of_conversions_does_not_pile_up_suffixes(self, tmp_path):
        a3m = _a3m(tmp_path)
        first = self._run(tmp_path, [{"id": "p1_msa", "sequences.id": "p1", "file": a3m}], "csv")
        assert first["id"].tolist() == ["p1_msa"]
        again = self._run(tmp_path, [{"id": "p1_msa", "file": first["file"][0]}], "a3m")
        assert again["id"].tolist() == ["p1_msa"]

    def test_a_provenance_the_id_does_not_spell_is_kept(self, tmp_path):
        a3m = _a3m(tmp_path)
        out = self._run(tmp_path, [{"id": "Panda_1", "sequences.id": "LID_1", "file": a3m}], "csv")
        assert out["sequences.id"].tolist() == ["LID_1"]
        assert MsaLookup.load(str(tmp_path / "out.csv")).path("LID_1") == out["file"][0]

    def test_the_declared_ids_are_the_input_ids(self):
        from biopipelines.msa import MSA

        class Ids:
            ids = ["p1", "p2_msa"]
            format = "a3m"

        tool = MSA.__new__(MSA)
        tool.msas_stream, tool.convert = Ids(), "csv"
        tool.output_msas_folder, tool.output_msas_csv, tool.output_folder = "/o/msas", "/o/msas.csv", "/o"
        tool.missing_csv = "/o/missing.csv"
        declared = tool.get_output_files()
        assert list(declared["msas"].ids) == ["p1", "p2_msa"]
        assert declared["tables"]["msas"].info.columns == ["id", "sequences.id", "sequence", "file"]


class TestMutagenesisFindsItsParentThroughTheLookup:

    @pytest.mark.parametrize("layout", list(LAYOUTS))
    def test_a_mutant_gets_its_parents_alignment_whatever_the_map_was_called(self, tmp_path, layout):
        a3m = _a3m(tmp_path)
        table = _write_table(tmp_path, [{**row, "file": a3m} for row in LAYOUTS[layout]], "in.csv")
        pd.DataFrame([{"id": "p1_5A", "original.id": "p1", "sequence": "A", "mutation_positions": "5",
                       "new_aa": "W"}]).to_csv(tmp_path / "mutants.csv", index=False)
        done = subprocess.run(
            [sys.executable, os.path.join(SCRIPTS, "pipe_mutagenesis_msa.py"),
             "--mutants", str(tmp_path / "mutants.csv"), "--input-msas", table,
             "--output-folder", str(tmp_path / "out"), "--output", str(tmp_path / "out.csv")],
            capture_output=True, text=True)
        assert done.returncode == 0, done.stdout + done.stderr
        out = pd.read_csv(tmp_path / "out.csv")
        assert out.columns.tolist() == ["id", "original.id", "sequence", "file"]
        assert out["id"].tolist() == ["p1_5A"]
        assert open(out["file"][0]).read().splitlines()[1] == PROTEIN[:4] + "W" + PROTEIN[5:]


# -- the rule, stated as source ------------------------------------------------------

def _tree(path):
    with open(path, encoding="utf-8") as handle:
        return ast.parse(handle.read())


def _python_files():
    for folder in (PACKAGE, SCRIPTS):
        for name in sorted(os.listdir(folder)):
            if name.endswith(".py"):
                yield os.path.join(folder, name)


class TestNoModuleBuildsTheSuffix:

    def test_no_module_builds_an_id_with_the_msa_suffix(self):
        """`<id>_msa` is read, as a legacy form, in id_map_utils and nowhere else."""
        offenders = []
        for path in _python_files():
            if os.path.basename(path) == "id_map_utils.py":
                continue
            for node in ast.walk(_tree(path)):
                built = (
                    (isinstance(node, ast.JoinedStr) and node.values
                     and isinstance(node.values[-1], ast.Constant) and str(node.values[-1].value).endswith("_msa"))
                    or (isinstance(node, ast.Constant) and node.value == "_msa")
                    or (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                        and node.func.attr in ("endswith", "removesuffix", "replace")
                        and node.args and isinstance(node.args[0], ast.Constant) and node.args[0].value == "_msa"))
                if built:
                    offenders.append(f"{os.path.relpath(path, REPO_ROOT)}:{node.lineno}")
        assert not offenders, f"ids built or stripped with the `_msa` suffix: {offenders}"

    def test_the_lookup_is_defined_once(self):
        homes = [os.path.relpath(path, REPO_ROOT) for path in _python_files()
                 if any(isinstance(node, ast.ClassDef) and node.name == "MsaLookup"
                        for node in ast.walk(_tree(path)))]
        assert homes == [os.path.join("biopipelines", "id_map_utils.py")]
