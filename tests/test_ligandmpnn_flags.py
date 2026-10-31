# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""LigandMPNN: the two newer switches are checked like their siblings, and `save_stats` documents the file upstream writes.

`pack_side_chains` and `pack_with_ligand_context` raise unless they are bools, but `save_stats` and `use_atom_context`, added in the same class, were not checked: `save_stats="yes"` silently emitted `--save_stats 1` and `use_atom_context=0` silently behaved like False.

Upstream `run.py` saves the stats of a structure as `<out_folder>/stats/<input file name without .pdb>.pt`, the name its FASTA carries, not the stream id; the docstring and the docs section promised `<id>.pt`. The two coincide only when the structure's file is `<id>.pdb`.
"""

import os
import re

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _build(new_pipeline, **kwargs):
    from biopipelines.ligand import Ligand
    from biopipelines.ligand_mpnn import LigandMPNN
    from biopipelines.mock import Mock

    pipeline = new_pipeline("lmpnn_flags")
    with pipeline:
        structures = Mock(ids=["s1"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
                          map_table_strategy="config")
        return LigandMPNN(structures=structures.streams.structures, ligand=Ligand(code="ATP"), **kwargs)


class TestTheSwitchesAreBools:

    @pytest.mark.parametrize("name", ["save_stats", "use_atom_context"])
    @pytest.mark.parametrize("value", ["yes", 0, 1, None])
    def test_a_value_that_is_not_a_bool_is_refused(self, name, value, local_config, isolated_cwd, new_pipeline):
        with pytest.raises(ValueError, match=f"{name} must be a bool"):
            _build(new_pipeline, **{name: value})

    @pytest.mark.parametrize("name", ["save_stats", "use_atom_context"])
    @pytest.mark.parametrize("value", [True, False])
    def test_a_bool_is_accepted(self, name, value, local_config, isolated_cwd, new_pipeline):
        _build(new_pipeline, **{name: value})


class TestSaveStatsDocumentsTheFileUpstreamWrites:

    @staticmethod
    def _save_stats_text():
        from biopipelines.ligand_mpnn import LigandMPNN

        docs = open(os.path.join(REPO_ROOT, "docs", "tool", "sequence_design.md"), encoding="utf-8").read()
        line = next(l for l in docs.splitlines() if l.startswith("- `save_stats`"))
        doc = LigandMPNN.__init__.__doc__
        docstring = doc[doc.index("save_stats:"):doc.index("use_atom_context:")]
        return {"docs": line, "docstring": " ".join(docstring.split())}

    @pytest.mark.parametrize("where", ["docs", "docstring"])
    def test_the_promised_name_is_the_input_file_not_the_id(self, where):
        text = self._save_stats_text()[where]
        assert "stats/<id>.pt" not in text
        assert re.search(r"stats/<input file name[^>]*>\.pt", text), text
        assert "<id>.pdb" in text
