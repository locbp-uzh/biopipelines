# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""TMalign tool: sequential and circularly permuted TM-align of structure pairs, with the TM-score difference.

Runs TM-align once per pair of a pair table, plain and with `-cp`, and reports TM-score, RMSD, aligned length, the residue alignment and the superposition of each, then the TM-score the permutation gained.

Reference: Zhang and Skolnick, Nucleic Acids Research 33, 2302-2309 (2005), doi:10.1093/nar/gki524. https://github.com/pylelab/USalign
"""

import os
from typing import Dict, List, Any, Union

try:
    from .base_config import StandardizedOutput, TableInfo
    from .datastream import DataStream
    from .pair_alignment import PairAligner, MODE_SEQUENTIAL, MODE_CIRCULAR_PERMUTATION
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import StandardizedOutput, TableInfo
    from datastream import DataStream
    from pair_alignment import PairAligner, MODE_SEQUENTIAL, MODE_CIRCULAR_PERMUTATION


ALIGNMENT_COLUMNS = [
    "id", "structure_1", "structure_2", "mode",
    "length_1", "length_2", "aligned_length", "rmsd", "tm_score_1", "tm_score_2",
    "rotation_11", "rotation_12", "rotation_13",
    "rotation_21", "rotation_22", "rotation_23",
    "rotation_31", "rotation_32", "rotation_33",
    "translation_1", "translation_2", "translation_3",
    "aligned_1", "aligned_2",
]
DELTA_COLUMNS = ["id", "structure_1", "structure_2", "delta_tm_1", "delta_tm_2"]


class TMalign(PairAligner):
    """
    TMalign: align every pair of a table of structure pairs in sequence order and with a circular permutation.

    For each pair (structure_1, structure_2) the tool runs TM-align plain and with
    `-cp`, and reports one row per pair and alignment plus, per pair, the TM-score
    the permutation gained. A positive `delta_tm` says the pair is better described
    as a circular permutant than as two folds that share an order. It answers the
    same question as USalign without the non-sequential alignment, and adds the
    TM-score difference; use USalign when the non-sequential alignment is wanted too.

    Structure_1 is the one that is rotated onto structure_2. Each TM-score is
    normalized by the length of one structure: `tm_score_1` by structure_1,
    `tm_score_2` by structure_2. TM-align recommends the one normalized by the
    reference structure, `tm_score_2`.

    Only the first chain of each structure is aligned. A pair whose ids are not both
    in `structures`, or whose alignment fails, keeps its rows with NaN values and is
    recorded in `missing`.

    Inputs:
        structures: PDB or mmCIF structures; every id used in the pair table must be one of them.
        structure_1: column of the pair table holding the id of the structure that is rotated,
            as `pairs.tables.data.<column>` or `(pairs.tables.data, "<column>")`.
        structure_2: column of the same table holding the id of the structure it is rotated onto.

    Outputs:
        Tables:
            alignments: id | structure_1 | structure_2 | mode | length_1 | length_2 |
                aligned_length | rmsd | tm_score_1 | tm_score_2 | rotation_11 ... rotation_33 |
                translation_1 ... translation_3 | aligned_1 | aligned_2
                One row per distinct pair and mode. `id` is `<structure_1>+<structure_2>_<mode>`,
                `mode` one of sequential, circular_permutation. `rotation_ij` and `translation_i`
                are the matrix u and vector t of TM-align, so a coordinate x of structure_1 lands
                at t + u x. `aligned_1` and `aligned_2` are the two aligned sequences, `-` for a
                gap; a `*` in `aligned_1` marks the circular permutation point of a permuted
                alignment.
            delta: id | structure_1 | structure_2 | delta_tm_1 | delta_tm_2
                One row per distinct pair, `id` being `<structure_1>+<structure_2>`. `delta_tm_1`
                is the TM-score of the circular_permutation alignment minus the sequential one,
                both normalized by structure_1; `delta_tm_2` the same normalized by structure_2.
                NaN where either alignment gave no TM-score. It is 0 when TM-align finds
                no circular permutation, since it then returns the sequential alignment.
            missing: id | removed_by | kind | cause

    Usage::

        pairs = Table("pairs.csv")        # one row per pair, columns "parent" and "variant"
        folds = ESMFold2(proteins=seqs)
        cp = TMalign(structures=folds,
                     structure_1=pairs.tables.data.variant,
                     structure_2=pairs.tables.data.parent)

        # pairs the permutation explains
        permuted = Panda(tables=cp.tables.delta,
                         operations=[Panda.filter("delta_tm_2 > 0.1")])

    See Also:
        USalign: the same two alignments plus the non-sequential alignment.
    """

    TOOL_NAME = "TMalign"
    TOOL_VERSION = "1.2"
    ENV_NAME = "usalign"
    BINARY = "TMalign"
    REPORTS_DELTA = True

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        # The usalign env ships TMalign beside USalign; delegate so a user installing TMalign need not know that.
        from .usalign import USalign
        return USalign._install_script(folders, env_manager=env_manager,
                                       force_reinstall=force_reinstall, **kwargs)

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 structure_1: Any,
                 structure_2: Any,
                 **kwargs):
        self._set_inputs(structures, structure_1, structure_2)
        super().__init__(**kwargs)

    def modes(self) -> List[Dict[str, Any]]:
        return [
            {"label": MODE_SEQUENTIAL, "args": []},
            {"label": MODE_CIRCULAR_PERMUTATION, "args": ["-cp"]},
        ]

    def get_output_files(self) -> Dict[str, Any]:
        tables = {
            "alignments": TableInfo(
                name="alignments",
                path=self.alignments_csv,
                columns=ALIGNMENT_COLUMNS,
                description="TM-score, RMSD, aligned length, residue alignment and superposition of every pair in both modes",
            ),
            "delta": TableInfo(
                name="delta",
                path=self.delta_csv,
                columns=DELTA_COLUMNS,
                description="TM-score gained by the circular permutation, per pair and normalization",
            ),
            "missing": self.missing_table_info(self.missing_csv),
        }
        return {
            "tables": tables,
            "output_folder": self.output_folder,
        }
