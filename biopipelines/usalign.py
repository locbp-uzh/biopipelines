# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""USalign tool: sequential, circularly permuted and non-sequential alignment of structure pairs.

Runs US-align once per pair of a pair table and per alignment mode, and reports TM-score, RMSD, aligned length, the residue alignment and the superposition of every one.

Reference: Zhang et al., Nature Methods 19, 1109-1115 (2022), doi:10.1038/s41592-022-01585-1. https://github.com/pylelab/USalign
"""

import os
from typing import Dict, List, Any, Union

try:
    from .base_config import StandardizedOutput, TableInfo
    from .datastream import DataStream
    from .pair_alignment import (
        PairAligner, MODE_SEQUENTIAL, MODE_CIRCULAR_PERMUTATION, MODE_NON_SEQUENTIAL,
    )
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import StandardizedOutput, TableInfo
    from datastream import DataStream
    from pair_alignment import (
        PairAligner, MODE_SEQUENTIAL, MODE_CIRCULAR_PERMUTATION, MODE_NON_SEQUENTIAL,
    )


# US-align `-mm` values: 5 is the fully non-sequential (fNS) alignment, 6 the semi-non-sequential (sNS) one.
NON_SEQUENTIAL_VARIANTS = {"fNS": 5, "sNS": 6}

ALIGNMENT_COLUMNS = [
    "id", "structure_1", "structure_2", "mode",
    "length_1", "length_2", "aligned_length", "rmsd", "tm_score_1", "tm_score_2",
    "rotation_11", "rotation_12", "rotation_13",
    "rotation_21", "rotation_22", "rotation_23",
    "rotation_31", "rotation_32", "rotation_33",
    "translation_1", "translation_2", "translation_3",
    "aligned_1", "aligned_2",
]


class USalign(PairAligner):
    """
    USalign: align every pair of a table of structure pairs three ways with US-align.

    For each pair (structure_1, structure_2) the tool runs the sequential alignment (`-mm 0`), the circularly permuted alignment (`-mm 3`) and the non-sequential alignment (`-mm 5`, or `-mm 6` with `non_sequential_variant="sNS"`), and reports one row per pair and alignment. Use it to ask whether a permuted or rearranged design keeps the fold of its parent: a circularly permuted protein aligns poorly in sequence order and well in the other two.

    Structure_1 is the one that is rotated onto structure_2. Each TM-score is normalized by the length of one structure: `tm_score_1` by structure_1, `tm_score_2` by structure_2. US-align recommends the one normalized by the reference structure, `tm_score_2`.

    Only the first chain of each structure is aligned (US-align's default for a pair of monomers). A pair whose ids are not both in `structures` keeps its rows with NaN values and is recorded in `missing` as `unmatched`; a pair whose alignment fails is recorded as `failure`.

    Inputs:
        structures: PDB or mmCIF structures; every id used in the pair table must be one of them.
        structure_1: column of the pair table holding the id of the structure that is rotated, as `pairs.tables.data.<column>` or `(pairs.tables.data, "<column>")`.
        structure_2: column of the same table holding the id of the structure it is rotated onto.
        non_sequential_variant: "fNS" (fully non-sequential, `-mm 5`) or "sNS" (semi-non-sequential, `-mm 6`, which keeps the sequence order within secondary-structure elements) (default "fNS").

    Outputs:
        Tables:
            alignments: id | structure_1 | structure_2 | mode | length_1 | length_2 | aligned_length | rmsd | tm_score_1 | tm_score_2 | rotation_11 ... rotation_33 | translation_1 ... translation_3 | aligned_1 | aligned_2
                One row per distinct pair and mode. `id` is `<structure_1>+<structure_2>_<mode>`, `mode` one of sequential, circular_permutation, non_sequential. `rotation_ij` and `translation_i` are the matrix u and vector t of US-align, so a coordinate x of structure_1 lands at t + u x. `aligned_1` and `aligned_2` are the two aligned sequences, `-` for a gap; a `*` in `aligned_1` marks the circular permutation point of a permuted alignment.
            missing: id | removed_by | kind | cause   (kind="unmatched" for a pair with an id absent from `structures`, kind="failure" for an alignment that failed)

    Usage::

        pairs = Table("pairs.csv")        # one row per pair, columns "parent" and "variant"
        folds = ESMFold2(proteins=seqs)
        cp = USalign(structures=folds,
                     structure_1=pairs.tables.data.variant,
                     structure_2=pairs.tables.data.parent)

        # variants whose permuted alignment to the parent is good
        permuted = Panda(tables=cp.tables.alignments,
                     operations=[Panda.filter("mode == 'circular_permutation' and tm_score_2 > 0.7")])

    See Also:
        TMalign: the sequential and circularly permuted alignments of TM-align, with their TM-score difference.
    """

    TOOL_NAME = "USalign"
    TOOL_VERSION = "1.5"
    ENV_NAME = "usalign"
    BINARY = "USalign"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        env = cls._install_env(env_manager)
        biopipelines = folders.get("biopipelines", "")
        env_check = cls._env_exists_check(env, env_manager)
        run = cls._env_run(env, env_manager)
        verified = f'{run}USalign -h >/dev/null 2>&1 && {run}TMalign -h >/dev/null 2>&1 && {run}python -c "import pandas, yaml" >/dev/null 2>&1'
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_check} && {verified}; then
    echo "USalign already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block(env, env_manager) if force_reinstall else ""
        env_block = cls._env_install_block(env, env_manager, biopipelines)
        return f"""echo "=== Installing USalign ==="
{skip}{remove_block}
if ! {env_check}; then
{env_block}
else
    echo "{env} environment already exists, skipping creation."
fi
{cls._link_conda_binaries(env, env_manager, ("USalign", "TMalign"))}
if {verified}; then
    touch "$INSTALL_SUCCESS"
    echo "=== USalign installation complete ==="
else
    echo "ERROR: USalign verification failed (USalign or TMalign not on PATH, or pandas/yaml not importable). Retry with force_reinstall=True to rebuild the env."
    exit 1
fi
"""

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 structure_1: Any,
                 structure_2: Any,
                 non_sequential_variant: str = "fNS",
                 **kwargs):
        self.non_sequential_variant = non_sequential_variant
        self._set_inputs(structures, structure_1, structure_2)
        super().__init__(**kwargs)

    def validate_params(self):
        super().validate_params()
        if self.non_sequential_variant not in NON_SEQUENTIAL_VARIANTS:
            raise ValueError(
                f"non_sequential_variant must be one of {tuple(NON_SEQUENTIAL_VARIANTS)}, "
                f"got {self.non_sequential_variant!r}")

    def modes(self) -> List[Dict[str, Any]]:
        return [
            {"label": MODE_SEQUENTIAL, "args": ["-mm", "0"]},
            {"label": MODE_CIRCULAR_PERMUTATION, "args": ["-mm", "3"]},
            {"label": MODE_NON_SEQUENTIAL,
             "args": ["-mm", str(NON_SEQUENTIAL_VARIANTS[self.non_sequential_variant])]},
        ]

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"NON-SEQUENTIAL VARIANT: {self.non_sequential_variant}")
        return lines

    def get_output_files(self) -> Dict[str, Any]:
        tables = {
            "alignments": TableInfo(
                name="alignments",
                path=self.alignments_csv,
                columns=ALIGNMENT_COLUMNS,
                description="TM-score, RMSD, aligned length, residue alignment and superposition of every pair in every mode",
            ),
            "missing": self.missing_table_info(self.missing_csv),
        }
        return {
            "tables": tables,
            "output_folder": self.output_folder,
        }
