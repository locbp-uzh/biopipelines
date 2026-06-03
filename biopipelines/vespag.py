# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""VespaG tool: zero-shot single-substitution fitness prediction from sequence.

VespaG is a small expert-guided model distilled from ESM-2 that predicts the
functional effect of single amino-acid substitutions, leaderboard-competitive on
ProteinGym while being MSA-free and fast. It complements ThermoMPNN: ThermoMPNN
scores fold *stability* (ddG) from structure, VespaG scores evolutionary/functional
*fitness* from sequence. Higher VespaG score = predicted more tolerated/beneficial.

This wrapper supports two modes:
  * saturation (default): score the whole single-site mutational landscape of
    each input sequence (all positions x 19 substitutions).
  * explicit mutation list: score only the substitutions named in `mutations`
    (a selection-style string or a table reference), translated into a VespaG
    mutation file at runtime.

GPU NOTE: VespaG generates ESM-2 (esm2_t36_3B_UR50D) embeddings internally on
first use — an ~11 GB weight download and impractically slow on CPU. This tool is
intended for a GPU runtime (Colab GPU / cluster GPU node), like P2Rank/NeuralPLexer.

Reference:
    Paper: https://academic.oup.com/bioinformatics/article/40/11/btae621/7907184
    Repo:  https://github.com/JSchlensok/VespaG
"""

import os
from typing import Dict, List, Any, Tuple, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream


def _mutations_arg(value) -> str:
    """Serialize the `mutations` selection for the pipe script (see ThermoMPNN).

    Plain string passes through; a (TableInfo, column) reference becomes a
    ``TABLE_REFERENCE:<path>:<column>`` token. Empty -> ``"-"`` (saturation).
    """
    if not value:
        return "-"
    return str(value)


class VespaG(BaseConfig):
    """
    VespaG: zero-shot single-substitution fitness scoring from sequence.

    Inputs:
        sequences:   amino-acid sequences (StandardizedOutput or DataStream;
                     uses the `sequences` stream / its content-bearing CSV).
        mutations:   optional. If empty (default), scores the full single-site
                     saturation landscape of each sequence. Otherwise scores only
                     the named point mutations. Accepts either a selection-style
                     string of wildtype-position-mutant tokens joined by '+'
                     (e.g. "A42G+L50V", 1-indexed), or a table column reference.

    Outputs:
        Tables:
            fitness:  id | sequences.id | position | wildtype | mutation | fitness
                      (one row per scored mutation; sorted by fitness descending,
                      i.e. most-tolerated/beneficial first)
            missing:  id | removed_by | cause (sequences VespaG could not score)
    """

    TOOL_NAME = "VespaG"
    TOOL_VERSION = "1.0"

    # ------------------------------------------------------------------
    # Install
    # ------------------------------------------------------------------

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        """Create the vespag env (pip-installs the package from GitHub source —
        VespaG is not on PyPI despite its README) and clone the repo for its
        committed model_weights/ (the small FNN heads; NOT included in the pip
        package, and `load_model` resolves them relative to the runtime cwd, so
        the pipe script runs `vespag predict` from this clone). The heavy ESM-2
        3B embedding weights download lazily on first prediction."""
        biopipelines = folders.get("biopipelines", "")
        repo_dir = folders.get("VespaG", "")
        parent_dir = os.path.dirname(repo_dir)
        weight_check = f'[ -f "{repo_dir}/model_weights/v2/esm2.pt" ]'
        env_check = cls._env_exists_check("vespag", env_manager)
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_check} && {weight_check} \\
   && {env_manager} run -n vespag python -c "import vespag, torch" >/dev/null 2>&1; then
    echo "vespag environment already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        clone_block = f"""mkdir -p "{parent_dir}"
if [ ! -d "{repo_dir}/.git" ]; then
    git clone https://github.com/JSchlensok/VespaG.git "{repo_dir}"
fi"""
        remove_block = cls._env_remove_block("vespag", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("vespag", env_manager, biopipelines)
        return f"""echo "=== Installing VespaG ==="
{skip}{remove_block}
{clone_block}

{env_block}
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create vespag environment."
    exit 1
fi

# Verify: package importable AND committed FNN weights present in the clone.
if {weight_check} \\
   && {env_manager} run -n vespag python -c "import vespag, torch" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== VespaG installation complete ==="
else
    echo "ERROR: VespaG verification failed (cannot import vespag/torch or model weights missing)"
    exit 1
fi
"""

    # ------------------------------------------------------------------
    # Lazy path descriptors
    # ------------------------------------------------------------------

    sequences_json = Path(lambda self: self.configuration_path("sequences.json"))
    fitness_csv = Path(lambda self: self.table_path("fitness"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_vespag.py"))

    # ------------------------------------------------------------------
    # Constructor
    # ------------------------------------------------------------------

    def __init__(self,
                 sequences: Union[DataStream, StandardizedOutput],
                 mutations: Union[str, Tuple['TableInfo', str]] = "",
                 **kwargs):
        if isinstance(sequences, StandardizedOutput):
            self.sequences_stream: DataStream = sequences.streams.sequences
        elif isinstance(sequences, DataStream):
            self.sequences_stream = sequences
        else:
            raise ValueError(
                f"sequences must be DataStream or StandardizedOutput, got {type(sequences).__name__}"
            )
        self.mutations = mutations
        super().__init__(**kwargs)

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def validate_params(self):
        if not self.sequences_stream or len(self.sequences_stream) == 0:
            raise ValueError("sequences parameter is required and must not be empty")
        if isinstance(self.mutations, str):
            _validate_freeform_string("mutations", self.mutations)

    # ------------------------------------------------------------------
    # Configure inputs
    # ------------------------------------------------------------------

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    # ------------------------------------------------------------------
    # Config display
    # ------------------------------------------------------------------

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.extend([
            f"INPUT SEQUENCES: {len(self.sequences_stream)}",
            f"MODE: {'explicit mutation list' if self.mutations else 'site-saturation'}",
        ])
        if self.mutations:
            lines.append(f"MUTATIONS: {self.mutations}")
        return lines

    # ------------------------------------------------------------------
    # Script generation
    # ------------------------------------------------------------------

    def generate_script(self, script_path: str) -> str:
        self.sequences_stream.save_json(self.sequences_json)

        # Content-bearing sequences CSV (id, sequence) — same source ESMFold uses.
        sequences_csv = self.sequences_stream.map_table
        repo_dir = self.folders.get("VespaG", "")

        script = "#!/bin/bash\n"
        script += "# VespaG fitness-scoring script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()  # vespag env

        mutations_arg = _mutations_arg(self.mutations)

        script += f"""echo "Running VespaG on {len(self.sequences_stream)} sequence(s)"

python "{self.helper_py}" \\
    --sequences-json "{self.sequences_json}" \\
    --sequences-csv "{sequences_csv}" \\
    --mutations "{mutations_arg}" \\
    --repo-dir "{repo_dir}" \\
    --fitness-csv "{self.fitness_csv}" \\
    --missing-csv "{self.missing_csv}"

if [ $? -ne 0 ]; then
    echo "Error: VespaG scoring failed"
    exit 1
fi

"""
        script += self.generate_completion_check_footer()
        return script

    # ------------------------------------------------------------------
    # Output prediction
    # ------------------------------------------------------------------

    def get_output_files(self) -> Dict[str, Any]:
        tables = {
            "fitness": TableInfo(
                name="fitness",
                path=self.fitness_csv,
                columns=["id", "sequences.id", "position",
                         "wildtype", "mutation", "fitness"],
                description="VespaG zero-shot single-substitution fitness score (higher = more tolerated/beneficial)",
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "kind", "cause"],
                description="Sequences VespaG could not score, with reason",
            ),
        }
        return {
            "tables": tables,
            "output_folder": self.output_folder,
        }
