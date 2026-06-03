# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""PLM_Sol tool: sequence-based protein solubility prediction.

PLM_Sol predicts whether a protein is soluble from its amino-acid sequence
alone, using ProtT5-XL embeddings fed to a biLSTM/TextCNN classifier trained on
the updated E. coli solubility dataset (UESolDS). It complements the structure-
based aggregation scorer Aggrescan3D: PLM_Sol needs no structure, only sequence.
Higher solubility = predicted more soluble; the binary call thresholds at 0.5.

GPU NOTE: the first stage generates ProtT5-XL (prottrans_t5_xl_u50) embeddings,
an ~5 GB weight download that is impractically slow on CPU. This tool is intended
for a GPU runtime (Colab GPU / cluster GPU node), like VespaG / ESMFold.

Reference:
    Paper: https://academic.oup.com/bib/article/25/5/bbae404/7739950
    Repo:  https://github.com/Violet969/PLM_Sol
"""

import os
from typing import Dict, List, Any, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream


class PLM_Sol(BaseConfig):
    """
    PLM_Sol: sequence-based protein solubility prediction.

    Scores each input sequence for solubility (ProtT5 embeddings + biLSTM/TextCNN).

    Inputs:
        sequences: amino-acid sequences (StandardizedOutput or DataStream;
                   uses the `sequences` stream / its content-bearing CSV).

    Outputs:
        Tables:
            solubility: id | sequences.id | solubility | soluble
                        (one row per sequence; solubility is the model
                        probability in [0,1], soluble is the >= 0.5 binary call;
                        sorted by solubility descending)
            missing:    id | removed_by | kind | cause (sequences PLM_Sol could
                        not score, with reason)
    """

    TOOL_NAME = "PLM_Sol"
    TOOL_VERSION = "1.0"

    # ------------------------------------------------------------------
    # Install
    # ------------------------------------------------------------------

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        """Create the plm_sol env (pip-installs the deps; PLM_Sol is not on PyPI)
        and clone the repo for its inference code and committed model checkpoint
        (model_param/model_param.t7). ProtT5 embeddings come from transformers
        (T5EncoderModel); the prot_t5_xl_uniref50 weights download lazily on the
        first prediction."""
        biopipelines = folders.get("biopipelines", "")
        repo_dir = folders.get("PLM_Sol", "")
        parent_dir = os.path.dirname(repo_dir)
        weight_check = f'[ -f "{repo_dir}/model_param/model_param.t7" ]'
        env_check = cls._env_exists_check("plm_sol", env_manager)
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_check} && {weight_check} \\
   && {env_manager} run -n plm_sol python -c "import torch, transformers" >/dev/null 2>&1; then
    echo "plm_sol environment already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        clone_block = f"""mkdir -p "{parent_dir}"
if [ ! -d "{repo_dir}/.git" ]; then
    git clone https://github.com/Violet969/PLM_Sol.git "{repo_dir}"
fi"""
        remove_block = cls._env_remove_block("plm_sol", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("plm_sol", env_manager, biopipelines)
        return f"""echo "=== Installing PLM_Sol ==="
{skip}{remove_block}
{clone_block}

{env_block}
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create plm_sol environment."
    exit 1
fi

# Verify: deps importable AND committed checkpoint present in the clone.
if {weight_check} \\
   && {env_manager} run -n plm_sol python -c "import torch, transformers" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== PLM_Sol installation complete ==="
else
    echo "ERROR: PLM_Sol verification failed (cannot import torch/transformers or model checkpoint missing)"
    exit 1
fi
"""

    # ------------------------------------------------------------------
    # Lazy path descriptors
    # ------------------------------------------------------------------

    sequences_json = Path(lambda self: self.configuration_path("sequences.json"))
    solubility_csv = Path(lambda self: self.table_path("solubility"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_plm_sol.py"))

    # ------------------------------------------------------------------
    # Constructor
    # ------------------------------------------------------------------

    def __init__(self,
                 sequences: Union[DataStream, StandardizedOutput],
                 **kwargs):
        if isinstance(sequences, StandardizedOutput):
            self.sequences_stream: DataStream = sequences.streams.sequences
        elif isinstance(sequences, DataStream):
            self.sequences_stream = sequences
        else:
            raise ValueError(
                f"sequences must be DataStream or StandardizedOutput, got {type(sequences).__name__}"
            )
        super().__init__(**kwargs)

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def validate_params(self):
        if not self.sequences_stream or len(self.sequences_stream) == 0:
            raise ValueError("sequences parameter is required and must not be empty")

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
        lines.append(f"INPUT SEQUENCES: {len(self.sequences_stream)}")
        return lines

    # ------------------------------------------------------------------
    # Script generation
    # ------------------------------------------------------------------

    def generate_script(self, script_path: str) -> str:
        self.sequences_stream.save_json(self.sequences_json)

        # Content-bearing sequences CSV (id, sequence) — same source ESMFold/VespaG use.
        sequences_csv = self.sequences_stream.map_table
        repo_dir = self.folders.get("PLM_Sol", "")

        script = "#!/bin/bash\n"
        script += "# PLM_Sol solubility-scoring script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()  # plm_sol env

        script += f"""echo "Running PLM_Sol on {len(self.sequences_stream)} sequence(s)"

python "{self.helper_py}" \\
    --sequences-json "{self.sequences_json}" \\
    --sequences-csv "{sequences_csv}" \\
    --repo-dir "{repo_dir}" \\
    --work-dir "{self.execution_folder}" \\
    --solubility-csv "{self.solubility_csv}" \\
    --missing-csv "{self.missing_csv}"

if [ $? -ne 0 ]; then
    echo "Error: PLM_Sol scoring failed"
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
            "solubility": TableInfo(
                name="solubility",
                path=self.solubility_csv,
                columns=["id", "sequences.id", "solubility", "soluble"],
                description="PLM_Sol sequence-based solubility (probability in [0,1]; "
                            "soluble = probability >= 0.5; higher = more soluble)",
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "kind", "cause"],
                description="Sequences PLM_Sol could not score, with reason",
            ),
        }
        return {
            "tables": tables,
            "output_folder": self.output_folder,
        }
