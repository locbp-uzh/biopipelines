# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""ThermoMPNN tool: predicted change in fold stability (ddG) upon point mutation.

ThermoMPNN is a lightweight GNN built on the ProteinMPNN backbone, trained on the
Megascale stability dataset to predict ddG for single point mutants. For each input
structure it runs a single forward pass and scores every position x 19 substitution
(site-saturation), which is materially faster than physics-based ddG estimation.

This wrapper supports two modes:
  * saturation (default): score all positions x 19 substitutions on `chain`.
  * explicit mutation list: score only the substitutions named in `mutations`
    (a selection-style string or a table reference). The native saturation pass
    is still run, then the pipe script filters its rows to the requested mutants.

Negative ddG_pred indicates a predicted stabilising mutation (lower folding free
energy), matching the upstream convention.

Reference:
    Paper: https://www.pnas.org/doi/10.1073/pnas.2314853121
    Repo:  https://github.com/Kuhlman-Lab/ThermoMPNN
"""

import os
from typing import Dict, List, Any, Optional, Tuple, Union

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
    """Serialize the `mutations` selection for the pipe script.

    A plain string passes through; a (TableInfo, column) table reference
    stringifies to a ``TABLE_REFERENCE:<path>:<column>`` token (mirrors the
    Frame2Seq/ProteinMPNN selection convention). Empty -> ``"-"`` (saturation).
    """
    if not value:
        return "-"
    return str(value)


class ThermoMPNN(BaseConfig):
    """
    ThermoMPNN: per-mutation fold-stability change (ddG) from structure.

    Inputs:
        structures:  backbone/complex PDBs (StandardizedOutput or DataStream).
        chain:       chain to score (default "A").
        mutations:   optional. If empty (default), runs site-saturation over all
                     positions of `chain`. Otherwise scores only the named point
                     mutations. Accepts either a selection-style string of
                     wildtype-position-mutant tokens joined by '+' (e.g.
                     "A42G+L50V"), or a table column reference.

    Outputs:
        Tables:
            ddg:      id | structures.id | chain | position | wildtype | mutation | ddG_pred
                      (one row per scored mutation; sorted by ddG_pred ascending,
                      i.e. most-stabilising first)
            missing:  id | removed_by | cause (structures ThermoMPNN could not score)
    """

    TOOL_NAME = "ThermoMPNN"
    TOOL_VERSION = "1.0"

    # ------------------------------------------------------------------
    # Install
    # ------------------------------------------------------------------

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        """Clone Kuhlman-Lab/ThermoMPNN (model weights ship in the repo's
        models/ dir) and create the env. Verification: repo cloned with the
        bundled checkpoint present, env exists, and torch importable."""
        biopipelines = folders.get("biopipelines", "")
        repo_dir = folders.get("ThermoMPNN", "")
        parent_dir = os.path.dirname(repo_dir)
        # The default checkpoint ships with the repo; custom_inference.py finds
        # it automatically under models/. We only assert one of the two
        # filenames the upstream loader looks for is present.
        ckpt_ckpt = os.path.join(repo_dir, "models", "thermoMPNN_default.ckpt")
        ckpt_pt = os.path.join(repo_dir, "models", "thermoMPNN_default.pt")

        env_check = cls._env_exists_check("thermompnn", env_manager)
        repo_check = f'[ -d "{repo_dir}/analysis" ]'
        ckpt_check = f'{{ [ -f "{ckpt_ckpt}" ] || [ -f "{ckpt_pt}" ]; }}'

        skip = "" if force_reinstall else f"""# Check if already installed
if {repo_check} && {ckpt_check} && {env_check} \\
   && {env_manager} run -n thermompnn python -c "import torch, wandb, pytorch_lightning, omegaconf" >/dev/null 2>&1; then
    echo "ThermoMPNN already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        clone_block = f"""mkdir -p {parent_dir}
cd {parent_dir}
if [ ! -d "{repo_dir}/.git" ]; then
    git clone https://github.com/Kuhlman-Lab/ThermoMPNN.git "{repo_dir}"
fi"""

        # The repo's local.yaml hardcodes the authors' cluster path
        # (thermompnn_dir: /proj/kuhl_lab/ThermoMPNN); get_protein_mpnn resolves
        # the ProteinMPNN vanilla weights (vanilla_model_weights/, shipped in the
        # repo) relative to it, so on any other machine the weight load 404s.
        # Repoint it at the actual clone dir.
        patch_local_yaml = (
            f'''sed -i 's|thermompnn_dir:.*|thermompnn_dir: "{repo_dir}"|' '''
            f'"{repo_dir}/local.yaml"'
        )

        remove_block = cls._env_remove_block("thermompnn", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("thermompnn", env_manager, biopipelines)

        return f"""echo "=== Installing ThermoMPNN ==="
{skip}{remove_block}
{clone_block}

# Repoint local.yaml's thermompnn_dir at the clone (vanilla ProteinMPNN weights).
{patch_local_yaml}

{env_block}
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create thermompnn environment."
    exit 1
fi

# Verify: bundled checkpoint + vanilla ProteinMPNN weights present, and the
# inference import chain works (custom_inference -> train_thermompnn -> wandb,
# pytorch_lightning, torch — wandb is a real inference dep, not training-only).
if {{ [ -f "{ckpt_ckpt}" ] || [ -f "{ckpt_pt}" ]; }} \\
   && [ -f "{repo_dir}/vanilla_model_weights/v_48_020.pt" ] \\
   && {env_manager} run -n thermompnn python -c "import torch, wandb, pytorch_lightning, omegaconf" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== ThermoMPNN installation complete ==="
else
    echo "ERROR: ThermoMPNN verification failed (missing checkpoint/weights or import chain)"
    exit 1
fi
"""

    # ------------------------------------------------------------------
    # Lazy path descriptors
    # ------------------------------------------------------------------

    structures_json = Path(lambda self: self.configuration_path("structures.json"))
    ddg_csv = Path(lambda self: self.table_path("ddg"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_thermompnn.py"))

    # ------------------------------------------------------------------
    # Constructor
    # ------------------------------------------------------------------

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 chain: str = "A",
                 mutations: Union[str, Tuple['TableInfo', str]] = "",
                 **kwargs):
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(
                f"structures must be DataStream or StandardizedOutput, got {type(structures).__name__}"
            )
        self.chain = chain
        self.mutations = mutations
        super().__init__(**kwargs)

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        _validate_freeform_string("chain", self.chain)
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
            f"INPUT STRUCTURES: {len(self.structures_stream)}",
            f"CHAIN: {self.chain}",
            f"MODE: {'explicit mutation list' if self.mutations else 'site-saturation'}",
        ])
        if self.mutations:
            lines.append(f"MUTATIONS: {self.mutations}")
        return lines

    # ------------------------------------------------------------------
    # Script generation
    # ------------------------------------------------------------------

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)

        repo_dir = self.folders.get("ThermoMPNN", "")

        script = "#!/bin/bash\n"
        script += "# ThermoMPNN ddG stability script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()  # thermompnn env

        mutations_arg = _mutations_arg(self.mutations)

        script += f"""echo "Running ThermoMPNN (chain={self.chain}) on {len(self.structures_stream)} structure(s)"

python "{self.helper_py}" \\
    --structures-json "{self.structures_json}" \\
    --repo-dir "{repo_dir}" \\
    --chain "{self.chain}" \\
    --mutations "{mutations_arg}" \\
    --ddg-csv "{self.ddg_csv}" \\
    --missing-csv "{self.missing_csv}"

if [ $? -ne 0 ]; then
    echo "Error: ThermoMPNN scoring failed"
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
            "ddg": TableInfo(
                name="ddg",
                path=self.ddg_csv,
                columns=["id", "structures.id", "chain", "position",
                         "wildtype", "mutation", "ddG_pred"],
                description="ThermoMPNN predicted fold-stability change per point mutation (kcal/mol; negative = stabilising)",
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "kind", "cause"],
                description="Structures ThermoMPNN could not score, with reason",
            ),
        }
        return {
            "tables": tables,
            "output_folder": self.output_folder,
        }
