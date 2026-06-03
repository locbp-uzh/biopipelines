# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""AF2BIND tool: small-molecule binding-site prediction from a protein structure.

AF2BIND probes a target protein with 20 "bait" amino acids and reads out
AlphaFold2's internal pair representation through a trained linear head,
yielding a per-residue binding probability ``p_bind`` without needing any
ligand. For each input structure the wrapper builds a ColabDesign ``binder``
model on the requested chain, runs a single AF2 pass, applies the AF2BIND
head, and writes per-residue scores plus a top-k binding-residue selection.

Reference: Gazizov, Lian, et al. (2023) AF2BIND. https://github.com/sokrypton/af2bind
"""

import os
from typing import Dict, List, Any, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from ._weights_cache import ensure_weights_block
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from _weights_cache import ensure_weights_block


# AF2BIND runs ColabDesign v1.1.1, which needs jax with jax.linear_util
# (removed in jax>=0.4.24). The cluster uses a dedicated `af2bind` conda env
# pinned to jax 0.4.23+cuda11 (environments/af2bind.yaml); on Colab, ColabDesign
# runs in the base Python's JAX/GPU stack.
#
# AlphaFold2 network params are SHARED across AF2-consuming tools via the
# `AlphaFoldParams` cache folder (its `params/` subdir holds params_model_*.npz,
# = ColabDesign's data_dir). The install ensures that folder exists and downloads
# the params only if absent, so a co-located AlphaFold install's weights are
# reused. The AF2BIND linear-head weights are tool-specific and live under the
# AF2BIND folder.
AF2_PARAMS_URL = "https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar"
AF2BIND_WEIGHTS_URL = "https://github.com/sokrypton/af2bind/raw/main/attempt_7_2k_lam0-03.zip"
# ColabDesign v1.1.1 (a branch in the upstream repo, matching the af2bind
# notebook's own recipe). git over https works through the cluster proxy via
# the http(s)_proxy env vars sourced in scheduler.init.
COLABDESIGN_SPEC = "git+https://github.com/sokrypton/ColabDesign.git@v1.1.1"


class AF2BIND(BaseConfig):
    """
    AF2BIND: predict ligand-binding residues from a protein structure.

    Inputs:
        structures: input PDB structures (one prediction per structure).
        chain: target chain to score (default "A").
        mask_sidechains: hide target side-chain atoms from AF2 (default True;
            selects AF2BIND's side-chain-masked weight set, the notebook
            default and the more transferable setting).
        mask_sequence: also hide the target sequence identity (default False).
        top_k: number of highest-p_bind residues reported as a chain-aware
            selection string in the summary table (default 15).

    Outputs:
        Streams:
            binding:  per-residue resi-csv (one <id>.csv per input) with columns
                      id | chain | resi | resn | p_bind. Consumable by the Selection
                      tool to retrieve predicted binding residues, e.g.
                      ``Selection.add(af2bind.streams.binding, include="p_bind>0.5")``.
        Tables:
            binding:  id | chain | resi | resn | p_bind
                      (one row per residue of the target chain)
            summary:  id | n_residues | top_resi | top_p_bind | binding_residues
                      (binding_residues = top_k residues by p_bind in
                      chain-aware format, e.g. "A45+A78-80", usable directly
                      as a downstream selection column)
            missing:  id | cause
    """

    TOOL_NAME = "AF2BIND"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        # AF2 network params live in the SHARED cache (its params/ subdir holds
        # params_model_*.npz). The AF2BIND linear head is tool-specific and
        # stays under the AF2BIND folder.
        af2_params_root = folders.get("AlphaFoldParams", "")
        params_dir = f"{af2_params_root}/params"
        head_root = folders.get("AF2BIND", "")
        head_dir = f"{head_root}/af2bind_params"
        # Verification: a representative AF2 ptm weight and the unzipped head exist.
        weights_check = (
            f'[ -f "{params_dir}/params_model_1_ptm.npz" ] && '
            f'[ -d "{head_dir}/attempt_7_2k_lam0-03" ]'
        )
        # AF2 params: shared cache, download-if-absent (cluster reuses
        # LocalColabFold's existing params and downloads nothing).
        af2_params_block = ensure_weights_block(
            dest_dir=params_dir,
            sentinel=f"{params_dir}/params_model_1_ptm.npz",
            download_cmds=f'curl -fsSL "{AF2_PARAMS_URL}" | tar x -C "{params_dir}"',
            label="AlphaFold2 params",
        )
        download_block = f"""{af2_params_block}
# AF2BIND linear-head weights are tool-specific (not shared).
mkdir -p "{head_dir}"
if [ ! -d "{head_dir}/attempt_7_2k_lam0-03" ]; then
    echo "Downloading AF2BIND head weights"
    curl -fsSL -o "{head_root}/attempt_7_2k_lam0-03.zip" "{AF2BIND_WEIGHTS_URL}"
    unzip -o -q "{head_root}/attempt_7_2k_lam0-03.zip" -d "{head_dir}"
    rm -f "{head_root}/attempt_7_2k_lam0-03.zip"
fi
"""

        # Both cluster and Colab use a dedicated af2bind env pinned to jax 0.4.23
        # (the last jax with jax.linear_util, which ColabDesign v1.1.1's haiku
        # imports; jax>=0.4.24 removed it). Colab's base Python now ships jax
        # 0.7.x, so the old "reuse base JAX" Colab path no longer works — the
        # dedicated env is the only viable route on both platforms. The conda
        # yaml is the python+numpy base; the JAX/ColabDesign stack is
        # pip-installed below as an explicit step. The jax[cuda11_pip] wheels
        # bundle their own CUDA 11 userspace, which runs on the T4 driver.
        biopipelines = folders.get("biopipelines", "")
        # MPLBACKEND=Agg: ColabDesign imports matplotlib, and Colab's inherited
        # MPLBACKEND=module://matplotlib_inline.backend_inline crashes the import
        # in this isolated env (no matplotlib_inline here).
        cd_check = f'MPLBACKEND=Agg {env_manager} run -n af2bind python -c "import colabdesign" >/dev/null 2>&1'
        skip = "" if force_reinstall else f"""# Check if already installed
if {cd_check} && {weights_check}; then
    echo "AF2BIND already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block("af2bind", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("af2bind", env_manager, biopipelines)
        # Explicit pip install (not a yaml pip: section): mamba env-create can
        # swallow a failing pip resolve and still exit 0, leaving jax/colabdesign
        # missing. Running it here surfaces failures and lets `set -e` abort.
        # --find-links provides the cuda11 jaxlib wheel for jax[cuda11_pip].
        # scipy<1.13: jax 0.4.23 imports scipy.linalg.tril, removed in scipy 1.13.
        pip_block = (
            f'{env_manager} run -n af2bind pip install '
            f'--find-links https://storage.googleapis.com/jax-releases/jax_cuda_releases.html '
            f'"jax[cuda11_pip]==0.4.23" "scipy<1.13" "dm-haiku==0.0.10" "chex==0.1.7" '
            f'"optax==0.1.7" "{COLABDESIGN_SPEC}"'
        )
        return f"""echo "=== Installing AF2BIND ==="
{skip}{remove_block}
{env_block}
echo "Installing JAX + ColabDesign into af2bind env"
{pip_block}
{download_block}
# Verify installation
if {cd_check} && {weights_check}; then
    touch "$INSTALL_SUCCESS"
    echo "=== AF2BIND installation complete ==="
else
    echo "ERROR: AF2BIND verification failed (colabdesign import or weights missing)"
    exit 1
fi
"""

    structures_json = Path(lambda self: self.configuration_path("structures.json"))
    binding_csv = Path(lambda self: self.table_path("binding"))
    binding_map = Path(lambda self: self.stream_map_path("binding"))
    summary_csv = Path(lambda self: self.table_path("summary"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_af2bind.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 chain: str = "A",
                 mask_sidechains: bool = True,
                 mask_sequence: bool = False,
                 top_k: int = 15,
                 **kwargs):
        self.structures = structures
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")
        self.chain = chain
        self.mask_sidechains = bool(mask_sidechains)
        self.mask_sequence = bool(mask_sequence)
        self.top_k = int(top_k)
        super().__init__(**kwargs)

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        if not self.chain:
            raise ValueError("chain parameter is required")
        if self.top_k < 1:
            raise ValueError("top_k must be >= 1")
        _validate_freeform_string("chain", self.chain)

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"STRUCTURES: {len(self.structures_stream)} files")
        lines.append(f"TARGET CHAIN: {self.chain}")
        lines.append(f"MASK SIDECHAINS: {self.mask_sidechains}")
        lines.append(f"MASK SEQUENCE: {self.mask_sequence}")
        lines.append(f"TOP K: {self.top_k}")
        return lines

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        # AF2 network params come from the shared cache (its params/ subdir);
        # the AF2BIND linear head lives under the tool's own folder. Forward-slash
        # join: this path runs in POSIX bash, so os.path.join would wrongly use
        # "\" when the pipeline is configured on Windows.
        af2_params_root = self.folders["AlphaFoldParams"]
        head_dir = f"{self.folders['AF2BIND']}/af2bind_params"

        sc_arg = " --mask-sidechains" if self.mask_sidechains else ""
        seq_arg = " --mask-sequence" if self.mask_sequence else ""
        upstream_missing_path = self._get_upstream_missing_table_path(self.structures)
        upstream_missing_flag = f' \\\n    --upstream-missing "{upstream_missing_path}"' if upstream_missing_path else ""
        script = "#!/bin/bash\n"
        script += "# AF2BIND binding-site prediction script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        # ColabDesign imports matplotlib; Colab's inherited MPLBACKEND crashes it.
        script += f"""export MPLBACKEND=Agg
echo "Running AF2BIND on {len(self.structures_stream)} structure(s) (chain={self.chain})"
python "{self.helper_py}" \\
    --structures-json "{self.structures_json}" \\
    --af2-params-dir "{af2_params_root}" \\
    --head-dir "{head_dir}" \\
    --chain "{self.chain}" \\
    --top-k {self.top_k} \\
    --binding-csv "{self.binding_csv}" \\
    --binding-dir "{self.stream_folder('binding')}" \\
    --binding-map-csv "{self.binding_map}" \\
    --summary-csv "{self.summary_csv}" \\
    --missing-csv "{self.missing_csv}"{sc_arg}{seq_arg}{upstream_missing_flag}
"""
        script += self.generate_completion_check_footer()
        return script

    def get_output_files(self) -> Dict[str, Any]:
        binding_stream = DataStream(
            name="binding",
            ids=self.structures_stream.ids,
            files=[self.stream_path("binding", "<id>.csv")],
            map_table=self.binding_map,
            format="resi-csv",
        )
        tables = {
            "binding": TableInfo(
                name="binding",
                path=self.binding_csv,
                columns=["id", "chain", "resi", "resn", "p_bind"],
                description="AF2BIND per-residue binding probability",
            ),
            "summary": TableInfo(
                name="summary",
                path=self.summary_csv,
                columns=["id", "n_residues", "top_resi", "top_p_bind", "binding_residues"],
                description="AF2BIND per-structure summary with top-k binding-residue selection",
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "kind", "cause"],
                description="IDs removed (upstream or local failure) with removal reason",
            ),
        }
        return {
            "binding": binding_stream,
            "tables": tables,
            "output_folder": self.output_folder,
        }
