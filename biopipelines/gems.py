# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""GEMS tool: protein-ligand binding affinity prediction via graph neural
network with protein/ligand language-model embeddings.

For each id, GEMS expects matched <id>.pdb (protein) and <id>.sdf (ligand)
files in a single directory. The wrapper stages the paired files into a
scratch directory, runs the upstream two-step workflow (dataprep -> inference),
and parses the predictions CSV into our standard tables/affinity.csv.

Reference: Holcomb et al. (2024) GEMS. https://github.com/camlab-ethz/GEMS
"""

import os
from typing import Dict, List, Any, Optional, Union

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


class GEMS(BaseConfig):
    """
    GEMS: GNN + protein/ligand language-model embeddings for binding affinity.

    Inputs:
        structures: per-id PDB protein structures.
        ligands:    a ligand-producing tool's StandardizedOutput. Coordinate
                    files are taken from its ``streams.structures`` and SMILES
                    metadata (for bond-order templating) from its
                    ``streams.compounds``. Mutually exclusive with
                    ``ligands_3d`` / ``ligands_smiles``.
        ligands_3d: a DataStream of per-id ligand coordinate files (sdf/pdb/
                    mol2), passed directly instead of via ``ligand``.
        ligands_smiles: an optional DataStream whose map_table carries a
                    ``smiles`` column for bond-order templating, paired with
                    ``ligands_3d``.
        skip_ligand_embedding: if True, use the GEMS18e variant with no
                               ChemBERTa ligand embedding (faster, less accurate).

    Pass either ``ligands`` (a StandardizedOutput) or ``ligands_3d``
    (+ optional ``ligands_smiles``) — not both.

    Outputs:
        Streams: (none)
        Tables:
            affinity: id | structures.id | ligands.id | pkd_pred
            missing:  id | cause
    """

    TOOL_NAME = "GEMS"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        biopipelines = folders.get("biopipelines", "")
        repo_dir = folders.get("GEMS", "")
        # Verify torch + ankh actually import (catches the MKL ABI breakage)
        # rather than just probing that the env directory exists.
        import_check = (f'{env_manager} run -n gems python -c '
                        f'"import torch, ankh, torch_geometric" >/dev/null 2>&1')
        repo_check = f'[ -d "{repo_dir}/model" ]'
        skip = "" if force_reinstall else f"""# Check if already installed
if {import_check} && {repo_check}; then
    echo "GEMS already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        clone_block = f"""mkdir -p "{os.path.dirname(repo_dir)}"
if [ ! -d "{repo_dir}" ]; then
    git clone https://github.com/camlab-ethz/GEMS.git "{repo_dir}"
fi"""
        # Upstream inference.py only calls model.eval() in the GPU (else) branch
        # of load_model_state; on a CPU device the models stay in train mode and
        # BatchNorm rejects a single-pair batch ("Expected more than 1 value per
        # channel"). Move the eval() out of the else so it always runs. Anchor on
        # the over-indented (8-space) eval line that sits inside the else block.
        patch_eval = (
            '''sed -i 's|^        model.eval()  # Set the model to evaluation mode|'''
            '''    model.eval()  # Set the model to evaluation mode (CPU too)|' '''
            f'"{repo_dir}/inference.py"'
        )
        remove_block = cls._env_remove_block("gems", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("gems", env_manager, biopipelines)
        return f"""echo "=== Installing GEMS ==="
{skip}{clone_block}

# Always eval() the model (upstream skips it on CPU -> BatchNorm train-mode crash).
{patch_eval}

{remove_block}
{env_block}

if {import_check} && [ -f "{repo_dir}/inference.py" ] && [ -d "{repo_dir}/model" ]; then
    touch "$INSTALL_SUCCESS"
    echo "=== GEMS installation complete ==="
else
    echo "ERROR: GEMS verification failed (env import, inference.py, or model weights missing)"
    exit 1
fi
"""

    structures_json = Path(lambda self: self.configuration_path("structures.json"))
    compounds_json = Path(lambda self: self.configuration_path("compounds.json"))
    smiles_json = Path(lambda self: self.configuration_path("smiles.json"))
    affinity_csv = Path(lambda self: self.table_path("affinity"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_gems.py"))
    gems_repo = Path(lambda self: self.folders["GEMS"])

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 ligands: Optional[StandardizedOutput] = None,
                 ligands_3d: Optional[DataStream] = None,
                 ligands_smiles: Optional[DataStream] = None,
                 skip_ligand_embedding: bool = False,
                 **kwargs):
        self.structures = structures
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        # Ligand inputs — two mutually exclusive forms (see validate_params):
        #   * ligands: a StandardizedOutput. Coordinate files come from its
        #     streams.structures; SMILES metadata (for bond-order templating of
        #     a coordinate-only PDB ligand) from its streams.compounds.
        #   * ligands_3d (+ optional ligands_smiles): the two streams directly.
        self._ligands_arg = ligands
        self._ligands_3d_arg = ligands_3d
        if ligands is not None and isinstance(ligands, StandardizedOutput):
            self.ligands_stream = getattr(ligands.streams, "structures", None)
            self.smiles_stream = getattr(ligands.streams, "compounds", None)
        else:
            self.ligands_stream = ligands_3d
            self.smiles_stream = ligands_smiles

        self.skip_ligand_embedding = bool(skip_ligand_embedding)
        super().__init__(**kwargs)

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")

        # Exactly one ligand input form: `ligands` (StandardizedOutput) XOR `ligands_3d`.
        if self._ligands_arg is not None and self._ligands_3d_arg is not None:
            raise ValueError("provide either `ligands` or `ligands_3d`, not both")
        if self._ligands_arg is None and self._ligands_3d_arg is None:
            raise ValueError("provide either `ligands` (a StandardizedOutput) or `ligands_3d` (a DataStream)")
        if self._ligands_arg is not None and not isinstance(self._ligands_arg, StandardizedOutput):
            raise ValueError(f"ligands must be a StandardizedOutput, got {type(self._ligands_arg).__name__}")
        if self._ligands_3d_arg is not None and not isinstance(self._ligands_3d_arg, DataStream):
            raise ValueError(f"ligands_3d must be a DataStream, got {type(self._ligands_3d_arg).__name__}")

        # Coordinate stream must actually carry files — a value-based compounds
        # stream (files=[], e.g. a raw Ligand) is not coordinates.
        files = getattr(self.ligands_stream, "files", None)
        has_files = bool(files if not isinstance(files, list) else any(files))
        if not has_files:
            raise ValueError(
                "ligand coordinate stream has no files (a value-based compounds "
                "stream is not coordinates). Run OpenBabel(compounds=..., "
                "convert_3d='sdf') first, or pass its output as `ligands`."
            )
        if len(self.ligands_stream) == 0:
            raise ValueError("ligands input must not be empty")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"STRUCTURES: {len(self.structures_stream)} proteins")
        lines.append(f"LIGANDS: {len(self.ligands_stream)} ligands")
        lines.append(f"SKIP LIGAND EMBEDDING: {self.skip_ligand_embedding}")
        return lines

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        self.ligands_stream.save_json(self.compounds_json)
        smiles_arg = ""
        if self.smiles_stream is not None:
            self.smiles_stream.save_json(self.smiles_json)
            smiles_arg = f' --smiles-json "{self.smiles_json}"'
        skip_arg = " --skip-ligand-embedding" if self.skip_ligand_embedding else ""
        upstream_missing_flag = self.upstream_missing_flag(self.structures, self._ligands_arg)
        script = "#!/bin/bash\n"
        script += "# GEMS protein-ligand affinity prediction\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += f"""echo "Predicting affinity for {len(self.structures_stream)} x {len(self.ligands_stream)} pair(s)"
python "{self.helper_py}" \\
    --structures-json "{self.structures_json}" \\
    --compounds-json "{self.compounds_json}"{smiles_arg} \\
    --gems-repo "{self.gems_repo}" \\
    --scratch-dir "{self.extras_path()}" \\
    --affinity-csv "{self.affinity_csv}" \\
    --missing-csv "{self.missing_csv}"{skip_arg}{upstream_missing_flag}
"""
        script += self.generate_completion_check_footer()
        return script

    def get_output_files(self) -> Dict[str, Any]:
        tables = {
            "affinity": TableInfo(
                name="affinity",
                path=self.affinity_csv,
                columns=["id", "structures.id", "ligands.id", "pkd_pred"],
                description="GEMS predicted pKd (higher = stronger binder)",
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "kind", "cause"],
                description="IDs removed (upstream or local failure) with removal reason",
            ),
        }
        return {"tables": tables, "output_folder": self.output_folder}
