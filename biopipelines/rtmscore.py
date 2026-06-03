# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""RTMScore tool: residue-atom distance-likelihood scoring with graph
transformers, for protein-ligand pose ranking and affinity-style ranking.

The wrapper takes a structures stream (complex PDBs) and a compounds stream
(ligand SDF files); for each id it runs the upstream rtmscore.py script
with `-gen_pocket` against the ligand SDF as the reference, producing one
score per ligand pose.

Reference: Shen et al. (2022) J. Med. Chem. https://github.com/sc8668/RTMScore
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


class RTMScore(BaseConfig):
    """
    RTMScore: pose / affinity-style ranking via residue-atom distance likelihood.

    Inputs:
        structures: complex PDB structures used as the protein source. The
                    pocket is auto-extracted within `cutoff` Å of the
                    ligand SDF reference.
        ligands:    a ligand-producing tool's StandardizedOutput. The ligand
                    SDF files (the references to be scored) are taken from its
                    ``streams.structures`` and SMILES metadata (for bond-order
                    templating) from its ``streams.compounds``. Mutually
                    exclusive with ``ligands_3d`` / ``ligands_smiles``.
        ligands_3d: a DataStream of per-id ligand coordinate files (sdf/pdb/
                    mol2), passed directly instead of via ``ligands``. Each
                    file may contain multiple poses; one row per pose is
                    written to the output table.
        ligands_smiles: an optional DataStream whose map_table carries a
                    ``smiles`` column for bond-order templating, paired with
                    ``ligands_3d``.
        cutoff:     pocket cutoff in Å around the reference ligand (default 10).
        model:      one of model1..model3 (default model1). Different
                    checkpoints shipped with the upstream repo.

    Pass either ``ligands`` (a StandardizedOutput) or ``ligands_3d``
    (+ optional ``ligands_smiles``) — not both.

    Outputs:
        Streams: (none)
        Tables:
            scores:  id | structures.id | ligands.id | pose | rtmscore
            missing: id | removed_by | kind | cause
    """

    TOOL_NAME = "RTMScore"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        biopipelines = folders.get("biopipelines", "")
        repo_dir = folders.get("RTMScore", "")
        # Verify the env can actually import the heavy deps (dgl/torch), not
        # just that the env directory exists — a half-built or ABI-broken env
        # must not be mistaken for a working install.
        import_check = (f'{env_manager} run -n rtmscore python -c '
                        f'"import torch, dgl, rdkit, torch_scatter" >/dev/null 2>&1')
        repo_check = f'[ -d "{repo_dir}/RTMScore" ]'
        model_check = f'[ -f "{repo_dir}/trained_models/rtmscore_model1.pth" ]'
        skip = "" if force_reinstall else f"""# Check if already installed
if {import_check} && {repo_check} && {model_check}; then
    echo "RTMScore already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        clone_block = f"""mkdir -p "{os.path.dirname(repo_dir)}"
if [ ! -d "{repo_dir}" ]; then
    git clone https://github.com/sc8668/RTMScore.git "{repo_dir}"
fi"""
        # The upstream example/rtmscore.py hardcodes the author's BABEL_LIBDIR
        # (os.environ["BABEL_LIBDIR"] = "/home/shenchao/.conda/.../openbabel/..."),
        # which clobbers any value we set and makes OpenBabel's plugin loader
        # fail on every other machine. Neutralize that line so the env var we
        # export at runtime (resolved from the env's own openbabel dir) wins.
        patch_block = (
            f'''sed -i 's|^os.environ\\["BABEL_LIBDIR"\\].*|'''
            f'''os.environ.setdefault("BABEL_LIBDIR", os.environ.get("BABEL_LIBDIR", ""))|' '''
            f'"{repo_dir}/example/rtmscore.py"'
        )
        # Upstream model2.py builds the batch-index tensor on CPU
        # (`th.tensor(range(B))`) then indexes it with C_mask, which is on the
        # model's device. On GPU that mismatches ("indices should be ... on the
        # same device"). Pin C_batch to C_mask.device so RTMScore runs on GPU.
        device_patch = (
            f'''sed -i 's|C_batch = th.tensor(range(B)).unsqueeze|'''
            f'''C_batch = th.tensor(range(B), device=C_mask.device).unsqueeze|' '''
            f'"{repo_dir}/RTMScore/model/model2.py"'
        )
        remove_block = cls._env_remove_block("rtmscore", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("rtmscore", env_manager, biopipelines)
        # torch_scatter comes from the conda yaml (pytorch_scatter pinned to a
        # cuda118 build) so it matches conda-forge pytorch's libtorch ABI — the
        # PyG pip wheel is built against PyPI torch and dies with an
        # undefined-symbol on conda-forge torch. torchdata still needs the pip
        # pin: the conda solver lands torchdata 0.11 (missing `datapipes` ->
        # breaks `import dgl`'s graphbolt), so force 0.9.0 (--no-deps leaves
        # torch untouched; 0.9.0 still ships datapipes).
        pip_block = (
            f'{env_manager} run -n rtmscore pip install --no-deps "torchdata==0.9.0"'
        )
        return f"""echo "=== Installing RTMScore ==="
{skip}{clone_block}

# Remove the hardcoded BABEL_LIBDIR so the runtime-exported one is honored.
{patch_block}
# Put the model's batch-index tensor on the same device as its mask (GPU fix).
{device_patch}

{remove_block}
{env_block}
{pip_block}

if {import_check} && {model_check}; then
    touch "$INSTALL_SUCCESS"
    echo "=== RTMScore installation complete ==="
else
    echo "ERROR: RTMScore verification failed (env or model weights missing)"
    exit 1
fi
"""

    structures_json = Path(lambda self: self.configuration_path("structures.json"))
    compounds_json = Path(lambda self: self.configuration_path("compounds.json"))
    smiles_json = Path(lambda self: self.configuration_path("smiles.json"))
    scores_csv = Path(lambda self: self.table_path("scores"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_rtmscore.py"))
    rtmscore_script = Path(
        lambda self: os.path.join(self.folders["RTMScore"], "example", "rtmscore.py")
    )
    model_path = Path(
        lambda self: os.path.join(self.folders["RTMScore"], "trained_models",
                                  f"rtmscore_{self.model}.pth")
    )

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 ligands: Optional[StandardizedOutput] = None,
                 ligands_3d: Optional[DataStream] = None,
                 ligands_smiles: Optional[DataStream] = None,
                 cutoff: float = 10.0,
                 model: str = "model1",
                 **kwargs):
        self.structures = structures
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        # Ligand inputs — two mutually exclusive forms (see validate_params):
        #   * ligands: a StandardizedOutput. RTMScore needs the per-id ligand
        #     coordinate files, taken from its streams.structures; SMILES
        #     metadata from its streams.compounds supplies bond-order templating
        #     for the pocket extraction (which needs a chemically valid
        #     reference ligand, not a coordinate-only PDB).
        #   * ligands_3d (+ optional ligands_smiles): the two streams directly.
        self._ligands_arg = ligands
        self._ligands_3d_arg = ligands_3d
        if ligands is not None and isinstance(ligands, StandardizedOutput):
            self.ligands_stream = getattr(ligands.streams, "structures", None)
            self.smiles_stream = getattr(ligands.streams, "compounds", None)
        else:
            self.ligands_stream = ligands_3d
            self.smiles_stream = ligands_smiles

        self.cutoff = float(cutoff)
        self.model = model
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
        if self.cutoff <= 0:
            raise ValueError("cutoff must be positive")
        if self.model not in {"model1", "model2", "model3"}:
            raise ValueError("model must be one of model1..model3")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"STRUCTURES: {len(self.structures_stream)} complexes")
        lines.append(f"LIGANDS: {len(self.ligands_stream)} ligand files")
        lines.append(f"CUTOFF: {self.cutoff} A")
        lines.append(f"MODEL: {self.model}")
        return lines

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        self.ligands_stream.save_json(self.compounds_json)
        smiles_arg = ""
        if self.smiles_stream is not None:
            self.smiles_stream.save_json(self.smiles_json)
            smiles_arg = f' --smiles-json "{self.smiles_json}"'
        upstream_missing_flag = self.upstream_missing_flag(self.structures, self._ligands_arg)
        script = "#!/bin/bash\n"
        script += "# RTMScore protein-ligand scoring script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += f"""echo "Scoring with RTMScore ({self.model}, cutoff={self.cutoff})"
python "{self.helper_py}" \\
    --structures-json "{self.structures_json}" \\
    --compounds-json "{self.compounds_json}"{smiles_arg} \\
    --rtmscore-script "{self.rtmscore_script}" \\
    --model-path "{self.model_path}" \\
    --cutoff {self.cutoff} \\
    --scratch-dir "{self.extras_path()}" \\
    --scores-csv "{self.scores_csv}" \\
    --missing-csv "{self.missing_csv}"{upstream_missing_flag}
"""
        script += self.generate_completion_check_footer()
        return script

    def get_output_files(self) -> Dict[str, Any]:
        tables = {
            "scores": TableInfo(
                name="scores",
                path=self.scores_csv,
                columns=["id", "structures.id", "ligands.id", "pose", "rtmscore"],
                description="RTMScore per-pose scores (higher = better)",
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "kind", "cause"],
                description="IDs removed (upstream or local failure) with removal reason",
            ),
        }
        return {"tables": tables, "output_folder": self.output_folder}
