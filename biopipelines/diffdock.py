# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
DiffDock tool for diffusion-based protein-ligand docking.

DiffDock samples K candidate poses per protein-ligand pair via a learned
diffusion process over translation, rotation, and torsion, then re-ranks them
with a confidence model. It is GPU-bound (CPU fallback exists but is much
slower) and works directly from PDB + SMILES/SDF without binding-box hints.

Reference:
    Corso et al. (2023) DiffDock: Diffusion Steps, Twists, and Turns for
    Molecular Docking. ICLR 2023. https://github.com/gcorso/DiffDock
"""

import os
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .combinatorics import predict_output_ids_with_provenance, generate_multiplied_ids_pattern
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from combinatorics import predict_output_ids_with_provenance, generate_multiplied_ids_pattern


class DiffDock(BaseConfig):
    """
    DiffDock: diffusion-based molecular docking with confidence re-ranking.

    For each (protein, ligand) pair, samples `samples_per_complex` poses and
    re-ranks them with the confidence model. Output IDs follow the multi-axis
    pattern <protein>+<ligand>_rank<N>; the rank-1 pose lives under the
    `structures` stream and the full ranked confidence is reported in a table.

    Inputs:
        structures:           protein backbones (PDB stream)
        compounds:            ligands as SMILES or SDF (compounds stream)
        samples_per_complex:  number of poses sampled per pair (default 10)
        inference_steps:      reverse-diffusion steps (default 20)
        actual_steps:         steps actually executed (<= inference_steps; None defaults to inference_steps - 1, the common DiffDock setting)
        no_final_step_noise:  disable noise on the last step (default True)
        batch_size:           DiffDock internal batch size (default 32, matching upstream)

    Outputs:
        Streams:
            structures:  rank-N pose SDF per (protein, ligand, rank) tuple
        Tables:
            confidence:  id | structures.id | compounds.id | rank | confidence
            missing:     id | removed_by | cause

    Usage::

        with Pipeline(project="Examples", job="DiffDock-demo"):
            Resources(gpu="A100", memory="32GB", time="4:00:00")
            target = PDB("4ufc", convert="pdb")
            lig = Ligand(smiles="CC(=O)Oc1ccccc1C(=O)O", ids="ASA")
            dock = DiffDock(structures=target, compounds=lig)
    """

    TOOL_NAME = "DiffDock"
    TOOL_VERSION = "1.0"

    # ------------------------------------------------------------------
    # Install
    # ------------------------------------------------------------------

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        """Clone gcorso/DiffDock + facebookresearch/esm, bring up its conda
        env, install PyG wheels against the active torch version.

        Recipe mirrors the canonical Colab notebook
        (https://colab.research.google.com/github/hgbrian/biocolabs/blob/master/DiffDock.ipynb).
        Pinned to DiffDock commit a6c5275 and esm commit ca8a710 per the
        notebook. openfold is intentionally NOT installed — ESM embeddings
        come from the facebookresearch/esm side-clone with PYTHONPATH pointed
        at it (set in generate_script when invoking inference).

        Verification: (a) repo cloned, (b) esm side-clone present,
        (c) torch + torch_geometric importable in the env.
        """
        biopipelines = folders.get("biopipelines", "")
        repo_dir = folders.get("DiffDock", "")
        parent_dir = os.path.dirname(repo_dir)
        esm_dir = os.path.join(repo_dir, "esm")  # side-clone, used via PYTHONPATH

        env_check = cls._env_exists_check("diffdock", env_manager)
        repo_check = f'[ -f "{repo_dir}/inference.py" ]'
        esm_check = f'[ -f "{esm_dir}/setup.py" ]'
        skip = "" if force_reinstall else f"""# Check if already installed
if {repo_check} && {esm_check} && {env_check} \\
   && {env_manager} run -n diffdock python -c "import torch, torch_geometric" >/dev/null 2>&1; then
    echo "DiffDock already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        # Clone DiffDock at the commit pinned by the canonical Colab notebook,
        # then clone facebookresearch/esm side-by-side and pin its commit too.
        clone_block = f"""mkdir -p {parent_dir}
cd {parent_dir}
if [ ! -d "{repo_dir}" ]; then
    git clone https://github.com/gcorso/DiffDock.git "{repo_dir}"
    cd "{repo_dir}"
    git checkout a6c5275
fi
if [ ! -d "{esm_dir}" ]; then
    git clone https://github.com/facebookresearch/esm "{esm_dir}"
    cd "{esm_dir}"
    git checkout ca8a710
fi"""

        remove_block = cls._env_remove_block("diffdock", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("diffdock", env_manager, biopipelines)

        # PyG wheels are pinned against whatever torch the env ends up with.
        # Following the notebook recipe: uninstall any prior PyG bits (no-op
        # on a fresh env), then reinstall against the active torch from the
        # data.pyg.org wheel index. `pyg` (the top-level package) is its own
        # PyPI dist; the notebook installs torch-geometric from git for the
        # newest features.
        pyg_block = f"""# Install PyG wheels against the active torch in the env.
TORCH_VER=$({env_manager} run -n diffdock python -c "import torch; print(torch.__version__)")
PYG_URL="https://data.pyg.org/whl/torch-${{TORCH_VER}}.html"
echo "PyG wheel index: $PYG_URL"
{env_manager} run -n diffdock pip uninstall -y torch-scatter torch-sparse torch-geometric torch-cluster torch-spline-conv >/dev/null 2>&1 || true
{env_manager} run -n diffdock pip install --quiet torch-scatter -f "$PYG_URL"
{env_manager} run -n diffdock pip install --quiet torch-sparse -f "$PYG_URL"
{env_manager} run -n diffdock pip install --quiet torch-cluster -f "$PYG_URL"
{env_manager} run -n diffdock pip install --quiet git+https://github.com/pyg-team/pytorch_geometric.git
"""

        return f"""echo "=== Installing DiffDock ==="
{skip}{clone_block}

{remove_block}
{env_block}
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create diffdock environment."
    exit 1
fi

{pyg_block}

# Verify installation
if {repo_check} && {esm_check} \\
   && {env_manager} run -n diffdock python -c "import torch, torch_geometric" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== DiffDock installation complete ==="
else
    echo "ERROR: DiffDock verification failed (repo, esm clone, or imports missing)"
    exit 1
fi
"""

    # ------------------------------------------------------------------
    # Lazy path descriptors — canonical sub-layout
    #   _configuration/  input_csv (protein_ligand_csv for DiffDock),
    #                    .input_structures.json, .input_compounds.json,
    #                    default_inference_args.yaml (per-run override)
    #   _execution/      raw DiffDock dump (<complex_name>/ folders)
    #   structures/      reorganised rank<N>.sdf per (prot+lig, rank) +
    #                    structures_map.csv
    #   tables/          confidence + missing
    # ------------------------------------------------------------------

    structures_json = Path(lambda self: self.configuration_path(".input_structures.json"))
    compounds_json = Path(lambda self: self.configuration_path(".input_compounds.json"))
    input_csv = Path(lambda self: self.configuration_path("protein_ligand.csv"))
    inference_args_yaml = Path(lambda self: self.configuration_path("inference_args.yaml"))
    raw_out_folder = Path(lambda self: self.execution_folder)
    structures_map = Path(lambda self: self.stream_map_path("structures"))
    confidence_csv = Path(lambda self: self.table_path("confidence"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    build_csv_py = Path(lambda self: self.pipe_script_path("pipe_diffdock_build_csv.py"))
    postprocess_py = Path(lambda self: self.pipe_script_path("pipe_diffdock_postprocess.py"))
    diffdock_inference_py = Path(lambda self: os.path.join(self.folders["DiffDock"], "inference.py"))

    # ------------------------------------------------------------------
    # Constructor
    # ------------------------------------------------------------------

    def __init__(
        self,
        structures: Union[DataStream, StandardizedOutput],
        compounds: Union[DataStream, StandardizedOutput],
        samples_per_complex: int = 10,
        inference_steps: int = 20,
        actual_steps: Optional[int] = None,
        no_final_step_noise: bool = True,
        batch_size: int = 32,
        **kwargs,
    ):
        # Keep the original inputs (StandardizedOutput or DataStream) for
        # predict_output_ids_with_provenance, which plucks the named stream
        # itself. Don't replace with the resolved streams.
        self.structures_input = structures
        self.compounds_input = compounds

        # Resolve protein input
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(
                f"structures must be DataStream or StandardizedOutput, got {type(structures).__name__}"
            )

        # Resolve compound input
        if isinstance(compounds, StandardizedOutput):
            self.compounds_stream: DataStream = compounds.streams.compounds
        elif isinstance(compounds, DataStream):
            self.compounds_stream = compounds
        else:
            raise ValueError(
                f"compounds must be DataStream or StandardizedOutput, got {type(compounds).__name__}"
            )

        self.samples_per_complex = samples_per_complex
        self.inference_steps = inference_steps
        # None defaults to inference_steps - 1, the common DiffDock setting.
        self.actual_steps = actual_steps if actual_steps is not None else inference_steps - 1
        self.no_final_step_noise = no_final_step_noise
        self.batch_size = batch_size
        super().__init__(**kwargs)

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        if not self.structures_stream.has_only_formats("pdb"):
            raise ValueError(
                f"DiffDock requires PDB-format structures, got '{self.structures_stream.format}'. "
                "Use convert='pdb' upstream (e.g. PDB(..., convert='pdb'))."
            )
        if not self.compounds_stream or len(self.compounds_stream) == 0:
            raise ValueError("compounds parameter is required and must not be empty")

        if self.samples_per_complex < 1:
            raise ValueError("samples_per_complex must be at least 1")
        if self.inference_steps < 1:
            raise ValueError("inference_steps must be at least 1")
        if not (1 <= self.actual_steps <= self.inference_steps):
            raise ValueError(
                f"actual_steps must satisfy 1 <= actual_steps <= inference_steps "
                f"(got {self.actual_steps} vs {self.inference_steps})"
            )
        if self.batch_size < 1:
            raise ValueError("batch_size must be at least 1")

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
            f"PROTEINS: {len(self.structures_stream)}",
            f"LIGANDS: {len(self.compounds_stream)}",
            f"PAIRS: {len(self.structures_stream) * len(self.compounds_stream)}",
            f"SAMPLES PER PAIR: {self.samples_per_complex}",
            f"INFERENCE STEPS: {self.inference_steps} (actual {self.actual_steps})",
            f"BATCH SIZE: {self.batch_size}",
        ])
        return lines

    # ------------------------------------------------------------------
    # Script generation
    # ------------------------------------------------------------------

    def generate_script(self, script_path: str) -> str:
        # Serialize DataStreams — pipe scripts iterate them at runtime to build
        # the protein_ligand.csv DiffDock consumes.
        self.structures_stream.save_json(self.structures_json)
        self.compounds_stream.save_json(self.compounds_json)

        script = "#!/bin/bash\n"
        script += "# DiffDock execution script\n"
        script += self.generate_completion_check_header()

        # Phase 1: build the protein_ligand.csv (biopipelines env: needs pandas
        # + biopipelines_io; lighter than the diffdock env).
        script += "# --- Build protein-ligand CSV (biopipelines env) ---\n"
        script += self.activate_environment(name="biopipelines")
        script += f"""python "{self.build_csv_py}" \\
    --structures-json "{self.structures_json}" \\
    --compounds-json "{self.compounds_json}" \\
    --output-csv "{self.input_csv}"

if [ $? -ne 0 ]; then
    echo "Error: building DiffDock input CSV failed"
    exit 1
fi

"""

        # Phase 2: inference under the diffdock env.
        script += "# --- DiffDock inference (diffdock env) ---\n"
        script += self.activate_environment()  # diffdock (primary)

        # At commit a6c5275 DiffDock has no default_inference_args.yaml; the
        # canonical hgbrian/biocolabs notebook just passes the inference knobs
        # via explicit flags. It also runs two preparation steps that this
        # version of DiffDock requires:
        #   1) datasets/esm_embedding_preparation.py — writes a FASTA derived
        #      from the protein_ligand csv.
        #   2) esm/scripts/extract.py — runs ESM2 over that FASTA and writes
        #      per-residue embeddings to data/esm2_output/, which the
        #      PDBBind dataset loader inside inference.py reads.
        # PYTHONPATH + HOME setup mirrors the notebook: facebookresearch/esm
        # lives as a side-clone under DiffDock/esm/, used as a
        # source-on-PYTHONPATH dependency (no openfold needed). HOME points
        # at esm's model_weights cache so first call doesn't re-download the
        # ESM2 checkpoint.
        no_noise_flag = " --no_final_step_noise" if self.no_final_step_noise else ""
        script += f"""mkdir -p "{self.raw_out_folder}"

cd "{self.folders["DiffDock"]}"
export PYTHONPATH="{self.folders["DiffDock"]}/esm:$PYTHONPATH"
# Cache the original HOME so we can restore it before postprocess (otherwise
# the next `micromamba activate` looks under ESM's weight cache and fails).
_BP_DIFFDOCK_ORIG_HOME="$HOME"
export HOME="{self.folders["DiffDock"]}/esm/model_weights"
mkdir -p "$HOME" data

# ESM embedding prep (FASTA from protein_ligand csv)
python datasets/esm_embedding_preparation.py \\
    --protein_ligand_csv "{self.input_csv}" \\
    --out_file data/prepared_for_esm.fasta
if [ $? -ne 0 ]; then
    echo "Error: ESM FASTA preparation failed"
    exit 1
fi

# ESM2 embedding extraction (writes data/esm2_output/<id>.pt files)
python esm/scripts/extract.py esm2_t33_650M_UR50D \\
    data/prepared_for_esm.fasta data/esm2_output \\
    --repr_layers 33 --include per_tok --truncation_seq_length 30000
if [ $? -ne 0 ]; then
    echo "Error: ESM2 embedding extraction failed"
    exit 1
fi

# DiffDock inference proper
python -m inference \\
    --protein_ligand_csv "{self.input_csv}" \\
    --out_dir "{self.raw_out_folder}" \\
    --samples_per_complex {self.samples_per_complex} \\
    --inference_steps {self.inference_steps} \\
    --actual_steps {self.actual_steps} \\
    --batch_size {self.batch_size}{no_noise_flag}

if [ $? -ne 0 ]; then
    echo "Error: DiffDock inference failed"
    exit 1
fi

# Restore HOME so subsequent env activations look in the right place.
export HOME="$_BP_DIFFDOCK_ORIG_HOME"
unset _BP_DIFFDOCK_ORIG_HOME

"""

        # Phase 3: post-processing under biopipelines env (parses confidence
        # from filenames, populates structures stream + tables).
        script += "# --- Post-processing (biopipelines env) ---\n"
        script += self.activate_environment(name="biopipelines")
        script += f"""python "{self.postprocess_py}" \\
    --raw-out "{self.raw_out_folder}" \\
    --input-csv "{self.input_csv}" \\
    --structures-folder "{self.stream_folder("structures")}" \\
    --structures-map "{self.structures_map}" \\
    --confidence-csv "{self.confidence_csv}" \\
    --missing-csv "{self.missing_csv}" \\
    --samples-per-complex {self.samples_per_complex}

if [ $? -ne 0 ]; then
    echo "Error: DiffDock post-processing failed"
    exit 1
fi

"""
        script += self.generate_completion_check_footer()
        return script

    # ------------------------------------------------------------------
    # Output prediction
    # ------------------------------------------------------------------

    def get_output_files(self) -> Dict[str, Any]:
        # Pass the original inputs (StandardizedOutput / DataStream), not the
        # resolved streams — the helper plucks the named stream from each
        # input itself.
        pair_ids, _ = predict_output_ids_with_provenance(
            structures=(self.structures_input, "structures"),
            compounds=(self.compounds_input, "compounds"),
        )

        ranked_ids: List[str] = []
        for pair_id in pair_ids:
            for n in range(1, self.samples_per_complex + 1):
                ranked_ids.append(f"{pair_id}_rank{n}")

        structures = DataStream(
            name="structures",
            ids=ranked_ids,
            files=[self.stream_path("structures", "<id>.sdf")],
            map_table=self.structures_map,
            format="sdf",
        )

        tables = {
            "confidence": TableInfo(
                name="confidence",
                path=self.confidence_csv,
                columns=["id", "structures.id", "compounds.id", "rank", "confidence"],
                description="DiffDock per-pose confidence scores",
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "kind", "cause"],
                description="Pairs that DiffDock could not dock, with reason",
            ),
        }

        return {
            "structures": structures,
            "tables": tables,
            "output_folder": self.output_folder,
        }
