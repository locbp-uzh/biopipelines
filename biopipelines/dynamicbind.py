# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
DynamicBind tool for flexible-backbone protein-ligand docking.

DynamicBind is a deep equivariant generative model that predicts ligand-specific
protein conformations, allowing the receptor backbone to flex toward the bound
state. It uses two conda envs at runtime: `dynamicbind` for inference and
`dynamicbind_relax` for OpenMM-based pose relaxation.

Reference:
    Lu et al. (2024) DynamicBind: predicting ligand-specific protein-ligand
    complex structure with a deep equivariant generative model.
    Nat. Commun. https://github.com/luwei0917/DynamicBind
"""

import os
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .combinatorics import predict_output_ids_with_provenance
    from .biopipelines_io import Resolve
    from .config_manager import ConfigManager
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from combinatorics import predict_output_ids_with_provenance
    from biopipelines_io import Resolve
    from config_manager import ConfigManager


# Zenodo v2 weights (matches upstream README's recommended pin).
DYNAMICBIND_WORKDIR_ZIP_URL = "https://zenodo.org/records/10183369/files/workdir.zip"


class DynamicBind(BaseConfig):
    """
    DynamicBind: deep equivariant flexible-backbone protein-ligand docking.

    For each (protein, ligand) pair: samples K candidate poses with flexed
    backbone, ranks them by predicted lDDT, optionally relaxes with OpenMM.
    Output IDs follow the multi-axis pattern <protein>+<ligand>_rank<N>.

    Inputs:
        structures:           protein PDBs (one entry per scaffold)
        compounds:            ligands as SMILES (compounds stream)
        num_samples:          candidate poses generated per pair (upstream
                              --samples_per_complex; default 40)
        num_saved:            ranked poses kept per pair (upstream
                              --savings_per_complex; None defaults to num_samples).
                              Must be <= num_samples.
        inference_steps:      reverse-diffusion steps (default 20)
        num_workers:          parallel relax workers (default 1)
        rigid_protein:        disable backbone flexing (default False — flex on)
        make_movie:           render an mp4/gif of each kept pose's reverse-
                              diffusion trajectory (default False; passes upstream
                              --movie and renders frames with PyMOL — needs
                              ProteinEnv and adds GPU/time cost)
        seed:                 random seed (default 42)

    Outputs:
        Streams:
            structures:  ranked pose SDF per (protein, ligand, rank)
            movies:      (only when make_movie=True) one mp4/gif per kept pose
        Tables:
            affinity:    id | structures.id | compounds.id | rank | lddt | affinity
            missing:     id | removed_by | cause

    Usage::

        with Pipeline(project="Examples", job="DynamicBind-demo"):
            Resources(gpu="A100", memory="32GB", time="6:00:00")
            target = PDB("4ufc", convert="pdb")
            lig = Ligand(smiles="CC(=O)Oc1ccccc1C(=O)O", ids="ASA")
            dock = DynamicBind(structures=target, compounds=lig)
    """

    TOOL_NAME = "DynamicBind"
    TOOL_VERSION = "1.1"

    # ------------------------------------------------------------------
    # Install
    # ------------------------------------------------------------------

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        """Clone luwei0917/DynamicBind, create both conda envs (dynamicbind +
        dynamicbind_relax), fetch the workdir.zip weights from Zenodo.

        Verification: both envs exist, torch importable in inference env,
        openmm importable in relax env, workdir/ directory present.
        """
        biopipelines = folders.get("biopipelines", "")
        repo_dir = folders.get("DynamicBind", "")
        parent_dir = os.path.dirname(repo_dir)
        workdir = os.path.join(repo_dir, "workdir")

        inf_check = cls._env_exists_check("dynamicbind", env_manager)
        relax_check = cls._env_exists_check("dynamicbind_relax", env_manager)
        repo_check = f'[ -f "{repo_dir}/run_single_protein_inference.py" ]'
        workdir_check = f'[ -d "{workdir}" ]'

        skip = "" if force_reinstall else f"""# Check if already installed
if {repo_check} && {inf_check} && {relax_check} && {workdir_check} \\
   && {env_manager} run -n dynamicbind python -c "import torch" >/dev/null 2>&1 \\
   && {env_manager} run -n dynamicbind_relax python -c "import openmm" >/dev/null 2>&1; then
    echo "DynamicBind already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        clone_block = f"""mkdir -p {parent_dir}
cd {parent_dir}
if [ ! -d "{repo_dir}" ]; then
    git clone https://github.com/luwei0917/DynamicBind.git "{repo_dir}"
fi"""

        inf_remove = cls._env_remove_block("dynamicbind", env_manager) if force_reinstall else ""
        relax_remove = cls._env_remove_block("dynamicbind_relax", env_manager) if force_reinstall else ""
        inf_env = cls._env_install_block("dynamicbind", env_manager, biopipelines)
        relax_env = cls._env_install_block("dynamicbind_relax", env_manager, biopipelines)

        weights_block = f"""# Fetch workdir.zip weights from Zenodo (~hundreds of MB).
if [ ! -d "{workdir}" ]; then
    cd "{repo_dir}"
    wget -q --show-progress "{DYNAMICBIND_WORKDIR_ZIP_URL}" -O workdir.zip
    unzip -q workdir.zip
    rm -f workdir.zip
fi"""

        return f"""echo "=== Installing DynamicBind ==="
{skip}{clone_block}

{inf_remove}
echo "--- Creating dynamicbind (inference) env ---"
{inf_env}
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create dynamicbind env."
    exit 1
fi

{relax_remove}
echo "--- Creating dynamicbind_relax env ---"
{relax_env}
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create dynamicbind_relax env."
    exit 1
fi

{weights_block}

# Verify installation
if {repo_check} && {workdir_check} \\
   && {env_manager} run -n dynamicbind python -c "import torch" >/dev/null 2>&1 \\
   && {env_manager} run -n dynamicbind_relax python -c "import openmm" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== DynamicBind installation complete ==="
else
    echo "ERROR: DynamicBind verification failed"
    exit 1
fi
"""

    # ------------------------------------------------------------------
    # Lazy path descriptors
    # ------------------------------------------------------------------

    structures_json = Path(lambda self: self.configuration_path(".input_structures.json"))
    compounds_json = Path(lambda self: self.configuration_path(".input_compounds.json"))
    # Per-protein ligand CSVs land under _configuration/ligand_csvs/<prot>.csv,
    # written by the staging pipe script.
    ligand_csvs_dir = Path(lambda self: self.configuration_path("ligand_csvs"))
    results_root = Path(lambda self: self.execution_folder)
    structures_map = Path(lambda self: self.stream_map_path("structures"))
    movies_map = Path(lambda self: self.stream_map_path("movies"))
    affinity_csv = Path(lambda self: self.table_path("affinity"))
    missing_csv = Path(lambda self: self.table_path("missing"))

    movie_jobs_csv = Path(lambda self: self.configuration_path("movie_jobs.csv"))
    build_csv_py = Path(lambda self: self.pipe_script_path("pipe_dynamicbind_build_csv.py"))
    postprocess_py = Path(lambda self: self.pipe_script_path("pipe_dynamicbind_postprocess.py"))
    movie_py = Path(lambda self: self.pipe_script_path("pipe_dynamicbind_movie.py"))
    db_inference_py = Path(lambda self: os.path.join(self.folders["DynamicBind"], "run_single_protein_inference.py"))

    # ------------------------------------------------------------------
    # Constructor
    # ------------------------------------------------------------------

    def __init__(
        self,
        structures: Union[DataStream, StandardizedOutput],
        compounds: Union[DataStream, StandardizedOutput],
        num_samples: int = 40,
        num_saved: Optional[int] = None,
        inference_steps: int = 20,
        num_workers: int = 1,
        rigid_protein: bool = False,
        make_movie: bool = False,
        seed: int = 42,
        **kwargs,
    ):
        # Keep original inputs for predict_output_ids_with_provenance.
        self.structures_input = structures
        self.compounds_input = compounds

        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(
                f"structures must be DataStream or StandardizedOutput, got {type(structures).__name__}"
            )

        if isinstance(compounds, StandardizedOutput):
            self.compounds_stream: DataStream = compounds.streams.compounds
        elif isinstance(compounds, DataStream):
            self.compounds_stream = compounds
        else:
            raise ValueError(
                f"compounds must be DataStream or StandardizedOutput, got {type(compounds).__name__}"
            )

        self.num_samples = num_samples
        # None keeps every generated pose.
        self.num_saved = num_saved if num_saved is not None else num_samples
        self.inference_steps = inference_steps
        self.num_workers = num_workers
        self.rigid_protein = rigid_protein
        self.make_movie = bool(make_movie)
        self.seed = seed
        super().__init__(**kwargs)

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        if not self.structures_stream.has_only_formats("pdb"):
            raise ValueError(
                f"DynamicBind requires PDB-format structures, got '{self.structures_stream.format}'."
            )
        if not self.compounds_stream or len(self.compounds_stream) == 0:
            raise ValueError("compounds parameter is required and must not be empty")

        if self.num_samples < 1:
            raise ValueError("num_samples must be at least 1")
        if self.num_saved < 1:
            raise ValueError("num_saved must be at least 1")
        if self.num_saved > self.num_samples:
            raise ValueError(
                f"num_saved ({self.num_saved}) cannot exceed num_samples ({self.num_samples})"
            )
        if self.inference_steps < 1:
            raise ValueError("inference_steps must be at least 1")
        if self.num_workers < 1:
            raise ValueError("num_workers must be at least 1")
        if self.seed < 0:
            raise ValueError("seed must be non-negative")

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
            f"SAMPLES GENERATED PER PAIR: {self.num_samples}",
            f"POSES SAVED PER PAIR: {self.num_saved}",
            f"INFERENCE STEPS: {self.inference_steps}",
            f"RIGID PROTEIN: {self.rigid_protein}",
            f"MAKE MOVIE: {self.make_movie}",
        ])
        return lines

    # ------------------------------------------------------------------
    # Script generation
    # ------------------------------------------------------------------

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        self.compounds_stream.save_json(self.compounds_json)

        script = "#!/bin/bash\n"
        script += "# DynamicBind execution script\n"
        script += self.generate_completion_check_header()
        script += self.warn_container_unsupported()

        # Phase 1: build per-protein ligand CSVs (biopipelines env).
        script += "# --- Build per-protein ligand CSVs (biopipelines env) ---\n"
        script += self.activate_environment(name="biopipelines")
        script += f"""python "{self.build_csv_py}" \\
    --structures-json "{self.structures_json}" \\
    --compounds-json "{self.compounds_json}" \\
    --output-dir "{self.ligand_csvs_dir}"

if [ $? -ne 0 ]; then
    echo "Error: DynamicBind staging failed"
    exit 1
fi

"""

        # Phase 2: inference loop over proteins (dynamicbind env).
        # run_single_protein_inference.py is itself a multi-env orchestrator
        # (it shells out to --python and --relax_python), so we resolve the
        # absolute env Python paths up front and let it dispatch.
        script += "# --- DynamicBind inference (dynamicbind env) ---\n"
        script += self.activate_environment()  # dynamicbind
        rigid_flag = " --rigid_protein" if self.rigid_protein else ""
        movie_flag = " --movie" if self.make_movie else ""
        relax_python_cmd = ConfigManager().get_env_python_command("dynamicbind_relax")
        script += f"""# Locate the two env Pythons via the same env_manager BioPipelines uses.
DB_PYTHON=$(python -c "import sys; print(sys.executable)")
RELAX_PYTHON={relax_python_cmd}

# rdkit (pulled as a transitive pip dep at v2026.03.2) links against a newer
# libstdc++ than Colab's system one (system caps at GLIBCXX_3.4.30; rdkit
# needs 3.4.31+). The env's lib/ has the right libstdc++.so.6.0.34, but
# `micromamba run` doesn't always inject it into LD_LIBRARY_PATH for child
# processes — and DynamicBind's wrapper spawns the inference and relax
# Pythons as subprocesses. Prepend both env lib dirs explicitly.
DB_ENV_LIB="$(dirname "$DB_PYTHON")/../lib"
RELAX_ENV_LIB="$(dirname "$RELAX_PYTHON")/../lib"
export LD_LIBRARY_PATH="$DB_ENV_LIB:$RELAX_ENV_LIB:$LD_LIBRARY_PATH"

mkdir -p "{self.results_root}"

cd "{self.folders["DynamicBind"]}"
for struct_id in {Resolve.stream_ids(self.structures_json)}; do
    PDB_FILE=$(resolve_stream_item "{self.structures_json}" "$struct_id")
    LIG_CSV="{self.ligand_csvs_dir}/${{struct_id}}.csv"

    if [ ! -f "$LIG_CSV" ]; then
        echo "  [skip] $struct_id: ligand CSV not found"
        continue
    fi

    echo "  [run]  DynamicBind on $struct_id"
    set +e
    python "{self.db_inference_py}" \\
        "$PDB_FILE" "$LIG_CSV" \\
        --header "$struct_id" \\
        --results "{self.results_root}" \\
        --samples_per_complex {self.num_samples} \\
        --savings_per_complex {self.num_saved} \\
        --inference_steps {self.inference_steps} \\
        --num_workers {self.num_workers} \\
        --seed {self.seed}{rigid_flag}{movie_flag} \\
        --python "$DB_PYTHON" \\
        --relax_python "$RELAX_PYTHON"
    rc=$?
    set -e
    if [ $rc -ne 0 ]; then
        echo "  [fail] $struct_id (exit $rc)"
    fi
done

"""

        # Phase 3: post-process (biopipelines env).
        script += "# --- Post-processing (biopipelines env) ---\n"
        script += self.activate_environment(name="biopipelines")
        movie_jobs_arg = f' \\\n    --movie-jobs-csv "{self.movie_jobs_csv}"' if self.make_movie else ""
        script += f"""python "{self.postprocess_py}" \\
    --results-root "{self.results_root}" \\
    --ligand-csvs-dir "{self.ligand_csvs_dir}" \\
    --structures-folder "{self.stream_folder("structures")}" \\
    --structures-map "{self.structures_map}" \\
    --affinity-csv "{self.affinity_csv}" \\
    --missing-csv "{self.missing_csv}" \\
    --num-saved {self.num_saved}{movie_jobs_arg}

if [ $? -ne 0 ]; then
    echo "Error: DynamicBind post-processing failed"
    exit 1
fi

"""

        # Phase 4 (optional): render a movie per kept pose with PyMOL (ProteinEnv).
        if self.make_movie:
            script += "# --- Movie rendering (ProteinEnv / PyMOL) ---\n"
            script += self.activate_environment(name="ProteinEnv")
            script += f"""python "{self.movie_py}" \\
    --movie-jobs-csv "{self.movie_jobs_csv}" \\
    --movies-folder "{self.stream_folder("movies")}" \\
    --movies-map "{self.movies_map}"

if [ $? -ne 0 ]; then
    echo "Warning: DynamicBind movie rendering failed (poses still produced)"
fi

"""
        script += self.generate_completion_check_footer()
        return script

    # ------------------------------------------------------------------
    # Output prediction
    # ------------------------------------------------------------------

    def get_output_files(self) -> Dict[str, Any]:
        pair_ids, _ = predict_output_ids_with_provenance(
            structures=(self.structures_input, "structures"),
            compounds=(self.compounds_input, "compounds"),
        )

        ranked_ids: List[str] = []
        for pair_id in pair_ids:
            for n in range(1, self.num_saved + 1):
                ranked_ids.append(f"{pair_id}_rank{n}")

        structures = DataStream(
            name="structures",
            ids=ranked_ids,
            files=[self.stream_path("structures", "<id>.sdf")],
            map_table=self.structures_map,
            format="sdf",
        )

        result_streams = {"structures": structures}
        if self.make_movie:
            result_streams["movies"] = DataStream(
                name="movies",
                ids=ranked_ids,
                files=[self.stream_path("movies", "<id>.mp4")],
                map_table=self.movies_map,
                format="mp4",
            )

        tables = {
            "affinity": TableInfo(
                name="affinity",
                path=self.affinity_csv,
                columns=["id", "structures.id", "compounds.id", "rank", "lddt", "affinity"],
                description="DynamicBind per-pose lDDT and predicted affinity",
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "kind", "cause"],
                description="Pairs DynamicBind could not dock, with reason",
            ),
        }

        return {
            **result_streams,
            "tables": tables,
            "output_folder": self.output_folder,
        }
