# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
NeuralPLexer tool for protein-ligand complex prediction.

NeuralPLexer is a physics-inspired flow-based generative model that predicts
state-specific protein-ligand complexes from a protein sequence/structure plus
a ligand SMILES/SDF, without requiring a binding-box hint.

Reference:
    Qiao et al. (2024) NeuralPLexer. https://github.com/zrqiao/NeuralPLexer
"""

import os
from typing import Dict, List, Any, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .combinatorics import predict_output_ids_with_provenance
    from .biopipelines_io import Resolve
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from combinatorics import predict_output_ids_with_provenance
    from biopipelines_io import Resolve


# Zenodo deposit for the pretrained checkpoints (CC BY-NC-SA 4.0).
NEURALPLEXER_WEIGHTS_URL = "https://zenodo.org/records/10373581/files/neuralplexermodels_downstream_datasets_predictions.zip"
# Top-level folder the zip expands into (under repo_dir/data/).
NEURALPLEXER_BUNDLE = "neuralplexermodels_downstream_datasets_predictions"


class NeuralPLexer(BaseConfig):
    """
    NeuralPLexer: protein-ligand complex prediction with state-specific sampling.

    For each (protein, ligand) pair, samples `n_samples` candidate complexes
    and (optionally) re-ranks them by confidence. Output IDs are the
    multi-axis cartesian product <protein>+<ligand>_rank<N>.

    Inputs:
        structures:   protein structures (PDB stream)
        compounds:    ligands as SMILES or SDF (compounds stream)
        n_samples:    candidate complexes per pair (default 16)
        num_steps:    diffusion sampling steps (default 40)
        chunk_size:   parallel batch size within a pair (default 4)
        sampler:      one of "langevin_simulated_annealing", "DDIM", "VPSDE"
                      (default: "langevin_simulated_annealing")
        cuda:         use GPU acceleration (default True; set False only for
                      CPU debug runs — sampling on CPU is impractically slow)
        use_template: pass a template PDB to bias the prediction (deferred)

    Outputs:
        Streams:
            structures:  per-rank predicted complex PDB (protein + ligand merged)
        Tables:
            confidence:  id | structures.id | compounds.id | rank | confidence
            missing:     id | removed_by | cause
    """

    TOOL_NAME = "NeuralPLexer"
    TOOL_VERSION = "1.0"

    # ------------------------------------------------------------------
    # Install
    # ------------------------------------------------------------------

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        """Clone zrqiao/NeuralPLexer, create env, install package, fetch
        Zenodo weights bundle. Verification: env exists, neuralplexer
        importable, checkpoint file present.
        """
        biopipelines = folders.get("biopipelines", "")
        repo_dir = folders.get("NeuralPLexer", "")
        parent_dir = os.path.dirname(repo_dir)
        # The Zenodo zip expands to a single top-level folder
        # (neuralplexermodels_downstream_datasets_predictions/) holding
        # models/, downstream_test_datasets/, predictions/ — it is NOT
        # flattened into data/. Keep the ckpt path in sync with NPLEXER_BUNDLE.
        models_dir = os.path.join(repo_dir, "data", NEURALPLEXER_BUNDLE, "models")
        ckpt_file = os.path.join(models_dir, "complex_structure_prediction.ckpt")

        env_check = cls._env_exists_check("neuralplexer", env_manager)
        repo_check = f'[ -d "{repo_dir}/neuralplexer" ]'
        ckpt_check = f'[ -f "{ckpt_file}" ]'

        skip = "" if force_reinstall else f"""# Check if already installed
if {repo_check} && {env_check} && {ckpt_check} \\
   && {env_manager} run -n neuralplexer python -c "import neuralplexer" >/dev/null 2>&1; then
    echo "NeuralPLexer already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        clone_block = f"""mkdir -p {parent_dir}
cd {parent_dir}
if [ ! -d "{repo_dir}" ]; then
    git clone https://github.com/zrqiao/NeuralPLexer.git "{repo_dir}"
fi"""

        remove_block = cls._env_remove_block("neuralplexer", env_manager) if force_reinstall else ""

        env_block = cls._env_install_block("neuralplexer", env_manager, biopipelines)

        # openfold builds a CUDA extension at install time. The build needs a
        # specific environment, all set *inside* the `micromamba run` subshell
        # (where $CONDA_PREFIX is the env and the activation hooks ran) — setting
        # these in the outer shell does not work, because micromamba run resets
        # the compiler vars on activation:
        #   sed c++14->c++17  — torch 2.1 headers require C++17; the pinned
        #     commit hardcodes nvcc -std=c++14 in setup.py.
        #   CC/CXX + NVCC_PREPEND_FLAGS=-ccbin — build with the env's gcc 10.
        #     nvcc 12.1 driving Colab's gcc 12 fails on torch 2.1's pybind11
        #     templates (`expected template-name`); gcc 10 compiles them cleanly.
        #   LDFLAGS -L$CONDA_PREFIX/lib — openfold links libcudart.so from the
        #     env lib/, which conda's compiler_compat ld does not search; rpath
        #     it too so the built .so resolves libcudart at load time.
        # Phase 3 has already downgraded setuptools<81 by the time this runs.
        openfold_dir = os.path.join(parent_dir, "openfold_build")
        openfold_block = f"""echo "--- build + install openfold (CUDA extension) ---"
{env_manager} run -n neuralplexer bash -c '
set -e
rm -rf "{openfold_dir}"
git clone --filter=blob:none --quiet https://github.com/aqlaboratory/openfold.git "{openfold_dir}"
cd "{openfold_dir}"
git checkout -q 103d0370ad9ce07579c20fa9c889a632f9b16618
sed -i "s/c++14/c++17/g" setup.py
# openfold setup.py hardcodes a compute-capability set incl. (3,7)=K80, emitting
# -gencode arch=compute_37 which CUDA 12.1 nvcc rejects ("Unsupported gpu
# architecture compute_37"). It narrows to the local GPU only if get_nvidia_cc()
# detects one — but installs run on a GPU-less node, so the full set is used.
# Drop the (3,7) entry (the only arch nvcc 12.1 refuses; 5.2/6.1/7.0/8.0 remain).
sed -i "/(3, 7),.*K80/d" setup.py
export CC="$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-gcc"
export CXX="$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-g++"
export NVCC_PREPEND_FLAGS="-ccbin $CXX"
export LDFLAGS="-L$CONDA_PREFIX/lib -Wl,-rpath,$CONDA_PREFIX/lib"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
pip install . --no-build-isolation --no-deps
'"""

        # NeuralPLexer is the package itself; install editable from the clone.
        pkg_install = f"""echo "--- pip install -e NeuralPLexer ---"
{env_manager} run -n neuralplexer pip install -e "{repo_dir}" --no-deps"""

        # The Zenodo bundle is ~8.7 GB; over a cluster compute-node proxy the
        # download is prone to truncation, and a partial zip then fails extract
        # with "End-of-central-directory signature not found". Harden it:
        #   * wget -c (resume a partial file) with built-in retries;
        #   * verify the archive with `unzip -t` BEFORE extracting; if it's
        #     incomplete, resume/re-download and re-check, up to N attempts;
        #   * only extract (and only delete the zip) once the archive is intact;
        #   * fail loudly with a clear message instead of leaving a half-unzipped
        #     data/ that the verification step would misreport as "ckpt missing".
        weights_block = f"""# Fetch + verify the Zenodo weights bundle if absent. Expands to data/.
if [ ! -f "{ckpt_file}" ]; then
    cd "{repo_dir}"
    mkdir -p data
    ZIP=nplexer_weights.zip
    ok=0
    for attempt in 1 2 3; do
        echo "Downloading NeuralPLexer weights (attempt $attempt)..."
        # -c resumes a partial $ZIP; --tries/--timeout ride out flaky proxies.
        wget -c --tries=5 --timeout=60 --show-progress "{NEURALPLEXER_WEIGHTS_URL}" -O "$ZIP" || true
        # Integrity check the whole archive before trusting it.
        if unzip -t "$ZIP" >/dev/null 2>&1; then
            ok=1
            break
        fi
        echo "Weights archive incomplete/corrupt after attempt $attempt; retrying download (resume)."
    done
    if [ "$ok" -ne 1 ]; then
        echo "ERROR: NeuralPLexer weights download could not be completed/verified after 3 attempts."
        echo "       (truncated zip from {NEURALPLEXER_WEIGHTS_URL}). Re-run install to resume."
        rm -f "$ZIP"
        exit 1
    fi
    # Extract ONLY the models/ subtree. The Zenodo bundle also ships a multi-GB
    # `predictions/` benchmark tree and `downstream_test_datasets/` that the tool
    # never uses at runtime — extracting everything wasted ~9 GB and filled the
    # home quota mid-unzip ("write error (disk full?)"), truncating the ckpt and
    # failing verification. The selective extract keeps only what's needed.
    # -o: overwrite without prompting (non-interactive installs must not hang).
    unzip -oq "$ZIP" "{NEURALPLEXER_BUNDLE}/models/*" -d data/
    rm -f "$ZIP"
fi"""

        return f"""echo "=== Installing NeuralPLexer ==="
{skip}{clone_block}

{remove_block}
{env_block}
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create neuralplexer environment."
    exit 1
fi

{openfold_block}
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to build/install openfold."
    exit 1
fi

{pkg_install}
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to install NeuralPLexer package."
    exit 1
fi

{weights_block}

# Verify installation. Importing neuralplexer pulls openfold.model.primitives,
# which loads the compiled attn_core_inplace_cuda — so this import also proves
# the CUDA extension built and loads.
if {repo_check} && {ckpt_check} \\
   && {env_manager} run -n neuralplexer python -c "import torch; import neuralplexer" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== NeuralPLexer installation complete ==="
else
    echo "ERROR: NeuralPLexer verification failed (repo, ckpt, or import missing)"
    exit 1
fi
"""

    # ------------------------------------------------------------------
    # Lazy path descriptors
    # ------------------------------------------------------------------

    structures_json = Path(lambda self: self.configuration_path(".input_structures.json"))
    compounds_json = Path(lambda self: self.configuration_path(".input_compounds.json"))
    staging_folder = Path(lambda self: self.configuration_path("staging"))
    raw_out_folder = Path(lambda self: self.execution_folder)
    structures_map = Path(lambda self: self.stream_map_path("structures"))
    confidence_csv = Path(lambda self: self.table_path("confidence"))
    missing_csv = Path(lambda self: self.table_path("missing"))

    stage_py = Path(lambda self: self.pipe_script_path("pipe_neuralplexer_stage.py"))
    postprocess_py = Path(lambda self: self.pipe_script_path("pipe_neuralplexer_postprocess.py"))
    nplexer_ckpt = Path(
        lambda self: os.path.join(self.folders["NeuralPLexer"], "data", NEURALPLEXER_BUNDLE, "models", "complex_structure_prediction.ckpt")
    )

    # ------------------------------------------------------------------
    # Constructor
    # ------------------------------------------------------------------

    def __init__(
        self,
        structures: Union[DataStream, StandardizedOutput],
        compounds: Union[DataStream, StandardizedOutput],
        n_samples: int = 16,
        num_steps: int = 40,
        chunk_size: int = 4,
        sampler: str = "langevin_simulated_annealing",
        cuda: bool = True,
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

        self.n_samples = n_samples
        self.num_steps = num_steps
        self.chunk_size = chunk_size
        self.sampler = sampler
        self.cuda = cuda
        super().__init__(**kwargs)

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        if not self.structures_stream.has_only_formats("pdb"):
            raise ValueError(
                f"NeuralPLexer requires PDB-format structures, got '{self.structures_stream.format}'."
            )
        if not self.compounds_stream or len(self.compounds_stream) == 0:
            raise ValueError("compounds parameter is required and must not be empty")

        if self.n_samples < 1:
            raise ValueError("n_samples must be at least 1")
        if self.num_steps < 1:
            raise ValueError("num_steps must be at least 1")
        if self.chunk_size < 1:
            raise ValueError("chunk_size must be at least 1")
        # Upstream samples in chunks: `for _ in range(n_samples // chunk_size)`.
        # If n_samples < chunk_size that loop runs zero times and produces no
        # structures (an opaque IndexError downstream). Require n_samples to be
        # a positive multiple of chunk_size so every requested sample is drawn.
        if self.n_samples % self.chunk_size != 0:
            raise ValueError(
                f"n_samples ({self.n_samples}) must be a multiple of chunk_size "
                f"({self.chunk_size}); upstream draws n_samples // chunk_size chunks, "
                f"so a smaller n_samples silently produces zero structures."
            )

        valid_samplers = {"langevin_simulated_annealing", "DDIM", "VPSDE"}
        if self.sampler not in valid_samplers:
            raise ValueError(f"sampler must be one of {sorted(valid_samplers)}; got {self.sampler!r}")

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
            f"N SAMPLES PER PAIR: {self.n_samples}",
            f"NUM STEPS: {self.num_steps}",
            f"CHUNK SIZE: {self.chunk_size}",
            f"SAMPLER: {self.sampler}",
            f"CUDA: {self.cuda}",
        ])
        return lines

    # ------------------------------------------------------------------
    # Script generation
    # ------------------------------------------------------------------

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        self.compounds_stream.save_json(self.compounds_json)

        script = "#!/bin/bash\n"
        script += "# NeuralPLexer execution script\n"
        script += self.generate_completion_check_header()

        # Phase 1: stage pairs (biopipelines env).
        script += "# --- Stage pairs (biopipelines env) ---\n"
        script += self.activate_environment(name="biopipelines")
        script += f"""python "{self.stage_py}" \\
    --structures-json "{self.structures_json}" \\
    --compounds-json "{self.compounds_json}" \\
    --staging-folder "{self.staging_folder}"

if [ $? -ne 0 ]; then
    echo "Error: NeuralPLexer staging failed"
    exit 1
fi

"""

        # Phase 2: per-pair inference loop (neuralplexer env). Upstream's
        # neuralplexer-inference is single-pair, so we loop ourselves.
        script += "# --- NeuralPLexer inference (neuralplexer env) ---\n"
        script += self.activate_environment()  # neuralplexer
        cuda_flag = " \\\n        --cuda" if self.cuda else ""
        # Colab exports MPLBACKEND=module://matplotlib_inline.backend_inline,
        # which leaks into this subprocess; matplotlib (pulled in transitively by
        # torchmetrics) rejects it as an invalid backend. Force a headless one.
        script += f"""export MPLBACKEND=Agg
mkdir -p "{self.raw_out_folder}"

for pair_dir in "{self.staging_folder}"/*/; do
    pair_id=$(basename "$pair_dir")
    receptor="${{pair_dir}}receptor.pdb"
    ligand="${{pair_dir}}ligand.sdf"
    out_path="{self.raw_out_folder}/${{pair_id}}"
    mkdir -p "$out_path"

    if [ ! -f "$receptor" ] || [ ! -f "$ligand" ]; then
        echo "  [skip] $pair_id: missing receptor or ligand"
        continue
    fi

    echo "  [run]  NeuralPLexer on $pair_id"
    set +e
    neuralplexer-inference --task=batched_structure_sampling \\
        --input-receptor "$receptor" \\
        --input-ligand "$ligand" \\
        --out-path "$out_path" \\
        --model-checkpoint "{self.nplexer_ckpt}" \\
        --n-samples {self.n_samples} \\
        --chunk-size {self.chunk_size} \\
        --num-steps={self.num_steps}{cuda_flag} \\
        --sampler={self.sampler} \\
        --separate-pdb \\
        --rank-outputs-by-confidence
    rc=$?
    set -e
    if [ $rc -ne 0 ]; then
        echo "  [fail] $pair_id (exit $rc)"
    fi
done

"""

        # Phase 3: post-process (biopipelines env).
        script += "# --- Post-processing (biopipelines env) ---\n"
        script += self.activate_environment(name="biopipelines")
        script += f"""python "{self.postprocess_py}" \\
    --raw-out "{self.raw_out_folder}" \\
    --structures-folder "{self.stream_folder("structures")}" \\
    --structures-map "{self.structures_map}" \\
    --confidence-csv "{self.confidence_csv}" \\
    --missing-csv "{self.missing_csv}" \\
    --n-samples {self.n_samples}

if [ $? -ne 0 ]; then
    echo "Error: NeuralPLexer post-processing failed"
    exit 1
fi

"""
        script += self.generate_completion_check_footer()
        return script

    # ------------------------------------------------------------------
    # Output prediction
    # ------------------------------------------------------------------

    def get_output_files(self) -> Dict[str, Any]:
        pair_ids, pair_provenance = predict_output_ids_with_provenance(
            structures=(self.structures_input, "structures"),
            compounds=(self.compounds_input, "compounds"),
        )

        ranked_ids: List[str] = []
        for pair_id in pair_ids:
            for n in range(1, self.n_samples + 1):
                ranked_ids.append(f"{pair_id}_rank{n}")

        structures = DataStream(
            name="structures",
            ids=ranked_ids,
            files=[self.stream_path("structures", "<id>.pdb")],
            map_table=self.structures_map,
            format="pdb",
        )

        tables = {
            "confidence": TableInfo(
                name="confidence",
                path=self.confidence_csv,
                columns=["id", "structures.id", "compounds.id", "rank", "confidence"],
                description="NeuralPLexer per-rank confidence (re-ranked outputs)",
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "kind", "cause"],
                description="Pairs NeuralPLexer could not predict, with reason",
            ),
        }

        return {
            "structures": structures,
            "tables": tables,
            "output_folder": self.output_folder,
        }
