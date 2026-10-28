# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""ESM2 tool: pseudo-log-likelihood and embeddings of protein sequences from the ESM-2 language model.

ESM-2 is a masked protein language model trained on UniRef50. This wrapper runs the Hugging Face `transformers` port of the released checkpoints (`facebook/esm2_t*_UR50D`, MIT license) on a `sequences` stream and writes two things per sequence:

  * a pseudo-log-likelihood (PLL): each residue is masked on its own, the model's log-probability of the true residue is read at that position, and the L values are summed (`pll`) or averaged (`pll_per_residue`). Higher (closer to 0) means the sequence is more probable under the model.
  * a mean-pooled embedding: the final-layer representation (after the final layer norm) averaged over the residue positions, one `.npy` vector per id.

Variant-effect scores (wild-type or masked marginals, Meier et al. 2021) are not offered: they are defined for a substitution relative to a reference sequence, which a single input sequence does not carry. VespaG scores the single-substitution landscape of a sequence.

Reference:
    Paper: Lin et al. (2023) Evolutionary-scale prediction of atomic-level protein structure with a language model. Science 379, 1123-1130. https://www.science.org/doi/10.1126/science.ade2574
    Repo:  https://github.com/facebookresearch/esm
    Weights: https://huggingface.co/facebook/esm2_t33_650M_UR50D (and the sibling checkpoints)
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


# Released ESM-2 checkpoints under facebook/ on the Hub, smallest first; the 15B one is left out (about 60 GB of fp32 weights).
ESM2_MODELS = (
    "esm2_t6_8M_UR50D",
    "esm2_t12_35M_UR50D",
    "esm2_t30_150M_UR50D",
    "esm2_t33_650M_UR50D",
    "esm2_t36_3B_UR50D",
)
DEFAULT_MODEL = "esm2_t33_650M_UR50D"
HUB_ORGANIZATION = "facebook"
DEVICES = ("auto", "cpu", "cuda")

# The 25 residue letters of the ESM-2 vocabulary: the 20 standard amino acids, then X, B, U, Z, O.
RESIDUE_ALPHABET = "LAGVSERTIDPKQNFYMHWCXBUZO"

SCORE_COLUMNS = ["id", "length", "pll", "pll_per_residue"]


class ESM2(BaseConfig):
    """
    ESM2: pseudo-log-likelihood and mean-pooled embeddings of protein sequences with ESM-2.

    Runs `EsmForMaskedLM` from Hugging Face `transformers` in fp32 on weights downloaded and pinned by `ESM2.install()`. Output ids equal input ids.

    Inputs:
        sequences:          protein sequences (one chain per id). Letters outside the ESM-2 alphabet (the 20 standard amino acids and X, B, U, Z, O, upper case) fail that id; sequences are never truncated.
        model:              ESM-2 checkpoint, one of esm2_t6_8M_UR50D, esm2_t12_35M_UR50D, esm2_t30_150M_UR50D, esm2_t33_650M_UR50D, esm2_t36_3B_UR50D (default: esm2_t33_650M_UR50D). The default is the size of the model card's own example: 2.6 GB of weights, which fits any current GPU, and the largest checkpoint the Hub ships as safetensors (the 3B one is two sharded pickle files, about 11 GB). The 15B checkpoint is not offered: 60 GB of fp32 weights do not fit one GPU without sharding or half precision, which this tool does not do.
        compute_pll:        write the `scores` table (default: True). It costs one forward pass per residue, so turn it off when only embeddings are wanted.
        compute_embeddings: write the `embeddings` stream (default: True). One forward pass per sequence.
        batch_size:         masked copies of one sequence scored per forward pass (default: 16). Lower it on out-of-memory; it changes memory and speed, not the result.
        device:             "auto" (cuda when available, else cpu), "cpu" or "cuda" (default: "auto"). "cuda" without a GPU fails the step rather than falling back.

    Pseudo-log-likelihood (Salazar et al. 2020; the same sum as the `pseudo-ppl` strategy of the variant-prediction example in facebookresearch/esm, which returns a sum of log-probabilities despite its name): for a sequence x_1..x_L tokenized as <cls> x_1..x_L <eos>, residue i is replaced by <mask> with every other token left as is, one forward pass per i, and log p(x_i | x without i) is the log-softmax over the model's full 33-token output vocabulary at that position, taken at the true token. `pll` is the sum over i, `pll_per_residue` the mean. The sum grows with length, so compare sequences of different lengths on `pll_per_residue`. Both are at most 0.

    Embeddings: one unmasked forward pass; the last hidden state of the encoder (after its final layer norm) at the L residue positions, `<cls>` and `<eos>` excluded, averaged over positions. This is the quantity `esm-extract --include mean` of fair-esm writes for the last layer, except that fair-esm truncates at 1022 residues by default and this tool does not. Shape (hidden_size,): 320, 480, 640, 1280, 2560 for the five checkpoints from smallest to largest, dtype float32.

    Outputs:
        Streams:
            embeddings: <id>.npy per sequence, the mean-pooled vector above (format npy). Empty when compute_embeddings is False.
        Tables:
            scores:  id | length | pll | pll_per_residue   (one row per input id; NaN for an id that failed. Absent when compute_pll is False.)
            missing: id | removed_by | kind | cause
    """

    TOOL_NAME = "ESM2"
    TOOL_VERSION = "1.0"
    ENV_NAME = "esm2"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, models=(DEFAULT_MODEL,), **kwargs):
        env = cls._install_env(env_manager)
        biopipelines = folders.get("biopipelines", "")
        weights_dir = folders.get("ESM2", "")
        if not weights_dir:
            raise ValueError(
                "ESM2: folders.repositories.ESM2 is not set in the config, so there is "
                "nowhere to put the weights. Add `ESM2: \"<data>/esm2\"`.")
        if isinstance(models, str):
            models = (models,)
        models = tuple(models)
        bad = [m for m in models if m not in ESM2_MODELS]
        if not models or bad:
            raise ValueError(f"ESM2.install: models must be a non-empty subset of {ESM2_MODELS}, got {list(models) or models}")
        env_check = cls._env_exists_check(env, env_manager)
        env_run = cls._env_run(env, env_manager)
        pipe_scripts = folders.get("pipe_scripts", f"{biopipelines}/pipe_scripts")
        weights_py = f"{pipe_scripts}/pipe_esm2_weights.py"
        imports_ok = f'{env_run}python -c "import torch, transformers, huggingface_hub" >/dev/null 2>&1'
        weights_ok = " && ".join(f'[ -f "{weights_dir}/{m}/REVISION" ]' for m in models)
        skip = "" if force_reinstall else f"""# Check if already fully installed (env + pinned weights)
if {env_check} && {imports_ok} && {weights_ok}; then
    echo "ESM2 already installed, skipping. Use force_reinstall=True to reinstall the env."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block(env, env_manager) if force_reinstall else ""
        env_block = cls._env_install_block(env, env_manager, biopipelines)
        return f"""echo "=== Installing ESM2 ({env}) ==="
{skip}{remove_block}
# Create {env} env (skip if it already exists)
if ! {env_check}; then
    {env_block}
else
    echo "{env} environment already exists, skipping creation."
fi

# Download each checkpoint at the Hub commit current now and record it in REVISION (skipped when present)
mkdir -p "{weights_dir}"
{env_run}python "{weights_py}" --dest "{weights_dir}" --models {" ".join(models)}

# Verify installation
if {imports_ok} && {weights_ok}; then
    touch "$INSTALL_SUCCESS"
    echo "=== ESM2 installation complete ==="
else
    echo "ERROR: ESM2 verification failed (torch/transformers missing from the env or a weights folder incomplete)"
    exit 1
fi
"""

    sequences_json = Path(lambda self: self.configuration_path("sequences.json"))
    config_yaml = Path(lambda self: self.configuration_path("config.yaml"))
    embeddings_map = Path(lambda self: self.stream_map_path("embeddings"))
    scores_csv = Path(lambda self: self.table_path("scores"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    local_missing_csv = Path(lambda self: self.execution_path("local_missing.csv"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_esm2.py"))

    def __init__(self,
                 sequences: Union[DataStream, StandardizedOutput],
                 model: str = DEFAULT_MODEL,
                 compute_pll: bool = True,
                 compute_embeddings: bool = True,
                 batch_size: int = 16,
                 device: str = "auto",
                 **kwargs):
        self.sequences = sequences
        if isinstance(sequences, StandardizedOutput):
            self.sequences_stream: DataStream = sequences.streams.sequences
        elif isinstance(sequences, DataStream):
            self.sequences_stream = sequences
        else:
            raise ValueError(f"sequences must be DataStream or StandardizedOutput, got {type(sequences)}")
        self.model = model
        self.compute_pll = compute_pll
        self.compute_embeddings = compute_embeddings
        self.batch_size = batch_size
        self.device = device
        super().__init__(**kwargs)

    def validate_params(self):
        if not self.sequences_stream or len(self.sequences_stream) == 0:
            raise ValueError("sequences parameter is required and must not be empty")
        if self.model not in ESM2_MODELS:
            raise ValueError(f"model must be one of {ESM2_MODELS}, got {self.model!r}")
        for name in ("compute_pll", "compute_embeddings"):
            if not isinstance(getattr(self, name), bool):
                raise ValueError(f"{name} must be True or False, got {getattr(self, name)!r}")
        if not (self.compute_pll or self.compute_embeddings):
            raise ValueError("compute_pll and compute_embeddings are both False, so there is nothing to compute")
        if isinstance(self.batch_size, bool) or not isinstance(self.batch_size, int) or self.batch_size < 1:
            raise ValueError(f"batch_size must be a positive integer, got {self.batch_size!r}")
        if self.device not in DEVICES:
            raise ValueError(f"device must be one of {DEVICES}, got {self.device!r}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"MODEL: {HUB_ORGANIZATION}/{self.model}")
        outputs = [name for name, on in (("pll", self.compute_pll), ("embeddings", self.compute_embeddings)) if on]
        lines.append(f"OUTPUTS: {', '.join(outputs)}")
        lines.append(f"DEVICE: {self.device}")
        return lines

    def _weights_folder(self) -> str:
        root = self.folders.get("ESM2", "")
        if not root:
            raise ValueError(
                "ESM2: no weights folder configured. Add `ESM2: \"<data>/esm2\"` under folders.repositories "
                "and run ESM2.install().")
        return os.path.join(root, self.model)

    def generate_script(self, script_path: str) -> str:
        weights_folder = self._weights_folder()
        self.sequences_stream.save_json(self.sequences_json)
        self._write_config_yaml()
        script = "#!/bin/bash\n"
        script += "# ESM2 script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += self.warn_container_unsupported()
        script += f"""echo "Running ESM-2 ({self.model}) on {len(self.sequences_stream)} sequence(s)"
python "{self.helper_py}" \\
    --config-yaml "{self.config_yaml}" \\
    --sequences-json "{self.sequences_json}" \\
    --weights-dir "{weights_folder}" \\
    --embeddings-dir "{self.stream_folder('embeddings')}" \\
    --embeddings-map-csv "{self.embeddings_map}" \\
    --scores-csv "{self.scores_csv}" \\
    --local-missing-csv "{self.local_missing_csv}"
"""
        script += self.generate_missing_propagation(
            self.sequences, local_missing=self.local_missing_csv, missing_csv=self.missing_csv
        )
        script += self.generate_completion_check_footer()
        return script

    def _write_config_yaml(self):
        import yaml
        with open(self.config_yaml, "w") as f:
            yaml.safe_dump({
                "model": self.model,
                "compute_pll": self.compute_pll,
                "compute_embeddings": self.compute_embeddings,
                "batch_size": self.batch_size,
                "device": self.device,
            }, f, sort_keys=False)

    def get_output_files(self) -> Dict[str, Any]:
        if self.compute_embeddings:
            embeddings = DataStream(
                name="embeddings",
                ids=self.sequences_stream.ids,
                files=[self.stream_path("embeddings", "<id>.npy")],
                map_table=self.embeddings_map,
                format="npy",
            )
        else:
            embeddings = DataStream.empty("embeddings", "npy")
        tables = {"missing": self.missing_table_info(self.missing_csv)}
        if self.compute_pll:
            tables["scores"] = TableInfo(
                name="scores",
                path=self.scores_csv,
                columns=SCORE_COLUMNS,
                description="ESM-2 pseudo-log-likelihood per sequence: sum and per-residue mean (higher, closer to 0 = more probable)",
            )
        return {
            "embeddings": embeddings,
            "tables": tables,
            "output_folder": self.output_folder,
        }
