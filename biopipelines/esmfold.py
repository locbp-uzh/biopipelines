# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
ESMFold tool for sequence-to-structure prediction.

ESMFold predicts protein 3D structure directly from amino acid sequence using
a language-model-based approach — no MSA required. Each input sequence yields
one predicted PDB structure with per-residue pLDDT confidence scores and a
global pTM score.

Reference:
    ESM: https://github.com/facebookresearch/esm
    Installation: https://github.com/mabr3112/One-command-install-ESMfold
"""

import os
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
except ImportError:
    import sys
    sys.path.insert(0, os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream


class ESMFold(BaseConfig):
    """
    ESMFold: language-model-based protein structure prediction.

    Predicts one PDB structure per input sequence. No MSA or templates needed.
    Confidence is reported as per-residue pLDDT (stored in B-factor column) and
    a global pTM score.

    Inputs:
        sequences:           amino acid sequences (StandardizedOutput or DataStream)
        chunk_size:          axial attention chunk size to reduce GPU memory
                             (default: None = no chunking, matching upstream ESMFold;
                             set 128, 64, or 32 for long sequences that OOM)
        num_recycles:        recycling iterations (default: 4, the upstream default used in training)
        max_tokens_per_batch: maximum tokens per GPU forward pass; lower to avoid OOM
                             on short-sequence batches (default: 1024, the upstream default)
        cpu_offload:         offload weights to CPU to reduce GPU memory (default: False)
        cpu_only:            run entirely on CPU — very slow (default: False)

    Outputs:
        Streams:
            structures:  predicted PDB files, one per input sequence
                         B-factor column = per-residue pLDDT (0–100)
        Tables:
            structures:  id | file | sequences.id
            confidence: id | file | plddt | ptm

    Usage::

        with Pipeline(project="Examples", job="ESMFold-demo"):
            Resources(memory="32GB", time="2:00:00", gpu=1)
            seqs = Sequence(["MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"],
                            ids=["p1"])
            esm = ESMFold(sequences=seqs)

            # Pass structures to downstream tools (e.g. CPred)
            cp = CPred(structures=esm)
    """

    TOOL_NAME = "ESMFold"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        """Create the esmfold env. Picks esmfold.<variant>.yaml automatically:
        the cluster yaml mirrors mabr3112/One-command-install-ESMfold, the
        colab yaml mirrors sokrypton/ColabFold's notebook recipe (sokrypton
        forks of esm + openfold, no nvcc compile). The Colab variant uses
        phased pip files (esmfold.pip.colab.1.txt, .2.txt) so torch is
        installed before openfold's setup.py imports it."""
        biopipelines = folders.get("biopipelines", "")
        env_check = cls._env_exists_check("esmfold", env_manager)
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_check}; then
    echo "esmfold environment already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block("esmfold", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("esmfold", env_manager, biopipelines)
        return f"""echo "=== Installing ESMFold ==="
{skip}{remove_block}
{env_block}
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create esmfold environment."
    exit 1
fi

# Verify installation
if {env_manager} run -n esmfold python -c "import esm" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== ESMFold installation complete ==="
else
    echo "ERROR: ESMFold verification failed (cannot import esm)"
    exit 1
fi
"""

    # ---------------------------------------------------------------------------
    # Lazy path descriptors
    # ---------------------------------------------------------------------------

    sequences_ds_json  = Path(lambda self: self.configuration_path(".input_sequences.json"))
    predictions_folder = Path(lambda self: self.stream_folder("structures"))
    structures_map     = Path(lambda self: self.stream_map_path("structures"))
    confidence_csv     = Path(lambda self: self.table_path("confidence"))
    inference_script   = Path(lambda self: os.path.join(self.folders["pipe_scripts"], "pipe_esmfold_inference.py"))
    helper_script      = Path(lambda self: os.path.join(self.folders["pipe_scripts"], "pipe_esmfold_postprocessing.py"))

    # ---------------------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------------------

    def __init__(
        self,
        sequences: Union[DataStream, StandardizedOutput],
        chunk_size: Optional[int] = None,
        num_recycles: int = 4,
        max_tokens_per_batch: int = 1024,
        cpu_offload: bool = False,
        cpu_only: bool = False,
        **kwargs
    ):
        """
        Configure ESMFold structure prediction.

        Args:
            sequences:            Input amino acid sequences as DataStream or
                                  StandardizedOutput (uses the sequences stream).
            chunk_size:           Axial attention chunk size. Lower values use less
                                  GPU memory at the cost of speed. Default None disables
                                  chunking (upstream ESMFold default); set 128/64/32 for
                                  long sequences that run out of memory.
            num_recycles:         Number of recycling iterations (default: 4).
            max_tokens_per_batch: Max tokens per forward pass for batching short
                                  sequences together (default: 1024).
            cpu_offload:          Offload model weights to CPU (default: False).
            cpu_only:             Run on CPU only — very slow (default: False).
            **kwargs:             Forwarded to BaseConfig.
        """
        if isinstance(sequences, StandardizedOutput):
            self.sequences_stream: DataStream = sequences.streams.sequences
        elif isinstance(sequences, DataStream):
            self.sequences_stream = sequences
        else:
            raise ValueError(
                f"sequences must be DataStream or StandardizedOutput, got {type(sequences).__name__}"
            )

        self.chunk_size = chunk_size
        self.num_recycles = num_recycles
        self.max_tokens_per_batch = max_tokens_per_batch
        self.cpu_offload = cpu_offload
        self.cpu_only = cpu_only
        super().__init__(**kwargs)

    # ---------------------------------------------------------------------------
    # Validation
    # ---------------------------------------------------------------------------

    def validate_params(self):
        if not self.sequences_stream:
            raise ValueError("sequences parameter is required and must not be empty")
        if self.num_recycles < 1:
            raise ValueError("num_recycles must be at least 1")
        if self.max_tokens_per_batch < 1:
            raise ValueError("max_tokens_per_batch must be at least 1")
        if self.chunk_size is not None and self.chunk_size < 1:
            raise ValueError("chunk_size must be a positive integer or None (no chunking)")

    # ---------------------------------------------------------------------------
    # Configure inputs
    # ---------------------------------------------------------------------------

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    # ---------------------------------------------------------------------------
    # Config display
    # ---------------------------------------------------------------------------

    def get_config_display(self) -> List[str]:
        config_lines = super().get_config_display()
        config_lines.extend([
            f"INPUT SEQUENCES: {len(self.sequences_stream)}",
            f"CHUNK SIZE: {self.chunk_size if self.chunk_size else 'disabled'}",
            f"NUM RECYCLES: {self.num_recycles}",
            f"MAX TOKENS PER BATCH: {self.max_tokens_per_batch}",
        ])
        if self.cpu_offload:
            config_lines.append("CPU OFFLOAD: enabled")
        if self.cpu_only:
            config_lines.append("CPU ONLY: enabled")
        return config_lines

    # ---------------------------------------------------------------------------
    # Script generation
    # ---------------------------------------------------------------------------

    def generate_script(self, script_path: str) -> str:
        # Serialize DataStream for use by post-processing helper at runtime
        self.sequences_stream.save_json(self.sequences_ds_json)

        # The sequences CSV (map_table) is the source of id,sequence pairs
        sequences_csv = self.sequences_stream.map_table

        script_content = "#!/bin/bash\n"
        script_content += "# ESMFold structure prediction script\n"
        script_content += self.generate_completion_check_header()

        # Phase 1: inference under esmfold env
        script_content += self.activate_environment()  # esmfold (primary env)

        n_sequences = len(self.sequences_stream)
        chunk_str = str(self.chunk_size) if self.chunk_size else "None"
        script_content += f'echo "Running ESMFold on {n_sequences} sequence(s)"\n'
        script_content += f'echo "Chunk size: {chunk_str} | Recycles: {self.num_recycles}"\n\n'

        # Build inference command
        optional_flags = ""
        if self.chunk_size:
            optional_flags += f" --chunk-size {self.chunk_size}"
        if self.cpu_offload:
            optional_flags += " --cpu-offload"
        if self.cpu_only:
            optional_flags += " --cpu-only"

        script_content += f"""python "{self.inference_script}" \\
    --sequences-csv "{sequences_csv}" \\
    --output-dir "{self.predictions_folder}" \\
    --num-recycles {self.num_recycles} \\
    --max-tokens-per-batch {self.max_tokens_per_batch}{optional_flags}

if [ $? -ne 0 ]; then
    echo "Error: ESMFold inference failed"
    exit 1
fi

"""

        # Phase 2: post-processing under biopipelines env
        script_content += "# --- Post-processing (biopipelines env) ---\n"
        script_content += self.activate_environment(name="biopipelines")

        script_content += f"""python "{self.helper_script}" \\
    --sequences-json "{self.sequences_ds_json}" \\
    --predictions-dir "{self.predictions_folder}" \\
    --map-csv "{self.structures_map}" \\
    --merged-csv "{self.confidence_csv}"

if [ $? -eq 0 ]; then
    echo "ESMFold post-processing completed successfully"
else
    echo "Error: ESMFold post-processing failed"
    exit 1
fi

"""
        script_content += self.generate_completion_check_footer()
        return script_content

    # ---------------------------------------------------------------------------
    # Output files
    # ---------------------------------------------------------------------------

    def get_output_files(self) -> Dict[str, Any]:
        structures = DataStream(
            name="structures",
            ids=self.sequences_stream.ids,
            files=[self.stream_path("structures", "<id>.pdb")],
            map_table=self.structures_map,
            format="pdb"
        )

        tables = {
            "structures": TableInfo(
                name="structures",
                path=self.structures_map,
                columns=["id", "file", "sequences.id"],
                description="ESMFold predicted structures map table",
            ),
            "confidence": TableInfo(
                name="confidence",
                path=self.confidence_csv,
                columns=["id", "file", "plddt", "ptm"],
                description="ESMFold confidence metrics (mean pLDDT and pTM) for all predicted structures",
            ),
        }

        return {
            "structures": structures,
            "sequences": DataStream.empty("sequences", "fasta"),
            "compounds": DataStream.empty("compounds", "sdf"),
            "tables": tables,
            "output_folder": self.output_folder,
            "rendering_parameters": {
                "structures": {
                    "color_by": "plddt",
                    "plddt_upper": 100,
                }
            }
        }

    # ---------------------------------------------------------------------------
    # Serialization
    # ---------------------------------------------------------------------------

    def to_dict(self) -> Dict[str, Any]:
        base_dict = super().to_dict()
        base_dict.update({
            "esmfold_params": {
                "chunk_size": self.chunk_size,
                "num_recycles": self.num_recycles,
                "max_tokens_per_batch": self.max_tokens_per_batch,
                "cpu_offload": self.cpu_offload,
                "cpu_only": self.cpu_only,
            }
        })
        return base_dict
