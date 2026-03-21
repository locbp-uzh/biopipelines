# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
ESMFold configuration for ultra-fast single-sequence protein structure prediction.

Uses Meta's ESMFold model (ESM-2 language model + structure module) for
single-sequence structure prediction without MSA generation.
"""

import os
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream, create_map_table
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream, create_map_table


class ESMFold(BaseConfig):
    """
    Configuration for ESMFold single-sequence structure prediction.

    Predicts protein structures from amino acid sequences using ESMFold,
    which leverages ESM-2 language model embeddings. Much faster than
    AlphaFold2 as it does not require MSA generation.

    Example:
        proteins = Sequence(["MKTVRQ...", "AETGFT..."], ids=["p1", "p2"])
        esm = ESMFold(proteins=proteins)

        # With more recycles for higher accuracy
        esm = ESMFold(proteins=mpnn_output, num_recycles=8)

        # With chunk size for long sequences
        esm = ESMFold(proteins=proteins, chunk_size=128)
    """

    TOOL_NAME = "ESMFold"

    # Lazy path descriptors
    queries_csv = Path(lambda self: os.path.join(self.output_folder, f"{self.pipeline_name}_queries.csv"))
    queries_fasta = Path(lambda self: os.path.join(self.output_folder, f"{self.pipeline_name}_queries.fasta"))
    confidence_csv = Path(lambda self: os.path.join(self.output_folder, "confidence.csv"))
    folding_folder = Path(lambda self: os.path.join(self.output_folder, "Folding"))
    structures_map = Path(lambda self: os.path.join(self.output_folder, "structures_map.csv"))
    missing_csv = Path(lambda self: os.path.join(self.output_folder, "missing.csv"))

    # Helper script paths
    fa_to_csv_fasta_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_fa_to_csv_fasta.py"))
    esm_fold_confidence_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_esm_fold_confidence.py"))
    propagate_missing_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_propagate_missing.py"))
    update_map_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_update_structures_map.py"))

    # HelpScript for Python API inference (Colab/pip)
    esm_fold_predict_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_esm_fold_predict.py"))

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        if env_manager in ("pip", "micromamba"):
            # Colab: use official ESM library (ColabFold approach)
            skip = "" if force_reinstall else """# Check if already installed
if [ -f "esmfold.model" ] && python -c "import esm" 2>/dev/null; then
    echo "ESMFold already installed, skipping. Use force_reinstall=True to reinstall."
    exit 0
fi
"""
            return f"""echo "=== Installing ESMFold (official ESM library) ==="
{skip}# Install dependencies
pip install omegaconf pytorch_lightning biopython ml_collections einops py3Dmol modelcif
pip install git+https://github.com/NVIDIA/dllogger.git
pip install git+https://github.com/sokrypton/openfold.git
pip install git+https://github.com/sokrypton/esm.git

# Download ESMFold model weights
if [ ! -f "esmfold.model" ]; then
    echo "Downloading ESMFold model weights..."
    apt-get install aria2 -qq 2>/dev/null || true
    if command -v aria2c &>/dev/null; then
        aria2c -x 16 https://colabfold.steineggerlab.workers.dev/esm/esmfold.model
    else
        wget https://colabfold.steineggerlab.workers.dev/esm/esmfold.model
    fi
fi

echo "=== ESMFold installation complete ==="
"""
        # Cluster: conda/mamba environment with fair-esm CLI
        skip = "" if force_reinstall else """# Check if already installed
if conda run -n ESMFoldEnv python -c "import esm" 2>/dev/null; then
    echo "ESMFold already installed, skipping. Use force_reinstall=True to reinstall."
    exit 0
fi
"""
        env_yaml = os.path.join(folders.get("biopipelines", ""), "Environments", "ESMFold.yaml")
        return f"""echo "=== Installing ESMFold ==="
{skip}{env_manager} env create -f {env_yaml} -y
{env_manager} run -n ESMFoldEnv pip install "fair-esm[esmfold]"
{env_manager} run -n ESMFoldEnv pip install 'dllogger @ git+https://github.com/NVIDIA/dllogger.git'
{env_manager} run -n ESMFoldEnv pip install 'openfold @ git+https://github.com/aqlaboratory/openfold.git@4b41059694619831a7db195b7e0988fc4ff3a307'

echo "=== ESMFold installation complete ==="
"""

    def __init__(self,
                 proteins: Union[DataStream, StandardizedOutput],
                 num_recycles: int = 4,
                 chunk_size: Optional[int] = None,
                 max_tokens_per_batch: Optional[int] = None,
                 cpu_offload: bool = False,
                 **kwargs):
        """
        Initialize ESMFold configuration.

        Args:
            proteins: Input protein sequences as DataStream or StandardizedOutput.
                      Use Sequence("MKTVRQ...") to create from raw sequence strings.
            num_recycles: Number of recycles to run (default 4, the training value)
            chunk_size: Chunk size for axial attention to reduce memory.
                        Recommended: 128, 64, or 32 for long sequences.
            max_tokens_per_batch: Maximum tokens per GPU batch. Set to 0 to disable
                                 batching. Lower values reduce memory usage.
            cpu_offload: If True, offload model parameters to CPU RAM.
                         Enables folding of longer sequences with less GPU memory.

        Output:
            Streams: structures (.pdb)
            Tables:
                structures: id | file
                confidence: id | structure | plddt
        """
        # Store original input for upstream missing table lookup
        self.proteins = proteins

        # Resolve input to DataStream
        if isinstance(proteins, StandardizedOutput):
            self.sequences_stream: DataStream = proteins.streams.sequences
        elif isinstance(proteins, DataStream):
            self.sequences_stream = proteins
        else:
            raise ValueError(f"proteins must be DataStream or StandardizedOutput, got {type(proteins)}")

        self.num_recycles = num_recycles
        self.chunk_size = chunk_size
        self.max_tokens_per_batch = max_tokens_per_batch
        self.cpu_offload = cpu_offload

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate ESMFold-specific parameters."""
        if not self.sequences_stream or len(self.sequences_stream) == 0:
            raise ValueError("proteins parameter is required and must not be empty")

        if self.num_recycles < 0:
            raise ValueError("num_recycles cannot be negative")

        if self.chunk_size is not None and self.chunk_size <= 0:
            raise ValueError("chunk_size must be positive")

        if self.max_tokens_per_batch is not None and self.max_tokens_per_batch < 0:
            raise ValueError("max_tokens_per_batch cannot be negative")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input sequences."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get ESMFold configuration display lines."""
        config_lines = super().get_config_display()
        config_lines.extend([
            f"INPUT PROTEINS: {len(self.sequences_stream)} sequences",
            f"NUM RECYCLES: {self.num_recycles}",
        ])
        if self.chunk_size is not None:
            config_lines.append(f"CHUNK SIZE: {self.chunk_size}")
        if self.max_tokens_per_batch is not None:
            config_lines.append(f"MAX TOKENS PER BATCH: {self.max_tokens_per_batch}")
        if self.cpu_offload:
            config_lines.append("CPU OFFLOAD: enabled")
        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate ESMFold execution script."""
        script_content = "#!/bin/bash\n"
        script_content += "# ESMFold execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self._generate_script_prepare_sequences()
        script_content += self._generate_script_run_esmfold()
        script_content += self._generate_script_collect_structures()
        script_content += self._generate_script_update_structures_map()
        script_content += self._generate_script_extract_confidence()
        script_content += self._generate_missing_table_propagation()
        script_content += self.generate_completion_check_footer()

        return script_content

    def _generate_script_prepare_sequences(self) -> str:
        """Generate the sequence preparation part of the script."""
        source_file = self.sequences_stream.map_table

        if source_file.endswith("seqs") and not source_file.endswith(".csv"):
            # ProteinMPNN output - convert .fa files to CSV and FASTA
            return f"""echo "Converting ProteinMPNN .fa files to queries CSV/FASTA"
python {self.fa_to_csv_fasta_py} {source_file} {self.queries_csv} {self.queries_fasta}

"""
        else:
            # Direct CSV - copy and convert to FASTA
            return f"""echo "Using sequences from: {source_file}"
if [ -f "{source_file}" ]; then
    cp "{source_file}" "{self.queries_csv}"
    echo "Successfully copied sequences file"
else
    echo "ERROR: Sequence file not found: {source_file}"
    echo "This usually means the previous step failed to generate the expected output"
    exit 1
fi

# Convert CSV to FASTA for ESMFold input
python -c "
import pandas as pd
df = pd.read_csv('{self.queries_csv}')
with open('{self.queries_fasta}', 'w') as f:
    for _, row in df.iterrows():
        f.write('>'+str(row['id'])+'\\n'+str(row['sequence'])+'\\n')
print('Wrote '+str(len(df))+' sequences to FASTA')
"

"""

    def _generate_script_run_esmfold(self) -> str:
        """Generate the ESMFold execution part of the script."""
        from .config_manager import ConfigManager
        env_manager = ConfigManager().get_env_manager()

        if env_manager in ("pip", "micromamba"):
            # Colab/pip: use HuggingFace Python API via HelpScript
            chunk_arg = f" --chunk-size {self.chunk_size}" if self.chunk_size is not None else ""
            return f"""echo "Running ESMFold (HuggingFace Transformers)"
echo "Output folder: {self.output_folder}"

# Create Folding subfolder for raw ESMFold outputs
mkdir -p "{self.folding_folder}"

# Run ESMFold via Python API
python {self.esm_fold_predict_py} \\
    --fasta "{self.queries_fasta}" \\
    --output "{self.folding_folder}" \\
    --num-recycles {self.num_recycles}{chunk_arg}

"""
        else:
            # Cluster: use fair-esm CLI
            options = f"--num-recycles {self.num_recycles}"
            if self.chunk_size is not None:
                options += f" --chunk-size {self.chunk_size}"
            if self.max_tokens_per_batch is not None:
                options += f" --max-tokens-per-batch {self.max_tokens_per_batch}"
            if self.cpu_offload:
                options += " --cpu-offload"

            return f"""echo "Running ESMFold"
echo "Options: {options}"
echo "Output folder: {self.output_folder}"

# Create Folding subfolder for raw ESMFold outputs
mkdir -p "{self.folding_folder}"

# Run ESMFold
esm-fold -i "{self.queries_fasta}" -o "{self.folding_folder}" {options}

"""

    def _generate_script_collect_structures(self) -> str:
        """Generate script to collect output PDB files from the Folding subfolder."""
        return f"""echo "Collecting predicted structures"
# ESMFold outputs PDB files named after FASTA headers: <id>.pdb
# Copy them to the main output folder

cd "{self.folding_folder}"
for file in *.pdb; do
    if [ -f "$file" ]; then
        cp "$file" "{self.output_folder}/$file"
        echo "Collected: $file"
    fi
done
cd - > /dev/null

"""

    def _generate_script_update_structures_map(self) -> str:
        """Generate script to update structures_map.csv with actual runtime output files."""
        return f"""echo "Updating structures map with actual output files"
python {self.update_map_py} --structures-map "{self.structures_map}" --output-folder "{self.output_folder}"

"""

    def _generate_script_extract_confidence(self) -> str:
        """Generate script to extract pLDDT confidence from PDB B-factor columns."""
        return f"""echo "Extracting confidence metrics from PDB B-factors"
python {self.esm_fold_confidence_py} "{self.output_folder}" "{self.confidence_csv}"

"""

    def _generate_missing_table_propagation(self) -> str:
        """Generate script section to propagate missing.csv from upstream tools."""
        upstream_missing_path = self._get_upstream_missing_table_path(
            self.proteins,
            self.sequences_stream
        )

        if not upstream_missing_path:
            return ""

        upstream_folder = os.path.dirname(upstream_missing_path)

        return f"""
# Propagate missing table from upstream tools
echo "Checking for upstream missing sequences..."
if [ -f "{upstream_missing_path}" ]; then
    echo "Found upstream missing.csv - propagating to current tool"
    python {self.propagate_missing_py} \\
        --upstream-folders "{upstream_folder}" \\
        --output-folder "{self.output_folder}"
else
    echo "No upstream missing.csv found - creating empty missing.csv"
    echo "id,removed_by,cause" > "{self.missing_csv}"
fi

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after ESMFold execution."""
        sequence_ids = list(self.sequences_stream.ids)

        # Structure files (one PDB per sequence)
        structure_files = [os.path.join(self.output_folder, "<id>.pdb")]

        # Create map_table for structures
        create_map_table(self.structures_map, sequence_ids, files=structure_files)

        structures = DataStream(
            name="structures",
            ids=sequence_ids,
            files=structure_files,
            map_table=self.structures_map,
            format="pdb"
        )

        # Tables
        tables = {
            "structures": TableInfo(
                name="structures",
                path=self.structures_map,
                columns=["id", "file"],
                description="ESMFold predicted structures",
                count=len(sequence_ids)
            ),
            "confidence": TableInfo(
                name="confidence",
                path=self.confidence_csv,
                columns=["id", "structure", "plddt"],
                description="ESMFold confidence metrics (mean pLDDT from B-factors)",
                count=len(sequence_ids)
            ),
        }

        # Check for upstream missing table
        upstream_missing_path = self._get_upstream_missing_table_path(
            self.proteins,
            self.sequences_stream
        )
        if upstream_missing_path:
            tables["missing"] = TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "cause"],
                description="IDs removed by upstream tools with removal reason",
                count="variable"
            )

        return {
            "structures": structures,
            "tables": tables,
            "output_folder": self.output_folder,
            "rendering_parameters": {
                "structures": {
                    "color_by": "plddt",
                    "plddt_upper": 100,
                }
            }
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including ESMFold-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "esmfold_params": {
                "num_recycles": self.num_recycles,
                "chunk_size": self.chunk_size,
                "max_tokens_per_batch": self.max_tokens_per_batch,
                "cpu_offload": self.cpu_offload
            }
        })
        return base_dict
