# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
MMseqs2 configuration for Multiple Sequence Alignment generation.

Provides both direct MMseqs2 search and server-based processing
with comprehensive parameter validation and script generation.
"""

import os
import json
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


class MMseqs2(BaseConfig):
    """
    Configuration for MMseqs2 client - processes sequences through server.

    Submits sequences to a running MMseqs2 server for MSA generation.

    ARCHITECTURE NOTE: MMseqs2 is an exception to the standard BioPipelines pattern.
    Unlike other tools that generate bash scripts and call pipe_<tool>.py at SLURM runtime,
    MMseqs2 calls existing bash scripts from HelpScripts (mmseqs2_client.sh) and helper
    python scripts (pipe_mmseqs2_sequences.py) directly. This is necessary because MMseqs2
    requires interaction with pre-existing server infrastructure rather than generating
    new computational workflows.
    """

    TOOL_NAME = "MMseqs2"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba"):
        return """echo "=== MMseqs2 ==="
echo "Uses biopipelines environment (no additional installation needed)."
echo "=== MMseqs2 ready ==="
"""

    # Lazy path descriptors
    output_msa_csv = Path(lambda self: os.path.join(self.output_folder, "msas.csv"))
    client_script = Path(lambda self: os.path.join(self.folders["HelpScripts"], "mmseqs2_client.sh"))
    helper_script = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_mmseqs2_sequences.py"))
    input_sequences_csv = Path(lambda self: self._get_input_sequences_path())

    def __init__(self, sequences: Union[str, List[str], DataStream, StandardizedOutput],
                 output_format: str = "csv",
                 timeout: int = 3600,
                 mask: Union[str, tuple] = "",
                 id_map: Dict[str, str] = None,
                 **kwargs):
        """
        Initialize MMseqs2 configuration.

        Args:
            sequences: Input sequences - can be sequence string, list, DataStream or StandardizedOutput
            output_format: Output format ("csv" or "a3m", default: csv)
            timeout: Timeout in seconds for server response
            mask: Positions to mask in MSA (excluding query sequence)
                  - String format: "10-20+30-40" (PyMOL selection style)
                  - Tuple format: (TableInfo, "column_name") for per-sequence masking
                  - Empty string: no masking (default)
            id_map: ID mapping pattern for matching sequence IDs to table IDs (default: {"*": "*_<N>"})
                  - Used when mask table IDs don't match sequence IDs
                  - Example: sequence ID "rifampicin_1_2" maps to table ID "rifampicin_1"
                  - Pattern {"*": "*_<N>"} strips last "_<number>" from sequence ID
            **kwargs: Additional parameters
        """
        if id_map is None:
            id_map = {"*": "*_<N>"}

        # Resolve sequences input
        self.sequences_source_file: Optional[str] = None
        self.sequences_stream: Optional[DataStream] = None
        self.raw_sequences: Optional[Union[str, List[str]]] = None

        if isinstance(sequences, StandardizedOutput):
            if sequences.streams.sequences and len(sequences.streams.sequences) > 0:
                self.sequences_stream = sequences.streams.sequences
                self.sequences_source_file = sequences.streams.sequences.map_table
            else:
                raise ValueError("No sequences found in StandardizedOutput")
        elif isinstance(sequences, DataStream):
            self.sequences_stream = sequences
            self.sequences_source_file = sequences.map_table
        elif isinstance(sequences, (str, list)):
            self.raw_sequences = sequences
        else:
            raise ValueError(f"sequences must be str, list, DataStream or StandardizedOutput, got {type(sequences)}")

        self.output_format = output_format
        self.timeout = timeout
        self.mask_positions = mask
        self.id_map = id_map

        super().__init__(**kwargs)

    def _get_input_sequences_path(self) -> str:
        """Get path to input sequences CSV file."""
        if self.sequences_source_file:
            return self.sequences_source_file
        return os.path.join(self.output_folder, "input_sequences.csv")

    def validate_params(self):
        """Validate MMseqs2-specific parameters."""
        if self.sequences_stream is None and self.raw_sequences is None:
            raise ValueError("sequences parameter is required for MMseqs2")

        if self.output_format not in ["csv", "a3m"]:
            raise ValueError("output_format must be 'csv' or 'a3m'")

        if self.timeout <= 0:
            raise ValueError("timeout must be positive")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input files and dependencies."""
        self.folders = pipeline_folders

        # Handle raw sequences input - write to CSV file
        if self.raw_sequences is not None:
            import pandas as pd

            sequences_data = []
            if isinstance(self.raw_sequences, str):
                sequences_data.append({"id": f"{self.pipeline_name}_1", "sequence": self.raw_sequences})
            elif isinstance(self.raw_sequences, list):
                for i, seq in enumerate(self.raw_sequences):
                    sequences_data.append({"id": f"{self.pipeline_name}_{i+1}", "sequence": seq})

            os.makedirs(self.output_folder, exist_ok=True)
            df = pd.DataFrame(sequences_data)
            df.to_csv(self.input_sequences_csv, index=False)

    def get_config_display(self) -> List[str]:
        """Get MMseqs2 configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"OUTPUT FORMAT: {self.output_format}",
            f"TIMEOUT: {self.timeout}s"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate MMseqs2 execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        script_content = "#!/bin/bash\n"
        script_content += "# MMseqs2 client script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self._generate_script_run_mmseqs2()
        script_content += self.generate_completion_check_footer()

        return script_content

    def _generate_script_run_mmseqs2(self) -> str:
        """
        Generate the MMseqs2 execution part of the script.

        NOTE: MMseqs2 is an exception to the standard architecture pattern.
        Unlike other tools that generate bash scripts and call pipe_<tool>.py,
        MMseqs2 calls existing bash scripts from HelpScripts (mmseqs2_client.sh)
        and helper python scripts (pipe_mmseqs2_sequences.py) directly.
        This is because MMseqs2 requires interaction with pre-existing server
        infrastructure rather than generating new computational workflows.
        """
        server_dir = self.folders.get("MMseqs2Server", "")

        return f"""echo "Starting MMseqs2 MSA generation"
echo "Input sequences: {self.input_sequences_csv}"
echo "Output format: {self.output_format}"
echo "Output MSA CSV: {self.output_msa_csv}"

# Process sequences through MMseqs2 server
# Note: Server checking and submission is handled by pipe_mmseqs2_sequences.py
echo "Processing sequences through MMseqs2 server..."
python {self.helper_script} \\
    "{self.input_sequences_csv}" \\
    "{self.output_msa_csv}" \\
    "{self.client_script}" \\
    --output_format {self.output_format} \\
    --server_dir "{server_dir}"{self._generate_mask_arguments()}

echo "MMseqs2 processing completed"

"""

    def _generate_mask_arguments(self) -> str:
        """Generate mask-related command line arguments for the helper script."""
        if not self.mask_positions:
            return ""

        id_map_json = json.dumps(self.id_map).replace('"', '\\"')

        # Handle tuple format: (TableInfo, "column_name")
        if isinstance(self.mask_positions, tuple):
            if len(self.mask_positions) == 2:
                table_info, column_name = self.mask_positions
                if hasattr(table_info, 'info'):
                    return f' \\\n    --mask_table "{table_info.info.path}" \\\n    --mask_column "{column_name}" \\\n    --id_map "{id_map_json}"'
                else:
                    raise ValueError(f"Invalid table reference in mask parameter: {self.mask_positions}")
            else:
                raise ValueError(f"Invalid tuple format for mask parameter: {self.mask_positions}")

        # Handle string format: direct selection like "10-20+30-40"
        elif isinstance(self.mask_positions, str):
            return f' \\\n    --mask_selection "{self.mask_positions}"'

        else:
            raise ValueError(f"Unsupported mask parameter type: {type(self.mask_positions)}")

    def _predict_sequence_ids(self) -> List[str]:
        """Predict sequence IDs from input sources."""
        # Get from sequences stream if available
        if self.sequences_stream and len(self.sequences_stream) > 0:
            return self.sequences_stream.ids

        # Handle raw sequences (string or list)
        if self.raw_sequences is not None:
            if isinstance(self.raw_sequences, str):
                return [f"{self.pipeline_name}_1"]
            elif isinstance(self.raw_sequences, list):
                return [f"{self.pipeline_name}_{i+1}" for i in range(len(self.raw_sequences))]

        raise ValueError("Could not determine sequence IDs - no valid input sequences found")

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after MMseqs2 execution."""
        sequence_ids = self._predict_sequence_ids()

        # Generate individual MSA file paths based on sequence IDs
        msa_files = []
        ext = "csv" if self.output_format == "csv" else "a3m"
        for seq_id in sequence_ids:
            msa_file = os.path.join(self.output_folder, f"{seq_id}.{ext}")
            msa_files.append(msa_file)

        msas = DataStream(
            name="msas",
            ids=sequence_ids,
            files=msa_files,
            map_table=self.output_msa_csv,
            format=ext
        )

        tables = {
            "msas": TableInfo(
                name="msas",
                path=self.output_msa_csv,
                columns=["id", "sequence_id", "sequence", "msa_file"],
                description="MSA files for sequence alignment",
                count=len(sequence_ids)
            )
        }

        return {
            "structures": DataStream.empty("structures", "pdb"),
            "sequences": DataStream.empty("sequences", "fasta"),
            "compounds": DataStream.empty("compounds", "sdf"),
            "msas": msas,
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including MMseqs2-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "mmseqs2_params": {
                "sequences_source": self.sequences_source_file if self.sequences_source_file else "raw_input",
                "output_format": self.output_format,
                "timeout": self.timeout,
                "mask_positions": str(self.mask_positions) if self.mask_positions else None,
                "id_map": self.id_map
            }
        })
        return base_dict


class MMseqs2Server(BaseConfig):
    """
    Configuration for MMseqs2 server management.

    Starts and manages MMseqs2 server processes with CPU or GPU modes.
    Does not process sequences - only manages server infrastructure.

    ARCHITECTURE NOTE: MMseqs2Server is an exception to the standard BioPipelines pattern.
    Unlike other tools that generate bash scripts and call pipe_<tool>.py at SLURM runtime,
    MMseqs2Server calls existing bash scripts from HelpScripts (mmseqs2_server_cpu.sh,
    mmseqs2_server_gpu.sh) directly. This is necessary because MMseqs2Server manages
    pre-existing server infrastructure rather than generating new computational workflows.
    """

    TOOL_NAME = "MMseqs2Server"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba"):
        return """echo "=== MMseqs2Server ==="
echo "Uses biopipelines environment (no additional installation needed)."
echo "=== MMseqs2Server ready ==="
"""

    # Lazy path descriptors
    cpu_server_script = Path(lambda self: os.path.join(self.folders["HelpScripts"], "mmseqs2_server_cpu.sh"))
    gpu_server_script = Path(lambda self: os.path.join(self.folders["HelpScripts"], "mmseqs2_server_gpu.sh"))
    shared_server_folder = Path(lambda self: self.folders.get("MMseqs2Server", ""))

    def __init__(self, mode: str = "cpu",
                 database: str = "uniref30_2302_db",
                 max_seqs: int = 10000,
                 threads: int = None,
                 poll_interval: int = 10,
                 **kwargs):
        """
        Initialize MMseqs2Server configuration.

        Args:
            mode: Server mode ("cpu" or "gpu")
            database: Database to use (default: uniref30_2302_db)
            max_seqs: Maximum sequences to return per query
            threads: Number of threads (auto-detect if None)
            poll_interval: Job polling interval in seconds
            **kwargs: Additional parameters
        """
        self.mode = mode
        self.database = database
        self.max_seqs = max_seqs
        self.threads = threads
        self.poll_interval = poll_interval

        # Set mode-specific default resources
        if mode == "cpu":
            mode_resources = {"gpu": "none", "memory": "16GB", "time": "24:00:00"}
        else:  # gpu mode
            mode_resources = {"gpu": "V100", "memory": "32GB", "time": "24:00:00"}

        # Merge mode-specific resources with any user-provided resources
        if 'resources' in kwargs:
            kwargs['resources'] = {**mode_resources, **kwargs['resources']}
        else:
            kwargs['resources'] = mode_resources

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate MMseqs2Server-specific parameters."""
        if self.mode not in ["cpu", "gpu"]:
            raise ValueError("mode must be 'cpu' or 'gpu'")

        if self.max_seqs <= 0:
            raise ValueError("max_seqs must be positive")

        if self.threads is not None and self.threads <= 0:
            raise ValueError("threads must be positive")

        if self.poll_interval <= 0:
            raise ValueError("poll_interval must be positive")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input files and dependencies."""
        self.folders = pipeline_folders

        # Create shared folders for queue and results
        os.makedirs(os.path.join(self.shared_server_folder, "job_queue"), exist_ok=True)
        os.makedirs(os.path.join(self.shared_server_folder, "results"), exist_ok=True)

    def get_config_display(self) -> List[str]:
        """Get MMseqs2Server configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"MODE: {self.mode.upper()}",
            f"DATABASE: {self.database}",
            f"MAX SEQUENCES: {self.max_seqs}",
            f"POLL INTERVAL: {self.poll_interval}s"
        ])

        if self.threads:
            config_lines.append(f"THREADS: {self.threads}")

        # Show dual folder structure
        if hasattr(self, 'shared_server_folder') and self.shared_server_folder:
            config_lines.extend([
                f"SHARED QUEUE/RESULTS: {self.shared_server_folder}",
                f"PIPELINE LOGS: {self.output_folder}"
            ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate MMseqs2Server execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        script_content = "#!/bin/bash\n"
        script_content += f"# MMseqs2Server {self.mode.upper()} script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self._generate_script_run_server()
        script_content += self.generate_completion_check_footer()

        return script_content

    def _generate_script_run_server(self) -> str:
        """
        Generate the MMseqs2Server execution part of the script.

        NOTE: MMseqs2Server is an exception to the standard architecture pattern.
        Unlike other tools that generate bash scripts and call pipe_<tool>.py,
        MMseqs2Server calls existing bash scripts from HelpScripts directly.
        This is because MMseqs2Server manages pre-existing server infrastructure
        rather than generating new computational workflows.
        """
        if self.mode == "gpu":
            return self._generate_gpu_server_script()
        else:
            return self._generate_cpu_server_script()

    def _generate_cpu_server_script(self) -> str:
        """Generate CPU server script by calling existing HelpScripts bash script."""
        env_vars = []
        if self.threads:
            env_vars.append(f"export OMP_NUM_THREADS={self.threads}")
        else:
            env_vars.append("export OMP_NUM_THREADS=$(nproc)")

        env_vars.extend([
            f"export MMSEQS2_DATABASE={self.database}",
            f"export MMSEQS2_MAX_SEQS={self.max_seqs}",
            f"export MMSEQS2_POLL_INTERVAL={self.poll_interval}",
            f"export MMSEQS2_SHARED_FOLDER={self.shared_server_folder}",
            f"export MMSEQS2_PIPELINE_LOG={self.output_folder}/server.log",
            f"export BIOPIPELINES_DATA_DIR={self.folders.get('data', '')}",
            f"export MMSEQS2_DIR={self.folders.get('MMseqs2', '')}"
        ])

        env_setup = "\n".join(env_vars)

        return f"""echo "Starting MMseqs2 CPU server"
echo "Database: {self.database}"
echo "Max sequences: {self.max_seqs}"
echo "Shared server folder: {self.shared_server_folder}"
echo "Pipeline log folder: {self.output_folder}"

# Set environment variables for server script
{env_setup}

# Call existing CPU server script from HelpScripts
echo "Executing MMseqs2 CPU server script..."
bash {self.cpu_server_script}

"""

    def _generate_gpu_server_script(self) -> str:
        """Generate GPU server script by calling existing HelpScripts bash script."""
        env_vars = []
        if self.threads:
            env_vars.append(f"export OMP_NUM_THREADS={self.threads}")
        else:
            env_vars.append("export OMP_NUM_THREADS=4")

        env_vars.extend([
            f"export MMSEQS2_DATABASE={self.database}",
            f"export MMSEQS2_MAX_SEQS={self.max_seqs}",
            f"export MMSEQS2_POLL_INTERVAL={self.poll_interval}",
            f"export MMSEQS2_SHARED_FOLDER={self.shared_server_folder}",
            f"export MMSEQS2_PIPELINE_LOG={self.output_folder}/server.log",
            f"export MMSEQS2_DB_DIR={self.folders.get('MMseqs2Databases', '')}",
            f"export BIOPIPELINES_DATA_DIR={self.folders.get('data', '')}",
            f"export MMSEQS2_DIR={self.folders.get('MMseqs2', '')}",
            "export CUDA_VISIBLE_DEVICES=0",
            "export CUDA_CACHE_MAXSIZE=2147483648",
            "export CUDA_CACHE_DISABLE=0"
        ])

        env_setup = "\n".join(env_vars)

        return f"""echo "Starting MMseqs2 GPU server"
echo "Database: {self.database}"
echo "Max sequences: {self.max_seqs}"
echo "Shared server folder: {self.shared_server_folder}"
echo "Pipeline log folder: {self.output_folder}"

# Check GPU availability
if ! command -v nvidia-smi &> /dev/null; then
    echo "Could not load GPU correctly: nvidia-smi could not be found"
    exit 1
fi

# Display GPU information
gpu_type=$(nvidia-smi --query-gpu=gpu_name --format=csv,noheader)
echo "GPU Type: $gpu_type"
nvidia-smi --query-gpu=memory.total,memory.used,memory.free --format=csv,noheader,nounits

# Set environment variables for server script
{env_setup}

# Call existing GPU server script from HelpScripts
echo "Executing MMseqs2 GPU server script..."
bash {self.gpu_server_script}

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after MMseqs2Server execution."""
        # Server doesn't produce output files - it just runs
        return {
            "structures": DataStream.empty("structures", "pdb"),
            "sequences": DataStream.empty("sequences", "fasta"),
            "compounds": DataStream.empty("compounds", "sdf"),
            "msas": DataStream.empty("msas", "a3m"),
            "tables": {},
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including MMseqs2Server-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "mmseqs2_server_params": {
                "mode": self.mode,
                "database": self.database,
                "max_seqs": self.max_seqs,
                "threads": self.threads,
                "poll_interval": self.poll_interval
            }
        })
        return base_dict