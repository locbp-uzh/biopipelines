"""
MMseqs2 configuration for Multiple Sequence Alignment generation.

Provides both direct MMseqs2 search and server-based processing
with comprehensive parameter validation and script generation.
"""

import os
import hashlib
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput


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
    DEFAULT_ENV = "ProteinEnv"  # Uses standard protein environment

    def __init__(self, sequences: Union[str, List[str], ToolOutput, StandardizedOutput],
                 output_format: str = "csv",
                 timeout: int = 3600,
                 **kwargs):
        """
        Initialize MMseqs2 configuration.

        Args:
            sequences: Input sequences - can be sequence string, list, or ToolOutput
            output_format: Output format ("csv" or "a3m", default: csv)
            timeout: Timeout in seconds for server response
            **kwargs: Additional parameters
        """
        # Store MMseqs2-specific parameters
        self.sequences = sequences
        self.sequences_is_tool_output = False
        self.sequences_source_file = None

        # Handle tool output for sequences input
        if isinstance(sequences, (ToolOutput, StandardizedOutput)):
            self.sequences_is_tool_output = True
            if isinstance(sequences, StandardizedOutput):
                # Get sequences from StandardizedOutput
                if hasattr(sequences, 'sequences') and sequences.sequences:
                    self.sequences_source_file = sequences.sequences[0] if isinstance(sequences.sequences, list) else sequences.sequences
                else:
                    raise ValueError("No sequences found in StandardizedOutput")
            else:  # ToolOutput
                # Get sequences from ToolOutput
                sequence_files = sequences.get_output_files("sequences")
                if sequence_files:
                    self.sequences_source_file = sequence_files[0]
                    # Add dependency
                    self.dependencies.append(sequences.config)
                else:
                    raise ValueError("No sequences found in ToolOutput")

        self.output_format = output_format
        self.timeout = timeout

        # Initialize base class
        super().__init__(**kwargs)

        # Initialize file paths
        self._initialize_file_paths()

    def validate_params(self):
        """Validate MMseqs2-specific parameters."""
        if not self.sequences:
            raise ValueError("sequences parameter is required for MMseqs2")

        if self.output_format not in ["csv", "a3m"]:
            raise ValueError("output_format must be 'csv' or 'a3m'")

        if self.timeout <= 0:
            raise ValueError("timeout must be positive")

    def _initialize_file_paths(self):
        """Initialize common file paths used throughout the class."""
        self.input_sequences_csv = None
        self.output_msa_csv = None
        self.pipeline_name = None
        self.client_script_path = None
        self.helper_script_path = None

    def _setup_file_paths(self):
        """Set up all file paths after output_folder is known."""
        # Extract pipeline name from folder structure
        self.pipeline_name = self._extract_pipeline_name()

        # Core output files
        self.output_msa_csv = os.path.join(self.output_folder, "msas.csv")

        # Helper script paths
        if hasattr(self, 'folders') and self.folders:
            self.client_script_path = os.path.join(self.folders["HelpScripts"], "mmseqs2_client.sh")
            self.helper_script_path = os.path.join(self.folders["HelpScripts"], "pipe_mmseqs2_sequences.py")

            # Input sequences file path
            if self.sequences_is_tool_output:
                # Use the actual file path from tool output
                self.input_sequences_csv = self.sequences_source_file
            else:
                # Create sequences CSV file
                self.input_sequences_csv = os.path.join(self.folders["runtime"], "sequences.csv")
        else:
            # Temporary placeholders when folders aren't available yet
            self.client_script_path = None
            self.helper_script_path = None
            self.input_sequences_csv = ""

    def _extract_pipeline_name(self) -> str:
        """Extract pipeline name from output folder structure."""
        folder_parts = self.output_folder.split(os.sep)
        for i, part in enumerate(folder_parts):
            if "MMseqs2" in part:
                if i > 0:
                    return folder_parts[i-1]
                break
        raise ValueError(f"Could not extract pipeline name from output folder: {self.output_folder}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input files and dependencies."""
        self.folders = pipeline_folders
        self._setup_file_paths()

        # Handle sequences input if not from tool output
        if not self.sequences_is_tool_output:
            # Create sequences CSV file from input
            import pandas as pd

            sequences_data = []
            if isinstance(self.sequences, str):
                sequences_data.append({"id": f"{self.pipeline_name}_1", "sequence": self.sequences})
            elif isinstance(self.sequences, list):
                for i, seq in enumerate(self.sequences):
                    sequences_data.append({"id": f"{self.pipeline_name}_{i+1}", "sequence": seq})

            # Create sequences CSV
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
        """
        Generate MMseqs2 execution script.

        Args:
            script_path: Path where script should be written

        Returns:
            Script content as string
        """
        runtime_folder = os.path.dirname(script_path)
        mmseqs_job_folder = self.output_folder
        os.makedirs(mmseqs_job_folder, exist_ok=True)

        # Generate script content
        script_content = "#!/bin/bash\n"
        script_content += "# MMseqs2 client script\n"
        script_content += "# Generated by ProteinNotebooks pipeline system\n\n"
        script_content += self.generate_completion_check_header()
        script_content += self.generate_script_run_mmseqs2()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_mmseqs2(self) -> str:
        """
        Generate the MMseqs2 execution part of the script.

        NOTE: MMseqs2 is an exception to the standard architecture pattern.
        Unlike other tools that generate bash scripts and call pipe_<tool>.py,
        MMseqs2 calls existing bash scripts from HelpScripts (mmseqs2_client.sh)
        and helper python scripts (pipe_mmseqs2_sequences.py) directly.
        This is because MMseqs2 requires interaction with pre-existing server
        infrastructure rather than generating new computational workflows.
        """

        return f"""echo "Starting MMseqs2 MSA generation"
echo "Input sequences: {self.input_sequences_csv}"
echo "Output format: {self.output_format}"
echo "Output MSA CSV: {self.output_msa_csv}"

# Check server status
echo "Checking MMseqs2 server status..."
{self.client_script_path} --status

# Process sequences through MMseqs2 server
echo "Processing sequences through MMseqs2 server..."
python {self.helper_script_path} \\
    "{self.input_sequences_csv}" \\
    "{self.output_msa_csv}" \\
    "{self.client_script_path}" \\
    --output_format {self.output_format}

echo "MMseqs2 processing completed"

"""

    def _predict_sequence_ids(self) -> List[str]:
        """
        Predict sequence IDs from input sources.

        Returns:
            List of expected sequence IDs that will have MSAs generated
        """
        sequence_ids = []

        # Case 1: ToolOutput input (from upstream tools like SDM)
        if hasattr(self.sequences, 'get_output_files'):
            upstream_tool = self.sequences
            output_files = upstream_tool.get_output_files()
            if 'sequence_ids' in output_files and output_files['sequence_ids']:
                sequence_ids = output_files['sequence_ids']

        # Case 2: StandardizedOutput input (input=tool)
        elif hasattr(self, 'standardized_input') and self.standardized_input:
            if hasattr(self.standardized_input, 'sequence_ids') and self.standardized_input.sequence_ids:
                sequence_ids = self.standardized_input.sequence_ids

        # Must have sequence IDs from input sources
        if not sequence_ids:
            raise ValueError("Could not determine sequence IDs - no valid input sequences found")

        return sequence_ids

    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after MMseqs2 execution.

        Returns:
            Dictionary mapping output types to file paths
        """
        # Ensure file paths are set up
        if not hasattr(self, 'pipeline_name') or self.pipeline_name is None:
            self._setup_file_paths()

        # Predict sequence IDs and individual MSA files
        sequence_ids = self._predict_sequence_ids()
        individual_msas = []

        # Generate individual MSA file paths based on sequence IDs
        for seq_id in sequence_ids:
            if self.output_format == "csv":
                msa_file = os.path.join(self.output_folder, f"{seq_id}.csv")
            else:  # a3m
                msa_file = os.path.join(self.output_folder, f"{seq_id}.a3m")
            individual_msas.append(msa_file)

        # Organize datasheets by content type
        datasheets = {
            "msas": {
                "path": self.output_msa_csv,
                "columns": ["id", "sequence_id", "sequence", "msa_file"],
                "description": "MSA files for sequence alignment",
                "count": "variable"
            }
        }

        return {
            "msas": individual_msas,  # Now returns individual MSA files like Boltz2 structures
            "msa_ids": sequence_ids,
            "sequences": [],
            "sequence_ids": sequence_ids,
            "structures": [],
            "structure_ids": [],
            "compounds": [],
            "compound_ids": [],
            "datasheets": datasheets,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including MMseqs2-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "mmseqs2_params": {
                "sequences": str(self.sequences) if not isinstance(self.sequences, (ToolOutput, StandardizedOutput)) else "tool_output",
                "output_format": self.output_format,
                "timeout": self.timeout
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
    DEFAULT_ENV = None  # MMseqs2Server doesn't require conda environment

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

        # Initialize base class
        super().__init__(**kwargs)

        # Initialize file paths
        self._initialize_file_paths()

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

    def _initialize_file_paths(self):
        """Initialize common file paths used throughout the class."""
        self.server_script_path = None
        self.pipeline_name = None
        self.shared_server_folder = None

    def _setup_file_paths(self):
        """Set up dual folder structure for MMseqs2Server."""
        # Extract user from output folder to create shared persistent location
        self.shared_server_folder = self._get_shared_server_folder()

        # Create shared folders for queue and results
        os.makedirs(os.path.join(self.shared_server_folder, "job_queue"), exist_ok=True)
        os.makedirs(os.path.join(self.shared_server_folder, "results"), exist_ok=True)

        # Server script path based on mode (in standard pipeline location for logs)
        if hasattr(self, 'folders') and self.folders:
            if self.mode == "gpu":
                # Create modified GPU server script
                self.server_script_path = os.path.join(self.output_folder, "mmseqs2_server_gpu_modified.sh")
            else:
                # Create modified CPU server script
                self.server_script_path = os.path.join(self.output_folder, "mmseqs2_server_cpu_modified.sh")

    def _get_shared_server_folder(self) -> str:
        """Get shared persistent folder for MMseqs2Server queue and results."""
        # Use global shared location that matches the client script
        # This is shared across all users
        return "/shares/locbp.chem.uzh/models/mmseqs2_server"

    def _extract_pipeline_name(self) -> str:
        """Extract pipeline name from output folder structure."""
        folder_parts = self.output_folder.split(os.sep)
        for i, part in enumerate(folder_parts):
            if "MMseqs2Server" in part:
                if i > 0:
                    return folder_parts[i-1]
                break
        raise ValueError(f"Could not extract pipeline name from output folder: {self.output_folder}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input files and dependencies."""
        self.folders = pipeline_folders
        self._setup_file_paths()

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
        """
        Generate MMseqs2Server execution script.

        Args:
            script_path: Path where script should be written

        Returns:
            Script content as string
        """
        runtime_folder = os.path.dirname(script_path)
        server_job_folder = self.output_folder
        os.makedirs(server_job_folder, exist_ok=True)

        # Generate script content
        script_content = "#!/bin/bash\n"
        script_content += f"# MMseqs2Server {self.mode.upper()} script\n"
        script_content += "# Generated by ProteinNotebooks pipeline system\n\n"
        script_content += self.generate_completion_check_header()
        script_content += self.generate_script_run_server()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_server(self) -> str:
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
        """
        Generate CPU server script by calling existing HelpScripts bash script.

        Rather than generating the entire server logic inline, this calls the
        pre-existing mmseqs2_server_cpu.sh script from HelpScripts with
        appropriate environment variables set.
        """
        # Set environment variables for the server script
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
            f"export MMSEQS2_PIPELINE_LOG={self.output_folder}/server.log"
        ])

        env_setup = "\n".join(env_vars)
        cpu_script_path = os.path.join(self.folders["HelpScripts"], "mmseqs2_server_cpu.sh")

        return f"""echo "Starting MMseqs2 CPU server"
echo "Database: {self.database}"
echo "Max sequences: {self.max_seqs}"
echo "Shared server folder: {self.shared_server_folder}"
echo "Pipeline log folder: {self.output_folder}"

# Set environment variables for server script
{env_setup}

# Call existing CPU server script from HelpScripts
echo "Executing MMseqs2 CPU server script..."
bash {cpu_script_path}

"""

    def _generate_gpu_server_script(self) -> str:
        """
        Generate GPU server script by calling existing HelpScripts bash script.

        Rather than generating the entire server logic inline, this calls the
        pre-existing mmseqs2_server_gpu.sh script from HelpScripts with
        appropriate environment variables set.
        """
        # Set environment variables for the server script
        env_vars = []
        if self.threads:
            env_vars.append(f"export OMP_NUM_THREADS={self.threads}")
        else:
            env_vars.append("export OMP_NUM_THREADS=16")

        env_vars.extend([
            f"export MMSEQS2_DATABASE={self.database}",
            f"export MMSEQS2_MAX_SEQS={self.max_seqs}",
            f"export MMSEQS2_POLL_INTERVAL={self.poll_interval}",
            f"export MMSEQS2_SHARED_FOLDER={self.shared_server_folder}",
            f"export MMSEQS2_PIPELINE_LOG={self.output_folder}/server.log",
            "export CUDA_VISIBLE_DEVICES=0",
            "export CUDA_CACHE_MAXSIZE=2147483648",
            "export CUDA_CACHE_DISABLE=0"
        ])

        env_setup = "\n".join(env_vars)
        gpu_script_path = os.path.join(self.folders["HelpScripts"], "mmseqs2_server_gpu.sh")

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
bash {gpu_script_path}

"""

    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after MMseqs2Server execution.

        Returns:
            Dictionary mapping output types to file paths
        """
        # Server doesn't produce output files - it just runs
        return {
            "msas": [],
            "msa_ids": [],
            "sequences": [],
            "sequence_ids": [],
            "structures": [],
            "structure_ids": [],
            "compounds": [],
            "compound_ids": [],
            "datasheets": {},
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