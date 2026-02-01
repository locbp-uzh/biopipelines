"""
MMseqs2 LCF (LocalColabFold) configuration for MSA generation using colabfold_search.

Provides GPU-accelerated MSA generation using colabfold_search, which searches
both UniRef30 and ColabFoldDB (environmental sequences) for richer alignments.

Key differences from standard MMseqs2:
- Uses colabfold_search command from localcolabfold installation
- Searches both UniRef30 and ColabFoldDB for more comprehensive MSAs
- Outputs A3M format (native to ColabFold)
- Requires localcolabfold installation with databases
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


class MMseqs2LCF(BaseConfig):
    """
    Configuration for MMseqs2 LCF client - processes sequences through LCF server.

    Submits sequences to a running MMseqs2ServerLCF for MSA generation using
    colabfold_search, which provides richer MSAs by searching both UniRef30
    and ColabFoldDB.

    ARCHITECTURE NOTE: MMseqs2LCF is an exception to the standard BioPipelines pattern.
    Unlike other tools that generate bash scripts and call pipe_<tool>.py at SLURM runtime,
    MMseqs2LCF calls existing bash scripts from HelpScripts (mmseqs2_lcf_client.sh) and helper
    python scripts (pipe_mmseqs2_lcf_sequences.py) directly. This is necessary because MMseqs2LCF
    requires interaction with pre-existing server infrastructure rather than generating
    new computational workflows.
    """

    TOOL_NAME = "MMseqs2LCF"

    # Lazy path descriptors
    output_msa_csv = Path(lambda self: os.path.join(self.output_folder, "msas.csv"))
    client_script = Path(lambda self: os.path.join(self.folders["HelpScripts"], "mmseqs2_lcf_client.sh"))
    helper_script = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_mmseqs2_lcf_sequences.py"))
    input_sequences_csv = Path(lambda self: self._get_input_sequences_path())

    def __init__(self, sequences: Union[str, List[str], DataStream, StandardizedOutput],
                 timeout: int = 3600,
                 mask: Union[str, tuple] = "",
                 id_map: Dict[str, str] = None,
                 **kwargs):
        """
        Initialize MMseqs2LCF configuration.

        Args:
            sequences: Input sequences - can be sequence string, list, DataStream or StandardizedOutput
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
            if sequences.sequences and len(sequences.sequences) > 0:
                self.sequences_stream = sequences.sequences
                self.sequences_source_file = sequences.sequences.map_table
            else:
                raise ValueError("No sequences found in StandardizedOutput")
        elif isinstance(sequences, DataStream):
            self.sequences_stream = sequences
            self.sequences_source_file = sequences.map_table
        elif isinstance(sequences, (str, list)):
            self.raw_sequences = sequences
        else:
            raise ValueError(f"sequences must be str, list, DataStream or StandardizedOutput, got {type(sequences)}")

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
        """Validate MMseqs2LCF-specific parameters."""
        if self.sequences_stream is None and self.raw_sequences is None:
            raise ValueError("sequences parameter is required for MMseqs2LCF")

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
        """Get MMseqs2LCF configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"OUTPUT FORMAT: a3m (native ColabFold)",
            f"DATABASES: UniRef30 + ColabFoldDB",
            f"TIMEOUT: {self.timeout}s"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate MMseqs2LCF execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        script_content = "#!/bin/bash\n"
        script_content += "# MMseqs2 LCF client script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self._generate_script_run_mmseqs2_lcf()
        script_content += self.generate_completion_check_footer()

        return script_content

    def _generate_script_run_mmseqs2_lcf(self) -> str:
        """
        Generate the MMseqs2LCF execution part of the script.

        NOTE: MMseqs2LCF is an exception to the standard architecture pattern.
        Unlike other tools that generate bash scripts and call pipe_<tool>.py,
        MMseqs2LCF calls existing bash scripts from HelpScripts (mmseqs2_lcf_client.sh)
        and helper python scripts (pipe_mmseqs2_lcf_sequences.py) directly.
        This is because MMseqs2LCF requires interaction with pre-existing server
        infrastructure rather than generating new computational workflows.
        """
        server_dir = self.folders.get("MMseqs2LCFServer", "")

        return f"""echo "Starting MMseqs2 LCF MSA generation"
echo "Input sequences: {self.input_sequences_csv}"
echo "Output format: a3m (native ColabFold)"
echo "Output MSA CSV: {self.output_msa_csv}"
echo "Databases: UniRef30 + ColabFoldDB"

# Process sequences through MMseqs2 LCF server
# Note: Server checking and submission is handled by pipe_mmseqs2_lcf_sequences.py
echo "Processing sequences through MMseqs2 LCF server..."
python {self.helper_script} \\
    "{self.input_sequences_csv}" \\
    "{self.output_msa_csv}" \\
    "{self.client_script}" \\
    --server_dir "{server_dir}"{self._generate_mask_arguments()}

echo "MMseqs2 LCF processing completed"

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
                if hasattr(table_info, 'path'):
                    return f' \\\n    --mask_table "{table_info.path}" \\\n    --mask_column "{column_name}" \\\n    --id_map "{id_map_json}"'
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
        """Get expected output files after MMseqs2LCF execution."""
        sequence_ids = self._predict_sequence_ids()

        # Generate individual MSA file paths based on sequence IDs (A3M format)
        msa_files = []
        for seq_id in sequence_ids:
            msa_file = os.path.join(self.output_folder, f"{seq_id}.a3m")
            msa_files.append(msa_file)

        msas = DataStream(
            name="msas",
            ids=sequence_ids,
            files=msa_files,
            map_table=self.output_msa_csv,
            format="a3m"
        )

        tables = {
            "msas": TableInfo(
                name="msas",
                path=self.output_msa_csv,
                columns=["id", "sequence_id", "sequence", "msa_file"],
                description="MSA files for sequence alignment (UniRef30 + ColabFoldDB)",
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
        """Serialize configuration including MMseqs2LCF-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "mmseqs2_lcf_params": {
                "sequences_source": self.sequences_source_file if self.sequences_source_file else "raw_input",
                "timeout": self.timeout,
                "mask_positions": str(self.mask_positions) if self.mask_positions else None,
                "id_map": self.id_map
            }
        })
        return base_dict


class MMseqs2ServerLCF(BaseConfig):
    """
    Configuration for MMseqs2 LCF server management.

    Starts and manages MMseqs2 LCF server processes using colabfold_search
    with GPU-accelerated searches against UniRef30 and ColabFoldDB.
    Does not process sequences - only manages server infrastructure.

    ARCHITECTURE NOTE: MMseqs2ServerLCF is an exception to the standard BioPipelines pattern.
    Unlike other tools that generate bash scripts and call pipe_<tool>.py at SLURM runtime,
    MMseqs2ServerLCF calls existing bash scripts from HelpScripts (mmseqs2_lcf_server_gpu.sh)
    directly. This is necessary because MMseqs2ServerLCF manages pre-existing server
    infrastructure rather than generating new computational workflows.
    """

    TOOL_NAME = "MMseqs2ServerLCF"

    # Lazy path descriptors
    server_script = Path(lambda self: os.path.join(self.folders["HelpScripts"], "mmseqs2_lcf_server_gpu.sh"))
    shared_server_folder = Path(lambda self: self._get_shared_server_folder())

    def __init__(self,
                 threads: int = None,
                 poll_interval: int = 10,
                 **kwargs):
        """
        Initialize MMseqs2ServerLCF configuration.

        Args:
            threads: Number of threads (auto-detect if None)
            poll_interval: Job polling interval in seconds
            **kwargs: Additional parameters
        """
        self.threads = threads
        self.poll_interval = poll_interval

        # Set GPU-specific default resources
        mode_resources = {"gpu": "80GB", "memory": "32GB", "time": "4:00:00"}

        # Merge mode-specific resources with any user-provided resources
        if 'resources' in kwargs:
            kwargs['resources'] = {**mode_resources, **kwargs['resources']}
        else:
            kwargs['resources'] = mode_resources

        super().__init__(**kwargs)

    def _get_shared_server_folder(self) -> str:
        """Get shared persistent folder for MMseqs2ServerLCF queue and results."""
        user = os.environ.get('USER', os.environ.get('USERNAME', 'unknown'))
        return f"/shares/locbp.chem.uzh/{user}/BioPipelines/MMseqs2LCFServer"

    def validate_params(self):
        """Validate MMseqs2ServerLCF-specific parameters."""
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
        """Get MMseqs2ServerLCF configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"MODE: GPU (colabfold_search)",
            f"DATABASES: UniRef30 + ColabFoldDB",
            f"POLL INTERVAL: {self.poll_interval}s"
        ])

        if self.threads:
            config_lines.append(f"THREADS: {self.threads}")

        config_lines.extend([
            f"SHARED QUEUE/RESULTS: {self.shared_server_folder}",
            f"PIPELINE LOGS: {self.output_folder}"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate MMseqs2ServerLCF execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        script_content = "#!/bin/bash\n"
        script_content += f"# MMseqs2ServerLCF GPU script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self._generate_script_run_server()
        script_content += self.generate_completion_check_footer()

        return script_content

    def _generate_script_run_server(self) -> str:
        """Generate the MMseqs2ServerLCF execution part of the script."""
        env_vars = []
        if self.threads:
            env_vars.append(f"export OMP_NUM_THREADS={self.threads}")
        else:
            env_vars.append("export OMP_NUM_THREADS=16")

        env_vars.extend([
            f"export MMSEQS2_POLL_INTERVAL={self.poll_interval}",
            f"export MMSEQS2_SHARED_FOLDER={self.shared_server_folder}",
            f"export MMSEQS2_PIPELINE_LOG={self.output_folder}/server.log",
            f"export COLABFOLD_DB_DIR={self.folders.get('ColabFoldDatabases', '')}",
            "export CUDA_VISIBLE_DEVICES=0",
            "export CUDA_CACHE_MAXSIZE=2147483648",
            "export CUDA_CACHE_DISABLE=0"
        ])

        env_setup = "\n".join(env_vars)

        return f"""echo "Starting MMseqs2 LCF GPU server"
echo "Databases: UniRef30 + ColabFoldDB"
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
echo "Executing MMseqs2 LCF GPU server script..."
bash {self.server_script}

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after MMseqs2ServerLCF execution."""
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
        """Serialize configuration including MMseqs2ServerLCF-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "mmseqs2_server_lcf_params": {
                "threads": self.threads,
                "poll_interval": self.poll_interval
            }
        })
        return base_dict
