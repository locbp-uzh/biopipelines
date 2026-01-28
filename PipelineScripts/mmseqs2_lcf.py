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
    

    def __init__(self, sequences: Union[str, List[str], ToolOutput, StandardizedOutput],
                 timeout: int = 3600,
                 mask: Union[str, tuple] = "",
                 id_map: Dict[str, str] = {"*": "*_<N>"},
                 **kwargs):
        """
        Initialize MMseqs2LCF configuration.

        Args:
            sequences: Input sequences - can be sequence string, list, or ToolOutput
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
        # Store MMseqs2LCF-specific parameters
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

        self.timeout = timeout
        self.mask_positions = mask
        self.id_map = id_map

        # Initialize base class
        super().__init__(**kwargs)

        # Initialize file paths
        self._initialize_file_paths()

    def validate_params(self):
        """Validate MMseqs2LCF-specific parameters."""
        if not self.sequences:
            raise ValueError("sequences parameter is required for MMseqs2LCF")

        if self.timeout <= 0:
            raise ValueError("timeout must be positive")

        # Validate mask parameter if provided
        if self.mask_positions:
            self.validate_table_reference(self.mask_positions)

    def _initialize_file_paths(self):
        """Initialize common file paths used throughout the class."""
        self.input_sequences_csv = None
        self.output_msa_csv = None
        self.pipeline_name = None
        self.client_script_path = None
        self.helper_script_path = None
        self.server_dir = None

    def _setup_file_paths(self):
        """Set up all file paths after output_folder is known."""
        # Extract pipeline name from folder structure
        self.pipeline_name = self._extract_pipeline_name()

        # Core output files
        self.output_msa_csv = os.path.join(self.output_folder, "msas.csv")

        # Helper script paths
        if hasattr(self, 'folders') and self.folders:
            self.client_script_path = os.path.join(self.folders["HelpScripts"], "mmseqs2_lcf_client.sh")
            self.helper_script_path = os.path.join(self.folders["HelpScripts"], "pipe_mmseqs2_lcf_sequences.py")

            # Server directory from folders configuration
            self.server_dir = self.folders["MMseqs2LCFServer"]

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
            self.server_dir = None

    def _extract_pipeline_name(self) -> str:
        """Extract pipeline name from output folder structure."""
        folder_parts = self.output_folder.split(os.sep)
        for i, part in enumerate(folder_parts):
            if "MMseqs2LCF" in part:
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
        """Get MMseqs2LCF configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"OUTPUT FORMAT: a3m (native ColabFold)",
            f"DATABASES: UniRef30 + ColabFoldDB",
            f"TIMEOUT: {self.timeout}s"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """
        Generate MMseqs2LCF execution script.

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
        script_content += "# MMseqs2 LCF client script\n"
        script_content += "# Generated by ProteinNotebooks pipeline system\n\n"
        script_content += self.generate_completion_check_header()
        script_content += self.generate_script_run_mmseqs2_lcf()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_mmseqs2_lcf(self) -> str:
        """
        Generate the MMseqs2LCF execution part of the script.

        NOTE: MMseqs2LCF is an exception to the standard architecture pattern.
        Unlike other tools that generate bash scripts and call pipe_<tool>.py,
        MMseqs2LCF calls existing bash scripts from HelpScripts (mmseqs2_lcf_client.sh)
        and helper python scripts (pipe_mmseqs2_lcf_sequences.py) directly.
        This is because MMseqs2LCF requires interaction with pre-existing server
        infrastructure rather than generating new computational workflows.
        """

        return f"""echo "Starting MMseqs2 LCF MSA generation"
echo "Input sequences: {self.input_sequences_csv}"
echo "Output format: a3m (native ColabFold)"
echo "Output MSA CSV: {self.output_msa_csv}"
echo "Databases: UniRef30 + ColabFoldDB"

# Process sequences through MMseqs2 LCF server
# Note: Server checking and submission is handled by pipe_mmseqs2_lcf_sequences.py
echo "Processing sequences through MMseqs2 LCF server..."
python {self.helper_script_path} \\
    "{self.input_sequences_csv}" \\
    "{self.output_msa_csv}" \\
    "{self.client_script_path}" \\
    --server_dir "{self.server_dir}"{self._generate_mask_arguments()}

echo "MMseqs2 LCF processing completed"

"""

    def _generate_mask_arguments(self) -> str:
        """
        Generate mask-related command line arguments for the helper script.

        Returns:
            String with mask arguments to append to python command
        """
        if not self.mask_positions:
            return ""

        # Convert id_map to JSON string for passing to script
        import json
        id_map_json = json.dumps(self.id_map).replace('"', '\\"')

        # Handle tuple format: (TableInfo, "column_name")
        if isinstance(self.mask_positions, tuple):
            if len(self.mask_positions) == 2:
                table_info, column_name = self.mask_positions
                if hasattr(table_info, 'path'):
                    # Per-sequence masking from table
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
        """
        Predict sequence IDs from input sources.

        Returns:
            List of expected sequence IDs that will have MSAs generated
        """
        # Try to get from standardized input first (highest priority)
        if hasattr(self, 'standardized_input') and self.standardized_input:
            if hasattr(self.standardized_input, 'sequence_ids') and self.standardized_input.sequence_ids:
                return self.standardized_input.sequence_ids

        # Try to get from sequences parameter if it's StandardizedOutput
        if hasattr(self.sequences, 'sequence_ids') and self.sequences.sequence_ids:
            return self.sequences.sequence_ids

        # Try to extract from dependencies (like SDM tool)
        if hasattr(self.sequences, 'config') and hasattr(self.sequences.config, '_predict_sequence_ids'):
            return self.sequences.config._predict_sequence_ids()

        # Handle raw sequences (string or list) - predict IDs based on pattern used in configure_inputs
        if not self.sequences_is_tool_output:
            # Need pipeline name to generate IDs
            if not hasattr(self, 'pipeline_name') or self.pipeline_name is None:
                # Try to extract from output_folder if available
                if hasattr(self, 'output_folder') and self.output_folder:
                    self.pipeline_name = self._extract_pipeline_name()
                else:
                    # Can't predict yet - return empty list (will be populated later)
                    return []

            if isinstance(self.sequences, str):
                return [f"{self.pipeline_name}_1"]
            elif isinstance(self.sequences, list):
                return [f"{self.pipeline_name}_{i+1}" for i in range(len(self.sequences))]

        # Must have sequence IDs from input sources
        raise ValueError("Could not determine sequence IDs - no valid input sequences found")

    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after MMseqs2LCF execution.

        Returns:
            Dictionary mapping output types to file paths
        """
        # Ensure file paths are set up
        if not hasattr(self, 'pipeline_name') or self.pipeline_name is None:
            self._setup_file_paths()

        # Predict sequence IDs and individual MSA files
        sequence_ids = self._predict_sequence_ids()
        individual_msas = []

        # Generate individual MSA file paths based on sequence IDs (A3M format)
        for seq_id in sequence_ids:
            msa_file = os.path.join(self.output_folder, f"{seq_id}.a3m")
            individual_msas.append(msa_file)

        # Organize tables by content type
        tables = {
            "msas": {
                "path": self.output_msa_csv,
                "columns": ["id", "sequence_id", "sequence", "msa_file"],
                "description": "MSA files for sequence alignment (UniRef30 + ColabFoldDB)",
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
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including MMseqs2LCF-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "mmseqs2_lcf_params": {
                "sequences": str(self.sequences) if not isinstance(self.sequences, (ToolOutput, StandardizedOutput)) else "tool_output",
                "timeout": self.timeout
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

        # Initialize base class
        super().__init__(**kwargs)

        # Initialize file paths
        self._initialize_file_paths()

    def validate_params(self):
        """Validate MMseqs2ServerLCF-specific parameters."""
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
        """Set up dual folder structure for MMseqs2ServerLCF."""
        # Extract user from output folder to create shared persistent location
        self.shared_server_folder = self._get_shared_server_folder()

        # Create shared folders for queue and results
        os.makedirs(os.path.join(self.shared_server_folder, "job_queue"), exist_ok=True)
        os.makedirs(os.path.join(self.shared_server_folder, "results"), exist_ok=True)

        # Server script path (in standard pipeline location for logs)
        if hasattr(self, 'folders') and self.folders:
            self.server_script_path = os.path.join(self.folders["HelpScripts"], "mmseqs2_lcf_server_gpu.sh")

    def _get_shared_server_folder(self) -> str:
        """Get shared persistent folder for MMseqs2ServerLCF queue and results."""
        # Use per-user location that matches the client script
        user = os.environ.get('USER', os.environ.get('USERNAME', 'unknown'))
        return f"/shares/locbp.chem.uzh/{user}/BioPipelines/MMseqs2LCFServer"

    def _extract_pipeline_name(self) -> str:
        """Extract pipeline name from output folder structure."""
        folder_parts = self.output_folder.split(os.sep)
        for i, part in enumerate(folder_parts):
            if "MMseqs2ServerLCF" in part:
                if i > 0:
                    return folder_parts[i-1]
                break
        raise ValueError(f"Could not extract pipeline name from output folder: {self.output_folder}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input files and dependencies."""
        self.folders = pipeline_folders
        self._setup_file_paths()

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

        # Show dual folder structure
        if hasattr(self, 'shared_server_folder') and self.shared_server_folder:
            config_lines.extend([
                f"SHARED QUEUE/RESULTS: {self.shared_server_folder}",
                f"PIPELINE LOGS: {self.output_folder}"
            ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """
        Generate MMseqs2ServerLCF execution script.

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
        script_content += f"# MMseqs2ServerLCF GPU script\n"
        script_content += "# Generated by ProteinNotebooks pipeline system\n\n"
        script_content += self.generate_completion_check_header()
        script_content += self.generate_script_run_server()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_server(self) -> str:
        """
        Generate the MMseqs2ServerLCF execution part of the script.

        NOTE: MMseqs2ServerLCF is an exception to the standard architecture pattern.
        Unlike other tools that generate bash scripts and call pipe_<tool>.py,
        MMseqs2ServerLCF calls existing bash scripts from HelpScripts directly.
        This is because MMseqs2ServerLCF manages pre-existing server infrastructure
        rather than generating new computational workflows.
        """
        # Set environment variables for the server script
        env_vars = []
        if self.threads:
            env_vars.append(f"export OMP_NUM_THREADS={self.threads}")
        else:
            env_vars.append("export OMP_NUM_THREADS=16")

        env_vars.extend([
            f"export MMSEQS2_POLL_INTERVAL={self.poll_interval}",
            f"export MMSEQS2_SHARED_FOLDER={self.shared_server_folder}",
            f"export MMSEQS2_PIPELINE_LOG={self.output_folder}/server.log",
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
bash {self.server_script_path}

"""

    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after MMseqs2ServerLCF execution.

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
