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
        """Generate the MMseqs2 execution part of the script."""

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

    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after MMseqs2 execution.

        Returns:
            Dictionary mapping output types to file paths
        """
        # Ensure file paths are set up
        if not hasattr(self, 'pipeline_name') or self.pipeline_name is None:
            self._setup_file_paths()

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
            "msas": [self.output_msa_csv],
            "msa_ids": [self.pipeline_name],
            "sequences": [],
            "sequence_ids": [],
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
        """Generate the MMseqs2Server execution part of the script."""
        if self.mode == "gpu":
            return self._generate_gpu_server_script()
        else:
            return self._generate_cpu_server_script()

    def _generate_cpu_server_script(self) -> str:
        """Generate CPU server script using dual folder structure."""
        threads_setting = f"export OMP_NUM_THREADS={self.threads}" if self.threads else "export OMP_NUM_THREADS=$(nproc)"

        return f"""echo "Starting MMseqs2 CPU server"
echo "Database: {self.database}"
echo "Max sequences: {self.max_seqs}"
echo "Shared server folder: {self.shared_server_folder}"
echo "Pipeline log folder: {self.output_folder}"

# Configuration - dual folder structure
USER=${{USER:-$(whoami)}}
JOB_QUEUE_DIR="{self.shared_server_folder}/job_queue"
RESULTS_DIR="{self.shared_server_folder}/results"
DB_DIR="/shares/locbp.chem.uzh/models/mmseqs2_databases/cpu"
DB_PATH="$DB_DIR/{self.database}"
POLL_INTERVAL={self.poll_interval}
MAX_SEQS={self.max_seqs}

# Pipeline-specific log in standard location
PIPELINE_LOG_FILE="{self.output_folder}/server.log"

mkdir -p "$JOB_QUEUE_DIR" "$RESULTS_DIR"
mkdir -p "{self.output_folder}"

# Logging setup - logs go to both shared and pipeline locations
LOG_FILE="$RESULTS_DIR/server.log"
: > "$LOG_FILE"
: > "$PIPELINE_LOG_FILE"
{threads_setting}

# Timestamped logging - write to both shared and pipeline logs
timestamp() {{ date '+%Y-%m-%d %H:%M:%S'; }}
log() {{ echo "[$(timestamp)] $*" | tee -a "$LOG_FILE" | tee -a "$PIPELINE_LOG_FILE"; }}

log "MMseqs2 CPU server starting (using $OMP_NUM_THREADS threads)"
log "Database path: $DB_PATH"

# Conversion functions
convert_to_a3m() {{
    local result_db=$1
    local query_db=$2
    local target_db=$3
    local output_file=$4
    local tmp_dir=$5

    mmseqs result2msa "$query_db" "$target_db" "$result_db" "$output_file" \\
        --msa-format-mode 5 \\
        --threads "$OMP_NUM_THREADS" \\
        2>&1 | tee -a "$LOG_FILE"

    # Remove null byte at end
    head -c -1 "$output_file" > "${{output_file}}.clean"
    mv "${{output_file}}.clean" "$output_file"

    local seq_count=$(grep -c "^>" "$output_file" || echo "0")
    log "A3M file created with $seq_count sequences"
}}

convert_a3m_to_csv() {{
    local a3m_file=$1
    local output_file=$2

    echo "key,sequence" > "$output_file"
    awk '/^[^>]/ {{print "-1," $0}}' "$a3m_file" >> "$output_file"

    local line_count=$(wc -l < "$output_file")
    log "CSV file created with $((line_count - 1)) data rows"
}}

handle_job() {{
    local job_meta="$1"
    local fname=$(basename "$job_meta")
    local job_id="${{fname%.job}}"
    local tmp="$RESULTS_DIR/tmp_$job_id"
    mkdir -p "$tmp"

    log "Picked up job $job_id"
    if [[ -f "$job_meta" ]]; then
        mv "$job_meta" "$tmp/"
    else
        log "Warning: job file $job_meta not found; skipping."
        return
    fi

    # Parse params
    output_format="csv"; fasta=""
    while IFS='=' read -r key val; do
        case $key in
            fasta)         fasta="$val"         ;;
            output_format) output_format="$val" ;;
        esac
    done < "$tmp/$fname"

    # Fallback FASTA
    [[ -n "$fasta" ]] || fasta="$JOB_QUEUE_DIR/${{job_id}}.fasta"
    if [[ ! -f "$fasta" ]]; then
        log "ERROR: FASTA not found at $fasta"
        echo -e "FAILED\\nInput FASTA not found" > "$RESULTS_DIR/$job_id.status"
        return
    fi

    cp "$fasta" "$tmp/query.fasta"
    local query_db="$tmp/queryDB"
    log "Creating queryDB"
    mmseqs createdb "$tmp/query.fasta" "$query_db" 2>&1 | tee -a "$LOG_FILE"

    local result_db="$tmp/resultDB"
    local out_prefix="$RESULTS_DIR/$job_id"

    log "Running search for $job_id"
    if ! mmseqs search \\
        "$query_db" "$DB_PATH" "$result_db" "$tmp/search_tmp" \\
        --db-load-mode 2 \\
        -s 7.5 \\
        -a 1 \\
        --alignment-mode 0 \\
        --threads "$OMP_NUM_THREADS" \\
        --max-seqs "$MAX_SEQS" \\
      2>&1 | tee -a "$LOG_FILE"; then
        log "Search failed for $job_id"
        echo -e "FAILED\\nSearch error" > "$RESULTS_DIR/$job_id.status"
        return
    fi

    log "Converting to $output_format format"
    case $output_format in
        a3m)
            convert_to_a3m "$result_db" "$query_db" "$DB_PATH" "$out_prefix.a3m" "$tmp"
            echo -e "SUCCESS\\noutput_file=$out_prefix.a3m" > "$RESULTS_DIR/$job_id.status"
            ;;
        csv)
            convert_to_a3m "$result_db" "$query_db" "$DB_PATH" "$tmp/temp.a3m" "$tmp"
            convert_a3m_to_csv "$tmp/temp.a3m" "$out_prefix.csv"
            rm "$tmp/temp.a3m"
            echo -e "SUCCESS\\noutput_file=$out_prefix.csv" > "$RESULTS_DIR/$job_id.status"
            ;;
        *)
            log "Unknown format $output_format"
            echo -e "FAILED\\nUnknown format" > "$RESULTS_DIR/$job_id.status"
            ;;
    esac

    # Cleanup
    rm -rf "$tmp"
    log "Job $job_id completed successfully"
}}

# Main processing loop
log "Entering job processing loop"
while true; do
  # Process current user's jobs first
  for job_meta in "$JOB_QUEUE_DIR"/*.job; do
    [[ -e "$job_meta" ]] || continue
    if [[ $(stat -c '%U' "$job_meta" 2>/dev/null || stat -f '%Su' "$job_meta") == "$USER" ]]; then
      handle_job "$job_meta"
    fi
  done

  # Then process other users' jobs
  for job_meta in "$JOB_QUEUE_DIR"/*.job; do
    [[ -e "$job_meta" ]] || continue
    if [[ $(stat -c '%U' "$job_meta" 2>/dev/null || stat -f '%Su' "$job_meta") != "$USER" ]]; then
      handle_job "$job_meta"
    fi
  done

  sleep $POLL_INTERVAL
done

"""

    def _generate_gpu_server_script(self) -> str:
        """Generate GPU server script using dual folder structure."""
        threads_setting = f"export OMP_NUM_THREADS={self.threads}" if self.threads else "export OMP_NUM_THREADS=16"

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

# Configuration - dual folder structure
USER=${{USER:-$(whoami)}}
JOB_QUEUE_DIR="{self.shared_server_folder}/job_queue"
RESULTS_DIR="{self.shared_server_folder}/results"
DB_DIR="/shares/locbp.chem.uzh/models/mmseqs2_databases/gpu"
DB_PATH="$DB_DIR/{self.database}"
THREADS={self.threads if self.threads else 32}
POLL_INTERVAL={self.poll_interval}
MAX_SEQS={self.max_seqs}

# Pipeline-specific log in standard location
PIPELINE_LOG_FILE="{self.output_folder}/server.log"

mkdir -p "$JOB_QUEUE_DIR" "$RESULTS_DIR"
mkdir -p "{self.output_folder}"

# Logging setup - logs go to both shared and pipeline locations
LOG_FILE="$RESULTS_DIR/server.log"
: > "$LOG_FILE"
: > "$PIPELINE_LOG_FILE"

# GPU Memory optimization
export MMSEQS_MAX_MEMORY=${{MMSEQS_MAX_MEMORY:-150G}}
{threads_setting}
export CUDA_VISIBLE_DEVICES=0
export CUDA_CACHE_MAXSIZE=2147483648
export CUDA_CACHE_DISABLE=0

# Logging function - write to both shared and pipeline logs
log() {{
  echo "$(date '+%Y-%m-%d %H:%M:%S') - $*" | tee -a "$LOG_FILE" | tee -a "$PIPELINE_LOG_FILE"
}}

log "Memory settings: MMSEQS_MAX_MEMORY=$MMSEQS_MAX_MEMORY, OMP_NUM_THREADS=$OMP_NUM_THREADS"
log "GPU Memory status:"
nvidia-smi --query-gpu=memory.total,memory.used,memory.free --format=csv,noheader,nounits | tee -a "$LOG_FILE" | tee -a "$PIPELINE_LOG_FILE"

# Start GPU server
log "Starting MMseqs2 GPU server for $DB_PATH"
CUDA_VISIBLE_DEVICES=0 mmseqs gpuserver "$DB_PATH" \\
  --max-seqs "$MAX_SEQS" \\
  --db-load-mode 0 \\
  --prefilter-mode 1 &
GPUSERVER_PID=$!
log "GPU server PID=$GPUSERVER_PID"

# Wait for GPU server to stabilize
log "Waiting for gpuserver to finish preloading DB into GPU memory"
prev_mem=-1
while true; do
  curr_mem=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits | head -1)
  if [[ "$curr_mem" -eq "$prev_mem" ]]; then
    log "GPU memory stabilized at ${{curr_mem}} MiB — assuming preload is done"
    break
  else
    log "GPU memory at ${{curr_mem}} MiB (loading…)"
    prev_mem=$curr_mem
    sleep 5
  fi
done

# Cleanup on exit
cleanup() {{
  log "Stopping GPU server PID=$GPUSERVER_PID"
  kill $GPUSERVER_PID || true
  sleep 5
  kill -9 $GPUSERVER_PID 2>/dev/null || true
  exit 0
}}
trap cleanup SIGINT SIGTERM

# GPU server processing loop (simplified - uses gpuserver)
log "Entering job processing loop"
while true; do
  for job_meta in "$JOB_QUEUE_DIR"/*.job; do
    [[ -e "$job_meta" ]] || continue

    fname=$(basename "$job_meta")
    job_id="${{fname%.job}}"
    tmp="$RESULTS_DIR/tmp_$job_id"
    mkdir -p "$tmp"

    log "Picked up job $job_id"
    mv "$job_meta" "$tmp/"

    # Parse params
    output_format="csv"; fasta=""
    while IFS='=' read -r key val; do
      case $key in
        fasta)         fasta="$val"         ;;
        output_format) output_format="$val" ;;
      esac
    done < "$tmp/$fname"

    [[ -n "$fasta" ]] || fasta="$JOB_QUEUE_DIR/${{job_id}}.fasta"
    if [[ ! -f "$fasta" ]]; then
      log "ERROR: FASTA not found at $fasta"
      echo -e "FAILED\\nInput FASTA not found" > "$RESULTS_DIR/$job_id.status"
      continue
    fi

    # Prepare query DB
    cp "$fasta" "$tmp/query.fasta"
    query_db="$tmp/queryDB"
    log "Creating queryDB"
    mmseqs createdb "$tmp/query.fasta" "$query_db" 2>&1 | tee -a "$LOG_FILE"

    result_db="$tmp/resultDB"
    out_prefix="$RESULTS_DIR/$job_id"

    gpu_mem_before=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits | head -1)
    log "GPU memory before search: ${{gpu_mem_before}}MB"

    if ! CUDA_VISIBLE_DEVICES=0 mmseqs search "$query_db" "$DB_PATH" "$result_db" "$tmp" \\
      --gpu 1 \\
      --gpu-server 1 \\
      --prefilter-mode 1 \\
      --db-load-mode 2 \\
      -a 1 \\
      --alignment-mode 0 \\
      --threads "$OMP_NUM_THREADS" \\
      --max-seqs "$MAX_SEQS" \\
      -s 7.5; then
      log "Search failed for job $job_id"
      gpu_mem_after=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits | head -1)
      log "GPU memory after failed search: ${{gpu_mem_after}}MB"
      echo "FAILED: Search failed" > "$RESULTS_DIR/$job_id.status"
      continue
    fi

    gpu_mem_after=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits | head -1)
    log "GPU memory after successful search: ${{gpu_mem_after}}MB"

    log "Converting to $output_format format"
    case $output_format in
      a3m)
        mmseqs result2msa "$query_db" "$DB_PATH" "$result_db" "$out_prefix.a3m" \\
            --msa-format-mode 5 \\
            --threads "$OMP_NUM_THREADS" \\
            2>&1 | tee -a "$LOG_FILE"
        head -c -1 "$out_prefix.a3m" > "${{out_prefix}}.a3m.clean"
        mv "${{out_prefix}}.a3m.clean" "$out_prefix.a3m"
        echo -e "SUCCESS\\noutput_file=$out_prefix.a3m" > "$RESULTS_DIR/$job_id.status"
        ;;
      csv)
        mmseqs result2msa "$query_db" "$DB_PATH" "$result_db" "$tmp/temp.a3m" \\
            --msa-format-mode 5 \\
            --threads "$OMP_NUM_THREADS" \\
            2>&1 | tee -a "$LOG_FILE"
        echo "key,sequence" > "$out_prefix.csv"
        awk '/^[^>]/ {{print "-1," $0}}' "$tmp/temp.a3m" >> "$out_prefix.csv"
        echo -e "SUCCESS\\noutput_file=$out_prefix.csv" > "$RESULTS_DIR/$job_id.status"
        ;;
      *)
        log "Unknown format $output_format"
        echo -e "FAILED\\nUnknown format" > "$RESULTS_DIR/$job_id.status"
        ;;
    esac

    # Cleanup
    rm -rf "$tmp"
    log "Job $job_id completed successfully"
  done
  sleep 5
done

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