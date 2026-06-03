# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
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
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .config_manager import ConfigManager
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from config_manager import ConfigManager


def _configured_server_mode() -> str:
    """Server mode the client would auto-start: tool_overrides.mmseqs2server.mode.

    Mirrors the resolution in pipe_mmseqs2_sequences.py (default 'cpu'). Used to
    reject GPU server mode at config time, since the GPU server path is currently
    unsupported.
    """
    try:
        overrides = ConfigManager()._config.get('tool_overrides', {}) or {}
        mode = (overrides.get('mmseqs2server', {}) or {}).get('mode', 'cpu')
    except Exception:
        return 'cpu'
    return str(mode).strip().lower()


class MMseqs2(BaseConfig):
    """
    Configuration for MMseqs2 client - processes sequences through server.

    Submits sequences to a running MMseqs2 server for MSA generation.

    ARCHITECTURE NOTE: MMseqs2 is an exception to the standard BioPipelines pattern.
    Unlike other tools that generate bash scripts and call pipe_<tool>.py at execution time,
    MMseqs2 calls existing bash scripts from pipe_scripts (mmseqs2_client.sh) and helper
    python scripts (pipe_mmseqs2_sequences.py) directly. This is necessary because MMseqs2
    requires interaction with pre-existing server infrastructure rather than generating
    new computational workflows.
    """

    TOOL_NAME = "MMseqs2"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== MMseqs2 ==="
echo "Uses biopipelines environment (no additional installation needed)."
touch "$INSTALL_SUCCESS"
echo "=== MMseqs2 ready ==="
"""

    # Lazy path descriptors — MSA CSV is the msas stream map_table.
    output_msa_csv = Path(lambda self: self.stream_map_path("msas"))
    client_script = Path(lambda self: self.pipe_script_path("mmseqs2_client.sh"))
    helper_script = Path(lambda self: self.pipe_script_path("pipe_mmseqs2_sequences.py"))
    input_sequences_csv = Path(lambda self: self._get_input_sequences_path())
    missing_csv = Path(lambda self: self.table_path("missing"))

    def __init__(self, sequences: Union[str, List[str], DataStream, StandardizedOutput],
                 output_format: str = "csv",
                 timeout: int = 3600,
                 mask: Union[str, tuple] = "",
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
            **kwargs: Additional parameters

        Output:
            Streams: msas (.csv/.a3m)
            Tables:
                msas: id | sequences.id | sequence | msa_file
        """
        # Resolve sequences input
        self.sequences_input = sequences  # retained for upstream missing-table lookup
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

        super().__init__(**kwargs)

    def _get_input_sequences_path(self) -> str:
        """Get path to input sequences CSV file."""
        if self.sequences_source_file:
            return self.sequences_source_file
        # Raw-sequences path: place the synthesized CSV under configuration/.
        return self.configuration_path("input_sequences.csv")

    def validate_params(self):
        """Validate MMseqs2-specific parameters."""
        if _configured_server_mode() == "gpu":
            raise NotImplementedError(
                "tool_overrides.mmseqs2server.mode is set to 'gpu', but the MMseqs2 "
                "GPU server is not currently supported. Set it to 'cpu' (the CPU "
                "server serves the RAM-resident ColabFold databases)."
            )

        if self.sequences_stream is None and self.raw_sequences is None:
            raise ValueError("sequences parameter is required for MMseqs2")

        if self.output_format not in ["csv", "a3m"]:
            raise ValueError("output_format must be 'csv' or 'a3m'")

        if self.timeout <= 0:
            raise ValueError("timeout must be positive")

        if isinstance(self.mask_positions, str):
            _validate_freeform_string("mask", self.mask_positions)

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

            # This runs in configure_inputs (before pipeline's layout step),
            # so mkdir the config folder explicitly here.
            os.makedirs(self.configuration_folder, exist_ok=True)
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
        MMseqs2 calls existing bash scripts from pipe_scripts (mmseqs2_client.sh)
        and helper python scripts (pipe_mmseqs2_sequences.py) directly.
        This is because MMseqs2 requires interaction with pre-existing server
        infrastructure rather than generating new computational workflows.
        """
        server_dir = self.folders.get("MMseqs2Server", "")

        # Sequences the upstream tool already dropped: forward its missing.csv so
        # the helper seeds those rows and excludes those ids before searching.
        upstream_missing_path = self._get_upstream_missing_table_path(
            self.sequences_input
        )
        upstream_flag = (
            f' \\\n    --upstream_missing "{upstream_missing_path}"'
            if upstream_missing_path else ""
        )

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
    --server_dir "{server_dir}" \\
    --missing_csv "{self.missing_csv}"{upstream_flag}{self._generate_mask_arguments()}

echo "MMseqs2 processing completed"

"""

    def _generate_mask_arguments(self) -> str:
        """Generate mask-related command line arguments for the helper script."""
        if not self.mask_positions:
            return ""

        # Handle TableReference format: table.column_name
        if hasattr(self.mask_positions, 'path') and hasattr(self.mask_positions, 'column'):
            return f' \\\n    --mask_table "{self.mask_positions.path}" \\\n    --mask_column "{self.mask_positions.column}"'

        # Handle string format: direct selection like "10-20+30-40"
        elif isinstance(self.mask_positions, str):
            return f' \\\n    --mask_selection "{self.mask_positions}"'

        else:
            raise ValueError(f"Unsupported mask parameter type: {type(self.mask_positions)}")

    def _predict_sequence_ids(self) -> List[str]:
        """Predict sequence IDs from input sources."""
        # Get from sequences stream if available
        if self.sequences_stream and len(self.sequences_stream) > 0:
            return list(self.sequences_stream.ids)

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
        ext = "csv" if self.output_format == "csv" else "a3m"

        # One <id> template (not a concrete path per id) so the completion check
        # expands it against the resolved ids at runtime — lazy ids stay lazy.
        msas = DataStream(
            name="msas",
            ids=sequence_ids,
            files=[self.stream_path("msas", f"<id>.{ext}")],
            map_table=self.output_msa_csv,
            format=ext
        )

        tables = {
            "msas": TableInfo(
                name="msas",
                path=self.output_msa_csv,
                columns=["id", "sequences.id", "sequence", "msa_file"],
                description="MSA files for sequence alignment"
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "kind", "cause"],
                description="Sequences with no MSA (search failure) plus any propagated from upstream"
            ),
        }

        return {
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
                "mask_positions": str(self.mask_positions) if self.mask_positions else None
            }
        })
        return base_dict


class MMseqs2Server(BaseConfig):
    """
    Configuration for MMseqs2 server management.

    Starts and manages MMseqs2 server processes with CPU or GPU modes.
    Does not process sequences - only manages server infrastructure.

    ARCHITECTURE NOTE: MMseqs2Server is an exception to the standard BioPipelines pattern.
    Unlike other tools that generate bash scripts and call pipe_<tool>.py at execution time,
    MMseqs2Server calls existing bash scripts from pipe_scripts (mmseqs2_server_cpu.sh,
    mmseqs2_server_gpu.sh) directly. This is necessary because MMseqs2Server manages
    pre-existing server infrastructure rather than generating new computational workflows.
    """

    TOOL_NAME = "MMseqs2Server"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False,
                        step=None, mode="gpu", **kwargs):
        """Install bash for the ColabFold MSA databases.

        Args:
            step: Which stage to run.
                None (default) — no-op: the server runs in the biopipelines env
                    and needs nothing installed, so just touch the marker.
                "databases" — download stage only: fetch the four ColabFold
                    tarballs + rsync the mmCIF snapshot into the ColabFoldDatabases
                    folder root (DOWNLOADS_ONLY=1). Mode-independent; shared by both
                    the cpu and gpu builds. No DB indexes are built.
                "build" — build the indexed databases from the already downloaded
                    tarballs into a mode subfolder (see `mode`). The downloads stay
                    at the ColabFoldDatabases root and are symlinked into the
                    subfolder so they aren't duplicated per mode.
            mode: Build flavour for step="build" (default "gpu").
                "gpu" — GPU-padded build (tsv2exprofiledb --gpu 1, makepaddedseqdb,
                    --index-subset 2) into <ColabFoldDatabases>/gpu/. Requires a GPU
                    node, so the user must call Resources(gpu=...) first; ensures a
                    GPU-capable MMseqs2 (release >=16, with gpuserver) is present.
                "cpu" — non-padded build (no GPU flags) into
                    <ColabFoldDatabases>/cpu/. Runs on a CPU node; uses the MMseqs2
                    already at the configured MMseqs2 folder (the AVX2 build the CPU
                    server installs is sufficient).

        The build is idempotent via setup_databases.sh's marker files
        (UNIREF30_READY, COLABDB_READY, PDB_READY, PDB100_READY) which live inside
        the mode subfolder, so re-running resumes where it left off. The download
        markers (DOWNLOADS_READY, PDB_MMCIF_READY) live at the root.
        """
        if step is None:
            # The server itself needs no install, but the CPU server runs
            # colabfold_search + the mmseqs bundled with LocalColabFold (under the
            # AlphaFold folder). Verify both are present; if not, point the user at
            # AlphaFold.install(). `folders` is {} during the override-probe in
            # base_config, so resolve defensively.
            cf_bin = os.path.join(folders.get("AlphaFold", ""), "colabfold-conda", "bin")
            return f"""echo "=== MMseqs2Server ==="
echo "The server needs no install of its own, but it runs colabfold_search and"
echo "the mmseqs bundled with LocalColabFold."
CF_BIN="{cf_bin}"
if [ -x "$CF_BIN/colabfold_search" ] && [ -x "$CF_BIN/mmseqs" ]; then
    echo "Found LocalColabFold (colabfold_search + mmseqs) at $CF_BIN."
    echo "Pass step=\\"databases\\" to download the ColabFold DBs, or step=\\"build\\" (mode=\\"cpu\\"|\\"gpu\\") to build the indexes."
    touch "$INSTALL_SUCCESS"
    echo "=== MMseqs2Server ready ==="
else
    echo "ERROR: LocalColabFold not found at $CF_BIN (need colabfold_search + mmseqs)."
    echo "Run AlphaFold.install() first - it installs LocalColabFold, which the"
    echo "MMseqs2 server uses for both colabfold_search and its CPU mmseqs binary."
    exit 1
fi
"""

        if step not in ("databases", "build"):
            raise ValueError(
                f"MMseqs2Server.install: step must be None, 'databases', or 'build', got {step!r}"
            )
        if step == "build" and mode not in ("cpu", "gpu"):
            raise ValueError(
                f"MMseqs2Server.install: mode must be 'cpu' or 'gpu', got {mode!r}"
            )

        # Trust the config: a missing folder key raises KeyError at config time,
        # which is the framework's intended "no fallbacks" behavior.
        db_dir = folders["ColabFoldDatabases"]
        mmseqs_dir = folders["MMseqs2"]
        setup_script = os.path.join(folders["pipe_scripts"], "colabfold_setup_databases.sh")

        if step == "databases":
            return f"""echo "=== MMseqs2Server: downloading ColabFold databases ==="
DB_DIR="{db_dir}"
mkdir -p "$DB_DIR"
echo "Target: $DB_DIR"
echo "Running download-only stage (tarballs + mmCIF rsync, no DB build)..."
DOWNLOADS_ONLY=1 bash "{setup_script}" "$DB_DIR"
echo "=== ColabFold database download complete ==="
touch "$INSTALL_SUCCESS"
"""

        # step == "build": build into a per-mode subfolder, sharing the downloads.
        # The downloads (tarballs, pdb/ mmCIF, download markers) stay at DB_DIR root
        # and are symlinked into DB_DIR/<mode>/ so setup_databases.sh — which builds
        # in its WORKDIR and reads the tarballs from there — finds them without
        # re-downloading. Only the built *_db* indexes end up under DB_DIR/<mode>/.
        link_downloads = f"""# Symlink the shared downloads into the mode subfolder so the upstream
# setup script (which builds in its WORKDIR) finds them without re-downloading.
# Tarballs (uniref30/envdb/pdb100-foldseek) and the pdb100 fasta (createdb input).
for f in "$DB_DIR"/*.tar.gz "$DB_DIR"/*.fasta.gz; do
    [ -e "$f" ] && ln -sf "$f" "$BUILD_DIR/$(basename "$f")"
done
for marker in DOWNLOADS_READY PDB_MMCIF_READY; do
    [ -e "$DB_DIR/$marker" ] && ln -sf "$DB_DIR/$marker" "$BUILD_DIR/$marker"
done
# mmCIF snapshot: symlink the whole pdb/ tree so PDB_MMCIF_READY stays valid.
[ -e "$DB_DIR/pdb" ] && ln -sfn "$DB_DIR/pdb" "$BUILD_DIR/pdb"
"""

        if mode == "cpu":
            return f"""echo "=== MMseqs2Server: building CPU (non-padded) ColabFold databases ==="
DB_DIR="{db_dir}"
MMSEQS_DIR="{mmseqs_dir}"
BUILD_DIR="$DB_DIR/cpu"
mkdir -p "$BUILD_DIR"

if [ ! -f "$DB_DIR/DOWNLOADS_READY" ]; then
    echo "ERROR: downloads not present at $DB_DIR (no DOWNLOADS_READY)."
    echo "Run MMseqs2Server.install(step=\\"databases\\") first."
    exit 1
fi

# Ensure an MMseqs2 binary is present (AVX2 CPU build is sufficient here).
MMSEQS_BIN="$MMSEQS_DIR/bin/mmseqs"
if [ ! -x "$MMSEQS_BIN" ]; then
    echo "MMseqs2 not found at $MMSEQS_BIN, downloading AVX2 build..."
    PARENT="$(dirname "$MMSEQS_DIR")"
    mkdir -p "$PARENT"
    cd "$PARENT"
    wget -q https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz
    tar xzf mmseqs-linux-avx2.tar.gz
    rm -f mmseqs-linux-avx2.tar.gz
fi
if [ ! -x "$MMSEQS_BIN" ]; then
    echo "ERROR: failed to obtain MMseqs2 at $MMSEQS_BIN"
    exit 1
fi
echo "Using MMseqs2: $MMSEQS_BIN ($("$MMSEQS_BIN" version 2>/dev/null))"
export PATH="$MMSEQS_DIR/bin:$PATH"

# createindex's indexdb child loads the whole DB into RAM (~273G for UniRef30
# alone) and OOM-kills an under-provisioned node. Cap it well below the node's
# RAM so indexdb splits into bounded passes instead. Derived from the SLURM
# memory allocation when available, else a safe default.
if [ -n "${{SLURM_MEM_PER_NODE:-}}" ]; then
    # SLURM_MEM_PER_NODE is in MB; use ~70% of it, floored at 32G.
    LIMIT_GB=$(( SLURM_MEM_PER_NODE * 70 / 100 / 1024 ))
    [ "$LIMIT_GB" -lt 32 ] && LIMIT_GB=32
    export MMSEQS_SPLIT_MEMORY_LIMIT="${{LIMIT_GB}}G"
else
    export MMSEQS_SPLIT_MEMORY_LIMIT="${{MMSEQS_SPLIT_MEMORY_LIMIT:-180G}}"
fi
echo "createindex split-memory-limit: $MMSEQS_SPLIT_MEMORY_LIMIT"

{link_downloads}
echo "Building non-padded indexes in $BUILD_DIR (resumes via marker files)..."
if ! bash "{setup_script}" "$BUILD_DIR"; then
    echo "ERROR: ColabFold CPU database build failed (see above)."
    exit 1
fi
echo "=== ColabFold CPU database build complete ==="
touch "$INSTALL_SUCCESS"
"""

        # mode == "gpu": GPU-padded index build into DB_DIR/gpu/.
        return f"""echo "=== MMseqs2Server: building GPU-padded ColabFold databases ==="
DB_DIR="{db_dir}"
MMSEQS_DIR="{mmseqs_dir}"
BUILD_DIR="$DB_DIR/gpu"
mkdir -p "$BUILD_DIR"

if [ ! -f "$DB_DIR/DOWNLOADS_READY" ]; then
    echo "ERROR: downloads not present at $DB_DIR (no DOWNLOADS_READY)."
    echo "Run MMseqs2Server.install(step=\\"databases\\") first."
    exit 1
fi

# A GPU build needs an actual GPU. The user must set Resources(gpu=...) before
# .install(step="build"); fail loudly rather than silently building CPU DBs.
if ! command -v nvidia-smi >/dev/null 2>&1; then
    echo "ERROR: mode=\\"gpu\\" build needs a GPU node, but nvidia-smi is not available."
    echo "Call Resources(gpu=\\"A100\\") before MMseqs2Server.install(step=\\"build\\", mode=\\"gpu\\")."
    exit 1
fi
nvidia-smi --query-gpu=gpu_name --format=csv,noheader || true

# Ensure a GPU-capable MMseqs2 (release >=16, with gpuserver) is present.
# The colabfold-conda mmseqs is v15 (no GPU), so we install the GPU build here.
MMSEQS_BIN="$MMSEQS_DIR/bin/mmseqs"
if [ ! -x "$MMSEQS_BIN" ] || ! "$MMSEQS_BIN" --help 2>/dev/null | grep -q gpuserver; then
    echo "GPU-capable MMseqs2 not found at $MMSEQS_BIN, downloading..."
    PARENT="$(dirname "$MMSEQS_DIR")"
    mkdir -p "$PARENT"
    cd "$PARENT"
    wget -q https://mmseqs.com/latest/mmseqs-linux-gpu.tar.gz
    tar xzf mmseqs-linux-gpu.tar.gz
    rm -f mmseqs-linux-gpu.tar.gz
fi
if [ ! -x "$MMSEQS_BIN" ] || ! "$MMSEQS_BIN" --help 2>/dev/null | grep -q gpuserver; then
    echo "ERROR: failed to obtain a GPU-capable MMseqs2 at $MMSEQS_BIN"
    exit 1
fi
echo "Using MMseqs2: $MMSEQS_BIN ($("$MMSEQS_BIN" version 2>/dev/null))"

# Put the GPU mmseqs first on PATH so setup_databases.sh's bare `mmseqs` calls
# hit it (the script invokes mmseqs unqualified).
export PATH="$MMSEQS_DIR/bin:$PATH"

{link_downloads}
echo "Building GPU-padded indexes in $BUILD_DIR (resumes via marker files)..."
if ! GPU=1 bash "{setup_script}" "$BUILD_DIR"; then
    echo "ERROR: ColabFold GPU database build failed (see above)."
    exit 1
fi
echo "=== ColabFold GPU database build complete ==="
touch "$INSTALL_SUCCESS"
"""

    # Lazy path descriptors
    cpu_server_script = Path(lambda self: self.pipe_script_path("mmseqs2_server_cpu.sh"))
    gpu_server_script = Path(lambda self: self.pipe_script_path("mmseqs2_server_gpu.sh"))
    shared_server_folder = Path(lambda self: self.folders.get("MMseqs2Server", ""))

    def __init__(self, mode: str = "cpu",
                 database: str = "uniref30_2302_db",
                 max_seqs: int = 10000,
                 threads: int = None,
                 poll_interval: int = 10,
                 gpus: int = 1,
                 idle_timeout: int = 600,
                 **kwargs):
        """
        Initialize MMseqs2Server configuration.

        Args:
            mode: Server mode ("cpu" or "gpu")
            database: Database to use (default: uniref30_2302_db)
            max_seqs: Maximum sequences to return per query
            threads: Number of threads (auto-detect if None)
            poll_interval: Job polling interval in seconds
            gpus: Number of GPUs for the GPU server (1 or 2). With 2, the UniRef30
                  and environmental gpuservers are pinned to separate GPUs
                  (CUDA_VISIBLE_DEVICES 0 and 1) so their prefilters don't share
                  one device; with 1, both run on GPU 0. Only meaningful for
                  mode="gpu".
            idle_timeout: Shut the server down after this many seconds with no
                  jobs in the queue (0 disables). A long-lived idle GPU server
                  squats the card and wrecks SLURM priority, so the server exits
                  once the queue has been empty for this long. Default 600 (10 min)
                  — re-initialising (locking the DB indices into RAM) takes ~4–6 min,
                  so 10 min idle is a small overhang over the re-warm cost.
            **kwargs: Additional parameters
        """
        self.mode = mode
        self.database = database
        self.max_seqs = max_seqs
        self.threads = threads
        self.poll_interval = poll_interval
        self.gpus = gpus
        self.idle_timeout = idle_timeout

        # Set mode-specific default resources. CPU mode runs the same
        # colabfold_search pipeline as GPU mode but keeps the DB indices resident
        # in the OS page cache (warmed at startup) instead of on the GPU. The CPU
        # indices are bigger than the GPU ones (~746GB vs ~410GB) because CPU mode
        # keeps the k-mer prefilter table the GPU build drops (--index-subset 2),
        # so it needs a ~1TB-class node to stay fully resident; the search itself
        # is CPU-bound and benefits from many cores.
        if mode == "cpu":
            mode_resources = {"gpu": "none", "memory": "800GB", "time": "24:00:00", "cpus": 32}
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

        if self.mode == "gpu":
            raise NotImplementedError(
                "MMseqs2Server GPU mode is not currently supported; use mode='cpu' "
                "(the CPU server serves the RAM-resident ColabFold databases). The "
                "GPU server path is retained but disabled."
            )

        if self.max_seqs <= 0:
            raise ValueError("max_seqs must be positive")

        if self.threads is not None and self.threads <= 0:
            raise ValueError("threads must be positive")

        if self.poll_interval <= 0:
            raise ValueError("poll_interval must be positive")

        if self.gpus not in (1, 2):
            raise ValueError("gpus must be 1 or 2")

        if self.idle_timeout < 0:
            raise ValueError("idle_timeout must be >= 0")

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
        MMseqs2Server calls existing bash scripts from pipe_scripts directly.
        This is because MMseqs2Server manages pre-existing server infrastructure
        rather than generating new computational workflows.
        """
        if self.mode == "gpu":
            return self._generate_gpu_server_script()
        else:
            return self._generate_cpu_server_script()

    def _generate_cpu_server_script(self) -> str:
        """Generate CPU server script by calling existing pipe_scripts bash script."""
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
            # Built databases live in the cpu/ subfolder of the DB root (see
            # MMseqs2Server.install(step="build", mode="cpu")).
            f"export MMSEQS2_DB_DIR={os.path.join(self.folders.get('MMseqs2Databases', ''), 'cpu')}",
            f"export BIOPIPELINES_DATA_DIR={self.folders.get('data', '')}",
            # LocalColabFold install: provides colabfold_search AND the CPU mmseqs
            # binary the server uses (colabfold-conda/bin). No separate MMseqs2
            # install needed for CPU mode — keeps it portable to CPU-only clusters.
            f"export COLABFOLD_DIR={self.folders.get('AlphaFold', '')}",
            # Auto-shutdown after this many idle seconds (0 disables) so an idle
            # large-RAM server doesn't squat the node and hurt SLURM priority.
            f"export MMSEQS2_IDLE_TIMEOUT={self.idle_timeout}"
        ])

        env_setup = "\n".join(env_vars)

        return f"""echo "Starting MMseqs2 CPU server"
echo "Database: {self.database}"
echo "Max sequences: {self.max_seqs}"
echo "Shared server folder: {self.shared_server_folder}"
echo "Pipeline log folder: {self.output_folder}"

# Set environment variables for server script
{env_setup}

# Call existing CPU server script from pipe_scripts
echo "Executing MMseqs2 CPU server script..."
bash {self.cpu_server_script}

"""

    def _generate_gpu_server_script(self) -> str:
        """Generate GPU server script by calling existing pipe_scripts bash script."""
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
            # Built databases live in the gpu/ subfolder of the DB root (see
            # MMseqs2Server.install(step="build", mode="gpu")).
            f"export MMSEQS2_DB_DIR={os.path.join(self.folders.get('MMseqs2Databases', ''), 'gpu')}",
            f"export BIOPIPELINES_DATA_DIR={self.folders.get('data', '')}",
            f"export MMSEQS2_DIR={self.folders.get('MMseqs2', '')}",
            # LocalColabFold install (ships colabfold_search); the server adds
            # its colabfold-conda/bin to PATH.
            f"export COLABFOLD_DIR={self.folders.get('AlphaFold', '')}",
            # Number of GPUs (1 or 2). With 2, the script pins the UniRef30 and
            # envdb gpuservers to separate devices.
            f"export MMSEQS2_GPUS={self.gpus}",
            f"export CUDA_VISIBLE_DEVICES={','.join(str(i) for i in range(self.gpus))}",
            # Auto-shutdown after this many idle seconds (0 disables) so an idle
            # GPU server doesn't squat the card and hurt SLURM priority.
            f"export MMSEQS2_IDLE_TIMEOUT={self.idle_timeout}",
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

# Call existing GPU server script from pipe_scripts
echo "Executing MMseqs2 GPU server script..."
bash {self.gpu_server_script}

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after MMseqs2Server execution."""
        # Server doesn't produce output files - it just runs
        return {
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