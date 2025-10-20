"""
ESMFold configuration for fast protein structure prediction.

Uses Meta's ESM-2 with ESMFold for efficient single-sequence structure prediction
without requiring MSAs. Models are loaded via torch.hub and cached in shared folder.
"""

import os
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo
except ImportError:
    # Fallback for direct execution
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo


class ESMFold(BaseConfig):
    """
    ESMFold configuration for fast protein structure prediction.

    Predicts protein structures from sequences using ESM-2 + ESMFold.
    No MSAs required - fast single-sequence prediction.
    """

    TOOL_NAME = "ESMFold"
    DEFAULT_ENV = "ProteinEnv"

    def __init__(self,
                 sequences: Union[str, List[str], ToolOutput, StandardizedOutput],
                 datasheets: Optional[List[str]] = None,
                 name: str = "",
                 num_recycle: int = 4,
                 chunk_size: Optional[int] = None,
                 **kwargs):
        """
        Initialize ESMFold configuration.

        Args:
            sequences: Input sequences - can be CSV file, list, ToolOutput, or StandardizedOutput
            datasheets: Input datasheet files for metadata
            name: Job name for output files
            num_recycle: Number of recycling iterations (default: 4)
            chunk_size: Chunk size for long sequences (default: None = auto)
            **kwargs: Additional parameters
        """
        # Handle different input formats
        if sequences is not None:
            if isinstance(sequences, StandardizedOutput):
                # StandardizedOutput object
                self.input_sequences = sequences.sequences
                self.input_datasheets = sequences.datasheets
                self.input_is_tool_output = False
                self.standardized_input = sequences
            elif isinstance(sequences, ToolOutput):
                # Direct ToolOutput object
                self.input_sequences = sequences
                self.input_datasheets = sequences.get_output_files("datasheets")
                self.input_is_tool_output = True
                self.standardized_input = None
            else:
                # Direct sequence(s) - string or list
                self.input_sequences = sequences
                self.input_datasheets = datasheets or {}
                self.input_is_tool_output = isinstance(sequences, ToolOutput)
                self.standardized_input = None
        else:
            raise ValueError("sequences parameter is required")

        # Store ESMFold-specific parameters
        self.name = name or kwargs.get('job_name', '')
        self.num_recycle = num_recycle
        self.chunk_size = chunk_size

        # Initialize base class
        super().__init__(**kwargs)

        # Initialize file paths (will be set in configure_inputs)
        self._initialize_file_paths()

    def validate_params(self):
        """Validate ESMFold-specific parameters."""
        if not self.input_sequences:
            raise ValueError("sequences parameter is required")

        if self.num_recycle < 0:
            raise ValueError("num_recycle cannot be negative")

        if self.chunk_size is not None and self.chunk_size < 1:
            raise ValueError("chunk_size must be at least 1")

    def _initialize_file_paths(self):
        """Initialize common file paths used throughout the class."""
        self.queries_csv = None
        self.structures_folder = None
        self.helper_script = None
        self.model_cache = None

    def _setup_file_paths(self):
        """Set up all file paths after output_folder is known."""
        # Extract pipeline/job name for file naming
        job_base = self.name or self._extract_job_name()

        # Core input/output files
        self.queries_csv = os.path.join(self.output_folder, f"{job_base}_queries.csv")
        self.structures_folder = os.path.join(self.output_folder, "structures")

        # Helper script path and model cache
        if hasattr(self, 'folders') and self.folders:
            self.helper_script = os.path.join(self.folders["HelpScripts"], "pipe_esmfold.py")
            # Use shared models folder for ESMFold cache
            self.model_cache = os.path.join(self.folders["group"], self.folders["shared"]["models"], "ESMFold")
        else:
            self.helper_script = None
            self.model_cache = None

    def _extract_job_name(self) -> str:
        """Extract job name from output folder structure."""
        # Get job name from parent folder
        # Structure: .../JobName_NNN/N_ESMFold -> JobName_NNN
        folder_parts = self.output_folder.split(os.sep)
        for i, part in enumerate(folder_parts):
            if "_ESMFold" in part:
                if i == 0:
                    raise ValueError(f"Invalid output folder structure: {self.output_folder}")
                return folder_parts[i-1]

        raise ValueError(f"Could not extract job name from output folder: {self.output_folder}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input sequences from various sources."""
        self.folders = pipeline_folders
        self._setup_file_paths()

        if self.input_is_tool_output:
            # Input from previous tool
            tool_output: ToolOutput = self.input_sequences

            # Try to get sequences - could be in various formats
            source_sequences = []

            # Try different output types
            for seq_type in ["sequences", "queries_csv", "fa_files"]:
                seq_files = tool_output.get_output_files(seq_type)
                if seq_files:
                    source_sequences = seq_files
                    break

            if not source_sequences:
                raise ValueError(f"No sequence outputs found from {tool_output.tool_type}")

            # Store source for script generation
            self.input_sources["sequences"] = source_sequences[0]

            # Add dependency
            self.dependencies.append(tool_output.config)

        elif isinstance(self.input_sequences, list):
            # Direct list of sequence file paths (from StandardizedOutput)
            if self.input_sequences:
                self.input_sources["sequences"] = self.input_sequences[0]
            else:
                raise ValueError("Empty sequence list provided")

        elif isinstance(self.input_sequences, str):
            # String input - could be file or direct sequence
            if self.input_sequences.endswith('.csv'):
                # CSV file with id,sequence columns
                csv_source = os.path.join(pipeline_folders.get("data", os.getcwd()), self.input_sequences)
                if not os.path.exists(csv_source):
                    raise ValueError(f"CSV file not found: {csv_source}")
                self.input_sources["sequences"] = csv_source
            else:
                # Direct sequence - will create CSV in generate_script
                if not self.name:
                    raise ValueError("name parameter required when providing direct sequence")
                self.input_sources["direct_sequence"] = self.input_sequences
        else:
            raise ValueError(f"Unsupported input type: {type(self.input_sequences)}")

    def get_config_display(self) -> List[str]:
        """Get ESMFold configuration display lines."""
        config_lines = super().get_config_display()

        # Input information
        if self.input_is_tool_output:
            config_lines.append(f"INPUT: {self.input_sequences.tool_type} output")
        else:
            config_lines.append(f"INPUT: {self.input_sequences}")

        config_lines.extend([
            f"NUM RECYCLE: {self.num_recycle}"
        ])

        if self.chunk_size:
            config_lines.append(f"CHUNK SIZE: {self.chunk_size}")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """
        Generate ESMFold execution script.

        Args:
            script_path: Path where script should be written

        Returns:
            Script content as string
        """
        script_content = "#!/bin/bash\n"
        script_content += "# ESMFold execution script\n"
        script_content += "# Generated by BioPipelines pipeline system\n\n"
        script_content += self.generate_completion_check_header()
        script_content += self._generate_prepare_sequences()
        script_content += self._generate_run_esmfold()
        script_content += self.generate_completion_check_footer()

        return script_content

    def _generate_prepare_sequences(self) -> str:
        """Generate the sequence preparation part of the script."""
        if "direct_sequence" in self.input_sources:
            # Create CSV from direct sequence
            sequence = self.input_sources["direct_sequence"]
            job_base = self.name or "sequence"

            return f"""echo "Creating queries CSV from direct sequence"
cat > {self.queries_csv} << EOF
id,sequence
{job_base},{sequence}
EOF

"""
        elif "sequences" in self.input_sources:
            source_file = self.input_sources["sequences"]

            return f"""echo "Using sequences from: {source_file}"
if [ -f "{source_file}" ]; then
    cp "{source_file}" "{self.queries_csv}"
    echo "Successfully copied sequences file"
else
    echo "ERROR: Sequence file not found: {source_file}"
    exit 1
fi

"""
        else:
            raise ValueError("No sequence input configured")

    def _generate_run_esmfold(self) -> str:
        """Generate the ESMFold execution part of the script."""
        # Build ESMFold options
        options = f"--num-recycles {self.num_recycle}"
        if self.chunk_size:
            options += f" --chunk-size {self.chunk_size}"

        return f"""echo "Running ESMFold"
echo "Model cache: {self.model_cache}"
echo "Options: {options}"
echo "Output folder: {self.structures_folder}"

# Create structures subfolder and model cache
mkdir -p {self.structures_folder}
mkdir -p {self.model_cache}

# Run ESMFold via helper script with custom cache location
python {self.helper_script} \\
    --input {self.queries_csv} \\
    --output {self.structures_folder} \\
    --model-cache {self.model_cache} \\
    {options}

if [ $? -eq 0 ]; then
    echo "ESMFold completed successfully"
else
    echo "ERROR: ESMFold failed"
    exit 1
fi

"""

    def _predict_structure_outputs(self) -> List[str]:
        """
        Predict structure files that ESMFold will generate.

        Based on input sequences.
        """
        structure_files = []
        sequence_ids = []

        # Try to get sequence IDs from standardized input
        if hasattr(self, 'standardized_input') and self.standardized_input:
            if hasattr(self.standardized_input, 'sequence_ids'):
                sequence_ids = self.standardized_input.sequence_ids
            elif 'sequence_ids' in self.standardized_input._data:
                sequence_ids = self.standardized_input._data['sequence_ids']

        # Try to get sequence info from direct ToolOutput
        if not sequence_ids and self.input_is_tool_output:
            tool_output: ToolOutput = self.input_sequences
            if hasattr(tool_output.config, '_predict_sequence_ids'):
                sequence_ids = tool_output.config._predict_sequence_ids()

        # If we have dependencies, try to get sequence IDs from them
        if not sequence_ids and self.dependencies:
            for dep in self.dependencies:
                if hasattr(dep, '_predict_sequence_ids'):
                    sequence_ids = dep._predict_sequence_ids()
                    break

        # Generate structure file paths from sequence IDs
        if sequence_ids:
            for seq_id in sequence_ids:
                pdb_path = os.path.join(self.structures_folder, f"{seq_id}.pdb")
                structure_files.append(pdb_path)

        if not structure_files:
            raise ValueError("Could not determine sequence IDs for structure prediction")

        return structure_files

    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after ESMFold execution.

        Returns:
            Dictionary mapping output types to file paths
        """
        # Ensure file paths are set up
        if not hasattr(self, 'queries_csv') or self.queries_csv is None:
            self._setup_file_paths()

        queries_csv = self.queries_csv

        # Predict structure outputs
        structure_files = self._predict_structure_outputs()

        # Extract structure IDs from predicted files
        structure_ids = []
        for struct_file in structure_files:
            struct_id = os.path.splitext(os.path.basename(struct_file))[0]
            structure_ids.append(struct_id)

        # Get sequence IDs from upstream input
        sequence_ids = []
        if hasattr(self, 'standardized_input') and self.standardized_input:
            if hasattr(self.standardized_input, 'sequence_ids'):
                sequence_ids = self.standardized_input.sequence_ids

        # Organize datasheets
        datasheets = {
            "structures": DatasheetInfo(
                name="structures",
                path=queries_csv,
                columns=["id", "sequence"],
                description="Input sequences used for ESMFold structure prediction",
                count=len(structure_ids)
            )
        }

        return {
            "structures": structure_files,
            "structure_ids": structure_ids,
            "compounds": [],
            "compound_ids": [],
            "sequences": [queries_csv],
            "sequence_ids": sequence_ids,
            "datasheets": datasheets,
            "output_folder": self.output_folder,
            "pdbs": structure_files,
            "queries_csv": [queries_csv]
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including ESMFold-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "esmfold_params": {
                "name": self.name,
                "num_recycle": self.num_recycle,
                "chunk_size": self.chunk_size,
                "input_type": "tool_output" if self.input_is_tool_output else "direct"
            }
        })
        return base_dict
