# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
ESMFold configuration for fast protein structure prediction.

Uses Meta's ESM-2 with ESMFold for efficient single-sequence structure prediction
without requiring MSAs. Models are loaded via torch.hub and cached in shared folder.
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
    ESMFold configuration for fast protein structure prediction.

    Predicts protein structures from sequences using ESM-2 + ESMFold.
    No MSAs required - fast single-sequence prediction.
    """

    TOOL_NAME = "ESMFold"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba"):
        return f"""echo "=== Installing ESMFold ==="
{env_manager} create -n esmfold python=3.9 -y
{env_manager} activate esmfold

{env_manager} install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia -y

pip install "fair-esm[esmfold]"
pip install 'dllogger @ git+https://github.com/NVIDIA/dllogger.git'
pip install 'openfold @ git+https://github.com/aqlaboratory/openfold.git@4b41059694619831a7db195b7e0988fc4ff3a307'

echo "=== ESMFold installation complete ==="
"""

    # Lazy path descriptors
    queries_csv = Path(lambda self: os.path.join(self.output_folder, f"{self.pipeline_name}_queries.csv"))
    structures_folder = Path(lambda self: os.path.join(self.output_folder, "structures"))
    structures_map = Path(lambda self: os.path.join(self.output_folder, "structures_map.csv"))

    # Helper script and model cache paths
    helper_script = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_esmfold.py"))
    model_cache = Path(lambda self: self.folders["ESMFold"])

    def __init__(self,
                 sequences: Union[DataStream, StandardizedOutput],
                 num_recycle: int = 4,
                 chunk_size: Optional[int] = None,
                 **kwargs):
        """
        Initialize ESMFold configuration.

        Args:
            sequences: Input sequences as DataStream or StandardizedOutput
            num_recycle: Number of recycling iterations (default: 4)
            chunk_size: Chunk size for long sequences (default: None = auto)
        """
        # Resolve input to DataStream
        if isinstance(sequences, StandardizedOutput):
            self.sequences_stream: DataStream = sequences.sequences
        elif isinstance(sequences, DataStream):
            self.sequences_stream = sequences
        else:
            raise ValueError(f"sequences must be DataStream or StandardizedOutput, got {type(sequences)}")

        # Store ESMFold-specific parameters
        self.num_recycle = num_recycle
        self.chunk_size = chunk_size

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate ESMFold-specific parameters."""
        if not self.sequences_stream or len(self.sequences_stream) == 0:
            raise ValueError("sequences parameter is required and must not be empty")

        if self.num_recycle < 0:
            raise ValueError("num_recycle cannot be negative")

        if self.chunk_size is not None and self.chunk_size < 1:
            raise ValueError("chunk_size must be at least 1")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input sequences."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get ESMFold configuration display lines."""
        config_lines = super().get_config_display()
        config_lines.extend([
            f"INPUT SEQUENCES: {len(self.sequences_stream)} sequences",
            f"NUM RECYCLE: {self.num_recycle}"
        ])

        if self.chunk_size:
            config_lines.append(f"CHUNK SIZE: {self.chunk_size}")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate ESMFold execution script."""
        script_content = "#!/bin/bash\n"
        script_content += "# ESMFold execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self._generate_script_prepare_sequences()
        script_content += self._generate_script_run_esmfold()
        script_content += self.generate_completion_check_footer()

        return script_content

    def _generate_script_prepare_sequences(self) -> str:
        """Generate the sequence preparation part of the script."""
        source_file = self.sequences_stream.map_table

        return f"""echo "Using sequences from: {source_file}"
if [ -f "{source_file}" ]; then
    cp "{source_file}" "{self.queries_csv}"
    echo "Successfully copied sequences file"
else
    echo "ERROR: Sequence file not found: {source_file}"
    echo "This usually means the previous step failed to generate the expected output"
    exit 1
fi

"""

    def _generate_script_run_esmfold(self) -> str:
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

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after ESMFold execution."""
        # Use sequence IDs from input to predict structure files
        sequence_ids = self.sequences_stream.ids

        # Generate structure file paths
        structure_files = []
        for seq_id in sequence_ids:
            pdb_path = os.path.join(self.structures_folder, f"{seq_id}.pdb")
            structure_files.append(pdb_path)

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
            )
        }

        return {
            "structures": structures,
            "sequences": DataStream.empty("sequences", "fasta"),
            "compounds": DataStream.empty("compounds", "sdf"),
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including ESMFold-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "esmfold_params": {
                "num_recycle": self.num_recycle,
                "chunk_size": self.chunk_size
            }
        })
        return base_dict
