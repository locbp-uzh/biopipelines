# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
SplitChains - Split single-chain sequences into multi-chain sequences.

This tool takes sequences where multiple protein chains are concatenated into
a single sequence and splits them into separate chains based on specified
split positions.
"""

import os
import json
from typing import Dict, List, Optional, Union, Any

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


class SplitChains(BaseConfig):
    """
    Split concatenated single-chain sequences into multi-chain sequences.

    Takes a sequences CSV where each row has one long sequence representing
    multiple chains, and splits it into multiple rows with separate chain sequences.
    """

    TOOL_NAME = "SplitChains"

    # Lazy path descriptors
    output_sequences_csv = Path(lambda self: os.path.join(self.output_folder, "split_sequences.csv"))
    config_file = Path(lambda self: os.path.join(self.output_folder, "split_config.json"))
    helper_script = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_split_chains.py"))

    def __init__(self,
                 sequences: Union[DataStream, StandardizedOutput],
                 split_positions: List[int],
                 chain_names: Optional[List[str]] = None,
                 **kwargs):
        """
        Initialize SplitChains configuration.

        Args:
            sequences: Input sequences as DataStream or StandardizedOutput
            split_positions: List of positions where to split sequences.
                           For example, [197] splits at position 197 creating 2 chains:
                           seq[0:197] and seq[197:]
            chain_names: Optional names for the chains (e.g., ["ChainA", "ChainB"])
                        If not provided, uses numeric suffixes (_1, _2, etc.)
            **kwargs: Additional parameters
        """
        # Resolve input to DataStream
        if isinstance(sequences, StandardizedOutput):
            self.sequences_stream: DataStream = sequences.sequences
        elif isinstance(sequences, DataStream):
            self.sequences_stream = sequences
        else:
            raise ValueError(f"sequences must be DataStream or StandardizedOutput, got {type(sequences)}")

        self.split_positions = split_positions
        self.chain_names = chain_names

        # Initialize base class
        super().__init__(**kwargs)

    def validate_params(self):
        """Validate SplitChains parameters."""
        if not self.sequences_stream or len(self.sequences_stream) == 0:
            raise ValueError("sequences parameter is required and must not be empty")

        # Validate split positions
        if not self.split_positions:
            raise ValueError("split_positions must be a non-empty list")

        if not all(isinstance(pos, int) and pos > 0 for pos in self.split_positions):
            raise ValueError("All split_positions must be positive integers")

        # Check that positions are in ascending order
        if self.split_positions != sorted(self.split_positions):
            raise ValueError("split_positions must be in ascending order")

        # Validate chain names if provided
        num_chains = len(self.split_positions) + 1
        if self.chain_names and len(self.chain_names) != num_chains:
            raise ValueError(f"Number of chain names ({len(self.chain_names)}) must match number of chains ({num_chains})")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure inputs."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        # Get input file for display
        input_file = self.sequences_stream.map_table or (self.sequences_stream.files[0] if self.sequences_stream.files else "N/A")

        config_lines.extend([
            f"INPUT SEQUENCES: {len(self.sequences_stream)} sequences",
            f"SPLIT POSITIONS: {self.split_positions}",
            f"NUMBER OF CHAINS: {len(self.split_positions) + 1}"
        ])
        if self.chain_names:
            config_lines.append(f"CHAIN NAMES: {', '.join(self.chain_names)}")
        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate the execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        # Get input sequences file
        input_file = self.sequences_stream.map_table or (self.sequences_stream.files[0] if self.sequences_stream.files else None)
        if not input_file:
            raise ValueError("No input sequences file available")

        # Create config file for the split operation
        config_data = {
            "input_csv": input_file,
            "output_csv": self.output_sequences_csv,
            "split_positions": self.split_positions,
            "chain_names": self.chain_names
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        script_content = "#!/bin/bash\n"
        script_content += "# SplitChains execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()

        script_content += f"""echo "Splitting sequences into multiple chains"
echo "Input: {input_file}"
echo "Split positions: {self.split_positions}"
echo "Output: {self.output_sequences_csv}"

python "{self.helper_script}" \\
  --config "{self.config_file}"

if [ $? -eq 0 ]; then
    echo "Results written to: {self.output_sequences_csv}"
else
    echo "Error: Splitting failed"
    exit 1
fi

"""
        script_content += self.generate_completion_check_footer()

        return script_content

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files."""
        # Build columns list
        columns = ['id', 'sequence', 'source_id', 'complex_id', 'chain_index', 'chain_name', 'chain_length']

        tables = {
            "sequences": TableInfo(
                name="sequences",
                path=self.output_sequences_csv,
                columns=columns,
                description="Split chain sequences with source tracking",
                count="variable"
            )
        }

        # Create sequences DataStream (IDs will be populated at runtime)
        sequences = DataStream(
            name="sequences",
            ids=[],  # Will be populated at runtime
            files=[self.output_sequences_csv],
            map_table=self.output_sequences_csv,
            format="csv"
        )

        return {
            "structures": DataStream.empty("structures", "pdb"),
            "sequences": sequences,
            "compounds": DataStream.empty("compounds", "sdf"),
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()

        # Get input file for serialization
        input_file = self.sequences_stream.map_table or (self.sequences_stream.files[0] if self.sequences_stream.files else None)

        base_dict.update({
            "tool_params": {
                "split_positions": self.split_positions,
                "chain_names": self.chain_names,
                "input_sequences_file": input_file
            }
        })
        return base_dict
