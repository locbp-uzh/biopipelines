"""
SplitChains - Split single-chain sequences into multi-chain sequences.

This tool takes sequences where multiple protein chains are concatenated into
a single sequence and splits them into separate chains based on specified
split positions.
"""

import os
import json
from typing import Dict, List, Optional, Union

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo


class SplitChains(BaseConfig):
    """
    Split concatenated single-chain sequences into multi-chain sequences.

    Takes a sequences CSV where each row has one long sequence representing
    multiple chains, and splits it into multiple rows with separate chain sequences.
    """

    TOOL_NAME = "SplitChains"
    

    def __init__(self,
                 sequences: Union[ToolOutput, StandardizedOutput],
                 split_positions: List[int],
                 chain_names: Optional[List[str]] = None,
                 **kwargs):
        """
        Initialize SplitChains configuration.

        Args:
            sequences: Input sequences (ToolOutput or StandardizedOutput)
            split_positions: List of positions where to split sequences.
                           For example, [197] splits at position 197 creating 2 chains:
                           seq[0:197] and seq[197:]
            chain_names: Optional names for the chains (e.g., ["ChainA", "ChainB"])
                        If not provided, uses numeric suffixes (_1, _2, etc.)
            **kwargs: Additional parameters
        """
        # Store parameters before calling super().__init__() which calls validate_params()
        self.sequences_input = sequences
        self.split_positions = split_positions
        self.chain_names = chain_names

        # Now call parent constructor (which will call validate_params)
        super().__init__(**kwargs)

        # Handle StandardizedOutput or ToolOutput
        if isinstance(sequences, StandardizedOutput):
            if hasattr(sequences, 'sequences') and sequences.sequences:
                self.input_sequences_file = sequences.sequences[0]
                if hasattr(sequences, 'tool_config'):
                    self.dependencies.append(sequences.tool_config)
            else:
                raise ValueError("No sequences found in StandardizedOutput")
        elif isinstance(sequences, ToolOutput):
            seq_files = sequences.get_output_files("sequences")
            if seq_files:
                self.input_sequences_file = seq_files[0]
                self.dependencies.append(sequences.config)
            else:
                raise ValueError("No sequences found in ToolOutput")
        else:
            raise TypeError("sequences must be StandardizedOutput or ToolOutput")

    def validate_params(self):
        """Validate SplitChains parameters."""
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
        """Configure inputs - called by pipeline during setup."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()
        config_lines.extend([
            f"SPLIT POSITIONS: {self.split_positions}",
            f"NUMBER OF CHAINS: {len(self.split_positions) + 1}",
            f"INPUT SEQUENCES: {os.path.basename(self.input_sequences_file)}"
        ])
        if self.chain_names:
            config_lines.append(f"CHAIN NAMES: {', '.join(self.chain_names)}")
        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate the execution script."""
        output_folder = self.output_folder
        os.makedirs(output_folder, exist_ok=True)

        # Output file
        output_sequences_file = os.path.join(output_folder, "split_sequences.csv")

        # Create config file for the split operation
        config_file = os.path.join(output_folder, "split_config.json")
        config_data = {
            "input_csv": self.input_sequences_file,
            "output_csv": output_sequences_file,
            "split_positions": self.split_positions,
            "chain_names": self.chain_names
        }

        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        # Build script
        script_content = f"""#!/bin/bash
# SplitChains execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Splitting sequences into multiple chains"
echo "Input: {self.input_sequences_file}"
echo "Split positions: {self.split_positions}"
echo "Output: {output_sequences_file}"

# Run Python splitting script
python "{os.path.join(self.folders['HelpScripts'], 'pipe_split_chains.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Results written to: {output_sequences_file}"
else
    echo "Error: Splitting failed"
    exit 1
fi

{self.generate_completion_check_footer()}
"""

        return script_content

    def get_output_files(self) -> Dict[str, List[str]]:
        """Get expected output files."""
        output_sequences_file = os.path.join(self.output_folder, "split_sequences.csv")

        # Build columns list
        columns = ['id', 'sequence', 'source_id', 'complex_id', 'chain_index', 'chain_name', 'chain_length']

        tables = {
            "sequences": TableInfo(
                name="sequences",
                path=output_sequences_file,
                columns=columns,
                description="Split chain sequences with source tracking",
                count="variable"
            )
        }

        return {
            "structures": [],
            "structure_ids": [],
            "compounds": [],
            "compound_ids": [],
            "sequences": [output_sequences_file],
            "sequence_ids": [],  # Will be populated at runtime
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                "split_positions": self.split_positions,
                "chain_names": self.chain_names,
                "input_sequences_file": self.input_sequences_file
            }
        })
        return base_dict
