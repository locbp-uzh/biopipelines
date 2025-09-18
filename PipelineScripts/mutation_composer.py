"""
MutationComposer tool for generating sequences based on mutation profiles.

Takes mutation frequency data from MutationProfiler and generates new sequences
using various composition strategies based on the observed mutation patterns.
"""

import os
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo


class MutationComposer(BaseConfig):
    """
    Pipeline tool for generating sequences based on mutation frequency profiles.
    
    Takes mutation frequency data (typically from MutationProfiler) and generates
    new sequences using various composition strategies informed by mutation patterns.
    
    Commonly used for:
    - Systematic single-point mutation design
    - Weighted random sequence generation
    - Hotspot-focused sequence design
    - Exploration of high-frequency mutation combinations
    """
    
    # Tool identification
    TOOL_NAME = "MutationComposer"
    DEFAULT_ENV = "MutationEnv" #with seaborn, pandas, matplotlib, logomaker
    DEFAULT_RESOURCES = {"gpu": "T4", "memory": "4GB", "time": "1:00:00"}
    
    def __init__(self,
                 frequencies: Union[List[Union[ToolOutput, StandardizedOutput, DatasheetInfo, str]],
                                  Union[ToolOutput, StandardizedOutput, DatasheetInfo, str]],
                 num_sequences: int = 10,
                 mode: str = "single_point",
                 min_frequency: float = 0.01,
                 max_mutations: int = None,
                 random_seed: int = None,
                 prefix: str = "",
                 hotspot_count: int = 10,
                 combination_strategy: str = "average",
                 **kwargs):
        """
        Initialize mutation composer tool.

        Args:
            frequencies: Input datasheet(s) with mutation frequencies. Can be:
                - Single datasheet: ToolOutput, StandardizedOutput, DatasheetInfo, or str
                - Multiple datasheets: List of the above types
            num_sequences: Number of sequences to generate
            mode: Generation strategy:
                - "single_point": One mutation per sequence, ordered by frequency
                - "weighted_random": Random sequences using frequency distributions
                - "hotspot_focused": Focus mutations on high-frequency positions
                - "top_mutations": Use top N most frequent mutations per position
            min_frequency: Minimum frequency threshold to consider mutations (default: 0.01)
            max_mutations: Maximum mutations per sequence (None = unlimited)
            random_seed: Random seed for reproducible generation (None = random)
            prefix: Prefix for sequence IDs (e.g., "single_point" -> "single_point_001")
            hotspot_count: Number of top hotspot positions to consider for hotspot_focused mode (default: 10)
            combination_strategy: Strategy for combining multiple frequency datasheets:
                - "average": Average frequencies across datasheets, then single mutations
                - "maximum": Take max frequencies across datasheets, then single mutations
                - "stack": One mutation from each datasheet per sequence (multi-point)
                - "round_robin": Alternate datasheets for single mutations
            **kwargs: Additional parameters

        Examples:
            # Generate single-point mutations from single datasheet
            composer = pipeline.add(MutationComposer(
                frequencies=profiler.datasheets.absolute_frequencies,
                num_sequences=20,
                mode="single_point"
            ))

            # Generate sequences with one mutation from each enantiomer
            composer = pipeline.add(MutationComposer(
                frequencies=[profiler_R.datasheets.absolute_frequencies,
                           profiler_S.datasheets.absolute_frequencies],
                num_sequences=50,
                mode="single_point",
                combination_strategy="stack"
            ))
        """
        # Normalize frequencies to list for consistent handling
        if isinstance(frequencies, list):
            self.mutation_inputs = frequencies
        else:
            self.mutation_inputs = [frequencies]

        self.combination_strategy = combination_strategy
        self.num_sequences = num_sequences
        self.mode = mode
        self.min_frequency = min_frequency
        self.max_mutations = max_mutations
        self.random_seed = random_seed
        self.prefix = prefix
        self.hotspot_count = hotspot_count
        
        # Validate mode
        valid_modes = ["single_point", "weighted_random", "hotspot_focused", "top_mutations"]
        if mode not in valid_modes:
            raise ValueError(f"Invalid mode: {mode}. Valid options: {', '.join(valid_modes)}")

        # Validate combination strategy
        valid_strategies = ["average", "maximum", "stack", "round_robin"]
        if combination_strategy not in valid_strategies:
            raise ValueError(f"Invalid combination_strategy: {combination_strategy}. Valid options: {', '.join(valid_strategies)}")
        
        # Initialize base class
        super().__init__(**kwargs)
        
        # Set up dependencies if inputs are from tools
        for mutation_input in self.mutation_inputs:
            if hasattr(mutation_input, 'config'):
                self.dependencies.append(mutation_input.config)
    
    def validate_params(self):
        """Validate MutationComposer parameters."""
        # Validate each mutation input
        for i, mutation_input in enumerate(self.mutation_inputs):
            if not isinstance(mutation_input, (ToolOutput, StandardizedOutput, DatasheetInfo, str)):
                raise ValueError(f"frequencies[{i}] must be a ToolOutput, StandardizedOutput, DatasheetInfo, or str object")
        
        if self.num_sequences <= 0:
            raise ValueError("num_sequences must be positive")
        
        if self.min_frequency < 0 or self.min_frequency > 1:
            raise ValueError("min_frequency must be between 0 and 1")
        
        if self.max_mutations is not None and self.max_mutations <= 0:
            raise ValueError("max_mutations must be positive or None")
    
    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input datasheets from MutationProfiler(s)."""
        self.folders = pipeline_folders

        # Extract mutation frequency datasheet paths
        self.frequency_datasheet_paths = []

        for i, mutation_input in enumerate(self.mutation_inputs):
            if hasattr(mutation_input, 'path'):
                # DatasheetInfo object
                frequency_path = mutation_input.path
            elif isinstance(mutation_input, str):
                # Direct file path
                frequency_path = mutation_input
            elif hasattr(mutation_input, 'datasheets'):
                # ToolOutput - look for absolute_frequencies or similar
                datasheets = mutation_input.datasheets

                if hasattr(datasheets, 'absolute_frequencies'):
                    frequency_path = datasheets.absolute_frequencies.path
                elif hasattr(datasheets, 'frequencies'):
                    frequency_path = datasheets.frequencies.path
                elif hasattr(datasheets, '_datasheets'):
                    # Look for frequency-related datasheet
                    frequency_path = None
                    for name, info in datasheets._datasheets.items():
                        if 'frequenc' in name.lower() or 'mutation' in name.lower():
                            frequency_path = info.path
                            break
                    if frequency_path is None:
                        raise ValueError(f"No frequency datasheet found in frequencies[{i}]")
                else:
                    raise ValueError(f"Could not find frequency datasheet in frequencies[{i}]")
            else:
                raise ValueError(f"Unsupported input type for frequencies[{i}]: {type(mutation_input)}")

            self.frequency_datasheet_paths.append(frequency_path)
    
    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()
        
        config_lines.extend([
            f"INPUTS: {len(self.mutation_inputs)} frequency datasheet(s)",
            f"COMBINATION_STRATEGY: {self.combination_strategy}",
            f"NUM_SEQUENCES: {self.num_sequences}",
            f"MODE: {self.mode}",
            f"MIN_FREQUENCY: {self.min_frequency}"
        ])
        
        if self.max_mutations:
            config_lines.append(f"MAX_MUTATIONS: {self.max_mutations}")
        
        return config_lines
    
    def generate_script(self, script_path: str) -> str:
        """
        Generate script to compose sequences based on mutation profiles.
        
        Args:
            script_path: Path where script should be written
            
        Returns:
            Script content as string
        """
        runtime_folder = os.path.dirname(script_path)
        output_folder = self.output_folder
        os.makedirs(output_folder, exist_ok=True)
        
        # Output files
        sequences_csv = os.path.join(output_folder, "sequences.csv")
        sequences_fasta = os.path.join(output_folder, "sequences.fasta")
        
        # Create config file for mutation composer
        config_file = os.path.join(output_folder, "mutation_composer_config.json")
        config_data = {
            "frequency_datasheets": self.frequency_datasheet_paths,
            "combination_strategy": self.combination_strategy,
            "num_sequences": self.num_sequences,
            "mode": self.mode,
            "min_frequency": self.min_frequency,
            "max_mutations": self.max_mutations,
            "random_seed": self.random_seed,
            "prefix": self.prefix,
            "hotspot_count": self.hotspot_count,
            "sequences_output": sequences_csv,
            "sequences_fasta": sequences_fasta
        }
        
        import json
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)
        
        # Generate script content
        script_content = f"""#!/bin/bash
# MutationComposer execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Generating sequences from mutation profiles"
echo "Input frequencies: {len(self.frequency_datasheet_paths)} datasheet(s)"
echo "Combination strategy: {self.combination_strategy}"
echo "Mode: {self.mode}"
echo "Sequences to generate: {self.num_sequences}"
echo "Output: {sequences_csv}"

# Run Python mutation composer script
python "{os.path.join(self.folders['HelpScripts'], 'pipe_mutation_composer.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Successfully generated {self.num_sequences} sequences"
    echo "Sequences CSV: {sequences_csv}"
    echo "Sequences FASTA: {sequences_fasta}"
else
    echo "Error: Failed to generate sequences"
    exit 1
fi

{self.generate_completion_check_footer()}
"""
        
        return script_content
    
    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after sequence generation.
        
        Returns:
            Dictionary with output file paths and datasheet information
        """
        sequences_csv = os.path.join(self.output_folder, "sequences.csv")
        sequences_fasta = os.path.join(self.output_folder, "sequences.fasta")
        
        # Generate sequence IDs based on prefix or mode
        if self.prefix:
            # Use provided prefix
            sequence_ids = [f"{self.prefix}_{i:03d}" for i in range(1, self.num_sequences + 1)]
        else:
            # Use default prefixes based on mode
            if self.mode == "single_point":
                sequence_ids = [f"mut_{i:03d}" for i in range(1, self.num_sequences + 1)]
            elif self.mode == "weighted_random":
                sequence_ids = [f"rand_{i:03d}" for i in range(1, self.num_sequences + 1)]
            elif self.mode == "hotspot_focused":
                sequence_ids = [f"hotspot_{i:03d}" for i in range(1, self.num_sequences + 1)]
            else:  # top_mutations
                sequence_ids = [f"top_{i:03d}" for i in range(1, self.num_sequences + 1)]
        
        # Define datasheets that will be created
        datasheets = {
            "sequences": DatasheetInfo(
                name="sequences",
                path=sequences_csv,
                columns=["id", "sequence", "mutations", "mutation_positions"],
                description=f"Generated sequences using {self.mode} mode",
                count=self.num_sequences
            )
        }
        
        return {
            "structures": [],
            "structure_ids": [],
            "compounds": [],
            "compound_ids": [],
            "sequences": [sequences_csv, sequences_fasta],
            "sequence_ids": sequence_ids,
            "datasheets": datasheets,
            "output_folder": self.output_folder
        }
    
    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                "num_sequences": self.num_sequences,
                "mode": self.mode,
                "combination_strategy": self.combination_strategy,
                "min_frequency": self.min_frequency,
                "max_mutations": self.max_mutations,
                "random_seed": self.random_seed,
                "prefix": self.prefix
            }
        })
        return base_dict