# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
MutationComposer tool for generating sequences based on mutation profiles.

Takes mutation frequency data from MutationProfiler and generates new sequences
using various composition strategies based on the observed mutation patterns.
"""

import os
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

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== MutationComposer ==="
echo "Requires MutationEnv (installed with MutationProfiler.install())"
echo "No additional installation needed."
echo "=== MutationComposer ready ==="
"""

    # Lazy path descriptors
    sequences_csv = Path(lambda self: os.path.join(self.output_folder, "sequences.csv"))
    sequences_fasta = Path(lambda self: os.path.join(self.output_folder, "sequences.fasta"))
    config_file = Path(lambda self: os.path.join(self.output_folder, "mutation_composer_config.json"))
    composer_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_mutation_composer.py"))

    def __init__(self,
                 frequencies: Union[List[Union[StandardizedOutput, TableInfo, str]],
                                  Union[StandardizedOutput, TableInfo, str]],
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
            frequencies: Input table(s) with mutation frequencies. Can be:
                - Single table: StandardizedOutput, TableInfo, or str
                - Multiple tables: List of the above types
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
            combination_strategy: Strategy for combining multiple frequency tables:
                - "average": Average frequencies across tables, then single mutations
                - "maximum": Take max frequencies across tables, then single mutations
                - "stack": One mutation from each table per sequence (multi-point)
                - "round_robin": Alternate tables for single mutations
            **kwargs: Additional parameters

        Output:
            Streams: sequences (.csv)
            Tables:
                sequences: id | sequence | mutations | mutation_positions
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

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate MutationComposer parameters."""
        # Validate mode
        valid_modes = ["single_point", "weighted_random", "hotspot_focused", "top_mutations"]
        if self.mode not in valid_modes:
            raise ValueError(f"Invalid mode: {self.mode}. Valid options: {', '.join(valid_modes)}")

        # Validate combination strategy
        valid_strategies = ["average", "maximum", "stack", "round_robin"]
        if self.combination_strategy not in valid_strategies:
            raise ValueError(f"Invalid combination_strategy: {self.combination_strategy}. Valid options: {', '.join(valid_strategies)}")

        # Validate each mutation input
        for i, mutation_input in enumerate(self.mutation_inputs):
            if not isinstance(mutation_input, (StandardizedOutput, TableInfo, str)):
                raise ValueError(f"frequencies[{i}] must be a StandardizedOutput, TableInfo, or str")

        if self.num_sequences <= 0:
            raise ValueError("num_sequences must be positive")

        if self.min_frequency < 0 or self.min_frequency > 1:
            raise ValueError("min_frequency must be between 0 and 1")

        if self.max_mutations is not None and self.max_mutations <= 0:
            raise ValueError("max_mutations must be positive or None")
    
    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input tables from MutationProfiler(s)."""
        self.folders = pipeline_folders

        # Extract mutation frequency table paths
        self.frequency_table_paths = []

        for i, mutation_input in enumerate(self.mutation_inputs):
            if isinstance(mutation_input, TableInfo):
                frequency_path = mutation_input.info.path
            elif isinstance(mutation_input, str):
                frequency_path = mutation_input
            elif isinstance(mutation_input, StandardizedOutput):
                # StandardizedOutput - look for absolute_frequencies or similar
                tables = mutation_input.tables
                if hasattr(tables, 'absolute_frequencies'):
                    frequency_path = tables.absolute_frequencies.info.path
                elif hasattr(tables, 'frequencies'):
                    frequency_path = tables.frequencies.info.path
                else:
                    raise ValueError(f"No frequency table found in frequencies[{i}]")
            else:
                raise ValueError(f"Unsupported input type for frequencies[{i}]: {type(mutation_input)}")

            self.frequency_table_paths.append(frequency_path)
    
    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()
        
        config_lines.extend([
            f"INPUTS: {len(self.mutation_inputs)} frequency table(s)",
            f"COMBINATION_STRATEGY: {self.combination_strategy}",
            f"NUM_SEQUENCES: {self.num_sequences}",
            f"MODE: {self.mode}",
            f"MIN_FREQUENCY: {self.min_frequency}"
        ])
        
        if self.max_mutations:
            config_lines.append(f"MAX_MUTATIONS: {self.max_mutations}")
        
        return config_lines
    
    def generate_script(self, script_path: str) -> str:
        """Generate MutationComposer execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        script_content = "#!/bin/bash\n"
        script_content += "# MutationComposer execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self.generate_script_run_composer()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_composer(self) -> str:
        """Generate the mutation composer execution part of the script."""
        import json

        # Create config file for mutation composer
        config_data = {
            "frequency_tables": self.frequency_table_paths,
            "combination_strategy": self.combination_strategy,
            "num_sequences": self.num_sequences,
            "mode": self.mode,
            "min_frequency": self.min_frequency,
            "max_mutations": self.max_mutations,
            "random_seed": self.random_seed,
            "prefix": self.prefix,
            "hotspot_count": self.hotspot_count,
            "sequences_output": self.sequences_csv,
            "sequences_fasta": self.sequences_fasta
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        return f"""echo "Generating sequences from mutation profiles"
echo "Input frequencies: {len(self.frequency_table_paths)} table(s)"
echo "Combination strategy: {self.combination_strategy}"
echo "Mode: {self.mode}"
echo "Sequences to generate: {self.num_sequences}"
echo "Output: {self.sequences_csv}"

python "{self.composer_py}" --config "{self.config_file}"

"""
    
    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after sequence generation."""
        # Generate sequence IDs based on prefix or mode
        if self.prefix:
            sequence_ids = [f"{self.prefix}_{i:03d}" for i in range(1, self.num_sequences + 1)]
        else:
            mode_prefixes = {
                "single_point": "mut",
                "weighted_random": "rand",
                "hotspot_focused": "hotspot",
                "top_mutations": "top"
            }
            prefix = mode_prefixes.get(self.mode, "seq")
            sequence_ids = [f"{prefix}_{i:03d}" for i in range(1, self.num_sequences + 1)]

        sequences = DataStream(
            name="sequences",
            ids=sequence_ids,
            files=[],
            map_table=self.sequences_csv,
            format="csv"
        )

        tables = {
            "sequences": TableInfo(
                name="sequences",
                path=self.sequences_csv,
                columns=["id", "sequence", "mutations", "mutation_positions"],
                description=f"Generated sequences using {self.mode} mode",
                count=self.num_sequences
            )
        }

        return {
            "sequences": sequences,
            "tables": tables,
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