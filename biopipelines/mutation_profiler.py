# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
MutationProfiler tool for analyzing mutation patterns between original and mutant sequences.

Takes original sequences and mutant sequences as input, analyzes mutation patterns
across all positions, and outputs comprehensive mutation statistics and frequency data.
Also generates sequence logo visualizations showing mutation patterns.
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

# Standard amino acids - guaranteed output structure
AMINO_ACIDS = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]


class MutationProfiler(BaseConfig):
    """
    Pipeline tool for analyzing mutation patterns between original and mutant sequences.
    
    Takes original and mutant sequences as input and generates comprehensive mutation
    analysis including position-wise statistics, amino acid frequencies, and sequence logos.
    
    Commonly used for:
    - LigandMPNN mutation pattern analysis
    - Protein design mutation profiling
    - Evolutionary analysis of sequence variants
    - Hotspot identification for further optimization
    """
    
    # Tool identification
    TOOL_NAME = "MutationProfiler"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        if env_manager == "pip":
            skip = "" if force_reinstall else """# Check if already installed
if python -c "import seaborn; import logomaker" 2>/dev/null; then
    echo "MutationEnv deps already installed, skipping. Use force_reinstall=True to reinstall."
    exit 0
fi
"""
            return f"""echo "=== Installing MutationEnv deps (pip) ==="
{skip}pip install seaborn matplotlib pandas logomaker scipy

echo "=== MutationEnv installation complete ==="
"""
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_manager} env list 2>/dev/null | grep -q "MutationEnv"; then
    echo "MutationEnv already installed, skipping. Use force_reinstall=True to reinstall."
    exit 0
fi
"""
        return f"""echo "=== Installing MutationEnv ==="
{skip}{env_manager} create -n MutationEnv seaborn matplotlib pandas logomaker scipy -y

echo "=== MutationEnv installation complete ==="
"""

    # Lazy path descriptors
    profile_csv = Path(lambda self: os.path.join(self.output_folder, "profile.csv"))
    mutations_csv = Path(lambda self: os.path.join(self.output_folder, "mutations.csv"))
    absolute_freq_csv = Path(lambda self: os.path.join(self.output_folder, "absolute_frequencies.csv"))
    relative_freq_csv = Path(lambda self: os.path.join(self.output_folder, "relative_frequencies.csv"))
    sequence_logo_relative_svg = Path(lambda self: os.path.join(self.output_folder, "sequence_logo_relative_frequencies.svg"))
    sequence_logo_relative_png = Path(lambda self: os.path.join(self.output_folder, "sequence_logo_relative_frequencies.png"))
    sequence_logo_absolute_svg = Path(lambda self: os.path.join(self.output_folder, "sequence_logo_absolute_frequencies.svg"))
    sequence_logo_absolute_png = Path(lambda self: os.path.join(self.output_folder, "sequence_logo_absolute_frequencies.png"))
    sequence_logo_counts_svg = Path(lambda self: os.path.join(self.output_folder, "sequence_logo_counts.svg"))
    sequence_logo_counts_png = Path(lambda self: os.path.join(self.output_folder, "sequence_logo_counts.png"))
    config_file = Path(lambda self: os.path.join(self.output_folder, "mutation_profiler_config.json"))
    profiler_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_mutation_profiler.py"))

    def __init__(self,
                 original: Union[DataStream, StandardizedOutput],
                 mutants: Union[DataStream, StandardizedOutput],
                 include_original: bool = True,
                 positions: Optional[str] = None,
                 **kwargs):
        """
        Initialize mutation profiler tool.

        Args:
            original: Original sequences as DataStream or StandardizedOutput
            mutants: Mutant sequences as DataStream or StandardizedOutput
            include_original: Whether to include original sequence in analysis (default: True)
            positions: PyMOL-style selection string for positions to display in plots (e.g., "141+143+145+147-149")
                      If None, shows all positions with mutations.
            **kwargs: Additional parameters

        Output:
            Streams: (none)
            Tables:
                profile: position | original | count | frequency
                mutations: position | original | A | C | D | ... | Y
                absolute_frequencies: position | original | A | C | D | ... | Y
                relative_frequencies: position | original | A | C | D | ... | Y
        """
        # Resolve original to DataStream
        if isinstance(original, StandardizedOutput):
            self.original_stream: DataStream = original.streams.sequences
        elif isinstance(original, DataStream):
            self.original_stream = original
        else:
            raise ValueError(f"original must be DataStream or StandardizedOutput, got {type(original)}")

        # Resolve mutants to DataStream
        if isinstance(mutants, StandardizedOutput):
            self.mutants_stream: DataStream = mutants.streams.sequences
        elif isinstance(mutants, DataStream):
            self.mutants_stream = mutants
        else:
            raise ValueError(f"mutants must be DataStream or StandardizedOutput, got {type(mutants)}")

        self.include_original = include_original
        self.positions = positions

        super().__init__(**kwargs)


    def validate_params(self):
        """Validate MutationProfiler parameters."""
        if not self.original_stream or len(self.original_stream) == 0:
            raise ValueError("original sequences cannot be empty")

        if not self.mutants_stream or len(self.mutants_stream) == 0:
            raise ValueError("mutants sequences cannot be empty")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input sequences."""
        self.folders = pipeline_folders
    
    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"ORIGINAL: {len(self.original_stream)} sequences",
            f"MUTANTS: {len(self.mutants_stream)} sequences",
            f"INCLUDE_ORIGINAL: {self.include_original}",
            f"POSITIONS: {self.positions if self.positions else 'auto (all mutations)'}"
        ])

        return config_lines
    
    def generate_script(self, script_path: str) -> str:
        """Generate MutationProfiler execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        script_content = "#!/bin/bash\n"
        script_content += "# MutationProfiler execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self.generate_script_run_profiler()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_profiler(self) -> str:
        """Generate the mutation profiler execution part of the script."""
        import json

        config_data = {
            "original_sequences": self.original_stream.map_table,
            "mutants_sequences": self.mutants_stream.map_table,
            "include_original": self.include_original,
            "positions": self.positions,
            "profile_output": self.profile_csv,
            "mutations_output": self.mutations_csv,
            "absolute_frequencies_output": self.absolute_freq_csv,
            "relative_frequencies_output": self.relative_freq_csv,
            "sequence_logo_relative_svg": self.sequence_logo_relative_svg,
            "sequence_logo_relative_png": self.sequence_logo_relative_png,
            "sequence_logo_absolute_svg": self.sequence_logo_absolute_svg,
            "sequence_logo_absolute_png": self.sequence_logo_absolute_png,
            "sequence_logo_counts_svg": self.sequence_logo_counts_svg,
            "sequence_logo_counts_png": self.sequence_logo_counts_png
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        return f"""echo "Analyzing mutation patterns"
echo "Original sequences: {self.original_stream.map_table}"
echo "Mutant sequences: {self.mutants_stream.map_table}"
echo "Output folder: {self.output_folder}"

python "{self.profiler_py}" --config "{self.config_file}"

"""
    
    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after mutation profiling."""
        aa_columns = ["position", "original"] + AMINO_ACIDS

        tables = {
            "profile": TableInfo(
                name="profile",
                path=self.profile_csv,
                columns=["position", "original", "count", "frequency"],
                description="Position-wise mutation statistics",
                count=0
            ),
            "mutations": TableInfo(
                name="mutations",
                path=self.mutations_csv,
                columns=aa_columns,
                description="Raw amino acid mutation counts per position",
                count=0
            ),
            "absolute_frequencies": TableInfo(
                name="absolute_frequencies",
                path=self.absolute_freq_csv,
                columns=aa_columns,
                description="Absolute amino acid frequencies per position (mutation_count/total_sequences)",
                count=0
            ),
            "relative_frequencies": TableInfo(
                name="relative_frequencies",
                path=self.relative_freq_csv,
                columns=aa_columns,
                description="Relative amino acid frequencies per position (mutation_count/mutations_at_position)",
                count=0
            )
        }

        return {
            "tables": tables,
            "visualizations": [
                self.sequence_logo_relative_svg, self.sequence_logo_relative_png,
                self.sequence_logo_absolute_svg, self.sequence_logo_absolute_png,
                self.sequence_logo_counts_svg, self.sequence_logo_counts_png
            ],
            "output_folder": self.output_folder
        }
    
    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                "include_original": self.include_original,
                "positions": self.positions
            }
        })
        return base_dict