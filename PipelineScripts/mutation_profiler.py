"""
MutationProfiler tool for analyzing mutation patterns between original and mutant sequences.

Takes original sequences and mutant sequences as input, analyzes mutation patterns
across all positions, and outputs comprehensive mutation statistics and frequency data.
Also generates sequence logo visualizations showing mutation patterns.
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
    DEFAULT_ENV = "MutationEnv"  #with seaborn pandas matplotlib logomaker
    COMPATIBLE_ENVS = ["MutationEnv"]
    DEFAULT_RESOURCES = {"gpu": "T4", "memory": "4GB", "time": "2:00:00"}
    
    def __init__(self,
                 original: Union[ToolOutput, StandardizedOutput],
                 mutants: Union[ToolOutput, StandardizedOutput],
                 include_original: bool = True,
                 **kwargs):
        """
        Initialize mutation profiler tool.
        
        Args:
            original: Original sequences (ToolOutput or StandardizedOutput)
            mutants: Mutant sequences to compare against original (ToolOutput or StandardizedOutput)  
            include_original: Whether to include original sequence in analysis (default: True)
            **kwargs: Additional parameters
            
        Examples:
            # Analyze LigandMPNN mutations against original
            profiler = pipeline.add(MutationProfiler(
                original=original_structure.output,
                mutants=lmpnn.output
            ))
            
            # Profile mutations without including original in statistics
            profiler = pipeline.add(MutationProfiler(
                original=wild_type.output,
                mutants=engineered_variants.output,
                include_original=False
            ))
        """
        self.original_input = original
        self.mutants_input = mutants
        self.include_original = include_original
        
        # Initialize base class
        super().__init__(**kwargs)
        
        # Set up dependencies
        if hasattr(original, 'config'):
            self.dependencies.append(original.config)
        if hasattr(mutants, 'config'):
            self.dependencies.append(mutants.config)
    
    def validate_params(self):
        """Validate MutationProfiler parameters."""
        if not isinstance(self.original_input, (ToolOutput, StandardizedOutput)):
            raise ValueError("original must be a ToolOutput or StandardizedOutput object")
        
        if not isinstance(self.mutants_input, (ToolOutput, StandardizedOutput)):
            raise ValueError("mutants must be a ToolOutput or StandardizedOutput object")
    
    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input sequences from previous tools."""
        self.folders = pipeline_folders
        
        # Extract original sequences
        self.original_sequences_path = self._extract_sequences_path(self.original_input, "original")
        
        # Extract mutant sequences  
        self.mutants_sequences_path = self._extract_sequences_path(self.mutants_input, "mutants")
    
    def _extract_sequences_path(self, input_obj: Union[ToolOutput, StandardizedOutput], input_type: str) -> str:
        """Extract sequences path from ToolOutput or StandardizedOutput."""
        # Try to get sequences datasheet from input
        if hasattr(input_obj, 'datasheets'):
            datasheets = input_obj.datasheets
            
            # Check for sequences datasheet
            if hasattr(datasheets, 'sequences'):
                if isinstance(datasheets.sequences,str):
                    return datasheets.sequences
                else:
                    if hasattr(datasheets.sequences,'path'):
                        return datasheets.sequences.path
                    else:
                        raise ValueError("Cannot extract path of sequences datasheet")
            elif hasattr(datasheets, '_datasheets'):
                # Standard BioPipelines format
                for name, info in datasheets._datasheets.items():
                    if 'sequence' in name.lower():
                        return info.path
            elif isinstance(datasheets, dict):
                # Dict format
                for name, info in datasheets.items():
                    if 'sequence' in name.lower():
                        if isinstance(info,str):
                            return info
                        elif isinstance(info, dict) and 'path' in info:
                            return info['path']
                        elif hasattr(info, 'path'):
                            return info.path
        
        # Fallback: try to predict sequences file in output folder
        if hasattr(input_obj, 'output_folder'):
            predicted_files = [
                os.path.join(input_obj.output_folder, 'sequences.csv'),
                os.path.join(input_obj.output_folder, 'datasheet.csv'),
                os.path.join(input_obj.output_folder, 'results.csv')
            ]
            return predicted_files[0]  # Return first prediction
        
        raise ValueError(f"Could not extract sequences path from {input_type} input")
    
    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()
        
        config_lines.extend([
            f"ORIGINAL: {type(self.original_input).__name__}",
            f"MUTANTS: {type(self.mutants_input).__name__}",
            f"INCLUDE_ORIGINAL: {self.include_original}"
        ])
        
        return config_lines
    
    def generate_script(self, script_path: str) -> str:
        """
        Generate script to analyze mutation patterns.
        
        Args:
            script_path: Path where script should be written
            
        Returns:
            Script content as string
        """
        runtime_folder = os.path.dirname(script_path)
        output_folder = self.output_folder
        os.makedirs(output_folder, exist_ok=True)
        
        # Output files
        profile_csv = os.path.join(output_folder, "profile.csv")
        mutations_csv = os.path.join(output_folder, "mutations.csv")
        absolute_freq_csv = os.path.join(output_folder, "absolute_frequencies.csv") 
        relative_freq_csv = os.path.join(output_folder, "relative_frequencies.csv")
        sequence_logo_relative_svg = os.path.join(output_folder, "sequence_logo_relative_frequencies.svg")
        sequence_logo_relative_png = os.path.join(output_folder, "sequence_logo_relative_frequencies.png")
        sequence_logo_absolute_svg = os.path.join(output_folder, "sequence_logo_absolute_frequencies.svg")
        sequence_logo_absolute_png = os.path.join(output_folder, "sequence_logo_absolute_frequencies.png")
        sequence_logo_counts_svg = os.path.join(output_folder, "sequence_logo_counts.svg")
        sequence_logo_counts_png = os.path.join(output_folder, "sequence_logo_counts.png")
        
        # Create config file for mutation profiler
        config_file = os.path.join(output_folder, "mutation_profiler_config.json")
        config_data = {
            "original_sequences": self.original_sequences_path,
            "mutants_sequences": self.mutants_sequences_path,
            "include_original": self.include_original,
            "profile_output": profile_csv,
            "mutations_output": mutations_csv,
            "absolute_frequencies_output": absolute_freq_csv,
            "relative_frequencies_output": relative_freq_csv,
            "sequence_logo_relative_svg": sequence_logo_relative_svg,
            "sequence_logo_relative_png": sequence_logo_relative_png,
            "sequence_logo_absolute_svg": sequence_logo_absolute_svg,
            "sequence_logo_absolute_png": sequence_logo_absolute_png,
            "sequence_logo_counts_svg": sequence_logo_counts_svg,
            "sequence_logo_counts_png": sequence_logo_counts_png
        }
        
        import json
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)
        
        # Generate script content
        script_content = f"""#!/bin/bash
# MutationProfiler execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Analyzing mutation patterns"
echo "Original sequences: {self.original_sequences_path}"
echo "Mutant sequences: {self.mutants_sequences_path}"
echo "Output folder: {output_folder}"

# Run Python mutation profiler script
python "{os.path.join(self.folders['HelpScripts'], 'pipe_mutation_profiler.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Successfully analyzed mutation patterns"
    echo "Profile data: {profile_csv}"
    echo "Mutations data: {mutations_csv}"
    echo "Absolute frequencies: {absolute_freq_csv}"
    echo "Relative frequencies: {relative_freq_csv}"
    echo "Sequence logos: {sequence_logo_relative_svg}, {sequence_logo_absolute_svg}, {sequence_logo_counts_svg}"
else
    echo "Error: Failed to analyze mutation patterns"
    exit 1
fi

{self.generate_completion_check_footer()}
"""
        
        return script_content
    
    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after mutation profiling.
        
        Returns:
            Dictionary with output file paths and datasheet information
        """
        profile_csv = os.path.join(self.output_folder, "profile.csv")
        mutations_csv = os.path.join(self.output_folder, "mutations.csv")
        absolute_freq_csv = os.path.join(self.output_folder, "absolute_frequencies.csv")
        relative_freq_csv = os.path.join(self.output_folder, "relative_frequencies.csv")
        sequence_logo_relative_svg = os.path.join(self.output_folder, "sequence_logo_relative_frequencies.svg")
        sequence_logo_relative_png = os.path.join(self.output_folder, "sequence_logo_relative_frequencies.png")
        sequence_logo_absolute_svg = os.path.join(self.output_folder, "sequence_logo_absolute_frequencies.svg")
        sequence_logo_absolute_png = os.path.join(self.output_folder, "sequence_logo_absolute_frequencies.png")
        sequence_logo_counts_svg = os.path.join(self.output_folder, "sequence_logo_counts.svg")
        sequence_logo_counts_png = os.path.join(self.output_folder, "sequence_logo_counts.png")
        
        # Define datasheets that will be created
        datasheets = {
            "profile": DatasheetInfo(
                name="profile",
                path=profile_csv,
                columns=["position", "original", "count", "frequency"],
                description="Position-wise mutation statistics",
                count=None  # Will be determined at runtime
            ),
            "mutations": DatasheetInfo(
                name="mutations",
                path=mutations_csv,
                columns=["position", "original", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"],
                description="Raw amino acid mutation counts per position",
                count=None
            ),
            "absolute_frequencies": DatasheetInfo(
                name="absolute_frequencies", 
                path=absolute_freq_csv,
                columns=["position", "original", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"],
                description="Absolute amino acid frequencies per position (mutation_count/total_sequences)",
                count=None
            ),
            "relative_frequencies": DatasheetInfo(
                name="relative_frequencies",
                path=relative_freq_csv, 
                columns=["position", "original", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"],
                description="Relative amino acid frequencies per position (mutation_count/mutations_at_position)",
                count=None
            )
        }
        
        return {
            "structures": [],
            "structure_ids": [],
            "compounds": [],
            "compound_ids": [],
            "sequences": [],  # This is an analysis tool, doesn't generate new sequences
            "sequence_ids": [],
            "datasheets": datasheets,
            "visualizations": [sequence_logo_relative_svg, sequence_logo_relative_png, 
                             sequence_logo_absolute_svg, sequence_logo_absolute_png,
                             sequence_logo_counts_svg, sequence_logo_counts_png],
            "output_folder": self.output_folder
        }
    
    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                "include_original": self.include_original
            }
        })
        return base_dict