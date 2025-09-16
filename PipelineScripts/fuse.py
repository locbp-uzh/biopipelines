"""
Fuse configuration for protein fusion sequence generation.

Creates fusion sequences by linking multiple proteins with flexible linkers
of varying lengths, generating all possible combinations for downstream
folding and analysis.
"""

import os
import json
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput


class Fuse(BaseConfig):
    """
    Configuration for protein fusion sequence generation.
    
    Generates all combinations of protein fusions with variable linker lengths,
    producing meaningful sequence IDs and structured output for downstream
    folding tools like AlphaFold.
    """
    
    TOOL_NAME = "Fuse"
    DEFAULT_ENV = "ProteinEnv"
    COMPATIBLE_ENVS = ["ProteinEnv"]
    DEFAULT_RESOURCES = {"gpu": "T4", "memory": "8GB", "time": "2:00:00"}
    
    def __init__(self, input: Union[str, List[str], ToolOutput, Dict[str, Any]] = None,
                 proteins: Union[List[str], str] = None,
                 sequences: Union[List[str], str] = None,
                 name: str = "",
                 linker: str = "GGGGSGGGGSGGGGSGGGGS",
                 linker_lengths: List[str] = None,
                 **kwargs):
        """
        Initialize Fuse configuration.
        
        Args:
            input: Complete standardized input dictionary with sequences, etc.
            proteins: List of protein sequences or PDB file paths
            sequences: Alias for proteins (for consistency with other tools)  
            name: Job name for output files
            linker: Linker sequence to use (will be truncated to specified lengths)
            linker_lengths: List of length ranges for each junction (e.g., ["1-6", "1-6"])
            **kwargs: Additional parameters
        """
        # Handle different input formats
        if input is not None:
            if isinstance(input, StandardizedOutput):
                # StandardizedOutput object (e.g., from upstream tool)
                self.input_proteins = self._extract_sequences_from_standardized(input)
                self.input_is_tool_output = False
                self.standardized_input = input
            elif isinstance(input, ToolOutput):
                # Direct ToolOutput object
                self.input_proteins = self._extract_sequences_from_tool_output(input)
                self.input_is_tool_output = True
                self.standardized_input = None
                self.dependencies.append(input.config)
            elif isinstance(input, dict):
                # Dictionary format with standardized keys
                self.input_proteins = input.get('sequences', [])
                self.input_is_tool_output = False
                self.standardized_input = None
            elif isinstance(input, list):
                # Direct list of sequences
                self.input_proteins = input
                self.input_is_tool_output = False
                self.standardized_input = None
            else:
                # Single sequence string
                self.input_proteins = [str(input)]
                self.input_is_tool_output = False
                self.standardized_input = None
        else:
            # Legacy format: proteins or sequences parameter
            if proteins is not None:
                self.input_proteins = proteins if isinstance(proteins, list) else [proteins]
            elif sequences is not None:
                self.input_proteins = sequences if isinstance(sequences, list) else [sequences]
            else:
                raise ValueError("Must provide either input, proteins, or sequences parameter")
            
            self.input_is_tool_output = False
            self.standardized_input = None
        
        # Store Fuse-specific parameters
        self.name = name or kwargs.get('job_name', '')
        self.linker = linker
        
        # Set default linker_lengths based on number of proteins if not provided
        if linker_lengths is None:
            expected_junctions = len(self.input_proteins) - 1
            self.linker_lengths = ["1-6"] * expected_junctions  # Default range for all junctions
        else:
            self.linker_lengths = linker_lengths
        
        # Validate we have enough proteins for fusion
        if len(self.input_proteins) < 2:
            raise ValueError("Fuse requires at least 2 proteins for fusion")
        
        # Validate linker_lengths matches protein count
        expected_junctions = len(self.input_proteins) - 1
        if len(self.linker_lengths) != expected_junctions:
            raise ValueError(f"linker_lengths must have {expected_junctions} entries for {len(self.input_proteins)} proteins")
        
        # Initialize base class
        super().__init__(**kwargs)
        
        # Initialize file paths (will be set in configure_inputs)
        self._initialize_file_paths()
    
    def _extract_sequences_from_standardized(self, standardized_output: StandardizedOutput) -> List[str]:
        """Extract protein sequences from StandardizedOutput."""
        # If it has sequence files, we'll need to read them at runtime
        # For now, return placeholder that will be handled in script generation
        if hasattr(standardized_output, 'sequences') and standardized_output.sequences:
            return standardized_output.sequences
        elif hasattr(standardized_output, 'structure_ids') and standardized_output.structure_ids:
            # Use structure IDs as protein identifiers
            return standardized_output.structure_ids
        else:
            raise ValueError("Cannot extract sequences from standardized output")
    
    def _extract_sequences_from_tool_output(self, tool_output: ToolOutput) -> List[str]:
        """Extract protein sequences from ToolOutput."""
        # Get sequences from tool output
        sequences = tool_output.get_output_files("sequences")
        if sequences:
            return sequences
        else:
            raise ValueError(f"No sequences found in {tool_output.tool_type} output")
    
    def validate_params(self):
        """Validate Fuse-specific parameters."""
        if not self.input_proteins:
            raise ValueError("input_proteins parameter is required")
        
        if len(self.input_proteins) < 2:
            raise ValueError("Fuse requires at least 2 proteins")
        
        if not self.linker:
            raise ValueError("linker sequence is required")
        
        if not self.linker_lengths:
            raise ValueError("linker_lengths is required")
        
        # Validate linker_lengths format (basic check)
        for length_spec in self.linker_lengths:
            if not isinstance(length_spec, str) or not any(c in length_spec for c in '0123456789'):
                raise ValueError(f"Invalid linker_lengths format: {length_spec}")
    
    def _initialize_file_paths(self):
        """Initialize common file paths used throughout the class."""
        self.queries_csv = None
        self.queries_fasta = None
        
        # Helper script paths
        self.fuse_queries_py = None
    
    def _setup_file_paths(self):
        """Set up all file paths after output_folder is known."""
        # Extract pipeline/job name for file naming
        job_base = self.name or self._extract_job_name()
        
        # Core output files
        self.queries_csv = os.path.join(self.output_folder, f"{job_base}_queries.csv")
        self.queries_fasta = os.path.join(self.output_folder, f"{job_base}_queries.fasta")
        
        # Helper script paths
        if hasattr(self, 'folders') and self.folders:
            self.fuse_queries_py = os.path.join(self.folders["HelpScripts"], "pipe_fuse_queries.py")
        else:
            self.fuse_queries_py = None
    
    def _extract_job_name(self) -> str:
        """Extract job name from output folder structure."""
        # Get job name from parent folder
        # Structure: .../JobName_NNN/N_Fuse -> JobName_NNN
        folder_parts = self.output_folder.split(os.sep)
        for i, part in enumerate(folder_parts):
            if "_Fuse" in part:
                if i == 0:
                    raise ValueError(f"Invalid output folder structure: {self.output_folder}")
                return folder_parts[i-1]
        
        raise ValueError(f"Could not extract job name from output folder: {self.output_folder}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input proteins from various sources."""
        self.folders = pipeline_folders
        self._setup_file_paths()  # Set up all file paths now that we have folders
        
        # For now, input_proteins are handled as-is
        # Runtime script will process them appropriately
        pass
    
    def _predict_sequence_ids(self) -> List[str]:
        """
        Predict sequence IDs that will be generated by fusion combinations.
        
        This allows downstream tools to predict their output file names.
        """
        # Parse linker lengths to get all combinations
        def parse_length_spec(spec: str) -> List[int]:
            """Parse length specification like '1-6' or '3+5-7'."""
            lengths = []
            if '+' in spec:
                for part in spec.split('+'):
                    lengths.extend(parse_length_spec(part.strip()))
            elif '-' in spec and not spec.startswith('-'):
                if spec.count('-') == 1:
                    start, end = map(int, spec.split('-'))
                    lengths.extend(range(start, end + 1))
                else:
                    # Handle negative values or complex ranges
                    lengths.append(int(spec))
            else:
                lengths.append(int(spec))
            return lengths
        
        # Get all length combinations
        length_ranges = [parse_length_spec(spec) for spec in self.linker_lengths]
        
        # Generate meaningful protein names
        protein_names = []
        for i, protein in enumerate(self.input_proteins):
            if isinstance(protein, str) and protein.endswith('.pdb'):
                # PDB file - use filename without extension
                protein_names.append(os.path.splitext(os.path.basename(protein))[0])
            elif isinstance(protein, str) and len(protein) < 50:
                # Short string - might be a name or short sequence
                protein_names.append(f"Protein{i+1}")
            else:
                # Long sequence or other
                protein_names.append(f"Protein{i+1}")
        
        # Generate all combinations
        sequence_ids = []
        
        def generate_combinations(current_name: str, junction_idx: int) -> None:
            if junction_idx >= len(length_ranges):
                # Base case - we've processed all junctions
                sequence_ids.append(current_name)
                return
            
            # Recursive case - add next protein with all possible linker lengths
            next_protein = protein_names[junction_idx + 1]
            for length in length_ranges[junction_idx]:
                if length >= 0:
                    new_name = f"{current_name}_Linker{length}_{next_protein}"
                else:
                    # Negative length - truncation
                    new_name = f"{current_name}{length}_{next_protein}"
                
                generate_combinations(new_name, junction_idx + 1)
        
        # Start with first protein
        generate_combinations(protein_names[0], 0)
        
        return sequence_ids
    
    def get_config_display(self) -> List[str]:
        """Get Fuse configuration display lines."""
        config_lines = super().get_config_display()
        
        config_lines.extend([
            f"PROTEINS: {len(self.input_proteins)} proteins",
            f"LINKER: {self.linker[:20]}{'...' if len(self.linker) > 20 else ''}",
            f"LINKER_LENGTHS: {', '.join(self.linker_lengths)}"
        ])
        
        return config_lines
    
    def generate_script(self, script_path: str) -> str:
        """
        Generate Fuse execution script.
        
        Args:
            script_path: Path where script should be written
            
        Returns:
            Script content as string
        """
        # Generate script content following modular pattern
        script_content = "#!/bin/bash\n"
        script_content += "# Fuse execution script\n"
        script_content += "# Generated by BioPipelines system\n\n"
        script_content += self.generate_completion_check_header()
        script_content += self.generate_script_prepare_proteins()
        script_content += self.generate_script_run_fuse()
        script_content += self.generate_completion_check_footer()
        
        return script_content
    
    def generate_script_prepare_proteins(self) -> str:
        """Generate the protein preparation part of the script."""
        # Build proteins string for fuse_queries.py
        if self.input_is_tool_output:
            # Need to extract sequences from upstream tool output at runtime
            return f"""echo "Extracting proteins from upstream tool output"
# TODO: Add logic to extract sequences from input files
PROTEINS_STR="{';'.join(self.input_proteins[:3])}"  # Placeholder
echo "Using proteins: $PROTEINS_STR"

"""
        else:
            # Direct protein sequences or file paths
            proteins_str = ";".join(self.input_proteins)
            return f"""echo "Using direct protein input"
PROTEINS_STR="{proteins_str}"
echo "Proteins: $PROTEINS_STR"

"""
    
    def generate_script_run_fuse(self) -> str:
        """Generate the fusion sequence generation part of the script."""
        # Format linker lengths for fuse_queries.py (needs L...L wrapping)
        linker_lengths_str = f"L{';'.join(self.linker_lengths)}L"
        
        job_base = self.name or "fuse"
        
        return f"""echo "Generating fusion sequence combinations"
echo "Linker: {self.linker}"
echo "Linker lengths: {linker_lengths_str}"

# Call fuse_queries.py to generate all fusion combinations
python {self.fuse_queries_py} "{job_base}" "$PROTEINS_STR" "{self.linker}" "{linker_lengths_str}" "{self.queries_csv}" "{self.queries_fasta}"

# Check if files were created successfully
if [ ! -f "{self.queries_csv}" ]; then
    echo "ERROR: Failed to generate queries CSV file"
    exit 1
fi

if [ ! -f "{self.queries_fasta}" ]; then
    echo "ERROR: Failed to generate queries FASTA file"
    exit 1
fi

echo "Successfully generated fusion sequences:"
echo "CSV file: {self.queries_csv}"
echo "FASTA file: {self.queries_fasta}"

# Show summary
NUM_SEQUENCES=$(tail -n +2 "{self.queries_csv}" | wc -l)
echo "Generated $NUM_SEQUENCES fusion sequence combinations"

"""
    
    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after Fuse execution.
        
        Uses pure path construction - no filesystem access.
        Returns expected paths based on fusion sequence generation.
        
        Returns:
            Dictionary mapping output types to file paths with standard keys:
            - structures: Empty (no structures from Fuse)
            - compounds: Empty (no compounds from Fuse) 
            - sequences: Main CSV and FASTA files
            - datasheets: CSV file with fusion metadata
            - output_folder: Tool's output directory
            - sequence_ids: Predicted sequence IDs for downstream tools
        """
        # Ensure file paths are set up
        if not hasattr(self, 'queries_csv') or self.queries_csv is None:
            self._setup_file_paths()
        
        queries_csv = self.queries_csv
        queries_fasta = self.queries_fasta
        
        # Predict sequence IDs for downstream tools
        sequence_ids = self._predict_sequence_ids()
        
        # Organize datasheets by content type
        datasheets = {
            "sequences": {
                "path": queries_csv,
                "columns": ["id", "sequence", "lengths"],
                "description": "Fusion sequences with linker length information",
                "count": len(sequence_ids)
            }
        }
        
        return {
            "structures": [],
            "structure_ids": [],
            "compounds": [],
            "compound_ids": [],
            "sequences": [queries_csv, queries_fasta],
            "sequence_ids": sequence_ids,
            "datasheets": datasheets,
            "output_folder": self.output_folder,
            # Keep legacy aliases for compatibility
            "queries_csv": [queries_csv],
            "queries_fasta": [queries_fasta]
        }
    
    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including Fuse-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "fuse_params": {
                "name": self.name,
                "linker": self.linker,
                "linker_lengths": self.linker_lengths,
                "num_proteins": len(self.input_proteins),
                "input_type": "tool_output" if self.input_is_tool_output else "direct"
            }
        })
        return base_dict