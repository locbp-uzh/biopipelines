"""
AlphaFold2/ColabFold configuration for protein structure prediction.

Handles sequence folding with AlphaFold2 using ColabFold implementation,
integrating with upstream sequence generation tools and providing comprehensive
ranking and analysis capabilities.
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


class AlphaFold(BaseConfig):
    """
    Configuration for AlphaFold2/ColabFold structure prediction.
    
    Predicts protein structures from amino acid sequences with support
    for various input sources, advanced parameters, and comprehensive analysis.
    """
    
    TOOL_NAME = "AlphaFold"
    DEFAULT_ENV = None  # Loaded from config.yaml
    
    def __init__(self, sequences: Union[str, List[str], ToolOutput, Dict[str, Any]] = None,
                 datasheets: Optional[List[str]] = None,
                 name: str = "",
                 num_relax: int = 0, num_recycle: int = 3, rand_seed: int = 0,
                 **kwargs):
        """
        Initialize AlphaFold configuration.

        Args:
            sequences: Input sequences - can be FASTA file, CSV file, list, ToolOutput, or dict with sequences
            datasheets: Input datasheet files for metadata
            name: Job name for output files
            rank: Rank structures by confidence metrics
            only_first: Only use first model for ranking comparisons
            num_relax: Number of best models to relax with AMBER
            num_recycle: Number of recycling iterations (default 3)
            rand_seed: Random seed for reproducible results (0 = random)
            **kwargs: Additional parameters
        """
        # Handle different input formats
        if sequences is not None:
            if isinstance(sequences, StandardizedOutput):
                # StandardizedOutput object (e.g., from pmpnn)
                self.input_sequences = sequences.sequences
                self.input_datasheets = sequences.datasheets
                self.input_is_tool_output = False  # Direct file paths now
                self.standardized_input = sequences  # Keep reference to get sequence_ids later
            elif isinstance(sequences, ToolOutput):
                # Direct ToolOutput object
                self.input_sequences = sequences
                self.input_datasheets = sequences.get_output_files("datasheets")
                self.input_is_tool_output = True
                self.standardized_input = None
            elif isinstance(sequences, dict):
                # Dictionary format with standardized keys
                self.input_sequences = sequences.get('sequences', [])
                self.input_datasheets = sequences.get('datasheets', {})
                self.input_is_tool_output = False  # Direct file paths
                self.standardized_input = None
            else:
                # Direct sequence(s) - string or list
                self.input_sequences = sequences
                self.input_datasheets = datasheets or {}
                self.input_is_tool_output = isinstance(sequences, ToolOutput)
                self.standardized_input = None
        else:
            raise ValueError("sequences parameter is required")
        
        # Store AlphaFold-specific parameters
        self.name = name or kwargs.get('job_name', '')
        self.num_relax = num_relax
        self.num_recycle = num_recycle
        self.rand_seed = rand_seed
        
        # Initialize base class
        super().__init__(**kwargs)
        
        # Initialize file paths (will be set in configure_inputs)
        self._initialize_file_paths()
    
    def validate_params(self):
        """Validate AlphaFold-specific parameters."""
        if not self.input_sequences:
            raise ValueError("input sequences parameter is required")
        
        if self.num_relax < 0:
            raise ValueError("num_relax cannot be negative")
        
        if self.num_recycle < 1:
            raise ValueError("num_recycle must be at least 1")
        
        if self.rand_seed < 0:
            raise ValueError("rand_seed cannot be negative")
    
    def _initialize_file_paths(self):
        """Initialize common file paths used throughout the class."""
        self.queries_csv = None
        self.queries_fasta = None
        
        # Helper script paths
        self.colabfold_batch = None
        self.fa_to_csv_fasta_py = None
    
    def _setup_file_paths(self):
        """Set up all file paths after output_folder is known."""
        # Extract pipeline/job name for file naming
        job_base = self.name or self._extract_job_name()
        
        # Core input/output files
        self.queries_csv = os.path.join(self.output_folder, f"{job_base}_queries.csv")
        self.queries_fasta = os.path.join(self.output_folder, f"{job_base}_queries.fasta")
        
        # Helper script paths (only set if folders are available)
        if hasattr(self, 'folders') and self.folders:
            self.colabfold_batch = os.path.join(self.folders["AlphaFold"], "colabfold-conda/bin/colabfold_batch")
            self.fa_to_csv_fasta_py = os.path.join(self.folders["HelpScripts"], "pipe_fa_to_csv_fasta.py")
        else:
            # Temporary placeholders when folders aren't available yet
            self.colabfold_batch = None
            self.fa_to_csv_fasta_py = None
    
    def _extract_job_name(self) -> str:
        """Extract job name from output folder structure."""
        # Get job name from parent folder
        # Structure: .../JobName_NNN/N_AlphaFold -> JobName_NNN
        folder_parts = self.output_folder.split(os.sep)
        for i, part in enumerate(folder_parts):
            if "_AlphaFold" in part:
                if i == 0:
                    raise ValueError(f"Invalid output folder structure: {self.output_folder}")
                return folder_parts[i-1]
        
        raise ValueError(f"Could not extract job name from output folder: {self.output_folder}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input sequences from various sources."""
        self.folders = pipeline_folders
        self._setup_file_paths()  # Set up all file paths now that we have folders
        
        if self.input_is_tool_output:
            # Input from previous tool (e.g., ProteinMPNN)
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
            self.input_sources["sequences"] = source_sequences[0]  # Use first sequence file
            
            # Add dependency
            self.dependencies.append(tool_output.config)
            
        elif isinstance(self.input_sequences, list):
            # Direct list of sequence file paths (from StandardizedOutput)
            if self.input_sequences:
                # Use first sequence file
                self.input_sources["sequences"] = self.input_sequences[0]
            else:
                raise ValueError("Empty sequence list provided")
                
        elif isinstance(self.input_sequences, str):
            # String input - could be file or direct sequence
            if self.input_sequences.endswith('.csv'):
                # CSV file with id,sequence columns
                csv_source = os.path.join(pipeline_folders.get("data", os.getcwd()), self.input_sequences)
                if os.path.exists(csv_source):
                    self.input_sources["sequences"] = csv_source
                else:
                    raise ValueError(f"CSV file not found: {csv_source}")
            elif self.input_sequences.endswith('.fasta') or self.input_sequences.endswith('.fa'):
                # FASTA file
                fasta_source = os.path.join(pipeline_folders.get("data", os.getcwd()), self.input_sequences)
                if os.path.exists(fasta_source):
                    self.input_sources["sequences"] = fasta_source
                else:
                    raise ValueError(f"FASTA file not found: {fasta_source}")
            else:
                # Direct sequence - create CSV file
                if not self.name:
                    raise ValueError("name parameter required when providing direct sequence")
                
                # Will create CSV in generate_script
                self.input_sources["direct_sequence"] = self.input_sequences
        else:
            raise ValueError(f"Unsupported input type: {type(self.input_sequences)}")
    
    def get_config_display(self) -> List[str]:
        """Get AlphaFold configuration display lines."""
        config_lines = super().get_config_display()
        
        # Input information
        if self.input_is_tool_output:
            config_lines.append(f"INPUT: {self.input_sequences.tool_type} output")
        else:
            config_lines.append(f"INPUT: {self.input_sequences}")
        
        config_lines.extend([
            f"NUM RELAX: {self.num_relax}",
            f"NUM RECYCLE: {self.num_recycle}"
        ])
        
        if self.rand_seed > 0:
            config_lines.append(f"RAND SEED: {self.rand_seed}")
        
        return config_lines

    def generate_script(self, script_path: str) -> str:
        """
        Generate AlphaFold execution script.
        
        Args:
            script_path: Path where script should be written
            
        Returns:
            Script content as string
        """
        # Generate script content following modular pattern
        script_content = "#!/bin/bash\n"
        script_content += "# AlphaFold execution script\n"
        script_content += "# Generated by ProteinNotebooks pipeline system\n\n"
        script_content += self.generate_completion_check_header()
        script_content += self.generate_script_prepare_sequences()
        script_content += self.generate_script_run_alphafold()
        script_content += self.generate_script_extract_best_rank()
        script_content += self.generate_completion_check_footer()
        
        return script_content
    
    def generate_script_prepare_sequences(self) -> str:
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
            
            # Check if we need to convert from FASTA folder to CSV
            if source_file.endswith("seqs") or "seqs" in source_file:
                # ProteinMPNN output - convert .fa files to CSV
                return f"""echo "Converting ProteinMPNN .fa files to queries CSV"
python {self.fa_to_csv_fasta_py} {source_file} {self.queries_csv} {self.queries_fasta}

"""
            else:
                # Direct CSV or FASTA file - copy/link with error checking
                return f"""echo "Using sequences from: {source_file}"
if [ -f "{source_file}" ]; then
    cp "{source_file}" "{self.queries_csv}"
    echo "Successfully copied sequences file"
else
    echo "ERROR: Sequence file not found: {source_file}"
    echo "This usually means the previous step (ProteinMPNN) failed to generate the expected output"
    exit 1
fi

"""
        else:
            raise ValueError("No sequence input configured")
    
    def generate_script_run_alphafold(self) -> str:
        """Generate the AlphaFold execution part of the script."""
        # Build AlphaFold options
        af_options = ""
        if self.num_relax > 0:
            af_options += f" --amber --use-gpu-relax --num-relax {self.num_relax}"
        if self.num_recycle != 3:
            af_options += f" --num-recycle {self.num_recycle}"
        if self.rand_seed != 0:
            af_options += f" --random-seed {self.rand_seed}"
        
        folding_folder = os.path.join(self.output_folder, "Folding")
        
        return f"""echo "Running AlphaFold2/ColabFold"
echo "Options: {af_options}"
echo "Output folder: {self.output_folder}"

# Create Folding subfolder for raw ColabFold outputs
mkdir -p "{folding_folder}"

# Run ColabFold batch - output to Folding subfolder
{self.colabfold_batch} {self.queries_csv} "{folding_folder}" {af_options}

"""
    def generate_script_extract_best_rank(self) -> str:
        folding_folder = os.path.join(self.output_folder, "Folding")
        
        return f"""echo "Extracting best structures from Folding subfolder"
# AlphaFold creates files like: 
# - sequenceid_unrelaxed_rank_001_alphafold2_ptm_model_N_seed_SSS.pdb
# - sequenceid_relaxed_rank_001_alphafold2_ptm_model_N_seed_SSS.pdb
# We want to copy the best ones to main directory as: sequenceid.pdb

cd "{folding_folder}"

# Handle relaxed format (preferred if both exist)
for file in *_relaxed_rank_001_*.pdb; do
    if [ -f "$file" ]; then
        # Extract sequence ID (everything before _relaxed_rank_001)
        base=$(echo "$file" | sed 's/_relaxed_rank_001_.*/.pdb/')
        cp "$file" "{self.output_folder}/$base"
        echo "Extracted: $file -> $base"
    fi
done

# Handle unrelaxed format (only if relaxed doesn't exist)
for file in *_unrelaxed_rank_001_*.pdb; do
    if [ -f "$file" ]; then
        # Extract sequence ID (everything before _unrelaxed_rank_001)
        base=$(echo "$file" | sed 's/_unrelaxed_rank_001_.*/.pdb/')
        # Only copy if the relaxed version doesn't already exist
        if [ ! -f "{self.output_folder}/$base" ]; then
            cp "$file" "{self.output_folder}/$base"
            echo "Extracted: $file -> $base"
        else
            echo "Skipped $file (relaxed version already exists as $base)"
        fi
    fi
done

cd - > /dev/null

"""

    
    def _predict_structure_outputs(self) -> List[str]:
        """
        Predict structure files that AlphaFold will generate.
        
        Based on input sequences and ColabFold naming patterns.
        """
        structure_files = []
        sequence_ids = []
        
        # AlphaFold/ColabFold generates files with specific patterns
        # Format: {sequence_id}_unrelaxed_rank_{rank}_alphafold2_ptm.pdb
        
        # Try to get sequence IDs from standardized input (pmpnn)
        if hasattr(self, 'standardized_input') and self.standardized_input:
            if hasattr(self.standardized_input, 'sequence_ids'):
                sequence_ids = self.standardized_input.sequence_ids
            elif 'sequence_ids' in self.standardized_input._data:
                sequence_ids = self.standardized_input._data['sequence_ids']
        
        # Try to get sequence info from direct ToolOutput
        if not sequence_ids and self.input_is_tool_output:
            # Get sequence IDs from upstream tool if possible
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
                # Use simple naming: just the sequence ID + .pdb
                pdb_path = os.path.join(self.output_folder, f"{seq_id}.pdb")
                structure_files.append(pdb_path)
        
        # Must have structure files from dependencies
        if not structure_files:
            raise ValueError("Could not determine sequence IDs for structure prediction - no upstream sequence data found")
        
        return structure_files
    
    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after AlphaFold execution.
        
        Uses pure path construction - no filesystem access.
        Returns expected paths based on AlphaFold/ColabFold output patterns.
        
        Returns:
            Dictionary mapping output types to file paths with standard keys:
            - structures: PDB structure files
            - compounds: Empty (no compounds from AlphaFold)
            - sequences: Input queries CSV file (main sequence source)
            - datasheets: Empty for now (ranking will be added later)
            - output_folder: Tool's output directory
        """
        # Ensure file paths are set up
        if not hasattr(self, 'queries_csv') or self.queries_csv is None:
            # Fallback if configure_inputs hasn't been called yet
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
        
        # Organize datasheets by content type (empty for now - ranking could be added later)
        datasheets = {
            "structures": {
                "path": queries_csv,  # Input sequences that generated the structures
                "columns": ["id", "source_id", "sequence"],
                "description": "Input sequences used for AlphaFold structure prediction",
                "count": len(structure_ids)
            }
        }
        
        return {
            "structures": structure_files,
            "structure_ids": structure_ids,
            "compounds": [],
            "compound_ids": [],
            "sequences": [queries_csv],  # Main sequence input/output
            "sequence_ids": sequence_ids,
            "datasheets": datasheets,
            "output_folder": self.output_folder,
            # Keep legacy aliases for compatibility
            "main": queries_csv,  # Legacy alias for backward compatibility
            "pdbs": structure_files,
            "queries_csv": [queries_csv]
        }
    
    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including AlphaFold-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "af_params": {
                "name": self.name,
                "num_relax": self.num_relax,
                "num_recycle": self.num_recycle,
                "rand_seed": self.rand_seed,
                "input_type": "tool_output" if self.input_is_tool_output else "direct"
            }
        })
        return base_dict