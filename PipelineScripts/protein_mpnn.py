"""
ProteinMPNN configuration for sequence design from protein structures.

Handles sequence generation with fixed position constraints, pLDDT thresholds,
and automatic integration with upstream structure generation tools.
"""

import os
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput


class ProteinMPNN(BaseConfig):
    """
    Configuration for ProteinMPNN sequence design.
    
    Generates protein sequences for given backbone structures with support
    for fixed positions, confidence-based constraints, and flexible input sources.
    """
    
    TOOL_NAME = "ProteinMPNN" 
    DEFAULT_ENV = "ProteinEnv"
    
    def __init__(self, input: Union[str, List[str], ToolOutput, Dict[str, Any]] = None,
                 structures: Union[str, List[str], ToolOutput] = "",
                 datasheets: Optional[List[str]] = None,
                 num_sequences: int = 1, fixed_positions: str = "",
                 designed_positions: str = "", fixed_chain: str = "A", 
                 plddt_threshold: float = 100.0, sampling_temp: float = 0.1, 
                 model_name: str = "v_48_020", soluble_model: bool = True, **kwargs):
        """
        Initialize ProteinMPNN configuration.
        
        Args:
            input: Complete standardized input dictionary with structures, datasheets, etc.
            structures: Input structures (PDB files, folder, or ToolOutput)
            datasheets: Input datasheet files for metadata
            num_sequences: Number of sequences to generate per structure
            fixed_positions: PyMOL-style selection or datasheet reference (e.g., "structures.fixed")
            designed_positions: PyMOL-style selection or datasheet reference (e.g., "structures.designed")
            fixed_chain: Chain to apply fixed positions to
            plddt_threshold: pLDDT threshold for automatic fixing (100 = no fixing)
            sampling_temp: Sampling temperature for sequence generation
            model_name: ProteinMPNN model variant to use
            soluble_model: Use soluble protein model
            **kwargs: Additional parameters
        """
        # Handle standardized input format
        if input is not None:
            if isinstance(input, StandardizedOutput):
                # StandardizedOutput object (e.g., rfd)
                self.input_structures = input.structures
                self.input_datasheets = input.datasheets
                self.input_is_tool_output = False  # Direct file paths now
                self.standardized_input = input  # Keep reference for metadata
            elif isinstance(input, ToolOutput):
                # Direct ToolOutput object
                self.input_structures = input
                self.input_datasheets = input.get_output_files("datasheets")
                self.input_is_tool_output = True
                self.standardized_input = None
            elif isinstance(input, dict):
                # Dictionary format with standardized keys
                self.input_structures = input.get('structures', [])
                self.input_datasheets = input.get('datasheets', {})
                self.input_is_tool_output = False  # Direct file paths
                self.standardized_input = None
            else:
                # Fallback to treating as input_structures
                self.input_structures = input
                self.input_datasheets = datasheets or {}
                self.input_is_tool_output = isinstance(input, ToolOutput)
                self.standardized_input = None
        else:
            # Legacy format: structures=previous_tool
            self.input_structures = structures
            self.input_datasheets = datasheets or {}
            self.input_is_tool_output = isinstance(structures, ToolOutput)
            self.standardized_input = None
        
        # Store ProteinMPNN-specific parameters
        self.num_sequences = num_sequences
        self.fixed_positions = fixed_positions
        self.designed_positions = designed_positions
        self.fixed_chain = fixed_chain
        self.plddt_threshold = plddt_threshold
        self.sampling_temp = sampling_temp
        self.model_name = model_name
        self.soluble_model = soluble_model
        
        # Track input source type
        self.input_is_tool_output = isinstance(structures, ToolOutput)
        self.input_is_folder = False
        self.input_pdb_files = []
        
        # Initialize base class
        super().__init__(**kwargs)

        # Initialize file paths (will be set in configure_inputs)
        self._initialize_file_paths()

        # Set up datasheets attribute for IDE autocompletion
        self._setup_datasheets_for_ide()

    def _setup_datasheets_for_ide(self):
        """Set up datasheets attribute with predefined columns for IDE autocompletion."""
        from .base_config import DatasheetContainer, DatasheetInfo

        # Create temporary DatasheetInfo objects with known columns for IDE support
        sequences_datasheet = DatasheetInfo(
            name="sequences",
            path="",  # Path will be set when output_folder is known
            columns=["id", "source_id", "source_pdb", "sequence", "score", "seq_recovery", "rmsd"],
            description="ProteinMPNN sequence generation results with scores and structure recovery metrics"
        )

        # Set up datasheets container for IDE autocompletion
        self.datasheets = DatasheetContainer({"sequences": sequences_datasheet})

        # CRITICAL: Explicitly set datasheet attributes for IDE autocompletion
        self.datasheets.sequences = sequences_datasheet

    def validate_params(self):
        """Validate ProteinMPNN-specific parameters."""
        if not self.input_structures:
            raise ValueError("input_structures parameter is required")
        
        if self.num_sequences <= 0:
            raise ValueError("num_sequences must be positive")
        
        if self.sampling_temp <= 0:
            raise ValueError("sampling_temp must be positive")
        
        if self.plddt_threshold < 0 or self.plddt_threshold > 100:
            raise ValueError("plddt_threshold must be between 0 and 100")
        
        # Validate datasheet references if provided
        if self.fixed_positions:
            self.validate_datasheet_reference(self.fixed_positions)
        if self.designed_positions:
            self.validate_datasheet_reference(self.designed_positions)
        
        # Validate model name
        valid_models = ["v_48_002", "v_48_010", "v_48_020", "v_48_030"]
        if self.model_name not in valid_models:
            raise ValueError(f"model_name must be one of: {valid_models}")
    
    def _initialize_file_paths(self):
        """Initialize common file paths used throughout the class."""
        self.parsed_pdbs_jsonl = None
        self.fixed_jsonl = None
        self.sele_csv = None
        self.seqs_folder = None
        self.main_datasheet = None
        self.queries_csv = None
        self.queries_fasta = None
        
        # Helper script paths
        self.parse_py = None
        self.fixed_py = None
        self.pmpnn_py = None
        self.datasheet_py = None
    
    def _setup_file_paths(self):
        """Set up all file paths after output_folder is known."""
        # Core output files
        self.parsed_pdbs_jsonl = os.path.join(self.output_folder, "parsed_pdbs.jsonl")
        self.fixed_jsonl = os.path.join(self.output_folder, "fixed_pos.jsonl")
        self.sele_csv = os.path.join(self.output_folder, "fixed_designed.csv")
        self.seqs_folder = os.path.join(self.output_folder, "seqs")
        self.main_datasheet = os.path.join(self.output_folder, "proteinmpnn_results.csv")
        
        # Pipeline-based files
        pipeline_name = self._extract_pipeline_name()
        self.queries_csv = os.path.join(self.output_folder, f"{pipeline_name}_queries.csv")
        self.queries_fasta = os.path.join(self.output_folder, f"{pipeline_name}_queries.fasta")
        
        # Helper script paths (only set if folders are available)
        if hasattr(self, 'folders') and self.folders:
            self.parse_py = os.path.join(self.folders["ProteinMPNN"], "helper_scripts", "parse_multiple_chains.py")
            self.fixed_py = os.path.join(self.folders["HelpScripts"], "pipe_pmpnn_fixed_positions.py")
            self.pmpnn_py = os.path.join(self.folders["ProteinMPNN"], "protein_mpnn_run.py")
            self.datasheet_py = os.path.join(self.folders["HelpScripts"], "pipe_pmpnn_datasheet.py")
            self.fa_to_csv_fasta_py = os.path.join(self.folders["HelpScripts"], "pipe_fa_to_csv_fasta.py")
        else:
            # Temporary placeholders when folders aren't available yet
            self.parse_py = None
            self.fixed_py = None
            self.pmpnn_py = None
            self.datasheet_py = None
            self.fa_to_csv_fasta_py = None
    
    def _extract_pipeline_name(self) -> str:
        """Extract pipeline name from output folder structure."""
        folder_parts = self.output_folder.split(os.sep)
        for i, part in enumerate(folder_parts):
            if "ProteinMPNN" in part:
                if i > 0:
                    return folder_parts[i-1]
                break
        raise ValueError(f"Could not extract pipeline name from output folder: {self.output_folder}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input structures from various sources."""
        self.folders = pipeline_folders
        self._setup_file_paths()  # Set up all file paths now that we have folders
        runtime_folder = os.path.join(self.output_folder, "RunTime")
        
        if self.input_is_tool_output:
            # Input from previous tool (e.g., RFdiffusion)
            tool_output: ToolOutput = self.input_structures
            source_pdbs = tool_output.get_output_files("pdbs")
            
            if not source_pdbs:
                raise ValueError(f"No PDB outputs found from {tool_output.tool_type}")
            
            # Add dependency
            self.dependencies.append(tool_output.config)
            
            # Copy or link to runtime folder
            self.input_pdb_files = []
            for pdb_path in source_pdbs:
                pdb_name = os.path.basename(pdb_path)
                runtime_pdb = os.path.join(runtime_folder, pdb_name)
                self.input_pdb_files.append(runtime_pdb)
                self.input_sources[pdb_name] = pdb_path
                
        elif isinstance(self.input_structures, list):
            # Direct list of PDB file paths (from rfd.structures)
            self.input_pdb_files = []
            for pdb_path in self.input_structures:
                pdb_name = os.path.basename(pdb_path)
                runtime_pdb = os.path.join(runtime_folder, pdb_name)
                self.input_pdb_files.append(runtime_pdb)
                
                # Only check existence for non-pipeline inputs (e.g., user files)
                # Pipeline-generated files don't exist yet during configuration
                if hasattr(self, 'is_pipeline_input') and self.is_pipeline_input:
                    # Pipeline input - don't check existence, files will be created during execution
                    self.input_sources[pdb_name] = pdb_path
                else:
                    # Direct file input - check existence
                    if os.path.exists(pdb_path):
                        self.input_sources[pdb_name] = pdb_path
                    else:
                        raise ValueError(f"PDB file not found: {pdb_path}")
        
        elif isinstance(self.input_structures, str):
            # String input - could be file or folder
            if self.input_structures.endswith('.pdb'):
                # Single PDB file
                pdb_source = os.path.join(pipeline_folders["PDBs"], self.input_structures)
                if os.path.exists(pdb_source):
                    runtime_pdb = os.path.join(runtime_folder, self.input_structures)
                    self.input_pdb_files = [runtime_pdb]
                    self.input_sources[self.input_structures] = pdb_source
                else:
                    raise ValueError(f"PDB file not found: {pdb_source}")
            else:
                # Folder of PDB files
                folder_path = os.path.join(os.getcwd(), self.input_structures)
                if os.path.exists(folder_path):
                    pdb_files = [f for f in os.listdir(folder_path) if f.endswith('.pdb')]
                    if not pdb_files:
                        raise ValueError(f"No PDB files found in folder: {folder_path}")
                    
                    self.input_is_folder = True
                    self.input_pdb_files = []
                    for pdb_file in pdb_files:
                        source_path = os.path.join(folder_path, pdb_file)
                        runtime_pdb = os.path.join(runtime_folder, pdb_file)
                        self.input_pdb_files.append(runtime_pdb)
                        self.input_sources[pdb_file] = source_path
                else:
                    raise ValueError(f"Input folder not found: {folder_path}")
    
    def get_config_display(self) -> List[str]:
        """Get ProteinMPNN configuration display lines.""" 
        config_lines = super().get_config_display()
        
        config_lines.extend([
            f"NUM SEQUENCES PER TARGET: {self.num_sequences}",
            f"FIXED: {self.fixed_positions or 'None'}",
            f"DESIGNED: {self.designed_positions or 'None'}",
            f"FIXED CHAIN: {self.fixed_chain}",
            f"pLDDT THR: {self.plddt_threshold}",
            f"SAMPLING T: {self.sampling_temp}",
            f"MODEL: {self.model_name}",
            f"SOLUBLE: {self.soluble_model}"
        ])
        
        return config_lines

    def generate_script(self, script_path: str) -> str:
        """
        Generate ProteinMPNN execution script.
        
        Args:
            script_path: Path where script should be written
            
        Returns:
            Script content as string
        """
        
        # Generate script content following modular pattern
        script_content = "#!/bin/bash\n"
        script_content += "# ProteinMPNN execution script\n"
        script_content += "# Generated by ProteinNotebooks pipeline system\n\n"
        script_content += self.generate_completion_check_header()
        script_content += self.generate_script_prepare_inputs()
        script_content += self.generate_script_run_proteinmpnn()
        script_content += self.generate_script_create_datasheet()
        script_content += self.generate_completion_check_footer()
        
        return script_content
    
    def generate_script_prepare_inputs(self) -> str:
        """Generate the input preparation part of the script."""
        # Determine input directory for parsing
        first_source = list(self.input_sources.values())[0]
        input_directory = os.path.dirname(first_source)
        
        # Determine input source and parameters for fixed positions script
        if (self.input_is_tool_output and hasattr(self, 'input_datasheets') and self.input_datasheets) or \
           (hasattr(self, 'is_pipeline_input') and self.is_pipeline_input and hasattr(self, 'input_datasheets') and self.input_datasheets):
            # Use datasheet from previous tool (e.g., RFdiffusion)
            input_source = "datasheet"
            if isinstance(self.input_datasheets, dict):
                # New named datasheet format - look for structures first (RFdiffusion), then main (legacy)
                if "structures" in self.input_datasheets:
                    input_datasheet = self.input_datasheets["structures"]["path"]
                elif "main" in self.input_datasheets:
                    input_datasheet = self.input_datasheets["main"]["path"]
                else:
                    # Use first available datasheet
                    first_key = next(iter(self.input_datasheets))
                    input_datasheet = self.input_datasheets[first_key]["path"]
            elif hasattr(self.input_datasheets, '_datasheets'):
                # DatasheetContainer object - get the path properly
                if 'structures' in self.input_datasheets._datasheets:
                    input_datasheet = self.input_datasheets._datasheets['structures'].path
                elif 'main' in self.input_datasheets._datasheets:
                    input_datasheet = self.input_datasheets._datasheets['main'].path
                else:
                    # Fallback to first available datasheet
                    first_name = list(self.input_datasheets._datasheets.keys())[0]
                    input_datasheet = self.input_datasheets._datasheets[first_name].path
            else:
                # Legacy format
                input_datasheet = self.input_datasheets[0] if isinstance(self.input_datasheets, list) else str(self.input_datasheets)
        elif self.fixed_positions and self.fixed_positions != "-":
            # Use direct fixed positions selection
            input_source = "selection"
            input_datasheet = "-"
        else:
            # Use pLDDT threshold method
            input_source = "plddt"
            input_datasheet = "-"
        
        # Resolve datasheet references in fixed/designed positions
        fixed_param = self.resolve_datasheet_reference(self.fixed_positions) if self.fixed_positions else "-"
        designed_param = self.resolve_datasheet_reference(self.designed_positions) if self.designed_positions else "-"
        
        return f"""echo "Determining fixed positions"
python {self.fixed_py} {input_directory} {input_source} {input_datasheet} {self.plddt_threshold} {fixed_param} {designed_param} {self.fixed_chain} {self.fixed_jsonl} {self.sele_csv}

echo "Parsing multiple PDBs" 
python {self.parse_py} --input_path {input_directory} --output_path {self.parsed_pdbs_jsonl}

"""
    
    def generate_script_run_proteinmpnn(self) -> str:
        """Generate the ProteinMPNN execution part of the script."""
        mpnn_job_folder = self.output_folder
        
        # Build ProteinMPNN options
        pmpnn_options = f"--num_seq_per_target {self.num_sequences}"
        pmpnn_options += f" --sampling_temp {self.sampling_temp}"
        pmpnn_options += f" --model_name {self.model_name}"
        
        if self.soluble_model:
            pmpnn_options += " --use_soluble_model"
        
        return f"""echo "Running model"
echo "Options: {pmpnn_options}"
echo "Output folder: {mpnn_job_folder}"

# Run ProteinMPNN
python {self.pmpnn_py} --jsonl_path {self.parsed_pdbs_jsonl} --fixed_positions_jsonl {self.fixed_jsonl} --out_folder {self.output_folder} {pmpnn_options}

"""

    def generate_script_create_datasheet(self) -> str:
        """Generate the datasheet creation part of the script."""
        pipeline_name = self._extract_pipeline_name()
        
        # Determine input datasheet to inherit columns from
        input_datasheet = "-"  # Default: no input datasheet
        if self.input_is_tool_output and hasattr(self, 'input_datasheets') and self.input_datasheets:
            if isinstance(self.input_datasheets, dict):
                # New named datasheet format - look for structures first (RFdiffusion), then main (legacy)
                if "structures" in self.input_datasheets:
                    input_datasheet = self.input_datasheets["structures"]["path"]
                elif "main" in self.input_datasheets:
                    input_datasheet = self.input_datasheets["main"]["path"]
                else:
                    # Use first available datasheet
                    first_key = next(iter(self.input_datasheets))
                    input_datasheet = self.input_datasheets[first_key]["path"]
            elif hasattr(self.input_datasheets, '_datasheets'):
                # DatasheetContainer object - get the path properly
                if 'structures' in self.input_datasheets._datasheets:
                    input_datasheet = self.input_datasheets._datasheets['structures'].path
                elif 'main' in self.input_datasheets._datasheets:
                    input_datasheet = self.input_datasheets._datasheets['main'].path
                else:
                    # Fallback to first available datasheet
                    first_name = list(self.input_datasheets._datasheets.keys())[0]
                    input_datasheet = self.input_datasheets._datasheets[first_name].path
            elif isinstance(self.input_datasheets, list):
                input_datasheet = self.input_datasheets[0]
            else:
                input_datasheet = str(self.input_datasheets)
        
        return f"""echo "Creating results datasheet and queries files"
# Create main datasheet with id, source_pdb, sequence and inherited columns
python {self.datasheet_py} {self.seqs_folder} {pipeline_name} {input_datasheet} {self.main_datasheet}

# Create queries CSV and FASTA from the main datasheet (needed for AlphaFold)
echo "Creating queries CSV and FASTA from results datasheet"
python {self.fa_to_csv_fasta_py} {self.seqs_folder} {self.queries_csv} {self.queries_fasta}

"""
    
    def _predict_sequence_ids(self) -> List[str]:
        """
        Predict the sequence IDs that ProteinMPNN will generate.
        
        Based on input structures and num_sequences parameter.
        Returns list of sequence IDs in the format: {pdb_base}_{seq_num}
        """
        sequence_ids = []
        
        # Use the original input parameter (set in constructor), not the processed one
        input_source = self.input if hasattr(self, 'input') else self.input_structures
        
        # Check for input data from upstream tools or direct file paths
        upstream_tool = None
        direct_file_paths = []
        
        # Case 1: ToolOutput input 
        if hasattr(input_source, 'get_output_files'):
            upstream_tool = input_source
        # Case 2: Direct file paths (from StandardizedOutput)
        elif isinstance(input_source, list):
            direct_file_paths = input_source
        
        if upstream_tool:
            # Get input PDB files from upstream tool
            input_pdbs = []
            
            # Get input PDB files from upstream tool
            input_pdbs = upstream_tool.get_output_files("pdbs")
            if not input_pdbs:
                input_pdbs = upstream_tool.get_output_files("structures")
            if not input_pdbs:
                raise ValueError(f"No PDB/structure files found in upstream tool {upstream_tool.tool_type}")
            
            # Process the PDB files if we got any
            if input_pdbs:
                for pdb_path in input_pdbs:
                    pdb_base = os.path.splitext(os.path.basename(pdb_path))[0]
                    # ProteinMPNN generates sequences numbered from 1
                    for seq_num in range(1, self.num_sequences + 1):
                        sequence_ids.append(f"{pdb_base}_{seq_num}")
        
        elif direct_file_paths:
            # Handle direct file paths from StandardizedOutput (input=rfd)
            for pdb_path in direct_file_paths:
                pdb_base = os.path.splitext(os.path.basename(pdb_path))[0]
                # ProteinMPNN generates sequences numbered from 1
                for seq_num in range(1, self.num_sequences + 1):
                    sequence_ids.append(f"{pdb_base}_{seq_num}")
        
        elif hasattr(self, 'input_sources') and self.input_sources:
            # Direct PDB file inputs
            for pdb_path in self.input_sources.values():
                pdb_base = os.path.splitext(os.path.basename(pdb_path))[0]
                for seq_num in range(1, self.num_sequences + 1):
                    sequence_ids.append(f"{pdb_base}_{seq_num}")
        
        # Must have sequence IDs from input sources
        if not sequence_ids:
            raise ValueError("Could not determine sequence IDs - no valid input structures found")
        
        return sequence_ids

    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after ProteinMPNN execution.
        
        Uses pure path construction - no filesystem access.
        Returns expected paths based on ProteinMPNN output patterns.
        
        Returns:
            Dictionary mapping output types to file paths with standard keys:
            - structures: Empty (no structures from ProteinMPNN)
            - compounds: Empty (no compounds from ProteinMPNN)
            - sequences: FASTA files
            - datasheets: Main datasheet CSV
            - output_folder: Tool's output directory
        """
        # Ensure file paths are set up
        if not hasattr(self, 'seqs_folder') or self.seqs_folder is None:
            # Fallback if configure_inputs hasn't been called yet
            self._setup_file_paths()
        
        seqs_folder = self.seqs_folder
        main_datasheet = self.main_datasheet
        queries_csv = self.queries_csv
        queries_fasta = self.queries_fasta
        
        # Expected FASTA files - ProteinMPNN generates one per input structure
        fasta_files = []
        
        # If we have dependency info, use it to predict filenames
        if self.input_is_tool_output:
            # Get expected input PDBs from dependency
            dependency_outputs = self.input_structures.get_output_files("pdbs")
            for pdb_path in dependency_outputs:
                pdb_base = os.path.splitext(os.path.basename(pdb_path))[0]
                fasta_path = os.path.join(seqs_folder, f"{pdb_base}.fa")
                fasta_files.append(fasta_path)
        else:
            # For direct file inputs, generate based on expected pattern
            # Default to pipeline-based naming
            pipeline_name = self._extract_pipeline_name()
            fasta_files = [os.path.join(seqs_folder, f"{pipeline_name}_sequences.fa")]
        
        # Predict sequence IDs for downstream tools
        sequence_ids = self._predict_sequence_ids()
        
        # Organize datasheets by content type with detailed metadata
        # Import DatasheetInfo
        from .base_config import DatasheetInfo

        datasheets = {
            "sequences": DatasheetInfo(
                name="sequences",
                path=main_datasheet,
                columns=["id", "source_id", "source_pdb", "sequence", "score", "seq_recovery", "rmsd"],
                description="ProteinMPNN sequence generation results with scores and structure recovery metrics",
                count=len(sequence_ids)  # Number of expected sequences
            )
        }
        
        return {
            "structures": [],
            "structure_ids": [],
            "compounds": [],
            "compound_ids": [],
            "sequences": [queries_csv],  # Main output is the datasheet with sequences and scores
            "sequence_ids": sequence_ids,
            "datasheets": datasheets,
            "output_folder": self.output_folder,
            # Keep legacy aliases for compatibility
            "main": main_datasheet,  # Legacy alias for backward compatibility
            "fa_files": fasta_files,  # Individual .fa files
            "queries_csv": [queries_csv],
            "queries_fasta": [queries_fasta],
            "seqs_folder": [seqs_folder]
        }
    
    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including ProteinMPNN-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "mpnn_params": {
                "num_sequences": self.num_sequences,
                "fixed_positions": self.fixed_positions,
                "designed_positions": self.designed_positions,
                "fixed_chain": self.fixed_chain,
                "plddt_threshold": self.plddt_threshold,
                "sampling_temp": self.sampling_temp,
                "model_name": self.model_name,
                "soluble_model": self.soluble_model,
                "input_type": "tool_output" if self.input_is_tool_output else "file/folder"
            }
        })
        return base_dict