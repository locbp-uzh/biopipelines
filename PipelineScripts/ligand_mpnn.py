"""
LigandMPNN configuration for ligand-aware sequence design.

Standardized to match ProteinMPNN structure with runtime fixed/redesigned selection
and ligand specification.
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


class LigandMPNN(BaseConfig):
    """
    LigandMPNN configuration for ligand-aware sequence design.
    
    Standardized to match ProteinMPNN: takes structures as input, outputs sequences,
    with runtime determination of fixed/redesigned positions and ligand specification.
    """
    
    TOOL_NAME = "LigandMPNN"
    DEFAULT_ENV = None  # Loaded from config.yaml
    
    def __init__(self,
                 structures: Union[str, List[str], ToolOutput],
                 ligand: str = "",
                 datasheets: Optional[List[str]] = None,
                 name: str = "",
                 num_sequences: int = 1,
                 fixed: str = "",
                 redesigned: str = "",
                 design_within: float = 5.0,
                 model: str = "v_32_010",
                 batch_size: int = 1,
                 **kwargs):
        """
        Initialize LigandMPNN configuration.

        Args:
            structures: Input structures (PDB files or ToolOutput from previous tool)
            ligand: Ligand identifier for binding site focus
            datasheets: Input datasheet files for metadata (from previous tool)
            name: Job name for output files
            num_sequences: Number of sequences to generate (batch_size * num_batches)
            fixed: Fixed positions (LigandMPNN format: "A3 A4 A5") or datasheet reference (e.g., (tool.datasheets.structures, "fixed"))
            redesigned: Designed positions (LigandMPNN format: "A3 A4 A5") or datasheet reference (e.g., (distance_analysis.datasheets.selections, "within"))  
            design_within: Distance in Å from ligand to redesign (fallback if positions not specified)
            model: LigandMPNN model version to use
            batch_size: Batch size for processing
            **kwargs: Additional parameters
        """
        # Store input parameters
        self.input_structures = structures
        self.input_datasheets = datasheets or {}
        self.input_is_tool_output = isinstance(structures, ToolOutput)
        
        # Store LigandMPNN-specific parameters
        self.ligand = ligand
        self.name = name or kwargs.get('job_name', '')
        self.num_sequences = num_sequences
        self.fixed_positions = fixed
        self.designed_positions = redesigned
        self.design_within = design_within
        self.model = model
        self.batch_size = batch_size
        
        # Calculate num_batches from num_sequences and batch_size
        self.num_batches = max(1, (num_sequences + batch_size - 1) // batch_size)
        
        # Initialize base class
        super().__init__(**kwargs)

        # Initialize file paths (will be set in configure_inputs)
        self._initialize_file_paths()

    def validate_params(self):
        """Validate LigandMPNN-specific parameters."""
        if not self.input_structures:
            raise ValueError("input structures parameter is required")
        
        if not self.ligand:
            raise ValueError("ligand parameter is required")
        
        if self.num_sequences <= 0:
            raise ValueError("num_sequences must be positive")
        
        if self.batch_size <= 0:
            raise ValueError("batch_size must be positive")
        
        if self.design_within <= 0:
            raise ValueError("design_within must be positive")
        
        # Validate datasheet references if provided
        if self.fixed_positions:
            self.validate_datasheet_reference(self.fixed_positions)
        if self.designed_positions:
            self.validate_datasheet_reference(self.designed_positions)
        
        # Validate model name
        valid_models = ["v_32_005", "v_32_010", "v_32_020", "v_32_025"]
        if self.model not in valid_models:
            raise ValueError(f"model must be one of: {valid_models}")
    
    
    def _initialize_file_paths(self):
        """Initialize common file paths used throughout the class."""
        self.seqs_folder = None
        self.queries_csv = None
        self.queries_fasta = None

        # Helper script paths
        self.fa_to_csv_fasta_py = None
        self.lmpnn_folder = None
        self.runtime_positions_py = None

        # Generated script paths
        self.commands_file = None
        self.replacement_script = None
    
    def _setup_file_paths(self):
        """Set up all file paths after output_folder is known."""
        # Extract pipeline/job name for file naming
        job_base = self.name or self._extract_job_name()

        # Core input/output files
        self.seqs_folder = os.path.join(self.output_folder, "seqs")
        self.queries_csv = os.path.join(self.output_folder, f"{job_base}_queries.csv")
        self.queries_fasta = os.path.join(self.output_folder, f"{job_base}_queries.fasta")

        # Generated script paths - defined once, used throughout
        self.commands_file = os.path.join(self.output_folder, "lmpnn_commands.sh")
        self.replacement_script = os.path.join(self.output_folder, "lmpnn_positions_replacement.sh")

        # Helper script paths (only set if folders are available)
        if hasattr(self, 'folders') and self.folders:
            self.fa_to_csv_fasta_py = os.path.join(self.folders["HelpScripts"], "pipe_fa_to_csv_fasta.py")
            self.lmpnn_folder = os.path.join(self.folders["data"], "LigandMPNN")
            self.runtime_positions_py = os.path.join(self.folders["HelpScripts"], "pipe_lmpnn_runtime_positions.py")
        else:
            # Temporary placeholders when folders aren't available yet
            self.fa_to_csv_fasta_py = None
            self.lmpnn_folder = None
            self.runtime_positions_py = None
    
    def _extract_job_name(self) -> str:
        """Extract job name from output folder structure."""
        # Get job name from parent folder
        # Structure: .../JobName_NNN/N_LigandMPNN -> JobName_NNN
        folder_parts = self.output_folder.split(os.sep)
        for i, part in enumerate(folder_parts):
            if "_LigandMPNN" in part:
                if i == 0:
                    raise ValueError(f"Invalid output folder structure: {self.output_folder}")
                return folder_parts[i-1]
        
        raise ValueError(f"Could not extract job name from output folder: {self.output_folder}")
    

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input structures from various sources."""
        self.folders = pipeline_folders
        self._setup_file_paths()  # Set up all file paths now that we have folders
        
        if self.input_is_tool_output:
            # Input from previous tool (e.g., RFdiffusion, Boltz2)
            tool_output: ToolOutput = self.input_structures
            
            # Try to get structures - could be in various formats
            source_structures = []
            
            # Try different output types
            for struct_type in ["structures", "pdbs"]:
                struct_files = tool_output.get_output_files(struct_type)
                if struct_files:
                    source_structures = struct_files
                    break
            
            if not source_structures:
                raise ValueError(f"No structure outputs found from {tool_output.tool_type}")
            
            # Store source for script generation
            self.input_sources["structures"] = source_structures
            
            # Add dependency
            self.dependencies.append(tool_output.config)
            
        elif isinstance(self.input_structures, list):
            # Direct list of structure file paths (from StandardizedOutput)
            if self.input_structures:
                self.input_sources["structures"] = self.input_structures
            else:
                raise ValueError("Empty structure list provided")
                
        elif isinstance(self.input_structures, StandardizedOutput):
            # StandardizedOutput object (from tool)
            if self.input_structures.structures:
                self.input_sources["structures"] = self.input_structures.structures
                # Also store structure IDs for proper tracking
                if self.input_structures.structure_ids:
                    self.input_sources["structure_ids"] = self.input_structures.structure_ids
                # Store datasheets for position references
                self.input_datasheets = self.input_structures.datasheets
                self.standardized_input = self.input_structures
            else:
                raise ValueError("No structures found in StandardizedOutput")
                
        elif isinstance(self.input_structures, str):
            # String input - single PDB file
            if self.input_structures.endswith('.pdb'):
                pdb_source = os.path.join(pipeline_folders["PDBs"], self.input_structures)
                if os.path.exists(pdb_source):
                    self.input_sources["structures"] = [pdb_source]
                else:
                    raise ValueError(f"PDB file not found: {pdb_source}")
            else:
                raise ValueError("String input must be a PDB file path")
        else:
            raise ValueError(f"Unsupported input type: {type(self.input_structures)}")
    
    def get_config_display(self) -> List[str]:
        """Get LigandMPNN configuration display lines."""
        config_lines = super().get_config_display()
        
        # Input information - show only structure count, not full details
        if self.input_is_tool_output:
            # Count structures from tool output
            structure_count = 0
            if hasattr(self.input_structures, 'structures') and self.input_structures.structures:
                structure_count = len(self.input_structures.structures)
            config_lines.append(f"INPUT: {self.input_structures.tool_type} output ({structure_count} structures)")
        else:
            config_lines.append(f"INPUT: {self.input_structures}")
        
        config_lines.extend([
            f"LIGAND: {self.ligand}",
            f"FIXED: {self.fixed_positions or 'Auto (from datasheet or ligand-based)'}",
            f"DESIGNED: {self.designed_positions or 'Auto (from datasheet or ligand-based)'}",
            f"DESIGN WITHIN: {self.design_within}Å",
            f"NUM SEQUENCES: {self.num_sequences}",
            f"BATCH SIZE: {self.batch_size}",
            f"MODEL: {self.model}"
        ])
        
        return config_lines

    def generate_script(self, script_path: str) -> str:
        """
        Generate LigandMPNN execution script.
        
        Args:
            script_path: Path where script should be written
            
        Returns:
            Script content as string
        """
        # Generate script content following modular pattern
        script_content = "#!/bin/bash\n"
        script_content += "# LigandMPNN execution script\n"
        script_content += "# Generated by ProteinNotebooks pipeline system\n\n"
        script_content += self.generate_completion_check_header()
        script_content += self.generate_script_setup_positions()
        script_content += self.generate_script_run_ligandmpnn()
        script_content += self.generate_script_convert_outputs()
        script_content += self.generate_completion_check_footer()
        
        return script_content
    
    def generate_script_setup_positions(self) -> str:
        """Generate the position setup and script modification part."""
        # Get specific structure files instead of directory
        if "structures" in self.input_sources:
            structure_files = self.input_sources["structures"]
            # Convert to comma-separated string for passing to script
            structure_files_str = ",".join(structure_files)
        else:
            raise ValueError("No structure sources found")
        
        # Determine fixed positions approach - prioritize explicit positions over datasheets
        if self.fixed_positions is not None or self.designed_positions is not None:
            # Use directly specified positions (highest priority) - includes empty strings
            input_source = "selection"  
            input_datasheet = "-"
        elif (self.input_is_tool_output and hasattr(self, 'input_datasheets') and self.input_datasheets) or \
             (hasattr(self, 'standardized_input') and self.standardized_input and hasattr(self.standardized_input, 'datasheets')):
            # Use datasheet from previous tool (e.g., RFdiffusion)
            input_source = "datasheet"
            if hasattr(self.standardized_input, 'datasheets') and self.standardized_input.datasheets:
                # StandardizedOutput format
                if hasattr(self.standardized_input.datasheets, '_datasheets') and 'structures' in self.standardized_input.datasheets._datasheets:
                    input_datasheet = self.standardized_input.datasheets._datasheets['structures'].path
                else:
                    input_datasheet = "-"
            elif isinstance(self.input_datasheets, dict):
                # Dictionary format
                if "structures" in self.input_datasheets:
                    input_datasheet = self.input_datasheets["structures"]["path"]
                elif "main" in self.input_datasheets:
                    input_datasheet = self.input_datasheets["main"]["path"]
                else:
                    first_key = next(iter(self.input_datasheets))
                    input_datasheet = self.input_datasheets[first_key]["path"]
            elif hasattr(self.input_datasheets, '_datasheets'):
                # DatasheetContainer format
                if 'structures' in self.input_datasheets._datasheets:
                    input_datasheet = self.input_datasheets._datasheets['structures'].path
                else:
                    first_name = list(self.input_datasheets._datasheets.keys())[0]
                    input_datasheet = self.input_datasheets._datasheets[first_name].path
            else:
                input_datasheet = "-"
        else:
            # Use ligand-based design with design_within cutoff
            input_source = "ligand"
            input_datasheet = "-"
        
        # Resolve datasheet references in fixed/designed positions
        resolved_fixed = self.resolve_datasheet_reference(self.fixed_positions) if self.fixed_positions else "-"
        resolved_designed = self.resolve_datasheet_reference(self.designed_positions) if self.designed_positions else "-"

        return f"""echo "Setting up LigandMPNN position constraints"
# Create a separate commands file that will be modified instead of this script
cat > {self.commands_file} << 'EOF'
#!/bin/bash
cd {self.lmpnn_folder}

# LigandMPNN commands with placeholders - will be replaced by position script
"""
    
    def generate_script_run_ligandmpnn(self) -> str:
        """Generate the LigandMPNN execution part of the script."""
        # Get specific structure files instead of directory
        if "structures" in self.input_sources:
            structure_files = self.input_sources["structures"]
            # Convert to comma-separated string for passing to script
            structure_files_str = ",".join(structure_files)
        else:
            raise ValueError("No structure sources found")
            
        # Determine fixed positions approach - prioritize explicit positions over datasheets
        if self.fixed_positions is not None or self.designed_positions is not None:
            # Use directly specified positions (highest priority) - includes empty strings
            input_source = "selection"  
            input_datasheet = "-"
        elif (self.input_is_tool_output and hasattr(self, 'input_datasheets') and self.input_datasheets) or \
             (hasattr(self, 'standardized_input') and self.standardized_input and hasattr(self.standardized_input, 'datasheets')):
            # Use datasheet from previous tool (e.g., RFdiffusion)
            input_source = "datasheet"
            if hasattr(self.standardized_input, 'datasheets') and self.standardized_input.datasheets:
                # StandardizedOutput format
                if hasattr(self.standardized_input.datasheets, '_datasheets') and 'structures' in self.standardized_input.datasheets._datasheets:
                    input_datasheet = self.standardized_input.datasheets._datasheets['structures'].path
                else:
                    input_datasheet = "-"
            elif isinstance(self.input_datasheets, dict):
                # Dictionary format
                if "structures" in self.input_datasheets:
                    input_datasheet = self.input_datasheets["structures"]["path"]
                elif "main" in self.input_datasheets:
                    input_datasheet = self.input_datasheets["main"]["path"]
                else:
                    first_key = next(iter(self.input_datasheets))
                    input_datasheet = self.input_datasheets[first_key]["path"]
            elif hasattr(self.input_datasheets, '_datasheets'):
                # DatasheetContainer format
                if 'structures' in self.input_datasheets._datasheets:
                    input_datasheet = self.input_datasheets._datasheets['structures'].path
                else:
                    first_name = list(self.input_datasheets._datasheets.keys())[0]
                    input_datasheet = self.input_datasheets._datasheets[first_name].path
            else:
                input_datasheet = "-"
        else:
            # Use ligand-based design with design_within cutoff
            input_source = "ligand"
            input_datasheet = "-"
        
        # Resolve datasheet references in fixed/designed positions
        resolved_fixed = self.resolve_datasheet_reference(self.fixed_positions) if self.fixed_positions else "-"
        resolved_designed = self.resolve_datasheet_reference(self.designed_positions) if self.designed_positions else "-"
        
        # Build base LigandMPNN options
        base_options = f'--model_type "ligand_mpnn"'
        base_options += f' --checkpoint_ligand_mpnn "./model_params/ligandmpnn_{self.model}_25.pt"'
        base_options += f' --batch_size {self.batch_size}'
        base_options += f' --number_of_batches {self.num_batches}'
        base_options += f' --ligand_mpnn_cutoff_for_score "{self.design_within}"'
        base_options += f' --out_folder "{self.output_folder}"'

        # Generate commands that will be written to the separate file
        commands = []
        
        # Use the input structure list to generate commands
        if "structures" in self.input_sources:
            structure_list = self.input_sources["structures"]
            for pdb_file in sorted(structure_list):  # Sort for consistent ordering
                pdb_name = os.path.basename(pdb_file)
                pdb_id = os.path.splitext(pdb_name)[0]
                
                commands.append(f"echo \"Processing {pdb_name} with positions for ID: {pdb_id}\"")
                commands.append(f"python run.py {base_options} --pdb_path \"{pdb_file}\" {pdb_id}_FIXED_OPTION_PLACEHOLDER {pdb_id}_REDESIGNED_OPTION_PLACEHOLDER")
                commands.append("")
        
        # Create the script content that writes commands and then executes them
        return f"""# Write LigandMPNN commands to separate file
{chr(10).join(commands)}
EOF

# Make the commands file executable
chmod +x {self.commands_file}

# Use existing HelpScript to create position replacement script
echo "Creating position replacement script..."
python {self.runtime_positions_py} "{structure_files_str}" "{input_source}" "{input_datasheet}" "{resolved_fixed}" "{resolved_designed}" "{self.ligand}" "{self.design_within}" "{self.replacement_script}"

# Run the replacement script on the commands file (not this script)
echo "Running position replacement script on commands file: {self.commands_file}"
bash {self.replacement_script} {self.commands_file}

# Now execute the modified commands file
echo "Executing LigandMPNN commands..."
bash {self.commands_file}

"""
    
    def generate_script_convert_outputs(self) -> str:
        """Generate the output conversion part of the script."""
        return f"""echo "Converting FASTA outputs to CSV format"
python {self.fa_to_csv_fasta_py} {self.seqs_folder} {self.queries_csv} {self.queries_fasta} --duplicates

"""
    
    def _predict_sequence_outputs(self) -> List[str]:
        """
        Predict sequence files that LigandMPNN will generate.
        
        Based on input structures and LigandMPNN output patterns.
        """
        sequence_files = []
        
        # LigandMPNN generates FASTA files in seqs folder
        # Use seqs folder as the main output (actual files determined at runtime)
        sequence_files = [self.seqs_folder]
        
        return sequence_files
    
    def _predict_sequence_ids(self) -> List[str]:
        """
        Predict the sequence IDs that LigandMPNN will generate.
        
        Based on input structures and num_sequences parameter.
        Returns list of sequence IDs in the format: {pdb_base}_{seq_num}
        """
        sequence_ids = []
        
        # Use structures parameter directly
        input_source = self.input_structures
        
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
                    # LigandMPNN generates sequences numbered from 1
                    for seq_num in range(1, self.num_sequences + 1):
                        sequence_ids.append(f"{pdb_base}_{seq_num}")
        
        elif direct_file_paths:
            # Handle direct file paths from StandardizedOutput (input=tool)
            for pdb_path in direct_file_paths:
                pdb_base = os.path.splitext(os.path.basename(pdb_path))[0]
                # LigandMPNN generates sequences numbered from 1
                for seq_num in range(1, self.num_sequences + 1):
                    sequence_ids.append(f"{pdb_base}_{seq_num}")
        
        elif hasattr(self, 'input_sources') and self.input_sources and "structures" in self.input_sources:
            # Direct PDB file inputs
            for pdb_path in self.input_sources["structures"]:
                pdb_base = os.path.splitext(os.path.basename(pdb_path))[0]
                for seq_num in range(1, self.num_sequences + 1):
                    sequence_ids.append(f"{pdb_base}_{seq_num}")
        
        # Must have sequence IDs from input sources
        if not sequence_ids:
            raise ValueError("Could not determine sequence IDs - no valid input structures found")
        
        return sequence_ids
    
    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after LigandMPNN execution.
        
        Uses pure path construction - no filesystem access.
        Returns expected paths based on LigandMPNN output patterns.
        
        Returns:
            Dictionary mapping output types to file paths with standard keys:
            - structures: Empty (no structures from LigandMPNN)
            - compounds: Empty (no compounds from LigandMPNN) 
            - sequences: FASTA files and converted CSV
            - datasheets: Main datasheet CSV
            - output_folder: Tool's output directory
        """
        # Ensure file paths are set up
        if not hasattr(self, 'seqs_folder') or self.seqs_folder is None:
            # Fallback if configure_inputs hasn't been called yet
            self._setup_file_paths()
        
        seqs_folder = self.seqs_folder
        queries_csv = self.queries_csv
        queries_fasta = self.queries_fasta
        
        # Predict sequence outputs
        sequence_files = self._predict_sequence_outputs()
        
        # Predict sequence IDs for downstream tools
        sequence_ids = self._predict_sequence_ids()
        
        # Organize datasheets by content type with detailed metadata
        datasheets = {
            "sequences": DatasheetInfo(
                name="sequences",
                path=queries_csv,
                columns=["id", "sequence", "sample", "T", "seed", "overall_confidence", "ligand_confidence", "seq_rec"],
                description="LigandMPNN ligand-aware sequence generation results with binding scores",
                count=len(sequence_ids)
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
            "main": queries_csv,  # Legacy alias for backward compatibility
            "fa_files": sequence_files,  # Individual .fa files
            "queries_csv": [queries_csv],
            "queries_fasta": [queries_fasta],
            "seqs_folder": [seqs_folder]
        }
    
    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including LigandMPNN-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "lmpnn_params": {
                "ligand": self.ligand,
                "name": self.name,
                "num_sequences": self.num_sequences,
                "fixed_positions": self.fixed_positions,
                "designed_positions": self.designed_positions,
                "design_within": self.design_within,
                "model": self.model,
                "batch_size": self.batch_size,
                "num_batches": self.num_batches,
                "input_type": "tool_output" if self.input_is_tool_output else "direct"
            }
        })
        return base_dict