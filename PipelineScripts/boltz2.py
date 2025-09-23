"""
Boltz2 configuration for protein-ligand complex prediction.

Handles apo and holo structure prediction with MSA caching,
ligand binding affinity calculation, and comprehensive analysis.
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


class Boltz2(BaseConfig):
    """
    Boltz2 configuration for protein-ligand complex prediction.
    
    Predicts both apo (protein-only) and holo (protein-ligand) structures
    with automatic MSA management and comprehensive scoring.
    """
    
    # Tool identification
    TOOL_NAME = "Boltz2"
    DEFAULT_ENV = "Boltz2Env"
    DEFAULT_RESOURCES = {"gpu": "V100", "memory": "16GB", "time": "24:00:00"}
    
    def __init__(self, 
                 # Standardized input parameter (like other pipeline tools)
                 input: Union[str, List[str], ToolOutput, Dict[str, Any]] = None,
                 # Primary input parameters (matching usage examples)
                 config: Optional[str] = None,
                 proteins: Union[str, List[str], ToolOutput] = None,
                 ligands: Union[str, ToolOutput, StandardizedOutput, None] = None,
                 msas: Optional[Union[str, ToolOutput]] = None,
                 # Legacy compatibility
                 sequences: Union[str, List[str], ToolOutput] = None,
                 ligand_smiles: Optional[str] = None,
                 # Library-based ligand inputs
                 ligand_library: Optional[str] = None,
                 primary_key: Optional[str] = None,
                 library_repr: str = "SMILES",
                 library_type: str = "noncovalent",
                 # Core prediction parameters
                 affinity: bool = True,
                 output_format: str = "pdb",
                 msa_server: str = "public",
                 global_msas_cache: bool = False,
                 # Advanced prediction parameters
                 recycling_steps: Optional[int] = None,
                 diffusion_samples: Optional[int] = None,
                 use_potentials: bool = False,
                 **kwargs):
        """
        Initialize Boltz2 configuration.
        
        Args:
            input: Complete standardized input with sequences, datasheets, etc. (e.g., lmpnn)
            config: Direct YAML configuration string (Example 1: Boltz2(config=yaml))
            proteins: Protein sequences - can be ToolOutput, file path, or direct sequence (Examples 2&3)
            ligands: Single ligand SMILES string, ToolOutput with compounds, or datasheet reference
            msas: MSA files for recycling (e.g., boltz2_apo.datasheets.msas)
            sequences: Legacy parameter, same as proteins (for backward compatibility)
            ligand_smiles: Legacy parameter, same as ligands (for backward compatibility)
            ligand_library: Path to CSV file with ligand library
            primary_key: Key column in library to filter by
            library_repr: Ligand representation ("SMILES" or "CCD")
            library_type: Binding type ("noncovalent" or "covalent")
            affinity: Whether to calculate binding affinity
            output_format: Output format ("pdb" or "mmcif")
            msa_server: MSA generation ("public" or "local")
            global_msas_cache: Enable global MSA caching to reuse MSAs across jobs
            recycling_steps: Number of recycling steps (default: Boltz2 default 3)
            diffusion_samples: Number of diffusion samples (default: Boltz2 default 1)
            use_potentials: Enable potentials for improved structure prediction
            **kwargs: Additional parameters
        """
        # Initialize default values
        self.config = config
        self.proteins = proteins or sequences  
        self.ligands = ligands or ligand_smiles
        self.msas = msas
        self.input_sequences = None
        self.input_compounds = []
        self.input_datasheets = {}
        self.input_is_tool_output = False
        self.standardized_input = None

        # Handle explicit parameter inputs from previous tools
        if isinstance(proteins, StandardizedOutput):
            # Explicit: proteins=tool
            self.input_sequences = proteins.sequences
            self.input_datasheets = getattr(proteins, 'datasheets', {})
            self.input_is_tool_output = False
            self.standardized_input = proteins
            # Keep proteins as StandardizedOutput reference
        elif isinstance(ligands, (StandardizedOutput, ToolOutput)):
            # Explicit: ligands=compounds from previous tool
            if isinstance(ligands, StandardizedOutput):
                self.input_compounds = getattr(ligands, 'compounds', [])
                self.input_datasheets = getattr(ligands, 'datasheets', {})
            else:  # ToolOutput
                self.input_compounds = ligands.get_output_files("compounds")
                self.input_datasheets = ligands.get_output_files("datasheets")
                # Add dependency for ToolOutput
                self.dependencies.append(ligands.config)
            # Keep ligands as tool reference
        elif isinstance(msas, StandardizedOutput):
            # Explicit: msas=previous_boltz
            self.input_datasheets = getattr(msas, 'datasheets', {})

        # Handle legacy standardized input format (less preferred)
        elif input is not None:
            # Legacy format: input=lmpnn (ambiguous)
            if isinstance(input, StandardizedOutput):
                # StandardizedOutput object - try to infer what to use
                self.input_sequences = getattr(input, 'sequences', [])
                self.input_compounds = getattr(input, 'compounds', [])
                self.input_datasheets = getattr(input, 'datasheets', {})
                self.input_is_tool_output = False  # Direct file paths now
                self.standardized_input = input  # Keep reference
                
                # Auto-assign based on what's available (proteins take priority)
                if self.input_sequences and not proteins:
                    self.proteins = self.input_sequences
                if self.input_compounds and not ligands and not ligand_smiles:
                    # Auto-use compounds as ligands if no explicit ligands provided
                    self.ligands = self.input_compounds[0] if self.input_compounds else None
            elif isinstance(input, ToolOutput):
                # Direct ToolOutput object
                self.input_sequences = input
                self.input_compounds = input.get_output_files("compounds")
                self.input_datasheets = input.get_output_files("datasheets")
                self.input_is_tool_output = True
                self.standardized_input = None
            elif isinstance(input, dict):
                # Dictionary format with standardized keys
                self.input_sequences = input.get('sequences', [])
                self.input_compounds = input.get('compounds', [])
                self.input_datasheets = input.get('datasheets', {})
                self.input_is_tool_output = False  # Direct file paths
                self.standardized_input = None
            else:
                # Fallback to treating as protein sequences
                self.input_sequences = input
                self.input_compounds = []
                self.input_datasheets = {}
                self.input_is_tool_output = isinstance(input, ToolOutput)
                self.standardized_input = None
            
            # When using input parameter, these should be None unless explicitly overridden
            self.config = config
            self.proteins = proteins or self.input_sequences
            self.ligands = ligands or ligand_smiles
            self.msas = msas
        else:
            # Legacy format: individual parameters
            self.config = config
            self.proteins = proteins or sequences  # Use proteins if provided, fallback to sequences
            self.ligands = ligands or ligand_smiles  # Use ligands if provided, fallback to ligand_smiles
            self.msas = msas
            
            # No standardized input
            self.input_sequences = self.proteins
            self.input_compounds = []
            self.input_datasheets = {}
            self.input_is_tool_output = isinstance(self.proteins, ToolOutput)
            self.standardized_input = None
        
        # Store other Boltz2-specific parameters
        self.ligand_library = ligand_library
        self.primary_key = primary_key
        self.library_repr = library_repr
        self.library_type = library_type
        # Override affinity to False if no ligands are present
        if self.ligands is None and affinity:
            print("Warning: No ligands detected, setting affinity=False")
            self.affinity = False
        else:
            self.affinity = affinity  # Renamed from calculate_affinity
        self.output_format = output_format
        self.msa_server = msa_server
        self.global_msas_cache = global_msas_cache
        self.recycling_steps = recycling_steps
        self.diffusion_samples = diffusion_samples
        self.use_potentials = use_potentials
        
        # Track input source type and files (tool-agnostic)
        self.input_yaml_entities = None
        self.input_fasta_files = []
        self.queries_csv_file = None
        
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
        affinity_datasheet = DatasheetInfo(
            name="affinity",
            path="",  # Path will be set when output_folder is known
            columns=["id", "affinity_pred_value", "affinity_confidence", "structure_confidence"],
            description="Boltz2 affinity predictions with confidence scores"
        )

        # Set up datasheets container for IDE autocompletion
        datasheets_dict = {"affinity": affinity_datasheet}
        self.datasheets = DatasheetContainer(datasheets_dict)

        # CRITICAL: Explicitly set datasheet attributes for IDE autocompletion
        self.datasheets.affinity = affinity_datasheet

    def _initialize_file_paths(self):
        """Initialize common file paths used throughout the class."""
        # Folder structure
        self.apo_config_folder = None
        self.bound_config_folder = None
        self.library_folder = None
        self.apo_prediction_folder = None
        self.bound_prediction_folder = None
        self.msa_cache_folder = None
        
        # Configuration files
        self.base_config_file = None
        self.bound_config_file = None
        self.queries_csv = None
        self.queries_fasta = None
        self.expanded_library_csv = None
        
        # Output files
        self.library_scores_csv = None
        self.config_txt_file = None
        self.results_zip = None
        
        # Helper script paths
        self.smiles_library_py = None
        self.boltz_config_py = None
        self.boltz_results_py = None
        self.mmseqs2_client_sh = None
        self.fa_to_csv_py = None
        
        # MSA settings
        self.msa_extension = None
        self.msa_option = None
    
    def _setup_file_paths(self):
        """Set up all file paths after output_folder is known."""
        # Folder structure in tool's output folder
        self.apo_config_folder = os.path.join(self.output_folder, "ApoConfig")
        self.bound_config_folder = os.path.join(self.output_folder, "BoundConfig")
        self.library_folder = os.path.join(self.output_folder, "Library")
        self.apo_prediction_folder = os.path.join(self.output_folder, "ApoPredictions")
        self.bound_prediction_folder = os.path.join(self.output_folder, "BoundPredictions")
        
        # MSA management - create MSAs folder in the tool directory
        self.msa_cache_folder = os.path.join(self.output_folder, "MSAs")
        self.msa_extension = ".a3m" if self.msa_server == "local" else ".csv"
        self.msa_option = "" if self.msa_server == "local" else " --use_msa_server"
        
        # Configuration files
        self.base_config_file = os.path.join(self.output_folder, "base_config.txt")
        self.bound_config_file = os.path.join(self.bound_config_folder, "bound_config.yaml")
        self.queries_csv = os.path.join(self.output_folder, f"{self.job_name}_queries.csv")
        self.queries_fasta = os.path.join(self.output_folder, f"{self.job_name}_queries.fasta")
        self.expanded_library_csv = os.path.join(self.output_folder, "expanded_smiles_library.csv")
        
        # Output files
        self.library_scores_csv = os.path.join(self.library_folder, "library_scores.csv")
        self.config_txt_file = os.path.join(self.output_folder, f"{self.job_name}_config.txt")
        self.results_zip = os.path.join(self.output_folder, f"{self.job_name}.zip")
        
        # Helper script paths (only set if folders are available)
        if hasattr(self, 'folders') and self.folders:
            self.smiles_library_py = os.path.join(self.folders["HelpScripts"], "pipe_smiles_library.py")
            self.boltz_config_py = os.path.join(self.folders["HelpScripts"], "pipe_boltz_config.py")
            self.boltz_protein_ligand_configs_py = os.path.join(self.folders["HelpScripts"], "pipe_boltz_protein_ligand_configs.py")
            self.boltz_results_py = os.path.join(self.folders["HelpScripts"], "pipe_boltz_results.py")
            self.mmseqs2_client_sh = os.path.join(self.folders.get("MMseqs2", "MMseqs2"), "mmseqs2_client.sh")
            self.fa_to_csv_py = os.path.join(self.folders["HelpScripts"], "pipe_fa_to_csv_fasta.py")
        else:
            # Temporary placeholders when folders aren't available yet
            self.smiles_library_py = None
            self.boltz_config_py = None
            self.boltz_protein_ligand_configs_py = None
            self.boltz_results_py = None
            self.mmseqs2_client_sh = None
            self.fa_to_csv_py = None

    def validate_params(self):
        """Validate Boltz2-specific parameters."""
        # Must have some form of input
        has_input = any([
            hasattr(self, 'input') and getattr(self, 'input', None) is not None,
            getattr(self, 'config', None) is not None,
            getattr(self, 'proteins', None) is not None,
            getattr(self, 'input_sequences', None) is not None
        ])
        if not has_input:
            raise ValueError("Either input, config, or proteins parameter is required")
        
        # Cannot specify multiple primary input methods (only check if they exist)
        primary_inputs = []
        if getattr(self, 'config', None) is not None:
            primary_inputs.append(self.config)
        if getattr(self, 'proteins', None) is not None:
            primary_inputs.append(self.proteins)
        
        if len(primary_inputs) > 1:
            raise ValueError("Cannot specify multiple primary input methods (config and proteins)")
        
        # For proteins input (not config), ligand information is typically required
        # But with standardized input, ligands might come from datasheets
        if (self.proteins and not self.config and 
            not self.ligands and not self.ligand_library and not self.msas and
            not self.standardized_input):
            # Only warn, don't error - ligands might be in datasheets or MSAs
            pass
        
        # Cannot specify both single ligand and library
        if self.ligands and self.ligand_library:
            raise ValueError("Cannot specify both ligands and ligand_library")
        
        # Validate datasheet references if provided (skip for StandardizedOutput)
        if self.ligands and isinstance(self.ligands, str):
            self.validate_datasheet_reference(self.ligands)
        
        # Validate enum values
        if self.library_repr not in ["SMILES", "CCD"]:
            raise ValueError("library_repr must be 'SMILES' or 'CCD'")
        
        if self.library_type not in ["noncovalent", "covalent"]:
            raise ValueError("library_type must be 'noncovalent' or 'covalent'")
        
        if self.output_format not in ["pdb", "mmcif"]:
            raise ValueError("output_format must be 'pdb' or 'mmcif'")
        
        if self.msa_server not in ["public", "local"]:
            raise ValueError("msa_server must be 'public' or 'local'")
        
        # Validate recycling_steps and diffusion_samples
        if self.recycling_steps is not None and (not isinstance(self.recycling_steps, int) or self.recycling_steps < 1):
            raise ValueError("recycling_steps must be a positive integer")
        
        if self.diffusion_samples is not None and (not isinstance(self.diffusion_samples, int) or self.diffusion_samples < 1):
            raise ValueError("diffusion_samples must be a positive integer")
    
    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input sequences from various sources (tool-agnostic)."""
        self.folders = pipeline_folders
        self._setup_file_paths()  # Set up all file paths now that we have folders
        
        # Determine primary input source (prioritize config, then proteins, then input_sequences)
        primary_input = self.config or self.proteins or self.input_sequences
        
        if self.input_is_tool_output:
            # Input from previous tool (any upstream tool via ToolOutput)
            tool_output: ToolOutput = primary_input
            
            # Look for queries.csv file or any compatible sequence format
            queries_csv = tool_output.get_output_files("queries_csv")
            if queries_csv:
                self.queries_csv_file = queries_csv[0]
            else:
                # Look for FASTA files or sequence files
                fasta_files = tool_output.get_output_files("sequences")
                if fasta_files:
                    self.input_fasta_files = fasta_files
                else:
                    raise ValueError(f"No sequence outputs found from upstream tool")
            
            # Add dependency
            self.dependencies.append(tool_output.config)
            
        elif self.config:
            # Direct YAML configuration (Example 1: config=yaml)
            if isinstance(self.config, str):
                if self.config.strip().startswith('sequences:') or '- protein:' in self.config:
                    # Direct YAML entities specification
                    self.input_yaml_entities = self.config
                else:
                    # File path to YAML config
                    config_path = os.path.join(pipeline_folders.get("data", "."), self.config)
                    if os.path.exists(config_path):
                        with open(config_path, 'r') as f:
                            self.input_yaml_entities = f.read()
                    else:
                        raise ValueError(f"Config file not found: {config_path}")
            else:
                raise ValueError(f"Invalid config type: {type(self.config)}")
        
        elif isinstance(primary_input, str):
            # String input - could be file path or direct sequence
            if primary_input.endswith('.csv'):
                # CSV queries file
                queries_path = os.path.join(pipeline_folders.get("data", "."), primary_input)
                if os.path.exists(queries_path):
                    self.queries_csv_file = queries_path
                else:
                    raise ValueError(f"Queries CSV file not found: {queries_path}")
            elif primary_input.endswith('.fasta') or primary_input.endswith('.fa'):
                # FASTA file
                fasta_path = os.path.join(pipeline_folders.get("data", "."), primary_input)
                if os.path.exists(fasta_path):
                    self.input_fasta_files = [fasta_path]
                else:
                    raise ValueError(f"FASTA file not found: {fasta_path}")
            else:
                # Direct protein sequence - will create YAML config from it
                self.input_direct_sequence = primary_input
        
        elif isinstance(primary_input, list):
            # List of FASTA files
            self.input_fasta_files = []
            for fasta_file in primary_input:
                fasta_path = os.path.join(pipeline_folders.get("data", "."), fasta_file)
                if os.path.exists(fasta_path):
                    self.input_fasta_files.append(fasta_path)
                else:
                    raise ValueError(f"FASTA file not found: {fasta_path}")
        
        elif hasattr(primary_input, '__class__') and 'StandardizedOutput' in str(type(primary_input)):
            # StandardizedOutput object (e.g., from lmpnn)
            if hasattr(primary_input, 'sequences') and primary_input.sequences:
                # Use the first sequence file from StandardizedOutput
                self.queries_csv_file = primary_input.sequences[0]
            else:
                raise ValueError("No sequences found in StandardizedOutput")
        
        else:
            raise ValueError(f"Invalid input type: {type(primary_input)}. Supported types: str, list, ToolOutput, StandardizedOutput")
        
        # Handle compound library input from various sources
        compounds_file = None
        
        # Priority 1: Explicit ligands parameter with StandardizedOutput
        if isinstance(getattr(self, 'ligands', None), StandardizedOutput):
            ligands_output = self.ligands
            if hasattr(ligands_output, 'compounds') and ligands_output.compounds:
                compounds_file = ligands_output.compounds[0]
                print(f"Using compounds from explicit ligands parameter: {compounds_file}")
        
        # Priority 2: input_compounds from legacy input parameter
        elif hasattr(self, 'input_compounds') and self.input_compounds:
            compounds_file = self.input_compounds[0]
            print(f"Using compounds from legacy input parameter: {compounds_file}")
        
        # Priority 3: ligands as direct file path
        elif isinstance(self.ligands, str) and self.ligands.endswith('.csv'):
            compounds_file = self.ligands
            print(f"Using compounds from ligands file path: {compounds_file}")
            
        # Set as ligand library if found
        if compounds_file:
            self.ligand_library = compounds_file
        
        # Handle existing ligand_library parameter (CSV file path)
        elif self.ligand_library:
            # Check if the ligand library file exists
            if not os.path.exists(self.ligand_library):
                # Try in project directory
                project_path = os.path.join(pipeline_folders["notebooks"], self.ligand_library)
                if os.path.exists(project_path):
                    self.ligand_library = project_path
                else:
                    raise ValueError(f"Ligand library file not found: {self.ligand_library}")
        
        # Handle MSA inputs - build YAML config with MSA paths
        if self.msas:
            self._integrate_msa_inputs()
        
        # Handle global MSA cache lookup
        if self.global_msas_cache:
            self._integrate_global_msa_cache()
    
    def _integrate_msa_inputs(self):
        """Build YAML config that includes MSA paths from previous tool."""
        if not self.msas:
            return
            
        # Get MSA datasheet from previous tool
        msa_paths = {}
        if hasattr(self.msas, 'datasheets'):
            if hasattr(self.msas, 'datasheets') and 'msas' in self.msas.datasheets:
                # MSAs should be provided in the datasheets.msas
                # Handle both DatasheetContainer (returns path string) and dict format
                msa_datasheet = self.msas.datasheets['msas']
                if isinstance(msa_datasheet, str):
                    # DatasheetContainer returns path directly
                    msa_datasheet_path = msa_datasheet
                elif isinstance(msa_datasheet, dict) and 'path' in msa_datasheet:
                    # Raw dict format
                    msa_datasheet_path = msa_datasheet['path']
                else:
                    raise ValueError(f"Invalid MSA datasheet format: {type(msa_datasheet)}")
                
                # This will be a CSV with columns: id, sequence_id, msa_file
                # We need to read this at runtime, not pipeline time
                # For now, just signal that MSAs will be integrated at runtime
                pass
            else:
                raise ValueError("MSA ToolOutput doesn't have datasheets.msas")
        elif isinstance(self.msas, str):
            # Direct path to MSA folder
            pass
        
        # Override config generation to include MSA paths in YAML
        # This ensures MSAs are specified in config rather than calculated
        self._build_yaml_config_with_msas()
    
    def _build_yaml_config_with_msas(self):
        """Build YAML configuration that specifies MSA paths instead of calculating them."""
        
        # For now, we'll modify the script generation to handle MSAs properly
        # The key insight is that we need to modify the YAML config at runtime
        # to include msa: path/to/file.csv entries for each protein
        
        # Set a flag that we need MSA integration
        self._needs_msa_integration = True
        
        # Clear any existing input yaml entities since we need to rebuild with MSAs
        if hasattr(self, 'input_yaml_entities'):
            self._original_yaml_entities = self.input_yaml_entities
            self.input_yaml_entities = None
    
    def _get_protein_sequence(self) -> str:
        """Get protein sequence for MSA config building."""
        if hasattr(self, 'input_direct_sequence') and self.input_direct_sequence:
            return self.input_direct_sequence
        elif hasattr(self, '_original_yaml_entities') and self._original_yaml_entities:
            # Extract sequence from original YAML - this is a simplified approach
            # In practice, we'd need a more robust YAML parser
            lines = self._original_yaml_entities.split('\n')
            for line in lines:
                if 'sequence:' in line:
                    return line.split('sequence:')[1].strip()
        return ""
    
    def _get_ligand_smiles(self) -> str:
        """Get ligand SMILES for MSA config building."""
        # Return the stored ligands value if it's a SMILES string
        if self.ligands and isinstance(self.ligands, str) and not self.ligands.endswith('.csv'):
            return self.ligands
        
        # If ligands is from YAML entities, parse it
        if hasattr(self, '_original_yaml_entities') and self._original_yaml_entities:
            lines = self._original_yaml_entities.split('\n')
            for line in lines:
                if 'smiles:' in line:
                    return line.split('smiles:')[1].strip()
        return ""
    
    def _generate_global_msa_cache_section(self, boltz_cache_folder: str) -> str:
        """
        Generate script section for global MSA cache lookup.
        
        Args:
            boltz_cache_folder: Path to Boltz cache folder
            
        Returns:
            Bash script content for global MSA cache handling
        """
        if not hasattr(self, '_cache_protein_sequence') or not self._cache_protein_sequence:
            return ""
        
        # Script to check global cache and potentially use cached MSA
        effective_job_name = self.get_effective_job_name()
        msa_filename = effective_job_name if effective_job_name is not None else "msa"
        
        return f'''
echo "Checking global MSA cache for protein sequence"
GLOBAL_MSA_CACHE_SCRIPT="{os.path.join(self.folders['HelpScripts'], 'pipe_manage_global_msa_cache.py')}"
CACHED_MSA_PATH="{os.path.join(self.msa_cache_folder, f'{msa_filename}.csv')}"

# Try to find MSA in global cache
python "$GLOBAL_MSA_CACHE_SCRIPT" lookup \\
    --cache-folder "{boltz_cache_folder}" \\
    --sequence "{self._cache_protein_sequence}" \\
    --output "$CACHED_MSA_PATH"

if [ $? -eq 0 ]; then
    echo "Found MSA in global cache, will use cached version"
    GLOBAL_MSA_FOUND=true
else
    echo "MSA not found in global cache, will calculate new MSA"
    GLOBAL_MSA_FOUND=false
fi

'''
    
    def _generate_global_msa_cache_storage_section(self, boltz_cache_folder: str) -> str:
        """
        Generate script section for storing newly calculated MSAs to global cache.
        
        This version extracts sequences from MSA files and caches them if not already present.
        
        Args:
            boltz_cache_folder: Path to Boltz cache folder
            
        Returns:
            Bash script content for MSA storage
        """
        # Script to extract sequences from MSA files and update global cache
        return f'''
# Update global MSA cache from organized MSAs folder
if [ "$GLOBAL_MSA_FOUND" = false ]; then
    echo "Extracting sequences from MSA files and updating global cache"
    
    MSA_FOLDER="{os.path.join(self.output_folder, 'MSAs')}"
    
    if [ -d "$MSA_FOLDER" ]; then
        python "{os.path.join(self.folders['HelpScripts'], 'pipe_update_global_msa_cache_from_folder.py')}" \\
            "$MSA_FOLDER" \\
            "{boltz_cache_folder}"
        
        if [ $? -eq 0 ]; then
            echo "Successfully updated global MSA cache from MSA folder"
        else
            echo "Warning: Failed to update global MSA cache"
        fi
    else
        echo "Warning: MSAs folder not found at $MSA_FOLDER"
    fi
else
    echo "Used cached MSA, skipping global cache update"
fi

'''
    
    def _integrate_global_msa_cache(self):
        """Set up global MSA cache integration."""
        # Mark that we need global MSA cache integration
        self._needs_global_msa_cache = True
        
        # Extract protein sequence for cache lookup
        self._cache_protein_sequence = self._extract_protein_sequence()
        
        if not self._cache_protein_sequence:
            print("Warning: Could not extract protein sequence for global MSA cache")
            self._needs_global_msa_cache = False
    
    def _extract_protein_sequence(self) -> str:
        """Extract protein sequence for global MSA cache lookup."""
        # Try different sources for protein sequence
        
        # First check self.proteins (most common case for direct sequence input)
        if self.proteins and isinstance(self.proteins, str):
            return self.proteins.strip()
        
        # Then check other sources
        if hasattr(self, 'input_direct_sequence') and self.input_direct_sequence:
            return self.input_direct_sequence.strip()
        
        elif hasattr(self, 'input_yaml_entities') and self.input_yaml_entities:
            # Extract from YAML entities
            lines = self.input_yaml_entities.split('\n')
            for line in lines:
                if 'sequence:' in line:
                    seq = line.split('sequence:')[1].strip()
                    if seq and not seq.startswith('id:'):
                        return seq
        
        elif hasattr(self, '_original_yaml_entities') and self._original_yaml_entities:
            # Extract from original YAML entities  
            lines = self._original_yaml_entities.split('\n')
            for line in lines:
                if 'sequence:' in line:
                    seq = line.split('sequence:')[1].strip()
                    if seq and not seq.startswith('id:'):
                        return seq
        
        return ""
    
    def generate_script(self, script_path: str) -> str:
        """
        Generate bash script for Boltz2 execution.
        
        Args:
            script_path: Path where script should be written
            
        Returns:
            Script content as string
        """
        # Boltz2 cache and configuration
        boltz_cache_folder = self.folders["BoltzCache"]
        
        # MSA configuration
        msa_option = "" if self.msa_server == "local" else " --use_msa_server"
        
        # Build script header
        script_content = "#!/bin/bash\n"
        script_content += "# Boltz2 execution script\n"
        script_content += "# Generated by ProteinNotebooks pipeline system\n\n"
        script_content += self.generate_completion_check_header()
        
        # Create basic folder structure
        script_content += f"""# Create output folders
mkdir -p {os.path.join(self.output_folder, "predictions")}
mkdir -p {self.msa_cache_folder}

"""
        
        # Global MSA cache lookup section
        if hasattr(self, '_needs_global_msa_cache') and self._needs_global_msa_cache:
            script_content += self._generate_global_msa_cache_section(boltz_cache_folder)
        
        # Handle MSA recycling if provided
        if self.msas:
            script_content += self._generate_msa_recycling_section()
        
        # Check if we need multiple config files execution (define early for use in conditions)
        has_multiple_ligands = False
        if self.ligands:
            if hasattr(self.ligands, 'compound_ids') and self.ligands.compound_ids:
                has_multiple_ligands = len(self.ligands.compound_ids) > 1
            elif hasattr(self.ligands, 'compounds') and self.ligands.compounds:
                has_multiple_ligands = len(self.ligands.compounds) > 1
        
        # Generate input configuration based on input type
        effective_job_name = self.get_effective_job_name()
        if effective_job_name is None:
            config_name = "prediction"
            #print("Tip: Pass pipeline=pipeline to Boltz2() to use meaningful pipeline job name.")
        else:
            config_name = effective_job_name
        config_file_path = os.path.join(self.output_folder, f"{config_name}.yaml")
        
        if self.config:
            # Direct YAML configuration (Example 1)
            script_content += f"""
echo "Using direct YAML configuration"
cat > {config_file_path} << 'EOF'
{self.config}
EOF

"""
        elif self.input_yaml_entities:
            # YAML entities from configure_inputs
            script_content += f"""
echo "Creating configuration from input entities"
cat > {config_file_path} << 'EOF'
{self.input_yaml_entities}
EOF

"""
        elif hasattr(self, 'input_direct_sequence') and self.input_direct_sequence and (has_multiple_ligands or hasattr(self, '_needs_global_msa_cache') or hasattr(self, '_needs_msa_integration')):
            # Direct protein sequence with multiple ligands OR MSA integration - use unified approach
            # Handle ligand parameter extraction
            ligand_param = "None"
            if self.ligands:
                if isinstance(self.ligands, StandardizedOutput):
                    # Extract compounds CSV path from StandardizedOutput
                    if hasattr(self.ligands, 'compounds') and self.ligands.compounds:
                        ligand_param = self.ligands.compounds[0]
                elif isinstance(self.ligands, ToolOutput):
                    # Extract compounds CSV path from ToolOutput
                    compounds_files = self.ligands.get_output_files("compounds")
                    if compounds_files:
                        ligand_param = compounds_files[0]
                elif isinstance(self.ligands, str):
                    # Direct file path or datasheet reference
                    ligand_param = self.resolve_datasheet_reference(self.ligands)
            affinity_flag = "--affinity" if self.affinity else ""
            config_files_dir = os.path.join(self.output_folder, "config_files")
            proteins_csv = os.path.join(self.output_folder, "proteins.csv")
            
            # Handle ligand parameter - create CSV if it's a direct string
            if ligand_param != "None":
                if ligand_param.endswith('.csv'):
                    # It's already a CSV file path
                    final_ligand_param = ligand_param
                else:
                    # It's a direct ligand string - create CSV
                    ligands_csv = os.path.join(self.output_folder, "ligands.csv")
                    final_ligand_param = ligands_csv
            else:
                final_ligand_param = "None"
            
            script_content += f"""
echo "Creating configurations for all protein-ligand combinations"
mkdir -p {config_files_dir}

# Create protein CSV from direct sequence
cat > {proteins_csv} << 'EOF'
id,sequence
protein,{self.input_direct_sequence}
EOF
"""
            
            # Add ligand CSV creation if needed
            if ligand_param != "None" and not ligand_param.endswith('.csv'):
                ligand_id = config_name if config_name != "prediction" else "ligand"
                script_content += f"""
# Create ligand CSV from direct ligand string
cat > {os.path.join(self.output_folder, "ligands.csv")} << 'EOF'
id,format,smiles,ccd
{ligand_id},smiles,{ligand_param},
EOF
"""
            
            # Prepare MSA parameter if available
            msa_param = ""
            if hasattr(self, '_needs_msa_integration') and self._needs_msa_integration:
                if hasattr(self.msas, 'datasheets') and 'msas' in self.msas.datasheets:
                    msa_ds = self.msas.datasheets['msas']
                    if isinstance(msa_ds, str):
                        msa_datasheet = msa_ds
                    elif isinstance(msa_ds, dict) and 'path' in msa_ds:
                        msa_datasheet = msa_ds['path']
                    else:
                        msa_datasheet = str(msa_ds)
                    msa_param = f'--msa-datasheet "{msa_datasheet}"'
            elif hasattr(self, '_needs_global_msa_cache') and self._needs_global_msa_cache:
                msa_param = f'--global-msa-cache --msa-folder "{self.msa_cache_folder}"'
            
            # Add job name parameter
            job_name_param = f'--job-name "{config_name}"' if config_name != "prediction" else ""
            
            script_content += f"""
# Generate config files for all protein-ligand combinations
python {self.boltz_protein_ligand_configs_py} {proteins_csv} "{final_ligand_param}" {config_files_dir} {affinity_flag} {msa_param} {job_name_param}

"""
        elif hasattr(self, '_needs_global_msa_cache') and self._needs_global_msa_cache:
            # Generate YAML config with global MSA cache integration (legacy single case)
            script_content += f"""
echo "Creating configuration with global MSA cache integration"
python {os.path.join(self.folders['HelpScripts'], 'pipe_build_boltz_config_with_global_cache.py')} \\
    --protein-sequence "{self._cache_protein_sequence}" \\
    --ligand-smiles "{self._get_ligand_smiles()}" \\
    --cached-msa-path "$CACHED_MSA_PATH" \\
    --output "{config_file_path}" \\
    --affinity "{self.affinity if self.affinity else 'false'}"

"""
        elif hasattr(self, '_needs_msa_integration') and self._needs_msa_integration and not self.queries_csv_file:
            # Generate YAML config with MSAs at runtime (single protein case)
            msa_datasheet = ""
            if hasattr(self.msas, 'datasheets') and 'msas' in self.msas.datasheets:
                # Handle both DatasheetContainer and dict formats
                msa_ds = self.msas.datasheets['msas']
                if isinstance(msa_ds, str):
                    msa_datasheet = msa_ds
                elif isinstance(msa_ds, dict) and 'path' in msa_ds:
                    msa_datasheet = msa_ds['path']
                else:
                    msa_datasheet = str(msa_ds)  # Fallback
            
            script_content += f"""
echo "Creating configuration with MSA integration"
python {os.path.join(self.folders['HelpScripts'], 'pipe_build_boltz_config_with_msas.py')} \\
    --protein-sequence "{self._get_protein_sequence()}" \\
    --ligand-smiles "{self._get_ligand_smiles()}" \\
    --msa-datasheet "{msa_datasheet}" \\
    --output "{config_file_path}" \\
    --affinity "{self.affinity if self.affinity else 'false'}"

"""
        elif self.queries_csv_file:
            # Convert CSV to separate YAML configs for each sequence
            config_files_dir = os.path.join(self.output_folder, "config_files")
            # Handle ligand parameter extraction
            ligand_param = "None"
            if self.ligands:
                if isinstance(self.ligands, StandardizedOutput):
                    # Extract compounds CSV path from StandardizedOutput
                    if hasattr(self.ligands, 'compounds') and self.ligands.compounds:
                        ligand_param = self.ligands.compounds[0]
                elif isinstance(self.ligands, ToolOutput):
                    # Extract compounds CSV path from ToolOutput
                    compounds_files = self.ligands.get_output_files("compounds")
                    if compounds_files:
                        ligand_param = compounds_files[0]
                elif isinstance(self.ligands, str):
                    # Direct file path or datasheet reference
                    ligand_param = self.resolve_datasheet_reference(self.ligands)
            affinity_flag = "--affinity" if self.affinity else ""
            
            # Add MSA parameter if MSAs are provided  
            msa_datasheet_flag = ""
            if self.msas:
                if hasattr(self.msas, 'datasheets'):
                    # ToolOutput or StandardizedOutput case - access MSA datasheet
                    if hasattr(self.msas.datasheets, '_datasheets') and 'msas' in self.msas.datasheets._datasheets:
                        msa_datasheet_path = self.msas.datasheets._datasheets['msas'].path
                        msa_datasheet_flag = f'--msa-datasheet "{msa_datasheet_path}"'
                    elif hasattr(self.msas.datasheets, 'msas'):
                        msa_datasheet_path = self.msas.datasheets.msas  # Direct path access
                        msa_datasheet_flag = f'--msa-datasheet "{msa_datasheet_path}"'
                elif isinstance(self.msas, str) and self.msas.endswith('.csv'):
                    # Direct CSV path
                    msa_datasheet_flag = f'--msa-datasheet "{self.msas}"'
            
            script_content += f"""
echo "Converting CSV to YAML configuration"
mkdir -p {config_files_dir}
python {self.boltz_protein_ligand_configs_py} {self.queries_csv_file} "{ligand_param}" {config_files_dir} {affinity_flag} {msa_datasheet_flag}

"""
        
        elif self.input_fasta_files:
            # Convert FASTA files to separate YAML configs
            config_files_dir = os.path.join(self.output_folder, "config_files")
            fasta_files_str = ",".join(self.input_fasta_files)
            # Handle ligand parameter extraction
            ligand_param = "None"
            if self.ligands:
                if isinstance(self.ligands, StandardizedOutput):
                    # Extract compounds CSV path from StandardizedOutput
                    if hasattr(self.ligands, 'compounds') and self.ligands.compounds:
                        ligand_param = self.ligands.compounds[0]
                elif isinstance(self.ligands, ToolOutput):
                    # Extract compounds CSV path from ToolOutput
                    compounds_files = self.ligands.get_output_files("compounds")
                    if compounds_files:
                        ligand_param = compounds_files[0]
                elif isinstance(self.ligands, str):
                    # Direct file path or datasheet reference
                    ligand_param = self.resolve_datasheet_reference(self.ligands)
            affinity_flag = "--affinity" if self.affinity else ""
            
            # Add MSA parameter if MSAs are provided  
            msa_datasheet_flag = ""
            if self.msas:
                if hasattr(self.msas, 'datasheets'):
                    # ToolOutput or StandardizedOutput case - access MSA datasheet
                    if hasattr(self.msas.datasheets, '_datasheets') and 'msas' in self.msas.datasheets._datasheets:
                        msa_datasheet_path = self.msas.datasheets._datasheets['msas'].path
                        msa_datasheet_flag = f'--msa-datasheet "{msa_datasheet_path}"'
                    elif hasattr(self.msas.datasheets, 'msas'):
                        msa_datasheet_path = self.msas.datasheets.msas  # Direct path access
                        msa_datasheet_flag = f'--msa-datasheet "{msa_datasheet_path}"'
                elif isinstance(self.msas, str) and self.msas.endswith('.csv'):
                    # Direct CSV path
                    msa_datasheet_flag = f'--msa-datasheet "{self.msas}"'
            
            script_content += f"""
echo "Converting FASTA files to YAML configuration"
mkdir -p {config_files_dir}
python {self.boltz_protein_ligand_configs_py} "{fasta_files_str}" "{ligand_param}" {config_files_dir} {affinity_flag} {msa_datasheet_flag}

"""
            base_config_file = self.queries_csv
        
        elif hasattr(self, 'input_direct_sequence') and self.input_direct_sequence:
            # Direct protein sequence - check if we need multiple configs for multiple ligands
            # Handle ligand parameter extraction
            ligand_param = "None"
            if self.ligands:
                if isinstance(self.ligands, StandardizedOutput):
                    # Extract compounds CSV path from StandardizedOutput
                    if hasattr(self.ligands, 'compounds') and self.ligands.compounds:
                        ligand_param = self.ligands.compounds[0]
                elif isinstance(self.ligands, ToolOutput):
                    # Extract compounds CSV path from ToolOutput
                    compounds_files = self.ligands.get_output_files("compounds")
                    if compounds_files:
                        ligand_param = compounds_files[0]
                elif isinstance(self.ligands, str):
                    # Direct file path or datasheet reference
                    ligand_param = self.resolve_datasheet_reference(self.ligands)
            affinity_flag = "--affinity" if self.affinity else ""
            
            # Check if we have multiple ligands
            has_multiple_ligands = False
            if self.ligands and hasattr(self.ligands, 'compound_ids') and self.ligands.compound_ids:
                has_multiple_ligands = len(self.ligands.compound_ids) > 1
            
            # Always generate separate config files for each protein-ligand combination
            config_files_dir = os.path.join(self.output_folder, "config_files")
            
            # Create protein CSV from the direct sequence  
            proteins_csv = os.path.join(self.output_folder, "proteins.csv")
            
            # Handle ligand parameter - create CSV if it's a direct string
            if ligand_param != "None":
                if ligand_param.endswith('.csv'):
                    # It's already a CSV file path
                    final_ligand_param = ligand_param
                else:
                    # It's a direct ligand string - create CSV
                    ligands_csv = os.path.join(self.output_folder, "ligands.csv")
                    final_ligand_param = ligands_csv
            else:
                final_ligand_param = "None"
            
            script_content += f"""
echo "Creating configurations for all protein-ligand combinations"
mkdir -p {config_files_dir}

# Create protein CSV from direct sequence
cat > {proteins_csv} << 'EOF'
id,sequence
protein,{self.input_direct_sequence}
EOF
"""
            
            # Add ligand CSV creation if needed
            if ligand_param != "None" and not ligand_param.endswith('.csv'):
                ligand_id = config_name if config_name != "prediction" else "ligand"
                script_content += f"""
# Create ligand CSV from direct ligand string
cat > {os.path.join(self.output_folder, "ligands.csv")} << 'EOF'
id,format,smiles,ccd
{ligand_id},smiles,{ligand_param},
EOF
"""
            
            # Prepare MSA parameter if available
            msa_param = ""
            if hasattr(self, '_needs_msa_integration') and self._needs_msa_integration:
                if hasattr(self.msas, 'datasheets') and 'msas' in self.msas.datasheets:
                    msa_ds = self.msas.datasheets['msas']
                    if isinstance(msa_ds, str):
                        msa_datasheet = msa_ds
                    elif isinstance(msa_ds, dict) and 'path' in msa_ds:
                        msa_datasheet = msa_ds['path']
                    else:
                        msa_datasheet = str(msa_ds)
                    msa_param = f'--msa-datasheet "{msa_datasheet}"'
            elif hasattr(self, '_needs_global_msa_cache') and self._needs_global_msa_cache:
                msa_param = f'--global-msa-cache --msa-folder "{self.msa_cache_folder}"'
            
            # Add job name parameter
            job_name_param = f'--job-name "{config_name}"' if config_name != "prediction" else ""
            
            script_content += f"""
# Generate config files for all protein-ligand combinations
python {self.boltz_protein_ligand_configs_py} {proteins_csv} "{final_ligand_param}" {config_files_dir} {affinity_flag} {msa_param} {job_name_param}

"""
        
        # Run Boltz2 prediction
        boltz_options = f"--cache {boltz_cache_folder} --out_dir {self.output_folder}{msa_option} --output_format {self.output_format}"
        
        # Add recycling_steps and diffusion_samples if specified
        if self.recycling_steps is not None:
            boltz_options += f" --recycling_steps {self.recycling_steps}"
        
        if self.diffusion_samples is not None:
            boltz_options += f" --diffusion_samples {self.diffusion_samples}"

        if self.use_potentials:
            boltz_options += " --use_potentials"
        
        if self.queries_csv_file or self.input_fasta_files or (hasattr(self, 'input_direct_sequence') and self.input_direct_sequence and (has_multiple_ligands or hasattr(self, '_needs_global_msa_cache') or hasattr(self, '_needs_msa_integration'))):
            # Multiple config files - run prediction for each one
            config_files_dir = os.path.join(self.output_folder, "config_files")
            script_content += f"""
echo "Running Boltz2 prediction on individual config files"
for config_file in {config_files_dir}/*.yaml; do
    if [ -f "$config_file" ]; then
        echo "Processing config: $config_file"
        boltz predict "$config_file" {boltz_options}
    fi
done

"""
        else:
            # Single config file
            script_content += f"""
echo "Running Boltz2 prediction"
boltz predict {config_file_path} {boltz_options}

"""
        
        # Post-process results using dedicated script
        script_content += self._generate_postprocess_with_script()
        
        # Add newly calculated MSAs to global cache after post-processing
        if hasattr(self, '_needs_global_msa_cache') and self._needs_global_msa_cache:
            script_content += self._generate_global_msa_cache_storage_section(boltz_cache_folder)
        
        script_content += self.generate_completion_check_footer()
        
        return script_content
    
    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after Boltz2 execution.
        
        Based on Boltz2 documentation, the output structure is:
        out_dir/predictions/[input_file]/
        ├── [input_file]_model_0.cif                    # Best structure
        ├── confidence_[input_file]_model_0.json        # Confidence scores
        ├── affinity_[input_file].json                  # Affinity scores (if enabled)
        └── [input_file]_model_[N].cif                  # Additional samples
        
        Returns:
            Dictionary mapping output types to file paths
        """
        # Predict input file names (will be config names or generated)
        input_file_names = self._predict_input_file_names()
        
        # Structure files (CIF or PDB based on output_format)
        structures = []
        structure_ids = []
        
        # Predict sequence IDs from input (these will be the final structure names after post-processing)
        predicted_sequence_ids = self._predict_sequence_ids()
        
        # After post-processing, structures will be named by sequence IDs in output folder
        for seq_id in predicted_sequence_ids:
            if self.output_format == "pdb":
                structure_file = os.path.join(self.output_folder, f"{seq_id}.pdb")
            else:  # mmcif
                structure_file = os.path.join(self.output_folder, f"{seq_id}.cif")
            
            structures.append(structure_file)
            structure_ids.append(seq_id)
        
        # Use the predicted sequence IDs from above
        sequence_ids = predicted_sequence_ids
        
        # MSA files for recycling (stored in tool MSAs folder)
        msa_files = []
        # Only predict MSA files if msa_cache_folder is set (after configure_inputs is called)
        if hasattr(self, 'msa_cache_folder') and self.msa_cache_folder:
            # MSA files will be created with config names and copied to proper names by post-processing
            for seq_id in sequence_ids:
                if self.msa_server == "local":
                    msa_file = os.path.join(self.msa_cache_folder, f"{seq_id}.a3m")
                else:
                    # Public server creates CSV MSA files that get renamed during post-processing
                    msa_file = os.path.join(self.msa_cache_folder, f"{seq_id}.csv")
                msa_files.append(msa_file)
        
        # Organize datasheets by content type with detailed metadata
        datasheets = {}
        
        # Confidence scores datasheet (converted from JSON to CSV for compatibility)
        confidence_csv = os.path.join(self.output_folder, "confidence_scores.csv")
        datasheets["confidence"] = DatasheetInfo(
            name="confidence",
            path=confidence_csv,
            columns=["id", "input_file", "confidence_score", "ptm", "iptm", "complex_plddt", "complex_iplddt"],
            description="Boltz2 confidence scores for structure predictions",
            count=len(input_file_names)
        )
        
        # Affinity scores datasheet (if affinity calculation enabled)
        if self.affinity:
            affinity_csv = os.path.join(self.output_folder, "affinity_scores.csv")
            datasheets["affinity"] = DatasheetInfo(
                name="affinity",
                path=affinity_csv,
                columns=["id", "input_file", "affinity_pred_value", "affinity_probability_binary"],
                description="Boltz2 binding affinity predictions",
                count=len(input_file_names)
            )
        
        # MSAs datasheet (for recycling between apo/holo predictions)
        if msa_files:
            msa_csv = os.path.join(self.output_folder, "msas.csv")
            datasheets["msas"] = DatasheetInfo(
                name="msas",
                path=msa_csv,
                columns=["id", "sequence_id", "sequence", "msa_file"],
                description="MSA files for sequence recycling between predictions",
                count=len(msa_files)
            )
        
        # Input sequences datasheet (for downstream tools)
        sequences_csv = os.path.join(self.output_folder, "sequences.csv")
        datasheets["sequences"] = DatasheetInfo(
            name="sequences",
            path=sequences_csv,
            columns=["id", "sequence"],
            description="Input protein sequences used for prediction",
            count=len(sequence_ids)
        )
        
        # Compounds handling
        compounds = []
        compound_ids = []

        # Handle compounds from various sources
        if hasattr(self, 'ligand_library') and self.ligand_library:
            compounds = [self.ligand_library]  # Point to the compound library CSV

            # Try to get compound_ids from standardized input
            if isinstance(getattr(self, 'ligands', None), StandardizedOutput):
                # Explicit ligands=compounds
                compound_ids = self.ligands.compound_ids
            elif hasattr(self, 'standardized_input') and self.standardized_input:
                # Legacy input=compounds
                compound_ids = getattr(self.standardized_input, 'compound_ids', ["library_compounds"])
            else:
                # Fallback placeholder
                compound_ids = ["library_compounds"]

            # Add compounds datasheet pointing to the ligand library
            datasheets["compounds"] = DatasheetInfo(
                name="compounds",
                path=self.ligand_library,
                columns=["id", "format", "smiles", "ccd"],
                description="Compound library used for ligand binding predictions",
                count=len(compound_ids) if isinstance(compound_ids, list) else 0
            )

        # Handle single SMILES string - create compounds table
        elif isinstance(self.ligands, str) and self.ligands and not self.ligands.endswith('.csv'):
            # Create a compounds table for the single SMILES string
            compounds_csv = os.path.join(self.output_folder, "ligands.csv")
            compounds = [compounds_csv]

            # Use same logic as script generation for ligand ID
            effective_job_name = self.get_effective_job_name()
            if effective_job_name is None:
                ligand_id = "ligand"
            else:
                ligand_id = effective_job_name
            compound_ids = [ligand_id]

            # Add compounds datasheet for single SMILES
            datasheets["compounds"] = DatasheetInfo(
                name="compounds",
                path=compounds_csv,
                columns=["id", "format", "smiles", "ccd"],
                description="Single ligand compound used for binding prediction",
                count=1
            )
        
        return {
            "structures": structures,
            "structure_ids": structure_ids,
            "compounds": compounds,
            "compound_ids": compound_ids,
            "sequences": [sequences_csv],  # Main sequence file
            "sequence_ids": sequence_ids,
            "msas": msa_files,  # Individual MSA files for recycling
            "datasheets": datasheets,
            "output_folder": self.output_folder,
            # Keep some legacy aliases for compatibility
            "pdbs": structures,
            "predictions_folder": [os.path.join(self.output_folder, "predictions")]
        }
    
    def _predict_input_file_names(self) -> List[str]:
        """
        Predict the input file names that Boltz2 will use for output folders.
        
        Returns:
            List of expected input file names
        """
        if self.config:
            # Direct YAML config - use job name or default
            return [self.job_name or "config"]
        elif hasattr(self, 'input_yaml_entities') and self.input_yaml_entities:
            # YAML entities - use job name or default  
            return [self.job_name or "sequences"]
        elif self.queries_csv_file:
            # CSV file - use basename without extension
            return [os.path.splitext(os.path.basename(self.queries_csv_file))[0]]
        elif self.input_fasta_files:
            # FASTA files - use basename without extension for each
            return [os.path.splitext(os.path.basename(f))[0] for f in self.input_fasta_files]
        else:
            # Default fallback
            return [self.job_name or "prediction"]
    
    def _predict_sequence_ids(self) -> List[str]:
        """
        Predict sequence IDs from input sources (tool-agnostic).
        
        Returns:
            List of expected sequence IDs
        """
        # Try to get from standardized input first (highest priority)
        if hasattr(self, 'standardized_input') and self.standardized_input:
            if hasattr(self.standardized_input, 'sequence_ids') and self.standardized_input.sequence_ids:
                return self.standardized_input.sequence_ids
        
        # Try to get from proteins parameter if it's StandardizedOutput
        if hasattr(self.proteins, 'sequence_ids') and self.proteins.sequence_ids:
            return self.proteins.sequence_ids
            
        # Try to get from input_sequences if it's StandardizedOutput  
        if hasattr(self.input_sequences, 'sequence_ids') and self.input_sequences.sequence_ids:
            return self.input_sequences.sequence_ids
        
        # Try to extract from dependencies
        for dep in self.dependencies:
            if hasattr(dep, '_predict_sequence_ids'):
                return dep._predict_sequence_ids()
        
        # Calculate actual expected structures based on proteins × compounds
        input_names = self._predict_input_file_names()
        
        # Count proteins
        if isinstance(self.proteins, str):
            # Single protein sequence
            num_proteins = 1
        elif isinstance(self.proteins, list):
            num_proteins = len(self.proteins)
        elif hasattr(self.proteins, 'sequence_ids'):
            num_proteins = len(self.proteins.sequence_ids)
        else:
            num_proteins = len(input_names) if input_names else 1
        
        # Count compounds - get from actual ligands input
        num_compounds = 1  # Default
        if hasattr(self, 'ligands') and self.ligands:
            if hasattr(self.ligands, 'compound_ids') and self.ligands.compound_ids:
                num_compounds = len(self.ligands.compound_ids)
            elif isinstance(self.ligands, str):
                # Direct SMILES string - single compound
                num_compounds = 1
            else:
                # Default to 1 for other cases (direct string, etc.)
                num_compounds = 1
        
        # Get actual ligand IDs if available
        ligand_ids = []
        if hasattr(self, 'ligands') and hasattr(self.ligands, 'compound_ids') and self.ligands.compound_ids:
            ligand_ids = self.ligands.compound_ids
        elif hasattr(self, 'standardized_input') and self.standardized_input and hasattr(self.standardized_input, 'compound_ids'):
            ligand_ids = self.standardized_input.compound_ids
        elif isinstance(self.ligands, str):
            # Direct SMILES string - use "ligand" as ID
            ligand_ids = ["ligand"]
        
        # Get actual protein IDs if available  
        protein_ids = []
        if hasattr(self, 'proteins') and hasattr(self.proteins, 'sequence_ids') and self.proteins.sequence_ids:
            protein_ids = self.proteins.sequence_ids
        
        # Generate IDs based on protein/ligand combinations
        if num_proteins == 1 and num_compounds > 1:
            # One protein, multiple ligands: use ligand IDs
            if ligand_ids:
                return ligand_ids[:num_compounds]
            else:
                return [f"ligand_{i:03d}" for i in range(1, num_compounds + 1)]
        elif num_proteins > 1 and num_compounds == 1:
            # Multiple proteins, one ligand: use protein IDs
            if protein_ids:
                return protein_ids[:num_proteins]
            else:
                return [f"protein_{i:03d}" for i in range(1, num_proteins + 1)]
        elif num_proteins > 1 and num_compounds > 1:
            # Multiple proteins, multiple ligands: use PROTEIN_LIGAND format
            sequence_ids = []
            for protein_idx in range(num_proteins):
                protein_id = protein_ids[protein_idx] if protein_ids and protein_idx < len(protein_ids) else f"protein_{protein_idx+1:03d}"
                for ligand_idx in range(num_compounds):
                    ligand_id = ligand_ids[ligand_idx] if ligand_ids and ligand_idx < len(ligand_ids) else f"ligand_{ligand_idx+1:03d}"
                    sequence_ids.append(f"{protein_id}_{ligand_id}")
            return sequence_ids
        else:
            # Single protein, single ligand: use effective job name or fallback
            if hasattr(self, '_needs_global_msa_cache') and self._needs_global_msa_cache:
                effective_job_name = self.get_effective_job_name()
                return [effective_job_name if effective_job_name is not None else "prediction"]
            base_name = input_names[0] if input_names else "structure"
            return [base_name]
    
    def _generate_msa_recycling_section(self) -> str:
        """
        Generate script section for MSA recycling.
        
        Returns:
            Bash script content for MSA handling
        """
        if hasattr(self.msas, 'datasheets'):
            # MSAs from previous Boltz2 prediction - use the MSA files directly
            if hasattr(self.msas, 'msas') and self.msas.msas:
                # Use the MSA files array from previous tool
                return f"""
echo "Recycling MSAs from previous prediction"
echo "Copying MSA files to cache folder"
mkdir -p {self.msa_cache_folder}
# Copy MSA files from previous Boltz2 output
cp {' '.join(self.msas.msas)} {self.msa_cache_folder}/

"""
            else:
                # Fallback: MSAs will be automatically found in shared cache folder
                return f"""
echo "Recycling MSAs from previous prediction"
# MSAs will be automatically found in shared cache folder
# Boltz2 will reuse existing MSAs if available

"""
        elif isinstance(self.msas, str):
            # Handle different string path types
            if self.msas.endswith('.csv'):
                # MSA datasheet CSV - read it and copy individual MSA files
                return f"""
echo "Using MSAs from datasheet: {self.msas}"
mkdir -p {self.msa_cache_folder}
# Read MSA datasheet and copy individual MSA files
python -c "
import pandas as pd
import shutil
import os

df = pd.read_csv('{self.msas}')
for _, row in df.iterrows():
    msa_file = row['msa_file']
    if os.path.exists(msa_file):
        dest_file = os.path.join('{self.msa_cache_folder}', os.path.basename(msa_file))
        shutil.copy2(msa_file, dest_file)
        print(f'Copied MSA: {{os.path.basename(msa_file)}}')
    else:
        print(f'Warning: MSA file not found: {{msa_file}}')
"

"""
            else:
                # Direct path to MSA folder
                msa_extension = ".a3m" if self.msa_server == "local" else ".csv"
                return f"""
echo "Using MSAs from directory: {self.msas}"
mkdir -p {self.msa_cache_folder}
cp {self.msas}/*{msa_extension} {self.msa_cache_folder}/

"""
        else:
            return ""
    
    
    def _generate_postprocess_section(self) -> str:
        """
        Generate script section for post-processing Boltz2 outputs.
        
        Returns:
            Bash script content for converting outputs to pipeline format
        """
        postprocess_script = f"""
echo "Post-processing Boltz2 outputs for pipeline compatibility"

# Create post-processing Python script
cat > {self.output_folder}/postprocess.py << 'PYEOF'
import json
import pandas as pd
import os
import glob

confidence_data = []
affinity_data = []
sequence_data = []

# Process all prediction folders
for pred_folder in glob.glob('{self.output_folder}/predictions/*'):
    input_name = os.path.basename(pred_folder)
    
    # Process confidence files
    for conf_file in glob.glob(f'{{pred_folder}}/confidence_*.json'):
        with open(conf_file, 'r') as f:
            conf_data = json.load(f)
        
        conf_row = {{
            'id': f'{{input_name}}_model_0',
            'input_file': input_name,
            'confidence_score': conf_data.get('confidence_score', 0),
            'ptm': conf_data.get('ptm', 0),
            'iptm': conf_data.get('iptm', 0),
            'complex_plddt': conf_data.get('complex_plddt', 0),
            'complex_iplddt': conf_data.get('complex_iplddt', 0)
        }}
        confidence_data.append(conf_row)
"""
        
        if self.affinity:
            postprocess_script += """
    
    # Process affinity files
    for aff_file in glob.glob(f'{pred_folder}/affinity_*.json'):
        with open(aff_file, 'r') as f:
            aff_data = json.load(f)
        
        aff_row = {
            'id': f'{input_name}_affinity',
            'input_file': input_name,
            'affinity_pred_value': aff_data.get('affinity_pred_value', 0),
            'affinity_probability_binary': aff_data.get('affinity_probability_binary', 0)
        }
        affinity_data.append(aff_row)
"""
        
        postprocess_script += f"""

# Save CSV files
if confidence_data:
    pd.DataFrame(confidence_data).to_csv('{self.output_folder}/confidence_scores.csv', index=False)
"""
        
        if self.affinity:
            postprocess_script += f"""
if affinity_data:
    pd.DataFrame(affinity_data).to_csv('{self.output_folder}/affinity_scores.csv', index=False)
"""
        
        postprocess_script += f"""

# Create sequences CSV (placeholder - would extract from input config)
sequence_data = [{{'id': 'seq_001', 'sequence': 'PLACEHOLDER'}}]
pd.DataFrame(sequence_data).to_csv('{self.output_folder}/sequences.csv', index=False)

print('Post-processing completed')
PYEOF

# Run post-processing
python {self.output_folder}/postprocess.py

"""
        return postprocess_script
    
    def _generate_postprocess_with_script(self) -> str:
        """
        Generate script section for post-processing using dedicated pipe_boltz_postprocessing.py script.
        
        Returns:
            Bash script content for post-processing with proper sequence ID handling
        """
        # Path to the post-processing script
        postprocess_script_path = os.path.join(self.folders["HelpScripts"], "pipe_boltz_postprocessing.py")
        
        # Create a sequences file for post-processing if we have direct sequence input
        if hasattr(self, '_cache_protein_sequence') and self._cache_protein_sequence:
            sequences_file = os.path.join(self.output_folder, "input_sequences.csv")
            effective_job_name = self.get_effective_job_name()
            if effective_job_name is None:
                sequence_id = "prediction"
            else:
                sequence_id = effective_job_name
                
            return f"""
echo "Post-processing Boltz2 results with sequence ID mapping"

# Create sequences file for post-processing
cat > {sequences_file} << 'EOF'
id,sequence
{sequence_id},{self._cache_protein_sequence}
EOF

# Run the post-processing script
python {postprocess_script_path} {self.output_folder} {self.output_folder} {sequences_file}

echo "Post-processing completed - structures renamed with sequence IDs"

"""
        else:
            # Use the original queries CSV file if available - prioritize input tool's sequences file
            sequence_ids_file = getattr(self, "queries_csv_file", None)
            if sequence_ids_file is None:
                # Fallback to Boltz2's own queries_csv if no input tool file
                sequence_ids_file = getattr(self, "queries_csv", None)
            if sequence_ids_file is None:
                # Try to extract sequence if available
                protein_sequence = self._extract_protein_sequence()
                sequences_file = os.path.join(self.output_folder, "input_sequences.csv")
                effective_job_name = self.get_effective_job_name()
                if effective_job_name is None:
                    sequence_id = "prediction"
                else:
                    sequence_id = effective_job_name
                
                return f"""
echo "Post-processing Boltz2 results with generated sequence ID"

# Create sequences file for post-processing
cat > {sequences_file} << 'EOF'
id,sequence
{sequence_id},{protein_sequence or ""}
EOF

# Run the post-processing script with generated CSV file
python {postprocess_script_path} {self.output_folder} {self.output_folder} {sequences_file}
"""
            else:
                return f"""
echo "Post-processing Boltz2 results with original sequence ID mapping"

# Run the post-processing script with original CSV file
python {postprocess_script_path} {self.output_folder} {self.output_folder} {sequence_ids_file}

echo "Post-processing completed - structures renamed with sequence IDs"

"""
    
    def get_config_display(self) -> List[str]:
        """Get configuration display lines for pipeline output."""
        config_lines = super().get_config_display()
        
        if self.ligands:
            config_lines.append(f"Ligands: {self.ligands}")
        elif hasattr(self, 'input_compounds') and self.input_compounds:
            config_lines.append(f"Compounds: {os.path.basename(self.input_compounds[0])} (from CompoundLibrary)")
        elif self.ligand_library:
            config_lines.append(f"Ligand library: {os.path.basename(self.ligand_library)}")
            if self.primary_key:
                config_lines.append(f"Primary key: {self.primary_key}")
        
        config_lines.extend([
            f"Library type: {self.library_type}",
            f"Output format: {self.output_format}",
            f"MSA server: {self.msa_server}",
            f"Affinity calculation: {self.affinity}"
        ])
        
        # Add advanced parameters if specified
        if self.recycling_steps is not None:
            config_lines.append(f"Recycling steps: {self.recycling_steps}")
        
        if self.diffusion_samples is not None:
            config_lines.append(f"Diffusion samples: {self.diffusion_samples}")

        if self.use_potentials:
            config_lines.append(f"Use potentials: {self.use_potentials}")

        return config_lines
    
    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including Boltz2-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "boltz2_params": {
                "config": self.config,
                "proteins": str(self.proteins) if not isinstance(self.proteins, (str, type(None))) else self.proteins,
                "ligands": self.ligands,
                "msas": str(self.msas) if not isinstance(self.msas, (str, type(None))) else self.msas,
                "ligand_library": self.ligand_library,
                "primary_key": self.primary_key,
                "library_repr": self.library_repr,
                "library_type": self.library_type,
                "affinity": self.affinity,
                "output_format": self.output_format,
                "msa_server": self.msa_server,
                "recycling_steps": self.recycling_steps,
                "diffusion_samples": self.diffusion_samples,
                "use_potentials": self.use_potentials
            }
        })
        return base_dict