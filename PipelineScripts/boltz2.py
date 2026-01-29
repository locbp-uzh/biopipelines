"""
Boltz2 configuration for protein-ligand complex prediction.

Handles apo and holo structure prediction with MSA caching,
ligand binding affinity calculation, and comprehensive analysis.
"""

import os
import json
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo
    from .combinatorics import generate_combinatorics_config, get_mode, CombinatoricsConfig, predict_output_ids
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo
    from combinatorics import generate_combinatorics_config, get_mode, CombinatoricsConfig, predict_output_ids


class Boltz2(BaseConfig):
    """
    Boltz2 configuration for protein-ligand complex prediction.
    
    Predicts both apo (protein-only) and holo (protein-ligand) structures
    with automatic MSA management and comprehensive scoring.
    """
    
    # Tool identification
    TOOL_NAME = "Boltz2"
    
    
    def __init__(self,
                 # Primary input parameters
                 config: Optional[str] = None,
                 proteins: Union[str, List[str], ToolOutput, StandardizedOutput] = None,
                 ligands: Union[str, ToolOutput, StandardizedOutput, None] = None,
                 msas: Optional[Union[str, ToolOutput]] = None,
                 # Library-based ligand inputs
                 ligand_library: Optional[str] = None,
                 primary_key: Optional[str] = None,
                 library_repr: str = "SMILES",
                 library_type: str = "noncovalent",
                 # Core prediction parameters
                 affinity: bool = True,
                 output_format: str = "pdb",
                 msa_server: str = "public",
                 # Advanced prediction parameters
                 recycling_steps: Optional[int] = None,
                 diffusion_samples: Optional[int] = None,
                 use_potentials: bool = False,
                 # Template parameters
                 template: Optional[str] = None,
                 template_chain_ids: Optional[List[str]] = None,
                 template_force: bool = True,
                 template_threshold: float = 5.0,
                 # Pocket constraint parameters
                 pocket_residues: Optional[List[int]] = None,
                 pocket_max_distance: float = 7.0,
                 pocket_force: bool = True,
                 # Glycosylation parameters
                 glycosylation: Optional[Dict[str, List[int]]] = None,
                 # Covalent linkage parameters
                 covalent_linkage: Optional[Dict[str, Any]] = None,
                 **kwargs):
        """
        Initialize Boltz2 configuration.

        Args:
            config: Direct YAML configuration string (Example 1: Boltz2(config=yaml))
            proteins: Protein sequences - can be ToolOutput, StandardizedOutput, file path, or direct sequence
            ligands: Single ligand SMILES string, ToolOutput with compounds, or table reference
            msas: MSA files for recycling (e.g., boltz2_apo.tables.msas)
            ligand_library: Path to CSV file with ligand library
            primary_key: Key column in library to filter by
            library_repr: Ligand representation ("SMILES" or "CCD")
            library_type: Binding type ("noncovalent" or "covalent")
            affinity: Whether to calculate binding affinity
            output_format: Output format ("pdb" or "mmcif")
            msa_server: MSA generation ("public" or "local")
            recycling_steps: Number of recycling steps (default: Boltz2 default 3)
            diffusion_samples: Number of diffusion samples (default: Boltz2 default 1)
            use_potentials: Enable potentials for improved structure prediction
            template: Path to PDB template file for structure guidance
            template_chain_ids: List of chain IDs to apply template to (e.g., ["A", "B"])
            template_force: Whether to force template usage (default: True)
            template_threshold: RMSD threshold for template matching (default: 5.0)
            pocket_residues: List of residue positions defining binding pocket (e.g., [50, 51, 52])
            pocket_max_distance: Maximum distance for pocket constraint (default: 7.0)
            pocket_force: Whether to force pocket constraint (default: True)
            glycosylation: Dict mapping chain IDs to Asn positions for N-glycosylation (e.g., {"A": [164]})
            covalent_linkage: Dict specifying covalent attachment (e.g., {"chain": "A", "position": 50, "protein_atom": "SG", "ligand_atom": "C1"})
            **kwargs: Additional parameters
        """
        # Initialize default values
        self.config = config
        self.proteins = proteins
        self.ligands = ligands
        self.msas = msas
        self.input_sequences = None
        self.input_compounds = []
        self.input_tables = {}
        self.input_is_tool_output = False
        self.standardized_input = None

        # Combinatorics modes (default is "each" for cartesian product)
        self._proteins_mode = get_mode(proteins)
        self._ligands_mode = get_mode(ligands)

        # Handle explicit parameter inputs from previous tools
        if isinstance(proteins, StandardizedOutput):
            # Explicit: proteins=tool
            self.input_sequences = proteins.sequences
            self.input_tables = getattr(proteins, 'tables', {})
            self.input_is_tool_output = False
            self.standardized_input = proteins
            # Keep proteins as StandardizedOutput reference
        elif isinstance(proteins, ToolOutput):
            # ToolOutput for proteins
            self.input_sequences = proteins
            self.input_tables = proteins.get_output_files("tables")
            self.input_is_tool_output = True
            self.standardized_input = None
            self.dependencies.append(proteins.config)
        elif isinstance(ligands, (StandardizedOutput, ToolOutput)):
            # Explicit: ligands=compounds from previous tool
            if isinstance(ligands, StandardizedOutput):
                self.input_compounds = getattr(ligands, 'compounds', [])
                self.input_tables = getattr(ligands, 'tables', {})
            else:  # ToolOutput
                self.input_compounds = ligands.get_output_files("compounds")
                self.input_tables = ligands.get_output_files("tables")
                # Add dependency for ToolOutput
                self.dependencies.append(ligands.config)
            # Keep ligands as tool reference
        elif isinstance(msas, StandardizedOutput):
            # Explicit: msas=previous_boltz
            self.input_tables = getattr(msas, 'tables', {})
        else:
            # Direct inputs (strings, lists, etc.)
            self.input_sequences = self.proteins
            self.input_compounds = []
            self.input_tables = {}
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
        self.recycling_steps = recycling_steps
        self.diffusion_samples = diffusion_samples
        self.use_potentials = use_potentials
        # Template parameters
        self.template = template
        self.template_chain_ids = template_chain_ids
        self.template_force = template_force
        self.template_threshold = template_threshold
        # Pocket constraint parameters
        self.pocket_residues = pocket_residues
        self.pocket_max_distance = pocket_max_distance
        self.pocket_force = pocket_force
        # Glycosylation parameters
        self.glycosylation = glycosylation
        # Covalent linkage parameters
        self.covalent_linkage = covalent_linkage
        
        # Track input source type and files (tool-agnostic)
        self.input_yaml_entities = None
        self.input_fasta_files = []
        self.queries_csv_file = None
        
        # Initialize base class
        super().__init__(**kwargs)

        # Initialize file paths (will be set in configure_inputs)
        self._initialize_file_paths()


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
        # Input list file (avoids "Argument list too long" with many fasta files)
        self.fasta_files_list_file = os.path.join(self.output_folder, ".input_fasta_files.txt")
        
        # Output files
        self.library_scores_csv = os.path.join(self.library_folder, "library_scores.csv")
        self.config_txt_file = os.path.join(self.output_folder, f"{self.job_name}_config.txt")
        self.results_zip = os.path.join(self.output_folder, f"{self.job_name}.zip")
        
        # Helper script paths (only set if folders are available)
        if hasattr(self, 'folders') and self.folders:
            self.smiles_library_py = os.path.join(self.folders["HelpScripts"], "pipe_smiles_library.py")
            self.boltz_config_unified_py = os.path.join(self.folders["HelpScripts"], "pipe_boltz_config_unified.py")
            self.boltz_results_py = os.path.join(self.folders["HelpScripts"], "pipe_boltz_results.py")
            self.mmseqs2_client_sh = os.path.join(self.folders.get("MMseqs2", "MMseqs2"), "mmseqs2_client.sh")
            self.fa_to_csv_py = os.path.join(self.folders["HelpScripts"], "pipe_fa_to_csv_fasta.py")
        else:
            # Temporary placeholders when folders aren't available yet
            self.smiles_library_py = None
            self.boltz_config_unified_py = None
            self.boltz_results_py = None
            self.mmseqs2_client_sh = None
            self.fa_to_csv_py = None

    def validate_params(self):
        """Validate Boltz2-specific parameters."""
        # Must have some form of input
        has_input = any([
            getattr(self, 'config', None) is not None,
            getattr(self, 'proteins', None) is not None,
            getattr(self, 'input_sequences', None) is not None
        ])
        if not has_input:
            raise ValueError("Either config or proteins parameter is required")
        
        # Cannot specify multiple primary input methods (only check if they exist)
        primary_inputs = []
        if getattr(self, 'config', None) is not None:
            primary_inputs.append(self.config)
        if getattr(self, 'proteins', None) is not None:
            primary_inputs.append(self.proteins)
        
        if len(primary_inputs) > 1:
            raise ValueError("Cannot specify multiple primary input methods (config and proteins)")
        
        # For proteins input (not config), ligand information is typically required
        # But with standardized input, ligands might come from tables
        if (self.proteins and not self.config and 
            not self.ligands and not self.ligand_library and not self.msas and
            not self.standardized_input):
            # Only warn, don't error - ligands might be in tables or MSAs
            pass
        
        # Cannot specify both single ligand and library
        if self.ligands and self.ligand_library:
            raise ValueError("Cannot specify both ligands and ligand_library")
        
        # Validate table references if provided (skip for StandardizedOutput)
        if self.ligands and isinstance(self.ligands, str):
            self.validate_table_reference(self.ligands)
        
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
                project_path = os.path.join(pipeline_folders["biopipelines"], self.ligand_library)
                if os.path.exists(project_path):
                    self.ligand_library = project_path
                else:
                    raise ValueError(f"Ligand library file not found: {self.ligand_library}")
        
        # Handle MSA inputs - build YAML config with MSA paths
        if self.msas:
            self._integrate_msa_inputs()
    
    def _integrate_msa_inputs(self):
        """Build YAML config that includes MSA paths from previous tool."""
        if not self.msas:
            return
            
        # Get MSA table from previous tool
        msa_paths = {}
        if hasattr(self.msas, 'tables'):
            if hasattr(self.msas, 'tables') and 'msas' in self.msas.tables:
                # MSAs should be provided in the tables.msas
                # Handle both TableContainer (returns path string) and dict format
                msa_table = self.msas.tables['msas']
                if isinstance(msa_table, str):
                    # TableContainer returns path directly
                    msa_table_path = msa_table
                elif isinstance(msa_table, dict) and 'path' in msa_table:
                    # Raw dict format
                    msa_table_path = msa_table['path']
                else:
                    raise ValueError(f"Invalid MSA table format: {type(msa_table)}")
                
                # This will be a CSV with columns: id, sequence_id, msa_file
                # We need to read this at runtime, not pipeline time
                # For now, just signal that MSAs will be integrated at runtime
                pass
            else:
                raise ValueError("MSA ToolOutput doesn't have tables.msas")
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
    
    def _generate_extra_config_params(self) -> str:
        """
        Generate extra CLI parameters for pipe_boltz_config_unified.py.

        Returns:
            String with extra CLI arguments for template, pocket, glycosylation, covalent linkage
        """
        extra_params = []

        # Template parameters
        if self.template:
            extra_params.append(f'--template "{self.template}"')
            if self.template_chain_ids:
                extra_params.append(f'--template-chains "{",".join(self.template_chain_ids)}"')
            if self.template_force:
                extra_params.append('--template-force')
            extra_params.append(f'--template-threshold {self.template_threshold}')

        # Pocket constraint parameters
        if self.pocket_residues:
            extra_params.append(f'--pocket-residues "{self.pocket_residues}"')
            extra_params.append(f'--pocket-max-distance {self.pocket_max_distance}')
            if self.pocket_force:
                extra_params.append('--pocket-force')

        # Glycosylation parameters
        if self.glycosylation:
            glyco_json = json.dumps(self.glycosylation)
            extra_params.append(f"--glycosylation '{glyco_json}'")

        # Covalent linkage parameters
        if self.covalent_linkage:
            covalent_json = json.dumps(self.covalent_linkage)
            extra_params.append(f"--covalent-linkage '{covalent_json}'")

        return " ".join(extra_params)

    def _write_combinatorics_config(self) -> str:
        """
        Write combinatorics config file at pipeline time.

        Returns:
            Path to the generated config file
        """
        config_path = os.path.join(self.output_folder, "combinatorics_config.json")
        generate_combinatorics_config(
            config_path,
            proteins=self.proteins,
            ligands=self.ligands
        )
        return config_path

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
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()

        # Create basic folder structure
        script_content += f"""# Create output folders
mkdir -p {os.path.join(self.output_folder, "predictions")}
mkdir -p {self.msa_cache_folder}

"""

        # Handle MSA recycling if provided
        if self.msas:
            script_content += self._generate_msa_recycling_section()

        # Generate input configuration based on input type
        effective_job_name = self.get_effective_job_name()
        if effective_job_name is None:
            config_name = "prediction"
        else:
            config_name = effective_job_name
        config_file_path = os.path.join(self.output_folder, f"{config_name}.yaml")

        # Write combinatorics config (describes proteins/ligands modes)
        combinatorics_config_path = self._write_combinatorics_config()

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
        else:
            # Use unified config generation via pipe_boltz_config_unified.py
            script_content += self._generate_unified_config_section(
                combinatorics_config_path, config_name
            )

        # Run Boltz2 prediction
        boltz_options = f"--cache {boltz_cache_folder} --out_dir {self.output_folder}{msa_option} --output_format {self.output_format}"

        # Add recycling_steps and diffusion_samples if specified
        if self.recycling_steps is not None:
            boltz_options += f" --recycling_steps {self.recycling_steps}"

        if self.diffusion_samples is not None:
            boltz_options += f" --diffusion_samples {self.diffusion_samples}"

        if self.use_potentials:
            boltz_options += " --use_potentials"

        # Determine if we use multiple config files or single
        uses_unified_config = not self.config and not self.input_yaml_entities

        if uses_unified_config:
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

        # Propagate missing table from upstream tools
        script_content += self._generate_missing_table_propagation()

        script_content += self._generate_boltz2_completion_check()

        return script_content

    def _generate_unified_config_section(self, combinatorics_config_path: str, config_name: str) -> str:
        """
        Generate script section that uses the unified pipe_boltz_config_unified.py.

        Args:
            combinatorics_config_path: Path to the combinatorics config JSON
            config_name: Job/config name for naming outputs

        Returns:
            Bash script content for config generation
        """
        config_files_dir = os.path.join(self.output_folder, "config_files")
        pipe_boltz_config_py = os.path.join(self.folders["HelpScripts"], "pipe_boltz_config_unified.py")

        script = f"""
echo "Generating Boltz2 configurations using unified config generator"
mkdir -p {config_files_dir}

"""
        # Generate proteins CSV if needed
        script += self._generate_proteins_csv_section()

        # Generate ligands CSV if needed
        script += self._generate_ligands_csv_section(config_name)

        # Build the command for unified config generator
        cmd_parts = [
            f'python {pipe_boltz_config_py}',
            f'--combinatorics-config "{combinatorics_config_path}"',
            f'--output-dir "{config_files_dir}"',
        ]

        # Add MSA table if available
        msa_table_flag = self._get_msa_table_flag()
        if msa_table_flag:
            cmd_parts.append(msa_table_flag)

        # Add affinity flag
        if self.affinity:
            cmd_parts.append('--affinity')

        # Add extra config parameters
        extra_params = self._generate_extra_config_params()
        if extra_params:
            cmd_parts.append(extra_params)

        script += '# Generate config files\n'
        script += ' \\\n    '.join(cmd_parts) + '\n\n'

        return script

    def _generate_proteins_csv_section(self) -> str:
        """Generate bash script section to create proteins CSV if needed."""
        # If we have a direct sequence, create a CSV from it
        if hasattr(self, 'input_direct_sequence') and self.input_direct_sequence:
            proteins_csv = os.path.join(self.output_folder, "proteins.csv")
            return f"""# Create proteins CSV from direct sequence
cat > {proteins_csv} << 'EOF'
id,sequence
protein,{self.input_direct_sequence}
EOF

"""
        # If we have queries_csv_file, it's already a CSV - no need to create
        # The combinatorics config already points to the source files
        return ""

    def _generate_ligands_csv_section(self, config_name: str) -> str:
        """Generate bash script section to create ligands CSV if needed."""
        # Check if ligands is a direct SMILES string (not a CSV path or tool output)
        if self.ligands and isinstance(self.ligands, str) and not self.ligands.endswith('.csv'):
            ligands_csv = os.path.join(self.output_folder, "ligands.csv")
            ligand_id = config_name if config_name != "prediction" else "ligand"
            return f"""# Create ligands CSV from direct SMILES string
cat > {ligands_csv} << 'EOF'
id,format,smiles,ccd
{ligand_id},smiles,{self.ligands},
EOF

"""
        return ""

    def _get_msa_table_flag(self) -> str:
        """Get the MSA table flag for the config generator."""
        if not self.msas:
            return ""

        if hasattr(self.msas, 'tables'):
            # ToolOutput or StandardizedOutput case
            if hasattr(self.msas.tables, '_tables') and 'msas' in self.msas.tables._tables:
                msa_table_path = self.msas.tables._tables['msas'].path
                return f'--msa-table "{msa_table_path}"'
            elif hasattr(self.msas.tables, 'msas'):
                msa_table_path = self.msas.tables.msas
                return f'--msa-table "{msa_table_path}"'
        elif isinstance(self.msas, str) and self.msas.endswith('.csv'):
            return f'--msa-table "{self.msas}"'

        return ""

    def _get_sequences_file_path(self) -> str:
        """Get the path to the sequences file used for input."""
        # Check if we're using protein queries (most common case)
        if hasattr(self, 'queries_csv_file') and self.queries_csv_file:
            return self.queries_csv_file

        # Check if we have tool input with sequences
        if hasattr(self, 'proteins') and hasattr(self.proteins, 'sequences') and self.proteins.sequences:
            return self.proteins.sequences[0]

        # No valid sequences file source found - fail explicitly
        raise ValueError("Cannot determine sequences file path: no valid source found (queries_csv_file, proteins.sequences)")

    def _generate_boltz2_completion_check(self) -> str:
        """
        Generate custom completion check that accounts for missing sequences from upstream tools.
        """
        # Check if we need to filter missing sequences
        has_missing_sequences = False
        missing_table_path = None

        # Try to find missing sequences table in proteins input
        if hasattr(self.proteins, 'tables'):
            tables = self.proteins.tables
            if hasattr(tables, '_tables'):
                # Standard BioPipelines format
                for name, info in tables._tables.items():
                    if 'missing' in name.lower():
                        has_missing_sequences = True
                        if hasattr(info, 'path'):
                            missing_table_path = info.path
                        elif isinstance(info, str):
                            missing_table_path = info
                        break
            elif isinstance(tables, dict):
                # Dict format
                for name, info in tables.items():
                    if 'missing' in name.lower():
                        has_missing_sequences = True
                        if isinstance(info, str):
                            missing_table_path = info
                        elif isinstance(info, dict) and 'path' in info:
                            missing_table_path = info['path']
                        break

        if has_missing_sequences and missing_table_path:
            return f"""
# Custom completion check that accounts for missing sequences
echo "Checking outputs and filtering missing sequences..."

# Create a temporary file to store expected sequence IDs
expected_ids_file=$(mktemp)

# Get sequence IDs excluding missing ones
python -c "
import pandas as pd
import sys

try:
    # Read main sequences
    sequences_df = pd.read_csv('{self._get_sequences_file_path()}')
    all_ids = set(sequences_df['id'].tolist())
    print(f'Found {{len(all_ids)}} total sequence IDs')

    # Read missing sequences
    try:
        missing_df = pd.read_csv('{missing_table_path}')
        missing_ids = set(missing_df['id'].tolist())
        print(f'Found {{len(missing_ids)}} missing sequence IDs: {{sorted(missing_ids)}}')
    except:
        missing_ids = set()
        print('No missing sequences found')

    # Expected IDs are all IDs minus missing IDs
    expected_ids = all_ids - missing_ids
    print(f'Expecting {{len(expected_ids)}} output structures')

    # Write expected IDs to file
    with open('$expected_ids_file', 'w') as f:
        for seq_id in sorted(expected_ids):
            f.write(seq_id + '\\n')

except Exception as e:
    print(f'Error processing sequences: {{e}}')
    sys.exit(1)
"

# Check if expected output files exist
missing_outputs=()
while IFS= read -r seq_id; do
    expected_file="{self.output_folder}/${{seq_id}}.{self.output_format.lower()}"
    if [ ! -f "$expected_file" ]; then
        missing_outputs+=("$expected_file")
    fi
done < "$expected_ids_file"

# Clean up temporary file
rm -f "$expected_ids_file"

# Report results
if [ ${{#missing_outputs[@]}} -eq 0 ]; then
    echo "All expected outputs found"
    echo "Boltz2 completed successfully"
else
    echo "Missing critical outputs for Boltz2:"
    echo "  structures:"
    printf '    - %s\\n' "${{missing_outputs[@]}}"
    echo "Boltz2 failed - some outputs missing"
    exit 1
fi
""" + self.generate_completion_check_footer()
        else:
            # No missing sequences, use standard completion check
            return self.generate_completion_check_footer()
    
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
        
        # Organize tables by content type with detailed metadata
        tables = {}
        
        # Confidence scores table (converted from JSON to CSV for compatibility)
        confidence_csv = os.path.join(self.output_folder, "confidence_scores.csv")
        tables["confidence"] = TableInfo(
            name="confidence",
            path=confidence_csv,
            columns=["id", "input_file", "confidence_score", "ptm", "iptm", "complex_plddt", "complex_iplddt"],
            description="Boltz2 confidence scores for structure predictions",
            count=len(input_file_names)
        )
        
        # Affinity scores table (if affinity calculation enabled)
        if self.affinity:
            affinity_csv = os.path.join(self.output_folder, "affinity_scores.csv")
            tables["affinity"] = TableInfo(
                name="affinity",
                path=affinity_csv,
                columns=["id", "input_file", "affinity_pred_value", "affinity_probability_binary"],
                description="Boltz2 binding affinity predictions",
                count=len(input_file_names)
            )
        
        # MSAs table (for recycling between apo/holo predictions)
        if msa_files:
            msa_csv = os.path.join(self.output_folder, "msas.csv")
            tables["msas"] = TableInfo(
                name="msas",
                path=msa_csv,
                columns=["id", "sequence_id", "sequence", "msa_file"],
                description="MSA files for sequence recycling between predictions",
                count=len(msa_files)
            )
        
        # Input sequences table (for downstream tools)
        sequences_csv = os.path.join(self.output_folder, "sequences.csv")
        tables["sequences"] = TableInfo(
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

            # Add compounds table pointing to the ligand library
            tables["compounds"] = TableInfo(
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

            # Add compounds table for single SMILES
            tables["compounds"] = TableInfo(
                name="compounds",
                path=compounds_csv,
                columns=["id", "format", "smiles", "ccd"],
                description="Single ligand compound used for binding prediction",
                count=1
            )

        # Add missing table if there's an upstream missing table to propagate
        upstream_missing_path = self._get_upstream_missing_table_path(
            self.proteins,
            self.input_sequences,
            getattr(self, 'standardized_input', None)
        )

        if upstream_missing_path:
            # Propagate missing table - IDs will be updated at SLURM runtime
            missing_csv = os.path.join(self.output_folder, "missing.csv")
            tables["missing"] = TableInfo(
                name="missing",
                path=missing_csv,
                columns=["id", "structure", "msa"],
                description="Sequences filtered out by upstream tools",
                count="variable"
            )

        return {
            "structures": structures,
            "structure_ids": structure_ids,
            "compounds": compounds,
            "compound_ids": compound_ids,
            "sequences": [sequences_csv],  # Main sequence file
            "sequence_ids": sequence_ids,
            "msas": msa_files,  # Individual MSA files for recycling
            "tables": tables,
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
    
    def _get_missing_sequence_ids(self) -> List[str]:
        """
        Get sequence IDs that are marked as missing in upstream tools.

        Returns:
            List of missing sequence IDs
        """
        missing_ids = []

        # Check proteins input for missing sequences table
        if hasattr(self.proteins, 'tables'):
            tables = self.proteins.tables
            if hasattr(tables, '_tables'):
                # Standard BioPipelines format
                for name, info in tables._tables.items():
                    if 'missing' in name.lower():
                        # This is a missing sequences table - we'll extract IDs at runtime
                        # For now, we can't predict them, so we'll handle this in the script
                        pass
            elif isinstance(tables, dict):
                # Dict format
                for name, info in tables.items():
                    if 'missing' in name.lower():
                        # This is a missing sequences table
                        pass

        return missing_ids

    def _predict_sequence_ids(self) -> List[str]:
        """
        Predict sequence IDs from input sources using combinatorics module.

        Returns:
            List of expected sequence IDs based on combinatorics modes
        """
        # Use the combinatorics module's predict_output_ids function
        # This ensures consistency between pipeline-time prediction and SLURM-time generation
        return predict_output_ids(
            bundled_name="bundled_complex",
            proteins=self.proteins,
            ligands=self.ligands
        )
    
    def _generate_msa_recycling_section(self) -> str:
        """
        Generate script section for MSA recycling.
        
        Returns:
            Bash script content for MSA handling
        """
        if hasattr(self.msas, 'tables'):
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
                # MSA table CSV - read it and copy individual MSA files
                return f"""
echo "Using MSAs from table: {self.msas}"
mkdir -p {self.msa_cache_folder}
# Read MSA table and copy individual MSA files
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
    
    def _generate_missing_table_propagation(self) -> str:
        """
        Generate script section to propagate missing.csv from upstream tools.

        Returns:
            Bash script content for missing table propagation
        """
        # Check if there's an upstream missing table to propagate
        upstream_missing_path = self._get_upstream_missing_table_path(
            self.proteins,
            self.input_sequences,
            getattr(self, 'standardized_input', None)
        )

        if not upstream_missing_path:
            # No upstream missing table - nothing to propagate
            return ""

        # Get the upstream output folder (parent of missing.csv)
        upstream_folder = os.path.dirname(upstream_missing_path)

        # Path to propagation script
        propagate_script = os.path.join(self.folders["HelpScripts"], "pipe_propagate_missing.py")

        # Determine file extensions based on output format
        structure_ext = ".pdb" if self.output_format == "pdb" else ".cif"
        msa_ext = ".csv" if self.msa_server == "public" else ".a3m"

        return f"""
# Propagate missing table from upstream tools
echo "Checking for upstream missing sequences..."
if [ -f "{upstream_missing_path}" ]; then
    echo "Found upstream missing.csv - propagating to current tool"
    python {propagate_script} \\
        --upstream-folders "{upstream_folder}" \\
        --output-folder "{self.output_folder}" \\
        --structure-ext "{structure_ext}" \\
        --msa-ext "{msa_ext}"
else
    echo "No upstream missing.csv found"
fi

"""

    def _generate_postprocess_with_script(self) -> str:
        """
        Generate script section for post-processing using dedicated pipe_boltz_postprocessing.py script.

        Returns:
            Bash script content for post-processing with proper sequence ID handling
        """
        # Path to the post-processing script
        postprocess_script_path = os.path.join(self.folders["HelpScripts"], "pipe_boltz_postprocessing.py")

        # Use the sequence_ids.csv file generated by pipe_boltz_config_unified.py
        # This file is written to the parent of config_files dir
        sequence_ids_file = os.path.join(self.output_folder, "sequence_ids.csv")

        # Fallback: check for queries_csv_file from input
        fallback_file = getattr(self, "queries_csv_file", None)
        if fallback_file is None:
            fallback_file = getattr(self, "queries_csv", None)

        return f"""
echo "Post-processing Boltz2 results"

# Use sequence_ids.csv generated by config generator, or fallback to input CSV
if [ -f "{sequence_ids_file}" ]; then
    python {postprocess_script_path} {self.output_folder} {self.output_folder} {sequence_ids_file}
elif [ -f "{fallback_file or ''}" ]; then
    python {postprocess_script_path} {self.output_folder} {self.output_folder} {fallback_file}
else
    echo "Warning: No sequence IDs file found, skipping post-processing"
fi

echo "Post-processing completed"

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

        # Template parameters
        if self.template:
            config_lines.append(f"Template: {os.path.basename(self.template)}")
            if self.template_chain_ids:
                config_lines.append(f"Template chains: {', '.join(self.template_chain_ids)}")

        # Pocket constraint parameters
        if self.pocket_residues:
            config_lines.append(f"Pocket residues: {self.pocket_residues}")

        # Glycosylation parameters
        if self.glycosylation:
            config_lines.append(f"Glycosylation: {self.glycosylation}")

        # Covalent linkage parameters
        if self.covalent_linkage:
            config_lines.append(f"Covalent linkage: {self.covalent_linkage}")

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
                "use_potentials": self.use_potentials,
                "template": self.template,
                "template_chain_ids": self.template_chain_ids,
                "template_force": self.template_force,
                "template_threshold": self.template_threshold,
                "pocket_residues": self.pocket_residues,
                "pocket_max_distance": self.pocket_max_distance,
                "pocket_force": self.pocket_force,
                "glycosylation": self.glycosylation,
                "covalent_linkage": self.covalent_linkage
            }
        })
        return base_dict