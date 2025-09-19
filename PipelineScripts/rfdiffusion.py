"""
RFdiffusion configuration for protein backbone generation.

Handles both standard RFdiffusion and RFdiffusion-AllAtom workflows
with proper parameter validation and script generation.
"""

import os
import shutil
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput


class RFdiffusion(BaseConfig):
    """
    Configuration for RFdiffusion protein backbone generation.
    
    Supports both standard RFdiffusion and RFdiffusion-AllAtom modes
    with comprehensive parameter validation and automatic script generation.
    """
    
    TOOL_NAME = "RFdiffusion"
    DEFAULT_ENV = "ProteinEnv"
    
    def __init__(self, pdb: Union[str, ToolOutput, StandardizedOutput] = "", contigs: str = "", inpaint: str = "",
                 num_designs: int = 1, active_site: bool = False, 
                 steps: int = 50, partial_steps: int = 0, 
                 reproducible: bool = False, reproducibility_number: int = 0,
                 **kwargs):
        """
        Initialize RFdiffusion configuration.
        
        Args:
            pdb: Input PDB file (optional template), ToolOutput, or StandardizedOutput
            contigs: Contig specification (e.g., "A1-100,10-20")
            inpaint: Inpainting specification (same format as contigs)
            num_designs: Number of designs to generate
            active_site: Use active site model for small motifs
            steps: Diffusion steps (default 50)
            partial_steps: Partial diffusion steps
            reproducible: Use deterministic sampling
            reproducibility_number: Seed for reproducible runs
            **kwargs: Additional parameters
        """
        # Store RFdiffusion-specific parameters
        self.pdb = pdb
        self.pdb_is_tool_output = False
        self.pdb_source_file = None

        # Handle tool output for PDB input
        if isinstance(pdb, (ToolOutput, StandardizedOutput)):
            self.pdb_is_tool_output = True
            if isinstance(pdb, StandardizedOutput):
                # Get first structure from StandardizedOutput
                if hasattr(pdb, 'structures') and pdb.structures:
                    self.pdb_source_file = pdb.structures[0]
                    if len(pdb.structures) > 1:
                        print(f"Warning: Multiple structures found in input, using first one: {os.path.basename(self.pdb_source_file)}")
                        print("Note: Multiple structure support not yet implemented")
                else:
                    raise ValueError("No structures found in StandardizedOutput for pdb parameter")
            else:  # ToolOutput
                # Get first structure from ToolOutput
                structures = pdb.get_output_files("structures")
                if structures:
                    self.pdb_source_file = structures[0]
                    if len(structures) > 1:
                        print(f"Warning: Multiple structures found in tool output, using first one: {os.path.basename(self.pdb_source_file)}")
                        print("Note: Multiple structure support not yet implemented")
                    # Add dependency
                    self.dependencies.append(pdb.config)
                else:
                    raise ValueError("No structures found in ToolOutput for pdb parameter")

        self.contigs = contigs
        self.inpaint = inpaint
        self.num_designs = num_designs
        self.active_site = active_site
        self.steps = steps
        self.partial_steps = partial_steps
        self.reproducible = reproducible
        self.reproducibility_number = reproducibility_number
        
        # Initialize base class
        super().__init__(**kwargs)
        
        # Initialize file paths (will be set in configure_inputs)
        self._initialize_file_paths()
    
    def validate_params(self):
        """Validate RFdiffusion-specific parameters."""
        if not self.contigs:
            raise ValueError("contigs parameter is required for RFdiffusion")
        
        if self.num_designs <= 0:
            raise ValueError("num_designs must be positive")
        
        if self.steps <= 0:
            raise ValueError("steps must be positive")
        
        if self.partial_steps < 0:
            raise ValueError("partial_steps cannot be negative")
        
        # Skip PDB validation for tool outputs (will be validated at runtime)
        if self.pdb and not self.pdb_is_tool_output:
            if not (isinstance(self.pdb, str) and (self.pdb.endswith('.pdb') or os.path.exists(self.pdb))):
                # Check if it exists in PDBs folder
                pdb_path = os.path.join(os.getcwd(), "PDBs", self.pdb + ".pdb" if not self.pdb.endswith('.pdb') else self.pdb)
                if not os.path.exists(pdb_path):
                    raise ValueError(f"PDB file not found: {self.pdb}")
    
    def _initialize_file_paths(self):
        """Initialize common file paths used throughout the class."""
        self.input_pdb_file = None
        self.main_datasheet = None
        self.rfd_log_file = None
        self.pipeline_name = None
        
        # Helper script paths
        self.inference_py_file = None
        self.datasheet_py_file = None
    
    def _setup_file_paths(self):
        """Set up all file paths after output_folder is known."""
        # Extract pipeline name from folder structure
        self.pipeline_name = self._extract_pipeline_name()
        
        # Core output files
        self.main_datasheet = os.path.join(self.output_folder, "rfdiffusion_results.csv")
        
        # Log file is created by pipeline in the main pipeline folder with pattern _{index}_{toolname}.log
        # Extract index from folder name (e.g., "1_RFdiffusion" -> "1")
        folder_name = os.path.basename(self.output_folder)
        pipeline_folder = os.path.dirname(self.output_folder)  # Get parent folder (DeNovoProtein_013)
        
        if '_' in folder_name and folder_name.split('_')[0].isdigit():
            index = folder_name.split('_')[0]
            tool_name = folder_name.split('_', 1)[1].lower()  # Get everything after first underscore, lowercase
            self.rfd_log_file = os.path.join(pipeline_folder, f"_{index:0>3}_{tool_name}.log")
        else:
            raise ValueError(f"Invalid output folder naming pattern: {folder_name}. Expected 'N_RFdiffusion' format.")
        
        # Helper script paths (only set if folders are available)
        if hasattr(self, 'folders') and self.folders:
            self.inference_py_file = os.path.join(self.folders["RFdiffusion"], "scripts", "run_inference.py")
            self.datasheet_py_file = os.path.join(self.folders["HelpScripts"], "pipe_rfdiffusion_datasheet.py")
            
            # Input PDB file path (if PDB is provided)
            if self.pdb_is_tool_output:
                # Use the actual file path from tool output
                self.input_pdb_file = self.pdb_source_file
            elif self.pdb:
                # Handle string PDB filename
                pdb_temp = self.pdb if self.pdb.endswith(".pdb") else self.pdb + ".pdb"
                self.input_pdb_file = os.path.join(self.folders["runtime"], pdb_temp)
            else:
                self.input_pdb_file = ""
        else:
            # Temporary placeholders when folders aren't available yet
            self.inference_py_file = None
            self.datasheet_py_file = None
            self.input_pdb_file = ""
    
    def _extract_pipeline_name(self) -> str:
        """Extract pipeline name from output folder structure."""
        folder_parts = self.output_folder.split(os.sep)
        for i, part in enumerate(folder_parts):
            if "RFdiffusion" in part:
                if i > 0:
                    return folder_parts[i-1]
                break
        raise ValueError(f"Could not extract job name from output folder: {self.output_folder}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input files and dependencies."""
        self.folders = pipeline_folders
        self._setup_file_paths()  # Set up all file paths now that we have folders
        
        # Handle PDB input if provided
        if self.pdb:
            pdb_temp = self.pdb if self.pdb.endswith(".pdb") else self.pdb + ".pdb"
            pdb_source = os.path.join(pipeline_folders["PDBs"], pdb_temp)
            
            if os.path.exists(pdb_source):
                # Will be copied in script generation
                self.input_sources["pdb"] = pdb_source
            else:
                raise ValueError(f"PDB file not found: {pdb_source}")
    
    def get_config_display(self) -> List[str]:
        """Get RFdiffusion configuration display lines."""
        config_lines = super().get_config_display()
        
        config_lines.extend([
            f"CONTIGS: {self.contigs}",
            f"NUM DESIGNS: {self.num_designs}",
            f"ACTIVE SITE: {self.active_site}",
            f"STEPS: {self.steps}"
        ])
        
        if self.pdb:
            config_lines.append(f"PDB: {self.pdb}")
        if self.inpaint:
            config_lines.append(f"INPAINT: {self.inpaint}")
        if self.partial_steps > 0:
            config_lines.append(f"PARTIAL STEPS: {self.partial_steps}")
        if self.reproducible:
            config_lines.append(f"REPRODUCIBLE: {self.reproducible} (seed: {self.reproducibility_number})")
        
        return config_lines

    def generate_script(self, script_path: str) -> str:
        """
        Generate RFdiffusion execution script.
        
        Args:
            script_path: Path where script should be written
            
        Returns:
            Script content as string
        """
        runtime_folder = os.path.dirname(script_path)
        rfd_job_folder = self.output_folder
        os.makedirs(rfd_job_folder, exist_ok=True)
        
        # Generate script content following modular pattern
        script_content = "#!/bin/bash\n"
        script_content += "# RFdiffusion execution script\n"
        script_content += "# Generated by ProteinNotebooks pipeline system\n\n"
        script_content += self.generate_completion_check_header()
        script_content += self.generate_script_run_rfdiffusion()
        script_content += self.generate_script_create_datasheet()
        script_content += self.generate_completion_check_footer()
        
        return script_content
    
    def generate_script_run_rfdiffusion(self) -> str:
        """Generate the RFdiffusion execution part of the script."""
        rfd_job_folder = self.output_folder
        
        # Copy input PDB if provided (only if source and destination are different)
        if self.pdb and "pdb" in self.input_sources:
            source_path = os.path.abspath(self.input_sources["pdb"])
            dest_path = os.path.abspath(self.input_pdb_file)
            if source_path != dest_path:
                shutil.copy(source_path, dest_path)
        
        # Build RFdiffusion options
        rfd_options = f"'contigmap.contigs=[{self.contigs}]'"
        
        if self.inpaint:
            rfd_options += f" 'contigmap.inpaint_seq=[{self.inpaint}]'"
        
        if self.input_pdb_file:
            rfd_options += f" inference.input_pdb={self.input_pdb_file}"
        
        # Output prefix for generated designs (use clean pipeline name for elegant consistency)
        prefix = os.path.join(rfd_job_folder, f"{self.pipeline_name}")
        rfd_options += f" inference.output_prefix={prefix}"
        rfd_options += f" inference.num_designs={self.num_designs}"
        rfd_options += f" inference.deterministic={self.reproducible}"
        rfd_options += f" inference.design_startnum={self.reproducibility_number}"
        
        if self.steps != 50:
            rfd_options += f" diffuser.T={self.steps}"
        
        if self.partial_steps > 0:
            rfd_options += f" diffuser.partial_T={self.partial_steps}"
        
        if self.active_site:
            rfd_options += " inference.ckpt_override_path=models/ActiveSite_ckpt.pt"
        
        return f"""echo "Starting RFdiffusion"
echo "Options: {rfd_options}"
echo "Output folder: {rfd_job_folder}"

# Run RFdiffusion
python {self.inference_py_file} {rfd_options}

"""

    def generate_script_create_datasheet(self) -> str:
        """Generate the datasheet creation part of the script."""
        rfd_job_folder = self.output_folder
        
        # Design character: '-' for RFdiffusion, '?' for RFdiffusion-AllAtom
        design_character = "-"
        
        return f"""echo "Creating results datasheet"
# Create main datasheet with id, pdb, fixed, designed columns by parsing RFdiffusion log
python {self.datasheet_py_file} "{rfd_job_folder}" "{self.rfd_log_file}" "{design_character}" "{self.pipeline_name}" {self.num_designs} "{self.main_datasheet}"

"""
    
    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after RFdiffusion execution.
        
        Uses pure path construction - no filesystem access.
        Returns expected paths based on naming patterns.
        
        Returns:
            Dictionary mapping output types to file paths with standard keys:
            - structures: PDB files
            - compounds: Empty (no compounds from RFdiffusion)
            - sequences: Empty (no sequences from RFdiffusion)  
            - datasheets: Main datasheet CSV
            - output_folder: Tool's output directory
        """
        # Ensure file paths are set up
        if not hasattr(self, 'pipeline_name') or self.pipeline_name is None:
            # Fallback if configure_inputs hasn't been called yet
            self._setup_file_paths()
        
        pipeline_name = self.pipeline_name
        main_datasheet = self.main_datasheet
        
        # Generate expected PDB file paths based on clean naming pattern
        design_pdbs = []
        structure_ids = []
        for i in range(self.num_designs):
            design_id = f"{pipeline_name}_{i}"
            design_path = os.path.join(self.output_folder, f"{design_id}.pdb")
            design_pdbs.append(design_path)
            structure_ids.append(design_id)
        
        # Organize datasheets by content type
        datasheets = {
            "structures": {
                "path": main_datasheet,
                "columns": ["id", "source_id", "pdb", "fixed", "designed", "contigs", "time", "status"],
                "description": "RFdiffusion structure generation results with fixed/designed regions",
                "count": self.num_designs
            }
        }
        
        return {
            "structures": design_pdbs,
            "structure_ids": structure_ids,
            "compounds": [],
            "compound_ids": [],
            "sequences": [],
            "sequence_ids": [],
            "datasheets": datasheets,
            "output_folder": self.output_folder,
            # Keep legacy aliases for compatibility
            "pdbs": design_pdbs,
            "main": main_datasheet  # Legacy alias for backward compatibility
        }
    
    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including RFdiffusion-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "rfd_params": {
                "pdb": self.pdb,
                "contigs": self.contigs,
                "inpaint": self.inpaint,
                "num_designs": self.num_designs,
                "active_site": self.active_site,
                "steps": self.steps,
                "partial_steps": self.partial_steps,
                "reproducible": self.reproducible,
                "reproducibility_number": self.reproducibility_number
            }
        })
        return base_dict


