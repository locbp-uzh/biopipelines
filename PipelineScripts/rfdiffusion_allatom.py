"""
RFdiffusion-AllAtom configuration for ligand-aware protein design.

Handles RFdiffusion-AllAtom workflows with ligand contexts, PPI design,
and all-atom generation capabilities.
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


class RFdiffusionAllAtom(BaseConfig):
    """
    RFdiffusion-AllAtom variant for ligand-aware protein design.
    
    Extends base RFdiffusion with support for ligand contexts and
    all-atom generation capabilities including PPI design.
    """
    
    TOOL_NAME = "RFdiffusionAllAtom"
    DEFAULT_ENV = None  # Loaded from config.yaml
    
    def __init__(self, ligand: str, pdb: Union[str, ToolOutput, StandardizedOutput] = "", contigs: str = "", inpaint: str = "",
                 num_designs: int = 1, active_site: bool = False, 
                 steps: int = 200, partial_steps: int = 0, 
                 reproducible: bool = False, design_startnum: int = 1,
                 ppi_design: bool = False, ppi_hotspot_residues: List[str] = None, 
                 ppi_binder_length: int = None, autogenerate_contigs: bool = False, 
                 model_only_neighbors: bool = False, num_recycles: int = 1, 
                 scaffold_guided: bool = False, align_motif: bool = True, 
                 deterministic: bool = False, inpaint_str: str = None, 
                 inpaint_seq: str = None, inpaint_length: int = None, 
                 guiding_potentials: str = None, **kwargs):
        """
        Initialize RFdiffusion-AllAtom configuration.
        
        Args:
            ligand: Ligand identifier (e.g., 'ZIT', 'RFP')
            pdb: Input PDB file (optional template), ToolOutput, or StandardizedOutput
            contigs: Contig specification (e.g., "A1-100,10-20")
            inpaint: Inpainting specification (same format as contigs)
            num_designs: Number of designs to generate
            active_site: Use active site model for small motifs
            steps: Diffusion steps (default 200 for AllAtom)
            partial_steps: Partial diffusion steps
            reproducible: Use deterministic sampling
            design_startnum: Starting number for design numbering (default: 1)
            ppi_design: Enable protein-protein interaction design
            ppi_hotspot_residues: List of hotspot residues for PPI (e.g., ["A116","A150"])
            ppi_binder_length: Length of PPI binder
            autogenerate_contigs: Auto-infer fixed contig segments from input PDB
            model_only_neighbors: Only remodel residues neighboring the scaffold
            num_recycles: Number of diffusion recycles (iterative refinement)
            scaffold_guided: Enforce strict adherence to input scaffold geometry
            align_motif: Pre-align any functional motif before diffusion
            deterministic: Use fixed RNG seeds for reproducible outputs
            inpaint_str: Secondary-structure pattern for inpainting (e.g., "HHHEE")
            inpaint_seq: Sequence pattern for inpainting (e.g., "ACDEFG")
            inpaint_length: Target length for each inpainted region
            guiding_potentials: JSON or filepath specifying custom external potentials
            **kwargs: Additional parameters
        """
        # In common with RFdiffusion
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
        self.design_startnum = design_startnum
        # Specific to RFdiffusion-AllAtom
        self.ligand = ligand
        self.ppi_design = ppi_design
        self.ppi_hotspot_residues = ppi_hotspot_residues or []
        self.ppi_binder_length = ppi_binder_length
        self.autogenerate_contigs = autogenerate_contigs
        self.model_only_neighbors = model_only_neighbors
        self.num_recycles = num_recycles
        self.scaffold_guided = scaffold_guided
        self.align_motif = align_motif
        self.deterministic = deterministic
        self.inpaint_str = inpaint_str
        self.inpaint_seq = inpaint_seq
        self.inpaint_length = inpaint_length
        self.guiding_potentials = guiding_potentials
        
        if contigs and '/' in contigs:
            print("Warning: Character '/' found in contigs. RFdiffusionAllAtom uses ','.")
        
        # Initialize base class
        super().__init__(**kwargs)

        # Initialize file paths (will be set in configure_inputs)
        self._initialize_file_paths()

    def _initialize_file_paths(self):
        """Initialize common file paths used throughout the class."""
        self.input_pdb_file = None
        self.main_table = None
        self.rfd_log_file = None
        self.pipeline_name = None
        
        # Helper script paths
        self.inference_py_file = None
        self.table_py_file = None
    
    def _extract_pipeline_name(self) -> str:
        """Extract pipeline name from output folder structure."""
        folder_parts = self.output_folder.split(os.sep)
        for i, part in enumerate(folder_parts):
            if "RFdiffusion" in part:
                if i > 0:
                    return folder_parts[i-1]
                break
        raise ValueError(f"Could not extract job name from output folder: {self.output_folder}")
    
    def validate_params(self):
        """Validate RFdiffusion-AllAtom parameters."""
        # Base RFdiffusion validation
        if not self.contigs:
            raise ValueError("contigs parameter is required for RFdiffusion-AllAtom")
        
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
        
        # Additional validation for PPI design
        if self.ppi_design and not self.ppi_hotspot_residues:
            raise ValueError("PPI design requires hotspot residues")
        
        if self.ppi_design and self.ppi_binder_length is None:
            raise ValueError("PPI design requires binder length")
        
        if self.num_recycles < 1:
            raise ValueError("num_recycles must be at least 1")
    
    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input files and dependencies."""
        self.folders = pipeline_folders
        self._setup_file_paths()  # Set up all file paths now that we have folders
        
        # Handle PDB input if provided
        if self.pdb_is_tool_output:
            # Tool output - file already exists at pdb_source_file path
            self.input_sources["pdb"] = self.pdb_source_file
        elif self.pdb:
            # String filename - look in PDBs folder
            pdb_temp = self.pdb if self.pdb.endswith(".pdb") else self.pdb + ".pdb"
            pdb_source = os.path.join(pipeline_folders["PDBs"], pdb_temp)

            if os.path.exists(pdb_source):
                # Will be copied in script generation
                self.input_sources["pdb"] = pdb_source
            else:
                raise ValueError(f"PDB file not found: {pdb_source}")
    
    def get_config_display(self) -> List[str]:
        """Get RFdiffusion-AllAtom configuration display lines."""
        config_lines = super().get_config_display()
        
        # Add RFdiffusion base parameters
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
            config_lines.append(f"REPRODUCIBLE: {self.reproducible}")
        
        # Add AllAtom-specific parameters
        if self.ligand:
            config_lines.append(f"LIGAND: {self.ligand}")
        if self.ppi_design:
            config_lines.append(f"PPI DESIGN: {self.ppi_design}")
            config_lines.append(f"PPI HOTSPOTS: {','.join(self.ppi_hotspot_residues)}")
            config_lines.append(f"PPI BINDER LENGTH: {self.ppi_binder_length}")
        if self.autogenerate_contigs:
            config_lines.append(f"AUTOGENERATE CONTIGS: {self.autogenerate_contigs}")
        if self.num_recycles > 1:
            config_lines.append(f"NUM RECYCLES: {self.num_recycles}")
        if self.scaffold_guided:
            config_lines.append(f"SCAFFOLD GUIDED: {self.scaffold_guided}")
        if not self.align_motif:
            config_lines.append(f"ALIGN MOTIF: {self.align_motif}")
        if self.deterministic:
            config_lines.append(f"DETERMINISTIC: {self.deterministic}")
        
        return config_lines

    def generate_script(self, script_path: str) -> str:
        """
        Generate RFdiffusion-AllAtom execution script.
        
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
        script_content += "# RFdiffusion-AllAtom execution script\n"
        script_content += "# Generated by ProteinNotebooks pipeline system\n\n"
        script_content += self.generate_completion_check_header()
        script_content += self.generate_script_run_rfdiffusion()
        script_content += self.generate_script_create_table()
        script_content += self.generate_completion_check_footer()
        
        return script_content
    
    def _setup_file_paths(self):
        """Set up all file paths after output_folder is known."""
        # Extract pipeline name from folder structure
        self.pipeline_name = self._extract_pipeline_name()
        
        # Core output files
        self.main_table = os.path.join(self.output_folder, "rfdiffusionAA_results.csv")
        
        # Log file is created by pipeline in Logs folder with pattern NNN_{toolname}.log
        # Extract index from folder name (e.g., "001_RFdiffusionAllAtom" -> "001")
        folder_name = os.path.basename(self.output_folder)
        pipeline_folder = os.path.dirname(self.output_folder)  # Get parent folder (DeNovoProtein_013)
        logs_folder = os.path.join(pipeline_folder, "Logs")

        if '_' in folder_name and folder_name.split('_')[0].isdigit():
            index = folder_name.split('_')[0]
            tool_name = folder_name.split('_', 1)[1].lower()  # Get everything after first underscore, lowercase
            self.rfd_log_file = os.path.join(logs_folder, f"{index}_{tool_name}.log")
        else:
            raise ValueError(f"Invalid output folder naming pattern: {folder_name}. Expected 'NNN_RFdiffusionAllAtom' format.")
        
        # Helper script paths (only set if folders are available)
        if hasattr(self, 'folders') and self.folders:
            # RFdiffusionAllAtom has script in root directory, not scripts/
            self.inference_py_file = "run_inference.py"
            self.table_py_file = os.path.join(self.folders["HelpScripts"], "pipe_rfdiffusion_table.py")
            
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
            self.inference_py_file = "run_inference.py"
            self.table_py_file = None
            self.input_pdb_file = ""
    
    def generate_script_run_rfdiffusion(self) -> str:
        """Generate the RFdiffusion-AllAtom execution part of the script."""
        rfd_job_folder = self.output_folder
        
        # Copy input PDB if provided (only if source and destination are different)
        if self.pdb and "pdb" in self.input_sources:
            source_path = os.path.abspath(self.input_sources["pdb"])
            dest_path = os.path.abspath(self.input_pdb_file)
            if source_path != dest_path:
                shutil.copy(source_path, dest_path)
        
        # Build RFdiffusion-AllAtom options (different format than regular RFdiffusion)
        aa_args = []
        
        # Core parameters
        aa_args.append(f"contigmap.contigs=[\\\'{self.contigs}\\\']")
        aa_args.append("inference.ckpt_path=RFDiffusionAA_paper_weights.pt")
        aa_args.append(f"diffuser.T={self.steps}")
        
        if self.input_pdb_file:
            aa_args.append(f"inference.input_pdb={self.input_pdb_file}")
        else:
            # Explicitly set empty input_pdb to override any default test file
            aa_args.append("inference.input_pdb=null")
        
        aa_args.append(f"inference.num_designs={self.num_designs}")
        aa_args.append(f"inference.output_prefix={os.path.join(rfd_job_folder, f'{self.pipeline_name}')}")
        aa_args.append(f"inference.design_startnum={self.design_startnum}")
        
        # Ligand parameter
        if self.ligand:
            aa_args.append(f"inference.ligand={self.ligand}")
        
        # PPI parameters
        if self.ppi_design:
            aa_args.append(f"inference.ppi_design={self.ppi_design}")
        if self.ppi_hotspot_residues:
            aa_args.append(f"ppi.hotspot_res=[\\\'{','.join(self.ppi_hotspot_residues)}\\\']")
        if self.ppi_binder_length is not None:
            aa_args.append(f"ppi.binderlen={self.ppi_binder_length}")
        
        # Advanced parameters
        if self.autogenerate_contigs:
            aa_args.append(f"inference.autogenerate_contigs={self.autogenerate_contigs}")
        if self.model_only_neighbors:
            aa_args.append(f"inference.model_only_neighbors={self.model_only_neighbors}")
        if self.num_recycles > 1:
            aa_args.append(f"inference.num_recycles={self.num_recycles}")
        if self.scaffold_guided:
            aa_args.append(f"inference.scaffold_guided={self.scaffold_guided}")
        if not self.align_motif:
            aa_args.append(f"inference.align_motif={self.align_motif}")
        if self.deterministic:
            aa_args.append(f"inference.deterministic={self.deterministic}")
        
        # Inpainting parameters
        if self.inpaint:
            aa_args.append(f"contigmap.inpaint_seq=[\\\'{self.inpaint}\\\']")
        if self.inpaint_str:
            aa_args.append(f"contigmap.inpaint_str={self.inpaint_str}")
        if self.inpaint_seq:
            aa_args.append(f"contigmap.inpaint_seq={self.inpaint_seq}")
        if self.inpaint_length is not None:
            aa_args.append(f"contigmap.length={self.inpaint_length}")
        
        # Partial steps
        if self.partial_steps > 0:
            aa_args.append(f"diffuser.partial_T={self.partial_steps}")
        
        # Guiding potentials
        if self.guiding_potentials:
            aa_args.append(f"potentials.guiding_potentials={self.guiding_potentials}")
        
        return f"""echo "Starting RFdiffusion-AllAtom"
echo "Arguments: {' '.join(aa_args)}"
echo "Output folder: {rfd_job_folder}"

# Run RFdiffusion-AllAtom
cd {self.folders["RFdiffusionAllAtom"]}
python {self.inference_py_file} {' '.join(aa_args)}

"""
    
    def generate_script_create_table(self) -> str:
        """Generate the table creation part of the script."""
        rfd_job_folder = self.output_folder
        
        # Design character: '?' for RFdiffusion-AllAtom (different from regular RFdiffusion)
        design_character = "?"
        
        return f"""echo "Creating results table"
# Create main table with id, pdb, fixed, designed columns by parsing RFdiffusion-AllAtom log
python {self.table_py_file} "{rfd_job_folder}" "{self.rfd_log_file}" "{design_character}" "{self.pipeline_name}" {self.num_designs} "{self.main_table}" {self.design_startnum}

"""
    
    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after RFdiffusion-AllAtom execution.
        
        Uses pure path construction - no filesystem access.
        Returns expected paths based on naming patterns.
        
        Returns:
            Dictionary mapping output types to file paths with standard keys:
            - structures: PDB files
            - compounds: Empty (no compounds from RFdiffusion-AllAtom)
            - sequences: Empty (no sequences from RFdiffusion-AllAtom)  
            - tables: Main table CSV
            - output_folder: Tool's output directory
        """
        # Ensure file paths are set up
        if not hasattr(self, 'pipeline_name') or self.pipeline_name is None:
            # Fallback if configure_inputs hasn't been called yet
            self._setup_file_paths()
        
        pipeline_name = self.pipeline_name
        main_table = self.main_table
        
        # Generate expected PDB file paths based on clean naming pattern
        design_pdbs = []
        structure_ids = []
        for i in range(self.num_designs):
            design_id = f"{pipeline_name}_{self.design_startnum + i}"
            design_path = os.path.join(self.output_folder, f"{design_id}.pdb")
            design_pdbs.append(design_path)
            structure_ids.append(design_id)
        
        # Import TableInfo
        from .base_config import TableInfo

        # Organize tables by content type
        tables = {
            "structures": TableInfo(
                name="structures",
                path=main_table,
                columns=["id", "source_id", "pdb", "fixed", "designed", "contigs", "time", "status"],
                description="RFdiffusion-AllAtom structure generation results with fixed/designed regions",
                count=self.num_designs
            )
        }
        
        return {
            "structures": design_pdbs,
            "structure_ids": structure_ids,
            "compounds": [],
            "compound_ids": [],
            "sequences": [],
            "sequence_ids": [],
            "tables": tables,
            "output_folder": self.output_folder,
            # Keep legacy aliases for compatibility
            "pdbs": design_pdbs,
            "main": main_table  # Legacy alias for backward compatibility
        }
    
    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration with all RFdiffusion-AllAtom parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                # Core diffusion parameters
                "pdb": self.pdb,
                "contigs": self.contigs,
                "inpaint": self.inpaint,
                "num_designs": self.num_designs,
                "active_site": self.active_site,
                "steps": self.steps,
                "partial_steps": self.partial_steps,
                "reproducible": self.reproducible,
                "design_startnum": self.design_startnum,
                # AllAtom-specific parameters
                "ligand": self.ligand,
                "ppi_design": self.ppi_design,
                "ppi_hotspot_residues": self.ppi_hotspot_residues,
                "ppi_binder_length": self.ppi_binder_length,
                "autogenerate_contigs": self.autogenerate_contigs,
                "model_only_neighbors": self.model_only_neighbors,
                "num_recycles": self.num_recycles,
                "scaffold_guided": self.scaffold_guided,
                "align_motif": self.align_motif,
                "deterministic": self.deterministic,
                "inpaint_str": self.inpaint_str,
                "inpaint_seq": self.inpaint_seq,
                "inpaint_length": self.inpaint_length,
                "guiding_potentials": self.guiding_potentials
            }
        })
        return base_dict