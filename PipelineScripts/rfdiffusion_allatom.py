"""
RFdiffusion-AllAtom configuration for ligand-aware protein design.

Handles RFdiffusion-AllAtom workflows with ligand contexts, PPI design,
and all-atom generation capabilities.
"""

import os
import shutil
from typing import Dict, List, Any, Optional, Union

try:
    from .rfdiffusion import RFdiffusion
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from rfdiffusion import RFdiffusion


class RFdiffusionAllAtom(RFdiffusion):
    """
    RFdiffusion-AllAtom variant for ligand-aware protein design.
    
    Extends base RFdiffusion with support for ligand contexts and
    all-atom generation capabilities including PPI design.
    """
    
    TOOL_NAME = "RFdiffusionAllAtom"
    DEFAULT_ENV = "ProteinEnv"
    COMPATIBLE_ENVS = ["ProteinEnv"]
    
    def __init__(self, ligand: str, ppi_design: bool = False,
                 ppi_hotspot_residues: List[str] = None, ppi_binder_length: int = None,
                 autogenerate_contigs: bool = False, model_only_neighbors: bool = False,
                 num_recycles: int = 1, scaffold_guided: bool = False,
                 align_motif: bool = True, deterministic: bool = False,
                 inpaint_str: str = None, inpaint_seq: str = None, 
                 inpaint_length: int = None, guiding_potentials: str = None,
                 **kwargs):
        """
        Initialize RFdiffusion-AllAtom configuration.
        
        Args:
            ligand: Ligand identifier (e.g., 'ZIT', 'RFP')
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
            **kwargs: RFdiffusion base parameters
        """
        # Store RFdiffusion-AllAtom specific parameters
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
        
        if 'contigs' in kwargs:
            if '/' in kwargs['contigs']:
                print("Warning: Character '/' found in contigs. RFdiffusionAllAtom uses ','.")

        # Override default steps for AllAtom (notebook uses denoising_steps=20 for debug, 100-200 normally)
        if 'steps' not in kwargs:
            kwargs['steps'] = 200
        
        super().__init__(**kwargs)
    
    def validate_params(self):
        """Validate RFdiffusion-AllAtom parameters."""
        super().validate_params()
        
        # Additional validation for PPI design
        if self.ppi_design and not self.ppi_hotspot_residues:
            raise ValueError("PPI design requires hotspot residues")
        
        if self.ppi_design and self.ppi_binder_length is None:
            raise ValueError("PPI design requires binder length")
        
        if self.num_recycles < 1:
            raise ValueError("num_recycles must be at least 1")
    
    def get_config_display(self) -> List[str]:
        """Get RFdiffusion-AllAtom configuration display lines."""
        config_lines = super().get_config_display()
        
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
    
    def _setup_file_paths(self):
        """Set up all file paths after output_folder is known."""
        # Extract pipeline name from folder structure
        self.pipeline_name = self._extract_pipeline_name()
        
        # Core output files
        self.main_datasheet = os.path.join(self.output_folder, "rfdiffusion_results.csv")
        
        # Log file is created by pipeline in the main pipeline folder with pattern _{index}_{toolname}.log
        # Extract index from folder name (e.g., "1_RFdiffusionAllAtom" -> "1")
        folder_name = os.path.basename(self.output_folder)
        pipeline_folder = os.path.dirname(self.output_folder)  # Get parent folder (DeNovoProtein_013)
        
        if '_' in folder_name and folder_name.split('_')[0].isdigit():
            index = folder_name.split('_')[0]
            tool_name = folder_name.split('_', 1)[1].lower()  # Get everything after first underscore, lowercase
            self.rfd_log_file = os.path.join(pipeline_folder, f"_{index}_{tool_name}.log")
        else:
            raise ValueError(f"Invalid output folder naming pattern: {folder_name}. Expected 'N_RFdiffusionAllAtom' format.")
        
        # Helper script paths (only set if folders are available)
        if hasattr(self, 'folders') and self.folders:
            # RFdiffusionAllAtom has script in root directory, not scripts/
            self.inference_py_file = "run_inference.py"
            self.datasheet_py_file = os.path.join(self.folders["HelpScripts"], "pipe_rfdiffusion_datasheet.py")
            
            # Input PDB file path (if PDB is provided)
            if self.pdb:
                pdb_temp = self.pdb if self.pdb.endswith(".pdb") else self.pdb + ".pdb"
                self.input_pdb_file = os.path.join(self.folders["runtime"], pdb_temp)
            else:
                self.input_pdb_file = ""
        else:
            # Temporary placeholders when folders aren't available yet
            self.inference_py_file = "run_inference.py"
            self.datasheet_py_file = None
            self.input_pdb_file = ""
    
    def generate_script_run_rfdiffusion(self) -> str:
        """Generate the RFdiffusion-AllAtom execution part of the script."""
        rfd_job_folder = self.output_folder
        
        # Copy input PDB if provided
        if self.pdb and "pdb" in self.input_sources:
            shutil.copy(self.input_sources["pdb"], self.input_pdb_file)
        
        # Build RFdiffusion-AllAtom options (different format than regular RFdiffusion)
        aa_args = []
        
        # Core parameters
        aa_args.append(f"contigmap.contigs=[\\\'{self.contigs}\\\']")
        aa_args.append("inference.ckpt_path=RFDiffusionAA_paper_weights.pt")
        aa_args.append(f"diffuser.T={self.steps}")
        
        if self.input_pdb_file:
            aa_args.append(f"inference.input_pdb={self.input_pdb_file}")
        
        aa_args.append(f"inference.num_designs={self.num_designs}")
        aa_args.append(f"inference.output_prefix={os.path.join(rfd_job_folder, f'{self.pipeline_name}')}")
        aa_args.append(f"inference.design_startnum={self.reproducibility_number}")
        
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
    
    def generate_script_create_datasheet(self) -> str:
        """Generate the datasheet creation part of the script."""
        rfd_job_folder = self.output_folder
        
        # Design character: '?' for RFdiffusion-AllAtom (different from regular RFdiffusion)
        design_character = "?"
        
        return f"""echo "Creating results datasheet"
# Create main datasheet with id, pdb, fixed, designed columns by parsing RFdiffusion-AllAtom log
python {self.datasheet_py_file} "{rfd_job_folder}" "{self.rfd_log_file}" "{design_character}" "{self.pipeline_name}" {self.num_designs} "{self.main_datasheet}"

"""
    
    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including RFdiffusion-AllAtom-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "rfd_aa_params": {
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