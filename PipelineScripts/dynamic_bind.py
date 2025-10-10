"""
DynamicBind configuration for ligand-specific protein-ligand complex structure prediction.

Predicts protein-ligand complex structures by recovering ligand-specific conformations
from unbound protein structures using diffusion-based generative models.
"""

import os
import shutil
import pandas as pd
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo


class DynamicBind(BaseConfig):
    """
    DynamicBind for ligand-specific protein-ligand complex structure prediction.

    Predicts protein-ligand binding conformations using equivariant diffusion models
    that recover ligand-specific protein conformations from unbound structures.
    """

    TOOL_NAME = "DynamicBind"
    DEFAULT_ENV = None  # Loaded from config.yaml

    def __init__(self,
                 proteins: Union[str, List[str], ToolOutput, StandardizedOutput],
                 ligands: Union[str, ToolOutput, StandardizedOutput],
                 samples_per_complex: int = 10,
                 savings_per_complex: int = 10,
                 inference_steps: int = 20,
                 no_relax: bool = False,
                 movie: bool = False,
                 no_clean: bool = False,
                 ligand_is_sdf: bool = False,
                 num_workers: int = 20,
                 paper: bool = False,
                 model: int = 1,
                 seed: int = 42,
                 rigid_protein: bool = False,
                 hts: bool = False,
                 **kwargs):
        """
        Initialize DynamicBind configuration.

        Args:
            proteins: Input protein - PDB filename(s), ToolOutput, or StandardizedOutput
                     If string/list: looks in PDBs folder
                     If ToolOutput/StandardizedOutput: takes all structures
            ligands: Input ligands - SMILES string, ToolOutput, or StandardizedOutput
                    If string: creates CSV with 'ligand' column
                    If ToolOutput/StandardizedOutput: uses compounds datasheet (must have 'smiles' column)
            samples_per_complex: Number of samples generated per complex (default: 10)
            savings_per_complex: Number of samples saved for movie generation (default: 10)
            inference_steps: Number of coordinate update steps (default: 20)
            no_relax: Skip relaxing last frame (default: False)
            movie: Generate movie of binding process (default: False)
            no_clean: Skip cleaning input protein file (default: False)
            ligand_is_sdf: Ligand file in SDF format instead of CSV (default: False)
            num_workers: Workers for relaxing final structure (default: 20)
            paper: Use paper version model (default: False)
            model: Model version to use (default: 1)
            seed: Random seed for reproducibility (default: 42)
            rigid_protein: No noise in reverse diffusion final step (default: False)
            hts: High-throughput screening mode (default: False)
            **kwargs: Additional parameters
        """
        self.proteins = proteins
        self.proteins_is_tool_output = False
        self.protein_source_files = []

        # Handle tool output for proteins input
        if isinstance(proteins, (ToolOutput, StandardizedOutput)):
            self.proteins_is_tool_output = True
            if isinstance(proteins, StandardizedOutput):
                # Get all structures from StandardizedOutput
                if hasattr(proteins, 'structures') and proteins.structures:
                    self.protein_source_files = proteins.structures
                else:
                    raise ValueError("No structures found in StandardizedOutput for proteins parameter")
            else:  # ToolOutput
                # Get all structures from ToolOutput
                structures = proteins.get_output_files("structures")
                if structures:
                    self.protein_source_files = structures
                    # Add dependency
                    self.dependencies.append(proteins.config)
                else:
                    raise ValueError("No structures found in ToolOutput for proteins parameter")
        elif isinstance(proteins, list):
            # List of PDB filenames
            self.protein_source_files = proteins
        else:
            # Single PDB filename
            self.protein_source_files = [proteins]

        self.ligands = ligands
        self.ligands_is_tool_output = False
        self.ligands_is_smiles = False
        self.ligands_source_file = None

        # Handle ligands input
        if isinstance(ligands, str):
            # Direct SMILES string
            self.ligands_is_smiles = True
            self.ligands_source_file = None  # Will create CSV in configure_inputs
        elif isinstance(ligands, (ToolOutput, StandardizedOutput)):
            self.ligands_is_tool_output = True
            if isinstance(ligands, StandardizedOutput):
                # Get compounds CSV from StandardizedOutput
                if hasattr(ligands, 'compounds') and ligands.compounds:
                    # compounds is a list, take the first CSV file
                    self.ligands_source_file = ligands.compounds[0] if isinstance(ligands.compounds, list) else ligands.compounds
                else:
                    raise ValueError("No compounds found in StandardizedOutput for ligands parameter")
            else:  # ToolOutput
                # Get compounds CSV from ToolOutput
                compounds = ligands.get_output_files("compounds")
                if compounds:
                    # compounds is a list, take the first CSV file
                    self.ligands_source_file = compounds[0] if isinstance(compounds, list) else compounds
                    # Add dependency
                    self.dependencies.append(ligands.config)
                else:
                    raise ValueError("No compounds found in ToolOutput for ligands parameter")
        else:
            raise ValueError("ligands must be a SMILES string, ToolOutput, or StandardizedOutput")

        # DynamicBind parameters
        self.samples_per_complex = samples_per_complex
        self.savings_per_complex = savings_per_complex
        self.inference_steps = inference_steps
        self.no_relax = no_relax
        self.movie = movie
        self.no_clean = no_clean
        self.ligand_is_sdf = ligand_is_sdf
        self.num_workers = num_workers
        self.paper = paper
        self.model = model
        self.seed = seed
        self.rigid_protein = rigid_protein
        self.hts = hts

        # Initialize base class
        super().__init__(**kwargs)

        # Initialize file paths (will be set in configure_inputs)
        self._initialize_file_paths()

    def _initialize_file_paths(self):
        """Initialize common file paths used throughout the class."""
        self.input_protein_files = []
        self.input_ligands_file = None
        self.output_compounds_datasheet = None
        self.main_datasheet = None
        self.pipeline_name = None
        self.inference_py_file = None

    def _extract_pipeline_name(self) -> str:
        """Extract pipeline name from output folder structure."""
        folder_parts = self.output_folder.split(os.sep)
        for i, part in enumerate(folder_parts):
            if "DynamicBind" in part:
                if i > 0:
                    return folder_parts[i-1]
                break
        raise ValueError(f"Could not extract job name from output folder: {self.output_folder}")

    def validate_params(self):
        """Validate DynamicBind parameters."""
        if not self.proteins:
            raise ValueError("proteins parameter is required for DynamicBind")

        if not self.ligands:
            raise ValueError("ligands parameter is required for DynamicBind")

        if self.samples_per_complex <= 0:
            raise ValueError("samples_per_complex must be positive")

        if self.savings_per_complex <= 0:
            raise ValueError("savings_per_complex must be positive")

        if self.inference_steps <= 0:
            raise ValueError("inference_steps must be positive")

        if self.num_workers <= 0:
            raise ValueError("num_workers must be positive")

        # Skip file validation for tool outputs (will be validated at runtime)
        if not self.proteins_is_tool_output:
            for protein_file in self.protein_source_files:
                if isinstance(protein_file, str):
                    if not (protein_file.endswith('.pdb') or os.path.exists(protein_file)):
                        # Check if it exists in PDBs folder
                        pdb_path = os.path.join(os.getcwd(), "PDBs", protein_file + ".pdb" if not protein_file.endswith('.pdb') else protein_file)
                        if not os.path.exists(pdb_path):
                            raise ValueError(f"Protein PDB file not found: {protein_file}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input files and dependencies."""
        self.folders = pipeline_folders
        self._setup_file_paths()

        # Handle proteins input
        if self.proteins_is_tool_output:
            # Tool output - files already exist at protein_source_files paths
            self.input_sources["proteins"] = self.protein_source_files
        else:
            # String filename(s) - look in PDBs folder
            resolved_proteins = []
            for protein in self.protein_source_files:
                pdb_temp = protein if protein.endswith(".pdb") else protein + ".pdb"
                pdb_source = os.path.join(pipeline_folders["PDBs"], pdb_temp)

                if os.path.exists(pdb_source):
                    resolved_proteins.append(pdb_source)
                else:
                    raise ValueError(f"Protein PDB file not found: {pdb_source}")
            self.input_sources["proteins"] = resolved_proteins

        # Handle ligands input
        if self.ligands_is_smiles:
            # Direct SMILES string - will create CSV in script generation
            self.input_sources["ligands"] = None  # Marker for SMILES input
        elif self.ligands_is_tool_output:
            # Tool output - file already exists at ligands_source_file path
            self.input_sources["ligands"] = self.ligands_source_file
        else:
            raise ValueError("Invalid ligands input type")

    def _setup_file_paths(self):
        """Set up all file paths after output_folder is known."""
        # Extract pipeline name from folder structure
        self.pipeline_name = self._extract_pipeline_name()

        # Core output files
        self.main_datasheet = os.path.join(self.output_folder, "dynamicbind_results.csv")
        self.output_compounds_datasheet = os.path.join(self.output_folder, "compounds.csv")

        # Helper script paths
        if hasattr(self, 'folders') and self.folders:
            self.inference_py_file = "run_single_protein_inference.py"

            # Input protein file paths - copy to runtime if needed
            if self.proteins_is_tool_output:
                self.input_protein_files = self.protein_source_files
            else:
                self.input_protein_files = []
                for protein in self.protein_source_files:
                    pdb_temp = protein if protein.endswith(".pdb") else protein + ".pdb"
                    self.input_protein_files.append(os.path.join(self.folders["runtime"], pdb_temp))

            # Input ligands file path
            if self.ligands_is_smiles:
                # Will create CSV from SMILES
                self.input_ligands_file = os.path.join(self.folders["runtime"], "ligands.csv")
            elif self.ligands_is_tool_output:
                # Will create CSV with 'ligand' column from tool output
                self.input_ligands_file = os.path.join(self.folders["runtime"], "ligands.csv")
            else:
                self.input_ligands_file = ""
        else:
            # Temporary placeholders when folders aren't available yet
            self.inference_py_file = "run_single_protein_inference.py"
            self.input_protein_files = []
            self.input_ligands_file = ""

    def get_config_display(self) -> List[str]:
        """Get DynamicBind configuration display lines."""
        config_lines = super().get_config_display()

        protein_display = f"{len(self.protein_source_files)} structures" if self.proteins_is_tool_output else str(self.proteins)
        ligand_display = "SMILES" if self.ligands_is_smiles else "Tool output"

        config_lines.extend([
            f"PROTEINS: {protein_display}",
            f"LIGANDS: {ligand_display}",
            f"SAMPLES PER COMPLEX: {self.samples_per_complex}",
            f"SAVINGS PER COMPLEX: {self.savings_per_complex}",
            f"INFERENCE STEPS: {self.inference_steps}",
            f"NO RELAX: {self.no_relax}",
            f"SEED: {self.seed}"
        ])

        if self.movie:
            config_lines.append(f"MOVIE: {self.movie}")
        if self.rigid_protein:
            config_lines.append(f"RIGID PROTEIN: {self.rigid_protein}")
        if self.hts:
            config_lines.append(f"HTS MODE: {self.hts}")
        if self.paper:
            config_lines.append(f"PAPER MODEL: {self.paper}")
        if self.model != 1:
            config_lines.append(f"MODEL VERSION: {self.model}")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """
        Generate DynamicBind execution script.

        Args:
            script_path: Path where script should be written

        Returns:
            Script content as string
        """
        runtime_folder = os.path.dirname(script_path)
        db_job_folder = self.output_folder
        os.makedirs(db_job_folder, exist_ok=True)

        # Generate script content
        script_content = "#!/bin/bash\n"
        script_content += "# DynamicBind execution script\n"
        script_content += "# Generated by BioPipelines system\n\n"
        script_content += self.generate_completion_check_header()
        script_content += self.generate_script_run_dynamicbind()
        script_content += self.generate_script_create_datasheet()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_dynamicbind(self) -> str:
        """Generate the DynamicBind execution part of the script."""
        db_job_folder = self.output_folder

        # Copy protein files if needed
        if "proteins" in self.input_sources:
            for source_path, dest_path in zip(self.input_sources["proteins"], self.input_protein_files):
                source_path = os.path.abspath(source_path)
                dest_path = os.path.abspath(dest_path)
                if source_path != dest_path:
                    shutil.copy(source_path, dest_path)

        # Handle ligands input
        if self.ligands_is_smiles:
            # Create CSV with 'ligand' column from SMILES string
            import pandas as pd
            ligands_df = pd.DataFrame({'ligand': [self.ligands]})
            ligands_df.to_csv(self.input_ligands_file, index=False)

            # Also create standardized compounds.csv
            compounds_df = pd.DataFrame({
                'id': ['ligand'],
                'format': ['smiles'],
                'smiles': [self.ligands],
                'ccd': ['']
            })
            compounds_df.to_csv(self.output_compounds_datasheet, index=False)
        elif self.ligands_is_tool_output and "ligands" in self.input_sources:
            # Read tool output compounds datasheet and create new CSV with 'ligand' column
            # This will be done in the bash script using the helper script
            pass

        # Build DynamicBind command arguments
        args = []
        args.append(f"--samples_per_complex {self.samples_per_complex}")
        args.append(f"--savings_per_complex {self.savings_per_complex}")
        args.append(f"--inference_steps {self.inference_steps}")
        args.append(f"--header {self.pipeline_name}")
        args.append(f"--results {db_job_folder}")
        args.append(f"--num_workers {self.num_workers}")
        args.append(f"--model {self.model}")
        args.append(f"--seed {self.seed}")

        if self.no_relax:
            args.append("--no_relax")
        if self.movie:
            args.append("--movie")
        if self.no_clean:
            args.append("--no_clean")
        if self.ligand_is_sdf:
            args.append("--ligand_is_sdf")
        if self.paper:
            args.append("--paper")
        if self.rigid_protein:
            args.append("--rigid_protein")
        if self.hts:
            args.append("--hts")

        # Generate bash script content
        script = f"""echo "Starting DynamicBind"
echo "Ligands: {self.input_ligands_file}"
echo "Output folder: {db_job_folder}"

"""

        # Add ligands CSV creation if from tool output
        if self.ligands_is_tool_output and "ligands" in self.input_sources:
            script += f"""# Create ligands CSV with 'ligand' column from tool output
python {self.folders["HelpScripts"]}/pipe_dynamicbind_prepare_ligands.py "{self.input_sources["ligands"]}" "{self.input_ligands_file}" "{self.output_compounds_datasheet}"

"""

        script += f"""# Activate dynamicbind environment and get python path
source $(conda info --base)/etc/profile.d/conda.sh
conda activate dynamicbind
DYNAMICBIND_PYTHON=$(which python)
echo "DynamicBind Python: $DYNAMICBIND_PYTHON"

# Activate relax environment and get python path
conda activate relax
RELAX_PYTHON=$(which python)
echo "Relax Python: $RELAX_PYTHON"

# Switch back to dynamicbind for main execution
conda activate dynamicbind

# Run DynamicBind for each protein
cd {self.folders["DynamicBind"]}
"""

        # Generate commands for each protein
        for i, protein_file in enumerate(self.input_protein_files):
            protein_name = os.path.basename(protein_file).replace('.pdb', '')
            script += f"""
echo "Processing protein {i+1}/{len(self.input_protein_files)}: {protein_name}"
python {self.inference_py_file} {protein_file} {self.input_ligands_file} {' '.join(args)} --python $DYNAMICBIND_PYTHON --relax_python $RELAX_PYTHON
"""

        return script + "\n"

    def generate_script_create_datasheet(self) -> str:
        """Generate the datasheet creation part of the script."""
        db_job_folder = self.output_folder

        # DynamicBind outputs affinity_prediction.csv
        # We'll create our standardized datasheet from it
        return f"""echo "Creating results datasheet"

# Parse DynamicBind outputs into standardized datasheet
python {self.folders["HelpScripts"]}/pipe_dynamicbind_datasheet.py "{db_job_folder}" "{self.main_datasheet}"

"""

    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after DynamicBind execution.

        Uses pure path construction - no filesystem access.
        Returns expected paths based on naming patterns.

        Returns:
            Dictionary mapping output types to file paths with standard keys:
            - structures: SDF files of predicted complexes
            - compounds: Empty (compounds are in input)
            - sequences: Empty (no sequences from DynamicBind)
            - datasheets: Main datasheet CSV
            - output_folder: Tool's output directory
        """
        # Ensure file paths are set up
        if not hasattr(self, 'pipeline_name') or self.pipeline_name is None:
            self._setup_file_paths()

        pipeline_name = self.pipeline_name
        main_datasheet = self.main_datasheet

        # DynamicBind generates structures but exact naming depends on scoring
        # We'll return the output folder for now and let the datasheet determine actual files
        structures = []
        structure_ids = []

        # Organize datasheets by content type
        datasheets = {
            "structures": DatasheetInfo(
                name="structures",
                path=main_datasheet,
                columns=["id", "ligand_id", "structure", "affinity", "lddt", "rank"],
                description="DynamicBind predicted protein-ligand complex structures with affinity and confidence scores",
                count=None  # Unknown until execution
            ),
            "compounds": DatasheetInfo(
                name="compounds",
                path=self.output_compounds_datasheet,
                columns=["id", "format", "smiles", "ccd"],  # Standardized format
                description="Compounds used for DynamicBind predictions",
                count=None
            )
        }

        return {
            "structures": structures,
            "structure_ids": structure_ids,
            "compounds": [self.output_compounds_datasheet],
            "compound_ids": [],
            "sequences": [],
            "sequence_ids": [],
            "datasheets": datasheets,
            "output_folder": self.output_folder,
            "main": main_datasheet
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration with all DynamicBind parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                "proteins": self.proteins,
                "ligands": self.ligands,
                "samples_per_complex": self.samples_per_complex,
                "savings_per_complex": self.savings_per_complex,
                "inference_steps": self.inference_steps,
                "no_relax": self.no_relax,
                "movie": self.movie,
                "no_clean": self.no_clean,
                "ligand_is_sdf": self.ligand_is_sdf,
                "num_workers": self.num_workers,
                "paper": self.paper,
                "model": self.model,
                "seed": self.seed,
                "rigid_protein": self.rigid_protein,
                "hts": self.hts
            }
        })
        return base_dict
