"""
PLIP (Protein-Ligand Interaction Profiler) configuration for protein-ligand interaction analysis.

Analyzes non-covalent protein-ligand interactions in 3D protein structures using the
PLIP singularity container. Generates interaction profiles, visualizations, and reports.
"""

import os
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo


class PLIP(BaseConfig):
    """
    PLIP configuration for protein-ligand interaction profiling.

    Uses the PLIP singularity container to analyze protein-ligand interactions,
    generating interaction profiles, XML/text reports, and PyMOL visualizations.
    """

    TOOL_NAME = "PLIP"
    DEFAULT_ENV = None  # Loaded from config.yaml

    def __init__(self,
                 structures: Union[str, List[str], ToolOutput, StandardizedOutput],
                 ligand: str = "",
                 output_format: List[str] = None,
                 create_pymol: bool = True,
                 create_images: bool = False,
                 analyze_peptides: bool = False,
                 analyze_intra: bool = False,
                 analyze_dna: bool = False,
                 max_threads: int = 4,
                 verbose: bool = True,
                 **kwargs):
        """
        Initialize PLIP configuration.

        Args:
            structures: Input structures (PDB files or ToolOutput from previous tool)
            ligand: Specific ligand identifier to analyze (empty = analyze all ligands)
            output_format: Output formats to generate ['xml', 'txt', 'pymol', 'images']
            create_pymol: Generate PyMOL session files (.pse)
            create_images: Generate ray-traced images
            analyze_peptides: Include protein-peptide interactions
            analyze_intra: Include intra-chain interactions
            analyze_dna: Include DNA/RNA interactions
            max_threads: Maximum threads for parallel processing
            verbose: Enable verbose output
            **kwargs: Additional parameters
        """
        # Store input parameters
        self.input_structures = structures
        self.ligand = ligand
        self.output_format = output_format or ['xml', 'txt']
        self.create_pymol = create_pymol
        self.create_images = create_images
        self.analyze_peptides = analyze_peptides
        self.analyze_intra = analyze_intra
        self.analyze_dna = analyze_dna
        self.max_threads = max_threads
        self.verbose = verbose

        # Determine input type
        self.input_is_tool_output = isinstance(structures, ToolOutput)

        # Initialize base class
        super().__init__(**kwargs)

        # Initialize file paths (will be set in configure_inputs)
        self._initialize_file_paths()

    def validate_params(self):
        """Validate PLIP-specific parameters."""
        if not self.input_structures:
            raise ValueError("structures parameter is required")

        valid_formats = ['xml', 'txt', 'pymol', 'images']
        for fmt in self.output_format:
            if fmt not in valid_formats:
                raise ValueError(f"Invalid output format '{fmt}'. Valid options: {valid_formats}")

        if self.max_threads < 1:
            raise ValueError("max_threads must be >= 1")

    def _initialize_file_paths(self):
        """Initialize common file paths used throughout the class."""
        self.results_csv = None
        self.summary_csv = None
        self.summary_txt = None
        self.plip_container = None
        self.helper_script = None

    def _setup_file_paths(self):
        """Set up all file paths after output_folder is known."""
        # Extract job name for file naming
        job_base = self._extract_job_name()

        # Core output files
        self.results_csv = os.path.join(self.output_folder, f"{job_base}_interactions.csv")
        self.summary_csv = os.path.join(self.output_folder, f"{job_base}_summary.csv")
        self.summary_txt = os.path.join(self.output_folder, f"{job_base}_summary.txt")

        # Tool paths (only set if folders are available)
        if hasattr(self, 'folders') and self.folders:
            self.plip_container = os.path.join(self.folders["containers"], "plip_3.0.0.simg")
            self.helper_script = os.path.join(self.folders["HelpScripts"], "pipe_plip_analysis.py")
        else:
            # Temporary placeholders when folders aren't available yet
            self.plip_container = None
            self.helper_script = None

    def _extract_job_name(self) -> str:
        """Extract job name from output folder structure."""
        folder_parts = self.output_folder.split(os.sep)
        for i, part in enumerate(folder_parts):
            if "PLIP" in part:
                if i > 0:
                    return folder_parts[i-1]
                break

        # Fallback
        return "plip"

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input structures from various sources."""
        self.folders = pipeline_folders
        self._setup_file_paths()

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

        elif hasattr(self.input_structures, 'structures'):
            # StandardizedOutput object (from pipeline.add)
            if self.input_structures.structures:
                self.input_sources["structures"] = self.input_structures.structures
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
        """Get PLIP configuration display lines."""
        config_lines = super().get_config_display()

        # Input information
        if self.input_is_tool_output:
            structure_count = len(getattr(self.input_structures, 'structures', []))
            config_lines.append(f"INPUT: {self.input_structures.tool_type} output ({structure_count} structures)")
        else:
            config_lines.append(f"INPUT: {self.input_structures}")

        config_lines.extend([
            f"LIGAND: {self.ligand or 'All ligands'}",
            f"OUTPUT FORMATS: {', '.join(self.output_format)}",
            f"CREATE PYMOL: {self.create_pymol}",
            f"CREATE IMAGES: {self.create_images}",
            f"ANALYZE PEPTIDES: {self.analyze_peptides}",
            f"ANALYZE INTRA: {self.analyze_intra}",
            f"ANALYZE DNA: {self.analyze_dna}",
            f"MAX THREADS: {self.max_threads}",
            f"VERBOSE: {self.verbose}"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """
        Generate PLIP execution script.

        Args:
            script_path: Path where script should be written

        Returns:
            Script content as string
        """
        script_content = "#!/bin/bash\n"
        script_content += "# PLIP execution script\n"
        script_content += "# Generated by BioPipelines pipeline system\n\n"
        script_content += self.generate_completion_check_header()
        script_content += self.generate_script_run_plip()
        script_content += self.generate_script_process_outputs()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_plip(self) -> str:
        """Generate the PLIP execution part of the script."""
        if "structures" not in self.input_sources:
            raise ValueError("No structure sources found")

        structure_files = self.input_sources["structures"]
        structure_files_str = ",".join(structure_files)

        # Build PLIP command options
        plip_options = []

        if self.verbose:
            plip_options.append("-v")

        if 'xml' in self.output_format:
            plip_options.append("-x")

        if 'txt' in self.output_format:
            plip_options.append("-t")

        if self.create_pymol or 'pymol' in self.output_format:
            plip_options.append("-y")

        if self.create_images or 'images' in self.output_format:
            plip_options.append("-p")

        if self.analyze_peptides:
            plip_options.append("--peptides")

        if self.analyze_intra:
            plip_options.append("--intra")

        if self.analyze_dna:
            plip_options.append("--dnareceptor")

        if self.max_threads > 1:
            plip_options.append(f"--maxthreads {self.max_threads}")

        plip_opts_str = " ".join(plip_options)

        return f"""echo "Running PLIP protein-ligand interaction profiler"
echo "Processing {len(structure_files)} structure(s)"

# Create output directory structure
mkdir -p {self.output_folder}/raw_outputs
mkdir -p {self.output_folder}/processed

# Process each structure with PLIP
IFS=',' read -ra PDB_FILES <<< "{structure_files_str}"
for pdb_file in "${{PDB_FILES[@]}}"; do
    pdb_name=$(basename "$pdb_file" .pdb)
    echo "Analyzing structure: $pdb_name"

    # Create individual output directory
    output_dir="{self.output_folder}/raw_outputs/$pdb_name"
    mkdir -p "$output_dir"

    # Run PLIP apptainer container
    apptainer exec {self.plip_container} plip -f "$pdb_file" {plip_opts_str} --outdir "$output_dir"

    if [ $? -eq 0 ]; then
        echo "PLIP analysis completed for $pdb_name"
    else
        echo "Error: PLIP analysis failed for $pdb_name"
        exit 1
    fi
done

echo "All PLIP analyses completed successfully"

"""

    def generate_script_process_outputs(self) -> str:
        """Generate the output processing part of the script."""
        structure_files = self.input_sources["structures"]
        structure_files_str = ",".join(structure_files)

        ligand_param = f'"{self.ligand}"' if self.ligand else '""'

        return f"""echo "Processing PLIP outputs into standardized format"
python {self.helper_script} \\
    --structures "{structure_files_str}" \\
    --raw_dir "{self.output_folder}/raw_outputs" \\
    --output_csv "{self.results_csv}" \\
    --summary_csv "{self.summary_csv}" \\
    --summary_txt "{self.summary_txt}" \\
    --ligand {ligand_param} \\
    --processed_dir "{self.output_folder}/processed"

echo "PLIP output processing completed"

"""

    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after PLIP execution.

        Returns:
            Dictionary mapping output types to file paths with standard keys
        """
        if not hasattr(self, 'results_csv') or self.results_csv is None:
            self._setup_file_paths()

        # Predict structure IDs for output tracking
        structure_ids = self._predict_structure_ids()

        # PLIP doesn't produce structures, compounds, or sequences
        output_files = {
            "structures": [],
            "structure_ids": [],
            "compounds": [],
            "compound_ids": [],
            "sequences": [],
            "sequence_ids": [],
            "tables": {
                "interactions": TableInfo(
                    name="interactions",
                    path=self.results_csv,
                    columns=["id", "ligand_id", "interaction_type", "residue", "distance", "angle", "energy"],
                    description="Protein-ligand interaction analysis results from PLIP",
                    count=len(structure_ids)  # Approximate
                ),
                "summary": TableInfo(
                    name="summary",
                    path=self.summary_csv,
                    columns=["id", "structure", "hbonds", "saltbridges", "hydrophobic", "pistacking", "pication", "halogen", "metal", "total_interactions"],
                    description="Aggregated interaction counts per structure",
                    count=len(structure_ids)
                )
            },
            "output_folder": self.output_folder,
            # Additional PLIP-specific outputs
            "interactions_csv": [self.results_csv],
            "summary_csv": [self.summary_csv],
            "summary_txt": [self.summary_txt],
            "raw_outputs": [os.path.join(self.output_folder, "raw_outputs")],
            "processed": [os.path.join(self.output_folder, "processed")]
        }

        return output_files

    def _predict_structure_ids(self) -> List[str]:
        """Predict structure IDs from input sources."""
        structure_ids = []

        if hasattr(self, 'input_sources') and "structures" in self.input_sources:
            for pdb_path in self.input_sources["structures"]:
                pdb_base = os.path.splitext(os.path.basename(pdb_path))[0]
                structure_ids.append(pdb_base)

        return structure_ids

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including PLIP-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "plip_params": {
                "ligand": self.ligand,
                "output_format": self.output_format,
                "create_pymol": self.create_pymol,
                "create_images": self.create_images,
                "analyze_peptides": self.analyze_peptides,
                "analyze_intra": self.analyze_intra,
                "analyze_dna": self.analyze_dna,
                "max_threads": self.max_threads,
                "verbose": self.verbose,
                "input_type": "tool_output" if self.input_is_tool_output else "direct"
            }
        })
        return base_dict