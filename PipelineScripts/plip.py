# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
PLIP (Protein-Ligand Interaction Profiler) configuration for protein-ligand interaction analysis.

Analyzes non-covalent protein-ligand interactions in 3D protein structures using the
PLIP singularity container. Generates interaction profiles, visualizations, and reports.
"""

import os
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
    from .config_manager import ConfigManager
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream
    from config_manager import ConfigManager


class PLIP(BaseConfig):
    """
    PLIP configuration for protein-ligand interaction profiling.

    Uses the PLIP singularity container to analyze protein-ligand interactions,
    generating interaction profiles, XML/text reports, and PyMOL visualizations.
    """

    TOOL_NAME = "PLIP"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba"):
        containers = folders.get("containers", "")
        return f"""echo "=== Installing PLIP ==="
echo "Downloading PLIP singularity container to {containers}"
mkdir -p {containers}
cd {containers}
singularity pull docker://pharmai/plip:3.0.0
mv plip_3.0.0.sif plip_3.0.0.simg

echo "=== PLIP installation complete ==="
"""

    # Lazy path descriptors
    results_csv = Path(lambda self: os.path.join(self.output_folder, f"{self.pipeline_name}_interactions.csv"))
    summary_csv = Path(lambda self: os.path.join(self.output_folder, f"{self.pipeline_name}_summary.csv"))
    summary_txt = Path(lambda self: os.path.join(self.output_folder, f"{self.pipeline_name}_summary.txt"))
    structures_ds_json = Path(lambda self: os.path.join(self.output_folder, "structures.json"))
    raw_outputs_folder = Path(lambda self: os.path.join(self.output_folder, "raw_outputs"))
    processed_folder = Path(lambda self: os.path.join(self.output_folder, "processed"))
    plip_container = Path(lambda self: os.path.join(self.folders["containers"], "plip_3.0.0.simg"))
    helper_script = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_plip_analysis.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
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
            structures: Input structures as DataStream or StandardizedOutput
            ligand: Specific ligand identifier to analyze (empty = analyze all ligands)
            output_format: Output formats to generate ['xml', 'txt', 'pymol', 'images']
            create_pymol: Generate PyMOL session files (.pse)
            create_images: Generate ray-traced images
            analyze_peptides: Include protein-peptide interactions
            analyze_intra: Include intra-chain interactions
            analyze_dna: Include DNA/RNA interactions
            max_threads: Maximum threads for parallel processing
            verbose: Enable verbose output
        """
        # Resolve input to DataStream
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        # PLIP parameters
        self.ligand = ligand
        self.output_format = output_format or ['xml', 'txt']
        self.create_pymol = create_pymol
        self.create_images = create_images
        self.analyze_peptides = analyze_peptides
        self.analyze_intra = analyze_intra
        self.analyze_dna = analyze_dna
        self.max_threads = max_threads
        self.verbose = verbose

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate PLIP-specific parameters."""
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")

        valid_formats = ['xml', 'txt', 'pymol', 'images']
        for fmt in self.output_format:
            if fmt not in valid_formats:
                raise ValueError(f"Invalid output format '{fmt}'. Valid options: {valid_formats}")

        if self.max_threads < 1:
            raise ValueError("max_threads must be >= 1")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input parameters."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get PLIP configuration display lines."""
        config_lines = super().get_config_display()
        config_lines.extend([
            f"STRUCTURES: {len(self.structures_stream)} structures",
            f"LIGAND: {self.ligand or 'All ligands'}",
            f"OUTPUT FORMATS: {', '.join(self.output_format)}",
            f"CREATE PYMOL: {self.create_pymol}",
            f"CREATE IMAGES: {self.create_images}",
            f"ANALYZE PEPTIDES: {self.analyze_peptides}",
            f"ANALYZE INTRA: {self.analyze_intra}",
            f"ANALYZE DNA: {self.analyze_dna}",
            f"MAX THREADS: {self.max_threads}"
        ])
        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate PLIP execution script."""
        import json
        # Serialize structures DataStream to JSON for HelpScript to load
        os.makedirs(self.output_folder, exist_ok=True)
        with open(self.structures_ds_json, 'w') as f:
            json.dump(self.structures_stream.to_dict(), f, indent=2)

        script_content = "#!/bin/bash\n"
        script_content += "# PLIP execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self._generate_script_run_plip()
        script_content += self._generate_script_process_outputs()
        script_content += self.generate_completion_check_footer()

        return script_content

    def _generate_script_run_plip(self) -> str:
        """Generate the PLIP execution part of the script."""
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
        container_exec = ConfigManager().get_container_executor()

        return f"""echo "Running PLIP protein-ligand interaction profiler"
echo "Processing {len(self.structures_stream)} structure(s)"

# Create output directory structure
mkdir -p {self.raw_outputs_folder}
mkdir -p {self.processed_folder}

# Process each structure with PLIP using Python to read DataStream
python -c "
import json
import sys
sys.path.append('{self.folders['HelpScripts']}')
from pipe_biopipelines_io import load_datastream, iterate_files

ds = load_datastream('{self.structures_ds_json}')
for struct_id, struct_path in iterate_files(ds):
    print(struct_id + '\\t' + struct_path)
" | while IFS=$'\\t' read -r struct_id pdb_file; do
    echo "Analyzing structure: $struct_id"

    # Create individual output directory
    output_dir="{self.raw_outputs_folder}/$struct_id"
    mkdir -p "$output_dir"

    # Run PLIP container
    {container_exec} exec {self.plip_container} plip -f "$pdb_file" {plip_opts_str} --outdir "$output_dir"

    if [ $? -eq 0 ]; then
        echo "PLIP analysis completed for $struct_id"
    else
        echo "Error: PLIP analysis failed for $struct_id"
        exit 1
    fi
done

echo "All PLIP analyses completed successfully"

"""

    def _generate_script_process_outputs(self) -> str:
        """Generate the output processing part of the script."""
        ligand_param = f'"{self.ligand}"' if self.ligand else '""'

        return f"""echo "Processing PLIP outputs into standardized format"
python {self.helper_script} \\
    --structures "{self.structures_ds_json}" \\
    --raw_dir "{self.raw_outputs_folder}" \\
    --output_csv "{self.results_csv}" \\
    --summary_csv "{self.summary_csv}" \\
    --summary_txt "{self.summary_txt}" \\
    --ligand {ligand_param} \\
    --processed_dir "{self.processed_folder}"

echo "PLIP output processing completed"

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after PLIP execution."""
        structure_ids = self.structures_stream.ids

        tables = {
            "interactions": TableInfo(
                name="interactions",
                path=self.results_csv,
                columns=["id", "ligand_id", "interaction_type", "residue", "distance", "angle", "energy"],
                description="Protein-ligand interaction analysis results from PLIP",
                count=len(structure_ids)
            ),
            "summary": TableInfo(
                name="summary",
                path=self.summary_csv,
                columns=["id", "structure", "hbonds", "saltbridges", "hydrophobic", "pistacking", "pication", "halogen", "metal", "total_interactions"],
                description="Aggregated interaction counts per structure",
                count=len(structure_ids)
            )
        }

        return {
            "structures": DataStream.empty("structures", "pdb"),
            "sequences": DataStream.empty("sequences", "fasta"),
            "compounds": DataStream.empty("compounds", "sdf"),
            "tables": tables,
            "output_folder": self.output_folder
        }

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
                "verbose": self.verbose
            }
        })
        return base_dict
