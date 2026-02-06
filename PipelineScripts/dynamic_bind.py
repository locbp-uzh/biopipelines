"""
DynamicBind configuration for ligand-specific protein-ligand complex structure prediction.

Predicts protein-ligand complex structures by recovering ligand-specific conformations
from unbound protein structures using diffusion-based generative models.
"""

import os
import pandas as pd
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream


class DynamicBind(BaseConfig):
    """
    DynamicBind for ligand-specific protein-ligand complex structure prediction.

    Predicts protein-ligand binding conformations using equivariant diffusion models
    that recover ligand-specific protein conformations from unbound structures.
    """

    TOOL_NAME = "DynamicBind"

    # Lazy path descriptors
    main_table = Path(lambda self: os.path.join(self.output_folder, "dynamicbind_results.csv"))
    output_compounds_table = Path(lambda self: os.path.join(self.output_folder, "compounds.csv"))
    prepare_ligands_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_dynamicbind_prepare_ligands.py"))
    table_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_dynamicbind_table.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 compounds: Union[DataStream, StandardizedOutput],
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
            structures: Input protein structures as DataStream or StandardizedOutput
            compounds: Input ligands as DataStream or StandardizedOutput
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
        """
        # Resolve structures input to DataStream
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        # Resolve compounds input to DataStream
        if isinstance(compounds, StandardizedOutput):
            self.compounds_stream: DataStream = compounds.streams.compounds
        elif isinstance(compounds, DataStream):
            self.compounds_stream = compounds
        else:
            raise ValueError(f"compounds must be DataStream or StandardizedOutput, got {type(compounds)}")

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

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate DynamicBind parameters."""
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")

        if not self.compounds_stream or len(self.compounds_stream) == 0:
            raise ValueError("compounds parameter is required and must not be empty")

        if self.samples_per_complex <= 0:
            raise ValueError("samples_per_complex must be positive")

        if self.savings_per_complex <= 0:
            raise ValueError("savings_per_complex must be positive")

        if self.inference_steps <= 0:
            raise ValueError("inference_steps must be positive")

        if self.num_workers <= 0:
            raise ValueError("num_workers must be positive")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input files and dependencies."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get DynamicBind configuration display lines."""
        config_lines = super().get_config_display()
        config_lines.extend([
            f"STRUCTURES: {len(self.structures_stream)} proteins",
            f"COMPOUNDS: {len(self.compounds_stream)} ligands",
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
        """Generate DynamicBind execution script."""
        script_content = "#!/bin/bash\n"
        script_content += "# DynamicBind execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self._generate_script_run_dynamicbind()
        script_content += self._generate_script_create_table()
        script_content += self.generate_completion_check_footer()

        return script_content

    def _generate_script_run_dynamicbind(self) -> str:
        """Generate the DynamicBind execution part of the script."""
        # Build DynamicBind command arguments
        args = []
        args.append(f"--samples_per_complex {self.samples_per_complex}")
        args.append(f"--savings_per_complex {self.savings_per_complex}")
        args.append(f"--inference_steps {self.inference_steps}")
        args.append(f"--results {self.output_folder}")
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

        # Get ligands source file
        ligands_source = self.compounds_stream.map_table

        script = f"""echo "Starting DynamicBind"
echo "Structures: {len(self.structures_stream.files)} proteins"
echo "Ligands source: {ligands_source}"
echo "Output folder: {self.output_folder}"

# Create compounds.csv with 'ligand' column from input
python {self.prepare_ligands_py} "{ligands_source}" "{self.output_compounds_table}"

# Get Python paths from both environments
"""
        script += self.activate_environment(index=0)
        script += """DYNAMICBIND_PYTHON=$(which python)
echo "DynamicBind Python: $DYNAMICBIND_PYTHON"

"""
        script += self.activate_environment(index=1)
        script += f"""RELAX_PYTHON=$(which python)
echo "Relax Python: $RELAX_PYTHON"

# Run DynamicBind for each protein
cd {self.folders["DynamicBind"]}
"""

        # Generate commands for each protein
        for i, (protein_id, protein_file) in enumerate(self.structures_stream):
            protein_args = args.copy()
            protein_args.append(f"--header {protein_id}")
            script += f"""
echo "Processing protein {i+1}/{len(self.structures_stream)}: {protein_id}"
"""
            script += self.activate_environment(index=0)
            script += f"""python run_single_protein_inference.py {protein_file} {self.output_compounds_table} {' '.join(protein_args)} --python $DYNAMICBIND_PYTHON --relax_python $RELAX_PYTHON
"""

        return script + "\n"

    def _generate_script_create_table(self) -> str:
        """Generate the table creation part of the script."""
        return f"""echo "Creating results table"

# Parse DynamicBind outputs into standardized table
python {self.table_py} "{self.output_folder}" "{self.main_table}"

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after DynamicBind execution."""
        # DynamicBind generates structures but exact files depend on scoring
        # We return empty structures DataStream as the actual files are in the table

        tables = {
            "structures": TableInfo(
                name="structures",
                path=self.main_table,
                columns=["id", "ligand_id", "structure", "affinity", "lddt", "rank"],
                description="DynamicBind predicted protein-ligand complex structures with affinity and confidence scores",
                count=len(self.structures_stream) * len(self.compounds_stream) * self.samples_per_complex
            ),
            "compounds": TableInfo(
                name="compounds",
                path=self.output_compounds_table,
                columns=["id", "code", "format", "smiles", "ligand", "ccd"],
                description="Compounds used for DynamicBind predictions",
                count=len(self.compounds_stream)
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
        """Serialize configuration with all DynamicBind parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "dynamicbind_params": {
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
