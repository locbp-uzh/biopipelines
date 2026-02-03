"""
RFdiffusion configuration for protein backbone generation.
"""

import os
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream, create_map_table
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream, create_map_table


class RFdiffusion(BaseConfig):
    """
    Configuration for RFdiffusion protein backbone generation.
    """

    TOOL_NAME = "RFdiffusion"

    # Lazy path descriptors
    main_table = Path(lambda self: os.path.join(self.output_folder, "rfdiffusion_results.csv"))
    table_py_file = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_rfdiffusion_table.py"))
    inference_py_file = Path(lambda self: os.path.join(self.folders["RFdiffusion"], "scripts", "run_inference.py"))

    def __init__(self,
                 contigs: str,
                 pdb: Optional[Union[DataStream, StandardizedOutput]] = None,
                 inpaint: str = "",
                 num_designs: int = 1,
                 active_site: bool = False,
                 steps: int = 50,
                 partial_steps: int = 0,
                 reproducible: bool = False,
                 design_startnum: int = 1,
                 **kwargs):
        """
        Initialize RFdiffusion configuration.

        Args:
            contigs: Contig specification (e.g., "A1-100,10-20")
            pdb: Optional input structure as DataStream or StandardizedOutput
            inpaint: Inpainting specification (same format as contigs)
            num_designs: Number of designs to generate
            active_site: Use active site model for small motifs
            steps: Diffusion steps (default 50)
            partial_steps: Partial diffusion steps
            reproducible: Use deterministic sampling
            design_startnum: Starting number for design numbering
        """
        # Resolve optional pdb input
        self.pdb_file: Optional[str] = None
        if pdb is not None:
            if isinstance(pdb, StandardizedOutput):
                self.pdb_file = pdb.streams.structures.files[0]
            elif isinstance(pdb, DataStream):
                self.pdb_file = pdb.files[0]
            else:
                raise ValueError(f"pdb must be DataStream or StandardizedOutput, got {type(pdb)}")

        self.contigs = contigs
        self.inpaint = inpaint
        self.num_designs = num_designs
        self.active_site = active_site
        self.steps = steps
        self.partial_steps = partial_steps
        self.reproducible = reproducible
        self.design_startnum = design_startnum

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate RFdiffusion-specific parameters."""
        if not self.contigs:
            raise ValueError("contigs parameter is required")

        if self.num_designs <= 0:
            raise ValueError("num_designs must be positive")

        if self.steps <= 0:
            raise ValueError("steps must be positive")

        if self.partial_steps < 0:
            raise ValueError("partial_steps cannot be negative")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input files."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get RFdiffusion configuration display lines."""
        config_lines = super().get_config_display()
        config_lines.extend([
            f"CONTIGS: {self.contigs}",
            f"NUM DESIGNS: {self.num_designs}",
            f"ACTIVE SITE: {self.active_site}",
            f"STEPS: {self.steps}"
        ])

        if self.pdb_file:
            config_lines.append(f"PDB: {self.pdb_file}")
        if self.inpaint:
            config_lines.append(f"INPAINT: {self.inpaint}")
        if self.partial_steps > 0:
            config_lines.append(f"PARTIAL STEPS: {self.partial_steps}")
        if self.reproducible:
            config_lines.append(f"REPRODUCIBLE: {self.reproducible}")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate RFdiffusion execution script."""
        script_content = "#!/bin/bash\n"
        script_content += "# RFdiffusion execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self._generate_script_run_rfdiffusion()
        script_content += self._generate_script_create_table()
        script_content += self.generate_completion_check_footer()
        return script_content

    def _generate_script_run_rfdiffusion(self) -> str:
        """Generate the RFdiffusion execution part of the script."""
        rfd_options = f"'contigmap.contigs=[{self.contigs}]'"

        if self.inpaint:
            rfd_options += f" 'contigmap.inpaint_seq=[{self.inpaint}]'"

        if self.pdb_file:
            rfd_options += f" inference.input_pdb={self.pdb_file}"

        prefix = os.path.join(self.output_folder, self.pipeline_name)
        rfd_options += f" inference.output_prefix={prefix}"
        rfd_options += f" inference.num_designs={self.num_designs}"
        rfd_options += f" inference.deterministic={self.reproducible}"
        rfd_options += f" inference.design_startnum={self.design_startnum}"

        if self.steps != 50:
            rfd_options += f" diffuser.T={self.steps}"

        if self.partial_steps > 0:
            rfd_options += f" diffuser.partial_T={self.partial_steps}"

        if self.active_site:
            rfd_options += " inference.ckpt_override_path=models/ActiveSite_ckpt.pt"

        return f"""echo "Starting RFdiffusion"
echo "Options: {rfd_options}"
echo "Output folder: {self.output_folder}"

cd {self.folders["RFdiffusion"]}
python {self.inference_py_file} {rfd_options}

"""

    def _generate_script_create_table(self) -> str:
        """Generate the table creation part of the script."""
        design_character = "-"

        return f"""echo "Creating results table"
python {self.table_py_file} "{self.output_folder}" "{self.log_file}" "{design_character}" "{self.pipeline_name}" {self.num_designs} "{self.main_table}" {self.design_startnum}

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after RFdiffusion execution."""
        # Generate expected PDB file paths and IDs
        design_pdbs = []
        structure_ids = []
        for i in range(self.num_designs):
            design_id = f"{self.pipeline_name}_{self.design_startnum + i}"
            design_path = os.path.join(self.output_folder, f"{design_id}.pdb")
            design_pdbs.append(design_path)
            structure_ids.append(design_id)

        # Create map_table for structures
        structures_map = os.path.join(self.output_folder, "structures_map.csv")
        create_map_table(structures_map, structure_ids, files=design_pdbs)

        structures = DataStream(
            name="structures",
            ids=structure_ids,
            files=design_pdbs,
            map_table=structures_map,
            format="pdb"
        )

        tables = {
            "structures": TableInfo(
                name="structures",
                path=self.main_table,
                columns=["id", "source_id", "pdb", "fixed", "designed", "contigs", "time", "status"],
                description="RFdiffusion structure generation results",
                count=self.num_designs
            )
        }

        return {
            "structures": structures,
            "sequences": DataStream.empty("sequences", "fasta"),
            "compounds": DataStream.empty("compounds", "sdf"),
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "rfd_params": {
                "pdb_file": self.pdb_file,
                "contigs": self.contigs,
                "inpaint": self.inpaint,
                "num_designs": self.num_designs,
                "active_site": self.active_site,
                "steps": self.steps,
                "partial_steps": self.partial_steps,
                "reproducible": self.reproducible,
                "design_startnum": self.design_startnum
            }
        })
        return base_dict
