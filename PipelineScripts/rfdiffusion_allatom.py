# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
RFdiffusion-AllAtom configuration for ligand-aware protein design.

Handles RFdiffusion-AllAtom workflows with ligand contexts, PPI design,
and all-atom generation capabilities.
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


class RFdiffusionAllAtom(BaseConfig):
    """
    RFdiffusion-AllAtom variant for ligand-aware protein design.

    Extends base RFdiffusion with support for ligand contexts and
    all-atom generation capabilities including PPI design.
    """

    TOOL_NAME = "RFdiffusionAllAtom"

    # Lazy path descriptors
    main_table = Path(lambda self: os.path.join(self.output_folder, "rfdiffusionAA_results.csv"))
    inference_py_file = Path(lambda self: "run_inference.py")
    table_py_file = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_rfdiffusion_table.py"))

    def __init__(self,
                 ligand: str,
                 pdb: Optional[Union[DataStream, StandardizedOutput]] = None,
                 contigs: str = "",
                 inpaint: str = "",
                 num_designs: int = 1,
                 active_site: bool = False,
                 steps: int = 200,
                 partial_steps: int = 0,
                 reproducible: bool = False,
                 design_startnum: int = 1,
                 ppi_design: bool = False,
                 ppi_hotspot_residues: List[str] = None,
                 ppi_binder_length: int = None,
                 autogenerate_contigs: bool = False,
                 model_only_neighbors: bool = False,
                 num_recycles: int = 1,
                 scaffold_guided: bool = False,
                 align_motif: bool = True,
                 deterministic: bool = False,
                 inpaint_str: str = None,
                 inpaint_seq: str = None,
                 inpaint_length: int = None,
                 guiding_potentials: str = None,
                 **kwargs):
        """
        Initialize RFdiffusion-AllAtom configuration.

        Args:
            ligand: Ligand identifier (e.g., 'ZIT', 'RFP')
            pdb: Input PDB structure as DataStream or StandardizedOutput
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
        # Resolve optional pdb input
        self.pdb_file: Optional[str] = None
        if pdb is not None:
            if isinstance(pdb, StandardizedOutput):
                self.pdb_file = pdb.streams.structures.files[0]
            elif isinstance(pdb, DataStream):
                self.pdb_file = pdb.files[0]
            else:
                raise ValueError(f"pdb must be DataStream or StandardizedOutput, got {type(pdb)}")

        # Core parameters
        self.contigs = contigs
        self.inpaint = inpaint
        self.num_designs = num_designs
        self.active_site = active_site
        self.steps = steps
        self.partial_steps = partial_steps
        self.reproducible = reproducible
        self.design_startnum = design_startnum

        # AllAtom-specific parameters
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

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate RFdiffusion-AllAtom parameters."""
        if not self.contigs:
            raise ValueError("contigs parameter is required for RFdiffusion-AllAtom")

        if self.num_designs <= 0:
            raise ValueError("num_designs must be positive")

        if self.steps <= 0:
            raise ValueError("steps must be positive")

        if self.partial_steps < 0:
            raise ValueError("partial_steps cannot be negative")

        if self.ppi_design and not self.ppi_hotspot_residues:
            raise ValueError("PPI design requires hotspot residues")

        if self.ppi_design and self.ppi_binder_length is None:
            raise ValueError("PPI design requires binder length")

        if self.num_recycles < 1:
            raise ValueError("num_recycles must be at least 1")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input files."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get RFdiffusion-AllAtom configuration display lines."""
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
        """Generate RFdiffusion-AllAtom execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        script_content = "#!/bin/bash\n"
        script_content += "# RFdiffusion-AllAtom execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self._generate_script_run_rfdiffusion()
        script_content += self._generate_script_create_table()
        script_content += self.generate_completion_check_footer()

        return script_content

    def _generate_script_run_rfdiffusion(self) -> str:
        """Generate the RFdiffusion-AllAtom execution part of the script."""
        aa_args = []

        # Core parameters
        aa_args.append(f"contigmap.contigs=[\\\'{self.contigs}\\\']")
        aa_args.append("inference.ckpt_path=RFDiffusionAA_paper_weights.pt")
        aa_args.append(f"diffuser.T={self.steps}")

        if self.pdb_file:
            aa_args.append(f"inference.input_pdb={self.pdb_file}")
        else:
            aa_args.append("inference.input_pdb=null")

        aa_args.append(f"inference.num_designs={self.num_designs}")
        aa_args.append(f"inference.output_prefix={os.path.join(self.output_folder, self.pipeline_name)}")
        aa_args.append(f"inference.design_startnum={self.design_startnum}")

        if self.ligand:
            aa_args.append(f"inference.ligand={self.ligand}")

        if self.ppi_design:
            aa_args.append(f"inference.ppi_design={self.ppi_design}")
        if self.ppi_hotspot_residues:
            aa_args.append(f"ppi.hotspot_res=[\\\'{','.join(self.ppi_hotspot_residues)}\\\']")
        if self.ppi_binder_length is not None:
            aa_args.append(f"ppi.binderlen={self.ppi_binder_length}")

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

        if self.inpaint:
            aa_args.append(f"contigmap.inpaint_seq=[\\\'{self.inpaint}\\\']")
        if self.inpaint_str:
            aa_args.append(f"contigmap.inpaint_str={self.inpaint_str}")
        if self.inpaint_seq:
            aa_args.append(f"contigmap.inpaint_seq={self.inpaint_seq}")
        if self.inpaint_length is not None:
            aa_args.append(f"contigmap.length={self.inpaint_length}")

        if self.partial_steps > 0:
            aa_args.append(f"diffuser.partial_T={self.partial_steps}")

        if self.guiding_potentials:
            aa_args.append(f"potentials.guiding_potentials={self.guiding_potentials}")

        return f"""echo "Starting RFdiffusion-AllAtom"
echo "Arguments: {' '.join(aa_args)}"
echo "Output folder: {self.output_folder}"

cd {self.folders["RFdiffusionAllAtom"]}
python {self.inference_py_file} {' '.join(aa_args)}

"""

    def _generate_script_create_table(self) -> str:
        """Generate the table creation part of the script."""
        design_character = "?"

        return f"""echo "Creating results table"
python {self.table_py_file} "{self.output_folder}" "{self.log_file}" "{design_character}" "{self.pipeline_name}" {self.num_designs} "{self.main_table}" {self.design_startnum}

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after RFdiffusion-AllAtom execution."""
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
                description="RFdiffusion-AllAtom structure generation results with fixed/designed regions",
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
        """Serialize configuration with all RFdiffusion-AllAtom parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "rfdaa_params": {
                "pdb_file": self.pdb_file,
                "contigs": self.contigs,
                "inpaint": self.inpaint,
                "num_designs": self.num_designs,
                "active_site": self.active_site,
                "steps": self.steps,
                "partial_steps": self.partial_steps,
                "reproducible": self.reproducible,
                "design_startnum": self.design_startnum,
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


class RFDAA_PrepareLigand(BaseConfig):
    """
    Preparation tool for RFdiffusion-AllAtom to add a dummy peptide to ligand-only PDB files.
    """

    TOOL_NAME = "RFDAA_PrepareLigand"

    # Lazy path descriptors
    prepared_pdb = Path(lambda self: os.path.join(self.output_folder, "prepared_ligand.pdb"))
    structures_csv = Path(lambda self: os.path.join(self.output_folder, "structures.csv"))
    helper_script = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_rfdaa_prepare_ligand.py"))

    def __init__(self,
                 ligand: Union[DataStream, StandardizedOutput],
                 **kwargs):
        """
        Initialize RFDAA_PrepareLigand tool.

        Args:
            ligand: Ligand structure as DataStream or StandardizedOutput
            **kwargs: Additional parameters
        """
        # Resolve ligand input
        if isinstance(ligand, StandardizedOutput):
            self.ligand_file = ligand.streams.structures.files[0]
        elif isinstance(ligand, DataStream):
            self.ligand_file = ligand.files[0]
        else:
            raise ValueError(f"ligand must be DataStream or StandardizedOutput, got {type(ligand)}")

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate parameters."""
        if not self.ligand_file:
            raise ValueError("ligand file is required")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure inputs."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display."""
        config_lines = super().get_config_display()
        config_lines.append(f"LIGAND_SOURCE: {os.path.basename(self.ligand_file)}")
        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate script to combine ligand with dummy peptide."""
        os.makedirs(self.output_folder, exist_ok=True)

        script_content = "#!/bin/bash\n"
        script_content += "# RFDAA_PrepareLigand execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += f"""
echo "Preparing ligand structure for RFdiffusion-AllAtom"
echo "Input ligand: {self.ligand_file}"
echo "Output: {self.prepared_pdb}"

python "{self.helper_script}" \\
  --ligand_pdb "{self.ligand_file}" \\
  --output_pdb "{self.prepared_pdb}" \\
  --output_csv "{self.structures_csv}" \\
  --pdbs_folder "{self.folders['PDBs']}"

if [ $? -eq 0 ]; then
    echo "Successfully prepared ligand structure"
else
    echo "Error: Failed to prepare ligand structure"
    exit 1
fi

"""
        script_content += self.generate_completion_check_footer()

        return script_content

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files."""
        structure_ids = ["prepared_ligand"]
        structure_files = [self.prepared_pdb]

        structures_map = os.path.join(self.output_folder, "structures_map.csv")
        create_map_table(structures_map, structure_ids, files=structure_files)

        structures = DataStream(
            name="structures",
            ids=structure_ids,
            files=structure_files,
            map_table=structures_map,
            format="pdb"
        )

        tables = {
            "structures": TableInfo(
                name="structures",
                path=self.structures_csv,
                columns=["id", "file_path"],
                description="Prepared ligand structure with dummy peptide",
                count=1
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
            "tool_params": {
                "ligand_source": self.ligand_file
            }
        })
        return base_dict
