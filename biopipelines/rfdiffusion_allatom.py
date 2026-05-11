# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
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
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream, create_map_table
    from .biopipelines_io import Resolve
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream, create_map_table
    from biopipelines_io import Resolve


class RFdiffusionAllAtom(BaseConfig):
    """
    RFdiffusion-AllAtom variant for ligand-aware protein design.

    Extends base RFdiffusion with support for ligand contexts and
    all-atom generation capabilities including PPI design.
    """

    TOOL_NAME = "RFdiffusionAllAtom"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        repo_dir = folders.get("RFdiffusionAllAtom", "")
        parent_dir = os.path.dirname(repo_dir)
        skip = "" if force_reinstall else f"""# Check if already installed
if [ -d "{repo_dir}" ] && [ -f "{repo_dir}/RFDiffusionAA_paper_weights.pt" ]; then
    echo "RFdiffusion-AllAtom already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        env_check = cls._env_exists_check("SE3nv", env_manager)
        return f"""echo "=== Installing RFdiffusion-AllAtom ==="
# SE3nv environment must exist — run RFdiffusion.install() first
if ! {env_check}; then
    echo "ERROR: SE3nv environment not found. Run RFdiffusion.install() before RFdiffusionAllAtom.install()."
    exit 1
fi
{skip}mkdir -p {parent_dir}
cd {parent_dir}
if [ ! -d "{repo_dir}" ]; then
    git clone https://github.com/baker-laboratory/rf_diffusion_all_atom.git
fi
cd {repo_dir}

# Download model weights
wget -nc http://files.ipd.uw.edu/pub/RF-All-Atom/weights/RFDiffusionAA_paper_weights.pt

# Initialize git submodules
git submodule init
git submodule update

# Additional dependencies into the shared SE3nv environment
{env_manager} run -n SE3nv pip install icecream openbabel-wheel assertpy

# Verify installation
if [ -f "{repo_dir}/RFDiffusionAA_paper_weights.pt" ]; then
    touch "$INSTALL_SUCCESS"
    echo "=== RFdiffusion-AllAtom installation complete ==="
    echo "Container mode: configure containers.RFdiffusionAllAtom in config.yaml"
else
    echo "ERROR: RFdiffusion-AllAtom verification failed (weights missing)"
    exit 1
fi
"""

    # Lazy path descriptors
    #   main_table — standalone TableInfo CSV (tables/structures.csv).
    #   pdb_ds_json — config-time input DataStream serialization.
    main_table = Path(lambda self: self.table_path("structures"))
    inference_py_file = Path(lambda self: "run_inference.py")
    table_py_file = Path(lambda self: self.pipe_script_path("pipe_rfdiffusion_table.py"))
    pdb_ds_json = Path(lambda self: self.configuration_path("input_structures.json"))
    update_map_py = Path(lambda self: self.pipe_script_path("pipe_update_structures_map.py"))

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

        Output:
            Streams: structures (.pdb)
            Tables:
                structures: id | pdb | fixed | designed | source_fixed | plddt_mean | status
        """
        # Resolve optional pdb input — store stream for runtime resolution
        self.pdb_stream: Optional[DataStream] = None
        self.pdb_input_id: Optional[str] = None
        if pdb is not None:
            if isinstance(pdb, StandardizedOutput):
                self.pdb_stream = pdb.streams.structures
            elif isinstance(pdb, DataStream):
                self.pdb_stream = pdb
            else:
                raise ValueError(f"pdb must be DataStream or StandardizedOutput, got {type(pdb)}")
            self.pdb_input_id = self.pdb_stream.ids[0]

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

        _validate_freeform_string("ligand", self.ligand)
        _validate_freeform_string("contigs", self.contigs)
        _validate_freeform_string("inpaint", self.inpaint)
        _validate_freeform_string("inpaint_str", self.inpaint_str)
        _validate_freeform_string("inpaint_seq", self.inpaint_seq)
        _validate_freeform_string("guiding_potentials", self.guiding_potentials)
        for i, res in enumerate(self.ppi_hotspot_residues):
            _validate_freeform_string(f"ppi_hotspot_residues[{i}]", res)

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input files."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get RFdiffusion-AllAtom configuration display lines."""
        config_lines = super().get_config_display()

        mode = "container" if self.uses_container() else "environment (SE3nv)"
        config_lines.extend([
            f"MODE: {mode}",
            f"CONTIGS: {self.contigs}",
            f"NUM DESIGNS: {self.num_designs}",
            f"ACTIVE SITE: {self.active_site}",
            f"STEPS: {self.steps}"
        ])

        if self.pdb_stream:
            config_lines.append(f"PDB: {self.pdb_input_id}")
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
        # Serialize input DataStream to JSON for runtime file resolution.
        # configuration/ is auto-created by the pipeline.
        if self.pdb_stream:
            self.pdb_stream.save_json(self.pdb_ds_json)

        script_content = "#!/bin/bash\n"
        script_content += "# RFdiffusion-AllAtom execution script\n"
        script_content += self.generate_completion_check_header()
        # e3nn 0.3.3 uses torch.load() without weights_only=False,
        # which fails on PyTorch 2.6+ where the default flipped to True
        script_content += "export TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1\n" 
        # activate_environment() activates the tool's configured env on the
        # host (SE3nv by default, or biopipelines fallback). Inference runs
        # under container_prefix when a container is set; host-side helpers
        # below run under the activated env either way.
        script_content += self.activate_environment()
        script_content += self._generate_script_run_rfdiffusion()
        script_content += self._generate_script_create_table()
        script_content += self._generate_script_update_structures_map()
        script_content += self.generate_completion_check_footer()

        return script_content

    def _build_inference_args(self) -> List[str]:
        """Build the inference argument list (shared between container and env modes)."""
        aa_args = []

        # Core parameters
        aa_args.append(f"contigmap.contigs=[\\\'{self.contigs}\\\']")
        aa_args.append("inference.ckpt_path=RFDiffusionAA_paper_weights.pt")
        aa_args.append(f"diffuser.T={self.steps}")

        if self.pdb_stream:
            aa_args.append("inference.input_pdb=$INPUT_PDB")
        else:
            aa_args.append("inference.input_pdb=null")

        aa_args.append(f"inference.num_designs={self.num_designs}")
        output_name = self.pdb_input_id if self.pdb_input_id else self.pipeline_name
        # Route PDB outputs into the structures/ stream folder.
        aa_args.append(f"inference.output_prefix={self.stream_path('structures', output_name)}")
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

        return aa_args

    def _generate_script_run_rfdiffusion(self) -> str:
        """Generate the RFdiffusion-AllAtom execution part of the script."""
        # Resolve input PDB at runtime if a DataStream is provided
        resolve_snippet = ""
        if self.pdb_stream:
            resolve_snippet = f"""INPUT_PDB_ID={Resolve.stream_ids(self.pdb_ds_json, index=0)}
INPUT_PDB={Resolve.stream_item(self.pdb_ds_json, '$INPUT_PDB_ID')}
"""

        aa_args = self._build_inference_args()
        args_str = ' '.join(aa_args)
        repo_dir = self.folders["RFdiffusionAllAtom"]
        mode = "container" if self.uses_container() else "environment"

        return f"""{resolve_snippet}echo "Starting RFdiffusion-AllAtom ({mode} mode)"
echo "Arguments: {args_str}"
echo "Output folder: {self.output_folder}"

cd {repo_dir}
{self.container_prefix()}python {self.inference_py_file} {args_str}

"""

    def _generate_script_create_table(self) -> str:
        """Generate the table creation part of the script."""
        output_name = self.pdb_input_id if self.pdb_input_id else self.pipeline_name
        # pipe_rfdiffusion_table.py scans the given folder for .pdb + .trb
        # files; those live in the structures/ stream folder.
        structures_dir = self.stream_folder("structures")
        return f"""echo "Creating results table"
python {self.table_py_file} "{structures_dir}" "{output_name}" {self.num_designs} "{self.main_table}" {self.design_startnum}

"""

    def _generate_script_update_structures_map(self) -> str:
        """Generate script to write structures_map.csv from the actual runtime PDBs."""
        structures_map = self.stream_map_path("structures")
        structures_dir = self.stream_folder("structures")
        # All designs share the same parent PDB (if a PDB input was given).
        prov_arg = (f' --set-provenance "structures.id={self.pdb_input_id}"'
                    if self.pdb_input_id else "")
        return f"""echo "Writing structures map from actual output files"
python {self.update_map_py} --structures-map "{structures_map}" --output-folder "{structures_dir}"{prov_arg}

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after RFdiffusion-AllAtom execution."""
        output_name = self.pdb_input_id if self.pdb_input_id else self.pipeline_name

        start = self.design_startnum
        end = self.design_startnum + self.num_designs - 1
        structure_ids = [f"{output_name}_<{start}..{end}>"]
        file_template = [self.stream_path("structures", "<id>.pdb")]

        # The per-design map_table is written at runtime by
        # _generate_script_update_structures_map(); here we only declare the
        # stream and its map_table path.
        structures_map = self.stream_map_path("structures")

        structures = DataStream(
            name="structures",
            ids=structure_ids,
            files=file_template,
            map_table=structures_map,
            format="pdb"
        )

        tables = {
            "structures": TableInfo(
                name="structures",
                path=self.main_table,
                columns=["id", "pdb", "fixed", "designed", "source_fixed", "plddt_mean", "status"],
                description="RFdiffusion-AllAtom structure generation results with fixed/designed regions"
            )
        }

        return {
            "structures": structures,
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration with all RFdiffusion-AllAtom parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "rfdaa_params": {
                "pdb_input_id": self.pdb_input_id,
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
    TOOL_VERSION = "1.0"

    # Lazy path descriptors
    #   prepared_pdb    — the single output PDB, lives in structures/.
    #   structures_csv  — standalone TableInfo (tables/structures.csv).
    #   ligand_ds_json  — config-time input DataStream serialization.
    prepared_pdb = Path(lambda self: self.stream_path("structures", "prepared_ligand.pdb"))
    structures_csv = Path(lambda self: self.table_path("structures"))
    helper_script = Path(lambda self: self.pipe_script_path("pipe_rfdaa_prepare_ligand.py"))
    ligand_ds_json = Path(lambda self: self.configuration_path("input_ligand.json"))

    def __init__(self,
                 ligand: Union[DataStream, StandardizedOutput],
                 **kwargs):
        """
        Initialize RFDAA_PrepareLigand tool.

        Args:
            ligand: Ligand structure as DataStream or StandardizedOutput
            **kwargs: Additional parameters

        Output:
            Streams: structures (.pdb)
            Tables:
                structures: id | file_path
        """
        # Resolve ligand input — store stream for runtime resolution
        if isinstance(ligand, StandardizedOutput):
            self.ligand_stream = ligand.streams.structures
        elif isinstance(ligand, DataStream):
            self.ligand_stream = ligand
        else:
            raise ValueError(f"ligand must be DataStream or StandardizedOutput, got {type(ligand)}")
        self.ligand_input_id = self.ligand_stream.ids[0]

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate parameters."""
        if not self.ligand_stream:
            raise ValueError("ligand input is required")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure inputs."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display."""
        config_lines = super().get_config_display()
        config_lines.append(f"LIGAND_SOURCE: {self.ligand_input_id}")
        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate script to combine ligand with dummy peptide."""
        # configuration/ + structures/ folders auto-created by the pipeline.
        self.ligand_stream.save_json(self.ligand_ds_json)

        script_content = "#!/bin/bash\n"
        script_content += "# RFDAA_PrepareLigand execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += f"""LIGAND_ID={Resolve.stream_ids(self.ligand_ds_json, index=0)}
LIGAND_FILE={Resolve.stream_item(self.ligand_ds_json, '$LIGAND_ID')}

echo "Preparing ligand structure for RFdiffusion-AllAtom"
echo "Input ligand: $LIGAND_FILE"
echo "Output: {self.prepared_pdb}"

python "{self.helper_script}" \\
  --ligand_pdb "$LIGAND_FILE" \\
  --output_pdb "{self.prepared_pdb}" \\
  --output_csv "{self.structures_csv}" \\
  --pdbs_folder "{self.folders['pdbs']}"

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

        structures_map = self.stream_map_path("structures")
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
                description="Prepared ligand structure with dummy peptide"
            )
        }

        return {
            "structures": structures,
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                "ligand_input_id": self.ligand_input_id
            }
        })
        return base_dict
