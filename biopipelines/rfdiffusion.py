# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
RFdiffusion configuration for protein backbone generation.

RFdiffusion is a diffusion-based generative model for designing protein backbones.
It supports unconditional generation, motif scaffolding, binder design, and partial
diffusion. See https://github.com/RosettaCommons/RFdiffusion for full documentation.
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

    RFdiffusion generates novel protein backbones using a denoising diffusion
    probabilistic model conditioned on structural motifs, partial structures,
    or protein-protein interaction targets.

    Main workflows:
        - Unconditional generation: produce backbones of a given length range
          without any input structure (contigs only, no pdb).
        - Motif scaffolding: hold a fixed structural motif in place and design
          the surrounding scaffold (pdb + contigs specifying the motif).
        - Binder design: design a new protein that binds to a target (pdb of
          the target + contigs describing the binder length).
        - Partial diffusion: apply limited noise to an existing structure and
          re-diffuse, producing near-neighbour variants (partial_steps).

    Model checkpoints (see WEIGHTS / DEFAULT_WEIGHTS):
        By default only Base and Complex_base are downloaded (~2 GB). Add others
        via RFdiffusion.install(weights=[...]) when needed.

    Installation:
        with Pipeline(...):
            RFdiffusion.install()                         # default weights
            RFdiffusion.install(weights=["Base",
                                         "InpaintSeq"])   # custom selection
            rfd = RFdiffusion(contigs="50-100", ...)

    Reference:
        Watson et al. (2023) De novo design of protein structure and function
        with RFdiffusion. Nature 620, 1089-1100.
        https://github.com/RosettaCommons/RFdiffusion
    """

    TOOL_NAME = "RFdiffusion"

    # Mapping of weight name -> (url_hash, filename)
    WEIGHTS = {
        "Base":              ("6f5902ac237024bdd0c176cb93063dc4", "Base_ckpt.pt"),
        "Complex_base":      ("e29311f6f1bf1af907f9ef9f44b8328b", "Complex_base_ckpt.pt"),
        "Complex_Fold_base": ("60f09a193fb5e5ccdc4980417708dbab", "Complex_Fold_base_ckpt.pt"),
        "InpaintSeq":        ("74f51cfb8b440f50d70878e05361d8f0", "InpaintSeq_ckpt.pt"),
        "InpaintSeq_Fold":   ("76d00716416567174cdb7ca96e208296", "InpaintSeq_Fold_ckpt.pt"),
        "ActiveSite":        ("5532d2e1f3a4738decd58b19d633b3c3", "ActiveSite_ckpt.pt"),
        "Base_epoch8":       ("12fc204edeae5b57713c5ad7dcb97d39", "Base_epoch8_ckpt.pt"),
        "RF_structure_prediction": ("1befcb9b28e2f778f53d47f18b7597fa", "RF_structure_prediction_weights.pt"),
    }
    # Weights downloaded by default — covers the two main workflows:
    #   Base          → unconditional generation, motif scaffolding
    #   Complex_base  → binder design / protein-protein interaction
    # Add more via RFdiffusion.install(weights=[...]):
    #   Complex_Fold_base  → binder design with fold/topology conditioning
    #   InpaintSeq         → auto-selected when contigmap.inpaint_seq is used
    #   InpaintSeq_Fold    → inpaint + fold conditioning together
    #   ActiveSite         → small active-site motifs (requires ckpt_override_path)
    #   Base_epoch8        → alternative base checkpoint
    #   RF_structure_prediction → legacy structure prediction (unrelated to diffusion)
    DEFAULT_WEIGHTS = ["Base", "Complex_base"]

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False,
                        weights=None, **kwargs):
        """
        Generate the bash installation script for RFdiffusion.

        Clones the RFdiffusion repository, downloads the requested model
        checkpoints, and creates (or reuses) the SE3nv conda/mamba environment.
        In pip mode (e.g. Google Colab) the environment creation step is skipped
        and dependencies are installed directly.

        Args:
            folders: Resolved pipeline folder paths (must contain "RFdiffusion"
                     and "biopipelines" keys).
            env_manager: "mamba", "conda", or "pip".
            force_reinstall: If True, skip the already-installed early-exit check.
            weights: List of checkpoint names to download. Defaults to
                     DEFAULT_WEIGHTS (["Base", "Complex_base"]).
                     Valid names are the keys of WEIGHTS.
        """
        biopipelines = folders.get("biopipelines", "")
        repo_dir = folders.get("RFdiffusion", "")
        parent_dir = os.path.dirname(repo_dir)

        if weights is None:
            weights = cls.DEFAULT_WEIGHTS
        invalid = [w for w in weights if w not in cls.WEIGHTS]
        if invalid:
            raise ValueError(
                f"Unknown weight(s): {invalid}. "
                f"Valid options are: {list(cls.WEIGHTS.keys())}"
            )
        wget_lines = "\n".join(
            f"wget -nc http://files.ipd.uw.edu/pub/RFdiffusion/{cls.WEIGHTS[w][0]}/{cls.WEIGHTS[w][1]}"
            f" || echo \"WARNING: failed to download {cls.WEIGHTS[w][1]}, skipping\""
            for w in weights
        )

        if env_manager == "pip":
            skip = "" if force_reinstall else f"""# Check if already installed
if [ -d "{repo_dir}/models" ] && [ -f "{repo_dir}/models/Base_ckpt.pt" ]; then
    echo "RFdiffusion already installed, skipping. Use force_reinstall=True to reinstall."
    exit 0
fi
"""
            return f"""echo "=== Installing RFdiffusion (pip) ==="
{skip}mkdir -p {parent_dir}
cd {parent_dir}
if [ ! -d "{repo_dir}" ]; then
    git clone https://github.com/RosettaCommons/RFdiffusion.git
fi
cd {repo_dir}

# Download model weights
mkdir -p models && cd models
{wget_lines}
cd ..

# Install dependencies via pip
pip install -r {biopipelines}/Environments/SE3nv_pip_requirements.txt

# Install SE3Transformer and RFdiffusion
cd env/SE3Transformer
pip install --no-cache-dir -r requirements.txt
python setup.py install
cd ../..
pip install -e .

echo "=== RFdiffusion installation complete ==="
"""
        skip = "" if force_reinstall else f"""# Check if already installed
if [ -d "{repo_dir}/models" ] && [ -f "{repo_dir}/models/Base_ckpt.pt" ] && {env_manager} env list 2>/dev/null | grep -q "SE3nv"; then
    echo "RFdiffusion already installed, skipping. Use force_reinstall=True to reinstall."
    exit 0
fi
"""
        return f"""echo "=== Installing RFdiffusion ==="
{skip}cd {parent_dir}
git clone https://github.com/RosettaCommons/RFdiffusion.git
cd {repo_dir}

# Download model weights
mkdir -p models && cd models
{wget_lines}
cd ..

# Create SE3nv environment
# Try BioPipelines SE3nv.yaml first (includes additional packages for the pipeline)
echo "Creating SE3nv environment from BioPipelines specification..."
{env_manager} env create -f {biopipelines}/Environments/SE3nv.yaml
if [ $? -ne 0 ]; then
    echo "WARNING: BioPipelines SE3nv.yaml failed. Trying official RFdiffusion environment..."
    {env_manager} env create -f env/SE3nv.yml
    if [ $? -ne 0 ]; then
        echo "ERROR: SE3nv environment creation failed with both methods."
        echo "This is likely a CUDA version mismatch for your system."
        echo "Options:"
        echo "  1. Edit {biopipelines}/Environments/SE3nv.yaml to match your CUDA version"
        echo "  2. Edit RFdiffusion/env/SE3nv.yml following https://github.com/RosettaCommons/RFdiffusion"
        exit 1
    fi
fi
{env_manager} activate SE3nv
pip install -r {biopipelines}/Environments/SE3nv_pip_requirements.txt

# Install SE3Transformer and RFdiffusion
cd env/SE3Transformer
pip install --no-cache-dir -r requirements.txt
python setup.py install
cd ../..
pip install -e .

echo "=== RFdiffusion installation complete ==="
"""

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
            contigs: Contig map describing which residues are fixed and which
                     are to be generated. Residues from the input PDB are
                     specified as chain+range (e.g. "A1-100"); new backbone
                     segments are specified as length ranges (e.g. "50-100").
                     Multiple segments are separated by "/" (e.g. "A1-50/30-50/A60-100").
                     For unconditional generation (no pdb) use a length range alone
                     (e.g. "100-200").
            pdb: Optional input structure. Required for motif scaffolding,
                 binder design, and partial diffusion. Accepts a DataStream
                 (single file) or a StandardizedOutput (first structure is used).
            inpaint: Residues whose sequence should be masked during diffusion
                     (same chain+range format as contigs). When set, the
                     InpaintSeq checkpoint is used automatically — ensure it
                     was downloaded at install time.
            num_designs: Number of independent backbone designs to generate.
            active_site: If True, use the ActiveSite checkpoint instead of
                         Base. Intended for scaffolding very small functional
                         motifs (< ~10 residues). Requires the ActiveSite
                         weight to have been downloaded at install time.
            steps: Number of denoising diffusion steps (default 50). Fewer
                   steps are faster but may reduce quality.
            partial_steps: Number of partial diffusion steps. When > 0, the
                           input structure is noised for this many steps and
                           then re-diffused, producing near-neighbour variants.
                           Must be less than steps.
            reproducible: If True, use deterministic (fixed-seed) sampling so
                          the same contigs always produce the same outputs.
            design_startnum: Integer appended to the pipeline name to number
                             output files (e.g. design_startnum=1 → name_1.pdb).
                             Useful when continuing a previous run.

        Output:
            Streams: structures (.pdb)
            Tables:
                structures: id | source_id | pdb | fixed | designed | contigs | time | status
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
