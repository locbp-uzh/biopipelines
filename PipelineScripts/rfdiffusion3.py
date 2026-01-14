"""
RFdiffusion3 configuration for all-atom protein design via foundry framework.

Third-generation diffusion model for fast, all-atom protein design with support
for hotspot-driven binder design, partial diffusion, and flexible structure control.
Approximately 10× faster than RFdiffusion2 with higher success rates.
"""

import os
import json
import re
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo
except ImportError:
    # Fallback for direct execution
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo


class RFdiffusion3(BaseConfig):
    """
    Configuration for RFdiffusion3 all-atom protein design.

    RFdiffusion3 is a fast all-atom diffusion model that operates at the atomic level
    (4 backbone + 10 sidechain atoms per residue) for precise protein design including
    hotspot-driven binder design, enzyme design, and symmetric assemblies.

    Requirements:
        - Python ≥3.12
        - foundry environment: pip install "rc-foundry[all]"
        - Checkpoints at: /home/$USER/data/rfdiffusion3/
        - Environment variable: FOUNDRY_CHECKPOINT_DIRS (optional override)

    Examples:
        # De novo design (no input PDB)
        designs = RFdiffusion3(length="100-120", num_designs=10)

        # Binder design with hotspots (requires input PDB)
        target = PDB(pdb="7KDL")
        binder = RFdiffusion3(
            pdb=target,
            contig="A50-100,80-100,\\0,A1-50",
            select_hotspots="A67,A89",
            num_designs=20
        )

        # Advanced: Full JSON control
        config = {
            "design_1": {
                "contig": "50-80,\\0,A1-100",
                "length": "150-200",
                "select_unfixed_sequence": "A20-35",
                "partial_t": 10.0
            }
        }
        designs = RFdiffusion3(json_config=config)

    Parameters:
        length (str or int): Length constraint for de novo design (no input PDB).
            Use "min-max" for range or int for exact length.
            Example: "100-150" or 120
        contig (str): Contig specification for motif-based design (requires input PDB).
            Use '\\0' for chain breaks. Chain letters reference input structure.
            Example: "A50-100,80-100,\\0,A1-50" (keep A50-100, design 80-100, break, keep A1-50)
        pdb (str or ToolOutput): Input PDB structure (required when using contig)
        ligand (str): Ligand selection by chemical component name
        num_designs (int): Number of designs to generate (default: 1)
        num_models (int): Number of models per design (default: 1).
            WARNING: RFdiffusion3's internal default is 8 models per design. Always
            explicitly pass this parameter to avoid unexpected behavior.
        prefix (str): Prefix for output file names (default: uses pipeline name)
        select_hotspots (str or dict): Hotspot residues for binder design
            String: "A67,A89" (all atoms) or "A67:CA,CB;A89:CA" (specific atoms)
            Dict: {"A67": "CA,CB", "A89": ""}
        json_config (str or dict): Override with full JSON configuration for advanced use
        design_startnum (int): Starting number for design numbering (default: 1)

    Outputs:
        structures: List of PDB files ({prefix}_d{D}_m{M}.pdb, ...)
        tables.structures: CSV with columns: id, design, model, pdb, contig, length, time, status

    Notes:
        - 10× faster than RFdiffusion/RFdiffusionAllAtom
        - All-atom model (4 backbone + 10 sidechain atoms)
        - Use 'length' for de novo design, 'contig' for motif-based design
        - 'contig' requires input PDB, even for numeric ranges
        - Chain breaks use '\\0' not '/' (different from RFdiffusion)
        - Advanced parameters available via json_config

    See Also:
        RFdiffusion, RFdiffusionAllAtom: Earlier versions
        ProteinMPNN: Sequence design for generated backbones

    References:
        Paper: https://www.biorxiv.org/content/10.1101/2024.11.13.623358v1
        GitHub: https://github.com/RosettaCommons/foundry
    """

    TOOL_NAME = "RFdiffusion3"
    DEFAULT_ENV = None  # Loaded from config.yaml

    def __init__(self,
                 contig: str = "",
                 length: Union[str, int] = None,
                 pdb: Union[str, ToolOutput, StandardizedOutput] = "",
                 ligand: str = "",
                 num_designs: int = 1,
                 num_models: int = 1,
                 prefix: str = None,
                 select_hotspots: Union[str, Dict[str, str]] = None,
                 json_config: Union[str, Dict] = None,
                 design_startnum: int = 1,
                 **kwargs):
        """
        Initialize RFdiffusion3 configuration.

        Args:
            contig: Contig specification (use '\\0' for chain breaks)
            length: Length constraint (str "min-max" or int)
            pdb: Input PDB structure (optional)
            ligand: Ligand selection by name
            num_designs: Number of designs to generate
            num_models: Number of models per design (default: 1). WARNING: RFdiffusion3's
                internal default is 8. Always pass explicitly.
            prefix: Prefix for output file names (defaults to pipeline name)
            select_hotspots: Hotspot residues specification
            json_config: Full JSON config override for advanced use
            design_startnum: Starting number for design IDs
            **kwargs: Additional parameters passed to BaseConfig
        """
        # Store parameters
        self.contig = contig
        self.length = length
        self.pdb = pdb
        self.ligand = ligand
        self.num_designs = num_designs
        self.num_models = num_models
        self.prefix = prefix
        self.select_hotspots = select_hotspots
        self.json_config = json_config
        self.design_startnum = design_startnum

        # Handle PDB input from tool outputs
        self.pdb_is_tool_output = False
        self.pdb_source_file = None

        if isinstance(pdb, StandardizedOutput):
            # StandardizedOutput from upstream tool
            self.pdb_is_tool_output = True
            if hasattr(pdb, 'structures') and pdb.structures:
                self.pdb_source_file = pdb.structures[0]
                if len(pdb.structures) > 1:
                    print(f"Warning: Multiple structures provided ({len(pdb.structures)}), using first: {pdb.structures[0]}")
        elif isinstance(pdb, ToolOutput):
            # Direct ToolOutput
            self.pdb_is_tool_output = True
            structures = pdb.get_output_files("structures")
            if structures:
                self.pdb_source_file = structures[0]
                self.dependencies.append(pdb.config)
                if len(structures) > 1:
                    print(f"Warning: Multiple structures from {pdb.tool_type}, using first")

        # Initialize base class
        super().__init__(**kwargs)

        # Initialize file paths (will be set properly in configure_inputs)
        self._initialize_file_paths()

    def validate_params(self):
        """Validate RFdiffusion3-specific parameters."""
        # Require either length, contig, or json_config
        if not self.length and not self.contig and not self.json_config:
            raise ValueError("Either length, contig, or json_config parameter is required")

        # Check for incorrect chain break syntax
        if self.contig and '/' in self.contig:
            raise ValueError(
                "RFdiffusion3 uses '\\0' for chain breaks, not '/'. "
                "Please update your contig specification. "
                "Example: '50-80,\\0,A1-100' instead of '50-80,/,A1-100'"
            )

        # Validate num_designs
        if self.num_designs <= 0:
            raise ValueError("num_designs must be positive")

        # Validate JSON config if provided
        if self.json_config:
            if isinstance(self.json_config, str):
                try:
                    json.loads(self.json_config)
                except json.JSONDecodeError as e:
                    raise ValueError(f"Invalid JSON config: {e}")

        # Validate hotspots format if provided as string
        if self.select_hotspots and isinstance(self.select_hotspots, str):
            # Basic validation - detailed parsing in _format_hotspots
            if not re.match(r'^[A-Z0-9,:; ]+$', self.select_hotspots):
                raise ValueError(
                    "Invalid hotspots format. Use 'A67,A89' or 'A67:CA,CB;A89:NE'"
                )

    def _initialize_file_paths(self):
        """Initialize file path placeholders."""
        self.pipeline_name = None
        self.json_file = None
        self.main_table = None
        self.metrics_csv = None
        self.specifications_csv = None
        self.rfd_log_file = None
        self.checkpoint_dir = None
        self.table_py_file = None
        self.postprocess_py_file = None
        self.input_pdb_file = None
        self.raw_output_folder = None

    def _extract_pipeline_name(self) -> str:
        """Extract pipeline/job name from output folder structure."""
        # Structure: .../JobName_NNN/NNN_RFdiffusion3 -> JobName_NNN
        folder_parts = self.output_folder.split(os.sep)
        for i, part in enumerate(folder_parts):
            if "_RFdiffusion3" in part or part.endswith("_RFdiffusion3"):
                if i == 0:
                    raise ValueError(f"Invalid output folder structure: {self.output_folder}")
                return folder_parts[i-1]

        raise ValueError(f"Could not extract pipeline name from: {self.output_folder}")

    def _setup_file_paths(self):
        """Set up all file paths after output_folder is known."""
        # Extract pipeline name
        self.pipeline_name = self._extract_pipeline_name()

        # Use provided prefix or default to pipeline name
        if self.prefix is None:
            self.prefix = self.pipeline_name

        # Core files
        self.json_file = os.path.join(
            self.output_folder,
            f"{self.prefix}_rfd3_input.json"
        )
        self.main_table = os.path.join(
            self.output_folder,
            "rfdiffusion3_results.csv"
        )
        self.metrics_csv = os.path.join(
            self.output_folder,
            "rfdiffusion3_metrics.csv"
        )
        self.specifications_csv = os.path.join(
            self.output_folder,
            "rfdiffusion3_specifications.csv"
        )

        # Raw output folder for CIF.gz files
        self.raw_output_folder = os.path.join(
            self.output_folder,
            "raw_output"
        )

        # Log file in parent Logs folder
        parent_dir = os.path.dirname(self.output_folder)
        folder_name = os.path.basename(self.output_folder)
        index = folder_name.split('_')[0] if '_' in folder_name else "000"
        self.rfd_log_file = os.path.join(parent_dir, "Logs", f"{index}_RFdiffusion3.log")

        # Helper script and checkpoint paths
        if hasattr(self, 'folders') and self.folders:
            # Checkpoint directory
            tool_data_folder = self.folders.get("tool_data", {})
            if isinstance(tool_data_folder, dict):
                checkpoint_base = tool_data_folder.get("RFdiffusion3", "rfdiffusion3")
            else:
                checkpoint_base = "rfdiffusion3"

            self.checkpoint_dir = os.path.join(
                os.path.expanduser("~"),
                "data",
                checkpoint_base
            )

            # HelpScript paths
            self.table_py_file = os.path.join(
                self.folders["HelpScripts"],
                "pipe_rfdiffusion3_table.py"
            )
            self.postprocess_py_file = os.path.join(
                self.folders["HelpScripts"],
                "pipe_rfdiffusion3_postprocess.py"
            )

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """
        Configure input sources from pipeline context.

        Args:
            pipeline_folders: Dictionary of pipeline folder paths
        """
        self.folders = pipeline_folders
        self._setup_file_paths()

        # Handle PDB input
        if self.pdb_is_tool_output:
            # File path already extracted in __init__
            self.input_sources["pdb"] = self.pdb_source_file
            self.input_pdb_file = self.pdb_source_file
        elif self.pdb:
            # String filename - look in PDBs folder
            pdb_filename = self.pdb if self.pdb.endswith(".pdb") else self.pdb + ".pdb"
            pdb_source = os.path.join(pipeline_folders["PDBs"], pdb_filename)

            if os.path.exists(pdb_source):
                self.input_sources["pdb"] = pdb_source
                self.input_pdb_file = pdb_source
            else:
                raise FileNotFoundError(f"PDB file not found: {pdb_source}")

    def _format_hotspots(self) -> Dict[str, str]:
        """
        Convert hotspots to JSON format.

        Formats:
            - Dict: {"A45": "CA,CB", "A67": ""} → return as-is
            - String: "A45,A67" → {"A45": "", "A67": ""}
            - String with atoms: "A45:CA,CB;A67:NE" → {"A45": "CA,CB", "A67": "NE"}

        Returns:
            Dictionary mapping residue IDs to atom names
        """
        if isinstance(self.select_hotspots, dict):
            return self.select_hotspots
        elif isinstance(self.select_hotspots, str):
            result = {}
            # Split by semicolon for different residues
            for spec in self.select_hotspots.split(';'):
                spec = spec.strip()
                if ':' in spec:
                    # Format: A45:CA,CB
                    residue, atoms = spec.split(':', 1)
                    result[residue.strip()] = atoms.strip()
                elif ',' in spec:
                    # Format: A45,A67 (comma-separated residues without atoms)
                    for res in spec.split(','):
                        res = res.strip()
                        if res:
                            result[res] = ""
                else:
                    # Single residue
                    if spec:
                        result[spec] = ""
            return result
        else:
            return {}

    def _build_json_config(self) -> Dict[str, Any]:
        """
        Build JSON configuration from parameters.

        Returns:
            Dictionary representing JSON config
        """
        # If full JSON config provided, use it
        if self.json_config:
            if isinstance(self.json_config, dict):
                return self.json_config
            else:
                return json.loads(self.json_config)

        # Build config from parameters
        # Use prefix (which defaults to pipeline_name if not set)
        prefix = self.prefix if self.prefix else self.pipeline_name

        # Create multiple design entries for num_designs
        config = {}
        for i in range(self.num_designs):
            design_key = f"{prefix}_design_{i}"
            config[design_key] = {}
            entry = config[design_key]

            # Required/common parameters
            if self.contig:
                entry["contig"] = self.contig

            if self.length is not None:
                entry["length"] = str(self.length)

            if self.input_pdb_file:
                entry["input"] = self.input_pdb_file

            if self.ligand:
                entry["ligand"] = self.ligand

            if self.select_hotspots:
                entry["select_hotspots"] = self._format_hotspots()

            # Set number of models per design
            if self.num_models > 1:
                entry["diffusion_batch_size"] = self.num_models

        return config

    def get_config_display(self) -> List[str]:
        """Get RFdiffusion3 configuration display lines."""
        config_lines = super().get_config_display()

        # Input information
        if self.pdb_is_tool_output:
            config_lines.append(f"INPUT PDB: {os.path.basename(self.pdb_source_file)}")
        elif self.pdb:
            config_lines.append(f"INPUT PDB: {self.pdb}")
        else:
            config_lines.append("INPUT: De novo design")

        if self.contig:
            config_lines.append(f"CONTIG: {self.contig}")

        if self.length:
            config_lines.append(f"LENGTH: {self.length}")

        if self.ligand:
            config_lines.append(f"LIGAND: {self.ligand}")

        if self.select_hotspots:
            hotspots_str = str(self.select_hotspots)
            if len(hotspots_str) > 50:
                hotspots_str = hotspots_str[:47] + "..."
            config_lines.append(f"HOTSPOTS: {hotspots_str}")

        config_lines.append(f"NUM DESIGNS: {self.num_designs}")
        config_lines.append(f"NUM MODELS: {self.num_models}")

        if self.prefix:
            config_lines.append(f"PREFIX: {self.prefix}")

        if self.json_config:
            config_lines.append("MODE: Advanced (JSON config)")

        return config_lines

    def _generate_json_section(self) -> str:
        """Generate bash section that creates JSON input file."""
        json_content = json.dumps(self._build_json_config(), indent=2)

        return f"""echo "Creating RFdiffusion3 JSON configuration"
cat > "{self.json_file}" << 'EOF'
{json_content}
EOF

"""

    def generate_script_run_rfdiffusion3(self) -> str:
        """Generate RFdiffusion3 execution bash code."""
        return f"""echo "Starting RFdiffusion3"
echo "JSON config: {self.json_file}"
echo "Output folder: {self.output_folder}"
echo "Raw output folder: {self.raw_output_folder}"

# Create raw output directory
mkdir -p "{self.raw_output_folder}"

# Set checkpoint directory
export FOUNDRY_CHECKPOINT_DIRS="{self.checkpoint_dir}"

# Check checkpoint directory exists
if [ ! -d "${{FOUNDRY_CHECKPOINT_DIRS}}" ]; then
    echo "ERROR: RFdiffusion3 checkpoints not found at ${{FOUNDRY_CHECKPOINT_DIRS}}"
    echo "Please ensure checkpoints are installed at /home/$USER/data/rfdiffusion3/"
    exit 1
fi

# Run RFdiffusion3 (outputs CIF.gz format to raw folder)
rfd3 design \\
    out_dir="{self.raw_output_folder}" \\
    inputs="{self.json_file}" \\
    global_prefix="{self.prefix}"

"""

    def _generate_postprocess_section(self) -> str:
        """Generate bash section to post-process RFdiffusion3 outputs.

        Converts CIF.gz outputs to PDB format and extracts metrics from JSON files.
        Runs in biopipelines environment (has BioPython, pandas).
        """
        return f"""echo "Post-processing RFdiffusion3 outputs"

# Process CIF.gz files: decompress, convert to PDB, extract metrics
mamba run -n biopipelines python "{self.postprocess_py_file}" \\
    --raw_folder "{self.raw_output_folder}" \\
    --output_folder "{self.output_folder}" \\
    --prefix "{self.prefix}" \\
    --num_designs {self.num_designs} \\
    --num_models {self.num_models} \\
    --design_startnum {self.design_startnum} \\
    --metrics_csv "{self.metrics_csv}" \\
    --specifications_csv "{self.specifications_csv}"

"""

    def generate_script_create_table(self) -> str:
        """Generate table creation bash code."""
        return f"""echo "Creating results table"
python "{self.table_py_file}" \\
    --output_folder "{self.output_folder}" \\
    --json_file "{self.json_file}" \\
    --pipeline_name "{self.prefix}" \\
    --num_designs {self.num_designs} \\
    --num_models {self.num_models} \\
    --table_path "{self.main_table}" \\
    --design_startnum {self.design_startnum}

"""

    def generate_script(self, script_path: str) -> str:
        """
        Generate RFdiffusion3 execution script.

        Args:
            script_path: Path where script should be written

        Returns:
            Script content as string
        """
        script_content = "#!/bin/bash\n"
        script_content += "# RFdiffusion3 execution script\n"
        script_content += "# Generated by BioPipelines\n\n"
        script_content += self.generate_completion_check_header()
        script_content += self._generate_json_section()
        script_content += self.generate_script_run_rfdiffusion3()
        script_content += self._generate_postprocess_section()
        script_content += self.generate_script_create_table()
        script_content += self.generate_completion_check_footer()

        return script_content

    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after RFdiffusion3 execution.

        Uses pure path construction - no filesystem access.

        Returns:
            Dictionary mapping output types to file paths with standard keys:
            - structures: PDB structure files
            - structure_ids: Structure identifiers
            - tables: TableInfo objects
            - output_folder: Tool's output directory
        """
        # Ensure file paths are set up
        if not hasattr(self, 'pipeline_name') or self.pipeline_name is None:
            self._setup_file_paths()

        # Generate expected structure paths
        design_pdbs = []
        structure_ids = []

        total_structures = self.num_designs * self.num_models
        for i in range(self.num_designs):
            for j in range(self.num_models):
                design_num = self.design_startnum + i
                model_num = self.design_startnum + j
                structure_id = f"{self.prefix}_d{design_num}_m{model_num}"
                structure_path = os.path.join(self.output_folder, f"{structure_id}.pdb")
                design_pdbs.append(structure_path)
                structure_ids.append(structure_id)

        # Define table structure
        tables = {
            "structures": TableInfo(
                name="structures",
                path=self.main_table,
                columns=["id", "design", "model", "pdb", "contig", "length", "time", "status"],
                description="RFdiffusion3 structure generation results",
                count=total_structures
            ),
            "metrics": TableInfo(
                name="metrics",
                path=self.metrics_csv,
                columns=[
                    "id", "design", "model", "max_ca_deviation", "n_chainbreaks",
                    "n_clashing_interresidue_w_sidechain", "n_clashing_interresidue_w_backbone",
                    "non_loop_fraction", "loop_fraction", "helix_fraction", "sheet_fraction",
                    "num_ss_elements", "radius_of_gyration", "alanine_content",
                    "glycine_content", "num_residues"
                ],
                description="RFdiffusion3 quality metrics extracted from JSON outputs",
                count=total_structures
            ),
            "specifications": TableInfo(
                name="specifications",
                path=self.specifications_csv,
                columns=[
                    "id", "design", "model", "sampled_contig", "num_tokens_in", "num_residues_in",
                    "num_chains", "num_atoms", "num_residues"
                ],
                description="RFdiffusion3 design specifications and statistics",
                count=total_structures
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
            # Legacy aliases for backward compatibility
            "pdbs": design_pdbs,
            "main": self.main_table
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including RFdiffusion3-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "rfd3_params": {
                "contig": self.contig,
                "length": self.length,
                "ligand": self.ligand,
                "num_designs": self.num_designs,
                "num_models": self.num_models,
                "prefix": self.prefix,
                "select_hotspots": self.select_hotspots,
                "has_json_config": self.json_config is not None,
                "design_startnum": self.design_startnum,
                "input_type": "tool_output" if self.pdb_is_tool_output else "direct"
            }
        })
        return base_dict
