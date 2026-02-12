# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
RFdiffusion3 configuration for all-atom protein design via foundry framework.

Third-generation diffusion model for fast, all-atom protein design with support
for hotspot-driven binder design, partial diffusion, and flexible structure control.
Approximately 10x faster than RFdiffusion2 with higher success rates.
"""

import os
import json
import re
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


class RFdiffusion3(BaseConfig):
    """
    Configuration for RFdiffusion3 all-atom protein design.

    RFdiffusion3 is a fast all-atom diffusion model that operates at the atomic level
    (4 backbone + 10 sidechain atoms per residue) for precise protein design including
    hotspot-driven binder design, enzyme design, and symmetric assemblies.

    Requirements:
        - Python >=3.12
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
        pdb (DataStream or StandardizedOutput): Input PDB structure (required when using contig)
        ligand_code (str): Ligand three-letter code identifying the molecule in input structure
        ligand_structure (DataStream or StandardizedOutput): Output from Ligand tool providing ligand PDB file.
            When provided, this PDB becomes the input structure for design.
        num_designs (int): Number of designs to generate (default: 1)
        num_models (int): Number of models per design (default: 1).
            WARNING: RFdiffusion3's internal default is 8 models per design. Always
            explicitly pass this parameter to avoid unexpected behavior.
        prefix (str): Prefix for output file names (default: uses pipeline name)
        select_hotspots (str or dict): Hotspot residues for binder design
            String: "A67,A89" (all atoms) or "A67:CA,CB;A89:CA" (specific atoms)
            Dict: {"A67": "CA,CB", "A89": ""}
        select_fixed_atoms (bool, str, or dict): Atoms with fixed 3D coordinates.
            True=all atoms fixed, ""=none fixed, dict=specific atoms per residue/ligand.
            Example: {"AXL": ""} to not fix any ligand atoms
        select_buried (str or dict): Atoms that should be buried in protein (RASA control).
            Example: {"AXL": "C1,C2,C3"} or "AXL" for all atoms
        select_exposed (str or dict): Atoms that should be solvent-exposed (RASA control).
            Example: {"AXL": "O1,O2"} or "AXL" for all atoms
        select_hbond_donor (dict): Hydrogen bond donor specification.
            Dict mapping residue/ligand to donor atoms. Example: {"AXL": "N1,N2"}
        select_hbond_acceptor (dict): Hydrogen bond acceptor specification.
            Dict mapping residue/ligand to acceptor atoms. Example: {"AXL": "O1,O2"}
        json_config (str or dict): Override with full JSON configuration for advanced use
        design_startnum (int): Starting number for design numbering (default: 1)

    Outputs:
        structures: DataStream of PDB files ({prefix}_d{D}_m{M}.pdb, ...)
        tables.structures: CSV with columns: id, design, model, pdb, contig, length, time, status

    Notes:
        - 10x faster than RFdiffusion/RFdiffusionAllAtom
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

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        repo_dir = folders.get("RFdiffusion3", "")
        skip = "" if force_reinstall else f"""# Check if already installed
if [ -d "{repo_dir}" ] && {env_manager} env list 2>/dev/null | grep -q "foundry"; then
    echo "RFdiffusion3 already installed, skipping. Use force_reinstall=True to reinstall."
    exit 0
fi
"""
        return f"""echo "=== Installing RFdiffusion3 (foundry) ==="
{skip}{env_manager} create -n foundry python=3.12 -y
{env_manager} activate foundry
pip install "rc-foundry[all]"

# Download model weights
mkdir -p {repo_dir}
foundry install rfd3 --checkpoint-dir {repo_dir}

echo "=== RFdiffusion3 installation complete ==="
"""

    # Lazy path descriptors
    json_file = Path(lambda self: os.path.join(self.output_folder, f"{self._get_prefix()}_rfd3_input.json"))
    main_table = Path(lambda self: os.path.join(self.output_folder, "rfdiffusion3_results.csv"))
    metrics_csv = Path(lambda self: os.path.join(self.output_folder, "rfdiffusion3_metrics.csv"))
    specifications_csv = Path(lambda self: os.path.join(self.output_folder, "rfdiffusion3_specifications.csv"))
    sequences_csv = Path(lambda self: os.path.join(self.output_folder, "rfdiffusion3_sequences.csv"))
    raw_output_folder = Path(lambda self: os.path.join(self.output_folder, "raw_output"))
    rfd_log_file = Path(lambda self: os.path.join(
        os.path.dirname(self.output_folder), "Logs",
        f"{os.path.basename(self.output_folder).split('_')[0] if '_' in os.path.basename(self.output_folder) else '000'}_RFdiffusion3.log"
    ))
    checkpoint_dir = Path(lambda self: self.folders["RFdiffusion3"])
    table_py_file = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_rfdiffusion3_table.py"))
    postprocess_py_file = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_rfdiffusion3_postprocess.py"))

    def __init__(self,
                 contig: str = "",
                 length: Union[str, int] = None,
                 pdb: Optional[Union[DataStream, StandardizedOutput]] = None,
                 ligand_code: str = "",
                 ligand_structure: Optional[Union[DataStream, StandardizedOutput]] = None,
                 num_designs: int = 1,
                 num_models: int = 1,
                 prefix: str = None,
                 select_hotspots: Union[str, Dict[str, str]] = None,
                 select_fixed_atoms: Union[bool, str, Dict[str, str]] = None,
                 select_buried: Union[str, Dict[str, str]] = None,
                 select_exposed: Union[str, Dict[str, str]] = None,
                 select_hbond_donor: Dict[str, str] = None,
                 select_hbond_acceptor: Dict[str, str] = None,
                 json_config: Union[str, Dict] = None,
                 design_startnum: int = 1,
                 **kwargs):
        """
        Initialize RFdiffusion3 configuration.

        Args:
            contig: Contig specification (use '\\0' for chain breaks)
            length: Length constraint (str "min-max" or int)
            pdb: Input PDB structure as DataStream or StandardizedOutput (optional)
            ligand_code: Ligand CCD code identifying the molecule
            ligand_structure: Output from Ligand tool providing ligand PDB file
            num_designs: Number of designs to generate
            num_models: Number of models per design (default: 1). WARNING: RFdiffusion3's
                internal default is 8. Always pass explicitly.
            prefix: Prefix for output file names (defaults to pipeline name)
            select_hotspots: Hotspot residues specification
            select_fixed_atoms: Atoms with fixed 3D coordinates (True/str/dict)
            select_buried: Atoms to bury in protein (RASA control)
            select_exposed: Atoms to expose to solvent (RASA control)
            select_hbond_donor: Hydrogen bond donor atoms (dict)
            select_hbond_acceptor: Hydrogen bond acceptor atoms (dict)
            json_config: Full JSON config override for advanced use
            design_startnum: Starting number for design IDs
            **kwargs: Additional parameters passed to BaseConfig
        """
        # Resolve PDB input
        self.input_pdb_file: Optional[str] = None

        # Handle ligand structure (takes precedence over pdb)
        if ligand_structure is not None:
            if isinstance(ligand_structure, StandardizedOutput):
                self.input_pdb_file = ligand_structure.streams.structures.files[0]
                # Auto-extract ligand code from compound_ids if not provided
                if not ligand_code and ligand_structure.streams.compounds:
                    ligand_code = ligand_structure.streams.compounds.ids[0] if ligand_structure.streams.compounds.ids else ""
            elif isinstance(ligand_structure, DataStream):
                self.input_pdb_file = ligand_structure.files[0]
            else:
                raise ValueError(f"ligand_structure must be DataStream or StandardizedOutput, got {type(ligand_structure)}")
        # Handle PDB input (only if ligand_structure not provided)
        elif pdb is not None:
            if isinstance(pdb, StandardizedOutput):
                self.input_pdb_file = pdb.streams.structures.files[0]
                if len(pdb.streams.structures.files) > 1:
                    print(f"Warning: Multiple structures provided ({len(pdb.streams.structures.files)}), using first: {pdb.streams.structures.files[0]}")
            elif isinstance(pdb, DataStream):
                self.input_pdb_file = pdb.files[0]
                if len(pdb.files) > 1:
                    print(f"Warning: Multiple structures provided ({len(pdb.files)}), using first: {pdb.files[0]}")
            else:
                raise ValueError(f"pdb must be DataStream or StandardizedOutput, got {type(pdb)}")

        # Store parameters
        self.contig = contig
        self.length = length
        self.ligand_code = ligand_code
        self.num_designs = num_designs
        self.num_models = num_models
        self.prefix = prefix
        self.select_hotspots = select_hotspots
        self.select_fixed_atoms = select_fixed_atoms
        self.select_buried = select_buried
        self.select_exposed = select_exposed
        self.select_hbond_donor = select_hbond_donor
        self.select_hbond_acceptor = select_hbond_acceptor
        self.json_config = json_config
        self.design_startnum = design_startnum

        # Initialize base class
        super().__init__(**kwargs)

    def _get_prefix(self) -> str:
        """Get the prefix to use for output files."""
        return self.prefix if self.prefix else self.pipeline_name

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

        # Validate constraint parameter types
        if self.select_fixed_atoms is not None:
            if not isinstance(self.select_fixed_atoms, (bool, str, dict)):
                raise ValueError("select_fixed_atoms must be bool, str, or dict")

        if self.select_buried is not None:
            if not isinstance(self.select_buried, (str, dict)):
                raise ValueError("select_buried must be str or dict")

        if self.select_exposed is not None:
            if not isinstance(self.select_exposed, (str, dict)):
                raise ValueError("select_exposed must be str or dict")

        if self.select_hbond_donor is not None:
            if not isinstance(self.select_hbond_donor, dict):
                raise ValueError("select_hbond_donor must be dict")

        if self.select_hbond_acceptor is not None:
            if not isinstance(self.select_hbond_acceptor, dict):
                raise ValueError("select_hbond_acceptor must be dict")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input files."""
        self.folders = pipeline_folders

    def _format_hotspots(self) -> Dict[str, str]:
        """
        Convert hotspots to JSON format.

        Formats:
            - Dict: {"A45": "CA,CB", "A67": ""} -> return as-is
            - String: "A45,A67" -> {"A45": "", "A67": ""}
            - String with atoms: "A45:CA,CB;A67:NE" -> {"A45": "CA,CB", "A67": "NE"}

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
        prefix = self._get_prefix()

        # Create single design entry (n_batches and diffusion_batch_size control multiplicity)
        config = {}
        design_key = f"{prefix}"
        config[design_key] = {}
        entry = config[design_key]

        # Input structure (ligand PDB or scaffold PDB)
        if self.input_pdb_file:
            entry["input"] = self.input_pdb_file

        # Ligand specification
        if self.ligand_code:
            entry["ligand"] = self.ligand_code

        # Contig specification (for motif scaffolding)
        if self.contig:
            entry["contig"] = self.contig

        # Length constraint
        if self.length is not None:
            entry["length"] = str(self.length)

        # Hotspot specification
        if self.select_hotspots:
            entry["select_hotspots"] = self._format_hotspots()

        # Constraint parameters for small molecule binder design
        if self.select_fixed_atoms is not None:
            entry["select_fixed_atoms"] = self.select_fixed_atoms

        if self.select_buried is not None:
            entry["select_buried"] = self.select_buried

        if self.select_exposed is not None:
            entry["select_exposed"] = self.select_exposed

        if self.select_hbond_donor is not None:
            entry["select_hbond_donor"] = self.select_hbond_donor

        if self.select_hbond_acceptor is not None:
            entry["select_hbond_acceptor"] = self.select_hbond_acceptor

        return config

    def get_config_display(self) -> List[str]:
        """Get RFdiffusion3 configuration display lines."""
        config_lines = super().get_config_display()

        # Input information
        if self.input_pdb_file:
            config_lines.append(f"INPUT PDB: {os.path.basename(self.input_pdb_file)}")
        else:
            config_lines.append("INPUT: De novo design")

        if self.contig:
            config_lines.append(f"CONTIG: {self.contig}")

        if self.length:
            config_lines.append(f"LENGTH: {self.length}")

        if self.ligand_code:
            config_lines.append(f"LIGAND CODE: {self.ligand_code}")

        if self.select_hotspots:
            hotspots_str = str(self.select_hotspots)
            if len(hotspots_str) > 50:
                hotspots_str = hotspots_str[:47] + "..."
            config_lines.append(f"HOTSPOTS: {hotspots_str}")

        # Display constraint parameters
        if self.select_fixed_atoms is not None:
            config_lines.append(f"FIXED ATOMS: {self._format_constraint_display(self.select_fixed_atoms)}")

        if self.select_buried is not None:
            config_lines.append(f"BURIAL CONSTRAINTS: {self._format_constraint_display(self.select_buried)}")

        if self.select_exposed is not None:
            config_lines.append(f"EXPOSURE CONSTRAINTS: {self._format_constraint_display(self.select_exposed)}")

        if self.select_hbond_donor is not None or self.select_hbond_acceptor is not None:
            config_lines.append("H-BOND CONSTRAINTS: defined")

        config_lines.append(f"NUM DESIGNS: {self.num_designs}")
        config_lines.append(f"NUM MODELS: {self.num_models}")

        if self.prefix:
            config_lines.append(f"PREFIX: {self.prefix}")

        if self.json_config:
            config_lines.append("MODE: Advanced (JSON config)")

        return config_lines

    def _format_constraint_display(self, constraint) -> str:
        """Format constraint parameter for display."""
        if isinstance(constraint, bool):
            return "all" if constraint else "none"
        elif isinstance(constraint, dict):
            if len(constraint) == 0:
                return "none"
            keys = list(constraint.keys())
            if len(keys) <= 2:
                return str(constraint)
            return f"{{{keys[0]}: ..., {keys[1]}: ...}} ({len(keys)} entries)"
        elif isinstance(constraint, str):
            if len(constraint) > 40:
                return constraint[:37] + "..."
            return constraint if constraint else "none"
        return str(constraint)

    def _generate_json_section(self) -> str:
        """Generate bash section that creates JSON input file."""
        json_config = self._build_json_config()

        # Write JSON config file at pipeline time (not SLURM time)
        os.makedirs(self.output_folder, exist_ok=True)
        with open(self.json_file, 'w') as f:
            json.dump(json_config, f, indent=2)

        return f"""echo "Using RFdiffusion3 JSON configuration: {self.json_file}"

"""

    def generate_script_run_rfdiffusion3(self) -> str:
        """Generate RFdiffusion3 execution bash code."""
        prefix = self._get_prefix()
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
    echo "Please ensure checkpoints are installed at ${{FOUNDRY_CHECKPOINT_DIRS}}"
    exit 1
fi

# Run RFdiffusion3 (outputs CIF.gz format to raw folder)
# n_batches: number of designs to generate (default: 1)
# diffusion_batch_size: number of models per design (default: 8)
rfd3 design \\
    out_dir="{self.raw_output_folder}" \\
    inputs="{self.json_file}" \\
    global_prefix="{prefix}" \\
    n_batches={self.num_designs} \\
    diffusion_batch_size={self.num_models}

"""

    def _generate_postprocess_section(self) -> str:
        """Generate bash section to post-process RFdiffusion3 outputs.

        Converts CIF.gz outputs to PDB format and extracts metrics from JSON files.
        Runs in biopipelines environment (has BioPython, pandas).
        """
        prefix = self._get_prefix()
        return f"""echo "Post-processing RFdiffusion3 outputs"

# Process CIF.gz files: decompress, convert to PDB, extract metrics and sequences
mamba run -n biopipelines python "{self.postprocess_py_file}" \\
    --raw_folder "{self.raw_output_folder}" \\
    --output_folder "{self.output_folder}" \\
    --prefix "{prefix}" \\
    --num_designs {self.num_designs} \\
    --num_models {self.num_models} \\
    --design_startnum {self.design_startnum} \\
    --metrics_csv "{self.metrics_csv}" \\
    --specifications_csv "{self.specifications_csv}" \\
    --sequences_csv "{self.sequences_csv}"

"""

    def generate_script_create_table(self) -> str:
        """Generate table creation bash code."""
        prefix = self._get_prefix()
        return f"""echo "Creating results table"
python "{self.table_py_file}" \\
    --output_folder "{self.output_folder}" \\
    --json_file "{self.json_file}" \\
    --pipeline_name "{prefix}" \\
    --num_designs {self.num_designs} \\
    --num_models {self.num_models} \\
    --table_path "{self.main_table}" \\
    --design_startnum {self.design_startnum} \\
    --specifications_csv "{self.specifications_csv}"

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
        script_content += self.activate_environment()
        script_content += self._generate_json_section()
        script_content += self.generate_script_run_rfdiffusion3()
        script_content += self._generate_postprocess_section()
        script_content += self.generate_script_create_table()
        script_content += self.generate_completion_check_footer()

        return script_content

    def get_output_files(self) -> Dict[str, Any]:
        """
        Get expected output files after RFdiffusion3 execution.

        Uses pure path construction - no filesystem access.

        Returns:
            Dictionary with DataStream objects:
            - structures: DataStream of PDB structure files
            - sequences: DataStream (references sequences CSV)
            - compounds: Empty DataStream
            - tables: Dict of TableInfo objects
            - output_folder: Tool's output directory
        """
        prefix = self._get_prefix()

        # Generate expected structure paths and IDs
        design_pdbs = []
        structure_ids = []

        total_structures = self.num_designs * self.num_models
        for i in range(self.num_designs):
            for j in range(self.num_models):
                design_num = self.design_startnum + i
                model_num = self.design_startnum + j

                # Conditional naming: include model suffix only if num_models > 1
                if self.num_models > 1:
                    structure_id = f"{prefix}_d{design_num}_m{model_num}"
                else:
                    structure_id = f"{prefix}_{design_num}"

                structure_path = os.path.join(self.output_folder, f"{structure_id}.pdb")
                design_pdbs.append(structure_path)
                structure_ids.append(structure_id)

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

        # Define table structures
        tables = {
            "structures": TableInfo(
                name="structures",
                path=self.main_table,
                columns=["id", "design", "model", "pdb", "fixed", "designed", "contig", "length", "time", "status"],
                description="RFdiffusion3 structure generation results with fixed/designed regions",
                count=total_structures
            ),
            "metrics": TableInfo(
                name="metrics",
                path=self.metrics_csv,
                columns=[
                    "id", "design", "model", "max_ca_deviation", "n_chainbreaks",
                    "n_clashing_interresidue_w_sidechain", "n_clashing_interresidue_w_backbone",
                    "ligand_clashes", "ligand_min_distance",
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
            ),
            "sequences": TableInfo(
                name="sequences",
                path=self.sequences_csv,
                columns=["id", "source_id", "source_pdb", "chain", "sequence", "length"],
                description="RFdiffusion3 designed protein sequences extracted from PDB",
                count=total_structures
            )
        }

        sequences = DataStream(
            name="sequences",
            ids=structure_ids,
            files=[],
            map_table=self.sequences_csv,
            format="csv"
        )

        return {
            "structures": structures,
            "sequences": sequences,
            "compounds": DataStream.empty("compounds", "sdf"),
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including RFdiffusion3-specific parameters."""
        base_dict = super().to_dict()

        base_dict.update({
            "rfd3_params": {
                "contig": self.contig,
                "length": self.length,
                "ligand_code": self.ligand_code,
                "num_designs": self.num_designs,
                "num_models": self.num_models,
                "prefix": self.prefix,
                "select_hotspots": self.select_hotspots,
                "select_fixed_atoms": self.select_fixed_atoms,
                "select_buried": self.select_buried,
                "select_exposed": self.select_exposed,
                "select_hbond_donor": self.select_hbond_donor,
                "select_hbond_acceptor": self.select_hbond_acceptor,
                "has_json_config": self.json_config is not None,
                "design_startnum": self.design_startnum,
                "input_pdb_file": self.input_pdb_file
            }
        })
        return base_dict
