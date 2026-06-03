# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
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
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .biopipelines_io import Resolve
    from .combinatorics import generate_multiplied_ids_pattern
    from .input_standardization import resolve_basic_input
    from .ligand import Ligand
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from biopipelines_io import Resolve
    from combinatorics import generate_multiplied_ids_pattern
    from input_standardization import resolve_basic_input
    from ligand import Ligand


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
        ligand (DataStream or StandardizedOutput): Ligand as a compounds stream
            (Ligand(code="LIG") or any compounds-producing tool). The residue
            `code` is read from the stream's `code` column at runtime. If the
            source also exposes a `structures` stream, that PDB becomes the
            input structure for design.
        num_designs (int): Number of designs to generate (default: 1)
        num_models (int): Number of models per design (default: 1).
            WARNING: RFdiffusion3's internal default is 8 models per design. Always
            explicitly pass this parameter to avoid unexpected behavior.
        prefix (str): Prefix for output file names (default: uses pipeline name)
        IMPORTANT — selection keys are <chain><resnum>, NOT residue names. Every
        per-residue/atom selector below (select_hotspots, select_fixed_atoms,
        select_buried, select_exposed, select_hbond_donor/acceptor) keys its dict by
        the chain id + residue number (e.g. "A84", and for a LIGAND its chain+resnum
        like "B1"), the same "contig index" form RFD3 uses internally — NOT the
        3-letter ligand code. Passing a ligand code key (e.g. {"LIG": "..."} or
        {"AXL": "..."}) silently selects zero atoms, which for RASA selectors fails
        with "could not broadcast ... (0,3)". To target a ligand, use its chain and
        residue number (find them in the input PDB; e.g. a dye on chain B residue 1
        is "B1").

        select_hotspots (str or dict): Hotspot residues for binder design
            String: "A67,A89" (all atoms) or "A67:CA,CB;A89:CA" (specific atoms)
            Dict: {"A67": "CA,CB", "A89": ""}
        select_fixed_atoms (bool, str, or dict): Atoms with fixed 3D coordinates.
            True=all atoms fixed, ""=none fixed, dict=specific atoms per <chain><resnum>.
            Example: {"B1": ""} to not fix any atoms of the ligand at chain B residue 1
        select_buried (str or dict): Atoms that should be buried in protein (RASA control).
            Example: {"B1": "C1,C2,C3"} (ligand at B1) or "B1" for all its atoms
        select_exposed (str or dict): Atoms that should be solvent-exposed (RASA control).
            Example: {"B1": "O1,O2"} or "B1" for all atoms
        select_hbond_donor (dict): Hydrogen bond donor specification.
            Dict mapping <chain><resnum> to donor atoms. Example: {"B1": "N1,N2"}
        select_hbond_acceptor (dict): Hydrogen bond acceptor specification.
            Dict mapping <chain><resnum> to acceptor atoms. Example: {"B1": "O1,O2"}
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
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        repo_dir = folders.get("RFdiffusion3", "")
        biopipelines = folders.get("biopipelines", "")
        env_check = cls._env_exists_check("foundry", env_manager)
        skip = "" if force_reinstall else f"""# Check if already fully installed (env + non-empty checkpoint dir)
if {env_check} && [ -d "{repo_dir}" ] && [ -n "$(ls -A "{repo_dir}" 2>/dev/null)" ]; then
    echo "RFdiffusion3 already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block("foundry", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("foundry", env_manager, biopipelines)
        return f"""echo "=== Installing RFdiffusion3 (foundry) ==="
{skip}{remove_block}
# Create foundry env (skip if it already exists)
if ! {env_check}; then
    {env_block}
else
    echo "foundry environment already exists, skipping creation."
fi

# Download model weights (skip if checkpoint dir already populated)
mkdir -p {repo_dir}
if [ -z "$(ls -A "{repo_dir}" 2>/dev/null)" ]; then
    {env_manager} run -n foundry foundry install rfd3 --checkpoint-dir {repo_dir}
else
    echo "Checkpoint dir {repo_dir} already populated, skipping weight download."
fi

# Verify installation
if [ -n "$(ls -A "{repo_dir}" 2>/dev/null)" ] && {env_manager} run -n foundry python -c "import foundry" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== RFdiffusion3 installation complete ==="
else
    echo "ERROR: RFdiffusion3 verification failed (checkpoints missing or foundry not importable)"
    exit 1
fi
"""

    # Lazy path descriptors — routed through the canonical sub-layout:
    #   configuration/  — input JSONs (tool config, input DataStreams)
    #   execution/      — raw_output/ goes here (foundry CIF.gz dumps)
    #   structures/     — <id>.pdb + structures_map.csv
    #   sequences/      — sequences.csv (content-bearing stream: map == content)
    #   tables/         — standalone TableInfo CSVs (main, metrics, specs)
    json_file = Path(lambda self: self.configuration_path(f"{self._get_prefix()}_rfd3_input.json"))
    # Config-time template entry (shared design params); the runtime builder
    # expands it into one keyed entry per input PDB, writing json_file.
    json_template = Path(lambda self: self.configuration_path(f"{self._get_prefix()}_rfd3_template.json"))
    build_inputs_py = Path(lambda self: self.pipe_script_path("pipe_rfdiffusion3_build_inputs.py"))
    main_table = Path(lambda self: self.table_path("structures"))
    metrics_csv = Path(lambda self: self.table_path("metrics"))
    specifications_csv = Path(lambda self: self.table_path("specifications"))
    sequences_csv = Path(lambda self: self.stream_path("sequences", "sequences.csv"))
    # Foundry dumps CIF.gz files directly into execution/ — that's what the
    # folder is for (raw model dumps). No extra subfolder needed.
    raw_output_folder = Path(lambda self: self.execution_folder)
    rfd_log_file = Path(lambda self: os.path.join(
        os.path.dirname(self.output_folder), "Logs",
        f"{os.path.basename(self.output_folder).split('_')[0] if '_' in os.path.basename(self.output_folder) else '000'}_RFdiffusion3.log"
    ))
    checkpoint_dir = Path(lambda self: self.folders["RFdiffusion3"])
    table_py_file = Path(lambda self: self.pipe_script_path("pipe_rfdiffusion3_table.py"))
    postprocess_py_file = Path(lambda self: self.pipe_script_path("pipe_rfdiffusion3_postprocess.py"))
    pdb_ds_json = Path(lambda self: self.configuration_path("input_structures.json"))
    ligand_json = Path(lambda self: self.configuration_path("input_ligand.json"))
    update_map_py = Path(lambda self: self.pipe_script_path("pipe_update_structures_map.py"))

    def __init__(self,
                 contig: str = "",
                 length: Union[str, int] = None,
                 pdb: Optional[Union[DataStream, StandardizedOutput]] = None,
                 ligand: Optional[Union[str, DataStream, StandardizedOutput]] = None,
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
            ligand: Ligand as a compounds stream (Ligand(code="LIG") or any
                    compounds-producing tool). The residue `code` is read from
                    the stream's `code` column at runtime. If the source also
                    exposes a `structures` stream (e.g. a Ligand with a bound
                    PDB), that structure is used as the input PDB.
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

        Output:
            Streams: structures (.pdb), sequences (.csv)
            Tables:
                structures: id | design | model | pdb | fixed | designed | contig | length | time | status
                metrics: id | design | model | max_ca_deviation | n_chainbreaks | ligand_clashes | ligand_min_distance | loop_fraction | helix_fraction | sheet_fraction | radius_of_gyration | ...
                specifications: id | design | model | sampled_contig | num_tokens_in | num_residues_in | num_chains | num_atoms | num_residues
                sequences: id | source_id | source_pdb | chain | sequence | length
        """
        # Resolve PDB input — store stream for runtime resolution. All input
        # PDBs are iterated at execution time (one foundry design entry per id);
        # never index ids[0] at config time (breaks under lazy IDs).
        self.pdb_stream: Optional[DataStream] = None

        # Ligand — a compounds stream supplying the residue `code` (resolved at
        # runtime). If the source also exposes a structures stream (a Ligand
        # with a bound PDB), those structures become the input PDBs.
        # A bare string is shorthand for an internal Ligand(code=...).
        self.ligand_stream: Optional[DataStream] = resolve_basic_input(
            ligand, Ligand, "compounds", "code")
        if ligand is not None and isinstance(ligand, StandardizedOutput):
            lig_structures = ligand.streams.structures
            if lig_structures and len(lig_structures) > 0:
                self.pdb_stream = lig_structures
        # Handle PDB input (only if a ligand structure wasn't provided)
        if self.pdb_stream is None and pdb is not None:
            if isinstance(pdb, StandardizedOutput):
                self.pdb_stream = pdb.streams.structures
            elif isinstance(pdb, DataStream):
                self.pdb_stream = pdb
            else:
                raise ValueError(f"pdb must be DataStream or StandardizedOutput, got {type(pdb)}")

        # Store parameters
        self.contig = contig
        self.length = length
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
        """Design key for the de-novo (no-PDB) single-entry case.

        With PDB input, each design key is its input pdb id (resolved at
        runtime), so this is only the fallback name for the no-PDB case and for
        config-time artifact filenames (e.g. the template JSON).
        """
        if self.prefix:
            return self.prefix
        return self.pipeline_name

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

        _validate_freeform_string("contig", self.contig)
        _validate_freeform_string("prefix", self.prefix)
        if isinstance(self.length, str):
            _validate_freeform_string("length", self.length)

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

    def _build_json_template(self) -> Dict[str, Any]:
        """Build the shared per-design template entry (one flat dict).

        This holds every design parameter EXCEPT the per-PDB `input` path and
        the `ligand` code, which the runtime builder fills in per entry (input
        resolved per input pdb id; ligand code broadcast). It is a single entry
        dict, not the keyed ``{key: {...}}`` config — the keying happens at
        runtime, one key per input PDB.
        """
        entry: Dict[str, Any] = {}

        if self.contig:
            entry["contig"] = self.contig
        if self.length is not None:
            entry["length"] = str(self.length)
        if self.select_hotspots:
            entry["select_hotspots"] = self._format_hotspots()
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

        return entry

    def get_config_display(self) -> List[str]:
        """Get RFdiffusion3 configuration display lines."""
        config_lines = super().get_config_display()

        # Input information
        if self.pdb_stream:
            config_lines.append(f"INPUT PDB: {', '.join(self.pdb_stream.ids)}")
        else:
            config_lines.append("INPUT: De novo design")

        if self.contig:
            config_lines.append(f"CONTIG: {self.contig}")

        if self.length:
            config_lines.append(f"LENGTH: {self.length}")

        if self.ligand_stream is not None:
            config_lines.append("LIGAND CODE: (resolved from compounds stream at runtime)")

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
        """Generate the bash section that materializes the inputs JSON.

        Advanced (json_config) mode writes the user's full keyed config verbatim
        at config time. Otherwise a config-time template entry is written and a
        runtime builder expands it into one keyed entry per input PDB (input
        paths resolved per id, ligand code broadcast) — replacing the old single
        -entry + jq-patch approach so multiple PDBs each get their own entry.
        """
        # Advanced override: a full {key: {...}} config, passed through as-is.
        if self.json_config:
            cfg = self.json_config if isinstance(self.json_config, dict) else json.loads(self.json_config)
            with open(self.json_file, 'w') as f:
                json.dump(cfg, f, indent=2)
            return f"""echo "Using RFdiffusion3 JSON configuration (advanced): {self.json_file}"

"""

        # Write the shared template entry at config time.
        with open(self.json_template, 'w') as f:
            json.dump(self._build_json_template(), f, indent=2)

        builder_args = f'--template "{self.json_template}" --output "{self.json_file}"'
        if self.pdb_stream:
            builder_args += f' --structures-json "{self.pdb_ds_json}"'
        else:
            builder_args += f' --denovo-prefix "{self._get_prefix()}"'
        if self.ligand_stream is not None:
            self.ligand_stream.save_json(self.ligand_json)
            builder_args += f' --ligand-json "{self.ligand_json}"'

        return f"""echo "Building RFdiffusion3 inputs JSON"
python "{self.build_inputs_py}" {builder_args}
echo "Using RFdiffusion3 JSON configuration: {self.json_file}"

"""

    def generate_script_run_rfdiffusion3(self) -> str:
        """Generate RFdiffusion3 execution bash code.

        No ``global_prefix`` is set: each design key in the inputs JSON is its
        own name (the input pdb id in the multi-PDB case), so foundry names
        outputs ``<key>_<design>_model_<model>``.
        """
        return f"""echo "Starting RFdiffusion3"
echo "JSON config: {self.json_file}"
echo "Output folder: {self.output_folder}"
echo "Raw (CIF.gz) output folder: {self.raw_output_folder}"

# execution/ folder already created by the pipeline's layout step.

# Set checkpoint directory
export FOUNDRY_CHECKPOINT_DIRS="{self.checkpoint_dir}"

# Check checkpoint directory exists
if [ ! -d "${{FOUNDRY_CHECKPOINT_DIRS}}" ]; then
    echo "ERROR: RFdiffusion3 checkpoints not found at ${{FOUNDRY_CHECKPOINT_DIRS}}"
    echo "Please ensure checkpoints are installed at ${{FOUNDRY_CHECKPOINT_DIRS}}"
    exit 1
fi

# Pre-import torch and torchvision to avoid circular import error on Colab
python -c "import torch; import torchvision" 2>/dev/null || true

# Run RFdiffusion3 (outputs CIF.gz format to raw folder)
# n_batches: number of designs to generate (default: 1)
# diffusion_batch_size: number of models per design (default: 8)
{self.container_prefix()}rfd3 design \\
    out_dir="{self.raw_output_folder}" \\
    inputs="{self.json_file}" \\
    n_batches={self.num_designs} \\
    diffusion_batch_size={self.num_models}

"""

    def _generate_postprocess_section(self) -> str:
        """Generate bash section to post-process RFdiffusion3 outputs.

        Converts CIF.gz outputs to PDB format and extracts metrics from JSON files.
        Runs in biopipelines environment (has BioPython, pandas).
        """
        structures_dir = self.stream_folder("structures")
        # Switch to the biopipelines env (BioPython + pandas) via the framework
        # helper — NOT a hardcoded `mamba run`, which doesn't exist on Colab
        # (micromamba) and left the CIF->PDB conversion unrun. On Colab this is a
        # no-op (biopipelines deps live in base Python); on cluster/local it
        # activates the env with whatever env_manager is configured.
        activate_biopipelines = self.activate_environment(name="biopipelines")
        # Output structure ids are keyed by the foundry design key (the input
        # pdb id per entry), so no single --prefix is passed; the postprocess
        # derives each id from the produced filename.
        return f"""echo "Post-processing RFdiffusion3 outputs"

{activate_biopipelines}
# Process CIF.gz files: decompress, convert to PDB, extract metrics and sequences.
# --output_folder is the PDB destination (structures/ stream folder).
python "{self.postprocess_py_file}" \\
    --raw_folder "{self.raw_output_folder}" \\
    --output_folder "{structures_dir}" \\
    --num_designs {self.num_designs} \\
    --num_models {self.num_models} \\
    --design_startnum {self.design_startnum} \\
    --metrics_csv "{self.metrics_csv}" \\
    --specifications_csv "{self.specifications_csv}" \\
    --sequences_csv "{self.sequences_csv}"

"""

    def generate_script_create_table(self) -> str:
        """Generate table creation bash code.

        The table builder scans the structures/ folder and matches each id back
        to its JSON design key, so it needs no prefix/num_designs.
        """
        return f"""echo "Creating results table"
python "{self.table_py_file}" \\
    --output_folder "{self.output_folder}" \\
    --json_file "{self.json_file}" \\
    --table_path "{self.main_table}" \\
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
        # Serialize input DataStream to JSON for runtime file resolution
        if self.pdb_stream:
            self.pdb_stream.save_json(self.pdb_ds_json)

        script_content = "#!/bin/bash\n"
        script_content += "# RFdiffusion3 execution script\n"
        script_content += "# Generated by BioPipelines\n\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self._generate_json_section()
        script_content += self.generate_script_run_rfdiffusion3()
        script_content += self._generate_postprocess_section()
        script_content += self.generate_script_create_table()
        script_content += self._generate_script_update_structures_map()
        script_content += self.generate_completion_check_footer()

        return script_content

    def _generate_script_update_structures_map(self) -> str:
        """Generate script to write structures_map.csv from the actual runtime PDBs."""
        structures_map = self.stream_map_path("structures")
        structures_dir = self.stream_folder("structures")
        # Each id is "<pdb_id>_d<D>_m<M>" (num_models>1) or "<pdb_id>_<D>";
        # strip that design/model suffix to recover the parent `structures.id`.
        # The regex is RFd3's id shape, passed to the tool-agnostic helper.
        if self.pdb_stream:
            suffix_re = r"_d\d+_m\d+$" if self.num_models > 1 else r"_\d+$"
            prov_arg = (f' --provenance-from-suffix "structures.id"'
                        f' --suffix-regex "{suffix_re}"')
        else:
            prov_arg = ""
        return f"""echo "Writing structures map from actual output files"
python {self.update_map_py} --structures-map "{structures_map}" --output-folder "{structures_dir}"{prov_arg}

"""

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
        # Generate pattern-based IDs. The base is each design key — the input
        # pdb id(s) in the multi-PDB case, or the de-novo prefix otherwise.
        d_start = self.design_startnum
        d_end = self.design_startnum + self.num_designs - 1
        m_start = self.design_startnum
        m_end = self.design_startnum + self.num_models - 1

        # Suffix appended after the helper's "_" separator. Keep parent ids
        # compact (no config-time expansion) via generate_multiplied_ids_pattern.
        if self.num_models > 1:
            suffix_pattern = f"d<{d_start}..{d_end}>_m<{m_start}..{m_end}>"
        else:
            suffix_pattern = f"<{d_start}..{d_end}>"

        base_keys = self.pdb_stream.ids if self.pdb_stream else [self._get_prefix()]
        structure_ids = generate_multiplied_ids_pattern(
            base_keys, suffix_pattern, input_stream_name="structures"
        )
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

        # Define table structures
        tables = {
            "structures": TableInfo(
                name="structures",
                path=self.main_table,
                columns=["id", "design", "model", "pdb", "fixed", "designed", "contig", "length", "time", "status"],
                description="RFdiffusion3 structure generation results with fixed/designed regions"
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
                description="RFdiffusion3 quality metrics extracted from JSON outputs"
            ),
            "specifications": TableInfo(
                name="specifications",
                path=self.specifications_csv,
                columns=[
                    "id", "design", "model", "sampled_contig", "num_tokens_in", "num_residues_in",
                    "num_chains", "num_atoms", "num_residues"
                ],
                description="RFdiffusion3 design specifications and statistics"
            ),
            "sequences": TableInfo(
                name="sequences",
                path=self.sequences_csv,
                columns=["id", "source_id", "source_pdb", "chain", "sequence", "length"],
                description="RFdiffusion3 designed protein sequences extracted from PDB"
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
                "ligand_ids": list(self.ligand_stream.ids) if self.ligand_stream else [],
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
                "pdb_input_ids": self.pdb_stream.ids if self.pdb_stream else None
            }
        })
        return base_dict
