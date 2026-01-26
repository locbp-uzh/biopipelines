"""
BoltzGen configuration for protein binder design.

Generates protein binders (proteins, peptides, or nanobodies) targeting specified molecules
using diffusion models, inverse folding, and iterative structure validation.
"""

import os
import json
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo


class BoltzGen(BaseConfig):
    """
    BoltzGen configuration for protein binder design.

    Generates and validates protein binders for target molecules using:
    - Diffusion-based backbone generation
    - Inverse folding for sequence design
    - Structure prediction for validation
    - Multi-metric filtering and ranking

    Commonly used for:
    - De novo binder design against protein targets
    - Peptide binder generation
    - Nanobody design
    - Protein-small molecule binder design
    """

    # Tool identification
    TOOL_NAME = "BoltzGen"
    DEFAULT_ENV = None  # Loaded from config.yaml

    def __init__(self,
                 # Design specification - Option 1: Manual YAML/dict
                 design_spec: Union[str, Dict[str, Any]] = None,
                 # Design specification - Option 2: Automatic building from BioPipelines inputs
                 # For protein-small_molecule: use ligand parameter
                 # For protein-anything: use target_structure parameter
                 ligand: Union[ToolOutput, StandardizedOutput] = None,
                 target_structure: Union[str, ToolOutput, StandardizedOutput] = None,
                 ligand_code: str = None,
                 binder_spec: Union[str, Dict[str, Any]] = None,
                 binding_region: str = None,
                 # Output configuration
                 protocol: str = "protein-anything",
                 num_designs: int = 10000,
                 budget: int = 100,
                 # Design generation parameters
                 design_checkpoints: Optional[List[str]] = None,
                 step_scale: Optional[float] = None,
                 noise_scale: Optional[float] = None,
                 diffusion_batch_size: Optional[int] = None,
                 # Inverse folding parameters
                 inverse_fold_num_sequences: int = 1,
                 skip_inverse_folding: bool = False,
                 # Filtering parameters
                 alpha: float = None,
                 filter_biased: bool = True,
                 refolding_rmsd_threshold: Optional[float] = None,
                 additional_filters: Optional[List[str]] = None,
                 metrics_override: Optional[Dict[str, float]] = None,
                 # Execution parameters
                 devices: Optional[int] = None,
                 reuse: Union[ToolOutput, StandardizedOutput, str, None] = None,
                 steps: Optional[List[str]] = None,
                 cache_dir: Optional[str] = None,
                 **kwargs):
        """
        Initialize BoltzGen configuration.

        You can provide EITHER design_spec (manual) OR automatic mode inputs.
        For automatic mode, use ONE of:
        - ligand + binder_spec: For protein-small_molecule protocol (binder against small molecule)
        - target_structure + binder_spec: For protein-anything protocol (binder against protein target)

        Args:
            design_spec: YAML configuration string/dict or path to YAML file defining:
                        - Target entities (proteins, ligands)
                        - Binder specification (sequence ranges)
                        - Binding site constraints
                        - Secondary structure specifications
            ligand: Ligand input from Ligand tool (ToolOutput/StandardizedOutput).
                   For protein-small_molecule protocol. The tool reads the compounds.csv
                   to determine format (smiles/ccd) automatically.
            target_structure: BioPipelines structure input (ToolOutput or file path) for
                            protein-anything protocol (binder against protein target)
            ligand_code: Ligand residue name in target structure (e.g., "LIG", "ATP")
                        for binding site detection when using target_structure
            binder_spec: Binder specification for automatic building:
                        - String format: "50-100" for length range, or "50" for fixed length
                        - Dict format: {"length": [50, 100]} or custom specification
            binding_region: Target binding region (optional, PyMOL-style: "A:9-140")
                          If not provided and ligand_code is given, will auto-detect from ligand proximity
            protocol: Design protocol - one of:
                     - "protein-anything": General protein binder design
                     - "peptide-anything": Peptide binder design
                     - "protein-small_molecule": Protein-small molecule binder design
                     - "nanobody-anything": Nanobody design
            num_designs: Total number of designs to generate (default: 10000)
            budget: Final number of designs after diversity filtering (default: 100)
            design_checkpoints: Model checkpoint paths for backbone generation
            step_scale: Diffusion sampling step scale parameter
            noise_scale: Diffusion sampling noise scale parameter
            diffusion_batch_size: Samples per trunk run (auto if None)
            inverse_fold_num_sequences: Sequences per backbone (default: 1)
            skip_inverse_folding: Skip sequence redesign step (default: False)
            alpha: Quality/diversity in sequence during filtering (0.0=quality, 1.0=diversity, default: 0.01 peptides and 0.001 others)
            filter_biased: Remove amino acid composition outliers (default: True)
            refolding_rmsd_threshold: RMSD threshold for filtering designs based on refolding (e.g., 2.5)
            additional_filters: Hard threshold expressions (e.g., ["design_ALA>0.3"])
            metrics_override: Per-metric ranking weights override
            devices: Number of GPUs to use (default: auto-detect)
            reuse: Previous BoltzGen run to resume or extend. Can be:
                  - ToolOutput/StandardizedOutput from a previous BoltzGen call
                  - String path to a previous BoltzGen output directory
                  When provided, operates on the existing output folder instead of
                  creating a new one. Use with `steps` to run specific steps only
                  (e.g., steps=["analysis"] or steps=["filtering"]).
            steps: Run only specific pipeline steps (default: all steps)
            cache_dir: Model/cache location (default: /home/$USER/data/boltzgen, auto-downloaded on first run)
            **kwargs: Additional parameters passed to BaseConfig
        """
        # Store design specification parameters
        self.design_spec = design_spec
        self.ligand = ligand
        self.target_structure = target_structure
        self.ligand_code = ligand_code
        self.binder_spec = binder_spec
        self.binding_region = binding_region

        # Store output configuration
        self.protocol = protocol
        self.num_designs = num_designs
        self.budget = budget

        # Design generation parameters
        self.design_checkpoints = design_checkpoints
        self.step_scale = step_scale
        self.noise_scale = noise_scale
        self.diffusion_batch_size = diffusion_batch_size

        # Inverse folding parameters
        self.inverse_fold_num_sequences = inverse_fold_num_sequences
        self.skip_inverse_folding = skip_inverse_folding

        # Filtering parameters
        self.alpha = alpha
        self.filter_biased = filter_biased
        self.refolding_rmsd_threshold = refolding_rmsd_threshold
        self.additional_filters = additional_filters or []
        self.metrics_override = metrics_override

        # Execution parameters
        self.devices = devices
        self.reuse = reuse
        self.reuse_path = None  # Will be set in configure_inputs if reuse is provided
        self.steps = steps or []
        self.cache_dir = cache_dir

        # Track inputs from previous tools
        self.ligand_compounds_csv = None  # Path to compounds.csv from Ligand tool
        self.target_structure_path = None  # Path to target protein structure

        # Determine specification mode
        if reuse is not None:
            # Reuse mode: operating on existing BoltzGen output
            self.spec_mode = "reuse"
            self.design_spec_is_dict = False
            self.design_spec_is_yaml_str = False
            self.design_spec_is_file = False
        elif design_spec is not None:
            # Manual mode: user provided design_spec
            self.spec_mode = "manual"
            self.design_spec_is_dict = isinstance(design_spec, dict)
            self.design_spec_is_yaml_str = isinstance(design_spec, str) and (
                'entities:' in design_spec or '- protein:' in design_spec
            )
            self.design_spec_is_file = isinstance(design_spec, str) and not self.design_spec_is_yaml_str
        elif ligand is not None and binder_spec is not None:
            # Automatic mode: ligand + binder_spec (protein-small_molecule)
            self.spec_mode = "ligand"
            self.design_spec_is_dict = False
            self.design_spec_is_yaml_str = False
            self.design_spec_is_file = False
        elif target_structure is not None and binder_spec is not None:
            # Automatic mode: target_structure + binder_spec (protein-anything)
            self.spec_mode = "target"
            self.design_spec_is_dict = False
            self.design_spec_is_yaml_str = False
            self.design_spec_is_file = False
        else:
            self.spec_mode = None  # Will fail validation

        # Initialize base class
        super().__init__(**kwargs)

        # Initialize file paths (will be set in configure_inputs)
        self._initialize_file_paths()

    def _initialize_file_paths(self):
        """Initialize common file paths used throughout the class."""
        # Design specification files
        self.design_spec_yaml_file = None

        # Output directories
        self.config_folder = None
        self.intermediate_designs_folder = None
        self.intermediate_inverse_folded_folder = None
        self.final_ranked_folder = None

        # Output files
        self.all_designs_metrics_csv = None
        self.final_designs_metrics_csv = None
        self.aggregate_metrics_csv = None
        self.per_target_metrics_csv = None
        self.results_overview_pdf = None

        # Helper script paths
        self.boltzgen_helper_py = None

    def _setup_file_paths(self):
        """Set up all file paths after output_folder is known."""
        # Design specification file
        self.design_spec_yaml_file = os.path.join(self.output_folder, "design_spec.yaml")

        # Output directory structure (mirrors BoltzGen output)
        self.config_folder = os.path.join(self.output_folder, "config")
        self.intermediate_designs_folder = os.path.join(self.output_folder, "intermediate_designs")
        self.intermediate_inverse_folded_folder = os.path.join(
            self.output_folder, "intermediate_designs_inverse_folded"
        )
        self.final_ranked_folder = os.path.join(self.output_folder, "final_ranked_designs")

        # Output files
        self.all_designs_metrics_csv = os.path.join(
            self.final_ranked_folder, "all_designs_metrics.csv"
        )
        self.final_designs_metrics_csv = os.path.join(
            self.final_ranked_folder, f"final_designs_metrics_{self.budget}.csv"
        )
        self.aggregate_metrics_csv = os.path.join(
            self.intermediate_inverse_folded_folder, "aggregate_metrics_analyze.csv"
        )
        self.per_target_metrics_csv = os.path.join(
            self.intermediate_inverse_folded_folder, "per_target_metrics_analyze.csv"
        )
        self.results_overview_pdf = os.path.join(
            self.final_ranked_folder, "results_overview.pdf"
        )

        # Helper script paths and cache directory (only set if folders are available)
        if hasattr(self, 'folders') and self.folders:
            self.boltzgen_helper_py = os.path.join(
                self.folders["HelpScripts"], "pipe_boltzgen.py"
            )
            # Use tool_data folder for cache if user didn't specify custom cache_dir
            if self.cache_dir is None:
                tool_data_folder = self.folders.get("tool_data", {})
                if isinstance(tool_data_folder, dict):
                    cache_base = tool_data_folder.get("BoltzGen", "boltzgen")
                else:
                    cache_base = "boltzgen"
                self.cache_dir = os.path.join(
                    os.path.expanduser("~"),
                    "data",
                    cache_base
                )
        else:
            self.boltzgen_helper_py = None

    def validate_params(self):
        """Validate BoltzGen parameters."""
        # Validate specification mode
        if self.spec_mode is None:
            raise ValueError(
                "Must provide EITHER design_spec (manual) OR "
                "(ligand + binder_spec) for protein-small_molecule OR "
                "(target_structure + binder_spec) for protein-anything OR "
                "reuse (to continue from existing run)"
            )

        if self.spec_mode == "reuse":
            # Validate reuse parameter
            if not isinstance(self.reuse, (ToolOutput, StandardizedOutput, str)):
                raise ValueError(
                    f"reuse must be ToolOutput, StandardizedOutput, or path string, "
                    f"got {type(self.reuse)}"
                )
            # Steps should typically be specified when reusing
            if not self.steps:
                # Allow full run with reuse (will skip completed steps)
                pass

        if self.spec_mode == "manual" and not self.design_spec:
            raise ValueError("design_spec cannot be empty in manual mode")

        if self.spec_mode == "ligand":
            if not self.ligand:
                raise ValueError("ligand is required for protein-small_molecule mode")
            if not self.binder_spec:
                raise ValueError("binder_spec is required for protein-small_molecule mode")
            # Validate ligand is ToolOutput or StandardizedOutput
            if not isinstance(self.ligand, (ToolOutput, StandardizedOutput)):
                raise ValueError(
                    f"ligand must be ToolOutput or StandardizedOutput from Ligand tool, "
                    f"got {type(self.ligand)}"
                )

        if self.spec_mode == "target":
            if not self.target_structure:
                raise ValueError("target_structure is required for protein-anything mode")
            if not self.binder_spec:
                raise ValueError("binder_spec is required for protein-anything mode")

        # Validate protocol
        valid_protocols = [
            "protein-anything", "peptide-anything",
            "protein-small_molecule", "nanobody-anything"
        ]
        if self.protocol not in valid_protocols:
            raise ValueError(
                f"protocol must be one of {valid_protocols}, got '{self.protocol}'"
            )

        # Validate numeric parameters
        if self.num_designs <= 0:
            raise ValueError("num_designs must be positive")

        if self.budget <= 0:
            raise ValueError("budget must be positive")

        if self.budget > self.num_designs:
            self.budget = self.num_designs
            print("[Warning] budget was set equal to num_designs.")

        if self.inverse_fold_num_sequences <= 0:
            raise ValueError("inverse_fold_num_sequences must be positive")

        if self.alpha and not (0.0 <= self.alpha <= 1.0):
            raise ValueError("alpha must be between 0.0 and 1.0")

        # Validate diffusion parameters if provided
        if self.step_scale is not None and self.step_scale <= 0:
            raise ValueError("step_scale must be positive")

        if self.noise_scale is not None and self.noise_scale <= 0:
            raise ValueError("noise_scale must be positive")

        if self.diffusion_batch_size is not None and self.diffusion_batch_size <= 0:
            raise ValueError("diffusion_batch_size must be positive")

        if self.devices is not None and self.devices <= 0:
            raise ValueError("devices must be positive")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input design specification and dependencies."""
        self.folders = pipeline_folders

        # Handle reuse mode first - extract path and override output_folder
        if self.spec_mode == "reuse":
            if isinstance(self.reuse, StandardizedOutput):
                if hasattr(self.reuse, 'output_folder') and self.reuse.output_folder:
                    self.reuse_path = self.reuse.output_folder
                elif hasattr(self.reuse, 'tool_folder') and self.reuse.tool_folder:
                    self.reuse_path = self.reuse.tool_folder
                else:
                    raise ValueError("StandardizedOutput has no output_folder or tool_folder")
            elif isinstance(self.reuse, ToolOutput):
                output_folder = self.reuse.get_output_files("output_folder")
                if output_folder:
                    self.reuse_path = output_folder[0] if isinstance(output_folder, list) else output_folder
                    self.dependencies.append(self.reuse.config)
                else:
                    raise ValueError("ToolOutput has no output_folder")
            elif isinstance(self.reuse, str):
                self.reuse_path = self.reuse
            else:
                raise ValueError(f"Invalid reuse type: {type(self.reuse)}")

            # Override output_folder with the reuse path
            self.output_folder = self.reuse_path

        self._setup_file_paths()

        if self.spec_mode == "manual":
            # Handle manual design specification input
            if self.design_spec_is_dict:
                # Dictionary specification - will convert to YAML in script generation
                pass
            elif self.design_spec_is_file:
                # File path - check if it exists
                if os.path.isabs(self.design_spec):
                    spec_path = self.design_spec
                else:
                    # Try relative to project folders
                    spec_path = os.path.join(pipeline_folders.get("notebooks", "."), self.design_spec)

                if not os.path.exists(spec_path):
                    raise ValueError(f"Design specification file not found: {self.design_spec}")

                self.design_spec = spec_path
            elif self.design_spec_is_yaml_str:
                # Direct YAML string - will write to file in script generation
                pass
            else:
                raise ValueError(
                    "design_spec must be a YAML string, dictionary, or path to YAML file"
                )

        elif self.spec_mode == "ligand":
            # Handle ligand mode: extract compounds CSV from Ligand ToolOutput
            if isinstance(self.ligand, StandardizedOutput):
                # Get compounds CSV path
                if hasattr(self.ligand, 'compounds') and self.ligand.compounds:
                    self.ligand_compounds_csv = self.ligand.compounds[0]
                elif hasattr(self.ligand, 'tables') and hasattr(self.ligand.tables, 'compounds'):
                    self.ligand_compounds_csv = self.ligand.tables.compounds
                else:
                    raise ValueError("Ligand StandardizedOutput has no compounds attribute")
            elif isinstance(self.ligand, ToolOutput):
                # Get compounds CSV path from ToolOutput
                compounds = self.ligand.get_output_files("compounds")
                if compounds:
                    self.ligand_compounds_csv = compounds[0]
                    self.dependencies.append(self.ligand.config)
                else:
                    raise ValueError("Ligand ToolOutput has no compounds output")
            else:
                raise ValueError(
                    f"ligand must be ToolOutput or StandardizedOutput, got {type(self.ligand)}"
                )

        elif self.spec_mode == "target":
            # Handle target mode: extract structure path from ToolOutput
            if isinstance(self.target_structure, (ToolOutput, StandardizedOutput)):
                # Extract first structure from tool output
                if hasattr(self.target_structure, 'structures') and self.target_structure.structures:
                    self.target_structure_path = self.target_structure.structures[0]
                elif hasattr(self.target_structure, 'get') and 'structures' in self.target_structure:
                    structures = self.target_structure.get('structures', [])
                    if structures:
                        self.target_structure_path = structures[0]
                    else:
                        raise ValueError("No structures found in target_structure ToolOutput")
                else:
                    raise ValueError("target_structure ToolOutput has no structures attribute")
            elif isinstance(self.target_structure, str):
                # Direct file path
                if not os.path.exists(self.target_structure):
                    raise ValueError(f"Target structure file not found: {self.target_structure}")
                self.target_structure_path = self.target_structure
            else:
                raise ValueError(
                    f"target_structure must be a file path or ToolOutput, got {type(self.target_structure)}"
                )

            # Will build design spec in script generation
            self.design_spec = None

    def get_config_display(self) -> List[str]:
        """Get BoltzGen configuration display lines."""
        config_lines = super().get_config_display()

        # Design specification display
        if self.spec_mode == "manual":
            if self.design_spec_is_dict:
                spec_display = f"Manual (Dictionary, {len(self.design_spec)} keys)"
            elif self.design_spec_is_file:
                spec_display = f"Manual (File: {os.path.basename(self.design_spec)})"
            else:
                spec_display = "Manual (YAML string)"
        elif self.spec_mode == "ligand":
            spec_display = "Automatic (from ligand + binder_spec)"
            if self.ligand_compounds_csv:
                spec_display += f"\n  Ligand CSV: {os.path.basename(self.ligand_compounds_csv)}"
            if self.binder_spec:
                spec_display += f"\n  Binder: {self.binder_spec}"
            if self.binding_region:
                spec_display += f"\n  Binding region: {self.binding_region}"
        elif self.spec_mode == "target":
            spec_display = "Automatic (from target_structure + binder_spec)"
            if self.target_structure_path:
                spec_display += f"\n  Target: {os.path.basename(self.target_structure_path)}"
            if self.ligand_code:
                spec_display += f"\n  Ligand code: {self.ligand_code}"
            if self.binder_spec:
                spec_display += f"\n  Binder: {self.binder_spec}"
            if self.binding_region:
                spec_display += f"\n  Binding region: {self.binding_region}"
        elif self.spec_mode == "reuse":
            spec_display = "Reuse (from previous BoltzGen run)"
            if self.reuse_path:
                spec_display += f"\n  Source: {os.path.basename(self.reuse_path)}"
            if self.steps:
                spec_display += f"\n  Steps: {', '.join(self.steps)}"
        else:
            spec_display = "Unknown mode"

        config_lines.extend([
            f"DESIGN SPEC: {spec_display}",
            f"PROTOCOL: {self.protocol}",
            f"NUM DESIGNS: {self.num_designs}",
            f"BUDGET: {self.budget}",
            f"INVERSE FOLD SEQUENCES: {self.inverse_fold_num_sequences}",
            f"ALPHA (quality/diversity): {self.alpha}",
            f"FILTER BIASED: {self.filter_biased}"
        ])

        if self.additional_filters:
            config_lines.append(f"ADDITIONAL FILTERS: {', '.join(self.additional_filters)}")

        if self.skip_inverse_folding:
            config_lines.append("SKIP INVERSE FOLDING: True")

        if self.steps:
            config_lines.append(f"STEPS: {', '.join(self.steps)}")

        if self.devices:
            config_lines.append(f"DEVICES: {self.devices}")

        return config_lines

    def _build_automatic_yaml_spec(self) -> Dict[str, Any]:
        """
        Build BoltzGen YAML specification from BioPipelines inputs.

        Returns:
            Dictionary representing the design specification
        """
        # Parse binder_spec
        if isinstance(self.binder_spec, dict):
            binder_length = self.binder_spec.get("length", [50, 100])
        elif isinstance(self.binder_spec, str):
            # Parse string format: "50-100" or "50"
            if "-" in self.binder_spec:
                parts = self.binder_spec.split("-")
                binder_length = [int(parts[0]), int(parts[1])]
            else:
                length_val = int(self.binder_spec)
                binder_length = [length_val, length_val]
        else:
            raise ValueError(f"Invalid binder_spec format: {self.binder_spec}")

        # Build entities section
        entities = []

        # Add target structure
        entities.append({
            "file": {
                "path": self.target_structure_path,
                "include": [{"chain": {"id": "A"}}]  # Default to chain A
            }
        })

        # Add binder
        entities.append({
            "protein": {
                "id": "binder",
                "sequence": f"{binder_length[0]}..{binder_length[1]}"
            }
        })

        # Build the full specification
        spec = {"entities": entities}

        # Add bindings section if binding region is specified
        if self.binding_region:
            # Parse binding region format: "A:9-140" -> target_chains: [A], target_region: "9-140"
            if ":" in self.binding_region:
                chain, region = self.binding_region.split(":", 1)
                spec["bindings"] = [{
                    "binder_chain": "binder",
                    "target_chains": [chain],
                    "target_region": region
                }]

        return spec

    def generate_script(self, script_path: str) -> str:
        """
        Generate BoltzGen execution script.

        Args:
            script_path: Path where script should be written

        Returns:
            Script content as string
        """
        # Create output directories
        os.makedirs(self.output_folder, exist_ok=True)
        os.makedirs(self.config_folder, exist_ok=True)

        # Parse binder spec
        binder_min, binder_max = self._parse_binder_spec()

        # Generate config generation section (all modes generate YAML at runtime)
        config_generation_section = self._generate_config_section(binder_min, binder_max)

        # Build BoltzGen command
        boltzgen_cmd = self._build_boltzgen_command()

        # Determine if postprocessing is needed (only for filtering step)
        # If no steps specified, all steps run including filtering
        has_filtering = "filtering" in self.steps if self.steps else True

        # Build postprocessing section
        if has_filtering:
            postprocessing_section = f"""
# Run postprocessing to extract sequences and validate outputs
echo "Running BoltzGen postprocessing..."
python "{self.boltzgen_helper_py}" \\
  --output_folder "{self.output_folder}" \\
  --budget {self.budget} \\
  --extract_sequences \\
  --validate

if [ $? -eq 0 ]; then
    echo "Postprocessing completed successfully"
else
    echo "Warning: Postprocessing had errors (continuing anyway)"
fi
"""
        else:
            postprocessing_section = ""

        # Generate script content
        script_content = f"""#!/bin/bash
# BoltzGen execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Running BoltzGen binder design"
echo "Protocol: {self.protocol}"
echo "Output folder: {self.output_folder}"

{config_generation_section}

echo "Design specification: {self.design_spec_yaml_file}"

# Run BoltzGen
{boltzgen_cmd}

if [ $? -eq 0 ]; then
    echo "BoltzGen completed successfully"
else
    echo "Error: BoltzGen failed"
    exit 1
fi
{postprocessing_section}
{self.generate_completion_check_footer()}
"""

        return script_content

    def _parse_binder_spec(self):
        """Parse binder_spec into min and max length values."""
        if isinstance(self.binder_spec, dict):
            # Dict format: {"length": [50, 100]}
            length = self.binder_spec.get("length", [100, 200])
            if isinstance(length, list) and len(length) == 2:
                return length[0], length[1]
            elif isinstance(length, int):
                return length, length
            else:
                return 100, 200
        elif isinstance(self.binder_spec, str):
            # String format: "50-100" or "50"
            if "-" in self.binder_spec:
                parts = self.binder_spec.split("-")
                return int(parts[0]), int(parts[1])
            else:
                length = int(self.binder_spec)
                return length, length
        else:
            return 100, 200

    def _generate_config_section(self, binder_min: int, binder_max: int) -> str:
        """Generate the config file creation section for the bash script."""
        config_py = os.path.join(self.folders["HelpScripts"], "pipe_boltzgen_config.py")

        if self.spec_mode == "reuse":
            # Reuse mode: design_spec.yaml already exists in the reuse folder
            return f'''# Using existing design specification from previous run
echo "Reusing design specification from: {self.reuse_path}"
echo "Design spec: {self.design_spec_yaml_file}"
'''

        if self.spec_mode == "ligand":
            # Ligand mode: generate YAML at runtime from compounds.csv
            target_args = ""
            if self.binding_region:
                target_args += f' --binding-region "{self.binding_region}"'

            return f'''# Generate design specification from ligand
echo "Generating design specification from ligand compounds..."
python "{config_py}" \\
  --ligand-csv "{self.ligand_compounds_csv}" \\
  --output-yaml "{self.design_spec_yaml_file}" \\
  --binder-length-min {binder_min} \\
  --binder-length-max {binder_max}{target_args}

if [ $? -ne 0 ]; then
    echo "Error: Failed to generate design specification"
    exit 1
fi
'''

        elif self.spec_mode == "target":
            # Target mode: generate YAML at runtime with target structure
            target_args = f' --target-structure "{self.target_structure_path}"'
            if self.ligand_code:
                target_args += f' --ligand-code "{self.ligand_code}"'
            if self.binding_region:
                target_args += f' --binding-region "{self.binding_region}"'

            # For target mode we need a ligand CSV - create a dummy one if not provided
            # Or require ligand input for target mode too
            return f'''# Generate design specification from target structure
echo "Generating design specification from target structure..."
# Note: Target mode with protein structure - YAML will be built at runtime
python "{config_py}" \\
  --ligand-csv "NONE" \\
  --output-yaml "{self.design_spec_yaml_file}" \\
  --binder-length-min {binder_min} \\
  --binder-length-max {binder_max}{target_args}

if [ $? -ne 0 ]; then
    echo "Error: Failed to generate design specification"
    exit 1
fi
'''

        elif self.spec_mode == "manual":
            # Manual mode: write YAML at pipeline time (user-provided spec)
            if self.design_spec_is_dict:
                # Convert dictionary to YAML and write to file at pipeline time
                import yaml
                with open(self.design_spec_yaml_file, 'w') as f:
                    yaml.dump(self.design_spec, f, default_flow_style=False)
                return f'# Using pre-generated design specification: {self.design_spec_yaml_file}\n'
            elif self.design_spec_is_yaml_str:
                # Write YAML string to file at pipeline time
                with open(self.design_spec_yaml_file, 'w') as f:
                    f.write(self.design_spec)
                return f'# Using pre-generated design specification: {self.design_spec_yaml_file}\n'
            elif self.design_spec_is_file:
                # Copy existing file to output folder at pipeline time
                import shutil
                shutil.copy(self.design_spec, self.design_spec_yaml_file)
                return f'# Using copied design specification: {self.design_spec_yaml_file}\n'

        return "# Unknown spec mode\n"

    def _build_boltzgen_command(self) -> str:
        """Build the boltzgen run command with all arguments."""
        # Base command
        cmd_parts = [
            "boltzgen run",
            f'"{self.design_spec_yaml_file}"',
            f'--output "{self.output_folder}"',
            f'--protocol {self.protocol}',
            f'--num_designs {self.num_designs}',
            f'--budget {self.budget}'
        ]

        # Design generation parameters
        if self.design_checkpoints:
            checkpoints_str = ' '.join(f'"{cp}"' for cp in self.design_checkpoints)
            cmd_parts.append(f'--design_checkpoints {checkpoints_str}')

        if self.step_scale is not None:
            cmd_parts.append(f'--step_scale {self.step_scale}')

        if self.noise_scale is not None:
            cmd_parts.append(f'--noise_scale {self.noise_scale}')

        if self.diffusion_batch_size is not None:
            cmd_parts.append(f'--diffusion_batch_size {self.diffusion_batch_size}')

        # Inverse folding parameters
        cmd_parts.append(f'--inverse_fold_num_sequences {self.inverse_fold_num_sequences}')

        if self.skip_inverse_folding:
            cmd_parts.append('--skip_inverse_folding')

        # Filtering parameters
        if self.alpha:
            cmd_parts.append(f'--alpha {self.alpha}')

        if not self.filter_biased:
            cmd_parts.append('--filter_biased false')

        if self.refolding_rmsd_threshold is not None:
            cmd_parts.append(f'--refolding_rmsd_threshold {self.refolding_rmsd_threshold}')

        if self.additional_filters:
            # Format: --additional_filters 'feature>threshold' 'feature<threshold'
            filters_str = ' '.join(f"'{f}'" for f in self.additional_filters)
            cmd_parts.append(f"--additional_filters {filters_str}")

        if self.metrics_override:
            # Format: --metrics_override metric_name=weight metric_name=weight
            # A larger weight down-weights that metric's rank
            metrics_str = ' '.join(f"{k}={v}" for k, v in self.metrics_override.items())
            cmd_parts.append(f"--metrics_override {metrics_str}")

        # Execution parameters
        if self.devices is not None:
            cmd_parts.append(f'--devices {self.devices}')

        if self.reuse:
            cmd_parts.append('--reuse')

        if self.steps:
            steps_str = ' '.join(self.steps)
            cmd_parts.append(f'--steps {steps_str}')

        if self.cache_dir:
            cmd_parts.append(f'--cache "{self.cache_dir}"')

        # Join with line continuations for readability
        return ' \\\n  '.join(cmd_parts)

    def get_output_files(self) -> Dict[str, Any]:
        """
        Get expected output files after BoltzGen execution.

        Predictions are step-dependent:
        - "analysis" in steps: aggregate_metrics_analyze.csv, per_target_metrics_analyze.csv
        - "filtering" in steps: structures in final_<budget>_designs/ (rankNNNN_*.cif),
                                all_designs_metrics.csv, final_designs_metrics_<budget>.csv,
                                results_overview.pdf
        - "folding" in steps (but not filtering): structures in refold_cif/ (design_spec_NNN.cif)
        - None of above: empty prediction

        Returns:
            Dictionary mapping output types to file paths with standard keys.
        """
        # Ensure file paths are set up
        if not hasattr(self, 'final_ranked_folder') or self.final_ranked_folder is None:
            self._setup_file_paths()

        # Determine which steps are being run
        # If no steps specified, BoltzGen runs all steps by default
        steps = self.steps if self.steps else [
            "design", "inverse_folding", "folding", "design_folding", "affinity", "analysis", "filtering"
        ]
        has_analysis = "analysis" in steps
        has_filtering = "filtering" in steps
        has_folding = "folding" in steps or "design_folding" in steps

        # Initialize empty outputs
        tables = {}
        structures = []
        structure_ids = []
        output = {
            "structures": structures,
            "structure_ids": structure_ids,
            "tables": tables,
            "output_folder": self.output_folder,
            "tool_folder": self.output_folder,
        }

        # If no relevant steps, return empty prediction
        if not has_analysis and not has_filtering and not has_folding:
            return output

        # Analysis step outputs
        if has_analysis:
            tables["aggregate_metrics"] = TableInfo(
                name="aggregate_metrics",
                path=self.aggregate_metrics_csv,
                columns=[
                    "id", "file_name", "designed_sequence", "designed_chain_sequence",
                    "num_prot_tokens", "num_lig_atoms", "num_resolved_tokens", "num_tokens",
                    "num_design", "UNK_fraction", "GLY_fraction", "ALA_fraction", "CYS_fraction",
                    "SER_fraction", "PRO_fraction", "THR_fraction", "VAL_fraction", "ILE_fraction",
                    "ASN_fraction", "ASP_fraction", "LEU_fraction", "MET_fraction", "GLN_fraction",
                    "GLU_fraction", "LYS_fraction", "HIS_fraction", "PHE_fraction", "ARG_fraction",
                    "TYR_fraction", "TRP_fraction", "loop", "helix", "sheet", "liability_score",
                    "liability_num_violations", "liability_high_severity_violations",
                    "liability_medium_severity_violations", "liability_low_severity_violations",
                    "liability_AspBridge_count", "liability_AspCleave_count", "liability_ProtTryp_count",
                    "liability_TrpOx_count", "liability_HydroPatch_count", "liability_UnpairedCys_count",
                    "native_rmsd", "native_rmsd_bb", "native_rmsd_refolded", "native_rmsd_bb_refolded",
                    "bb_rmsd", "bb_rmsd_design", "bb_rmsd_target", "design_ptm", "design_iptm",
                    "design_to_target_iptm", "min_design_to_target_pae", "min_interaction_pae",
                    "affinity_pred_value", "affinity_probability_binary1"
                ],
                description="Aggregate metrics for all designs from analysis step",
                count=None
            )
            tables["per_target_metrics"] = TableInfo(
                name="per_target_metrics",
                path=self.per_target_metrics_csv,
                columns=[
                    "target_id", "num_prot_tokens", "num_lig_atoms", "num_resolved_tokens",
                    "num_tokens", "num_design", "UNK_fraction", "GLY_fraction", "ALA_fraction",
                    "CYS_fraction", "SER_fraction", "PRO_fraction", "THR_fraction", "VAL_fraction",
                    "ILE_fraction", "ASN_fraction", "ASP_fraction", "LEU_fraction", "MET_fraction",
                    "GLN_fraction", "GLU_fraction", "LYS_fraction", "HIS_fraction", "PHE_fraction",
                    "ARG_fraction", "TYR_fraction", "TRP_fraction", "loop", "helix", "sheet",
                    "liability_score", "liability_num_violations", "bb_rmsd", "bb_rmsd_design",
                    "design_ptm", "design_iptm", "design_to_target_iptm", "min_design_to_target_pae",
                    "affinity_pred_value", "affinity_probability_binary1"
                ],
                description="Per-target metrics from analysis step",
                count=None
            )

        # Filtering step outputs
        if has_filtering:
            # Final designs folder
            final_designs_folder = os.path.join(self.final_ranked_folder, f"final_{self.budget}_designs")
            output["final_designs_folder"] = final_designs_folder

            # Structure pattern: rankNNNN_*.cif
            # We predict the pattern, actual files will be rankNNNN_<design_id>.cif
            structures = [os.path.join(final_designs_folder, "rankNNNN_*.cif")]
            structure_ids = [f"rank{i:04d}" for i in range(1, self.budget + 1)]
            output["structures"] = structures
            output["structure_ids"] = structure_ids

            # All designs metrics
            tables["all_designs_metrics"] = TableInfo(
                name="all_designs_metrics",
                path=self.all_designs_metrics_csv,
                columns=[
                    "id", "final_rank", "designed_sequence", "designed_chain_sequence",
                    "num_design", "affinity_probability_binary1", "design_to_target_iptm",
                    "min_design_to_target_pae", "design_ptm", "filter_rmsd",
                    "designfolding-filter_rmsd", "plip_hbonds_refolded", "delta_sasa_refolded",
                    "design_largest_hydrophobic_patch_refolded", "design_chain_hydrophobicity",
                    "design_hydrophobicity", "loop", "helix", "sheet", "file_name",
                    "num_prot_tokens", "num_lig_atoms", "num_resolved_tokens", "num_tokens",
                    "UNK_fraction", "GLY_fraction", "ALA_fraction", "CYS_fraction",
                    "SER_fraction", "PRO_fraction", "THR_fraction", "VAL_fraction",
                    "ILE_fraction", "ASN_fraction", "ASP_fraction", "LEU_fraction",
                    "MET_fraction", "GLN_fraction", "GLU_fraction", "LYS_fraction",
                    "HIS_fraction", "PHE_fraction", "ARG_fraction", "TYR_fraction",
                    "TRP_fraction", "liability_score", "liability_num_violations",
                    "liability_high_severity_violations", "liability_medium_severity_violations",
                    "liability_low_severity_violations", "liability_AspCleave_count",
                    "liability_ProtTryp_count", "liability_MetOx_count",
                    "liability_HydroPatch_count", "liability_AspBridge_count",
                    "liability_AspBridge_position", "liability_AspBridge_length",
                    "liability_AspBridge_severity", "liability_AspBridge_details",
                    "liability_AspBridge_positions", "liability_AspBridge_num_positions",
                    "liability_AspBridge_global_details", "liability_AspBridge_avg_severity",
                    "liability_DPP4_count", "liability_DPP4_position", "liability_DPP4_length",
                    "liability_DPP4_severity", "liability_DPP4_details", "liability_DPP4_positions",
                    "liability_DPP4_num_positions", "liability_DPP4_global_details",
                    "liability_DPP4_avg_severity", "liability_MetOx_position",
                    "liability_MetOx_length", "liability_MetOx_severity", "liability_MetOx_details",
                    "liability_MetOx_positions", "liability_MetOx_num_positions",
                    "liability_MetOx_global_details", "liability_MetOx_avg_severity",
                    "liability_ProtTryp_position", "liability_ProtTryp_length",
                    "liability_ProtTryp_severity", "liability_ProtTryp_details",
                    "liability_ProtTryp_positions", "liability_ProtTryp_num_positions",
                    "liability_ProtTryp_global_details", "liability_ProtTryp_avg_severity",
                    "liability_HydroPatch_position", "liability_HydroPatch_length",
                    "liability_HydroPatch_severity", "liability_HydroPatch_details",
                    "liability_HydroPatch_positions", "liability_HydroPatch_num_positions",
                    "liability_HydroPatch_global_details", "liability_HydroPatch_avg_severity",
                    "liability_NTCycl_count", "liability_NTCycl_position",
                    "liability_NTCycl_length", "liability_NTCycl_severity",
                    "liability_NTCycl_details", "liability_NTCycl_positions",
                    "liability_NTCycl_num_positions", "liability_NTCycl_global_details",
                    "liability_NTCycl_avg_severity", "liability_AspCleave_position",
                    "liability_AspCleave_length", "liability_AspCleave_severity",
                    "liability_AspCleave_details", "liability_AspCleave_positions",
                    "liability_AspCleave_num_positions", "liability_AspCleave_global_details",
                    "liability_AspCleave_avg_severity", "liability_TrpOx_count",
                    "liability_TrpOx_position", "liability_TrpOx_length",
                    "liability_TrpOx_severity", "liability_TrpOx_details",
                    "liability_TrpOx_positions", "liability_TrpOx_num_positions",
                    "liability_TrpOx_global_details", "liability_TrpOx_avg_severity",
                    "liability_violations_summary", "native_rmsd", "native_rmsd_bb",
                    "native_rmsd_refolded", "native_rmsd_bb_refolded",
                    "designfolding-bb_rmsd", "designfolding-bb_rmsd_design",
                    "designfolding-bb_rmsd_target", "designfolding-bb_rmsd_design_target",
                    "designfolding-bb_target_aligned_rmsd_design", "designfolding-bb_rmsd<2.5",
                    "designfolding-bb_target_aligned<2.5", "designfolding-bb_designability_rmsd_2",
                    "designfolding-bb_designability_rmsd_4", "designfolding-ligand_iptm",
                    "designfolding-interaction_pae", "designfolding-min_interaction_pae",
                    "designfolding-min_design_to_target_pae", "designfolding-iptm",
                    "designfolding-ptm", "designfolding-protein_iptm", "designfolding-design_iptm",
                    "designfolding-design_iiptm", "designfolding-design_to_target_iptm",
                    "designfolding-target_ptm", "designfolding-design_ptm",
                    "designfolding-min_interaction_pae<1.5", "designfolding-min_interaction_pae<2",
                    "designfolding-min_interaction_pae<2.5", "designfolding-min_interaction_pae<3",
                    "designfolding-min_interaction_pae<4", "designfolding-min_interaction_pae<5",
                    "designfolding-design_ptm>80", "designfolding-design_ptm>75",
                    "designfolding-design_iptm>80", "designfolding-design_iptm>70",
                    "designfolding-design_iptm>60", "designfolding-design_iptm>50",
                    "bb_rmsd", "bb_rmsd_design", "bb_rmsd_target", "bb_rmsd_design_target",
                    "bb_target_aligned_rmsd_design", "bb_rmsd<2.5", "bb_target_aligned<2.5",
                    "bb_designability_rmsd_2", "bb_designability_rmsd_4", "ligand_iptm",
                    "interaction_pae", "min_interaction_pae", "iptm", "ptm", "protein_iptm",
                    "design_iptm", "design_iiptm", "target_ptm", "min_interaction_pae<1.5",
                    "min_interaction_pae<2", "min_interaction_pae<2.5", "min_interaction_pae<3",
                    "min_interaction_pae<4", "min_interaction_pae<5", "design_ptm>80",
                    "design_ptm>75", "design_iptm>80", "design_iptm>70", "design_iptm>60",
                    "design_iptm>50", "design_sasa_unbound_refolded", "design_sasa_bound_refolded",
                    "plip_saltbridge_refolded", "affinity_pred_value", "affinity_probability_binary",
                    "affinity_probability_binary2", "affinity_pred_value2", "affinity_pred_value1",
                    "affinity_probability_binary1>50", "affinity_probability_binary1>75",
                    "liability_UnpairedCys_count", "liability_UnpairedCys_position",
                    "liability_UnpairedCys_length", "liability_UnpairedCys_severity",
                    "liability_UnpairedCys_details", "liability_UnpairedCys_positions",
                    "liability_UnpairedCys_num_positions", "liability_UnpairedCys_global_details",
                    "liability_UnpairedCys_avg_severity", "liability_details",
                    "filter_rmsd_design", "neg_min_design_to_target_pae",
                    "neg_design_hydrophobicity", "neg_design_largest_hydrophobic_patch_refolded",
                    "neg_min_interaction_pae", "has_x", "num_filters_passed", "pass_has_x_filter",
                    "pass_filters", "pass_filter_rmsd_filter", "pass_filter_rmsd_design_filter",
                    "pass_designfolding-filter_rmsd_filter", "pass_ALA_fraction_filter",
                    "pass_GLY_fraction_filter", "pass_GLU_fraction_filter",
                    "pass_LEU_fraction_filter", "pass_VAL_fraction_filter",
                    "rank_design_to_target_iptm", "rank_design_ptm",
                    "rank_neg_min_design_to_target_pae", "rank_affinity_probability_binary1",
                    "rank_plip_hbonds_refolded", "rank_plip_saltbridge_refolded",
                    "rank_delta_sasa_refolded", "rank_design_iiptm", "max_rank",
                    "secondary_rank", "quality_score"
                ],
                description="Metrics for all designs considered by filtering",
                count=None
            )

            # Final designs metrics
            tables["final_designs_metrics"] = TableInfo(
                name="final_designs_metrics",
                path=self.final_designs_metrics_csv,
                columns=[
                    "id", "final_rank", "designed_sequence", "designed_chain_sequence",
                    "num_design", "affinity_probability_binary1", "design_to_target_iptm",
                    "min_design_to_target_pae", "design_ptm", "filter_rmsd",
                    "designfolding-filter_rmsd", "plip_hbonds_refolded", "delta_sasa_refolded",
                    "design_largest_hydrophobic_patch_refolded", "design_chain_hydrophobicity",
                    "design_hydrophobicity", "loop", "helix", "sheet", "file_name",
                    "num_prot_tokens", "num_lig_atoms", "num_resolved_tokens", "num_tokens",
                    "UNK_fraction", "GLY_fraction", "ALA_fraction", "CYS_fraction",
                    "SER_fraction", "PRO_fraction", "THR_fraction", "VAL_fraction",
                    "ILE_fraction", "ASN_fraction", "ASP_fraction", "LEU_fraction",
                    "MET_fraction", "GLN_fraction", "GLU_fraction", "LYS_fraction",
                    "HIS_fraction", "PHE_fraction", "ARG_fraction", "TYR_fraction",
                    "TRP_fraction", "liability_score", "liability_num_violations",
                    "liability_high_severity_violations", "liability_medium_severity_violations",
                    "liability_low_severity_violations", "liability_AspCleave_count",
                    "liability_ProtTryp_count", "liability_MetOx_count",
                    "liability_HydroPatch_count", "liability_AspBridge_count",
                    "liability_AspBridge_position", "liability_AspBridge_length",
                    "liability_AspBridge_severity", "liability_AspBridge_details",
                    "liability_AspBridge_positions", "liability_AspBridge_num_positions",
                    "liability_AspBridge_global_details", "liability_AspBridge_avg_severity",
                    "liability_DPP4_count", "liability_DPP4_position", "liability_DPP4_length",
                    "liability_DPP4_severity", "liability_DPP4_details", "liability_DPP4_positions",
                    "liability_DPP4_num_positions", "liability_DPP4_global_details",
                    "liability_DPP4_avg_severity", "liability_MetOx_position",
                    "liability_MetOx_length", "liability_MetOx_severity", "liability_MetOx_details",
                    "liability_MetOx_positions", "liability_MetOx_num_positions",
                    "liability_MetOx_global_details", "liability_MetOx_avg_severity",
                    "liability_ProtTryp_position", "liability_ProtTryp_length",
                    "liability_ProtTryp_severity", "liability_ProtTryp_details",
                    "liability_ProtTryp_positions", "liability_ProtTryp_num_positions",
                    "liability_ProtTryp_global_details", "liability_ProtTryp_avg_severity",
                    "liability_HydroPatch_position", "liability_HydroPatch_length",
                    "liability_HydroPatch_severity", "liability_HydroPatch_details",
                    "liability_HydroPatch_positions", "liability_HydroPatch_num_positions",
                    "liability_HydroPatch_global_details", "liability_HydroPatch_avg_severity",
                    "liability_NTCycl_count", "liability_NTCycl_position",
                    "liability_NTCycl_length", "liability_NTCycl_severity",
                    "liability_NTCycl_details", "liability_NTCycl_positions",
                    "liability_NTCycl_num_positions", "liability_NTCycl_global_details",
                    "liability_NTCycl_avg_severity", "liability_AspCleave_position",
                    "liability_AspCleave_length", "liability_AspCleave_severity",
                    "liability_AspCleave_details", "liability_AspCleave_positions",
                    "liability_AspCleave_num_positions", "liability_AspCleave_global_details",
                    "liability_AspCleave_avg_severity", "liability_TrpOx_count",
                    "liability_TrpOx_position", "liability_TrpOx_length",
                    "liability_TrpOx_severity", "liability_TrpOx_details",
                    "liability_TrpOx_positions", "liability_TrpOx_num_positions",
                    "liability_TrpOx_global_details", "liability_TrpOx_avg_severity",
                    "liability_violations_summary", "native_rmsd", "native_rmsd_bb",
                    "native_rmsd_refolded", "native_rmsd_bb_refolded",
                    "designfolding-bb_rmsd", "designfolding-bb_rmsd_design",
                    "designfolding-bb_rmsd_target", "designfolding-bb_rmsd_design_target",
                    "designfolding-bb_target_aligned_rmsd_design", "designfolding-bb_rmsd<2.5",
                    "designfolding-bb_target_aligned<2.5", "designfolding-bb_designability_rmsd_2",
                    "designfolding-bb_designability_rmsd_4", "designfolding-ligand_iptm",
                    "designfolding-interaction_pae", "designfolding-min_interaction_pae",
                    "designfolding-min_design_to_target_pae", "designfolding-iptm",
                    "designfolding-ptm", "designfolding-protein_iptm", "designfolding-design_iptm",
                    "designfolding-design_iiptm", "designfolding-design_to_target_iptm",
                    "designfolding-target_ptm", "designfolding-design_ptm",
                    "designfolding-min_interaction_pae<1.5", "designfolding-min_interaction_pae<2",
                    "designfolding-min_interaction_pae<2.5", "designfolding-min_interaction_pae<3",
                    "designfolding-min_interaction_pae<4", "designfolding-min_interaction_pae<5",
                    "designfolding-design_ptm>80", "designfolding-design_ptm>75",
                    "designfolding-design_iptm>80", "designfolding-design_iptm>70",
                    "designfolding-design_iptm>60", "designfolding-design_iptm>50",
                    "bb_rmsd", "bb_rmsd_design", "bb_rmsd_target", "bb_rmsd_design_target",
                    "bb_target_aligned_rmsd_design", "bb_rmsd<2.5", "bb_target_aligned<2.5",
                    "bb_designability_rmsd_2", "bb_designability_rmsd_4", "ligand_iptm",
                    "interaction_pae", "min_interaction_pae", "iptm", "ptm", "protein_iptm",
                    "design_iptm", "design_iiptm", "target_ptm", "min_interaction_pae<1.5",
                    "min_interaction_pae<2", "min_interaction_pae<2.5", "min_interaction_pae<3",
                    "min_interaction_pae<4", "min_interaction_pae<5", "design_ptm>80",
                    "design_ptm>75", "design_iptm>80", "design_iptm>70", "design_iptm>60",
                    "design_iptm>50", "design_sasa_unbound_refolded", "design_sasa_bound_refolded",
                    "plip_saltbridge_refolded", "affinity_pred_value", "affinity_probability_binary",
                    "affinity_probability_binary2", "affinity_pred_value2", "affinity_pred_value1",
                    "affinity_probability_binary1>50", "affinity_probability_binary1>75",
                    "liability_UnpairedCys_count", "liability_UnpairedCys_position",
                    "liability_UnpairedCys_length", "liability_UnpairedCys_severity",
                    "liability_UnpairedCys_details", "liability_UnpairedCys_positions",
                    "liability_UnpairedCys_num_positions", "liability_UnpairedCys_global_details",
                    "liability_UnpairedCys_avg_severity", "liability_details",
                    "filter_rmsd_design", "neg_min_design_to_target_pae",
                    "neg_design_hydrophobicity", "neg_design_largest_hydrophobic_patch_refolded",
                    "neg_min_interaction_pae", "has_x", "num_filters_passed", "pass_has_x_filter",
                    "pass_filters", "pass_filter_rmsd_filter", "pass_filter_rmsd_design_filter",
                    "pass_designfolding-filter_rmsd_filter", "pass_ALA_fraction_filter",
                    "pass_GLY_fraction_filter", "pass_GLU_fraction_filter",
                    "pass_LEU_fraction_filter", "pass_VAL_fraction_filter",
                    "rank_design_to_target_iptm", "rank_design_ptm",
                    "rank_neg_min_design_to_target_pae", "rank_affinity_probability_binary1",
                    "rank_plip_hbonds_refolded", "rank_plip_saltbridge_refolded",
                    "rank_delta_sasa_refolded", "rank_design_iiptm", "max_rank",
                    "secondary_rank", "quality_score"
                ],
                description=f"Metrics for top {self.budget} designs after filtering",
                count=self.budget
            )

            # Results overview PDF
            output["results_pdf"] = self.results_overview_pdf
            output["main"] = self.final_designs_metrics_csv

        # Folding step outputs (only if filtering not present)
        elif has_folding:
            # Structures in refold_cif folder
            refold_cif_folder = os.path.join(self.intermediate_inverse_folded_folder, "refold_cif")
            output["refold_cif_folder"] = refold_cif_folder

            # Structure pattern: design_spec_NNN.cif
            structures = [os.path.join(refold_cif_folder, "design_spec_*.cif")]
            structure_ids = [f"design_spec_{i}" for i in range(self.num_designs)]
            output["structures"] = structures
            output["structure_ids"] = structure_ids

        output["tables"] = tables
        return output

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration with all BoltzGen parameters."""
        base_dict = super().to_dict()

        # Prepare design spec for serialization
        if self.spec_mode == "automatic":
            design_spec_value = f"<automatic from {self.target_structure_path}>"
        elif self.design_spec_is_dict:
            design_spec_value = "<dict>"
        else:
            design_spec_value = self.design_spec

        base_dict.update({
            "tool_params": {
                "spec_mode": self.spec_mode,
                "design_spec": design_spec_value,
                "target_structure": str(self.target_structure) if self.target_structure else None,
                "ligand": self.ligand,
                "binder_spec": self.binder_spec,
                "binding_region": self.binding_region,
                "protocol": self.protocol,
                "num_designs": self.num_designs,
                "budget": self.budget,
                "design_checkpoints": self.design_checkpoints,
                "step_scale": self.step_scale,
                "noise_scale": self.noise_scale,
                "diffusion_batch_size": self.diffusion_batch_size,
                "inverse_fold_num_sequences": self.inverse_fold_num_sequences,
                "skip_inverse_folding": self.skip_inverse_folding,
                "alpha": self.alpha,
                "filter_biased": self.filter_biased,
                "additional_filters": self.additional_filters,
                "metrics_override": self.metrics_override,
                "devices": self.devices,
                "reuse": self.reuse,
                "steps": self.steps,
                "cache_dir": self.cache_dir
            }
        })
        return base_dict


class BoltzGenMerge(BaseConfig):
    """
    Merge multiple BoltzGen output directories into one.

    Use this to combine results from parallel BoltzGen runs before
    re-running filtering on the combined set.

    Example:
        run1 = BoltzGen(ligand=..., binder_spec="140-180", ...)
        run2 = BoltzGen(ligand=..., binder_spec="140-180", ...)
        merged = BoltzGenMerge(sources=[run1, run2])
        filtered = BoltzGen(reuse=merged, steps=["filtering"], budget=100, ...)

    With ID renaming to avoid collisions:
        merged = BoltzGenMerge(sources=[run1, run2, ...],
                              id_template="batch{i:03d}_")
        # This renames design_spec_001 -> batch000_design_spec_001, etc.
    """

    TOOL_NAME = "BoltzGenMerge"
    DEFAULT_ENV = None  # Uses same env as BoltzGen

    def __init__(self,
                 sources: List[Union[ToolOutput, StandardizedOutput, str]] = None,
                 id_template: str = None,
                 **kwargs):
        """
        Initialize BoltzGenMerge configuration.

        Args:
            sources: List of BoltzGen outputs to merge. Can be ToolOutput,
                    StandardizedOutput, or direct paths to output directories.
            id_template: Template for renaming design IDs to avoid collisions.
                        If provided, uses custom merge logic instead of
                        `boltzgen merge` command.
                        Supports:
                        - {i} for source index (0-based)
                        - {name} for source folder name
                        Examples:
                        - "batch{i:03d}_" -> batch000_design_spec_001
                        - "{name}_" -> BoltzGen_Run1_design_spec_001
                        Default: None (uses standard boltzgen merge)
            **kwargs: Additional parameters passed to BaseConfig
        """
        self.sources = sources or []
        self.source_paths = []
        self.id_template = id_template

        # Initialize base class
        super().__init__(**kwargs)

    def validate_params(self):
        """Validate BoltzGenMerge parameters."""
        if not self.sources or len(self.sources) < 2:
            raise ValueError("BoltzGenMerge requires at least 2 source directories to merge")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input sources."""
        self.folders = pipeline_folders

        # Store HelpScripts path for custom merge
        self.merge_helper_py = os.path.join(
            pipeline_folders.get("HelpScripts", "HelpScripts"),
            "pipe_boltzgen_merge.py"
        )

        # Reset source_paths to avoid accumulation if configure_inputs is called multiple times
        self.source_paths = []

        # Extract paths from sources
        for source in self.sources:
            if isinstance(source, StandardizedOutput):
                if hasattr(source, 'output_folder') and source.output_folder:
                    self.source_paths.append(source.output_folder)
                else:
                    raise ValueError("StandardizedOutput has no output_folder")
            elif isinstance(source, ToolOutput):
                # Get output folder from ToolOutput
                output_folder = source.get_output_files("output_folder")
                if output_folder:
                    self.source_paths.append(output_folder[0] if isinstance(output_folder, list) else output_folder)
                    self.dependencies.append(source.config)
                else:
                    raise ValueError("ToolOutput has no output_folder")
            elif isinstance(source, str):
                self.source_paths.append(source)
            else:
                raise ValueError(f"Invalid source type: {type(source)}")

    def generate_script(self, script_path: str) -> str:
        """Generate BoltzGenMerge execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        # Build source paths string
        sources_str = ' '.join(f'"{p}"' for p in self.source_paths)

        # Use first source for design_spec.yaml copy
        first_source = self.source_paths[0]
        design_spec_dest = os.path.join(self.output_folder, "design_spec.yaml")

        if self.id_template:
            # Use custom merge with ID renaming
            script_content = f"""#!/bin/bash
# BoltzGenMerge execution script (with ID renaming)
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Merging BoltzGen outputs with ID renaming"
echo "Sources: {len(self.source_paths)} directories"
echo "Output: {self.output_folder}"
echo "ID template: {self.id_template}"

# Run custom merge with ID renaming
python "{self.merge_helper_py}" \\
    --sources {sources_str} \\
    --output "{self.output_folder}" \\
    --id_template "{self.id_template}"

if [ $? -eq 0 ]; then
    echo "BoltzGenMerge with ID renaming completed successfully"
else
    echo "Error: BoltzGenMerge failed"
    exit 1
fi

# Copy design_spec.yaml from first source for reuse compatibility
if [ -f "{first_source}/design_spec.yaml" ]; then
    cp "{first_source}/design_spec.yaml" "{design_spec_dest}"
    echo "Copied design_spec.yaml from first source"
else
    echo "Warning: No design_spec.yaml found in first source directory"
fi

{self.generate_completion_check_footer()}
"""
        else:
            # Use standard boltzgen merge command
            script_content = f"""#!/bin/bash
# BoltzGenMerge execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Merging BoltzGen outputs"
echo "Sources: {', '.join(os.path.basename(p) for p in self.source_paths)}"
echo "Output: {self.output_folder}"

# Run BoltzGen merge
boltzgen merge {sources_str} --output "{self.output_folder}"

if [ $? -eq 0 ]; then
    echo "BoltzGenMerge completed successfully"
else
    echo "Error: BoltzGenMerge failed"
    exit 1
fi

# Copy design_spec.yaml from first source for reuse compatibility
if [ -f "{first_source}/design_spec.yaml" ]; then
    cp "{first_source}/design_spec.yaml" "{design_spec_dest}"
    echo "Copied design_spec.yaml from first source"
else
    echo "Warning: No design_spec.yaml found in first source directory"
fi

{self.generate_completion_check_footer()}
"""
        return script_content

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after merge."""
        return {
            "output_folder": self.output_folder,
            "tool_folder": self.output_folder,
            # Merged outputs will have same structure as BoltzGen
            "structures": [],  # Will be populated after merge
            "structure_ids": [],
            "tables": {}
        }

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()
        config_lines.append(f"SOURCES: {len(self.sources)} directories")
        for i, path in enumerate(self.source_paths):
            config_lines.append(f"  [{i+1}] {os.path.basename(path)}")
        if self.id_template:
            config_lines.append(f"ID TEMPLATE: {self.id_template}")
            config_lines.append("  (Using custom merge with ID renaming)")
        else:
            config_lines.append("ID TEMPLATE: None (using standard boltzgen merge)")
        return config_lines


class BoltzGenImport(BaseConfig):
    """
    Import external structures into BoltzGen filesystem format.

    This tool converts structures (e.g., from RFdiffusion) into the BoltzGen
    directory structure, enabling use of downstream BoltzGen steps like
    folding, analysis, and filtering.

    Two modes are supported:
    1. Design only: Import backbone structures (simulates "design" step output)
       BoltzGenImport(designs=rfd_output, ligand=lig, binder_spec="140-180")

    2. Design + Inverse Folding: Import structures with sequences
       (simulates "inverse_folding" step output)
       BoltzGenImport(designs=rfd_output, sequences=lmpnn_output, ligand=lig, binder_spec="140-180")

    After import, use BoltzGen(reuse=imported, steps=["folding", ...]) to
    continue the pipeline.

    Notes:
    - Chain reassignment: Input structures (e.g., from RFdiffusion) may have
      protein+ligand in chain A. This tool separates them into chain A (protein)
      and chain B (ligand) as BoltzGen expects.
    - When sequences are provided, sidechain atom coordinates are set to (0,0,0)
      as BoltzGen expects - the folding step will predict them.
    - NPZ files are not generated (BoltzGen should handle their absence).
    """

    TOOL_NAME = "BoltzGenImport"
    DEFAULT_ENV = None

    def __init__(self,
                 designs: Union[ToolOutput, StandardizedOutput] = None,
                 sequences: Union[ToolOutput, StandardizedOutput] = None,
                 # Design spec parameters (same as BoltzGen)
                 ligand: Union[ToolOutput, StandardizedOutput] = None,
                 binder_spec: Union[str, Dict[str, Any]] = None,
                 protocol: str = "protein-small_molecule",
                 ligand_chain: str = "B",
                 protein_chain: str = "A",
                 id_map: Dict[str, str] = {"*": "*_<N>"},
                 ligand_name: str = "LIG1",
                 num_designs: int = None,
                 **kwargs):
        """
        Initialize BoltzGenImport configuration.

        Args:
            designs: Structure input from a design tool (RFdiffusion, etc.).
                    Must have 'structures' output with PDB/CIF files.
            sequences: Optional sequence input from inverse folding tool
                      (ProteinMPNN, LigandMPNN, etc.). Must have a table
                      with 'id' and 'sequence' columns that map to designs.
                      If provided, simulates inverse_folding step output.
            ligand: Ligand input from Ligand tool (for design_spec generation).
                   Used to get SMILES/CCD for the design specification.
            binder_spec: Binder length specification (e.g., "140-180" or {"length": [140, 180]}).
                        Used for design_spec generation.
            protocol: BoltzGen protocol (default: "protein-small_molecule").
                     Used when running minimal design step to generate molecule pickle.
            ligand_chain: Chain ID for ligand in output (default: "B")
            protein_chain: Chain ID for protein in output (default: "A")
            id_map: ID mapping pattern for matching structure IDs to sequence IDs.
                   Default {"*": "*_<N>"} handles recursive numeric suffixes:
                   - Structure "design_1" matches sequence "design_1_1", "design_1_1_1", etc.
                   - Set to {"*": "*"} for exact ID matching.
            ligand_name: Residue name for ligand in output (default: "LIG1").
                        BoltzGen expects ligands to match pattern ^LIG\\d+ (e.g., LIG1, LIG2).
            num_designs: Number of designs being imported (auto-detected if None).
                        Used for consistent naming with BoltzGen conventions.
            **kwargs: Additional parameters passed to BaseConfig
        """
        self.designs = designs
        self.sequences = sequences
        self.ligand = ligand
        self.binder_spec = binder_spec
        self.protocol = protocol
        self.ligand_chain = ligand_chain
        self.protein_chain = protein_chain
        self.id_map = id_map
        self.ligand_name = ligand_name
        self.num_designs = num_designs

        # Paths to be resolved in configure_inputs
        self.design_structures = []  # List of PDB/CIF paths
        self.sequences_csv = None    # Path to sequences CSV
        self.ligand_compounds_csv = None  # Path to ligand compounds CSV

        # Determine mode
        self.mode = "inverse_folding" if sequences is not None else "design"

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate BoltzGenImport parameters."""
        if self.designs is None:
            raise ValueError("designs parameter is required")

        if not isinstance(self.designs, (ToolOutput, StandardizedOutput)):
            raise ValueError(
                f"designs must be ToolOutput or StandardizedOutput, "
                f"got {type(self.designs)}"
            )

        if self.sequences is not None:
            if not isinstance(self.sequences, (ToolOutput, StandardizedOutput)):
                raise ValueError(
                    f"sequences must be ToolOutput or StandardizedOutput, "
                    f"got {type(self.sequences)}"
                )

        if self.binder_spec is None:
            raise ValueError("binder_spec is required for design_spec generation")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input sources."""
        self.folders = pipeline_folders

        # Helper script path
        self.import_helper_py = os.path.join(
            pipeline_folders.get("HelpScripts", "HelpScripts"),
            "pipe_boltzgen_import.py"
        )

        # Extract structure paths from designs
        if isinstance(self.designs, StandardizedOutput):
            if hasattr(self.designs, 'structures') and self.designs.structures:
                self.design_structures = self.designs.structures
            else:
                raise ValueError("designs StandardizedOutput has no structures")
        elif isinstance(self.designs, ToolOutput):
            structures = self.designs.get_output_files("structures")
            if structures:
                self.design_structures = structures
                self.dependencies.append(self.designs.config)
            else:
                raise ValueError("designs ToolOutput has no structures")

        # Extract sequences CSV if provided
        if self.sequences is not None:
            if isinstance(self.sequences, StandardizedOutput):
                # Look for sequences table or CSV
                if hasattr(self.sequences, 'tables'):
                    if hasattr(self.sequences.tables, 'sequences'):
                        self.sequences_csv = self.sequences.tables.sequences.path
                    elif hasattr(self.sequences.tables, '_tables'):
                        for name, info in self.sequences.tables._tables.items():
                            if 'sequence' in name.lower() or name == 'main':
                                self.sequences_csv = info.path
                                break
                if not self.sequences_csv and hasattr(self.sequences, 'sequences'):
                    seqs = self.sequences.sequences
                    if seqs and len(seqs) > 0:
                        self.sequences_csv = seqs[0]
            elif isinstance(self.sequences, ToolOutput):
                # Try to get sequences table
                output_files = self.sequences.get_output_files()
                if 'tables' in output_files:
                    tables = output_files['tables']
                    if isinstance(tables, dict):
                        for name, info in tables.items():
                            if 'sequence' in name.lower() or name == 'main':
                                if hasattr(info, 'path'):
                                    self.sequences_csv = info.path
                                elif isinstance(info, dict) and 'path' in info:
                                    self.sequences_csv = info['path']
                                break
                self.dependencies.append(self.sequences.config)

            if not self.sequences_csv:
                raise ValueError(
                    "Could not find sequences CSV from sequences input. "
                    "Ensure it has a table with 'id' and 'sequence' columns."
                )

        # Extract ligand compounds CSV if provided
        if self.ligand is not None:
            if isinstance(self.ligand, StandardizedOutput):
                if hasattr(self.ligand, 'compounds') and self.ligand.compounds:
                    self.ligand_compounds_csv = self.ligand.compounds[0]
                elif hasattr(self.ligand, 'tables') and hasattr(self.ligand.tables, 'compounds'):
                    self.ligand_compounds_csv = self.ligand.tables.compounds.path
            elif isinstance(self.ligand, ToolOutput):
                compounds = self.ligand.get_output_files("compounds")
                if compounds:
                    self.ligand_compounds_csv = compounds[0]
                    self.dependencies.append(self.ligand.config)

    def _parse_binder_spec(self):
        """Parse binder_spec into min and max length values."""
        if isinstance(self.binder_spec, dict):
            length = self.binder_spec.get("length", [100, 200])
            if isinstance(length, list) and len(length) == 2:
                return length[0], length[1]
            elif isinstance(length, int):
                return length, length
            else:
                return 100, 200
        elif isinstance(self.binder_spec, str):
            if "-" in self.binder_spec:
                parts = self.binder_spec.split("-")
                return int(parts[0]), int(parts[1])
            else:
                length = int(self.binder_spec)
                return length, length
        else:
            return 100, 200

    def generate_script(self, script_path: str) -> str:
        """Generate BoltzGenImport execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        # Create subdirectories
        intermediate_designs = os.path.join(self.output_folder, "intermediate_designs")
        intermediate_inverse_folded = os.path.join(
            self.output_folder, "intermediate_designs_inverse_folded"
        )
        config_dir = os.path.join(self.output_folder, "config")

        os.makedirs(intermediate_designs, exist_ok=True)
        os.makedirs(config_dir, exist_ok=True)
        if self.mode == "inverse_folding":
            os.makedirs(intermediate_inverse_folded, exist_ok=True)

        # Parse binder spec
        binder_min, binder_max = self._parse_binder_spec()

        # Build structures list file
        structures_list_file = os.path.join(self.output_folder, ".input_structures.txt")
        with open(structures_list_file, 'w') as f:
            for struct in self.design_structures:
                f.write(f"{struct}\n")

        # Write steps.yaml
        steps_yaml_file = os.path.join(self.output_folder, "steps.yaml")
        if self.mode == "design":
            steps_content = """steps:
- name: design
  config_file: config/design.yaml
"""
        else:  # inverse_folding
            steps_content = """steps:
- name: design
  config_file: config/design.yaml
- name: inverse_folding
  config_file: config/inverse_folding.yaml
"""
        with open(steps_yaml_file, 'w') as f:
            f.write(steps_content)

        # Serialize id_map as JSON string for passing to script
        import json
        id_map_json = json.dumps(self.id_map)

        # Build command arguments
        cmd_args = [
            f'--structures "{structures_list_file}"',
            f'--output "{self.output_folder}"',
            f'--mode {self.mode}',
            f'--protein-chain {self.protein_chain}',
            f'--ligand-chain {self.ligand_chain}',
            f'--binder-min {binder_min}',
            f'--binder-max {binder_max}',
            f'--ligand-name {self.ligand_name}',
            f"--id-map '{id_map_json}'",
        ]

        if self.sequences_csv:
            cmd_args.append(f'--sequences "{self.sequences_csv}"')

        if self.ligand_compounds_csv:
            cmd_args.append(f'--ligand-csv "{self.ligand_compounds_csv}"')

        cmd_str = ' \\\n    '.join(cmd_args)

        # Determine target directory for molecules
        if self.mode == "inverse_folding":
            molecules_target_dir = os.path.join(intermediate_inverse_folded, "molecules_out_dir")
        else:
            molecules_target_dir = os.path.join(intermediate_designs, "molecules_out_dir")

        # Build molecule generation command if ligand is provided
        molecule_gen_section = ""
        if self.ligand_compounds_csv:
            temp_design_dir = os.path.join(self.output_folder, ".temp_molecule_gen")
            design_spec_yaml = os.path.join(self.output_folder, "design_spec.yaml")
            molecule_gen_section = f'''
# Step 1: Generate design_spec.yaml for molecule generation
echo "Generating design_spec.yaml for molecule pickle generation..."
python "{self.import_helper_py}" \\
    --generate-design-spec \\
    --ligand-csv "{self.ligand_compounds_csv}" \\
    --output "{design_spec_yaml}" \\
    --protein-chain {self.protein_chain} \\
    --ligand-chain {self.ligand_chain} \\
    --binder-min {binder_min} \\
    --binder-max {binder_max}

if [ $? -ne 0 ]; then
    echo "Error: Failed to generate design_spec.yaml"
    exit 1
fi

# Step 2: Generate molecule pickle by running minimal BoltzGen design
echo "Generating molecule pickle via BoltzGen design step..."
TEMP_MOL_DIR="{temp_design_dir}"
mkdir -p "$TEMP_MOL_DIR"

# Run BoltzGen with 1 design to generate the molecule pickle
boltzgen run "{design_spec_yaml}" \\
    --output "$TEMP_MOL_DIR" \\
    --protocol {self.protocol} \\
    --num_designs 1 \\
    --steps design

if [ $? -ne 0 ]; then
    echo "Error: Failed to generate molecule pickle"
    exit 1
fi

# Copy molecule pickle to target directory
MOLECULES_SOURCE="$TEMP_MOL_DIR/intermediate_designs/molecules_out_dir"
MOLECULES_TARGET="{molecules_target_dir}"
mkdir -p "$MOLECULES_TARGET"

if [ -d "$MOLECULES_SOURCE" ]; then
    echo "Copying molecule pickles from $MOLECULES_SOURCE to $MOLECULES_TARGET"
    cp -r "$MOLECULES_SOURCE"/* "$MOLECULES_TARGET"/
else
    echo "Warning: No molecules_out_dir found in BoltzGen output"
fi

# Keep temp directory for debugging (contains molecule pickle generation output)
echo "Molecule pickle generation complete (temp files preserved in $TEMP_MOL_DIR)"
'''

        script_content = f"""#!/bin/bash
# BoltzGenImport execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Importing structures into BoltzGen format"
echo "Mode: {self.mode}"
echo "Input structures: {len(self.design_structures)} files"
echo "Output: {self.output_folder}"
echo "Binder spec: {binder_min}-{binder_max}"
{"echo 'Sequences CSV: " + self.sequences_csv + "'" if self.sequences_csv else ""}
{"echo 'Ligand CSV: " + self.ligand_compounds_csv + "'" if self.ligand_compounds_csv else ""}
{molecule_gen_section}
# Step 3: Run import helper to convert structures
echo "Converting structures to BoltzGen format..."
python "{self.import_helper_py}" \\
    {cmd_str}

if [ $? -eq 0 ]; then
    echo "BoltzGenImport completed successfully"
else
    echo "Error: BoltzGenImport failed"
    exit 1
fi

{self.generate_completion_check_footer()}
"""
        return script_content

    def get_output_files(self) -> Dict[str, Any]:
        """
        Get expected output files.

        BoltzGenImport does not predict outputs - it creates the filesystem
        structure for subsequent BoltzGen steps.
        """
        return {
            "output_folder": self.output_folder,
            "tool_folder": self.output_folder,
            "structures": [],
            "structure_ids": [],
            "tables": {}
        }

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()
        config_lines.append(f"MODE: {self.mode}")
        config_lines.append(f"PROTOCOL: {self.protocol}")
        config_lines.append(f"INPUT STRUCTURES: {len(self.design_structures)} files")
        config_lines.append(f"BINDER SPEC: {self.binder_spec}")
        config_lines.append(f"CHAINS: protein={self.protein_chain}, ligand={self.ligand_chain}")
        config_lines.append(f"LIGAND NAME: {self.ligand_name}")
        config_lines.append(f"ID MAP: {self.id_map}")
        if self.sequences_csv:
            config_lines.append(f"SEQUENCES CSV: {os.path.basename(self.sequences_csv)}")
        if self.ligand_compounds_csv:
            config_lines.append(f"LIGAND CSV: {os.path.basename(self.ligand_compounds_csv)}")
        return config_lines


# =============================================================================
# BOLTZGEN INTERNAL FORMAT DOCUMENTATION
# =============================================================================
# This section documents the internal file formats used by BoltzGen, extracted
# from the boltzgen package source code (v0.1.4). This information is needed
# when importing external structures (e.g., from RFdiffusion) into BoltzGen
# filesystem format.
#
# Source: boltzgen/task/predict/data_from_generated.py
#         boltzgen/data/const.py
#
# -----------------------------------------------------------------------------
# NPZ METADATA FILES
# -----------------------------------------------------------------------------
# BoltzGen expects each CIF structure file to have a corresponding NPZ metadata
# file with the same base name (e.g., design_spec_0.cif -> design_spec_0.npz).
#
# The NPZ file contains per-token arrays where N_tokens = N_protein_residues + N_ligand_atoms
#
# REQUIRED FIELDS:
# ----------------
# design_mask : np.ndarray, dtype=float32, shape=(N_tokens,)
#     Indicates which tokens are designable (part of the binder).
#     Values: 1.0 = designable (binder protein residues)
#             0.0 = fixed (target/ligand atoms)
#
# OPTIONAL FIELDS (may cause issues if missing):
# ----------------------------------------------
# inverse_fold_design_mask : np.ndarray or object
#     Same as design_mask, or can be overridden for inverse folding.
#     Used when use_new_design_mask=True in BoltzGen.
#     Can be object type for complex specifications.
#
# ss_type : np.ndarray, dtype=int64, shape=(N_tokens,)
#     Secondary structure type per token.
#     Values (from boltzgen/data/const.py):
#         0 = UNSPECIFIED
#         1 = LOOP
#         2 = HELIX
#         3 = SHEET
#
# binding_type : np.ndarray, dtype=int64 or float32, shape=(N_tokens,)
#     Binding type per token.
#     Values (from boltzgen/data/const.py):
#         0 = UNSPECIFIED
#         1 = BINDING
#         2 = NOT_BINDING
#
# token_resolved_mask : np.ndarray, dtype=float32, shape=(N_tokens,)
#     Indicates which tokens have resolved coordinates.
#     Values: 1.0 = resolved (has coordinates)
#             0.0 = unresolved
#
# NOTE: mol_type is NOT loaded from NPZ - it's computed by the tokenizer from
# the structure. Chain types are:
#     0 = PROTEIN
#     1 = DNA
#     2 = RNA
#     3 = NONPOLYMER (ligands)
#
# -----------------------------------------------------------------------------
# EXAMPLE NPZ GENERATION (for BoltzGenImport)
# -----------------------------------------------------------------------------
# For a structure with N protein residues and M ligand atoms:
#
#     import numpy as np
#
#     n_tokens = n_protein_residues + n_ligand_atoms
#
#     # Design mask: 1.0 for protein (designed), 0.0 for ligand (fixed)
#     design_mask = np.concatenate([
#         np.ones(n_protein_residues, dtype=np.float32),
#         np.zeros(n_ligand_atoms, dtype=np.float32)
#     ])
#
#     # Secondary structure: all unspecified
#     ss_type = np.zeros(n_tokens, dtype=np.int64)
#
#     # Binding type: all unspecified
#     binding_type = np.zeros(n_tokens, dtype=np.int64)
#
#     # Token resolved mask: all resolved
#     token_resolved_mask = np.ones(n_tokens, dtype=np.float32)
#
#     np.savez(
#         output_npz_path,
#         design_mask=design_mask,
#         inverse_fold_design_mask=design_mask,  # Same as design_mask
#         ss_type=ss_type,
#         binding_type=binding_type,
#         token_resolved_mask=token_resolved_mask
#     )
#
# -----------------------------------------------------------------------------
# FILESYSTEM STRUCTURE
# -----------------------------------------------------------------------------
# BoltzGen expects the following directory structure:
#
# <output_folder>/
# ├── design_spec.yaml              # Design specification (entities, bindings)
# ├── steps.yaml                    # Completed pipeline steps
# ├── config/                       # Step configuration files
# │   ├── design.yaml
# │   └── inverse_folding.yaml
# ├── intermediate_designs/         # Output from "design" step
# │   ├── design_spec_0.cif
# │   ├── design_spec_0.npz         # <-- Required metadata
# │   ├── design_spec_1.cif
# │   ├── design_spec_1.npz
# │   └── ...
# └── intermediate_designs_inverse_folded/  # Output from "inverse_folding" step
#     ├── design_spec_0.cif
#     ├── design_spec_0.npz         # <-- Required metadata
#     ├── design_spec_1.cif
#     ├── design_spec_1.npz
#     └── ...
#
# =============================================================================

