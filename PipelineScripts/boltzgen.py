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
                 target_structure: Union[str, ToolOutput, StandardizedOutput] = None,
                 ligand: str = None,
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
                 alpha: float = 0.5,
                 filter_biased: bool = True,
                 additional_filters: Optional[List[str]] = None,
                 metrics_override: Optional[Dict[str, float]] = None,
                 # Execution parameters
                 devices: Optional[int] = None,
                 reuse: bool = False,
                 steps: Optional[List[str]] = None,
                 cache_dir: Optional[str] = None,
                 **kwargs):
        """
        Initialize BoltzGen configuration.

        You can provide EITHER design_spec (manual) OR target_structure + binder_spec (automatic).

        Args:
            design_spec: YAML configuration string/dict or path to YAML file defining:
                        - Target entities (proteins, ligands)
                        - Binder specification (sequence ranges)
                        - Binding site constraints
                        - Secondary structure specifications
            target_structure: BioPipelines structure input (ToolOutput or file path) for automatic spec building
            ligand: Ligand residue name in target structure (e.g., "LIG", "ATP") for automatic binding site detection
            binder_spec: Binder specification for automatic building:
                        - String format: "50-100" for length range, or "50" for fixed length
                        - Dict format: {"length": [50, 100]} or custom specification
            binding_region: Target binding region (optional, PyMOL-style: "A:9-140")
                          If not provided and ligand is given, will auto-detect from ligand proximity
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
            alpha: Quality/diversity trade-off (0.0=quality, 1.0=diversity, default: 0.5)
            filter_biased: Remove amino acid composition outliers (default: True)
            additional_filters: Hard threshold expressions (e.g., ["design_ALA>0.3"])
            metrics_override: Per-metric ranking weights override
            devices: Number of GPUs to use (default: auto-detect)
            reuse: Resume interrupted runs (default: False)
            steps: Run only specific pipeline steps (default: all steps)
            cache_dir: Model/cache location (default: /home/$USER/data/boltzgen, auto-downloaded on first run)
            **kwargs: Additional parameters passed to BaseConfig
        """
        # Store design specification parameters
        self.design_spec = design_spec
        self.target_structure = target_structure
        self.ligand = ligand
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
        self.additional_filters = additional_filters or []
        self.metrics_override = metrics_override

        # Execution parameters
        self.devices = devices
        self.reuse = reuse
        self.steps = steps or []
        self.cache_dir = cache_dir

        # Determine specification mode
        if design_spec is not None:
            # Manual mode: user provided design_spec
            self.spec_mode = "manual"
            self.design_spec_is_dict = isinstance(design_spec, dict)
            self.design_spec_is_yaml_str = isinstance(design_spec, str) and (
                'entities:' in design_spec or '- protein:' in design_spec
            )
            self.design_spec_is_file = isinstance(design_spec, str) and not self.design_spec_is_yaml_str
        elif target_structure is not None and binder_spec is not None:
            # Automatic mode: build from BioPipelines inputs
            self.spec_mode = "automatic"
            self.design_spec_is_dict = False
            self.design_spec_is_yaml_str = False
            self.design_spec_is_file = False
        else:
            self.spec_mode = None  # Will fail validation

        # Track inputs from previous tools
        self.input_structures = None
        self.input_compounds = None
        self.target_structure_path = None

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
                "(target_structure + binder_spec) for automatic building"
            )

        if self.spec_mode == "manual" and not self.design_spec:
            raise ValueError("design_spec cannot be empty in manual mode")

        if self.spec_mode == "automatic":
            if not self.target_structure:
                raise ValueError("target_structure is required in automatic mode")
            if not self.binder_spec:
                raise ValueError("binder_spec is required in automatic mode")

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
            raise ValueError(
                f"budget ({self.budget}) cannot exceed num_designs ({self.num_designs})"
            )

        if self.inverse_fold_num_sequences <= 0:
            raise ValueError("inverse_fold_num_sequences must be positive")

        if not (0.0 <= self.alpha <= 1.0):
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

        elif self.spec_mode == "automatic":
            # Handle automatic mode: extract structure path from ToolOutput
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
        else:  # automatic
            spec_display = "Automatic (from target_structure + binder_spec)"
            if self.target_structure_path:
                spec_display += f"\n  Target: {os.path.basename(self.target_structure_path)}"
            if self.ligand:
                spec_display += f"\n  Ligand: {self.ligand}"
            if self.binder_spec:
                spec_display += f"\n  Binder: {self.binder_spec}"
            if self.binding_region:
                spec_display += f"\n  Binding region: {self.binding_region}"

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

        # Prepare design specification file
        if self.spec_mode == "automatic":
            # Build YAML spec automatically
            import yaml
            auto_spec = self._build_automatic_yaml_spec()
            with open(self.design_spec_yaml_file, 'w') as f:
                yaml.dump(auto_spec, f, default_flow_style=False, sort_keys=False)
        elif self.design_spec_is_dict:
            # Convert dictionary to YAML and write to file
            import yaml
            with open(self.design_spec_yaml_file, 'w') as f:
                yaml.dump(self.design_spec, f, default_flow_style=False)
        elif self.design_spec_is_yaml_str:
            # Write YAML string to file
            with open(self.design_spec_yaml_file, 'w') as f:
                f.write(self.design_spec)
        elif self.design_spec_is_file:
            # Copy existing file to output folder
            import shutil
            shutil.copy(self.design_spec, self.design_spec_yaml_file)

        # Build BoltzGen command
        boltzgen_cmd = self._build_boltzgen_command()

        # Generate script content
        script_content = f"""#!/bin/bash
# BoltzGen execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Running BoltzGen binder design"
echo "Design specification: {self.design_spec_yaml_file}"
echo "Protocol: {self.protocol}"
echo "Output folder: {self.output_folder}"

# Run BoltzGen
{boltzgen_cmd}

if [ $? -eq 0 ]; then
    echo "BoltzGen completed successfully"
else
    echo "Error: BoltzGen failed"
    exit 1
fi

# Run postprocessing to extract sequences and format tables
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

{self.generate_completion_check_footer()}
"""

        return script_content

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
        cmd_parts.append(f'--alpha {self.alpha}')

        if not self.filter_biased:
            cmd_parts.append('--filter_biased false')

        if self.additional_filters:
            for filter_expr in self.additional_filters:
                cmd_parts.append(f"--additional_filters '{filter_expr}'")

        if self.metrics_override:
            # Convert dict to JSON string for command line
            metrics_json = json.dumps(self.metrics_override)
            cmd_parts.append(f"--metrics_override '{metrics_json}'")

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

        Returns:
            Dictionary mapping output types to file paths with standard keys:
            - structures: List of final designed structure files
            - sequences: Path to sequences CSV
            - tables: Analysis metrics and results
            - output_folder: Tool's output directory
        """
        # Ensure file paths are set up
        if not hasattr(self, 'final_ranked_folder') or self.final_ranked_folder is None:
            self._setup_file_paths()

        # BoltzGen generates structures in final_ranked_designs/final_<budget>_designs/
        final_designs_folder = os.path.join(self.final_ranked_folder, f"final_{self.budget}_designs")

        # Structures list file (created by postprocessing script)
        structures_list_file = os.path.join(self.output_folder, "final_structures.txt")

        # Sequences CSV (created by postprocessing script)
        sequences_csv = os.path.join(self.output_folder, "sequences.csv")

        # Structures CSV (created by postprocessing script)
        structures_csv = os.path.join(self.output_folder, "structures.csv")

        # Organize tables by content type
        tables = {
            "all_metrics": TableInfo(
                name="all_metrics",
                path=self.all_designs_metrics_csv,
                columns=[
                    "id", "RMSD", "hydrogen_bonds", "packing_quality",
                    "interface_contacts", "binding_energy", "design_plddt"
                ],
                description="Comprehensive metrics for all generated designs",
                count=None  # Unknown until execution
            ),
            "final_metrics": TableInfo(
                name="final_metrics",
                path=self.final_designs_metrics_csv,
                columns=[
                    "id", "RMSD", "hydrogen_bonds", "packing_quality",
                    "interface_contacts", "binding_energy", "design_plddt", "rank"
                ],
                description=f"Metrics for top {self.budget} designs after filtering",
                count=self.budget
            ),
            "aggregate_metrics": TableInfo(
                name="aggregate_metrics",
                path=self.aggregate_metrics_csv,
                columns=["metric", "mean", "std", "min", "max"],
                description="Aggregate statistics across all designs",
                count=None
            ),
            "per_target_metrics": TableInfo(
                name="per_target_metrics",
                path=self.per_target_metrics_csv,
                columns=["target_id", "num_designs", "avg_rmsd", "avg_plddt"],
                description="Per-target performance metrics",
                count=None
            ),
            "sequences": TableInfo(
                name="sequences",
                path=sequences_csv,
                columns=["id", "sequence"],
                description="Extracted sequences from designed structures",
                count=self.budget
            ),
            "structures": TableInfo(
                name="structures",
                path=structures_csv,
                columns=["id", "pdb"],
                description="List of designed structure files",
                count=self.budget
            )
        }

        return {
            "structures": structures_list_file,  # Will be populated at runtime
            "structure_ids": [],  # Will be populated at runtime
            "sequences": sequences_csv,
            "sequence_ids": [],  # Will be populated at runtime
            "tables": tables,
            "output_folder": self.output_folder,
            "final_designs_folder": final_designs_folder,
            "main": self.final_designs_metrics_csv
        }

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
