"""
RF3 (RoseTTAFold3) configuration for biomolecular structure prediction.

Handles protein and protein-ligand complex prediction with batch processing
and MSA caching support.
"""

import os
import json
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo


class RF3(BaseConfig):
    """
    RF3 (RoseTTAFold3) configuration for biomolecular structure prediction.

    Supports protein and protein-ligand complex prediction with batch processing.
    """

    # Tool identification
    TOOL_NAME = "RF3"
    DEFAULT_ENV = None  # Loaded from config.yaml

    def __init__(self,
                 # Primary input parameters
                 proteins: Union[str, List[str], ToolOutput, StandardizedOutput] = None,
                 ligands: Union[str, ToolOutput, StandardizedOutput, None] = None,
                 msas: Optional[Union[str, ToolOutput]] = None,
                 # Core prediction parameters
                 output_format: str = "pdb",
                 checkpoint_path: Optional[str] = None,
                 early_stopping_plddt: Optional[float] = None,
                 # Advanced parameters
                 use_templates: bool = False,
                 **kwargs):
        """
        Initialize RF3 configuration.

        Args:
            proteins: Protein sequences - can be ToolOutput, StandardizedOutput, file path, or direct sequence
            ligands: Single ligand SMILES string, ToolOutput with compounds, or datasheet reference
            msas: MSA files for recycling (optional)
            output_format: Output format ("pdb" or "cif")
            checkpoint_path: Path to RF3 checkpoint file
            early_stopping_plddt: pLDDT threshold for early stopping
            use_templates: Enable template-based prediction
            **kwargs: Additional parameters
        """
        # Store RF3-specific parameters
        self.proteins = proteins
        self.ligands = ligands
        self.msas = msas
        self.output_format = output_format.lower()
        self.checkpoint_path = checkpoint_path
        self.early_stopping_plddt = early_stopping_plddt
        self.use_templates = use_templates

        # Input tracking
        self.input_sequences = None
        self.input_compounds = []
        self.input_datasheets = {}
        self.input_is_tool_output = False
        self.standardized_input = None

        # Initialize base class
        super().__init__(**kwargs)

        # Process inputs
        self._process_inputs()

        # Validate output format
        if self.output_format not in ["pdb", "cif"]:
            raise ValueError(f"Invalid output_format: {self.output_format}. Must be 'pdb' or 'cif'")

    def _process_inputs(self):
        """Process protein and ligand inputs."""
        # Handle protein inputs
        if isinstance(self.proteins, StandardizedOutput):
            self.input_sequences = self.proteins.sequences
            self.input_datasheets = getattr(self.proteins, 'datasheets', {})
            self.standardized_input = self.proteins
        elif isinstance(self.proteins, ToolOutput):
            self.input_sequences = self.proteins
            self.input_datasheets = self.proteins.get_output_files("datasheets")
            self.input_is_tool_output = True
            self.dependencies.append(self.proteins.config)
        elif isinstance(self.proteins, str):
            # Single sequence or file path
            self.input_sequences = self.proteins
        elif isinstance(self.proteins, list):
            # List of sequences or file paths
            self.input_sequences = self.proteins

        # Handle ligand inputs
        if self.ligands:
            if isinstance(self.ligands, StandardizedOutput):
                # Get compounds from StandardizedOutput
                self.input_compounds = getattr(self.ligands, 'compounds', [])
                self.input_datasheets = getattr(self.ligands, 'datasheets', {})
            elif isinstance(self.ligands, ToolOutput):
                # Get compounds from ToolOutput
                self.input_compounds = self.ligands.get_output_files("compounds")
                self.input_datasheets = self.ligands.get_output_files("datasheets")
                self.dependencies.append(self.ligands.config)
            elif isinstance(self.ligands, str):
                # Single SMILES string or file path
                self.input_compounds.append(self.ligands)

    def validate_params(self):
        """Validate RF3 parameters."""
        if not self.proteins:
            raise ValueError("proteins parameter is required")

        if self.output_format not in ["pdb", "cif"]:
            raise ValueError("output_format must be 'pdb' or 'cif'")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input parameters."""
        self.folders = pipeline_folders

        # Set default checkpoint path if not provided
        if self.checkpoint_path is None:
            self.checkpoint_path = os.path.join(pipeline_folders["ModelForge"], "rf3_latest.pt")

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"PROTEINS: {self._format_input_display(self.proteins)}",
            f"LIGANDS: {self._format_input_display(self.ligands) if self.ligands else 'None (apo prediction)'}",
            f"OUTPUT_FORMAT: {self.output_format}",
            f"USE_TEMPLATES: {self.use_templates}"
        ])

        if self.early_stopping_plddt:
            config_lines.append(f"EARLY_STOPPING_PLDDT: {self.early_stopping_plddt}")

        return config_lines

    def _format_input_display(self, input_val):
        """Format input value for display."""
        if isinstance(input_val, (ToolOutput, StandardizedOutput)):
            return f"<{type(input_val).__name__}>"
        elif isinstance(input_val, list):
            return f"List[{len(input_val)} items]"
        elif isinstance(input_val, str):
            return input_val if len(input_val) < 50 else f"{input_val[:47]}..."
        return str(input_val)

    def generate_script(self, script_path: str) -> str:
        """
        Generate RF3 execution script.

        Args:
            script_path: Path where script should be written

        Returns:
            Script content as string
        """
        runtime_folder = os.path.dirname(script_path)
        output_folder = self.output_folder
        os.makedirs(output_folder, exist_ok=True)

        # Create configuration file for RF3 processing
        config_file = os.path.join(output_folder, "rf3_config.json")
        config_data = {
            "proteins": self._serialize_input(self.proteins),
            "ligands": self._serialize_input(self.ligands) if self.ligands else None,
            "msas": self._serialize_input(self.msas) if self.msas else None,
            "output_format": self.output_format,
            "checkpoint_path": self.checkpoint_path,
            "early_stopping_plddt": self.early_stopping_plddt,
            "use_templates": self.use_templates,
            "output_folder": output_folder,
            "input_datasheets": {k: v.path if hasattr(v, 'path') else str(v)
                                for k, v in self.input_datasheets.items()}
        }

        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        # Generate script
        script_content = f"""#!/bin/bash
# RF3 execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Running RF3 structure prediction"
echo "Output folder: {output_folder}"

# Run RF3 processing script
python "{os.path.join(self.folders['HelpScripts'], 'pipe_rf3.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "RF3 prediction completed successfully"
else
    echo "Error: RF3 prediction failed"
    exit 1
fi

{self.generate_completion_check_footer()}
"""

        return script_content

    def _serialize_input(self, input_val):
        """Serialize input value for JSON."""
        if isinstance(input_val, (ToolOutput, StandardizedOutput)):
            if hasattr(input_val, 'datasheets'):
                return {"type": "tool_output", "datasheets": list(input_val.datasheets.keys())}
            return {"type": "tool_output"}
        elif isinstance(input_val, list):
            return {"type": "list", "values": input_val}
        elif isinstance(input_val, str):
            return {"type": "string", "value": input_val}
        return None

    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after RF3 prediction.

        Returns:
            Dictionary mapping output types to file paths
        """
        # Predict structure file names based on inputs
        structure_files = []

        # For now, return generic structure pattern
        # This will be populated by the runtime script
        structures_csv = os.path.join(self.output_folder, "structures.csv")
        confidence_csv = os.path.join(self.output_folder, "confidence_scores.csv")

        datasheets = {}

        if os.path.exists(structures_csv):
            datasheets["structures"] = DatasheetInfo(
                name="structures",
                path=structures_csv,
                columns=["id", "model_id", "file_path", "plddt_score"],
                description="Predicted structures with confidence scores"
            )

        if os.path.exists(confidence_csv):
            datasheets["confidence"] = DatasheetInfo(
                name="confidence",
                path=confidence_csv,
                columns=["id", "model_id", "plddt_score", "ptm_score"],
                description="Per-model confidence metrics"
            )

        outputs = {
            "structures": [structures_csv],
            "datasheets": datasheets,
            "output_folder": self.output_folder
        }

        return outputs

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including RF3-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "rf3_params": {
                "output_format": self.output_format,
                "early_stopping_plddt": self.early_stopping_plddt,
                "use_templates": self.use_templates
            }
        })
        return base_dict
