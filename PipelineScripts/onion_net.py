"""
OnionNet and OnionNet-2 configurations for protein-ligand binding affinity prediction.

Predicts binding affinities from protein-ligand complex structures using CNN-based models.
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


class OnionNet(BaseConfig):
    """
    OnionNet configuration for protein-ligand binding affinity prediction.

    Uses CNN-based approach with rotation-free element-pair-specific contacts.
    """

    # Tool identification
    TOOL_NAME = "OnionNet"
    DEFAULT_ENV = None  # Loaded from config.yaml

    def __init__(self,
                 # Primary input parameters
                 structures: Union[str, ToolOutput, StandardizedOutput] = None,
                 # Model parameters
                 model_weights: Optional[str] = None,
                 scaler_model: Optional[str] = None,
                 # Output parameters
                 output_format: str = "csv",
                 **kwargs):
        """
        Initialize OnionNet configuration.

        Args:
            structures: Protein-ligand complex structures (PDB files or ToolOutput)
            model_weights: Path to model weights file (.h5)
            scaler_model: Path to scaler model file
            output_format: Output format for predictions (default: "csv")
            **kwargs: Additional parameters
        """
        # Store OnionNet-specific parameters
        self.structures = structures
        self.model_weights = model_weights
        self.scaler_model = scaler_model
        self.output_format = output_format.lower()

        # Input tracking
        self.input_structures = []
        self.input_tables = {}
        self.input_is_tool_output = False
        self.standardized_input = None

        # Initialize base class
        super().__init__(**kwargs)

        # Process inputs
        self._process_inputs()

    def _process_inputs(self):
        """Process structure inputs."""
        if isinstance(self.structures, StandardizedOutput):
            self.input_structures = getattr(self.structures, 'structures', [])
            self.input_tables = getattr(self.structures, 'tables', {})
            self.standardized_input = self.structures
        elif isinstance(self.structures, ToolOutput):
            self.input_structures = self.structures.get_output_files("structures")
            self.input_tables = self.structures.get_output_files("tables")
            self.input_is_tool_output = True
            self.dependencies.append(self.structures.config)
        elif isinstance(self.structures, str):
            # Single structure file or directory
            self.input_structures = [self.structures]
        elif isinstance(self.structures, list):
            # List of structure files
            self.input_structures = self.structures

    def validate_params(self):
        """Validate OnionNet parameters."""
        if not self.structures:
            raise ValueError("structures parameter is required")

        if self.output_format not in ["csv", "json"]:
            raise ValueError("output_format must be 'csv' or 'json'")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input parameters."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"STRUCTURES: {self._format_input_display(self.structures)}",
            f"MODEL_WEIGHTS: {self.model_weights if self.model_weights else 'Default'}",
            f"OUTPUT_FORMAT: {self.output_format}"
        ])

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
        Generate OnionNet execution script.

        Args:
            script_path: Path where script should be written

        Returns:
            Script content as string
        """
        runtime_folder = os.path.dirname(script_path)
        output_folder = self.output_folder
        os.makedirs(output_folder, exist_ok=True)

        # Create configuration file for OnionNet processing
        config_file = os.path.join(output_folder, "onionnet_config.json")
        config_data = {
            "structures": self._serialize_input(self.structures),
            "model_weights": self.model_weights,
            "scaler_model": self.scaler_model,
            "output_format": self.output_format,
            "output_folder": output_folder,
            "input_tables": {k: v.path if hasattr(v, 'path') else str(v)
                                for k, v in self.input_tables.items()}
        }

        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        # Generate script
        script_content = f"""#!/bin/bash
# OnionNet execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Running OnionNet binding affinity prediction"
echo "Output folder: {output_folder}"

# Run OnionNet processing script
python "{os.path.join(self.folders['HelpScripts'], 'pipe_onion_net.py')}" \\
  --config "{config_file}" \\
  --version 1

if [ $? -eq 0 ]; then
    echo "OnionNet prediction completed successfully"
else
    echo "Error: OnionNet prediction failed"
    exit 1
fi

{self.generate_completion_check_footer()}
"""

        return script_content

    def _serialize_input(self, input_val):
        """Serialize input value for JSON."""
        if isinstance(input_val, (ToolOutput, StandardizedOutput)):
            if hasattr(input_val, 'tables'):
                return {"type": "tool_output", "tables": list(input_val.tables.keys())}
            return {"type": "tool_output"}
        elif isinstance(input_val, list):
            return {"type": "list", "values": input_val}
        elif isinstance(input_val, str):
            return {"type": "string", "value": input_val}
        return None

    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after OnionNet prediction.

        Returns:
            Dictionary mapping output types to file paths
        """
        predictions_csv = os.path.join(self.output_folder, "affinity_predictions.csv")

        tables = {}

        if os.path.exists(predictions_csv):
            tables["affinities"] = TableInfo(
                name="affinities",
                path=predictions_csv,
                columns=["id", "structure_path", "predicted_affinity_pKa"],
                description="Predicted protein-ligand binding affinities in pKa scale"
            )

        outputs = {
            "predictions": [predictions_csv],
            "tables": tables,
            "output_folder": self.output_folder
        }

        return outputs

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including OnionNet-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "onionnet_params": {
                "model_weights": self.model_weights,
                "scaler_model": self.scaler_model,
                "output_format": self.output_format
            }
        })
        return base_dict


class OnionNet2(BaseConfig):
    """
    OnionNet-2 configuration for protein-ligand binding affinity prediction.

    Improved version with higher accuracy and lower computational cost.
    Uses residue-atom contacting shells in CNN architecture.
    """

    # Tool identification
    TOOL_NAME = "OnionNet2"
    DEFAULT_ENV = None  # Loaded from config.yaml

    def __init__(self,
                 # Primary input parameters
                 structures: Union[str, ToolOutput, StandardizedOutput] = None,
                 # Model parameters
                 model_path: Optional[str] = None,
                 scaler_path: Optional[str] = None,
                 shells: int = 62,
                 # Output parameters
                 output_format: str = "csv",
                 **kwargs):
        """
        Initialize OnionNet-2 configuration.

        Args:
            structures: Protein-ligand complex structures (PDB files or ToolOutput)
            model_path: Path to trained model file
            scaler_path: Path to scaler file
            shells: Number of contacting shells (default: 62)
            output_format: Output format for predictions (default: "csv")
            **kwargs: Additional parameters
        """
        # Store OnionNet2-specific parameters
        self.structures = structures
        self.model_path = model_path
        self.scaler_path = scaler_path
        self.shells = shells
        self.output_format = output_format.lower()

        # Input tracking
        self.input_structures = []
        self.input_tables = {}
        self.input_is_tool_output = False
        self.standardized_input = None

        # Initialize base class
        super().__init__(**kwargs)

        # Process inputs
        self._process_inputs()

    def _process_inputs(self):
        """Process structure inputs."""
        if isinstance(self.structures, StandardizedOutput):
            self.input_structures = getattr(self.structures, 'structures', [])
            self.input_tables = getattr(self.structures, 'tables', {})
            self.standardized_input = self.structures
        elif isinstance(self.structures, ToolOutput):
            self.input_structures = self.structures.get_output_files("structures")
            self.input_tables = self.structures.get_output_files("tables")
            self.input_is_tool_output = True
            self.dependencies.append(self.structures.config)
        elif isinstance(self.structures, str):
            # Single structure file or directory
            self.input_structures = [self.structures]
        elif isinstance(self.structures, list):
            # List of structure files
            self.input_structures = self.structures

    def validate_params(self):
        """Validate OnionNet-2 parameters."""
        if not self.structures:
            raise ValueError("structures parameter is required")

        if self.shells < 1:
            raise ValueError("shells must be >= 1")

        if self.output_format not in ["csv", "json"]:
            raise ValueError("output_format must be 'csv' or 'json'")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input parameters."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"STRUCTURES: {self._format_input_display(self.structures)}",
            f"MODEL_PATH: {self.model_path if self.model_path else 'Default'}",
            f"SHELLS: {self.shells}",
            f"OUTPUT_FORMAT: {self.output_format}"
        ])

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
        Generate OnionNet-2 execution script.

        Args:
            script_path: Path where script should be written

        Returns:
            Script content as string
        """
        runtime_folder = os.path.dirname(script_path)
        output_folder = self.output_folder
        os.makedirs(output_folder, exist_ok=True)

        # Create configuration file for OnionNet-2 processing
        config_file = os.path.join(output_folder, "onionnet2_config.json")
        config_data = {
            "structures": self._serialize_input(self.structures),
            "model_path": self.model_path,
            "scaler_path": self.scaler_path,
            "shells": self.shells,
            "output_format": self.output_format,
            "output_folder": output_folder,
            "input_tables": {k: v.path if hasattr(v, 'path') else str(v)
                                for k, v in self.input_tables.items()}
        }

        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        # Generate script
        script_content = f"""#!/bin/bash
# OnionNet-2 execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Running OnionNet-2 binding affinity prediction"
echo "Output folder: {output_folder}"
echo "Number of shells: {self.shells}"

# Run OnionNet-2 processing script
python "{os.path.join(self.folders['HelpScripts'], 'pipe_onion_net.py')}" \\
  --config "{config_file}" \\
  --version 2

if [ $? -eq 0 ]; then
    echo "OnionNet-2 prediction completed successfully"
else
    echo "Error: OnionNet-2 prediction failed"
    exit 1
fi

{self.generate_completion_check_footer()}
"""

        return script_content

    def _serialize_input(self, input_val):
        """Serialize input value for JSON."""
        if isinstance(input_val, (ToolOutput, StandardizedOutput)):
            if hasattr(input_val, 'tables'):
                return {"type": "tool_output", "tables": list(input_val.tables.keys())}
            return {"type": "tool_output"}
        elif isinstance(input_val, list):
            return {"type": "list", "values": input_val}
        elif isinstance(input_val, str):
            return {"type": "string", "value": input_val}
        return None

    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after OnionNet-2 prediction.

        Returns:
            Dictionary mapping output types to file paths
        """
        predictions_csv = os.path.join(self.output_folder, "affinity_predictions.csv")

        tables = {}

        if os.path.exists(predictions_csv):
            tables["affinities"] = TableInfo(
                name="affinities",
                path=predictions_csv,
                columns=["id", "structure_path", "predicted_affinity_pKa"],
                description="Predicted protein-ligand binding affinities in pKa scale (OnionNet-2)"
            )

        outputs = {
            "predictions": [predictions_csv],
            "tables": tables,
            "output_folder": self.output_folder
        }

        return outputs

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including OnionNet-2-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "onionnet2_params": {
                "model_path": self.model_path,
                "scaler_path": self.scaler_path,
                "shells": self.shells,
                "output_format": self.output_format
            }
        })
        return base_dict
