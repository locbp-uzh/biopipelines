# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
OnionNet and OnionNet-2 configurations for protein-ligand binding affinity prediction.

Predicts binding affinities from protein-ligand complex structures using CNN-based models.
"""

import os
import json
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream

class OnionNet2(BaseConfig):
    """
    OnionNet-2 configuration for protein-ligand binding affinity prediction.

    Improved version with higher accuracy and lower computational cost.
    Uses residue-atom contacting shells in CNN architecture.
    """

    TOOL_NAME = "OnionNet2"

    # Lazy path descriptors
    config_file = Path(lambda self: os.path.join(self.output_folder, "onionnet2_config.json"))
    predictions_csv = Path(lambda self: os.path.join(self.output_folder, "affinity_predictions.csv"))
    helper_script = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_onion_net.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 model_path: Optional[str] = None,
                 scaler_path: Optional[str] = None,
                 shells: int = 62,
                 output_format: str = "csv",
                 **kwargs):
        """
        Initialize OnionNet-2 configuration.

        Args:
            structures: Protein-ligand complex structures as DataStream or StandardizedOutput
            model_path: Path to trained model file
            scaler_path: Path to scaler file
            shells: Number of contacting shells (default: 62)
            output_format: Output format for predictions (default: "csv")
        """
        # Resolve input to DataStream
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        # Store OnionNet2-specific parameters
        self.model_path = model_path
        self.scaler_path = scaler_path
        self.shells = shells
        self.output_format = output_format.lower()

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate OnionNet-2 parameters."""
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")

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
            f"STRUCTURES: {len(self.structures_stream)} complexes",
            f"MODEL_PATH: {self.model_path or 'Default'}",
            f"SHELLS: {self.shells}",
            f"OUTPUT_FORMAT: {self.output_format}"
        ])
        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate OnionNet-2 execution script."""
        script_content = "#!/bin/bash\n"
        script_content += "# OnionNet-2 execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self._generate_script_run_onionnet2()
        script_content += self.generate_completion_check_footer()

        return script_content

    def _generate_script_run_onionnet2(self) -> str:
        """Generate the OnionNet-2 execution part of the script."""
        # Create configuration data
        config_data = {
            "structure_files": self.structures_stream.files,
            "structure_ids": self.structures_stream.ids,
            "model_path": self.model_path,
            "scaler_path": self.scaler_path,
            "shells": self.shells,
            "output_format": self.output_format,
            "output_folder": self.output_folder
        }

        # Write config file at pipeline time
        os.makedirs(self.output_folder, exist_ok=True)
        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        return f"""echo "Running OnionNet-2 binding affinity prediction"
echo "Output folder: {self.output_folder}"
echo "Number of shells: {self.shells}"

# Run OnionNet-2 processing script
python "{self.helper_script}" \\
  --config "{self.config_file}" \\
  --version 2

if [ $? -eq 0 ]; then
    echo "OnionNet-2 prediction completed successfully"
else
    echo "Error: OnionNet-2 prediction failed"
    exit 1
fi

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after OnionNet-2 prediction."""
        structure_ids = self.structures_stream.ids

        tables = {
            "affinities": TableInfo(
                name="affinities",
                path=self.predictions_csv,
                columns=["id", "structure_path", "predicted_affinity_pKa"],
                description="Predicted protein-ligand binding affinities in pKa scale (OnionNet-2)",
                count=len(structure_ids)
            )
        }

        return {
            "structures": DataStream.empty("structures", "pdb"),
            "sequences": DataStream.empty("sequences", "fasta"),
            "compounds": DataStream.empty("compounds", "sdf"),
            "tables": tables,
            "output_folder": self.output_folder
        }

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
