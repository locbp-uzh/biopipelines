"""
MobileDistance analysis for calculating distances between mobile regions in APO and HOLO structures.

Analyzes protein structures to calculate sum of Cα distances between APO and HOLO forms
over mobile residues, normalized by √(number of residues).
Outputs CSV with mobile distance metrics for all structure pairs.
"""

import numpy as np
from typing import Dict, List, Any, Optional, Union, Tuple

import os

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo


class MobileDistance(BaseConfig):
    """
    Pipeline tool for analyzing mobile region distances between APO and HOLO structures.

    Takes APO structures and HOLO structures as input and calculates the sum of Cα distances
    between corresponding residues in mobile regions, normalized by √(number of residues).

    Generates CSV with mobile distance metrics for all structure pairs.

    Commonly used for:
    - Conformational change analysis
    - Mobility assessment between APO/HOLO states
    - Binding-induced structural changes
    - Flexibility region identification
    """

    # Tool identification
    TOOL_NAME = "MobileDistance"
    DEFAULT_ENV = "ProteinEnv"
    DEFAULT_RESOURCES = {"gpu": "T4", "memory": "4GB", "time": "1:00:00"}

    def __init__(self,
                 apo_structures: Union[ToolOutput, StandardizedOutput],
                 holo_structures: Union[ToolOutput, StandardizedOutput],
                 mobile: Union[str, ToolOutput],
                 metric_name: str = None,
                 **kwargs):
        """
        Initialize mobile distance analysis tool.

        Args:
            apo_structures: APO structures from previous tool (ToolOutput or StandardizedOutput)
            holo_structures: HOLO structures from previous tool (ToolOutput or StandardizedOutput)
            mobile: Mobile region specification - string or datasheet column reference
                   String format: '10-20+30-40' (residue ranges)
                   Datasheet format: <tool_output>.datasheets.<datasheet_name>.<column_name>
            metric_name: Custom name for the distance column (default: "mobile_distance")
            **kwargs: Additional parameters

        Mobile Selection Syntax:
            String format:
            - '10-20' → residues 10 to 20
            - '10-20+30-40' → residues 10-20 and 30-40
            - '145+147+150' → specific residues 145, 147, and 150

            Datasheet format:
            - boltz_results.datasheets.structures.designed → use 'designed' column values

        Examples:
            # Analyze mobile distance with fixed mobile region
            mobile_analysis = pipeline.add(MobileDistance(
                apo_structures=boltz_apo,
                holo_structures=boltz_holo,
                mobile='10-20+30-40',
                metric_name='conformational_change'
            ))

            # Use mobile regions from datasheet
            mobile_analysis = pipeline.add(MobileDistance(
                apo_structures=boltz_apo,
                holo_structures=boltz_holo,
                mobile=rfdaa.datasheets.structures.designed,
                metric_name='design_mobility'
            ))
        """
        self.apo_input = apo_structures
        self.holo_input = holo_structures
        self.mobile_regions = mobile
        self.custom_metric_name = metric_name

        # Initialize base class
        super().__init__(**kwargs)

        # Set up dependencies
        if hasattr(apo_structures, 'config'):
            self.dependencies.append(apo_structures.config)
        if hasattr(holo_structures, 'config'):
            self.dependencies.append(holo_structures.config)
        if hasattr(mobile, 'config'):
            self.dependencies.append(mobile.config)

    def get_metric_name(self) -> str:
        """Get the default metric name."""
        if self.custom_metric_name:
            return self.custom_metric_name
        return "mobile_distance"

    def get_analysis_csv_path(self) -> str:
        """Get the path for the analysis CSV file - defined once, used everywhere."""
        return os.path.join(self.output_folder, "mobile_distance_analysis.csv")

    def validate_params(self):
        """Validate MobileDistance parameters."""
        if not isinstance(self.apo_input, (ToolOutput, StandardizedOutput)):
            raise ValueError("APO structures input must be a ToolOutput or StandardizedOutput object")

        if not isinstance(self.holo_input, (ToolOutput, StandardizedOutput)):
            raise ValueError("HOLO structures input must be a ToolOutput or StandardizedOutput object")

        if not self.mobile_regions:
            raise ValueError("Mobile regions specification cannot be empty")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input structures from previous tools."""
        self.folders = pipeline_folders

        # Predict APO input structures paths
        self.apo_structures = []
        if hasattr(self.apo_input, 'structures'):
            if isinstance(self.apo_input.structures, list):
                self.apo_structures = self.apo_input.structures
            else:
                self.apo_structures = [self.apo_input.structures]
        elif hasattr(self.apo_input, 'output_folder'):
            # Predict structure files in output folder
            output_folder = self.apo_input.output_folder
            predicted_structures = [
                os.path.join(output_folder, "predicted_structures.pdb"),
                os.path.join(output_folder, "structures.pdb"),
                os.path.join(output_folder, "output.pdb")
            ]
            self.apo_structures = predicted_structures

        # Predict HOLO input structures paths
        self.holo_structures = []
        if hasattr(self.holo_input, 'structures'):
            if isinstance(self.holo_input.structures, list):
                self.holo_structures = self.holo_input.structures
            else:
                self.holo_structures = [self.holo_input.structures]
        elif hasattr(self.holo_input, 'output_folder'):
            # Predict structure files in output folder
            output_folder = self.holo_input.output_folder
            predicted_structures = [
                os.path.join(output_folder, "predicted_structures.pdb"),
                os.path.join(output_folder, "structures.pdb"),
                os.path.join(output_folder, "output.pdb")
            ]
            self.holo_structures = predicted_structures

        if not self.apo_structures:
            raise ValueError(f"Could not predict APO structure paths from: {self.apo_input}")
        if not self.holo_structures:
            raise ValueError(f"Could not predict HOLO structure paths from: {self.holo_input}")

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        mobile_display = str(self.mobile_regions)
        if hasattr(self.mobile_regions, 'datasheet_path'):
            mobile_display = f"From datasheet: {self.mobile_regions.datasheet_path}"

        config_lines.extend([
            f"APO STRUCTURES: {len(getattr(self, 'apo_structures', []))} files",
            f"HOLO STRUCTURES: {len(getattr(self, 'holo_structures', []))} files",
            f"MOBILE REGIONS: {mobile_display}",
            f"OUTPUT METRIC: {self.get_metric_name()}"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """
        Generate mobile distance analysis execution script.

        Args:
            script_path: Path where script should be written

        Returns:
            Script content as string
        """
        runtime_folder = os.path.dirname(script_path)
        output_folder = self.output_folder
        os.makedirs(output_folder, exist_ok=True)

        # Output CSV path - defined once in get_analysis_csv_path()
        analysis_csv = self.get_analysis_csv_path()

        # Create config file for mobile distance calculation
        config_file = os.path.join(output_folder, "mobile_distance_config.json")

        # Handle mobile regions input
        mobile_config = None
        if isinstance(self.mobile_regions, str):
            mobile_config = {"type": "fixed", "value": self.mobile_regions}
        else:
            # Datasheet reference
            mobile_config = {
                "type": "datasheet",
                "datasheet_path": getattr(self.mobile_regions, 'datasheet_path', ''),
                "column_name": getattr(self.mobile_regions, 'column_name', 'mobile')
            }

        config_data = {
            "apo_structures": self.apo_structures,
            "holo_structures": self.holo_structures,
            "mobile_regions": mobile_config,
            "metric_name": self.get_metric_name(),
            "output_csv": analysis_csv
        }

        import json
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        # Generate script content
        script_content = f"""#!/bin/bash
# MobileDistance execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Running mobile distance analysis"
echo "APO structures: {len(self.apo_structures)}"
echo "HOLO structures: {len(self.holo_structures)}"
echo "Mobile regions: {self.mobile_regions}"
echo "Output: {analysis_csv}"

# Run Python analysis script
python "{os.path.join(self.folders['HelpScripts'], 'pipe_mobile_distance.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Mobile distance analysis completed successfully"
    echo "Results written to: {analysis_csv}"
else
    echo "Error: Mobile distance analysis failed"
    exit 1
fi

{self.generate_completion_check_footer()}
"""

        return script_content

    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after mobile distance analysis.

        Returns:
            Dictionary with output file paths
        """
        analysis_csv = self.get_analysis_csv_path()

        datasheets = {
            "mobile_analysis": DatasheetInfo(
                name="mobile_analysis",
                path=analysis_csv,
                columns=["id", "apo_structure", "holo_structure", "mobile", self.get_metric_name()],
                description=f"Mobile distance analysis between APO and HOLO structures",
                count=len(getattr(self, 'apo_structures', []))
            )
        }

        return {
            "structures": [],
            "structure_ids": [],
            "compounds": [],
            "compound_ids": [],
            "sequences": [],
            "sequence_ids": [],
            "datasheets": datasheets,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                "mobile_regions": str(self.mobile_regions),
                "metric_name": self.get_metric_name()
            }
        })
        return base_dict