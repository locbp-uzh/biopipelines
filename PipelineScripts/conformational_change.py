"""
ConformationalChange analysis for measuring structural changes between reference and target structures.

Analyzes protein structures to quantify conformational changes through multiple metrics including
RMSD, maximum distance, mean distance, and sum over square root normalization.
Supports structural alignment using PyMOL's align, super, or cealign methods.
"""

import numpy as np
from typing import Dict, List, Any, Optional, Union, Tuple

import os

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo


class ConformationalChange(BaseConfig):
    """
    Pipeline tool for analyzing conformational changes between reference and target structures.

    Takes reference structures and target structures as input, aligns them using specified method,
    and calculates multiple metrics to quantify conformational changes in selected regions.

    Generates CSV with comprehensive conformational change analysis.

    Commonly used for:
    - Conformational change analysis between different states
    - Binding-induced structural changes
    - Flexibility region identification
    - Comparison of predicted vs experimental structures
    """

    # Tool identification
    TOOL_NAME = "ConformationalChange"
    DEFAULT_ENV = None  # Loaded from config.yaml

    def __init__(self,
                 reference_structures: Union[str, ToolOutput, StandardizedOutput],
                 target_structures: Union[ToolOutput, StandardizedOutput],
                 selection: Union[str, ToolOutput, None] = None,
                 alignment: str = "align",
                 **kwargs):
        """
        Initialize conformational change analysis tool.

        Args:
            reference_structures: Reference structures - can be:
                                - Single PDB file path (string)
                                - ToolOutput/StandardizedOutput (single or multiple structures)
            target_structures: Target structures from previous tool (ToolOutput or StandardizedOutput)
            selection: Region specification - string, table column reference, or None for all residues
                      String format: '10-20+30-40' (residue ranges)
                      Table format: <tool_output>.tables.<table_name>.<column_name>
                      None: Compare all CA atoms (whole structure RMSD)
            alignment: Alignment method - "align", "super", or "cealign" (default: "align")
            **kwargs: Additional parameters

        Selection Syntax:
            String format:
            - '10-20' → residues 10 to 20
            - '10-20+30-40' → residues 10-20 and 30-40
            - '145+147+150' → specific residues 145, 147, and 150
            - None → all CA atoms

            Table format:
            - boltz_results.tables.structures.designed → use 'designed' column values

        Alignment Methods:
            - "align": PyMOL align (sequence-dependent, fast)
            - "super": PyMOL super (structure-based superposition)
            - "cealign": PyMOL cealign (combinatorial extension alignment)

        Examples:
            # Whole structure RMSD (no selection)
            conf_analysis = pipeline.add(ConformationalChange(
                reference_structures=design,
                target_structures=refolded,
                alignment='cealign'
            ))

            # Analyze conformational change with fixed selection
            conf_analysis = pipeline.add(ConformationalChange(
                reference_structures="reference.pdb",
                target_structures=boltz_results,
                selection='10-20+30-40',
                alignment='super'
            ))

            # Use selection from table
            conf_analysis = pipeline.add(ConformationalChange(
                reference_structures=boltz_ref,
                target_structures=boltz_targets,
                selection=rfdaa.tables.structures.designed,
                alignment='cealign'
            ))
        """
        self.reference_input = reference_structures
        self.target_input = target_structures
        self.selection_spec = selection
        self.alignment_method = alignment

        # Initialize base class
        super().__init__(**kwargs)

        # Set up dependencies
        if hasattr(reference_structures, 'config'):
            self.dependencies.append(reference_structures.config)
        if hasattr(target_structures, 'config'):
            self.dependencies.append(target_structures.config)
        if hasattr(selection, 'config'):
            self.dependencies.append(selection.config)

    def get_analysis_csv_path(self) -> str:
        """Get the path for the analysis CSV file - defined once, used everywhere."""
        return os.path.join(self.output_folder, "conformational_change_analysis.csv")

    def validate_params(self):
        """Validate ConformationalChange parameters."""
        if not self.target_input:
            raise ValueError("Target structures input cannot be empty")

        if not isinstance(self.target_input, (ToolOutput, StandardizedOutput)):
            raise ValueError("Target structures input must be a ToolOutput or StandardizedOutput object")

        # selection can be None (whole structure RMSD), string, or table reference

        if self.alignment_method not in ["align", "super", "cealign"]:
            raise ValueError(f"Alignment method must be 'align', 'super', or 'cealign', got: {self.alignment_method}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input structures from previous tools."""
        self.folders = pipeline_folders

        # Predict reference input structures paths
        self.reference_structures = []
        if isinstance(self.reference_input, str):
            # Single PDB file path
            self.reference_structures = [self.reference_input]
        elif hasattr(self.reference_input, 'structures'):
            if isinstance(self.reference_input.structures, list):
                self.reference_structures = self.reference_input.structures
            else:
                self.reference_structures = [self.reference_input.structures]
        elif hasattr(self.reference_input, 'output_folder'):
            # Predict structure files in output folder
            output_folder = self.reference_input.output_folder
            predicted_structures = [
                os.path.join(output_folder, "predicted_structures.pdb"),
                os.path.join(output_folder, "structures.pdb"),
                os.path.join(output_folder, "output.pdb")
            ]
            self.reference_structures = predicted_structures

        # Predict target input structures paths
        self.target_structures = []
        if hasattr(self.target_input, 'structures'):
            if isinstance(self.target_input.structures, list):
                self.target_structures = self.target_input.structures
            else:
                self.target_structures = [self.target_input.structures]
        elif hasattr(self.target_input, 'output_folder'):
            # Predict structure files in output folder
            output_folder = self.target_input.output_folder
            predicted_structures = [
                os.path.join(output_folder, "predicted_structures.pdb"),
                os.path.join(output_folder, "structures.pdb"),
                os.path.join(output_folder, "output.pdb")
            ]
            self.target_structures = predicted_structures

        if not self.reference_structures:
            raise ValueError(f"Could not predict reference structure paths from: {self.reference_input}")
        if not self.target_structures:
            raise ValueError(f"Could not predict target structure paths from: {self.target_input}")

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        if self.selection_spec is None:
            selection_display = "All CA atoms (whole structure)"
        elif hasattr(self.selection_spec, 'table_path'):
            selection_display = f"From table: {self.selection_spec.table_path}"
        else:
            selection_display = str(self.selection_spec)

        config_lines.extend([
            f"REFERENCE STRUCTURES: {len(getattr(self, 'reference_structures', []))} files",
            f"TARGET STRUCTURES: {len(getattr(self, 'target_structures', []))} files",
            f"SELECTION: {selection_display}",
            f"ALIGNMENT METHOD: {self.alignment_method}",
            f"METRICS: RMSD, max_distance, mean_distance, sum_over_square_root"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """
        Generate conformational change analysis execution script.

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

        # Create config file for conformational change calculation
        config_file = os.path.join(output_folder, "conformational_change_config.json")

        # Handle selection input
        selection_config = None
        if self.selection_spec is None:
            # No selection = compare all CA atoms
            selection_config = {"type": "all"}
        elif isinstance(self.selection_spec, str):
            selection_config = {"type": "fixed", "value": self.selection_spec}
        elif isinstance(self.selection_spec, tuple) and len(self.selection_spec) == 2:
            # Tuple format: (TableInfo, column_name)
            table_info, column_name = self.selection_spec
            selection_config = {
                "type": "table",
                "table_path": table_info.path if hasattr(table_info, 'path') else '',
                "column_name": column_name
            }
        else:
            # Table reference with attributes
            selection_config = {
                "type": "table",
                "table_path": getattr(self.selection_spec, 'table_path', ''),
                "column_name": getattr(self.selection_spec, 'column_name', 'selection')
            }

        config_data = {
            "reference_structures": self.reference_structures,
            "target_structures": self.target_structures,
            "selection": selection_config,
            "alignment_method": self.alignment_method,
            "output_csv": analysis_csv
        }

        import json
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        # Generate script content
        script_content = f"""#!/bin/bash
# ConformationalChange execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Running conformational change analysis"
echo "Reference structures: {len(self.reference_structures)}"
echo "Target structures: {len(self.target_structures)}"
echo "Selection: {self.selection_spec}"
echo "Alignment method: {self.alignment_method}"
echo "Output: {analysis_csv}"

# Run Python analysis script
python "{os.path.join(self.folders['HelpScripts'], 'pipe_conformational_change.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Conformational change analysis completed successfully"
    echo "Results written to: {analysis_csv}"
else
    echo "Error: Conformational change analysis failed"
    exit 1
fi

{self.generate_completion_check_footer()}
"""

        return script_content

    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after conformational change analysis.

        Returns:
            Dictionary with output file paths
        """
        analysis_csv = self.get_analysis_csv_path()

        tables = {
            "conformational_analysis": TableInfo(
                name="conformational_analysis",
                path=analysis_csv,
                columns=["id", "reference_structure", "target_structure", "selection",
                        "num_residues", "RMSD", "max_distance", "mean_distance", "sum_over_square_root"],
                description=f"Conformational change analysis between reference and target structures",
                count=len(getattr(self, 'target_structures', []))
            )
        }

        return {
            "structures": [],
            "structure_ids": [],
            "compounds": [],
            "compound_ids": [],
            "sequences": [],
            "sequence_ids": [],
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                "selection": str(self.selection_spec),
                "alignment_method": self.alignment_method
            }
        })
        return base_dict
