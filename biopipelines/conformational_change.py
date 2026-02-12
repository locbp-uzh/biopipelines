# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
ConformationalChange analysis for measuring structural changes between reference and target structures.

Analyzes protein structures to quantify conformational changes via RMSD computed
by PyMOL's align, super, or cealign methods.
"""

import os
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

    @classmethod
    def _install_script(cls, folders, env_manager="mamba"):
        return """echo "=== ConformationalChange ==="
echo "Requires ProteinEnv (installed with PyMOL.install())"
echo "No additional installation needed."
echo "=== ConformationalChange ready ==="
"""

    # Lazy path descriptors
    analysis_csv = Path(lambda self: os.path.join(self.output_folder, "conformational_change_analysis.csv"))
    config_file = Path(lambda self: os.path.join(self.output_folder, "conformational_change_config.json"))
    reference_ds_json = Path(lambda self: os.path.join(self.output_folder, "reference_structures.json"))
    target_ds_json = Path(lambda self: os.path.join(self.output_folder, "target_structures.json"))
    analysis_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_conformational_change.py"))

    def __init__(self,
                 reference_structures: Union[DataStream, StandardizedOutput],
                 target_structures: Union[DataStream, StandardizedOutput],
                 selection: Optional[str] = None,
                 alignment: str = "align",
                 **kwargs):
        """
        Initialize conformational change analysis tool.

        Args:
            reference_structures: Reference structures as DataStream or StandardizedOutput
            target_structures: Target structures as DataStream or StandardizedOutput
            selection: Region specification - string or None for all residues
                      String format: '10-20+30-40' (residue ranges)
                      None: Compare all CA atoms (whole structure RMSD)
            alignment: Alignment method - "align", "super", or "cealign" (default: "align")
            **kwargs: Additional parameters

        Selection Syntax:
            - '10-20' → residues 10 to 20
            - '10-20+30-40' → residues 10-20 and 30-40
            - '145+147+150' → specific residues 145, 147, and 150
            - None → all CA atoms

        Alignment Methods:
            - "align": PyMOL align (sequence-dependent, fast)
            - "super": PyMOL super (structure-based superposition)
            - "cealign": PyMOL cealign (combinatorial extension alignment)

        Examples:
            # Whole structure RMSD (no selection)
            conf_analysis = ConformationalChange(
                reference_structures=design,
                target_structures=refolded,
                alignment='cealign'
            )

            # Analyze conformational change with fixed selection
            conf_analysis = ConformationalChange(
                reference_structures=rfdaa,
                target_structures=boltz_results,
                selection='10-20+30-40',
                alignment='super'
            )
        """
        # Resolve reference structures to DataStream
        if isinstance(reference_structures, StandardizedOutput):
            self.reference_stream: DataStream = reference_structures.streams.structures
        elif isinstance(reference_structures, DataStream):
            self.reference_stream = reference_structures
        else:
            raise ValueError(f"reference_structures must be DataStream or StandardizedOutput, got {type(reference_structures)}")

        # Resolve target structures to DataStream
        if isinstance(target_structures, StandardizedOutput):
            self.target_stream: DataStream = target_structures.streams.structures
        elif isinstance(target_structures, DataStream):
            self.target_stream = target_structures
        else:
            raise ValueError(f"target_structures must be DataStream or StandardizedOutput, got {type(target_structures)}")

        self.selection_spec = selection
        self.alignment_method = alignment

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate ConformationalChange parameters."""
        if not self.reference_stream or len(self.reference_stream) == 0:
            raise ValueError("reference_structures cannot be empty")

        if not self.target_stream or len(self.target_stream) == 0:
            raise ValueError("target_structures cannot be empty")

        if self.alignment_method not in ["align", "super", "cealign"]:
            raise ValueError(f"Alignment method must be 'align', 'super', or 'cealign', got: {self.alignment_method}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input structures."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        selection_display = self.selection_spec if self.selection_spec else "All CA atoms (whole structure)"

        config_lines.extend([
            f"REFERENCE STRUCTURES: {len(self.reference_stream)} files",
            f"TARGET STRUCTURES: {len(self.target_stream)} files",
            f"SELECTION: {selection_display}",
            f"ALIGNMENT METHOD: {self.alignment_method}",
            f"METRICS: RMSD (from PyMOL {self.alignment_method})"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate conformational change analysis execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        script_content = "#!/bin/bash\n"
        script_content += "# ConformationalChange execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self.generate_script_run_analysis()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_analysis(self) -> str:
        """Generate the conformational change analysis part of the script."""
        import json

        # Serialize DataStreams to JSON for HelpScript to load
        with open(self.reference_ds_json, 'w') as f:
            json.dump(self.reference_stream.to_dict(), f, indent=2)

        with open(self.target_ds_json, 'w') as f:
            json.dump(self.target_stream.to_dict(), f, indent=2)

        # Handle selection input
        if self.selection_spec is None:
            selection_config = {"type": "all"}
        else:
            selection_config = {"type": "fixed", "value": self.selection_spec}

        config_data = {
            "reference_structures_json": self.reference_ds_json,
            "target_structures_json": self.target_ds_json,
            "selection": selection_config,
            "alignment_method": self.alignment_method,
            "output_csv": self.analysis_csv
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        return f"""echo "Running conformational change analysis"
echo "Reference structures: {len(self.reference_stream)}"
echo "Target structures: {len(self.target_stream)}"
echo "Selection: {self.selection_spec}"
echo "Alignment method: {self.alignment_method}"
echo "Output: {self.analysis_csv}"

python "{self.analysis_py}" --config "{self.config_file}"

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after conformational change analysis."""
        tables = {
            "changes": TableInfo(
                name="changes",
                path=self.analysis_csv,
                columns=["id", "reference_structure", "target_structure", "selection",
                        "num_aligned_atoms", "RMSD"],
                description="Conformational change analysis between reference and target structures",
                count=len(self.target_stream)
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
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                "selection": str(self.selection_spec),
                "alignment_method": self.alignment_method
            }
        })
        return base_dict
