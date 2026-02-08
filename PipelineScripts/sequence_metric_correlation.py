# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# This software is freely available for use, modification, and redistribution.
# If you use this software or any derivative of it in published work,
# you must cite the original author and this repository.
"""
SequenceMetricCorrelation tool for analyzing correlations between sequence mutations and metrics.

Computes correlation signals c(i) and c(i,aa) where:
c(i) = (mean_metric_mutated - mean_metric_wt) / sqrt(var_mutated + var_wt)
This measures the statistical correlation between mutations at position i and metric performance.
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

# Standard amino acids - guaranteed output structure
AMINO_ACIDS = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]


class SequenceMetricCorrelation(BaseConfig):
    """
    Pipeline tool for analyzing correlations between mutations and a metric.

    Computes correlation signals c(i) and c(i,aa) where:
    c(i) = (mean_metric_mutated - mean_metric_wt) / sqrt(var_mutated + var_wt)

    The formula quantifies how mutations at position i correlate with metric changes.
    Positive values indicate mutations improve the metric, negative values indicate
    mutations worsen the metric.

    Commonly used for:
    - Identifying positions where mutations correlate with metric improvements
    - Understanding mutation-metric relationships in iterative design
    - Generating correlation-based sequence logos
    """

    # Tool identification
    TOOL_NAME = "SequenceMetricCorrelation"

    # Lazy path descriptors
    correlation_1d_csv = Path(lambda self: os.path.join(self.output_folder, "correlation_1d.csv"))
    correlation_2d_csv = Path(lambda self: os.path.join(self.output_folder, "correlation_2d.csv"))
    logo_svg = Path(lambda self: os.path.join(self.output_folder, "correlation_logo.svg"))
    logo_png = Path(lambda self: os.path.join(self.output_folder, "correlation_logo.png"))
    config_file = Path(lambda self: os.path.join(self.output_folder, "analysis_config.json"))
    correlation_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_sequence_metric_correlation.py"))

    def __init__(self,
                 mutants: Union[DataStream, StandardizedOutput, TableInfo, str,
                               List[Union[DataStream, StandardizedOutput, TableInfo, str]]],
                 data: Union[TableInfo, str, List[Union[TableInfo, str]]],
                 original: Union[str, DataStream, StandardizedOutput],
                 metric: str,
                 positions: Optional[str] = None,
                 **kwargs):
        """
        Initialize sequence-metric correlation analysis.

        Args:
            mutants: Mutant sequences with 'id' and 'sequence' columns
                     Can be single table or list of tables (for multiple cycles)
            data: Table(s) with metric values (must have matching 'id' column and metric column)
                  Can be single table or list of tables (for multiple cycles)
            original: Original/reference sequence for mutation calling
                      Can be string (sequence) or DataStream/StandardizedOutput
            metric: Column name of the metric to analyze
            positions: PyMOL-style selection string for positions to display in plots
            **kwargs: Additional parameters

        Examples:
            # Single cycle analysis
            correlation = SequenceMetricCorrelation(
                mutants=filtered.tables.sequences,
                data=filtered.tables.merged,
                original=original_holo,
                metric="affinity_pred_value"
            )
        """
        # Handle list inputs
        self.mutants_input = mutants if isinstance(mutants, list) else [mutants]
        self.data_input = data if isinstance(data, list) else [data]
        self.original_input = original
        self.metric = metric
        self.positions = positions

        super().__init__(**kwargs)


    def validate_params(self):
        """Validate SequenceMetricCorrelation parameters."""
        # Validate mutants inputs
        for mutant_input in self.mutants_input:
            if not isinstance(mutant_input, (DataStream, StandardizedOutput, TableInfo, str)):
                raise ValueError("Each mutants input must be DataStream, StandardizedOutput, TableInfo, or string path")

        # Validate data inputs
        for data_input in self.data_input:
            if not isinstance(data_input, (TableInfo, str)):
                raise ValueError("Each data input must be TableInfo or string path")

        # Validate same number of mutants and data inputs
        if len(self.mutants_input) != len(self.data_input):
            raise ValueError(f"Number of mutants ({len(self.mutants_input)}) must match number of data inputs ({len(self.data_input)})")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input sequences and metrics."""
        self.folders = pipeline_folders

        # Extract paths for all mutants tables
        self.mutants_paths = []
        for mutant_input in self.mutants_input:
            path = self._extract_path(mutant_input, "sequences")
            self.mutants_paths.append(path)

        # Extract paths for all data tables
        self.data_paths = []
        for data_input in self.data_input:
            path = self._extract_path(data_input, "data")
            self.data_paths.append(path)

        # Extract reference sequence
        self.original_sequence = self._extract_reference_sequence(self.original_input)

    def _extract_path(self, input_obj: Union[DataStream, StandardizedOutput, TableInfo, str], input_type: str) -> str:
        """Extract table path from various input types."""
        if isinstance(input_obj, str):
            return input_obj
        elif isinstance(input_obj, TableInfo):
            return input_obj.path
        elif isinstance(input_obj, DataStream):
            return input_obj.map_table
        elif isinstance(input_obj, StandardizedOutput):
            if input_type == "sequences" and hasattr(input_obj.tables, 'sequences'):
                return input_obj.tables.sequences.path
            raise ValueError(f"Could not extract {input_type} path from StandardizedOutput")
        else:
            raise ValueError(f"Unsupported input type: {type(input_obj)}")

    def _extract_reference_sequence(self, ref_input: Union[str, DataStream, StandardizedOutput]) -> str:
        """Extract reference sequence string from input."""
        if isinstance(ref_input, str):
            # Check if it's a sequence (contains only amino acids) or a path
            if all(c in 'ACDEFGHIKLMNPQRSTVWY' for c in ref_input.upper()):
                return ref_input
            return f"@{ref_input}"  # @ prefix indicates path

        if isinstance(ref_input, DataStream):
            return f"@{ref_input.map_table}"

        if isinstance(ref_input, StandardizedOutput):
            if hasattr(ref_input.tables, 'sequences'):
                return f"@{ref_input.tables.sequences.path}"

        raise ValueError("original must be a sequence string, DataStream, or StandardizedOutput")

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"METRIC: {self.metric}",
            f"MUTANTS INPUTS: {len(self.mutants_input)}",
            f"DATA INPUTS: {len(self.data_input)}",
            f"POSITIONS: {self.positions if self.positions else 'auto (all correlations)'}"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate SequenceMetricCorrelation execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        script_content = "#!/bin/bash\n"
        script_content += "# SequenceMetricCorrelation execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self.generate_script_run_correlation_analysis()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_correlation_analysis(self) -> str:
        """Generate the correlation analysis execution part of the script."""
        import json

        config_data = {
            "mutants_paths": self.mutants_paths,
            "data_paths": self.data_paths,
            "original_sequence": self.original_sequence,
            "metric": self.metric,
            "positions": self.positions,
            "correlation_1d_output": self.correlation_1d_csv,
            "correlation_2d_output": self.correlation_2d_csv,
            "logo_svg_output": self.logo_svg,
            "logo_png_output": self.logo_png
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        return f"""echo "Running sequence-metric correlation analysis"
echo "Metric: {self.metric}"
echo "Output folder: {self.output_folder}"

python "{self.correlation_py}" --config "{self.config_file}"

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after analysis."""
        aa_columns = ["position", "wt_aa"] + AMINO_ACIDS

        tables = {
            "correlation_1d": TableInfo(
                name="correlation_1d",
                path=self.correlation_1d_csv,
                columns=["position", "wt_aa", "correlation", "mean_mutated", "mean_wt", "var_mutated", "var_wt", "n_mutated", "n_wt"],
                description="1D correlation signal c(i) for position i",
                count=0
            ),
            "correlation_2d": TableInfo(
                name="correlation_2d",
                path=self.correlation_2d_csv,
                columns=aa_columns,
                description="2D correlation signal c(i,aa) for each position and amino acid",
                count=0
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
                "metric": self.metric,
                "num_mutants_inputs": len(self.mutants_input),
                "num_data_inputs": len(self.data_input),
                "positions": self.positions
            }
        })
        return base_dict
