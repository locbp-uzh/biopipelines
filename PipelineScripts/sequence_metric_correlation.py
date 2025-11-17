"""
SequenceMetricCorrelation tool for analyzing correlations between sequence mutations and metrics.

Computes correlation signals c(i) and c(i,aa) where:
c(i) = (mean_metric_mutated - mean_metric_wt) / sqrt(var_mutated + var_wt)
This measures the statistical correlation between mutations at position i and metric performance.
"""

import os
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo


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
    DEFAULT_ENV = None  # Loaded from config.yaml

    def __init__(self,
                 mutants: Union[ToolOutput, StandardizedOutput, List[Union[ToolOutput, StandardizedOutput, TableInfo, str]]],
                 data: Union[ToolOutput, StandardizedOutput, TableInfo, str, List[Union[ToolOutput, StandardizedOutput, TableInfo, str]]],
                 original: Union[str, ToolOutput, StandardizedOutput],
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
                      Can be:
                      - String: direct sequence
                      - ToolOutput/StandardizedOutput: extracts first sequence
            metric: Column name of the metric to analyze
            positions: PyMOL-style selection string for positions to display in plots (e.g., "141+143+145+147-149")
                      If None, shows all positions with correlations. This ensures consistent x-axis across tools.
            **kwargs: Additional parameters

        Examples:
            # Single cycle analysis
            correlation = SequenceMetricCorrelation(
                mutants=filtered.tables.sequences,
                data=filtered.tables.merged,
                original=original_holo,
                metric="affinity_pred_value"
            )

            # Multi-cycle analysis (accumulate data) with position filter
            correlation = SequenceMetricCorrelation(
                mutants=[cycle1.tables.sequences, cycle2.tables.sequences],
                data=[cycle1.tables.merged, cycle2.tables.merged],
                original=original_holo,
                metric="affinity_pred_value",
                positions="141+143+145+147-149+151-152"
            )
        """
        # Handle list inputs
        if isinstance(mutants, list):
            self.mutants_input = mutants
        else:
            self.mutants_input = [mutants]

        if isinstance(data, list):
            self.data_input = data
        else:
            self.data_input = [data]

        self.original_input = original
        self.metric = metric
        self.positions = positions

        # Initialize base class
        super().__init__(**kwargs)

        # Set up dependencies
        for mutant_input in self.mutants_input:
            if hasattr(mutant_input, 'config'):
                self.dependencies.append(mutant_input.config)

        for data_input in self.data_input:
            if hasattr(data_input, 'config'):
                self.dependencies.append(data_input.config)

        if hasattr(original, 'config'):
            self.dependencies.append(original.config)


    def validate_params(self):
        """Validate SequenceMetricCorrelation parameters."""
        # Validate mutants inputs
        for mutant_input in self.mutants_input:
            if not isinstance(mutant_input, (ToolOutput, StandardizedOutput, TableInfo, str)):
                raise ValueError("Each mutants input must be a ToolOutput, StandardizedOutput, TableInfo, or string path")

        # Validate data inputs
        for data_input in self.data_input:
            if not isinstance(data_input, (ToolOutput, StandardizedOutput, TableInfo, str)):
                raise ValueError("Each data input must be a ToolOutput, StandardizedOutput, TableInfo, or string path")

        # Validate same number of mutants and data inputs
        if len(self.mutants_input) != len(self.data_input):
            raise ValueError(f"Number of mutants ({len(self.mutants_input)}) must match number of data inputs ({len(self.data_input)})")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input sequences and metrics from previous tools."""
        self.folders = pipeline_folders

        # Extract paths for all mutants tables
        self.mutants_paths = []
        for mutant_input in self.mutants_input:
            path = self._extract_sequences_path(mutant_input)
            self.mutants_paths.append(path)

        # Extract paths for all data tables
        self.data_paths = []
        for data_input in self.data_input:
            path = self._extract_metrics_path(data_input)
            self.data_paths.append(path)

        # Extract reference sequence
        self.original_sequence = self._extract_reference_sequence(self.original_input)

    def _extract_sequences_path(self, input_obj: Union[ToolOutput, StandardizedOutput, TableInfo, str]) -> str:
        """Extract sequences table path from input."""
        # Handle direct string path
        if isinstance(input_obj, str):
            return input_obj

        # Handle TableInfo object directly
        if isinstance(input_obj, TableInfo):
            if hasattr(input_obj, 'path'):
                return input_obj.path
            raise ValueError("TableInfo object does not have a 'path' attribute")

        # Handle ToolOutput/StandardizedOutput
        if hasattr(input_obj, 'tables'):
            tables = input_obj.tables

            # Check for sequences table
            if hasattr(tables, 'sequences'):
                if isinstance(tables.sequences, str):
                    return tables.sequences
                elif hasattr(tables.sequences, 'path'):
                    return tables.sequences.path

            # Check for merged table (from Filter output)
            if hasattr(tables, 'merged'):
                if isinstance(tables.merged, str):
                    return tables.merged
                elif hasattr(tables.merged, 'path'):
                    return tables.merged.path

            # Try _tables dict
            if hasattr(tables, '_tables'):
                for name, info in tables._tables.items():
                    if 'sequence' in name.lower():
                        return info.path

        # Fallback: predict in output folder
        if hasattr(input_obj, 'output_folder'):
            return os.path.join(input_obj.output_folder, 'sequences.csv')

        raise ValueError("Could not extract sequences path from input")

    def _extract_metrics_path(self, input_obj: Union[ToolOutput, StandardizedOutput, TableInfo, str]) -> str:
        """Extract metrics table path from input."""
        # Direct path string
        if isinstance(input_obj, str):
            return input_obj

        # TableInfo object
        if isinstance(input_obj, TableInfo):
            return input_obj.path

        # ToolOutput or StandardizedOutput
        if hasattr(input_obj, 'tables'):
            tables = input_obj.tables

            # Check for common table names
            for table_name in ['merged', 'affinity', 'confidence', 'analysis', 'results']:
                if hasattr(tables, table_name):
                    table = getattr(tables, table_name)
                    if isinstance(table, str):
                        return table
                    elif hasattr(table, 'path'):
                        return table.path

            # Try _tables dict
            if hasattr(tables, '_tables'):
                # Return first table
                for name, info in tables._tables.items():
                    return info.path

        # Fallback: predict in output folder
        if hasattr(input_obj, 'output_folder'):
            return os.path.join(input_obj.output_folder, 'metrics.csv')

        raise ValueError("Could not extract metrics path from input")

    def _extract_reference_sequence(self, ref_input: Union[str, ToolOutput, StandardizedOutput]) -> str:
        """Extract reference sequence string from input."""
        # Direct string
        if isinstance(ref_input, str):
            # Check if it's a sequence (contains only amino acids) or a path
            if all(c in 'ACDEFGHIKLMNPQRSTVWY' for c in ref_input.upper()):
                return ref_input
            # Otherwise assume it's a path - will be read at runtime
            return f"@{ref_input}"  # @ prefix indicates path

        # ToolOutput or StandardizedOutput - extract sequences table
        if hasattr(ref_input, 'tables') and hasattr(ref_input.tables, 'sequences'):
            seq_table = ref_input.tables.sequences
            if isinstance(seq_table, str):
                return f"@{seq_table}"  # Path to sequences file
            elif hasattr(seq_table, 'path'):
                return f"@{seq_table.path}"

        raise ValueError(
            "reference_sequence must be a sequence string or ToolOutput with sequences table"
        )

    def _extract_history_path(self, history_obj: Union[ToolOutput, StandardizedOutput]) -> str:
        """Extract mutation_statistics table path from history."""
        if hasattr(history_obj, 'tables'):
            tables = history_obj.tables

            # Look for mutation_statistics table
            if hasattr(tables, 'mutation_statistics'):
                table = tables.mutation_statistics
                if isinstance(table, str):
                    return table
                elif hasattr(table, 'path'):
                    return table.path

        raise ValueError("Could not extract mutation_statistics table from history")

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
        """
        Generate sequence-metric correlation analysis execution script.

        Args:
            script_path: Path where script should be written

        Returns:
            Script content as string
        """
        runtime_folder = os.path.dirname(script_path)
        output_folder = self.output_folder
        os.makedirs(output_folder, exist_ok=True)

        # Output files
        correlation_1d_csv = os.path.join(output_folder, "correlation_1d.csv")
        correlation_2d_csv = os.path.join(output_folder, "correlation_2d.csv")
        logo_svg = os.path.join(output_folder, "correlation_logo.svg")
        logo_png = os.path.join(output_folder, "correlation_logo.png")

        # Create config file
        config_file = os.path.join(output_folder, "analysis_config.json")
        config_data = {
            "mutants_paths": self.mutants_paths,
            "data_paths": self.data_paths,
            "original_sequence": self.original_sequence,
            "metric": self.metric,
            "positions": self.positions,
            "correlation_1d_output": correlation_1d_csv,
            "correlation_2d_output": correlation_2d_csv,
            "logo_svg_output": logo_svg,
            "logo_png_output": logo_png
        }

        import json
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        # Generate script content
        script_content = f"""#!/bin/bash
# SequenceMetricCorrelation execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Running sequence-metric correlation analysis"
echo "Metric: {self.metric}"
echo "Output folder: {output_folder}"

# Run Python analysis script
python "{os.path.join(self.folders['HelpScripts'], 'pipe_sequence_metric_correlation.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Successfully computed correlation signals"
    echo "1D correlations: {correlation_1d_csv}"
    echo "2D correlations: {correlation_2d_csv}"
    echo "Logo plot: {logo_svg}"
else
    echo "Error: Failed to compute correlation signals"
    exit 1
fi

{self.generate_completion_check_footer()}
"""

        return script_content

    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after analysis.

        Returns:
            Dictionary with output file paths and table information
        """
        correlation_1d_csv = os.path.join(self.output_folder, "correlation_1d.csv")
        correlation_2d_csv = os.path.join(self.output_folder, "correlation_2d.csv")

        # Define tables
        tables = {
            "correlation_1d": TableInfo(
                name="correlation_1d",
                path=correlation_1d_csv,
                columns=["position", "wt_aa", "correlation", "mean_mutated", "mean_wt", "var_mutated", "var_wt", "n_mutated", "n_wt"],
                description=f"1D correlation signal c(i) for position i",
                count=None  # Determined at runtime
            ),
            "correlation_2d": TableInfo(
                name="correlation_2d",
                path=correlation_2d_csv,
                columns=["position", "wt_aa", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"],
                description=f"2D correlation signal c(i,aa) for each position and amino acid",
                count=None
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
                "metric": self.metric,
                "num_mutants_inputs": len(self.mutants_input),
                "num_data_inputs": len(self.data_input),
                "positions": self.positions
            }
        })
        return base_dict
