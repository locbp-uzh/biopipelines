"""
SequenceMetricsCorrelation tool for analyzing correlations between sequence mutations and metrics.

Analyzes mutation-metric correlations across sequence pools to identify beneficial/detrimental
mutations for data-driven sequence optimization. Outputs statistical analysis and scored
mutation tables compatible with MutationComposer.
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


class SequenceMetricsCorrelation(BaseConfig):
    """
    Pipeline tool for analyzing correlations between mutations and metrics.

    Analyzes sequence pools with associated metrics to identify which mutations
    correlate with better/worse performance. Generates scored mutation tables
    for data-driven sequence optimization.

    Commonly used for:
    - Learning mutation-metric correlations in iterative design cycles
    - Identifying beneficial mutations for binding affinity optimization
    - Data-driven guidance for sequence composition
    - Position-specific fitness landscape analysis
    """

    # Tool identification
    TOOL_NAME = "SequenceMetricsCorrelation"
    DEFAULT_ENV = None  # Loaded from config.yaml

    def __init__(self,
                 sequences: Union[ToolOutput, StandardizedOutput],
                 metrics: Union[ToolOutput, StandardizedOutput, TableInfo, str],
                 reference_sequence: Union[str, ToolOutput, StandardizedOutput],
                 metric_columns: Union[str, List[str]],
                 primary_metric: str,
                 mode: str = "minimize",
                 min_observations: int = 3,
                 history: Optional[Union[ToolOutput, StandardizedOutput]] = None,
                 **kwargs):
        """
        Initialize sequence-metric correlation analysis.

        Args:
            sequences: Sequence pool with 'id' and 'sequence' columns
            metrics: Table with metric values (must have matching 'id' column)
            reference_sequence: Reference sequence for mutation calling
                               Can be:
                               - String: direct sequence
                               - ToolOutput/StandardizedOutput: extracts first sequence
            metric_columns: Column name(s) to analyze from metrics table
                           String for single metric or list for multiple
            primary_metric: Which metric to use for scoring (mutation_deltas, mutation_zscores)
                           Must be in metric_columns
            mode: Optimization direction for primary_metric
                  "minimize" → lower is better (e.g., affinity, RMSD)
                  "maximize" → higher is better (e.g., pLDDT, contacts)
            min_observations: Minimum count to assign non-zero scores in mutation_deltas/zscores
            history: Previous SequenceMetricAnalysis results to accumulate with
            **kwargs: Additional parameters

        Examples:
            # Single metric analysis (first cycle)
            analysis = SequenceMetricAnalysis(
                sequences=filtered.tables.sequences,
                metrics=filtered.tables.merged,
                reference_sequence=original_holo,
                metric_columns="affinity_pred_value",
                primary_metric="affinity_pred_value",
                mode="minimize"
            )

            # Multi-metric analysis
            analysis = SequenceMetricAnalysis(
                sequences=filtered.tables.sequences,
                metrics=filtered.tables.merged,
                reference_sequence=original_holo,
                metric_columns=["affinity_pred_value", "complex_plddt", "contacts"],
                primary_metric="affinity_pred_value",
                mode="minimize",
                min_observations=5
            )

            # Accumulate with previous cycle
            analysis = SequenceMetricAnalysis(
                sequences=filtered.tables.sequences,
                metrics=filtered.tables.merged,
                reference_sequence=original_holo,
                metric_columns=["affinity_pred_value", "complex_plddt"],
                primary_metric="affinity_pred_value",
                mode="minimize",
                history=previous_analysis
            )
        """
        self.sequences_input = sequences
        self.metrics_input = metrics
        self.reference_sequence_input = reference_sequence

        # Handle single metric or list of metrics
        if isinstance(metric_columns, str):
            self.metric_columns = [metric_columns]
        else:
            self.metric_columns = metric_columns

        self.primary_metric = primary_metric
        self.mode = mode
        self.min_observations = min_observations
        self.history_input = history

        # Validate mode
        if mode not in ["minimize", "maximize"]:
            raise ValueError(f"mode must be 'minimize' or 'maximize', got '{mode}'")

        # Validate primary_metric is in metric_columns
        if primary_metric not in self.metric_columns:
            raise ValueError(
                f"primary_metric '{primary_metric}' must be in metric_columns {self.metric_columns}"
            )

        # Initialize base class
        super().__init__(**kwargs)

        # Set up dependencies
        if hasattr(sequences, 'config'):
            self.dependencies.append(sequences.config)
        if hasattr(metrics, 'config'):
            self.dependencies.append(metrics.config)
        if hasattr(reference_sequence, 'config'):
            self.dependencies.append(reference_sequence.config)
        if history and hasattr(history, 'config'):
            self.dependencies.append(history.config)


    def validate_params(self):
        """Validate SequenceMetricAnalysis parameters."""
        if not isinstance(self.sequences_input, (ToolOutput, StandardizedOutput, TableInfo)):
            raise ValueError("sequences must be a ToolOutput, StandardizedOutput, or TableInfo object")

        if not isinstance(self.metrics_input, (ToolOutput, StandardizedOutput, TableInfo, str)):
            raise ValueError("metrics must be a ToolOutput, StandardizedOutput, TableInfo, or string path")

        if self.history_input is not None and not isinstance(self.history_input, (ToolOutput, StandardizedOutput)):
            raise ValueError("history must be a ToolOutput or StandardizedOutput object if provided")

        if self.min_observations < 1:
            raise ValueError(f"min_observations must be >= 1, got {self.min_observations}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input sequences and metrics from previous tools."""
        self.folders = pipeline_folders

        # Extract sequences path
        self.sequences_path = self._extract_sequences_path(self.sequences_input)

        # Extract metrics path
        self.metrics_path = self._extract_metrics_path(self.metrics_input)

        # Extract reference sequence
        self.reference_sequence = self._extract_reference_sequence(self.reference_sequence_input)

        # Extract history path if provided
        if self.history_input is not None:
            self.history_path = self._extract_history_path(self.history_input)
        else:
            self.history_path = None

    def _extract_sequences_path(self, input_obj: Union[ToolOutput, StandardizedOutput, TableInfo]) -> str:
        """Extract sequences table path from input."""
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
            f"METRIC COLUMNS: {', '.join(self.metric_columns)}",
            f"PRIMARY METRIC: {self.primary_metric}",
            f"MODE: {self.mode}",
            f"MIN OBSERVATIONS: {self.min_observations}",
            f"HISTORY: {'Yes' if self.history_input else 'No'}"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """
        Generate sequence-metric analysis execution script.

        Args:
            script_path: Path where script should be written

        Returns:
            Script content as string
        """
        runtime_folder = os.path.dirname(script_path)
        output_folder = self.output_folder
        os.makedirs(output_folder, exist_ok=True)

        # Output files
        mutation_stats_csv = os.path.join(output_folder, "mutation_statistics.csv")
        mutation_deltas_csv = os.path.join(output_folder, "mutation_deltas.csv")
        mutation_zscores_csv = os.path.join(output_folder, "mutation_zscores.csv")
        top_mutations_csv = os.path.join(output_folder, "top_mutations.csv")
        coverage_csv = os.path.join(output_folder, "coverage.csv")

        # Create config file
        config_file = os.path.join(output_folder, "analysis_config.json")
        config_data = {
            "sequences_path": self.sequences_path,
            "metrics_path": self.metrics_path,
            "reference_sequence": self.reference_sequence,
            "metric_columns": self.metric_columns,
            "primary_metric": self.primary_metric,
            "mode": self.mode,
            "min_observations": self.min_observations,
            "history_path": self.history_path,
            "mutation_statistics_output": mutation_stats_csv,
            "mutation_deltas_output": mutation_deltas_csv,
            "mutation_zscores_output": mutation_zscores_csv,
            "top_mutations_output": top_mutations_csv,
            "coverage_output": coverage_csv
        }

        import json
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        # Generate script content
        script_content = f"""#!/bin/bash
# SequenceMetricAnalysis execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Running sequence-metric correlation analysis"
echo "Sequences: {self.sequences_path}"
echo "Metrics: {self.metrics_path}"
echo "Primary metric: {self.primary_metric} ({self.mode})"
echo "Output folder: {output_folder}"

# Run Python analysis script
python "{os.path.join(self.folders['HelpScripts'], 'pipe_sequence_metrics_correlation.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Successfully analyzed sequence-metric correlations"
    echo "Mutation statistics: {mutation_stats_csv}"
    echo "Mutation deltas: {mutation_deltas_csv}"
    echo "Mutation z-scores: {mutation_zscores_csv}"
    echo "Top mutations: {top_mutations_csv}"
    echo "Coverage: {coverage_csv}"
else
    echo "Error: Failed to analyze sequence-metric correlations"
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
        mutation_stats_csv = os.path.join(self.output_folder, "mutation_statistics.csv")
        mutation_deltas_csv = os.path.join(self.output_folder, "mutation_deltas.csv")
        mutation_zscores_csv = os.path.join(self.output_folder, "mutation_zscores.csv")
        top_mutations_csv = os.path.join(self.output_folder, "top_mutations.csv")
        coverage_csv = os.path.join(self.output_folder, "coverage.csv")

        # Build dynamic columns for mutation_statistics based on metrics
        stats_columns = ["position", "wt_aa", "mut_aa", "count"]
        for metric in self.metric_columns:
            stats_columns.extend([
                f"{metric}_mean",
                f"{metric}_std",
                f"{metric}_min",
                f"{metric}_max"
            ])

        # Top mutations columns
        top_columns = ["position", "wt_aa", "best_mutation", "count", "delta_score", "zscore"]
        for metric in self.metric_columns:
            top_columns.append(f"{metric}_mean")

        # Define tables
        tables = {
            "mutation_statistics": TableInfo(
                name="mutation_statistics",
                path=mutation_stats_csv,
                columns=stats_columns,
                description=f"Detailed mutation-metric statistics for all observed mutations",
                count=None  # Determined at runtime
            ),
            "mutation_deltas": TableInfo(
                name="mutation_deltas",
                path=mutation_deltas_csv,
                columns=["position", "wt_aa", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"],
                description=f"Delta scores (improvement vs WT) for {self.primary_metric} - compatible with MutationComposer",
                count=None
            ),
            "mutation_zscores": TableInfo(
                name="mutation_zscores",
                path=mutation_zscores_csv,
                columns=["position", "wt_aa", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"],
                description=f"Z-scores (standardized) for {self.primary_metric} - compatible with MutationComposer",
                count=None
            ),
            "top_mutations": TableInfo(
                name="top_mutations",
                path=top_mutations_csv,
                columns=top_columns,
                description="Best mutation at each position",
                count=None
            ),
            "coverage": TableInfo(
                name="coverage",
                path=coverage_csv,
                columns=["position", "wt_aa", "n_observations", "n_mutations_tested", "coverage_fraction", "max_count", "min_count", "mean_count"],
                description="Position exploration statistics",
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
                "metric_columns": self.metric_columns,
                "primary_metric": self.primary_metric,
                "mode": self.mode,
                "min_observations": self.min_observations,
                "has_history": self.history_input is not None
            }
        })
        return base_dict
