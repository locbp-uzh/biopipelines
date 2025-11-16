"""
BayesianAdjuster tool for adjusting mutation frequencies based on correlation signals.

Applies Bayesian log-odds updates to mutation frequency tables using correlation signals
from SequenceMetricCorrelation. Implements the feedback loop for iterative sequence
optimization by adjusting probabilities to favor mutations that correlate with metric improvement.

Mathematical formula:
    p(i,aa|c) = σ(σ⁻¹(p₀(i,aa)) + γ·c(i,aa))

where:
    - p₀(i,aa) = prior probability (from MutationProfiler)
    - c(i,aa) = correlation signal (from SequenceMetricCorrelation)
    - γ = strength hyperparameter
    - σ(x) = 1/(1+e⁻ˣ) (sigmoid function)
    - σ⁻¹(p) = log(p/(1-p)) (logit function)
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


class BayesianAdjuster(BaseConfig):
    """
    Pipeline tool for Bayesian adjustment of mutation frequencies.

    Takes frequency tables e.g. from MutationProfiler and correlation signals e.g. from
    SequenceMetricCorrelation, applies Bayesian log-odds updates to generate
    adjusted probabilities that favor beneficial mutations.

    Commonly used for:
    - Iterative sequence optimization with feedback
    - Adapting mutation sampling based on performance data
    - Creating probability distributions informed by correlation signals
    """

    # Tool identification
    TOOL_NAME = "BayesianAdjuster"
    DEFAULT_ENV = None  # Loaded from config.yaml

    def __init__(self,
                 frequencies: Union[TableInfo, str],
                 correlations: Union[TableInfo, str],
                 mode: str = "min",
                 gamma: float = 3.0,
                 kappa: float = 10.0,
                 **kwargs):
        """
        Initialize Bayesian frequency adjuster.

        Args:
            frequencies: Frequency table from MutationProfiler (typically absolute_frequencies)
                        Should have columns: position, original, A, C, D, ..., Y
            correlations: Correlation table from SequenceMetricCorrelation (correlation_2d)
                         Should have columns: position, wt_aa, A, C, D, ..., Y
            mode: Optimization direction - "min" or "max"
                  "min" = lower metric values are better (e.g., binding affinity)
                  "max" = higher metric values are better (e.g., activity)
            gamma: Strength hyperparameter for Bayesian update (default: 3.0)
                   Higher values = more aggressive adjustment based on correlations
                   Lower values = more conservative, stays closer to prior
            kappa: Pseudo-observations for sample size shrinkage (default: 10.0)
                   Used to down-weight correlations from small sample sizes
            **kwargs: Additional parameters

        Examples:
            # Adjust frequencies based on correlation signals
            adjuster = BayesianAdjuster(
                frequencies=profiler.tables.absolute_frequencies,
                correlations=correlation_analysis.tables.correlation_2d,
                mode="min",  # Minimize affinity
                gamma=3.0
            )

            # Use adjusted probabilities in MutationComposer
            composer = MutationComposer(
                frequencies=adjuster.tables.absolute_probabilities,
                num_sequences=50,
                mode="weighted_random"
            )
        """
        self.frequencies_input = frequencies
        self.correlations_input = correlations
        self.mode = mode
        self.gamma = gamma
        self.kappa = kappa

        # Initialize base class
        super().__init__(**kwargs)

        # Set up dependencies
        if hasattr(frequencies, 'config'):
            self.dependencies.append(frequencies.config)
        if hasattr(correlations, 'config'):
            self.dependencies.append(correlations.config)

    def validate_params(self):
        """Validate BayesianAdjuster parameters."""
        # Validate frequencies input
        if not isinstance(self.frequencies_input, (TableInfo, str)):
            raise ValueError("frequencies must be a TableInfo object or string path")

        # Validate correlations input
        if not isinstance(self.correlations_input, (TableInfo, str)):
            raise ValueError("correlations must be a TableInfo object or string path")

        # Validate mode
        if self.mode not in ["min", "max"]:
            raise ValueError(f"mode must be 'min' or 'max', got: {self.mode}")

        # Validate hyperparameters
        if self.gamma <= 0:
            raise ValueError(f"gamma must be positive, got: {self.gamma}")
        if self.kappa < 0:
            raise ValueError(f"kappa must be non-negative, got: {self.kappa}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input tables from previous tools."""
        self.folders = pipeline_folders

        # Extract frequencies table path
        self.frequencies_path = self._extract_table_path(self.frequencies_input, "frequencies")

        # Extract correlations table path
        self.correlations_path = self._extract_table_path(self.correlations_input, "correlations")

    def _extract_table_path(self, input_obj: Union[TableInfo, str], name: str) -> str:
        """Extract table path from input."""
        # Direct path string
        if isinstance(input_obj, str):
            return input_obj

        # TableInfo object
        if isinstance(input_obj, TableInfo):
            return input_obj.path

        raise ValueError(f"Could not extract {name} table path from input")

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"MODE: {self.mode}",
            f"GAMMA: {self.gamma}",
            f"KAPPA: {self.kappa}",
            f"FREQUENCIES: {self.frequencies_path if hasattr(self, 'frequencies_path') else 'not configured'}",
            f"CORRELATIONS: {self.correlations_path if hasattr(self, 'correlations_path') else 'not configured'}"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """
        Generate Bayesian adjustment execution script.

        Args:
            script_path: Path where script should be written

        Returns:
            Script content as string
        """
        runtime_folder = os.path.dirname(script_path)
        output_folder = self.output_folder
        os.makedirs(output_folder, exist_ok=True)

        # Output files
        adjusted_probs_csv = os.path.join(output_folder, "adjusted_probabilities.csv")
        absolute_probs_csv = os.path.join(output_folder, "absolute_probabilities.csv")
        relative_probs_csv = os.path.join(output_folder, "relative_probabilities.csv")
        adjustment_log_csv = os.path.join(output_folder, "adjustment_log.csv")

        # Logo plots
        adjusted_logo_svg = os.path.join(output_folder, "adjusted_probabilities_logo.svg")
        adjusted_logo_png = os.path.join(output_folder, "adjusted_probabilities_logo.png")
        absolute_logo_svg = os.path.join(output_folder, "absolute_probabilities_logo.svg")
        absolute_logo_png = os.path.join(output_folder, "absolute_probabilities_logo.png")
        relative_logo_svg = os.path.join(output_folder, "relative_probabilities_logo.svg")
        relative_logo_png = os.path.join(output_folder, "relative_probabilities_logo.png")

        # Create config file
        config_file = os.path.join(output_folder, "adjustment_config.json")
        config_data = {
            "frequencies_path": self.frequencies_path,
            "correlations_path": self.correlations_path,
            "mode": self.mode,
            "gamma": self.gamma,
            "kappa": self.kappa,
            "adjusted_probabilities_output": adjusted_probs_csv,
            "absolute_probabilities_output": absolute_probs_csv,
            "relative_probabilities_output": relative_probs_csv,
            "adjustment_log_output": adjustment_log_csv,
            "adjusted_logo_svg": adjusted_logo_svg,
            "adjusted_logo_png": adjusted_logo_png,
            "absolute_logo_svg": absolute_logo_svg,
            "absolute_logo_png": absolute_logo_png,
            "relative_logo_svg": relative_logo_svg,
            "relative_logo_png": relative_logo_png
        }

        import json
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        # Generate script content
        script_content = f"""#!/bin/bash
# BayesianAdjuster execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Running Bayesian frequency adjustment"
echo "Mode: {self.mode}"
echo "Gamma: {self.gamma}"
echo "Kappa: {self.kappa}"
echo "Output folder: {output_folder}"

# Run Python adjustment script
python "{os.path.join(self.folders['HelpScripts'], 'pipe_bayesian_adjuster.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Successfully adjusted frequencies"
    echo "Adjusted probabilities: {adjusted_probs_csv}"
    echo "Absolute probabilities: {absolute_probs_csv}"
    echo "Relative probabilities: {relative_probs_csv}"
    echo "Logo plots generated"
else
    echo "Error: Failed to adjust frequencies"
    exit 1
fi

{self.generate_completion_check_footer()}
"""

        return script_content

    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after adjustment.

        Returns:
            Dictionary with output file paths and table information
        """
        adjusted_probs_csv = os.path.join(self.output_folder, "adjusted_probabilities.csv")
        absolute_probs_csv = os.path.join(self.output_folder, "absolute_probabilities.csv")
        relative_probs_csv = os.path.join(self.output_folder, "relative_probabilities.csv")
        adjustment_log_csv = os.path.join(self.output_folder, "adjustment_log.csv")

        # Define tables
        tables = {
            "adjusted_probabilities": TableInfo(
                name="adjusted_probabilities",
                path=adjusted_probs_csv,
                columns=["position", "original", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"],
                description="Raw Bayesian-adjusted probabilities before normalization",
                count=None
            ),
            "absolute_probabilities": TableInfo(
                name="absolute_probabilities",
                path=absolute_probs_csv,
                columns=["position", "original", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"],
                description="Normalized absolute probabilities (comparable to MutationProfiler absolute_frequencies)",
                count=None
            ),
            "relative_probabilities": TableInfo(
                name="relative_probabilities",
                path=relative_probs_csv,
                columns=["position", "original", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"],
                description="Normalized relative probabilities (comparable to MutationProfiler relative_frequencies)",
                count=None
            ),
            "adjustment_log": TableInfo(
                name="adjustment_log",
                path=adjustment_log_csv,
                columns=["position", "wt_aa", "aa", "prior_freq", "correlation", "adjusted_prob", "change"],
                description="Log of Bayesian adjustments for debugging",
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
                "mode": self.mode,
                "gamma": self.gamma,
                "kappa": self.kappa
            }
        })
        return base_dict
