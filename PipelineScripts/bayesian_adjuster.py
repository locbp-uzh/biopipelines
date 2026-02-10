# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

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
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path

# Standard amino acids - guaranteed output structure
AMINO_ACIDS = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    from datastream import DataStream


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

    # Lazy path descriptors
    adjusted_probs_csv = Path(lambda self: os.path.join(self.output_folder, "adjusted_probabilities.csv"))
    absolute_probs_csv = Path(lambda self: os.path.join(self.output_folder, "absolute_probabilities.csv"))
    relative_probs_csv = Path(lambda self: os.path.join(self.output_folder, "relative_probabilities.csv"))
    adjustment_log_csv = Path(lambda self: os.path.join(self.output_folder, "adjustment_log.csv"))
    adjusted_logo_svg = Path(lambda self: os.path.join(self.output_folder, "adjusted_probabilities_logo.svg"))
    adjusted_logo_png = Path(lambda self: os.path.join(self.output_folder, "adjusted_probabilities_logo.png"))
    absolute_logo_svg = Path(lambda self: os.path.join(self.output_folder, "absolute_probabilities_logo.svg"))
    absolute_logo_png = Path(lambda self: os.path.join(self.output_folder, "absolute_probabilities_logo.png"))
    relative_logo_svg = Path(lambda self: os.path.join(self.output_folder, "relative_probabilities_logo.svg"))
    relative_logo_png = Path(lambda self: os.path.join(self.output_folder, "relative_probabilities_logo.png"))
    config_file = Path(lambda self: os.path.join(self.output_folder, "adjustment_config.json"))
    adjuster_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_bayesian_adjuster.py"))

    def __init__(self,
                 frequencies: Union[TableInfo, str],
                 correlations: Union[TableInfo, str],
                 mode: str = "min",
                 gamma: float = 3.0,
                 kappa: float = 10.0,
                 pseudocount: float = 0.01,
                 positions: Optional[str] = None,
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
            pseudocount: Pseudocount added only to amino acids with good correlations (default: 0.01)
                        Only amino acids where c(i,aa) is nonzero and good (negative for min mode,
                        positive for max mode) receive pseudocounts. This allows correlation signals
                        to resurrect unseen but beneficial mutations while preventing all amino acids
                        from having nonzero probability. After adding, frequencies are normalized to
                        preserve the original sum at each position.
            positions: PyMOL-style selection string for positions to display in plots (e.g., "141+143+145+147-149")
                      If None, shows all positions with adjustments. This ensures consistent x-axis across tools.
            **kwargs: Additional parameters

        Examples:
            # Adjust frequencies based on correlation signals
            adjuster = BayesianAdjuster(
                frequencies=profiler.tables.absolute_frequencies,
                correlations=correlation_analysis.tables.correlation_2d,
                mode="min",  # Minimize affinity
                gamma=3.0,
                pseudocount=0.01
            )

            # With position filter for plot consistency
            adjuster = BayesianAdjuster(
                frequencies=profiler.tables.absolute_frequencies,
                correlations=correlation_analysis.tables.correlation_2d,
                mode="min",
                gamma=3.0,
                pseudocount=0.01,
                positions="141+143+145+147-149+151-152"
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
        self.pseudocount = pseudocount
        self.positions = positions

        super().__init__(**kwargs)

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
        #if self.gamma <= 0:
        #    raise ValueError(f"gamma must be positive, got: {self.gamma}")
        #if self.kappa < 0:
        #    raise ValueError(f"kappa must be non-negative, got: {self.kappa}")
        if self.pseudocount < 0:
            raise ValueError(f"pseudocount must be non-negative, got: {self.pseudocount}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input tables."""
        self.folders = pipeline_folders

        # Extract table paths
        self.frequencies_path = self.frequencies_input.info.path if isinstance(self.frequencies_input, TableInfo) else self.frequencies_input
        self.correlations_path = self.correlations_input.info.path if isinstance(self.correlations_input, TableInfo) else self.correlations_input

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"MODE: {self.mode}",
            f"GAMMA: {self.gamma}",
            f"KAPPA: {self.kappa}",
            f"PSEUDOCOUNT: {self.pseudocount}",
            f"POSITIONS: {self.positions if self.positions else 'auto (all adjustments)'}",
            f"FREQUENCIES: {self.frequencies_path if hasattr(self, 'frequencies_path') else 'not configured'}",
            f"CORRELATIONS: {self.correlations_path if hasattr(self, 'correlations_path') else 'not configured'}"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate Bayesian adjustment execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        script_content = "#!/bin/bash\n"
        script_content += "# BayesianAdjuster execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self.generate_script_run_adjustment()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_adjustment(self) -> str:
        """Generate the Bayesian adjustment part of the script."""
        import json

        config_data = {
            "frequencies_path": self.frequencies_path,
            "correlations_path": self.correlations_path,
            "mode": self.mode,
            "gamma": self.gamma,
            "kappa": self.kappa,
            "pseudocount": self.pseudocount,
            "positions": self.positions,
            "adjusted_probabilities_output": self.adjusted_probs_csv,
            "absolute_probabilities_output": self.absolute_probs_csv,
            "relative_probabilities_output": self.relative_probs_csv,
            "adjustment_log_output": self.adjustment_log_csv,
            "adjusted_logo_svg": self.adjusted_logo_svg,
            "adjusted_logo_png": self.adjusted_logo_png,
            "absolute_logo_svg": self.absolute_logo_svg,
            "absolute_logo_png": self.absolute_logo_png,
            "relative_logo_svg": self.relative_logo_svg,
            "relative_logo_png": self.relative_logo_png
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        return f"""echo "Running Bayesian frequency adjustment"
echo "Mode: {self.mode}"
echo "Gamma: {self.gamma}"
echo "Kappa: {self.kappa}"
echo "Pseudocount: {self.pseudocount}"
echo "Output folder: {self.output_folder}"

python "{self.adjuster_py}" --config "{self.config_file}"

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after adjustment."""
        aa_columns = ["position", "original"] + AMINO_ACIDS

        tables = {
            "adjusted_probabilities": TableInfo(
                name="adjusted_probabilities",
                path=self.adjusted_probs_csv,
                columns=aa_columns,
                description="Raw Bayesian-adjusted probabilities before normalization",
                count=0
            ),
            "absolute_probabilities": TableInfo(
                name="absolute_probabilities",
                path=self.absolute_probs_csv,
                columns=aa_columns,
                description="Normalized absolute probabilities (comparable to MutationProfiler absolute_frequencies)",
                count=0
            ),
            "relative_probabilities": TableInfo(
                name="relative_probabilities",
                path=self.relative_probs_csv,
                columns=aa_columns,
                description="Normalized relative probabilities (comparable to MutationProfiler relative_frequencies)",
                count=0
            ),
            "adjustment_log": TableInfo(
                name="adjustment_log",
                path=self.adjustment_log_csv,
                columns=["position", "wt_aa", "aa", "prior_freq", "correlation", "adjusted_prob", "change"],
                description="Log of Bayesian adjustments for debugging",
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
                "mode": self.mode,
                "gamma": self.gamma,
                "kappa": self.kappa,
                "pseudocount": self.pseudocount,
                "positions": self.positions
            }
        })
        return base_dict
