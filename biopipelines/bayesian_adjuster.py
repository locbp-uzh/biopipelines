# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
BayesianAdjuster tool for adjusting mutation frequencies based on correlation signals.

Applies Bayesian log-odds updates to amino-acid frequency tables using correlation
signals. Implements the feedback loop for iterative sequence optimization by
adjusting probabilities to favor mutations that correlate with metric improvement.

Mathematical formula:
    p(i,aa|c) = σ(σ⁻¹(p₀(i,aa)) + γ·c̃(i,aa))
    c̃(i,aa) = c(i,aa)·n(i,aa)/(n(i,aa)+κ)   (sample-size shrinkage; κ=0 ⟹ c̃=c)

where:
    - p₀(i,aa) = prior probability from the input frequency table
    - c(i,aa) = signal from the input correlation table
    - n(i,aa) = sample count for (i,aa) from the input sample-counts table
      (e.g. MutationProfiler.tables.mutations); only used when κ is provided
    - γ = strength hyperparameter
    - κ = shrinkage pseudo-observations; down-weights correlations with low
      sample support. κ=None/0 leaves the correlation unchanged.
    - σ(x) = 1/(1+e⁻ˣ) (sigmoid function)
    - σ⁻¹(p) = log(p/(1-p)) (logit function)
"""

import os
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream

# Standard amino acids - guaranteed output structure
AMINO_ACIDS = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]


class BayesianAdjuster(BaseConfig):
    """
    Pipeline tool for Bayesian adjustment of mutation frequencies.

    Takes amino-acid frequency and correlation/evidence tables, applies Bayesian
    log-odds updates, and generates adjusted probabilities that favor beneficial
    mutations.

    Commonly used for:
    - Iterative sequence optimization with feedback
    - Adapting mutation sampling based on performance data
    - Creating probability distributions informed by correlation signals
    """

    # Tool identification
    TOOL_NAME = "BayesianAdjuster"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== BayesianAdjuster ==="
echo "Requires MutationEnv (installed with MutationProfiler.install())"
echo "No additional installation needed."
touch "$INSTALL_SUCCESS"
echo "=== BayesianAdjuster ready ==="
"""

    # Lazy path descriptors
    adjusted_probs_csv = Path(lambda self: self.table_path("adjusted_probabilities"))
    absolute_probs_csv = Path(lambda self: self.table_path("absolute_probabilities"))
    relative_probs_csv = Path(lambda self: self.table_path("relative_probabilities"))
    adjustment_log_csv = Path(lambda self: self.table_path("adjustment_log"))
    adjusted_logo_svg = Path(lambda self: self.stream_path("images", "adjusted_probabilities_logo.svg"))
    adjusted_logo_png = Path(lambda self: self.stream_path("images", "adjusted_probabilities_logo.png"))
    absolute_logo_svg = Path(lambda self: self.stream_path("images", "absolute_probabilities_logo.svg"))
    absolute_logo_png = Path(lambda self: self.stream_path("images", "absolute_probabilities_logo.png"))
    relative_logo_svg = Path(lambda self: self.stream_path("images", "relative_probabilities_logo.svg"))
    relative_logo_png = Path(lambda self: self.stream_path("images", "relative_probabilities_logo.png"))
    config_file = Path(lambda self: self.configuration_path("adjustment_config.json"))
    adjuster_py = Path(lambda self: self.pipe_script_path("pipe_bayesian_adjuster.py"))

    def __init__(self,
                 frequencies: Union[TableInfo, str],
                 correlations: Union[TableInfo, str],
                 mode: str = "min",
                 gamma: float = 3.0,
                 kappa: Optional[float] = None,
                 sample_counts: Optional[Union[TableInfo, str]] = None,
                 pseudocount: float = 0.01,
                 positions: Optional[str] = None,
                 **kwargs):
        """
        Initialize Bayesian frequency adjuster.

        Args:
            frequencies: Frequency table with columns: position, original, A, C, D, ..., Y
            correlations: Correlation/evidence table with columns: position, original, A, C, D, ..., Y
            sample_counts: Sample-count table with columns: position, original, A, C, D, ..., Y.
                           Required when kappa is provided.
            mode: Optimization direction - "min" or "max"
                  "min" = lower metric values are better (e.g., binding affinity)
                  "max" = higher metric values are better (e.g., activity)
            gamma: Strength hyperparameter for Bayesian update (default: 3.0)
                   Higher values = more aggressive adjustment based on correlations
                   Lower values = more conservative, stays closer to prior
            kappa: Pseudo-observations for sample size shrinkage. If provided,
                   sample_counts is required.
            pseudocount: Pseudocount added only to amino acids with good correlations (default: 0.01)
                        Only amino acids where c(i,aa) is nonzero and good (negative for min mode,
                        positive for max mode) receive pseudocounts. This allows correlation signals
                        to resurrect unseen but beneficial mutations while preventing all amino acids
                        from having nonzero probability. After adding, frequencies are normalized to
                        preserve the original sum at each position.
            positions: PyMOL-style selection string for positions to display in plots (e.g., "141+143+145+147-149")
                      If None, shows all positions with adjustments. This ensures consistent x-axis across tools.
            **kwargs: Additional parameters

        Output:
            Streams: (none)
            Tables:
                adjusted_probabilities: position | original | A | C | D | ... | Y
                absolute_probabilities: position | original | A | C | D | ... | Y
                relative_probabilities: position | original | A | C | D | ... | Y
                adjustment_log: position | original | aa | prior_freq | correlation | raw_correlation | sample_count | shrinkage_weight | adjusted_prob | change
        """
        self.frequencies_input = frequencies
        self.correlations_input = correlations
        self.sample_counts_input = sample_counts
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

        # Validate sample counts input
        if self.sample_counts_input is not None and not isinstance(self.sample_counts_input, (TableInfo, str)):
            raise ValueError("sample_counts must be a TableInfo object or string path")

        # Validate mode
        if self.mode not in ["min", "max"]:
            raise ValueError(f"mode must be 'min' or 'max', got: {self.mode}")

        # Validate hyperparameters
        if self.gamma < 0:
            raise ValueError(f"gamma must be non-negative, got: {self.gamma}")
        if self.kappa is not None and self.kappa < 0:
            raise ValueError(f"kappa must be non-negative, got: {self.kappa}")
        if self.kappa is not None and self.sample_counts_input is None:
            raise ValueError("sample_counts is required when kappa is provided")
        if self.pseudocount < 0:
            raise ValueError(f"pseudocount must be non-negative, got: {self.pseudocount}")

        # Shell-safety: the table inputs reach bash as paths in the config JSON
        # and in get_config_display() echo lines. Validate the string form (a
        # TableInfo path is framework-generated and trusted).
        if isinstance(self.frequencies_input, str):
            _validate_freeform_string("frequencies", self.frequencies_input)
        if isinstance(self.correlations_input, str):
            _validate_freeform_string("correlations", self.correlations_input)
        if isinstance(self.sample_counts_input, str):
            _validate_freeform_string("sample_counts", self.sample_counts_input)
        _validate_freeform_string("positions", self.positions)

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input tables."""
        self.folders = pipeline_folders

        # Extract table paths
        self.frequencies_path = self.frequencies_input.info.path if isinstance(self.frequencies_input, TableInfo) else self.frequencies_input
        self.correlations_path = self.correlations_input.info.path if isinstance(self.correlations_input, TableInfo) else self.correlations_input
        self.sample_counts_path = None
        if self.sample_counts_input is not None:
            self.sample_counts_path = self.sample_counts_input.info.path if isinstance(self.sample_counts_input, TableInfo) else self.sample_counts_input

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
            f"CORRELATIONS: {self.correlations_path if hasattr(self, 'correlations_path') else 'not configured'}",
            f"SAMPLE_COUNTS: {self.sample_counts_path if getattr(self, 'sample_counts_path', None) else 'not provided'}"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate Bayesian adjustment execution script."""
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
            "sample_counts_path": self.sample_counts_path,
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

        # Colab exports MPLBACKEND=module://matplotlib_inline.backend_inline,
        # which leaks into this subprocess; logomaker/matplotlib reject it as an
        # invalid backend at import. Force a headless backend (correct on HPC too).
        return f"""export MPLBACKEND=Agg
echo "Running Bayesian frequency adjustment"
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
                description="Raw Bayesian-adjusted probabilities before normalization"
            ),
            "absolute_probabilities": TableInfo(
                name="absolute_probabilities",
                path=self.absolute_probs_csv,
                columns=aa_columns,
                description="Normalized probabilities across all amino acids at each position"
            ),
            "relative_probabilities": TableInfo(
                name="relative_probabilities",
                path=self.relative_probs_csv,
                columns=aa_columns,
                description="Normalized probabilities across non-original amino acids at each position"
            ),
            "adjustment_log": TableInfo(
                name="adjustment_log",
                path=self.adjustment_log_csv,
                columns=[
                    "position", "original", "aa", "prior_freq", "correlation",
                    "raw_correlation", "sample_count", "shrinkage_weight",
                    "adjusted_prob", "change"
                ],
                description="Log of Bayesian adjustments for debugging"
            )
        }

        # Declaring the images stream is what makes the framework create the
        # images/ stream folder; without it the pipe script's plt.savefig into
        # images/ fails with FileNotFoundError (the dir never exists).
        images = DataStream(
            name="images",
            ids=["adjusted_probabilities_logo", "absolute_probabilities_logo", "relative_probabilities_logo"],
            files=[self.adjusted_logo_png, self.absolute_logo_png, self.relative_logo_png],
            format="png"
        )

        return {
            "images": images,
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
                "sample_counts": self.sample_counts_path if hasattr(self, "sample_counts_path") else None,
                "pseudocount": self.pseudocount,
                "positions": self.positions
            }
        })
        return base_dict
