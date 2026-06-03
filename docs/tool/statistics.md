# Sequence Statistics

[← Back to Tool Reference](../tool_reference.md)

---

### BayesianAdjuster

Adjusts mutation frequencies using correlation signals via Bayesian log-odds updates. Boosts beneficial mutations and suppresses detrimental ones based on correlation evidence.

**Bayesian Update Formula**:
```
p(i,aa|c) = σ(σ⁻¹(p₀(i,aa)) + γ·c(i,aa))
```
Where:
- p₀(i,aa) = prior probability from the input frequency table
- c(i,aa) = signal from the input correlation table
- γ = strength hyperparameter (higher = more aggressive adjustment)
- σ(x) = sigmoid function = 1/(1+e⁻ˣ)
- σ⁻¹(p) = logit function = log(p/(1-p))

**Installation**: Requires the standard mutation-analysis Python environment with pandas, NumPy, matplotlib, and logomaker.

**Parameters**:
- `frequencies`: Union[TableInfo, str] (required) - Frequency table with columns `position`, `original`, and one column per amino acid: `A`, `C`, `D`, ..., `Y`.
- `correlations`: Union[TableInfo, str] (required) - Correlation/evidence table with columns `position`, `original`, and one column per amino acid: `A`, `C`, `D`, ..., `Y`.
- `sample_counts`: Optional[Union[TableInfo, str]] = None - Sample-count table with columns `position`, `original`, and one column per amino acid: `A`, `C`, `D`, ..., `Y`. Required when `kappa` is provided.
- `mode`: str = "min" - Optimization direction:
  - "min" = lower metric values are better (e.g., binding affinity)
  - "max" = higher metric values are better (e.g., activity)
- `gamma`: float = 3.0 - Strength hyperparameter for Bayesian update (higher = more aggressive)
- `kappa`: Optional[float] = None - Pseudo-observations for sample size shrinkage. When provided, correlation signals are weighted by `n / (n + kappa)` using `sample_counts`.
- `pseudocount`: float = 0.01 - Pseudocount added to all amino acids before adjustment. Ensures no amino acid has zero probability, allowing correlation signals to resurrect beneficial mutations. After adding, frequencies are normalized to preserve original sum.
- `positions`: Optional[str] = None - PyMOL-style position filter (e.g., "141+143+145+147-149")

**Outputs**:
- `tables.adjusted_probabilities`: Raw Bayesian-adjusted probabilities
- `tables.absolute_probabilities`: Normalized as absolute probabilities
- `tables.relative_probabilities`: Normalized as relative probabilities (original AA = 0)
- `tables.adjustment_log`: Detailed log of all adjustments

**Example**:
```python
from biopipelines.bayesian_adjuster import BayesianAdjuster
from biopipelines.mutation_composer import MutationComposer

# Basic usage
adjuster = BayesianAdjuster(
    frequencies=frequency_table,
    correlations=correlation_table,
    mode="min",  # Minimize affinity
    gamma=3.0,
    pseudocount=0.01
)

# More aggressive adjustment with position filter
adjuster = BayesianAdjuster(
    frequencies=frequency_table,
    correlations=correlation_table,
    sample_counts=sample_counts_table,
    mode="min",
    gamma=5.0,  # More aggressive
    kappa=10.0,  # Shrink low-sample correlation signals
    pseudocount=0.005,  # Lower threshold for resurrected mutations
    positions="141+143+145+147-149"
)

# Use adjusted frequencies in MutationComposer
composer = MutationComposer(
    frequencies=adjuster.tables.absolute_probabilities,
    num_sequences=50,
    mode="weighted_random"
)
```

**Key Features**:
- **Pseudocount mechanism**: Allows resurrection of beneficial mutations that weren't initially present (probability = 0) but have strong positive correlations
- **Position filtering**: Restricts sequence-logo plots to selected positions for easier visual comparison
- **Dual normalization**: Provides both all-amino-acid normalization and relative probabilities with the original amino acid excluded

---

### MutationProfiler

Analyzes mutation patterns across sequence sets. Calculates position-specific amino acid frequencies for understanding sequence diversity.

**Installation**: This tool only needs the a small environment:
```bash
mamba create -n MutationEnv seaborn matplotlib pandas logomaker scipy
```

**Parameters**:
- `original`: Union[ToolOutput, StandardizedOutput] (required) - Original/reference sequences
- `mutants`: Union[ToolOutput, StandardizedOutput] (required) - Mutant sequences to analyze
- `include_original`: bool = True - Include original sequence in frequency analysis
- `positions`: Optional[str] = None - PyMOL-style position filter restricting the analysis/logos to selected positions (e.g., "141+143+145+147-149")
- `color_palette`: str = "okabe-ito" - Color palette for sequence logos. Options: "okabe-ito" (colorblind-friendly), "standard" (property-based shading)

**Tables**:
- `profile`:

  | position | original | count | frequency |
  |----------|----------|-------|-----------|

- `mutations`:

  | position | original | A | C | D | E | F | G | H | I | K | L | M | N | P | Q | R | S | T | V | W | Y |
  |----------|----------|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|

- `absolute_frequencies`:

  | position | original | A | C | ... | Y |
  |----------|----------|---|---|-----|---|

- `relative_frequencies`:

  | position | original | A | C | ... | Y |
  |----------|----------|---|---|-----|---|

**Example**:
```python
from biopipelines.mutation_profiler import MutationProfiler

profiler = MutationProfiler(
    original=template,
    mutants=lmpnn
)
```
---

### SequenceMetricCorrelation

Computes correlation signals between sequence mutations and performance metrics. Quantifies how specific mutations affect metrics using standardized effect sizes.

**Correlation Formulas**:

1D correlation (position-level):
```
c(i) = (m̄ᵢ₋ₘᵤₜ - m̄ᵢ₋wₜ) / √(s²ᵢ₋ₘᵤₜ + s²ᵢ₋wₜ)
```
Where:
- m̄ᵢ₋ₘᵤₜ = mean metric for sequences with mutation at position i
- m̄ᵢ₋wₜ = mean metric for sequences with wildtype at position i
- s²ᵢ₋ₘᵤₜ, s²ᵢ₋wₜ = unbiased variances (ddof=1)

2D correlation (amino acid-specific):
```
c(i,aa) = (m̄ᵢ,ₐₐ - m̄ᵢ,¬ₐₐ) / √(s²ᵢ,ₐₐ + s²ᵢ,¬ₐₐ)
```
Where:
- m̄ᵢ,ₐₐ = mean metric for sequences with amino acid 'aa' at position i
- m̄ᵢ,¬ₐₐ = mean metric for sequences with other amino acids at position i
- s²ᵢ,ₐₐ, s²ᵢ,¬ₐₐ = unbiased variances (ddof=1)

Note: When n ≤ 1 for either group, correlation is set to 0.

**Installation**: Requires the standard mutation-analysis Python environment with pandas, NumPy, matplotlib, and logomaker.

**Parameters**:
- `mutants`: Union[ToolOutput, StandardizedOutput, List[...]] (required) - Mutant sequences
- `data`: Union[ToolOutput, StandardizedOutput, TableInfo, str, List[...]] (required) - Table(s) with metric values
- `original`: Union[str, ToolOutput, StandardizedOutput] (required) - Reference sequence
- `metric`: str (required) - Metric column name to analyze
- `positions`: Optional[str] = None - PyMOL-style position filter (e.g., "141+143+145+147-149")

**Outputs**:
- `tables.correlation_1d`:

  | position | original | correlation | mean_mutated | mean_wt | var_mutated | var_wt | n_mutated | n_wt |
  |----------|-------|-------------|--------------|---------|-------------|--------|-----------|------|

- `tables.correlation_2d`:

  | position | original | A | C | ... | Y |
  |----------|-------|---|---|-----|---|

- `tables.sample_counts_2d`:

  | position | original | A | C | ... | Y |
  |----------|-------|---|---|-----|---|

**Example**:
```python
from biopipelines.sequence_metric_correlation import SequenceMetricCorrelation

# Single cycle
correlation = SequenceMetricCorrelation(
    mutants=filtered.tables.sequences,
    data=filtered.tables.merged,
    original=original_holo,
    metric="affinity_pred_value"
)

# Multi-cycle accumulation
correlation = SequenceMetricCorrelation(
    mutants=[cycle1.tables.sequences, cycle2.tables.sequences],
    data=[cycle1.tables.merged, cycle2.tables.merged],
    original=original_holo,
    metric="affinity_pred_value",
    positions="141+143+145+147-149"  # Filter plot to these positions
)
```

---
