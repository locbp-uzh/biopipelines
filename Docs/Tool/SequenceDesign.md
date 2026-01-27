# Sequence Design

[← Back to Tool Reference](../ToolReference.md)

---

### ProteinMPNN

Designs protein sequences for given backbone structures. Uses graph neural networks to optimize sequences for structure stability while respecting fixed/designed region constraints.

**References**: https://www.science.org/doi/10.1126/science.add2187.

**Installation**: Go to your data folder and clone the official repository (https://github.com/dauparas/ProteinMPNN). The model will then work in the same environment as RFdiffusion.
```bash
git clone https://github.com/dauparas/ProteinMPN
```

**Parameters**:
- `structures`: Union[str, List[str], ToolOutput] (required) - Input structures
- `tables`: Optional[List[str]] = None - Input table files
- `num_sequences`: int = 1 - Number of sequences per structure
- `fixed`: str = "" - Fixed positions (PyMOL selection or table reference)
- `redesigned`: str = "" - Redesigned positions (PyMOL selection or table reference)
- `fixed_chain`: str = "A" - Chain to apply fixed positions
- `plddt_threshold`: float = 100.0 - pLDDT threshold for automatic fixing (residues above threshold are fixed)
- `sampling_temp`: float = 0.1 - Sampling temperature
- `model_name`: str = "v_48_020" - ProteinMPNN model variant
- `soluble_model`: bool = True - Use soluble protein model

**Outputs**:
- `sequences`: CSV file with generated sequences
- `tables.sequences`:

  | id | source_id | source_pdb | sequence | score | seq_recovery | rmsd |
  |----|-----------|------------|----------|-------|--------------|------|

**Note**: Sample 0 is the original/template sequence, samples 1+ are designs.

**Example**:
```python
from PipelineScripts.protein_mpnn import ProteinMPNN

pmpnn = ProteinMPNN(
    structures=rfd,
    num_sequences=10,
    fixed="1-10+50-60",
    redesigned="20-40"
)
```

---

### LigandMPNN

Designs protein sequences optimized for ligand binding. Specialized version of ProteinMPNN that considers protein-ligand interactions during sequence design.

**References**: https://www.nature.com/articles/s41592-025-02626-1.

**Installation**: As from the official repository (https://github.com/dauparas/LigandMPNN), go to your data folder then run:
```bash
git clone https://github.com/dauparas/LigandMPNN.git
cd LigandMPNN
bash get_model_params.sh "./model_params"
mamba create -n ligandmpnn_env python=3.11
pip3 install -r requirements.txt
```

**Parameters**:
- `structures`: Union[str, List[str], ToolOutput] (required) - Input structures
- `ligand`: str (required) - Ligand identifier for binding site focus
- `tables`: Optional[List[str]] = None - Input table files
- `name`: str = "" - Job name for output files
- `num_sequences`: int = 1 - Number of sequences to generate per structure
- `fixed`: str = "" - Fixed positions (LigandMPNN format "A3 A4 A5" or table reference)
- `redesigned`: str = "" - Designed positions (LigandMPNN format or table reference)
- `design_within`: float = 5.0 - Distance in Angstroms from ligand for post-generation analysis only (does not control design). For actually designing residues within a distance, use [DistanceSelector](Analysis.md#distanceselector) to select positions first.
- `model`: str = "v_32_010" - LigandMPNN model version (v_32_005, v_32_010, v_32_020, v_32_025)
- `batch_size`: int = 1 - Batch size for processing

**Outputs**:
- `sequences`: CSV file with generated sequences
- `tables.sequences`:

  | id | sequence | sample | T | seed | overall_confidence | ligand_confidence | seq_rec |
  |----|----------|--------|---|------|-------------------|-------------------|---------|

**Example**:
```python
from PipelineScripts.ligand_mpnn import LigandMPNN

lmpnn = LigandMPNN(
    structures=rfdaa,
    ligand="LIG",
    num_sequences=5,
    redesigned=rfdaa.tables.structures.designed
)
```

---

### MutationComposer

Generates new protein sequences by composing mutations based on frequency analysis. Creates combinatorial mutants from mutation profiles with different sampling strategies.

**Installation**: Same environment as MutationProfiler.

**Parameters**:
- `frequencies`: Union[List, ToolOutput, StandardizedOutput, TableInfo, str] (required) - Mutation frequency table(s) from MutationProfiler
- `num_sequences`: int = 10 - Number of sequences to generate
- `mode`: str = "single_point" - Generation strategy:
  - "single_point": One mutation per sequence
  - "weighted_random": Random mutations weighted by frequency
  - "hotspot_focused": Focus on high-frequency positions
  - "top_mutations": Use only top N mutations
- `min_frequency`: float = 0.01 - Minimum frequency threshold for mutations
- `max_mutations`: int = None - Maximum mutations per sequence
- `random_seed`: int = None - Random seed for reproducibility
- `prefix`: str = "" - Prefix for sequence IDs
- `hotspot_count`: int = 10 - Number of top hotspot positions (for hotspot_focused mode)
- `combination_strategy`: str = "average" - Strategy for combining multiple tables (average, maximum, stack, round_robin)

**Outputs**:
- `sequences`: CSV file with composed sequences
- `tables.sequences`:

  | id | sequence | mutations | mutation_positions |
  |----|----------|-----------|-------------------|

**Example**:
```python
from PipelineScripts.mutation_composer import MutationComposer
from PipelineScripts.mutation_profiler import MutationProfiler

profiler = MutationProfiler(original=ref, mutants=variants)
composer = MutationComposer(
    frequencies=profiler.tables.relative_frequencies,
    num_sequences=50,
    mode="weighted_random",
    max_mutations=5
)
```

---

### SDM (SiteDirectedMutagenesis)

Performs site-directed mutagenesis at specified positions. Generates systematic amino acid substitutions for experimental library design or computational scanning.

**Environment**: `ProteinEnv`

**Parameters**:
- `original`: Union[str, ToolOutput, StandardizedOutput] (required) - Input structure/sequence
- `position`: int (required) - Target position for mutagenesis (1-indexed)
- `mode`: str = "saturation" - Mutagenesis strategy:
  - "saturation": All 20 amino acids
  - "hydrophobic": Hydrophobic residues only
  - "hydrophilic": Hydrophilic residues only
  - "charged": Charged residues only
  - "polar": Polar residues only
  - "nonpolar": Nonpolar residues only
  - "aromatic": Aromatic residues only
  - "aliphatic": Aliphatic residues only
  - "positive": Positively charged residues only
  - "negative": Negatively charged residues only
- `include_original`: bool = False - Include original amino acid in output
- `exclude`: str = "" - Amino acids to exclude (single letter codes as string, e.g., "CP")
- `prefix`: str = "" - Prefix for sequence IDs

**Outputs**:
- `sequences`: CSV file with mutant sequences
- `tables.sequences`:

  | id | sequence | mutation | position | original_aa | new_aa |
  |----|----------|----------|----------|-------------|--------|

- `tables.missing_sequences`:

  | id | sequence | reason |
  |----|----------|--------|

**Example**:
```python
from PipelineScripts.site_directed_mutagenesis import SDM

sdm = SDM(
    original=template,
    position=42,
    mode="saturation",
    exclude="CP"
)
```

---

### Fuse

Concatenates multiple protein sequences with flexible linkers. Creates fusion proteins with customizable linker lengths for domain engineering. Outputs include domain/linker position columns in PyMOL selection format for easy visualization.

**Environment**: `ProteinEnv`

**Parameters**:
- `proteins`: Union[List[str], str] (required) - List of protein sequences or PDB file paths
- `sequences`: Union[List[str], str] = None - Alias for proteins
- `name`: str = "" - Job name for output files
- `linker`: str = "GGGGSGGGGSGGGGSGGGGS" - Linker sequence that will be cut based on `linker_lengths` if specified
- `linker_lengths`: List[str] = None - List of length ranges for each junction to generate multiple variants by cutting the linker (e.g., ["1-6", "1-6"])

**Outputs**:
- `sequences`: CSV file with fused sequences
- `tables.sequences`:

  | id | sequence | lengths | D1 | L1 | D2 | L2 | D3 | ... |
  |----|----------|---------|----|----|----|----|----| --- |

  - `lengths`: Shortname of the lengths e.e. 2_4, 5_2_4, ...
  - `D1`, `D2`, `D3`, ...: Domain positions in PyMOL selection format (e.g., "1-73", "76-237")
  - `L1`, `L2`, ...: Linker positions in PyMOL selection format (e.g., "74-75", "238-240")
  - Number of columns depends on number of input proteins: n proteins → n domain columns (D1...Dn) and n-1 linker columns (L1...Ln-1)

**Example**:
```python
from PipelineScripts.fuse import Fuse
from PipelineScripts.pdb import PDB

N="GNH..."
mid=PDB("...")
C="EFT..."
fused = Fuse(
    proteins=[N, mid, C],
    linker="GSGAG",
    linker_lengths=["2-4", "2-4"],
    name="protein_fusion"
)
```

---

### StitchSequences

Combines a template sequence with two types of modifications: **substitutions** (position-to-position copying from equal-length sequences) and **indels** (segment replacement that can change sequence length). Generates all Cartesian product combinations.

**Environment**: `ProteinEnv`

**Parameters**:
- `template`: Union[str, ToolOutput, StandardizedOutput] - Base sequence (raw string or tool output). Optional if using concatenation mode.
- `substitutions`: Dict[str, Union[List[str], ToolOutput]] = None - Position-to-position substitutions from equal-length sequences. For each position in the selection, the residue at that position in the substitution sequence replaces the residue at that position in the template.
  - Keys: Position strings like `"11-19"` or `"11-19+31-44"`, or table references
  - Values: ToolOutput with sequences (must be same length as template)
- `indels`: Dict[str, Union[List[str], ToolOutput]] = None - Segment replacements where each contiguous segment is replaced with the given sequence. Can change sequence length.
  - Keys: Position strings like `"50-55"` or `"6-7+9-10+17-18"`, or integers for concatenation mode
  - Values: List of raw sequences (each segment replaced with full sequence)
- `id_map`: Dict[str, str] = {"*": "*_<N>"} - ID mapping pattern for matching sequences

**Position Syntax**:
- `"10-20"` → positions 10 to 20 (inclusive, 1-indexed)
- `"10-20+30-40"` → positions 10-20 and 30-40
- `"145+147+150"` → specific positions 145, 147, and 150

**Processing Order**:
1. Substitutions are applied first (position-to-position, same length)
2. Indels are applied second (segment replacement, can change length)

**Indel Segment Behavior**: For discontinuous selections like `"6-7+9-10+17-18"`, each contiguous segment is replaced with the full replacement sequence. So `"6-7+9-10": "GP"` replaces segment 6-7 with "GP" AND segment 9-10 with "GP".

**Concatenation Mode**: When `template` is omitted and `indels` keys are integers (1, 2, 3...), sequences are concatenated in order.

**Outputs**:
- `sequences`: CSV file with stitched sequences
- `tables.sequences`:

  | id | sequence |
  |----|----------|

**Examples**:
```python
from PipelineScripts.stitch_sequences import StitchSequences

# Position-to-position substitution from ToolOutput
# Both template and substitution sequences are 180 residues
stitched = StitchSequences(
    template=pmpnn,
    substitutions={
        "11-19+31-44": lmpnn  # Copy residues at these positions from lmpnn
    }
)

# Segment replacement with indels
stitched = StitchSequences(
    template="MKTAYIAKQRQISFVKSHFS...",
    indels={
        "11-15": ["AAAAA", "GGGGG"],  # Replace segment with 5-char options
        "20-22": ["XX", "YYY", "ZZZZ"]  # Can change length
    }
)
# Output: 2 × 3 = 6 combinations

# Combined: substitutions then indels
stitched = StitchSequences(
    template=pmpnn,
    substitutions={
        "6-12+19+21": lmpnn  # Position-to-position from lmpnn
    },
    indels={
        "50-55": ["LINKER", "GGG"]  # Replace segment 50-55
    }
)

# Discontinuous indel: replace multiple segments with same sequence
stitched = StitchSequences(
    template="ABCDEFGHIJKLMNOPQRSTUVWXYZ",
    indels={
        "3-4+7-8+11-12": ["XX", "YY"]  # Each segment replaced with "XX" or "YY"
    }
)
# "3-4+7-8+11-12": "XX" -> "ABXXEFXXIJXXMNOPQRSTUVWXYZ"

# ToolOutput with table-based positions
stitched = StitchSequences(
    template=pmpnn,
    substitutions={
        distances.tables.selections.within: lmpnn
    }
)

# Concatenation mode (no template, integer keys in indels)
stitched = StitchSequences(
    indels={
        1: ["AAAA", "BBBB"],         # First segment
        2: ["CCCC"],                 # Second segment
        3: ["DDDD", "EEEE", "FFFF"]  # Third segment
    }
)
# Output: 2 × 1 × 3 = 6 concatenated sequences
```

---

### SplitChains

Splits concatenated single-chain sequences into multi-chain sequences. Takes sequences where multiple protein chains have been fused into a single sequence and separates them into individual chains based on specified split positions. Useful for preparing multi-chain inputs for structure prediction tools.

**Environment**: None (pure Python)

**Parameters**:
- `sequences`: Union[ToolOutput, StandardizedOutput] (required) - Input sequences containing concatenated chains
- `split_positions`: List[int] (required) - Positions where to split sequences (1-indexed). For example, [197] splits at position 197 creating chains seq[0:197] and seq[197:]
- `chain_names`: Optional[List[str]] = None - Names for the chains (e.g., ["ChainA", "ChainB"]). If not provided, uses numeric suffixes (_1, _2, etc.)

**Outputs**:
- `sequences`: CSV file with split chain sequences
- `tables.sequences`:

  | id | sequence | source_id | complex_id | chain_index | chain_name | chain_length |
  |----|----------|-----------|------------|-------------|------------|--------------|

Note: `complex_id` is an alias for `source_id` and is used by Boltz2 to group chains belonging to the same multi-chain complex.

**Example**:
```python
from PipelineScripts.split_chains import SplitChains

# Split a 400-residue fusion into two chains at position 200
split = SplitChains(
    sequences=fused_sequences,
    split_positions=[200],
    chain_names=["Heavy", "Light"]
)

# Split into three chains
split = SplitChains(
    sequences=lmpnn,
    split_positions=[150, 300],
    chain_names=["A", "B", "C"]
)
```

---

### DNAEncoder

Reverse-translates protein sequences to DNA with organism-specific codon optimization. Uses thresholded weighted codon sampling based on CoCoPUTs genome frequency tables.

**Environment**: `biopipelines`

**Parameters**:
- `sequences`: Union[ToolOutput, StandardizedOutput] (required) - Input protein sequences
- `organism`: str = "EC" - Target organism for codon optimization:
  - "EC" (Escherichia coli)
  - "SC" (Saccharomyces cerevisiae)
  - "HS" (Homo sapiens)
  - Combinations: "EC&HS", "EC&SC", "HS&SC", "EC&HS&SC"

**Outputs**:
- `tables.dna`:

  | id | protein_sequence | dna_sequence | organism | method |
  |----|------------------|--------------|----------|--------|

- Excel file with color-coded codons (red <5‰, orange 5-10‰, black ≥10‰)

**Example**:
```python
from PipelineScripts.dna_encoder import DNAEncoder

dna = DNAEncoder(
    sequences=lmpnn,
    organism="EC&HS"  # Conservative optimization for both E. coli and human
)
```

**Note**: Uses thresholded weighted sampling (codons ≥10‰, fallback to ≥5‰). For multi-organism optimization, uses minimum frequency across organisms. Please cite CoCoPUTs (HIVE) when using.

---
