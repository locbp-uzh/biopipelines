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
- `chain`: str = "auto" - Chain to apply fixed positions ("auto" detects from input structure)
- `plddt_threshold`: float = 100.0 - pLDDT threshold for automatic fixing (residues above threshold are fixed)
- `sampling_temp`: float = 0.1 - Sampling temperature
- `model_name`: str = "v_48_020" - ProteinMPNN model variant
- `soluble_model`: bool = True - Use soluble protein model
- `fill_gaps`: str = "G" - Fill gaps in the protein with an amino acid (default glycine). 

**Streams**: `sequences`

**Tables**:
- `sequences`:

  | id | structures.id | source_pdb | sequence | score | seq_recovery | rmsd | gaps |
  |----|---------------|------------|----------|-------|--------------|------|------|

**Note**: Sample 0 is the original/template sequence, samples 1+ are designs.

**Example**:
```python
from biopipelines.protein_mpnn import ProteinMPNN

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
- `num_sequences`: int = 1 - Number of sequences per batch
- `fixed`: str = "" - Fixed positions (LigandMPNN format "A3 A4 A5" or table reference)
- `redesigned`: str = "" - Designed positions (LigandMPNN format or table reference)
- `design_within`: float = 5.0 - Distance in Angstroms from ligand for post-generation analysis only (does not control design). For actually designing residues within a distance, use [DistanceSelector](Analysis.md#distanceselector) to select positions first.
- `model`: str = "v_32_010" - LigandMPNN model version (v_32_005, v_32_010, v_32_020, v_32_025)
- `num_batches`: int = 1 - Number of batches to run. Total sequences = num_sequences × num_batches
- `fill_gaps`: str = "G" - Fill gaps in the protein with an amino acid (default glycine). 

**Streams**: `sequences`

**Tables**:
- `sequences`:

  | id | sequence | sample | T | seed | overall_confidence | ligand_confidence | seq_rec | gaps |
  |----|----------|--------|---|------|-------------------|-------------------|---------|------|

**Example**:
```python
from biopipelines.ligand_mpnn import LigandMPNN

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

**Streams**: `sequences`

**Tables**:
- `sequences`:

  | id | sequence | mutations | mutation_positions |
  |----|----------|-----------|-------------------|

**Example**:
```python
from biopipelines.mutation_composer import MutationComposer
from biopipelines.mutation_profiler import MutationProfiler

profiler = MutationProfiler(original=ref, mutants=variants)
composer = MutationComposer(
    frequencies=profiler.tables.relative_frequencies,
    num_sequences=50,
    mode="weighted_random",
    max_mutations=5
)
```

---

### Mutagenesis

Performs mutagenesis at specified positions. Generates systematic amino acid substitutions for experimental library design or computational scanning.

**Environment**: `MutationEnv`

**Parameters**:
- `original`: Union[str, ToolOutput, StandardizedOutput] (required) - Input structure/sequence
- `position`: int (required) - Target position for mutagenesis (1-indexed)
- `mutate_to`: str = "" - Target amino acid(s) for "specific" mode (e.g., "A" for alanine, "AV" for alanine and valine). Required when mode is "specific".
- `mode`: str = "specific" - Mutagenesis strategy:
  - "specific": Only the amino acid(s) given in `mutate_to` (default)
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

**Streams**: `sequences`

**Tables**:
- `sequences`:

  | id | sequences.id | sequence | mutations | mutation_positions | original_aa | new_aa |
  |----|--------------|----------|-----------|--------------------|-------------|--------|

  When chaining multiple Mutagenesis steps, `mutations` accumulates (e.g., `A42V,G50L`) and `mutation_positions` uses PyMOL selection format (e.g., `42+50`).

- `missing`:

  | id | removed_by | cause |
  |----|------------|-------|

**Example**:
```python
from biopipelines.mutagenesis import Mutagenesis

# Convert position 42 to alanine
sdm = Mutagenesis(original=template, position=42, mutate_to="A")

# Saturation mutagenesis at position 42 (excluding cysteine and proline)
sdm = Mutagenesis(original=template, position=42, mode="saturation", exclude="CP")
```

---

### Fuse

Concatenates multiple sequences with flexible linkers. Creates fusion sequences with customizable linker lengths for domain engineering. Works with both protein and DNA sequences. Outputs include sequence/linker position columns in PyMOL selection format for easy visualization.

**Environment**: `biopipelines`

**Parameters**:
- `sequences`: Union[List[str], str] (required) - List of sequences or PDB file paths
- `name`: str = "" - Job name for output files
- `linker`: str = "GGGGSGGGGSGGGGSGGGGS" - Linker sequence that will be cut based on `linker_lengths` if specified
- `linker_lengths`: List[str] = None - List of length ranges for each junction to generate multiple variants by cutting the linker (e.g., ["1-6", "1-6"])

**Streams**: `sequences`

**Tables**:
- `sequences`:

  | id | sequence | lengths | S1 | L1 | S2 | L2 | S3 | ... |
  |----|----------|---------|----|----|----|----|----| --- |

  - `lengths`: Shortname of the lengths e.g. 2-4, 5-2-4, ...
  - `S1`, `S2`, `S3`, ...: Sequence positions in PyMOL selection format (e.g., "1-73", "76-237")
  - `L1`, `L2`, ...: Linker positions in PyMOL selection format (e.g., "74-75", "238-240")
  - Number of columns depends on number of input sequences: n sequences → n sequence columns (S1...Sn) and n-1 linker columns (L1...Ln-1)

**Example**:
```python
from biopipelines.fuse import Fuse
from biopipelines.pdb import PDB

N="GNH..."
mid=PDB("...")
C="EFT..."
fused = Fuse(
    sequences=[N, mid, C],
    linker="GSGAG",
    linker_lengths=["2-4", "2-4"],
    name="protein_fusion"
)
```

---

### StitchSequences

Combines a template sequence with two types of modifications: **substitutions** (position-to-position copying from equal-length sequences) and **indels** (segment replacement that can change sequence length). Generates all Cartesian product combinations.

**Environment**: `biopipelines`

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

**Streams**: `sequences`

**Tables**:
- `sequences`:

  | id | sequence |
  |----|----------|

**Examples**:
```python
from biopipelines.stitch_sequences import StitchSequences

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

**Environment**: `biopipelines`

**Parameters**:
- `sequences`: Union[ToolOutput, StandardizedOutput] (required) - Input sequences containing concatenated chains
- `split_positions`: List[int] (required) - Positions where to split sequences (1-indexed). For example, [197] splits at position 197 creating chains seq[0:197] and seq[197:]
- `chain_names`: Optional[List[str]] = None - Names for the chains (e.g., ["ChainA", "ChainB"]). If not provided, uses numeric suffixes (_1, _2, etc.)

**Streams**: `sequences`

**Tables**:
- `sequences`:

  | id | sequence | source_id | complex_id | chain_index | chain_name | chain_length |
  |----|----------|-----------|------------|-------------|------------|--------------|

Note: `complex_id` is an alias for `source_id` and is used by Boltz2 to group chains belonging to the same multi-chain complex.

**Example**:
```python
from biopipelines.split_chains import SplitChains

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

**Tables**:
- `dna`:

  | id | protein_sequence | dna_sequence | organism | method |
  |----|------------------|--------------|----------|--------|

- Excel file with color-coded codons (red <5‰, orange 5-10‰, black ≥10‰)

**Example**:
```python
from biopipelines.dna_encoder import DNAEncoder

dna = DNAEncoder(
    sequences=lmpnn,
    organism="EC&HS"  # Conservative optimization for both E. coli and human
)
```

**Note**: Uses thresholded weighted sampling (codons ≥10‰, fallback to ≥5‰). For multi-organism optimization, uses minimum frequency across organisms. Please cite CoCoPUTs (HIVE) when using.

---

### RBSDesigner

Designs synthetic ribosome binding sites (RBS) to control protein expression in bacteria. Uses the Salis thermodynamic model to predict translation initiation rates and a simulated annealing optimizer to design RBS sequences matching a target expression level. Requires ViennaRNA for RNA free energy calculations.

**Reference**: Salis, Mirsky & Voigt, *Nat. Biotechnol.* **27**, 946–950 (2009). doi:10.1038/nbt.1568

**Environment**: `rbs_designer` (ViennaRNA from bioconda, requires flexible channel priority)

**Installation**:
```python
RBSDesigner.install()
```

**Parameters**:
- `sequences`: Union[ToolOutput, StandardizedOutput] (required) — Input DNA sequences (typically from DNAEncoder)
- `tir`: Union[str, int, float] = "medium" — Target translation initiation rate:
  - "low" (100 au)
  - "medium" (1000 au)
  - "high" (10000 au)
  - "maximum" (100000 au)
  - Or any numeric value (au)
- `pre_sequence`: str = "" — Optional fixed 5'UTR DNA to prepend before the designed RBS

**Tables**:
- `rbs`:

  | id | dna_sequence | rbs_sequence | full_gene | dg_total | tir_predicted | target_tir | target_dg | spacing | dg_mrna_rrna | dg_start | dg_spacing | dg_mrna | dg_standby |
  |----|-------------|--------------|-----------|----------|---------------|------------|-----------|---------|-------------|----------|------------|---------|------------|

  - `full_gene` = `pre_sequence` + `rbs_sequence` + `dna_sequence` (complete DNA ready for synthesis)

**Example**:
```python
from biopipelines.dna_encoder import DNAEncoder
from biopipelines.rbs_designer import RBSDesigner

# Codon-optimize then design RBS for high expression in E. coli
dna = DNAEncoder(sequences=proteins, organism="EC")
rbs = RBSDesigner(sequences=dna, tir="high")

# Design RBS for specific TIR with 5'UTR prefix
rbs = RBSDesigner(sequences=dna, tir=5000, pre_sequence="AATTAA")
```

**Thermodynamic model** (Equation 2):
```
dG_tot = dG_mRNA:rRNA + dG_start + dG_spacing - dG_standby - dG_mRNA
```
- `dG_mRNA:rRNA`: SD / anti-SD hybridization energy (ViennaRNA duplexfold)
- `dG_start`: Start codon identity (AUG = −1.194, GUG = −0.075 kcal/mol)
- `dG_spacing`: Penalty for non-optimal SD-to-start-codon distance
- `dG_standby`: Energy to unfold the 4-nt standby site upstream of SD
- `dG_mRNA`: Local mRNA folding energy (70 nt window around start codon)

**Note**: RBS design uses simulated annealing with adaptive temperature control (5–20% acceptance rate). Each sequence typically requires thousands of energy evaluations. Computation time scales with the number of input sequences. Please cite Salis et al. 2009 when using.

---
