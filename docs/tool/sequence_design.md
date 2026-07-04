# Sequence Design

[← Back to Tool Reference](../tool_reference.md)

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
- `sequences`:

  | id | sequence | protein_sequence | organism | method |
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

### Frame2Seq

Fast structure-conditioned inverse folding. Frame2Seq is a non-autoregressive masked-language model that generates multiple sequences per backbone in a single forward pass — same role as ProteinMPNN (structure → sequence), but materially faster, with slightly higher native-sequence recovery on CATH 4.2. Output IDs follow the ProteinMPNN multiplier convention `<structure_id>_<n>`.

**References**: https://github.com/dakpinaroglu/Frame2seq · https://arxiv.org/abs/2312.02447

**Environment**: `frame2seq`

**Installation**: `Frame2Seq.install()` creates the env and pip-installs the package (weights ship with it). Runs on CPU; a GPU speeds up large/many inputs.

**Parameters**:
- `structures`: DataStream | StandardizedOutput (required) — Input backbones.
- `num_sequences`: int = 1 — Sequences to sample per structure.
- `temperature`: float = 1.0 — Sampling temperature (>0).
- `chain`: str = "A" — Chain to redesign.
- `omit_aa`: str = "" — Single-letter codes to exclude from sampling (e.g. `"CM"` to omit Cys and Met).
- `fixed`: str | (TableInfo, column) = "" — Residues to keep at their input identity. Chain-aware selection (`"A10-20+A30"`) or a table column reference. Mutually exclusive with `redesigned`.
- `redesigned`: str | (TableInfo, column) = "" — Residues to redesign; everything else on `chain` is held fixed.

**Streams**: `sequences`, `fasta`

**Tables**:
- `sequences`:

  | id | sequence | score | recovery | structures.id |
  |----|----------|-------|----------|---------------|

- `missing`: | id | removed_by | cause |

**Example**:
```python
from biopipelines.frame2seq import Frame2Seq

seqs = Frame2Seq(structures=rfd, num_sequences=10, temperature=0.5)
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

### LASErMPNN

Ligand-conditioned inverse folding with all-atom sidechain packing. LASErMPNN redesigns a protein sequence around a bound ligand and packs the sidechains in a single pass, emitting full-atom complex PDBs (designed sequence + rotamers + ligand). Unlike LigandMPNN it reads the ligand straight from the input PDB's HETATM records — there is no ligand-code argument — and it controls fixed vs designed positions through the input B-factor column (`fixed`/`redesigned` are stamped as B-factors at runtime and run with `--fix_beta`). Output IDs follow the multiplier convention `<structure_id>_<n>`.

**References**: https://github.com/polizzilab/LASErMPNN

**Environment**: `lasermpnn`

**Installation**: `LASErMPNN.install()` clones the repo (weights ship in-repo, no download) and creates the `lasermpnn` env from the vendored spec (`environments/lasermpnn.<device>.yaml` + `.pip.txt`: torch, the matching PyG extension wheels, ProDy, pydssp). The repo has no `setup.py`; it is run as a module (`python -m LASErMPNN.run_batch_inference`) from the clone's parent directory. Pass `device="cpu"` for a CPU-only install (default `"gpu"`). Requires a GPU for practical runtimes; falls back to CPU when no GPU is present.

**Parameters** (defaults match the upstream `run_batch_inference` defaults):
- `structures`: DataStream | StandardizedOutput (required) — Input complexes. PDB-only; the ligand is read from HETATM.
- `num_sequences`: int = 1 — Designs to generate per input (upstream `designs_per_input`).
- `temperature`: float | None = None — Sequence sampling temperature (`--sequence_temp`). None uses the model's built-in temperature.
- `first_shell_temperature`: float | None = None — Temperature for binding-site residues (`--first_shell_sequence_temp`).
- `chi_temperature`: float | None = None — Rotamer sampling temperature (`--chi_temp`).
- `model`: str = "default" — Weights: `"default"` (`laser_weights_0p1A_nothing_heldout`), `"ligandmpnn_split"`, or `"soluble"`.
- `designs_per_batch`: int = 30 — Designs per GPU batch (tune for GPU memory).
- `inputs_per_pass`: int = 5 — Input files processed per GPU pass.
- `disabled_residues`: str = "X,C" — Comma-separated one-letter codes to forbid (upstream default omits Cys).
- `fixed`: str | (TableInfo, column) = "" — Positions to hold at their input identity. Stamped as B-factor 1.0 with `--fix_beta`. Mutually exclusive with `redesigned`.
- `redesigned`: str | (TableInfo, column) = "" — Positions to design; everything else on `chain` is held fixed.
- `chain`: str = "A" — Default chain for chainless position input.
- `repack_only`: bool = False — Keep the input sequence, only repack sidechains (`--repack_only_input_sequence`).
- `repack_all`: bool = False — Repack all residues, including fixed ones (`--repack_all`).
- `ignore_ligand`: bool = False — Ignore the ligand during design (`--ignore_ligand`).
- `constrain_ala_gly`: bool = False — Cap Ala/Gly over-sampling in exposed non-secondary-structure regions (`-c`). Upstream flags this as generally helpful for realistic designs.
- `ala_budget`: int = 4 — Max Ala in the constrained region.
- `gly_budget`: int = 0 — Max Gly in the constrained region.

**Streams**: `structures` (full-atom designed complexes), `sequences`

**Tables**:
- `sequences`:

  | id | structures.id | sequence | score |
  |----|---------------|----------|-------|

  `score` is LASErMPNN's per-design mean log-probability.

- `missing`: | id | removed_by | kind | cause |

**Example**:
```python
from biopipelines.lasermpnn import LASErMPNN

# Design the pocket around the bound ligand, everything else fixed
laser = LASErMPNN(
    structures=rfdaa,
    num_sequences=8,
    temperature=0.1,
    redesigned=rfdaa.tables.structures.designed,
    constrain_ala_gly=True,
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
- `structures`: Union[DataStream, StandardizedOutput] (required) - Input structures
- `ligand`: Optional[Union[str, DataStream, StandardizedOutput]] = None - Compounds stream (`Ligand(code="LIG")` or any compounds-producing tool) or a 3-letter code naming the bound ligand for binding-site focus; the residue `code` is read from the stream at runtime
- `num_sequences`: int = 1 - Number of sequences per batch
- `fixed`: str | (TableInfo, column) = "" - Fixed positions (LigandMPNN format "A3 A4 A5" or table reference)
- `redesigned`: str | (TableInfo, column) = "" - Designed positions (LigandMPNN format or table reference)
- `design_within`: float = 5.0 - Distance in Angstroms from ligand for post-generation analysis only (does not control design). For actually designing residues within a distance, use [DistanceSelector](analysis.md#distanceselector) to select positions first.
- `chain`: str = "A" - Default chain ID applied to chainless position input (e.g. when positions are given as "10-20" without chain prefix)
- `model`: str = "v_32_010" - LigandMPNN model version (v_32_005, v_32_010, v_32_020, v_32_025)
- `num_batches`: int = 1 - Number of batches to run. Total sequences = num_sequences × num_batches
- `remove_duplicates`: bool = True - Drop duplicate sequences from the output
- `fill_gaps`: str = "G" - Fill gaps in the protein with an amino acid (default glycine).
- `temperature`: float = 0.0 - Sampling temperature (0.0 = argmax / deterministic)
- `bias_AA_per_residue`: str = "" - Per-residue amino-acid bias (LigandMPNN JSONL path or spec)
- `seed`: int = 0 - Random seed (0 = random)

**Streams**: `sequences`

**Tables**:
- `sequences`:

  | id | sequence | sample | T | seed | overall_confidence | ligand_confidence | seq_rec | gaps |
  |----|----------|--------|---|------|-------------------|-------------------|---------|------|

**Example**:
```python
from biopipelines.ligand_mpnn import LigandMPNN
from biopipelines.ligand import Ligand

lmpnn = LigandMPNN(
    structures=rfdaa,
    ligand=Ligand(code="LIG"),
    num_sequences=5,
    redesigned=rfdaa.tables.structures.designed
)
```

---

### Mutagenesis

Performs mutagenesis at specified positions. Generates systematic amino acid substitutions for experimental library design or computational scanning.

**Environment**: `MutationEnv`

**Parameters**:
- `original`: Union[DataStream, StandardizedOutput] (required) - Input structure/sequence
- `position`: Union[int, str, TableReference, StandardizedOutput] = None (required in practice) - Target position(s) for mutagenesis:
  - `int`: Fixed position (1-indexed) for all sequences
  - `str`: PyMOL-style selection (e.g., `"141+143+145-149"`)
  - `TableReference`: Per-row position lookup (e.g., `fuse.tables.sequences.L1`)
  - `StandardizedOutput`: From Selection tool (extracts `selections.selection` column)
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
- `combinatorial`: bool = False - When multiple positions are given, generate the full Cartesian product of substitutions across positions instead of mutating each position independently
- `msas`: Union[DataStream, StandardizedOutput] = None - Optional precomputed MSAs for the original protein(s) (e.g. from AlphaFold, MMseqs2, or the MSA tool). When provided, a synthetic per-mutant MSA is derived by copying the parent's MSA (matched on the original protein id) and substituting **only the query (first) row** at the mutated position(s); homolog rows pass through unchanged and the format (a3m/csv) is preserved. The emitted `msas` stream is keyed by the mutant ids, so it feeds straight into a downstream folding tool alongside the mutant sequences. No realignment is performed — point substitutions introduce no gaps, so alignment columns are unchanged. A mutant whose parent has no MSA is skipped (warned, not fabricated). When None (default), no `msas` stream is emitted.

**Streams**: `sequences`; `msas` (only when `msas=` is given)

**Tables**:
- `sequences`:

  | id | sequences.id | sequence | mutations | mutation_positions | original_aa | new_aa |
  |----|--------------|----------|-----------|--------------------|-------------|--------|

  When chaining multiple Mutagenesis steps, `mutations` accumulates (e.g., `A42V,G50L`) and `mutation_positions` uses PyMOL selection format (e.g., `42+50`).

- `missing`:

  | id | removed_by | cause |
  |----|------------|-------|

- `msas` (only when `msas=` is given):

  | id | sequences.id | original.id | sequence | msa_file |
  |----|--------------|-------------|----------|----------|

  `id` and `sequences.id` are both the mutant id (so downstream folding tools match the mutant query); `original.id` records the parent protein the MSA was derived from; `sequence` is the mutated query sequence.

**Example**:
```python
from biopipelines.mutagenesis import Mutagenesis

# Convert position 42 to alanine
sdm = Mutagenesis(original=template, position=42, mutate_to="A")

# Saturation mutagenesis at position 42 (excluding cysteine and proline)
sdm = Mutagenesis(original=template, position=42, mode="saturation", exclude="CP")

# Multiple positions
sdm = Mutagenesis(original=template, position="42+50+55-60", mode="saturation")

# Per-row positions from a table column (e.g., linker positions from Fuse)
sdm = Mutagenesis(original=fused, position=fused.tables.sequences.L1, mode="saturation")

# Positions from Selection tool
sdm = Mutagenesis(original=template, position=selection_output, mode="saturation")

# Reuse the wild-type MSA for all mutants: derive one synthetic MSA per mutant
# (query row substituted, homolog rows kept) and fold without re-querying.
af = AlphaFold(proteins=template)            # builds the original MSA
sdm = Mutagenesis(original=template, position=42, mode="saturation", msas=af)
folded = AlphaFold(proteins=sdm, msas=sdm)   # per-mutant MSAs, keyed by mutant id
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

### ProteinMPNN

Designs protein sequences for given backbone structures. Uses graph neural networks to optimize sequences for structure stability while respecting fixed/designed region constraints.

**References**: https://www.science.org/doi/10.1126/science.add2187.

**Installation**: Go to your data folder and clone the official repository (https://github.com/dauparas/ProteinMPNN). The model will then work in the same environment as RFdiffusion.
```bash
git clone https://github.com/dauparas/ProteinMPN
```

**Parameters**:
- `structures`: Union[DataStream, StandardizedOutput] (required) - Input structures
- `num_sequences`: int = 1 - Number of sequences per structure
- `fixed`: str | (TableInfo, column) = "" - Fixed positions (PyMOL selection or table reference)
- `redesigned`: str | (TableInfo, column) = "" - Redesigned positions (PyMOL selection or table reference)
- `chain`: str = "auto" - Chain to apply fixed positions ("auto" detects from input structure)
- `sampling_temp`: float = 0.1 - Sampling temperature
- `model_name`: str = "v_48_020" - ProteinMPNN model variant
- `soluble_model`: bool = False - Use the soluble protein model (see `SolubleMPNN` for a convenience wrapper that locks this on)
- `remove_duplicates`: bool = True - Drop duplicate sequences from the output
- `fill_gaps`: str = "G" - Fill gaps in the protein with an amino acid (default glycine).
- `bias_AA_jsonl`: str = "" - Path to a ProteinMPNN amino-acid bias JSONL
- `omit_AA_jsonl`: str = "" - Path to a ProteinMPNN per-position omit-AA JSONL
- `seed`: int = 0 - Random seed (0 = random)
- `ca_noise_std`: float = 0.0 - Std. dev. of Gaussian noise added to Cα coordinates before design

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

### SolubleMPNN

`ProteinMPNN` with the soluble model locked on. Identical to `ProteinMPNN` in every other respect — same parameters (minus `soluble_model`), same environment and output. Use it instead of `ProteinMPNN(..., soluble_model=True)` when designing for soluble expression.

**Example**:
```python
from biopipelines import SolubleMPNN

seqs = SolubleMPNN(structures=rfd, num_sequences=10)
```

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
- `add_start_codon`: bool = False — Prepend an ATG start codon to the gene if absent

**Tables**:
- `rbs`:

  | id | sequence | rbs_sequence | full_gene | dg_total | tir_predicted | target_tir | target_dg | spacing | dg_mrna_rrna | dg_start | dg_spacing | dg_mrna | dg_standby |
  |----|-------------|--------------|-----------|----------|---------------|------------|-----------|---------|-------------|----------|------------|---------|------------|

  - `full_gene` = `pre_sequence` + `rbs_sequence` + `sequence` (complete DNA ready for synthesis)

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

### StitchSequences

Combines a template sequence with two types of modifications: **substitutions** (position-to-position copying from equal-length sequences) and **indels** (segment replacement that can change sequence length). Generates all Cartesian product combinations.

**Environment**: `biopipelines`

**Parameters**:
- `template`: Union[str, DataStream, StandardizedOutput] - Base sequence (raw string or tool output). Optional if using concatenation mode.
- `substitutions`: Dict[str, Union[List[str], DataStream, StandardizedOutput]] = None - Position-to-position substitutions from equal-length sequences. For each position in the selection, the residue at that position in the substitution sequence replaces the residue at that position in the template.
  - Keys: Position strings like `"11-19"` or `"11-19+31-44"`, table references, OR a **marker key** with no digits (e.g. `"X"`) — fills every template position whose residue is one of the marker chars from the source at the same index. Use this to repair gap-marked sequences (e.g. the `X` padding `get_protein_sequence` writes where a structure has missing residues) against a reference; a marked position the source can't fill (its residue is also a marker) is left as-is — the sequence is still kept and the unfilled positions are logged (not removed). The reference is matched by the template's own id.
  - Values: tool output with sequences (must be same length as template)
- `indels`: Dict[str, Union[List[str], DataStream, StandardizedOutput]] = None - Segment replacements where each contiguous segment is replaced with the given sequence. Can change sequence length.
  - Keys: Position strings like `"50-55"` or `"6-7+9-10+17-18"`, or integers for concatenation mode
  - Values: List of raw sequences (each segment replaced with full sequence)
- `remove_duplicates`: bool = True - Drop duplicate sequences from the generated combinations

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

# Marker fill: repair X-gapped sequences against a reference (matched by id)
repaired = StitchSequences(
    template=gapped,                 # sequences with 'X' where residues were missing
    substitutions={"X": reference},  # fill each X from the reference at the same index
)
```

---
