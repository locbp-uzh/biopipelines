# DataStream I/O Standardization Refactor - COMPLETE

This document tracks the refactoring of HelpScripts to use `pipe_biopipelines_io.py` for standardized DataStream and table handling at SLURM runtime.

## Status: COMPLETE

All HelpScripts that iterate over file-based DataStreams have been refactored.

---

## Overview

The `pipe_biopipelines_io.py` module provides standardized functions for scripts that **consume file-based DataStreams** (iterate over structure files with ID matching):

### DataStream Utilities
- `load_datastream(source)` - Load DataStream from JSON file or dict
- `iterate_files(ds)` - Iterate `(id, resolved_file_path)` for file-based streams
- `iterate_values(ds, columns)` - Iterate `(id, value_dict)` for value-based streams
- `resolve_file(ds, item_id)` - Resolve single file for an ID (handles wildcards)
- `get_value(ds, item_id, column)` - Get single value from map_table
- `get_all_values(ds, item_id)` - Get all values from map_table

### Table Reference Utilities
- `load_table(reference)` - Load table from path or `DATASHEET_REFERENCE:path:column` string
- `lookup_table_value(table, item_id, column)` - Look up value for specific ID
- `iterate_table_values(table, item_ids, column)` - Iterate over `(id, value)` pairs

---

## Refactored HelpScripts (10 total)

These scripts iterate over file-based DataStreams and now use `pipe_biopipelines_io`:

| HelpScript | PipelineScript | Description |
|------------|----------------|-------------|
| `pipe_distance_selector.py` | `distance_selector.py` | Distance-based residue selection |
| `pipe_lmpnn_runtime_positions.py` | `ligand_mpnn.py` | LigandMPNN position handling |
| `pipe_pmpnn_fixed_positions.py` | `protein_mpnn.py` | ProteinMPNN fixed positions |
| `pipe_conformational_change.py` | `conformational_change.py` | RMSD calculation between structures |
| `pipe_pose_distance.py` | `pose_distance.py` | Pose comparison metrics |
| `pipe_protein_ligand_contacts.py` | `protein_ligand_contacts.py` | Contact analysis |
| `pipe_residue_atom_distance.py` | `residue_atom_distance.py` | Distance calculations |
| `pipe_plip_analysis.py` | `plip.py` | PLIP interaction analysis |
| `pipe_sasa.py` | `sasa.py` | Solvent accessible surface area |
| `pipe_selection_editor.py` | `selection_editor.py` | Structure selection editing |

---

## HelpScripts NOT Requiring Refactoring

The following scripts do NOT iterate over file-based DataStreams and therefore do not need `pipe_biopipelines_io`:

### Producer Scripts (Create DataStreams/Tables from tool output)
These parse tool output and CREATE result tables, not consume DataStreams:

| Script | Reason |
|--------|--------|
| `pipe_pdb.py` | Creates structures DataStream from input PDB files |
| `pipe_sequence.py` | Creates sequences DataStream from input |
| `pipe_ligand.py` | Creates ligands DataStream from input |
| `pipe_boltz_postprocessing.py` | Parses Boltz2 output, creates result tables |
| `pipe_boltz_results.py` | Extracts confidence metrics from Boltz2 output |
| `pipe_boltzgen.py` | Parses BoltzGen output structures |
| `pipe_boltzgen_config.py` | Generates BoltzGen config files |
| `pipe_boltzgen_import.py` | Imports BoltzGen structures |
| `pipe_boltzgen_merge.py` | Merges BoltzGen results |
| `pipe_esmfold.py` | Runs ESMFold on sequences CSV |
| `pipe_alphafold_confidence.py` | Parses AlphaFold confidence metrics |
| `pipe_dynamicbind_table.py` | Parses DynamicBind SDF output |
| `pipe_rfdiffusion_table.py` | Parses RFdiffusion logs |
| `pipe_rfdiffusion3_table.py` | Parses RFdiffusion3 logs |
| `pipe_rfdiffusion3_postprocess.py` | Postprocesses RFdiffusion3 output |
| `pipe_pmpnn_table.py` | Parses ProteinMPNN FASTA output |
| `pipe_compound_library.py` | Creates compound library DataStream |
| `pipe_smiles_library.py` | Validates SMILES input |
| `pipe_smiles_properties.py` | Calculates compound properties |

### Table-to-Table Transformations (CSV in/out)
These read/write tables directly without iterating structure files:

| Script | Reason |
|--------|--------|
| `pipe_split_chains.py` | Reads CSV, splits sequences, writes CSV |
| `pipe_dna_encoder.py` | Encodes sequences to DNA codons |
| `pipe_extract_metrics.py` | Aggregates metrics from multiple tables |
| `pipe_csv_to_fasta.py` | Format conversion |
| `pipe_fa_to_csv_fasta.py` | Format conversion |
| `pipe_stitch_sequences.py` | Sequence stitching from table |
| `pipe_mmseqs2_sequences.py` | Processes sequences for MMseqs2 |
| `pipe_mmseqs2_lcf_sequences.py` | MMseqs2 LCF processing |
| `pipe_panda.py` | General table operations |
| `pipe_remove_duplicates.py` | Table deduplication |
| `pipe_plot.py` | Table visualization |
| `pipe_bayesian_adjuster.py` | Bayesian score adjustment |
| `pipe_sequence_metric_correlation.py` | Correlation analysis |
| `pipe_mutation_composer.py` | Mutation generation from table |
| `pipe_mutation_profiler.py` | Mutation profiling |
| `pipe_site_directed_mutagenesis.py` | SDM operations |
| `pipe_fuse_queries.py` | Fusion sequence generation |
| `pipe_dynamicbind_prepare_ligands.py` | Ligand table preparation |
| `pipe_rfdaa_prepare_ligand.py` | RFdiffusionAA ligand prep |

### Sequence-Based Matching Scripts (MSA lookup by sequence, not ID)
These use sequence-based matching which is fundamentally different from ID-based DataStream iteration:

| Script | Reason |
|--------|--------|
| `pipe_boltz_config_unified.py` | MSA lookup by sequence (same sequence = same MSA) |
| `pipe_boltz_msa_copy.py` | MSA handling based on sequence identity |
| `pipe_alphafold_msas.py` | MSA sequence matching |

### Specialized Workflow Scripts
These have custom patterns that don't fit the DataStream iteration model:

| Script | Reason |
|--------|--------|
| `pipe_pymol.py` | Uses operation-based config with `structures_ref` |
| `pipe_onion_net.py` | Loads structures from custom table format |
| `pipe_boltz_completion_check.py` | Validates file existence, not iteration |
| `pipe_boltz_direct_sequence_config.py` | Config generator from CLI args |

### Utility Scripts
These handle pipeline orchestration, not structure iteration:

| Script | Reason |
|--------|--------|
| `pipe_check_completion.py` | Completion validation |
| `pipe_filter_execution.py` | Execution filtering logic |
| `pipe_load_output_filter.py` | Output loading |
| `pipe_propagate_missing.py` | Missing data handling |

---

## Refactoring Pattern

When a HelpScript iterates over structure files, use this pattern:

### Before (manual file handling)
```python
import json
import glob

with open(structures_json) as f:
    structures = json.load(f)

for struct_id in structures['ids']:
    pattern = f"{structures_dir}/{struct_id}*.pdb"
    matches = glob.glob(pattern)
    struct_file = matches[0]  # Fragile!
    process(struct_id, struct_file)
```

### After (using pipe_biopipelines_io)
```python
from pipe_biopipelines_io import load_datastream, iterate_files

structures_ds = load_datastream(structures_json)

for struct_id, struct_file in iterate_files(structures_ds):
    process(struct_id, struct_file)
```

---

## Key Insight

`pipe_biopipelines_io` is specifically for scripts that **consume file-based DataStreams** by iterating over structure files and matching them to IDs. Most HelpScripts in this repository are either:

1. **Producers** - Parse tool output and create result tables
2. **Table processors** - Transform CSV/table data without file iteration
3. **Sequence-based** - Match by sequence content, not ID
4. **Specialized** - Have custom workflows that don't fit the pattern

Only 10 scripts actually iterate over structure files with ID matching, and all 10 have been refactored.
