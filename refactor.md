# DataStream I/O Standardization Refactor Checklist

This checklist tracks the refactoring of HelpScripts to use `pipe_biopipelines_io.py` for standardized DataStream and table handling at SLURM runtime.

## Overview

The `pipe_biopipelines_io.py` module provides standardized functions for:

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

## HelpScripts Refactoring Status

### Already Using `pipe_biopipelines_io`
- [x] `pipe_distance_selector.py`
- [x] `pipe_lmpnn_runtime_positions.py`
- [x] `pipe_pmpnn_fixed_positions.py`

### Entity HelpScripts (Priority: High)
These scripts process basic entities and should standardize early.

| HelpScript | Associated PipelineScript | Status | Notes |
|------------|---------------------------|--------|-------|
| `pipe_pdb.py` | `pdb.py` | [ ] | Structure loading, sequence extraction |
| `pipe_sequence.py` | `sequence.py` | [ ] | Sequence file creation |
| `pipe_ligand.py` | `ligand.py` | [ ] | Ligand SMILES/SDF handling |
| `pipe_compound_library.py` | `compound_library.py` | [ ] | Compound expansion |
| `pipe_smiles_library.py` | `compound_library.py` | [ ] | SMILES validation |
| `pipe_smiles_properties.py` | `compound_library.py` | [ ] | Property calculation |

### Structure Prediction HelpScripts (Priority: High)
Core prediction tools that process DataStreams extensively.

| HelpScript | Associated PipelineScript | Status | Notes |
|------------|---------------------------|--------|-------|
| `pipe_boltz_config_unified.py` | `boltz2.py` | N/A | Uses sequence-based MSA matching |
| `pipe_boltz_postprocessing.py` | `boltz2.py` | N/A | Producer script |
| `pipe_boltz_results.py` | `boltz2.py` | [ ] | Confidence metrics |
| `pipe_boltz_completion_check.py` | `boltz2.py` | [ ] | Job completion |
| `pipe_boltz_msa_copy.py` | `boltz2.py` | N/A | Uses sequence-based MSA matching |
| `pipe_boltz_direct_sequence_config.py` | `boltz2.py` | [ ] | Direct sequence mode |
| `pipe_boltzgen.py` | `boltzgen.py` | [ ] | BoltzGen processing |
| `pipe_boltzgen_config.py` | `boltzgen.py` | [ ] | BoltzGen config |
| `pipe_boltzgen_import.py` | `boltzgen.py` | [ ] | BoltzGen imports |
| `pipe_boltzgen_merge.py` | `boltzgen.py` | [ ] | BoltzGen merging |
| `pipe_esmfold.py` | `esmfold.py` | [ ] | ESMFold processing |
| `pipe_alphafold_msas.py` | `alphafold.py` | N/A | Uses sequence-based MSA matching |
| `pipe_alphafold_confidence.py` | `alphafold.py` | [ ] | AlphaFold confidence |
| `pipe_dynamicbind_table.py` | `dynamic_bind.py` | [ ] | DynamicBind results |
| `pipe_dynamicbind_prepare_ligands.py` | `dynamic_bind.py` | [ ] | Ligand preparation |

### Structure Design HelpScripts (Priority: High)
Design tools that read structures and output sequences.

| HelpScript | Associated PipelineScript | Status | Notes |
|------------|---------------------------|--------|-------|
| `pipe_rfdiffusion_table.py` | `rfdiffusion.py` | [ ] | RFdiffusion results |
| `pipe_rfdiffusion3_table.py` | `rfdiffusion3.py` | [ ] | RFdiffusion3 results |
| `pipe_rfdiffusion3_postprocess.py` | `rfdiffusion3.py` | [ ] | RFdiffusion3 postprocess |
| `pipe_rfdaa_prepare_ligand.py` | `rfdiffusion_allatom.py` | [ ] | RFdiffusionAA ligand prep |
| `pipe_pmpnn_table.py` | `protein_mpnn.py` | [ ] | ProteinMPNN results |

### Sequence/MSA Processing HelpScripts (Priority: Medium)

| HelpScript | Associated PipelineScript | Status | Notes |
|------------|---------------------------|--------|-------|
| `pipe_mmseqs2_sequences.py` | `mmseqs2.py` | [ ] | MMseqs2 sequence output |
| `pipe_mmseqs2_lcf_sequences.py` | `mmseqs2_lcf.py` | [ ] | MMseqs2 LCF sequences |
| `pipe_csv_to_fasta.py` | Multiple | [ ] | CSV to FASTA conversion |
| `pipe_fa_to_csv_fasta.py` | Multiple | [ ] | FASTA processing |
| `pipe_stitch_sequences.py` | `stitch_sequences.py` | [ ] | Sequence stitching |
| `pipe_dna_encoder.py` | `dna_encoder.py` | [ ] | DNA encoding |

### Analysis/Metrics HelpScripts (Priority: Medium)

| HelpScript | Associated PipelineScript | Status | Notes |
|------------|---------------------------|--------|-------|
| `pipe_extract_metrics.py` | `extract_metrics.py` | [ ] | Metric extraction |
| `pipe_conformational_change.py` | `conformational_change.py` | [x] | RMSD calculation |
| `pipe_pose_distance.py` | `pose_distance.py` | [x] | Pose comparison |
| `pipe_protein_ligand_contacts.py` | `protein_ligand_contacts.py` | [x] | Contact analysis |
| `pipe_residue_atom_distance.py` | `residue_atom_distance.py` | [x] | Distance calculation |
| `pipe_plip_analysis.py` | `plip.py` | [x] | PLIP interactions |
| `pipe_sasa.py` | `sasa.py` | [x] | Solvent accessibility |
| `pipe_onion_net.py` | `onion_net.py` | [ ] | OnionNet scoring |

### Table/Data Processing HelpScripts (Priority: Medium)

| HelpScript | Associated PipelineScript | Status | Notes |
|------------|---------------------------|--------|-------|
| `pipe_panda.py` | `panda.py` | [ ] | Table operations |
| `pipe_remove_duplicates.py` | `remove_duplicates.py` | [ ] | Deduplication |
| `pipe_plot.py` | `plot.py` | [ ] | Plotting |
| `pipe_bayesian_adjuster.py` | `bayesian_adjuster.py` | [ ] | Bayesian scoring |
| `pipe_sequence_metric_correlation.py` | `sequence_metric_correlation.py` | [ ] | Correlation analysis |

### Mutation HelpScripts (Priority: Medium)

| HelpScript | Associated PipelineScript | Status | Notes |
|------------|---------------------------|--------|-------|
| `pipe_mutation_composer.py` | `mutation_composer.py` | [ ] | Mutation composition |
| `pipe_mutation_profiler.py` | `mutation_profiler.py` | [ ] | Mutation profiling |
| `pipe_site_directed_mutagenesis.py` | `site_directed_mutagenesis.py` | [ ] | SDM |

### Structure Manipulation HelpScripts (Priority: Low)

| HelpScript | Associated PipelineScript | Status | Notes |
|------------|---------------------------|--------|-------|
| `pipe_split_chains.py` | `split_chains.py` | [ ] | Chain splitting |
| `pipe_selection_editor.py` | `selection_editor.py` | [x] | Selection editing |
| `pipe_fuse_queries.py` | `fuse.py` | [ ] | Structure fusion |
| `pipe_pymol.py` | `pymol.py` | [ ] | PyMOL rendering |

### Utility HelpScripts (Priority: Low)

| HelpScript | Associated PipelineScript | Status | Notes |
|------------|---------------------------|--------|-------|
| `pipe_check_completion.py` | Multiple | [ ] | Completion checking |
| `pipe_filter_execution.py` | Multiple | [ ] | Execution filtering |
| `pipe_load_output_filter.py` | `load.py` | [ ] | Output loading |
| `pipe_propagate_missing.py` | Multiple | [ ] | Missing data handling |

---

## PipelineScripts Without Dedicated HelpScripts

These PipelineScripts either don't require SLURM-time Python execution or use inline bash:

- `base_config.py` - Base class (no HelpScript needed)
- `combinatorics.py` - Bundle/Each logic (pipeline time only)
- `converters.py` - Utility converters
- `datastream.py` - DataStream class definition
- `datastream_resolver.py` - Resolver logic
- `entities.py` - Entity imports
- `file_paths.py` - Path utilities
- `folders.py` - Folder management
- `pipeline.py` - Pipeline orchestration
- `table.py` - Table entity
- `table_utils.py` - Table utilities
- `config_manager.py` - Configuration

---

## Refactoring Pattern

When refactoring a HelpScript to use `pipe_biopipelines_io`, follow this pattern:

### Before (manual file handling)
```python
import pandas as pd
import glob
import json

# Manual structure iteration
with open(structures_json) as f:
    structures = json.load(f)

for struct_id in structures['ids']:
    # Manual wildcard resolution
    pattern = f"{structures_dir}/{struct_id}*.pdb"
    matches = glob.glob(pattern)
    struct_file = matches[0]  # Fragile!

    # Manual table lookup
    table = pd.read_csv(table_path)
    row = table[table['id'] == struct_id]
    value = row['column'].iloc[0]
```

### After (using pipe_biopipelines_io)
```python
from pipe_biopipelines_io import (
    load_datastream, iterate_files,
    load_table, lookup_table_value
)

# Load DataStream
structures_ds = load_datastream(structures_json)

# Iterate with automatic wildcard resolution
for struct_id, struct_file in iterate_files(structures_ds):
    # Clean table lookup with proper ID matching
    table, column = load_table(table_reference)
    value = lookup_table_value(table, struct_id, column)
```

---

## Testing Checklist

After refactoring each HelpScript:

- [ ] Verify DataStream loading works with JSON files
- [ ] Verify wildcard patterns resolve correctly
- [ ] Verify table lookups match by id/pdb columns
- [ ] Test with real pipeline output
- [ ] Run associated example pipeline

---

## HelpScripts NOT Requiring Refactoring

These scripts don't iterate over DataStream files at SLURM time:

### Producer Scripts (Create DataStreams/Tables)
These scripts are "producers" that create DataStreams from raw input rather than consuming them:
- `pipe_pdb.py` - Creates structures DataStream from input PDB files
- `pipe_sequence.py` - Creates sequences DataStream from input
- `pipe_boltz_postprocessing.py` - Processes Boltz2 output, creates result tables
- `pipe_rfdiffusion_table.py` - Parses RFdiffusion logs, creates result tables
- `pipe_rfdiffusion3_table.py` - Similar for RFdiffusion3
- `pipe_pmpnn_table.py` - Creates ProteinMPNN result tables

### Table-to-Table Transformations
These scripts read/write tables directly without iterating files:
- `pipe_split_chains.py` - Reads CSV, splits sequences, writes CSV
- `pipe_dna_encoder.py` - Reads sequences CSV, encodes to DNA, writes CSV
- `pipe_extract_metrics.py` - Extracts metrics from multiple tables
- `pipe_fuse_queries.py` - Generates fusion sequences from command-line input

### Sequence-Based Matching Scripts
These scripts use sequence-based matching rather than ID-based matching (by design):
- `pipe_boltz_config_unified.py` - MSA lookup by sequence (same sequence = same MSA)
- `pipe_boltz_msa_copy.py` - MSA handling based on sequence identity
- `pipe_alphafold_msas.py` - Similar MSA sequence matching

### Complex Operation Scripts
These scripts have specialized workflows that don't fit the standard pattern:
- `pipe_pymol.py` - Uses operation-based config with `structures_ref` pattern
- `pipe_onion_net.py` - Uses table-based structure loading

---

## Notes

- Priority is based on frequency of use and complexity of current file handling
- Some HelpScripts may only need partial adoption (e.g., just table lookups)
- The `pipe_biopipelines_io` module handles caching of loaded tables
- Error messages from the standardized functions are more informative
