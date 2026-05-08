# BioPipelines Tool Index (SI)

Single-source-of-truth listing for every public-API tool. Generated from
`grep TOOL_NAME = biopipelines/*.py` and aligned with the categories in
[`tool_reference.md`](tool_reference.md).

**Public-API count: 45** (51 `TOOL_NAME` entries in the codebase, minus 4 internal
helpers — `BoltzGenMerge`, `BoltzGenImport`, `RFDAA_PrepareLigand`, `Mock` —
and 2 base-class scaffolding entries: `base`, `install`).

| # | Tool | Category | Source file |
|---|------|----------|-------------|
| 1  | ADMETAI                   | Analysis             | `biopipelines/admet_ai.py` |
| 2  | AlphaFold                 | Structure Prediction | `biopipelines/alphafold.py` |
| 3  | Angle                     | Analysis             | `biopipelines/angle.py` |
| 4  | BayesianAdjuster          | Statistics           | `biopipelines/bayesian_adjuster.py` |
| 5  | Boltz2                    | Structure Prediction | `biopipelines/boltz2.py` |
| 6  | BoltzGen                  | Structure Generation | `biopipelines/boltzgen.py` |
| 7  | CABSflex                  | Analysis             | `biopipelines/cabsflex.py` |
| 8  | CompoundLibrary           | Utilities            | `biopipelines/compound_library.py` |
| 9  | ConformationalChange      | Analysis             | `biopipelines/conformational_change.py` |
| 10 | Contacts                  | Analysis             | `biopipelines/contacts.py` |
| 11 | Distance                  | Analysis             | `biopipelines/distance.py` |
| 12 | DistanceSelector          | Analysis             | `biopipelines/distance_selector.py` |
| 13 | DNAEncoder                | Sequence Design      | `biopipelines/dna_encoder.py` |
| 14 | ExtractMetrics            | Data Management      | `biopipelines/extract_metrics.py` |
| 15 | Fuse                      | Sequence Design      | `biopipelines/fuse.py` |
| 16 | Gnina                     | Structure Prediction | `biopipelines/gnina.py` |
| 17 | Ligand                    | Utilities            | `biopipelines/ligand.py` |
| 18 | LigandMPNN                | Sequence Design      | `biopipelines/ligand_mpnn.py` |
| 19 | Load                      | Utilities            | `biopipelines/load.py` |
| 20 | MMseqs2                   | Utilities            | `biopipelines/mmseqs2.py` |
| 21 | MMseqs2Server             | Utilities            | `biopipelines/mmseqs2.py` |
| 22 | MSA                       | Utilities            | `biopipelines/msa.py` |
| 23 | Mutagenesis               | Sequence Design      | `biopipelines/mutagenesis.py` |
| 24 | MutationComposer          | Sequence Design      | `biopipelines/mutation_composer.py` |
| 25 | MutationProfiler          | Statistics           | `biopipelines/mutation_profiler.py` |
| 26 | Panda                     | Data Management      | `biopipelines/panda.py` |
| 27 | PDB                       | Utilities            | `biopipelines/pdb.py` |
| 28 | Plot                      | Utilities            | `biopipelines/plot.py` |
| 29 | Pool                      | Data Management      | `biopipelines/pool.py` |
| 30 | PoseBusters               | Analysis             | `biopipelines/posebusters.py` |
| 31 | PoseChange                | Analysis             | `biopipelines/pose_change.py` |
| 32 | ProteinMPNN               | Sequence Design      | `biopipelines/protein_mpnn.py` |
| 33 | PyMOL                     | Utilities            | `biopipelines/pymol.py` |
| 34 | RBSDesigner               | Sequence Design      | `biopipelines/rbs_designer.py` |
| 35 | RCSB                      | Utilities            | `biopipelines/rcsb.py` |
| 36 | ReMap                     | Data Management      | `biopipelines/remap.py` |
| 37 | RFdiffusion               | Structure Generation | `biopipelines/rfdiffusion.py` |
| 38 | RFdiffusion3              | Structure Generation | `biopipelines/rfdiffusion3.py` |
| 39 | RFdiffusionAllAtom        | Structure Generation | `biopipelines/rfdiffusion_allatom.py` |
| 40 | SASA                      | Analysis             | `biopipelines/sasa.py` |
| 41 | Selection                 | Data Management      | `biopipelines/selection.py` |
| 42 | Sequence                  | Utilities            | `biopipelines/sequence.py` |
| 43 | SequenceMetricCorrelation | Statistics           | `biopipelines/sequence_metric_correlation.py` |
| 44 | StitchSequences           | Sequence Design      | `biopipelines/stitch_sequences.py` |
| 45 | Table                     | Data Management      | `biopipelines/table.py` |

## Internal / auxiliary classes (not user-facing)

These appear as `TOOL_NAME = ...` for code-reuse but are **not** part of the
advertised public API:

| Tool | Role | Source file |
|------|------|-------------|
| Mock              | Test-only stub-output generator (used by the pytest suite) | `biopipelines/mock.py` |
| BoltzGenMerge     | Internal post-processing step of `BoltzGen` workflow        | `biopipelines/boltzgen.py` |
| BoltzGenImport    | Internal import step of `BoltzGen` workflow                  | `biopipelines/boltzgen.py` |
| RFDAA_PrepareLigand | Ligand preparation helper for `RFdiffusionAllAtom`        | `biopipelines/rfdiffusion_allatom.py` |
| `base`            | `BaseConfig` scaffolding marker                              | `biopipelines/base_config.py` |
| `install`         | Dynamic per-tool installation marker                         | `biopipelines/base_config.py` |
