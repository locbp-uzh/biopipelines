# BioPipelines Tool Index (SI)

Single-source-of-truth listing for every public-API tool. Generated from
`TOOL_NAME` / `TOOL_VERSION` in `biopipelines/*.py` and grouped by the categories
in [`tool_reference.md`](tool_reference.md).

**Public-API count: 75** (80 `TOOL_NAME` entries in the codebase, minus 4 internal
helpers — `BoltzGenMerge`, `BoltzGenImport`, `RFDAA_PrepareLigand`, `Mock` —
and 2 base-class scaffolding entries: `base`, `install`).

| # | Tool | Category | Version | Source file |
|---|------|----------|---------|-------------|
| 1 | BoltzGen                  | Structure Generation           | 1.0 | `biopipelines/boltzgen.py` |
| 2 | PocketGen                 | Structure Generation           | 1.0 | `biopipelines/pocketgen.py` |
| 3 | RFdiffusion               | Structure Generation           | 1.0 | `biopipelines/rfdiffusion.py` |
| 4 | RFdiffusion3              | Structure Generation           | 1.0 | `biopipelines/rfdiffusion3.py` |
| 5 | RFdiffusionAllAtom        | Structure Generation           | 1.0 | `biopipelines/rfdiffusion_allatom.py` |
| 6 | DNAEncoder                | Sequence Design                | 1.0 | `biopipelines/dna_encoder.py` |
| 7 | Frame2Seq                 | Sequence Design                | 1.0 | `biopipelines/frame2seq.py` |
| 8 | Fuse                      | Sequence Design                | 1.0 | `biopipelines/fuse.py` |
| 9 | LigandMPNN                | Sequence Design                | 1.0 | `biopipelines/ligand_mpnn.py` |
| 10 | Mutagenesis               | Sequence Design                | 1.0 | `biopipelines/mutagenesis.py` |
| 11 | MutationComposer          | Sequence Design                | 1.0 | `biopipelines/mutation_composer.py` |
| 12 | ProteinMPNN               | Sequence Design                | 1.0 | `biopipelines/protein_mpnn.py` |
| 13 | RBSDesigner               | Sequence Design                | 1.0 | `biopipelines/rbs_designer.py` |
| 14 | StitchSequences           | Sequence Design                | 1.0 | `biopipelines/stitch_sequences.py` |
| 15 | AlphaFold                 | Structure Prediction & Docking | 1.0 | `biopipelines/alphafold.py` |
| 16 | Boltz2                    | Structure Prediction & Docking | 1.0 | `biopipelines/boltz2.py` |
| 17 | DiffDock                  | Structure Prediction & Docking | 1.0 | `biopipelines/diffdock.py` |
| 18 | DynamicBind               | Structure Prediction & Docking | 1.0 | `biopipelines/dynamicbind.py` |
| 19 | ESMFold                   | Structure Prediction & Docking | 1.0 | `biopipelines/esmfold.py` |
| 20 | Gnina                     | Structure Prediction & Docking | 1.0 | `biopipelines/gnina.py` |
| 21 | NeuralPLexer              | Structure Prediction & Docking | 1.0 | `biopipelines/neuralplexer.py` |
| 22 | PLACER                    | Structure Prediction & Docking | 1.0 | `biopipelines/placer.py` |
| 23 | ADMETAI                   | Analysis                       | 1.0 | `biopipelines/admet_ai.py` |
| 24 | AF2BIND                   | Analysis                       | 1.0 | `biopipelines/af2bind.py` |
| 25 | Aggrescan3D               | Analysis                       | 1.0 | `biopipelines/aggrescan3d.py` |
| 26 | Angle                     | Analysis                       | 1.0 | `biopipelines/angle.py` |
| 27 | APBS                      | Analysis                       | 1.0 | `biopipelines/apbs.py` |
| 28 | BioEmu                    | Analysis                       | 1.0 | `biopipelines/bioemu.py` |
| 29 | CABSflex                  | Analysis                       | 1.0 | `biopipelines/cabsflex.py` |
| 30 | ConformationalChange      | Analysis                       | 1.0 | `biopipelines/conformational_change.py` |
| 31 | Consensus                 | Analysis                       | 1.0 | `biopipelines/consensus.py` |
| 32 | Contacts                  | Analysis                       | 1.0 | `biopipelines/contacts.py` |
| 33 | Distance                  | Analysis                       | 1.0 | `biopipelines/distance.py` |
| 34 | DistanceSelector          | Analysis                       | 1.0 | `biopipelines/distance_selector.py` |
| 35 | DSSP                      | Analysis                       | 1.0 | `biopipelines/dssp.py` |
| 36 | EnsembleAnalysis          | Analysis                       | 1.0 | `biopipelines/ensemble_analysis.py` |
| 37 | FPocket                   | Analysis                       | 1.0 | `biopipelines/fpocket.py` |
| 38 | GEMS                      | Analysis                       | 1.0 | `biopipelines/gems.py` |
| 39 | OpenMM                    | Analysis                       | 1.0 | `biopipelines/openmm.py` |
| 40 | P2Rank                    | Analysis                       | 1.0 | `biopipelines/p2rank.py` |
| 41 | PLIP                      | Analysis                       | 1.0 | `biopipelines/plip.py` |
| 42 | PLM_Sol                   | Analysis                       | 1.0 | `biopipelines/plm_sol.py` |
| 43 | PoseBusters               | Analysis                       | 1.0 | `biopipelines/posebusters.py` |
| 44 | PoseChange                | Analysis                       | 1.0 | `biopipelines/pose_change.py` |
| 45 | Prodigy                   | Analysis                       | 1.0 | `biopipelines/prodigy.py` |
| 46 | ProLIF                    | Analysis                       | 1.0 | `biopipelines/prolif.py` |
| 47 | Reduce                    | Analysis                       | 1.0 | `biopipelines/reduce.py` |
| 48 | RTMScore                  | Analysis                       | 1.0 | `biopipelines/rtmscore.py` |
| 49 | SASA                      | Analysis                       | 1.0 | `biopipelines/sasa.py` |
| 50 | ThermoMPNN                | Analysis                       | 1.0 | `biopipelines/thermompnn.py` |
| 51 | VespaG                    | Analysis                       | 1.0 | `biopipelines/vespag.py` |
| 52 | XTB                       | Analysis                       | 1.0 | `biopipelines/xtb.py` |
| 53 | OpenBabel                 | Cheminformatics                | 1.1 | `biopipelines/openbabel.py` |
| 54 | RDKit                     | Cheminformatics                | 1.0 | `biopipelines/rdkit_descriptors.py` |
| 55 | BayesianAdjuster          | Sequence Statistics            | 1.0 | `biopipelines/bayesian_adjuster.py` |
| 56 | MutationProfiler          | Sequence Statistics            | 1.0 | `biopipelines/mutation_profiler.py` |
| 57 | SequenceMetricCorrelation | Sequence Statistics            | 1.0 | `biopipelines/sequence_metric_correlation.py` |
| 58 | ExtractMetrics            | Data Management                | 1.0 | `biopipelines/extract_metrics.py` |
| 59 | Panda                     | Data Management                | 1.0 | `biopipelines/panda.py` |
| 60 | Pool                      | Data Management                | 1.0 | `biopipelines/pool.py` |
| 61 | ReMap                     | Data Management                | 1.0 | `biopipelines/remap.py` |
| 62 | Selection                 | Data Management                | 1.0 | `biopipelines/selection.py` |
| 63 | MMseqs2                   | MSAs                           | 1.0 | `biopipelines/mmseqs2.py` |
| 64 | MMseqs2Server             | MSAs                           | 1.0 | `biopipelines/mmseqs2.py` |
| 65 | MSA                       | MSAs                           | 1.0 | `biopipelines/msa.py` |
| 66 | CompoundLibrary           | Inputs & I/O                   | 1.0 | `biopipelines/compound_library.py` |
| 67 | Ligand                    | Inputs & I/O                   | 1.0 | `biopipelines/ligand.py` |
| 68 | Load                      | Inputs & I/O                   | 1.0 | `biopipelines/load.py` |
| 69 | PDB                       | Inputs & I/O                   | 1.0 | `biopipelines/pdb.py` |
| 70 | Plot                      | Inputs & I/O                   | 1.0 | `biopipelines/plot.py` |
| 71 | PyMOL                     | Inputs & I/O                   | 1.0 | `biopipelines/pymol.py` |
| 72 | RCSB                      | Inputs & I/O                   | 1.0 | `biopipelines/rcsb.py` |
| 73 | Sequence                  | Inputs & I/O                   | 1.0 | `biopipelines/sequence.py` |
| 74 | Table                     | Inputs & I/O                   | 1.0 | `biopipelines/table.py` |
| 75 | UniProt                   | Inputs & I/O                   | 1.0 | `biopipelines/uniprot.py` |

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
