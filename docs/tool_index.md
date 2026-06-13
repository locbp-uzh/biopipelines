# BioPipelines Tool Index (SI)

Single-source-of-truth listing for every public-API tool. Generated from
`TOOL_NAME` / `TOOL_VERSION` in `biopipelines/*.py` and grouped by the categories
in [`tool_reference.md`](tool_reference.md).

**Public-API count: 76** (81 `TOOL_NAME` entries in the codebase, minus 4 internal
helpers — `BoltzGenMerge`, `BoltzGenImport`, `RFDAA_PrepareLigand`, `Mock` —
and 2 base-class scaffolding entries: `base`, `install`).

| # | Tool | Category | Version | Source file |
|---|------|----------|---------|-------------|
| 1 | BoltzGen                  | Structure Generation           | 1.0 | `biopipelines/boltzgen.py` |
| 2 | PocketGen                 | Structure Generation           | 1.0 | `biopipelines/pocketgen.py` |
| 3 | RFdiffusion               | Structure Generation           | 1.0 | `biopipelines/rfdiffusion.py` |
| 4 | RFdiffusion3              | Structure Generation           | 1.0 | `biopipelines/rfdiffusion3.py` |
| 5 | RFdiffusionAllAtom        | Structure Generation           | 1.0 | `biopipelines/rfdiffusion_allatom.py` |
| 6 | HBDesigner                | Structure Generation           | 1.0 | `biopipelines/hbdesigner.py` |
| 7 | DNAEncoder                | Sequence Design                | 1.0 | `biopipelines/dna_encoder.py` |
| 8 | Frame2Seq                 | Sequence Design                | 1.0 | `biopipelines/frame2seq.py` |
| 9 | Fuse                      | Sequence Design                | 1.0 | `biopipelines/fuse.py` |
| 10 | LigandMPNN                | Sequence Design                | 1.0 | `biopipelines/ligand_mpnn.py` |
| 11 | Mutagenesis               | Sequence Design                | 1.0 | `biopipelines/mutagenesis.py` |
| 12 | MutationComposer          | Sequence Design                | 1.0 | `biopipelines/mutation_composer.py` |
| 13 | ProteinMPNN               | Sequence Design                | 1.0 | `biopipelines/protein_mpnn.py` |
| 14 | RBSDesigner               | Sequence Design                | 1.0 | `biopipelines/rbs_designer.py` |
| 15 | StitchSequences           | Sequence Design                | 1.0 | `biopipelines/stitch_sequences.py` |
| 16 | AlphaFold                 | Structure Prediction & Docking | 1.0 | `biopipelines/alphafold.py` |
| 17 | Boltz2                    | Structure Prediction & Docking | 1.0 | `biopipelines/boltz2.py` |
| 18 | DiffDock                  | Structure Prediction & Docking | 1.0 | `biopipelines/diffdock.py` |
| 19 | DynamicBind               | Structure Prediction & Docking | 1.0 | `biopipelines/dynamicbind.py` |
| 20 | ESMFold                   | Structure Prediction & Docking | 1.0 | `biopipelines/esmfold.py` |
| 21 | Gnina                     | Structure Prediction & Docking | 1.0 | `biopipelines/gnina.py` |
| 22 | NeuralPLexer              | Structure Prediction & Docking | 1.0 | `biopipelines/neuralplexer.py` |
| 23 | PLACER                    | Structure Prediction & Docking | 1.0 | `biopipelines/placer.py` |
| 24 | ADMETAI                   | Analysis                       | 1.0 | `biopipelines/admet_ai.py` |
| 25 | AF2BIND                   | Analysis                       | 1.0 | `biopipelines/af2bind.py` |
| 26 | Aggrescan3D               | Analysis                       | 1.0 | `biopipelines/aggrescan3d.py` |
| 27 | Angle                     | Analysis                       | 1.0 | `biopipelines/angle.py` |
| 28 | APBS                      | Analysis                       | 1.0 | `biopipelines/apbs.py` |
| 29 | BioEmu                    | Analysis                       | 1.0 | `biopipelines/bioemu.py` |
| 30 | CABSflex                  | Analysis                       | 1.0 | `biopipelines/cabsflex.py` |
| 31 | ConformationalChange      | Analysis                       | 1.0 | `biopipelines/conformational_change.py` |
| 32 | Consensus                 | Analysis                       | 1.0 | `biopipelines/consensus.py` |
| 33 | Contacts                  | Analysis                       | 1.0 | `biopipelines/contacts.py` |
| 34 | Distance                  | Analysis                       | 1.0 | `biopipelines/distance.py` |
| 35 | DistanceSelector          | Analysis                       | 1.0 | `biopipelines/distance_selector.py` |
| 36 | DSSP                      | Analysis                       | 1.0 | `biopipelines/dssp.py` |
| 37 | EnsembleAnalysis          | Analysis                       | 1.0 | `biopipelines/ensemble_analysis.py` |
| 38 | FPocket                   | Analysis                       | 1.0 | `biopipelines/fpocket.py` |
| 39 | GEMS                      | Analysis                       | 1.0 | `biopipelines/gems.py` |
| 40 | OpenMM                    | Analysis                       | 1.0 | `biopipelines/openmm.py` |
| 41 | P2Rank                    | Analysis                       | 1.0 | `biopipelines/p2rank.py` |
| 42 | PLIP                      | Analysis                       | 1.0 | `biopipelines/plip.py` |
| 43 | PLM_Sol                   | Analysis                       | 1.0 | `biopipelines/plm_sol.py` |
| 44 | PoseBusters               | Analysis                       | 1.0 | `biopipelines/posebusters.py` |
| 45 | PoseChange                | Analysis                       | 1.0 | `biopipelines/pose_change.py` |
| 46 | Prodigy                   | Analysis                       | 1.0 | `biopipelines/prodigy.py` |
| 47 | ProLIF                    | Analysis                       | 1.0 | `biopipelines/prolif.py` |
| 48 | Reduce                    | Analysis                       | 1.0 | `biopipelines/reduce.py` |
| 49 | RTMScore                  | Analysis                       | 1.0 | `biopipelines/rtmscore.py` |
| 50 | SASA                      | Analysis                       | 1.0 | `biopipelines/sasa.py` |
| 51 | ThermoMPNN                | Analysis                       | 1.0 | `biopipelines/thermompnn.py` |
| 52 | VespaG                    | Analysis                       | 1.0 | `biopipelines/vespag.py` |
| 53 | XTB                       | Analysis                       | 1.0 | `biopipelines/xtb.py` |
| 54 | OpenBabel                 | Cheminformatics                | 1.1 | `biopipelines/openbabel.py` |
| 55 | RDKit                     | Cheminformatics                | 1.0 | `biopipelines/rdkit_descriptors.py` |
| 56 | BayesianAdjuster          | Sequence Statistics            | 1.0 | `biopipelines/bayesian_adjuster.py` |
| 57 | MutationProfiler          | Sequence Statistics            | 1.0 | `biopipelines/mutation_profiler.py` |
| 58 | SequenceMetricCorrelation | Sequence Statistics            | 1.0 | `biopipelines/sequence_metric_correlation.py` |
| 59 | ExtractMetrics            | Data Management                | 1.0 | `biopipelines/extract_metrics.py` |
| 60 | Panda                     | Data Management                | 1.0 | `biopipelines/panda.py` |
| 61 | Pool                      | Data Management                | 1.0 | `biopipelines/pool.py` |
| 62 | ReMap                     | Data Management                | 1.0 | `biopipelines/remap.py` |
| 63 | Selection                 | Data Management                | 1.0 | `biopipelines/selection.py` |
| 64 | MMseqs2                   | MSAs                           | 1.0 | `biopipelines/mmseqs2.py` |
| 65 | MMseqs2Server             | MSAs                           | 1.0 | `biopipelines/mmseqs2.py` |
| 66 | MSA                       | MSAs                           | 1.0 | `biopipelines/msa.py` |
| 67 | CompoundLibrary           | Inputs & I/O                   | 1.0 | `biopipelines/compound_library.py` |
| 68 | Ligand                    | Inputs & I/O                   | 1.0 | `biopipelines/ligand.py` |
| 69 | Load                      | Inputs & I/O                   | 1.0 | `biopipelines/load.py` |
| 70 | PDB                       | Inputs & I/O                   | 1.0 | `biopipelines/pdb.py` |
| 71 | Plot                      | Inputs & I/O                   | 1.0 | `biopipelines/plot.py` |
| 72 | PyMOL                     | Inputs & I/O                   | 1.0 | `biopipelines/pymol.py` |
| 73 | RCSB                      | Inputs & I/O                   | 1.0 | `biopipelines/rcsb.py` |
| 74 | Sequence                  | Inputs & I/O                   | 1.0 | `biopipelines/sequence.py` |
| 75 | Table                     | Inputs & I/O                   | 1.0 | `biopipelines/table.py` |
| 76 | UniProt                   | Inputs & I/O                   | 1.0 | `biopipelines/uniprot.py` |

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
