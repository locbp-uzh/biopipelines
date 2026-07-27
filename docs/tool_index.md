# BioPipelines Tool Index (SI)

Single-source-of-truth listing for every public-API tool. Generated from
`TOOL_NAME` / `TOOL_VERSION` in `biopipelines/*.py` and grouped by the categories
in [`tool_reference.md`](tool_reference.md).

**Public-API count: 81** (listed below). Internal helpers — `BoltzGenMerge`,
`BoltzGenImport`, `RFDAA_PrepareLigand`, `Mock` — the `TemplateTool` skeleton,
and the 2 base-class scaffolding entries (`base`, `install`) are excluded.

| # | Tool | Category | Version | Source file |
|---|------|----------|---------|-------------|
| 1 | BoltzGen                  | Structure Generation           | 1.0 | `biopipelines/boltzgen.py` |
| 2 | PocketGen                 | Structure Generation           | 1.1 | `biopipelines/pocketgen.py` |
| 3 | RFdiffusion               | Structure Generation           | 1.1 | `biopipelines/rfdiffusion.py` |
| 4 | RFdiffusion3              | Structure Generation           | 2.0 | `biopipelines/rfdiffusion3.py` |
| 5 | RFdiffusionAllAtom        | Structure Generation           | 1.1 | `biopipelines/rfdiffusion_allatom.py` |
| 6 | RFdiffusion2              | Structure Generation           | 1.1 | `biopipelines/rfdiffusion2.py` |
| 7 | HBDesigner                | Structure Generation           | 1.1 | `biopipelines/hbdesigner.py` |
| 8 | DNAEncoder                | Sequence Design                | 1.0 | `biopipelines/dna_encoder.py` |
| 9 | Frame2Seq                 | Sequence Design                | 1.1 | `biopipelines/frame2seq.py` |
| 10 | Fuse                      | Sequence Design                | 1.0 | `biopipelines/fuse.py` |
| 11 | LigandMPNN                | Sequence Design                | 1.1 | `biopipelines/ligand_mpnn.py` |
| 12 | LASErMPNN                 | Sequence Design                | 1.0 | `biopipelines/lasermpnn.py` |
| 13 | Mutagenesis               | Sequence Design                | 1.0 | `biopipelines/mutagenesis.py` |
| 14 | MutationComposer          | Sequence Design                | 1.0 | `biopipelines/mutation_composer.py` |
| 15 | ProteinMPNN               | Sequence Design                | 1.1 | `biopipelines/protein_mpnn.py` |
| 16 | RBSDesigner               | Sequence Design                | 1.0 | `biopipelines/rbs_designer.py` |
| 17 | StitchSequences           | Sequence Design                | 1.0 | `biopipelines/stitch_sequences.py` |
| 18 | AlphaFold                 | Structure Prediction & Docking | 1.0 | `biopipelines/alphafold.py` |
| 19 | Boltz2                    | Structure Prediction & Docking | 1.0 | `biopipelines/boltz2.py` |
| 20 | DiffDock                  | Structure Prediction & Docking | 1.1 | `biopipelines/diffdock.py` |
| 21 | DynamicBind               | Structure Prediction & Docking | 1.1 | `biopipelines/dynamicbind.py` |
| 22 | ESMFold                   | Structure Prediction & Docking | 1.2 | `biopipelines/esmfold.py` |
| 23 | ESMFold2                  | Structure Prediction & Docking | 1.0 | `biopipelines/esmfold2.py` |
| 24 | Gnina                     | Structure Prediction & Docking | 1.1 | `biopipelines/gnina.py` |
| 25 | NeuralPLexer              | Structure Prediction & Docking | 1.1 | `biopipelines/neuralplexer.py` |
| 26 | PLACER                    | Structure Prediction & Docking | 1.1 | `biopipelines/placer.py` |
| 27 | ADMETAI                   | Analysis                       | 1.0 | `biopipelines/admet_ai.py` |
| 28 | AF2BIND                   | Analysis                       | 1.1 | `biopipelines/af2bind.py` |
| 29 | Aggrescan3D               | Analysis                       | 1.0 | `biopipelines/aggrescan3d.py` |
| 30 | Angle                     | Analysis                       | 1.0 | `biopipelines/angle.py` |
| 31 | APBS                      | Analysis                       | 1.1 | `biopipelines/apbs.py` |
| 32 | BioEmu                    | Analysis                       | 1.1 | `biopipelines/bioemu.py` |
| 33 | CABSflex                  | Analysis                       | 1.0 | `biopipelines/cabsflex.py` |
| 34 | ConformationalChange      | Analysis                       | 1.0 | `biopipelines/conformational_change.py` |
| 35 | Consensus                 | Analysis                       | 1.1 | `biopipelines/consensus.py` |
| 36 | Contacts                  | Analysis                       | 1.0 | `biopipelines/contacts.py` |
| 37 | Distance                  | Analysis                       | 1.0 | `biopipelines/distance.py` |
| 38 | DistanceSelector          | Analysis                       | 2.0 | `biopipelines/distance_selector.py` |
| 39 | DSSP                      | Analysis                       | 1.1 | `biopipelines/dssp.py` |
| 40 | EnsembleAnalysis          | Analysis                       | 1.0 | `biopipelines/ensemble_analysis.py` |
| 41 | FPocket                   | Analysis                       | 1.1 | `biopipelines/fpocket.py` |
| 42 | GEMS                      | Analysis                       | 1.1 | `biopipelines/gems.py` |
| 43 | OpenMM                    | Analysis                       | 1.1 | `biopipelines/openmm.py` |
| 44 | P2Rank                    | Analysis                       | 1.1 | `biopipelines/p2rank.py` |
| 45 | PLIP                      | Analysis                       | 1.0 | `biopipelines/plip.py` |
| 46 | PLM_Sol                   | Analysis                       | 1.2 | `biopipelines/plm_sol.py` |
| 47 | PoseBusters               | Analysis                       | 1.0 | `biopipelines/posebusters.py` |
| 48 | PoseChange                | Analysis                       | 1.0 | `biopipelines/pose_change.py` |
| 49 | Prodigy                   | Analysis                       | 1.1 | `biopipelines/prodigy.py` |
| 50 | ProLIF                    | Analysis                       | 1.1 | `biopipelines/prolif.py` |
| 51 | Reduce                    | Analysis                       | 1.1 | `biopipelines/reduce.py` |
| 52 | RTMScore                  | Analysis                       | 1.1 | `biopipelines/rtmscore.py` |
| 53 | SASA                      | Analysis                       | 1.0 | `biopipelines/sasa.py` |
| 54 | ThermoMPNN                | Analysis                       | 1.1 | `biopipelines/thermompnn.py` |
| 55 | VespaG                    | Analysis                       | 1.1 | `biopipelines/vespag.py` |
| 56 | XTB                       | Analysis                       | 1.1 | `biopipelines/xtb.py` |
| 57 | AiZynthFinder             | Cheminformatics                | 1.1 | `biopipelines/aizynthfinder.py` |
| 58 | OpenBabel                 | Cheminformatics                | 1.1 | `biopipelines/openbabel.py` |
| 59 | RDKit                     | Cheminformatics                | 1.1 | `biopipelines/rdkit_descriptors.py` |
| 60 | BayesianAdjuster          | Sequence Statistics            | 1.0 | `biopipelines/bayesian_adjuster.py` |
| 61 | MutationProfiler          | Sequence Statistics            | 1.0 | `biopipelines/mutation_profiler.py` |
| 62 | SequenceMetricCorrelation | Sequence Statistics            | 1.0 | `biopipelines/sequence_metric_correlation.py` |
| 63 | ExtractMetrics            | Data Management                | 1.0 | `biopipelines/extract_metrics.py` |
| 64 | Panda                     | Data Management                | 1.1 | `biopipelines/panda.py` |
| 65 | Pool                      | Data Management                | 1.0 | `biopipelines/pool.py` |
| 66 | ReMap                     | Data Management                | 1.0 | `biopipelines/remap.py` |
| 67 | Selection                 | Data Management                | 1.1 | `biopipelines/selection.py` |
| 68 | MMseqs2                   | MSAs                           | 1.0 | `biopipelines/mmseqs2.py` |
| 69 | MMseqs2Server             | MSAs                           | 1.0 | `biopipelines/mmseqs2.py` |
| 70 | MSA                       | MSAs                           | 1.0 | `biopipelines/msa.py` |
| 71 | CompoundLibrary           | Inputs & I/O                   | 1.0 | `biopipelines/compound_library.py` |
| 72 | Ligand                    | Inputs & I/O                   | 1.0 | `biopipelines/ligand.py` |
| 73 | Load                      | Inputs & I/O                   | 1.0 | `biopipelines/load.py` |
| 74 | PDB                       | Inputs & I/O                   | 1.0 | `biopipelines/pdb.py` |
| 75 | Plot                      | Inputs & I/O                   | 1.0 | `biopipelines/plot.py` |
| 76 | PyMOL                     | Inputs & I/O                   | 1.0 | `biopipelines/pymol.py` |
| 77 | Scripting                 | Inputs & I/O                   | 1.0 | `biopipelines/scripting.py` |
| 78 | RCSB                      | Inputs & I/O                   | 1.0 | `biopipelines/rcsb.py` |
| 79 | Sequence                  | Inputs & I/O                   | 1.0 | `biopipelines/sequence.py` |
| 80 | Table                     | Inputs & I/O                   | 1.0 | `biopipelines/table.py` |
| 81 | UniProt                   | Inputs & I/O                   | 1.0 | `biopipelines/uniprot.py` |

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
