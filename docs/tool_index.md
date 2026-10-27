# BioPipelines Tool Index (SI)

Single-source-of-truth listing for every public-API tool. Generated from
`TOOL_NAME` / `TOOL_VERSION` in `biopipelines/*.py` and grouped by the categories
in [`tool_reference.md`](tool_reference.md).

**Public-API count: 96** (listed below). Internal helpers — `BoltzGenMerge`,
`BoltzGenImport`, `RFDAA_PrepareLigand`, `Mock` — the `TemplateTool` skeleton,
and the 2 base-class scaffolding entries (`base`, `install`) are excluded.

| # | Tool | Category | Version | Source file |
|---|------|----------|---------|-------------|
| 1 | BoltzGen | Structure Generation           | 2.6 | `biopipelines/boltzgen.py` |
| 2 | PocketGen | Structure Generation           | 2.0 | `biopipelines/pocketgen.py` |
| 3 | RFdiffusion | Structure Generation           | 2.4 | `biopipelines/rfdiffusion.py` |
| 4 | RFdiffusion3 | Structure Generation           | 3.7 | `biopipelines/rfdiffusion3.py` |
| 5 | RFD4Proteina | Structure Generation           | 1.0 | `biopipelines/rfd4proteina.py` |
| 6 | RFdiffusionAllAtom | Structure Generation           | 2.3 | `biopipelines/rfdiffusion_allatom.py` |
| 7 | RFdiffusion2              | Structure Generation           | 1.4 | `biopipelines/rfdiffusion2.py` |
| 8 | HBDesigner | Structure Generation           | 2.1 | `biopipelines/hbdesigner.py` |
| 9 | DNAEncoder                | Sequence Design                | 1.3 | `biopipelines/dna_encoder.py` |
| 10 | Frame2Seq | Sequence Design                | 2.0 | `biopipelines/frame2seq.py` |
| 11 | Fuse | Sequence Design                | 1.3 | `biopipelines/fuse.py` |
| 12 | LigandMPNN | Sequence Design                | 3.1 | `biopipelines/ligand_mpnn.py` |
| 13 | LASErMPNN | Sequence Design                | 2.2 | `biopipelines/lasermpnn.py` |
| 14 | Mutagenesis               | Sequence Design                | 1.3 | `biopipelines/mutagenesis.py` |
| 15 | MutationComposer          | Sequence Design                | 1.0 | `biopipelines/mutation_composer.py` |
| 16 | ProteinMPNN | Sequence Design                | 2.10 | `biopipelines/protein_mpnn.py` |
| 17 | RBSDesigner | Sequence Design                | 2.0 | `biopipelines/rbs_designer.py` |
| 18 | StitchSequences           | Sequence Design                | 1.0 | `biopipelines/stitch_sequences.py` |
| 19 | AlphaFold                 | Structure Prediction & Docking | 1.13 | `biopipelines/alphafold.py` |
| 20 | Boltz2 | Structure Prediction & Docking | 2.10 | `biopipelines/boltz2.py` |
| 21 | DiffDock | Structure Prediction & Docking | 2.1 | `biopipelines/diffdock.py` |
| 22 | DynamicBind | Structure Prediction & Docking | 2.2 | `biopipelines/dynamicbind.py` |
| 23 | ESMFold | Structure Prediction & Docking | 2.1 | `biopipelines/esmfold.py` |
| 24 | ESMFold2 | Structure Prediction & Docking | 2.6 | `biopipelines/esmfold2.py` |
| 25 | Gnina | Structure Prediction & Docking | 2.2 | `biopipelines/gnina.py` |
| 26 | Vina | Structure Prediction & Docking | 2.2 | `biopipelines/gnina.py` |
| 27 | NeuralPLexer | Structure Prediction & Docking | 2.1 | `biopipelines/neuralplexer.py` |
| 28 | OpenFold3 | Structure Prediction & Docking | 1.5 | `biopipelines/openfold3.py` |
| 29 | PLACER | Structure Prediction & Docking | 2.1 | `biopipelines/placer.py` |
| 30 | ADMETAI | Analysis                       | 2.0 | `biopipelines/admet_ai.py` |
| 31 | AF2BIND | Analysis                       | 2.0 | `biopipelines/af2bind.py` |
| 32 | Aggrescan3D | Analysis                       | 2.2 | `biopipelines/aggrescan3d.py` |
| 33 | Angle                     | Analysis                       | 1.0 | `biopipelines/angle.py` |
| 34 | APBS | Analysis                       | 2.0 | `biopipelines/apbs.py` |
| 35 | BindingData               | Analysis                       | 1.2 | `biopipelines/binding_data.py` |
| 36 | BioEmu | Analysis                       | 2.0 | `biopipelines/bioemu.py` |
| 37 | BioEmu2 | Analysis                       | 1.0 | `biopipelines/bioemu2.py` |
| 38 | CABSflex | Analysis                       | 2.1 | `biopipelines/cabsflex.py` |
| 39 | ConformationalChange      | Analysis                       | 1.2 | `biopipelines/conformational_change.py` |
| 40 | Consensus                 | Analysis                       | 1.1 | `biopipelines/consensus.py` |
| 41 | Contacts                  | Analysis                       | 1.1 | `biopipelines/contacts.py` |
| 42 | Distance                  | Analysis                       | 1.0 | `biopipelines/distance.py` |
| 43 | DistanceSelector          | Analysis                       | 2.0 | `biopipelines/distance_selector.py` |
| 44 | DSSP | Analysis                       | 2.2 | `biopipelines/dssp.py` |
| 45 | EnsembleAnalysis          | Analysis                       | 1.0 | `biopipelines/ensemble_analysis.py` |
| 46 | Foldseek | Analysis                       | 1.0 | `biopipelines/foldseek.py` |
| 47 | FPocket | Analysis                       | 2.1 | `biopipelines/fpocket.py` |
| 48 | GEMS | Analysis                       | 2.0 | `biopipelines/gems.py` |
| 49 | IPSAE                     | Analysis                       | 1.0 | `biopipelines/ipsae.py` |
| 50 | LigandAtomSelector        | Analysis                       | 1.2 | `biopipelines/ligand_atom_selector.py` |
| 51 | LigandDdG                 | Analysis                       | 1.0 | `biopipelines/ligand_ddg.py` |
| 52 | OpenMM | Analysis                       | 3.3 | `biopipelines/openmm.py` |
| 53 | P2Rank | Analysis                       | 2.0 | `biopipelines/p2rank.py` |
| 54 | PPIDdG                    | Analysis                       | 1.0 | `biopipelines/ppi_ddg.py` |
| 55 | PLIP | Analysis                       | 2.2 | `biopipelines/plip.py` |
| 56 | PLM_Sol | Analysis                       | 2.0 | `biopipelines/plm_sol.py` |
| 57 | PoseBusters | Analysis                       | 2.2 | `biopipelines/posebusters.py` |
| 58 | PoseChange                | Analysis                       | 1.1 | `biopipelines/pose_change.py` |
| 59 | Peppr                     | Analysis                       | 1.0 | `biopipelines/peppr.py` |
| 60 | Prodigy | Analysis                       | 2.0 | `biopipelines/prodigy.py` |
| 61 | ProLIF | Analysis                       | 2.1 | `biopipelines/prolif.py` |
| 62 | Reduce | Analysis                       | 2.0 | `biopipelines/reduce.py` |
| 63 | RTMScore | Analysis                       | 2.1 | `biopipelines/rtmscore.py` |
| 64 | SASA                      | Analysis                       | 1.0 | `biopipelines/sasa.py` |
| 65 | SequenceIdentityGate      | Analysis                       | 1.0 | `biopipelines/sequence_identity_gate.py` |
| 66 | StructureCluster          | Analysis                       | 1.1 | `biopipelines/structure_cluster.py` |
| 67 | ThermoMPNN | Analysis                       | 2.0 | `biopipelines/thermompnn.py` |
| 68 | VespaG | Analysis                       | 2.0 | `biopipelines/vespag.py` |
| 69 | XTB | Analysis                       | 2.1 | `biopipelines/xtb.py` |
| 70 | BFactor                   | Analysis                       | 1.2 | `biopipelines/bfactor.py` |
| 71 | AiZynthFinder | Cheminformatics                | 2.2 | `biopipelines/aizynthfinder.py` |
| 72 | OpenBabel                 | Cheminformatics                | 1.1 | `biopipelines/openbabel.py` |
| 73 | RDKit                     | Cheminformatics                | 1.2 | `biopipelines/rdkit_descriptors.py` |
| 74 | BayesianAdjuster          | Sequence Statistics            | 1.0 | `biopipelines/bayesian_adjuster.py` |
| 75 | MutationProfiler | Sequence Statistics            | 2.0 | `biopipelines/mutation_profiler.py` |
| 76 | SequenceMetricCorrelation | Sequence Statistics            | 1.0 | `biopipelines/sequence_metric_correlation.py` |
| 77 | ExtractMetrics            | Data Management                | 1.1 | `biopipelines/extract_metrics.py` |
| 78 | Panda                     | Data Management                | 1.5 | `biopipelines/panda.py` |
| 79 | Pool                      | Data Management                | 1.2 | `biopipelines/pool.py` |
| 80 | ReMap | Data Management                | 1.2 | `biopipelines/remap.py` |
| 81 | Selection                 | Data Management                | 1.2 | `biopipelines/selection.py` |
| 82 | SelectionTransfer | Data Management | 1.0 | `biopipelines/selection_transfer.py` |
| 83 | MMseqs2                   | MSAs                           | 1.10 | `biopipelines/mmseqs2.py` |
| 84 | MMseqs2Server             | MSAs                           | 1.9 | `biopipelines/mmseqs2.py` |
| 85 | MSA                       | MSAs                           | 1.4 | `biopipelines/msa.py` |
| 86 | CompoundLibrary           | Inputs & I/O                   | 1.0 | `biopipelines/compound_library.py` |
| 87 | Ligand                    | Inputs & I/O                   | 1.5 | `biopipelines/ligand.py` |
| 88 | Load | Inputs & I/O                   | 1.6 | `biopipelines/load.py` |
| 89 | PDB | Inputs & I/O                   | 1.6 | `biopipelines/pdb.py` |
| 90 | Plot                      | Inputs & I/O                   | 1.0 | `biopipelines/plot.py` |
| 91 | PyMOL | Inputs & I/O                   | 2.1 | `biopipelines/pymol.py` |
| 92 | Scripting | Inputs & I/O                   | 1.3 | `biopipelines/scripting.py` |
| 93 | RCSB                      | Inputs & I/O                   | 1.5 | `biopipelines/rcsb.py` |
| 94 | Sequence                  | Inputs & I/O                   | 1.0 | `biopipelines/sequence.py` |
| 95 | Table                     | Inputs & I/O                   | 1.1 | `biopipelines/table.py` |
| 96 | UniProt | Inputs & I/O                   | 1.1 | `biopipelines/uniprot.py` |

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
