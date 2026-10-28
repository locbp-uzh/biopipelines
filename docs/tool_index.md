# BioPipelines Tool Index (SI)

Single-source-of-truth listing for every public-API tool. Generated from
`TOOL_NAME` / `TOOL_VERSION` in `biopipelines/*.py` and grouped by the categories
in [`tool_reference.md`](tool_reference.md).

**Public-API count: 101** (listed below). Internal helpers — `BoltzGenMerge`,
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
| 12 | LigandMPNN | Sequence Design                | 3.2 | `biopipelines/ligand_mpnn.py` |
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
| 30 | Protenix | Structure Prediction & Docking | 1.0 | `biopipelines/protenix.py` |
| 31 | ADMETAI | Analysis                       | 2.0 | `biopipelines/admet_ai.py` |
| 32 | AF2BIND | Analysis                       | 2.0 | `biopipelines/af2bind.py` |
| 33 | Aggrescan3D | Analysis                       | 2.2 | `biopipelines/aggrescan3d.py` |
| 34 | Angle                     | Analysis                       | 1.0 | `biopipelines/angle.py` |
| 35 | APBS | Analysis                       | 2.0 | `biopipelines/apbs.py` |
| 36 | BindingData               | Analysis                       | 1.2 | `biopipelines/binding_data.py` |
| 37 | BioEmu | Analysis                       | 2.0 | `biopipelines/bioemu.py` |
| 38 | BioEmu2 | Analysis                       | 1.0 | `biopipelines/bioemu2.py` |
| 39 | CABSflex | Analysis                       | 2.1 | `biopipelines/cabsflex.py` |
| 40 | ConformationalChange      | Analysis                       | 1.2 | `biopipelines/conformational_change.py` |
| 41 | Consensus                 | Analysis                       | 1.1 | `biopipelines/consensus.py` |
| 42 | Contacts                  | Analysis                       | 1.1 | `biopipelines/contacts.py` |
| 43 | Distance                  | Analysis                       | 1.0 | `biopipelines/distance.py` |
| 44 | DistanceSelector          | Analysis                       | 2.0 | `biopipelines/distance_selector.py` |
| 45 | DSSP | Analysis                       | 2.2 | `biopipelines/dssp.py` |
| 46 | EnsembleAnalysis          | Analysis                       | 1.0 | `biopipelines/ensemble_analysis.py` |
| 47 | Foldseek | Analysis                       | 1.0 | `biopipelines/foldseek.py` |
| 48 | FPocket | Analysis                       | 2.1 | `biopipelines/fpocket.py` |
| 49 | GEMS | Analysis                       | 2.0 | `biopipelines/gems.py` |
| 50 | IPSAE                     | Analysis                       | 1.0 | `biopipelines/ipsae.py` |
| 51 | LigandAtomSelector        | Analysis                       | 1.2 | `biopipelines/ligand_atom_selector.py` |
| 52 | LigandDdG                 | Analysis                       | 1.0 | `biopipelines/ligand_ddg.py` |
| 53 | OpenMM | Analysis                       | 3.3 | `biopipelines/openmm.py` |
| 54 | P2Rank | Analysis                       | 2.0 | `biopipelines/p2rank.py` |
| 55 | PPIDdG                    | Analysis                       | 1.0 | `biopipelines/ppi_ddg.py` |
| 56 | PLIP | Analysis                       | 2.2 | `biopipelines/plip.py` |
| 57 | PLM_Sol | Analysis                       | 2.0 | `biopipelines/plm_sol.py` |
| 58 | PoseBusters | Analysis                       | 2.2 | `biopipelines/posebusters.py` |
| 59 | PoseChange                | Analysis                       | 1.1 | `biopipelines/pose_change.py` |
| 60 | Peppr                     | Analysis                       | 1.0 | `biopipelines/peppr.py` |
| 61 | Prodigy | Analysis                       | 2.0 | `biopipelines/prodigy.py` |
| 62 | ProLIF | Analysis                       | 2.1 | `biopipelines/prolif.py` |
| 63 | Reduce | Analysis                       | 2.0 | `biopipelines/reduce.py` |
| 64 | RTMScore | Analysis                       | 2.1 | `biopipelines/rtmscore.py` |
| 65 | SASA                      | Analysis                       | 1.0 | `biopipelines/sasa.py` |
| 66 | SequenceIdentityGate      | Analysis                       | 1.0 | `biopipelines/sequence_identity_gate.py` |
| 67 | StructureCluster          | Analysis                       | 1.1 | `biopipelines/structure_cluster.py` |
| 68 | ThermoMPNN | Analysis                       | 2.0 | `biopipelines/thermompnn.py` |
| 69 | VespaG | Analysis                       | 2.0 | `biopipelines/vespag.py` |
| 70 | XTB | Analysis                       | 2.1 | `biopipelines/xtb.py` |
| 71 | BFactor                   | Analysis                       | 1.2 | `biopipelines/bfactor.py` |
| 72 | ESM2                      | Analysis                       | 1.0 | `biopipelines/esm2.py` |
| 73 | AiZynthFinder | Cheminformatics                | 2.2 | `biopipelines/aizynthfinder.py` |
| 74 | OpenBabel                 | Cheminformatics                | 1.1 | `biopipelines/openbabel.py` |
| 75 | RDKit                     | Cheminformatics                | 1.2 | `biopipelines/rdkit_descriptors.py` |
| 76 | REINVENT4                 | Cheminformatics                | 1.0 | `biopipelines/reinvent4.py` |
| 77 | BayesianAdjuster          | Sequence Statistics            | 1.0 | `biopipelines/bayesian_adjuster.py` |
| 78 | MutationProfiler | Sequence Statistics            | 2.0 | `biopipelines/mutation_profiler.py` |
| 79 | SequenceMetricCorrelation | Sequence Statistics            | 1.0 | `biopipelines/sequence_metric_correlation.py` |
| 80 | ExtractMetrics            | Data Management                | 1.1 | `biopipelines/extract_metrics.py` |
| 81 | Panda                     | Data Management                | 1.5 | `biopipelines/panda.py` |
| 82 | Pool                      | Data Management                | 1.3 | `biopipelines/pool.py` |
| 83 | ReMap | Data Management                | 1.2 | `biopipelines/remap.py` |
| 84 | Selection                 | Data Management                | 1.2 | `biopipelines/selection.py` |
| 85 | SelectionTransfer | Data Management | 1.0 | `biopipelines/selection_transfer.py` |
| 86 | MMseqs2                   | MSAs                           | 1.10 | `biopipelines/mmseqs2.py` |
| 87 | MMseqs2Server             | MSAs                           | 1.9 | `biopipelines/mmseqs2.py` |
| 88 | MSA                       | MSAs                           | 1.4 | `biopipelines/msa.py` |
| 89 | CompoundLibrary           | Inputs & I/O                   | 1.0 | `biopipelines/compound_library.py` |
| 90 | Ligand                    | Inputs & I/O                   | 1.5 | `biopipelines/ligand.py` |
| 91 | Load | Inputs & I/O                   | 1.6 | `biopipelines/load.py` |
| 92 | PDB | Inputs & I/O                   | 1.6 | `biopipelines/pdb.py` |
| 93 | Plot                      | Inputs & I/O                   | 1.0 | `biopipelines/plot.py` |
| 94 | PyMOL | Inputs & I/O                   | 2.1 | `biopipelines/pymol.py` |
| 95 | Scripting | Inputs & I/O                   | 1.4 | `biopipelines/scripting.py` |
| 96 | RCSB                      | Inputs & I/O                   | 1.5 | `biopipelines/rcsb.py` |
| 97 | Sequence                  | Inputs & I/O                   | 1.0 | `biopipelines/sequence.py` |
| 98 | Table                     | Inputs & I/O                   | 1.1 | `biopipelines/table.py` |
| 99 | UniProt | Inputs & I/O                   | 1.1 | `biopipelines/uniprot.py` |
| 100 | RF3 | Structure Prediction & Docking | 1.0 | `biopipelines/rf3.py` |
| 101 | FoundryMPNN | Sequence Design                | 1.0 | `biopipelines/foundry_mpnn.py` |

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
