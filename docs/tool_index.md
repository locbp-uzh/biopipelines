# BioPipelines Tool Index (SI)

Single-source-of-truth listing for every public-API tool. Generated from
`TOOL_NAME` / `TOOL_VERSION` in `biopipelines/*.py` and grouped by the categories
in [`tool_reference.md`](tool_reference.md).

**Public-API count: 106** (listed below). Internal helpers — `BoltzGenMerge`,
`BoltzGenImport`, `RFDAA_PrepareLigand`, `Mock` — the `TemplateTool` skeleton,
and the 2 base-class scaffolding entries (`base`, `install`) are excluded.

| # | Tool | Category | Version | Source file |
|---|------|----------|---------|-------------|
| 1 | BoltzGen | Structure Generation           | 2.6 | `biopipelines/boltzgen.py` |
| 2 | PocketGen | Structure Generation           | 2.0 | `biopipelines/pocketgen.py` |
| 3 | RFdiffusion | Structure Generation           | 2.4 | `biopipelines/rfdiffusion.py` |
| 4 | RFdiffusion3 | Structure Generation           | 3.8 | `biopipelines/rfdiffusion3.py` |
| 5 | RFD4Proteina | Structure Generation           | 1.3 | `biopipelines/rfd4proteina.py` |
| 6 | RFdiffusionAllAtom | Structure Generation           | 2.3 | `biopipelines/rfdiffusion_allatom.py` |
| 7 | RFdiffusion2              | Structure Generation           | 1.4 | `biopipelines/rfdiffusion2.py` |
| 8 | HBDesigner | Structure Generation           | 2.1 | `biopipelines/hbdesigner.py` |
| 9 | FLOWR | Structure Generation           | 1.2 | `biopipelines/flowr.py` |
| 10 | DNAEncoder                | Sequence Design                | 1.3 | `biopipelines/dna_encoder.py` |
| 11 | Frame2Seq | Sequence Design                | 2.0 | `biopipelines/frame2seq.py` |
| 12 | Fuse | Sequence Design                | 1.3 | `biopipelines/fuse.py` |
| 13 | LigandMPNN | Sequence Design                | 3.4 | `biopipelines/ligand_mpnn.py` |
| 14 | LASErMPNN | Sequence Design                | 2.2 | `biopipelines/lasermpnn.py` |
| 15 | Mutagenesis               | Sequence Design                | 1.3 | `biopipelines/mutagenesis.py` |
| 16 | MutationComposer          | Sequence Design                | 1.0 | `biopipelines/mutation_composer.py` |
| 17 | ProteinMPNN | Sequence Design                | 2.10 | `biopipelines/protein_mpnn.py` |
| 18 | RBSDesigner | Sequence Design                | 2.0 | `biopipelines/rbs_designer.py` |
| 19 | StitchSequences           | Sequence Design                | 1.0 | `biopipelines/stitch_sequences.py` |
| 20 | AlphaFold                 | Structure Prediction & Docking | 1.13 | `biopipelines/alphafold.py` |
| 21 | Boltz2 | Structure Prediction & Docking | 2.11 | `biopipelines/boltz2.py` |
| 22 | DiffDock | Structure Prediction & Docking | 2.1 | `biopipelines/diffdock.py` |
| 23 | DynamicBind | Structure Prediction & Docking | 2.2 | `biopipelines/dynamicbind.py` |
| 24 | ESMFold | Structure Prediction & Docking | 2.1 | `biopipelines/esmfold.py` |
| 25 | ESMFold2 | Structure Prediction & Docking | 2.8 | `biopipelines/esmfold2.py` |
| 26 | Gnina | Structure Prediction & Docking | 2.2 | `biopipelines/gnina.py` |
| 27 | Vina | Structure Prediction & Docking | 2.2 | `biopipelines/gnina.py` |
| 28 | NeuralPLexer | Structure Prediction & Docking | 2.1 | `biopipelines/neuralplexer.py` |
| 29 | OpenFold3 | Structure Prediction & Docking | 1.5 | `biopipelines/openfold3.py` |
| 30 | PLACER | Structure Prediction & Docking | 2.1 | `biopipelines/placer.py` |
| 31 | Protenix | Structure Prediction & Docking | 1.4 | `biopipelines/protenix.py` |
| 32 | ADMETAI | Analysis                       | 2.0 | `biopipelines/admet_ai.py` |
| 33 | AF2BIND | Analysis                       | 2.0 | `biopipelines/af2bind.py` |
| 34 | Aggrescan3D | Analysis                       | 2.2 | `biopipelines/aggrescan3d.py` |
| 35 | Angle                     | Analysis                       | 1.0 | `biopipelines/angle.py` |
| 36 | APBS | Analysis                       | 2.0 | `biopipelines/apbs.py` |
| 37 | BindingData               | Analysis                       | 1.2 | `biopipelines/binding_data.py` |
| 38 | BioEmu | Analysis                       | 2.0 | `biopipelines/bioemu.py` |
| 39 | BioEmu2 | Analysis                       | 1.4 | `biopipelines/bioemu2.py` |
| 40 | CABSflex | Analysis                       | 2.1 | `biopipelines/cabsflex.py` |
| 41 | ConformationalChange      | Analysis                       | 1.2 | `biopipelines/conformational_change.py` |
| 42 | Consensus                 | Analysis                       | 1.1 | `biopipelines/consensus.py` |
| 43 | Contacts                  | Analysis                       | 1.1 | `biopipelines/contacts.py` |
| 44 | Distance                  | Analysis                       | 1.0 | `biopipelines/distance.py` |
| 45 | DistanceSelector          | Analysis                       | 2.0 | `biopipelines/distance_selector.py` |
| 46 | DSSP | Analysis                       | 2.2 | `biopipelines/dssp.py` |
| 47 | EnsembleAnalysis          | Analysis                       | 1.0 | `biopipelines/ensemble_analysis.py` |
| 48 | Foldseek | Analysis                       | 1.4 | `biopipelines/foldseek.py` |
| 49 | FPocket | Analysis                       | 2.1 | `biopipelines/fpocket.py` |
| 50 | GEMS | Analysis                       | 2.0 | `biopipelines/gems.py` |
| 51 | IPSAE                     | Analysis                       | 1.4 | `biopipelines/ipsae.py` |
| 52 | LigandAtomSelector        | Analysis                       | 1.2 | `biopipelines/ligand_atom_selector.py` |
| 53 | LigandDdG                 | Analysis                       | 1.1 | `biopipelines/ligand_ddg.py` |
| 54 | OpenMM | Analysis                       | 3.3 | `biopipelines/openmm.py` |
| 55 | P2Rank | Analysis                       | 2.0 | `biopipelines/p2rank.py` |
| 56 | PPIDdG                    | Analysis                       | 1.0 | `biopipelines/ppi_ddg.py` |
| 57 | PLIP | Analysis                       | 2.2 | `biopipelines/plip.py` |
| 58 | PLM_Sol | Analysis                       | 2.0 | `biopipelines/plm_sol.py` |
| 59 | PoseBusters | Analysis                       | 2.2 | `biopipelines/posebusters.py` |
| 60 | PoseChange                | Analysis                       | 1.1 | `biopipelines/pose_change.py` |
| 61 | Peppr                     | Analysis                       | 1.2 | `biopipelines/peppr.py` |
| 62 | Prodigy | Analysis                       | 2.0 | `biopipelines/prodigy.py` |
| 63 | ProLIF | Analysis                       | 2.1 | `biopipelines/prolif.py` |
| 64 | Reduce | Analysis                       | 2.0 | `biopipelines/reduce.py` |
| 65 | RTMScore | Analysis                       | 2.1 | `biopipelines/rtmscore.py` |
| 66 | SASA                      | Analysis                       | 1.0 | `biopipelines/sasa.py` |
| 67 | SequenceIdentityGate      | Analysis                       | 1.0 | `biopipelines/sequence_identity_gate.py` |
| 68 | StructureCluster          | Analysis                       | 1.1 | `biopipelines/structure_cluster.py` |
| 69 | ThermoMPNN | Analysis                       | 2.0 | `biopipelines/thermompnn.py` |
| 70 | VespaG | Analysis                       | 2.0 | `biopipelines/vespag.py` |
| 71 | XTB | Analysis                       | 2.1 | `biopipelines/xtb.py` |
| 72 | TMalign                   | Analysis                       | 1.1 | `biopipelines/tmalign.py` |
| 73 | USalign                   | Analysis                       | 1.1 | `biopipelines/usalign.py` |
| 74 | BFactor                   | Analysis                       | 1.2 | `biopipelines/bfactor.py` |
| 75 | ESM2                      | Analysis                       | 1.1 | `biopipelines/esm2.py` |
| 76 | AiZynthFinder | Cheminformatics                | 2.5 | `biopipelines/aizynthfinder.py` |
| 77 | OpenBabel                 | Cheminformatics                | 1.1 | `biopipelines/openbabel.py` |
| 78 | RDKit                     | Cheminformatics                | 1.2 | `biopipelines/rdkit_descriptors.py` |
| 79 | REINVENT4                 | Cheminformatics                | 1.1 | `biopipelines/reinvent4.py` |
| 80 | BayesianAdjuster          | Sequence Statistics            | 1.0 | `biopipelines/bayesian_adjuster.py` |
| 81 | MutationProfiler | Sequence Statistics            | 2.0 | `biopipelines/mutation_profiler.py` |
| 82 | SequenceMetricCorrelation | Sequence Statistics            | 1.0 | `biopipelines/sequence_metric_correlation.py` |
| 83 | ExtractMetrics            | Data Management                | 1.1 | `biopipelines/extract_metrics.py` |
| 84 | Panda                     | Data Management                | 1.7 | `biopipelines/panda.py` |
| 85 | Pool                      | Data Management                | 1.3 | `biopipelines/pool.py` |
| 86 | ReMap | Data Management                | 1.2 | `biopipelines/remap.py` |
| 87 | Selection                 | Data Management                | 1.2 | `biopipelines/selection.py` |
| 88 | SelectionTransfer | Data Management | 1.0 | `biopipelines/selection_transfer.py` |
| 89 | MMseqs2                   | MSAs                           | 1.10 | `biopipelines/mmseqs2.py` |
| 90 | MMseqs2Server             | MSAs                           | 1.9 | `biopipelines/mmseqs2.py` |
| 91 | MSA                       | MSAs                           | 1.4 | `biopipelines/msa.py` |
| 92 | CompoundLibrary           | Inputs & I/O                   | 1.0 | `biopipelines/compound_library.py` |
| 93 | Ligand                    | Inputs & I/O                   | 1.5 | `biopipelines/ligand.py` |
| 94 | Load | Inputs & I/O                   | 1.7 | `biopipelines/load.py` |
| 95 | PDB | Inputs & I/O                   | 1.6 | `biopipelines/pdb.py` |
| 96 | Plot                      | Inputs & I/O                   | 1.0 | `biopipelines/plot.py` |
| 97 | PyMOL | Inputs & I/O                   | 2.1 | `biopipelines/pymol.py` |
| 98 | Scripting | Inputs & I/O                   | 1.4 | `biopipelines/scripting.py` |
| 99 | RCSB                      | Inputs & I/O                   | 1.5 | `biopipelines/rcsb.py` |
| 100 | Sequence                  | Inputs & I/O                   | 1.1 | `biopipelines/sequence.py` |
| 101 | Table                     | Inputs & I/O                   | 1.1 | `biopipelines/table.py` |
| 102 | UniProt | Inputs & I/O                   | 1.1 | `biopipelines/uniprot.py` |
| 103 | RF3 | Structure Prediction & Docking | 1.1 | `biopipelines/rf3.py` |
| 104 | FoundryMPNN | Sequence Design                | 1.1 | `biopipelines/foundry_mpnn.py` |
| 105 | Foldcomp | Inputs & I/O                   | 1.1 | `biopipelines/foldcomp.py` |
| 106 | Embedding                 | Analysis                       | 1.1 | `biopipelines/embedding.py` |

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
