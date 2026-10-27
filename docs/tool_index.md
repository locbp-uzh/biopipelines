# BioPipelines Tool Index (SI)

Single-source-of-truth listing for every public-API tool. Generated from
`TOOL_NAME` / `TOOL_VERSION` in `biopipelines/*.py` and grouped by the categories
in [`tool_reference.md`](tool_reference.md).

**Public-API count: 88** (listed below). Internal helpers — `BoltzGenMerge`,
`BoltzGenImport`, `RFDAA_PrepareLigand`, `Mock` — the `TemplateTool` skeleton,
and the 2 base-class scaffolding entries (`base`, `install`) are excluded.

| # | Tool | Category | Version | Source file |
|---|------|----------|---------|-------------|
| 1 | BoltzGen | Structure Generation           | 2.6 | `biopipelines/boltzgen.py` |
| 2 | PocketGen | Structure Generation           | 2.0 | `biopipelines/pocketgen.py` |
| 3 | RFdiffusion | Structure Generation           | 2.4 | `biopipelines/rfdiffusion.py` |
| 4 | RFdiffusion3 | Structure Generation           | 3.7 | `biopipelines/rfdiffusion3.py` |
| 5 | RFdiffusionAllAtom | Structure Generation           | 2.3 | `biopipelines/rfdiffusion_allatom.py` |
| 6 | RFdiffusion2              | Structure Generation           | 1.4 | `biopipelines/rfdiffusion2.py` |
| 7 | HBDesigner | Structure Generation           | 2.1 | `biopipelines/hbdesigner.py` |
| 8 | DNAEncoder                | Sequence Design                | 1.3 | `biopipelines/dna_encoder.py` |
| 9 | Frame2Seq | Sequence Design                | 2.0 | `biopipelines/frame2seq.py` |
| 10 | Fuse | Sequence Design                | 1.3 | `biopipelines/fuse.py` |
| 11 | LigandMPNN | Sequence Design                | 3.0 | `biopipelines/ligand_mpnn.py` |
| 12 | LASErMPNN | Sequence Design                | 2.2 | `biopipelines/lasermpnn.py` |
| 13 | Mutagenesis               | Sequence Design                | 1.3 | `biopipelines/mutagenesis.py` |
| 14 | MutationComposer          | Sequence Design                | 1.0 | `biopipelines/mutation_composer.py` |
| 15 | ProteinMPNN | Sequence Design                | 2.10 | `biopipelines/protein_mpnn.py` |
| 16 | RBSDesigner | Sequence Design                | 2.0 | `biopipelines/rbs_designer.py` |
| 17 | StitchSequences           | Sequence Design                | 1.0 | `biopipelines/stitch_sequences.py` |
| 18 | AlphaFold                 | Structure Prediction & Docking | 1.11 | `biopipelines/alphafold.py` |
| 19 | Boltz2 | Structure Prediction & Docking | 2.9 | `biopipelines/boltz2.py` |
| 20 | DiffDock | Structure Prediction & Docking | 2.1 | `biopipelines/diffdock.py` |
| 21 | DynamicBind | Structure Prediction & Docking | 2.2 | `biopipelines/dynamicbind.py` |
| 22 | ESMFold | Structure Prediction & Docking | 2.1 | `biopipelines/esmfold.py` |
| 23 | ESMFold2 | Structure Prediction & Docking | 2.5 | `biopipelines/esmfold2.py` |
| 24 | Gnina | Structure Prediction & Docking | 2.2 | `biopipelines/gnina.py` |
| 25 | Vina | Structure Prediction & Docking | 2.2 | `biopipelines/gnina.py` |
| 26 | NeuralPLexer | Structure Prediction & Docking | 2.1 | `biopipelines/neuralplexer.py` |
| 27 | OpenFold3 | Structure Prediction & Docking | 1.4 | `biopipelines/openfold3.py` |
| 28 | PLACER | Structure Prediction & Docking | 2.1 | `biopipelines/placer.py` |
| 29 | ADMETAI | Analysis                       | 2.0 | `biopipelines/admet_ai.py` |
| 30 | AF2BIND | Analysis                       | 2.0 | `biopipelines/af2bind.py` |
| 31 | Aggrescan3D | Analysis                       | 2.2 | `biopipelines/aggrescan3d.py` |
| 32 | Angle                     | Analysis                       | 1.0 | `biopipelines/angle.py` |
| 33 | APBS | Analysis                       | 2.0 | `biopipelines/apbs.py` |
| 34 | BindingData               | Analysis                       | 1.2 | `biopipelines/binding_data.py` |
| 35 | BioEmu | Analysis                       | 2.0 | `biopipelines/bioemu.py` |
| 36 | CABSflex | Analysis                       | 2.1 | `biopipelines/cabsflex.py` |
| 37 | ConformationalChange      | Analysis                       | 1.2 | `biopipelines/conformational_change.py` |
| 38 | Consensus                 | Analysis                       | 1.1 | `biopipelines/consensus.py` |
| 39 | Contacts                  | Analysis                       | 1.1 | `biopipelines/contacts.py` |
| 40 | Distance                  | Analysis                       | 1.0 | `biopipelines/distance.py` |
| 41 | DistanceSelector          | Analysis                       | 2.0 | `biopipelines/distance_selector.py` |
| 42 | DSSP | Analysis                       | 2.2 | `biopipelines/dssp.py` |
| 43 | EnsembleAnalysis          | Analysis                       | 1.0 | `biopipelines/ensemble_analysis.py` |
| 44 | FPocket | Analysis                       | 2.1 | `biopipelines/fpocket.py` |
| 45 | GEMS | Analysis                       | 2.0 | `biopipelines/gems.py` |
| 46 | LigandAtomSelector        | Analysis                       | 1.2 | `biopipelines/ligand_atom_selector.py` |
| 47 | OpenMM | Analysis                       | 3.3 | `biopipelines/openmm.py` |
| 48 | P2Rank | Analysis                       | 2.0 | `biopipelines/p2rank.py` |
| 49 | PLIP | Analysis                       | 2.2 | `biopipelines/plip.py` |
| 50 | PLM_Sol | Analysis                       | 2.0 | `biopipelines/plm_sol.py` |
| 51 | PoseBusters | Analysis                       | 2.2 | `biopipelines/posebusters.py` |
| 52 | PoseChange                | Analysis                       | 1.1 | `biopipelines/pose_change.py` |
| 53 | Prodigy | Analysis                       | 2.0 | `biopipelines/prodigy.py` |
| 54 | ProLIF | Analysis                       | 2.1 | `biopipelines/prolif.py` |
| 55 | Reduce | Analysis                       | 2.0 | `biopipelines/reduce.py` |
| 56 | RTMScore | Analysis                       | 2.1 | `biopipelines/rtmscore.py` |
| 57 | SASA                      | Analysis                       | 1.0 | `biopipelines/sasa.py` |
| 58 | StructureCluster          | Analysis                       | 1.1 | `biopipelines/structure_cluster.py` |
| 59 | ThermoMPNN | Analysis                       | 2.0 | `biopipelines/thermompnn.py` |
| 60 | VespaG | Analysis                       | 2.0 | `biopipelines/vespag.py` |
| 61 | XTB | Analysis                       | 2.1 | `biopipelines/xtb.py` |
| 62 | BFactor                   | Analysis                       | 1.2 | `biopipelines/bfactor.py` |
| 63 | AiZynthFinder | Cheminformatics                | 2.1 | `biopipelines/aizynthfinder.py` |
| 64 | OpenBabel                 | Cheminformatics                | 1.1 | `biopipelines/openbabel.py` |
| 65 | RDKit                     | Cheminformatics                | 1.2 | `biopipelines/rdkit_descriptors.py` |
| 66 | BayesianAdjuster          | Sequence Statistics            | 1.0 | `biopipelines/bayesian_adjuster.py` |
| 67 | MutationProfiler | Sequence Statistics            | 2.0 | `biopipelines/mutation_profiler.py` |
| 68 | SequenceMetricCorrelation | Sequence Statistics            | 1.0 | `biopipelines/sequence_metric_correlation.py` |
| 69 | ExtractMetrics            | Data Management                | 1.1 | `biopipelines/extract_metrics.py` |
| 70 | Panda                     | Data Management                | 1.5 | `biopipelines/panda.py` |
| 71 | Pool                      | Data Management                | 1.2 | `biopipelines/pool.py` |
| 72 | ReMap | Data Management                | 1.2 | `biopipelines/remap.py` |
| 73 | Selection                 | Data Management                | 1.2 | `biopipelines/selection.py` |
| 74 | SelectionTransfer | Data Management | 1.0 | `biopipelines/selection_transfer.py` |
| 75 | MMseqs2                   | MSAs                           | 1.10 | `biopipelines/mmseqs2.py` |
| 76 | MMseqs2Server             | MSAs                           | 1.9 | `biopipelines/mmseqs2.py` |
| 77 | MSA                       | MSAs                           | 1.4 | `biopipelines/msa.py` |
| 78 | CompoundLibrary           | Inputs & I/O                   | 1.0 | `biopipelines/compound_library.py` |
| 79 | Ligand                    | Inputs & I/O                   | 1.5 | `biopipelines/ligand.py` |
| 80 | Load | Inputs & I/O                   | 1.6 | `biopipelines/load.py` |
| 81 | PDB | Inputs & I/O                   | 1.6 | `biopipelines/pdb.py` |
| 82 | Plot                      | Inputs & I/O                   | 1.0 | `biopipelines/plot.py` |
| 83 | PyMOL | Inputs & I/O                   | 2.1 | `biopipelines/pymol.py` |
| 84 | Scripting | Inputs & I/O                   | 1.3 | `biopipelines/scripting.py` |
| 85 | RCSB                      | Inputs & I/O                   | 1.5 | `biopipelines/rcsb.py` |
| 86 | Sequence                  | Inputs & I/O                   | 1.0 | `biopipelines/sequence.py` |
| 87 | Table                     | Inputs & I/O                   | 1.1 | `biopipelines/table.py` |
| 88 | UniProt | Inputs & I/O                   | 1.1 | `biopipelines/uniprot.py` |

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
