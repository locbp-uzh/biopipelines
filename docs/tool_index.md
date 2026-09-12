# BioPipelines Tool Index (SI)

Single-source-of-truth listing for every public-API tool. Generated from
`TOOL_NAME` / `TOOL_VERSION` in `biopipelines/*.py` and grouped by the categories
in [`tool_reference.md`](tool_reference.md).

**Public-API count: 85** (listed below). Internal helpers — `BoltzGenMerge`,
`BoltzGenImport`, `RFDAA_PrepareLigand`, `Mock` — the `TemplateTool` skeleton,
and the 2 base-class scaffolding entries (`base`, `install`) are excluded.

| # | Tool | Category | Version | Source file |
|---|------|----------|---------|-------------|
| 1 | BoltzGen | Structure Generation           | 2.1 | `biopipelines/boltzgen.py` |
| 2 | PocketGen | Structure Generation           | 2.0 | `biopipelines/pocketgen.py` |
| 3 | RFdiffusion | Structure Generation           | 2.2 | `biopipelines/rfdiffusion.py` |
| 4 | RFdiffusion3 | Structure Generation           | 3.4 | `biopipelines/rfdiffusion3.py` |
| 5 | RFdiffusionAllAtom | Structure Generation           | 2.3 | `biopipelines/rfdiffusion_allatom.py` |
| 6 | RFdiffusion2              | Structure Generation           | 1.4 | `biopipelines/rfdiffusion2.py` |
| 7 | HBDesigner | Structure Generation           | 2.1 | `biopipelines/hbdesigner.py` |
| 8 | DNAEncoder                | Sequence Design                | 1.2 | `biopipelines/dna_encoder.py` |
| 9 | Frame2Seq | Sequence Design                | 2.0 | `biopipelines/frame2seq.py` |
| 10 | Fuse | Sequence Design                | 1.3 | `biopipelines/fuse.py` |
| 11 | LigandMPNN | Sequence Design                | 2.5 | `biopipelines/ligand_mpnn.py` |
| 12 | LASErMPNN | Sequence Design                | 2.1 | `biopipelines/lasermpnn.py` |
| 13 | Mutagenesis               | Sequence Design                | 1.2 | `biopipelines/mutagenesis.py` |
| 14 | MutationComposer          | Sequence Design                | 1.0 | `biopipelines/mutation_composer.py` |
| 15 | ProteinMPNN | Sequence Design                | 2.3 | `biopipelines/protein_mpnn.py` |
| 16 | RBSDesigner | Sequence Design                | 2.0 | `biopipelines/rbs_designer.py` |
| 17 | StitchSequences           | Sequence Design                | 1.0 | `biopipelines/stitch_sequences.py` |
| 18 | AlphaFold                 | Structure Prediction & Docking | 1.2 | `biopipelines/alphafold.py` |
| 19 | Boltz2 | Structure Prediction & Docking | 2.4 | `biopipelines/boltz2.py` |
| 20 | DiffDock | Structure Prediction & Docking | 2.1 | `biopipelines/diffdock.py` |
| 21 | DynamicBind | Structure Prediction & Docking | 2.2 | `biopipelines/dynamicbind.py` |
| 22 | ESMFold | Structure Prediction & Docking | 2.0 | `biopipelines/esmfold.py` |
| 23 | ESMFold2 | Structure Prediction & Docking | 2.1 | `biopipelines/esmfold2.py` |
| 24 | Gnina | Structure Prediction & Docking | 2.2 | `biopipelines/gnina.py` |
| 25 | Vina | Structure Prediction & Docking | 2.2 | `biopipelines/gnina.py` |
| 26 | NeuralPLexer | Structure Prediction & Docking | 2.1 | `biopipelines/neuralplexer.py` |
| 27 | PLACER | Structure Prediction & Docking | 2.1 | `biopipelines/placer.py` |
| 28 | ADMETAI | Analysis                       | 2.0 | `biopipelines/admet_ai.py` |
| 29 | AF2BIND | Analysis                       | 2.0 | `biopipelines/af2bind.py` |
| 30 | Aggrescan3D | Analysis                       | 2.1 | `biopipelines/aggrescan3d.py` |
| 31 | Angle                     | Analysis                       | 1.0 | `biopipelines/angle.py` |
| 32 | APBS | Analysis                       | 2.0 | `biopipelines/apbs.py` |
| 33 | BindingData               | Analysis                       | 1.2 | `biopipelines/binding_data.py` |
| 34 | BioEmu | Analysis                       | 2.0 | `biopipelines/bioemu.py` |
| 35 | CABSflex | Analysis                       | 2.1 | `biopipelines/cabsflex.py` |
| 36 | ConformationalChange      | Analysis                       | 1.2 | `biopipelines/conformational_change.py` |
| 37 | Consensus                 | Analysis                       | 1.1 | `biopipelines/consensus.py` |
| 38 | Contacts                  | Analysis                       | 1.0 | `biopipelines/contacts.py` |
| 39 | Distance                  | Analysis                       | 1.0 | `biopipelines/distance.py` |
| 40 | DistanceSelector          | Analysis                       | 2.0 | `biopipelines/distance_selector.py` |
| 41 | DSSP | Analysis                       | 2.0 | `biopipelines/dssp.py` |
| 42 | EnsembleAnalysis          | Analysis                       | 1.0 | `biopipelines/ensemble_analysis.py` |
| 43 | FPocket | Analysis                       | 2.1 | `biopipelines/fpocket.py` |
| 44 | GEMS | Analysis                       | 2.0 | `biopipelines/gems.py` |
| 45 | LigandAtomSelector        | Analysis                       | 1.1 | `biopipelines/ligand_atom_selector.py` |
| 46 | OpenMM | Analysis                       | 3.2 | `biopipelines/openmm.py` |
| 47 | P2Rank | Analysis                       | 2.0 | `biopipelines/p2rank.py` |
| 48 | PLIP | Analysis                       | 2.2 | `biopipelines/plip.py` |
| 49 | PLM_Sol | Analysis                       | 2.0 | `biopipelines/plm_sol.py` |
| 50 | PoseBusters | Analysis                       | 2.2 | `biopipelines/posebusters.py` |
| 51 | PoseChange                | Analysis                       | 1.1 | `biopipelines/pose_change.py` |
| 52 | Prodigy | Analysis                       | 2.0 | `biopipelines/prodigy.py` |
| 53 | ProLIF | Analysis                       | 2.1 | `biopipelines/prolif.py` |
| 54 | Reduce | Analysis                       | 2.0 | `biopipelines/reduce.py` |
| 55 | RTMScore | Analysis                       | 2.1 | `biopipelines/rtmscore.py` |
| 56 | SASA                      | Analysis                       | 1.0 | `biopipelines/sasa.py` |
| 57 | StructureCluster          | Analysis                       | 1.1 | `biopipelines/structure_cluster.py` |
| 58 | ThermoMPNN | Analysis                       | 2.0 | `biopipelines/thermompnn.py` |
| 59 | VespaG | Analysis                       | 2.0 | `biopipelines/vespag.py` |
| 60 | XTB | Analysis                       | 2.1 | `biopipelines/xtb.py` |
| 61 | AiZynthFinder | Cheminformatics                | 2.1 | `biopipelines/aizynthfinder.py` |
| 62 | OpenBabel                 | Cheminformatics                | 1.1 | `biopipelines/openbabel.py` |
| 63 | RDKit                     | Cheminformatics                | 1.2 | `biopipelines/rdkit_descriptors.py` |
| 64 | BayesianAdjuster          | Sequence Statistics            | 1.0 | `biopipelines/bayesian_adjuster.py` |
| 65 | MutationProfiler | Sequence Statistics            | 2.0 | `biopipelines/mutation_profiler.py` |
| 66 | SequenceMetricCorrelation | Sequence Statistics            | 1.0 | `biopipelines/sequence_metric_correlation.py` |
| 67 | ExtractMetrics            | Data Management                | 1.0 | `biopipelines/extract_metrics.py` |
| 68 | Panda                     | Data Management                | 1.4 | `biopipelines/panda.py` |
| 69 | Pool                      | Data Management                | 1.2 | `biopipelines/pool.py` |
| 70 | ReMap | Data Management                | 1.2 | `biopipelines/remap.py` |
| 71 | Selection                 | Data Management                | 1.1 | `biopipelines/selection.py` |
| 72 | MMseqs2                   | MSAs                           | 1.6 | `biopipelines/mmseqs2.py` |
| 73 | MMseqs2Server             | MSAs                           | 1.5 | `biopipelines/mmseqs2.py` |
| 74 | MSA                       | MSAs                           | 1.1 | `biopipelines/msa.py` |
| 75 | CompoundLibrary           | Inputs & I/O                   | 1.0 | `biopipelines/compound_library.py` |
| 76 | Ligand                    | Inputs & I/O                   | 1.4 | `biopipelines/ligand.py` |
| 77 | Load | Inputs & I/O                   | 1.4 | `biopipelines/load.py` |
| 78 | PDB | Inputs & I/O                   | 1.4 | `biopipelines/pdb.py` |
| 79 | Plot                      | Inputs & I/O                   | 1.0 | `biopipelines/plot.py` |
| 80 | PyMOL | Inputs & I/O                   | 2.1 | `biopipelines/pymol.py` |
| 81 | Scripting | Inputs & I/O                   | 1.1 | `biopipelines/scripting.py` |
| 82 | RCSB                      | Inputs & I/O                   | 1.5 | `biopipelines/rcsb.py` |
| 83 | Sequence                  | Inputs & I/O                   | 1.0 | `biopipelines/sequence.py` |
| 84 | Table                     | Inputs & I/O                   | 1.1 | `biopipelines/table.py` |
| 85 | UniProt | Inputs & I/O                   | 1.1 | `biopipelines/uniprot.py` |

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
