# BioPipelines Tool Reference

Complete tool reference, organized into the same categories as the project
[README](https://github.com/locbp-uzh/biopipelines#tools).

For a single flat list of every public-API tool (with category and source-file
mapping), see [tool_index.md](tool_index.md).

---

## [Structure Generation](tool/structure_generation.md)

Generate novel protein structures and pockets.

- [BoltzGen](tool/structure_generation.md#boltzgen) — Protein/peptide/nanobody binder design
- [HBDesigner](tool/structure_generation.md#hbdesigner) — Hydrogen-bond network design
- [PocketGen](tool/structure_generation.md#pocketgen) — Ligand-tailored pocket co-design
- [RFdiffusion](tool/structure_generation.md#rfdiffusion) — Backbone generation
- [RFdiffusion2](tool/structure_generation.md#rfdiffusion2) — Atomic active-site scaffolding
- [RFdiffusion3](tool/structure_generation.md#rfdiffusion3) — Third-generation design
- [RFdiffusionAllAtom](tool/structure_generation.md#rfdiffusionallatom) — All-atom with ligands

---

## [Sequence Design](tool/sequence_design.md)

Design and manipulate sequences.

- [DNAEncoder](tool/sequence_design.md#dnaencoder) — Codon optimization
- [Frame2Seq](tool/sequence_design.md#frame2seq) — Fast inverse folding
- [Fuse](tool/sequence_design.md#fuse) — Fusion proteins with linkers
- [LigandMPNN](tool/sequence_design.md#ligandmpnn) — Ligand-aware sequence design
- [Mutagenesis](tool/sequence_design.md#mutagenesis) — Systematic substitutions
- [MutationComposer](tool/sequence_design.md#mutationcomposer) — Combinatorial mutants
- [ProteinMPNN](tool/sequence_design.md#proteinmpnn) — Sequence design for backbones
- [RBSDesigner](tool/sequence_design.md#rbsdesigner) — Ribosome binding site design
- [StitchSequences](tool/sequence_design.md#stitchsequences) — Sequence combinations

---

## [Structure Prediction & Docking](tool/structure_prediction.md)

Predict protein and complex structures, and dock ligands.

- [AlphaFold](tool/structure_prediction.md#alphafold) — AlphaFold2 prediction
- [Boltz2](tool/structure_prediction.md#boltz2) — Biomolecular complex prediction
- [DiffDock](tool/structure_prediction.md#diffdock) — Blind diffusion docking
- [DynamicBind](tool/structure_prediction.md#dynamicbind) — Flexible-backbone docking
- [ESMFold](tool/structure_prediction.md#esmfold) — Single-sequence prediction (no MSA)
- [Gnina](tool/structure_prediction.md#gnina) — Docking or no-search pose scoring (score/minimize) with CNN
- [NeuralPLexer](tool/structure_prediction.md#neuralplexer) — Protein–ligand complex prediction
- [PLACER](tool/structure_prediction.md#placer) — Ligand-pose ensemble in a bound pocket

---

## [Analysis](tool/analysis.md)

Analyze structures, interactions, stability, and fitness.

- [ADMETAI](tool/analysis.md#admetai) — ADMET endpoint predictions
- [AF2BIND](tool/analysis.md#af2bind) — AF2-based binding-residue prediction
- [Aggrescan3D](tool/analysis.md#aggrescan3d) — Structure-based aggregation propensity
- [Angle](tool/analysis.md#angle) — Bond, torsional, and vector angles
- [APBS](tool/analysis.md#apbs) — Electrostatic surface potential
- [BioEmu](tool/analysis.md#bioemu) — Equilibrium ensemble sampling
- [CABSflex](tool/analysis.md#cabsflex) — Fast flexibility simulation
- [ConformationalChange](tool/analysis.md#conformationalchange) — Structural changes
- [Consensus](tool/analysis.md#consensus) — Per-group aggregation of a resi-csv stream
- [Contacts](tool/analysis.md#contacts) — Contact analysis
- [Distance](tool/analysis.md#distance) — Distance measurements
- [DistanceSelector](tool/analysis.md#distanceselector) — Proximity-based selection
- [DSSP](tool/analysis.md#dssp) — Secondary-structure assignment
- [EnsembleAnalysis](tool/analysis.md#ensembleanalysis) — Per-residue RMSF from any ensemble
- [FPocket](tool/analysis.md#fpocket) — Alpha-sphere pocket detection
- [GEMS](tool/analysis.md#gems) — Protein–ligand affinity (GNN)
- [OpenMM](tool/analysis.md#openmm) — Energy minimization
- [P2Rank](tool/analysis.md#p2rank) — Ligand binding-site prediction
- [PLIP](tool/analysis.md#plip) — Interaction profiling
- [PLM_Sol](tool/analysis.md#plm_sol) — Sequence-based solubility prediction
- [PoseBusters](tool/analysis.md#posebusters) — Pose validation
- [PoseChange](tool/analysis.md#posechange) — Ligand pose comparison
- [Prodigy](tool/analysis.md#prodigy) — Protein–protein affinity
- [ProLIF](tool/analysis.md#prolif) — Interaction fingerprints
- [Reduce](tool/analysis.md#reduce) — Add explicit hydrogens
- [RTMScore](tool/analysis.md#rtmscore) — Pose ranking
- [SASA](tool/analysis.md#sasa) — Solvent-accessible surface area
- [ThermoMPNN](tool/analysis.md#thermompnn) — Fold-stability change (ddG) from structure
- [VespaG](tool/analysis.md#vespag) — Zero-shot fitness from sequence
- [XTB](tool/analysis.md#xtb) — Quantum interaction energy

---

## [Cheminformatics](tool/cheminformatics.md)

Convert molecules and compute chemical descriptors.

- [OpenBabel](tool/cheminformatics.md#openbabel) — Format conversion, 3-D coordinates, protonation
- [RDKit](tool/cheminformatics.md#rdkit) — Molecular descriptors

---

## [Sequence Statistics](tool/statistics.md)

Mutation analysis and frequency optimization.

- [BayesianAdjuster](tool/statistics.md#bayesianadjuster) — Correlation-based frequency adjustment
- [MutationProfiler](tool/statistics.md#mutationprofiler) — Mutation patterns
- [SequenceMetricCorrelation](tool/statistics.md#sequencemetriccorrelation) — Mutation–metric correlations

---

## [MSAs](tool/msas.md)

Generate and convert multiple-sequence alignments.

- [MMseqs2](tool/msas.md#mmseqs2) — MSA generation
- [MMseqs2Server](tool/msas.md#mmseqs2server) — Local MSA server management
- [MSA](tool/msas.md#msa) — MSA format conversion (CSV/A3M)

---

## [Data Management](tool/data_management.md)

Filter, transform, and route tables and streams.

- [ExtractMetrics](tool/data_management.md#extractmetrics) — Metric extraction for Prism
- [Panda](tool/data_management.md#panda) — Unified table transformations
- [Pool](tool/data_management.md#pool) — Gather parallel-run outputs
- [ReMap](tool/data_management.md#remap) — Rename IDs across streams and tables
- [Selection](tool/data_management.md#selection) — Selection string manipulation

---

## [Inputs & I/O](tool/inputs_io.md)

Bring data into a pipeline and take results out.

- [CompoundLibrary](tool/inputs_io.md#compoundlibrary) — Create compound collections
- [Ligand](tool/inputs_io.md#ligand) — Fetch small molecules
- [Load / LoadMultiple](tool/inputs_io.md#load) — Reload previous outputs
- [PDB](tool/inputs_io.md#pdb) — Fetch protein structures
- [Plot](tool/inputs_io.md#plot) — Plots
- [PyMOL](tool/inputs_io.md#pymol) — Session creation and rendering
- [RCSB](tool/inputs_io.md#rcsb) — Search RCSB PDB and download
- [Scripting](tool/inputs_io.md#scripting) — Run a custom two-phase script as a typed step
- [Sequence](tool/inputs_io.md#sequence) — Create sequences from strings
- [Table](tool/inputs_io.md#table) — Load an existing CSV/Excel table
- [UniProt](tool/inputs_io.md#uniprot) — Fetch sequences and annotations
