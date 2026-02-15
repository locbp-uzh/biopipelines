# BioPipelines Tool Reference

Complete tool reference organized by category.

---

## [Utilities](Tool/Utilities.md)

Input entities and utility tools.

**Basic Input Types**:
- [PDB](Tool/Utilities.md#pdb) - Fetch protein structures
- [Sequence](Tool/Utilities.md#sequence) - Create sequences from strings
- [Ligand](Tool/Utilities.md#ligand) - Fetch small molecules
- [CompoundLibrary](Tool/Utilities.md#compoundlibrary) - Create compound collections

**Loading & MSA**:
- [LoadOutput / LoadOutputs](Tool/Utilities.md#loadoutput) - Load previous outputs
- [MMseqs2](Tool/Utilities.md#mmseqs2) - MSA generation

**Visualization**:
- [PyMOL](Tool/Utilities.md#pymol) - Session creation
- [Plot](Tool/Utilities.md#plot) - plots

---

## [Structure Generation](Tool/StructureGeneration.md)

Generate novel protein structures.

- [RFdiffusion](Tool/StructureGeneration.md#rfdiffusion) - Backbone generation
- [RFdiffusionAllAtom](Tool/StructureGeneration.md#rfdiffusionallatom) - All-atom with ligands
- [RFdiffusion3](Tool/StructureGeneration.md#rfdiffusion3) - Third-generation design
- [BoltzGen](Tool/StructureGeneration.md#boltzgen) - End-to-end binder design

---

## [Sequence Design](Tool/SequenceDesign.md)

Design and manipulate sequences.

- [ProteinMPNN](Tool/SequenceDesign.md#proteinmpnn) - Sequence design for backbones
- [LigandMPNN](Tool/SequenceDesign.md#ligandmpnn) - Ligand-aware sequence design
- [MutationComposer](Tool/SequenceDesign.md#mutationcomposer) - Combinatorial mutants
- [Mutagenesis](Tool/SequenceDesign.md#sdm-sitedirectedmutagenesis) - Systematic substitutions
- [Fuse](Tool/SequenceDesign.md#fuse) - Fusion proteins with linkers
- [StitchSequences](Tool/SequenceDesign.md#stitchsequences) - Sequence combinations
- [SplitChains](Tool/SequenceDesign.md#splitchains) - Split concatenated sequences
- [DNAEncoder](Tool/SequenceDesign.md#dnaencoder) - Codon optimization
- [RBSDesigner](Tool/SequenceDesign.md#rbsdesigner) - Ribosome binding site design

---

## [Structure Prediction](Tool/StructurePrediction.md)

Predict protein and complex structures.

- [AlphaFold](Tool/StructurePrediction.md#alphafold) - AlphaFold2 prediction
- [Boltz2](Tool/StructurePrediction.md#boltz2) - Biomolecular complex prediction

---

## [Analysis](Tool/Analysis.md)

Analyze structures and interactions.

- [Distance](Tool/Analysis.md#distance) - Distance measurements
- [Angle](Tool/Analysis.md#angle) - Bond and torsional angles
- [DistanceSelector](Tool/Analysis.md#distanceselector) - Proximity-based selection
- [ConformationalChange](Tool/Analysis.md#conformationalchange) - Structural changes
- [Contacts](Tool/Analysis.md#contacts) - Contact analysis
- [PoseChange](Tool/Analysis.md#posechange) - Ligand pose comparison (to be unified with ConformationalChange)

---

## [Statistics](Tool/Statistics.md)

Mutation analysis and frequency optimization.

- [MutationProfiler](Tool/Statistics.md#mutationprofiler) - Mutation patterns
- [SequenceMetricCorrelation](Tool/Statistics.md#sequencemetriccorrelation) - Mutation-metric correlations
- [BayesianAdjuster](Tool/Statistics.md#bayesianadjuster) - Correlation-based frequency adjustment

---

## [Data Management](Tool/DataManagement.md)

Filter, transform, and manipulate tables.

- [Panda](Tool/DataManagement.md#panda) - Unified table transformations (filter, sort, merge, concat, rank, ...)
- [ExtractMetrics](Tool/DataManagement.md#extractmetrics) - Metric extraction for Prism
- [SelectionEditor](Tool/DataManagement.md#selectioneditor) - Selection string manipulation

