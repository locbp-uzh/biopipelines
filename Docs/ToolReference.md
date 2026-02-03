# BioPipelines Tool Reference

Complete tool reference organized by category.

---

## [Utilities](Tool/Utilities.md)

Input entities and utility tools.

**Entity Types** (from `PipelineScripts.entities`):
- [PDB](Tool/Utilities.md#pdb) - Fetch protein structures
- [Sequence](Tool/Utilities.md#sequence) - Create sequences from strings
- [Ligand](Tool/Utilities.md#ligand) - Fetch small molecules
- [CompoundLibrary](Tool/Utilities.md#compoundlibrary) - Create compound collections

**Loading & MSA**:
- [LoadOutput / LoadOutputs](Tool/Utilities.md#loadoutput) - Load previous outputs
- [MMseqs2](Tool/Utilities.md#mmseqs2) - MSA generation

**Visualization**:
- [PyMOL](Tool/Utilities.md#pymol) - Session creation
- [Plot](Tool/Utilities.md#plot) - Publication-ready plots

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
- [SDM](Tool/SequenceDesign.md#sdm-sitedirectedmutagenesis) - Systematic substitutions
- [Fuse](Tool/SequenceDesign.md#fuse) - Fusion proteins with linkers
- [StitchSequences](Tool/SequenceDesign.md#stitchsequences) - Sequence combinations
- [SplitChains](Tool/SequenceDesign.md#splitchains) - Split concatenated sequences
- [DNAEncoder](Tool/SequenceDesign.md#dnaencoder) - Codon optimization

---

## [Structure Prediction](Tool/StructurePrediction.md)

Predict protein and complex structures.

- [AlphaFold](Tool/StructurePrediction.md#alphafold) - AlphaFold2 prediction
- [Boltz2](Tool/StructurePrediction.md#boltz2) - Biomolecular complex prediction

---

## [Analysis](Tool/Analysis.md)

Analyze structures and interactions.

- [ResidueAtomDistance](Tool/Analysis.md#residueatomdistance) - Distance measurements
- [DistanceSelector](Tool/Analysis.md#distanceselector) - Proximity-based selection
- [ConformationalChange](Tool/Analysis.md#conformationalchange) - Structural changes
- [MutationProfiler](Tool/Analysis.md#mutationprofiler) - Mutation patterns
- [ProteinLigandContacts](Tool/Analysis.md#proteinligandcontacts) - Contact analysis
- [PoseDistance](Tool/Analysis.md#posedistance) - Ligand pose comparison

---

## [Data Management](Tool/DataManagement.md)

Filter, transform, and manipulate tables.

- [Panda](Tool/DataManagement.md#panda) - Unified table transformations (filter, sort, merge, concat, rank, ...)
- [RemoveDuplicates](Tool/DataManagement.md#removeduplicates) - Cross-table deduplication
- [ExtractMetrics](Tool/DataManagement.md#extractmetrics) - Metric extraction for Prism
- [SelectionEditor](Tool/DataManagement.md#selectioneditor) - Selection string manipulation

## Under development
- PLIP (Protein-Ligand Interaction Profiler) - Interaction analysis 
- DynamicBind- Ligand-specific conformations
- SequenceMetricCorrelation - Mutation-metric correlations
- BayesianAdjuster - Correlation-based frequency adjustment
- SequenceMetricAnalysis - Multi-metric mutation analysis
- ESMFold - Fast single-sequence prediction (under development)
