# BioPipelines Tool Reference

This document contains the complete, verified tool reference with all parameters, default environments, and output specifications.

## Index

### [Structure Generation](Tool/StructureGeneration.md)

Tools for generating novel protein backbone structures and binders.

- [RFdiffusion](Tool/StructureGeneration.md#rfdiffusion) - Diffusion-based backbone generation
- [RFdiffusionAllAtom](Tool/StructureGeneration.md#rfdiffusionallatom) - All-atom diffusion with ligand modeling
- [RFdiffusion3](Tool/StructureGeneration.md#rfdiffusion3) - Third-generation atomic-level design
- [BoltzGen](Tool/StructureGeneration.md#boltzgen) - End-to-end binder design pipeline

### [Sequence Design](Tool/SequenceDesign.md)

Tools for designing and manipulating protein sequences.

- [ProteinMPNN](Tool/SequenceDesign.md#proteinmpnn) - GNN-based sequence design for backbones
- [LigandMPNN](Tool/SequenceDesign.md#ligandmpnn) - Sequence design optimized for ligand binding
- [MutationComposer](Tool/SequenceDesign.md#mutationcomposer) - Combinatorial mutant generation
- [SDM (SiteDirectedMutagenesis)](Tool/SequenceDesign.md#sdm-sitedirectedmutagenesis) - Systematic amino acid substitutions
- [Fuse](Tool/SequenceDesign.md#fuse) - Fusion protein creation with linkers
- [StitchSequences](Tool/SequenceDesign.md#stitchsequences) - Cartesian product sequence combinations
- [DNAEncoder](Tool/SequenceDesign.md#dnaencoder) - Codon-optimized reverse translation

### [Structure Prediction](Tool/StructurePrediction.md)

Tools for predicting protein and complex structures.

- [AlphaFold](Tool/StructurePrediction.md#alphafold) - AlphaFold2 structure prediction
- [GhostFold](Tool/StructurePrediction.md#ghostfold) - Database-free prediction with synthetic MSAs (under development)
- [ESMFold](Tool/StructurePrediction.md#esmfold) - Fast single-sequence prediction (under development)
- [Boltz2](Tool/StructurePrediction.md#boltz2) - Biomolecular complex prediction
- [RF3](Tool/StructurePrediction.md#rf3) - RoseTTAFold3 prediction (under development)
- [OnionNet](Tool/StructurePrediction.md#onionnet) - Binding affinity prediction (under development)
- [OnionNet2](Tool/StructurePrediction.md#onionnet2) - Improved affinity prediction (under development)

### [Analysis](Tool/Analysis.md)

Tools for analyzing structures, interactions, and mutations.

- [DynamicBind](Tool/Analysis.md#dynamicbind) - Ligand-specific conformations (under development)
- [ResidueAtomDistance](Tool/Analysis.md#residueatomdistance) - Distance measurements
- [PLIP (Protein-Ligand Interaction Profiler)](Tool/Analysis.md#plip-protein-ligand-interaction-profiler) - Interaction analysis (under development)
- [DistanceSelector](Tool/Analysis.md#distanceselector) - Proximity-based residue selection
- [SelectionEditor](Tool/Analysis.md#selectioneditor) - Selection string manipulation
- [ConformationalChange](Tool/Analysis.md#conformationalchange) - Structural change quantification
- [MutationProfiler](Tool/Analysis.md#mutationprofiler) - Mutation pattern analysis
- [SequenceMetricCorrelation](Tool/Analysis.md#sequencemetriccorrelation) - Mutation-metric correlations
- [BayesianAdjuster](Tool/Analysis.md#bayesianadjuster) - Correlation-based frequency adjustment
- [SequenceMetricAnalysis](Tool/Analysis.md#sequencemetricanalysis) - Multi-metric mutation analysis
- [ProteinLigandContacts](Tool/Analysis.md#proteinligandcontacts) - Contact analysis
- [PoseDistance](Tool/Analysis.md#posedistance) - Ligand pose comparison

### [Data Management](Tool/DataManagement.md)

Tools for filtering, ranking, and manipulating data tables.

- [Filter](Tool/DataManagement.md#filter) - Expression-based filtering
- [Rank](Tool/DataManagement.md#rank) - Metric-based ranking with ID renaming
- [SelectBest](Tool/DataManagement.md#selectbest) - Single/multi-objective selection
- [RemoveDuplicates](Tool/DataManagement.md#removeduplicates) - Deduplication
- [MergeTables](Tool/DataManagement.md#mergetables) - Table joining
- [ConcatenateTables](Tool/DataManagement.md#concatenatetables) - Vertical table stacking
- [SliceTable](Tool/DataManagement.md#slicetable) - Row/column extraction
- [ExtractMetrics](Tool/DataManagement.md#extractmetrics) - Metric aggregation
- [AverageByTable](Tool/DataManagement.md#averagebytable) - Grouped averaging

### [Utilities](Tool/Utilities.md)

Tools for loading data, fetching structures, and visualization.

- [LoadOutput / LoadOutputs](Tool/Utilities.md#loadoutput--loadoutputs) - Load previous pipeline outputs
- [MMseqs2](Tool/Utilities.md#mmseqs2) - MSA generation
- [MMseqs2Server](Tool/Utilities.md#mmseqs2server) - Local MSA server
- [CompoundLibrary](Tool/Utilities.md#compoundlibrary) - Chemical compound management
- [PDB](Tool/Utilities.md#pdb) - Protein structure fetching
- [Ligand](Tool/Utilities.md#ligand) - Small molecule fetching
- [PyMOL](Tool/Utilities.md#pymol) - Session creation and visualization
