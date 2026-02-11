# BioPipelines Paper — ACS Chemical Biology

---

## Part 1: Publication Strategy

### Target journal

**ACS Chemical Biology** — covers the interface of chemistry and biology, including computational tools that enable chemical biology research. The audience consists of experimental chemists and chemical biologists who use computational methods but are not primarily modelers.

### Framing

The paper should be framed not as a software architecture contribution, but as an **enabling tool for chemical biology research**. The emphasis should be on:

- What scientific questions BioPipelines lets chemical biologists answer
- How it lowers the barrier to using computational protein/ligand design
- Concrete applications: enzyme engineering, drug-target screening, ligand-binding optimization
- Accessibility: no programming beyond a few lines of Python, no infrastructure expertise needed

Technical details (DataStreams, auto-registration, `__new__` override) should be minimized or moved to Supporting Information. The architecture section should focus on *what the user experiences*, not how the internals work.

### What strengthens the paper for this audience

- **Real experimental validation** — even a single experimentally characterized design would dramatically strengthen the manuscript
- **Accessible language** — avoid software engineering jargon; explain in terms of the scientific workflow
- **Figures showing molecular outputs** — structures, binding poses, affinity distributions, not software architecture diagrams
- **Comparison framed as accessibility** — not "fewer lines of code" but "accessible to non-specialists"

---

## Part 2: Paper Draft

---

# BioPipelines: accessible computational protein and ligand design for chemical biology

**Gianluca Quargnali**, **Pablo Rivera-Fuentes**

Laboratory of Computational Biophysics and Protein Design (LOCBP), Department of Chemistry, University of Zurich, Switzerland

## Abstract

Deep learning methods for protein structure generation, sequence design, and structure prediction have created unprecedented opportunities for protein engineering and drug discovery. However, using these tools requires navigating incompatible software environments, diverse input/output formats, and high-performance computing infrastructure — barriers that limit adoption by chemical biology laboratories. Here we present BioPipelines, an open-source Python framework that allows researchers to define multi-step computational design workflows in a few lines of code. BioPipelines handles environment management, data conversion between tools, and job scheduling on computing clusters, enabling researchers to focus on the scientific question rather than computational logistics. The framework integrates over 30 tools for structure generation, sequence design, structure prediction, compound screening, and analysis. The same workflow code can be prototyped interactively in a Jupyter notebook and then submitted for production-scale runs without modification. We demonstrate applications in de novo protein design, binding site optimization, compound library screening, and iterative sequence engineering. BioPipelines is available under the MIT license at https://github.com/[TODO]/biopipelines.

**Keywords:** protein design, protein engineering, computational workflow, drug discovery, deep learning

## Introduction

Computational protein design has progressed from an expert-only discipline to a broadly useful tool for chemical biology. Structure generation methods such as RFdiffusion^1^ can create entirely new protein backbones or redesign portions of existing proteins. Sequence design tools including ProteinMPNN^2^ and LigandMPNN^3^ determine amino acid sequences that will fold into a target structure while accommodating small-molecule binding partners. Structure prediction with AlphaFold2^4^ or Boltz-2^5^ validates whether a designed sequence is likely to adopt the intended fold and can predict protein-ligand complex geometries with binding affinity estimates.

These advances are especially relevant to chemical biology, where researchers routinely need to engineer proteins with altered binding specificity, design enzyme variants for new substrates, or screen compound libraries against protein targets. A typical computational campaign might involve generating protein scaffolds around a ligand binding site, designing sequences compatible with that scaffold, predicting the structures of the resulting designs, and ranking them by predicted binding affinity — all before any experimental work begins.

In practice, however, chaining these tools together presents significant challenges for laboratories without dedicated computational support. Each tool requires its own software environment, uses different input and output file formats, and must be configured for the computing cluster. Running a multi-tool workflow manually involves writing and debugging shell scripts, tracking intermediate files across tools, and managing job dependencies on the cluster scheduler. These logistical hurdles — rather than the underlying science — are often the rate-limiting step in adopting computational design.

Several workflow frameworks have been developed to address this problem (Table 1). ProtFlow^6^ provides Python wrappers around design tools with cluster job management but requires writing verbose configuration code and maintaining a running Python process throughout execution. Ovo^7^ offers a web interface built on the Nextflow workflow engine but requires a database server and containerization infrastructure. ProteinDJ^8^ achieves efficient multi-GPU parallelism but restricts users to nine predefined pipeline configurations without support for custom workflows or iterative optimization.

Here we present BioPipelines, a framework designed to make computational protein and ligand design accessible to research groups with minimal computational expertise. The key design goals are conciseness (workflows defined in a few lines of Python), transparency (every computation generates human-readable scripts that can be inspected), and flexibility (arbitrary tool combinations including iterative multi-cycle optimization).

## Results and Discussion

### Workflow definition

A BioPipelines workflow reads like a description of the scientific experiment. The following example redesigns a portion of lysozyme, generates sequences for the new backbone, and predicts whether those sequences will fold correctly:

```python
from biopipelines.pipeline import *
from biopipelines.entities import *
from biopipelines.rfdiffusion import RFdiffusion
from biopipelines.protein_mpnn import ProteinMPNN
from biopipelines.alphafold import AlphaFold

with Pipeline(project="Lysozyme", job="Redesign"):
    Resources(gpu="A100", time="4:00:00", memory="16GB")

    lysozyme = PDB("168L")

    backbones = RFdiffusion(pdb=lysozyme,
                            contigs='50-70/A81-140',
                            num_designs=10)

    sequences = ProteinMPNN(structures=backbones,
                            num_sequences=5,
                            redesigned=backbones.tables.structures.designed)

    predictions = AlphaFold(proteins=sequences)
```

This pipeline (i) fetches the lysozyme crystal structure from the PDB, (ii) generates 10 backbone designs where residues 1-80 are replaced with 50-70 new residues while keeping the C-terminal domain fixed, (iii) designs 5 sequences per backbone with ProteinMPNN — automatically restricting redesign to only the residues that RFdiffusion generated — and (iv) folds all 50 sequences with AlphaFold2 to assess structural plausibility. The framework handles all file conversions, environment switching, and job scheduling.

The `redesigned=backbones.tables.structures.designed` argument illustrates how information flows between tools: each backbone has different designed positions (determined by RFdiffusion at runtime), and this per-structure information is automatically passed to ProteinMPNN without any manual file parsing.

### Interactive prototyping and production runs

A key feature for adoption in experimental laboratories is that the same workflow code works in two modes without any modification. When run in a Jupyter notebook — a familiar interface for many experimentalists — each tool executes immediately, allowing researchers to inspect intermediate results, adjust parameters, and build the workflow step by step. When the workflow is ready for production, the same script is submitted to the cluster with `biopipelines-submit pipeline.py`, where it generates optimized batch jobs.

This means a researcher can prototype a design campaign interactively on a single GPU, verify that the workflow produces sensible intermediate results, and then submit the identical script for a production run with hundreds of designs — all without rewriting any code or learning a separate batch submission workflow.

### Compound library screening

Computational screening of compound libraries against protein targets is a common need in chemical biology, particularly for identifying lead compounds or optimizing binding interactions. BioPipelines provides tools for building combinatorial compound libraries and screening them with structure prediction:

```python
with Pipeline(project="Screening", job="FragmentScreen"):
    Resources(gpu="A100", time="8:00:00", memory="32GB")

    target = Sequence("MDPLNLS...", ids="DRD2")
    cofactor = Ligand("ATP")

    library = CompoundLibrary(
        library={
            "candidate": "<aryl><linker><amide>",
            "aryl": ["<methoxyphenyl>", "<fluorophenyl>"],
            "methoxyphenyl": r"c1cc(OC)cc(c1)",
            "fluorophenyl": r"c1cc(F)ccc1",
            "linker": ["CC", "CCC", "CCCC"],
            "amide": ["C(=O)N", "C(=O)NC"],
        },
        primary_key="candidate"
    )
    # Generates 2 x 3 x 2 = 12 compounds with fragment tracking

    boltz = Boltz2(
        proteins=target,
        ligands=Bundle(Each(library), cofactor),
        affinity=True
    )
    # 12 predictions, each compound co-predicted with the ATP cofactor
```

The `CompoundLibrary` tool expands a combinatorial template defined in SMILES notation into all possible compounds, tracking which fragment was used at each position. The `Bundle(Each(library), cofactor)` expression ensures that each library compound is predicted in the presence of ATP as a cofactor — a common experimental scenario. The framework records which fragments contributed to each prediction, enabling structure-activity analysis at the fragment level.

### Iterative sequence optimization

Many protein engineering campaigns involve multiple rounds of design and evaluation. BioPipelines supports this directly within a single workflow definition:

```python
with Pipeline(project="Optimization", job="BindingOptimization"):
    Resources(gpu="A100", time="4:00:00", memory="16GB")

    protein = PDB("6U32")
    ligand = Ligand(smiles="O=C(...)NCCOCCCCCCCl", ids="HALOTAG_LIGAND")

    current_best = Boltz2(proteins=protein, ligands=ligand)
    current_best = Panda(table=current_best.tables.affinity,
                         operations=[], pool=current_best)

    for cycle in range(3):
        Suffix(f"Cycle{cycle+1}")

        # Identify residues near the ligand binding site
        pocket = DistanceSelector(structures=current_best,
                                  ligand="LIG", distance=5)

        # Generate 1000 sequence variants near the binding site
        variants = LigandMPNN(structures=current_best,
                              ligand="LIG",
                              num_sequences=1000,
                              redesigned=pocket.tables.selections.within)

        # Analyze mutation frequencies and compose new candidates
        profile = MutationProfiler(original=current_best, mutants=variants)
        candidates = MutationComposer(
            frequencies=profile.tables.absolute_frequencies,
            num_sequences=3, mode="weighted_random", max_mutations=3)

        # Predict structures and select the best across all cycles
        predicted = Boltz2(proteins=candidates, ligands=ligand)
        current_best = Panda(
            tables=[current_best.tables.result, predicted.tables.affinity],
            operations=[
                Panda.concat(add_source=True),
                Panda.sort("affinity_pred_value", ascending=True),
                Panda.head(1)
            ],
            pool=[current_best, predicted]
        )
```

Each cycle generates sequence variants near the binding pocket, profiles mutation frequencies, composes new candidates by weighted sampling from the observed distribution, predicts their structures, and selects the overall best design. The `current_best` variable accumulates the best result across cycles, carrying forward the corresponding structure files. This iterative refinement — analogous to directed evolution but performed computationally — can identify beneficial mutations in specific positions before any experimental testing.

### Data management and analysis

After structure prediction, researchers typically need to filter results by quality metrics, merge data from multiple sources, and generate visualizations. The `Panda` tool provides declarative table operations (filtering, sorting, merging, grouping, and calculated columns) that also manage the associated structure files:

```python
# Filter by confidence, sort, and keep top designs with their structure files
best_designs = Panda(
    table=predictions.tables.confidence,
    operations=[
        Panda.filter("plddt > 80"),
        Panda.sort("ptm", ascending=False),
        Panda.head(20)
    ],
    pool=predictions,
    rename="top"
)

# Generate publication-ready plots
Plot(
    Plot.Scatter(data=predictions.tables.confidence,
                 x="plddt", y="ptm",
                 title="AlphaFold Confidence Metrics"),
    Plot.Histogram(data=predictions.tables.confidence,
                   x="plddt", bins=20,
                   title="pLDDT Distribution")
)

# Create a PyMOL session with aligned and colored structures
PyMOL(
    PyMOL.Load(best_designs),
    PyMOL.ColorAF(best_designs),
    PyMOL.Align(),
    session="Top designs"
)
```

The `pool` argument links table operations to their associated files: when filtering reduces a table to 20 rows, the corresponding 20 structure files are automatically copied to the output folder. This keeps data and structures synchronized without manual file management.

### Supported tools and installation

BioPipelines integrates over 30 tools organized by function (Table 2). Each external tool — such as RFdiffusion, AlphaFold2, or Boltz-2 — requires its own software environment with specific dependencies. Installing these tools manually is often a significant barrier for non-computational groups, as it involves cloning repositories, creating conda environments, downloading model weights, and resolving version conflicts.

BioPipelines addresses this with a built-in installation system. Each tool defines its own installation procedure, and users can install all required software through a single pipeline script:

```python
from biopipelines.pipeline import *
from biopipelines.rfdiffusion import RFdiffusion
from biopipelines.protein_mpnn import ProteinMPNN
from biopipelines.boltz2 import Boltz2

with Pipeline(project="Setup", job="InstallTools"):
    Resources(time="8:00:00", memory="32GB")

    RFdiffusion.install()   # Clones repository, downloads weights, creates environment
    ProteinMPNN.install()   # Clones repository (shares environment with RFdiffusion)
    Boltz2.install()        # Creates environment, installs via pip
```

Each `install()` call adds an installation step to the pipeline, just like any other tool. When submitted to the cluster, each step runs the tool's installation script — cloning repositories, creating conda environments, downloading model weights — as a scheduled job. This means installation uses the same familiar workflow as running a design campaign, and the generated scripts can be inspected or rerun if a step fails. Tools that share an environment (e.g., ProteinMPNN uses the same environment as RFdiffusion) indicate this dependency clearly, and tools that require only the base BioPipelines environment need no separate installation.

*Table 2. Tools integrated in BioPipelines, organized by function.*

| Function | Tools | Description |
|----------|-------|-------------|
| Structure generation | RFdiffusion, RFdiffusion All-Atom, RFdiffusion3, BoltzGen | Generate protein backbones, including with bound ligands |
| Sequence design | ProteinMPNN, LigandMPNN, MutationComposer, Mutagenesis | Design amino acid sequences for target structures |
| Structure prediction | AlphaFold2, Boltz-2 | Predict 3D structures and binding affinities |
| Compound handling | Ligand, CompoundLibrary | Fetch and build small-molecule libraries |
| Structure analysis | Distance, Angle, Contacts, ConformationalChange, PoseChange | Measure geometric and interaction properties |
| Sequence statistics | MutationProfiler, SequenceMetricCorrelation, BayesianAdjuster | Analyze mutation patterns and correlations |
| Data management | Panda, ExtractMetrics, Table, LoadOutput | Transform, filter, and manage tabular data |
| Visualization | Plot, PyMOL | Generate plots and molecular visualization sessions |
| Sequence utilities | Fuse, StitchSequences, SplitChains, DNAEncoder | Construct and modify protein sequences |
| MSA generation | MMseqs2 | Generate multiple sequence alignments (basic uniref30 implementation) |

New tools can be added by implementing four methods (parameter validation, input configuration, script generation, and output prediction), and the framework provides base classes that handle environment activation, completion tracking, and standardized output formatting.

### Writing workflows and extending the framework

**IDE-assisted pipeline authoring.** Because BioPipelines is a pure Python package with typed interfaces, users benefit from modern IDE features when writing workflows. We recommend Visual Studio Code with the Python and Pylance extensions, which provide real-time autocompletion, parameter hints, and inline documentation for every tool and parameter. For example, typing `RFdiffusion(` immediately displays all available arguments with their types and descriptions, and passing an incorrect type triggers an underline warning before the script is ever run. This guided authoring experience means that researchers can explore available tool options without consulting the documentation, significantly lowering the learning curve for new users.

**AI-assisted tool development.** The standardized tool interface also makes it straightforward to extend BioPipelines with new tools using large language model coding assistants. In practice, we have found that tools corresponding to published GitHub repositories can be implemented by providing Claude Code (Anthropic, Opus 4.6 model) with the repository URL and a single instruction: *"Implement a BioPipelines tool for this repository, conforming to the existing tool standard."* The assistant reads the repository's installation instructions, input/output formats, and command-line interface, then generates a complete tool module — including the installation script, parameter validation, bash script generation, and output parsing — that integrates directly into the framework. A representative example of this process is provided in the Supporting Information (Appendix A). This approach dramatically reduces the effort required to keep the framework current with the rapidly evolving landscape of computational biology tools, and enables individual laboratories to add specialized or in-house tools without deep familiarity with the framework internals.

### Comparison with existing frameworks

Table 1 compares BioPipelines with other available frameworks. The key distinctions for chemical biology users are: (i) pipeline definitions require minimal code and no infrastructure configuration beyond a single settings file; (ii) the same code works interactively in notebooks and in production on the cluster; (iii) compound library screening with combinatorial control is built in; and (iv) iterative multi-cycle optimization workflows are supported.

*Table 1. Comparison of protein design workflow frameworks.*

| | BioPipelines | ProtFlow | Ovo | ProteinDJ |
|--|-------------|----------|-----|-----------|
| User interface | Python (few lines) | Python (verbose) | Web UI + Nextflow | Nextflow config |
| Interactive prototyping | Jupyter (same code as production) | Jupyter (different from batch) | Web interface | Not supported |
| Custom workflows | Unlimited combinations | Unlimited combinations | Nextflow pipelines | 9 predefined modes |
| Compound screening | Built-in combinatorics | Not supported | Not supported | Not supported |
| Iterative optimization | Built-in (Python loops) | Manual implementation | Not supported | Not supported |
| Infrastructure needed | Python + SLURM | Python + SLURM | Nextflow + database + web server | Nextflow + Apptainer |
| Inspectable outputs | Human-readable bash scripts | None | Nextflow DAG | Nextflow DAG |
| Tool count | 30+ | ~15 | ~10 | 9 |

## Conclusion

BioPipelines makes computational protein and ligand design accessible to chemical biology laboratories by handling the computational logistics — environment management, file format conversion, cluster job scheduling, and data tracking — that typically require dedicated bioinformatics support. Researchers define workflows as concise Python scripts that read as descriptions of the scientific experiment, prototype interactively in Jupyter notebooks, and submit the same code for production runs. The framework's support for compound library screening, iterative optimization, and fragment-level analysis addresses common needs in chemical biology that are not covered by existing frameworks.

BioPipelines is freely available under the MIT license at https://github.com/[TODO]/biopipelines with documentation at https://[TODO].readthedocs.io.

## Associated Content

**Supporting Information.** Technical architecture details, DataStream specification, auto-registration mechanism, full tool parameter reference, and additional pipeline examples. **Appendix A:** Transcript of an AI-assisted tool implementation session, showing how a new tool was added to BioPipelines from a GitHub repository URL using a single natural-language instruction to Claude Code.

## Author Information

**Corresponding Author:** Pablo Rivera-Fuentes — Laboratory of Computational Biophysics and Protein Design (LOCBP), Department of Chemistry, University of Zurich, Switzerland

## Acknowledgments

[TODO]

## Funding

[TODO]

## References

(1) Watson, J. L.; et al. De novo design of protein structure and function with RFdiffusion. *Nature* **2023**, *620*, 1089-1100.

(2) Dauparas, J.; et al. Robust deep learning-based protein sequence design using ProteinMPNN. *Science* **2022**, *378*, 49-56.

(3) Dauparas, J.; et al. Atomic context-conditioned protein sequence design using LigandMPNN. *bioRxiv* **2024**, 2023.12.22.573103.

(4) Jumper, J.; et al. Highly accurate protein structure prediction with AlphaFold. *Nature* **2021**, *596*, 583-589.

(5) Wohlwend, J.; et al. Boltz-2: Exploring the frontiers of biomolecular prediction and generation. *bioRxiv* **2025**, 2025.03.28.646001.

(6) Breitwieser, M. ProtFlow. GitHub repository, 2024. https://github.com/mabr3112/ProtFlow

(7) Harteveld, Z.; et al. Ovo, an open-source ecosystem for de novo protein design. *bioRxiv* **2025**, 2025.11.27.691041.

(8) Silke, A. C.; et al. ProteinDJ: a high-performance and modular protein design pipeline. *Protein Science* **2025**, *35*, e70464.
