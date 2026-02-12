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

### From journal page:

Articles. Concise, yet comprehensive, reports (usually ≤12 printed pages) of original research presenting an advance of immediate, broad, and lasting impact. Articles are not intended to be follow-up papers, unless they contain new and extensive information that will advance the understanding of the system or biological process. Articles contain an unreferenced abstract of ~250 words. Abstracts should not contain abbreviations or acronyms unless essential. An unheaded but referenced introduction of ≤1000 words should expand on the background of the work. Articles include the following headed sections: combined Results and Discussion, and Methods. In general, Articles include ≤8 display items (figures/tables/schemes) and ~50 references. Supporting Information may be included. Articles must be <6500 words in length, including the abstract, body text, methods, references, and figure/scheme legends. Articles include a graphical Table of Contents entry.

Letter. A cover letter must accompany every manuscript submission. During the submission process, you may type it or paste it into the submission system, or you may attach it as a file.

The cover letter must contain the following elements:

the manuscript title,
the name of the corresponding author,
the name(s) of any other author(s),
a paragraph explaining why the paper is appropriate for ACS Chemical Biology, and
a description of any Supporting Information and/or Review-Only Material.
Additionally, authors should note whether the manuscript was discussed with an ACS Chemical Biology editor before submission.


---

## Part 2: Paper Draft

---

# BioPipelines: accessible computational protein and ligand design for chemical biology

**Gianluca Quargnali**, **Pablo Rivera-Fuentes**

Laboratory of Computational Biophysics and Protein Design (LOCBP), Department of Chemistry, University of Zurich, Switzerland

## Abstract

Deep learning methods for protein structure generation, sequence design, and structure prediction have created unprecedented opportunities for protein engineering and drug discovery. However, using these tools requires navigating incompatible software environments, diverse input/output formats, and high-performance computing infrastructure — barriers that limit adoption by chemical biology laboratories. Here we present BioPipelines, an open-source Python framework that allows researchers to define multi-step computational design workflows in a few lines of code. BioPipelines handles environment management, data conversion between tools, and job scheduling on computing clusters, enabling researchers to focus on the scientific question rather than computational logistics. The framework integrates over 30 tools for structure generation, sequence design, structure prediction, compound screening, and analysis. The same workflow code can be prototyped interactively in a Jupyter notebook and then submitted for production-scale runs without modification. We demonstrate applications in de novo protein design, binding site optimization, compound library screening, iterative sequence engineering, fusion protein linker optimization, ligand pose sensitivity analysis, and gene synthesis preparation. BioPipelines is available under the MIT license at https://github.com/[TODO]/biopipelines.

**Keywords:** protein design, protein engineering, computational workflow, drug discovery, deep learning

## Introduction

Computational protein design has progressed from an expert-only discipline to a broadly useful tool for chemical biology. Deep network hallucination first demonstrated that neural networks trained for structure prediction could be repurposed to generate novel protein folds,^10^ and diffusion-based methods including RFdiffusion,^1^ Chroma,^11^ and FrameDiff^12^ have since made backbone generation routine, producing designable scaffolds with programmable structural features. In parallel, inverse folding models — notably ProteinMPNN^2^ and ESM-IF^13^ — have enabled rapid sequence design for target structures, with extensions such as LigandMPNN^3^ accommodating small-molecule binding partners. On the prediction side, AlphaFold2^4^ transformed structure prediction, and subsequent methods including RoseTTAFold,^14^ ESMFold,^15^ AlphaFold3,^16^ and Boltz-2^5^ have extended prediction to biomolecular complexes encompassing proteins, nucleic acids, and small molecules, with binding affinity estimation. Generative protein language models such as ProGen^17^ have further expanded the design toolkit by generating functional sequences directly from learned evolutionary distributions, without requiring an explicit structural template. Together, these advances provide a comprehensive computational pipeline — from backbone generation through sequence design, structure prediction, and experimental gene synthesis — that can dramatically accelerate protein engineering campaigns.

These advances are especially relevant to chemical biology, where researchers routinely need to engineer proteins with altered binding specificity, design enzyme variants for new substrates, or screen compound libraries against protein targets. A typical computational campaign might involve generating protein scaffolds around a ligand binding site, designing sequences compatible with that scaffold, predicting the structures of the resulting designs, and ranking them by predicted binding affinity — all before any experimental work begins. In practice, however, chaining these tools together presents significant challenges for laboratories without dedicated computational support. Each tool requires its own software environment, uses different input and output file formats, and must be configured for the computing cluster. Running a multi-tool workflow manually involves writing and debugging shell scripts, tracking intermediate files across tools, and managing job dependencies on the cluster scheduler. These logistical hurdles — rather than the underlying science — are often the rate-limiting step in adopting computational design.

Several workflow frameworks have been developed to address this problem (Table 1). ProtFlow^6^  (under development) provides Python wrappers around design tools with cluster job management but the absence of typed data for common entities in chemical biology (ligands, proteins, ...) impairs its modularity. Moreover, it requires writing verbose configuration code and maintaining a running Python process throughout execution. Ovo^7^ offers a web interface built on the Nextflow workflow engine but requires a database server and containerization infrastructure. ProteinDJ^8^ achieves efficient multi-GPU parallelism but restricts users to nine predefined pipeline configurations without support for custom workflows or iterative optimization.

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

    refolded = AlphaFold(proteins=sequences)
```

This pipeline (i) fetches the lysozyme crystal structure from the PDB, (ii) generates 10 backbone designs where residues 1-80 are replaced with 50-70 new residues while keeping the C-terminal domain fixed, (iii) designs 5 sequences per backbone with ProteinMPNN — automatically restricting redesign to only the residues that RFdiffusion generated — and (iv) folds all 50 sequences with AlphaFold2 to assess structural plausibility. The framework handles all file conversions, environment switching, and job scheduling.

The `redesigned=backbones.tables.structures.designed` argument illustrates how information flows between tools: each backbone has different designed positions (determined by RFdiffusion at runtime), and this per-structure information is automatically passed to ProteinMPNN without any manual file parsing.

### Interactive prototyping and production runs

A key feature for adoption in experimental laboratories is that the same workflow code can be tested and manually inspected prior to scaling up This feature is built in the architecture of BioPipelines tools. When run in a Jupyter notebook — a familiar interface for many experimentalists — each tool executes immediately, and provides intermediate visualization output, allowing researchers to inspect intermediate results, adjust parameters, and build the workflow step by step. When the workflow is ready for production, the same script is submitted to the cluster with a `biopipelines-submit` command. This means a researcher can prototype a design campaign interactively on a single GPU, verify that the workflow produces sensible intermediate results, and then submit the identical script for a production run with hundreds of designs — all without rewriting any code or learning a separate batch submission workflow.

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

### From designed protein to DNA: gene synthesis preparation

A common endpoint of any computational design campaign is ordering synthetic genes for experimental characterization. This final step — reverse-translating a protein sequence into codon-optimized DNA — is typically performed manually through vendor websites, one sequence at a time. BioPipelines integrates this step directly into the workflow, so the top designs from a structure prediction run can be converted to synthesis-ready DNA sequences without leaving the pipeline:

```python
from biopipelines.pipeline import *
from biopipelines.entities import *
from biopipelines.rfdiffusion import RFdiffusion
from biopipelines.protein_mpnn import ProteinMPNN
from biopipelines.alphafold import AlphaFold
from biopipelines.panda import Panda
from biopipelines.dna_encoder import DNAEncoder

with Pipeline(project="GFPSensor", job="GeneSynthesis"):
    Resources(gpu="A100", time="4:00:00", memory="16GB")

    gfp = PDB("1GFL")

    backbones = RFdiffusion(pdb=gfp,
                            contigs='A1-142/10-20/A157-230',
                            num_designs=10)

    sequences = ProteinMPNN(structures=backbones,
                            num_sequences=10,
                            redesigned=backbones.tables.structures.designed)

    predictions = AlphaFold(proteins=sequences)

    best = Panda(
        table=predictions.tables.confidence,
        operations=[
            Panda.filter("plddt > 85"),
            Panda.sort("ptm", ascending=False),
            Panda.head(5)
        ],
        pool=predictions,
        rename="top"
    )

    # Reverse-translate to E. coli-optimized DNA for gene synthesis
    dna = DNAEncoder(sequences=best, organism="EC")
```

This pipeline redesigns the loop connecting beta-strands 7 and 8 in GFP (a common target for engineering split-GFP complementation sensors), designs and validates sequences computationally, selects the top 5 by AlphaFold confidence, and produces codon-optimized DNA sequences ready for synthesis ordering. The `DNAEncoder` uses organism-specific codon usage tables from CoCoPUTs (HIVE, April 2024), applying thresholded weighted sampling to avoid rare codons. Multi-organism optimization (e.g., `organism="EC&HS"` for dual expression in *E. coli* and human cells) is supported for constructs that must be expressed in multiple hosts.

### Fusion protein linker optimization

Fusion proteins are central to chemical biology — FRET biosensors, targeted protein degraders (PROTACs/molecular glues), and reporter systems all depend on connecting two functional domains with a linker of appropriate length and flexibility. The critical design question is which linker length keeps both domains properly folded and correctly oriented. BioPipelines enables systematic screening of linker variants with structural validation:

```python
with Pipeline(project="Biosensor", job="LinkerScreen"):
    Resources(gpu="A100", time="8:00:00", memory="16GB")

    donor = PDB("5WJ2")     # mClover3
    acceptor = PDB("5WJ4")  # mRuby3

    fusions = Fuse(
        proteins=[donor, acceptor],
        name="FRET_sensor",
        linker="GGGGSGGGGSGGGGSGGGGSGGGGSGGGGS",
        linker_lengths=["5-25"]
    )

    folded = AlphaFold(proteins=fusions)

    chromophore_dist = Distance(
        structures=folded,
        residue=["66", "-66"],
        method="min",
        metric_name="FRET_distance"
    )

    orientation = Angle(
        structures=folded,
        atoms=["66.CA", "1.CA", "-1.CA", "-66.CA"],
        metric_name="domain_orientation"
    )

    combined = Panda(
        tables=[chromophore_dist.tables.distances, orientation.tables.angles],
        operations=[
            Panda.merge(on="id"),
            Panda.sort("FRET_distance", ascending=True),
        ],
        pool=folded
    )

    Plot(
        Plot.Scatter(
            data=combined.tables.result,
            x="FRET_distance", y="domain_orientation",
            title="Linker Length vs FRET Geometry",
            xlabel="Chromophore Distance (A)",
            ylabel="Domain Orientation (deg)", grid=True
        ),
    )

    PyMOL(
        PyMOL.Load(folded),
        PyMOL.ColorAF(folded),
        PyMOL.Color(folded,
                    selection=fusions.tables.sequences.L1,
                    color="yellow"),
        PyMOL.Align(),
        session="FRET_linker_screen"
    )
```

The `Fuse` tool generates all linker-length variants from 5 to 25 residues (21 constructs), concatenating the donor and acceptor sequences with truncated segments of a (GGGGS)~n~ repeat. AlphaFold predicts the structure of each fusion, and the `Distance` and `Angle` tools measure the chromophore separation and inter-domain orientation — the two geometric parameters that determine FRET efficiency. The results are merged into a single table, enabling direct comparison of how linker length affects sensor geometry. In the PyMOL session, the linker region is highlighted in yellow using the per-construct positional information that `Fuse` records automatically, while the two domains are colored by AlphaFold confidence.

### Ligand pose sensitivity analysis

Understanding how protein mutations affect small-molecule binding is fundamental to medicinal chemistry and resistance biology. A common question is: which residue positions tolerate substitution without disrupting the drug binding pose, and which mutations cause the ligand to shift or dissociate? BioPipelines addresses this with a saturation mutagenesis-to-pose analysis pipeline:

```python
with Pipeline(project="Imatinib", job="PoseSensitivity"):
    Resources(gpu="A100", time="8:00:00", memory="32GB")

    abl1 = PDB("1IEP")
    imatinib = Ligand("STI")

    mutants = Mutagenesis(original=abl1,
                          position=315,
                          mode="saturation",
                          include_original=True)

    predictions = Boltz2(proteins=mutants,
                         ligands=imatinib,
                         affinity=True)

    pose = PoseChange(reference_structure=abl1,
                      sample_structures=predictions,
                      reference_ligand="STI",
                      calculate_centroid=True,
                      calculate_orientation=True)

    contacts = Contacts(structures=predictions,
                        ligand="STI",
                        contact_threshold=4.0)

    analysis = Panda(
        tables=[
            pose.tables.changes,
            contacts.tables.contacts,
            predictions.tables.affinity,
        ],
        operations=[
            Panda.merge(on="id"),
            Panda.sort("ligand_rmsd", ascending=True),
        ],
        pool=predictions
    )

    Plot(
        Plot.Scatter(
            data=analysis.tables.result,
            x="ligand_rmsd", y="contacts",
            title="Pose Deviation vs Binding Contacts",
            xlabel="Ligand RMSD (A)",
            ylabel="Contacts (<4 A)", grid=True
        ),
        Plot.Scatter(
            data=analysis.tables.result,
            x="ligand_rmsd", y="affinity_pred_value",
            title="Pose Deviation vs Predicted Affinity",
            xlabel="Ligand RMSD (A)",
            ylabel="Predicted Affinity", grid=True
        ),
    )

    PyMOL(
        PyMOL.Load(predictions),
        PyMOL.ColorAF(predictions),
        PyMOL.Align(),
        session="abl1_imatinib_mutants"
    )
```

This pipeline performs computational saturation mutagenesis at the gatekeeper position (Thr315) of ABL1 kinase — the clinically important residue whose T315I mutation confers resistance to imatinib in chronic myeloid leukemia.^9^ All 20 amino acid variants are generated by `Mutagenesis`, and each mutant–imatinib complex is predicted with Boltz-2 including binding affinity estimation. `PoseChange` then superimposes each predicted complex onto the wild-type crystal structure and measures how far the imatinib pose has shifted (RMSD, centroid distance, and orientation angle), while `Contacts` counts the number of binding site interactions preserved. By merging these metrics, the pipeline produces a comprehensive structure-activity landscape: mutations that maintain low RMSD and high contact counts are predicted to preserve drug sensitivity, while those with large pose deviations suggest resistance mechanisms — directly connecting structural predictions to clinically relevant phenotypes.

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

### Cluster configuration and installation

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

### Supported tools

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

BioPipelines makes computational protein and ligand design accessible to chemical biology laboratories by handling the computational logistics — environment management, file format conversion, cluster job scheduling, and data tracking — that typically require dedicated bioinformatics support. Researchers define workflows as concise Python scripts that read as descriptions of the scientific experiment, prototype interactively in Jupyter notebooks, and submit the same code for production runs. The framework's support for compound library screening, iterative optimization, fusion protein engineering, ligand pose sensitivity analysis, and end-to-end gene synthesis preparation addresses common needs in chemical biology that are not covered by existing frameworks.

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

(9) Gorre, M. E.; et al. Clinical resistance to STI-571 cancer therapy caused by BCR-ABL gene mutation or amplification. *Science* **2001**, *293*, 876-880.

(10) Anishchenko, I.; et al. De novo protein design by deep network hallucination. *Nature* **2021**, *600*, 547-552.

(11) Ingraham, J. B.; et al. Illuminating protein space with a programmable generative model. *Nature* **2023**, *623*, 1070-1078.

(12) Yim, J.; et al. SE(3) diffusion model with application to protein backbone generation. *Proc. 40th Int. Conf. Mach. Learn. (ICML)*, PMLR 202, **2023**.

(13) Hsu, C.; et al. Learning inverse folding from millions of predicted structures. *Proc. 39th Int. Conf. Mach. Learn. (ICML)*, PMLR 162, **2022**, 8946-8970.

(14) Baek, M.; et al. Accurate prediction of protein structures and interactions using a three-track neural network. *Science* **2021**, *373*, 871-876.

(15) Lin, Z.; et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science* **2023**, *379*, 1123-1130.

(16) Abramson, J.; et al. Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature* **2024**, *630*, 493-500.

(17) Madani, A.; et al. Large language models generate functional protein sequences across diverse families. *Nat. Biotechnol.* **2023**, *41*, 1099-1106.
