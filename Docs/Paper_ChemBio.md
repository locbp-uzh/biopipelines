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

Technical details should be minimized or moved to Supporting Information. The architecture section should focus on *what the user experiences*, not how the internals work.

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


### To do

- **Figures showing molecular outputs** — structures, binding poses, affinity distributions, not software architecture diagrams
- **Proper references** — not "fewer lines of code" but "accessible to non-specialists"
- **Debug** — especially google colab should be tested

---

## Part 2: Paper Draft

---

# BioPipelines: accessible computational protein and ligand design for chemical biology

**Gianluca Quargnali**, **Pablo Rivera-Fuentes**

Laboratory of Chemical and Biological Probes (LOCBP), Department of Chemistry, University of Zurich, Switzerland

## Abstract

Deep learning methods for protein structure generation, sequence design, and structure prediction have created unprecedented opportunities for protein engineering and drug discovery. However, using these tools requires navigating incompatible software environments, diverse input/output formats, and high-performance computing infrastructure — barriers that limit adoption by chemical biology laboratories. Here we present BioPipelines, an open-source Python framework that allows researchers to define multi-step computational design workflows in a few lines of code. Further, its robust yet modular architecture provides a straightforward way to expand the toolkit with different functionalities with little effort. The framework currently integrates over 30 tools for structure generation, sequence design, structure prediction, compound screening, and analysis. The same workflow code can be prototyped interactively in a Jupyter notebook and then submitted for production-scale runs without modification. We demonstrate applications in de novo protein design, binding site optimization, compound library screening, iterative sequence engineering, fusion protein linker optimization, ligand pose sensitivity analysis, and gene synthesis preparation. We hope this framework will empower researchers allowing them to focus on the scientific question rather than computational logistics. BioPipelines is available under the MIT license at https://github.com/[TODO]/biopipelines.

**Keywords:** protein design, protein engineering, computational workflow, drug discovery, deep learning

## Introduction

Computational protein design has progressed from an expert-only discipline to a broadly useful tool for chemical biology. Deep network hallucination first demonstrated that neural networks trained for structure prediction could be repurposed to generate novel protein folds,^10^ and diffusion-based methods including RFdiffusion,^1^ Chroma,^11^ and FrameDiff^12^ have since made backbone generation routine, producing designable scaffolds with programmable structural features. In parallel, inverse folding models — notably ProteinMPNN^2^ and ESM-IF^13^ — have enabled rapid sequence design for target structures, with extensions such as LigandMPNN^3^ accommodating small-molecule binding partners. On the prediction side, AlphaFold2^4^ transformed the field of structural biology, and subsequent methods including RoseTTAFold,^14^ ESMFold,^15^ AlphaFold3,^16^ and Boltz-2^5^ have extended prediction to biomolecular complexes encompassing proteins, nucleic acids, and small molecules, with binding affinity estimation.

These advances are especially relevant to chemical biology, where researchers routinely need to engineer proteins with altered binding specificity, design enzyme variants for new substrates, or screen compound libraries against protein targets. A typical computational campaign might involve generating protein scaffolds around a ligand binding site, designing sequences compatible with that scaffold, predicting the structures of the resulting designs, and ranking them by predicted binding affinity. In practice, however, chaining these tools together presents significant challenges for laboratories without dedicated computational support. Each tool requires its own software environment, uses different input and output file formats, and must be configured for the computing cluster. Running a multi-tool workflow manually involves writing and debugging shell scripts, tracking intermediate files across tools, and managing job dependencies on the cluster scheduler. These logistical hurdles — rather than the underlying science — are often the rate-limiting step in adopting computational design.

Some workflow frameworks have been developed to address this problem (Table 1). ProtFlow^6^  (unpublished) provides Python wrappers around design tools with cluster job management but the absence of typed data for common entities in chemical biology (ligands, proteins, ...) impairs its modularity. Moreover, it requires writing verbose configuration code and maintaining a running Python process throughout execution. Ovo^7^ offers a web interface built on the Nextflow workflow engine but requires a database server and containerization infrastructure. ProteinDJ^8^ achieves efficient multi-GPU parallelism but restricts users to nine predefined pipeline configurations without support for custom workflows or iterative optimization.

Here we present BioPipelines, a framework designed to make computational protein and ligand design accessible to research groups with minimal computational expertise. The key design goals are clarity and conciseness (workflows defined in a few lines of Python), transparency (every computation generates human-readable scripts that can be inspected), and flexibility (arbitrary tool combinations including iterative multi-cycle optimization). A BioPipeline workflow is meant to read like a description of a scientific experiment. We showcase how this framework can be employed in simple and self-explanatory pipelines.

## Results and Discussion

### Enhancing stability and solubility of a protein

In the following pipeline, the protein of interest Human Carbonic Anhydrase II (PDB: 3KS3) is loaded from the Protein DataBank. The sequence, excluding residues from the active site, is redesigned with the ProteinMPNN model trained on soluble proteins, and the resulting sequences are folded with AlphaFold2 for inspection. Finally, it produces codon-optimized DNA sequences ready for synthesis. The `DNAEncoder` uses organism-specific codon usage tables from CoCoPUTs (HIVE, April 2024), applying as a default a thresholded weighted sampling to avoid rare codons while reducing gene repeats. Multi-organism optimization (e.g., `organism="EC&HS"` for dual expression in *E. coli* and human cells) is supported for constructs that must be expressed in multiple hosts.

```python
from biopipelines.pipeline import *
from biopipelines.entities import *
from biopipelines.protein_mpnn import ProteinMPNN
from biopipelines.alphafold import AlphaFold
from biopipelines.dna_encoder import DNAEncoder

with Pipeline(project="GFPSensor", job="GeneSynthesis"):
    Resources(gpu="A100", time="4:00:00", memory="16GB")

    caII = PDB("3KS3")

    sequences = ProteinMPNN(structures=caII,
                            fixed="94+96+119",
                            num_sequences=10,
                            soluble_model=True)

    folded = AlphaFold(proteins=sequences)

    dna = DNAEncoder(sequences=sequences, organism="EC")
```

Instead of relying on manual inspection, analysis and filtering within the pipeline is also possible. Commonly, one would discard folds that don't resemble the original, or that AlphaFold2 predicts with low confidence. This can easily be done with integrated tools in BioPipelines:

```python
# imports omitted
    conf_change = ConformationalChange(reference_structures = caII,
                                       target_structures = folded)

    filtered = Panda(
        tables=[folded.tables.confidence,conf_change.tables.changes],
        operations=[
            Panda.merge(),
            Panda.filter("RMSD < 1.5 and plddt > 80")
        ],
        pool=folded
    )

    dna = DNAEncoder(sequences=filtered, organism="EC")
```

### Redesign of a protein domain

The following example redesigns the LID, non-essential domain of the protein adenylate kinase (PDB: 4AKE). RFdiffusion is employed to generate ten new protein backbones, replacing the segment A118-160 with a new backbone of length between 50 and 70. ProteinMPNN is used for the inverse-fold of the new backbones, generating two sequences that can fold with the geometry suggested by RFdiffusion. Finally, AlphaFold2 is used for validation. 

```python
# imports omitted
with Pipeline(project="AdenylateKinase", job="LID_Redesign"):
    Resources(gpu="A100", time="4:00:00", memory="16GB")

    kinase = PDB("4AKE")

    backbones = RFdiffusion(pdb=kinase,
                            contigs='A1-117/50-70/A161-214', 
                            num_designs=10)

    sequences = ProteinMPNN(structures=backbones,
                            num_sequences=2,
                            redesigned=backbones.tables.structures.designed)

    refolded = AlphaFold(proteins=sequences)
```

The `redesigned=backbones.tables.structures.designed` argument illustrates how information flows between tools: each backbone has different designed positions (determined by RFdiffusion at runtime), and this per-structure information is automatically passed to ProteinMPNN without any manual file parsing. Pipelines involving RFdiffusion usually involve hundreds if not thousands of designs, hence here instead of only filtering based on the conformational change of the non-designed portion of the protein, we sort them and obtain the best designed based on the model confidence.

```python
    # imports omitted
    conf_change = ConformationalChange(reference_structures = kinase,
                                       target_structures = refolded)

    top3 = Panda(
        tables=[folded.tables.confidence,conf_change.tables.changes],
        operations=[
            Panda.merge(),
            Panda.filter("RMSD < 1.5 and plddt > 80"),
            Panda.sort("plddt"),
            Panda.head(3)
        ],
        pool=folded
    )
```


### Screening a library of compounds against a target protein

Computational screening of compound libraries against protein targets is a common need in chemical biology, particularly for identifying lead compounds or optimizing binding interactions. BioPipelines provides tools for building combinatorial compound libraries and screening them with structure prediction. In the following pipeline, different modification of salicylic acid are screened against homodimer.

```python
# imports omitted
with Pipeline(project="Screening", job="FragmentScreen"):
    Resources(gpu="A100", time="8:00:00", memory="32GB")

    target = Sequence("LFNEIIPLGRLIHMVNQKKDRLLNEYLSPLDITAAQFKVLCSIRCAACITPVELKKVLSVDLGALTRMLDRLVCKGWVERLPNPNDKRGVLVKLTTGGAAICEQCHQLVGQDLHQELTKNLTADEVATLEYLLKKVLP", ids="MarR")

    library = CompoundLibrary(
        library={
            "candidate": "<aryl><carboxylate>",
            "aryl": ["<o-hydroxyphenyl>", "<m-hydroxyphenyl>","<p-hydroxyphenyl>"],
            "o-hydroxyphenyl": r"c1ccc(O)c(c1)",
            "m-hydroxyphenyl": r"c1cc(O)cc(c1)",
            "p-hydroxyphenyl": r"c1c(O)ccc(c1)",
            "carboxylate": r"C(=O)[O-]",
        },
        primary_key="candidate"
    )

    boltz = Boltz2(
        proteins=Bundle(target, target), # both at once -> homodimer
        ligands=Each(library),
        affinity=True
    )

    merged = Panda(
        tables=[boltz.tables.affinity, library.tables.compounds],
        operations=[Panda.merge(on="id")]
    )

    Plot(
        Plot.Scatter(
            data=merged.tables.result,
            x="aryl", y="affinity_pred_value",
            title="Predicted Affinity by Aryl Substituent",
            xlabel="Aryl Group", ylabel="Predicted Affinity", grid=True
        ),
    )
```

### Measuring the influence of mutations on the binding affinity

Understanding how protein mutations affect small-molecule binding is fundamental to medicinal chemistry and resistance biology. A common question is: which residue positions tolerate substitution without disrupting the drug binding pose, and which mutations cause the ligand to shift or dissociate? BioPipelines tools can easily be combined to address this with a saturation mutagenesis-to-pose analysis pipeline. The following pipeline performs computational saturation mutagenesis at the gatekeeper position (Thr315) of ABL1 kinase — the clinically important residue whose T315I mutation confers resistance to imatinib in chronic myeloid leukemia.^9^ All 20 amino acid variants are generated by `Mutagenesis`, and each mutant–imatinib complex is predicted with Boltz-2 including binding affinity estimation. `PoseChange` then superimposes each predicted complex onto the wild-type crystal structure and measures how far the imatinib pose has shifted (RMSD, centroid distance, and orientation angle). By merging these metrics, the pipeline produces a comprehensive structure-activity landscape: mutations that maintain low RMSD and low affinity are predicted to preserve drug sensitivity, while those with large pose deviations suggest resistance mechanisms — directly connecting structural predictions to clinically relevant phenotypes.

```python
# imports omitted
with Pipeline(project="Imatinib", job="PoseSensitivity"):
    Resources(gpu="A100", time="8:00:00", memory="32GB")

    abl1 = PDB("3QRK")
    imatinib = Ligand("STI")

    original = Boltz2(proteins=abl1,
                      ligands=imatinib)

    mutants_sequences = Mutagenesis(original=abl1,
                                    position=93, #AUTH 315
                                    mode="saturation")

    mutants = Boltz2(proteins=mutants_sequences,
                    ligands=imatinib)

    pose = PoseChange(reference_structure=original,
                      sample_structures=mutants,
                      reference_ligand="ATP",
                      sample_ligand="LIG")

    analysis = Panda(
        tables=[pose.tables.changes, mutants.tables.affinity,
                mutants_sequences.tables.sequences],
        operations=[Panda.merge()])

    Plot(
        Plot.Scatter(
            data=analysis.tables.result,
            x="ligand_rmsd", y="affinity_pred_value",
            title="Pose Deviation vs Predicted Affinity",
            xlabel="Ligand RMSD (A)",
            ylabel="Predicted Affinity", grid=True
        ),
        Plot.Bar(
            data=analysis.tables.result,
            x="mutation", y="affinity_pred_value",
            title="Predicted Affinity per Mutant",
            xlabel="Mutation", ylabel="Predicted Affinity",
            x_tick_rotation=45, grid=True
        ),
        Plot.Bar(
            data=analysis.tables.result,
            x="mutation", y="ligand_rmsd",
            title="Ligand RMSD per Mutant",
            xlabel="Mutation", ylabel="Ligand RMSD (A)",
            x_tick_rotation=45, grid=True
        ),
    )
```


### Designing a FRET-based calcium sensor

Genetically encoded FRET sensors typically sandwich a conformationally responsive sensing domain between a donor and acceptor fluorescent protein. A critical design question is which linker lengths keep both fluorescent domains properly folded while maximizing the change in FRET efficiency upon ligand binding. BioPipelines enables systematic screening of linker variants with structural validation under both apo and holo conditions. The following pipeline constructs a calmodulin-based calcium sensor by fusing ECFP (donor), calmodulin (sensing domain), and EYFP (acceptor) with variable-length flexible linkers. Both apo and Ca2+-bound forms are predicted with AlphaFold, and the chromophore distance and inter-domain orientation are compared between states.

```python
# imports omitted
with Pipeline(project="Biosensor", job="CaFRET"):
    Resources(gpu="A100", time="8:00:00", memory="16GB")

    donor = Sequence("MVSKGEELFTGV...")     # ECFP
    cam = PDB("1CFD")                       # Calmodulin
    acceptor = Sequence("MVSKGEELFTG...")   # EYFP
    calcium = Ligand("CA")

    fusions = Fuse(proteins=[donor, cam, acceptor],
                   name="CaFRET",
                   linker="GSGAG",
                   linker_lengths=["3-5", "3-5"])

    apo = Boltz2(proteins=fusions)
    holo = Boltz2(proteins=fusions, 
                  ligands=calcium,
                  msas=apo) # recycle msas

    dist_apo = Distance(structures=apo,
                        residue=["66", "-173"], method="min",
                        metric_name="FRET_distance_apo")
    dist_holo = Distance(structures=holo,
                         residue=["66", "-173"], method="min",
                         metric_name="FRET_distance_holo")

    # Measure inter-domain orientation
    dihedral_apo = Angle(structures=apo,
                         atoms=["64.CA", "66.CA", "-175.CA", "-173.CA"],
                         metric_name="domain_orientation_apo")

    dihedral_holo = Angle(structures=holo,
                          atoms=["64.CA", "66.CA", "-175.CA", "-173.CA"],
                          metric_name="domain_orientation_holo")

    # Merge apo and holo metrics, compute FRET distance change per linker variant
    analysis = Panda(
        tables=[dist_apo.tables.distances, dist_holo.tables.distances,
                dihedral_apo.tables.angles, dihedral_holo.tables.angles],
        operations=[Panda.merge(on="id"),
                    Panda.calculate({"delta_distance": "FRET_distance_holo - FRET_distance_apo",
                                     "delta_orientation": "domain_orientation_holo - domain_orientation_apo"})])

    Plot(
        Plot.Scatter(
            data=analysis.tables.result,
            x="delta_distance", y="delta_orientation",
            title="Calcium-Induced FRET Geometry Change",
            xlabel="Distance Change, Apo to Holo (A)",
            ylabel="Orientation Change (deg)", grid=True
        )
    )
```

The `Fuse` tool generates all linker-length combinations (3 x 3 = 9 constructs), concatenating the three domains with short flexible linkers. Both apo and holo forms are predicted for each construct, and the `Distance` and `Angle` tools measure the chromophore separation (residue 66 in the donor, residue -173 counting from the C-terminus in the acceptor) and inter-domain orientation — the two geometric parameters that determine FRET efficiency. The scatter plot reveals which linker combinations produce the largest calcium-dependent change in FRET geometry, directly guiding sensor optimization.


### Iterative sequence optimization

Many protein engineering campaigns involve multiple rounds of design and evaluation. BioPipelines supports this directly within a single workflow definition:

```python
# imports omitted
with Pipeline(project="Optimization", job="BindingOptimization"):
    Resources(gpu="A100", time="4:00:00", memory="16GB")

    protein = PDB("6U32")
    ligand = Ligand("JF646-HaloTag ligand", ids="JF646")

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

After structure prediction, researchers typically need to filter results by quality metrics, merge data from multiple sources, and generate visualizations. The `Panda` tool in BioPipelines, based on python package pandas, provides declarative table operations (filtering, sorting, merging, grouping, and calculated columns) that also manage the associated structure files. The `pool` argument links table operations to their associated files: when filtering reduces a table to 20 rows, the corresponding 20 structure files are automatically copied to the output folder. This keeps data and structures synchronized without manual file management.

```python
# imports omitted
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
```

### Interactive prototyping and production runs

A key feature for adoption in experimental laboratories is that the same workflow code can be tested and manually inspected prior to scaling up. This feature is built in the architecture of BioPipelines tools. When run in a Jupyter notebook each tool executes immediately, allowing researchers to inspect intermediate results, adjust parameters, and build the workflow step by step. When running in a Jupyter notebook, all tools that produce structures (e.g., AlphaFold, Boltz-2, RFdiffusion) automatically render an interactive py3Dmol 3D viewer inline, allowing researchers to rotate and inspect predicted structures without leaving the notebook. The `Plot` tool generates publication-ready figures as PNG files using matplotlib, which are displayed inline in Jupyter and saved to the output folder during cluster execution. The `PyMOL` tool generates a `.pse` session file when executed on the cluster; in Jupyter, it instead renders an interactive py3Dmol viewer showing the loaded structures, providing immediate visual feedback without requiring a local PyMOL installation. This consistent notebook integration means that every step of a pipeline produces inspectable output during prototyping.

```python
# imports omitted
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

When the workflow is ready for production, the same script is submitted to the cluster with a `biopipelines-submit` command. This means a researcher can prototype a design campaign interactively on a single GPU, verify that the workflow produces sensible intermediate results, and then submit the identical script for a production run with hundreds of designs — all without rewriting any code or learning a separate batch submission workflow.

### Cluster configuration and installation

BioPipelines integrates over 30 tools organized by function (Table 2). Each external tool — such as RFdiffusion, AlphaFold2, or Boltz-2 — requires its own software environment with specific dependencies. Installing these tools manually is often a significant barrier for non-computational groups, as it involves cloning repositories, creating conda environments, downloading model weights, and resolving version conflicts.

BioPipelines addresses this with a built-in installation system. Each tool defines its own installation procedure, and users can install all required software through a single pipeline script:

```python
from biopipelines.pipeline import *
from biopipelines.boltz2 import Boltz2
from biopipelines.pymol import PyMOL

with Pipeline(project="Setup", job="InstallTools"):
    Resources(time="8:00:00", memory="32GB")

    RFdiffusion.install()   # Clones repository, downloads weights, creates environment
    PyMOL.install()         # Creates environment
```

Each `install()` call adds an installation step to the pipeline, just like any other tool. This is particularly important on Google Colab, when often tools have to be installed at every new session. When submitted to the cluster, each step runs the tool's installation script — cloning repositories, creating conda environments, downloading model weights — as a scheduled job. This means installation uses the same familiar workflow as running a design campaign, and the generated scripts can be inspected or rerun if a step fails. Tools that share an environment (e.g., ProteinMPNN uses the same environment as RFdiffusion) indicate this dependency clearly, and tools that require only the base BioPipelines environment need no separate installation. 

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

## Conclusion

BioPipelines makes computational protein and ligand design accessible to chemical biology laboratories by handling the computational logistics — environment management, file format conversion, cluster job scheduling, and data tracking — that typically require dedicated bioinformatics support. Researchers define workflows as concise Python scripts that read as descriptions of the scientific experiment, prototype interactively in Jupyter notebooks, and submit the same code for production runs. The framework's support for compound library screening, iterative optimization, fusion protein engineering, ligand pose sensitivity analysis, and end-to-end gene synthesis preparation addresses common needs in chemical biology that are not covered by existing frameworks.

BioPipelines is freely available under the MIT license at https://github.com/[TODO]/biopipelines with documentation at https://[TODO].readthedocs.io.

## Associated Content

**Supporting Information.** Technical architecture details, DataStream specification, auto-registration mechanism, full tool parameter reference, and additional pipeline examples. **Appendix A:** Transcript of an AI-assisted tool implementation session, showing how a new tool was added to BioPipelines from a GitHub repository URL using a Claude Code Opus 4.6.

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
