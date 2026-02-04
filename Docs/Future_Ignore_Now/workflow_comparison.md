# Workflow Comparison: BioPipelines vs Other Systems

This compares the `rfdaa_pmpnn_ligandmpnn_mmseqs_boltz2.py` pipeline across workflow systems.

## Original BioPipelines (55 lines of logic)

```python
from PipelineScripts.pipeline import *
from PipelineScripts.entities import *
from PipelineScripts.rfdiffusion_allatom import RFdiffusionAllAtom
from PipelineScripts.distance_selector import DistanceSelector
from PipelineScripts.protein_mpnn import ProteinMPNN
from PipelineScripts.ligand_mpnn import LigandMPNN
from PipelineScripts.stitch_sequences import StitchSequences
from PipelineScripts.mmseqs2 import MMseqs2
from PipelineScripts.boltz2 import Boltz2
from PipelineScripts.panda import Panda
from PipelineScripts.plot import Plot
from PipelineScripts.pymol import PyMOL

with Pipeline(project="Examples", job="RFDAA-Pipeline", description="redesign of N terminus"):

    Resources(gpu="any", time="4:00:00", memory="16GB")

    protein = PDB("9IAF")
    rifampicin = Ligand("rifampicin")
    quinolinone = Ligand("A1I1V")

    rfdaa = RFdiffusionAllAtom(pdb=protein, ligand='RFP', contigs='10-20,A6-140', num_designs=3)

    distances = DistanceSelector(structures=rfdaa, ligand="RFP", distance=5,
                                  restrict_to=rfdaa.tables.structures.designed)

    pmpnn = ProteinMPNN(structures=rfdaa, num_sequences=2,
                        redesigned=distances.tables.selections.beyond)

    lmpnn = LigandMPNN(structures=rfdaa, ligand="RFP", num_sequences=2,
                       redesigned=distances.tables.selections.within)

    sequences = StitchSequences(template=pmpnn,
                                substitutions={distances.tables.selections.within: lmpnn})

    Resources(gpu=None, time="4:00:00", memory="16GB")
    msas = MMseqs2(sequences=sequences)

    Resources(gpu="any", time="4:00:00", memory="16GB")
    boltz_quino = Boltz2(proteins=sequences, ligands=quinolinone, msas=msas)
    boltz_rif = Boltz2(proteins=sequences, ligands=rifampicin, msas=msas)

    merged = Panda(
        tables=[boltz_quino.tables.affinity, boltz_rif.tables.affinity],
        operations=[
            Panda.merge(on="id", prefixes=["quino_", "rif_"]),
            Panda.calculate({"affinity_delta": "rif_affinity_pred_value - quino_affinity_pred_value"}),
            Panda.sort("affinity_delta", ascending=True)
        ]
    )

    Plot(Plot.Scatter(data=merged.tables.result, x="rif_affinity_pred_value",
                      y="quino_affinity_pred_value", title="Rif vs Quino"))

    PyMOL(PyMOL.RenderEach(structures=boltz_rif, color_protein="plddt"))
```

---

## Snakemake (~200+ lines)

```python
# Snakefile

configfile: "config.yaml"

# Must manually define all expected outputs upfront
# Snakemake requires knowing the DAG at parse time

DESIGNS = [f"design_{i}" for i in range(config["num_designs"])]
SEQUENCES = expand("{design}_seq{s}", design=DESIGNS, s=range(config["num_sequences"]))

rule all:
    input:
        "results/merged_affinity.csv",
        "results/scatter_plot.png",
        "results/pymol_renders/"

rule fetch_pdb:
    output:
        "inputs/9IAF.pdb"
    shell:
        "wget -O {output} https://files.rcsb.org/download/9IAF.pdb"

rule fetch_ligand_rifampicin:
    output:
        "inputs/rifampicin.sdf"
    shell:
        """
        python scripts/fetch_pubchem.py --name rifampicin --output {output}
        """

rule fetch_ligand_quinolinone:
    output:
        "inputs/quinolinone.sdf"
    shell:
        """
        python scripts/fetch_ccd.py --code A1I1V --output {output}
        """

rule rfdiffusion_allatom:
    input:
        pdb="inputs/9IAF.pdb"
    output:
        structures=expand("rfdaa/design_{i}.pdb", i=range(config["num_designs"])),
        table="rfdaa/structures.csv"
    params:
        contigs="10-20,A6-140",
        ligand="RFP"
    resources:
        gpu=1,
        mem_mb=16000,
        time="4:00:00"
    conda:
        "envs/ProteinEnv.yaml"
    shell:
        """
        python {config[rfdaa_path]}/scripts/run_inference.py \
            inference.input_pdb={input.pdb} \
            inference.ligand={params.ligand} \
            contigmap.contigs=[{params.contigs}] \
            inference.num_designs={config[num_designs]} \
            inference.output_prefix=rfdaa/design
        python scripts/generate_rfdaa_table.py --input_dir rfdaa --output {output.table}
        """

rule distance_selector:
    input:
        structure="rfdaa/{design}.pdb",
        designed_table="rfdaa/structures.csv"
    output:
        selections="distance_selector/{design}_selections.csv"
    params:
        ligand="RFP",
        distance=5
    conda:
        "envs/biopipelines.yaml"
    shell:
        """
        python scripts/distance_selector.py \
            --structure {input.structure} \
            --ligand {params.ligand} \
            --distance {params.distance} \
            --designed_col $(python scripts/get_designed.py {input.designed_table} {wildcards.design}) \
            --output {output.selections}
        """

rule protein_mpnn:
    input:
        structure="rfdaa/{design}.pdb",
        selections="distance_selector/{design}_selections.csv"
    output:
        sequences="pmpnn/{design}_sequences.fa"
    params:
        num_sequences=config["num_sequences"]
    resources:
        gpu=1
    conda:
        "envs/ProteinEnv.yaml"
    shell:
        """
        # Get "beyond" selection from CSV
        REDESIGNED=$(python scripts/get_column.py {input.selections} beyond)

        python {config[pmpnn_path]}/protein_mpnn_run.py \
            --pdb_path {input.structure} \
            --out_folder pmpnn/ \
            --num_seq_per_target {params.num_sequences} \
            --redesigned_residues "$REDESIGNED"
        """

rule ligand_mpnn:
    input:
        structure="rfdaa/{design}.pdb",
        selections="distance_selector/{design}_selections.csv"
    output:
        sequences="lmpnn/{design}_sequences.fa"
    params:
        num_sequences=config["num_sequences"],
        ligand="RFP"
    resources:
        gpu=1
    conda:
        "envs/ligandmpnn_env.yaml"
    shell:
        """
        REDESIGNED=$(python scripts/get_column.py {input.selections} within)

        python {config[lmpnn_path]}/run.py \
            --model_type ligand_mpnn \
            --pdb_path {input.structure} \
            --out_folder lmpnn/ \
            --num_seq_per_target {params.num_sequences} \
            --redesigned_residues "$REDESIGNED" \
            --ligand_name {params.ligand}
        """

rule stitch_sequences:
    input:
        pmpnn="pmpnn/{design}_sequences.fa",
        lmpnn="lmpnn/{design}_sequences.fa",
        selections="distance_selector/{design}_selections.csv"
    output:
        sequences="stitched/{design}_sequences.fa"
    conda:
        "envs/biopipelines.yaml"
    shell:
        """
        python scripts/stitch_sequences.py \
            --template {input.pmpnn} \
            --substitution {input.lmpnn} \
            --positions $(python scripts/get_column.py {input.selections} within) \
            --output {output.sequences}
        """

rule combine_sequences:
    input:
        expand("stitched/{design}_sequences.fa", design=DESIGNS)
    output:
        "stitched/all_sequences.fa"
    shell:
        "cat {input} > {output}"

rule mmseqs2:
    input:
        sequences="stitched/all_sequences.fa"
    output:
        msas=directory("msas/")
    resources:
        mem_mb=16000,
        time="4:00:00"
    conda:
        "envs/biopipelines.yaml"
    shell:
        """
        python scripts/run_mmseqs2.py \
            --sequences {input.sequences} \
            --output_dir {output.msas}
        """

rule boltz2_quinolinone:
    input:
        sequences="stitched/all_sequences.fa",
        msas="msas/",
        ligand="inputs/quinolinone.sdf"
    output:
        structures=directory("boltz_quino/structures/"),
        confidence="boltz_quino/confidence.csv",
        affinity="boltz_quino/affinity.csv"
    resources:
        gpu=1,
        mem_mb=16000,
        time="4:00:00"
    conda:
        "envs/Boltz2Env.yaml"
    shell:
        """
        # Must write config YAML for each sequence
        python scripts/generate_boltz_configs.py \
            --sequences {input.sequences} \
            --ligand {input.ligand} \
            --msas {input.msas} \
            --output_dir boltz_quino/configs/

        # Run Boltz2 for each config
        for config in boltz_quino/configs/*.yaml; do
            boltz predict $config --output_dir {output.structures}
        done

        # Aggregate results
        python scripts/aggregate_boltz_results.py \
            --input_dir {output.structures} \
            --confidence_out {output.confidence} \
            --affinity_out {output.affinity}
        """

rule boltz2_rifampicin:
    input:
        sequences="stitched/all_sequences.fa",
        msas="msas/",
        ligand="inputs/rifampicin.sdf"
    output:
        structures=directory("boltz_rif/structures/"),
        confidence="boltz_rif/confidence.csv",
        affinity="boltz_rif/affinity.csv"
    resources:
        gpu=1
    conda:
        "envs/Boltz2Env.yaml"
    shell:
        """
        python scripts/generate_boltz_configs.py \
            --sequences {input.sequences} \
            --ligand {input.ligand} \
            --msas {input.msas} \
            --output_dir boltz_rif/configs/

        for config in boltz_rif/configs/*.yaml; do
            boltz predict $config --output_dir {output.structures}
        done

        python scripts/aggregate_boltz_results.py \
            --input_dir {output.structures} \
            --confidence_out {output.confidence} \
            --affinity_out {output.affinity}
        """

rule merge_affinities:
    input:
        quino="boltz_quino/affinity.csv",
        rif="boltz_rif/affinity.csv"
    output:
        "results/merged_affinity.csv"
    conda:
        "envs/biopipelines.yaml"
    shell:
        """
        python scripts/merge_tables.py \
            --tables {input.quino} {input.rif} \
            --on id \
            --prefixes quino_ rif_ \
            --calculate "affinity_delta=rif_affinity_pred_value-quino_affinity_pred_value" \
            --sort affinity_delta \
            --output {output}
        """

rule plot:
    input:
        "results/merged_affinity.csv"
    output:
        "results/scatter_plot.png"
    conda:
        "envs/biopipelines.yaml"
    shell:
        """
        python scripts/plot_scatter.py \
            --input {input} \
            --x rif_affinity_pred_value \
            --y quino_affinity_pred_value \
            --title "Rif vs Quino" \
            --output {output}
        """

rule pymol_render:
    input:
        structures="boltz_rif/structures/",
        affinity="boltz_rif/affinity.csv"
    output:
        directory("results/pymol_renders/")
    conda:
        "envs/ProteinEnv.yaml"
    shell:
        """
        python scripts/pymol_render.py \
            --structures {input.structures} \
            --color_protein plddt \
            --data_table {input.affinity} \
            --output_dir {output}
        """
```

**Additional files needed for Snakemake:**
- `config.yaml` - configuration
- `envs/*.yaml` - conda environment definitions (one per environment)
- `scripts/fetch_pubchem.py`
- `scripts/fetch_ccd.py`
- `scripts/generate_rfdaa_table.py`
- `scripts/distance_selector.py`
- `scripts/get_column.py`
- `scripts/get_designed.py`
- `scripts/stitch_sequences.py`
- `scripts/run_mmseqs2.py`
- `scripts/generate_boltz_configs.py`
- `scripts/aggregate_boltz_results.py`
- `scripts/merge_tables.py`
- `scripts/plot_scatter.py`
- `scripts/pymol_render.py`

---

## Nextflow (~250+ lines)

```groovy
#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.pdb_id = "9IAF"
params.ligand_name = "rifampicin"
params.ligand_ccd = "A1I1V"
params.contigs = "10-20,A6-140"
params.num_designs = 3
params.num_sequences = 2
params.distance = 5
params.outdir = "results"

process FETCH_PDB {
    conda 'envs/biopipelines.yaml'

    output:
    path "${params.pdb_id}.pdb"

    script:
    """
    wget -O ${params.pdb_id}.pdb https://files.rcsb.org/download/${params.pdb_id}.pdb
    """
}

process FETCH_LIGAND_PUBCHEM {
    conda 'envs/biopipelines.yaml'

    input:
    val ligand_name

    output:
    path "${ligand_name}.sdf"

    script:
    """
    python $projectDir/scripts/fetch_pubchem.py --name ${ligand_name} --output ${ligand_name}.sdf
    """
}

process FETCH_LIGAND_CCD {
    conda 'envs/biopipelines.yaml'

    input:
    val ccd_code

    output:
    path "${ccd_code}.sdf"

    script:
    """
    python $projectDir/scripts/fetch_ccd.py --code ${ccd_code} --output ${ccd_code}.sdf
    """
}

process RFDIFFUSION_ALLATOM {
    conda 'envs/ProteinEnv.yaml'
    label 'gpu'
    memory '16 GB'
    time '4h'

    input:
    path pdb

    output:
    path "designs/*.pdb", emit: structures
    path "structures.csv", emit: table

    script:
    """
    mkdir -p designs
    python ${params.rfdaa_path}/scripts/run_inference.py \
        inference.input_pdb=${pdb} \
        inference.ligand=RFP \
        contigmap.contigs=[${params.contigs}] \
        inference.num_designs=${params.num_designs} \
        inference.output_prefix=designs/design

    python $projectDir/scripts/generate_rfdaa_table.py --input_dir designs --output structures.csv
    """
}

process DISTANCE_SELECTOR {
    conda 'envs/biopipelines.yaml'

    input:
    path structure
    path designed_table

    output:
    tuple val(design_id), path("${design_id}_selections.csv")

    script:
    design_id = structure.baseName
    """
    python $projectDir/scripts/distance_selector.py \
        --structure ${structure} \
        --ligand RFP \
        --distance ${params.distance} \
        --designed_table ${designed_table} \
        --design_id ${design_id} \
        --output ${design_id}_selections.csv
    """
}

process PROTEIN_MPNN {
    conda 'envs/ProteinEnv.yaml'
    label 'gpu'

    input:
    path structure
    tuple val(design_id), path(selections)

    output:
    tuple val(design_id), path("${design_id}_pmpnn.fa")

    script:
    """
    REDESIGNED=\$(python $projectDir/scripts/get_column.py ${selections} beyond)

    python ${params.pmpnn_path}/protein_mpnn_run.py \
        --pdb_path ${structure} \
        --out_folder ./ \
        --num_seq_per_target ${params.num_sequences} \
        --redesigned_residues "\$REDESIGNED"

    mv seqs/*.fa ${design_id}_pmpnn.fa
    """
}

process LIGAND_MPNN {
    conda 'envs/ligandmpnn_env.yaml'
    label 'gpu'

    input:
    path structure
    tuple val(design_id), path(selections)

    output:
    tuple val(design_id), path("${design_id}_lmpnn.fa")

    script:
    """
    REDESIGNED=\$(python $projectDir/scripts/get_column.py ${selections} within)

    python ${params.lmpnn_path}/run.py \
        --model_type ligand_mpnn \
        --pdb_path ${structure} \
        --out_folder ./ \
        --num_seq_per_target ${params.num_sequences} \
        --redesigned_residues "\$REDESIGNED" \
        --ligand_name RFP

    mv seqs/*.fa ${design_id}_lmpnn.fa
    """
}

process STITCH_SEQUENCES {
    conda 'envs/biopipelines.yaml'

    input:
    tuple val(design_id), path(pmpnn_seqs), path(lmpnn_seqs), path(selections)

    output:
    path "${design_id}_stitched.fa"

    script:
    """
    python $projectDir/scripts/stitch_sequences.py \
        --template ${pmpnn_seqs} \
        --substitution ${lmpnn_seqs} \
        --selections ${selections} \
        --output ${design_id}_stitched.fa
    """
}

process MMSEQS2 {
    conda 'envs/biopipelines.yaml'
    memory '16 GB'
    time '4h'

    input:
    path sequences

    output:
    path "msas/*", emit: msas

    script:
    """
    cat ${sequences} > all_sequences.fa
    mkdir -p msas
    python $projectDir/scripts/run_mmseqs2.py \
        --sequences all_sequences.fa \
        --output_dir msas/
    """
}

process BOLTZ2 {
    conda 'envs/Boltz2Env.yaml'
    label 'gpu'
    memory '16 GB'
    time '4h'

    input:
    path sequences
    path ligand
    path msas
    val ligand_id

    output:
    path "structures/*.pdb", emit: structures
    path "confidence.csv", emit: confidence
    path "affinity.csv", emit: affinity

    script:
    """
    mkdir -p structures configs

    python $projectDir/scripts/generate_boltz_configs.py \
        --sequences ${sequences} \
        --ligand ${ligand} \
        --msas ${msas} \
        --output_dir configs/

    for config in configs/*.yaml; do
        boltz predict \$config --output_dir structures/
    done

    python $projectDir/scripts/aggregate_boltz_results.py \
        --input_dir structures/ \
        --confidence_out confidence.csv \
        --affinity_out affinity.csv
    """
}

process MERGE_AFFINITIES {
    conda 'envs/biopipelines.yaml'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path quino_affinity
    path rif_affinity

    output:
    path "merged_affinity.csv"

    script:
    """
    python $projectDir/scripts/merge_tables.py \
        --tables ${quino_affinity} ${rif_affinity} \
        --on id \
        --prefixes quino_ rif_ \
        --calculate "affinity_delta=rif_affinity_pred_value-quino_affinity_pred_value" \
        --sort affinity_delta \
        --output merged_affinity.csv
    """
}

process PLOT {
    conda 'envs/biopipelines.yaml'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path merged

    output:
    path "scatter_plot.png"

    script:
    """
    python $projectDir/scripts/plot_scatter.py \
        --input ${merged} \
        --x rif_affinity_pred_value \
        --y quino_affinity_pred_value \
        --title "Rif vs Quino" \
        --output scatter_plot.png
    """
}

process PYMOL_RENDER {
    conda 'envs/ProteinEnv.yaml'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path structures
    path affinity

    output:
    path "renders/*"

    script:
    """
    mkdir -p renders
    python $projectDir/scripts/pymol_render.py \
        --structures ${structures} \
        --color_protein plddt \
        --data_table ${affinity} \
        --output_dir renders/
    """
}

workflow {
    // Fetch inputs
    pdb = FETCH_PDB()
    rifampicin = FETCH_LIGAND_PUBCHEM(params.ligand_name)
    quinolinone = FETCH_LIGAND_CCD(params.ligand_ccd)

    // RFdiffusion
    rfdaa = RFDIFFUSION_ALLATOM(pdb)

    // Distance selection - need to scatter over structures
    structures_ch = rfdaa.structures.flatten()
    distances = DISTANCE_SELECTOR(structures_ch, rfdaa.table)

    // MPNN steps - need to join channels properly
    pmpnn = PROTEIN_MPNN(structures_ch, distances)
    lmpnn = LIGAND_MPNN(structures_ch, distances)

    // Stitch - need to join pmpnn, lmpnn, and selections by design_id
    stitched_input = pmpnn.join(lmpnn).join(distances)
    stitched = STITCH_SEQUENCES(stitched_input)

    // MSAs
    msas = MMSEQS2(stitched.collect())

    // Boltz2 for both ligands
    all_seqs = stitched.collect()
    boltz_quino = BOLTZ2(all_seqs, quinolinone, msas.msas, "quino")
    boltz_rif = BOLTZ2(all_seqs, rifampicin, msas.msas, "rif")

    // Merge and visualize
    merged = MERGE_AFFINITIES(boltz_quino.affinity, boltz_rif.affinity)
    PLOT(merged)
    PYMOL_RENDER(boltz_rif.structures, boltz_rif.affinity)
}
```

**Additional files needed:** Same helper scripts as Snakemake, plus `nextflow.config` for executor settings.

---

## CWL (Common Workflow Language) (~400+ lines across multiple files)

CWL requires separate files for each tool, plus a workflow file. Here's a subset:

**workflow.cwl:**
```yaml
cwlVersion: v1.2
class: Workflow

requirements:
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  MultipleInputFeatureRequirement: {}

inputs:
  pdb_id: string
  ligand_name: string
  ligand_ccd: string
  contigs: string
  num_designs: int
  num_sequences: int
  distance: float

outputs:
  merged_affinity:
    type: File
    outputSource: merge_affinities/merged
  scatter_plot:
    type: File
    outputSource: plot/image
  pymol_renders:
    type: Directory
    outputSource: pymol_render/renders

steps:
  fetch_pdb:
    run: tools/fetch_pdb.cwl
    in:
      pdb_id: pdb_id
    out: [pdb_file]

  fetch_rifampicin:
    run: tools/fetch_ligand_pubchem.cwl
    in:
      ligand_name: ligand_name
    out: [ligand_file]

  fetch_quinolinone:
    run: tools/fetch_ligand_ccd.cwl
    in:
      ccd_code: ligand_ccd
    out: [ligand_file]

  rfdiffusion:
    run: tools/rfdiffusion_allatom.cwl
    in:
      pdb: fetch_pdb/pdb_file
      contigs: contigs
      num_designs: num_designs
    out: [structures, table]

  distance_selector:
    run: tools/distance_selector.cwl
    scatter: structure
    in:
      structure: rfdiffusion/structures
      designed_table: rfdiffusion/table
      distance: distance
    out: [selections]

  protein_mpnn:
    run: tools/protein_mpnn.cwl
    scatter: [structure, selections]
    scatterMethod: dotproduct
    in:
      structure: rfdiffusion/structures
      selections: distance_selector/selections
      num_sequences: num_sequences
    out: [sequences]

  ligand_mpnn:
    run: tools/ligand_mpnn.cwl
    scatter: [structure, selections]
    scatterMethod: dotproduct
    in:
      structure: rfdiffusion/structures
      selections: distance_selector/selections
      num_sequences: num_sequences
    out: [sequences]

  stitch_sequences:
    run: tools/stitch_sequences.cwl
    scatter: [pmpnn_seqs, lmpnn_seqs, selections]
    scatterMethod: dotproduct
    in:
      pmpnn_seqs: protein_mpnn/sequences
      lmpnn_seqs: ligand_mpnn/sequences
      selections: distance_selector/selections
    out: [stitched]

  combine_sequences:
    run: tools/combine_fasta.cwl
    in:
      fasta_files: stitch_sequences/stitched
    out: [combined]

  mmseqs2:
    run: tools/mmseqs2.cwl
    in:
      sequences: combine_sequences/combined
    out: [msas]

  boltz2_quinolinone:
    run: tools/boltz2.cwl
    in:
      sequences: combine_sequences/combined
      ligand: fetch_quinolinone/ligand_file
      msas: mmseqs2/msas
    out: [structures, confidence, affinity]

  boltz2_rifampicin:
    run: tools/boltz2.cwl
    in:
      sequences: combine_sequences/combined
      ligand: fetch_rifampicin/ligand_file
      msas: mmseqs2/msas
    out: [structures, confidence, affinity]

  merge_affinities:
    run: tools/merge_tables.cwl
    in:
      table1: boltz2_quinolinone/affinity
      table2: boltz2_rifampicin/affinity
    out: [merged]

  plot:
    run: tools/plot_scatter.cwl
    in:
      data: merge_affinities/merged
    out: [image]

  pymol_render:
    run: tools/pymol_render.cwl
    in:
      structures: boltz2_rifampicin/structures
      affinity: boltz2_rifampicin/affinity
    out: [renders]
```

**tools/rfdiffusion_allatom.cwl:**
```yaml
cwlVersion: v1.2
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: rfdiffusion:latest
  ResourceRequirement:
    coresMin: 4
    ramMin: 16000
    cudaDeviceCount: 1

baseCommand: [python, /app/scripts/run_inference.py]

inputs:
  pdb:
    type: File
    inputBinding:
      prefix: inference.input_pdb=
      separate: false
  contigs:
    type: string
    inputBinding:
      prefix: contigmap.contigs=
      separate: false
  num_designs:
    type: int
    inputBinding:
      prefix: inference.num_designs=
      separate: false

outputs:
  structures:
    type: File[]
    outputBinding:
      glob: "designs/*.pdb"
  table:
    type: File
    outputBinding:
      glob: "structures.csv"
```

**tools/boltz2.cwl:**
```yaml
cwlVersion: v1.2
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: boltz2:latest
  ResourceRequirement:
    coresMin: 4
    ramMin: 16000
    cudaDeviceCount: 1
  InitialWorkDirRequirement:
    listing:
      - $(inputs.sequences)
      - $(inputs.ligand)
      - $(inputs.msas)

baseCommand: [python, /scripts/run_boltz2_workflow.py]

inputs:
  sequences:
    type: File
    inputBinding:
      prefix: --sequences
  ligand:
    type: File
    inputBinding:
      prefix: --ligand
  msas:
    type: Directory
    inputBinding:
      prefix: --msas

outputs:
  structures:
    type: File[]
    outputBinding:
      glob: "structures/*.pdb"
  confidence:
    type: File
    outputBinding:
      glob: "confidence.csv"
  affinity:
    type: File
    outputBinding:
      glob: "affinity.csv"
```

**Additional files needed for CWL:**
- ~15 separate `.cwl` tool definitions
- Same helper scripts
- Docker images or conda environments for each tool
- `inputs.yaml` for parameters

---

## Galaxy (~500+ lines XML)

Galaxy requires XML tool wrappers for each tool. Here's a partial example:

**tools/rfdiffusion_allatom.xml:**
```xml
<tool id="rfdiffusion_allatom" name="RFdiffusion All-Atom" version="1.0.0">
    <description>Generate protein structures with ligand awareness</description>

    <requirements>
        <container type="docker">rfdiffusion:latest</container>
    </requirements>

    <command detect_errors="exit_code"><![CDATA[
        python '$__tool_directory__/scripts/run_inference.py'
            inference.input_pdb='$input_pdb'
            inference.ligand='$ligand_name'
            contigmap.contigs='[$contigs]'
            inference.num_designs=$num_designs
            inference.output_prefix=output/design
        &&
        python '$__tool_directory__/scripts/generate_table.py'
            --input_dir output
            --output '$structures_table'
    ]]></command>

    <inputs>
        <param name="input_pdb" type="data" format="pdb" label="Input PDB structure"/>
        <param name="ligand_name" type="text" value="LIG" label="Ligand residue name"/>
        <param name="contigs" type="text" value="10-20,A6-140" label="Contig specification"/>
        <param name="num_designs" type="integer" value="3" min="1" label="Number of designs"/>
    </inputs>

    <outputs>
        <collection name="structures" type="list" label="Generated structures">
            <discover_datasets pattern="output/(?P&lt;designation&gt;.+)\.pdb" format="pdb"/>
        </collection>
        <data name="structures_table" format="csv" label="Structures table"/>
    </outputs>

    <tests>
        <test>
            <param name="input_pdb" value="test.pdb"/>
            <param name="num_designs" value="1"/>
            <output_collection name="structures" type="list" count="1"/>
        </test>
    </tests>

    <help><![CDATA[
        RFdiffusion All-Atom generates protein backbone structures
        while being aware of ligand positions.
    ]]></help>

    <citations>
        <citation type="doi">10.1038/s41586-023-06415-8</citation>
    </citations>
</tool>
```

**tools/boltz2.xml:**
```xml
<tool id="boltz2" name="Boltz2 Structure Prediction" version="1.0.0">
    <description>Predict protein-ligand complex structures</description>

    <requirements>
        <container type="docker">boltz2:latest</container>
    </requirements>

    <command detect_errors="exit_code"><![CDATA[
        mkdir -p configs structures
        &&
        python '$__tool_directory__/scripts/generate_boltz_configs.py'
            --sequences '$sequences'
            --ligand '$ligand'
            #if $msas
                --msas '$msas'
            #end if
            --output_dir configs/
        &&
        for config in configs/*.yaml; do
            boltz predict \$config --output_dir structures/
        done
        &&
        python '$__tool_directory__/scripts/aggregate_results.py'
            --input_dir structures/
            --confidence_out '$confidence'
            --affinity_out '$affinity'
    ]]></command>

    <inputs>
        <param name="sequences" type="data" format="fasta" label="Input sequences"/>
        <param name="ligand" type="data" format="sdf,mol2" label="Ligand structure"/>
        <param name="msas" type="data" format="a3m" optional="true" label="Pre-computed MSAs"/>
    </inputs>

    <outputs>
        <collection name="structures" type="list" label="Predicted structures">
            <discover_datasets pattern="structures/(?P&lt;designation&gt;.+)\.pdb" format="pdb"/>
        </collection>
        <data name="confidence" format="csv" label="Confidence scores"/>
        <data name="affinity" format="csv" label="Affinity predictions"/>
    </outputs>

    <help><![CDATA[
        Boltz2 predicts protein-ligand complex structures from sequences.
    ]]></help>
</tool>
```

**workflow/rfdaa_pipeline.ga (Galaxy workflow JSON):**
```json
{
    "a]_]galaxy_workflow": "true",
    "format-version": "0.1",
    "name": "RFDAA-PMPNN-LMPNN-MMseqs-Boltz2",
    "steps": {
        "0": {
            "tool_id": "fetch_pdb",
            "inputs": [{"name": "pdb_id", "value": "9IAF"}]
        },
        "1": {
            "tool_id": "fetch_ligand_pubchem",
            "inputs": [{"name": "ligand_name", "value": "rifampicin"}]
        },
        "2": {
            "tool_id": "fetch_ligand_ccd",
            "inputs": [{"name": "ccd_code", "value": "A1I1V"}]
        },
        "3": {
            "tool_id": "rfdiffusion_allatom",
            "input_connections": {"input_pdb": {"id": 0, "output_name": "pdb_file"}},
            "inputs": [
                {"name": "contigs", "value": "10-20,A6-140"},
                {"name": "num_designs", "value": "3"}
            ]
        }
        // ... 15+ more steps with connections
    }
}
```

**Additional files needed for Galaxy:**
- ~15 XML tool wrappers
- Tool shed installation configs
- Test data files
- Helper scripts
- Docker images or conda environments

---

## Summary Comparison

| Aspect | BioPipelines | Snakemake | Nextflow | CWL | Galaxy |
|--------|-------------|-----------|----------|-----|--------|
| **Lines of code** | ~55 | ~200+ | ~250+ | ~400+ | ~500+ |
| **Files needed** | 1 | ~18 | ~18 | ~20 | ~20 |
| **Learning curve** | Python only | Python + DSL | Groovy DSL | YAML + JSON | XML + GUI |
| **Combinatorics** | `Bundle`/`Each` | Manual expand() | scatter + join | scatter + dotproduct | Collection mapping |
| **Table operations** | `Panda.merge()` | Custom script | Custom script | Custom script | Join tool |
| **Column references** | `tool.tables.x.col` | Manual parsing | Manual parsing | Manual parsing | Manual parsing |
| **Resource batching** | `Resources()` | `resources:` per rule | `label` per process | `ResourceRequirement` | Job config |
| **Output prediction** | Automatic | Manual glob patterns | Manual emit | Manual outputBinding | Manual discover |
| **ID tracking** | Automatic | Manual wildcards | Manual channels | Manual scatter | Manual collections |

### Key Differentiators of BioPipelines

1. **Declarative combinatorics**: `Bundle`/`Each` abstracts away the complexity of cartesian products vs grouping
2. **Automatic output prediction**: Tools predict their outputs without manual glob patterns
3. **Table-aware**: Column references like `tool.tables.selections.within` are first-class
4. **Resource batching**: Multiple `Resources()` calls create SLURM job dependencies automatically
5. **Single-file workflows**: No need for separate tool definitions

### What BioPipelines Lacks

1. **Portability**: Tied to specific cluster infrastructure
2. **Reproducibility features**: No container hashes, no workflow provenance
3. **Community ecosystem**: No shared tool repository
4. **Resumability**: Snakemake/Nextflow can resume from failures
5. **Parallelization control**: Other systems offer finer-grained control
