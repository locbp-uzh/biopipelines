# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

# tested:


"""
Minimal debug pipeline exercising every tool with the smallest possible inputs.
Use this to verify that all tools are installed and wired correctly end-to-end.

Structure:
  - RFdiffusion          : 2 de-novo backbones
  - RFdiffusionAllAtom   : 2 ligand-aware backbones (1AKE with AMP)
  - RFdiffusion3         : 2 de-novo backbones
  - ProteinMPNN          : 2 sequences per RFdiffusion backbone
  - LigandMPNN           : 2 sequences per RFdiffusionAllAtom backbone
  - MutationComposer     : 2 composed sequences from ProteinMPNN frequencies
  - Mutagenesis          : saturation at position 1 of a short Sequence
  - Fuse                 : join two short sequences with a linker
  - StitchSequences      : stitch pmpnn / lmpnn sequences
  - DNAEncoder           : codon-optimize pmpnn sequences for E. coli
  - RBSDesigner          : design RBS for the codon-optimized sequences
  - AlphaFold            : fold pmpnn sequences (single_sequence mode, fast)
  - Boltz2               : predict protein-ligand complex for lmpnn sequences
  - Gnina                : dock imatinib into the Boltz2-predicted structures
  - Distance             : measure Cl-to-residue distance in Boltz2 structures
  - Angle                : N-CA-C bond angle in Boltz2 structures
  - DistanceSelector     : select residues within 5 Å of ligand
  - ConformationalChange : RMSD between RFdiffusion backbones and AlphaFold models
  - Contacts             : protein-ligand contacts in Boltz2 structures
  - PoseBusters          : validate Boltz2 ligand poses
  - PoseChange           : compare Gnina poses to crystal reference
  - CABSflex             : flexibility of a single input PDB (2 models)
  - MMseqs2              : generate MSAs for pmpnn sequences
  - MutationProfiler     : profile pmpnn mutations vs original
  - SequenceMetricCorrelation : correlate mutations with Boltz2 affinity
  - BayesianAdjuster     : adjust mutation frequencies by correlation
  - Panda                : filter / merge / sort tables
  - ReMap                : rename IDs
  - ExtractMetrics       : extract metrics for Prism
  - Selection            : expand DistanceSelector selection
  - PyMOL                : visualize Boltz2 structures
  - Plot                 : scatter / histogram / column plots
"""

from biopipelines.pipeline import *

from biopipelines.rfdiffusion import RFdiffusion
from biopipelines.rfdiffusion_allatom import RFdiffusionAllAtom
from biopipelines.rfdiffusion3 import RFdiffusion3
from biopipelines.protein_mpnn import ProteinMPNN
from biopipelines.ligand_mpnn import LigandMPNN
from biopipelines.mutation_composer import MutationComposer
from biopipelines.mutagenesis import Mutagenesis
from biopipelines.fuse import Fuse
from biopipelines.stitch_sequences import StitchSequences
from biopipelines.dna_encoder import DNAEncoder
from biopipelines.rbs_designer import RBSDesigner
from biopipelines.alphafold import AlphaFold
from biopipelines.boltz2 import Boltz2
from biopipelines.gnina import Gnina
from biopipelines.distance import Distance
from biopipelines.angle import Angle
from biopipelines.distance_selector import DistanceSelector
from biopipelines.conformational_change import ConformationalChange
from biopipelines.contacts import Contacts
from biopipelines.posebusters import PoseBusters
from biopipelines.pose_change import PoseChange
from biopipelines.cabsflex import CABSflex
from biopipelines.mmseqs2 import MMseqs2
from biopipelines.mutation_profiler import MutationProfiler
from biopipelines.sequence_metric_correlation import SequenceMetricCorrelation
from biopipelines.bayesian_adjuster import BayesianAdjuster
from biopipelines.panda import Panda
from biopipelines.remap import ReMap
from biopipelines.extract_metrics import ExtractMetrics
from biopipelines.selection import Selection
from biopipelines.pymol import PyMOL
from biopipelines.plot import Plot

with Pipeline(project="Debug",
              job="TestAll",
              description="Minimal smoke-test for every tool"):

    Resources(gpu="A100", time="24:00:00", memory="32GB")

    # ------------------------------------------------------------------
    # Input entities
    # ------------------------------------------------------------------
    abl1 = PDB("3QRK")                          # protein with a ligand (imatinib STI)
    imatinib = Ligand("STI")
    short_seq = Sequence("MKTVRQERLKSIVRILERSKEPVSGAQ", ids="short")
    library = CompoundLibrary({
        "imatinib": "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",
        "analog":   "Cc1ccc(NC(=O)c2ccc(CN3CCN(CC(F)(F)(F))CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",
    })

    # ------------------------------------------------------------------
    # Structure generation
    # ------------------------------------------------------------------
    rfd = RFdiffusion(
        contigs="50-70",
        num_designs=2,
        steps=20,           # fewer steps for speed
    )

    rfdaa = RFdiffusionAllAtom(
        pdb=abl1,
        ligand="STI",
        contigs="10-20,A1-267",
        num_designs=2,
        steps=20,
    )

    rfd3 = RFdiffusion3(
        contig="50-70",
        num_designs=2,
    )

    # ------------------------------------------------------------------
    # Sequence design
    # ------------------------------------------------------------------
    pmpnn = ProteinMPNN(
        structures=rfd3,
        num_sequences=2,
    )

    lmpnn = LigandMPNN(
        structures=rfdaa,
        ligand="STI",
        num_sequences=2,
        redesigned=rfdaa.tables.structures.designed,
    )

    profiler = MutationProfiler(
        original=rfd3,
        mutants=pmpnn,
    )

    composer = MutationComposer(
        frequencies=profiler.tables.absolute_frequencies,
        num_sequences=2,
        mode="weighted_random",
        max_mutations=3,
    )

    sdm = Mutagenesis(
        original=short_seq,
        position=1,
        mode="saturation",
    )

    fused = Fuse(
        sequences=[short_seq, short_seq],
        linker="GSGSG",
        linker_lengths=["1-3"],
    )

    distances_rfd = DistanceSelector(
        structures=rfdaa,
        ligand="STI",
        distance=5.0,
    )

    stitched = StitchSequences(
        template=rfdaa,
        substitutions={
            distances_rfd.tables.selections.beyond: pmpnn,
            distances_rfd.tables.selections.within: lmpnn,
        },
    )

    dna = DNAEncoder(
        sequences=pmpnn,
        organism="EC",
    )

    rbs = RBSDesigner(
        sequences=dna,
        tir="high",
    )

    # ------------------------------------------------------------------
    # Structure prediction
    # ------------------------------------------------------------------
    af = AlphaFold(
        sequences=pmpnn,
        msa_mode="single_sequence",
        num_recycle=1
    )

    boltz_apo = Boltz2(
        proteins=pmpnn
    )
    boltz_holo = Boltz2(
        proteins=pmpnn,
        ligands=imatinib,
        msas=boltz_apo,
        affinity=True
    )

    docking = Gnina(
        structures=boltz_holo,
        compounds=imatinib,
        exhaustiveness=8,
        num_runs=2,
        num_modes=3
    )

    # ------------------------------------------------------------------
    # Analysis
    # ------------------------------------------------------------------
    dist = Distance(
        structures=boltz_holo,
        atom="LIG.Cl",
        residue="145",
        method="min",
        metric_name="cl_dist"
    )

    ang = Angle(
        structures=boltz_holo,
        atoms=("1.N", "1.CA", "1.C"),
        metric_name="nca_angle"
    )

    selector = DistanceSelector(
        structures=boltz_holo,
        ligand="LIG",
        distance=5.0
    )

    expanded = Selection(
        Selection.add(selector.tables.selections.within),
        Selection.expand(2),
        structures=boltz_holo
    )

    conf_change = ConformationalChange(
        reference_structures=rfd,
        target_structures=af,
        atoms="CA"
    )

    contacts = Contacts(
        structures=boltz_holo,
        ligand="LIG",
        contact_threshold=5.0
    )

    pb = PoseBusters(
        structures=boltz_holo,
        ligand="LIG"
    )

    pc = PoseChange(
        reference_structure=abl1,
        sample_structures=docking,
        reference_ligand="STI",
        sample_ligand="LIG"
    )

    flex = CABSflex(
        structures=abl1,
        num_models=2,
        mc_cycles=5,
        mc_steps=5,
        mc_annealing=5,
        aa_rebuild=False
    )

    # ------------------------------------------------------------------
    # MSA generation
    # ------------------------------------------------------------------
    #msas = MMseqs2(sequences=pmpnn)

    # ------------------------------------------------------------------
    # Statistics
    # ------------------------------------------------------------------
    correlation = SequenceMetricCorrelation(
        mutants=pmpnn,
        data=boltz_holo.tables.affinity,
        original=rfd,
        metric="affinity_pred_value",
    )

    adjuster = BayesianAdjuster(
        frequencies=profiler.tables.absolute_frequencies,
        correlations=correlation.tables.correlation_2d,
        mode="min",
        gamma=3.0,
    )

    # ------------------------------------------------------------------
    # Data management
    # ------------------------------------------------------------------
    best = Panda(
        tables=boltz_holo.tables.affinity,
        operations=[
            Panda.sort("affinity_pred_value", ascending=True),
            Panda.head(1),
        ],
        pool=boltz_holo,
        rename="best",
    )

    merged = Panda(
        tables=[boltz_holo.tables.affinity, boltz_holo.tables.confidence],
        operations=[
            Panda.merge(on="id", prefixes=["aff_", "conf_"]),
        ]
    )

    remapped = ReMap(source=pmpnn, onto="design")

    extracted = ExtractMetrics(
        tables=[boltz_holo.tables.affinity],
        metrics=["affinity_pred_value"],
        table_names=["Boltz_holo"],
    )

    # ------------------------------------------------------------------
    # Visualization
    # ------------------------------------------------------------------
    PyMOL(
        PyMOL.Load(boltz_holo),
        PyMOL.ColorAF(boltz_holo, upper=1),
        PyMOL.Align("align"),
        session="test_all_boltz",
    )

    Plot(
        Plot.Scatter(
            data=merged.tables.result,
            x="aff_affinity_pred_value",
            y="conf_confidence_score",
            title="Affinity vs Confidence",
            xlabel="Affinity",
            ylabel="Confidence",
        ),
        Plot.Histogram(
            data=boltz_holo.tables.confidence,
            x="confidence_score",
            bins=5,
            title="Confidence Score Distribution",
        ),
        Plot.Column(
            data=[boltz_apo.tables.confidence, boltz_holo.tables.confidence],
            y="complex_plddt",
            labels=["Apo", "Holo"],
            style="box",
        ),
    )
