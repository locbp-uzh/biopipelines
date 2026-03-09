# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

# tested:


"""
Lazy-ID test pipeline: starts from a PDB and Mutagenesis with lazy IDs,
then chains through every compatible tool to verify lazy ID propagation.

Structure:
  - PDB                  : fetch 3QRK (ABL1 kinase with imatinib STI)
  - Sequence             : extract sequence from 3QRK
  - Mutagenesis          : multi-position saturation ("1+2") → lazy IDs
  - AlphaFold            : fold mutant sequences
  - Boltz2 (apo)         : predict apo structures
  - Boltz2 (holo)        : predict protein-ligand complexes with imatinib
  - Gnina                : dock imatinib into Boltz2 structures
  - Distance             : measure Cl-to-residue distance in Boltz2 structures
  - Angle                : N-CA-C bond angle in Boltz2 structures
  - DistanceSelector     : select residues within 5 Å of ligand
  - Selection            : expand DistanceSelector selection
  - Mutagenesis (2nd)    : selection-based positions from DistanceSelector → lazy IDs
  - ConformationalChange : RMSD between PDB and AlphaFold models
  - Contacts             : protein-ligand contacts in Boltz2 structures
  - PoseBusters          : validate Boltz2 ligand poses
  - PoseChange           : compare Gnina poses to crystal reference
  - MutationProfiler     : profile mutations vs original sequence
  - SequenceMetricCorrelation : correlate mutations with Boltz2 affinity
  - BayesianAdjuster     : adjust mutation frequencies by correlation
  - MutationComposer     : compose new sequences from adjusted frequencies
  - DNAEncoder           : codon-optimize mutant sequences for E. coli
  - RBSDesigner          : design RBS for DNA sequences
  - Fuse                 : fuse mutant sequences with a linker
  - StitchSequences      : stitch mutant sequences into template
  - Panda                : filter / merge tables
  - ReMap                : rename IDs
  - ExtractMetrics       : extract metrics for Prism
  - PyMOL                : visualize Boltz2 structures
  - Plot                 : scatter / histogram / column plots
"""

from biopipelines.pipeline import *

from biopipelines.mutagenesis import Mutagenesis
from biopipelines.alphafold import AlphaFold
from biopipelines.boltz2 import Boltz2
from biopipelines.gnina import Gnina
from biopipelines.distance import Distance
from biopipelines.angle import Angle
from biopipelines.distance_selector import DistanceSelector
from biopipelines.selection import Selection
from biopipelines.conformational_change import ConformationalChange
from biopipelines.contacts import Contacts
from biopipelines.posebusters import PoseBusters
from biopipelines.pose_change import PoseChange
from biopipelines.mutation_profiler import MutationProfiler
from biopipelines.sequence_metric_correlation import SequenceMetricCorrelation
from biopipelines.bayesian_adjuster import BayesianAdjuster
from biopipelines.mutation_composer import MutationComposer
from biopipelines.dna_encoder import DNAEncoder
from biopipelines.rbs_designer import RBSDesigner
from biopipelines.fuse import Fuse
from biopipelines.stitch_sequences import StitchSequences
from biopipelines.panda import Panda
from biopipelines.remap import ReMap
from biopipelines.extract_metrics import ExtractMetrics
from biopipelines.pymol import PyMOL
from biopipelines.plot import Plot

with Pipeline(project="Debug",
              job="TestAllLazy",
              description="Lazy-ID smoke-test: Mutagenesis → all downstream tools"):

    Resources(gpu="A100", time="24:00:00", memory="32GB")

    # ------------------------------------------------------------------
    # Input entities
    # ------------------------------------------------------------------
    abl1 = PDB("3QRK")                          # ABL1 kinase with imatinib
    abl1_seq = Sequence("3QRK")                  # extract sequence from RCSB
    imatinib = Ligand("9DP")

    # ------------------------------------------------------------------
    # Mutagenesis with lazy IDs (multi-position string → bracket pattern)
    # ------------------------------------------------------------------
    sdm = Mutagenesis(
        original=abl1_seq,
        position="1+2",               # multi-position → lazy IDs
        mode="aliphatic",              # A, G, I, L, V → small set for testing
        exclude="G",                   # exclude glycine
    )

    # ------------------------------------------------------------------
    # Structure prediction from lazy-ID sequences
    # ------------------------------------------------------------------
    af = AlphaFold(
        proteins=sdm,
        num_recycle=4,
    )

    boltz_apo = Boltz2(
        proteins=sdm,
    )

    boltz_holo = Boltz2(
        proteins=sdm,
        ligands=imatinib,
        msas=boltz_apo,
    )

    # ------------------------------------------------------------------
    # Docking
    # ------------------------------------------------------------------
    docking = Gnina(
        structures=boltz_holo,
        compounds=imatinib,
        exhaustiveness=8,
        num_runs=2,
        num_modes=3,
    )

    # ------------------------------------------------------------------
    # Structural analysis (all fed by lazy-ID structures)
    # ------------------------------------------------------------------
    dist = Distance(
        structures=boltz_holo,
        atom="LIG.Cl",
        residue="145",
        method="min",
        metric_name="cl_dist",
    )

    ang = Angle(
        structures=boltz_holo,
        atoms=("1.N", "1.CA", "1.C"),
        metric_name="nca_angle",
    )

    selector = DistanceSelector(
        structures=boltz_holo,
        ligand="LIG",
        distance=5.0,
    )

    expanded = Selection(
        Selection.add(selector.tables.selections.within),
        Selection.expand(2),
        structures=boltz_holo,
    )

    conf_change = ConformationalChange(
        reference_structures=abl1,
        target_structures=af,
        atoms="CA",
    )

    contacts = Contacts(
        structures=boltz_holo,
        ligand="LIG",
        contact_threshold=5.0,
    )

    pb = PoseBusters(
        structures=boltz_holo,
        ligand="LIG",
    )

    pc = PoseChange(
        reference_structure=abl1,
        sample_structures=docking,
        reference_ligand="9DP",
        sample_ligand="LIG",
    )

    # ------------------------------------------------------------------
    # Second Mutagenesis: selection-based positions → lazy IDs
    # ------------------------------------------------------------------
    sdm_selection = Mutagenesis(
        original=abl1_seq,
        position=expanded,             # Selection output → lazy IDs
        mode="charged",                # D, E, H, K, R
    )

    # ------------------------------------------------------------------
    # Mutation analysis (lazy-ID mutants vs original)
    # ------------------------------------------------------------------
    profiler = MutationProfiler(
        original=abl1_seq,
        mutants=sdm,
    )

    correlation = SequenceMetricCorrelation(
        mutants=sdm,
        data=boltz_holo.tables.affinity,
        original=abl1_seq,
        metric="affinity_pred_value",
    )

    adjuster = BayesianAdjuster(
        frequencies=profiler.tables.absolute_frequencies,
        correlations=correlation.tables.correlation_2d,
        mode="min",
        gamma=3.0,
    )

    composer = MutationComposer(
        frequencies=profiler.tables.absolute_frequencies,
        num_sequences=2,
        mode="weighted_random",
        max_mutations=3,
    )

    # ------------------------------------------------------------------
    # Sequence processing
    # ------------------------------------------------------------------
    dna = DNAEncoder(
        sequences=sdm,
        organism="EC",
    )

    rbs = RBSDesigner(
        sequences=dna,
        tir="high",
    )

    fused = Fuse(
        sequences=[sdm, sdm],
        linker="GSGSG",
        linker_lengths=["1-3"],
    )

    stitched = StitchSequences(
        template=abl1_seq,
        substitutions={
            "1-10": sdm,
        },
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
        ],
    )

    remapped = ReMap(source=sdm, onto="mutant")

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
        session="test_lazy_boltz",
    )

    Plot(
        Plot.Scatter(
            data=merged.tables.result,
            x="aff_affinity_pred_value",
            y="conf_confidence_score",
            title="Affinity vs Confidence (lazy IDs)",
            xlabel="Affinity",
            ylabel="Confidence",
        ),
        Plot.Histogram(
            data=boltz_holo.tables.confidence,
            x="confidence_score",
            bins=5,
            title="Confidence Score Distribution (lazy IDs)",
        ),
        Plot.Column(
            data=[boltz_apo.tables.confidence, boltz_holo.tables.confidence],
            y="complex_plddt",
            labels=["Apo", "Holo"],
            style="box",
        ),
    )
