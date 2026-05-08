# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Exhaustive BioPipelines-native input audit (Reviewer 3A follow-up).

For every wrapper kwarg whose type union accepts a BioPipelines-native
form (DataStream, StandardizedOutput, TableInfo, TableReference,
Tuple[TableInfo, str]), this pipeline drives one chained step that
feeds the native form rather than a raw value. Audited downstream by
docs/reviewer_evidence/c3_parameter_coverage/audit/audit_cluster_run.py.

Inputs are kept minimal (num_designs=2, steps=20, short contigs,
single-sequence target) so the cluster job stays cheap. Each Suffix
block isolates a single native-input form so the audit can grep the
emitted artefacts per kwarg.
"""

from biopipelines.pipeline import *
from biopipelines.alphafold import AlphaFold
from biopipelines.boltz2 import Boltz2
from biopipelines.boltzgen import BoltzGen
from biopipelines.cabsflex import CABSflex
from biopipelines.conformational_change import ConformationalChange
from biopipelines.dna_encoder import DNAEncoder
from biopipelines.gnina import Gnina
from biopipelines.ligand_mpnn import LigandMPNN
from biopipelines.mmseqs2 import MMseqs2
from biopipelines.mutagenesis import Mutagenesis
from biopipelines.panda import Panda
from biopipelines.posebusters import PoseBusters
from biopipelines.protein_mpnn import ProteinMPNN
from biopipelines.pymol import PyMOL
from biopipelines.rbs_designer import RBSDesigner
from biopipelines.rfdiffusion import RFdiffusion
from biopipelines.rfdiffusion3 import RFdiffusion3
from biopipelines.rfdiffusion_allatom import RFdiffusionAllAtom
from biopipelines.sasa import SASA
from biopipelines.selection import Selection
from biopipelines.distance_selector import DistanceSelector


with Pipeline(project="ToolParameters",
              job="ChainNativeInputs",
              description="Native-input forms (DataStream / StandardizedOutput / TableReference) across every wrapper"):

    Resources(gpu="A100", time="12:00:00", memory="32GB")

    # ──────────────────────────────────────────────────────────────────────
    # Seed inputs
    # ──────────────────────────────────────────────────────────────────────
    abl1 = PDB("3QRK", ids="ABL1", convert="pdb")
    short_seq = Sequence("MKTVRQERLKSIVRILERSKEPVSGAQ", ids="short")
    imatinib = Ligand("9DP", ids="imatinib", codes="LIG")

    # ──────────────────────────────────────────────────────────────────────
    # 1. RFdiffusion → ProteinMPNN: structures as StandardizedOutput
    # ──────────────────────────────────────────────────────────────────────
    Suffix("rfd_to_pmpnn")
    rfd = RFdiffusion(contigs="50-70", num_designs=2, steps=20)
    pmpnn_from_rfd = ProteinMPNN(structures=rfd, num_sequences=2)

    # ──────────────────────────────────────────────────────────────────────
    # 2. RFdiffusionAllAtom → LigandMPNN with redesigned=TableReference
    #    Tests Tuple[TableInfo, str] / TableReference form for `redesigned`.
    # ──────────────────────────────────────────────────────────────────────
    Suffix("rfdaa_to_lmpnn_tableref")
    rfdaa = RFdiffusionAllAtom(
        pdb=abl1,
        ligand="9DP",
        contigs="A227-377,10-20",
        num_designs=2,
        steps=20,
    )
    lmpnn_tableref = LigandMPNN(
        structures=rfdaa,
        ligand="9DP",
        num_sequences=2,
        redesigned=rfdaa.tables.structures.designed,  # TableReference
    )

    # ──────────────────────────────────────────────────────────────────────
    # 3. ProteinMPNN with fixed=TableReference (the second flavour)
    # ──────────────────────────────────────────────────────────────────────
    Suffix("pmpnn_fixed_tableref")
    pmpnn_fixed = ProteinMPNN(
        structures=rfdaa,
        num_sequences=2,
        fixed=rfdaa.tables.structures.fixed,  # TableReference
    )

    # ──────────────────────────────────────────────────────────────────────
    # 4. Sequence chaining: ProteinMPNN.streams.sequences → AlphaFold,
    #    DNAEncoder, RBSDesigner, MMseqs2 (all consume sequence StandardizedOutput).
    # ──────────────────────────────────────────────────────────────────────
    Suffix("seq_chain_af")
    af_chained = AlphaFold(proteins=pmpnn_from_rfd, num_recycle=1)

    Suffix("seq_chain_dna")
    dna_chained = DNAEncoder(sequences=pmpnn_from_rfd, organism="EC")

    Suffix("seq_chain_rbs")
    rbs_chained = RBSDesigner(sequences=dna_chained, tir="medium")

    Suffix("seq_chain_mmseqs")
    mmseqs_chained = MMseqs2(sequences=pmpnn_from_rfd, output_format="csv")

    # ──────────────────────────────────────────────────────────────────────
    # 5. Boltz2 with msas=upstream (StandardizedOutput chaining for msas)
    #    and ligands as Ligand StandardizedOutput.
    # ──────────────────────────────────────────────────────────────────────
    Suffix("boltz_apo")
    boltz_apo = Boltz2(proteins=short_seq)

    Suffix("boltz_holo_msas_chain")
    boltz_holo = Boltz2(
        proteins=short_seq,
        ligands=imatinib,
        msas=boltz_apo,  # StandardizedOutput chaining for msas
    )

    # ──────────────────────────────────────────────────────────────────────
    # 6. Gnina with structures=Boltz2.StandardizedOutput,
    #    compounds=Ligand StandardizedOutput.
    # ──────────────────────────────────────────────────────────────────────
    Suffix("gnina_chain")
    gnina_chained = Gnina(
        structures=boltz_holo,
        compounds=imatinib,
        exhaustiveness=4,
        num_runs=1,
        num_modes=3,
    )

    # ──────────────────────────────────────────────────────────────────────
    # 7. PoseBusters chained from Boltz2 (structures=upstream).
    # ──────────────────────────────────────────────────────────────────────
    Suffix("posebusters_chain")
    pb_chained = PoseBusters(structures=boltz_holo, ligand="LIG", mode="dock")

    # ──────────────────────────────────────────────────────────────────────
    # 8. SASA chained from Boltz2 (structures=upstream).
    # ──────────────────────────────────────────────────────────────────────
    Suffix("sasa_chain")
    sasa_chained = SASA(structures=boltz_holo, ligand="LIG")

    # ──────────────────────────────────────────────────────────────────────
    # 9. ConformationalChange across two Boltz2 outputs.
    # ──────────────────────────────────────────────────────────────────────
    Suffix("confchange_chain")
    cc_chained = ConformationalChange(
        reference_structures=boltz_apo,
        target_structures=boltz_holo,
        atoms="CA",
    )

    # ──────────────────────────────────────────────────────────────────────
    # 10. CABSflex chained from Boltz2 (structures=upstream).
    # ──────────────────────────────────────────────────────────────────────
    Suffix("cabsflex_chain")
    cabsflex_chained = CABSflex(structures=boltz_holo, num_models=3, mc_cycles=10, mc_steps=10)

    # ──────────────────────────────────────────────────────────────────────
    # 11. RFdiffusion3 with pdb=PDB StandardizedOutput.
    # ──────────────────────────────────────────────────────────────────────
    Suffix("rfd3_pdb_input")
    rfd3 = RFdiffusion3(
        pdb=abl1,
        ligand_code="9DP",
        contig="A227-377,10-20",
        num_designs=2,
        num_models=1,
    )

    # ──────────────────────────────────────────────────────────────────────
    # 12. DistanceSelector + Selection(Selection.add(TableReference)).
    #     Covers TableReference-based selection.
    # ──────────────────────────────────────────────────────────────────────
    Suffix("selection_tableref")
    selector = DistanceSelector(structures=boltz_holo, ligand="LIG", distance=5.0)
    expanded = Selection(
        Selection.add(selector.tables.selections.within),
        Selection.expand(2),
        structures=boltz_holo,
    )

    # ──────────────────────────────────────────────────────────────────────
    # 13. PyMOL with selection=TableReference
    # ──────────────────────────────────────────────────────────────────────
    Suffix("pymol_selection_tableref")
    PyMOL(
        PyMOL.Load(structures=boltz_holo),
        PyMOL.Color(color="cyan", selection=selector.tables.selections.within),
        session="native_inputs",
    )

    # ──────────────────────────────────────────────────────────────────────
    # 14. Panda merging multiple TableInfo refs from upstream.
    # ──────────────────────────────────────────────────────────────────────
    Suffix("panda_table_chain")
    merged = Panda(
        tables=[boltz_holo.tables.confidence, boltz_holo.tables.affinity],
        operations=[Panda.merge()],
    )

    # ──────────────────────────────────────────────────────────────────────
    # 15. Mutagenesis with position=TableReference (the column is parsed
    #     into per-row residue indices).
    # ──────────────────────────────────────────────────────────────────────
    Suffix("mutagenesis_tableref")
    sdm_chained = Mutagenesis(
        original=short_seq,
        position=expanded.tables.selections.expanded,  # TableReference column
        mode="saturation",
    )
