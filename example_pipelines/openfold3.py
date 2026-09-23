# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

# tested:

"""OpenFold3 with the same inputs as the Boltz2 example, so the two can be compared directly.

The parallel is the point: OpenFold3 takes the same axes and the same combinatorics, so a
campaign written for one runs on the other by changing the tool name. Where a case from
`boltz2.py` is missing here, it is missing because OpenFold3 does not do it — there is no
affinity head and no covalent linkage, so glycosylation, covalent bonds and contact constraints
have no counterpart. Those stay Boltz2's.
"""

from biopipelines.pipeline import *
from biopipelines import OpenFold3, MMseqs2, Bundle, Each

with Pipeline(project="Examples",
              job="OpenFold3",
              description="OpenFold3 across input shapes, combinatorics and options"):

    Resources(gpu="A100",
              time="24:00:00",
              memory="32GB")

    # =========================================================================
    # Basic usage
    # =========================================================================

    # 1: direct sequence
    Suffix("1")
    of3_seq = OpenFold3(
        proteins=Sequence("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH")
    )

    # 2: PDB fetch as the protein input
    Suffix("2")
    lysozyme = PDB("1AKI", ids="LYZ")
    of3_pdb = OpenFold3(
        proteins=lysozyme
    )

    # 3: PDB output format -- per-atom pLDDT lands in the B-factor column
    Suffix("3")
    of3_pdb_format = OpenFold3(
        proteins=lysozyme,
        output_format="pdb"
    )

    # 4: protein + ligand
    Suffix("4")
    of3_ligand = OpenFold3(
        proteins=Sequence("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH"),
        ligands=Ligand("ethanol")
    )

    # 5: CompoundLibrary with several ligands
    Suffix("5")
    compounds = CompoundLibrary({
        'ethanol': 'CCO',
        'methanol': 'CO',
        'propanol': 'CCCO'
    })
    of3_multi = OpenFold3(
        proteins=lysozyme,
        ligands=compounds
    )

    # =========================================================================
    # Combinatorics -- identical semantics to the Boltz2 example
    # =========================================================================

    protein_a = Sequence("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH", ids="ProteinA")
    protein_b = Sequence("MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAIGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAALINMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAKRVITTFRTGTWDAYKNL", ids="ProteinB")
    ligand_library = Ligand(['aspirin', 'caffeine', 'ibuprofen'])

    # 6: Each x Each -- 2 proteins x 3 ligands = 6 predictions
    Suffix("6")
    of3_each = OpenFold3(
        proteins=Each(protein_a, protein_b),
        ligands=ligand_library
    )

    # 7: bundled ligands -- each protein folded with all three ligands present, 2 predictions
    Suffix("7")
    of3_bundle_ligands = OpenFold3(
        proteins=Each(protein_a, protein_b),
        ligands=Bundle(ligand_library)
    )

    # 8: bundled proteins -- both chains with each ligand, 3 predictions
    Suffix("8")
    of3_bundle_proteins = OpenFold3(
        proteins=Bundle(protein_a, protein_b),
        ligands=ligand_library
    )

    # 9: nested -- for each library ligand, bundle it with a common cofactor
    Suffix("9")
    atp = Ligand("ATP")
    of3_nested = OpenFold3(
        proteins=protein_a,
        ligands=Bundle(Each(ligand_library), atp)
    )

    # =========================================================================
    # Nucleic acids
    # =========================================================================

    # 10: dsDNA -- the reverse-complement strand is generated as a second chain
    Suffix("10")
    dna_strand = Sequence("GATTACAGATTACA", type="dna", ids="DNA_strand")
    of3_dna = OpenFold3(
        dsDNA=dna_strand
    )

    # 11: ssDNA -- one chain only
    Suffix("11")
    of3_ssdna = OpenFold3(
        ssDNA=dna_strand
    )

    # 12: dsDNA + ligand
    Suffix("12")
    dna_target = Sequence("AATTAATTAATTAATT", type="dna", ids="DNA_target")
    of3_dna_ligand = OpenFold3(
        dsDNA=dna_target,
        ligands=Ligand("daunorubicin")
    )

    # 13: three axes -- protein x dsDNA x ligand
    Suffix("13")
    of3_three_axis = OpenFold3(
        proteins=Each(protein_a, protein_b),
        dsDNA=dna_target,
        ligands=Each(Ligand("aspirin"), Ligand("caffeine"))
    )

    # 14: RNA
    Suffix("14")
    of3_rna = OpenFold3(
        ssRNA=Sequence("GGCACGUAGCUAGCUAGCUAGCUUGCC", type="rna", ids="RNA_hairpin")
    )

    # =========================================================================
    # Precomputed MSAs
    # =========================================================================

    # 15: OpenFold3 reads alignments as .a3m/.sto (or pre-parsed .npz) -- NOT the `key,sequence`
    # csv that is the Boltz2 convention and MMseqs2's default. So the MSA step has to be asked
    # for a3m, and the server has to be turned off, since a run takes one or the other. Passing
    # a csv msas stream is refused by name rather than folded single-sequence in silence.
    Suffix("15")
    msa_a3m = MMseqs2(sequences=protein_a, output_format="a3m")
    of3_msa = OpenFold3(
        proteins=protein_a,
        msas=msa_a3m,
        use_msa_server=False
    )

    # =========================================================================
    # Sampling, seeds and memory
    # =========================================================================

    # 16: every diffusion sample surfaced as its own structure (<id>_1..N)
    Suffix("16")
    of3_samples = OpenFold3(
        proteins=protein_a,
        num_diffusion_samples=3,
        top_only=False
    )

    # 17: explicit seeds -- the way to make a run reproducible by seed. Passing
    # num_model_seeds as well with a different count is refused rather than silently resolved.
    Suffix("17")
    of3_seeded = OpenFold3(
        proteins=protein_a,
        seeds=[100, 101]
    )

    # 18: low_mem -- the upstream preset, for a GPU that cannot hold the default
    Suffix("18")
    of3_low_mem = OpenFold3(
        proteins=protein_b,
        low_mem=True
    )
