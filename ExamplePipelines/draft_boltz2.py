# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

from biopipelines.pipeline import *
from biopipelines.boltz2 import Boltz2
from biopipelines.combinatorics import Bundle, Each

with Pipeline(project="Examples",
              job="Boltz2",
              description="Boltz2 with various inputs"):

    Resources(gpu="any",
              time="4:00:00",
              memory="16GB")

    # =========================================================================
    # Basic usage examples
    # =========================================================================

    # 1: Boltz2 with direct sequence
    Suffix("1") # Outout folder is <Execution order>_Boltz_<Suffix>
    boltz_seq = Boltz2(
        proteins=Sequence("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH")
    )

    # 2: PDB fetch and Boltz2 with structure input
    Suffix("2")
    lysozyme = PDB("1AKI", ids="LYZ")
    boltz_pdb = Boltz2(
        proteins=lysozyme
    )

    # 3: Boltz2 with recycled MSAs
    Suffix("3")
    boltz_msa = Boltz2(
        proteins=lysozyme,
        msas=boltz_pdb
    )

    # 4: Boltz2 with ligand (tests compound input)
    Suffix("4")
    boltz_ligand = Boltz2(
        proteins=Sequence("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH"),
        ligands=Ligand("ethanol")
    )

    # 5: CompoundLibrary with multiple ligands
    Suffix("5")
    compounds = CompoundLibrary({
        'ethanol': 'CCO',
        'methanol': 'CO',
        'propanol': 'CCCO'
    })
    boltz_multi = Boltz2(
        proteins=lysozyme,
        ligands=compounds
    )

    # =========================================================================
    # Combinatorics
    # =========================================================================

    # Setup: multiple proteins and ligands for combinatorics tests
    protein_a = Sequence("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH", ids="ProteinA")
    protein_b = Sequence("MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAIGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAALINMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAKRVITTFRTGTWDAYKNL", ids="ProteinB")

    ligand_library = Ligand(['aspirin',
                             'caffeine',
                             'ibuprofen'])

    # 6: Default behavior (Each) - Cartesian product
    Suffix("6")
    # With 2 proteins and 3 ligands: generates 6 predictions
    # ProteinA_aspirin, ProteinA_caffeine, ProteinA_ibuprofen,
    # ProteinB_aspirin, ProteinB_caffeine, ProteinB_ibuprofen
    boltz_each = Boltz2(
        proteins=Each(protein_a, protein_b),  # Each is default, explicit here for clarity
        ligands=ligand_library
    )

    # 7: Bundle ligands - Each protein gets all ligands together
    Suffix("7")
    # With 2 proteins and 3 ligands bundled: generates 2 predictions
    # ProteinA (with aspirin+caffeine+ibuprofen), ProteinB (with aspirin+caffeine+ibuprofen)
    boltz_bundle_ligands = Boltz2(
        proteins=Each(protein_a, protein_b),
        ligands=Bundle(ligand_library)
    )

    # 8: Bundle proteins - All proteins together with each ligand
    Suffix("8")
    # With 2 proteins bundled and 3 ligands: generates 3 predictions
    # (ProteinA+ProteinB)_aspirin, (ProteinA+ProteinB)_caffeine, (ProteinA+ProteinB)_ibuprofen
    boltz_bundle_proteins = Boltz2(
        proteins=Bundle(protein_a, protein_b),
        ligands=ligand_library
    )

    # 9: Bundle both - All proteins with all ligands in one prediction
    Suffix("9")
    # With 2 proteins and 3 ligands both bundled: generates 1 prediction
    # (ProteinA+ProteinB) with (aspirin+caffeine+ibuprofen)
    boltz_bundle_all = Boltz2(
        proteins=Bundle(protein_a, protein_b),
        ligands=Bundle(ligand_library)
    )

    # 10: Nested combinatorics - Bundle containing 
    Suffix("10")
    # This is useful when you want to predict each protein-ligand pair with a common cofactor (ATP).
    # Bundle(Each(library), ATP) means: for each ligand in the library, bundle it with ATP.
    # With 1 protein and 3 ligands: generates 3 predictions
    # ProteinA + (aspirin + ATP), ProteinA + (caffeine + ATP), ProteinA + (ibuprofen + ATP)
    #
    # IMPORTANT: Boltz2 calculates affinity only for the FIRST ligand in a bundle.
    # Here, affinity is calculated for aspirin/caffeine/ibuprofen (in presence of ATP).
    # Also consider: if two ligands are identical, set affinity=False (calculation not supported by Boltz)
    atp = Ligand("ATP")
    boltz_nested = Boltz2(
        proteins=protein_a,
        ligands=Bundle(Each(ligand_library), atp)
    )

    # 11: Same pattern with multiple proteins - full cartesian product
    Suffix("11")
    # With 2 proteins and Bundle(Each(3 ligands), ATP): generates 6 predictions
    # ProteinA + (aspirin + ATP), ProteinA + (caffeine + ATP), ProteinA + (ibuprofen + ATP),
    # ProteinB + (aspirin + ATP), ProteinB + (caffeine + ATP), ProteinB + (ibuprofen + ATP)
    boltz_nested_multi = Boltz2(
        proteins=Each(protein_a, protein_b),
        ligands=Bundle(Each(ligand_library), atp)
    )

    # 12: Reversed order - affinity for cofactor instead of library ligands
    Suffix("12")
    # Bundle(ATP, Each(library)) means: ATP is first, so affinity is calculated for ATP
    # in the presence of each library ligand.
    # With 1 protein and 3 ligands: generates 3 predictions
    # ProteinA + (ATP + aspirin), ProteinA + (ATP + caffeine), ProteinA + (ATP + ibuprofen)
    # Affinity is calculated for ATP in each case.
    boltz_atp_affinity = Boltz2(
        proteins=protein_a,
        ligands=Bundle(atp, Each(ligand_library))
    )

    # =========================================================================
    # Covalent modifications and constraints
    # =========================================================================

    # 13: N-glycosylation on IgG1 Fc (Asn-297) - with NAG (N-acetylglucosamine)
    Suffix("13")
    fc_fragment = Sequence(
        "TCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPG",
        ids="IgG1_Fc"
    )
    boltz_glyco = Boltz2(
        proteins=fc_fragment,
        glycosylation={"A": [73]}
    )

    # 14: Benzylated SNAP-Tag - Cys-145 is the reactive cysteine
    Suffix("14")
    SnapTag = Sequence("MDKDCEMKRTTLDSPLGKLELSGCEQGLHRIIFLGKGTSAADAVEVPAPAAVLGGPEPLMQATAWLNAYFHQPEAIEEFPVPALHHPVFQQESFTRQVLWKLLKVVKFGEVISYSHLAALAGNPAATAAVKTALSGNPVPILIPCHRVVQGDLDVGGYEGGLAVKEWLLAHEGHRLGKPGLG")
    Bn = Ligand(smiles="CC1=CC=CC=C1")

    # Run without the covalent linkage and inspect to verify the atom number
    boltz_inspection = Boltz2(
        proteins=SnapTag,
        ligands=Bn
    )

    boltz_covalent = Boltz2(
        proteins=SnapTag,
        ligands=Bn,
        msas=boltz_inspection,
        covalent_linkage={
            "chain": "A",
            "position": 145,
            "protein_atom": "SG",
            "ligand_atom": "C15"
        }
    )

    # 15: Glycosylated protein with a covalent ligand
    Suffix("15")
    boltz_glyco_covalent = Boltz2(
        proteins=SnapTag,
        ligands=Bn,
        msas=boltz_covalent,
        glycosylation={"A": [67]},
        covalent_linkage={
            "chain": "A",
            "position": 145,
            "protein_atom": "SG",
            "ligand_atom": "C15" # Run without the covalent linkage and inspect to verify the atom number
        }
    )

    # 16: Contact constraints - specific residues to be near the ligand
    # token1/token2 use [chain, residue_index] for proteins and [chain, atom_name] for ligands.
    Suffix("16")
    boltz_contact = Boltz2(
        proteins=SnapTag,
        ligands=Bn,
        msas=boltz_inspection,
        contacts=[
            {"token1": ["A", 145], "token2": ["B", "C15"], "max_distance": 8.0, "force": False}
        ]
    )

    # =========================================================================
    # DNA examples
    # =========================================================================

    # 17: DNA only - predict structure of a DNA duplex
    Suffix("17")
    dna_strand = Sequence("ACGTACGTACGTACGT", type="dna", ids="DNA_strand")
    boltz_dna = Boltz2(
        dna=dna_strand
    )

    # 18: DNA + ligand - predict binding of a small molecule to DNA
    Suffix("18")
    dna_target = Sequence("AATTAATTAATTAATT", type="dna", ids="DNA_target")
    daunorubicin = Ligand("daunorubicin")
    boltz_dna_ligand = Boltz2(
        dna=dna_target,
        ligands=daunorubicin
    )

    # 19: Protein + DNA + ligand - three-axis combinatorics
    Suffix("19")
    # With 2 proteins, 1 DNA, and 2 ligands (all Each): generates 4 predictions
    # ProteinA_DNA_target_aspirin, ProteinA_DNA_target_caffeine,
    # ProteinB_DNA_target_aspirin, ProteinB_DNA_target_caffeine
    boltz_three_axis = Boltz2(
        proteins=Each(protein_a, protein_b),
        dna=dna_target,
        ligands=Each(Ligand("aspirin"), Ligand("caffeine"))
    )
