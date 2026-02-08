# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# This software is freely available for use, modification, and redistribution.
# If you use this software or any derivative of it in published work,
# you must cite the original author and this repository.
"""
Test pipeline for refactored tools - Boltz2 tests
Tests: PDB, CompoundLibrary, Boltz2, MMseqs2, Bundle, Each

Combinatorics (Bundle/Each) examples:
- Each: Iterate over elements (cartesian product) - DEFAULT behavior
- Bundle: Group elements together as one entity

Example with 2 proteins (P1, P2) and 3 ligands (L1, L2, L3):
- Default (Each proteins, Each ligands): 6 predictions (P1_L1, P1_L2, P1_L3, P2_L1, P2_L2, P2_L3)
- Bundle(ligands): 2 predictions (P1 with all ligands, P2 with all ligands)
- Bundle(proteins): 3 predictions (all proteins with L1, all proteins with L2, all proteins with L3)
- Bundle(proteins), Bundle(ligands): 1 prediction (all proteins with all ligands)

Nested combinatorics - Bundle containing Each:
- Bundle(Each(library), cofactor): Each ligand from library is bundled with cofactor
  Example with 2 proteins, 3 ligands, 1 cofactor (ATP):
    - 6 predictions: (P1 + L1 + ATP), (P1 + L2 + ATP), (P1 + L3 + ATP),
                     (P2 + L1 + ATP), (P2 + L2 + ATP), (P2 + L3 + ATP)
  Useful for: predicting protein-ligand complexes where a cofactor is always present

IMPORTANT - Affinity calculation order:
  Boltz2 calculates affinity only for the FIRST ligand in a bundle.
  - Bundle(Each(library), ATP): affinity calculated for each library ligand (in presence of ATP)
  - Bundle(ATP, Each(library)): affinity calculated for ATP (in presence of each library ligand)
"""

from PipelineScripts.pipeline import *
from PipelineScripts.entities import *
from PipelineScripts.boltz2 import Boltz2
from PipelineScripts.combinatorics import Bundle, Each

with Pipeline(project="Examples",
              job="Boltz2",
              description="Test Boltz2 with various inputs"):

    Resources(gpu="any",
              time="4:00:00",
              memory="16GB")

    # Test 1: Boltz2 with direct sequence
    Suffix("1") # Outout folder is <Execution order>_Boltz_<Suffix>
    boltz_seq = Boltz2(
        proteins=Sequence("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH")
    )

    # Test 2: PDB fetch and Boltz2 with structure input
    Suffix("2")
    lysozyme = PDB("1AKI", ids="LYZ")
    boltz_pdb = Boltz2(
        proteins=lysozyme
    )

    # Test 3: Boltz2 with ligand (tests compound input)
    Suffix("3")
    boltz_ligand = Boltz2(
        proteins=Sequence("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH"),
        ligands=Ligand("ethanol")
    )

    # Test 4: CompoundLibrary with multiple ligands
    Suffix("4")
    compounds = CompoundLibrary({
        'ethanol': 'CCO',
        'methanol': 'CO',
        'propanol': 'CCCO'
    })
    boltz_multi = Boltz2(
        proteins=lysozyme,
        ligands=compounds
    )

    # Test 5: Boltz2 with MSAs
    Suffix("5")
    boltz_msa = Boltz2(
        proteins=lysozyme,
        msas=boltz_pdb
    )

    # =========================================================================
    # Combinatorics Examples with Bundle and Each
    # =========================================================================

    # Setup: multiple proteins and ligands for combinatorics tests
    protein_a = Sequence("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH", ids="ProteinA")
    protein_b = Sequence("MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAIGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAALINMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAKRVITTFRTGTWDAYKNL", ids="ProteinB")

    ligand_library = CompoundLibrary({
        'aspirin': 'CC(=O)OC1=CC=CC=C1C(=O)O',
        'caffeine': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
        'ibuprofen': 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O'
    })

    # Test 6: Default behavior (Each) - Cartesian product
    Suffix("6")
    # With 2 proteins and 3 ligands: generates 6 predictions
    # ProteinA_aspirin, ProteinA_caffeine, ProteinA_ibuprofen,
    # ProteinB_aspirin, ProteinB_caffeine, ProteinB_ibuprofen
    boltz_each = Boltz2(
        proteins=Each(protein_a, protein_b),  # Each is default, explicit here for clarity
        ligands=ligand_library
    )

    # Test 7: Bundle ligands - Each protein gets all ligands together
    Suffix("7")
    # With 2 proteins and 3 ligands bundled: generates 2 predictions
    # ProteinA (with aspirin+caffeine+ibuprofen), ProteinB (with aspirin+caffeine+ibuprofen)
    boltz_bundle_ligands = Boltz2(
        proteins=Each(protein_a, protein_b),
        ligands=Bundle(ligand_library)
    )

    # Test 8: Bundle proteins - All proteins together with each ligand
    Suffix("8")
    # With 2 proteins bundled and 3 ligands: generates 3 predictions
    # (ProteinA+ProteinB)_aspirin, (ProteinA+ProteinB)_caffeine, (ProteinA+ProteinB)_ibuprofen
    boltz_bundle_proteins = Boltz2(
        proteins=Bundle(protein_a, protein_b),
        ligands=ligand_library
    )

    # Test 9: Bundle both - All proteins with all ligands in one prediction
    Suffix("9")
    # With 2 proteins and 3 ligands both bundled: generates 1 prediction
    # (ProteinA+ProteinB) with (aspirin+caffeine+ibuprofen)
    boltz_bundle_all = Boltz2(
        proteins=Bundle(protein_a, protein_b),
        ligands=Bundle(ligand_library)
    )

    # Test 10: Nested combinatorics - Bundle containing 
    Suffix("10")
    # This is useful when you want to predict each protein-ligand pair with a common cofactor (ATP).
    # Bundle(Each(library), ATP) means: for each ligand in the library, bundle it with ATP.
    # With 1 protein and 3 ligands: generates 3 predictions
    # ProteinA + (aspirin + ATP), ProteinA + (caffeine + ATP), ProteinA + (ibuprofen + ATP)
    #
    # IMPORTANT: Boltz2 calculates affinity only for the FIRST ligand in a bundle.
    # Here, affinity is calculated for aspirin/caffeine/ibuprofen (in presence of ATP).
    atp = Ligand("ATP")
    boltz_nested = Boltz2(
        proteins=protein_a,
        ligands=Bundle(Each(ligand_library), atp)
    )

    # Test 11: Same pattern with multiple proteins - full cartesian product
    Suffix("11")
    # With 2 proteins and Bundle(Each(3 ligands), ATP): generates 6 predictions
    # ProteinA + (aspirin + ATP), ProteinA + (caffeine + ATP), ProteinA + (ibuprofen + ATP),
    # ProteinB + (aspirin + ATP), ProteinB + (caffeine + ATP), ProteinB + (ibuprofen + ATP)
    boltz_nested_multi = Boltz2(
        proteins=Each(protein_a, protein_b),
        ligands=Bundle(Each(ligand_library), atp)
    )

    # Test 12: Reversed order - affinity for cofactor instead of library ligands
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

    print("=== Boltz2 Test Summary ===")
    print(f"boltz_seq structures: {len(boltz_seq.streams.structures.ids)}")
    print(f"boltz_pdb structures: {len(boltz_pdb.streams.structures.ids)}")
    print(f"boltz_ligand structures: {len(boltz_ligand.streams.structures.ids)}")
    print(f"boltz_multi structures: {len(boltz_multi.streams.structures.ids)}")
    print(f"boltz_msa structures: {len(boltz_msa.streams.structures.ids)}")

    print("\n=== Combinatorics Examples ===")
    print(f"boltz_each (Each proteins x Each ligands): {len(boltz_each.streams.structures.ids)} structures")
    print(f"  IDs: {boltz_each.streams.structures.ids}")
    print(f"boltz_bundle_ligands (Each proteins x Bundle ligands): {len(boltz_bundle_ligands.streams.structures.ids)} structures")
    print(f"  IDs: {boltz_bundle_ligands.streams.structures.ids}")
    print(f"boltz_bundle_proteins (Bundle proteins x Each ligands): {len(boltz_bundle_proteins.streams.structures.ids)} structures")
    print(f"  IDs: {boltz_bundle_proteins.streams.structures.ids}")
    print(f"boltz_bundle_all (Bundle proteins x Bundle ligands): {len(boltz_bundle_all.streams.structures.ids)} structures")
    print(f"  IDs: {boltz_bundle_all.streams.structures.ids}")

    print("\n=== Nested Combinatorics ===")
    print(f"boltz_nested (1 protein x Bundle(Each(library), ATP)): {len(boltz_nested.streams.structures.ids)} structures (expected: 3)")
    print(f"  IDs: {boltz_nested.streams.structures.ids}")
    print(f"  Affinity calculated for: library ligands (aspirin, caffeine, ibuprofen)")
    print(f"boltz_nested_multi (2 proteins x Bundle(Each(library), ATP)): {len(boltz_nested_multi.streams.structures.ids)} structures (expected: 6)")
    print(f"  IDs: {boltz_nested_multi.streams.structures.ids}")
    print(f"  Affinity calculated for: library ligands (aspirin, caffeine, ibuprofen)")
    print(f"boltz_atp_affinity (1 protein x Bundle(ATP, Each(library))): {len(boltz_atp_affinity.streams.structures.ids)} structures (expected: 3)")
    print(f"  IDs: {boltz_atp_affinity.streams.structures.ids}")
    print(f"  Affinity calculated for: ATP (in presence of each library ligand)")
