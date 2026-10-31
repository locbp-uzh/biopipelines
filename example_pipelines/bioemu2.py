# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

# tested:

"""BioEmu2 across input shapes, combinatorics and constraints.

The axes and the combinatorics are the ones `Boltz2` and `Protenix` take: a run is one system, `Each` (the default) makes one run per sequence and ligand, and `Bundle` puts several rows in one system.
The conformers of a run are `<run>_1` to `<run>_<num_samples>`, and the `physicality` table says which of them upstream's checks rejected.
"""

from biopipelines.pipeline import *
from biopipelines import BioEmu2, Boltz2, Bundle, Each

with Pipeline(project="Examples",
              job="BioEmu2",
              description="BioEmu2 equilibrium ensembles: apo, holo, recycled alignments, combinatorics and constraints"):

    Resources(gpu="A100",
              time="04:00:00",
              memory="48GB")

    hbb = Sequence("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH", ids="HBB_fragment")

    # =========================================================================
    # Basic usage
    # =========================================================================

    # 1: an apo ensemble from a sequence. Without msas, bioemu2 searches ColabFold for the alignment
    Suffix("1")
    bioemu_apo = BioEmu2(
        sequences=hbb,
        num_samples=50
    )

    # 2: a holo ensemble with ATP. A ligand with a CCD code is passed to bioemu2 by that code, and placed by it from scratch
    Suffix("2")
    bioemu_holo = BioEmu2(
        sequences=hbb,
        ligands=Ligand("ATP"),
        num_samples=50,
        seed=1
    )

    # 3: reusing an alignment. Boltz2 writes a csv alignment per sequence, which bioemu2 is given in place of a search
    Suffix("3")
    boltz_alignment = Boltz2(proteins=hbb)
    bioemu_recycled = BioEmu2(
        sequences=hbb,
        msas=boltz_alignment,
        num_samples=50
    )

    # =========================================================================
    # Combinatorics
    # =========================================================================

    hbb_variant = Sequence("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSA", ids="HBB_variant")
    cofactors = Ligand(["ATP", "ADP"])

    # 4: Each x Each. Two sequences and two ligands are four runs: HBB_fragment+ATP, HBB_fragment+ADP, HBB_variant+ATP, HBB_variant+ADP
    Suffix("4")
    bioemu_each = BioEmu2(
        sequences=Each(hbb, hbb_variant),
        ligands=cofactors,
        num_samples=20
    )

    # 5: Bundle puts both ligands in one system with each sequence: two runs, ATP+ADP+HBB_fragment and ATP+ADP+HBB_variant, with chain A the protein and B and C the ligands
    Suffix("5")
    bioemu_bundle = BioEmu2(
        sequences=Each(hbb, hbb_variant),
        ligands=Bundle(cofactors),
        num_samples=20
    )

    # =========================================================================
    # Constraints
    # =========================================================================

    # 6: a covalent bond between a protein atom and an atom of the first ligand, here the reactive cysteine of SnapTag
    # and the anomeric carbon of NAG. The bond names the ligand's residue, so the ligand has to be given by CCD code.
    # The pair is illustrative: a real adduct needs its own code and atom names, so check them on an unconstrained run first.
    Suffix("6")
    snaptag = Sequence(
        "MDKDCEMKRTTLDSPLGKLELSGCEQGLHRIIFLGKGTSAADAVEVPAPAAVLGGPEPLMQATAWLNAYFHQPEAIEEFPVPALHHPVFQQESFTRQVLWKLLKVVKFGEVISYSHLAALAGNPAATAAVKTALSGNPVPILIPCHRVVQGDLDVGGYEGGLAVKEWLLAHEGHRLGKPGLG",
        ids="SnapTag"
    )
    bioemu_covalent = BioEmu2(
        sequences=snaptag,
        ligands=Ligand("NAG"),
        num_samples=50,
        covalent_linkage={
            "chain": "A",
            "position": 145,
            "protein_atom": "SG",
            "ligand_atom": "C1"
        }
    )

    # 7: a restraint under a steered sampler, between residue 3 of the protein (chain A) and the gamma phosphorus PG of ATP (chain B)
    Suffix("7")
    bioemu_contact = BioEmu2(
        sequences=hbb,
        ligands=Ligand("ATP"),
        sampler="edm_gd",
        num_samples=50,
        contacts=[{"token1": {"asym_id": "A", "res_id": 3},
                   "token2": {"asym_id": "B", "atom_name": "PG"},
                   "max_distance": 6.0}]
    )
