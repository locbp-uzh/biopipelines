# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

# tested:

"""Protenix with the same inputs as the Boltz2 and OpenFold3 examples, so the three can be compared directly.

The axes and the combinatorics are the ones `Boltz2` takes, so a campaign written for one runs on the others by changing the tool name.
Where a case from `boltz2.py` is missing here it is because Protenix does not do it: there is no affinity head, and the contact and pocket constraints exist only for an older checkpoint.
"""

from biopipelines.pipeline import *
from biopipelines import Protenix, MMseqs2, IPSAE, Bundle, Each

with Pipeline(project="Examples",
              job="Protenix",
              description="Protenix across input shapes, combinatorics and options"):

    Resources(gpu="A100",
              time="24:00:00",
              memory="32GB")

    protein_a = Sequence("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH", ids="ProteinA")
    protein_b = Sequence("MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAIGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAALINMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAKRVITTFRTGTWDAYKNL", ids="ProteinB")
    ligand_library = Ligand(['aspirin', 'caffeine', 'ibuprofen'])

    # 1: a protein alone. Without msas the alignment is searched through upstream's MSA service
    Suffix("1")
    px_seq = Protenix(proteins=protein_a)

    # 2: protein + ligand
    Suffix("2")
    px_ligand = Protenix(proteins=protein_a, ligands=Ligand("ethanol"))

    # 3: Each x Each -- 2 proteins x 3 ligands = 6 predictions
    Suffix("3")
    px_each = Protenix(proteins=Each(protein_a, protein_b), ligands=ligand_library)

    # 4: bundled proteins -- both chains with each ligand, 3 predictions
    Suffix("4")
    px_bundle = Protenix(proteins=Bundle(protein_a, protein_b), ligands=ligand_library)

    # 5: protein + double-stranded DNA; the reverse-complement strand is generated as a second chain
    Suffix("5")
    px_dna = Protenix(proteins=protein_a,
                      dsDNA=Sequence("AATTAATTAATTAATT", type="dna", ids="DNA_target"))

    # 6: precomputed alignments. Protenix reads .a3m, so the MSA step is asked for a3m; the service is
    # then off, and a protein chain with no alignment skips its query instead of reaching the server
    Suffix("6")
    msa_a3m = MMseqs2(sequences=protein_a, output_format="a3m")
    px_msa = Protenix(proteins=protein_a, msas=msa_a3m)

    # 7: no alignment at all, every chain folded from its single sequence
    Suffix("7")
    px_single = Protenix(proteins=protein_a, use_msa_server=False)

    # 8: five seeds of three samples, every structure surfaced as <id>_1..15 instead of only the best
    Suffix("8")
    px_samples = Protenix(proteins=protein_a,
                          seeds=[1, 2, 3, 4, 5],
                          num_diffusion_samples=3,
                          top_only=False)

    # 9: another checkpoint, here upstream's previous generation base model
    Suffix("9")
    px_v1 = Protenix(proteins=protein_a, model_name="protenix_base_default_v1.0.0")

    # 10: a covalent bond between entity 1 (the protein, residue 5, atom SG) and entity 2 (the ligand, CCD atom C1)
    Suffix("10")
    px_covalent = Protenix(proteins=protein_a,
                           ligands=Ligand("NAG"),
                           covalent_bonds=[{"entity1": 1, "position1": 5, "atom1": "SG",
                                            "entity2": 2, "position2": 1, "atom2": "C1"}])

    # 11: the pae stream has the Boltz2 layout, which IPSAE reads
    Suffix("11")
    px_dimer = Protenix(proteins=Bundle(protein_a, protein_b))
    interface = IPSAE(structures=px_dimer, pae=px_dimer)
