"""
Example pipeline: Boltz2 with glycosylation and covalent bonding.

Demonstrates:
1. N-glycosylation on an antibody Fc fragment (NAG at Asn-297)
2. Covalent inhibitor bound to a cysteine protease (e.g., to Cys-25)
3. Combined: glycosylated protein with a covalent ligand
4. Contact constraints to guide ligand placement
"""

import os, sys
sys.path.insert(0, os.getcwd())

from PipelineScripts.pipeline import *
from PipelineScripts.entities import *
from PipelineScripts.boltz2 import Boltz2

with Pipeline(project="Examples",
              job="BoltzCovalent",
              description="Boltz2 glycosylation and covalent bonding examples"):

    Resources(gpu="any",
              time="4:00:00",
              memory="16GB")

    # =========================================================================
    # Example 1: N-glycosylation on IgG1 Fc (Asn-297)
    # =========================================================================
    # The IgG1 Fc fragment is glycosylated at Asn-297 on each chain.
    # Boltz2 adds an NAG (N-acetylglucosamine) moiety at the specified positions.

    Suffix("glyco")
    fc_fragment = Sequence(
        "TCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPG",
        ids="IgG1_Fc"
    )

    boltz_glyco = Boltz2(
        proteins=fc_fragment,
        glycosylation={"A": [73]}
    )

    # =========================================================================
    # Example 2: Benzylated SNAP-Tag
    # =========================================================================
    # Cys-145 is the reactive cysteine

    Suffix("covalent")
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

    # =========================================================================
    # Example 3: Glycosylated protein with a covalent ligand
    # =========================================================================

    Suffix("glyco_covalent")
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

    # =========================================================================
    # Example 4: Contact constraints
    # =========================================================================
    # Guide ligand placement by constraining specific residues to be near
    # the ligand. token1/token2 use [chain, residue] or [chain, residue, atom].
    # max_distance defaults to 6A, force defaults to False.

    Suffix("contact")
    boltz_contact = Boltz2(
        proteins=SnapTag,
        ligands=Bn,
        msas=boltz_inspection,
        contacts=[
            {"token1": ["A", 145], "token2": ["B", 1], "max_distance": 8.0, "force": False}
        ]
    )

    print("=== Boltz2 Covalent & Glycosylation Examples ===")
    print(f"Example 1 (glycosylation): {len(boltz_glyco.streams.structures.ids)} structures")
    print(f"  IDs: {boltz_glyco.streams.structures.ids}")
    print(f"Example 2 (covalent):      {len(boltz_covalent.streams.structures.ids)} structures")
    print(f"  IDs: {boltz_covalent.streams.structures.ids}")
    print(f"Example 3 (both):          {len(boltz_glyco_covalent.streams.structures.ids)} structures")
    print(f"  IDs: {boltz_glyco_covalent.streams.structures.ids}")
    print(f"Example 4 (contacts):      {len(boltz_contact.streams.structures.ids)} structures")
    print(f"  IDs: {boltz_contact.streams.structures.ids}")
