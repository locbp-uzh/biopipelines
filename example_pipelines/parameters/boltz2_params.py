# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Minimal Boltz2 runs exercising every kwarg covered by
tests/tool_parameters/test_boltz2_params.py. Each Suffix block picks one
constraint type and applies it to a small system so the cluster job stays
short. Reviewer 3A: pocket / contacts / glycosylation / covalent / disulfide /
metal_coord are the constraint families surfaced by the wrapper.
"""

from biopipelines.pipeline import *
from biopipelines import Boltz2

with Pipeline(project="ToolParameters",
              job="Boltz2",
              description="Boltz2 parameter coverage — constraints, sampling knobs, output format"):

    Resources(gpu="A100", time="6:00:00", memory="16GB")

    short_seq = Sequence("MKTAYIAKQRQISFVKSHFSRQLEERLGL", ids="short")
    # IgG1 Fc fragment with native N-glycosylation site at Asn-297 (sequence pos 73 here).
    fc = Sequence(
        "TCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPG",
        ids="IgG1_Fc",
    )
    # SnapTag: Cys at residues 5, 24, 145; Asn at 67, 123, 137.
    snap = Sequence(
        "MDKDCEMKRTTLDSPLGKLELSGCEQGLHRIIFLGKGTSAADAVEVPAPAAVLGGPEPLMQATAWLNAYFHQPEAIEEFPVPALHHPVFQQESFTRQVLWKLLKVVKFGEVISYSHLAALAGNPAATAAVKTALSGNPVPILIPCHRVVQGDLDVGGYEGGLAVKEWLLAHEGHRLGKPGLG",
        ids="SnapTag",
    )

    # 1: recycling / diffusion / potentials sampling knobs
    Suffix("sampling")
    Boltz2(
        proteins=short_seq,
        recycling_steps=5,
        diffusion_samples=2,
        use_potentials=True,
    )

    # 2: output format mmcif (msa_server defaults to public)
    Suffix("output_mmcif")
    Boltz2(proteins=short_seq, output_format="mmcif")

    # 3: glycosylation at IgG1 Fc Asn-297 (residue 73 of this fragment)
    Suffix("glycosylation")
    Boltz2(proteins=fc, glycosylation={"A": [73]})

    # 4: pocket constraints around SnapTag's reactive cysteine 145
    Suffix("pocket")
    Boltz2(
        proteins=snap,
        pocket_residues=[143, 144, 145],
        pocket_max_distance=8.0,
        pocket_force=True,
    )

    # 5: contact constraints between two SnapTag residues
    Suffix("contacts")
    Boltz2(
        proteins=snap,
        contacts=[{"token1": ["A", 5], "token2": ["A", 145], "max_distance": 6.0, "force": False}],
    )

    # 6: disulfide bond between SnapTag Cys-5 and Cys-24
    Suffix("disulfide")
    Boltz2(
        proteins=snap,
        disulfide_bonds=[{"token1": ["A", 5], "token2": ["A", 24]}],
    )

    # 7: metal coordination — SnapTag His-29 ND1 coordinating a zinc ion
    Suffix("metal_coord")
    zn = Ligand("ZN", ids="ZN")
    Boltz2(
        proteins=snap,
        ligands=zn,
        metal_coord=[{"atom1": ["A", 29, "ND1"], "atom2": ["B", 1, "ZN"]}],
    )
