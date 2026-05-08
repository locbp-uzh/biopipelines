# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Minimal BoltzGen runs exercising every kwarg covered by
tests/tool_parameters/test_boltzgen_params.py. Each Suffix block keeps
num_designs and budget tiny so the cluster job stays cheap.
"""

from biopipelines.pipeline import *
from biopipelines.boltzgen import BoltzGen

with Pipeline(project="ToolParameters",
              job="BoltzGen",
              description="BoltzGen parameter coverage — protocol/num_designs/budget/sampling/filter knobs"):

    Resources(gpu="80GB|96GB", time="6:00:00", memory="16GB")

    target = PDB("1AKI", ids="LYZ")  # small protein target

    # 1: protocol + sampling knobs (step_scale, noise_scale, batch size)
    Suffix("sampling")
    BoltzGen(
        target_structure=target,
        binder_spec="80",
        protocol="protein-anything",
        num_designs=20,
        budget=4,
        step_scale=1.25,
        noise_scale=0.6,
        diffusion_batch_size=8,
    )

    # 2: inverse-folding knobs
    Suffix("inverse_fold")
    BoltzGen(
        target_structure=target,
        binder_spec="80",
        protocol="protein-anything",
        num_designs=20,
        budget=4,
        inverse_fold_num_sequences=4,
    )

    # 3: skip inverse folding entirely
    Suffix("skip_if")
    BoltzGen(
        target_structure=target,
        binder_spec="80",
        protocol="protein-anything",
        num_designs=20,
        budget=4,
        skip_inverse_folding=True,
    )

    # 4: filter knobs (alpha, refolding RMSD, filter_biased off)
    Suffix("filters")
    BoltzGen(
        target_structure=target,
        binder_spec="80",
        protocol="protein-anything",
        num_designs=20,
        budget=4,
        alpha=0.7,
        refolding_rmsd_threshold=2.5,
        filter_biased=False,
    )
