# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
De novo antibiotic binding protein design using BoltzGen.
"""

from biopipelines.pipeline import *
from biopipelines.load import LoadOutput
from biopipelines.boltzgen import BoltzGenMerge, BoltzGen

with Pipeline(project="Examples",
            job="Dopamine_BoltzGen_10x1000designs_AnalysisFiltering",
            description="BoltzGen-based de novo protein binder design against dopamine - analysis and filtering"):   
    Dependencies([str(n) for n in range(FIRST_JOB_ID,FIRST_JOB_ID+10)])
    Resources(gpu="any",time="24:00:00", memory="64GB")
    # Load results from 10 parallel batches
    boltzgens=[]
    for i in range(1,11):
        jsonpath=f"/path/to/Dopamine_BoltzGen_1000designs_{i:03d}/ToolOutputs/002_BoltzGen.json"
        boltzgens.append(LoadOutput(jsonpath,validate_files=False))
    # Run analysis
    analyzed = [BoltzGen(reuse=bg1000,
            steps=["analysis"],
            binder_spec="140-180",  # Binder length range per preprint
            protocol="protein-small_molecule",
            num_designs=1000,
            budget=100,  # Top 100 highest-confidence designs This will matter in the filtering step.
            # RMSD < 2.5Å filter relative to Boltz-2 refolded models (per preprint)
            refolding_rmsd_threshold=2.5,
            # Modified metric weights per preprint
            metrics_override={
                "design_iiptm": 1.1,
                "design_ptm": 1.1,
                "neg_min_design_to_target_pae": 1.1,
                "plip_hbonds_refolded": 2.0,
                "plip_saltbridge_refolded": 2.0,
                "delta_sasa_refolded": 2.0
            }
    ) for bg1000 in boltzgens]
    # We merge and rename design_spec_001 -> batch00_design_spec_001, batch01_design_spec_001, etc.
    merge=BoltzGenMerge(analyzed, id_template="batch{i:02d}_") # this step only works with id_template if analysis has been conducted before (it will merge the analysis results)
    # Final filtering.
    BoltzGen(reuse=merge, # After this is done, one can run the filtering again by loading the merge json.
            steps=["filtering"],
            binder_spec="140-180",  # Binder length range per preprint
            protocol="protein-small_molecule",
            num_designs=1000, # This is irrelevant, it will count how many designs are there
            budget=100,  # Top 100 highest-confidence designs
            # RMSD < 2.5Å filter relative to Boltz-2 refolded models (per preprint)
            refolding_rmsd_threshold=2.5,
            # Modified metric weights per preprint
            metrics_override={
                "design_iiptm": 1.1,
                "design_ptm": 1.1,
                "neg_min_design_to_target_pae": 1.1,
                "plip_hbonds_refolded": 2.0,
                "plip_saltbridge_refolded": 2.0,
                "delta_sasa_refolded": 2.0
            }
    )
