"""
This pipeline shows how to use RF3 (RoseTTAFold3)
"""

import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import Pipeline
from PipelineScripts.rf3 import RF3
from PipelineScripts.pdb import PDB
from PipelineScripts.compound_library import CompoundLibrary

pipeline = Pipeline(
    pipeline_name="Examples",
    job_name="RF3",
    job_description="RF3 structure prediction examples")

pipeline.resources(
    gpu="A100",
    time="24:00:00",
    memory="16GB"
)

"""
RF3 run with direct protein sequence and ligand SMILES
"""
rf3_1 = pipeline.add(RF3(
    proteins="MAEIGTGFPFDPHYVEVLGERMHYVDVGPRDGTPVLFLHGNPTSSYVWRNIIPHVAPTHRCIAPDLIGMGKSDKPDLGYFFDDHVRFMDAFIEALGLEEVVLVIHDWGSALGFHWAKRNPERVKGIAFMEFIRPIPTWDEWPEFARETFQAFRTTDVGRKLIIDQNVFIEGTLPMGVVRPLTEVEMDHYREPFLNPVDREPLWRFPNELPIAGEPANIVALVEEYMDWLHQSPVPKLLFWGTPGVLIPPAEAARLAKSLPNCKAVDIGPGLNLLQEDNPDLIGSEIARWLSTLEISG",
    ligands=r"CC/1(C)C2=C(C=CC=C2)N(C)\C1=C\C=C\C=C\C=C\[C@@]34[C@@](CC5=CN(CCOCCOCCCCCCCl)N=N5)(CC(=O)N3C)C6=C(C=CC=C6)N4C",
    num_models=3
))

"""
RF3 apo prediction from PDB structure
"""
# Load HaloTag from PDB
HaloTag = pipeline.add(PDB("6U32", "HT"))

rf3_2 = pipeline.add(RF3(
    proteins=HaloTag,
    num_models=5
))

"""
RF3 holo prediction using compounds from previous RF3 run
"""
rf3_3 = pipeline.add(RF3(
    proteins=HaloTag,
    ligands=rf3_1,  # Output of RF3 contains CSV with ligand SMILES
    num_models=3
))

"""
RF3 run with multiple ligands from compound library
"""
compounds = pipeline.add(CompoundLibrary(
    {
        'S_Cy7-C-RCG+': r'CC/1(C)C2=C(C=CC=C2)N(C)\C1=C\C=C\C=C\C=C\C3=[N+](C)C4=C(C=CC=C4)[C@]3(CC5=CN(CCOCCOCCCCCCCl)N=N5)CC(=O)NC',
        'S_Cy7-C-CAP+': r'CC1(C)C2=C(C=CC=C2)[N+](=C1/C=C/C=C/C=C/C=C\3/[C@@](CC4=CN(CCOCCOCCCCCCCl)N=N4)(CC(=O)NC)C5=C(C=CC=C5)N3C)C',
        'R_Cy7-C-RCG+': r'CC/1(C)C2=C(C=CC=C2)N(C)\C1=C\C=C\C=C\C=C\C3=[N+](C)C4=C(C=CC=C4)[C@@]3(CC5=CN(CCOCCOCCCCCCCl)N=N5)CC(=O)NC',
        'R_Cy7-C-CAP+': r'CC1(C)C2=C(C=CC=C2)[N+](=C1/C=C/C=C/C=C/C=C\3/[C@](CC4=CN(CCOCCOCCCCCCCl)N=N4)(CC(=O)NC)C5=C(C=CC=C5)N3C)C'
    }
))

rf3_4 = pipeline.add(RF3(
    proteins=HaloTag,
    ligands=compounds,
    num_models=3,
    early_stopping_plddt=85.0
))

"""
RF3 with expandable ligand library
"""
expanded_library = pipeline.add(CompoundLibrary(
    {
        'CyC': r'<Linker><Cy7MethylAmide>',
        'Linker': ['<OriginalLinker>', '<AmideLinker>'],
        'OriginalLinker': 'ClCCCCCCOCCOCC',
        'AmideLinker': '<ChloroCarbonyl><AminoPEG>',
        'ChloroCarbonyl': ['ClCCC(=O)', 'ClCCCC(=O)', 'ClCCCCC(=O)', 'ClCCCCCC(=O)', 'ClCCCCCCC(=O)'],
        'AminoPEG': ['NCCOCC', 'NCCOCCOCC'],
        'Cy7MethylAmide': r'N%91C=C(C[C<Configuration>]%92(CC(NC)=O)C(/C=C/C=C/C=C/C=C(N(C%95=C%96C=CC=C%95)C)C%96(C)C)=[N+](C)C%93=CC=CC=C%93%92)N=N%91',
        'Configuration': '@'  # Only R enantiomer
    },
    primary_key='CyC'
))

rf3_5 = pipeline.add(RF3(
    proteins=HaloTag,
    ligands=expanded_library,
    num_models=2
))

pipeline.slurm()
