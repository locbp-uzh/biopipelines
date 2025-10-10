"""
This pipeline shows how to use boltz2
"""

import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import Pipeline
from PipelineScripts.boltz2 import Boltz2
from PipelineScripts.pdb import PDB
from PipelineScripts.compound_library import CompoundLibrary
from PipelineScripts.mmseqs2 import MMseqs2

pipeline = Pipeline(
    pipeline_name="Examples",
    job_name="Boltz",
    job_description="Folding of HaloTag7 with Cy7 methyl amide close enantiomer SS")

pipeline.resources(
    gpu="V100",
    time="24:00:00",
    memory="16GB"
)

"""
Boltz run with yaml config
"""
config = """
sequences:
- protein:
    id: A
    sequence: MAEIGTGFPFDPHYVEVLGERMHYVDVGPRDGTPVLFLHGNPTSSYVWRNIIPHVAPTHRCIAPDLIGMGKSDKPDLGYFFDDHVRFMDAFIEALGLEEVVLVIHDWGSALGFHWAKRNPERVKGIAFMEFIRPIPTWDEWPEFARETFQAFRTTDVGRKLIIDQNVFIEGTLPMGVVRPLTEVEMDHYREPFLNPVDREPLWRFPNELPIAGEPANIVALVEEYMDWLHQSPVPKLLFWGTPGVLIPPAEAARLAKSLPNCKAVDIGPGLNLLQEDNPDLIGSEIARWLSTLEISG
- ligand:
    id: B
    smiles: C[N+](C1=C2C=CC=C1)=C([C@]2(CC3=CN(N=N3)CCOCCOCCCCCCCl)CC(NC)=O)/C=C/C=C/C=C(N(C4=C5C=CC=C4)C)\C5(C)C
properties:
- affinity:
    binder: B
"""
boltz_0 = pipeline.add(Boltz2(config=config))

"""
Boltz run with provided sequence and ligand smiles
"""

boltz_1 = pipeline.add(Boltz2(proteins="MAEIGTGFPFDPHYVEVLGERMHYVDVGPRDGTPVLFLHGNPTSSYVWRNIIPHVAPTHRCIAPDLIGMGKSDKPDLGYFFDDHVRFMDAFIEALGLEEVVLVIHDWGSALGFHWAKRNPERVKGIAFMEFIRPIPTWDEWPEFARETFQAFRTTDVGRKLIIDQNVFIEGTLPMGVVRPLTEVEMDHYREPFLNPVDREPLWRFPNELPIAGEPANIVALVEEYMDWLHQSPVPKLLFWGTPGVLIPPAEAARLAKSLPNCKAVDIGPGLNLLQEDNPDLIGSEIARWLSTLEISG",
ligands=r"CC/1(C)C2=C(C=CC=C2)N(C)\C1=C\C=C\C=C\C=C\[C@@]34[C@@](CC5=CN(CCOCCOCCCCCCCl)N=N5)(CC(=O)N3C)C6=C(C=CC=C6)N4C", #use r prefix to avoid misinterpretation of \
affinity=True, # this is set to True as default when there is a ligand, and false otherwise. So it only makes sense to override when there is a ligand, and you don't want the affinity
global_msas_cache=True,
pipeline=pipeline))# We add pipeline=pipeline only so the tool can access the job name. This is not done when folding sequences from previous tools.


"""
Boltz run with protein sequence fetched from pdb
"""

# We load pdb id 6U32 with structure id HT
HaloTag = pipeline.add(PDB("6U32","HT"))

boltz_2 = pipeline.add(Boltz2(proteins=HaloTag,
ligands=boltz_1, # output of Boltz contains a csv file with the smiles of the ligands used.
global_msas_cache=True,
pipeline=pipeline)) 

"""
Boltz run with multiple ligands
"""
compounds = pipeline.add(CompoundLibrary(
    {
        'S_Cy7-C-RCG+':r'CC/1(C)C2=C(C=CC=C2)N(C)\C1=C\C=C\C=C\C=C\C3=[N+](C)C4=C(C=CC=C4)[C@]3(CC5=CN(CCOCCOCCCCCCCl)N=N5)CC(=O)NC',
        'S_Cy7-C-CAP+':r'CC1(C)C2=C(C=CC=C2)[N+](=C1/C=C/C=C/C=C/C=C\3/[C@@](CC4=CN(CCOCCOCCCCCCCl)N=N4)(CC(=O)NC)C5=C(C=CC=C5)N3C)C',
        'R_Cy7-C-RCG+':r'CC/1(C)C2=C(C=CC=C2)N(C)\C1=C\C=C\C=C\C=C\C3=[N+](C)C4=C(C=CC=C4)[C@@]3(CC5=CN(CCOCCOCCCCCCCl)N=N5)CC(=O)NC',
        'R_Cy7-C-CAP+':r'CC1(C)C2=C(C=CC=C2)[N+](=C1/C=C/C=C/C=C/C=C\3/[C@](CC4=CN(CCOCCOCCCCCCCl)N=N4)(CC(=O)NC)C5=C(C=CC=C5)N3C)C'
    }
))

boltz_3 = pipeline.add(Boltz2(proteins=HaloTag,
ligands=compounds,
global_msas_cache=True,
pipeline=pipeline))

"""
Boltz with an expandable ligand library
"""
expanded_library = pipeline.add(CompoundLibrary(
    {
        'CyC':r'<Linker><Cy7MethylAmide>',
        'Linker': ['<OriginalLinker>','<AmideLinker>'],
        'OriginalLinker': 'ClCCCCCCOCCOCC',
        'AmideLinker': '<ChloroCarbonyl><AminoPEG>',
        'ChloroCarbonyl': ['ClCCC(=O)','ClCCCC(=O)','ClCCCCC(=O)','ClCCCCCC(=O)','ClCCCCCCC(=O)'],
        'AminoPEG': ['NCCOCC','NCCOCCOCC'],
        'Cy7MethylAmide': r'N%91C=C(C[C<Configuration>]%92(CC(NC)=O)C(/C=C/C=C/C=C/C=C(N(C%95=C%96C=CC=C%95)C)C%96(C)C)=[N+](C)C%93=CC=CC=C%93%92)N=N%91',
        'Configuration': '@' #Only R enantiomer
    },
    primary_key='CyC'
))

boltz_4 = pipeline.add(Boltz2(proteins=HaloTag,
ligands=expanded_library,
global_msas_cache=True,
pipeline=pipeline))

"""
Boltz with msas from MMseqs2
"""

msa = pipeline.add(MMseqs2(sequences=HaloTag))
boltz_5 = pipeline.add(Boltz2(proteins=HaloTag,msas=msa))


pipeline.slurm()