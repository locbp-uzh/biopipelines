"""
This pipeline shows how to improve the difference in predicted binding affinity between open and close form of a carbocyanine 7 chloride and halotag7 starting from a Boltz model of the open form.
"""

#!/usr/bin/env python3
import os, sys
sys.path.insert(0, os.getcwd()) #to see scripts in current folder

from PipelineScripts.pipeline import Pipeline
from PipelineScripts.load_output import LoadOutput
from PipelineScripts.distance_selector import DistanceSelector
from PipelineScripts.ligand_mpnn import LigandMPNN
from PipelineScripts.boltz2 import Boltz2
from PipelineScripts.residue_atom_distance import ResidueAtomDistance
from PipelineScripts.merge_datasheets import MergeDatasheets
from PipelineScripts.filter import Filter

pipeline = Pipeline(
    pipeline_name="LigandMPNN-Boltz", #Will create a folder in /shares/USER/<pipeline_name>
    job_name="HT7_Cy7_ChlorineFilter", #Unique job folder in /shares/USER/<pipeline_name>/job_name_NNN
    job_description="Test on filter based on distance between chlorine and aspartate")

pipeline.resources(
    gpu="V100",
    time="24:00:00",
    memory="16GB"
)

"""
It is best practise to start from a Boltz2 output with the open form, to have a benchmark affinity.
One can then load it with the LoadOutput tool, which will contain the same structures (pdbs), ids, and datasheets as the Boltz2 tool of the past pipeline.  
"""
original = pipeline.add(LoadOutput(
    '/shares/locbp.chem.uzh/public/BioPipelines/Boltz/HT7_Cy7_C_R_001/ToolOutputs/1_Boltz2_output.json'
    #'path/to/job/ToolOutputs/<Job>_Boltz2_output.json'
))

"""
Generate distance-based residue selections
This new approach replaces the previous design_within post-processing with explicit
distance analysis that generates PyMOL-formatted selections for each structure.
"""
distance_analysis = pipeline.add(DistanceSelector(structures=original,
                                                       ligand="LIG",
                                                       distance=4.0))

"""
Diversify with LigandMPNN using distance-based selections
"""
lmpnn = pipeline.add(LigandMPNN(structures=original, #this is equivalent to boltz2
                                ligand="LIG", #in ligand mpnn you should always specify the ligand name, which is LIG if from Boltz
                                num_sequences=3,
                                redesigned=(distance_analysis.datasheets.selections, "within") #use residues within 4Å of ligand
))

"""
We run the Apo version first. One can also extract confidence parameters from it, and in general here is where we calculate the MSAs, which will be recycled later on with the msas input parameter.
"""
boltz_apo = pipeline.add(Boltz2(proteins=lmpnn))

"""
Run with open form and closed form.
Important: msas are passed with <tool>, not with <tool>.msas
"""
boltz_holo_open = pipeline.add(Boltz2(proteins=lmpnn,
                                ligands=original,
                                msas=boltz_apo,
                                affinity=True))

"""
Now we want to calculate the distance between the active aspartate and the chlorine in the ligand. This can be done with an analysis tool called ResidueAtomDistance. The residue can be specified with it position e.g. residue="106" or with its context e.g. "D in IHDWG". The context is useful if the analysis follows e.g. RFdiffusion with non-fixed length contigs, but make sure LigandMPNN is not mutating it. The ligand name has to be specified, which is always LIG when using Boltz for folding. The output will be a datasheet containing columns id, <metric_name>, where id comes from boltz and therefore is associated to a specific protein sequence.
"""
open_chlorine_aspartate_distance = pipeline.add(ResidueAtomDistance(input=boltz_holo_open,
                                                                    residue='D in IHDWG',
                                                                    atom='LIG.Cl',
                                                                    metric_name='open_chlorine_distance'))

"""
Now we merge the datasheets based on the id (i.e. same protein sequence) and have some prefixes to distinguish between apo and holo affinity parameters. Here we can also define new columns based on previous ones.
Analysis tools always output a datasheet called analysis.
"""
analysis = pipeline.add(MergeDatasheets(datasheets=[boltz_holo_open.datasheets.affinity,
                                                    open_chlorine_aspartate_distance.datasheets.analysis]))
"""
Now we want to make sure we only take poses in which the chlorine atom is close to the aspartate.
By defining the pool, we make sure we will find in the output folder of the filter tool all the structures with open form respecting the expression. If pool is not defined, we will only have a datasheet.
"""
filtered  = pipeline.add(Filter(pool=boltz_holo_open,
                                data=analysis,
                                expression="open_chlorine_distance < 5.0"))

#Prints
pipeline.save()
pipeline.slurm() 
