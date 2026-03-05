# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

# tested: 

from biopipelines.pipeline import *
from biopipelines.compound_library import CompoundLibrary
from biopipelines.boltz2 import Boltz2
from biopipelines.panda import Panda
from biopipelines.plot import Plot

with Pipeline(project="CarbonicAnhydrase", job="FragmentScreen"):
    Resources(gpu="A100", time="8:00:00", memory="32GB")
    CaII = Sequence("3KS3", ids="CaII")
    library = CompoundLibrary(library={"candidate": "<sulfonamide><ring><acetamide>",
                                       "sulfonamide": r"NS(=O)(=O)",
                                       "ring": ["<thiadiazole>", 
                                                "<furan>",
                                                "<cyclopentane>"],
                                       "thiadiazole": r"C1=NN=C(S1)",
                                       "furan": r"C1=CC=C(O1)",
                                       "cyclopentane": r"C1CCC(C1)",
                                       "acetamide": r"NC(=O)C"},
                              primary_key="candidate")
    zn = Ligand("ZN")
    cofolded = Boltz2(proteins=CaII,
                      ligands=Bundle(Each(library),zn))
    merged = Panda(tables=[cofolded.tables.affinity, 
                           library.tables.compounds],
                   operations=[Panda.merge(on="id"),
                               Panda.calculate({"aff_uM":"10**affinity_pred_value"})])
    Plot(Plot.Scatter(data=merged.tables.result,
                      x="aryl", 
                      y="aff_uM",
                      title="Predicted Affinity by Aryl Substituent",
                      xlabel="Aryl Group", 
                      ylabel="Predicted Affinity [uM]", 
                      grid=True))

    


    

