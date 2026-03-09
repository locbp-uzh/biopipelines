# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

# tested: 

from biopipelines.pipeline import *
from biopipelines.compound_library import CompoundLibrary
from biopipelines.boltz2 import Boltz2
from biopipelines.panda import Panda
from biopipelines.plot import Plot

with Pipeline(project="TrpRepressor", job="CompoundLibraryScreen"):
    Resources(gpu="A100", time="8:00:00", memory="32GB")
    TrpR = Sequence("MAQQSPYSAAMAEERHQEWLRFVDLLKNAYQNDLHLPLLNLMLTPDEREALGTRVRIVEELLRGEMSQRELKNELGAGIATITRGSNSLKAAPVELRQWLEEVLLKSD",
                    ids="TrpR")
    DNA = Sequence("TGTACTAGTTAACTAGTAC",
                   ids="dna")
    library = CompoundLibrary("./ExamplePipelines/compound_library_TrpR.cdxml")
    cofolded = Boltz2(proteins=Bundle(TrpR,TrpR),
                      dsDNA=DNA,
                      ligands=Each(library))
    merged = Panda(tables=[library.tables.compounds, cofolded.tables.affinity],
                    operations=[Panda.merge(),
                                Panda.calculate({"aff_uM":"10**affinity_pred_value"})])
    Plot(Plot.Scatter(data=merged.tables.result,
                        x="R1",
                        y="aff_uM", # affinity_probability_binary
                        title="Predicted Affinity by R1 group",
                        xlabel="R1",
                        ylabel="Predicted Affinity [uM]" # Predicted Binding Probability [0-1]
                     ),
         Plot.Scatter(data=merged.tables.result,
                        x="R1",
                        y="aff_uM", # affinity_probability_binary
                        title="Predicted Affinity by R1 group",
                        xlabel="R1",
                        ylabel="Predicted Affinity [uM]" # Predicted Binding Probability [0-1]
                     ))


    


    

