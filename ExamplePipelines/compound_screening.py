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
    TrpR = Sequence("MAQQSPYSAAMAEERHQEWLRFVDLLKNAYQNDLHLPLLNLMLTPDEREALGTRVRIVEELLRGEMSQRELKNELGAGIATITRGSNSLKAAPVELRQWLEEVLLKSD")
    DNA = Sequence("TGTACTAGTTAACTAGTAC")
    library = CompoundLibrary("/path/to/library.cdxml")
    cofolded = Boltz2(proteins=Bundle(TrpR,TrpR),
                      dna=DNA,
                      ligands=Each(library))
    merged = Panda(tables=cofolded.tables.affinity,
                   operations=[Panda.calculate({"aff_uM":"10**affinity_pred_value"})])


    


    

