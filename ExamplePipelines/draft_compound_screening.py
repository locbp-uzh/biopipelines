# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

# tested: 13.2.2026

from biopipelines.pipeline import *
from biopipelines.compound_library import CompoundLibrary
from biopipelines.boltz2 import Boltz2
from biopipelines.panda import Panda
from biopipelines.plot import Plot

with Pipeline(project="MarR", job="FragmentScreen"):
    Resources(gpu="A100", time="8:00:00", memory="32GB")
    MarR = Sequence("LFNEIIPLGRLIHMVNQKKDRLLNEYLSPLDITAAQFKVLCSIRCAACITPVELKKVLSVDLGALTRMLDRLVCKGWVERLPNPNDKRGVLVKLTTGGAAICEQCHQLVGQDLHQELTKNLTADEVATLEYLLKKVLP")
    library = CompoundLibrary(library={"candidate": "<aryl><carboxylate>",
                                       "aryl": ["<o-hydroxyphenyl>", 
                                                "<m-hydroxyphenyl>",
                                                "<p-hydroxyphenyl>"],
                                       "o-hydroxyphenyl": r"c1ccc(O)c(c1)",
                                       "m-hydroxyphenyl": r"c1cc(O)cc(c1)",
                                       "p-hydroxyphenyl": r"c1c(O)ccc(c1)",
                                       "carboxylate": r"C(=O)[O-]"},
                              primary_key="candidate")
    boltz = Boltz2(proteins=Bundle(MarR, MarR),
                   ligands=Each(library))
    merged = Panda(tables=[boltz.tables.affinity, 
                           library.tables.compounds],
                   operations=[Panda.merge(on="id")])
    Plot(Plot.Scatter(data=merged.tables.result,
                      x="aryl", 
                      y="affinity_pred_value",
                      title="Predicted Affinity by Aryl Substituent",
                      xlabel="Aryl Group", 
                      ylabel="Predicted Affinity", 
                      grid=True))
    

