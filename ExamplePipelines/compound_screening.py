# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

# tested: 13.2.2026

from biopipelines.pipeline import *
from biopipelines.compound_library import CompoundLibrary
from biopipelines.boltz2 import Boltz2
from biopipelines.panda import Panda
from biopipelines.plot import Plot
#    MarR = Sequence("LFNEIIPLG...")
#    MarR = Sequence("LFNEIIPLGRLIHMVNQKKDRLLNEYLSPLDITAAQFKVLCSIRCAACITPVELKKVLSVDLGALTRMLDRLVCKGWVERLPNPNDKRGVLVKLTTGGAAICEQCHQLVGQDLHQELTKNLTADEVATLEYLLKKVLP")

with Pipeline(project="MarR", job="FragmentScreen"):
    Resources(gpu="A100", time="8:00:00", memory="32GB")
    MarR = Sequence("LFNEIIPLG...")
    library = CompoundLibrary(library={"candidate": "<aryl><carboxylate>",
                                       "aryl": ["<o-hydroxyphenyl>", 
                                                "<m-hydroxyphenyl>",
                                                "<p-hydroxyphenyl>"],
                                       "o-hydroxyphenyl": r"c1ccc(O)c(c1)",
                                       "m-hydroxyphenyl": r"c1cc(O)cc(c1)",
                                       "p-hydroxyphenyl": r"c1c(O)ccc(c1)",
                                       "carboxylate": r"C(=O)[O-]"},
                              primary_key="candidate")
    cofolded = Boltz2(proteins=Bundle(MarR, MarR),
                   ligands=Each(library))
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
    
    atp=Ligand("ATP")
    boltz_ligand_affinity = Boltz2(proteins=Bundle(MarR, MarR),
                                   ligands=Bundle(Each(library),atp))
    boltz_atp_affinity = Boltz2(proteins=Bundle(MarR, MarR),
                                   ligands=Bundle(atp,Each(library)))
    

    


    

