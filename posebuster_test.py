from biopipelines.pipeline import *
from biopipelines.load import Load
from biopipelines.posebusters import PoseBusters

with Pipeline("Test","PoseBusters"):
    boltz_predictions = Load("/shares/locbp.chem.uzh/gquarg/BioPipelines/Imatinib/PoseSensitivityGnina_002/005_Boltz2")
    posebusters = PoseBusters(structures=boltz_predictions,
                              ligand="LIG")