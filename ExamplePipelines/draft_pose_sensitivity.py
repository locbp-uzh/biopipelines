# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

# debug 13.02: pose change has to align based on proteins only, FAILED

from biopipelines.pipeline import *
from biopipelines.mutagenesis import Mutagenesis
from biopipelines.boltz2 import Boltz2
from biopipelines.pose_change import PoseChange
from biopipelines.panda import Panda
from biopipelines.plot import Plot

with Pipeline(project="Imatinib", job="PoseSensitivity"):
    Resources(gpu="A100", time="8:00:00", memory="32GB")
    abl1 = PDB("3QRK")
    imatinib = Ligand("STI")
    original = Boltz2(proteins=abl1,
                      ligands=imatinib)
    mutants_sequences = Mutagenesis(original=abl1,
                                    position=93, #AUTH 315
                                    mode="saturation")
    mutants = Boltz2(proteins=mutants_sequences,
                    ligands=imatinib)
    pose = PoseChange(reference_structure=original,
                      sample_structures=mutants,
                      reference_ligand="LIG")
    analysis = Panda(tables=[pose.tables.changes, 
                             mutants.tables.affinity,
                             mutants_sequences.tables.sequences],
                     operations=[Panda.merge()])
    Plot(Plot.Bar(data=analysis.tables.result,
                  x="mutation", 
                  y="ligand_rmsd",
                  title="Ligand RMSD per Mutant",
                  xlabel="Mutation", 
                  ylabel="Ligand RMSD (A)",
                  x_tick_rotation=45, 
                  grid=True))


