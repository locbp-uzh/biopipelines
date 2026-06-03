# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Minimal LigandMPNN runs exercising every kwarg covered by
tests/tool_parameters/test_ligand_mpnn_params.py. Hemoglobin α-chain (1A3N)
is used as a heme-binding backbone so the `ligand` and `design_within`
flags have a meaningful target.

LigandMPNN's per-position bias flag is --bias_AA_per_residue (JSON, not
JSONL): {"<chain><residue>": {"<AA>": <bias>, ...}, ...}.
"""

import json
from pathlib import Path

from biopipelines.pipeline import *
from biopipelines import LigandMPNN, Ligand

_HERE = Path(__file__).resolve().parent
BIAS_PATH = str(_HERE / "_lmpnn_bias_per_residue.json")
Path(BIAS_PATH).write_text(json.dumps({"A50": {"H": 2.0}}))


with Pipeline(project="ToolParameters",
              job="LigandMPNN",
              description="LigandMPNN parameter coverage — bias_AA_per_residue, temperature, seed, model"):

    Resources(gpu="A100", time="2:00:00", memory="16GB")

    backbone = PDB("1A3N", ids="HBA")
    hem = Ligand(code="HEM")
    atp = Ligand(code="ATP")

    # 1: defaults with HEM ligand
    Suffix("hem_default")
    LigandMPNN(structures=backbone, ligand=hem)

    # 2: alternative ligand + design_within
    Suffix("atp_within")
    LigandMPNN(structures=backbone, ligand=atp, design_within=8.0)

    # 3: batch / batches knobs
    Suffix("batches")
    LigandMPNN(structures=backbone, ligand=hem, num_sequences=4, num_batches=3)

    # 4: model variant
    Suffix("model_v32_020")
    LigandMPNN(structures=backbone, ligand=hem, model="v_32_020")

    # 5: temperature
    Suffix("temperature")
    LigandMPNN(structures=backbone, ligand=hem, temperature=0.3, num_sequences=2)

    # 6: per-residue bias (LigandMPNN-native --bias_AA_per_residue)
    Suffix("bias")
    LigandMPNN(structures=backbone, ligand=hem, bias_AA_per_residue=BIAS_PATH, num_sequences=2)

    # 7: deterministic seed
    Suffix("seed")
    LigandMPNN(structures=backbone, ligand=hem, seed=4242, num_sequences=2)
