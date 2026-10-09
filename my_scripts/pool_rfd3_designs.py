# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Pool a packed RFdiffusion3 + ProteinMPNN run into one sequence table and one design table.

A packed run holds one RFdiffusion3 step per opening, and `prefix=` does not reach RFdiffusion3's output ids: every step emits `<pdb>_1 ... <pdb>_N`, so pooling them naively collapses N openings onto N ids. Ids are rewritten here to `s<rfd3 step>_<mpnn id>`, which is unique because the step number is.

The design table carries each sequence's contig and `fixed` selection, read from the RFdiffusion3 step that produced its backbone, so downstream scoring needs no cross-run id matching -- which matters when the folding half runs on a different machine.

Usage: python pool_rfd3_designs.py <run_dir> <sequences_csv> <designs_csv>
"""

import csv
import re
import sys
from pathlib import Path


def pool(run_dir):
    run = Path(run_dir)
    steps = sorted(p.name for p in run.iterdir() if re.fullmatch(r"\d+_ProteinMPNN", p.name))
    sequences, designs = [], []

    for mpnn in steps:
        number = int(mpnn.split("_")[0])
        # The RFdiffusion3 step is the one immediately before its ProteinMPNN.
        rfd3 = next(p.name for p in run.iterdir()
                    if p.name.startswith(f"{number - 1:03d}_RFdiffusion3"))

        with open(run / rfd3 / "tables" / "structures.csv", encoding="utf-8") as handle:
            backbones = {row["id"]: row for row in csv.DictReader(handle)}

        with open(run / mpnn / "sequences" / "sequences.csv", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                backbone = backbones[row["structures.id"]]
                pooled = f"s{number - 1:03d}_{row['id']}"
                sequences.append({"id": pooled, "sequence": row["sequence"]})
                designs.append({"id": pooled,
                                "step": f"{number - 1:03d}",
                                "backbone": row["structures.id"],
                                "contigs": backbone["contig"],
                                "fixed": backbone["fixed"]})

    return sequences, designs


def write(rows, path, fields):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    run_dir, sequences_csv, designs_csv = sys.argv[1:4]
    sequences, designs = pool(run_dir)
    write(sequences, sequences_csv, ["id", "sequence"])
    write(designs, designs_csv, ["id", "step", "backbone", "contigs", "fixed"])
    print(f"{len(sequences)} sequences over {len({d['step'] for d in designs})} openings")
