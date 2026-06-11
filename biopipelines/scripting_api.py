# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Public surface a user ``Scripting`` script imports.

A Scripting script implements two functions and uses only the declarations
here, never the framework internals:

    from biopipelines.scripting_api import Stream, Table

    def configuration(inputs):
        ids = inputs["structures"].ids
        return {"structures": Stream("pdb", ids),
                "scores": Table(columns=["id", "rmsd"])}

    def execution(inputs, outputs):
        for sid, path in inputs["structures"].iterate():
            out = outputs.file("structures", f"{sid}.pdb")
            ...
            outputs.row("scores", {"id": sid, "rmsd": ...})
        return {"structures": [{"id": sid, "file": out} for ...]}

``configuration`` runs in the pipeline's own Python at config time; it must be
pure shape-prediction (no disk reads). ``execution`` runs at execution time in
the configured env, against real files.
"""

from typing import List, Optional


class Stream:
    """Declared output stream from ``configuration``.

    ``format`` is the stream format (``pdb``, ``cif``, ``sdf``, ``csv``, ...).
    A value-based stream is just ``format="csv"`` whose rows carry no ``file``
    column; a file-based stream carries per-id files. ``ids`` are the predicted
    ids, derived from the input streams' ids (lazy/pattern ids pass through).
    """

    def __init__(self, format: str, ids: List[str]):
        self.format = format
        self.ids = list(ids)


class Table:
    """Declared standalone output table from ``configuration``."""

    def __init__(self, columns: List[str], description: str = ""):
        self.columns = list(columns)
        self.description = description
