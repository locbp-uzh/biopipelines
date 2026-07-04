# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Public surface a user ``Scripting`` script imports.

A Scripting script implements two functions:

    def configuration(inputs):
        from biopipelines.scripting_api import Stream, Table
        ids = inputs["structures"].ids
        return {"structures": Stream("pdb", ids),
                "scores": Table(columns=["id", "rmsd"])}

    def execution(inputs, outputs):
        for sid, path in inputs["structures"].iterate():
            out = outputs["structures"].file(sid, f"{sid}.pdb")
            ...
            outputs["scores"].row({"id": sid, "rmsd": ...})

``configuration`` runs in the pipeline's own Python at config time; it must be pure shape-prediction (no disk reads). ``execution`` runs at execution time in the configured env, against real files, and returns nothing.

Imports: keep the ``scripting_api`` import inside ``configuration`` (only it needs Stream/Table) and put tool deps (numpy, rdkit, ...) inside ``execution``. A script may import biopipelines from ``execution`` (e.g. ``get_mapped_ids``, ``DataStream``) when its env carries biopipelines' deps — the default biopipelines env does, and a custom env does after ``pip install -e ".[scripting]"``.

Outputs are handles reached by ``outputs[name]``, where ``name`` is a key the ``configuration`` dict declared:

    outputs[name].file(id, filename)               # file stream: alloc a path under the stream folder + record the {id, file} row
    outputs[name].add(id=, file=, value=, **cols)  # explicit row (any columns)
    outputs[name].row({"id": ..., ...})            # append one row mapping
    outputs[name].dataframe(df)                     # hand over a whole DataFrame
    outputs[name] = "path/to.csv"                   # adopt an existing CSV

``file()`` is file-stream-only; ``row``/``add``/``dataframe`` fill value-based streams and standalone tables alike (both are id-keyed CSVs).

When ``execution`` intentionally drops an input id (a data-dependent filter), report it with ``outputs.drop(id, cause="...")`` so the completion check excuses that id's absent outputs instead of flagging the step as failed.
"""

from typing import List, Optional


class Stream:
    """Declared output stream from ``configuration``.

    ``format`` is the stream format (``pdb``, ``cif``, ``sdf``, ``csv``, ...). A value-based stream is just ``format="csv"`` whose rows carry no ``file`` column; a file-based stream carries per-id files. ``ids`` are the predicted ids, derived from the input streams' ids (lazy/pattern ids pass through).
    """

    def __init__(self, format: str, ids: List[str]):
        self.format = format
        self.ids = list(ids)


class Table:
    """Declared standalone output table from ``configuration``."""

    def __init__(self, columns: List[str], description: str = ""):
        self.columns = list(columns)
        self.description = description
