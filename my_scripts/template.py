# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Template for a Scripting tool script.

A Scripting script is an escape hatch for one-off glue logic that does not
justify a full tool wrapper. It implements two functions:

  configuration(inputs)            — runs in the pipeline's own Python at config
                                     time. PURE shape prediction: derive output
                                     ids from input ids, declare output Streams /
                                     Tables. No disk reads, no real files exist.

  execution(inputs, outputs)       — runs at execution time in the configured
                                     env, against real files. Read each input,
                                     write per-id files via outputs.file(...),
                                     accumulate table rows via outputs.row(...),
                                     and return the stream map rows.

Use it from a pipeline (Scripting searches the configured scripts folder, so a
bare filename works):

    Scripting("template.py",
              inputs={"structures": some_tool},   # StandardizedOutput / DataStream
              env="biopipelines")                 # optional; defaults to biopipelines

The keys of `inputs` are the names this script reads back from `inputs[...]`.
"""

from biopipelines.scripting_api import Stream, Table


def configuration(inputs):
    """Declare the output shape. Runs at config time — prediction only.

    `inputs[name]` exposes `.ids` (the declared input ids). Derive output ids
    from them; never touch the filesystem here.
    """
    ids = inputs["structures"].ids

    # Return a dict mapping each output name to a Stream or a Table.
    #   Stream(format, ids) — per-id files (pdb/cif/sdf/...) or, with
    #                         format="csv", a value-based stream (no files).
    #   Table(columns)      — a standalone results CSV (no id tracking).
    return {
        "structures": Stream("pdb", ids),
        "scores": Table(columns=["id", "n_atoms"]),
    }


def execution(inputs, outputs):
    """Do the work. Runs at execution time against real files.

    `inputs["structures"].iterate()` yields (id, file_path) for a file-based
    stream, or (id, value_dict) for a value-based one. `outputs.file(stream,
    name)` returns a path under that stream's folder to write into;
    `outputs.row(table, mapping)` accumulates one results row.

    Return a dict of {stream_name: [row, ...]} where each row has at least
    `id` and (for file streams) `file`. The runner writes the map_tables.
    """
    structure_rows = []

    for sid, path in inputs["structures"].iterate():
        # Example payload: copy the input through and count its atoms. Replace
        # this block with the actual transformation.
        out_path = outputs.file("structures", f"{sid}.pdb")
        n_atoms = 0
        with open(path) as src, open(out_path, "w") as dst:
            for line in src:
                if line.startswith(("ATOM", "HETATM")):
                    n_atoms += 1
                dst.write(line)

        structure_rows.append({"id": sid, "file": out_path})
        outputs.row("scores", {"id": sid, "n_atoms": n_atoms})

    return {"structures": structure_rows}
