"""PDB terminator records must be written at full width, not bare.

`END` and `TER` are 3-character records, and the PDB spec allows that. But
third-party parsers routinely index fixed columns before checking the record
type — Aggrescan3D's reads `line[21]` (the chain column) on every line, so a
bare "END\n" raises `IndexError: string index out of range` and the whole run
dies on the first structure.

Every structure this framework writes is an input to some other tool, so the
terminators are padded to 80 columns on write. Cheap here, and it removes a
failure mode that is invisible until a downstream tool crashes.
"""

import glob
import io

# Files that write PDB structures downstream tools consume.
WRITERS = [
    "pipe_scripts/pipe_gnina.py",
    "pipe_scripts/pipe_neuralplexer_postprocess.py",
    "pipe_scripts/pipe_openmm.py",
    "pipe_scripts/pipe_posebusters.py",
    "pipe_scripts/pipe_rfdaa_prepare_ligand.py",
]

BARE_END = 'write("END' + chr(92) + 'n")'
BARE_TER = 'write("TER' + chr(92) + 'n")'


def test_no_bare_terminator_writes_in_known_writers():
    offenders = []
    for path in WRITERS:
        src = io.open(path, encoding="utf-8").read()
        if BARE_END in src or BARE_TER in src:
            offenders.append(path)
    assert not offenders, (
        "bare END/TER writes found in %s — pad to 80 columns, a parser that "
        "indexes line[21] before checking the record type crashes on them"
        % offenders)


def test_no_bare_terminator_writes_anywhere():
    """New writers must not reintroduce it."""
    offenders = []
    for path in glob.glob("pipe_scripts/*.py") + glob.glob("biopipelines/*.py"):
        src = io.open(path, encoding="utf-8", errors="replace").read()
        if BARE_END in src or BARE_TER in src:
            offenders.append(path)
    assert not offenders, "bare END/TER writes found in %s" % offenders


def test_writers_pad_terminators():
    """The known writers actually emit the padded form."""
    for path in WRITERS:
        src = io.open(path, encoding="utf-8").read()
        assert 'ljust(80)' in src, "%s no longer pads its terminator records" % path
