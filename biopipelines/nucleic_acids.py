# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""How the framework's polymer axes map onto a co-folding model's chains.

Every co-folding wrapper faces the same two questions: what does the axis name `dsDNA` mean to
the model, and what does "double-stranded" expand to. The answer is the same for all of them —
a double-stranded axis is two chains, the second being the reverse complement — so it lives
here rather than once per config generator. Boltz2 and OpenFold3 both read it.

Kept dependency-free so a `pipe_*` script can import it out of the biopipelines env.
"""

# The axis names a user writes, mapped to the entity type a config generator switches on.
# `sequences` and `compounds` are the stream names, accepted so a config written before the
# axis names were recorded still resolves.
AXIS_NAME_TO_ENTITY_TYPE = {
    "proteins": "protein",
    "sequences": "protein",
    "ssDNA": "ssdna",
    "dsDNA": "dsdna",
    "ssRNA": "ssrna",
    "dsRNA": "dsrna",
    "ligands": "ligand",
    "compounds": "ligand",
}

# These emit two chains; the second is the reverse complement of the first.
DOUBLE_STRANDED_ENTITY_TYPES = {"dsdna", "dsrna"}

RNA_ENTITY_TYPES = ("rna", "ssrna", "dsrna")
DNA_ENTITY_TYPES = ("dna", "ssdna", "dsdna")

_DNA_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")
_RNA_COMPLEMENT = str.maketrans("ACGUacgu", "UGCAugca")


def reverse_complement(sequence: str, entity_type: str) -> str:
    """Return the reverse complement of a DNA or RNA sequence."""
    table = _RNA_COMPLEMENT if entity_type in RNA_ENTITY_TYPES else _DNA_COMPLEMENT
    return sequence.translate(table)[::-1]
