"""Filesystem-layout constants shared across the framework.

Kept in a tiny standalone module so both ``pipeline.py`` and the test suite
can import them without pulling in the heavier pipeline machinery.
"""

# Framework-owned subfolder that holds auto-generated internal tools (e.g. a
# Ligand created from the ``ligand="LIG"`` shorthand). Hidden from a plain
# ``ls`` by the leading dot, and excluded from the public tool numbering.
INTERNAL_FOLDER = ".internal"
