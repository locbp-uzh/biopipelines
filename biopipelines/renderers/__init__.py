"""Output renderers, loaded by path rather than imported.

`pipeline.py` and `outputs.py` resolve these against this package's directory,
which is the same spelling from a clone and from site-packages — so the paths
in `config.*.yaml` stay `renderers/<name>.py` either way.
"""
