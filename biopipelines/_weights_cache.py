# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Shared bash-emitting helpers for populating and reusing weight caches.

Several tools download large model weights (AlphaFold2 params, ESM embeddings,
…) into a cache. When the same weights are needed by more than one tool — or
already exist on the host — re-downloading them per tool wastes bandwidth and
disk. These helpers let an install script (a) populate a *shared* cache only if
it isn't already there (order-independent: whichever tool installs first fills
it), and (b) point a tool's hardcoded cache path at the shared copy without
duplicating the bytes.

Both functions return bash text for use inside a tool's ``_install_script``.
They are deliberately weight-agnostic — the caller supplies the directory, the
sentinel that proves the weights are present, and the download commands — so
the same logic serves AlphaFold params, ESM weights, or anything similar.
"""

from typing import Optional


def ensure_weights_block(dest_dir: str, sentinel: str, download_cmds: str,
                         label: str = "weights") -> str:
    """Bash: ensure ``dest_dir`` exists and is populated, downloading only if absent.

    Args:
        dest_dir: the shared cache directory to create and populate.
        sentinel: a path (typically under ``dest_dir``) whose existence proves
            the weights are already present. The download is skipped iff this
            file/dir exists, so it must be something the download produces last
            or a representative final artefact (e.g. ``params_model_1_ptm.npz``).
        download_cmds: bash run inside the ``if`` body to fetch the weights into
            ``dest_dir`` when the sentinel is missing.
        label: human-readable name for log lines.

    The ``mkdir -p`` runs unconditionally so the folder always exists for any
    other tool to populate later, even when this install downloads nothing.
    """
    return f"""# Ensure shared {label} cache exists, then download only if absent.
mkdir -p "{dest_dir}"
if [ ! -e "{sentinel}" ]; then
    echo "{label} not found in shared cache; downloading to {dest_dir}"
{_indent(download_cmds)}
else
    echo "Reusing existing {label} in shared cache: {dest_dir}"
fi"""


def link_weights_block(target_dir: str, link_path: str,
                       marker: Optional[str] = None,
                       label: str = "weights") -> str:
    """Bash: point a tool's hardcoded cache ``link_path`` at the shared ``target_dir``.

    Tries to create ``link_path`` as a symlink to ``target_dir`` (no byte
    duplication). If symlinking fails — some network/overlay filesystems don't
    support it — it falls back to copying the contents. An existing real
    directory at ``link_path`` is left untouched so we never clobber real data.
    Optionally touches ``marker`` afterwards — a sentinel some upstreams check
    to decide the cache is "ready" and skip their own download (e.g. BioEmu's
    inlined ColabFold checks for ``download_finished.txt``).

    Args:
        target_dir: the shared cache the tool should use.
        link_path: the path the tool's code hardcodes (e.g. ``~/.cache/colabfold/params``).
        marker: optional sentinel to touch after linking/copying.
        label: human-readable name for log lines.
    """
    marker_line = f'\ntouch "{marker}"' if marker else ""
    return f"""# Point {label} cache ({link_path}) at the shared copy, symlink if possible,
# else copy (some filesystems don't support symlinks).
mkdir -p "$(dirname "{link_path}")"
if [ -e "{link_path}" ] && [ ! -L "{link_path}" ]; then
    echo "WARNING: {link_path} is a real path, not relinking {label} cache"
elif ln -sfn "{target_dir}" "{link_path}" 2>/dev/null && [ -e "{link_path}/." ]; then
    echo "Linked {label} cache {link_path} -> {target_dir}"
else
    echo "Symlink unsupported; copying {label} into {link_path}"
    rm -f "{link_path}"
    mkdir -p "{link_path}"
    cp -a "{target_dir}/." "{link_path}/"
fi{marker_line}"""


def _indent(text: str, spaces: int = 4) -> str:
    pad = " " * spaces
    return "\n".join(pad + line if line else line for line in text.splitlines())
