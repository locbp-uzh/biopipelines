# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Format-aware slicers for shared-file DataStreams.

When a tool (typically Panda or Pool) filters a stream whose ``files`` is a
single shared artifact (e.g. a multi-record FASTA), the artifact must be
re-emitted containing only the kept IDs — and optionally with those IDs
renamed so the file agrees with the downstream stream's renamed ids. The
slicing rule depends on the file format — text vs binary, line-oriented
vs record-oriented — so each format owns a small function that knows how
to walk its own structure.

Contract: ``slicer(src_path, dest_path, kept_ids, rename_map=None) -> None``
    - Read ``src_path``, keep only records whose ID is in ``kept_ids``
      (preserving original order), write the subset to ``dest_path``.
    - If ``rename_map`` is provided, rewrite each kept record's ID from
      ``orig`` to ``rename_map[orig]`` before writing. Records with no
      rename entry keep their original ID.
    - Create the parent directory of ``dest_path`` as needed.
    - Idempotent: re-running with the same arguments yields the same output.
    - IDs that appear in ``kept_ids`` but not in the source file are silently
      ignored — the caller (e.g. Panda) decides whether that's an error via
      its own ignore_missing semantics.

A slicer registered without a ``rename_map`` parameter is callable too:
``get_slicer`` wraps it so a ``rename_map=None`` call is dropped at the
boundary. To register a format that should support renaming, accept
``rename_map`` in the signature.

Registering a new format:

    from biopipelines.stream_slicers import register

    @register("sdf")
    def _slice_sdf(src, dest, kept_ids, rename_map=None):
        ...
"""

import inspect
import os
import re
from typing import Callable, Dict, Iterable, List, Optional, Set, Tuple


SlicerFn = Callable[..., None]
MergerFn = Callable[[List[str], str], None]

_REGISTRY: Dict[str, SlicerFn] = {}
_MERGERS: Dict[str, MergerFn] = {}


def register(fmt: str):
    """Decorator that registers a slicer under a format key.

    Format keys are case-insensitive. A function may be registered under
    multiple keys via stacked decorators (e.g. @register("fasta") and
    @register("fa") for the same parser).
    """
    def deco(fn: SlicerFn) -> SlicerFn:
        _REGISTRY[fmt.lower()] = fn
        return fn
    return deco


def register_merger(fmt: str):
    """Decorator that registers a merger for a format key.

    Mergers combine several already-sliced files into a single output.
    Used by Panda's multi-pool path: slice each pool independently, then
    merge the results so the downstream stream has one shared artifact
    instead of N pool-specific ones.
    """
    def deco(fn: MergerFn) -> MergerFn:
        _MERGERS[fmt.lower()] = fn
        return fn
    return deco


def _accepts_rename_map(fn: SlicerFn) -> bool:
    """True if ``fn`` declares a ``rename_map`` parameter or **kwargs."""
    try:
        sig = inspect.signature(fn)
    except (TypeError, ValueError):
        return False
    for p in sig.parameters.values():
        if p.kind is inspect.Parameter.VAR_KEYWORD:
            return True
        if p.name == "rename_map":
            return True
    return False


def get_slicer(fmt: str) -> SlicerFn:
    """Return the slicer registered for ``fmt``.

    Raises ValueError if no slicer is registered — no silent fallback to
    "copy whole file", because that would let stale records leak through
    a filter step. Register a slicer in this module instead.

    Slicers whose signature omits ``rename_map`` get a thin wrapper that
    drops the keyword at the boundary, so callers always invoke with the
    full contract regardless of how the underlying function was declared.
    """
    key = (fmt or "").lower()
    if key not in _REGISTRY:
        raise ValueError(
            f"No shared-file slicer registered for format '{fmt}'. "
            f"Available: {sorted(_REGISTRY)}. "
            f"Register one in biopipelines/stream_slicers.py."
        )
    fn = _REGISTRY[key]
    if _accepts_rename_map(fn):
        return fn

    def _wrapped(src, dest, kept_ids, rename_map=None, **kwargs):
        if rename_map:
            raise ValueError(
                f"Slicer for format '{fmt}' does not support rename_map; "
                f"update its signature to accept rename_map=None."
            )
        return fn(src, dest, kept_ids, **kwargs)

    return _wrapped


def get_merger(fmt: str) -> MergerFn:
    """Return the merger registered for ``fmt`` (combines sliced files)."""
    key = (fmt or "").lower()
    if key not in _MERGERS:
        raise ValueError(
            f"No shared-file merger registered for format '{fmt}'. "
            f"Available: {sorted(_MERGERS)}. "
            f"Register one in biopipelines/stream_slicers.py."
        )
    return _MERGERS[key]


def available_formats() -> List[str]:
    """Sorted list of format keys with a registered slicer."""
    return sorted(_REGISTRY)


# Header tokenizer for FASTA records. The id is the longest prefix of the
# header (after the leading '>') that contains no whitespace, comma, or
# semicolon — `,` and `;` are common annotation-field separators in raw
# upstream dumps (e.g. ProteinMPNN: `>4LCD, score=...`). The pipe `|` is
# left as part of the id because NCBI-style compound ids (`sp|P12345|FOO`)
# embed it as a structural separator inside the id itself.
_FASTA_ID_TERMINATORS = re.compile(r"[\s,;]")


def _parse_fasta_header(line: str) -> Tuple[str, str]:
    """Return (id, suffix) from a FASTA header line.

    ``id`` is the canonical record id; ``suffix`` is the remainder of the
    line including the separator character that terminated the id (so the
    rest of the line can be re-emitted verbatim around a renamed id).
    """
    body = line[1:]  # strip leading '>'
    m = _FASTA_ID_TERMINATORS.search(body)
    if m is None:
        return body.rstrip("\n"), ""
    return body[:m.start()], body[m.start():]


@register("fasta")
@register("fa")
def _slice_fasta(src: str, dest: str, kept_ids: Iterable[str],
                 rename_map: Optional[Dict[str, str]] = None) -> None:
    """Slice a FASTA file by record ID, optionally renaming kept headers.

    Each FASTA record is ``>ID[,;\\s]...\\n<sequence lines>``. The record
    ID is the longest prefix of the header (after ``>``) that contains
    no whitespace, comma, or semicolon — pipes are preserved (NCBI-style
    compound ids stay intact). Records whose ID is in ``kept_ids`` are
    written to ``dest`` in their original order; all others are dropped.
    If ``rename_map`` is provided, kept headers are rewritten so the FASTA
    agrees with the downstream stream's ids; everything after the id
    terminator is preserved verbatim (descriptions, score fields, etc.).
    """
    keep: Set[str] = {str(i) for i in kept_ids}
    rmap = rename_map or {}
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)

    cur_id = None
    cur_lines: List[str] = []

    with open(src, "r") as f_in, open(dest, "w") as f_out:
        def flush():
            if cur_id is not None and cur_id in keep:
                new_id = rmap.get(cur_id, cur_id)
                if new_id != cur_id and cur_lines:
                    header = cur_lines[0]
                    _, suffix = _parse_fasta_header(header.rstrip("\n"))
                    nl = "\n" if header.endswith("\n") else ""
                    cur_lines[0] = f">{new_id}{suffix}{nl}"
                f_out.writelines(cur_lines)

        for line in f_in:
            if line.startswith(">"):
                flush()
                cur_id, _ = _parse_fasta_header(line.rstrip("\n"))
                cur_lines = [line]
            else:
                cur_lines.append(line)
        flush()


@register("csv")
def _slice_csv(src: str, dest: str, kept_ids: Iterable[str],
               rename_map: Optional[Dict[str, str]] = None) -> None:
    """Slice a CSV by ``id`` column, optionally renaming kept rows.

    Reads the CSV with pandas, keeps only rows whose ``id`` value is in
    ``kept_ids``, and writes the subset to ``dest``. Requires an ``id``
    column — raises ValueError otherwise. If ``rename_map`` is provided,
    the kept rows' ``id`` values are rewritten before write so the CSV
    agrees with the downstream stream's ids.
    """
    import pandas as pd

    df = pd.read_csv(src)
    if "id" not in df.columns:
        raise ValueError(f"CSV slicer requires an 'id' column in {src}")
    keep = {str(i) for i in kept_ids}
    out = df[df["id"].astype(str).isin(keep)].copy()
    if rename_map:
        out["id"] = out["id"].astype(str).map(lambda v: rename_map.get(v, v))
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
    out.to_csv(dest, index=False)


# ── mergers ───────────────────────────────────────────────────────────────────

@register_merger("fasta")
@register_merger("fa")
def _merge_fasta(parts: List[str], dest: str) -> None:
    """Concatenate FASTA parts byte-for-byte, ensuring record separation."""
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
    with open(dest, "w") as f_out:
        for i, part in enumerate(parts):
            if not part or not os.path.exists(part):
                continue
            with open(part, "r") as f_in:
                text = f_in.read()
            if not text:
                continue
            # Make sure consecutive parts don't fuse a sequence line into the
            # next header — add a newline if the previous chunk lacked one.
            if i > 0 and not text.startswith(">"):
                # Defensive: if a part starts mid-sequence, prepend a newline.
                f_out.write("\n")
            f_out.write(text)
            if not text.endswith("\n"):
                f_out.write("\n")


@register_merger("csv")
def _merge_csv(parts: List[str], dest: str) -> None:
    """Vertical-concat CSV parts; column union, NaN-fill missing columns."""
    import pandas as pd
    frames = []
    for part in parts:
        if part and os.path.exists(part):
            frames.append(pd.read_csv(part))
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
    if not frames:
        # Touch an empty file so downstream existence checks pass.
        open(dest, "w").close()
        return
    merged = pd.concat(frames, ignore_index=True, sort=False)
    merged.to_csv(dest, index=False)
