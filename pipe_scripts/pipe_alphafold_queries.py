#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Build the ColabFold queries CSV from an AlphaFold combinatorics config.

AlphaFold has a single ``proteins`` axis. Each iterated element of that axis
becomes one ColabFold query row (``id,sequence``). For a Bundle, the bundled
sequences are collapsed into ONE query whose ``sequence`` field is the chains
colon-joined (``SEQ_A:SEQ_B``) — ColabFold's native multimer representation,
which auto-selects the multimer model and runs its default paired+unpaired MSA
pipeline.

Axis semantics (mirrors combinatorics.py):
- ``Each(a, b)``      → one monomer query per element: ids ``a``, ``b``.
- ``Bundle(a, b)``    → one complex query: id ``a+b``, sequence ``SEQ_a:SEQ_b``.
- ``Bundle(x, Each(y))`` → one complex per iterated ``y``: id ``x+y_i`` (or
  ``y_i+x`` depending on source order), with ``x`` held static in every row.

Output row ids are taken from the config's pre-computed ``predicted_ids`` so
they agree exactly with what ``get_output_files()`` declared at config time.
The colon-join order within a row follows the source order recorded in the
axis (static-vs-iterated placement honoured the same way Boltz2 does).
"""

import argparse
import os
import sys

import pandas as pd

# Resolve combinatorics IDs identically to configuration time.
_biopipelines_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'biopipelines')
sys.path.insert(0, _biopipelines_dir)
from combinatorics import CombinatoricsConfig, predict_single_output_id  # noqa: E402


def _load_id_to_seq(path):
    """Load an {id: sequence} map from a sequences CSV (columns: id, sequence)."""
    df = pd.read_csv(path)
    if 'id' not in df.columns or 'sequence' not in df.columns:
        raise ValueError(f"sequences CSV {path} must have 'id' and 'sequence' columns")
    return {str(r['id']): str(r['sequence']).strip() for _, r in df.iterrows()}


def _ordered_sources(axis):
    """Return axis sources sorted by their recorded 'order' (stable for ties)."""
    sources = axis.sources or []
    return sorted(
        (s if isinstance(s, dict) else {"path": s, "iterate": True, "order": 0} for s in sources),
        key=lambda s: s.get("order", 0),
    )


def build_queries(config_path, output_csv):
    cfg = CombinatoricsConfig.load(config_path)
    axis = cfg.get_axis("proteins")
    if axis is None:
        raise ValueError("combinatorics config has no 'proteins' axis")

    sources = _ordered_sources(axis)
    if not sources:
        raise ValueError("'proteins' axis has no sources")

    # Per source: (id_to_seq, ids_in_file_order, iterate_flag), preserving the
    # axis source order so colon-join order is deterministic and matches IDs.
    loaded = []
    for src in sources:
        path = src.get("path")
        if not path or not os.path.exists(path):
            raise ValueError(f"sequences source not found: {path}")
        id_to_seq = _load_id_to_seq(path)
        loaded.append({
            "id_to_seq": id_to_seq,
            "ids": list(id_to_seq.keys()),
            "iterate": bool(src.get("iterate", True)),
        })

    static_ids = [i for src in loaded if not src["iterate"] for i in src["ids"]]
    iterated_sources = [src for src in loaded if src["iterate"]]
    iterated_ids = [i for src in iterated_sources for i in src["ids"]]

    # Lookup any id -> sequence across all sources.
    seq_of = {}
    for src in loaded:
        seq_of.update(src["id_to_seq"])

    rows = []
    if axis.mode == "bundle" and not iterated_sources:
        # Pure bundle: one complex from all (unique, order-preserving) ids.
        chain_ids = []
        for src in loaded:
            for i in src["ids"]:
                if i not in chain_ids:
                    chain_ids.append(i)
        row_id = predict_single_output_id(proteins=("bundle", chain_ids, None, [], False))
        rows.append((row_id, [seq_of[i] for i in chain_ids], chain_ids))
    else:
        # Each (no static) or Bundle(static, Each(...)): one complex per iterated id.
        # static_first matches the recorded source order.
        static_first = False
        if static_ids:
            first_static_order = min(s.get("order", 0) for s in sources if not s.get("iterate", True))
            first_iter_order = min(s.get("order", 0) for s in sources if s.get("iterate", True))
            static_first = first_static_order < first_iter_order
        for iid in iterated_ids:
            row_id = predict_single_output_id(
                proteins=("each", iterated_ids, iterated_ids.index(iid), static_ids, static_first)
            )
            chain_ids = ([*static_ids, iid] if static_first else [iid, *static_ids])
            rows.append((row_id, [seq_of[i] for i in chain_ids], chain_ids))

    out = pd.DataFrame(
        [{"id": rid, "sequence": ":".join(seqs)} for rid, seqs, _ in rows],
        columns=["id", "sequence"],
    )
    os.makedirs(os.path.dirname(output_csv) or ".", exist_ok=True)
    out.to_csv(output_csv, index=False)

    for rid, seqs, chain_ids in rows:
        kind = "complex" if len(seqs) > 1 else "monomer"
        print(f"  query {rid}: {kind} ({len(seqs)} chain(s): {'+'.join(chain_ids)})")
    print(f"Wrote {len(rows)} ColabFold query row(s) to {output_csv}")


def main():
    parser = argparse.ArgumentParser(description="Build ColabFold queries CSV from a combinatorics config")
    parser.add_argument("--combinatorics-config", required=True, help="Path to combinatorics config JSON")
    parser.add_argument("--output-csv", required=True, help="Path to write the queries CSV (id,sequence)")
    args = parser.parse_args()
    build_queries(args.combinatorics_config, args.output_csv)


if __name__ == "__main__":
    main()
