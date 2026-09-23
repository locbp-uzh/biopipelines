#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Reading a combinatorics axis's rows at execution time.

Every multi-axis tool -- Boltz2, OpenFold3 -- faces the same job before it can write a config:
turn an axis's source map_tables into the records that axis contributes, honouring `Each` vs
`Bundle`, the id filter a narrowed stream carries, and `group_by`. That is what lives here.

What deliberately does NOT live here is the walk that pairs the axes into configs. It assigns
chain ids as it goes, and each model interleaves its own entry building and constraints with
that assignment, so a shared walk would have to be parameterised into something less readable
than the two loops it replaced. The loading is identical; the walk is not.

Imported by pipe scripts that run out of the biopipelines env, so it keeps to pandas and the
framework's own id helpers.
"""

import os
import sys
from typing import Dict, List, Optional

import pandas as pd

_biopipelines_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'biopipelines')
sys.path.insert(0, _biopipelines_dir)
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.id_map_utils import get_mapped_ids
import id_patterns


def load_axis_data(axis_config: Dict) -> tuple:
    """
    Load data from an axis's source CSV files.

    Args:
        axis_config: Dict with 'name', 'mode', 'sources' keys
                    sources can be list of strings or list of dicts with 'path', 'iterate', and 'order'

    Returns:
        Tuple of (iterated_data, static_data, static_first):
        - iterated_data: List of dicts from sources with iterate=True
        - static_data: List of dicts from sources with iterate=False
        - static_first: True if static sources should be added before iterated (based on order)
    """
    iterated_data = []
    static_data = []
    sources = axis_config.get('sources', [])

    # Track order to determine if static should come first
    min_iterated_order = float('inf')
    min_static_order = float('inf')

    for source in sources:
        # Handle both old string format and new dict format
        keep_ids = None
        if isinstance(source, dict):
            source_path = source.get('path')
            is_iterate = source.get('iterate', True)
            order = source.get('order', 0)
            keep_ids = source.get('ids')
        else:
            source_path = source
            is_iterate = True
            order = 0

        if not source_path or not os.path.exists(source_path):
            print(f"Warning: Source file not found: {source_path}")
            continue
        try:
            df = pd.read_csv(source_path, dtype={'id': str})
            records = df.to_dict('records')
            if keep_ids:
                by_id = {str(r['id']): r for r in records if 'id' in r}
                selected = id_patterns.select_ids(
                    [str(p) for p in keep_ids], list(by_id.keys())
                )
                records = [by_id[i] for i in selected]
            if is_iterate:
                if isinstance(source, dict) and source.get('group_by'):
                    records = group_records(records, source)
                iterated_data.extend(records)
                min_iterated_order = min(min_iterated_order, order)
            else:
                static_data.extend(records)
                min_static_order = min(min_static_order, order)
        except Exception as e:
            print(f"Error loading {source_path}: {e}")
            sys.exit(1)

    # Static comes first if its minimum order is less than iterated's minimum order
    static_first = min_static_order < min_iterated_order

    return iterated_data, static_data, static_first


def group_records(records: List[Dict], source: Dict) -> List[Dict]:
    """Collapse a source's rows into one iteration element per group key.

    A grouped element carries its members under ``__members__`` and takes the group's own
    id, so downstream id prediction and provenance see one item where the config will hold
    several chains. Membership is the framework's id matching, which is what lets a chain
    row find its design however the ids upstream were renamed.
    """
    group_path = source['group_by']
    if not os.path.exists(group_path):
        print(f"Error: group source not found: {group_path}")
        sys.exit(1)
    group_rows = pd.read_csv(group_path, dtype={'id': str}).to_dict('records')
    group_ids = [str(r['id']) for r in group_rows if 'id' in r]
    patterns = [str(p) for p in (source.get('group_ids') or [])]
    if patterns:
        group_ids = id_patterns.select_ids(patterns, group_ids)

    by_id = {str(r['id']): r for r in records if 'id' in r}
    members = get_mapped_ids(group_ids, list(by_id.keys()), unique=False)

    grouped = []
    for gid in group_ids:
        member_ids = members.get(gid) or []
        if not member_ids:
            print(f"Warning: group '{gid}' has no member rows; skipping")
            continue
        grouped.append({'id': gid, '__members__': [by_id[m] for m in member_ids]})
    print(f"Grouped {len(records)} row(s) into {len(grouped)} group(s)")
    return grouped
