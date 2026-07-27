# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Resolve a single value from a TABLE_REFERENCE table column.

Usage:
    python resolve_table_column.py "TABLE_REFERENCE:/path/table.csv:column" "item_id" [map_table.csv]

The optional map_table is the consumed stream's map_table; pass it so an id an
upstream tool renamed still resolves against the table's original id space.

Prints the resolved value to stdout.
"""

import sys
from biopipelines.biopipelines_io import load_table, lookup_table_value


def main():
    if len(sys.argv) not in (3, 4):
        print("Usage: python resolve_table_column.py <TABLE_REFERENCE> <item_id> [map_table]",
              file=sys.stderr)
        sys.exit(1)

    reference = sys.argv[1]
    item_id = sys.argv[2]
    map_table_paths = [sys.argv[3]] if len(sys.argv) == 4 else None

    table, column = load_table(reference)
    value = lookup_table_value(table, item_id, column, map_table_paths=map_table_paths)
    print(value)


if __name__ == "__main__":
    main()
