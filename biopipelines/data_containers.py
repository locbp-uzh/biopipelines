# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Data containers: the typed views a tool's declared outputs present to a pipeline.

TableMetadata holds a table's declaration (name, path, columns, description); TableInfo wraps one so that attribute access yields a TableReference to a column; TableContainer, IndexedTableContainer and StreamContainer group them by name, by id and by stream.

Named for the containers rather than for tables alone, because StreamContainer lives here too; not tables.py, because table.py (the Table tool) and table_utils.py already exist.

Layered below outputs.py and base_config.py, both of which re-export these names, so importing them from base_config keeps working.
"""

from typing import Any, Dict, List, Optional

try:
    from .biopipelines_io import TableReference
except ImportError:
    import os
    import sys
    sys.path.append(os.path.dirname(__file__))
    from biopipelines_io import TableReference


class TableMetadata:
    """Holds table metadata (name, path, columns, description)."""

    def __init__(self, name: str, path: str, columns: List[str] = None,
                 description: str = ""):
        self.name = name
        self.path = path
        self.columns = columns or []
        self.description = description

    def __repr__(self) -> str:
        return f"TableMetadata(name='{self.name}', path='{self.path}', columns={self.columns})"


class TableInfo:
    """Information about a table including name, path, and expected columns.

    Metadata is accessed via the .info property:
        table.info.path, table.info.columns, etc.

    Any other public attribute access returns a TableReference:
        table.fixed -> TableReference(path, "fixed")

    Names starting with ``_`` are refused with AttributeError, because they are never column names and answering them made every ``hasattr(table, ...)`` duck-type probe elsewhere in the codebase come back true.
    """

    def __init__(self, name: str, path: str, columns: List[str] = None,
                 description: str = ""):
        self._info = TableMetadata(name, path, columns, description)

        # Set column attributes for IDE autocompletion
        for column in self._info.columns:
            setattr(self, column, self._create_column_reference(column))

    @property
    def info(self) -> TableMetadata:
        """Access table metadata (name, path, columns, description)."""
        return self._info

    def _create_column_reference(self, column_name: str):
        """Create a TableReference for column access."""
        return TableReference(self._info.path, column_name)

    def __getattr__(self, column_name: str):
        """Return a TableReference for any public column name, refusing private ones.

        Declared columns never reach here (``__init__`` sets them as real instance attributes), so this only serves tables whose schema is not declared up front -- and, before the guard, every ``hasattr`` probe, dunder lookup and copy/pickle protocol hook, which is why a ``hasattr(table, '_entries')`` check in ``table_utils`` fired for every ordinary table.
        """
        if column_name.startswith("_"):
            raise AttributeError(
                f"{type(self).__name__!r} object has no attribute {column_name!r} "
                f"(names starting with '_' are not treated as column references)"
            )
        return self._create_column_reference(column_name)

    def __str__(self) -> str:
        if self._info.columns:
            col_display = ', '.join(self._info.columns)
            return f"${self._info.name} ({col_display})"
        return f"${self._info.name}"

    def __repr__(self) -> str:
        return f"TableInfo(name='{self._info.name}', path='{self._info.path}', columns={self._info.columns})"

    def to_dict(self) -> Dict[str, Any]:
        """Convert TableInfo to dictionary for JSON serialization."""
        return {
            "name": self._info.name,
            "path": self._info.path,
            "columns": self._info.columns.copy() if self._info.columns else [],
            "description": self._info.description,
        }


def resolve_table_reference(value, parameter: str = "value"):
    """Normalize a documented column reference to a str or TableReference.

    Accepts a literal string, a TableReference (``tool.tables.x.col``), a
    ``(TableInfo, "col")`` tuple, or a ``(path, "col")`` pair. ``None`` passes
    through. A TableReference serializes to TABLE_REFERENCE:path:column, which
    is what pipe scripts resolve per id.

    A ``{chain: selection}`` dict says which residues belong to which chain, which is the
    only unambiguous way to express a selection over a multi-chain structure. Each value
    is resolved by the same rules, so a chain may carry a literal or a column reference.
    """
    if value is None or isinstance(value, (str, TableReference)):
        return value
    if isinstance(value, dict):
        resolved = {}
        for chain, selection in value.items():
            if not isinstance(chain, str) or len(chain) != 1 or not chain.isalnum():
                raise ValueError(
                    f"{parameter} dict keys must be single-character chain ids, got {chain!r}")
            resolved[chain] = resolve_table_reference(selection, f"{parameter}[{chain!r}]")
        return resolved
    if isinstance(value, tuple) and len(value) == 2:
        table, column = value
        if not isinstance(column, str):
            raise ValueError(
                f"{parameter} tuple must be (TableInfo|path, column_name), "
                f"got column of type {type(column).__name__}")
        # TableInfo routes bare attribute access to TableReference, so the path is only on .info.
        if isinstance(table, TableInfo):
            return TableReference(table.info.path, column)
        if isinstance(table, str):
            return TableReference(table, column)
        raise ValueError(
            f"{parameter} tuple's first element must be a TableInfo or path string, "
            f"got {type(table).__name__}")
    raise ValueError(
        f"{parameter} must be a string, a TableReference "
        f"(e.g. tool.tables.structures.designed), or a (TableInfo, column) tuple; "
        f"got {type(value).__name__}")


class IndexedTableContainer:
    """A collection of TableInfo objects indexed by ID, all sharing the same schema.

    Used when a tool produces one table per input ID (e.g., per-structure RMSF).
    Integrates into TableContainer as a single named entry.

    Usage:
        rmsf = IndexedTableContainer(
            name="rmsf",
            columns=["id", "chain", "resi", "rmsf"],
            description="Per-residue RMSF"
        )
        rmsf.add("1a2j", path="/output/1a2j_RMSF.csv")
        rmsf.add("1gfl", path="/output/1gfl_RMSF.csv")

        rmsf["1a2j"]              # -> TableInfo
        rmsf["1a2j"].info.path    # -> path string

        for struct_id, table_info in rmsf:
            print(struct_id, table_info.info.path)
    """

    def __init__(self, name: str, columns: List[str] = None,
                 description: str = ""):
        self.name = name
        self.columns = columns or []
        self.description = description
        self._entries: Dict[str, TableInfo] = {}

    def add(self, entry_id: str, path: str) -> 'IndexedTableContainer':
        """Add a per-ID table entry. Returns self for chaining."""
        self._entries[entry_id] = TableInfo(
            name=f"{self.name}_{entry_id}",
            path=path,
            columns=self.columns.copy(),
            description=f"{self.description} ({entry_id})",
        )
        return self

    def __getitem__(self, entry_id: str) -> TableInfo:
        """Get TableInfo by ID."""
        if entry_id in self._entries:
            return self._entries[entry_id]
        raise KeyError(f"No entry for ID '{entry_id}' in indexed table '{self.name}'. "
                       f"Available IDs: {list(self._entries.keys())}")

    def get(self, entry_id: str, default=None) -> Optional[TableInfo]:
        """Get TableInfo by ID with default."""
        return self._entries.get(entry_id, default)

    def __contains__(self, entry_id: str) -> bool:
        return entry_id in self._entries

    def __iter__(self):
        """Iterate over (id, TableInfo) pairs."""
        return iter(self._entries.items())

    def __len__(self) -> int:
        return len(self._entries)

    @property
    def ids(self) -> List[str]:
        """All entry IDs."""
        return list(self._entries.keys())

    def keys(self):
        """All entry IDs."""
        return self._entries.keys()

    def values(self):
        """All TableInfo objects."""
        return self._entries.values()

    def items(self):
        """All (id, TableInfo) pairs."""
        return self._entries.items()

    def to_dict(self) -> Dict[str, Any]:
        """Serialize for JSON persistence."""
        return {
            "_type": "IndexedTableContainer",
            "name": self.name,
            "columns": self.columns.copy(),
            "description": self.description,
            "entries": {
                entry_id: table_info.to_dict()
                for entry_id, table_info in self._entries.items()
            }
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'IndexedTableContainer':
        """Deserialize from JSON dict."""
        container = cls(
            name=data["name"],
            columns=data.get("columns", []),
            description=data.get("description", "")
        )
        for entry_id, info_dict in data.get("entries", {}).items():
            container._entries[entry_id] = TableInfo(
                name=info_dict.get("name", f"{data['name']}_{entry_id}"),
                path=info_dict.get("path", ""),
                columns=info_dict.get("columns", []),
                description=info_dict.get("description", ""),
            )
        return container

    def __str__(self) -> str:
        ids_display = ', '.join(self.ids[:5])
        if len(self.ids) > 5:
            ids_display += f', ... ({len(self.ids)} total)'
        col_display = ', '.join(self.columns) if self.columns else 'unknown'
        return f"${self.name}[{ids_display}] ({col_display})"

    def __repr__(self) -> str:
        return f"IndexedTableContainer(name='{self.name}', ids={self.ids}, columns={self.columns})"


class TableContainer:
    """Container for named tables with dot-notation access."""
    
    def __init__(self, tables: Dict[str, TableInfo]):
        self._tables = tables

        # Set attributes for dot notation access to TableInfo objects (not just paths)
        # Note: For IDE autocompletion, attributes should also be explicitly set in calling code
        for name, info in tables.items():
            setattr(self, name, info)
    
    def __getitem__(self, key: str):
        """Get table path by name, or IndexedTableContainer directly."""
        if key in self._tables:
            entry = self._tables[key]
            if isinstance(entry, IndexedTableContainer):
                return entry
            return entry.info.path

        raise KeyError(f"No table named '{key}' in tables")
    
    def __getattr__(self, name: str) -> TableInfo:
        """Get TableInfo object by name via dot notation."""
        if name in self._tables:
            return self._tables[name]
        raise AttributeError(f"No table named '{name}'")
    
    def keys(self):
        """Get all table names."""
        return self._tables.keys()

    def items(self):
        """Get all name, path/container pairs."""
        result = []
        for name, entry in self._tables.items():
            if isinstance(entry, IndexedTableContainer):
                result.append((name, entry))
            else:
                result.append((name, entry.info.path))
        return result
    
    def __contains__(self, key: str) -> bool:
        """Support 'in' operator: 'table_name' in tables"""
        return key in self._tables
    
    def get(self, key: str, default: str = "") -> str:
        """Get table path with default (like dict.get())."""
        try:
            return self.__getitem__(key)
        except KeyError:
            return default
    
    def __str__(self) -> str:
        if not self._tables:
            return "{}"

        lines = []
        for name, entry in self._tables.items():
            if isinstance(entry, IndexedTableContainer):
                lines.append(f"    {str(entry)}:")
                for eid, tinfo in entry:
                    filename = tinfo.info.path.split('/')[-1]
                    lines.append(f"        [{eid}] '<output_folder>/{filename}'")
            else:
                # Format: table name with columns, then indented file path
                # Extract just the filename from the full path
                filename = entry.info.path.split('/')[-1]
                path_display = f"<output_folder>/{filename}"
                lines.append(f"    {str(entry)}:")
                lines.append(f"        – '{path_display}'")
        
        return "\n".join(lines)
    
    def __repr__(self) -> str:
        return f"TableContainer({list(self._tables.keys())})"


class StreamContainer:
    """Container for named DataStreams with dot-notation access."""

    def __init__(self, streams: Dict[str, Any]):
        self._streams = streams

        # Set attributes for dot notation access to DataStream objects
        for name, stream in streams.items():
            setattr(self, name, stream)

    def __getitem__(self, key: str):
        """Get DataStream by name."""
        if key in self._streams:
            return self._streams[key]
        raise KeyError(f"No stream named '{key}' in streams")

    def __getattr__(self, name: str):
        """Get DataStream by name via dot notation. Returns None if not present."""
        if '_streams' in self.__dict__ and name in self._streams:
            return self._streams[name]
        return None

    def keys(self):
        """Get all stream names."""
        return self._streams.keys()

    def items(self):
        """Get all name, stream pairs."""
        return self._streams.items()

    def __contains__(self, key: str) -> bool:
        """Support 'in' operator: 'stream_name' in streams"""
        return key in self._streams

    def get(self, key: str, default=None):
        """Get stream with default (like dict.get())."""
        return self._streams.get(key, default)

    def __len__(self) -> int:
        """Return number of streams."""
        return len(self._streams)

    def __str__(self) -> str:
        if not self._streams:
            return "{}"
        return f"StreamContainer({list(self._streams.keys())})"

    def __repr__(self) -> str:
        return f"StreamContainer({list(self._streams.keys())})"
