# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
DataStream: Unified container for tool I/O data.

Provides a consistent interface for structures, sequences, compounds, and MSAs
regardless of their underlying representation (files, SMILES strings, etc.).

Supports pattern-based IDs for compact representation:
    ids=["5HG6_<0..49>"]  →  lazily expands to 50 IDs
    files=["<id>.pdb"]    →  lazily expands using expanded IDs
"""

import os
import pandas as pd
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Set, Iterator, Tuple, Union

from . import id_patterns


@dataclass
class DataStream:
    """
    Unified container for any data type flowing between tools.

    A DataStream represents a collection of items (structures, sequences, compounds, etc.)
    with a consistent interface regardless of storage format.

    Attributes:
        name: Name of this data stream (e.g., "structures", "sequences", "compounds")
        ids: List of identifiers, may contain patterns (e.g., ["5HG6_<0..49>"])
        files: List of file paths, may contain "<id>" template or "*" globs
        map_table: Path to CSV file mapping ids to files/values and additional metadata
        format: Data format ("pdb", "cif", "fasta", "csv", "sdf", "smiles", "ccd", etc.)
        metadata: Additional tool-specific metadata

    The map_table CSV always contains at minimum:
        - id: The item identifier
        - file: Path to file (empty string if format is inline like "smiles")
        - value: Inline value (empty string if format is file-based)

    Examples:
        # Structure data (file-based)
        structures = DataStream(
            name="structures",
            ids=["protein_1", "protein_2"],
            files=["/path/to/protein_1.pdb", "/path/to/protein_2.pdb"],
            map_table="/path/to/structures_map.csv",
            format="pdb"
        )

        # Pattern-based IDs (compact representation)
        designs = DataStream(
            name="structures",
            ids=["5HG6_<0..49>"],
            files=[os.path.join(folder, "<id>.pdb")],
            format="pdb"
        )
        # len(designs) == 50, but ids stores only 1 pattern string

        # Iteration yields single-item DataStreams that can be passed to tools
        for structure in structures:
            print(f"Processing {structure.ids[0]}: {structure.files[0]}")
            tool = SomeTool(structures=structure)
    """

    name: str = ""
    ids: List[str] = field(default_factory=list)
    files: List[str] = field(default_factory=list)
    map_table: str = ""
    format: str = "pdb"
    metadata: Dict[str, Any] = field(default_factory=dict)

    # Internal caches (not serialized, not compared)
    _ids_expanded: Optional[List[str]] = field(default=None, repr=False, compare=False)
    _files_expanded: Optional[List[str]] = field(default=None, repr=False, compare=False)
    _map_data: Optional[pd.DataFrame] = field(default=None, repr=False, compare=False)

    def __post_init__(self):
        """Validate DataStream after initialization."""
        # Skip validation when patterns or <id> template are present
        if self.has_patterns() or self._has_file_template():
            return
        # Validation rules for files:
        # - Empty files list: valid (value-based, data in map_table)
        # - Single file: valid (one file contains all records, e.g., CSV)
        # - Multiple files: must match number of ids (one file per record)
        if len(self.files) > 1 and len(self.files) != len(self.ids):
            raise ValueError(
                f"Length mismatch: {len(self.ids)} ids but {len(self.files)} files. "
                f"Use empty files list or single file for table-based data, "
                f"or provide one file per id."
            )

    def _invalidate_caches(self):
        """Clear expansion caches (call when ids or files change)."""
        self._ids_expanded = None
        self._files_expanded = None

    def has_patterns(self) -> bool:
        """True if any ID contains a pattern slot."""
        return any(id_patterns.contains_pattern(s) for s in self.ids)

    @property
    def is_lazy(self) -> bool:
        """True if any ID contains [...] bracket segments (runtime-dependent)."""
        return any(id_patterns.is_lazy(s) for s in self.ids)

    def _has_file_template(self) -> bool:
        """True if files list uses <id> template substitution."""
        return len(self.files) == 1 and '<id>' in self.files[0]

    @property
    def ids_expanded(self) -> List[str]:
        """Lazily expand patterns in ids. Cached after first call.

        For lazy patterns (with brackets), expands only the deterministic prefix.
        """
        if self._ids_expanded is None:
            if self.is_lazy:
                self._ids_expanded, _ = id_patterns.try_expand_ids(self.ids)
            elif self.has_patterns():
                self._ids_expanded = id_patterns.expand_ids(self.ids)
            else:
                self._ids_expanded = self.ids  # No copy — same list reference
        return self._ids_expanded

    @property
    def files_expanded(self) -> List[str]:
        """Lazily expand file patterns using expanded IDs."""
        if self._files_expanded is None:
            if not self.files:
                self._files_expanded = []
            elif self._has_file_template():
                template = self.files[0]
                self._files_expanded = [
                    id_patterns.expand_file_pattern(template, eid)
                    for eid in self.ids_expanded
                ]
            else:
                self._files_expanded = self.files  # Already explicit
        return self._files_expanded

    def __len__(self) -> int:
        """
        Return number of items in the stream.

        Uses pattern counting when possible to avoid full expansion.
        """
        if self._ids_expanded is not None:
            return len(self._ids_expanded)
        return id_patterns.count_ids(self.ids)

    def __iter__(self) -> Iterator['DataStream']:
        """
        Iterate over single-item DataStream objects.

        Each yielded DataStream has the same name and format as the parent,
        with ids containing a single ID and files containing the corresponding
        file (for file-based streams) or empty (for value-based streams).

        Yields:
            Single-item DataStream for each item
        """
        expanded_ids = self.ids_expanded
        expanded_files = self.files_expanded
        if expanded_files:
            for item_id, item_file in zip(expanded_ids, expanded_files):
                yield DataStream(
                    name=self.name,
                    ids=[item_id],
                    files=[item_file],
                    format=self.format
                )
        else:
            for item_id in expanded_ids:
                yield DataStream(
                    name=self.name,
                    ids=[item_id],
                    files=[],
                    format=self.format
                )

    def __getitem__(self, index):
        """Get item by index or slice.

        For integer index: returns single-item DataStream.
        For slice: returns new DataStream with sliced items.
        """
        if isinstance(index, slice):
            sliced_ids = self.ids_expanded[index]
            sliced_files = self.files_expanded[index] if self.files_expanded else []
            return DataStream(
                name=self.name,
                ids=sliced_ids,
                files=sliced_files,
                map_table=self.map_table,
                format=self.format
            )

        # Integer indexing — try O(1) expand_at for single-pattern case
        n = len(self)
        if index < 0 or index >= n:
            raise IndexError(f"Index {index} out of range for DataStream with {n} items")

        if len(self.ids) == 1 and self.has_patterns():
            item_id = id_patterns.expand_at(self.ids[0], index)
        else:
            item_id = self.ids_expanded[index]

        if self.files_expanded:
            return DataStream(
                name=self.name,
                ids=[item_id],
                files=[self.files_expanded[index]],
                format=self.format
            )
        else:
            return DataStream(
                name=self.name,
                ids=[item_id],
                files=[],
                format=self.format
            )

    def __bool__(self) -> bool:
        """DataStream is truthy if it has any items."""
        return len(self.ids) > 0

    def _get_map_data(self) -> Optional[pd.DataFrame]:
        """Lazily load map_table data."""
        if self._map_data is None and self.map_table and os.path.exists(self.map_table):
            self._map_data = pd.read_csv(self.map_table)
        return self._map_data

    def get_file(self, item_id: str) -> Optional[str]:
        """
        Get file path for a specific ID.

        Args:
            item_id: The item identifier

        Returns:
            File path if found, None otherwise
        """
        expanded_ids = self.ids_expanded
        expanded_files = self.files_expanded
        if item_id in expanded_ids and expanded_files:
            idx = expanded_ids.index(item_id)
            return expanded_files[idx]
        return None

    def get_value(self, item_id: str, column: str = 'value') -> Optional[str]:
        """
        Get a value from the map_table for a specific ID.

        Args:
            item_id: The item identifier
            column: Column name to retrieve (default: 'value')

        Returns:
            Value if found, None otherwise
        """
        map_data = self._get_map_data()
        if map_data is not None and column in map_data.columns:
            row = map_data[map_data['id'] == item_id]
            if not row.empty:
                return row.iloc[0][column]
        return None

    def filter_by_ids(self, keep_ids: Union[Set[str], List[str]]) -> 'DataStream':
        """
        Create a new DataStream containing only specified IDs.

        Triggers expansion (unavoidable — must check each ID).
        """
        if isinstance(keep_ids, list):
            keep_ids = set(keep_ids)

        expanded_ids = self.ids_expanded
        expanded_files = self.files_expanded

        new_ids = []
        new_files = []

        for i, item_id in enumerate(expanded_ids):
            if item_id in keep_ids:
                new_ids.append(item_id)
                if expanded_files:
                    new_files.append(expanded_files[i])

        return DataStream(
            name=self.name,
            ids=new_ids,
            files=new_files,
            map_table=self.map_table,
            format=self.format,
            metadata={**self.metadata, '_filtered': True, '_original_count': len(self)}
        )

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert to dictionary for serialization.

        Stores pattern-based ids/files (compact representation).
        """
        return {
            'name': self.name,
            'ids': self.ids.copy(),
            'files': self.files.copy(),
            'map_table': self.map_table,
            'format': self.format,
            'metadata': self.metadata.copy()
        }

    def save_json(self, path: str) -> str:
        """Save DataStream to JSON file with compact patterns.

        Patterns stay compact; pipe scripts expand at runtime via DataStreamRuntime.
        """
        import json
        os.makedirs(os.path.dirname(path), exist_ok=True)
        data = {
            'name': self.name,
            'ids': self.ids.copy(),
            'files': self.files.copy(),
            'map_table': self.map_table,
            'format': self.format,
            'metadata': self.metadata.copy()
        }
        with open(path, 'w') as f:
            json.dump(data, f, indent=2)
        return path

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'DataStream':
        """
        Create DataStream from dictionary.

        Accepts both pattern and explicit ids. Silently ignores
        legacy 'files_contain_wildcards' key for backward compat.
        """
        return cls(
            name=data.get('name', ''),
            ids=data.get('ids', []),
            files=data.get('files', []),
            map_table=data.get('map_table', ''),
            format=data.get('format', 'pdb'),
            metadata=data.get('metadata', {})
        )

    @classmethod
    def empty(cls, name: str, format: str) -> 'DataStream':
        """Create an empty DataStream."""
        return cls(name=name, ids=[], files=[], map_table="", format=format, metadata={})

    def is_file_based(self) -> bool:
        """
        Check if this DataStream uses file-based storage.

        Returns:
            True if format uses files (pdb, cif, sdf, fasta, etc.)
            False if format uses inline values (smiles, ccd)
        """
        value_based_formats = {'smiles', 'ccd', 'sequence'}
        return self.format.lower() not in value_based_formats

    def __repr__(self) -> str:
        """Detailed string representation."""
        return (
            f"DataStream(name='{self.name}', format='{self.format}', "
            f"items={len(self)}, "
            f"files={len(self.files)}, "
            f"map_table={'set' if self.map_table else 'unset'})"
        )

    def __str__(self) -> str:
        """Human-readable string representation."""
        name_part = self.name or "unnamed"
        if len(self) == 0:
            return f"DataStream({name_part}, {self.format}): empty"

        if self.has_patterns():
            id_preview = ", ".join(self.ids)
            return f"DataStream({name_part}, {self.format}): [{id_preview}] ({len(self)} total)"

        if len(self) <= 3:
            id_preview = ", ".join(self.ids)
        else:
            id_preview = f"{self.ids[0]}, {self.ids[1]}, ... ({len(self)} total)"

        return f"DataStream({name_part}, {self.format}): [{id_preview}]"


def create_map_table(
    output_path: str,
    ids: List[str],
    files: Optional[List[str]] = None,
    values: Optional[List[str]] = None,
    additional_columns: Optional[Dict[str, List[Any]]] = None,
    provenance: Optional[Dict[str, List[str]]] = None,
    parent_ids: Optional[List[str]] = None,
    suffix_pattern: Optional[str] = None,
    input_stream_name: Optional[str] = None
) -> str:
    """
    Create a map_table CSV file for a DataStream.

    Accepts pattern-based ids/files and expands them internally.

    Args:
        output_path: Path where CSV will be written
        ids: List of item identifiers (may contain patterns)
        files: List of file paths (optional, may contain "<id>" template)
        values: List of inline values (optional, for value-based formats like SMILES)
        additional_columns: Dict of column_name -> values to include
        provenance: Dict of {stream_name: [input_id_per_output, ...]}
        parent_ids: Parent IDs (may contain patterns) for auto-provenance
        suffix_pattern: Suffix pattern (e.g., "<1..3>") — when set along with
            parent_ids and input_stream_name, auto-generates provenance by
            stripping the suffix from each expanded ID.
        input_stream_name: Name of the input stream for auto-provenance column

    Returns:
        Path to created CSV file
    """
    # Expand pattern-based ids (skip if lazy — lazy IDs are expanded at runtime)
    if any(id_patterns.is_lazy(s) for s in ids):
        expanded_ids, _ = id_patterns.try_expand_ids(ids)
    elif any(id_patterns.contains_pattern(s) for s in ids):
        expanded_ids = id_patterns.expand_ids(ids)
    else:
        expanded_ids = ids

    # Expand file template if needed
    if files and len(files) == 1 and '<id>' in files[0]:
        template = files[0]
        expanded_files = [id_patterns.expand_file_pattern(template, eid) for eid in expanded_ids]
    elif files and any(id_patterns.contains_pattern(s) for s in files):
        expanded_files = id_patterns.expand_ids(files)
    else:
        expanded_files = files

    data = {'id': expanded_ids}

    if expanded_files:
        data['file'] = expanded_files
    else:
        data['file'] = [''] * len(expanded_ids)

    if values:
        data['value'] = values
    else:
        data['value'] = [''] * len(expanded_ids)

    if additional_columns:
        for col_name, col_values in additional_columns.items():
            if len(col_values) != len(expanded_ids):
                raise ValueError(
                    f"Column '{col_name}' has {len(col_values)} values "
                    f"but expected {len(expanded_ids)}"
                )
            data[col_name] = col_values

    # Auto-generate provenance from parent_ids + suffix_pattern
    if parent_ids is not None and suffix_pattern and input_stream_name:
        suffix_values = id_patterns.expand_pattern(suffix_pattern) if id_patterns.contains_pattern(suffix_pattern) else [suffix_pattern]
        n_suffixes = len(suffix_values)
        expanded_parents = id_patterns.expand_ids(parent_ids) if any(id_patterns.contains_pattern(s) for s in parent_ids) else parent_ids
        prov_list = []
        for pid in expanded_parents:
            prov_list.extend([pid] * n_suffixes)
        col_name = f"{input_stream_name}.id"
        data[col_name] = prov_list

    if provenance:
        for stream_name, prov_ids in provenance.items():
            col_name = f"{stream_name}.id"
            if len(prov_ids) != len(expanded_ids):
                raise ValueError(
                    f"Provenance column '{col_name}' has {len(prov_ids)} values "
                    f"but expected {len(expanded_ids)}"
                )
            data[col_name] = prov_ids

    df = pd.DataFrame(data)

    # Ensure parent directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    df.to_csv(output_path, index=False)
    return output_path
