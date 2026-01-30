"""
DataStream: Unified container for tool I/O data.

Provides a consistent interface for structures, sequences, compounds, and MSAs
regardless of their underlying representation (files, SMILES strings, etc.).
"""

import os
import pandas as pd
from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Set, Iterator, Tuple, Union


@dataclass
class DataStream:
    """
    Unified container for any data type flowing between tools.

    A DataStream represents a collection of items (structures, sequences, compounds, etc.)
    with a consistent interface regardless of storage format.

    Attributes:
        name: Name of this data stream (e.g., "structures", "sequences", "compounds")
        ids: List of unique identifiers for each item
        files: List of file paths (may be empty for inline formats like SMILES)
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

        # Compound data with SMILES (value-based, no files)
        compounds = DataStream(
            name="compounds",
            ids=["ligand_1", "ligand_2"],
            files=[],  # Empty for SMILES format
            map_table="/path/to/compounds_map.csv",  # Contains id,smiles columns
            format="smiles"
        )

        # Iteration works uniformly
        for item_id, item_path in structures:
            print(f"Processing {item_id}: {item_path}")
    """

    name: str = ""
    ids: List[str] = field(default_factory=list)
    files: List[str] = field(default_factory=list)
    map_table: str = ""
    format: str = "pdb"
    metadata: Dict[str, Any] = field(default_factory=dict)

    # Internal cache for map_table data (loaded lazily)
    _map_data: Optional[pd.DataFrame] = field(default=None, repr=False, compare=False)

    def __post_init__(self):
        """Validate DataStream after initialization."""
        # For file-based formats, files and ids should match in length
        if self.files and len(self.files) != len(self.ids):
            raise ValueError(
                f"Length mismatch: {len(self.ids)} ids but {len(self.files)} files. "
                f"For value-based formats (like 'smiles'), leave files empty."
            )

    def __len__(self) -> int:
        """
        Return number of items in the stream.

        This is always the number of IDs, not files, because:
        - File-based formats: one file per ID
        - Value-based formats (SMILES): zero files, values in map_table
        """
        return len(self.ids)

    def __iter__(self) -> Iterator[Tuple[str, str]]:
        """
        Iterate over (id, file_or_value) pairs.

        For file-based formats, yields (id, file_path).
        For value-based formats (like SMILES), yields (id, value_from_map_table).

        Yields:
            Tuple of (item_id, item_path_or_value)
        """
        if self.files:
            # File-based format
            yield from zip(self.ids, self.files)
        else:
            # Value-based format - get values from map_table
            map_data = self._get_map_data()
            if map_data is not None and 'value' in map_data.columns:
                for item_id in self.ids:
                    row = map_data[map_data['id'] == item_id]
                    if not row.empty:
                        yield (item_id, row.iloc[0]['value'])
                    else:
                        yield (item_id, "")
            else:
                # Fallback: yield empty values
                for item_id in self.ids:
                    yield (item_id, "")

    def __getitem__(self, index: int) -> Tuple[str, str]:
        """Get item by index as (id, file_or_value) tuple."""
        if index < 0 or index >= len(self.ids):
            raise IndexError(f"Index {index} out of range for DataStream with {len(self.ids)} items")

        item_id = self.ids[index]
        if self.files:
            return (item_id, self.files[index])
        else:
            # Value-based format
            map_data = self._get_map_data()
            if map_data is not None and 'value' in map_data.columns:
                row = map_data[map_data['id'] == item_id]
                if not row.empty:
                    return (item_id, row.iloc[0]['value'])
            return (item_id, "")

    def __bool__(self) -> bool:
        """DataStream is truthy if it has any items."""
        return len(self.ids) > 0

    def _get_map_data(self) -> Optional[pd.DataFrame]:
        """Lazily load map_table data."""
        if self._map_data is None and self.map_table and os.path.exists(self.map_table):
            try:
                self._map_data = pd.read_csv(self.map_table)
            except Exception:
                self._map_data = None
        return self._map_data

    def get_file(self, item_id: str) -> Optional[str]:
        """
        Get file path for a specific ID.

        Args:
            item_id: The item identifier

        Returns:
            File path if found, None otherwise
        """
        if item_id in self.ids and self.files:
            idx = self.ids.index(item_id)
            return self.files[idx]
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

        Args:
            keep_ids: Set or list of IDs to keep

        Returns:
            New DataStream with filtered items

        Example:
            passing_ids = {"protein_1", "protein_3"}
            filtered = structures.filter_by_ids(passing_ids)
        """
        if isinstance(keep_ids, list):
            keep_ids = set(keep_ids)

        # Filter ids and files together
        new_ids = []
        new_files = []

        for i, item_id in enumerate(self.ids):
            if item_id in keep_ids:
                new_ids.append(item_id)
                if self.files:
                    new_files.append(self.files[i])

        # Note: map_table path stays the same but represents different subset
        # A filtered map_table should be created by the tool if needed
        return DataStream(
            name=self.name,
            ids=new_ids,
            files=new_files,
            map_table=self.map_table,  # Original map still valid for lookups
            format=self.format,
            metadata={**self.metadata, '_filtered': True, '_original_count': len(self.ids)}
        )

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert to dictionary for serialization.

        Returns:
            Dictionary representation suitable for JSON serialization
        """
        return {
            'name': self.name,
            'ids': self.ids.copy(),
            'files': self.files.copy(),
            'map_table': self.map_table,
            'format': self.format,
            'metadata': self.metadata.copy()
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'DataStream':
        """
        Create DataStream from dictionary.

        Args:
            data: Dictionary with DataStream fields

        Returns:
            New DataStream instance
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
        """
        Create an empty DataStream.

        Args:
            name: Name for the empty stream
            format: Data format for the empty stream

        Returns:
            Empty DataStream instance
        """
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
    additional_columns: Optional[Dict[str, List[Any]]] = None
) -> str:
    """
    Create a map_table CSV file for a DataStream.

    Args:
        output_path: Path where CSV will be written
        ids: List of item identifiers
        files: List of file paths (optional, for file-based formats)
        values: List of inline values (optional, for value-based formats like SMILES)
        additional_columns: Dict of column_name -> values to include

    Returns:
        Path to created CSV file

    Example:
        # For structures (file-based)
        create_map_table(
            "/path/to/structures_map.csv",
            ids=["prot_1", "prot_2"],
            files=["/path/prot_1.pdb", "/path/prot_2.pdb"]
        )

        # For SMILES (value-based)
        create_map_table(
            "/path/to/compounds_map.csv",
            ids=["lig_1", "lig_2"],
            values=["CCO", "CC(=O)O"],
            additional_columns={"name": ["ethanol", "acetic acid"]}
        )
    """
    data = {'id': ids}

    if files:
        data['file'] = files
    else:
        data['file'] = [''] * len(ids)

    if values:
        data['value'] = values
    else:
        data['value'] = [''] * len(ids)

    if additional_columns:
        for col_name, col_values in additional_columns.items():
            if len(col_values) != len(ids):
                raise ValueError(
                    f"Column '{col_name}' has {len(col_values)} values "
                    f"but expected {len(ids)}"
                )
            data[col_name] = col_values

    df = pd.DataFrame(data)

    # Ensure parent directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    df.to_csv(output_path, index=False)
    return output_path
