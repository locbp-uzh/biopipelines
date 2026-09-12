# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Tool output views: what a tool hands downstream once its files are declared.

StandardizedOutput gives dot-notation access to a tool's output keys, indexing by id, chunking,
rendering and download; ToolOutput is the per-tool wrapper around a config's declared outputs.
Both carry the producer back-reference the pipeline walks to recover dataflow edges.

Layered between data_containers.py and base_config.py, which re-exports these names, so importing them from base_config keeps working.
"""

import os
from typing import Any, Dict, List, Optional, Union, TYPE_CHECKING

try:
    from .data_containers import (
        TableInfo,
        IndexedTableContainer,
        TableContainer,
        StreamContainer,
    )
    from .biopipelines_io import TableReference
    from .config_manager import ConfigManager
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from data_containers import (
        TableInfo,
        IndexedTableContainer,
        TableContainer,
        StreamContainer,
    )
    from biopipelines_io import TableReference
    from config_manager import ConfigManager

if TYPE_CHECKING:
    from .base_config import BaseConfig


class StandardizedOutput:
    """
    Provides dot-notation access to standardized output keys.

    DataStreams are accessed via the streams container:

        for structure in output.streams.structures:
            print(f"Processing {structure.ids[0]}: {structure.files[0]}")

        print(f"Generated {len(output.streams.structures)} structures")

    If all streams share the same IDs, the output itself is iterable,
    yielding single-item StandardizedOutput objects:

        for item in output:
            tool = SomeTool(structures=item)  # item.streams has one id per stream
    """

    def __init__(self, output_files: Dict[str, Any]):
        """Initialize with output files dictionary."""
        self._data = output_files.copy()

        # Set before the trailing setattr loop so a stray data key of the same name is skipped.
        self._producer = None

        # Handle tables - convert to TableInfo objects if needed
        tables_raw = output_files.get('tables', [])
        self.tables = self._process_tables(tables_raw)

        # Build streams container from all DataStream objects in output_files
        from .datastream import DataStream
        streams_dict = {}
        for key, value in output_files.items():
            if isinstance(value, DataStream):
                streams_dict[key] = value
        self.streams = StreamContainer(streams_dict)

        self.output_folder = output_files.get('output_folder', '')

        # Handle filtering metadata
        self.filter_metadata = output_files.get('filter_metadata', {})
        self.is_filtered = self.filter_metadata.get('is_filtered', False)

        # Store additional attributes (non-DataStream, non-reserved)
        reserved_keys = {'tables', 'output_folder', 'filter_metadata', 'streams'}
        for key, value in output_files.items():
            if key not in reserved_keys and not isinstance(value, DataStream) and not hasattr(self, key):
                setattr(self, key, value)
    
    @property
    def producer(self) -> Optional['BaseConfig']:
        """The tool config that produced this output, or None if unknown."""
        return self._producer

    def _stamp_producer(self, config: 'BaseConfig') -> 'StandardizedOutput':
        """Record ``config`` as this output's producer and stamp its streams.

        Streams are stamped only when unclaimed: a tool that re-exports the very DataStream object it was handed produces nothing new, and the first producer is the right answer for it.
        """
        self._producer = config
        self.streams.__dict__["_producer"] = config
        from .datastream import DataStream
        for stream in self.streams._streams.values():
            if isinstance(stream, DataStream) and getattr(stream, "_producer", None) is None:
                stream._producer = config
        return self

    def _inherit_producer(self, parent: 'StandardizedOutput') -> 'StandardizedOutput':
        """Carry ``parent``'s producer onto this derived output and its fresh streams."""
        if parent._producer is not None:
            self._stamp_producer(parent._producer)
        return self

    def _process_tables(self, tables_raw: Any) -> 'TableContainer':
        """Process tables into TableContainer with named access."""
        if isinstance(tables_raw, dict):
            # Already in named format
            table_infos = {}
            for name, info in tables_raw.items():
                if isinstance(info, IndexedTableContainer):
                    table_infos[name] = info
                elif isinstance(info, dict):
                    # Check for serialized IndexedTableContainer
                    if info.get('_type') == 'IndexedTableContainer':
                        table_infos[name] = IndexedTableContainer.from_dict(info)
                    else:
                        # Handle columns - ensure it's always a list
                        columns = info.get('columns', [])
                        if isinstance(columns, str):
                            # No placeholders - if it's a string, it should be an actual column name
                            columns = [columns]  # Single string becomes list

                        table_infos[name] = TableInfo(
                            name=name,
                            path=info.get('path', ''),
                            columns=columns,
                            description=info.get('description', ''),
                        )
                elif isinstance(info, TableInfo):
                    table_infos[name] = info
                else:
                    # Legacy format - just a path
                    table_infos[name] = TableInfo(name=name, path=str(info))
            return TableContainer(table_infos)
        elif isinstance(tables_raw, list):
            # Legacy format - convert first item to "main"
            table_infos = {}
            for i, path in enumerate(tables_raw):
                name = "main" if i == 0 else f"sheet_{i}"
                if isinstance(path, str):
                    table_infos[name] = TableInfo(name=name, path=path)
                else:
                    # Assume it's already a TableInfo or dict
                    table_infos[name] = path
            return TableContainer(table_infos)
        else:
            return TableContainer({})

    def __matmul__(self, column: str) -> 'TableReference':
        """Column lookup sugar: ``output @ "col"`` finds the table declaring
        ``col`` and returns its ``TableReference``, same object as
        ``output.tables.<name>.col``.

        Errors fast: raises if no declared table has the column, if more than
        one does (caller must use ``output.tables.<name>.col``), or if the only
        match is a per-ID ``IndexedTableContainer`` (caller must index by ID
        first, then access the column).
        """
        if not isinstance(column, str):
            raise TypeError(
                f"Column lookup expects a string, got {type(column).__name__}"
            )

        matches = []
        indexed_only = []
        for name, entry in self.tables._tables.items():
            if isinstance(entry, IndexedTableContainer):
                if column in entry.columns:
                    indexed_only.append(name)
                continue
            if column in entry.info.columns:
                matches.append(name)

        if len(matches) == 1:
            return getattr(self.tables._tables[matches[0]], column)
        if len(matches) > 1:
            raise ValueError(
                f"Column '{column}' is declared in tables {matches}; "
                f"use output.tables.<name>.{column} to disambiguate."
            )
        if indexed_only:
            raise ValueError(
                f"Column '{column}' belongs to per-ID indexed table(s) {indexed_only}; "
                f"index by ID first, e.g. output.tables.{indexed_only[0]}[id].{column}."
            )
        declared = {
            name: (entry.columns if isinstance(entry, IndexedTableContainer)
                   else entry.info.columns)
            for name, entry in self.tables._tables.items()
        }
        raise ValueError(
            f"No declared table has a column '{column}'. Declared columns: {declared}"
        )

    def __getitem__(self, key: str):
        """
        ID-based selection: ``output["CP1"]`` returns a new ``StandardizedOutput``
        whose streams carry ``ids=[key]`` and the original streams' ``files``,
        ``map_table``, and ``format`` unchanged.

        ``key`` must appear in every non-empty stream's ``ids_expanded``.
        ``ids_expanded`` is mode-aware: deterministic streams enumerate
        fully at config time (typos fail fast); lazy streams at config
        time enumerate only the deterministic prefix (so only prefix-shaped
        keys pass); lazy streams at runtime (``_runtime_mode=True``) read
        from the materialized map_table. See developer_manual.md "ID Patterns".
        """
        from .datastream import DataStream

        if not isinstance(key, str):
            raise TypeError(
                f"StandardizedOutput key must be a string ID, got {type(key).__name__}"
            )

        streams = [
            (name, ds) for name, ds in self.streams.items()
            if isinstance(ds, DataStream) and len(ds) > 0
        ]

        if not streams:
            raise KeyError(f"No non-empty streams to select ID '{key}' from")

        for name, ds in streams:
            if key not in ds.ids_expanded:
                raise KeyError(
                    f"ID '{key}' not found in stream '{name}' "
                    f"(ids: {list(ds.ids_expanded)})"
                )

        single_data = {}
        for name, ds in streams:
            # Pick the right files shape for the sub-stream so the
            # constructor's "len(files) must match len(ids)" rule holds.
            # Templates and shared-file form already cover every id, so they
            # pass through unchanged; explicit per-id lists need slicing.
            if ds.is_shared_file or ds._has_file_template() or not ds.files:
                sub_files = ds.files
            else:
                # `files` is positional against whichever id list it was built from. Indexing an
                # expanded position into a pattern-space list returns a neighbor's file, silently.
                expanded = list(ds.ids_expanded)
                if len(ds.files) == len(expanded):
                    sub_files = [ds.files[expanded.index(key)]]
                else:
                    raise KeyError(
                        f"Cannot resolve a file for id '{key}' in stream '{name}': it declares "
                        f"{len(ds.files)} explicit file(s) against {len(expanded)} expanded id(s) "
                        f"({len(ds.ids)} before expansion), so no file corresponds to one id. "
                        "Give the stream a '<id>' file template, a shared file, or one explicit "
                        "file per expanded id."
                    )
            single_data[name] = DataStream(
                name=ds.name,
                ids=[key],
                files=sub_files,
                map_table=ds.map_table,
                format=ds.format
            )
        single_data["output_folder"] = self.output_folder
        single_data["tables"] = self._data.get("tables", [])
        return StandardizedOutput(single_data)._inherit_producer(self)
    
    def __iter__(self):
        """
        Iterate over single-item StandardizedOutput objects.

        Requires at least one stream and that all streams share the same IDs.
        Each yielded item is a StandardizedOutput with single-item DataStreams
        (one ID/file per stream), preserving the same stream names.

        Raises:
            ValueError: If there are no streams or streams have mismatched IDs.
        """
        from .datastream import DataStream

        # Collect all streams
        streams = [
            (name, ds) for name, ds in self.streams.items()
            if isinstance(ds, DataStream) and len(ds) > 0
        ]

        if not streams:
            raise ValueError(
                "Cannot iterate over StandardizedOutput: no non-empty streams. "
                "Iterate over a specific stream instead (e.g., output.streams.structures)."
            )

        # Verify all streams share the same IDs
        reference_name, reference_ds = streams[0]
        reference_ids = list(reference_ds.ids)
        for name, ds in streams[1:]:
            if list(ds.ids) != reference_ids:
                raise ValueError(
                    f"Cannot iterate over StandardizedOutput: streams have mismatched IDs. "
                    f"'{reference_name}' has {len(reference_ids)} IDs, "
                    f"'{name}' has {len(ds)} IDs. "
                    f"Iterate over a specific stream instead (e.g., output.streams.{reference_name})."
                )

        # Yield one StandardizedOutput per ID
        for idx in range(len(reference_ids)):
            single_streams = {}
            for name, ds in streams:
                if ds.is_shared_file:
                    sub_files = ds.files  # str — shared across ids
                elif isinstance(ds.files, list) and len(ds.files) > idx:
                    sub_files = [ds.files[idx]]
                else:
                    sub_files = []
                single_streams[name] = DataStream(
                    name=ds.name,
                    ids=[ds.ids[idx]],
                    files=sub_files,
                    map_table=ds.map_table,
                    format=ds.format
                )
            single_streams["output_folder"] = self.output_folder
            yield StandardizedOutput(single_streams)._inherit_producer(self)

    def chunks(self, count: Optional[int] = None, *, size: Optional[int] = None
               ) -> List['StandardizedOutput']:
        """Split aligned streams into groups of ``StandardizedOutput`` objects.

        Every stream must carry the same expanded IDs in the same order. Use a
        specific ``DataStream.chunks()`` when an output contains independent ID
        axes.
        """
        from .datastream import DataStream

        streams = [
            (name, ds) for name, ds in self.streams.items()
            if isinstance(ds, DataStream)
        ]
        if not streams:
            raise ValueError(
                "Cannot chunk StandardizedOutput: it contains no streams."
            )

        reference_name, reference_ds = streams[0]
        reference_ids = list(reference_ds.ids_expanded)
        for name, ds in streams[1:]:
            ids = list(ds.ids_expanded)
            if ids != reference_ids:
                raise ValueError(
                    "Cannot chunk StandardizedOutput: streams have mismatched IDs. "
                    f"'{reference_name}' has {len(reference_ids)} IDs, "
                    f"'{name}' has {len(ids)} IDs. "
                    f"Chunk a specific stream instead (e.g., "
                    f"output.streams.{reference_name}.chunks(...))."
                )

        chunked = [
            (name, ds.chunks(count, size=size)) for name, ds in streams
        ]
        n_chunks = len(chunked[0][1])
        result = []
        for index in range(n_chunks):
            chunk_data = self._data.copy()
            for name, chunks in chunked:
                chunk_data[name] = chunks[index]
            result.append(StandardizedOutput(chunk_data)._inherit_producer(self))
        return result
    
    def keys(self):
        """Get all available keys."""
        return self._data.keys()
    
    def items(self):
        """Get all key-value pairs."""
        return self._data.items()
    
    def get(self, key: str, default=None):
        """Get value with default."""
        return self._data.get(key, default)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert back to dictionary."""
        return self._data.copy()
    
    def __contains__(self, key: str) -> bool:
        """Support 'in' operator: 'structures' in output"""
        return key in self._data
    
    def get_legacy_table(self, name: str = "main") -> str:
        """
        Get table path using legacy naming conventions.
        
        Provides backward compatibility for code that expects:
        - output.tables["main"]
        - output.get_legacy_table("main")
        
        Args:
            name: Legacy table name (default "main")
            
        Returns:
            Path to table file, or empty string if not found
        """
        # Check if legacy "main" alias exists in _data
        if name == "main" and "main" in self._data:
            return self._data["main"]
        
        # Check in tables container
        if hasattr(self.tables, name):
            table_info = getattr(self.tables, name)
            return table_info.info.path if hasattr(table_info, 'info') else str(table_info)

        # For "main", try common table types in order of preference
        if name == "main":
            for fallback_name in ['sequences', 'structures', 'compounds']:
                if hasattr(self.tables, fallback_name):
                    table_info = getattr(self.tables, fallback_name)
                    return table_info.info.path if hasattr(table_info, 'info') else str(table_info)
        
        return ""
    
    def pretty(self) -> str:
        """Pretty formatted representation of the output."""
        from .datastream import DataStream
        import pandas as pd

        lines = []

        # Helper function to make paths relative to output_folder
        def make_relative_path(file_path: str) -> str:
            if self.output_folder and file_path.startswith(self.output_folder):
                relative_path = os.path.relpath(file_path, self.output_folder)
                return f"<output_folder>/{relative_path}"
            return file_path

        def format_datastream(stream_key: str, ds: DataStream) -> List[str]:
            """Format a DataStream for display."""
            import pandas as pd

            if not ds or len(ds) == 0:
                return []

            result = [f"{stream_key}:"]

            # Metadata fields
            result.append(f"    name: {ds.name}")
            result.append(f"    format: {ds.format}")
            result.append(f"    items: {len(ds)}")
            result.append(f"    map_table: '{make_relative_path(ds.map_table) if ds.map_table else ''}'")
            if ds.has_patterns():
                result.append(f"    has_patterns: True")
            if ds.metadata:
                meta_str = ", ".join(f"{k}={v}" for k, v in ds.metadata.items())
                result.append(f"    metadata: {meta_str}")

            # Data: try map_table first, then id/file fallback
            map_data = ds._get_map_data()
            if map_data is not None and len(map_data) > 0:
                columns = list(map_data.columns)
                n_rows = len(map_data)
                if n_rows <= 4:
                    row_indices = list(range(n_rows))
                    ellipsis_after = None
                else:
                    row_indices = list(range(2)) + list(range(n_rows - 2, n_rows))
                    ellipsis_after = 2

                for row_i, idx in enumerate(row_indices):
                    if ellipsis_after is not None and row_i == ellipsis_after:
                        result.append(f"    – ... ({n_rows - 4} more) ...")
                    row = map_data.iloc[idx]
                    parts = []
                    for col in columns:
                        val = str(row[col]) if pd.notna(row[col]) else ''
                        if len(val) > 40:
                            val = val[:37] + '...'
                        parts.append(f"{col}={val}")
                    result.append(f"    – {', '.join(parts)}")
            else:
                if ds.is_shared_file:
                    items = [(iid, ds.files) for iid in ds.ids]
                elif ds.files:
                    items = list(zip(ds.ids, ds.files))
                else:
                    items = [(iid, "") for iid in ds.ids]
                n_items = len(items)
                if n_items <= 4:
                    display_items = items
                    ellipsis_after = None
                else:
                    display_items = items[:2] + items[-2:]
                    ellipsis_after = 2

                for item_i, (item_id, item_file) in enumerate(display_items):
                    if ellipsis_after is not None and item_i == ellipsis_after:
                        result.append(f"    – ... ({n_items - 4} more) ...")
                    if item_file:
                        rel_path = make_relative_path(item_file)
                        result.append(f"    – {item_id}: '{rel_path}'")
                    else:
                        result.append(f"    – {item_id}")

            return result

        # Streams summary
        active_streams = [
            (k, ds) for k, ds in self.streams.items()
            if isinstance(ds, DataStream) and len(ds) > 0
        ]
        if active_streams:
            lines.append("streams:")
            for attr_name, ds in active_streams:
                mt = make_relative_path(ds.map_table) if ds.map_table else ''
                lines.append(f"    {attr_name}: format={ds.format}, items={len(ds)}, map_table='{mt}'")

        # Per-stream detail
        for attr_name, ds in active_streams:
            lines.extend(format_datastream(attr_name, ds))

        # Tables summary
        if 'tables' in self._data and hasattr(self.tables, '_tables') and self.tables._tables:
            lines.append("tables:")
            for name, info in self.tables._tables.items():
                col_display = ', '.join(info.info.columns[:]) if info.info.columns else ''
                lines.append(f"    {name} ({col_display}):")
                relative_path = make_relative_path(info.info.path)
                lines.append(f"        – '{relative_path}'")

            # Per-table detail
            for name, info in self.tables._tables.items():
                t_meta = info.info
                lines.append(f"{name}:")
                lines.append(f"    name: {t_meta.name}")
                lines.append(f"    path: '{make_relative_path(t_meta.path) if t_meta.path else ''}'")
                lines.append(f"    columns: {', '.join(t_meta.columns) if t_meta.columns else ''}")
                lines.append(f"    description: {t_meta.description}")
                if t_meta.path and os.path.exists(t_meta.path):
                    try:
                        table_df = pd.read_csv(t_meta.path)
                        if len(table_df) > 0:
                            columns = list(table_df.columns)
                            n_rows = len(table_df)
                            if n_rows <= 4:
                                row_indices = list(range(n_rows))
                                ellipsis_after = None
                            else:
                                row_indices = list(range(2)) + list(range(n_rows - 2, n_rows))
                                ellipsis_after = 2
                            for row_i, idx in enumerate(row_indices):
                                if ellipsis_after is not None and row_i == ellipsis_after:
                                    lines.append(f"    – ... ({n_rows - 4} more) ...")
                                row = table_df.iloc[idx]
                                parts = []
                                for col in columns:
                                    val = str(row[col]) if pd.notna(row[col]) else ''
                                    if len(val) > 40:
                                        val = val[:37] + '...'
                                    parts.append(f"{col}={val}")
                                lines.append(f"    – {', '.join(parts)}")
                    except Exception:
                        pass

        # Output folder
        if self.output_folder:
            lines.append("output_folder:")
            lines.append(f"    – '{self.output_folder}'")

        # Additional attributes (not DataStreams, tables, or output_folder)
        processed_keys = {'structures', 'sequences', 'compounds', 'msas',
                         'tables', 'output_folder', 'filter_metadata'}
        for key in sorted(self._data.keys()):
            if key in processed_keys:
                continue
            value = self._data[key]
            if isinstance(value, DataStream):
                continue  # Already handled
            if isinstance(value, str):
                lines.append(f"{key}:")
                lines.append(f"    – '{make_relative_path(value)}'")
            elif isinstance(value, list) and value:
                lines.append(f"{key}:")
                for item in value[:6]:
                    if isinstance(item, str):
                        lines.append(f"    – '{make_relative_path(item)}'")
                    else:
                        lines.append(f"    – {item}")
                if len(value) > 6:
                    lines.append(f"    – ... ({len(value) - 6} more)")

        return "\n".join(lines)
    
    def __str__(self) -> str:
        """String representation with improved formatting."""
        return self.pretty()
    
    def __repr__(self) -> str:
        """Detailed representation."""
        return f"StandardizedOutput({dict(self._data)})"

    @staticmethod
    def _is_colab() -> bool:
        """Check if currently running on Google Colab."""
        try:
            import google.colab  # noqa: F401
            return True
        except ImportError:
            return False

    def download(self):
        """
        Download the tool's output folder as a zip.
        On Colab, triggers a browser download. Elsewhere, prints the path.
        """
        if not self.output_folder or not os.path.isdir(self.output_folder):
            print("No output folder found.")
            return

        import zipfile
        folder_name = os.path.basename(self.output_folder)
        zip_path = os.path.join(os.path.dirname(self.output_folder), f"{folder_name}.zip")
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
            for root, _dirs, files in os.walk(self.output_folder):
                for fname in files:
                    fpath = os.path.join(root, fname)
                    arcname = os.path.relpath(fpath, self.output_folder)
                    zf.write(fpath, arcname)

        if self._is_colab():
            from google.colab import files as colab_files
            size_mb = os.path.getsize(zip_path) / (1024 * 1024)
            print(f"Downloading: {os.path.basename(zip_path)} ({size_mb:.1f} MB)")
            print("If no download starts, check that your browser allows popups from colab.research.google.com")
            colab_files.download(zip_path)
        else:
            print(f"Output zipped to: {zip_path}")

    @staticmethod
    def _is_notebook() -> bool:
        """Check if currently running in a Jupyter/Colab notebook."""
        try:
            from IPython import get_ipython
            shell = get_ipython()
            if shell is not None:
                if (shell.__class__.__name__ == "ZMQInteractiveShell"
                        or getattr(shell.__class__, "__module__", "") == "google.colab._shell"):
                    return True
        except (ImportError, NameError):
            pass
        return False

    # Cache for loaded renderer modules
    _renderer_cache = {}

    @staticmethod
    def _load_render_fn(script_path):
        """Load and cache a renderer's render() function from a script path."""
        if script_path not in StandardizedOutput._renderer_cache:
            import importlib.util
            spec = importlib.util.spec_from_file_location("renderer", script_path)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            StandardizedOutput._renderer_cache[script_path] = mod.render
        return StandardizedOutput._renderer_cache[script_path]

    def _get_renderers_config(self):
        """
        Load renderers config from config.yaml, resolving script paths
        relative to the repository root.
        """
        try:
            config = ConfigManager().get_renderers_config()
        except Exception:
            return {}
        if not config:
            return {}

        # Config names them as `renderers/<name>.py`, resolved against the package
        # so the same spelling works from a clone and from site-packages.
        package_dir = os.path.dirname(os.path.abspath(__file__))
        resolved = {}
        for section in ("streams", "tables"):
            if section in config:
                resolved[section] = {}
                for key, path in config[section].items():
                    resolved[section][key] = os.path.join(package_dir, path)
        return resolved

    _CSS = """<style>
.bp-table { border-collapse: collapse; font-family: monospace; font-size: 0.9em; margin: 4px 0 12px 0; }
.bp-table th { background: #f0f0f0; padding: 4px 10px; border: 1px solid #ddd; text-align: left; }
.bp-table td { padding: 4px 10px; border: 1px solid #ddd; }
.bp-table tr:nth-child(even) { background: #fafafa; }
.bp-table .bp-ellipsis td { color: #888; font-style: italic; text-align: center; }
.bp-section { margin-top: 8px; font-family: monospace; }
.bp-section-title { font-weight: bold; margin-bottom: 4px; }
.bp-table-toggle summary { cursor: pointer; font-family: monospace; color: #555; margin-bottom: 4px; }
.bp-table-toggle[open] + .bp-table-collapsed { display: none; }
</style>"""

    def _repr_html_(self) -> Optional[str]:
        """
        Rich HTML display for Jupyter notebooks.

        Dispatches rendering to scripts configured in config.yaml renderers section.
        Lookup order for streams: stream name -> format -> _default.
        Lookup order for tables: table name -> _default.
        """
        if not self._is_notebook():
            return None

        from .datastream import DataStream
        import html as html_module

        renderers_config = self._get_renderers_config()
        stream_renderers = renderers_config.get("streams", {})
        table_renderers = renderers_config.get("tables", {})
        metadata_script = stream_renderers.get("_metadata")

        html_parts = [self._CSS]

        def _rel(path):
            if self.output_folder and path.startswith(self.output_folder):
                return "<output_folder>/" + os.path.relpath(path, self.output_folder)
            return path

        # --- Tool-specific visualization (backwards compat) ---
        # Run first: if the tool provides its own rendering, skip stream dispatch.
        tool_html = ""
        if hasattr(self, "tool") and hasattr(self.tool, "_repr_notebook_html"):
            try:
                tool_html = self.tool._repr_notebook_html(self) or ""
            except Exception:
                pass

        # --- Streams summary ---
        active_streams = [
            (sn, s) for sn, s in self.streams.items()
            if isinstance(s, DataStream) and len(s) > 0
        ]

        if active_streams:
            html_parts.append(
                '<div class="bp-section">'
                '<div class="bp-section-title">STREAMS</div>'
            )
            html_parts.append(
                '<table class="bp-table">'
                '<tr><th>name</th><th>format</th><th>items</th><th>map_table</th></tr>'
            )
            for stream_name, stream in active_streams:
                html_parts.append(
                    f'<tr><td>{html_module.escape(stream_name)}</td>'
                    f'<td>{html_module.escape(stream.format)}</td>'
                    f'<td>{len(stream)}</td>'
                    f'<td>{html_module.escape(_rel(stream.map_table) if stream.map_table else "")}</td></tr>'
                )
            html_parts.append('</table></div>')

        # --- Render each stream via config ---
        if not tool_html:
            for stream_name, stream in self.streams.items():
                if not isinstance(stream, DataStream) or len(stream) == 0:
                    continue
                # Lookup: stream name first, then format, then _default
                script = (stream_renderers.get(stream_name)
                          or stream_renderers.get(stream.format)
                          or stream_renderers.get(stream.format.lower())
                          or stream_renderers.get("_default"))
                if not script:
                    continue
                try:
                    # If the matched renderer is not the metadata script itself,
                    # render metadata first, then the specialized renderer.
                    if metadata_script and script != metadata_script:
                        meta_fn = self._load_render_fn(metadata_script)
                        meta_result = meta_fn(stream, self)
                        if meta_result:
                            html_parts.append(meta_result)
                    render_fn = self._load_render_fn(script)
                    result = render_fn(stream, self)
                    if result:
                        html_parts.append(result)
                except Exception:
                    pass
        else:
            html_parts.append(tool_html)

        # --- Tables summary ---
        if "tables" in self._data and hasattr(self.tables, "_tables") and self.tables._tables:
            html_parts.append(
                '<div class="bp-section">'
                '<div class="bp-section-title">TABLES</div>'
            )
            html_parts.append(
                '<table class="bp-table">'
                '<tr><th>name</th><th>columns</th><th>path</th></tr>'
            )
            for name, info in self.tables._tables.items():
                cols = ", ".join(info.info.columns) if info.info.columns else ""
                path = _rel(info.info.path) if info.info.path else ""
                html_parts.append(
                    f"<tr><td>{html_module.escape(name)}</td>"
                    f"<td>{html_module.escape(cols)}</td>"
                    f"<td>{html_module.escape(path)}</td></tr>"
                )
            html_parts.append('</table></div>')

            # Render each table via config
            for name, info in self.tables._tables.items():
                script = (table_renderers.get(name)
                          or table_renderers.get("_default"))
                if script:
                    try:
                        render_fn = self._load_render_fn(script)
                        result = render_fn(info, self)
                        if result:
                            html_parts.append(result)
                    except Exception:
                        pass

        if len(html_parts) <= 1:  # only CSS
            return None
        return "\n".join(html_parts)

    def get_filter_info(self) -> Dict[str, Any]:
        """
        Get filtering information if this is filtered output.
        
        Returns:
            Dictionary with filter information or empty dict if not filtered
        """
        return self.filter_metadata.copy() if self.filter_metadata else {}
    
    def get_kept_items_count(self) -> int:
        """
        Get count of items that passed filtering.

        Returns:
            Number of kept items (structures + sequences + compounds)
        """
        count = 0
        for attr in (self.streams.structures, self.streams.sequences, self.streams.compounds):
            if attr is not None:
                count += len(attr)
        return count
    
    def get_original_items_count(self) -> Optional[int]:
        """
        Get count of original items before filtering (if available).
        
        Returns:
            Original item count or None if not available
        """
        if self.is_filtered and 'input_count' in self.filter_metadata:
            return self.filter_metadata['input_count']
        return None
    
    def get_filter_pass_rate(self) -> Optional[float]:
        """
        Calculate filter pass rate if filtering information is available.

        Returns:
            Pass rate (0.0 to 1.0) or None if not applicable
        """
        original_count = self.get_original_items_count()
        if original_count is not None and original_count > 0:
            kept_count = self.get_kept_items_count()
            return kept_count / original_count
        return None


class ToolOutput:
    """
    Container for tool output information.

    Returned by pipeline.add() to provide rich metadata about tool outputs
    and enable flexible chaining between tools.
    """
    
    def __init__(self, config: 'BaseConfig'):
        """Initialize with reference to tool configuration."""
        self.config = config
        self.tool_type = config.TOOL_NAME
        self.environments = config.environments
        self.output_folder = config.output_folder
        self.execution_order = config.execution_order
        
        # Will be populated after tool configuration
        self._output_files = {}
        self._metadata = {}
    
    def update_outputs(self, output_files: Dict[str, List[str]], metadata: Dict[str, Any] = None):
        """Update output files and metadata after tool configuration."""
        self._output_files = output_files
        self._metadata = metadata or {}
    
    def get_output_files(self, output_type: str = None) -> Union[Dict[str, List[str]], List[str]]:
        """
        Get output files, optionally filtered by type.
        
        Args:
            output_type: Specific output type to retrieve (e.g., 'pdbs', 'sequences')
            
        Returns:
            All outputs if output_type is None, otherwise specific output list
        """
        # If no cached outputs, delegate to config
        if not self._output_files and hasattr(self.config, 'get_output_files'):
            config_outputs = self.config.get_output_files()
            if output_type is None:
                return config_outputs
            return config_outputs.get(output_type, [])

        # Use cached outputs
        if output_type is None:
            return self._output_files
        return self._output_files.get(output_type, [])
    
    @property
    def output_pdbs(self) -> List[str]:
        """Convenience property for PDB outputs."""
        return self.get_output_files('pdbs')
    
    @property
    def output_sequences(self) -> List[str]:
        """Convenience property for sequence outputs."""
        return self.get_output_files('sequences')
    
    @property
    def output_structures(self) -> List[str]:
        """Convenience property for structure outputs (PDbs)."""
        return self.output_pdbs
    
    @property
    def output_tables(self) -> List[str]:
        """Convenience property for table outputs (JSON, CSV, etc.)."""
        return self.get_output_files('tables')

    @property
    def tables(self):
        """
        Access tool's tables for IDE autocompletion and column access.

        Enables usage like: tool_output.tables.structures.fixed
        Returns the tool's tables container with full column access.
        """
        if hasattr(self.config, 'tables'):
            return self.config.tables

        # Create TableContainer from get_output_files if available
        if hasattr(self.config, 'get_output_files'):
            output_files = self.config.get_output_files()
            tables_dict = output_files.get('tables', {})
            if tables_dict:
                return TableContainer(tables_dict)

        # Return empty container if nothing available
        return TableContainer({})
    
    @property
    def job_name(self) -> str:
        """Get job name from configuration."""
        return self.config.job_name
    
    @property
    def dependencies(self) -> List:
        """Get dependencies from configuration.""" 
        return self.config.dependencies
    
    @property
    def output(self) -> StandardizedOutput:
        """
        Get standardized output with dot notation access.

        Allows usage like:
        - rfd.streams.structures
        - rfd.tables
        - rfd['streams'] (dict-style)
        """
        # Get current output files from the config
        if hasattr(self.config, 'get_output_files'):
            output_files = self.config.get_output_files()
        else:
            output_files = self._output_files

        standardized_output = StandardizedOutput(output_files)
        # Add tool reference for user access
        standardized_output.tool = self.config
        standardized_output._stamp_producer(self.config)

        return standardized_output

    @property
    def o(self) -> StandardizedOutput:
        """
        Shorthand for .output - elegant access to tool outputs.

        Allows usage like:
        - tool.o.streams.structures
        - tool.o.tables.sequences
        - tool.o.tables.concatenated
        """
        return self

    def _repr_html_(self):
        """Delegate rich HTML display to StandardizedOutput."""
        return self.output._repr_html_()

    def __str__(self) -> str:
        return f"ToolOutput({self.tool_type}, {len(self._output_files)} output types)"

    def __repr__(self) -> str:
        return self.__str__()

