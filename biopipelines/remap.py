# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
ReMap tool for renaming IDs across all streams and tables from a source tool output.

Provides a lightweight way to rename IDs without re-running any computation.
At execution time, files are symlinked (or copied on Windows) and CSV tables
are rewritten with the new IDs.

Usage:
    # String basename: appends _1, _2, ...
    remapped = ReMap(source=tool_a, onto="design")

    # List of new IDs (matched per-stream where lengths agree)
    remapped = ReMap(source=tool_a, onto=["kinase_apo", "kinase_holo"])

    # Dict mapping old->new (for selective remapping)
    remapped = ReMap(source=tool_a, onto={"prot1_lig1": "complex_A"})

    # Align onto another tool's IDs (per-stream zip where lengths match,
    # or provenance-based matching when they don't)
    remapped = ReMap(source=tool_a, onto=tool_b)
    remapped = ReMap(source=tool_a, onto=tool_b.streams.structures)

    # Use an intermediate tool as bridge for provenance-based mapping
    remapped = ReMap(source=tool_a, onto=tool_c, map=tool_b)
"""

import os
import json
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream, create_map_table
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream, create_map_table


class ReMap(BaseConfig):
    """
    ReMap tool for renaming IDs across all streams and tables.

    Takes a source tool output and renames its IDs according to the provided
    mapping specification. At execution time, creates symlinks for files and
    rewrites CSV tables with remapped IDs.
    """

    TOOL_NAME = "ReMap"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== ReMap ==="
echo "Uses biopipelines environment (no additional installation needed)."
echo "=== ReMap ready ==="
"""

    # Lazy path descriptors
    remap_config_json = Path(lambda self: os.path.join(self.output_folder, "remap_config.json"))
    remap_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_remap.py"))

    def __init__(self,
                 source: StandardizedOutput,
                 onto: Union[str, List[str], Dict[str, str], StandardizedOutput, DataStream],
                 map: Optional[Union[StandardizedOutput, DataStream]] = None,
                 **kwargs):
        """
        Initialize ReMap configuration.

        Args:
            source: StandardizedOutput from a previous tool whose IDs will be renamed.
            onto: Target ID specification:
                - str: basename for auto-numbered IDs (e.g., "design" -> design_1, design_2, ...)
                - list: explicit new IDs (matched per-stream where lengths agree)
                - dict: old_id -> new_id mapping (selective remapping)
                - StandardizedOutput: align source IDs onto this tool's IDs (per-stream zip)
                - DataStream: align source IDs onto this stream's IDs
            map: Optional intermediate tool whose provenance bridges source and onto IDs.
                For example, if source=tool_a and onto=tool_c, passing map=tool_b
                uses tool_b's provenance columns to build the mapping from tool_a
                IDs to tool_c IDs.
            **kwargs: Additional parameters passed to BaseConfig
        """
        self.source = source
        self.onto_spec = onto
        self.map_source = map

        # Resolve onto IDs
        self.onto_ids = self._resolve_onto_ids()

        # Collect source stream info
        self.source_streams = self._collect_source_streams()

        # Compute the per-stream old->new mappings
        self.id_mapping = self._compute_id_mapping()

        # Collect source table info
        self.source_tables = self._collect_source_tables()

        super().__init__(**kwargs)

    def _resolve_onto_ids(self) -> Union[str, List[str], Dict[str, str], Dict[str, List[str]]]:
        """
        Resolve the onto specification into a usable form.

        Returns:
            - str: basename for auto-numbering
            - list: explicit ID list
            - dict (str->str): explicit old->new mapping
            - dict (str->list): per-stream-name -> IDs (from StandardizedOutput/DataStream)
        """
        onto = self.onto_spec

        if isinstance(onto, str):
            return onto

        if isinstance(onto, dict):
            return onto

        if isinstance(onto, list):
            return onto

        if isinstance(onto, DataStream):
            # Single stream: return its IDs keyed by stream name
            return {onto.name: list(onto.ids)}

        if isinstance(onto, StandardizedOutput):
            # Collect IDs from all non-empty streams
            stream_ids = {}
            for stream_name, stream in onto.streams.items():
                if isinstance(stream, DataStream) and len(stream) > 0:
                    stream_ids[stream_name] = list(stream.ids)
            return stream_ids

        raise ValueError(f"onto must be str, list, dict, DataStream, or StandardizedOutput, got {type(onto)}")

    def _collect_source_streams(self) -> List[Dict[str, Any]]:
        """Collect info about source streams for config generation."""
        streams_info = []
        for stream_name, stream in self.source.streams.items():
            if isinstance(stream, DataStream) and len(stream) > 0:
                info = {
                    "name": stream.name,
                    "stream_key": stream_name,
                    "format": stream.format,
                    "ids": list(stream.ids),
                    "files": list(stream.files),
                    "map_table": stream.map_table,
                    "files_contain_wildcards": stream.files_contain_wildcards,
                }
                streams_info.append(info)
        return streams_info

    def _compute_id_mapping(self) -> Dict[str, str]:
        """
        Compute a single old->new ID mapping that covers all streams.

        Strategy depends on the resolved onto_ids type:
        - str: auto-number based on each stream's IDs
        - list: match against streams with the same length
        - dict (str->str): direct mapping
        - dict (str->list): per-stream matching by name or length
        """
        onto = self.onto_ids
        mapping = {}

        if isinstance(onto, str):
            # String basename: generate numbered IDs from all source streams
            all_source_ids = []
            seen = set()
            for stream_info in self.source_streams:
                for sid in stream_info["ids"]:
                    if sid not in seen:
                        all_source_ids.append(sid)
                        seen.add(sid)
            for i, old_id in enumerate(all_source_ids, start=1):
                mapping[old_id] = f"{onto}_{i}"
            return mapping

        if isinstance(onto, list):
            # List: find streams where length matches, zip them
            for stream_info in self.source_streams:
                source_ids = stream_info["ids"]
                if len(source_ids) == len(onto):
                    for old_id, new_id in zip(source_ids, onto):
                        mapping[old_id] = new_id
            if not mapping:
                raise ValueError(
                    f"onto list has {len(onto)} IDs but no source stream has "
                    f"matching length. Stream lengths: "
                    f"{[(s['stream_key'], len(s['ids'])) for s in self.source_streams]}"
                )
            return mapping

        if isinstance(onto, dict):
            # Could be str->str (direct mapping) or str->list (per-stream)
            if not onto:
                return mapping

            first_val = next(iter(onto.values()))

            if isinstance(first_val, str):
                # Direct old->new mapping
                for stream_info in self.source_streams:
                    for old_id in stream_info["ids"]:
                        if old_id in onto:
                            mapping[old_id] = onto[old_id]
                return mapping

            if isinstance(first_val, list):
                # Per-stream mapping from StandardizedOutput/DataStream
                if self.map_source is not None:
                    return self._compute_mapping_via_bridge(onto)
                return self._compute_mapping_via_zip(onto)

        raise ValueError(f"Unexpected onto_ids type: {type(onto)}")

    def _resolve_via_provenance(self, source_ids: List[str], onto_stream: DataStream) -> Dict[str, str]:
        """
        Build mapping using provenance columns from the onto stream's map_table.

        Map_tables include provenance columns ({stream_name}.id) that track which
        input ID produced each output. This method reads those columns to find
        which source IDs map to which onto IDs.

        Args:
            source_ids: IDs from the source tool to map from
            onto_stream: DataStream whose map_table contains provenance columns

        Returns:
            Dict mapping source_id -> onto_id for matched IDs
        """
        mapping = {}
        if not onto_stream.map_table:
            return mapping

        map_data = onto_stream._get_map_data()
        if map_data is None or 'id' not in map_data.columns:
            return mapping

        source_id_set = set(str(s) for s in source_ids)
        provenance_cols = [c for c in map_data.columns if c.endswith('.id') and c != 'id']

        for prov_col in provenance_cols:
            prov_values = set(map_data[prov_col].astype(str))
            if not prov_values.intersection(source_id_set):
                continue

            # This provenance column links source IDs to onto IDs
            for _, row in map_data.iterrows():
                src_id = str(row[prov_col])
                onto_id = str(row['id'])
                if src_id in source_id_set and src_id not in mapping:
                    mapping[src_id] = onto_id

            if mapping:
                return mapping

        return mapping

    def _get_onto_streams(self) -> List[DataStream]:
        """Get the original DataStream objects from the onto specification."""
        onto = self.onto_spec
        if isinstance(onto, DataStream):
            return [onto]
        if isinstance(onto, StandardizedOutput):
            streams = []
            for stream_name, stream in onto.streams.items():
                if isinstance(stream, DataStream) and len(stream) > 0:
                    streams.append(stream)
            return streams
        return []

    def _compute_mapping_via_zip(self, onto_stream_ids: Dict[str, List[str]]) -> Dict[str, str]:
        """
        Build mapping by zipping source and onto IDs where stream lengths match.

        For each source stream, looks for an onto stream with matching length
        (first by name, then by length). When neither matches, falls back to
        provenance-based resolution using map_table provenance columns.
        """
        mapping = {}
        onto_by_length = {}
        for name, ids in onto_stream_ids.items():
            onto_by_length.setdefault(len(ids), []).append((name, ids))

        unmapped_source_streams = []

        for stream_info in self.source_streams:
            source_ids = stream_info["ids"]
            stream_key = stream_info["stream_key"]
            n = len(source_ids)

            # Try matching by stream name first
            if stream_key in onto_stream_ids and len(onto_stream_ids[stream_key]) == n:
                for old_id, new_id in zip(source_ids, onto_stream_ids[stream_key]):
                    mapping[old_id] = new_id
                continue

            # Fall back to matching by length
            if n in onto_by_length and onto_by_length[n]:
                _, onto_ids = onto_by_length[n][0]
                for old_id, new_id in zip(source_ids, onto_ids):
                    mapping[old_id] = new_id
                continue

            unmapped_source_streams.append(stream_info)

        # Provenance-based fallback for streams that couldn't be matched by name/length
        if unmapped_source_streams:
            onto_streams = self._get_onto_streams()
            for stream_info in unmapped_source_streams:
                source_ids = stream_info["ids"]
                for onto_stream in onto_streams:
                    prov_mapping = self._resolve_via_provenance(source_ids, onto_stream)
                    if prov_mapping:
                        mapping.update(prov_mapping)
                        break

        return mapping

    def _compute_mapping_via_bridge(self, onto_stream_ids: Dict[str, List[str]]) -> Dict[str, str]:
        """
        Build mapping using an intermediate tool's provenance to bridge source and onto IDs.

        The bridge tool (self.map_source) connects source IDs to onto IDs through its
        provenance columns. For example, if source=tool_a and onto=tool_c,
        bridge=tool_b:
        - tool_b's map_table has provenance columns (e.g., sequences.id) linking
          tool_a input IDs to tool_b output IDs
        - tool_b output IDs match tool_c's IDs
        - So: tool_a_id -> tool_b_id -> tool_c_id

        The bridge's map_tables contain provenance columns ({stream}.id) that link
        input IDs to output IDs.
        """
        mapping = {}
        bridge = self.map_source

        # Collect all onto IDs into a flat set for matching
        all_onto_ids = set()
        for ids in onto_stream_ids.values():
            all_onto_ids.update(ids)

        # Collect all source IDs for provenance matching
        all_source_ids = set()
        for stream_info in self.source_streams:
            all_source_ids.update(stream_info["ids"])

        # Look through bridge's streams for provenance that connects source to onto
        for stream_name, stream in bridge.streams.items():
            if not isinstance(stream, DataStream) or not stream.map_table:
                continue

            bridge_ids = list(stream.ids)

            # Check if bridge output IDs overlap with onto IDs
            if not all_onto_ids.intersection(bridge_ids):
                continue

            # This bridge stream's IDs match onto IDs.
            # Read its map_table to find provenance columns linking back to source IDs.
            map_data = stream._get_map_data()
            if map_data is None or 'id' not in map_data.columns:
                continue

            provenance_cols = [c for c in map_data.columns if c.endswith('.id') and c != 'id']

            for prov_col in provenance_cols:
                prov_values = set(map_data[prov_col].astype(str))
                if not prov_values.intersection(all_source_ids):
                    continue

                # This provenance column links source IDs to bridge IDs
                for _, row in map_data.iterrows():
                    src_id = str(row[prov_col])
                    bridge_id = str(row['id'])
                    if src_id in all_source_ids and bridge_id in all_onto_ids:
                        mapping[src_id] = bridge_id

                if mapping:
                    return mapping

        # Fallback: try zip approach if bridge mapping didn't work
        if not mapping:
            return self._compute_mapping_via_zip(onto_stream_ids)

        return mapping

    def _collect_source_tables(self) -> List[Dict[str, str]]:
        """Collect info about source tables."""
        tables_info = []
        if hasattr(self.source, 'tables') and hasattr(self.source.tables, '_tables'):
            for name, table_info in self.source.tables._tables.items():
                tables_info.append({
                    "name": name,
                    "path": table_info.info.path,
                    "columns": table_info.info.columns,
                })
        return tables_info

    def validate_params(self):
        """Validate ReMap parameters."""
        pass

    def configure_inputs(self, pipeline_folders):
        """Configure inputs from pipeline context."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()
        config_lines.append(f"MAPPINGS: {len(self.id_mapping)} IDs")
        for old_id, new_id in list(self.id_mapping.items())[:5]:
            config_lines.append(f"  {old_id} -> {new_id}")
        if len(self.id_mapping) > 5:
            config_lines.append(f"  ... and {len(self.id_mapping) - 5} more")
        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate bash script for ReMap execution."""
        config = {
            "id_mapping": self.id_mapping,
            "output_folder": self.output_folder,
            "streams": self.source_streams,
            "tables": self.source_tables,
        }

        os.makedirs(self.output_folder, exist_ok=True)
        with open(self.remap_config_json, 'w') as f:
            json.dump(config, f, indent=2)

        script = "#!/bin/bash\n"
        script += f"# ReMap: Rename IDs in streams and tables\n\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()

        script += f'echo "=== ReMap: Renaming IDs ==="\n'
        script += f'python "{self.remap_py}" "{self.remap_config_json}"\n'
        script += f'echo "=== ReMap complete ==="\n\n'

        script += self.generate_completion_check_footer()
        return script

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files with remapped IDs."""
        output = {}

        for stream_info in self.source_streams:
            stream_key = stream_info["stream_key"]
            fmt = stream_info["format"]
            old_ids = stream_info["ids"]
            old_files = stream_info["files"]

            new_ids = [self.id_mapping.get(oid, oid) for oid in old_ids]

            # Build new file paths
            new_files = []
            if old_files and len(old_files) == len(old_ids):
                # One file per ID: create new paths
                for new_id, old_file in zip(new_ids, old_files):
                    ext = os.path.splitext(old_file)[1]
                    new_file = os.path.join(self.output_folder, f"{new_id}{ext}")
                    new_files.append(new_file)
            elif old_files and len(old_files) == 1:
                # Single shared file (e.g., sequences CSV): will be rewritten
                new_files = [os.path.join(self.output_folder,
                             os.path.basename(old_files[0]))]
            else:
                new_files = []

            # Create map_table for this stream
            map_table_path = os.path.join(self.output_folder, f"{stream_info['name']}_map.csv")

            # Build provenance: trace new IDs back to old IDs
            provenance = {stream_key: old_ids}

            create_map_table(
                map_table_path,
                ids=new_ids,
                files=new_files if new_files else None,
                provenance=provenance,
            )

            output[stream_key] = DataStream(
                name=stream_info["name"],
                ids=new_ids,
                files=new_files,
                map_table=map_table_path,
                format=fmt,
                files_contain_wildcards=stream_info.get("files_contain_wildcards", False),
            )

        # Build remapped tables
        tables = {}
        for table_info in self.source_tables:
            table_name = table_info["name"]
            new_table_path = os.path.join(self.output_folder, f"{table_name}.csv")
            tables[table_name] = {
                "path": new_table_path,
                "columns": table_info.get("columns", []),
                "description": f"Remapped {table_name} table",
            }

        output["tables"] = tables
        output["output_folder"] = self.output_folder

        return output
