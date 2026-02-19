# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
ReMap tool for renaming IDs across all streams and tables from a source tool output.

Provides a lightweight way to rename auto-generated IDs to user-friendly names
without re-running any computation. At execution time, files are symlinked
(or copied on Windows) and CSV tables are rewritten with the new IDs.

Usage:
    # String basename: appends _1, _2, ...
    remapped = ReMap(source=boltz, ids="design")

    # List of new IDs (must match length)
    remapped = ReMap(source=boltz, ids=["kinase_apo", "kinase_holo"])

    # Dict mapping old->new (for selective remapping)
    remapped = ReMap(source=boltz, ids={"prot1_lig1": "complex_A", "prot1_lig2": "complex_B"})
"""

import os
import json
from typing import Dict, List, Any, Union

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
                 ids: Union[str, List[str], Dict[str, str]],
                 **kwargs):
        """
        Initialize ReMap configuration.

        Args:
            source: StandardizedOutput from a previous tool
            ids: Remapping specification:
                - str: basename for auto-numbered IDs (e.g., "design" -> design_1, design_2, ...)
                - list: explicit new IDs (must match count of source IDs)
                - dict: old_id -> new_id mapping (all keys must exist in source)
            **kwargs: Additional parameters passed to BaseConfig
        """
        self.source = source
        self.ids_spec = ids

        # Resolve source IDs from the first non-empty stream
        self.source_ids = self._get_source_ids()
        if not self.source_ids:
            raise ValueError("ReMap source has no IDs in any stream")

        # Compute the old->new mapping
        self.id_mapping = self._compute_id_mapping()

        # Collect source stream info for config generation
        self.source_streams = self._collect_source_streams()

        # Collect source table info
        self.source_tables = self._collect_source_tables()

        # Collect source map_table paths
        self.source_map_tables = self._collect_source_map_tables()

        super().__init__(**kwargs)

    def _get_source_ids(self) -> List[str]:
        """Get IDs from the first non-empty stream in source."""
        for stream_name, stream in self.source.streams.items():
            if isinstance(stream, DataStream) and len(stream) > 0:
                return list(stream.ids)
        return []

    def _compute_id_mapping(self) -> Dict[str, str]:
        """Compute old->new ID mapping from the ids specification."""
        if isinstance(self.ids_spec, str):
            # String basename: generate numbered IDs
            mapping = {}
            for i, old_id in enumerate(self.source_ids, start=1):
                mapping[old_id] = f"{self.ids_spec}_{i}"
            return mapping

        elif isinstance(self.ids_spec, list):
            # List of new IDs: must match length
            if len(self.ids_spec) != len(self.source_ids):
                raise ValueError(
                    f"ids list length ({len(self.ids_spec)}) does not match "
                    f"source ID count ({len(self.source_ids)})"
                )
            return dict(zip(self.source_ids, self.ids_spec))

        elif isinstance(self.ids_spec, dict):
            # Dict mapping: validate all keys exist
            missing = set(self.ids_spec.keys()) - set(self.source_ids)
            if missing:
                raise ValueError(f"ids dict contains unknown source IDs: {missing}")
            # For IDs not in the dict, keep the original name
            mapping = {}
            for old_id in self.source_ids:
                mapping[old_id] = self.ids_spec.get(old_id, old_id)
            return mapping

        else:
            raise ValueError(f"ids must be str, list, or dict, got {type(self.ids_spec)}")

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

    def _collect_source_map_tables(self) -> List[str]:
        """Collect all map_table paths from source streams."""
        paths = []
        for stream_info in self.source_streams:
            if stream_info["map_table"]:
                paths.append(stream_info["map_table"])
        return paths

    def validate_params(self):
        """Validate ReMap parameters."""
        pass

    def configure_inputs(self, pipeline_folders):
        """Configure inputs from pipeline context."""
        self.folders = pipeline_folders

    def generate_script(self, script_path: str) -> str:
        """Generate bash script for ReMap execution."""
        # Write the remap config JSON at configuration time
        config = {
            "id_mapping": self.id_mapping,
            "output_folder": self.output_folder,
            "streams": self.source_streams,
            "tables": self.source_tables,
            "map_tables": self.source_map_tables,
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

        # Build remapped streams
        for stream_info in self.source_streams:
            stream_key = stream_info["stream_key"]
            fmt = stream_info["format"]
            old_ids = stream_info["ids"]
            old_files = stream_info["files"]

            new_ids = [self.id_mapping.get(oid, oid) for oid in old_ids]

            # Build new file paths (symlinks will be created at execution time)
            new_files = []
            if old_files:
                for old_id, new_id, old_file in zip(old_ids, new_ids, old_files):
                    ext = os.path.splitext(old_file)[1]
                    new_file = os.path.join(self.output_folder, f"{new_id}{ext}")
                    new_files.append(new_file)

            # Create map_table for this stream
            map_table_path = os.path.join(self.output_folder, f"{stream_info['name']}_map.csv")

            # Build provenance: trace new IDs back to old IDs via stream_key
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
