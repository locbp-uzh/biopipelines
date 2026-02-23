# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Load tool for loading previously saved pipeline tool outputs.

Allows reusing results from previous pipeline runs without re-execution,
enabling incremental pipeline development and efficient result reuse.
"""

import os
import re
import json
from typing import Dict, List, Any, Optional, Union
from datetime import datetime

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream


class Load(BaseConfig):
    """
    Load results from a previously executed pipeline tool, with optional filtering.

    This tool allows reusing outputs from previous pipeline runs without
    re-execution, enabling incremental pipeline development and result reuse.

    The Load tool:
    - Loads .expected_outputs.json from a tool's output folder
    - Rebases file paths when the folder has been moved
    - Optionally applies filters to loaded data (expression string or Filter tool output)
    - Provides filtered results through the standard pipeline interface
    - Enables subsequent tools to work with the filtered subset

    This is the ONLY tool that checks file existence at configuration time (not execution time).
    """

    # Tool identification
    TOOL_NAME = "Load"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== Load ==="
echo "Uses biopipelines environment (no additional installation needed)."
echo "=== Load ready ==="
"""

    def __init__(self, path: str, filter = None, validate_files: bool = True, **kwargs):
        """
        Initialize Load tool.

        Args:
            path: Path to the tool's output folder (containing .expected_outputs.json)
            filter: Optional filter result (ToolOutput from Filter tool) to only load specific IDs
            validate_files: Whether to validate that all referenced files exist
            **kwargs: Additional parameters for BaseConfig

        Output:
            Streams: inherits from the loaded tool (structures, sequences, compounds, etc.)
            Tables: inherits all tables from the loaded tool
        """
        self.tool_folder = os.path.abspath(path)
        self.result_file = os.path.join(self.tool_folder, ".expected_outputs.json")
        self.filter_input = filter
        self.validate_files = validate_files
        self.loaded_result = None
        self.missing_files = []
        self.original_tool_name = None
        self.filtered_ids = None

        # Load and validate the result file
        self._load_and_validate_result()

        # Process filter if provided
        if self.filter_input:
            self._process_filter()

        # Set up job name based on loaded result
        if not kwargs.get('job_name'):
            original_job_name = self.loaded_result.get('job_name')

            if not original_job_name or original_job_name == 'unknown':
                if ('configuration' in self.loaded_result and
                    'pipeline_context' in self.loaded_result['configuration']):
                    pipeline_context = self.loaded_result['configuration']['pipeline_context']
                    original_job_name = pipeline_context.get('pipeline_job_name', 'unknown')

            if not original_job_name:
                original_job_name = 'unknown'

            kwargs['job_name'] = f"load_{original_job_name}"

        # Initialize base class
        super().__init__(**kwargs)

        # Set up dependency on filter input if provided
        if self.filter_input and hasattr(self.filter_input, 'config'):
            self.dependencies.append(self.filter_input.config)

        # Initialize folders to prevent AttributeError during script generation
        self.folders = {"HelpScripts": "HelpScripts"}  # Will be updated in configure_inputs

    def _load_and_validate_result(self):
        """Load and validate the result file, rebasing paths if the folder was moved."""
        if not os.path.exists(self.result_file):
            raise ValueError(f".expected_outputs.json not found in: {self.tool_folder}")

        try:
            with open(self.result_file, 'r') as f:
                self.loaded_result = json.load(f)
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON in {self.result_file}: {e}")

        # Validate required fields
        required_fields = ['tool_name', 'tool_class', 'output_structure']
        for field in required_fields:
            if field not in self.loaded_result:
                raise ValueError(f"Missing required field '{field}' in {self.result_file}")

        self.original_tool_name = self.loaded_result['tool_name']

        # Rebase paths if the folder has been moved
        output_structure = self.loaded_result['output_structure']
        old_output_folder = output_structure.get('output_folder', '')
        if old_output_folder and os.path.normpath(old_output_folder) != os.path.normpath(self.tool_folder):
            # The folder was moved — rebase all paths
            old_prefix = os.path.dirname(old_output_folder)  # parent of old tool folder
            new_prefix = os.path.dirname(self.tool_folder)    # parent of new tool folder
            print(f"Load: Rebasing paths from {old_prefix} -> {new_prefix}")
            self._rebase_paths(output_structure, old_prefix, new_prefix)

        # Validate file existence if requested
        if self.validate_files:
            self._validate_file_existence()

    def _rebase_paths(self, obj, old_prefix: str, new_prefix: str):
        """
        Recursively rebase all file paths in the output structure.

        Replaces old_prefix with new_prefix in all string values that
        start with old_prefix. Modifies the structure in place.
        """
        old_norm = os.path.normpath(old_prefix)

        if isinstance(obj, dict):
            for key in obj:
                if isinstance(obj[key], str) and os.path.normpath(obj[key]).startswith(old_norm):
                    obj[key] = obj[key].replace(old_prefix, new_prefix, 1)
                elif isinstance(obj[key], (dict, list)):
                    self._rebase_paths(obj[key], old_prefix, new_prefix)
        elif isinstance(obj, list):
            for i, item in enumerate(obj):
                if isinstance(item, str) and os.path.normpath(item).startswith(old_norm):
                    obj[i] = item.replace(old_prefix, new_prefix, 1)
                elif isinstance(item, (dict, list)):
                    self._rebase_paths(item, old_prefix, new_prefix)

    def _validate_file_existence(self):
        """Validate that all referenced files in the output structure exist."""
        self.missing_files = []
        output_structure = self.loaded_result['output_structure']

        # Check DataStream files
        for stream_name in ['structures', 'sequences', 'compounds']:
            if stream_name in output_structure and isinstance(output_structure[stream_name], dict):
                ds_dict = output_structure[stream_name]
                for file_path in ds_dict.get('files', []):
                    if isinstance(file_path, str) and not os.path.exists(file_path):
                        self.missing_files.append(file_path)

        # Check table files
        if 'tables' in output_structure:
            tables = output_structure['tables']
            if isinstance(tables, dict):
                for table_name, table_info in tables.items():
                    if isinstance(table_info, dict) and 'path' in table_info:
                        path = table_info['path']
                        if not os.path.exists(path):
                            self.missing_files.append(path)

        # Check output folder
        if 'output_folder' in output_structure:
            output_folder = output_structure['output_folder']
            if not os.path.exists(output_folder):
                self.missing_files.append(output_folder)

        # Report missing files but don't fail (files might be on a different machine)
        if self.missing_files:
            print(f"Warning: {len(self.missing_files)} files referenced in {self.result_file} are missing:")
            for missing_file in self.missing_files[:5]:
                print(f"  - {missing_file}")
            if len(self.missing_files) > 5:
                print(f"  ... and {len(self.missing_files) - 5} more")

    def _resolve_file_paths(self, output_structure: Dict[str, Any]) -> Dict[str, Any]:
        """
        Resolve glob patterns and match IDs to actual file paths.

        Called at configuration time when validate_files=True.

        Args:
            output_structure: The output structure to resolve

        Returns:
            Updated output structure with resolved file paths
        """
        import glob as glob_module

        resolved = output_structure.copy()

        for file_type in ['structures', 'sequences', 'compounds']:
            if file_type not in resolved or not isinstance(resolved[file_type], dict):
                continue

            ds_dict = resolved[file_type]
            files = ds_dict.get('files', [])
            ids = ds_dict.get('ids', [])

            if not files or not ids:
                continue

            # Already matched (same length, no globs)
            if len(files) == len(ids):
                has_glob = any('*' in f or '?' in f for f in files if isinstance(f, str))
                if not has_glob:
                    continue

            # Need to resolve
            print(f"Load: Resolving {file_type} paths ({len(files)} patterns -> {len(ids)} IDs)")

            # Collect all potential files
            all_files = []
            for file_pattern in files:
                if not isinstance(file_pattern, str):
                    continue

                if '*' in file_pattern or '?' in file_pattern:
                    expanded = glob_module.glob(file_pattern)
                    all_files.extend(expanded)
                    print(f"  Expanded glob: {len(expanded)} files")
                elif os.path.isfile(file_pattern):
                    all_files.append(file_pattern)
                elif os.path.isdir(file_pattern):
                    extensions = ['.pdb', '.cif', '.mmcif', '.fasta', '.fa', '.csv', '.sdf', '.mol2']
                    for ext in extensions:
                        all_files.extend(glob_module.glob(os.path.join(file_pattern, f"*{ext}")))

            # Also check output_folder if we don't have enough files
            if len(all_files) < len(ids) and 'output_folder' in resolved:
                output_folder = resolved['output_folder']
                if os.path.isdir(output_folder):
                    extensions = ['.pdb', '.cif', '.mmcif', '.fasta', '.fa', '.csv', '.sdf', '.mol2']
                    for ext in extensions:
                        for f in glob_module.glob(os.path.join(output_folder, f"**/*{ext}"), recursive=True):
                            if f not in all_files:
                                all_files.append(f)

            # Match IDs to files
            resolved_files = []
            resolved_ids = []

            for item_id in ids:
                matched_file = None

                for file_path in all_files:
                    basename = os.path.splitext(os.path.basename(file_path))[0]

                    if basename == item_id:
                        matched_file = file_path
                        break
                    if item_id in basename:
                        matched_file = file_path
                        break
                    if basename in item_id:
                        matched_file = file_path
                        break
                    if item_id.startswith("rank") and basename.startswith("rank"):
                        id_num = ''.join(filter(str.isdigit, item_id))
                        file_rank_part = basename.split('_')[0]
                        file_num = ''.join(filter(str.isdigit, file_rank_part))
                        if id_num and file_num and id_num == file_num:
                            matched_file = file_path
                            break

                if matched_file:
                    resolved_files.append(matched_file)
                    actual_id = os.path.splitext(os.path.basename(matched_file))[0]
                    resolved_ids.append(actual_id)
                else:
                    print(f"  Warning: No file found for ID '{item_id}'")

            resolved[file_type]['files'] = resolved_files
            resolved[file_type]['ids'] = resolved_ids
            resolved[file_type]['files_contain_wildcards'] = False

            print(f"  Resolved: {len(resolved_files)}/{len(ids)} {file_type}")

        return resolved

    def _process_filter(self):
        """Process filter input to determine which IDs to keep."""
        import pandas as pd

        # Find the main table to filter
        output_structure = self.loaded_result['output_structure']
        main_table_path = None

        if 'tables' in output_structure:
            tables = output_structure['tables']
            if isinstance(tables, dict):
                priority_names = ['structures', 'analysis', 'combined', 'results']
                for name in priority_names:
                    if name in tables:
                        ds_info = tables[name]
                        if isinstance(ds_info, dict) and 'path' in ds_info:
                            main_table_path = ds_info['path']
                            break

                if not main_table_path:
                    first_ds = next(iter(tables.values()))
                    if isinstance(first_ds, dict) and 'path' in first_ds:
                        main_table_path = first_ds['path']

        if not main_table_path or not os.path.exists(main_table_path):
            raise ValueError(f"Cannot find main table to filter in loaded output: {main_table_path}")

        print(f"Load: Applying filter to table: {main_table_path}")

        try:
            df = pd.read_csv(main_table_path)
            print(f"  - Loaded {len(df)} rows")
        except Exception as e:
            raise ValueError(f"Error loading table for filtering: {e}")

        # Apply filter
        if isinstance(self.filter_input, str):
            try:
                filtered_df = df.query(self.filter_input)
                print(f"  - Applied expression: {self.filter_input}")
            except Exception as e:
                raise ValueError(f"Error applying filter expression '{self.filter_input}': {e}")

        elif hasattr(self.filter_input, 'tables'):
            filter_tables = self.filter_input.tables
            if isinstance(filter_tables, dict):
                if 'filtered' in filter_tables:
                    filter_ds_info = filter_tables['filtered']
                else:
                    filter_ds_info = next(iter(filter_tables.values()))

                if isinstance(filter_ds_info, dict) and 'path' in filter_ds_info:
                    filter_csv_path = filter_ds_info['path']
                elif isinstance(filter_ds_info, str):
                    filter_csv_path = filter_ds_info
                else:
                    raise ValueError("Cannot determine filter CSV path")

                if not os.path.exists(filter_csv_path):
                    raise ValueError(f"Filter CSV file not found: {filter_csv_path}")

                try:
                    filter_df = pd.read_csv(filter_csv_path)
                    if 'id' in filter_df.columns and 'id' in df.columns:
                        filtered_df = df[df['id'].isin(filter_df['id'])]
                        print(f"  - Filtered using precomputed results: {len(filter_df)} IDs")
                    else:
                        raise ValueError("Cannot match filtered IDs - missing 'id' columns")
                except Exception as e:
                    raise ValueError(f"Error loading filter results from {filter_csv_path}: {e}")
            else:
                raise ValueError("Invalid filter input format")
        else:
            raise ValueError(f"Invalid filter input type: {type(self.filter_input)}")

        print(f"  - Result: {len(filtered_df)}/{len(df)} items kept")

        if 'id' in filtered_df.columns:
            self.filtered_ids = set(filtered_df['id'].tolist())
        else:
            self.filtered_ids = set(filtered_df.index.tolist())

        print(f"  - Filtered IDs: {sorted(list(self.filtered_ids)) if len(self.filtered_ids) <= 10 else f'{len(self.filtered_ids)} items'}")

    def _apply_filter_to_output_structure(self, output_structure: Dict[str, Any]) -> Dict[str, Any]:
        """Apply filtering to the output structure based on filtered IDs."""

        print(f"Load: Filtering output structure to {len(self.filtered_ids)} items")

        filtered_structure = output_structure.copy()

        # Filter DataStream dicts
        for stream_name in ['structures', 'sequences', 'compounds']:
            if stream_name not in output_structure or not isinstance(output_structure[stream_name], dict):
                continue

            ds_dict = output_structure[stream_name]
            ids = ds_dict.get('ids', [])
            files = ds_dict.get('files', [])

            filtered_ids = []
            filtered_files = []

            if len(files) == len(ids):
                for item_id, file_path in zip(ids, files):
                    if item_id in self.filtered_ids:
                        filtered_ids.append(item_id)
                        filtered_files.append(file_path)
            elif len(files) <= 1:
                # Single file or no files — just filter IDs
                filtered_ids = [i for i in ids if i in self.filtered_ids]
                # For single CSV files, update path to Load's output folder
                filtered_files = []
                for file_path in files:
                    if isinstance(file_path, str) and file_path.endswith('.csv'):
                        filtered_files.append(os.path.join(self.output_folder, os.path.basename(file_path)))
                    else:
                        filtered_files.append(file_path)

            filtered_structure[stream_name] = dict(ds_dict)
            filtered_structure[stream_name]['ids'] = filtered_ids
            filtered_structure[stream_name]['files'] = filtered_files
            print(f"  - Filtered {stream_name}: {len(filtered_ids)}/{len(ids)}")

        # Update tables paths to point to Load's output folder and update counts
        if 'tables' in filtered_structure:
            for ds_name, ds_info in filtered_structure['tables'].items():
                if isinstance(ds_info, dict):
                    if 'path' in ds_info:
                        original_path = ds_info['path']
                        if isinstance(original_path, str) and original_path.endswith('.csv'):
                            ds_info['path'] = os.path.join(self.output_folder, f"{ds_name}.csv")
                    if 'count' in ds_info and isinstance(ds_info['count'], int):
                        ds_info['count'] = len(self.filtered_ids)

        return filtered_structure

    def _exclude_missing_ids(self, output_structure: Dict[str, Any]) -> set:
        """
        Check for 'missing' table and return set of IDs to exclude.
        """
        import pandas as pd

        if 'tables' not in output_structure:
            return set()

        tables = output_structure['tables']
        if not isinstance(tables, dict) or 'missing' not in tables:
            return set()

        missing_info = tables['missing']
        if isinstance(missing_info, dict) and 'path' in missing_info:
            missing_path = missing_info['path']
        else:
            return set()

        if not os.path.exists(missing_path):
            raise FileNotFoundError(f"Missing table referenced but file not found: {missing_path}")

        missing_df = pd.read_csv(missing_path)
        if 'id' not in missing_df.columns:
            raise ValueError(f"missing.csv does not have required 'id' column: {missing_path}")

        missing_ids = set(missing_df['id'].tolist())
        print(f"Load: Excluding {len(missing_ids)} IDs from missing.csv")
        return missing_ids

    def validate_params(self):
        """Validate Load parameters."""
        if not self.loaded_result:
            raise ValueError("No result loaded - initialization failed")

        if self.validate_files and len(self.missing_files) > 0:
            print(f"Warning: {len(self.missing_files)} referenced files are missing")
            print("Consider setting validate_files=False if files have been moved")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure Load - no inputs needed since we're loading existing results."""
        self.folders = pipeline_folders

    def get_output_files(self) -> Dict[str, Any]:
        """
        Return the loaded output structure converted to DataStream format.

        Returns:
            Dict with DataStream objects for structures, sequences, compounds,
            plus tables dict and output_folder string.
        """
        if not self.loaded_result:
            raise RuntimeError("No result loaded")

        output_structure = self.loaded_result['output_structure'].copy()

        # Resolve glob patterns and match IDs to files when validate_files=True
        if self.validate_files:
            output_structure = self._resolve_file_paths(output_structure)

        # Apply filtering if filter was provided
        if self.filtered_ids is not None:
            output_structure = self._apply_filter_to_output_structure(output_structure)
        else:
            # Check for missing table and exclude those IDs
            missing_ids = self._exclude_missing_ids(output_structure)
            if missing_ids:
                self._filter_missing_ids(output_structure, missing_ids)

        return self._convert_to_datastream_format(output_structure)

    def _filter_missing_ids(self, output_structure: Dict[str, Any], missing_ids: set):
        """Filter out missing IDs from output structure in place."""
        for stream_name in ['structures', 'sequences', 'compounds']:
            if stream_name not in output_structure or not isinstance(output_structure[stream_name], dict):
                continue
            ds_dict = output_structure[stream_name]
            ids = ds_dict.get('ids', [])
            files = ds_dict.get('files', [])
            if len(files) == len(ids):
                filtered_files = []
                filtered_ids = []
                for f, i in zip(files, ids):
                    if i not in missing_ids:
                        filtered_files.append(f)
                        filtered_ids.append(i)
                ds_dict['files'] = filtered_files
                ds_dict['ids'] = filtered_ids
                if len(filtered_files) < len(files):
                    print(f"  - Filtered {stream_name}: {len(filtered_files)}/{len(files)} kept")

    def _convert_to_datastream_format(self, output_structure: Dict[str, Any]) -> Dict[str, Any]:
        """
        Convert output structure to DataStream format.

        Args:
            output_structure: Dict with DataStream dicts

        Returns:
            Dict with DataStream objects for structures, sequences, compounds
        """
        from .datastream import DataStream

        output_folder = output_structure.get('output_folder', '')

        def convert_data_type(data_key: str, default_format: str) -> DataStream:
            raw_data = output_structure.get(data_key)

            if isinstance(raw_data, dict) and 'ids' in raw_data and 'files' in raw_data:
                return DataStream.from_dict(raw_data)

            return None

        structures = convert_data_type('structures', 'pdb')
        sequences = convert_data_type('sequences', 'fasta')
        compounds = convert_data_type('compounds', 'sdf')

        # Convert tables to TableInfo objects
        tables_dict = {}
        raw_tables = output_structure.get('tables', {})
        if isinstance(raw_tables, dict):
            for table_name, table_info in raw_tables.items():
                if isinstance(table_info, dict) and 'path' in table_info:
                    tables_dict[table_name] = TableInfo(
                        name=table_name,
                        path=table_info['path'],
                        columns=table_info.get('columns', []),
                        description=table_info.get('description', ''),
                        count=table_info.get('count', 0)
                    )
                elif isinstance(table_info, TableInfo):
                    tables_dict[table_name] = table_info

        return {
            "structures": structures,
            "sequences": sequences,
            "compounds": compounds,
            "tables": tables_dict,
            "output_folder": output_folder
        }

    def generate_script(self, script_path: str) -> str:
        """
        Generate no-op script since files already exist.
        """
        original_tool = self.loaded_result.get('tool_name', 'Unknown')
        original_job = self.loaded_result.get('job_name')

        if not original_job or original_job == 'unknown':
            if ('configuration' in self.loaded_result and
                'pipeline_context' in self.loaded_result['configuration']):
                pipeline_context = self.loaded_result['configuration']['pipeline_context']
                original_job = pipeline_context.get('pipeline_job_name', 'unknown')

        if not original_job:
            original_job = 'unknown'

        script_content = "#!/bin/bash\n"
        script_content += f"# Load script - loading results from {original_tool}\n"
        script_content += f"# Original job: {original_job}\n"
        script_content += f"# Tool folder: {self.tool_folder}\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += f"""
echo "Loading output from previous {original_tool} execution"
echo "Original job: {original_job}"
echo "Tool folder: {self.tool_folder}"

# Validate that key files exist
echo "Validating loaded files..."

# Check output folder
if [ ! -d "{self.loaded_result['output_structure'].get('output_folder', '')}" ]; then
    echo "Warning: Output folder missing: {self.loaded_result['output_structure'].get('output_folder', '')}"
fi

# Check some key files (first few from each category)
"""

        output_structure = self.loaded_result['output_structure']

        for stream_name in ['structures', 'sequences', 'compounds']:
            if stream_name in output_structure and isinstance(output_structure[stream_name], dict):
                files_list = output_structure[stream_name].get('files', [])
                files_to_check = files_list[:3]
                for file_path in files_to_check:
                    if isinstance(file_path, str):
                        script_content += f"""
if [ ! -f "{file_path}" ]; then
    echo "Warning: {stream_name} file missing: {file_path}"
fi"""

        # Add filtering section if filter was applied
        if self.filter_input and self.filtered_ids:
            filter_config = {
                "filtered_ids": list(self.filtered_ids),
                "output_structure": self.loaded_result['output_structure'],
                "output_folder": self.output_folder
            }

            filter_config_path = os.path.join(self.output_folder, "filter_config.json")
            with open(filter_config_path, 'w') as f:
                json.dump(filter_config, f, indent=2)

            script_content += f"""

# Apply filter and create filtered tables
echo "Creating filtered copies of tables..."
echo "Filter: {self.filter_input if isinstance(self.filter_input, str) else 'ToolOutput filter'}"
echo "Filtered IDs: {len(self.filtered_ids)} items"
echo "Output folder: {self.output_folder}"

python "{os.path.join(self.folders.get('HelpScripts', 'HelpScripts'), 'pipe_load_output_filter.py')}" \\
  --config "{filter_config_path}"

if [ $? -ne 0 ]; then
    echo "Error: Failed to create filtered tables"
    exit 1
fi

echo "Filtered tables created successfully"
"""

        script_content += f"""

echo "Load complete"
echo "Files loaded from {original_tool} are ready for use"

"""
        script_content += self.generate_completion_check_footer()

        return script_content

    def get_config_display(self) -> List[str]:
        """Get Load configuration display."""
        config_lines = super().get_config_display()

        if self.loaded_result:
            original_tool = self.loaded_result.get('tool_name', 'Unknown')
            original_job = self.loaded_result.get('job_name')

            if not original_job or original_job == 'unknown':
                if ('configuration' in self.loaded_result and
                    'pipeline_context' in self.loaded_result['configuration']):
                    pipeline_context = self.loaded_result['configuration']['pipeline_context']
                    original_job = pipeline_context.get('pipeline_job_name', 'unknown')

            if not original_job:
                original_job = 'unknown'

            execution_order = self.loaded_result.get('execution_order', 'unknown')

            config_lines.extend([
                f"Loading from: {original_tool}",
                f"Original job: {original_job}",
                f"Original order: {execution_order}",
                f"Tool folder: {self.tool_folder}",
                f"File validation: {'enabled' if self.validate_files else 'disabled'}"
            ])

            output_structure = self.loaded_result['output_structure']
            for stream_name in ['structures', 'sequences', 'compounds']:
                if stream_name in output_structure and isinstance(output_structure[stream_name], dict):
                    count = len(output_structure[stream_name].get('ids', []))
                    if count > 0:
                        config_lines.append(f"Loaded {stream_name}: {count}")

            if 'tables' in output_structure:
                table_count = len(output_structure['tables'])
                config_lines.append(f"Loaded tables: {table_count}")

            if self.missing_files:
                config_lines.append(f"Missing files: {len(self.missing_files)}")

        return config_lines

    def get_loaded_metadata(self) -> Dict[str, Any]:
        """
        Get metadata about the loaded result.
        """
        if not self.loaded_result:
            return {}

        metadata = {
            'original_tool_name': self.loaded_result.get('tool_name'),
            'original_tool_class': self.loaded_result.get('tool_class'),
            'original_job_name': self.loaded_result.get('job_name'),
            'execution_order': self.loaded_result.get('execution_order'),
            'tool_folder': self.tool_folder,
            'missing_files_count': len(self.missing_files),
            'validation_enabled': self.validate_files
        }

        if 'execution_metadata' in self.loaded_result:
            metadata['execution_metadata'] = self.loaded_result['execution_metadata']

        if 'configuration' in self.loaded_result:
            metadata['original_configuration'] = self.loaded_result['configuration']

        return metadata

    def to_dict(self) -> Dict[str, Any]:
        """Serialize Load configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "load_params": {
                "tool_folder": self.tool_folder,
                "validate_files": self.validate_files,
                "original_tool_name": self.original_tool_name,
                "missing_files_count": len(self.missing_files) if self.missing_files else 0
            }
        })
        return base_dict

    def __str__(self) -> str:
        """String representation."""
        original_tool = self.original_tool_name or 'Unknown'
        return f"Load(from {original_tool}: {os.path.basename(self.tool_folder)})"


def LoadMultiple(path: str,
                tool: Optional[str] = None,
                suffix: Optional[str] = None,
                in_suffix: Optional[Union[str, List[str]]] = None,
                not_in_suffix: Optional[Union[str, List[str]]] = None,
                ascending: bool = True,
                **load_output_kwargs) -> Dict[str, Load]:
    """
    Load multiple tool outputs from a job folder by scanning tool subfolders.

    Scans for subfolders matching the NNN_ToolName[_Suffix] pattern and loads
    .expected_outputs.json from each.

    Args:
        path: Path to the job output folder (containing NNN_ToolName/ subfolders)
        tool: Filter by tool name (e.g., "Boltz2", "RFdiffusion")
        suffix: Filter by exact folder suffix
        in_suffix: Filter requiring suffix to contain string(s)
        not_in_suffix: Filter excluding suffix containing string(s)
        ascending: Sort order (True = ascending by name, False = descending)
        **load_output_kwargs: Additional parameters passed to Load constructor

    Returns:
        Dictionary mapping folder names to Load objects

    Examples:
        # Load all tool outputs from a job
        >>> data = LoadMultiple("/path/to/Job_001")

        # Load only Boltz2 outputs
        >>> boltz = LoadMultiple("/path/to/Job_001", tool="Boltz2")

        # Load outputs with a specific suffix
        >>> cycle10 = LoadMultiple("/path/to/Job_001", suffix="Cycle10")
    """
    job_folder = os.path.abspath(path)

    if not os.path.exists(job_folder):
        raise ValueError(f"Job folder not found: {job_folder}")

    if not os.path.isdir(job_folder):
        raise ValueError(f"Path is not a directory: {job_folder}")

    # Scan for tool subfolders matching NNN_ToolName[_Suffix] pattern
    tool_folder_pattern = re.compile(r'^(\d{3})_(.+)$')

    candidates = []
    for entry in os.listdir(job_folder):
        entry_path = os.path.join(job_folder, entry)
        if not os.path.isdir(entry_path):
            continue
        match = tool_folder_pattern.match(entry)
        if not match:
            continue
        # Check that .expected_outputs.json exists
        expected_json = os.path.join(entry_path, ".expected_outputs.json")
        if not os.path.exists(expected_json):
            continue
        candidates.append((entry, match.group(1), match.group(2)))  # (folder_name, index, rest)

    if not candidates:
        raise ValueError(f"No tool folders with .expected_outputs.json found in: {job_folder}")

    # Filter and load outputs
    loaded_outputs = {}

    for folder_name, index, rest in candidates:
        # Parse tool name and suffix from the rest (e.g., "Boltz2" or "Boltz2_Cycle10")
        parts = rest.split('_', 1)
        folder_tool_name = parts[0]
        folder_suffix = parts[1] if len(parts) > 1 else None

        # Apply tool name filter
        if tool is not None and folder_tool_name != tool:
            continue

        # Apply suffix filters
        if suffix is not None:
            if folder_suffix != suffix:
                continue

        if in_suffix is not None:
            if folder_suffix is None:
                continue
            checks = [in_suffix] if isinstance(in_suffix, str) else in_suffix
            if not all(c in folder_suffix for c in checks):
                continue

        if not_in_suffix is not None:
            if folder_suffix is not None:
                checks = [not_in_suffix] if isinstance(not_in_suffix, str) else not_in_suffix
                if any(c in folder_suffix for c in checks):
                    continue

        # Passed all filters — create Load
        folder_path = os.path.join(job_folder, folder_name)
        try:
            load_output = Load(folder_path, **load_output_kwargs)
            loaded_outputs[folder_name] = load_output
        except Exception as e:
            print(f"Warning: Could not create Load for {folder_name}: {e}")
            continue

    if not loaded_outputs:
        raise ValueError(
            f"No outputs matched the specified filters:\n"
            f"  Path: {job_folder}\n"
            f"  Tool: {tool if tool else 'any'}\n"
            f"  Suffix: {suffix if suffix else 'any'}\n"
            f"  Found {len(candidates)} tool folders, but none matched filters"
        )

    # Sort outputs by name
    sorted_outputs = dict(sorted(loaded_outputs.items(), reverse=not ascending))

    print(f"LoadMultiple: Loaded {len(sorted_outputs)} outputs from {job_folder}")
    if tool:
        print(f"  - Tool filter: {tool}")
    if suffix:
        print(f"  - Suffix filter: {suffix}")
    print(f"  - Sort order: {'ascending' if ascending else 'descending'}")
    print(f"  - Keys: {list(sorted_outputs.keys())}")

    return sorted_outputs
