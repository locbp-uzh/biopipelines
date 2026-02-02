"""
LoadOutput tool for loading previously saved pipeline tool outputs.

Allows reusing results from previous pipeline runs without re-execution,
enabling incremental pipeline development and efficient result reuse.
"""

import os
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


class LoadOutput(BaseConfig):
    """
    Load results from a previously executed pipeline tool, with optional filtering.
    
    This tool allows reusing outputs from previous pipeline runs without
    re-execution, enabling incremental pipeline development and result reuse.
    
    The LoadOutput tool:
    - Validates file existence and loads output metadata at pipeline runtime
    - Optionally applies filters to loaded data (expression string or Filter tool output)  
    - Provides filtered results through the standard pipeline interface
    - Enables subsequent tools to work with the filtered subset
    
    This is the ONLY tool that checks file existence at pipeline runtime (not SLURM runtime).
    """
    
    # Tool identification
    TOOL_NAME = "LoadOutput"
    
    
    def __init__(self, output_json: str, filter = None, validate_files: bool = True, **kwargs):
        """
        Initialize LoadOutput tool.
        
        Args:
            output_json: Path to saved tool result JSON file
            filter: Optional filter result (ToolOutput from Filter tool) to only load specific IDs
            validate_files: Whether to validate that all referenced files exist
            **kwargs: Additional parameters for BaseConfig
        """
        self.result_file = output_json
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
            # Extract job name from loaded result if available (will be processed below)
            
            # Try multiple locations for the job name
            original_job_name = self.loaded_result.get('job_name')
            
            if not original_job_name or original_job_name == 'unknown':
                # Try pipeline_job_name from pipeline context
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
        """Load and validate the result file."""
        if not os.path.exists(self.result_file):
            raise ValueError(f"Result file not found: {self.result_file}")

        try:
            with open(self.result_file, 'r') as f:
                self.loaded_result = json.load(f)
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON in result file {self.result_file}: {e}")
        except Exception as e:
            raise ValueError(f"Error loading result file {self.result_file}: {e}")

        # Validate required fields
        required_fields = ['tool_name', 'tool_class', 'output_structure']
        for field in required_fields:
            if field not in self.loaded_result:
                raise ValueError(f"Missing required field '{field}' in result file")

        self.original_tool_name = self.loaded_result['tool_name']

        # Validate file existence if requested
        if self.validate_files:
            self._validate_file_existence()
    
    def _validate_file_existence(self):
        """Validate that all referenced files in the output structure exist."""
        self.missing_files = []
        output_structure = self.loaded_result['output_structure']
        
        # Check structure files
        for file_list in ['structures', 'sequences', 'compounds']:
            if file_list in output_structure:
                for file_path in output_structure[file_list]:
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
                    elif isinstance(table_info, str) and not os.path.exists(table_info):
                        self.missing_files.append(table_info)
        
        # Check output folder
        if 'output_folder' in output_structure:
            output_folder = output_structure['output_folder']
            if not os.path.exists(output_folder):
                self.missing_files.append(output_folder)
        
        # Report missing files but don't fail (files might be moved/renamed)
        if self.missing_files:
            print(f"Warning: {len(self.missing_files)} files referenced in {self.result_file} are missing:")
            for missing_file in self.missing_files[:5]:  # Show first 5
                print(f"  - {missing_file}")
            if len(self.missing_files) > 5:
                print(f"  ... and {len(self.missing_files) - 5} more")

    def _resolve_file_paths(self, output_structure: Dict[str, Any]) -> Dict[str, Any]:
        """
        Resolve glob patterns and match IDs to actual file paths.

        This handles both legacy and new output formats:
        - Legacy: structures=[...], structure_ids=[...]
        - New: structures={ids:[...], files:[...], ...}

        Called at pipeline runtime when validate_files=True.

        Args:
            output_structure: The output structure to resolve

        Returns:
            Updated output structure with resolved file paths
        """
        import glob as glob_module

        resolved = output_structure.copy()

        for file_type, id_type in [('structures', 'structure_ids'),
                                    ('sequences', 'sequence_ids'),
                                    ('compounds', 'compound_ids')]:

            # Handle new DataStream dict format
            if file_type in resolved and isinstance(resolved[file_type], dict):
                ds_dict = resolved[file_type]
                files = ds_dict.get('files', [])
                ids = ds_dict.get('ids', [])
            # Handle legacy format
            elif file_type in resolved and id_type in resolved:
                files = resolved[file_type]
                ids = resolved[id_type]
            else:
                continue

            if not files or not ids:
                continue

            # Case 1: Already matched (same length, no globs)
            if len(files) == len(ids):
                has_glob = any('*' in f or '?' in f for f in files if isinstance(f, str))
                if not has_glob:
                    continue  # Already resolved

            # Case 2: Need to resolve - either glob pattern or mismatched lengths
            print(f"LoadOutput: Resolving {file_type} paths ({len(files)} patterns -> {len(ids)} IDs)")

            # Collect all potential files
            all_files = []
            for file_pattern in files:
                if not isinstance(file_pattern, str):
                    continue

                if '*' in file_pattern or '?' in file_pattern:
                    # Expand glob pattern
                    expanded = glob_module.glob(file_pattern)
                    all_files.extend(expanded)
                    print(f"  Expanded glob: {len(expanded)} files")
                elif os.path.isfile(file_pattern):
                    all_files.append(file_pattern)
                elif os.path.isdir(file_pattern):
                    # Search directory for structure files
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

                    # Direct match
                    if basename == item_id:
                        matched_file = file_path
                        break

                    # ID contained in filename
                    if item_id in basename:
                        matched_file = file_path
                        break

                    # Filename contained in ID
                    if basename in item_id:
                        matched_file = file_path
                        break

                    # Rank number matching (e.g., rank0001 -> rank0001_something.cif)
                    if item_id.startswith("rank") and basename.startswith("rank"):
                        id_num = ''.join(filter(str.isdigit, item_id))
                        # Extract rank number from basename (before first underscore)
                        file_rank_part = basename.split('_')[0]
                        file_num = ''.join(filter(str.isdigit, file_rank_part))
                        if id_num and file_num and id_num == file_num:
                            matched_file = file_path
                            break

                if matched_file:
                    resolved_files.append(matched_file)
                    # Update ID to match the actual filename (use basename without extension)
                    actual_id = os.path.splitext(os.path.basename(matched_file))[0]
                    resolved_ids.append(actual_id)
                else:
                    print(f"  Warning: No file found for ID '{item_id}'")

            # Store back in appropriate format
            if isinstance(resolved.get(file_type), dict):
                # New DataStream dict format - update in place
                resolved[file_type]['files'] = resolved_files
                resolved[file_type]['ids'] = resolved_ids
            else:
                # Legacy format - store as separate lists
                resolved[file_type] = resolved_files
                resolved[id_type] = resolved_ids

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
                # Look for main table (priority order)
                priority_names = ['structures', 'analysis', 'combined', 'results']
                for name in priority_names:
                    if name in tables:
                        ds_info = tables[name]
                        if isinstance(ds_info, dict) and 'path' in ds_info:
                            main_table_path = ds_info['path']
                            break
                        elif isinstance(ds_info, str):
                            main_table_path = ds_info
                            break
                
                # If no priority match, take first available
                if not main_table_path:
                    first_ds = next(iter(tables.values()))
                    if isinstance(first_ds, dict) and 'path' in first_ds:
                        main_table_path = first_ds['path']
                    elif isinstance(first_ds, str):
                        main_table_path = first_ds
        
        if not main_table_path or not os.path.exists(main_table_path):
            raise ValueError(f"Cannot find main table to filter in loaded output: {main_table_path}")
        
        print(f"LoadOutput: Applying filter to table: {main_table_path}")
        
        # Load the table
        try:
            df = pd.read_csv(main_table_path)
            print(f"  - Loaded {len(df)} rows")
        except Exception as e:
            raise ValueError(f"Error loading table for filtering: {e}")
        
        # Apply filter
        if isinstance(self.filter_input, str):
            # Simple expression filter
            try:
                filtered_df = df.query(self.filter_input)
                print(f"  - Applied expression: {self.filter_input}")
            except Exception as e:
                raise ValueError(f"Error applying filter expression '{self.filter_input}': {e}")
        
        elif hasattr(self.filter_input, 'tables'):
            # Precomputed filter result - load the filtered IDs
            filter_tables = self.filter_input.tables
            if isinstance(filter_tables, dict):
                # Find filtered results
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
                    # Use the 'id' column from filtered results to filter the original data
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

        # Store filtered IDs
        if 'id' in filtered_df.columns:
            self.filtered_ids = set(filtered_df['id'].tolist())
        else:
            # Fall back to row indices
            self.filtered_ids = set(filtered_df.index.tolist())

        print(f"  - Filtered IDs: {sorted(list(self.filtered_ids)) if len(self.filtered_ids) <= 10 else f'{len(self.filtered_ids)} items'}")

    def _apply_filter_to_output_structure(self, output_structure: Dict[str, Any]) -> Dict[str, Any]:
        """Apply filtering to the output structure based on filtered IDs."""
        
        print(f"LoadOutput: Filtering output structure to {len(self.filtered_ids)} items")
        
        filtered_structure = output_structure.copy()
        
        # Filter structures and structure_ids
        if 'structures' in output_structure and 'structure_ids' in output_structure:
            original_structures = output_structure['structures']
            original_structure_ids = output_structure['structure_ids']
            
            if len(original_structures) == len(original_structure_ids):
                filtered_structures = []
                filtered_structure_ids = []
                
                for structure_path, structure_id in zip(original_structures, original_structure_ids):
                    if structure_id in self.filtered_ids:
                        filtered_structures.append(structure_path)
                        filtered_structure_ids.append(structure_id)
                
                filtered_structure['structures'] = filtered_structures
                filtered_structure['structure_ids'] = filtered_structure_ids
                print(f"  - Filtered structures: {len(filtered_structures)}/{len(original_structures)}")
        
        # Filter compounds and compound_ids
        if 'compounds' in output_structure and 'compound_ids' in output_structure:
            original_compounds = output_structure['compounds']
            original_compound_ids = output_structure['compound_ids']

            # Filter compound_ids to only those in filtered set
            filtered_compound_ids = [cid for cid in original_compound_ids if cid in self.filtered_ids]

            # Handle compounds: update CSV paths to point to LoadOutput's filtered outputs
            if len(original_compounds) == len(original_compound_ids):
                # 1:1 mapping of files to IDs (e.g., one file per compound)
                filtered_compounds = []
                for compound_path, compound_id in zip(original_compounds, original_compound_ids):
                    if compound_id in self.filtered_ids:
                        # Update path to point to LoadOutput's output folder where filtered CSV will be
                        if isinstance(compound_path, str) and compound_path.endswith('.csv'):
                            new_path = os.path.join(self.output_folder, os.path.basename(compound_path))
                            filtered_compounds.append(new_path)
                        else:
                            filtered_compounds.append(compound_path)
            else:
                # Single CSV file containing multiple compounds (e.g., CompoundLibrary)
                # Update path to point to LoadOutput's output folder where filtered CSV will be created
                filtered_compounds = []
                for compound_path in original_compounds:
                    if isinstance(compound_path, str) and compound_path.endswith('.csv'):
                        new_path = os.path.join(self.output_folder, os.path.basename(compound_path))
                        filtered_compounds.append(new_path)
                    else:
                        filtered_compounds.append(compound_path)

            filtered_structure['compounds'] = filtered_compounds
            filtered_structure['compound_ids'] = filtered_compound_ids
            print(f"  - Filtered compounds: {len(filtered_compound_ids)}/{len(original_compound_ids)} IDs, {len(filtered_compounds)} file(s)")
        
        # Filter sequences and sequence_ids
        if 'sequences' in output_structure and 'sequence_ids' in output_structure:
            original_sequences = output_structure['sequences']
            original_sequence_ids = output_structure['sequence_ids']

            # Filter sequence_ids to only those in filtered set
            filtered_sequence_ids = [sid for sid in original_sequence_ids if sid in self.filtered_ids]

            # Handle sequences: update CSV paths to point to LoadOutput's filtered outputs
            if len(original_sequences) == len(original_sequence_ids):
                # 1:1 mapping of files to IDs (e.g., one file per sequence)
                filtered_sequences = []
                for sequence_path, sequence_id in zip(original_sequences, original_sequence_ids):
                    if sequence_id in self.filtered_ids:
                        # Update path to point to LoadOutput's output folder where filtered CSV will be
                        if isinstance(sequence_path, str) and sequence_path.endswith('.csv'):
                            new_path = os.path.join(self.output_folder, os.path.basename(sequence_path))
                            filtered_sequences.append(new_path)
                        else:
                            filtered_sequences.append(sequence_path)
            else:
                # Single CSV file containing multiple sequences
                # Update path to point to LoadOutput's output folder where filtered CSV will be created
                filtered_sequences = []
                for sequence_path in original_sequences:
                    if isinstance(sequence_path, str) and sequence_path.endswith('.csv'):
                        new_path = os.path.join(self.output_folder, os.path.basename(sequence_path))
                        filtered_sequences.append(new_path)
                    else:
                        filtered_sequences.append(sequence_path)

            filtered_structure['sequences'] = filtered_sequences
            filtered_structure['sequence_ids'] = filtered_sequence_ids
            print(f"  - Filtered sequences: {len(filtered_sequence_ids)}/{len(original_sequence_ids)} IDs, {len(filtered_sequences)} file(s)")

        # Update tables paths to point to LoadOutput's output folder and update counts
        if 'tables' in filtered_structure:
            for ds_name, ds_info in filtered_structure['tables'].items():
                if isinstance(ds_info, dict):
                    # Update path to LoadOutput's output folder
                    if 'path' in ds_info:
                        original_path = ds_info['path']
                        if isinstance(original_path, str) and original_path.endswith('.csv'):
                            ds_info['path'] = os.path.join(self.output_folder, f"{ds_name}.csv")

                    # Update count to reflect filtering
                    if 'count' in ds_info and isinstance(ds_info['count'], int):
                        ds_info['count'] = len(self.filtered_ids)
        
        return filtered_structure

    def _update_ids_from_csv(self, output_structure: Dict[str, Any]):
        """Update compound_ids and sequence_ids with actual IDs from CSV files."""
        import pandas as pd

        # Update compound_ids from compounds CSV files
        if 'compounds' in output_structure and output_structure['compounds']:
            updated_compound_ids = []
            for compound_file in output_structure['compounds']:
                if isinstance(compound_file, str) and os.path.exists(compound_file) and compound_file.endswith('.csv'):
                    df = pd.read_csv(compound_file)
                    if 'id' in df.columns:
                        updated_compound_ids.extend(df['id'].tolist())

            if updated_compound_ids:
                output_structure['compound_ids'] = updated_compound_ids

        # Update sequence_ids from sequences CSV files
        if 'sequences' in output_structure and output_structure['sequences']:
            updated_sequence_ids = []
            for sequence_file in output_structure['sequences']:
                if isinstance(sequence_file, str) and os.path.exists(sequence_file) and sequence_file.endswith('.csv'):
                    df = pd.read_csv(sequence_file)
                    if 'id' in df.columns:
                        updated_sequence_ids.extend(df['id'].tolist())

            if updated_sequence_ids:
                output_structure['sequence_ids'] = updated_sequence_ids

    def validate_params(self):
        """Validate LoadOutput parameters."""
        if not self.loaded_result:
            raise ValueError("No result loaded - initialization failed")
        
        if self.validate_files and len(self.missing_files) > 0:
            print(f"Warning: {len(self.missing_files)} referenced files are missing")
            print("Consider setting validate_files=False if files have been moved")
    
    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure LoadOutput - no inputs needed since we're loading existing results."""
        self.folders = pipeline_folders

        # No input configuration needed - we're loading pre-existing results
        pass

    def _exclude_missing_ids(self, output_structure: Dict[str, Any]) -> set:
        """
        Check for 'missing' table and return set of IDs to exclude.

        This handles cases where Filter tool has created a missing.csv listing
        filtered-out IDs that don't have corresponding files.

        Args:
            output_structure: The output structure to check

        Returns:
            Set of IDs that should be excluded (empty set if no missing table)
        """
        import pandas as pd

        # Check if there's a 'missing' table
        if 'tables' not in output_structure:
            return set()

        tables = output_structure['tables']
        if not isinstance(tables, dict) or 'missing' not in tables:
            return set()

        # Get the missing table path
        missing_info = tables['missing']
        if isinstance(missing_info, dict) and 'path' in missing_info:
            missing_path = missing_info['path']
        elif isinstance(missing_info, str):
            missing_path = missing_info
        else:
            return set()

        # Read the missing.csv file
        if not os.path.exists(missing_path):
            raise FileNotFoundError(f"Missing table referenced but file not found: {missing_path}")

        missing_df = pd.read_csv(missing_path)
        if 'id' not in missing_df.columns:
            raise ValueError(f"missing.csv does not have required 'id' column: {missing_path}")

        missing_ids = set(missing_df['id'].tolist())
        print(f"LoadOutput: Excluding {len(missing_ids)} IDs from missing.csv")
        return missing_ids

    def get_output_files(self) -> Dict[str, Any]:
        """
        Return the loaded output structure converted to new DataStream format.

        Handles legacy outputs (raw lists) and converts them to DataStream objects.
        When validate_files=True, resolves glob patterns to actual file paths.

        Returns:
            Dict with DataStream objects for structures, sequences, compounds,
            plus tables dict and output_folder string.
        """
        if not self.loaded_result:
            raise RuntimeError("No result loaded")

        # Get the original output structure
        output_structure = self.loaded_result['output_structure'].copy()

        # Backward compatibility: rename "datasheets" to "tables"
        if 'datasheets' in output_structure and 'tables' not in output_structure:
            print(f"LoadOutput: Converting legacy 'datasheets' to 'tables' for compatibility")
            output_structure['tables'] = output_structure.pop('datasheets')

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

            # Update compound_ids with actual IDs from CSV files (only when not filtering)
            self._update_ids_from_csv(output_structure)

        # Convert legacy format to DataStream format
        return self._convert_to_datastream_format(output_structure)

    def _filter_missing_ids(self, output_structure: Dict[str, Any], missing_ids: set):
        """Filter out missing IDs from output structure in place."""
        for file_type, id_type in [('structures', 'structure_ids'),
                                    ('compounds', 'compound_ids'),
                                    ('sequences', 'sequence_ids')]:
            if file_type in output_structure and id_type in output_structure:
                files = output_structure[file_type]
                ids = output_structure[id_type]
                if len(files) == len(ids):
                    filtered_files = []
                    filtered_ids = []
                    for f, i in zip(files, ids):
                        if i not in missing_ids:
                            filtered_files.append(f)
                            filtered_ids.append(i)
                    output_structure[file_type] = filtered_files
                    output_structure[id_type] = filtered_ids
                    if len(filtered_files) < len(files):
                        print(f"  - Filtered {file_type}: {len(filtered_files)}/{len(files)} kept")

    def _convert_to_datastream_format(self, output_structure: Dict[str, Any]) -> Dict[str, Any]:
        """
        Convert legacy output format (raw lists) to new DataStream format.

        Handles both:
        - Legacy format: structures=[...], structure_ids=[...]
        - New format: structures={name:..., ids:..., files:..., format:...}

        Args:
            output_structure: Dict with structures/structure_ids lists or DataStream dicts

        Returns:
            Dict with DataStream objects for structures, sequences, compounds
        """
        from .datastream import DataStream, create_map_table

        output_folder = output_structure.get('output_folder', '')

        # Helper to detect format from file extension
        def detect_format(files: List[str], default: str) -> str:
            for f in files:
                if isinstance(f, str):
                    ext = os.path.splitext(f)[1].lower()
                    if ext in ['.pdb']:
                        return 'pdb'
                    elif ext in ['.cif', '.mmcif']:
                        return 'cif'
                    elif ext in ['.fasta', '.fa']:
                        return 'fasta'
                    elif ext in ['.csv']:
                        return 'csv'
                    elif ext in ['.sdf', '.mol2']:
                        return 'sdf'
            return default

        # Helper to convert a single data type (structures, sequences, or compounds)
        def convert_data_type(data_key: str, ids_key: str, default_format: str) -> DataStream:
            raw_data = output_structure.get(data_key)

            # Case 1: Already a DataStream dict (new format)
            if isinstance(raw_data, dict) and 'ids' in raw_data and 'files' in raw_data:
                return DataStream.from_dict(raw_data)

            # Case 2: Legacy format - separate lists
            files = raw_data if isinstance(raw_data, list) else []
            ids = output_structure.get(ids_key, [])

            if files and ids and len(files) == len(ids):
                fmt = detect_format(files, default_format)
                return DataStream(
                    name=data_key,
                    ids=ids,
                    files=files,
                    map_table="",
                    format=fmt
                )
            elif ids and not files:
                # IDs but no files (value-based like SMILES)
                return DataStream(
                    name=data_key,
                    ids=ids,
                    files=[],
                    map_table="",
                    format=default_format
                )
            else:
                return DataStream.empty(data_key, default_format)

        # Convert all data types using the helper
        structures = convert_data_type('structures', 'structure_ids', 'pdb')
        sequences = convert_data_type('sequences', 'sequence_ids', 'fasta')
        compounds = convert_data_type('compounds', 'compound_ids', 'sdf')

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
                elif isinstance(table_info, str):
                    # Legacy format: just a path string
                    tables_dict[table_name] = TableInfo(
                        name=table_name,
                        path=table_info,
                        columns=[],
                        description='',
                        count=0
                    )
                elif isinstance(table_info, TableInfo):
                    # Already a TableInfo object
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
        
        Args:
            script_path: Path where script should be written
            
        Returns:
            Script content as string
        """
        original_tool = self.loaded_result.get('tool_name', 'Unknown')
        # Try multiple locations for the job name
        original_job = self.loaded_result.get('job_name')
        
        if not original_job or original_job == 'unknown':
            # Try pipeline_job_name from pipeline context
            if ('configuration' in self.loaded_result and 
                'pipeline_context' in self.loaded_result['configuration']):
                pipeline_context = self.loaded_result['configuration']['pipeline_context']
                original_job = pipeline_context.get('pipeline_job_name', 'unknown')
        
        if not original_job:
            original_job = 'unknown'
        
        script_content = "#!/bin/bash\n"
        script_content += f"# LoadOutput script - loading results from {original_tool}\n"
        script_content += f"# Original job: {original_job}\n"
        script_content += f"# Result file: {self.result_file}\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += f"""
echo "Loading output from previous {original_tool} execution"
echo "Original job: {original_job}"
echo "Result file: {self.result_file}"

# Validate that key files exist
echo "Validating loaded files..."

# Check output folder
if [ ! -d "{self.loaded_result['output_structure'].get('output_folder', '')}" ]; then
    echo "Warning: Output folder missing: {self.loaded_result['output_structure'].get('output_folder', '')}"
fi

# Check some key files (first few from each category)
"""
        
        # Add file existence checks for first few files in each category
        output_structure = self.loaded_result['output_structure']
        
        for file_type in ['structures', 'sequences', 'compounds']:
            if file_type in output_structure and output_structure[file_type]:
                files_to_check = output_structure[file_type][:3]  # Check first 3
                for file_path in files_to_check:
                    if isinstance(file_path, str):
                        script_content += f"""
if [ ! -f "{file_path}" ]; then
    echo "Warning: {file_type} file missing: {file_path}"
fi"""
        
        # Add filtering section if filter was applied
        if self.filter_input and self.filtered_ids:
            # Create filter configuration file
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

echo "LoadOutput complete"
echo "Files loaded from {original_tool} are ready for use"

"""
        script_content += self.generate_completion_check_footer()

        return script_content
    
    def get_config_display(self) -> List[str]:
        """Get LoadOutput configuration display."""
        config_lines = super().get_config_display()
        
        if self.loaded_result:
            original_tool = self.loaded_result.get('tool_name', 'Unknown')
            # Try multiple locations for the job name
        original_job = self.loaded_result.get('job_name')
        
        if not original_job or original_job == 'unknown':
            # Try pipeline_job_name from pipeline context
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
                f"Result file: {os.path.basename(self.result_file)}",
                f"File validation: {'enabled' if self.validate_files else 'disabled'}"
            ])
            
            # Show output summary
            output_structure = self.loaded_result['output_structure']
            for data_type in ['structures', 'sequences', 'compounds']:
                if data_type in output_structure and output_structure[data_type]:
                    count = len(output_structure[data_type])
                    config_lines.append(f"Loaded {data_type}: {count}")
            
            if 'tables' in output_structure:
                table_count = len(output_structure['tables'])
                config_lines.append(f"Loaded tables: {table_count}")
            
            # Warn about missing files
            if self.missing_files:
                config_lines.append(f"Missing files: {len(self.missing_files)}")
        
        return config_lines
    
    def get_loaded_metadata(self) -> Dict[str, Any]:
        """
        Get metadata about the loaded result.
        
        Returns:
            Dictionary with metadata about the original tool execution
        """
        if not self.loaded_result:
            return {}
        
        metadata = {
            'original_tool_name': self.loaded_result.get('tool_name'),
            'original_tool_class': self.loaded_result.get('tool_class'),
            'original_job_name': self.loaded_result.get('job_name'),
            'execution_order': self.loaded_result.get('execution_order'),
            'result_file': self.result_file,
            'missing_files_count': len(self.missing_files),
            'validation_enabled': self.validate_files
        }
        
        # Add execution metadata if available
        if 'execution_metadata' in self.loaded_result:
            metadata['execution_metadata'] = self.loaded_result['execution_metadata']
        
        # Add configuration metadata if available
        if 'configuration' in self.loaded_result:
            metadata['original_configuration'] = self.loaded_result['configuration']
        
        return metadata
    
    def to_dict(self) -> Dict[str, Any]:
        """Serialize LoadOutput configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "load_params": {
                "result_file": self.result_file,
                "validate_files": self.validate_files,
                "original_tool_name": self.original_tool_name,
                "missing_files_count": len(self.missing_files) if self.missing_files else 0
            }
        })
        return base_dict
    
    def __str__(self) -> str:
        """String representation."""
        original_tool = self.original_tool_name or 'Unknown'
        return f"LoadOutput(from {original_tool}: {os.path.basename(self.result_file)})"


def LoadOutputs(path: str,
                tool: Optional[str] = None,
                suffix: Optional[str] = None,
                in_suffix: Optional[Union[str, List[str]]] = None,
                not_in_suffix: Optional[Union[str, List[str]]] = None,
                ascending: bool = True,
                **load_output_kwargs) -> Dict[str, LoadOutput]:
    """
    Load multiple tool outputs from a ToolOutputs folder with filtering and sorting.

    This helper function scans a ToolOutputs folder and loads multiple outputs
    matching specified criteria, returning them as a dictionary for easy access.

    Args:
        path: Path to the job folder or ToolOutputs folder
              - If job folder: automatically appends "ToolOutputs"
              - If ToolOutputs folder: uses directly
        tool: Filter by tool name (e.g., "MergeTables", "Filter", "Boltz2")
              - Matches the tool name in JSON files
        suffix: Filter by filename suffix (e.g., "Cycle10", "Affinity")
                - Matches pattern: NNN_ToolName_Suffix.json
        [not]_in_suffix: Filter by filename suffix, whether a string or list of strings is contained in it or not
        ascending: Sort order (True = ascending by name, False = descending by name)
        **load_output_kwargs: Additional parameters passed to LoadOutput constructor
                              (e.g., filter, validate_files)

    Returns:
        Dictionary mapping identifiers to LoadOutput objects
        - Keys are generated from: "{execution_order}_{tool_name}_{suffix}" (if suffix)
                               or: "{execution_order}_{tool_name}" (otherwise)
        - Values are LoadOutput instances

    Examples:
        # Load all MergeTables outputs from a job
        >>> data = LoadOutputs(
        ...     path="/shares/user/BioPipelines/Project/Job_001",
        ...     tool="MergeTables"
        ... )
        >>> # Access outputs: data["003_MergeTables"], data["007_MergeTables"], etc.

        # Load all outputs from Cycle10 with validation disabled
        >>> cycle10_data = LoadOutputs(
        ...     path="/shares/user/BioPipelines/Project/Job_001/ToolOutputs",
        ...     suffix="Cycle10",
        ...     validate_files=False
        ... )

        # Load all Filter outputs in descending order
        >>> filters = LoadOutputs(
        ...     path="/shares/user/BioPipelines/Project/Job_001",
        ...     tool="Filter",
        ...     ascending=False
        ... )

        # Combine filters: specific tool with suffix
        >>> merged_cycle5 = LoadOutputs(
        ...     path="/shares/user/BioPipelines/Project/Job_001",
        ...     tool="MergeTables",
        ...     suffix="Cycle5"
        ... )
    """
    # Normalize path - check if it's ToolOutputs folder or job folder
    if os.path.basename(path) == "ToolOutputs":
        tool_outputs_folder = path
    else:
        # Assume it's a job folder - append ToolOutputs
        tool_outputs_folder = os.path.join(path, "ToolOutputs")

    # Validate folder exists
    if not os.path.exists(tool_outputs_folder):
        raise ValueError(f"ToolOutputs folder not found: {tool_outputs_folder}")

    if not os.path.isdir(tool_outputs_folder):
        raise ValueError(f"Path is not a directory: {tool_outputs_folder}")

    # Scan for JSON files
    json_files = [f for f in os.listdir(tool_outputs_folder) if f.endswith('.json')]

    if not json_files:
        raise ValueError(f"No JSON files found in: {tool_outputs_folder}")

    # Filter and load outputs
    loaded_outputs = {}

    for json_filename in json_files:
        json_path = os.path.join(tool_outputs_folder, json_filename)

        # Load JSON to check filters
        try:
            with open(json_path, 'r') as f:
                json_data = json.load(f)
        except (json.JSONDecodeError, Exception) as e:
            print(f"Warning: Could not load {json_filename}: {e}")
            continue

        # Apply filters
        # 1. Tool name filter
        if tool is not None:
            tool_name = json_data.get('tool_name', '')
            if tool_name != tool:
                continue

        # 2. Suffix filter
        # Parse filename: NNN_ToolName.json or NNN_ToolName_Suffix.json
        if suffix is not None or in_suffix is not None or not_in_suffix is not None:
            # Extract suffix from filename (everything between last _ and .json)
            # Pattern: NNN_ToolName_Suffix.json
            basename = os.path.splitext(json_filename)[0]  # Remove .json
            parts = basename.split('_')

            # Check if there's a suffix (more than 2 parts: index, toolname, suffix...)
            if len(parts) > 2:
                # Everything after tool name is considered suffix
                file_suffix = '_'.join(parts[2:])
                if suffix is not None:
                    if file_suffix != suffix:
                        continue
                if in_suffix is not None:
                    if isinstance(in_suffix,str):
                        if not in_suffix in suffix:
                            continue
                    else: #list of strings
                        for cond in in_suffix:
                            if not cond in suffix:
                                continue
                if not_in_suffix is not None:
                    if isinstance(not_in_suffix,str):
                        if not not_in_suffix in suffix:
                            continue
                    else: #list of strings
                        for cond in not_in_suffix:
                            if cond in suffix:
                                continue
            else:
                # No suffix in filename, skip this file
                continue

        # File passed all filters - create LoadOutput
        try:
            load_output = LoadOutput(json_path, **load_output_kwargs)

            # Generate key from filename components
            basename = os.path.splitext(json_filename)[0]  # Remove .json

            # Use basename as key (includes index, tool name, and suffix if present)
            loaded_outputs[basename] = load_output

        except Exception as e:
            print(f"Warning: Could not create LoadOutput for {json_filename}: {e}")
            continue

    if not loaded_outputs:
        raise ValueError(
            f"No outputs matched the specified filters:\n"
            f"  Path: {tool_outputs_folder}\n"
            f"  Tool: {tool if tool else 'any'}\n"
            f"  Suffix: {suffix if suffix else 'any'}\n"
            f"  Found {len(json_files)} JSON files, but none matched filters"
        )

    # Sort outputs by name (ascending or descending)
    sorted_outputs = dict(sorted(loaded_outputs.items(), reverse=not ascending))

    print(f"LoadOutputs: Loaded {len(sorted_outputs)} outputs from {tool_outputs_folder}")
    if tool:
        print(f"  - Tool filter: {tool}")
    if suffix:
        print(f"  - Suffix filter: {suffix}")
    print(f"  - Sort order: {'ascending' if ascending else 'descending'}")
    print(f"  - Keys: {list(sorted_outputs.keys())}")

    return sorted_outputs