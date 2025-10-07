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
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo


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
    DEFAULT_ENV = None  # Loaded from config.yaml
    
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
        
        # Check datasheet files
        if 'datasheets' in output_structure:
            datasheets = output_structure['datasheets']
            if isinstance(datasheets, dict):
                for datasheet_name, datasheet_info in datasheets.items():
                    if isinstance(datasheet_info, dict) and 'path' in datasheet_info:
                        path = datasheet_info['path']
                        if not os.path.exists(path):
                            self.missing_files.append(path)
                    elif isinstance(datasheet_info, str) and not os.path.exists(datasheet_info):
                        self.missing_files.append(datasheet_info)
        
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
    
    def _process_filter(self):
        """Process filter input to determine which IDs to keep."""
        import pandas as pd
        
        # Find the main datasheet to filter
        output_structure = self.loaded_result['output_structure']
        main_datasheet_path = None
        
        if 'datasheets' in output_structure:
            datasheets = output_structure['datasheets']
            if isinstance(datasheets, dict):
                # Look for main datasheet (priority order)
                priority_names = ['structures', 'analysis', 'combined', 'results']
                for name in priority_names:
                    if name in datasheets:
                        ds_info = datasheets[name]
                        if isinstance(ds_info, dict) and 'path' in ds_info:
                            main_datasheet_path = ds_info['path']
                            break
                        elif isinstance(ds_info, str):
                            main_datasheet_path = ds_info
                            break
                
                # If no priority match, take first available
                if not main_datasheet_path:
                    first_ds = next(iter(datasheets.values()))
                    if isinstance(first_ds, dict) and 'path' in first_ds:
                        main_datasheet_path = first_ds['path']
                    elif isinstance(first_ds, str):
                        main_datasheet_path = first_ds
        
        if not main_datasheet_path or not os.path.exists(main_datasheet_path):
            raise ValueError(f"Cannot find main datasheet to filter in loaded output: {main_datasheet_path}")
        
        print(f"LoadOutput: Applying filter to datasheet: {main_datasheet_path}")
        
        # Load the datasheet
        try:
            df = pd.read_csv(main_datasheet_path)
            print(f"  - Loaded {len(df)} rows")
        except Exception as e:
            raise ValueError(f"Error loading datasheet for filtering: {e}")
        
        # Apply filter
        if isinstance(self.filter_input, str):
            # Simple expression filter
            try:
                filtered_df = df.query(self.filter_input)
                print(f"  - Applied expression: {self.filter_input}")
            except Exception as e:
                raise ValueError(f"Error applying filter expression '{self.filter_input}': {e}")
        
        elif hasattr(self.filter_input, 'datasheets'):
            # Precomputed filter result - load the filtered IDs
            filter_datasheets = self.filter_input.datasheets
            if isinstance(filter_datasheets, dict):
                # Find filtered results
                if 'filtered' in filter_datasheets:
                    filter_ds_info = filter_datasheets['filtered']
                else:
                    filter_ds_info = next(iter(filter_datasheets.values()))
                
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

        # IMPORTANT: Create filtered copies of ALL datasheets in LoadOutput's output folder
        # This is necessary so downstream tools receive the filtered data, not the original
        self._create_filtered_datasheets(filtered_df)

    def _create_filtered_datasheets(self, filtered_main_df):
        """
        Create filtered copies of all datasheets in LoadOutput's output folder.

        This ensures downstream tools receive filtered data, not the original unfiltered files.
        """
        import pandas as pd

        # Ensure output folder exists (it should be set during pipeline add)
        if not hasattr(self, 'output_folder') or not self.output_folder:
            raise ValueError("LoadOutput output_folder not set. Cannot create filtered datasheets.")

        os.makedirs(self.output_folder, exist_ok=True)
        print(f"LoadOutput: Creating filtered datasheets in {self.output_folder}")

        output_structure = self.loaded_result['output_structure']

        # Process all datasheets
        if 'datasheets' in output_structure:
            datasheets = output_structure['datasheets']
            if isinstance(datasheets, dict):
                for ds_name, ds_info in datasheets.items():
                    # Get original datasheet path
                    original_path = None
                    if isinstance(ds_info, dict) and 'path' in ds_info:
                        original_path = ds_info['path']
                    elif isinstance(ds_info, str):
                        original_path = ds_info

                    if not original_path or not os.path.exists(original_path):
                        print(f"  Warning: Skipping datasheet '{ds_name}' - file not found: {original_path}")
                        continue

                    try:
                        # Load original datasheet
                        df = pd.read_csv(original_path)

                        # Filter by IDs
                        if 'id' in df.columns:
                            filtered_df = df[df['id'].isin(self.filtered_ids)]
                        else:
                            # If no 'id' column, skip filtering for this datasheet
                            print(f"  Warning: Datasheet '{ds_name}' has no 'id' column, copying as-is")
                            filtered_df = df

                        # Save filtered datasheet to LoadOutput's output folder
                        filtered_path = os.path.join(self.output_folder, f"{ds_name}.csv")
                        filtered_df.to_csv(filtered_path, index=False)
                        print(f"  ✓ {ds_name}.csv: {len(filtered_df)}/{len(df)} rows")

                        # Update the datasheet info to point to the new filtered file
                        if isinstance(ds_info, dict):
                            ds_info['path'] = filtered_path
                            ds_info['count'] = len(filtered_df)
                        else:
                            # If it was just a string, replace it with updated path
                            datasheets[ds_name] = filtered_path

                    except Exception as e:
                        print(f"  Error filtering datasheet '{ds_name}': {e}")
                        continue

        # Also handle special compound/sequence lists if they exist as CSV files
        for list_name in ['compounds', 'sequences']:
            if list_name in output_structure:
                file_list = output_structure[list_name]
                if isinstance(file_list, list):
                    filtered_files = []
                    for file_path in file_list:
                        if isinstance(file_path, str) and file_path.endswith('.csv') and os.path.exists(file_path):
                            try:
                                df = pd.read_csv(file_path)
                                if 'id' in df.columns:
                                    filtered_df = df[df['id'].isin(self.filtered_ids)]
                                    # Save filtered version
                                    filtered_path = os.path.join(self.output_folder, os.path.basename(file_path))
                                    filtered_df.to_csv(filtered_path, index=False)
                                    filtered_files.append(filtered_path)
                                    print(f"  ✓ {os.path.basename(file_path)}: {len(filtered_df)}/{len(df)} rows")
                                else:
                                    # No filtering possible, keep original
                                    filtered_files.append(file_path)
                            except Exception as e:
                                print(f"  Error filtering {file_path}: {e}")
                                filtered_files.append(file_path)
                        else:
                            # Not a CSV or doesn't exist, keep original
                            filtered_files.append(file_path)

                    # Update the list with filtered file paths
                    output_structure[list_name] = filtered_files

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
            
            if len(original_compounds) == len(original_compound_ids):
                filtered_compounds = []
                filtered_compound_ids = []
                
                for compound_path, compound_id in zip(original_compounds, original_compound_ids):
                    if compound_id in self.filtered_ids:
                        filtered_compounds.append(compound_path)
                        filtered_compound_ids.append(compound_id)
                
                filtered_structure['compounds'] = filtered_compounds
                filtered_structure['compound_ids'] = filtered_compound_ids
        
        # Filter sequences and sequence_ids
        if 'sequences' in output_structure and 'sequence_ids' in output_structure:
            original_sequences = output_structure['sequences']
            original_sequence_ids = output_structure['sequence_ids']
            
            if len(original_sequences) == len(original_sequence_ids):
                filtered_sequences = []
                filtered_sequence_ids = []
                
                for sequence_path, sequence_id in zip(original_sequences, original_sequence_ids):
                    if sequence_id in self.filtered_ids:
                        filtered_sequences.append(sequence_path)
                        filtered_sequence_ids.append(sequence_id)
                
                filtered_structure['sequences'] = filtered_sequences
                filtered_structure['sequence_ids'] = filtered_sequence_ids
        
        # Update datasheets count information
        if 'datasheets' in filtered_structure:
            for ds_name, ds_info in filtered_structure['datasheets'].items():
                if isinstance(ds_info, dict) and 'count' in ds_info:
                    # Update count to reflect filtering
                    if isinstance(ds_info['count'], int):
                        # Rough estimate - could be more precise if we tracked per-datasheet filtering
                        filtered_structure['datasheets'][ds_name]['count'] = len(self.filtered_ids)
        
        return filtered_structure

    def _update_ids_from_csv(self, output_structure: Dict[str, Any]):
        """Update compound_ids and sequence_ids with actual IDs from CSV files."""
        import pandas as pd

        # Update compound_ids from compounds CSV files
        if 'compounds' in output_structure and output_structure['compounds']:
            updated_compound_ids = []
            for compound_file in output_structure['compounds']:
                if isinstance(compound_file, str) and os.path.exists(compound_file) and compound_file.endswith('.csv'):
                    try:
                        df = pd.read_csv(compound_file)
                        if 'id' in df.columns:
                            updated_compound_ids.extend(df['id'].tolist())
                    except Exception as e:
                        # If we can't read the file, keep the original IDs
                        print(f"Warning: Could not read compound file {compound_file}: {e}")
                        continue

            if updated_compound_ids:
                output_structure['compound_ids'] = updated_compound_ids

        # Update sequence_ids from sequences CSV files
        if 'sequences' in output_structure and output_structure['sequences']:
            updated_sequence_ids = []
            for sequence_file in output_structure['sequences']:
                if isinstance(sequence_file, str) and os.path.exists(sequence_file) and sequence_file.endswith('.csv'):
                    try:
                        df = pd.read_csv(sequence_file)
                        if 'id' in df.columns:
                            updated_sequence_ids.extend(df['id'].tolist())
                    except Exception as e:
                        # If we can't read the file, keep the original IDs
                        print(f"Warning: Could not read sequence file {sequence_file}: {e}")
                        continue

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
    
    def get_output_files(self) -> Dict[str, Any]:
        """
        Return the loaded output structure, optionally filtered.

        Returns:
            The output structure from the saved result, filtered if filter was applied
        """
        if not self.loaded_result:
            raise RuntimeError("No result loaded")

        # Get the original output structure
        output_structure = self.loaded_result['output_structure'].copy()

        # Apply filtering if filter was provided
        if self.filtered_ids is not None:
            output_structure = self._apply_filter_to_output_structure(output_structure)

        # Update compound_ids with actual IDs from CSV files
        self._update_ids_from_csv(output_structure)

        # Ensure we have all expected keys for compatibility
        default_structure = {
            "structures": [],
            "structure_ids": [],
            "compounds": [],
            "compound_ids": [],
            "sequences": [],
            "sequence_ids": [],
            "datasheets": {},
            "output_folder": ""
        }

        # Merge with defaults
        for key, default_value in default_structure.items():
            if key not in output_structure:
                output_structure[key] = default_value

        return output_structure
    
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
        
        script_content = f"""#!/bin/bash
# LoadOutput script - loading results from {original_tool}
# Original job: {original_job}
# Result file: {self.result_file}

{self.generate_completion_check_header()}

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
        
        script_content += f"""

echo "LoadOutput validation complete"
echo "Files loaded from {original_tool} are ready for use"

{self.generate_completion_check_footer()}
"""
        
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
            
            if 'datasheets' in output_structure:
                datasheet_count = len(output_structure['datasheets'])
                config_lines.append(f"Loaded datasheets: {datasheet_count}")
            
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