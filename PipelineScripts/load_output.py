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
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput


class LoadOutput(BaseConfig):
    """
    Load results from a previously executed pipeline tool.
    
    This tool allows reusing outputs from previous pipeline runs without
    re-execution, enabling incremental pipeline development and result reuse.
    
    The LoadOutput tool validates file existence, loads output metadata,
    and provides the loaded results through the standard pipeline interface.
    """
    
    # Tool identification
    TOOL_NAME = "LoadOutput"
    DEFAULT_ENV = "ProteinEnv"  # Can run in any environment since it just loads files
    COMPATIBLE_ENVS = ["ProteinEnv", "Boltz2Env", "ligandmpnn_env"]
    DEFAULT_RESOURCES = {"memory": "1GB", "time": "0:10:00"}  # Minimal resources needed
    
    def __init__(self, result_file: str, validate_files: bool = True, **kwargs):
        """
        Initialize LoadOutput tool.
        
        Args:
            result_file: Path to saved tool result JSON file
            validate_files: Whether to validate that all referenced files exist
            **kwargs: Additional parameters for BaseConfig
        """
        self.result_file = result_file
        self.validate_files = validate_files
        self.loaded_result = None
        self.missing_files = []
        self.original_tool_name = None
        
        # Load and validate the result file
        self._load_and_validate_result()
        
        # Set up job name based on loaded result
        if not kwargs.get('job_name'):
            # Debug: Print available keys to understand the structure
            if self.loaded_result:
                print(f"Debug: Available keys in loaded result: {list(self.loaded_result.keys())}")
                if 'job_name' in self.loaded_result:
                    print(f"Debug: job_name found: {self.loaded_result['job_name']}")
                else:
                    print("Debug: job_name not found at top level")
                    # Check if job_name is nested under configuration
                    if 'configuration' in self.loaded_result and 'pipeline_context' in self.loaded_result['configuration']:
                        pipeline_context = self.loaded_result['configuration']['pipeline_context']
                        if 'pipeline_job_name' in pipeline_context:
                            print(f"Debug: Found pipeline_job_name: {pipeline_context['pipeline_job_name']}")
            
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
        Return the loaded output structure.
        
        Returns:
            The complete output structure from the saved result
        """
        if not self.loaded_result:
            raise RuntimeError("No result loaded")
        
        # Return the saved output structure
        output_structure = self.loaded_result['output_structure'].copy()
        
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