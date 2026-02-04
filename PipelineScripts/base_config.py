"""
Base configuration class for all modeling tools.

Provides common interface and functionality for tool configuration,
validation, and integration with the pipeline system.
"""

import pandas as pd
import os
import json
from abc import ABC, abstractmethod
from typing import Dict, List, Any, Optional, Union, TypeVar, overload, Type

try:
    from .config_manager import ConfigManager
    from .file_paths import Path
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from config_manager import ConfigManager
    from file_paths import Path


# TypeVar for preserving concrete tool types in type hints
T = TypeVar('T', bound='BaseConfig')


class BaseConfig(ABC):
    """
    Abstract base class for all tool configurations.

    Provides common functionality for parameter validation,
    environment management, and integration with Pipeline.
    """

    # Tool-specific defaults (override in subclasses)
    TOOL_NAME = "base"

    # Common path descriptors available to all tools
    pipeline_name = Path(lambda self: self._extract_pipeline_name())
    log_file = Path(lambda self: self._compute_log_file_path())
    

    @overload
    def __new__(cls: Type[T], *args, **kwargs) -> T:
        """
        Type hint overload for IDE autocomplete support.

        This tells type checkers that Boltz(...) returns a Boltz instance,
        enabling IDE parameter suggestions for tool constructors even though
        the actual runtime returns StandardizedOutput when in Pipeline context.
        """
        ...

    def __new__(cls, *args, **kwargs):
        """
        Create a new tool instance with optional auto-registration.

        If called within a Pipeline context manager, this automatically registers
        the tool and returns a ToolOutput object instead of the tool instance itself.
        This enables clean syntax like:

            with Pipeline(...) as pipeline:
                rfdaa = RFdiffusionAllAtom(...)  # Returns ToolOutput
                distances = DistanceSelector(structures=rfdaa, ...)  # Also ToolOutput

        Outside a Pipeline context, returns the tool instance normally for
        backward compatibility with explicit pipeline.add() usage.

        Args:
            *args: Positional arguments for tool initialization
            **kwargs: Keyword arguments for tool initialization

        Returns:
            ToolOutput if within Pipeline context, tool instance otherwise
        """
        # Create the tool instance normally
        instance = super(BaseConfig, cls).__new__(cls)

        # Check for active pipeline context
        # Import here to avoid circular dependency
        from .pipeline import Pipeline
        active_pipeline = Pipeline.get_active_pipeline()

        if active_pipeline is not None:
            # We're in a Pipeline context - auto-register this tool
            # Initialize the instance first so it's properly configured
            instance.__init__(*args, **kwargs)

            # Auto-register with the pipeline and get ToolOutput
            tool_output = active_pipeline._auto_register(instance)

            # Return the StandardizedOutput from ToolOutput for chaining
            return tool_output.output
        else:
            # No active pipeline - return the tool instance normally
            # (will be initialized by Python calling __init__ automatically)
            return instance

    def __init__(self, **kwargs):
        """Initialize base configuration with common parameters."""
        # Core identification
        self.tool_name = self.TOOL_NAME
        self.job_name = kwargs.get('name', '')

        # Pipeline reference for getting job name
        self.pipeline = kwargs.get('pipeline', None)

        # Environment(s) from config.yaml
        self._load_environments()
        self.resources = kwargs.get('resources', {})
        
        # Pipeline integration
        self.dependencies = kwargs.get('dependencies', [])
        self.pipeline_ref = None  # Set by Pipeline when added
        self.execution_order = 0   # Set by Pipeline
        
        # I/O tracking
        self.input_sources = {}    # What this tool takes as input
        self.output_files = {}     # What this tool produces
        self.output_folder = ""    # Set when added to pipeline
        
        # Execution state
        self.configured = False
        self.executed = False
        
        # Store all parameters for validation and serialization
        self.params = kwargs

        # Validate configuration
        self.validate_params()

    def _load_environments(self):
        """
        Load environment(s) for this tool from config.yaml.

        The environments section in config.yaml can specify either:
        - A single environment string: "ProteinEnv"
        - A list of environments: ["dynamicbind", "relax"]

        Tools that need multiple environments (like DynamicBind) can access
        them via activate_environment(index=N).
        """
        config_manager = ConfigManager()
        env_config = config_manager.get_environment(self.TOOL_NAME)

        if env_config is None:
            # Default to biopipelines if not configured
            self.environments = ["biopipelines"]
        elif isinstance(env_config, list):
            # Multiple environments specified
            self.environments = env_config
        else:
            # Single environment string
            self.environments = [env_config]

    def activate_environment(self, index: int = 0, name: Optional[str] = None) -> str:
        """
        Generate bash script snippet to activate a conda environment.

        Args:
            index: Index of the environment to activate (default: 0, the first/primary environment)
            name: Optional explicit environment name (overrides index-based lookup)

        Returns:
            Bash script content for activating the environment with diagnostics

        Raises:
            IndexError: If index is out of range for configured environments
        """
        if name is not None:
            env_name = name
        else:
            if index >= len(self.environments):
                raise IndexError(
                    f"Environment index {index} out of range for {self.TOOL_NAME}. "
                    f"Available environments: {self.environments}"
                )
            env_name = self.environments[index]

        return f"""# Activate environment: {env_name}
echo "=== Activating Environment ==="
echo "Requested: {env_name}"
eval "$(mamba shell hook --shell bash)"
mamba activate {env_name}
echo "Environment: $CONDA_DEFAULT_ENV"
echo "Location: $CONDA_PREFIX"
echo "Python: $(which python)"
echo "Python version: $(python --version 2>&1)"
echo "=============================="

"""

    @abstractmethod
    def validate_params(self):
        """Validate tool-specific parameters. Override in subclasses."""
        pass

    def _extract_pipeline_name(self) -> str:
        """
        Extract pipeline name from output folder structure.

        The output folder follows the pattern: .../PipelineName_NNN/NNN_ToolName/
        This method extracts PipelineName from that structure using TOOL_NAME.

        Returns:
            Pipeline name string

        Raises:
            ValueError: If pipeline name cannot be extracted
        """
        folder_parts = self.output_folder.split(os.sep)
        for i, part in enumerate(folder_parts):
            if self.TOOL_NAME in part:
                if i > 0:
                    return folder_parts[i - 1]
        raise ValueError(f"Could not extract pipeline name from output folder: {self.output_folder}")

    def _compute_log_file_path(self) -> str:
        """
        Compute log file path from output folder naming pattern.

        Log files are stored in the Logs folder with pattern NNN_ToolName.log

        Returns:
            Path to log file

        Raises:
            ValueError: If folder naming pattern is invalid
        """
        folder_name = os.path.basename(self.output_folder)
        pipeline_folder = os.path.dirname(self.output_folder)
        logs_folder = os.path.join(pipeline_folder, "Logs")

        if '_' in folder_name and folder_name.split('_')[0].isdigit():
            index = folder_name.split('_')[0]
            tool_name = folder_name.split('_', 1)[1]
            return os.path.join(logs_folder, f"{index}_{tool_name}.log")
        raise ValueError(f"Invalid output folder naming pattern: {folder_name}. Expected 'NNN_ToolName' format.")

    def get_effective_job_name(self) -> str:
        """
        Get the effective job name, preferring pipeline job name over tool name.

        Returns:
            Job name from pipeline if available, otherwise tool job name, or None if neither available
        """
        if self.pipeline and hasattr(self.pipeline, 'job_name') and self.pipeline.job_name:
            return self.pipeline.job_name
        return self.job_name if self.job_name else None
    
    @abstractmethod
    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """
        Configure input sources from pipeline context and dependencies.
        
        Args:
            pipeline_folders: Dictionary of pipeline folder paths
        """
        pass
    
    @abstractmethod
    def generate_script(self, script_path: str) -> str:
        """
        Generate bash script for tool execution.
        
        Args:
            script_path: Path where script should be written
            
        Returns:
            Script content as string
        """
        pass
    
    @abstractmethod
    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after execution.
        
        Returns:
            Dictionary mapping output type to file paths
        """
        pass
    
    def get_expected_output_paths(self) -> Dict[str, List[str]]:
        """
        Get expected output file paths without validating existence.
        
        Default implementation just calls get_output_files().
        Override if file existence validation is problematic.
        
        Returns:
            Dictionary mapping output type to expected file paths
        """
        return self.get_output_files()
    
    def set_pipeline_context(self, pipeline_ref, execution_order: int, output_folder: str, suffix: str = ""):
        """Set context when added to pipeline."""
        self.pipeline_ref = pipeline_ref
        self.execution_order = execution_order
        self.output_folder = output_folder
        self.suffix = suffix
        self.configured = True
    
    def resolve_dependency_outputs(self, dependency):
        """
        Resolve outputs from a dependency tool.
        
        Args:
            dependency: Another BaseConfig instance or output specification
            
        Returns:
            Dictionary of available outputs from dependency
        """
        if hasattr(dependency, 'get_output_files'):
            return dependency.get_output_files()
        elif isinstance(dependency, dict):
            return dependency
        else:
            raise ValueError(f"Cannot resolve outputs from dependency: {dependency}")
    
    def get_resource_requirements(self) -> Dict[str, str]:
        """Get SLURM resource requirements for this tool."""
        return self.resources.copy()
    
    def get_config_display(self) -> List[str]:
        """
        Get configuration display lines for pipeline config output.
        Override in subclasses to show tool-specific parameters.
        
        Returns:
            List of configuration strings for display
        """
        config_lines = []
        
        # Add common parameters
        if hasattr(self, 'job_name') and self.job_name:
            config_lines.append(f"Job: {self.job_name}")
        
        # Subclasses should override to add tool-specific parameters
        return config_lines
    
    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration to dictionary."""
        return {
            'tool_name': self.tool_name,
            'job_name': self.job_name,
            'environments': self.environments,
            'resources': self.resources,
            'dependencies': len(self.dependencies),
            'execution_order': self.execution_order
        }
    
    def save_config(self, config_file: str):
        """Save configuration to JSON file."""
        with open(config_file, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)
    
    def generate_completion_check_header(self) -> str:
        """
        Generate bash script header that checks for completion status.

        Returns:
            Bash script content for checking completion
        """
        # Get step number from output folder (e.g., "1_RFdiffusion")
        folder_name = os.path.basename(self.output_folder)
        if '_' in folder_name and folder_name.split('_')[0].isdigit():
            step_number = folder_name.split('_')[0]
            completed_file = f"{step_number}_{self.TOOL_NAME}_COMPLETED"
        else:
            completed_file = f"{self.TOOL_NAME}_COMPLETED"

        parent_dir = os.path.dirname(self.output_folder)

        return f"""# Check if already completed
if [ -f "{parent_dir}/{completed_file}" ]; then
    echo "{self.TOOL_NAME} already completed, skipping..."
    exit 0
fi

"""

    def generate_completion_check_footer(self) -> str:
        """
        Generate bash script footer that checks outputs and creates status files.
        
        Returns:
            Bash script content for final completion check
        """
        # Get expected outputs as JSON string
        expected_outputs = self.get_expected_output_paths()
        
        # Convert TableInfo objects to dictionaries for JSON serialization
        def make_json_safe(obj, seen=None):
            """Recursively convert TableInfo objects to dictionaries."""
            if seen is None:
                seen = set()

            # Avoid infinite recursion from circular references
            obj_id = id(obj)
            if obj_id in seen:
                return f"<circular reference to {type(obj).__name__}>"

            if isinstance(obj, TableInfo):
                # Explicitly convert TableInfo to plain dict, avoiding all attributes
                return {
                    "name": obj.name,
                    "path": obj.path,
                    "columns": list(obj.columns) if obj.columns else [],
                    "description": obj.description,
                    "count": obj.count
                }
            elif isinstance(obj, dict):
                seen.add(obj_id)
                result = {k: make_json_safe(v, seen) for k, v in obj.items()}
                seen.remove(obj_id)
                return result
            elif isinstance(obj, (list, tuple)):
                seen.add(obj_id)
                result = [make_json_safe(item, seen) for item in obj]
                seen.remove(obj_id)
                return result
            elif hasattr(obj, '__dict__'):
                # Handle custom objects - check for TableInfo first
                if hasattr(obj, 'name') and hasattr(obj, 'path') and hasattr(obj, 'columns'):
                    # This looks like a TableInfo object that wasn't caught above
                    return {
                        "name": getattr(obj, 'name', ''),
                        "path": getattr(obj, 'path', ''),
                        "columns": list(getattr(obj, 'columns', [])),
                        "description": getattr(obj, 'description', ''),
                        "count": getattr(obj, 'count', None)
                    }
                # Handle DataStream objects - extract files and map_table for completion checking
                elif hasattr(obj, 'files') and hasattr(obj, 'ids') and hasattr(obj, 'map_table'):
                    # Return list of all paths that need to exist: files + map_table
                    paths = list(getattr(obj, 'files', []))
                    map_table = getattr(obj, 'map_table', '')
                    if map_table:
                        paths.append(map_table)
                    return paths
                else:
                    # Convert other custom objects to string
                    return str(obj)
            else:
                return obj

        json_safe_outputs = make_json_safe(expected_outputs)

        # Debug: try to serialize and catch the exact error
        try:
            expected_outputs_json = json.dumps(json_safe_outputs).replace('"', '\\"')
        except TypeError as e:
            print(f"JSON serialization error: {e}")
            print(f"Expected outputs keys: {list(expected_outputs.keys()) if isinstance(expected_outputs, dict) else 'Not a dict'}")
            print(f"JSON safe outputs keys: {list(json_safe_outputs.keys()) if isinstance(json_safe_outputs, dict) else 'Not a dict'}")

            # Try to find the problematic object
            for key, value in json_safe_outputs.items():
                try:
                    json.dumps(value)
                    print(f"✓ Key '{key}' serializes OK")
                except Exception as sub_e:
                    print(f"✗ Key '{key}' failed: {sub_e}")
                    print(f"  Type: {type(value)}")
                    if isinstance(value, dict):
                        for subkey, subvalue in value.items():
                            try:
                                json.dumps(subvalue)
                                print(f"  ✓ Subkey '{subkey}' OK")
                            except Exception as subsub_e:
                                print(f"  ✗ Subkey '{subkey}' failed: {subsub_e} (type: {type(subvalue)})")
                                # Drill down further if it's a dict
                                if isinstance(subvalue, dict):
                                    print(f"    Contents: {list(subvalue.keys())}")
                                    for k, v in subvalue.items():
                                        print(f"    '{k}': {type(v)} = {repr(v)[:100]}...")
                                        if hasattr(v, '__dict__'):
                                            print(f"      __dict__ keys: {list(v.__dict__.keys()) if hasattr(v, '__dict__') else 'None'}")

            # Fallback: convert everything to strings
            def stringify_all(obj):
                if isinstance(obj, dict):
                    return {k: stringify_all(v) for k, v in obj.items()}
                elif isinstance(obj, (list, tuple)):
                    return [stringify_all(item) for item in obj]
                else:
                    return str(obj)

            print("Using fallback string conversion for all objects")
            json_safe_outputs = stringify_all(expected_outputs)
            expected_outputs_json = json.dumps(json_safe_outputs).replace('"', '\\"')
        
        pipe_check_completion = os.path.join(
            self.folders.get("HelpScripts", "HelpScripts"),
            "pipe_check_completion.py"
        )

        # Write expected outputs to a temporary JSON file to avoid "Argument list too long" error
        expected_outputs_file = os.path.join(self.output_folder, ".expected_outputs.json")

        return f"""
# Check completion and create status files
echo "Checking outputs and creating completion status..."

# Write expected outputs to JSON file to avoid argument list length limits
cat > "{expected_outputs_file}" << 'EXPECTED_OUTPUTS_EOF'
{json.dumps(json_safe_outputs, indent=2)}
EXPECTED_OUTPUTS_EOF

python {pipe_check_completion} "{self.output_folder}" "{self.TOOL_NAME}" "{expected_outputs_file}"

if [ $? -eq 0 ]; then
    echo "{self.TOOL_NAME} completed successfully"
else
    echo "{self.TOOL_NAME} failed - some outputs missing"
    exit 1
fi
"""
    
    def __str__(self) -> str:
        """String representation of configuration."""
        return f"{self.TOOL_NAME}(envs={self.environments}, order={self.execution_order})"
    
    def __repr__(self) -> str:
        return self.__str__()
    
    def resolve_table_reference(self, reference) -> str:
        """
        Resolve table column references to actual values.

        Supports two formats:
        - String: 'table_name.column_name' (e.g., "structures.fixed", "selections.within")
        - Tuple: (table_object, "column_name") (e.g., (tool.tables.selections, "within"))

        This provides a general way to reference any column from any input table across all tools.

        Args:
            reference: Reference string or tuple (e.g., "structures.fixed" or (tool.tables.selections, "within"))

        Returns:
            Resolved value from table or original reference if not a table reference
        """
        if not reference:
            return reference

        # Handle tuple format: (table_object, "column_name")
        if isinstance(reference, tuple) and len(reference) == 2:
            table_object, column_name = reference

            # Get the table path from the object
            if hasattr(table_object, 'path'):
                table_path = table_object.path
            else:
                raise ValueError(f"Invalid table object in tuple reference: {table_object}")

            # Read and resolve the column value
            return self._resolve_column_from_table(table_path, column_name)

        # Handle string format: "table.column"
        if isinstance(reference, str) and "." in reference:
            parts = reference.split(".")
            if len(parts) == 2:
                # Table reference format: "structures.fixed"
                table_name = parts[0]
                column_name = parts[1]

                # Get the table path
                table_path = self._find_table_path(table_name)

                if not table_path:
                    raise ValueError(f"Table '{table_name}' not found in input tables")

                # Read and resolve the column value
                return self._resolve_column_from_table(table_path, column_name)
            else:
                # Not a table reference
                return reference
        else:
            # Not a table reference
            return reference

    def _resolve_column_from_table(self, table_path: str, column_name: str) -> str:
        """
        Resolve a column value from a table file.

        Args:
            table_path: Path to the table CSV file
            column_name: Name of the column to resolve

        Returns:
            Placeholder string for script generation
        """
        # For script generation, return a placeholder that indicates this is a table reference
        # The actual resolution will happen in the script generation where we have access to pandas
        return f"DATASHEET_REFERENCE:{table_path}:{column_name}"
    
    def _find_table_path(self, table_name: str) -> Optional[str]:
        """
        Find the path to a named table from various input sources.
        
        Args:
            table_name: Name of the table to find
            
        Returns:
            Path to the table file, or None if not found
        """
        # Check StandardizedOutput format
        if hasattr(self, 'standardized_input') and self.standardized_input and hasattr(self.standardized_input, 'tables'):
            if hasattr(self.standardized_input.tables, '_tables') and table_name in self.standardized_input.tables._tables:
                return self.standardized_input.tables._tables[table_name].path
        
        # Check dictionary format
        if hasattr(self, 'input_tables') and isinstance(self.input_tables, dict) and table_name in self.input_tables:
            return self.input_tables[table_name]["path"]
        
        # Check TableContainer format  
        if hasattr(self, 'input_tables') and hasattr(self.input_tables, '_tables') and table_name in self.input_tables._tables:
            return self.input_tables._tables[table_name].path
        
        return None
    
    def validate_table_reference(self, reference):
        """Validate table reference format."""
        if not reference:
            return

        # Handle tuple format: (table_object, "column_name")
        if isinstance(reference, tuple):
            if len(reference) == 2:
                table_object, column_name = reference
                if hasattr(table_object, 'path') and isinstance(column_name, str):
                    # Valid tuple format
                    return
                else:
                    raise ValueError(f"Invalid tuple format: expected (table_object, 'column_name'), got {reference}")
            else:
                raise ValueError(f"Invalid tuple format: expected (table_object, 'column_name'), got {reference}")

        # Handle string format: table.column
        if isinstance(reference, str) and "." in reference:
            parts = reference.split(".")
            if len(parts) == 2:
                # Valid table reference format
                return
            else:
                # Multiple dots - not a valid table reference
                return

        # No dots - regular value, not a table reference

    def _get_upstream_missing_table_path(self, *input_sources) -> Optional[str]:
        """
        Get the path to the missing table from upstream tool outputs.

        This should be called during get_output_files() to check if there's an upstream
        missing table that needs to be propagated.

        Args:
            *input_sources: Variable number of StandardizedOutput or ToolOutput objects

        Returns:
            Path to upstream missing.csv, or None if no missing table exists
        """
        for input_source in input_sources:
            if input_source is None:
                continue

            # Try to find missing table from input source
            if hasattr(input_source, 'tables'):
                tables = input_source.tables
                if hasattr(tables, '_tables') and 'missing' in tables._tables:
                    missing_info = tables._tables['missing']
                    return missing_info.path if hasattr(missing_info, 'path') else str(missing_info)
                elif isinstance(tables, dict) and 'missing' in tables:
                    missing_info = tables['missing']
                    if isinstance(missing_info, dict) and 'path' in missing_info:
                        return missing_info['path']
                    elif hasattr(missing_info, 'path'):
                        return missing_info.path
                    else:
                        return str(missing_info)

        return None


class TableInfo:
    """Information about a table including name, path, and expected columns."""

    def __init__(self, name: str, path: str, columns: List[str] = None,
                 description: str = "", count: int = 0):
        self.name = name
        self.path = path
        self.columns = columns or []
        self.description = description
        self.count = count

        # Set column attributes for IDE autocompletion
        # Avoid overwriting existing attributes (name, path, columns, description, count)
        reserved_attributes = {'name', 'path', 'columns', 'description', 'count'}
        for column in self.columns:
            if column not in reserved_attributes:
                setattr(self, column, self._create_column_reference(column))

    def _create_column_reference(self, column_name: str):
        """Create a column reference tuple for table access."""
        return (self, column_name)

    def __getattr__(self, column_name: str):
        """
        Handle access to any column name.
        This enables both IDE autocompletion for known columns and dynamic access.
        """
        # Return column reference tuple for any column name
        return self._create_column_reference(column_name)

    def __str__(self) -> str:
        if self.columns:
            # Show all columns
            col_display = ', '.join(self.columns)
            return f"${self.name} ({col_display})"
        return f"${self.name}"

    def __repr__(self) -> str:
        return f"TableInfo(name='{self.name}', path='{self.path}', columns={self.columns}, count={self.count})"

    def to_dict(self) -> Dict[str, Any]:
        """Convert TableInfo to dictionary for JSON serialization."""
        # Only return core data, not the dynamically set column attributes which contain circular references
        return {
            "name": self.name,
            "path": self.path,
            "columns": self.columns.copy() if self.columns else [],
            "description": self.description,
            "count": self.count
        }


class TableContainer:
    """Container for named tables with dot-notation access."""
    
    def __init__(self, tables: Dict[str, TableInfo]):
        self._tables = tables

        # Set attributes for dot notation access to TableInfo objects (not just paths)
        # Note: For IDE autocompletion, attributes should also be explicitly set in calling code
        for name, info in tables.items():
            setattr(self, name, info)
    
    def __getitem__(self, key: str) -> str:
        """Get table path by name with legacy 'main' support."""
        if key in self._tables:
            return self._tables[key].path

        raise KeyError(f"No table named '{key}' in tables")
    
    def __getattr__(self, name: str) -> TableInfo:
        """Get TableInfo object by name via dot notation."""
        if name in self._tables:
            return self._tables[name]
        raise AttributeError(f"No table named '{name}'")
    
    def keys(self):
        """Get all table names."""
        return self._tables.keys()

    def info(self, name: str) -> TableInfo:
        """Get TableInfo object by name."""
        if name in self._tables:
            return self._tables[name]
        raise ValueError(f"No table named '{name}'")
    
    def items(self):
        """Get all name, path pairs."""
        return [(name, info.path) for name, info in self._tables.items()]
    
    def info(self, name: str) -> TableInfo:
        """Get full TableInfo object."""
        return self._tables.get(name)
    
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
        for name, info in self._tables.items():
            # Format: table name with columns, then indented file path
            # Extract just the filename from the full path
            filename = info.path.split('/')[-1]
            path_display = f"<output_folder>/{filename}"
            lines.append(f"    {str(info)}:")
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
        """Get DataStream by name via dot notation."""
        if '_streams' in self.__dict__ and name in self._streams:
            return self._streams[name]
        raise AttributeError(f"No stream named '{name}'")

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


class StandardizedOutput:
    """
    Provides dot-notation access to standardized output keys.

    DataStreams are accessed via the streams container:

        for struct_id, pdb_path in output.streams.structures:
            print(f"Processing {struct_id}: {pdb_path}")

        print(f"Generated {len(output.streams.structures)} structures")
    """

    def __init__(self, output_files: Dict[str, Any]):
        """Initialize with output files dictionary."""
        self._data = output_files.copy()

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
    
    def _process_tables(self, tables_raw: Any) -> 'TableContainer':
        """Process tables into TableContainer with named access."""
        if isinstance(tables_raw, dict):
            # Already in named format
            table_infos = {}
            for name, info in tables_raw.items():
                if isinstance(info, dict):
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
                        count=info.get('count', 0)
                    )
                else:
                    # Could be a TableInfo object or just a path
                    if isinstance(info, TableInfo):
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
    
    def __getitem__(self, key: str):
        """Dictionary-style access: output['structures']"""
        return self._data.get(key, [])
    
    def __iter__(self):
        """Allow iteration over keys."""
        return iter(self._data)
    
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
            return table_info.path if hasattr(table_info, 'path') else str(table_info)
        
        # For "main", try common table types in order of preference
        if name == "main":
            for fallback_name in ['sequences', 'structures', 'compounds']:
                if hasattr(self.tables, fallback_name):
                    table_info = getattr(self.tables, fallback_name)
                    return table_info.path if hasattr(table_info, 'path') else str(table_info)
        
        return ""
    
    def pretty(self) -> str:
        """Pretty formatted representation of the output."""
        from .datastream import DataStream

        lines = []

        # Helper function to make paths relative to output_folder
        def make_relative_path(file_path: str) -> str:
            if self.output_folder and file_path.startswith(self.output_folder):
                relative_path = os.path.relpath(file_path, self.output_folder)
                return f"<output_folder>/{relative_path}"
            return file_path

        def format_datastream(name: str, ds: DataStream) -> List[str]:
            """Format a DataStream for display."""
            if not ds or len(ds) == 0:
                return []

            result = [f"{name}: ({ds.format}, {len(ds)} items)"]

            # Show (id, file) pairs
            items = list(ds)
            if len(items) <= 6:
                for item_id, item_file in items:
                    if item_file:
                        rel_path = make_relative_path(item_file)
                        result.append(f"    – {item_id}: '{rel_path}'")
                    else:
                        result.append(f"    – {item_id}")
            else:
                # Truncated: first 3, ..., last 3
                for item_id, item_file in items[:3]:
                    if item_file:
                        rel_path = make_relative_path(item_file)
                        result.append(f"    – {item_id}: '{rel_path}'")
                    else:
                        result.append(f"    – {item_id}")
                result.append(f"    – ... ({len(items) - 6} more) ...")
                for item_id, item_file in items[-3:]:
                    if item_file:
                        rel_path = make_relative_path(item_file)
                        result.append(f"    – {item_id}: '{rel_path}'")
                    else:
                        result.append(f"    – {item_id}")

            return result

        # Main DataStream attributes
        for attr_name in ('structures', 'sequences', 'compounds', 'msas'):
            ds = getattr(self, attr_name, None)
            if isinstance(ds, DataStream) and len(ds) > 0:
                lines.extend(format_datastream(attr_name, ds))

        # Tables
        if 'tables' in self._data and hasattr(self.tables, '_tables'):
            lines.append("tables:")
            for name, info in self.tables._tables.items():
                col_display = ', '.join(info.columns[:]) if info.columns else ''
                lines.append(f"    ${name} ({col_display}):")
                relative_path = make_relative_path(info.path)
                lines.append(f"        – '{relative_path}'")

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
    
    def __init__(self, config: BaseConfig):
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

    def __str__(self) -> str:
        return f"ToolOutput({self.tool_type}, {len(self._output_files)} output types)"
    
    def __repr__(self) -> str:
        return self.__str__()