# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

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

    @classmethod
    def _install_script(cls, folders: Dict[str, str], env_manager: str = "mamba",
                        force_reinstall: bool = False, **kwargs) -> Optional[str]:
        """Override in subclasses to provide installation bash commands.

        Args:
            folders: Resolved pipeline folder paths from config.yaml
                     (e.g., folders["data"], folders["RFdiffusion"], etc.)
            env_manager: Environment manager command from config.yaml
                         (e.g., "mamba" or "conda")
            force_reinstall: If False, skip installation when the tool is
                             already installed (e.g. repo already cloned,
                             environment already created). If True, always run
                             the full installation.
            **kwargs: Additional tool-specific arguments.

        Returns:
            Bash script content for installing this tool, or None if not defined.
        """
        return None

    @classmethod
    def install(cls, force_reinstall: bool = False, **kwargs):
        """Add an installation step for this tool to the active pipeline.

        Must be called within a Pipeline context. Tools must override
        _install_script() to provide installation bash commands.

        Usage:
            with Pipeline(...):
                Resources(...)
                RFdiffusion.install()                    # skip if already installed
                RFdiffusion.install(force_reinstall=True) # always reinstall
                rfd = RFdiffusion(...)

        Args:
            force_reinstall: If False (default), skip installation when the tool
                is already installed. If True, always run the full installation.
            **kwargs: Additional tool-specific arguments forwarded to
                _install_script().

        Returns:
            StandardizedOutput (via auto-registration, just like tool instantiation)

        Raises:
            NotImplementedError: If the tool hasn't defined _install_script()
        """
        # Verify _install_script is overridden (call with empty dict/default to check)
        if cls._install_script({}, "mamba") is None:
            raise NotImplementedError(
                f"{cls.__name__} does not define installation steps. "
                f"Override _install_script() to provide bash commands."
            )
        return _Installer(parent_tool_cls=cls, force_reinstall=force_reinstall,
                          install_kwargs=kwargs)

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

        In pip mode (e.g. Google Colab), no environments are needed since
        everything runs in the pre-existing Python environment.
        """
        config_manager = ConfigManager()
        env_config = config_manager.get_environment(self.TOOL_NAME)

        if env_config is None and config_manager.get_env_manager() == "pip":
            # pip mode with no environments configured: nothing to activate
            self.environments = []
            return

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
        Generate bash script snippet to activate an environment.

        For pip mode (e.g. Google Colab), returns a no-op comment since
        everything runs in the pre-existing Python environment.

        Args:
            index: Index of the environment to activate (default: 0, the first/primary environment)
            name: Optional explicit environment name (overrides index-based lookup)

        Returns:
            Bash script content for activating the environment with diagnostics

        Raises:
            IndexError: If index is out of range for configured environments
        """
        config_manager = ConfigManager()

        # pip mode: no environment activation needed
        if config_manager.get_env_manager() == "pip":
            return "# pip mode: no environment activation needed\n\n"

        if name is not None:
            env_name = name
        else:
            if index >= len(self.environments):
                raise IndexError(
                    f"Environment index {index} out of range for {self.TOOL_NAME}. "
                    f"Available environments: {self.environments}"
                )
            env_name = self.environments[index]

        shell_hook = config_manager.get_shell_hook_command()
        activate_cmd = config_manager.get_activate_command(env_name)

        return f"""# Activate environment: {env_name}
echo "=== Activating Environment ==="
echo "Requested: {env_name}"
{shell_hook}
{activate_cmd}
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

    def _repr_notebook_html(self, output: 'StandardizedOutput') -> str:
        """
        Generate custom HTML visualization for notebook display.

        Override in subclasses to provide tool-specific visualization.
        The default implementation renders an interactive py3Dmol viewer
        if the tool outputs a non-empty structures stream.

        Args:
            output: The StandardizedOutput being displayed

        Returns:
            HTML string for notebook display, or empty string if nothing to show
        """
        from .datastream import DataStream

        structures_ds = output.streams.get("structures")
        if not isinstance(structures_ds, DataStream) or len(structures_ds) == 0:
            return ""

        return self._build_py3dmol_html(structures_ds)

    @staticmethod
    def _build_py3dmol_html(structures_ds, max_structures: int = 50) -> str:
        """
        Build py3Dmol interactive 3D viewer HTML with prev/next navigation.

        Displays one structure at a time with left/right buttons and an ID label.

        Args:
            structures_ds: DataStream containing structure files (pdb/cif)
            max_structures: Maximum number of structures to load

        Returns:
            HTML string with the navigable 3D viewer
        """
        import json as json_module

        if structures_ds.format not in ("pdb", "cif"):
            return ""

        pdb_data = []
        for struct_id, file_path in structures_ds:
            if file_path and os.path.isfile(file_path):
                try:
                    with open(file_path, "r") as f:
                        pdb_data.append((struct_id, f.read()))
                except Exception:
                    pass
            if len(pdb_data) >= max_structures:
                break

        if not pdb_data:
            return ""

        fmt = "pdb" if structures_ds.format == "pdb" else "cif"

        # Unique ID for this viewer instance to avoid conflicts with multiple viewers
        import random
        viewer_id = f"bp3d_{random.randint(100000, 999999)}"

        # Serialize structure data for JavaScript
        struct_ids_json = json_module.dumps([sid for sid, _ in pdb_data])
        struct_data_json = json_module.dumps([content for _, content in pdb_data])

        truncated = len(structures_ds) > max_structures
        total_label = f"{len(pdb_data)} structure{'s' if len(pdb_data) != 1 else ''}"
        if truncated:
            total_label += f" (of {len(structures_ds)} total)"

        colors = [
            "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
            "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
        ]
        colors_json = json_module.dumps(colors)

        html = f"""
<script src="https://cdn.jsdelivr.net/npm/3dmol@2.5.2/build/3Dmol-min.js"></script>
<div style="margin-top: 12px;">
  <strong>3D Structure Viewer</strong> ({total_label})
</div>
<div id="{viewer_id}_container" style="position: relative; width: 800px;">
  <div id="{viewer_id}_viewer" style="width: 800px; height: 500px; position: relative;"></div>
  <div style="display: flex; align-items: center; justify-content: center; gap: 12px; margin-top: 6px; font-family: monospace;">
    <button id="{viewer_id}_prev" onclick="{viewer_id}_navigate(-1)"
            style="padding: 4px 14px; font-size: 1.1em; cursor: pointer; border: 1px solid #ccc; border-radius: 4px; background: #f5f5f5;">&#9664;</button>
    <span id="{viewer_id}_label" style="min-width: 200px; text-align: center; font-size: 0.95em;"></span>
    <button id="{viewer_id}_next" onclick="{viewer_id}_navigate(1)"
            style="padding: 4px 14px; font-size: 1.1em; cursor: pointer; border: 1px solid #ccc; border-radius: 4px; background: #f5f5f5;">&#9654;</button>
  </div>
</div>
<script>
(function() {{
  var ids = {struct_ids_json};
  var data = {struct_data_json};
  var fmt = "{fmt}";
  var colors = {colors_json};
  var idx = 0;
  var viewer = null;

  function initViewer() {{
    if (typeof $3Dmol === "undefined") {{
      setTimeout(initViewer, 200);
      return;
    }}
    var el = document.getElementById("{viewer_id}_viewer");
    viewer = $3Dmol.createViewer(el, {{backgroundColor: "white"}});
    showStructure(0);
  }}

  function showStructure(i) {{
    if (!viewer) return;
    idx = i;
    if (idx < 0) idx = ids.length - 1;
    if (idx >= ids.length) idx = 0;
    viewer.removeAllModels();
    viewer.addModel(data[idx], fmt);
    var color = colors[idx % colors.length];
    viewer.setStyle({{}}, {{"cartoon": {{"color": color}}}});
    viewer.zoomTo();
    viewer.render();
    document.getElementById("{viewer_id}_label").innerHTML =
      '<span style="display:inline-block;width:12px;height:12px;background:' + color +
      ';border-radius:2px;vertical-align:middle;margin-right:6px;"></span>' +
      ids[idx] + '  <span style="color:#888;">(' + (idx+1) + '/' + ids.length + ')</span>';
  }}

  window.{viewer_id}_navigate = function(delta) {{
    showStructure(idx + delta);
  }};

  initViewer();
}})();
</script>
"""
        return html

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
                return obj.to_dict()
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
                if hasattr(obj, 'info') and hasattr(obj, 'to_dict') and isinstance(obj, TableInfo):
                    # This looks like a TableInfo object that wasn't caught above
                    return obj.to_dict()
                # Handle DataStream objects - use to_dict() for consistent serialization
                elif hasattr(obj, 'files') and hasattr(obj, 'ids') and hasattr(obj, 'map_table'):
                    return obj.to_dict()
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

        # Write expected outputs to JSON file at pipeline time (not SLURM time)
        # Wrap in metadata envelope so LoadOutput can load directly from the tool folder
        expected_outputs_wrapped = {
            "tool_name": self.TOOL_NAME,
            "tool_class": self.__class__.__name__,
            "output_structure": json_safe_outputs,
        }
        expected_outputs_file = os.path.join(self.output_folder, ".expected_outputs.json")
        os.makedirs(self.output_folder, exist_ok=True)
        with open(expected_outputs_file, 'w') as f:
            json.dump(expected_outputs_wrapped, f, indent=2)

        return f"""
# Check completion and create status files
echo "Checking outputs and creating completion status..."

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
            if hasattr(table_object, 'info'):
                table_path = table_object.info.path
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
                return self.standardized_input.tables._tables[table_name].info.path
        
        # Check dictionary format
        if hasattr(self, 'input_tables') and isinstance(self.input_tables, dict) and table_name in self.input_tables:
            return self.input_tables[table_name]["path"]
        
        # Check TableContainer format  
        if hasattr(self, 'input_tables') and hasattr(self.input_tables, '_tables') and table_name in self.input_tables._tables:
            return self.input_tables._tables[table_name].info.path
        
        return None
    
    def validate_table_reference(self, reference):
        """Validate table reference format."""
        if not reference:
            return

        # Handle tuple format: (table_object, "column_name")
        if isinstance(reference, tuple):
            if len(reference) == 2:
                table_object, column_name = reference
                if hasattr(table_object, 'info') and isinstance(column_name, str):
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
                    return missing_info.info.path if hasattr(missing_info, 'info') else str(missing_info)
                elif isinstance(tables, dict) and 'missing' in tables:
                    missing_info = tables['missing']
                    if isinstance(missing_info, dict) and 'path' in missing_info:
                        return missing_info['path']
                    elif hasattr(missing_info, 'info'):
                        return missing_info.info.path
                    else:
                        return str(missing_info)

        return None


class _Installer(BaseConfig):
    """Generic installer step added to the pipeline by BaseConfig.install()."""

    TOOL_NAME = "install"  # Overridden dynamically in __init__

    def __init__(self, parent_tool_cls: type, force_reinstall: bool = False,
                 install_kwargs: Dict[str, Any] = None, **kwargs):
        self._parent_tool_cls = parent_tool_cls
        self.TOOL_NAME = f"{parent_tool_cls.TOOL_NAME}_installation"
        self._parent_tool_name = parent_tool_cls.TOOL_NAME
        self._force_reinstall = force_reinstall
        self._install_kwargs = install_kwargs or {}
        super().__init__(**kwargs)

    def _load_environments(self):
        """Use biopipelines environment since install scripts create their own environments."""
        self.environments = ["biopipelines"]

    def validate_params(self):
        pass

    def configure_inputs(self, pipeline_folders):
        self.folders = pipeline_folders

    def generate_script(self, script_path: str) -> str:
        config_manager = ConfigManager()
        env_manager = config_manager.get_env_manager()
        install_commands = self._parent_tool_cls._install_script(
            self.folders, env_manager,
            force_reinstall=self._force_reinstall, **self._install_kwargs)
        script = "#!/bin/bash\n"
        script += f"# Installation: {self._parent_tool_name}\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += install_commands + "\n"
        script += self.generate_completion_check_footer()
        return script

    def get_output_files(self):
        from .datastream import DataStream
        return {
            "structures": DataStream.empty("structures", "pdb"),
            "sequences": DataStream.empty("sequences", "fasta"),
            "compounds": DataStream.empty("compounds", "sdf"),
            "tables": {},
            "output_folder": self.output_folder
        }


class TableMetadata:
    """Holds table metadata (name, path, columns, description, count)."""

    def __init__(self, name: str, path: str, columns: List[str] = None,
                 description: str = "", count: int = 0):
        self.name = name
        self.path = path
        self.columns = columns or []
        self.description = description
        self.count = count

    def __repr__(self) -> str:
        return f"TableMetadata(name='{self.name}', path='{self.path}', columns={self.columns}, count={self.count})"


class TableInfo:
    """Information about a table including name, path, and expected columns.

    Metadata is accessed via the .info property:
        table.info.path, table.info.columns, table.info.count, etc.

    All other attribute access returns column reference tuples (TableInfo, column_name).
    """

    def __init__(self, name: str, path: str, columns: List[str] = None,
                 description: str = "", count: int = 0):
        self._info = TableMetadata(name, path, columns, description, count)

        # Set column attributes for IDE autocompletion
        for column in self._info.columns:
            setattr(self, column, self._create_column_reference(column))

    @property
    def info(self) -> TableMetadata:
        """Access table metadata (name, path, columns, description, count)."""
        return self._info

    def _create_column_reference(self, column_name: str):
        """Create a column reference tuple for table access."""
        return (self, column_name)

    def __getattr__(self, column_name: str):
        """Return column reference tuple for any attribute access."""
        return self._create_column_reference(column_name)

    def __str__(self) -> str:
        if self._info.columns:
            col_display = ', '.join(self._info.columns)
            return f"${self._info.name} ({col_display})"
        return f"${self._info.name}"

    def __repr__(self) -> str:
        return f"TableInfo(name='{self._info.name}', path='{self._info.path}', columns={self._info.columns}, count={self._info.count})"

    def to_dict(self) -> Dict[str, Any]:
        """Convert TableInfo to dictionary for JSON serialization."""
        return {
            "name": self._info.name,
            "path": self._info.path,
            "columns": self._info.columns.copy() if self._info.columns else [],
            "description": self._info.description,
            "count": self._info.count
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
            return self._tables[key].info.path

        raise KeyError(f"No table named '{key}' in tables")
    
    def __getattr__(self, name: str) -> TableInfo:
        """Get TableInfo object by name via dot notation."""
        if name in self._tables:
            return self._tables[name]
        raise AttributeError(f"No table named '{name}'")
    
    def keys(self):
        """Get all table names."""
        return self._tables.keys()

    def items(self):
        """Get all name, path pairs."""
        return [(name, info.info.path) for name, info in self._tables.items()]
    
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
            filename = info.info.path.split('/')[-1]
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
            result.append(f"    files_contain_wildcards: {ds.files_contain_wildcards}")
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
                items = list(ds)
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
                lines.append(f"    count: {t_meta.count}")
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

    # Image formats recognized for inline display in notebooks
    _IMAGE_FORMATS = {"png", "jpg", "jpeg", "svg", "gif", "bmp", "tiff", "webp"}

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

    def _repr_html_(self) -> Optional[str]:
        """
        Rich HTML display for Jupyter notebooks.

        Renders streams as HTML tables, displays images inline,
        and delegates tool-specific visualization (e.g. 3D viewer)
        to the tool's _repr_notebook_html() method.
        """
        if not self._is_notebook():
            return None

        from .datastream import DataStream
        import base64
        import html as html_module
        import pandas as pd

        html_parts = []
        css = """<style>
.bp-table { border-collapse: collapse; font-family: monospace; font-size: 0.9em; margin: 4px 0 12px 0; }
.bp-table th { background: #f0f0f0; padding: 4px 10px; border: 1px solid #ddd; text-align: left; }
.bp-table td { padding: 4px 10px; border: 1px solid #ddd; }
.bp-table tr:nth-child(even) { background: #fafafa; }
.bp-table .bp-ellipsis td { color: #888; font-style: italic; text-align: center; }
.bp-section { margin-top: 8px; font-family: monospace; }
.bp-section-title { font-weight: bold; margin-bottom: 4px; }
</style>"""

        # Helper: make path relative to output_folder
        def rel(path: str) -> str:
            if self.output_folder and path.startswith(self.output_folder):
                return "<output_folder>/" + os.path.relpath(path, self.output_folder)
            return path

        # Helper: render a DataFrame as an HTML table with 2+...+2 truncation
        def _render_dataframe(df, html_out):
            columns = list(df.columns)
            html_out.append('<table class="bp-table"><tr>')
            for col in columns:
                html_out.append(f'<th>{html_module.escape(str(col))}</th>')
            html_out.append('</tr>')

            n_rows = len(df)
            if n_rows <= 4:
                display_rows = list(range(n_rows))
                ellipsis_after = None
            else:
                display_rows = list(range(2)) + list(range(n_rows - 2, n_rows))
                ellipsis_after = 2

            for row_i, idx in enumerate(display_rows):
                if ellipsis_after is not None and row_i == ellipsis_after:
                    html_out.append(
                        f'<tr class="bp-ellipsis"><td colspan="{len(columns)}">'
                        f'... {n_rows - 4} more ...</td></tr>'
                    )
                row = df.iloc[idx]
                html_out.append('<tr>')
                for col in columns:
                    val = str(row[col]) if pd.notna(row[col]) else ''
                    display_val = val if len(val) <= 60 else val[:57] + '...'
                    html_out.append(f'<td>{html_module.escape(display_val)}</td>')
                html_out.append('</tr>')
            html_out.append('</table>')

        # --- Streams ---
        # Collect non-empty, non-image streams
        active_streams = [
            (sn, s) for sn, s in self.streams.items()
            if isinstance(s, DataStream) and len(s) > 0
            and s.format.lower() not in self._IMAGE_FORMATS
        ]

        if active_streams:
            # Streams summary table
            html_parts.append(
                '<div class="bp-section">'
                '<div class="bp-section-title">streams</div>'
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
                    f'<td>{html_module.escape(rel(stream.map_table) if stream.map_table else "")}</td></tr>'
                )
            html_parts.append('</table></div>')

            # Per-stream detail
            for stream_name, stream in active_streams:
                html_parts.append(
                    f'<div class="bp-section">'
                    f'<div class="bp-section-title">{stream_name}</div>'
                )

                # Metadata table (horizontal)
                meta_headers = ['name', 'format', 'items', 'map_table', 'files_contain_wildcards']
                meta_values = [
                    html_module.escape(stream.name),
                    html_module.escape(stream.format),
                    str(len(stream)),
                    html_module.escape(rel(stream.map_table) if stream.map_table else ''),
                    str(stream.files_contain_wildcards),
                ]
                if stream.metadata:
                    meta_headers.append('metadata')
                    meta_values.append(html_module.escape(
                        ", ".join(f"{k}={v}" for k, v in stream.metadata.items())
                    ))
                html_parts.append('<table class="bp-table"><tr>')
                for h in meta_headers:
                    html_parts.append(f'<th>{h}</th>')
                html_parts.append('</tr><tr>')
                for v in meta_values:
                    html_parts.append(f'<td>{v}</td>')
                html_parts.append('</tr></table>')

                # Data: map_table or id/file fallback
                map_data = stream._get_map_data()
                if map_data is not None and len(map_data) > 0:
                    _render_dataframe(map_data, html_parts)
                else:
                    # Fallback: render id/file pairs from iterator
                    items = list(stream)
                    n_items = len(items)
                    html_parts.append('<table class="bp-table"><tr><th>id</th><th>file</th></tr>')
                    if n_items <= 4:
                        display_items = items
                        ellipsis_after = None
                    else:
                        display_items = items[:2] + items[-2:]
                        ellipsis_after = 2

                    for item_i, (item_id, item_file) in enumerate(display_items):
                        if ellipsis_after is not None and item_i == ellipsis_after:
                            html_parts.append(
                                f'<tr class="bp-ellipsis"><td colspan="2">'
                                f'... {n_items - 4} more ...</td></tr>'
                            )
                        html_parts.append(
                            f"<tr><td>{html_module.escape(str(item_id))}</td>"
                            f"<td>{html_module.escape(rel(item_file) if item_file else '')}</td></tr>"
                        )
                    html_parts.append('</table>')

                html_parts.append('</div>')

        # --- Tables ---
        if "tables" in self._data and hasattr(self.tables, "_tables") and self.tables._tables:
            # Tables summary
            html_parts.append(
                '<div class="bp-section">'
                '<div class="bp-section-title">tables</div>'
            )
            html_parts.append(
                '<table class="bp-table">'
                '<tr><th>name</th><th>columns</th><th>path</th></tr>'
            )
            for name, info in self.tables._tables.items():
                cols = ", ".join(info.info.columns) if info.info.columns else ""
                path = rel(info.info.path) if info.info.path else ""
                html_parts.append(
                    f"<tr><td>{html_module.escape(name)}</td>"
                    f"<td>{html_module.escape(cols)}</td>"
                    f"<td>{html_module.escape(path)}</td></tr>"
                )
            html_parts.append('</table></div>')

            # Per-table detail
            for name, info in self.tables._tables.items():
                html_parts.append(
                    f'<div class="bp-section">'
                    f'<div class="bp-section-title">{html_module.escape(name)}</div>'
                )

                # Table metadata (horizontal)
                t_meta = info.info
                t_headers = ['name', 'path', 'columns', 'description', 'count']
                t_values = [
                    html_module.escape(t_meta.name),
                    html_module.escape(rel(t_meta.path) if t_meta.path else ''),
                    html_module.escape(", ".join(t_meta.columns) if t_meta.columns else ''),
                    html_module.escape(t_meta.description),
                    str(t_meta.count),
                ]
                html_parts.append('<table class="bp-table"><tr>')
                for h in t_headers:
                    html_parts.append(f'<th>{h}</th>')
                html_parts.append('</tr><tr>')
                for v in t_values:
                    html_parts.append(f'<td>{v}</td>')
                html_parts.append('</tr></table>')

                # Try to load and display CSV content
                if t_meta.path and os.path.exists(t_meta.path):
                    try:
                        table_df = pd.read_csv(t_meta.path)
                        if len(table_df) > 0:
                            _render_dataframe(table_df, html_parts)
                    except Exception:
                        pass

                html_parts.append('</div>')

        # --- Image streams: display inline ---
        for stream_name, stream in self.streams.items():
            if not isinstance(stream, DataStream) or len(stream) == 0:
                continue
            if stream.format.lower() not in self._IMAGE_FORMATS:
                continue

            html_parts.append(
                f'<div class="bp-section">'
                f'<div class="bp-section-title">{stream_name} '
                f'<span style="font-weight: normal; color: #666;">({stream.format}, {len(stream)} items)</span></div>'
            )

            for item_id, file_path in stream:
                if not file_path or not os.path.isfile(file_path):
                    continue

                if stream.format.lower() == "svg":
                    try:
                        with open(file_path, "r") as f:
                            svg_content = f.read()
                        html_parts.append(
                            f'<div style="margin: 4px 0;">'
                            f'<div style="color: #666; font-size: 0.85em;">{html_module.escape(str(item_id))}</div>'
                            f"{svg_content}</div>"
                        )
                    except Exception:
                        pass
                else:
                    try:
                        with open(file_path, "rb") as f:
                            img_data = base64.b64encode(f.read()).decode("utf-8")
                        mime = f"image/{stream.format.lower()}"
                        if stream.format.lower() in ("jpg", "jpeg"):
                            mime = "image/jpeg"
                        html_parts.append(
                            f'<div style="margin: 4px 0;">'
                            f'<div style="color: #666; font-size: 0.85em;">{html_module.escape(str(item_id))}</div>'
                            f'<img src="data:{mime};base64,{img_data}" '
                            f'style="max-width: 100%; height: auto;" /></div>'
                        )
                    except Exception:
                        pass

            html_parts.append("</div>")

        # --- Tool-specific visualization (e.g. 3D viewer) ---
        if hasattr(self, "tool") and hasattr(self.tool, "_repr_notebook_html"):
            try:
                tool_html = self.tool._repr_notebook_html(self)
                if tool_html:
                    html_parts.append(tool_html)
            except Exception:
                pass

        if not html_parts:
            return None
        return css + "\n" + "\n".join(html_parts)

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