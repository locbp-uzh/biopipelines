"""
File Path Manager for automatic path management in tools.

Provides descriptor-based automatic file path generation,
eliminating manual _initialize_file_paths() and _setup_file_paths() methods.
"""

import os
from typing import Callable, Union, Optional


class FilePathDescriptor:
    """
    Descriptor for auto-managed file paths.

    Automatically generates file paths based on output_folder and templates.
    Paths are lazily evaluated when first accessed.

    Examples:
        >>> class MyTool(BaseConfig):
        ...     # Simple filename
        ...     results_csv = FilePathDescriptor("results.csv")
        ...
        ...     # Template with placeholders
        ...     queries_csv = FilePathDescriptor("{pipeline_name}_queries.csv")
        ...
        ...     # Callable template for complex logic
        ...     output_pdb = FilePathDescriptor(
        ...         lambda self: f"{self.name}_{self.execution_order}.pdb"
        ...     )
        ...
        ...     # Custom folder
        ...     helper_script = FilePathDescriptor(
        ...         "parse.py",
        ...         folder_key="HelpScripts"
        ...     )
    """

    def __init__(self,
                 filename_template: Union[str, Callable],
                 folder_key: str = "output_folder"):
        """
        Initialize file path descriptor.

        Args:
            filename_template: Either:
                - String with optional {placeholders} (e.g., "{tool_name}_results.csv")
                - Callable taking instance and returning filename string
            folder_key: Attribute name to get folder from (default: "output_folder")

        Supported placeholders in string templates:
            - {tool_name}: Tool's TOOL_NAME attribute
            - {pipeline_name}: Extracted pipeline name
            - {job_name}: Job name
            - {execution_order}: Tool's execution order in pipeline
        """
        self.filename_template = filename_template
        self.folder_key = folder_key
        self._attr_name = None

    def __set_name__(self, owner, name):
        """Called when descriptor is assigned to class attribute."""
        self._attr_name = f"_{name}_cached"

    def __get__(self, instance, owner):
        """Get the file path, computing it lazily on first access."""
        if instance is None:
            return self

        # Check if already computed and cached
        if hasattr(instance, self._attr_name):
            return getattr(instance, self._attr_name)

        # Compute the path
        path = self._compute_path(instance)

        # Cache the result
        setattr(instance, self._attr_name, path)

        return path

    def __set__(self, instance, value):
        """Allow manual override of the path if needed."""
        setattr(instance, self._attr_name, value)

    def _compute_path(self, instance) -> str:
        """Compute the full file path."""
        # Get the folder
        if self.folder_key == "output_folder":
            folder = getattr(instance, 'output_folder', '')
        elif hasattr(instance, 'folders') and isinstance(instance.folders, dict):
            folder = instance.folders.get(self.folder_key, '')
        else:
            folder = getattr(instance, self.folder_key, '')

        # Get the filename
        if callable(self.filename_template):
            # Callable template - call it with instance
            filename = self.filename_template(instance)
        else:
            # String template - format with available attributes
            filename = self._format_template(instance, self.filename_template)

        # Combine folder and filename
        if folder:
            return os.path.join(folder, filename)
        return filename

    def _format_template(self, instance, template: str) -> str:
        """Format a string template with instance attributes."""
        # Build format dict with available attributes
        format_dict = {}

        # Tool name
        if hasattr(instance, 'TOOL_NAME'):
            format_dict['tool_name'] = instance.TOOL_NAME

        # Job name
        if hasattr(instance, 'name'):
            format_dict['job_name'] = instance.name
        elif hasattr(instance, 'job_name'):
            format_dict['job_name'] = instance.job_name

        # Execution order
        if hasattr(instance, 'execution_order'):
            format_dict['execution_order'] = instance.execution_order

        # Pipeline name (try to extract if method exists)
        if hasattr(instance, '_extract_pipeline_name'):
            try:
                format_dict['pipeline_name'] = instance._extract_pipeline_name()
            except:
                format_dict['pipeline_name'] = 'pipeline'
        else:
            # Fallback: try to extract from output_folder
            if hasattr(instance, 'output_folder'):
                folder_parts = instance.output_folder.split(os.sep)
                # Look for pattern like "JobName_NNN/N_ToolName"
                for i, part in enumerate(folder_parts):
                    if '_' in part and i > 0:
                        format_dict['pipeline_name'] = folder_parts[i-1]
                        break
                else:
                    format_dict['pipeline_name'] = 'pipeline'
            else:
                format_dict['pipeline_name'] = 'pipeline'

        # Format the template
        try:
            return template.format(**format_dict)
        except KeyError as e:
            raise ValueError(
                f"Template '{template}' contains unknown placeholder: {e}. "
                f"Available: {list(format_dict.keys())}"
            )

    def invalidate_cache(self, instance):
        """
        Invalidate cached path (useful if instance attributes change).

        Args:
            instance: Tool instance

        Example:
            >>> tool.output_folder = "/new/path"
            >>> FilePathDescriptor.invalidate_cache(tool)
            >>> # Next access will recompute path
        """
        if hasattr(instance, self._attr_name):
            delattr(instance, self._attr_name)


class SubfolderDescriptor(FilePathDescriptor):
    """
    Descriptor for subfolders within output_folder.

    Similar to FilePathDescriptor but creates directories on first access.

    Example:
        >>> class MyTool(BaseConfig):
        ...     seqs_folder = SubfolderDescriptor("seqs")
        ...     results_folder = SubfolderDescriptor("Results")
    """

    def __get__(self, instance, owner):
        """Get folder path and ensure it exists."""
        path = super().__get__(instance, owner)

        # Create directory if it doesn't exist (only when accessed)
        if instance is not None and path and not os.path.exists(path):
            try:
                os.makedirs(path, exist_ok=True)
            except OSError:
                pass  # May fail during pipeline setup before folders exist

        return path


class HelperScriptDescriptor(FilePathDescriptor):
    """
    Descriptor specifically for helper script paths.

    Automatically uses HelpScripts folder and adds 'pipe_' prefix convention.

    Example:
        >>> class MyTool(BaseConfig):
        ...     filter_script = HelperScriptDescriptor("filter.py")
        ...     # Resolves to: {HelpScripts}/pipe_filter.py
        ...
        ...     # Custom name (no auto-prefix)
        ...     custom_script = HelperScriptDescriptor("my_script.py", add_prefix=False)
    """

    def __init__(self, script_name: str, add_prefix: bool = True):
        """
        Initialize helper script descriptor.

        Args:
            script_name: Name of the script file
            add_prefix: Whether to add 'pipe_' prefix (default: True)
        """
        if add_prefix and not script_name.startswith('pipe_'):
            script_name = f"pipe_{script_name}"

        super().__init__(script_name, folder_key="HelpScripts")
