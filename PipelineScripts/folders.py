# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Folder management utilities for pipeline execution.

Handles automatic creation and organization of all pipeline-related directories
following the established conventions from the biopipelines.
"""

import getpass
import os
import re
from typing import Dict

try:
    from .config_manager import ConfigManager
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from config_manager import ConfigManager


class FolderResolutionError(Exception):
    """Raised when a folder key cannot be resolved."""
    pass


# Regex pattern to match placeholders like <key>
PLACEHOLDER_PATTERN = re.compile(r'<([a-zA-Z_][a-zA-Z0-9_]*)>')

# Order in which config sections are resolved
SECTION_ORDER = ['base', 'infrastructure', 'repositories', 'cache', 'derived', 'server']


class FolderManager:
    """
    Manages folder structure for pipeline execution.

    Creates and maintains the standard directory hierarchy used by all tools,
    following conventions established in the biopipelines.
    """

    def __init__(self, project: str, job: str, local_output: bool=False):
        """
        Initialize folder manager for a pipeline.

        Args:
            project: Name of the folder (used for output folders)
            job: Name of the specific job (a unique numeric id NNN will be appended to it) (used for output folders)
            local_output: If True, write output to ./tests/ (current working
                        directory) instead of the config-configured path.
        """
        self._folders: Dict[str, str] = {}
        self._local_output = local_output

        # Load configuration
        config_manager = ConfigManager()
        folder_config = config_manager.get_folder_config()

        # Inject runtime values
        self._inject_runtime_values()

        # Resolve all config sections in order
        self._resolve_config_sections(folder_config)

        # Override biopipelines_output for local output mode
        if local_output:
            self._folders["biopipelines_output"] = os.path.join(os.getcwd(), "tests")

        # Setup runtime paths (project, output, runtime, logs)
        self._setup_runtime_paths(project, job)

    def _inject_runtime_values(self):
        """Inject runtime values that are computed from the environment."""
        # Resolve from package location (parent of PipelineScripts/)
        # This works regardless of the current working directory,
        biopipelines_folder = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        user_name = getpass.getuser()
        self.user_name = user_name

        # These are the runtime placeholders available in config
        # - username: actual username (e.g., gquarg)
        # - cwd: biopipelines root directory
        self._folders["username"] = user_name
        self._folders["cwd"] = biopipelines_folder

    def _resolve_config_sections(self, folder_config: Dict):
        """Resolve all config sections in the defined order."""
        for section_name in SECTION_ORDER:
            section = folder_config.get(section_name, {})
            for key, path_template in section.items():
                resolved_path = self._resolve_path(path_template, key)
                self._folders[key] = resolved_path

    def _resolve_path(self, path_template: str, key_being_resolved: str) -> str:
        """
        Resolve a path template by substituting all placeholders.

        Args:
            path_template: Path string with <placeholder> syntax
            key_being_resolved: The key being resolved (for error messages)

        Returns:
            Fully resolved path string

        Raises:
            FolderResolutionError: If a placeholder cannot be resolved
        """
        result = path_template
        max_iterations = 10  # Prevent infinite loops from circular references

        for _ in range(max_iterations):
            matches = PLACEHOLDER_PATTERN.findall(result)
            if not matches:
                break

            for placeholder in matches:
                if placeholder not in self._folders:
                    raise FolderResolutionError(
                        f"Cannot resolve folder key '{key_being_resolved}': "
                        f"placeholder '<{placeholder}>' not found. "
                        f"Available keys: {list(self._folders.keys())}"
                    )
                result = result.replace(f"<{placeholder}>", self._folders[placeholder])
        else:
            # If we've hit max iterations, there's likely a circular reference
            raise FolderResolutionError(
                f"Cannot resolve folder key '{key_being_resolved}': "
                f"possible circular reference in path '{path_template}'"
            )

        return result

    def _setup_runtime_paths(self, project: str, job: str):
        """Setup runtime paths that depend on project/job names."""
        # Project path
        self._folders["project"] = os.path.join(self._folders["biopipelines_output"], project)

        # Create necessary directories
        folders_to_create = ["biopipelines_output", "PDBs", "Ligands"]
        if not self._local_output:
            folders_to_create.append("user")
        for folder_key in folders_to_create:
            os.makedirs(self._folders[folder_key], exist_ok=True)
        os.makedirs(self._folders["project"], exist_ok=True)

        # Output path with unique name
        self._folders["output"] = FolderManager.unique_name(
            self._folders["project"], job, full_path=True
        )
        os.makedirs(self._folders["output"], exist_ok=True)

        # Extract job ID from the output folder name
        output_basename = os.path.basename(self._folders["output"])
        if "_" in output_basename:
            self.job_id = output_basename.split("_")[-1]
        else:
            raise ValueError(f"Could not extract job ID from output folder name: {output_basename}")

        # Runtime and logs folders
        self._folders["runtime"] = os.path.join(self._folders["output"], "RunTime")
        os.makedirs(self._folders["runtime"], exist_ok=True)

        self._folders["logs"] = os.path.join(self._folders["output"], "Logs")
        os.makedirs(self._folders["logs"], exist_ok=True)

    def get_folders(self) -> Dict[str, str]:
        """
        Get dictionary of all folder paths.

        Returns:
            Dictionary mapping folder names to absolute paths
        """
        return self._folders.copy()

    def unique_name(directory: str, root: str, ext: str = "", w: int = 3, full_path: bool = False) -> str:
        """
        Generate unique filename/foldername to avoid conflicts.

        Based on the unique_name function from biopipelines.

        Args:
            directory: Directory to check for existing files
            root: Base name for file/folder
            ext: File extension (optional)
            w: Width of numeric suffix (default 3)

        Returns:
            Unique name that doesn't exist in directory
        """
        i = 1
        while True:
            u_name = root + "_" + "{:0>{width}}".format(i, width=w)
            if not os.path.exists(os.path.join(directory, u_name + ext)):
                if full_path:   return os.path.join(directory, u_name + ext)
                else:           return u_name + ext
            i += 1

    def __getitem__(self, key: str) -> str:
        """Allow dictionary-like access to folders."""
        if key not in self._folders:
            raise FolderResolutionError(
                f"Folder key '{key}' not found. "
                f"Available keys: {list(self._folders.keys())}"
            )
        return self._folders[key]

    def __contains__(self, key: str) -> bool:
        """Check if folder key exists."""
        return key in self._folders
