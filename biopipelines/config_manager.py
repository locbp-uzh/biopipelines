# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Configuration management for BioPipelines.

Handles loading, validating, and managing the config.yaml file,
including pull/repull functionality from the repository.
"""

import os
import yaml
import shutil
import urllib.request
from typing import Dict, Any, List, Optional


class ConfigManager:
    """
    Manages BioPipelines configuration from config.yaml.

    Provides methods to load, validate, pull, and manage configuration
    for folder paths and default environments.
    """

    _instance = None
    _config = None

    def __new__(cls):
        """Singleton pattern to ensure only one config manager exists."""
        if cls._instance is None:
            cls._instance = super(ConfigManager, cls).__new__(cls)
        return cls._instance

    def __init__(self):
        """Initialize configuration manager."""
        if self._config is None:
            self._config = self._load_config()

    @staticmethod
    def _get_config_path() -> str:
        """
        Get the path to the config file.

        On Google Colab, looks for colab.yaml first; falls back to config.yaml.

        Returns:
            Absolute path to the config file in the repository root
        """
        # Get repository root (assumes this file is in biopipelines/)
        script_dir = os.path.dirname(os.path.abspath(__file__))
        repo_root = os.path.dirname(script_dir)

        # Auto-detect Google Colab
        try:
            import google.colab  # noqa: F401
            colab_path = os.path.join(repo_root, "colab.yaml")
            if os.path.exists(colab_path):
                return colab_path
        except ImportError:
            pass

        return os.path.join(repo_root, "config.yaml")

    def _load_config(self) -> Dict[str, Any]:
        """
        Load configuration from config.yaml.

        Returns:
            Dictionary containing configuration

        Raises:
            FileNotFoundError: If config.yaml doesn't exist
            ValueError: If config.yaml is invalid
        """
        config_path = self._get_config_path()

        if not os.path.exists(config_path):
            raise FileNotFoundError(
                f"config.yaml not found at {config_path}. "
                f"Run ConfigManager.pull_config() to download the default configuration."
            )

        try:
            with open(config_path, 'r') as f:
                config = yaml.safe_load(f)

            # Validate required sections
            required_sections = ['folders']
            for section in required_sections:
                if section not in config:
                    raise ValueError(f"config.yaml missing required section: {section}")

            # Validate folder subsections
            folders = config.get('folders', {})
            required_folder_sections = ['base', 'infrastructure', 'repositories']
            for folder_section in required_folder_sections:
                if folder_section not in folders:
                    raise ValueError(
                        f"config.yaml folders section missing required subsection: {folder_section}"
                    )

            return config

        except yaml.YAMLError as e:
            raise ValueError(f"Invalid YAML in config.yaml: {e}")

    def reload(self):
        """Reload configuration from disk."""
        self._config = self._load_config()

    def get_folder_config(self) -> Dict[str, Any]:
        """
        Get folder configuration section.

        Returns:
            Dictionary containing folder paths and settings
        """
        return self._config.get('folders', {})

    def get_environment(self, tool_name: str) -> Optional[str]:
        """
        Get default environment for a tool.

        Args:
            tool_name: Name of the tool (e.g., "Boltz2", "RFdiffusion")

        Returns:
            Environment name or None if not configured
        """
        environments = self._config.get('environments', {})
        return environments.get(tool_name)

    def get_all_environments(self) -> Dict[str, Optional[str]]:
        """
        Get all configured environments.

        Returns:
            Dictionary mapping tool names to environment names
        """
        return self._config.get('environments', {})

    def _get_cluster_config(self) -> Dict[str, Any]:
        """Get cluster configuration section. Raises if missing."""
        cluster = self._config.get('cluster')
        if cluster is None:
            raise KeyError("'cluster' section not found in config.yaml. Please add it (see config.yaml template).")
        return cluster

    def get_emails(self) -> Dict[str, str]:
        """Get email mapping (Unix username -> email address)."""
        return self._get_cluster_config().get('emails', {})

    def get_env_manager(self) -> str:
        """Get configured environment manager ("mamba" or "conda")."""
        val = self._get_cluster_config().get('env_manager')
        if not val:
            raise KeyError("'env_manager' not set in config.yaml cluster section.")
        return val

    def get_container_executor(self) -> str:
        """Get configured container executor ("apptainer" or "singularity")."""
        val = self._get_cluster_config().get('container_executor')
        if not val:
            raise KeyError("'container_executor' not set in config.yaml cluster section.")
        return val

    def get_container(self, tool_name: str) -> Optional[str]:
        """Get container path for a tool, or None if not configured.

        Args:
            tool_name: Name of the tool (e.g., "RFdiffusionAllAtom", "PLIP")

        Returns:
            Container path string (may contain folder placeholders), or None
        """
        containers = self._config.get('containers', {})
        if containers is None:
            return None
        return containers.get(tool_name)

    def get_slurm_modules(self) -> List[str]:
        """Get list of modules to load in SLURM wrapper scripts."""
        val = self._get_cluster_config().get('slurm_modules')
        if val is None:
            raise KeyError("'slurm_modules' not set in config.yaml cluster section. Use [] for no modules.")
        return val

    def get_shell_hook_command(self) -> str:
        """Get the shell hook initialization command derived from env_manager."""
        mgr = self.get_env_manager()
        if mgr == "pip":
            return ""
        if mgr == "conda":
            return 'eval "$(conda shell.bash hook)"'
        return f'eval "$({mgr} shell hook --shell bash)"'

    def get_activate_command(self, env_name: str) -> str:
        """Get the activate command for a specific environment."""
        mgr = self.get_env_manager()
        if mgr == "pip":
            return ""
        return f'{mgr} activate {env_name}'

    def get_module_load_line(self) -> str:
        """Get the full 'module load ...' line for SLURM scripts, or empty string if none."""
        modules = self.get_slurm_modules()
        if not modules:
            return ""
        return "module load " + " ".join(modules)

    def get_repository_url(self) -> str:
        """
        Get repository URL for config updates.

        Returns:
            URL to fetch config.yaml from repository
        """
        repo_config = self._config.get('repository', {})
        return repo_config.get('url', '')

    @classmethod
    def pull_config(cls, force: bool = False) -> str:
        """
        Pull default config.yaml from repository.

        Args:
            force: If True, overwrites existing config (default: False)

        Returns:
            Path to downloaded config file

        Raises:
            FileExistsError: If config exists and force=False
            urllib.error.URLError: If download fails
        """
        config_path = cls._get_config_path()

        # Check if config already exists
        if os.path.exists(config_path) and not force:
            raise FileExistsError(
                f"config.yaml already exists at {config_path}. "
                f"Use force=True to overwrite, or use repull_config() to backup and update."
            )

        # Try to get URL from existing config, otherwise use default
        default_url = "https://raw.githubusercontent.com/gquargnali/biopipelines/main/config.yaml"

        try:
            # If config exists, try to read repository URL from it
            if os.path.exists(config_path):
                with open(config_path, 'r') as f:
                    existing_config = yaml.safe_load(f)
                    repo_config = existing_config.get('repository', {})
                    url = repo_config.get('url', default_url)
            else:
                url = default_url
        except Exception:
            url = default_url

        # Download config
        try:
            print(f"Downloading config.yaml from {url}...")
            urllib.request.urlretrieve(url, config_path)
            print(f"Config saved to {config_path}")
            return config_path
        except Exception as e:
            raise urllib.error.URLError(f"Failed to download config.yaml: {e}")

    @classmethod
    def repull_config(cls, backup: bool = True) -> str:
        """
        Re-download config.yaml from repository, optionally backing up existing one.

        Args:
            backup: If True, creates a backup of existing config (default: True)

        Returns:
            Path to downloaded config file
        """
        config_path = cls._get_config_path()

        # Create backup if requested and config exists
        if backup and os.path.exists(config_path):
            backup_path = f"{config_path}.backup"
            counter = 1
            while os.path.exists(backup_path):
                backup_path = f"{config_path}.backup{counter}"
                counter += 1

            shutil.copy2(config_path, backup_path)
            print(f"Backed up existing config to {backup_path}")

        # Pull new config with force=True
        return cls.pull_config(force=True)

    @classmethod
    def create_user_config(cls) -> str:
        """
        Create a user config.yaml if it doesn't exist, using the default template.

        This is called automatically by ConfigManager when no config is found.

        Returns:
            Path to created config file
        """
        config_path = cls._get_config_path()

        if os.path.exists(config_path):
            return config_path

        # Try to pull from repository
        try:
            return cls.pull_config(force=False)
        except Exception as e:
            print(f"Warning: Could not pull config from repository: {e}")
            print("Please ensure config.yaml exists in the repository root.")
            raise FileNotFoundError(
                f"config.yaml not found and could not be downloaded. "
                f"Expected location: {config_path}"
            )

    def __repr__(self) -> str:
        """String representation of config manager."""
        return f"ConfigManager(config_path={self._get_config_path()})"


# Convenience functions for direct access
def get_environment(tool_name: str) -> Optional[str]:
    """
    Get default environment for a tool.

    Args:
        tool_name: Name of the tool

    Returns:
        Environment name or None
    """
    manager = ConfigManager()
    return manager.get_environment(tool_name)


def get_folder_config() -> Dict[str, Any]:
    """
    Get folder configuration.

    Returns:
        Dictionary containing folder settings
    """
    manager = ConfigManager()
    return manager.get_folder_config()


def reload_config():
    """Reload configuration from disk."""
    manager = ConfigManager()
    manager.reload()
