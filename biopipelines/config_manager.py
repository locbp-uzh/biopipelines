# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Configuration management for BioPipelines.

Handles loading, validating, and managing the config.<variant>.yaml files,
including pull/repull functionality from the repository.

Variants:
  - cluster (default outside Colab) -> config.cluster.yaml
  - colab   (default inside Colab)  -> config.colab.yaml
  - <other> (user-defined site)     -> config.<other>.yaml

Set the variant via `Pipeline(config="<variant>")` (recommended) or by
constructing `ConfigManager(variant="<variant>")` before any other code reads
the config. Without an explicit variant, the loader auto-detects in this
order: colab (if google.colab imports), then any config.<variant>.yaml
whose `machine.username` matches the current Unix user, otherwise cluster.
"""

import datetime
import getpass
import glob
import os
import shutil
import urllib.request
from typing import Dict, Any, List, Optional


def backup_file(path) -> str:
    """Copy ``path`` to a timestamped ``<path>.<YYYYMMDD-HHMMSS>.bak`` and
    return the backup path.

    Timestamped so successive edits never clobber an earlier backup (the old
    fixed ``<path>.bak`` lost the previous copy on every save). If two saves
    land in the same second, a numeric suffix is appended to stay unique.
    """
    import shutil as _shutil
    from pathlib import Path
    path = Path(path)
    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    bak = path.with_name(f"{path.name}.{stamp}.bak")
    n = 1
    while bak.exists():
        bak = path.with_name(f"{path.name}.{stamp}-{n}.bak")
        n += 1
    _shutil.copy2(path, bak)
    return str(bak)


def _autodetect_variant() -> str:
    """Pick a variant name based on the runtime environment.

    Selection order:
      1. If running inside Google Colab, return ``"colab"``.
      2. Otherwise, scan every repo-root ``config.<variant>.yaml`` (excluding
         ``config.colab.yaml``) and classify each by its ``machine.username``:
           - **match**: ``username`` equals the current Unix user.
           - **wildcard**: ``username`` is absent or empty.
         If exactly one match exists, return it. If no match exists but at
         least one wildcard does, return the first wildcard (alphabetical).
      3. Fall back to ``"cluster"``.

    Raises:
        RuntimeError: If two or more variants claim the same ``username``
            (ambiguous — the user must select one explicitly).
    """
    try:
        import google.colab  # noqa: F401
        return "colab"
    except ImportError:
        pass

    try:
        current_user = getpass.getuser()
    except Exception:
        current_user = ""

    try:
        import yaml
    except ImportError:
        return "cluster"

    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    matches: List[str] = []
    wildcards: List[str] = []

    for path in sorted(glob.glob(os.path.join(repo_root, "config.*.yaml"))):
        base = os.path.basename(path)
        variant = base[len("config."):-len(".yaml")]
        if variant == "colab":
            continue
        try:
            with open(path, "r") as f:
                data = yaml.safe_load(f) or {}
        except Exception:
            continue
        username = (data.get("machine") or {}).get("username") or ""
        if username and current_user and username == current_user:
            matches.append(variant)
        elif not username:
            wildcards.append(variant)

    if len(matches) > 1:
        raise RuntimeError(
            f"Ambiguous variant auto-detection: multiple config files claim "
            f"username={current_user!r} ({matches}). Select one explicitly via "
            f"Pipeline(config=...)."
        )
    if matches:
        return matches[0]
    if wildcards:
        return wildcards[0]
    return "cluster"


class ConfigManager:
    """
    Manages BioPipelines configuration from config.<variant>.yaml.

    Provides methods to load, validate, pull, and manage configuration
    for folder paths and default environments.
    """

    _instance = None
    _config = None
    _variant: Optional[str] = None

    def __new__(cls, variant: Optional[str] = None):
        """Singleton pattern to ensure only one config manager exists."""
        if cls._instance is None:
            cls._instance = super(ConfigManager, cls).__new__(cls)
        return cls._instance

    def __init__(self, variant: Optional[str] = None):
        """Initialize configuration manager.

        Args:
            variant: Config variant to load (e.g. "cluster", "colab", or any
                user-defined name matching `config.<variant>.yaml`). If None
                and no variant has been set previously, auto-detects from the
                runtime environment. Re-initializing with a different variant
                resets the cached config so the new one is loaded.
        """
        if variant is not None and variant != type(self)._variant:
            type(self)._variant = variant
            type(self)._config = None
        elif type(self)._variant is None:
            type(self)._variant = _autodetect_variant()

        if self._config is None:
            type(self)._config = self._load_config()

    @classmethod
    def get_variant(cls) -> str:
        """Return the active variant name (e.g. 'cluster', 'colab')."""
        if cls._variant is None:
            cls._variant = _autodetect_variant()
        return cls._variant

    @classmethod
    def _get_config_path(cls, variant: Optional[str] = None) -> str:
        """
        Get the path to the config file for the given (or active) variant.

        Args:
            variant: Variant name. If None, uses the active variant
                     (auto-detected on first call).

        Returns:
            Absolute path to `config.<variant>.yaml` in the repository root.
        """
        # Get repository root (assumes this file is in biopipelines/)
        script_dir = os.path.dirname(os.path.abspath(__file__))
        repo_root = os.path.dirname(script_dir)

        if variant is None:
            variant = cls._variant or _autodetect_variant()

        return os.path.join(repo_root, f"config.{variant}.yaml")

    def _load_config(self) -> Dict[str, Any]:
        """
        Load configuration from `config.<variant>.yaml`.

        Returns:
            Dictionary containing configuration

        Raises:
            FileNotFoundError: If the config file doesn't exist
            ValueError: If the config file is invalid
        """
        config_path = self._get_config_path()
        config_name = os.path.basename(config_path)

        if not os.path.exists(config_path):
            raise FileNotFoundError(
                f"{config_name} not found at {config_path}. "
                f"Run ConfigManager.pull_config() to download the default configuration."
            )

        try:
            import yaml
            with open(config_path, 'r') as f:
                config = yaml.safe_load(f)

            # Validate required sections
            required_sections = ['folders']
            for section in required_sections:
                if section not in config:
                    raise ValueError(f"{config_name} missing required section: {section}")

            # Validate folder subsections
            folders = config.get('folders', {})
            required_folder_sections = ['base', 'infrastructure', 'repositories']
            for folder_section in required_folder_sections:
                if folder_section not in folders:
                    raise ValueError(
                        f"{config_name} folders section missing required subsection: {folder_section}"
                    )

            self._validate_shell_safety(config, config_name)

            return config

        except yaml.YAMLError as e:
            raise ValueError(f"Invalid YAML in {config_name}: {e}")

    # Characters that would break double-quoted bash interpolation.
    _SHELL_UNSAFE_CHARS = ('"', '`', '$', '\\')

    @classmethod
    def _reject_shell_unsafe(cls, field: str, value: Any, config_name: str) -> None:
        """Reject characters that break double-quoted bash interpolation.

        `field` is the dotted config path (e.g. "folders.base.home"); used in
        error messages so the user can locate the offending line.
        """
        if value is None:
            return
        if not isinstance(value, str):
            raise ValueError(
                f"{config_name}: {field!r} must be a string, got "
                f"{type(value).__name__}"
            )
        for ch in cls._SHELL_UNSAFE_CHARS:
            if ch in value:
                raise ValueError(
                    f"{config_name}: {field}={value!r} contains {ch!r}, which "
                    'would break the generated shell script. Remove the four '
                    'characters " ` $ \\ from config values that are used as '
                    'paths or module names.'
                )

    @classmethod
    def _validate_shell_safety(cls, config: Dict[str, Any], config_name: str) -> None:
        """Validate config fields that are interpolated into generated bash.

        Narrow by design: only covers fields whose values end up in a shell
        script (paths, module names, a small handful of enum-ish settings).
        Does not attempt schema validation beyond that.
        """
        # All folder path values across every subsection.
        folders = config.get('folders') or {}
        if isinstance(folders, dict):
            for subsection, entries in folders.items():
                if not isinstance(entries, dict):
                    continue
                for key, value in entries.items():
                    cls._reject_shell_unsafe(
                        f"folders.{subsection}.{key}", value, config_name
                    )

        # Container paths.
        containers = config.get('containers') or {}
        if isinstance(containers, dict):
            for key, value in containers.items():
                cls._reject_shell_unsafe(
                    f"containers.{key}", value, config_name
                )

        # Machine section. env_manager and scheduler are dicts of the form
        # {name, init, [modules]}. Their `name` fields are interpolated
        # unquoted into shell expressions (e.g. `eval "$(<mgr> shell hook
        # ...)"` is no longer auto-derived but the value still flows into
        # `<mgr> activate <env>`), so they must be restricted to a fixed
        # allowlist. The `init` lists carry raw bash and are emitted
        # verbatim — we only sanity-check against NUL bytes (paste-error
        # protection); arbitrary code in this field is intentional and
        # documented (the field is a trust boundary, not user input).
        machine = config.get('machine') or config.get('cluster') or {}
        if isinstance(machine, dict):
            _NAME_ENUMS = {
                'env_manager': ('mamba', 'conda', 'micromamba', 'pip'),
                'scheduler':   ('slurm', 'colab', 'none'),
            }
            for key, allowed in _NAME_ENUMS.items():
                block = machine.get(key)
                if block is None:
                    continue
                if not isinstance(block, dict):
                    raise ValueError(
                        f"{config_name}: machine.{key} must be a dict "
                        f"of the form {{name: ..., init: [...]}}, got "
                        f"{type(block).__name__}."
                    )
                name = block.get('name')
                if name is None:
                    continue  # _machine_block() will raise on use; leave to caller.
                if name not in allowed:
                    raise ValueError(
                        f"{config_name}: machine.{key}.name={name!r} is "
                        f"not allowed. Must be one of {allowed}."
                    )
                init_lines = block.get('init') or []
                if not isinstance(init_lines, list):
                    raise ValueError(
                        f"{config_name}: machine.{key}.init must be a "
                        f"list of bash strings, got {type(init_lines).__name__}."
                    )
                for i, line in enumerate(init_lines):
                    if not isinstance(line, str):
                        raise ValueError(
                            f"{config_name}: machine.{key}.init[{i}] must "
                            f"be a string, got {type(line).__name__}."
                        )
                    if "\x00" in line:
                        raise ValueError(
                            f"{config_name}: machine.{key}.init[{i}] "
                            "contains a NUL byte."
                        )

            # container_executor is still a flat string (no init / no modules
            # attached today); same allowlist as before.
            container_executor = machine.get('container_executor')
            if container_executor is not None:
                if container_executor not in ('apptainer', 'singularity', 'none'):
                    raise ValueError(
                        f"{config_name}: machine.container_executor="
                        f"{container_executor!r} is not allowed. "
                        "Must be one of ('apptainer', 'singularity', 'none')."
                    )

            # username is free-form (a Unix login) but still flows into paths
            # and the email-lookup key; sanity-check for shell metas.
            cls._reject_shell_unsafe(
                "machine.username", machine.get('username'), config_name
            )

            # scheduler.modules is interpolated as `module load <m1> <m2> ...`
            # so each entry must be a tight identifier with no shell metas.
            scheduler_block = machine.get('scheduler')
            if isinstance(scheduler_block, dict):
                modules = scheduler_block.get('modules') or []
                if isinstance(modules, list):
                    for i, mod in enumerate(modules):
                        cls._reject_shell_unsafe(
                            f"machine.scheduler.modules[{i}]", mod, config_name
                        )

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

    def _get_machine_config(self) -> Dict[str, Any]:
        """Get machine configuration section. Raises if missing.

        Accepts the legacy 'cluster' key as a fallback so older config files
        keep loading until they are migrated.
        """
        machine = self._config.get('machine')
        if machine is None:
            machine = self._config.get('cluster')
        if machine is None:
            config_name = os.path.basename(self._get_config_path())
            raise KeyError(
                f"'machine' section not found in {config_name}. "
                f"Please add it (see config.cluster.yaml template)."
            )
        return machine

    def get_email(self) -> str:
        """Get the SLURM-notification email address (empty string if unset).

        Each config serves a single user, so this is a flat ``machine.email``
        string. The legacy ``machine.emails`` dict ({username: address}) is
        still read for back-compat: the current user's entry wins, else the
        first address present.
        """
        machine = self._get_machine_config()
        email = machine.get('email')
        if email is not None:
            return email
        legacy = machine.get('emails')
        if isinstance(legacy, dict) and legacy:
            import getpass
            try:
                user = getpass.getuser()
            except Exception:
                user = None
            if user in legacy:
                return legacy[user]
            return next(iter(legacy.values()))
        return ""

    def _machine_block(self, key: str) -> Dict[str, Any]:
        """Return the dict at ``machine.<key>``; raise if missing or wrong shape.

        Both ``env_manager`` and ``scheduler`` are required to be dicts of
        the form ``{name: ..., init: [...], (modules: [...])}``. This is
        the single source of truth — neither the flat-string nor an
        implicit default form is supported.
        """
        machine = self._get_machine_config()
        block = machine.get(key)
        if block is None:
            config_name = os.path.basename(self._get_config_path())
            raise KeyError(
                f"'machine.{key}' is missing in {config_name}. "
                f"Use the dict form: {key}: {{name: ..., init: [...]}}"
            )
        if not isinstance(block, dict):
            config_name = os.path.basename(self._get_config_path())
            raise ValueError(
                f"{config_name}: machine.{key} must be a dict "
                f"(got {type(block).__name__}). Use the dict form: "
                f"{key}: {{name: ..., init: [...]}}"
            )
        if not block.get('name'):
            config_name = os.path.basename(self._get_config_path())
            raise KeyError(
                f"{config_name}: machine.{key}.name is missing or empty."
            )
        return block

    def get_env_manager(self) -> str:
        """Get configured environment manager ("mamba", "conda", "micromamba", "pip")."""
        return self._machine_block('env_manager')['name']

    def get_env_manager_init(self) -> List[str]:
        """Raw bash lines that activate the env-manager shell hook.

        Emitted verbatim before any per-tool ``<mgr> activate <env>``
        line. The user supplies this in the YAML — the framework no
        longer derives it from the manager name.
        """
        block = self._machine_block('env_manager')
        return list(block.get('init', []) or [])

    def get_scheduler(self) -> str:
        """Get configured scheduler ("slurm", "colab", or "none")."""
        return self._machine_block('scheduler')['name']

    def get_scheduler_init(self) -> List[str]:
        """Raw bash lines that initialise the scheduler's host environment.

        Emitted verbatim at the very top of every generated batch script,
        before any GPU setup, module load, or per-tool activation. Use
        this to e.g. source ``/etc/profile.d/lmod.sh`` so the ``module``
        builtin works in non-login SLURM shells.
        """
        block = self._machine_block('scheduler')
        return list(block.get('init', []) or [])

    def get_container_executor(self) -> str:
        """Get configured container executor ("apptainer" or "singularity")."""
        val = self._get_machine_config().get('container_executor')
        if not val:
            config_name = os.path.basename(self._get_config_path())
            raise KeyError(f"'container_executor' not set in {config_name} machine section.")
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

    def get_container_bind_sections(self) -> List[str]:
        """Folder sections whose resolved paths are auto-bound into containers.

        Apptainer/Singularity only bind $HOME, /tmp, $PWD by default, so
        pipeline outputs and shared caches aren't visible inside the container
        unless we pass them explicitly with -B. These sections are the safe,
        generally-useful ones.
        """
        return ['base', 'infrastructure', 'cache']

    def get_slurm_modules(self) -> List[str]:
        """Modules to load in SLURM batch scripts (machine.scheduler.modules)."""
        block = self._machine_block('scheduler')
        # 'modules' is optional; an empty list is fine.
        modules = block.get('modules')
        if modules is None:
            return []
        if not isinstance(modules, list):
            config_name = os.path.basename(self._get_config_path())
            raise ValueError(
                f"{config_name}: machine.scheduler.modules must be a list, "
                f"got {type(modules).__name__}."
            )
        return list(modules)

    def get_shell_hook_command(self) -> str:
        """Bash to source the env-manager into the current shell.

        Returns the user-supplied ``machine.env_manager.init`` lines
        joined with newlines (so the existing call sites that interpolate
        a single string still work). Empty string for env managers like
        ``pip`` that don't need a shell hook (the user expresses this by
        leaving init empty).
        """
        return "\n".join(self.get_env_manager_init())

    def get_scheduler_init_block(self) -> str:
        """Bash to initialise the scheduler's host environment.

        Returns the user-supplied ``machine.scheduler.init`` lines
        joined with newlines, ready to be interpolated into the top of
        every emitted SLURM batch script. Empty string when init is empty.
        """
        return "\n".join(self.get_scheduler_init())

    def get_activate_command(self, env_name: str) -> str:
        """Get the activate command for a specific environment."""
        mgr = self.get_env_manager()
        if mgr == "pip":
            return ""
        return f'{mgr} activate {env_name}'

    def get_env_python_command(self, env_name: str) -> str:
        """Bash that prints the Python executable of ``env_name`` via the
        configured env manager.

        Multi-env tools (a tool whose runtime dispatches to a *second* env's
        Python without activating it) must resolve that Python through the
        same manager the rest of the pipeline uses — never a hard-coded
        ``micromamba || mamba || conda`` fallback chain, which can pick a
        different manager than the active config. In pip mode there is no
        second env, so this resolves to the current interpreter.
        """
        mgr = self.get_env_manager()
        py = 'python -c "import sys; print(sys.executable)"'
        if mgr == "pip":
            return f'$({py})'
        return f'$({mgr} run -n {env_name} {py})'

    def get_module_load_line(self) -> str:
        """Get the full 'module load ...' line for SLURM scripts, or empty string if none."""
        modules = self.get_slurm_modules()
        if not modules:
            return ""
        return "module load " + " ".join(modules)

    def get_renderers_config(self) -> Dict[str, Any]:
        """
        Get renderers configuration section.

        Returns:
            Dictionary with 'streams' and 'tables' sub-dicts mapping
            format/name to renderer script path, or empty dict if not configured.
        """
        return self._config.get('renderers', {})

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

        config_name = os.path.basename(config_path)

        # Check if config already exists
        if os.path.exists(config_path) and not force:
            raise FileExistsError(
                f"{config_name} already exists at {config_path}. "
                f"Use force=True to overwrite, or use repull_config() to backup and update."
            )

        # Try to get URL from existing config, otherwise use default for the active variant
        variant = cls._variant or _autodetect_variant()
        default_url = f"https://raw.githubusercontent.com/gquargnali/biopipelines/main/config.{variant}.yaml"

        try:
            # If config exists, try to read repository URL from it
            if os.path.exists(config_path):
                import yaml
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
            print(f"Downloading {config_name} from {url}...")
            urllib.request.urlretrieve(url, config_path)
            print(f"Config saved to {config_path}")
            return config_path
        except Exception as e:
            raise urllib.error.URLError(f"Failed to download {config_name}: {e}")

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

        config_name = os.path.basename(config_path)

        # Try to pull from repository
        try:
            return cls.pull_config(force=False)
        except Exception as e:
            print(f"Warning: Could not pull config from repository: {e}")
            print(f"Please ensure {config_name} exists in the repository root.")
            raise FileNotFoundError(
                f"{config_name} not found and could not be downloaded. "
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
