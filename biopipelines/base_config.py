# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Base configuration class for all modeling tools.

Provides common interface and functionality for tool configuration,
validation, and integration with the pipeline system.
"""

import pandas as pd
import functools
import inspect
import os
import sys
import re
import json
from abc import ABC, abstractmethod
from typing import Dict, List, Any, Optional, TypeVar, overload, Type

try:
    from ._layout import INTERNAL_FOLDER
except ImportError:
    from _layout import INTERNAL_FOLDER


_SAFE_FREEFORM_RE = re.compile(r'^[^"`$\\]*$')


def _validate_freeform_string(field: str, value: Optional[str]) -> None:
    """
    Reject characters that would break double-quoted bash interpolation.

    Call from a tool's validate_params() on any user-supplied string that is
    later interpolated into a generated shell script (echo lines, CLI args,
    filenames built from the value, etc.). Characters " ` $ \\ would either
    terminate the surrounding quote context or trigger shell expansion.

    Values that are legitimately shell-expression-shaped (e.g. Panda filter
    expressions) must NOT use this helper — they need their own validator.
    """
    if value is None:
        return
    if not isinstance(value, str):
        raise ValueError(
            f"{field!r} must be a string, got {type(value).__name__}"
        )
    if not _SAFE_FREEFORM_RE.match(value):
        raise ValueError(
            f"{field!r}={value!r} contains one of \" ` $ \\ , which would "
            "break the generated shell script. Please remove these characters."
        )


def _escape_for_double_quotes(value: str) -> str:
    """Escape characters that retain special meaning inside a bash "..." context.

    Use this for filesystem paths and other framework-generated strings that
    are emitted into bash but would legitimately contain ``\\`` on Windows or
    ``$`` in some environments. User-supplied free-form strings should still
    go through ``_validate_freeform_string`` at construction time; this helper
    is only for emission-time escaping.

    Order matters: backslash must be escaped first so its replacement isn't
    re-processed by the subsequent rules.
    """
    return (
        value.replace("\\", "\\\\")
             .replace('"', '\\"')
             .replace("$", "\\$")
             .replace("`", "\\`")
    )


# CLI dialects `render_extra_args` can emit: "argparse" is `--key value`, "hydra" is `key=value`, which the whole RFdiffusion family needs because its entry points are hydra apps rather than argparse ones.
EXTRA_ARGS_DIALECTS = ("argparse", "hydra")

# The env the framework itself runs in, and the fallback when a config declares none.
FRAMEWORK_ENV = "biopipelines"

# A forwarded key becomes an upstream flag verbatim, so it is allowlisted rather than denylisted. Dotted paths and a leading "+" are here because a hydra override is spelled `denoiser.noise_scale_ca=1` or `+potentials.substrate=LIG`; those reach a constructor as `Tool(**{"denoiser.noise_scale_ca": 1})`, which Python accepts even though the key is not an identifier. An inner "-" is allowed too, because many upstream parsers spell flags `--msa-mode`.
_SAFE_EXTRA_KEY_RE = re.compile(r"^\+?[A-Za-z_][A-Za-z0-9_-]*(\.[A-Za-z_][A-Za-z0-9_-]*)*$")

# Stricter than _SAFE_FREEFORM_RE because LigandMPNN emits its command through `eval`, which re-parses
# the line and so strips one layer of quoting: a `;` in a forwarded value would start a new command.
_SHELL_METACHARACTERS_RE = re.compile(r"[;|&<>()\n\r]")


def _validate_extra_arg_scalar(field: str, value: Any) -> None:
    if value is None or isinstance(value, (bool, int, float)):
        return
    if isinstance(value, str):
        _validate_freeform_string(field, value)
        found = _SHELL_METACHARACTERS_RE.search(value)
        if found:
            raise ValueError(
                f"{field}={value!r} contains the shell metacharacter {found.group()!r}, which is never "
                "part of an upstream command-line value and would be re-interpreted as shell syntax."
            )
        return
    raise ValueError(
        f"{field}={value!r} of type {type(value).__name__} cannot be forwarded as an "
        "upstream command-line argument; only strings, numbers, booleans, None and "
        "flat lists of those can."
    )


def _validate_extra_arg(key: str, value: Any) -> None:
    """Shell-safety layer 4 for one forwarded kwarg, applied to the key and the value.

    Runs at construction time so an unusable value raises at the user's desk instead of reaching bash on a compute node.
    """
    if not _SAFE_EXTRA_KEY_RE.match(key):
        raise ValueError(
            f"forwarded parameter name {key!r} is not an identifier or a dotted "
            "path, so it cannot be rendered as an upstream flag."
        )
    if isinstance(value, (list, tuple)):
        for i, item in enumerate(value):
            if isinstance(item, (list, tuple)):
                raise ValueError(
                    f"{key}[{i}]: a nested list cannot be rendered as a flag value."
                )
            _validate_extra_arg_scalar(f"{key}[{i}]", item)
        return
    _validate_extra_arg_scalar(key, value)


def render_extra_args(extras: Dict[str, Any], dialect: str = "argparse") -> List[str]:
    """Render forwarded kwargs as raw upstream argv tokens.

    The single renderer for every forwarding tool, so 73 pipe scripts and their wrappers cannot each invent their own spelling. Underscores are preserved: both dialects in use take them verbatim.

    ``True`` becomes a bare flag under argparse and an explicit ``key=True`` under hydra, where a bare token is not a valid override. ``False`` and ``None`` are omitted entirely, since "off" is the upstream default the wrapper never touched. A list repeats an argparse flag once per element, but renders as a single hydra ``key=[a,b]``: a repeated hydra override keeps only the last value, which would drop elements silently.
    """
    if dialect not in EXTRA_ARGS_DIALECTS:
        raise ValueError(
            f"unknown extra-args dialect {dialect!r}; expected one of {EXTRA_ARGS_DIALECTS}"
        )
    tokens: List[str] = []
    for key, value in extras.items():
        if value is False or value is None:
            continue
        if dialect == "hydra":
            if isinstance(value, (list, tuple)):
                tokens.append(f"{key}=[{','.join(str(v) for v in value)}]")
            else:
                tokens.append(f"{key}={value}")
            continue
        if value is True:
            tokens.append(f"--{key}")
        elif isinstance(value, (list, tuple)):
            for item in value:
                tokens += [f"--{key}", str(item)]
        else:
            tokens += [f"--{key}", str(value)]
    return tokens


def _constructor_parameter_names(cls: type) -> List[str]:
    """The named parameters of a tool's own constructor, for the misspelling check."""
    try:
        parameters = inspect.signature(cls.__init__).parameters
    except (TypeError, ValueError):
        return []
    skip = (inspect.Parameter.VAR_KEYWORD, inspect.Parameter.VAR_POSITIONAL)
    return [n for n, p in parameters.items() if n != "self" and p.kind not in skip]


def _resolve_parameter_aliases(cls: type, kwargs: Dict[str, Any]) -> Dict[str, Any]:
    """Rewrite an alternative parameter spelling to the parameter's current name.

    Applied before the tool's own ``__init__`` binds anything, so an alias binds the real parameter instead of falling through to the leftover-kwargs logic — which would otherwise report it as a probable misspelling and, on a forwarding tool, put it on the upstream command line.

    A key the framework itself reads (`name`) is copied rather than moved: the tool gets the value under its new parameter name AND `BaseConfig` still sees `name`, so a tool that used to capture the job name keeps its old behaviour while the job name finally lands.
    """
    aliases = cls.PARAMETER_ALIASES
    if not aliases:
        return kwargs
    resolved = dict(kwargs)
    for old, new in aliases.items():
        if old not in resolved:
            continue
        reserved = old in contract_enforcement.RESERVED_KWARGS
        if new in resolved:
            # `name` alongside the tool's own parameter is the uniform spelling, not a clash: it is the job name and nothing else.
            if reserved:
                continue
            raise ValueError(
                f"{cls.TOOL_NAME}: {old}= is another spelling of {new}=; "
                f"pass one or the other, not both."
            )
        value = resolved[old] if reserved else resolved.pop(old)
        resolved[new] = value
        if old in cls.DEPRECATED_ALIASES:
            contract_enforcement.report(contract_enforcement.Violation(
                check="deprecated_alias",
                message=f"{cls.TOOL_NAME}: {old}= is a synonym for {new}= and will soon be deprecated.",
                hint=(f"The value was bound to {new}=, so this pipeline still runs as "
                      f"written. Silence this line with "
                      f"BIOPIPELINES_ENFORCE_DEPRECATED_ALIAS=off."),
            ))
    return resolved


def _alias_resolving_init(init):
    """Wrap a tool's ``__init__`` so `PARAMETER_ALIASES` is applied to its kwargs first."""
    @functools.wraps(init)
    def wrapper(self, *args, **kwargs):
        return init(self, *args, **_resolve_parameter_aliases(type(self), kwargs))

    wrapper._resolves_aliases = True
    return wrapper


try:
    from .config_manager import ConfigManager
    from .file_paths import Path
    from . import contract_enforcement
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from config_manager import ConfigManager
    from file_paths import Path
    import contract_enforcement



# Re-exported so the modules importing these from base_config keep working after the split.
try:
    from .data_containers import (  # noqa: F401
        TableMetadata,
        TableInfo,
        IndexedTableContainer,
        TableContainer,
        StreamContainer,
        resolve_table_reference,
    )
    from .outputs import StandardizedOutput, ToolOutput  # noqa: F401
except ImportError:
    from data_containers import (  # noqa: F401
        TableMetadata,
        TableInfo,
        IndexedTableContainer,
        TableContainer,
        StreamContainer,
        resolve_table_reference,
    )
    from outputs import StandardizedOutput, ToolOutput  # noqa: F401

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
    # Wrapper version. Bump in the subclass when the wrapper's external contract
    # changes (CLI flag mapping, generated bash, validate_params, output streams).
    # Independent of the upstream tool's version. The pre-commit hook enforces
    # that this is bumped whenever the wrapper or its pipe scripts change.
    TOOL_VERSION = "1.0"
    # The environment this tool builds and runs in. Declared here so no install script hides an env name in its body; only env_manager "pip" ever falls back to it (see _install_env).
    ENV_NAME = ""

    # True when the tool's user names its streams (Mock(streams={...}), Scripting), so no registry can know them.
    USER_STREAM_NAMES = False

    # "" = off. Set to a dialect from EXTRA_ARGS_DIALECTS to declare that this tool renders constructor kwargs it does not type onto its upstream command line, so an advanced user can reach an upstream option the wrapper never curated.
    FORWARD_UNKNOWN_KWARGS = ""

    # {alternative spelling: current parameter name}. A renamed parameter keeps its old spelling here so a pipeline written against it still binds the new parameter; a synonym a user reasonably expects (RFdiffusion3's `contigs`, the plural its three sibling tools take) is declared the same way. Resolved before the constructor binds its parameters — see _resolve_parameter_aliases.
    PARAMETER_ALIASES: Dict[str, str] = {}

    # The PARAMETER_ALIASES keys that are a retired spelling rather than a first-class synonym: only these report a deprecation, so a synonym stays silent.
    DEPRECATED_ALIASES: tuple = ()

    # Common path descriptors available to all tools
    pipeline_name = Path(lambda self: self._extract_pipeline_name())
    log_file = Path(lambda self: self._compute_log_file_path())

    # Sub-layout inside each tool's output_folder. Tools should route their
    # files through these descriptors instead of writing into output_folder
    # directly:
    #   configuration/ — config-time artifacts (input JSONs, .expected_outputs.json)
    #   execution/     — execution-time artifacts (logs, raw model dumps)
    #   tables/        — TableInfo CSVs and map_table CSVs
    #   _extras/       — catch-all for anything that doesn't fit the above
    # Per-stream folders live at <output_folder>/<stream_name>/ and are
    # obtained via self.stream_folder(name).
    configuration_folder = Path(lambda self: os.path.join(self.output_folder, "_configuration"))
    execution_folder = Path(lambda self: os.path.join(self.output_folder, "_execution"))
    tables_folder = Path(lambda self: os.path.join(self.output_folder, "tables"))
    extras_folder = Path(lambda self: os.path.join(self.output_folder, "_extras"))
    

    def __init_subclass__(cls, **kwargs):
        """Attribute every contract check raised inside get_output_files() to this tool.

        Wrapping here rather than at each of the ~90 call sites is what lets a DataStream stay constructible standalone while still naming its tool.
        """
        super().__init_subclass__(**kwargs)
        if cls.FORWARD_UNKNOWN_KWARGS:
            if cls.FORWARD_UNKNOWN_KWARGS not in EXTRA_ARGS_DIALECTS:
                raise ValueError(
                    f"{cls.TOOL_NAME}: FORWARD_UNKNOWN_KWARGS must be one of "
                    f"{EXTRA_ARGS_DIALECTS} or '' (off), got "
                    f"{cls.FORWARD_UNKNOWN_KWARGS!r}."
                )
            contract_enforcement.assert_can_forward(cls)
        if cls.PARAMETER_ALIASES or cls.DEPRECATED_ALIASES:
            known = set(_constructor_parameter_names(cls)) | contract_enforcement.RESERVED_KWARGS
            unbound = sorted(set(cls.PARAMETER_ALIASES.values()) - known)
            if unbound:
                raise ValueError(
                    f"{cls.TOOL_NAME}: PARAMETER_ALIASES points at {unbound}, which "
                    f"{cls.__name__}.__init__ does not name, so the alias would land "
                    f"back in **kwargs."
                )
            undeclared = sorted(set(cls.DEPRECATED_ALIASES) - set(cls.PARAMETER_ALIASES))
            if undeclared:
                raise ValueError(
                    f"{cls.TOOL_NAME}: DEPRECATED_ALIASES names {undeclared}, which "
                    f"PARAMETER_ALIASES does not map to a current parameter."
                )
            if not getattr(cls.__init__, "_resolves_aliases", False):
                cls.__init__ = _alias_resolving_init(cls.__init__)
        produce_outputs = cls.__dict__.get("get_output_files")
        if produce_outputs is None or getattr(produce_outputs, "_contract_context", False):
            return

        @functools.wraps(produce_outputs)
        def wrapper(self, *args, **kwargs):
            with contract_enforcement.tool_context(
                self.TOOL_NAME, user_stream_names=type(self).USER_STREAM_NAMES
            ):
                return produce_outputs(self, *args, **kwargs)

        wrapper._contract_context = True
        cls.get_output_files = wrapper

    @overload
    def __new__(cls: Type[T], *args, **kwargs) -> T:
        """
        Type hint overload for IDE autocomplete support.

        This tells type checkers that Tool(...) returns a Tool instance,
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
    def _install_env(cls, env_manager: str, env_name: Optional[str] = None) -> str:
        """Return the env name this tool's install script must build.

        Resolution order: an explicit ``env_name`` from the caller, then the tool's ``environments:`` entry in the active config, then raise. The config entry is what ``_load_environments`` activates at run time, so reading it here is what keeps install and runtime pointed at the same env: an install script that hard-codes its own name builds one env while the pipeline activates another, and the mismatch only surfaces later on a compute node.

        No fallback to ``ENV_NAME`` outside pip mode is deliberate. A tool with no entry under the active variant is not set up there, and raising names the missing ``environments:`` entry instead of leaving a job to fail on an env that was never created.

        ``env_manager: pip`` is exempt from the raise, because a laptop config declares no ``environments:`` at all: it falls back to the tool's own name as a label, and :meth:`_env_install_block` refuses to build a conda-style env there rather than emitting ``pip env create`` and ``pip run -n``, which are not pip commands and could only ever fail.
        """
        if env_name:
            return env_name

        config_manager = ConfigManager()
        configured = config_manager.get_environment(cls.TOOL_NAME)
        if isinstance(configured, list):
            # A multi-env tool lists its primary env first; the rest are its own business.
            configured = configured[0] if configured else None
        if configured:
            return configured

        if env_manager == "pip":
            # Exempt from the raise, not from the config: a laptop variant declares no environments:.
            return cls.ENV_NAME

        config_name = os.path.basename(config_manager._get_config_path())
        variant = config_manager.get_variant()
        if configured == "":
            raise ValueError(
                f"{cls.TOOL_NAME}: environments.{cls.TOOL_NAME} is empty in "
                f"{config_name} (variant '{variant}'), which means an image supplies "
                f"the interpreter and there is no environment to build. Installing "
                f"{cls.TOOL_NAME} under this variant is not supported: use the "
                f"configured image, or give environments.{cls.TOOL_NAME} a real "
                f"environment name."
            )
        raise ValueError(
            f"{cls.TOOL_NAME}: no environment is configured for it in {config_name} "
            f"(variant '{variant}'), so there is no name to install into. Add an "
            f"environments: entry -- environments.{cls.TOOL_NAME}: \"<env name>\" -- "
            f"or install {cls.TOOL_NAME} under a variant that configures it."
        )

    @staticmethod
    def _env_exists_check(env_name: str, env_manager: str) -> str:
        """Return a bash test expression that's true iff the env exists.

        Probes by trying to launch python inside the env rather than parsing
        ``env list`` output — the latter's column layout (leading whitespace,
        active-env ``*`` prefix, name-vs-path columns) varies across mamba,
        micromamba, and conda, and across hosts. The probe form is uniform.

        Use as: ``if {check}; then ...`` or ``if ! {check}; then ...``.
        """
        if env_manager == "venv":
            return f'"{ConfigManager().get_venv_path(env_name)}/bin/python" -c "" >/dev/null 2>&1'
        return f'{env_manager} run -n {env_name} python -c "" >/dev/null 2>&1'

    @staticmethod
    def _env_run(env_name: str, env_manager: str) -> str:
        """Prefix that runs the next command inside ``env_name``.

        ``<mgr> run -n <env>`` for the conda family; venv has no such
        subcommand, so its env is addressed by path instead.

        Use as: ``f"{cls._env_run(env, mgr)}pip install ..."`` — note no space,
        the venv form resolves to ``<path>/bin/pip``. Every command the tools
        wrap this way (python, pip, mkdssp, xtb, …) lives in that bin/.
        """
        if env_manager == "venv":
            return f'{ConfigManager().get_venv_path(env_name)}/bin/'
        return f'{env_manager} run -n {env_name} '

    @classmethod
    def _env_remove_block(cls, env_name: str, env_manager: str) -> str:
        """Return bash that removes an existing env if present (idempotent).

        Pair with ``_env_install_block`` under ``force_reinstall=True``: the
        helper emits ``<env_manager> env create -f ...``, which fails when the
        env already exists. Calling this beforehand makes the reinstall path
        actually work. The existence guard avoids non-zero exit codes from
        ``env remove`` on a non-existent env aborting the script under set -e.
        """
        # A tool mapped to the framework's own env (config.cluster.yaml does this for PLIP and Gnina)
        # would otherwise delete the env every other step activates.
        if env_name == FRAMEWORK_ENV:
            raise ValueError(
                f"{cls.TOOL_NAME}.install(force_reinstall=True) resolves to the framework env "
                f"{FRAMEWORK_ENV!r}, which every other step activates, so removing it would break the "
                f"whole pipeline. Its `environments:` entry points there deliberately (the tool needs "
                f"only the shared deps), so there is nothing tool-specific to rebuild. Reinstall "
                f"without force, or rebuild {FRAMEWORK_ENV!r} yourself from "
                f"environments/{FRAMEWORK_ENV}.yaml."
            )
        check = cls._env_exists_check(env_name, env_manager)
        # venv has no registry and no `env remove`; the env IS its directory.
        if env_manager == "venv":
            path = ConfigManager().get_venv_path(env_name)
            return (
                f'if [ -d "{path}" ]; then\n'
                f'    echo "Removing existing {env_name} env for reinstall"\n'
                f'    rm -rf "{path}"\n'
                f'fi'
            )
        return (
            f'if {check}; then\n'
            f'    echo "Removing existing {env_name} env for reinstall"\n'
            f'    {env_manager} env remove -n {env_name} -y\n'
            f'fi'
        )

    @staticmethod
    def _link_conda_binaries(env_name: str, env_manager: str, binaries) -> str:
        """Symlink conda-provided executables into a venv shim's bin/.

        A venv built by a conda env's python inherits that env's *Python
        packages* through ``--system-site-packages``, but not its executables:
        ``<venv>/bin`` holds only python and pip. A tool whose conda package
        ships a binary (vina, obabel, mkdssp, …) must link it in, or the binary
        is missing from PATH once the shim is activated.

        No-op unless the env manager is venv — the conda family puts binaries on
        PATH itself.
        """
        if env_manager != "venv":
            return ""
        cm = ConfigManager()
        venv_bin = f"{cm.get_venv_path(env_name)}/bin"
        conda_bin = f"{cm.get_conda_env_root()}/{env_name}/bin"
        names = " ".join(binaries)
        return (
            f'# venv shims inherit packages, not executables.\n'
            f'for _bin in {names}; do\n'
            f'    if [ -x "{conda_bin}/$_bin" ]; then\n'
            f'        ln -sf "{conda_bin}/$_bin" "{venv_bin}/$_bin"\n'
            f'    fi\n'
            f'done\n'
        )

    @staticmethod
    def _env_create_name_flag(env_name: str, env_manager: str) -> str:
        """Return the ``-n <env>`` flag for an ``env create``, or "" under pip.

        A spec yaml carries its own ``name:`` key, which is what ``env create`` would name the env; passing ``-n`` overrides it so a redirected ``environments:`` entry creates the env the pipeline actually activates. Empty under pip, whose emitted commands are inert and are left exactly as they were.
        """
        if env_manager == "pip":
            return ""
        return f" -n {env_name}"

    @staticmethod
    def _env_install_block(env_name: str, env_manager: str, biopipelines: str) -> str:
        """Return bash that creates an env from environments/<env>.<variant>.yaml
        and (if present) installs the matching pip requirements on top.

        For both the conda YAML and the pip txt, falls back to the variant-less
        filename (``<env>.yaml`` / ``<env>.pip.txt``) when the variant-specific
        one is absent — envs that resolve identically across variants are
        stored under a single file.

        Phased pip installs: if numbered pip files exist
        (``<env>.pip.<variant>.<N>.txt`` for N=1,2,...), they are installed in
        numeric order as separate ``pip install -r`` calls. Use this when one
        package's setup.py imports another at module level (e.g. openfold
        importing torch) so build isolation can't see it via metadata. A pip
        file whose first non-empty line is ``# bp:no-build-isolation`` is
        installed with ``--no-build-isolation``, letting it pick up packages
        already installed in the env from an earlier phase. When no numbered
        files exist, falls back to the single-file form
        ``<env>.pip.<variant>.txt`` / ``<env>.pip.txt``.

        The variant is read from the active ConfigManager so each install script
        only needs to call this helper with the env name.
        """
        import glob, re
        if env_manager == "pip":
            # `pip env create` and `pip run -n` are not pip commands; emitting them only failed later.
            message = (
                f"Cannot build the {env_name} environment: this variant sets machine.env_manager to"
                " pip, which has no environment manager. Install the tool's dependencies into the"
                " active interpreter, or use a variant whose env_manager is mamba, micromamba,"
                " conda or venv."
            )
            return f'echo "{message}" >&2' + chr(10) + "exit 1" + chr(10)
        variant = ConfigManager().get_variant()
        env_dir = f"{biopipelines}/environments"
        env_yaml_variant = f"{env_dir}/{env_name}.{variant}.yaml"
        env_yaml_shared = f"{env_dir}/{env_name}.yaml"
        if env_manager == "venv":
            import yaml as _yaml
            src = env_yaml_variant if os.path.isfile(env_yaml_variant) else env_yaml_shared
            try:
                with open(src, encoding="utf-8") as f:
                    _y = _yaml.safe_load(f) or {}
                deps = _y.get("dependencies", []) or []
                channels = _y.get("channels", []) or []
            except OSError:
                deps, channels = [], []
            extra = [d for d in deps if isinstance(d, str)
                     and re.split(r"[=<>\s]", d, 1)[0] not in ("python", "pip")]
            cm = ConfigManager()
            conda_for_venv = cm.get_conda_for_venv()
            if extra and not conda_for_venv:
                raise ValueError(
                    f"{env_name}: env_manager 'venv' cannot install conda packages "
                    f"{extra} from {os.path.basename(src)}. Move them to a pip file, "
                    f"or supply them via the tool's EDF image."
                )
            venv_path = cm.get_venv_path(env_name)
            # A package with no wheel for this platform compiles, and a login node
            # may ship no python3-dev. Point CPATH at borrowed headers if the
            # config provides them; harmless when every dependency is a wheel.
            hdrs = (cm.get_folder_config().get("infrastructure") or {}).get("python_headers", "")
            cpath = f'export CPATH="{cm._resolve_folder_template(hdrs)}"\n' if hdrs else ""
            if extra:
                # Conda packages with no wheel (pymol, mkdssp, …). Build a conda
                # env for them, then a venv created BY its python: that inherits
                # them via --system-site-packages and supplies the bin/activate
                # the framework sources. Binaries are not inherited — a tool
                # needing one symlinks it in its own install script.
                mamba = cm.get_conda_for_venv()
                prefix = f"{cm.get_conda_env_root()}/{env_name}"
                chan = " ".join(f"-c {c}" for c in channels) or "-c conda-forge"
                pyspec = next((d for d in deps if isinstance(d, str)
                               and re.split(r"[=<>\s]", d, 1)[0] == "python"), "python")
                env_block = (
                    cpath +
                    f'if [ ! -x "{prefix}/bin/python" ]; then\n'
                    f'    "{mamba}" create -y -q -p "{prefix}" {chan} "{pyspec}" '
                    f'{" ".join(chr(34) + d + chr(34) for d in extra)}\n'
                    f'fi\n'
                    f'if [ ! -d "{venv_path}" ]; then\n'
                    f'    "{prefix}/bin/python" -m venv --system-site-packages "{venv_path}"\n'
                    f'fi\n'
                    f'source "{venv_path}/bin/activate"\n'
                )
            else:
                env_block = (
                    cpath +
                    f'if [ ! -d "{venv_path}" ]; then\n'
                    f'    python -m venv --system-site-packages "{venv_path}"\n'
                    f'fi\n'
                    f'source "{venv_path}/bin/activate"\n'
                )
        else:
            named = BaseConfig._env_create_name_flag(env_name, env_manager)
            # An env may be built by an upstream installer rather than from a spec (RF3's 'modelforge'), and `env create -f` on an absent file only fails.
            env_block = (
                f'if [ -f "{env_yaml_variant}" ]; then\n'
                f'    {env_manager} env create -f "{env_yaml_variant}"{named} -y\n'
                f'elif [ -f "{env_yaml_shared}" ]; then\n'
                f'    {env_manager} env create -f "{env_yaml_shared}"{named} -y\n'
                f'else\n'
                f'    echo "No environments/{env_name}.{variant}.yaml or {env_name}.yaml;'
                f' assuming the {env_name} environment already exists." >&2\n'
                f'fi\n'
            )

        # Discover numbered pip files in phase order. Variant-specific
        # (<env>.pip.<variant>.<N>.txt) takes precedence; if none exist we
        # fall back to the shared form (<env>.pip.<N>.txt) so environments
        # that resolve identically across variants don't need duplicated files.
        def _numbered(pattern_glob, pattern_re):
            return sorted(
                ((int(m.group(1)), p.replace("\\", "/")) for p, m in (
                    (p, re.search(pattern_re, p))
                    for p in glob.glob(pattern_glob)
                ) if m),
                key=lambda t: t[0],
            )

        numbered = _numbered(
            f"{env_dir}/{env_name}.pip.{variant}.*.txt",
            rf"{re.escape(env_name)}\.pip\.{re.escape(variant)}\.(\d+)\.txt$",
        )
        if not numbered:
            numbered = _numbered(
                f"{env_dir}/{env_name}.pip.*.txt",
                rf"{re.escape(env_name)}\.pip\.(\d+)\.txt$",
            )

        # venv's env_block already activated the env; the conda family has not.
        pip = "pip" if env_manager == "venv" else f'{env_manager} run -n {env_name} pip'

        def _pip_line(path: str) -> str:
            # Honor leading `# bp:*` markers. `# bp:no-build-isolation` opts
            # the phase out of PEP 517 isolation (so it can see packages
            # already installed in the env by an earlier phase or the conda
            # layer). `# bp:no-deps` tells pip not to resolve transitive
            # dependencies (used when the conda layer already provides them
            # and pip would otherwise try to install from sdist). Both markers
            # must appear at the top of the file, in any order, before the
            # first requirement; the first non-comment line ends the header.
            flags = []
            try:
                with open(path, encoding="utf-8") as f:
                    for line in f:
                        s = line.strip()
                        if not s:
                            continue
                        if s == "# bp:no-build-isolation":
                            flags.append("--no-build-isolation")
                            continue
                        if s == "# bp:no-deps":
                            flags.append("--no-deps")
                            continue
                        if s.startswith("#"):
                            continue
                        break
            except OSError as exc:
                # Unreadable marker header just means no flags; warn so a botched
                # install isn't silently missing --no-deps/--no-build-isolation.
                print(f"Warning: could not read pip marker header {path!r}: {exc}",
                      file=sys.stderr)
            flag_str = (" " + " ".join(flags)) if flags else ""
            return f'{pip} install{flag_str} -r "{path}"'

        if numbered:
            pip_block = "\n".join(_pip_line(p) for _, p in numbered)
            return env_block + pip_block

        pip_txt_variant = f"{env_dir}/{env_name}.pip.{variant}.txt"
        pip_txt_shared = f"{env_dir}/{env_name}.pip.txt"
        return env_block + (
            f'if [ -f "{pip_txt_variant}" ]; then\n'
            f'    {pip} install -r "{pip_txt_variant}"\n'
            f'elif [ -f "{pip_txt_shared}" ]; then\n'
            f'    {pip} install -r "{pip_txt_shared}"\n'
            f'fi'
        )

    @classmethod
    def _container_pull_block(cls, folders: Dict[str, str], force_reinstall: bool = False) -> str:
        """Bash that pulls the tool's container image to the configured path, if a container is configured.

        Returns "" (nothing to pull) unless BOTH sides are present: a DESTINATION ``.sif`` path configured for this tool (``container:<TOOL_NAME>`` in folders) AND a SOURCE for ``cls.TOOL_NAME`` in ``environments/_containers.yaml``. Two source forms are supported: a registry URI (``docker://...``), fetched with ``<executor> pull``; or a direct ``http(s)://`` link to a ``.sif``, fetched with ``wget`` (it is already an image, nothing to convert). Skips when the destination exists unless ``force_reinstall``. A container tool calls this from its ``_install_script`` alongside (or instead of) the conda-env block; the same destination path then drives execution via ``container_prefix()``.
        """
        image = folders.get(f"container:{cls.TOOL_NAME}")
        if not image:
            return ""
        executor = ConfigManager().get_container_executor()
        if executor in (None, "none", ""):
            return ""
        source = cls._container_source(folders)
        if not source:
            return f'echo "WARNING: {cls.TOOL_NAME} has a container path configured but no source in environments/_containers.yaml; cannot pull."\n'
        # A registry URI is converted by the executor; a plain http(s) .sif is already an image — just fetch it.
        # `pull` refuses to overwrite an existing file, so force_reinstall needs --force (wget -O already overwrites).
        if source.startswith(("http://", "https://")):
            fetch = f'wget -O "{image}" "{source}"\n'
        else:
            force_flag = "--force " if force_reinstall else ""
            fetch = f'{executor} pull {force_flag}"{image}" "{source}"\n'
        pull = f'mkdir -p "$(dirname "{image}")"\n{fetch}'
        if force_reinstall:
            body = pull
        else:
            body = (
                f'if [ -f "{image}" ]; then\n'
                f'    echo "Container image already present: {image}"\n'
                f'else\n'
                + "".join(f"    {line}\n" for line in pull.splitlines())
                + 'fi\n'
            )
        return f'echo "=== Pulling {cls.TOOL_NAME} container image ==="\n{body}'

    @classmethod
    def _container_source(cls, folders: Dict[str, str]) -> Optional[str]:
        """The pull source URI for this tool from environments/_containers.yaml (keyed by TOOL_NAME), or None."""
        biopipelines = folders.get("biopipelines", "")
        path = os.path.join(biopipelines, "environments", "_containers.yaml")
        if not os.path.isfile(path):
            return None
        import yaml
        with open(path) as f:
            data = yaml.safe_load(f) or {}
        return data.get(cls.TOOL_NAME)

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
        # Probed with the configured manager, not a fixed "mamba", because _install_env resolves per manager.
        if cls._install_script({}, ConfigManager().get_env_manager()) is None:
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
        # Internal tools auto-register and execute but are hidden from public numbering/layout.
        self.internal = bool(kwargs.pop("_internal", False))

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

        # Whatever is left after the subclass bound its named parameters is a framework key, a typo, or — on a forwarding tool — an upstream option the wrapper does not type.
        dialect = type(self).FORWARD_UNKNOWN_KWARGS
        self.extra_args = {
            k: v for k, v in kwargs.items()
            if k not in contract_enforcement.RESERVED_KWARGS
        } if dialect else {}

        contract_enforcement.check_kwargs(
            self.TOOL_NAME, kwargs,
            known=_constructor_parameter_names(type(self)),
            forwards=bool(dialect),
        )
        for key, value in self.extra_args.items():
            _validate_extra_arg(key, value)

        # Validate configuration
        self.validate_params()

    def extra_args_tokens(self) -> List[str]:
        """Raw upstream argv tokens for this tool's forwarded kwargs; empty when it does not forward."""
        extras = getattr(self, "extra_args", None)
        if not extras:
            return []
        return render_extra_args(extras, type(self).FORWARD_UNKNOWN_KWARGS)

    def extra_args_bash_tokens(self) -> List[str]:
        """`extra_args_tokens`, each wrapped in double quotes for a bash array element or command line."""
        return ['"' + t + '"' for t in self.extra_args_tokens()]

    def extra_args_bash(self) -> str:
        """The forwarded arguments as one bash-ready fragment, or '' when there are none."""
        return " ".join(self.extra_args_bash_tokens())

    def extra_args_echo(self) -> str:
        """A bash `echo` naming the forwarded arguments, so they are visible in the log rather than silent.

        The tokens carry no `" ` $ \\` — `_validate_extra_arg` rejected those at construction — so they are safe unquoted inside the echo's own double quotes.
        """
        tokens = self.extra_args_tokens()
        if not tokens:
            return ""
        return (f'echo "Forwarding unrecognized argument(s) to {self.TOOL_NAME}: '
                f'{" ".join(tokens)}"\n')

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

        if env_config == "":
            # Explicit empty: the tool's EDF image supplies the interpreter.
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

    def uses_container(self) -> bool:
        """True iff a container image is configured for this tool."""
        folders = getattr(self, 'folders', {}) or {}
        return f"container:{self.TOOL_NAME}" in folders

    def container_prefix(self) -> str:
        """Shell prefix that runs the next command inside this tool's container.

        Returns '<executor> exec --nv -B <path1>,<path2>,... <image> ' when a
        container is configured, else ''. Bind mounts come from the pipeline's
        resolved folders so outputs and caches are visible inside the image.
        Tools build commands as f"{self.container_prefix()}python foo.py ...".
        """
        if not self.uses_container():
            return ""
        folders = self.folders
        image = folders[f"container:{self.TOOL_NAME}"]
        binds = folders.get("__container_binds__") or []
        executor = ConfigManager().get_container_executor()
        bind_arg = f"-B {','.join(binds)} " if binds else ""
        return f"{executor} exec --nv {bind_arg}{image} "

    def warn_container_unsupported(self) -> str:
        """Bash that warns when a container is configured for a tool that can't use it.

        Some tools run their model in-process alongside biopipelines imports, or
        dispatch host-env sub-processes, so they have no single binary boundary to
        wrap and ignore a configured container. They call this in generate_script
        so a user who sets one gets feedback instead of silent env-mode execution.
        Returns '' when no container is configured.
        """
        if not self.uses_container():
            return ""
        return (
            f'echo "WARNING: a container is configured for {self.TOOL_NAME}, but this '
            f'tool does not support container execution; running in environment mode."\n'
        )

    def _resolve_env_placeholders(self, env_name: str) -> str:
        """Substitute ``<folder>`` placeholders in an env value against the
        resolved folder map.

        Lets the config target a path-based (conda ``-p``) env by reusing a
        folder key, e.g. ``AF2BIND: "<AlphaFold>/colabfold-conda"`` ->
        ``/.../colabfold/colabfold-conda``. Bare names pass through untouched.
        Raises if a placeholder names an unknown folder key — silently leaving
        an unresolved ``<...>`` would produce a broken ``activate`` at runtime.
        """
        if not env_name or "<" not in env_name:
            return env_name
        folders = getattr(self, "folders", {}) or {}
        import re as _re
        def _sub(m):
            key = m.group(1)
            if key not in folders:
                raise ValueError(
                    f"Environment '{env_name}' for {self.TOOL_NAME} references "
                    f"unknown folder placeholder '<{key}>'. Known: {sorted(folders)}"
                )
            return folders[key]
        return _re.sub(r"<([^<>]+)>", _sub, env_name)

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

        # Source resolve_stream_item.sh for runtime file resolution
        resolve_source = ""
        pipe_scripts = getattr(self, 'folders', {}).get("pipe_scripts", "")
        if pipe_scripts:
            resolve_sh = os.path.join(pipe_scripts, "resolve_stream_item.sh")
            resolve_source = f'source "{resolve_sh}"\n'

        # pip mode: no environment activation needed
        if config_manager.get_env_manager() == "pip":
            return f"# pip mode: no environment activation needed\n{resolve_source}\n"

        # An empty `environments:` entry alongside an `edf:` one means the image
        # supplies the interpreter on PATH, so there is nothing to activate.
        if not self.environments and config_manager.get_edf(self.TOOL_NAME):
            return (f"# environment provided by the tool's EDF image\n"
                    f'echo "Python: $(which python)"\n{resolve_source}\n')

        # Container mode still activates a host env: the tool's configured env
        # (or the biopipelines fallback from _load_environments) handles any
        # host-side helper scripts that run outside container_prefix.
        if name is not None:
            env_name = name
        else:
            if index >= len(self.environments):
                raise IndexError(
                    f"Environment index {index} out of range for {self.TOOL_NAME}. "
                    f"Available environments: {self.environments}"
                )
            env_name = self.environments[index]

        # Resolve <folder> placeholders in the env value against the resolved
        # folder map, so the config can point a tool at a path-based (conda
        # ``-p``) env — e.g. ``AF2BIND: "<AlphaFold>/colabfold-conda"``. A bare
        # name (no ``<...>``) passes through unchanged. mamba/conda accept both
        # a registered name and an absolute prefix path after ``activate``;
        # prefix envs are NOT name-resolvable, so the placeholder form is the
        # only way to target them from the env map.
        env_name = self._resolve_env_placeholders(env_name)

        # On Colab the `biopipelines` env is never created — its deps live in
        # base Python (pip install -e). We don't `activate` it, but if another
        # per-tool env is currently active (a multi-phase tool that ran, e.g.,
        # foundry/diffdock inference then switches back to biopipelines for
        # post-processing) its python is still on PATH and lacks the base deps
        # (BioPython, pandas). Deactivate to fall back to base Python before
        # running biopipelines-env work. Harmless when nothing is active.
        if config_manager.get_scheduler() == "colab" and env_name == "biopipelines":
            shell_hook = config_manager.get_shell_hook_command()
            return (
                f"# Colab: biopipelines deps in base Python; deactivate any active env\n"
                f"{shell_hook}\n"
                f"micromamba deactivate 2>/dev/null || true\n"
                f"{resolve_source}\n"
            )

        shell_hook = config_manager.get_shell_hook_command()
        activate_cmd = config_manager.get_activate_command(env_name)

        # venv sets VIRTUAL_ENV, not the CONDA_* pair.
        if config_manager.get_env_manager() == "venv":
            env_report = 'echo "Environment: $VIRTUAL_ENV"'
        else:
            env_report = 'echo "Environment: $CONDA_DEFAULT_ENV"\necho "Location: $CONDA_PREFIX"'

        return f"""# Activate environment: {env_name}
echo "=== Activating Environment ==="
echo "Requested: {env_name}"
{shell_hook}
{activate_cmd}
{env_report}
echo "Python: $(which python)"
echo "Python version: $(python --version 2>&1)"
echo "=============================="
{resolve_source}
"""

    def generate_filtered_map_table_block(
        self,
        ds_json: str,
        out_csv: str,
        columns: Optional[List[str]] = None,
        required_columns: Optional[List[str]] = None,
    ) -> str:
        """
        Bash block that materializes a filtered map_table under the biopipelines env.

        Projects the stream's map_table to its expanded ids (honoring an upstream
        id filter) and writes ``out_csv``. The materializer imports biopipelines,
        so it must run in the biopipelines env — not the tool env. Emit this
        before activating the tool env for the consumer step.
        """
        from .biopipelines_io import Resolve

        block = "# Materialize filtered input CSV (biopipelines env)\n"
        block += self.activate_environment(name="biopipelines")
        block += Resolve.filtered_map_table(
            ds_json, out_csv, columns=columns, required_columns=required_columns
        ) + "\n"
        block += """if [ $? -ne 0 ]; then
    echo "Error: could not materialize filtered input CSV"
    exit 1
fi

"""
        return block

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

    def stream_folder(self, name: str) -> str:
        """Path to the subfolder that holds files for a declared output stream.

        Tools should emit per-ID stream files (e.g. ``<id>.pdb``) into
        ``self.stream_folder(stream_name)`` rather than directly into
        ``self.output_folder``. The folder is not created here — the tool
        (or its pipe script) is responsible for calling ``os.makedirs`` when
        it actually writes files.
        """
        if not name or not isinstance(name, str):
            raise ValueError(f"stream name must be a non-empty string, got {name!r}")
        return os.path.join(self.output_folder, name)

    def stream_path(self, name: str, *parts: str) -> str:
        """Join ``parts`` under ``self.stream_folder(name)``.

        Prefer this over ad-hoc ``os.path.join(self.stream_folder(name), ...)``
        for files that live inside a stream folder (content tables, per-stream
        artefacts like sequence logos).
        """
        return os.path.join(self.stream_folder(name), *parts)

    def stream_map_path(self, name: str) -> str:
        """Canonical path for a stream's ``map_table`` CSV.

        The map_table is part of the stream, so it lives inside the stream
        folder (next to the per-ID files), not under ``tables/``. Tables
        folder is reserved for standalone ``TableInfo`` outputs that are
        not tied to a specific stream (e.g. Boltz2's ``confidence.csv``).

        For content-bearing streams whose map_table IS the content table
        (Sequence, Ligand, CompoundLibrary), prefer
        ``self.stream_path(name, f"{name}.csv")`` and point both the
        ``DataStream.map_table`` and the ``TableInfo.path`` at it.
        """
        if not name or not isinstance(name, str):
            raise ValueError(f"stream name must be a non-empty string, got {name!r}")
        return self.stream_path(name, f"{name}_map.csv")

    def table_path(self, name: str) -> str:
        """Canonical path for a ``TableInfo`` CSV (``tables/<name>.csv``).

        Prefer this over ad-hoc ``os.path.join(self.output_folder, ...)``
        so the layout stays consistent across tools.
        """
        if not name or not isinstance(name, str):
            raise ValueError(f"table name must be a non-empty string, got {name!r}")
        return os.path.join(self.tables_folder, f"{name}.csv")

    def configuration_path(self, *parts: str) -> str:
        """Join ``parts`` under ``self.configuration_folder``.

        Prefer this over ad-hoc ``os.path.join(self.configuration_folder, ...)``
        for config-time artifacts (input JSONs/YAMLs, ``.expected_outputs.json``).
        """
        return os.path.join(self.configuration_folder, *parts)

    def execution_path(self, *parts: str) -> str:
        """Join ``parts`` under ``self.execution_folder``.

        Prefer this over ad-hoc ``os.path.join(self.execution_folder, ...)``
        for execution-time artifacts (raw model dumps, per-run working dirs).
        """
        return os.path.join(self.execution_folder, *parts)

    def extras_path(self, *parts: str) -> str:
        """Join ``parts`` under ``self.extras_folder``.

        Prefer this over ad-hoc ``os.path.join(self.extras_folder, ...)``
        for ancillary files (session dumps, debug aids) that don't fit under
        configuration/, execution/, tables/, or a stream folder.
        """
        return os.path.join(self.extras_folder, *parts)

    def pipe_script_path(self, *parts: str) -> str:
        """Join ``parts`` under the repo's ``pipe_scripts/`` directory.

        Prefer this over ad-hoc ``os.path.join(self.folders["pipe_scripts"], ...)``
        for resolving the path of a runtime helper (``pipe_*.py``).
        """
        return os.path.join(self.folders["pipe_scripts"], *parts)

    def _materialize_output_layout(self, output_files: Dict[str, Any]) -> None:
        """Create the full sub-layout under ``output_folder`` at config time.

        Called by the pipeline once the tool's ``get_output_files()`` has
        returned, so tool authors never need to ``mkdir`` anything — the
        four canonical sub-dirs (``configuration/``, ``execution/``,
        ``tables/``, ``_extras/``) plus one folder per declared stream are
        materialized in a single pass. Idempotent.

        Streams are identified by duck-typing (``.map_table`` attribute),
        same rule used by ``get_id_provenance``.
        """
        if not self.output_folder:
            return
        canonical = (self.configuration_folder, self.execution_folder,
                     self.tables_folder, self.extras_folder)
        for folder in canonical:
            os.makedirs(folder, exist_ok=True)

        if not isinstance(output_files, dict):
            return
        reserved = {"tables", "output_folder", "filter_metadata", "streams"}
        for name, value in output_files.items():
            if name in reserved:
                continue
            if hasattr(value, "map_table"):  # DataStream duck-type
                os.makedirs(self.stream_folder(name), exist_ok=True)

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

    def get_id_provenance(self) -> Dict[str, Dict[str, str]]:
        """
        Get provenance mapping for this tool's output IDs.

        Reads the map_table CSV of every output stream and extracts
        provenance columns ({stream_name}.id) to build child->parent
        mappings. Skips axes already found in an earlier stream.

        Returns:
            Dict mapping axis_name -> {output_id: parent_id}.
            Empty dict if no provenance data is available.
        """
        from .datastream import DataStream

        output_files = self.get_output_files()

        result = {}
        seen_tables = set()

        for key, val in output_files.items():
            if not isinstance(val, DataStream):
                continue
            if not val.map_table or not os.path.exists(val.map_table):
                continue
            if val.map_table in seen_tables:
                continue
            seen_tables.add(val.map_table)

            try:
                df = pd.read_csv(val.map_table)
            except Exception:
                continue
            if "id" not in df.columns:
                continue

            for col in df.columns:
                if col.endswith(".id") and col != "id":
                    axis_name = col[:-3]  # "structures.id" -> "structures"
                    if axis_name in result:
                        continue  # already found from another stream
                    mapping = {}
                    for _, row in df.iterrows():
                        output_id = str(row["id"])
                        parent_id = str(row[col])
                        if parent_id and parent_id != "nan":
                            mapping[output_id] = parent_id
                    if mapping:
                        result[axis_name] = mapping

        return result

    def get_expected_output_paths(self) -> Dict[str, List[str]]:
        """
        Get expected output file paths without validating existence.
        
        Default implementation just calls get_output_files().
        Override if file existence validation is problematic.
        
        Returns:
            Dictionary mapping output type to expected file paths
        """
        return self.get_output_files()
    
    def set_pipeline_context(self, pipeline_ref, execution_order: int, output_folder: str,
                             suffix: str = "", public_step: int = None, internal_order: int = None):
        """Set context when added to pipeline."""
        self.pipeline_ref = pipeline_ref
        self.execution_order = execution_order
        self.output_folder = output_folder
        self.suffix = suffix
        self.public_step = public_step
        self.internal_order = internal_order
        # Single source of truth for this tool's RunTime/<basename>.sh and Logs/<basename>.log.
        # Internal tools nest under a .internal/ subdir of RunTime and Logs.
        if self.internal:
            self.script_basename = os.path.join(INTERNAL_FOLDER, f"{internal_order:03d}_{self.TOOL_NAME}")
        elif suffix:
            self.script_basename = f"{public_step:03d}_{self.TOOL_NAME}_{suffix}"
        else:
            self.script_basename = f"{public_step:03d}_{self.TOOL_NAME}"
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

            if isinstance(obj, IndexedTableContainer):
                return obj.to_dict()
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
            self.folders.get("pipe_scripts", "pipe_scripts"),
            "pipe_check_completion.py"
        )

        # Write expected outputs to JSON file at configuration time (not execution time)
        # Wrap in metadata envelope so Load can load directly from the tool folder
        expected_outputs_wrapped = {
            "tool_name": self.TOOL_NAME,
            "tool_class": self.__class__.__name__,
            "output_structure": json_safe_outputs,
        }
        # .expected_outputs.json is a top-level status/manifest file; it
        # stays at the tool-folder root alongside .log so users and Load()
        # can find it with a single ls.
        expected_outputs_file = os.path.join(self.output_folder, ".expected_outputs.json")
        os.makedirs(self.output_folder, exist_ok=True)
        with open(expected_outputs_file, 'w') as f:
            json.dump(expected_outputs_wrapped, f, indent=2)

        completion_environment = self.activate_environment(name="biopipelines")

        return f"""# Capture the main command's status if nothing already did. A tool that emits
# the missing-propagation block has it set there, right after its own command;
# one that does not reaches this footer directly, so $? is still its own.
BP_MAIN_RC=${{BP_MAIN_RC:-$?}}
{completion_environment}
# Check completion and create status files
echo "Checking outputs and creating completion status..."

python "{pipe_check_completion}" "{self.output_folder}" "{self.TOOL_NAME}" "{expected_outputs_file}"

if [ $? -eq 0 ]; then
    echo "{self.TOOL_NAME} completed successfully"
else
    echo "{self.TOOL_NAME} - some outputs missing"
fi

# The step's own status is what the failure guard reads. Trailing blocks run
# whatever happened -- a partly-successful tool still owes downstream a missing
# manifest -- but they must not overwrite the verdict.
if [ "${{BP_MAIN_RC:-0}}" -ne 0 ]; then
    echo "ERROR: {self.TOOL_NAME} exited ${{BP_MAIN_RC}}"
    exit "${{BP_MAIN_RC}}"
fi
"""
    
    def __str__(self) -> str:
        """String representation of configuration."""
        return f"{self.TOOL_NAME}(envs={self.environments}, order={self.execution_order})"
    
    def __repr__(self) -> str:
        return self.__str__()
    
    @staticmethod
    def _missing_path_of(input_source) -> Optional[str]:
        """Extract one input source's `missing` table path, or None."""
        if input_source is None or not hasattr(input_source, 'tables'):
            return None
        tables = input_source.tables
        if hasattr(tables, '_tables') and 'missing' in tables._tables:
            missing_info = tables._tables['missing']
            return missing_info.info.path if hasattr(missing_info, 'info') else str(missing_info)
        if isinstance(tables, dict) and 'missing' in tables:
            missing_info = tables['missing']
            if isinstance(missing_info, dict) and 'path' in missing_info:
                return missing_info['path']
            if hasattr(missing_info, 'info'):
                return missing_info.info.path
            return str(missing_info)
        return None

    def _collect_upstream_missing_paths(self, *input_sources) -> List[str]:
        """All distinct upstream `missing` table paths across input sources.

        Unlike ``_get_upstream_missing_table_path`` (first match only), this
        returns every input axis's manifest so a filter on more than one axis
        (e.g. both proteins and ligands) propagates fully. Order-preserving,
        de-duplicated.
        """
        paths: List[str] = []
        seen = set()
        for src in input_sources:
            p = self._missing_path_of(src)
            if p and p not in seen:
                seen.add(p)
                paths.append(p)
        return paths

    def _get_upstream_missing_table_path(self, *input_sources) -> Optional[str]:
        """
        Get the path to the missing table from upstream tool outputs.

        This should be called during get_output_files() to check if there's an upstream
        missing table that needs to be propagated. Returns the FIRST input source's
        missing table; use ``_collect_upstream_missing_paths`` to merge several axes.

        Args:
            *input_sources: Variable number of StandardizedOutput or ToolOutput objects

        Returns:
            Path to upstream missing.csv, or None if no missing table exists
        """
        paths = self._collect_upstream_missing_paths(*input_sources)
        return paths[0] if paths else None

    def upstream_missing_flag(self, *input_sources, flag: str = "--upstream-missing") -> str:
        """A CLI flag carrying EVERY upstream `missing` path, for pipe scripts.

        For tools that propagate inside their own pipe script (read upstream
        rows, merge with local failures, write their `missing.csv`): emits
        ``" --upstream-missing \\"p1\\" \\"p2\\""`` across all id-bearing input
        axes, or ``""`` when none exist. The pipe script reads them with
        ``biopipelines_io.read_upstream_missing`` (accepts the list). Passing
        every axis fixes the multi-input case where only the first manifest
        would otherwise propagate.
        """
        paths = self._collect_upstream_missing_paths(*input_sources)
        if not paths:
            return ""
        return " " + flag + "".join(f' "{p}"' for p in paths)

    def missing_table_info(self, missing_csv: Optional[str] = None) -> "TableInfo":
        """The standard `missing` TableInfo (schema id|removed_by|kind|cause).

        Tools that consume ids declare this in ``get_output_files()`` whenever an
        upstream manifest exists, so ``pipe_check_completion`` excuses the dropped
        ids instead of flagging them FAILED. Defaults the path to
        ``self.table_path("missing")``.
        """
        return TableInfo(
            name="missing",
            path=missing_csv or self.table_path("missing"),
            columns=["id", "removed_by", "kind", "cause"],
            description="IDs removed by upstream tools (or local failures) with removal reason",
        )

    def generate_missing_propagation(self, *input_sources, missing_csv: Optional[str] = None,
                                     local_missing: Optional[str] = None) -> str:
        """Bash that propagates the union of upstream `missing` manifests.

        The single inheritable propagation block: collects every input axis's
        missing-table folder and hands them to ``pipe_propagate_missing.py``,
        which merges them (de-duping by id) into this tool's own
        ``tables/missing.csv``. The tool OWNS the resulting file rather than
        pointing downstream at an upstream path, so it survives the upstream
        folder being cleaned or moved. Returns "" when no upstream manifest
        exists AND no ``local_missing`` is given (the tool then declares no
        `missing` table).

        ``local_missing`` is a path to rows this tool produced itself (e.g. a Scripting step's ``outputs.drop()`` ids); it is merged in alongside the upstream manifests, so a tool can excuse ids it filtered internally even with no upstream manifest.

        Rows are propagated in the UPSTREAM (input-axis) id space; mapping them
        into this tool's output id space (for products / ``_N`` multipliers /
        group keys) is done once, centrally, by ``pipe_check_completion`` when it
        excuses missing files — so every tool propagates raw rows uniformly.
        """
        paths = self._collect_upstream_missing_paths(*input_sources)
        if not paths and not local_missing:
            return ""
        out_csv = missing_csv or self.table_path("missing")
        propagate_py = self.pipe_script_path("pipe_propagate_missing.py")
        folders = " ".join(f'"{os.path.dirname(p)}"' for p in paths) or '""'
        local_arg = f' --local-missing "{local_missing}"' if local_missing else ""
        return f"""
# Capture the main command's status BEFORE anything else runs: this block and
# the completion check that follows both succeed on their own, and without this
# the script would exit 0 and the failure guard would see a clean step.
BP_MAIN_RC=$?

# Propagate `missing` manifest(s) from upstream tools — writes tables/missing.csv.
echo "Propagating upstream missing manifest(s)"
python "{propagate_py}" \\
    --upstream-folders {folders} \\
    --output-folder "{self.output_folder}" \\
    --missing-csv "{out_csv}"{local_arg}

"""


class _Installer(BaseConfig):
    """Generic installer step added to the pipeline by BaseConfig.install()."""

    TOOL_NAME = "install"  # Overridden dynamically in __init__

    # Sentinel file each install script must touch on success. Surfaces as
    # a tracked output so the missing-check reports COMPLETED vs FAILED.
    install_success_file = Path(lambda self: self.execution_path("install.success"))

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
        # Each _install_script must touch "$INSTALL_SUCCESS" only after a
        # successful verification (e.g. an import). The missing-check then
        # reports COMPLETED only when the file exists.
        script += f'export INSTALL_SUCCESS="{self.install_success_file}"\n'
        script += f'rm -f "$INSTALL_SUCCESS"\n'
        # Run the install body in a subshell so any ``exit`` inside it
        # (e.g. the "already installed" early return) only exits the
        # subshell. The footer below always runs and writes COMPLETED or
        # FAILED based on whether $INSTALL_SUCCESS was touched.
        script += "(\n"
        script += install_commands + "\n"
        script += ")\n"
        script += self.generate_completion_check_footer()
        return script

    def get_output_files(self):
        return {
            "tables": {
                "install_success": TableInfo(
                    name="install_success",
                    path=self.install_success_file,
                    columns=[],
                    description="Sentinel file touched by the install script on successful verification.",
                )
            },
            "output_folder": self.output_folder
        }


# ── dataflow edge recovery ───────────────────────────────────────────

# Never a dataflow input, and walking them would reach the pipeline or this tool's own outputs.
_EDGE_SKIP_ATTRS = frozenset({
    "pipeline", "pipeline_ref", "folders", "output_files", "input_sources",
    "tables", "streams", "extra_args", "resources", "environments",
    "dependencies", "filter_metadata",
})

_EDGE_WALK_MAX_DEPTH = 5
_EDGE_WALK_MAX_ITEMS = 5000


def producer_of(obj) -> Optional['BaseConfig']:
    """The tool config that produced ``obj``, or None when it carries no back-reference."""
    from .datastream import DataStream

    if isinstance(obj, StandardizedOutput):
        return obj._producer
    if isinstance(obj, ToolOutput):
        return obj.config
    if isinstance(obj, DataStream):
        return getattr(obj, "_producer", None)
    return None


def _walk_producer_refs(value, found: List[Any], depth: int = 0, seen=None) -> None:
    """Collect producer-carrying objects reachable from ``value``.

    Recurses through the shapes tools actually store an input in: Bundle/Each wrappers, lists of inputs, and dicts. Stops at anything that carries its own back-reference, and never descends into a tool config or a Pipeline.
    """
    from .datastream import DataStream
    from .combinatorics import Bundle, Each

    if value is None or depth > _EDGE_WALK_MAX_DEPTH:
        return
    if isinstance(value, (str, bytes, bool, int, float)):
        return
    if seen is None:
        seen = set()
    marker = id(value)
    if marker in seen:
        return
    seen.add(marker)

    if isinstance(value, (StandardizedOutput, ToolOutput, DataStream)):
        found.append(value)
        return
    if isinstance(value, (Bundle, Each)):
        for source in value.sources:
            _walk_producer_refs(source, found, depth + 1, seen)
        return
    if isinstance(value, (list, tuple, set, frozenset)):
        if len(value) > _EDGE_WALK_MAX_ITEMS:
            return
        for item in value:
            _walk_producer_refs(item, found, depth + 1, seen)
        return
    if isinstance(value, dict):
        if len(value) > _EDGE_WALK_MAX_ITEMS:
            return
        for item in value.values():
            _walk_producer_refs(item, found, depth + 1, seen)


def recover_input_edges(tool_config) -> List[Dict[str, Any]]:
    """Recover ``tool_config``'s dataflow inputs as a list of edge records.

    Each record is ``{"producer", "stream", "shape", "arguments"}`` — the producing tool config, the stream name when the stored input was a bare DataStream (None when a whole StandardizedOutput was stored), the class name of what was stored, and the attribute/kwarg names it was found under.
    """
    from .datastream import DataStream

    channels: List[tuple] = []
    params = getattr(tool_config, "params", None)
    if isinstance(params, dict):
        for key, value in params.items():
            if key not in _EDGE_SKIP_ATTRS:
                channels.append((key, value))
    for key, value in list(vars(tool_config).items()):
        if key in _EDGE_SKIP_ATTRS or key == "params":
            continue
        channels.append((key.lstrip("_"), value))

    edges: Dict[tuple, Dict[str, Any]] = {}
    for label, value in channels:
        found: List[Any] = []
        try:
            _walk_producer_refs(value, found)
        except Exception:
            continue
        for obj in found:
            producer = producer_of(obj)
            if producer is None or producer is tool_config:
                continue
            stream = obj.name if isinstance(obj, DataStream) else None
            key = (id(producer), stream)
            entry = edges.get(key)
            if entry is None:
                edges[key] = {
                    "producer": producer,
                    "stream": stream,
                    "shape": type(obj).__name__,
                    "arguments": [label],
                }
            elif label not in entry["arguments"]:
                entry["arguments"].append(label)

    # Storing both the whole output and a stream out of it yields two edges to one producer.
    named_producers = {k[0] for k, e in edges.items() if e["stream"]}
    return [e for k, e in edges.items() if not (k[1] is None and k[0] in named_producers)]
