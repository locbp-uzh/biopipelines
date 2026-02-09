# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Command-line entry points for BioPipelines.

Provides biopipelines-submit, biopipelines-run, and biopipelines-config
commands that work from any directory by locating the repository root automatically.
"""

import os
import re
import sys
import getpass
import subprocess


def _get_repo_root():
    """Get the biopipelines repository root (parent of PipelineScripts/)."""
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _run_script(script_name):
    """Run a bash script from the repo root, passing through all CLI arguments."""
    repo_root = _get_repo_root()
    script_path = os.path.join(repo_root, script_name)

    if not os.path.exists(script_path):
        print(f"ERROR: {script_name} not found at {script_path}")
        sys.exit(1)

    # The first non-flag argument is the pipeline script path.
    # Resolve it to an absolute path so the bash script finds it
    # regardless of its cwd.
    args = sys.argv[1:]
    for i, arg in enumerate(args):
        if not arg.startswith("-"):
            args[i] = os.path.abspath(arg)
            break

    result = subprocess.run(
        ["bash", script_path] + args,
        cwd=repo_root
    )
    sys.exit(result.returncode)


def submit():
    """Entry point for biopipelines-submit."""
    _run_script("submit")


def run():
    """Entry point for biopipelines-run."""
    _run_script("run")


def _resolve_folders():
    """Resolve all folder paths from config.yaml without creating directories."""
    from .config_manager import ConfigManager

    config_manager = ConfigManager()
    folder_config = config_manager.get_folder_config()

    # Runtime placeholders
    repo_root = _get_repo_root()
    folders = {
        "username": getpass.getuser(),
        "cwd": repo_root,
    }

    # Resolve sections in order
    placeholder_re = re.compile(r'<([a-zA-Z_][a-zA-Z0-9_]*)>')
    section_order = ['base', 'infrastructure', 'repositories', 'cache', 'derived', 'server']

    for section_name in section_order:
        section = folder_config.get(section_name, {})
        for key, path_template in section.items():
            result = path_template
            for _ in range(10):
                matches = placeholder_re.findall(result)
                if not matches:
                    break
                for ph in matches:
                    if ph in folders:
                        result = result.replace(f"<{ph}>", folders[ph])
                    else:
                        break
                else:
                    continue
                break
            folders[key] = result

    return folders


def config():
    """Entry point for biopipelines-config."""
    from .config_manager import ConfigManager

    args = sys.argv[1:]

    if not args or args[0] in ("-h", "--help"):
        print("Usage: biopipelines-config <command> [args]")
        print("")
        print("Commands:")
        print("  show                  Show full resolved configuration")
        print("  folder <key>          Get resolved folder path")
        print("  folder                List all folder keys")
        print("  env <tool>            Get conda environment for a tool")
        print("  env                   List all tool environments")
        print("  cluster <key>         Get cluster setting (env_manager, container_executor, ...)")
        print("  cluster               List all cluster settings")
        print("")
        print("Examples:")
        print("  biopipelines-config folder containers")
        print("  biopipelines-config env Boltz2")
        print("  biopipelines-config cluster env_manager")
        print("  biopipelines-config show")
        return

    command = args[0]
    rest = args[1:]

    if command == "show":
        _cmd_show()
    elif command == "folder":
        _cmd_folder(rest)
    elif command == "env":
        _cmd_env(rest)
    elif command == "cluster":
        _cmd_cluster(rest)
    else:
        print(f"Unknown command: {command}")
        print("Run 'biopipelines-config --help' for usage.")
        sys.exit(1)


def _cmd_show():
    """Print full resolved configuration."""
    from .config_manager import ConfigManager

    config_manager = ConfigManager()

    print("=== Folders ===")
    folders = _resolve_folders()
    for key, path in folders.items():
        print(f"  {key}: {path}")

    print("")
    print("=== Environments ===")
    environments = config_manager.get_all_environments()
    for tool, env in sorted(environments.items()):
        print(f"  {tool}: {env}")

    print("")
    print("=== Cluster ===")
    print(f"  env_manager: {config_manager.get_env_manager()}")
    print(f"  container_executor: {config_manager.get_container_executor()}")
    print(f"  slurm_modules: {config_manager.get_slurm_modules()}")
    print(f"  shell_hook: {config_manager.get_shell_hook_command()}")
    print(f"  module_load: {config_manager.get_module_load_line()}")

    emails = config_manager.get_emails()
    if emails:
        print(f"  emails: {len(emails)} configured")


def _cmd_folder(args):
    """Handle 'folder' subcommand."""
    folders = _resolve_folders()

    if not args:
        # List all folder keys
        for key, path in folders.items():
            print(f"{key}: {path}")
        return

    key = args[0]
    if key in folders:
        print(folders[key])
    else:
        print(f"Unknown folder key: {key}", file=sys.stderr)
        print(f"Available keys: {', '.join(folders.keys())}", file=sys.stderr)
        sys.exit(1)


def _cmd_env(args):
    """Handle 'env' subcommand."""
    from .config_manager import ConfigManager

    config_manager = ConfigManager()

    if not args:
        # List all environments
        environments = config_manager.get_all_environments()
        for tool, env in sorted(environments.items()):
            print(f"{tool}: {env}")
        return

    tool = args[0]
    env = config_manager.get_environment(tool)
    if env is not None:
        print(env)
    else:
        print(f"No environment configured for: {tool}", file=sys.stderr)
        sys.exit(1)


def _cmd_cluster(args):
    """Handle 'cluster' subcommand."""
    from .config_manager import ConfigManager

    config_manager = ConfigManager()

    cluster_values = {
        "env_manager": config_manager.get_env_manager(),
        "container_executor": config_manager.get_container_executor(),
        "slurm_modules": " ".join(config_manager.get_slurm_modules()),
        "shell_hook": config_manager.get_shell_hook_command(),
        "module_load": config_manager.get_module_load_line(),
        "activate": config_manager.get_activate_command("ENV_NAME"),
    }

    if not args:
        # List all cluster settings
        for key, value in cluster_values.items():
            print(f"{key}: {value}")
        return

    key = args[0]
    if key == "emails":
        emails = config_manager.get_emails()
        for user, email in emails.items():
            print(f"{user}: {email}")
    elif key in cluster_values:
        print(cluster_values[key])
    else:
        print(f"Unknown cluster key: {key}", file=sys.stderr)
        print(f"Available keys: {', '.join(list(cluster_values.keys()) + ['emails'])}", file=sys.stderr)
        sys.exit(1)
