# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Command-line entry points for BioPipelines.

Provides biopipelines-submit, biopipelines-run, biopipelines-config, and biopipelines-otf
commands that work from any directory by locating the repository root automatically.
"""

import json
import os
import re
import sys
import getpass
import subprocess
import tempfile


def _get_repo_root():
    """Get the biopipelines repository root (parent of biopipelines/)."""
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _notebook_to_script(notebook_path):
    """Convert a Jupyter notebook to a temporary Python script.

    Extracts code cells, skipping IPython magics and shell commands,
    and writes them to a temporary .py file in the same directory
    so that relative paths in the pipeline still work.
    """
    with open(notebook_path, 'r', encoding='utf-8') as f:
        nb = json.load(f)

    code_lines = []
    for cell in nb.get('cells', []):
        if cell.get('cell_type') != 'code':
            continue
        source = cell.get('source', [])
        # Join source lines (they may or may not end with \n)
        cell_code = ''.join(source).strip()
        if not cell_code:
            continue
        # Skip IPython magics and shell commands
        first_line = cell_code.split('\n')[0].strip()
        if first_line.startswith(('%', '!')):
            continue
        code_lines.append(cell_code)

    script_content = '\n\n'.join(code_lines) + '\n'

    # Write to temp file in same directory (preserves relative paths)
    notebook_dir = os.path.dirname(notebook_path)
    stem = os.path.splitext(os.path.basename(notebook_path))[0]
    tmp_path = os.path.join(notebook_dir, f'._biopipelines_tmp_{stem}.py')

    with open(tmp_path, 'w', encoding='utf-8') as f:
        f.write(script_content)

    return tmp_path


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
    tmp_script = None
    for i, arg in enumerate(args):
        if not arg.startswith("-"):
            args[i] = os.path.abspath(arg)
            # Convert notebook to temporary Python script if needed
            if args[i].endswith('.ipynb'):
                tmp_script = _notebook_to_script(args[i])
                args[i] = tmp_script
            break

    try:
        result = subprocess.run(
            ["bash", script_path] + args,
            cwd=repo_root
        )
    finally:
        if tmp_script and os.path.exists(tmp_script):
            os.remove(tmp_script)

    sys.exit(result.returncode)


def submit():
    """Entry point for biopipelines-submit."""
    _run_script("submit")


def run():
    """Entry point for biopipelines-run."""
    _run_script("run")


def _generate_gpu_line(gpu_spec):
    """Generate SBATCH GPU line from a GPU specification string."""
    if gpu_spec is None or gpu_spec in ("none", ""):
        return ""
    elif gpu_spec == "high-memory":
        return '#SBATCH --gpus=1\n#SBATCH --constraint="GPUMEM32GB|GPUMEM80GB|GPUMEM96GB"'
    elif gpu_spec in ("gpu", "any"):
        return "#SBATCH --gpus=1"
    elif gpu_spec.startswith("!"):
        excluded = gpu_spec[1:]
        if excluded.upper() == "L4":
            return '#SBATCH --gpus=1\n#SBATCH --constraint="GPUMEM32GB|GPUMEM80GB|GPUMEM96GB"'
        return f'#SBATCH --gpus=1\n#SBATCH --constraint="~GPU{excluded}"'
    elif gpu_spec in ("24GB", "32GB", "80GB", "96GB") or "|" in gpu_spec:
        if "|" in gpu_spec:
            parts = [f"GPUMEM{m}" for m in gpu_spec.split("|")]
            constraint = "|".join(parts)
        else:
            constraint = f"GPUMEM{gpu_spec}"
        return f'#SBATCH --gpus=1\n#SBATCH --constraint="{constraint}"'
    else:
        return f"#SBATCH --gpus={gpu_spec}:1"


def otf():
    """Entry point for bp-otf: submit a pipeline for on-the-fly execution on SLURM.

    Uses sbatch so the job survives terminal disconnects, then tails the output
    file to stream results back to the terminal.  Ctrl+C stops tailing but the
    job keeps running.  Use --detach to skip tailing.
    """
    from .config_manager import ConfigManager

    args = sys.argv[1:]

    # Parse arguments
    pipeline_path = None
    mem = "16GB"
    time_limit = "24:00:00"
    gpu = None
    output_path = None
    detach = False
    extra_sbatch = []

    i = 0
    while i < len(args):
        arg = args[i]
        if arg in ("-h", "--help"):
            print("Usage: bp-otf [options] <pipeline.py|notebook.ipynb>")
            print("")
            print("Submit a pipeline for on-the-fly execution on a SLURM compute node.")
            print("The pipeline runs as a live Python process (on_the_fly mode).")
            print("Output is streamed to the terminal. Ctrl+C stops tailing but the")
            print("job keeps running on the compute node.")
            print("")
            print("Options:")
            print("  --mem=<size>      Memory allocation (default: 16GB)")
            print("  --time=<time>     Wall time limit (default: 24:00:00)")
            print("  --gpu=<spec>      GPU specification (e.g., any, A100, 80GB)")
            print("  --output=<path>   SLURM output file (default: <pipeline_dir>/<name>_otf.out)")
            print("  --detach          Submit and exit without tailing output")
            print("  --<key>=<value>   Any additional SBATCH parameter")
            print("")
            print("Examples:")
            print("  bp-otf pipeline.py")
            print("  bp-otf pipeline.py --mem=32GB --time=4:00:00 --gpu=any")
            print("  bp-otf --detach notebook.ipynb --gpu=A100 --partition=gpu")
            return
        elif arg == "--detach":
            detach = True
        elif arg.startswith("--mem="):
            mem = arg.split("=", 1)[1]
        elif arg.startswith("--time="):
            time_limit = arg.split("=", 1)[1]
        elif arg.startswith("--gpu="):
            gpu = arg.split("=", 1)[1]
        elif arg.startswith("--output="):
            output_path = arg.split("=", 1)[1]
        elif arg.startswith("--"):
            extra_sbatch.append(arg)
        else:
            pipeline_path = os.path.abspath(arg)
        i += 1

    if pipeline_path is None:
        print("ERROR: No pipeline script specified.")
        print("Run 'bp-otf --help' for usage.")
        sys.exit(1)

    if not os.path.exists(pipeline_path):
        print(f"ERROR: Pipeline script not found: {pipeline_path}")
        sys.exit(1)

    # Convert notebook to script if needed
    tmp_script = None
    if pipeline_path.endswith('.ipynb'):
        tmp_script = _notebook_to_script(pipeline_path)
        pipeline_path = tmp_script

    # Resolve output path
    pipeline_dir = os.path.dirname(pipeline_path)
    pipeline_stem = os.path.splitext(os.path.basename(pipeline_path))[0]
    if tmp_script:
        # Use the original notebook name for the output file
        pipeline_stem = pipeline_stem.replace("._biopipelines_tmp_", "")
    if output_path is None:
        output_path = os.path.join(pipeline_dir, f"{pipeline_stem}_otf.out")

    # Get cluster config
    config_manager = ConfigManager()
    module_load = config_manager.get_module_load_line()
    shell_hook = config_manager.get_shell_hook_command()
    activate = config_manager.get_activate_command("biopipelines")

    # Build GPU line
    gpu_line = _generate_gpu_line(gpu)

    # Build extra SBATCH lines
    extra_lines = ""
    for opt in extra_sbatch:
        # --partition=gpu -> #SBATCH --partition=gpu
        extra_lines += f"\n#SBATCH {opt}"

    # Generate sbatch script
    sbatch_content = f"""#!/usr/bin/bash
#SBATCH --mem={mem}
#SBATCH --time={time_limit}
#SBATCH --output={output_path}
#SBATCH --job-name={pipeline_stem}_otf"""

    if gpu_line:
        sbatch_content += f"\n{gpu_line}"

    if extra_lines:
        sbatch_content += extra_lines

    sbatch_content += f"""

# Make all files group-writable by default
umask 002
{module_load}
{shell_hook}
{activate}

# Force on-the-fly mode
export BIOPIPELINES_OTF=1

python "{pipeline_path}"
"""

    # Write temp sbatch script and submit
    sbatch_path = None
    try:
        fd, sbatch_path = tempfile.mkstemp(suffix=".sh", prefix="bp_otf_")
        with os.fdopen(fd, 'w') as f:
            f.write(sbatch_content)

        result = subprocess.run(
            ["sbatch", "--parsable", sbatch_path],
            capture_output=True, text=True
        )

        if result.returncode != 0:
            print(f"ERROR: sbatch failed: {result.stderr.strip()}")
            sys.exit(1)

        job_id = result.stdout.strip()
        print(f"Submitted job {job_id}")
        print(f"Output: {output_path}")

        if not detach:
            # Tail the output file so the user sees live output.
            # Poll squeue to stop automatically when the job finishes.
            # Ctrl+C stops tailing but the SLURM job keeps running.
            print(f"Tailing output (Ctrl+C to detach, job {job_id} keeps running)...")
            print("---")
            try:
                subprocess.run(["bash", "-c",
                    f'touch "{output_path}";'
                    f' tail -f "{output_path}" &'
                    f' TAIL_PID=$!;'
                    f' while squeue -j {job_id} -h 2>/dev/null | grep -q {job_id}; do sleep 5; done;'
                    f' sleep 1; kill $TAIL_PID 2>/dev/null; wait $TAIL_PID 2>/dev/null;'
                    f' echo "---"; echo "Job {job_id} finished."'
                ])
            except KeyboardInterrupt:
                print(f"\n---\nDetached. Job {job_id} is still running.")
                print(f"Re-attach with: tail -f {output_path}")
                print(f"Cancel with:   scancel {job_id}")

    finally:
        if sbatch_path and os.path.exists(sbatch_path):
            os.remove(sbatch_path)
        # Clean up notebook tmp script if created
        if tmp_script and os.path.exists(tmp_script):
            os.remove(tmp_script)


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
