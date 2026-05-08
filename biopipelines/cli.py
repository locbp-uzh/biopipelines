# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Command-line entry points for BioPipelines.

Provides biopipelines-submit and biopipelines-config commands that
work from any directory by locating the repository root automatically.
"""

import json
import os
import re
import sys
import getpass
import subprocess


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
    """Entry point for bp-run.

    Activates the biopipelines environment, exports BIOPIPELINES_OTF=1
    so Pipeline(...) treats this as on-the-fly, and runs the given
    pipeline script inline (no SLURM).
    """
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
        print("  list                  List every config.<variant>.yaml in the repo")
        print("  path [--variant V]    Print absolute path of active config YAML")
        print("  edit [--variant V]    Open the config in an interactive TUI editor")
        print("                        (no --variant: pop a variant picker first)")
        print("  auto [--variant V]    Probe the host (username, env manager, scheduler,")
        print("                        modules, lmod profile, container runtime, git email)")
        print("                        and write the discovered values into the chosen")
        print("                        config (a .bak is created first). Without --variant,")
        print("                        a picker shows every variant + a [dry-run] preview.")
        print("  folder <key>          Get resolved folder path")
        print("  folder                List all folder keys")
        print("  env <tool>            Get conda environment for a tool")
        print("  env                   List all tool environments")
        print("  machine <key>         Get machine setting (env_manager, scheduler, ...)")
        print("  machine               List all machine settings")
        print("")
        print("Examples:")
        print("  biopipelines-config path")
        print("  $EDITOR \"$(biopipelines-config path)\"")
        print("  biopipelines-config list")
        print("  biopipelines-config edit")
        print("  biopipelines-config edit --variant cluster")
        print("  biopipelines-config auto                       # picker (variants + dry-run)")
        print("  biopipelines-config auto --variant cluster     # write straight to cluster.yaml")
        print("  biopipelines-config folder containers")
        print("  biopipelines-config env Boltz2")
        print("  biopipelines-config machine env_manager")
        print("  biopipelines-config show")
        return

    command = args[0]
    rest = args[1:]

    if command == "show":
        _cmd_show()
    elif command == "list":
        _cmd_list()
    elif command == "path":
        _cmd_path(rest)
    elif command == "edit":
        _cmd_edit(rest)
    elif command == "auto":
        _cmd_auto(rest)
    elif command == "folder":
        _cmd_folder(rest)
    elif command == "env":
        _cmd_env(rest)
    elif command in ("machine", "cluster"):
        _cmd_machine(rest)
    else:
        print(f"Unknown command: {command}")
        print("Run 'biopipelines-config --help' for usage.")
        sys.exit(1)


def _cmd_show():
    """Print the active config as YAML.

    A header line names the active variant and the absolute path of the
    loaded file; the rest is the raw resolved config dumped through
    ``yaml.safe_dump`` so the output is a faithful reflection of what
    ``ConfigManager`` actually loaded — no synthesised keys, no
    hand-picked subset.
    """
    import yaml
    from .config_manager import ConfigManager

    cm = ConfigManager()
    print(f"config: {cm.get_variant()}   ({cm._get_config_path()})")
    print()
    yaml.safe_dump(
        cm._config,
        sys.stdout,
        sort_keys=False,
        default_flow_style=False,
    )


def _parse_variant_flag(args):
    """Pop ``--variant <name>`` from args (in place); return variant or None."""
    variant = None
    if "--variant" in args:
        i = args.index("--variant")
        if i + 1 >= len(args):
            print("--variant requires an argument", file=sys.stderr)
            sys.exit(2)
        variant = args[i + 1]
        del args[i:i + 2]
    return variant


def _cmd_path(args):
    """Print the absolute path of the active config YAML.

    Designed to be used in shells:
      $EDITOR "$(biopipelines-config path)"
    """
    from .config_manager import ConfigManager

    args = list(args)
    variant = _parse_variant_flag(args)
    if args:
        print(f"Unexpected argument: {args[0]}", file=sys.stderr)
        sys.exit(2)

    path = ConfigManager._get_config_path(variant=variant)
    print(path)


def _discover_variants():
    """Return the sorted list of variant names found in the repo root.

    A variant is the ``<v>`` of every ``config.<v>.yaml`` file at the
    repo root. Order: alphabetical, matching ``ls`` output.
    """
    import glob
    repo_root = _get_repo_root()
    variants = []
    for path in sorted(glob.glob(os.path.join(repo_root, "config.*.yaml"))):
        base = os.path.basename(path)
        v = base[len("config."):-len(".yaml")]
        if v:
            variants.append(v)
    return variants


def _cmd_list():
    """Print every config.<variant>.yaml at the repo root, marking the active one."""
    from .config_manager import ConfigManager

    variants = _discover_variants()
    if not variants:
        print("No config.<variant>.yaml files found in the repo root.", file=sys.stderr)
        sys.exit(1)

    try:
        active = ConfigManager().get_variant()
    except Exception:
        active = None

    repo_root = _get_repo_root()
    for v in variants:
        marker = " (active)" if v == active else ""
        print(f"  {v}{marker}\t{os.path.join(repo_root, f'config.{v}.yaml')}")


def _pick_variant_interactive(active, extra_options=None, prepend_options=None,
                              title=None, default_key=None):
    """Arrow-key picker for the available config variants.

    The cursor starts on the active variant when it appears in the list,
    otherwise on the first row. ↑/↓ moves; Enter selects; q / Ctrl+C /
    Esc aborts. Returns the chosen value (a variant name or one of the
    extra-option keys), or None on abort.

    prepend_options / extra_options: optional lists of (key, label, style)
    tuples placed before / after the variants respectively. Style is a
    prompt_toolkit class name applied to the row when the cursor isn't
    on it (e.g. "muted" for a dry-run entry rendered in italic gray).

    Skips the picker entirely when only one variant exists and no
    extra rows are given. Falls back to the active variant when stdin
    is not a TTY, or exits with a clear message when there is no
    active variant in non-interactive mode.
    """
    variants = _discover_variants()
    if not variants:
        print("No config.<variant>.yaml files found in the repo root.", file=sys.stderr)
        sys.exit(1)

    prepend_options = list(prepend_options or [])
    extra_options = list(extra_options or [])
    if len(variants) == 1 and not prepend_options and not extra_options:
        return variants[0]

    # Build the unified list of (key, label, style) for the picker.
    # Variant entries get an empty style; extras carry their own.
    items = list(prepend_options)
    items.extend((v, v + (" (active)" if v == active else ""), "") for v in variants)
    items.extend(extra_options)

    if not sys.stdin.isatty():
        if active and active in variants:
            return active
        print(
            "Cannot pick a variant non-interactively; pass --variant <name>.",
            file=sys.stderr,
        )
        sys.exit(2)

    from prompt_toolkit import Application
    from prompt_toolkit.key_binding import KeyBindings
    from prompt_toolkit.layout import Layout, Window, HSplit, FormattedTextControl
    from prompt_toolkit.layout.dimension import Dimension
    from prompt_toolkit.formatted_text import FormattedText
    from prompt_toolkit.styles import Style

    # State held in a 1-item list so the closures can mutate it.
    # Cursor starts on default_key if given, else on the active variant
    # if present, else on the first item.
    if default_key is not None:
        start = next((i for i, it in enumerate(items) if it[0] == default_key), 0)
    else:
        start = next((i for i, it in enumerate(items) if it[0] == active), 0)
    cursor = [start]
    chosen = [None]
    title_text = title or "Choose config variant"

    def _render():
        lines = [("class:title", f"{title_text} (↑/↓ to move, Enter to select, q to cancel)\n\n")]
        for i, (_key, label, item_style) in enumerate(items):
            if i == cursor[0]:
                lines.append(("class:cursor", f"  > {label}\n"))
            else:
                style_class = f"class:{item_style}" if item_style else ""
                lines.append((style_class, f"    {label}\n"))
        return FormattedText(lines)

    body = Window(
        content=FormattedTextControl(_render),
        height=Dimension(min=len(items) + 2, preferred=len(items) + 3),
        always_hide_cursor=True,
    )

    kb = KeyBindings()

    @kb.add("up")
    def _(event):
        cursor[0] = (cursor[0] - 1) % len(items)

    @kb.add("down")
    def _(event):
        cursor[0] = (cursor[0] + 1) % len(items)

    @kb.add("enter")
    def _(event):
        chosen[0] = items[cursor[0]][0]
        event.app.exit()

    @kb.add("q")
    @kb.add("escape")
    @kb.add("c-c")
    def _(event):
        event.app.exit()

    style = Style.from_dict({
        "title": "bold",
        "cursor": "reverse",
        "muted": "fg:ansigray italic",
    })

    app = Application(
        layout=Layout(body),
        key_bindings=kb,
        full_screen=False,
        style=style,
    )
    app.run()
    return chosen[0]


def _cmd_edit(args):
    """Open a config YAML in the interactive TUI editor.

    With ``--variant <name>``: edit that variant directly.
    Without: pop a variant picker first; the picker marks the currently-
    active variant (resolved from the username mapping) so the user knows
    which one is the default.
    """
    from .config_manager import ConfigManager
    from .config_editor import run_editor

    args = list(args)
    variant = _parse_variant_flag(args)
    if args:
        print(f"Unexpected argument: {args[0]}", file=sys.stderr)
        sys.exit(2)

    if variant is None:
        try:
            active = ConfigManager().get_variant()
        except Exception:
            active = None
        variant = _pick_variant_interactive(active)
        if variant is None:
            sys.exit(1)

    path = ConfigManager._get_config_path(variant=variant)
    rc = run_editor(path)
    # Tell the user which file they were editing on the way out, so the
    # variant choice is visible above any of the editor's own messages.
    print(f"config has been set to {variant}: {path}")
    sys.exit(rc)


# ── bp-config auto ────────────────────────────────────────────────────────────
#
# The probes below favour experiment over hardcoded names. The basic
# loop is: source candidate shell-init files in a subshell, then `module
# load <candidate>` candidate module names, then check whether the
# binary we care about appeared on PATH. The first combination that
# produces a working binary wins.
#
# Per-call state is bundled into a small ProbeContext namespace so the
# helpers can share the verbose flag and a warnings list without global
# variables. The cli.config() dispatcher creates a fresh context per
# command invocation.


class _ProbeContext:
    """Per-invocation state shared across probe helpers."""

    def __init__(self, verbose: bool = False) -> None:
        self.verbose = verbose
        self.warnings: list = []
        # Cached results so a second probe doesn't re-run the same shell.
        self._scheduler_init_files = None  # type: Optional[list]
        self._module_avail = None  # type: Optional[str]

    def log(self, msg: str) -> None:
        if self.verbose:
            print(f"  [probe] {msg}", file=sys.stderr)

    def warn(self, msg: str) -> None:
        self.warnings.append(msg)


def _which(*names):
    """Return the first executable from `names` found on PATH, or None."""
    import shutil
    for n in names:
        p = shutil.which(n)
        if p:
            return p, n
    return None, None


def _probe_username():
    """The active Unix username."""
    import getpass
    try:
        return getpass.getuser()
    except Exception:
        return os.environ.get("USER") or os.environ.get("USERNAME") or ""


def _probe_email():
    """git config user.email if set; empty otherwise."""
    import subprocess
    try:
        r = subprocess.run(
            ["git", "config", "--get", "user.email"],
            capture_output=True, text=True, timeout=2,
        )
        if r.returncode == 0:
            return r.stdout.strip()
    except Exception:
        pass
    return ""


def _probe_emails():
    """Build an emails dict from `{<username>: <git user.email>}`. Empty if either is unknown."""
    user = _probe_username()
    email = _probe_email()
    if user and email:
        return {user: email}
    return {}


# ── shell-init discovery (lmod / site-specific profile.d files) ──────────────

def _probe_scheduler_init_files(ctx: _ProbeContext) -> list:
    """Find which /etc/profile.d/*.sh files actually contribute lmod state.

    Algorithm: take a baseline subshell with a clean environment, record
    its `MODULEPATH` and whether `module` is a function. For each
    /etc/profile.d/*.sh, source it on top of the baseline and check
    whether either changed. Keep the files that contributed.

    This generalises beyond the S3IT-specific `lmod.sh + z*_lmodenv*.sh`
    pattern to any cluster that puts module setup in profile.d under
    arbitrary names (`spack.sh`, `cluster-init.sh`, `cc_modules.sh`,
    etc.). Returns the list of contributing absolute paths in
    alphabetical order (matching the order /etc/profile would source
    them in a normal login).
    """
    if ctx._scheduler_init_files is not None:
        return ctx._scheduler_init_files
    import glob
    import subprocess

    profile_d = "/etc/profile.d"
    if not os.path.isdir(profile_d):
        ctx.log(f"{profile_d} not present; no scheduler.init")
        ctx._scheduler_init_files = []
        return []

    all_files = sorted(glob.glob(os.path.join(profile_d, "*.sh")))
    if not all_files:
        ctx.log(f"no *.sh under {profile_d}; no scheduler.init")
        ctx._scheduler_init_files = []
        return []

    def _probe_state(extra_lines: list) -> tuple:
        """Run a subshell that sources `extra_lines`, return
        (module_is_function: bool, modulepath: str)."""
        prelude = "; ".join(extra_lines)
        if prelude:
            prelude += "; "
        cmd = (
            f"{prelude}"
            "if type module >/dev/null 2>&1 && [ \"$(type -t module)\" = function ]; "
            "then echo MODULE_FUNC=1; else echo MODULE_FUNC=0; fi; "
            "echo MODULEPATH=${MODULEPATH:-}"
        )
        try:
            r = subprocess.run(
                ["bash", "--noprofile", "--norc", "-c", cmd],
                capture_output=True, text=True, timeout=10, env={"PATH": os.environ.get("PATH", "")},
            )
            txt = r.stdout
        except Exception as e:
            ctx.log(f"subshell failed for state probe: {e}")
            return False, ""
        is_func = "MODULE_FUNC=1" in txt
        mp = ""
        for line in txt.splitlines():
            if line.startswith("MODULEPATH="):
                mp = line[len("MODULEPATH="):]
                break
        return is_func, mp

    base_func, base_mp = _probe_state([])
    ctx.log(f"baseline: module-is-function={base_func}, MODULEPATH={base_mp!r}")

    contributing = []
    for f in all_files:
        # Source baseline (nothing) + this single file.
        line = f'[ -f "{f}" ] && source "{f}"'
        is_func, mp = _probe_state([line])
        if (is_func and not base_func) or (mp != base_mp):
            ctx.log(f"  {f}: contributes (func={is_func}, MODULEPATH delta={mp != base_mp})")
            contributing.append(f)
        else:
            ctx.log(f"  {f}: no effect; skipped")

    ctx._scheduler_init_files = contributing
    return contributing


def _scheduler_init_lines(ctx: _ProbeContext) -> list:
    """Format the discovered profile.d files as one bash for-loop line.

    Empty list when nothing contributes (cluster doesn't use lmod /
    profile.d, or the host has no module system). The `[ -f "$f" ]`
    guard makes the line robust against later sysadmin churn.
    """
    files = _probe_scheduler_init_files(ctx)
    if not files:
        return []
    return [
        'for f in ' + ' '.join(files) + '; do [ -f "$f" ] && source "$f"; done'
    ]


# ── module discovery ─────────────────────────────────────────────────────────

def _module_avail(ctx: _ProbeContext) -> str:
    """Cached `module --terse avail` output after sourcing scheduler.init."""
    if ctx._module_avail is not None:
        return ctx._module_avail
    import subprocess
    init = _scheduler_init_lines(ctx)
    init_block = "; ".join(init) + "; " if init else ""
    try:
        r = subprocess.run(
            ["bash", "-c", f"{init_block}module --terse avail 2>&1 || true"],
            capture_output=True, text=True, timeout=15,
        )
        ctx._module_avail = r.stdout + r.stderr
    except Exception as e:
        ctx.log(f"module avail failed: {e}")
        ctx._module_avail = ""
    return ctx._module_avail


def _module_candidates(ctx: _ProbeContext, regex: str) -> list:
    """All module names whose bare-alias matches `regex` (case-insensitive).

    Returns a list of (preferred_name, display_label) where preferred_name
    is the bare alias if present, else the highest-versioned variant.
    Used by the env-manager and container-runtime probes to enumerate
    cluster modules worth trying.
    """
    import re
    avail = _module_avail(ctx)
    pat = re.compile(regex, re.IGNORECASE)
    bare: set = set()
    versioned: dict = {}
    for ln in avail.splitlines():
        ln = ln.strip()
        if not ln or ln.endswith(":"):
            continue
        if ln.endswith("/"):
            base = ln.rstrip("/")
            bare.add(base)
            continue
        if "/" in ln:
            base, _, _ver = ln.partition("/")
            versioned.setdefault(base, []).append(ln)
        else:
            bare.add(ln)

    out = []
    seen = set()
    for name in sorted(bare | set(versioned.keys())):
        if not pat.search(name) or name in seen:
            continue
        seen.add(name)
        if name in bare:
            out.append(name)
        else:
            out.append(sorted(versioned[name])[-1])
    return out


def _probe_module_for_binary(
    ctx: _ProbeContext,
    binary_names: list,
    module_regex: str,
) -> tuple:
    """Find the module name + the binary it exposes by experiment.

    For each cluster module matching `module_regex`, run
    `module load <name>` in a subshell and check whether any of
    `binary_names` appears on PATH. Returns (module_name, binary_name,
    binary_path) for the first hit, or (None, None, None) if nothing
    works.

    Honours the cached scheduler.init so module commands actually
    function in the subshell.
    """
    import subprocess
    init = _scheduler_init_lines(ctx)
    init_block = "; ".join(init) + "; " if init else ""
    candidates = _module_candidates(ctx, module_regex)
    ctx.log(f"  module candidates for /{module_regex}/: {candidates}")
    for mod in candidates:
        for bin_name in binary_names:
            cmd = (
                f"{init_block}module load {mod} 2>/dev/null; "
                f"command -v {bin_name} 2>/dev/null || true"
            )
            try:
                r = subprocess.run(
                    ["bash", "-c", cmd],
                    capture_output=True, text=True, timeout=10,
                )
                p = r.stdout.strip()
            except Exception as e:
                ctx.log(f"  module load {mod} → subshell failed ({e})")
                continue
            if p:
                ctx.log(f"  module {mod} → {bin_name} at {p}")
                return mod, bin_name, p
    return None, None, None


def _probe_binary_on_clean_path(ctx: _ProbeContext, binary_names: list) -> tuple:
    """Find a binary on the user's existing PATH (no module load).

    Useful for laptop installs where the env manager is on PATH from
    the user's `.bashrc` and not behind a module system.
    """
    for bin_name in binary_names:
        p, _ = _which(bin_name)
        if p:
            ctx.log(f"  {bin_name} on PATH at {p}")
            return bin_name, p
    return None, None


# ── env manager + container executor probes ────────────────────────────────

def _shell_hook_line(mgr_name: str) -> str:
    """The eval-able hook line for a given env-manager binary name."""
    if mgr_name == "conda":
        return 'eval "$(conda shell.bash hook)"'
    return f'eval "$({mgr_name} shell hook --shell bash)"'


def _probe_env_manager(ctx: _ProbeContext) -> tuple:
    """Discover the env manager. Returns ({name, init}, source_label).

    Order of attempts:
      1. mamba / micromamba / conda already on the user's PATH.
      2. mamba / micromamba / conda via a cluster module.
      3. pixi (PATH or module) — known but not supported as an
         env manager; we warn and leave the user to decide.
      4. pip on PATH — degraded mode, no env activation.

    Emits warnings when the result is degraded (pixi-only, pip-only,
    nothing).
    """
    SUPPORTED = ("mamba", "micromamba", "conda")
    UNSUPPORTED = ("pixi",)

    # 1. PATH probe.
    bin_name, bin_path = _probe_binary_on_clean_path(ctx, list(SUPPORTED))
    if bin_name:
        return {"name": bin_name, "init": [_shell_hook_line(bin_name)]}, f"PATH ({bin_path})"

    # 2. Module probe.
    mod, bin_name, bin_path = _probe_module_for_binary(
        ctx, list(SUPPORTED),
        r"^(miniforge|mambaforge|miniconda|anaconda|mamba|conda|micromamba|python(/|$))",
    )
    if bin_name:
        return {"name": bin_name, "init": [_shell_hook_line(bin_name)]}, f"module {mod} ({bin_path})"

    # 3. Pixi check — recognised but unsupported.
    px_name, px_path = _probe_binary_on_clean_path(ctx, list(UNSUPPORTED))
    if not px_name:
        _, px_name, px_path = _probe_module_for_binary(
            ctx, list(UNSUPPORTED), r"^pixi",
        )
    if px_name:
        ctx.warn(
            f"found pixi at {px_path} but pixi is not currently supported "
            "as an env_manager; leaving env_manager unchanged. "
            "(install mamba / micromamba / conda or hand-edit env_manager.)"
        )
        return None, f"pixi only ({px_path}; unsupported)"

    # 4. pip fallback.
    pip_name, pip_path = _probe_binary_on_clean_path(ctx, ["pip"])
    if pip_name:
        ctx.warn(
            f"only pip found (at {pip_path}) — no conda-style env manager available. "
            "BioPipelines pipelines will skip env activation; install mamba / "
            "micromamba / conda for the standard cluster setup."
        )
        return {"name": "pip", "init": []}, f"pip fallback ({pip_path})"

    # Nothing found.
    ctx.warn(
        "no env manager found (mamba / micromamba / conda / pip not on PATH "
        "or available as a module). Leaving env_manager unchanged."
    )
    return None, "not found"


def _probe_container_executor(ctx: _ProbeContext) -> tuple:
    """Discover the container runtime. Returns (name, source_label).

    Same two-stage approach as the env manager: PATH first, then
    cluster modules. If nothing is found, returns ('none', '').
    """
    NAMES = ("apptainer", "singularity")

    bin_name, bin_path = _probe_binary_on_clean_path(ctx, list(NAMES))
    if bin_name:
        return bin_name, f"PATH ({bin_path})"

    mod, bin_name, bin_path = _probe_module_for_binary(
        ctx, list(NAMES),
        r"^(apptainer|singularity|container)",
    )
    if bin_name:
        return bin_name, f"module {mod} ({bin_path})"

    return "none", "not found"


# ── scheduler probe ──────────────────────────────────────────────────────────

def _probe_scheduler(ctx: _ProbeContext) -> tuple:
    """Return ({name, init, modules}, source_label) for machine.scheduler.

    `modules` is the list of cluster modules to `module load` in
    SLURM batch scripts. We populate it from whatever modules the
    env-manager and container probes resolved to (since those are the
    binaries each batch script needs available).
    """
    sinfo, _ = _which("sinfo")
    if not sinfo:
        return {"name": "none", "init": [], "modules": []}, "no sinfo on PATH"

    init = _scheduler_init_lines(ctx)

    # Re-run the env-mgr / container module probes specifically to pick
    # up the module *names* (already cached in ctx so this is cheap).
    em_block, _em_src = _probe_env_manager(ctx)
    ce_name, _ce_src = _probe_container_executor(ctx)

    modules: list = []
    # Re-discover which module exposed mamba / conda / etc.
    if em_block and em_block.get("name") in ("mamba", "micromamba", "conda"):
        m, _, _ = _probe_module_for_binary(
            ctx, [em_block["name"]],
            r"^(miniforge|mambaforge|miniconda|anaconda|mamba|conda|micromamba|python(/|$))",
        )
        if m:
            modules.append(m)
    if ce_name in ("apptainer", "singularity"):
        m, _, _ = _probe_module_for_binary(
            ctx, [ce_name], r"^(apptainer|singularity|container)",
        )
        if m:
            modules.append(m)

    return {
        "name": "slurm",
        "init": init,
        "modules": modules,
    }, "sinfo present"


def _build_machine_block(ctx: _ProbeContext) -> tuple:
    """Run all the probes and return (machine_dict, [(field, source), ...]).

    machine_dict matches the on-disk schema. The source list explains
    where each field's value came from (PATH, module name, "not
    found"), and is printed in the auto report.
    """
    em_block, em_src = _probe_env_manager(ctx)
    sched_block, sched_src = _probe_scheduler(ctx)
    ce_name, ce_src = _probe_container_executor(ctx)

    machine = {
        "username": _probe_username(),
        "emails": _probe_emails(),
        "env_manager": em_block,  # may be None on miss
        "scheduler": sched_block,
        "container_executor": ce_name,
    }
    sources = [
        ("env_manager", em_src),
        ("scheduler", sched_src),
        ("container_executor", ce_src),
    ]
    return machine, sources


def _diff_machine(old, new):
    """Yield (path, old_value, new_value) tuples for keys that differ."""
    if not isinstance(old, dict):
        old = {}
    keys = sorted(set(old.keys()) | set(new.keys()))
    for k in keys:
        ov, nv = old.get(k), new.get(k)
        if ov != nv:
            yield k, ov, nv


def _machine_for_write(new_machine, old_machine):
    """Drop env_manager from the write payload when the probe failed.

    Falls back to whatever the target variant already had for
    env_manager so the existing curated value is preserved. Anywhere
    else we treat `None` as a probe miss the user explicitly does not
    want overwritten.
    """
    if new_machine.get('env_manager') is None:
        out = {k: v for k, v in new_machine.items() if k != 'env_manager'}
        if isinstance(old_machine, dict) and 'env_manager' in old_machine:
            out['env_manager'] = old_machine['env_manager']
        return out
    return new_machine


def _read_machine_block(variant):
    """Return the existing `machine` block for a variant (dict; {} if absent)."""
    from .config_manager import ConfigManager
    from ruamel.yaml import YAML
    path = ConfigManager._get_config_path(variant=variant)
    if not os.path.isfile(path):
        return {}, path
    yaml = YAML()
    yaml.preserve_quotes = True
    try:
        with open(path, "r", encoding="utf-8") as f:
            doc = yaml.load(f)
    except Exception:
        return {}, path
    if doc is None:
        return {}, path
    return (doc.get("machine") or {}), path


def _yaml_value(value, indent=4):
    """Format a value as YAML block-style text, indented by `indent` spaces.

    Scalars render on a single line ("mamba"); mappings and sequences
    render as multi-line YAML. Long string values (e.g. the lmod
    for-loop in scheduler.init) are kept on a single line by setting
    a very wide line-wrap; the goal is to mirror what the on-disk
    config will actually look like, not auto-wrap.
    """
    from ruamel.yaml import YAML
    import io
    pad = " " * indent
    if value is None:
        return pad + "null\n"
    if isinstance(value, bool):
        return pad + ("true" if value else "false") + "\n"
    if isinstance(value, (int, float)):
        return pad + repr(value) + "\n"
    if isinstance(value, str):
        # Quote when ruamel would (special chars, leading whitespace);
        # otherwise emit bare. Easiest: round-trip through ruamel and
        # strip the document-end '...' marker that gets appended for
        # a top-level scalar dump.
        y = YAML()
        y.width = 10**6
        buf = io.StringIO()
        y.dump(value, buf)
        text = buf.getvalue()
        if text.endswith("...\n"):
            text = text[:-4]
        return pad + text
    if isinstance(value, dict) and not value:
        return pad + "{}\n"
    if isinstance(value, list) and not value:
        return pad + "[]\n"
    y = YAML()
    y.default_flow_style = False
    y.width = 10**6
    buf = io.StringIO()
    y.dump(value, buf)
    out = buf.getvalue()
    return "".join(pad + line for line in out.splitlines(keepends=True))


def _pick_apply_target(new_machine, active, warnings=None):
    """Apply-to picker with live diff pane against the highlighted variant.

    Layout (top-down): the picker on top, a divider, the live diff
    against the currently-highlighted target below. Returns the chosen
    variant name to write to, or None if the user picked
    "Exit without saving" / hit Esc / Ctrl+C.

    `warnings` (optional list of strings) is rendered above the picker
    so it survives the full-screen-clear that prompt_toolkit performs
    on startup; the discovered-values block printed before this picker
    is gone by the time the TUI is up, so anything the user needs to
    see while choosing belongs here.
    """
    variants = _discover_variants()
    if not variants:
        print("No config.<variant>.yaml files found in the repo root.", file=sys.stderr)
        sys.exit(1)

    EXIT_KEY = "__exit__"
    items = []
    for v in variants:
        label = v + (" (active)" if v == active else "")
        items.append((v, label, ""))
    items.append((EXIT_KEY, "Exit without saving", "muted"))

    if not sys.stdin.isatty():
        # Non-interactive: write to active variant if it exists.
        if active and active in variants:
            return active
        print(
            "Cannot pick a variant non-interactively; pass --variant <name>.",
            file=sys.stderr,
        )
        sys.exit(2)

    from prompt_toolkit import Application
    from prompt_toolkit.key_binding import KeyBindings
    from prompt_toolkit.layout import Layout, Window, HSplit, FormattedTextControl
    from prompt_toolkit.layout.dimension import Dimension
    from prompt_toolkit.formatted_text import FormattedText
    from prompt_toolkit.styles import Style

    # Cache the per-variant `machine` block so we don't re-read YAML on
    # every keystroke — variants don't change during the picker.
    machine_cache = {}

    def _get_machine(variant_key):
        if variant_key in machine_cache:
            return machine_cache[variant_key]
        m, _path = _read_machine_block(variant_key)
        # ruamel CommentedMap → plain dict for diffing.
        m_plain = _ruamel_to_plain(m)
        machine_cache[variant_key] = m_plain
        return m_plain

    start = next((i for i, it in enumerate(items) if it[0] == active), 0)
    cursor = [start]
    chosen = [None]
    warnings = list(warnings or [])
    # Diff scroll offset (lines from the top), reset to 0 whenever
    # the picker cursor moves to a different variant.
    diff_top = [0]
    last_cursor = [cursor[0]]
    # Filled by _render_diff each frame so the scroll-key bindings
    # know how many lines exist + how many are visible.
    diff_total_lines = [0]
    diff_visible_lines = [10]

    def _render_warnings():
        if not warnings:
            return FormattedText([])
        out = [("class:warn_title", "Warnings:\n")]
        for w in warnings:
            out.append(("class:warn", f"  ! {w}\n"))
        out.append(("", "\n"))
        return FormattedText(out)

    def _render_picker():
        lines = [(
            "class:title",
            "Apply discovered values to: (↑/↓ to move, Enter to write, "
            "q / Esc to cancel; PgUp/PgDn or j/k or Ctrl-U/Ctrl-D scrolls diff)\n\n",
        )]
        for i, (_key, label, item_style) in enumerate(items):
            if i == cursor[0]:
                lines.append(("class:cursor", f"  > {label}\n"))
            else:
                style_class = f"class:{item_style}" if item_style else ""
                lines.append((style_class, f"    {label}\n"))
        return FormattedText(lines)

    def _build_diff_lines():
        """Return the full list of (style, line) tuples for the current target."""
        key = items[cursor[0]][0]
        if key == EXIT_KEY:
            return [("class:muted", "(no changes will be written)\n")]
        old = _get_machine(key)
        new = _machine_for_write(new_machine, old)
        diffs = list(_diff_machine(old, new))
        if not diffs:
            return [
                ("class:title", f"Diff against {key}\n"),
                ("", "\n"),
                ("class:muted",
                 f"  config.{key}.yaml already matches the discovered "
                 "values — nothing to write.\n"),
            ]
        out = [("class:title", f"Diff against {key}\n"), ("", "\n")]
        for k, ov, nv in diffs:
            out.append(("class:diff_key", f"  machine.{k}:\n"))
            out.append(("class:diff_label", "    old:\n"))
            for line in _yaml_value(ov, indent=6).splitlines(keepends=True):
                out.append(("class:diff_old", line))
            out.append(("class:diff_label", "    new:\n"))
            for line in _yaml_value(nv, indent=6).splitlines(keepends=True):
                out.append(("class:diff_new", line))
        return out

    def _render_diff():
        # Reset the scroll offset whenever the picker moves to a new
        # variant; otherwise the user would land mid-diff after
        # scrolling on the previous one.
        if last_cursor[0] != cursor[0]:
            last_cursor[0] = cursor[0]
            diff_top[0] = 0
        full = _build_diff_lines()
        diff_total_lines[0] = len(full)
        height = max(1, diff_visible_lines[0])

        # Body height = visible rows minus reserved hint rows. A hint
        # row appears when there's content off-screen on that side.
        # Compute body_height *first* (anchored at top=0 not scrolling
        # → no hints; otherwise reserve 1-2 rows), then clamp `top` so
        # the last body_height-sized window can actually reach the
        # bottom of `full`. The previous version clamped top to
        # `total - height`, which left up to 2 lines unreachable.
        body_height = height
        # Initial pass to figure out how many hint rows we'll have at
        # the eventual top. Pessimistic: assume both hints, then refine.
        if len(full) > height:
            body_height = max(1, height - 2)
        max_top = max(0, len(full) - body_height)
        top = max(0, min(diff_top[0], max_top))
        diff_top[0] = top
        # Now compute the actual hints based on the clamped top.
        has_above = top > 0
        has_below = (top + body_height) < len(full)
        # Recover an extra row when one side has no hint.
        body_height = height - (1 if has_above else 0) - (1 if has_below else 0)
        body_height = max(1, body_height)
        # Re-clamp now that body_height may have grown.
        max_top = max(0, len(full) - body_height)
        if top > max_top:
            top = max_top
            diff_top[0] = top
            has_above = top > 0
            has_below = (top + body_height) < len(full)

        out = []
        if has_above:
            out.append((
                "class:scroll_hint",
                f"  ↑ {top} more above (PgUp / k / Ctrl-U)\n",
            ))
        out.extend(full[top:top + body_height])
        if has_below:
            remaining = len(full) - (top + body_height)
            out.append((
                "class:scroll_hint",
                f"  ↓ {remaining} more below (PgDn / j / Ctrl-D)\n",
            ))
        return FormattedText(out)

    # Height of the warnings pane: 0 if no warnings, else title + N + blank.
    warn_height = 0 if not warnings else len(warnings) + 2

    warn_window = Window(
        content=FormattedTextControl(_render_warnings),
        height=Dimension.exact(warn_height) if warn_height else Dimension.exact(0),
        always_hide_cursor=True,
        wrap_lines=True,
    )
    picker_window = Window(
        content=FormattedTextControl(_render_picker),
        height=Dimension(min=len(items) + 2, preferred=len(items) + 3),
        always_hide_cursor=True,
    )
    divider = Window(
        content=FormattedTextControl(lambda: FormattedText([("class:divider", "─" * 60 + "\n")])),
        height=1,
    )
    diff_window = Window(
        content=FormattedTextControl(_render_diff),
        always_hide_cursor=True,
        wrap_lines=False,
    )

    def _refresh_visible_height():
        """Read the diff window's actual rendered height after each frame.

        prompt_toolkit only knows the window's pixel-row count after
        the layout has been measured. Pulling it from `render_info`
        each frame keeps the scroll math correct when the terminal is
        resized.
        """
        ri = diff_window.render_info
        if ri is not None:
            diff_visible_lines[0] = max(1, ri.window_height)

    kb = KeyBindings()

    @kb.add("up")
    def _(event):
        cursor[0] = (cursor[0] - 1) % len(items)

    @kb.add("down")
    def _(event):
        cursor[0] = (cursor[0] + 1) % len(items)

    def _scroll(delta):
        # Bump the offset; _render_diff is the single source of truth
        # for clamping (it knows the body height after reserving hint
        # rows). Set a generous upper bound here so we don't need to
        # duplicate that logic.
        _refresh_visible_height()
        diff_top[0] = max(0, diff_top[0] + delta)

    @kb.add("pageup")
    @kb.add("c-u")
    def _(event):
        _refresh_visible_height()
        _scroll(-max(1, diff_visible_lines[0] - 1))

    @kb.add("pagedown")
    @kb.add("c-d")
    def _(event):
        _refresh_visible_height()
        _scroll(+max(1, diff_visible_lines[0] - 1))

    @kb.add("k")
    def _(event):
        _scroll(-1)

    @kb.add("j")
    def _(event):
        _scroll(+1)

    @kb.add("home")
    def _(event):
        diff_top[0] = 0

    @kb.add("end")
    def _(event):
        # Big number → _render_diff clamps to the true bottom on next paint.
        _refresh_visible_height()
        diff_top[0] = diff_total_lines[0]

    @kb.add("enter")
    def _(event):
        key = items[cursor[0]][0]
        chosen[0] = None if key == EXIT_KEY else key
        event.app.exit()

    @kb.add("q")
    @kb.add("escape")
    @kb.add("c-c")
    def _(event):
        chosen[0] = None
        event.app.exit()

    style = Style.from_dict({
        "title": "bold",
        "cursor": "reverse",
        "muted": "fg:ansigray italic",
        "divider": "fg:ansigray",
        "diff_key": "bold",
        "diff_label": "fg:ansibrightblack",
        "diff_old": "fg:ansired",
        "diff_new": "fg:ansigreen",
        "warn_title": "bold fg:ansiyellow",
        "warn": "fg:ansiyellow",
        "scroll_hint": "fg:ansibrightblack italic",
    })

    layout_panes = []
    if warn_height:
        layout_panes.append(warn_window)
    layout_panes.extend([picker_window, divider, diff_window])

    app = Application(
        layout=Layout(HSplit(layout_panes)),
        key_bindings=kb,
        full_screen=True,
        style=style,
    )
    app.run()
    return chosen[0]


def _ruamel_to_plain(node):
    """Recursively coerce ruamel CommentedMap/CommentedSeq to plain dict/list."""
    from ruamel.yaml.comments import CommentedMap, CommentedSeq
    if isinstance(node, (CommentedMap, dict)):
        return {k: _ruamel_to_plain(v) for k, v in node.items()}
    if isinstance(node, (CommentedSeq, list)):
        return [_ruamel_to_plain(v) for v in node]
    return node


def _cmd_auto(args):
    """biopipelines-config auto — probe host + apply to a variant.

    Probes the host first, prints the discovered values, then pops a
    full-screen picker that lets the user choose which variant to
    write to (or exit without saving). The picker shows a live diff
    of the discovered values against the highlighted target as the
    cursor moves between variants.

    With ``--variant <name>`` the picker is skipped — values are
    written straight into that variant (with a .bak). With
    ``--verbose`` / ``-v`` each probe step is streamed to stderr.
    """
    from .config_manager import ConfigManager
    from ruamel.yaml import YAML

    args = list(args)
    verbose = False
    if "--verbose" in args:
        args.remove("--verbose")
        verbose = True
    if "-v" in args:
        args.remove("-v")
        verbose = True
    variant = _parse_variant_flag(args)
    if args:
        print(f"Unexpected argument: {args[0]}", file=sys.stderr)
        sys.exit(2)

    ctx = _ProbeContext(verbose=verbose)
    new_machine, sources = _build_machine_block(ctx)

    print("Discovered values:")
    print(f"  username           = {new_machine['username']!r}")
    print(f"  emails             = {dict(new_machine['emails'])}")
    em = new_machine['env_manager']
    if em is None:
        print(f"  env_manager        = <unchanged — probe failed>")
    else:
        print(f"  env_manager.name   = {em['name']!r}")
        for line in em['init']:
            print(f"  env_manager.init   = {line!r}")
    print(f"  scheduler.name     = {new_machine['scheduler']['name']!r}")
    for line in new_machine['scheduler']['init']:
        print(f"  scheduler.init     = {line!r}")
    print(f"  scheduler.modules  = {list(new_machine['scheduler']['modules'])}")
    print(f"  container_executor = {new_machine['container_executor']!r}")
    print()
    print("Probe sources:")
    for field, src in sources:
        print(f"  {field:<18} ← {src}")
    print()

    if ctx.warnings:
        print("Warnings:")
        for w in ctx.warnings:
            print(f"  ! {w}")
        print()

    if variant is None:
        try:
            active = ConfigManager().get_variant()
        except Exception:
            active = None
        variant = _pick_apply_target(new_machine, active, warnings=ctx.warnings)
        if variant is None:
            print("Exited without saving.")
            return

    path = ConfigManager._get_config_path(variant=variant)
    if not os.path.isfile(path):
        print(f"Config file not found: {path}", file=sys.stderr)
        sys.exit(1)

    yaml = YAML()
    yaml.preserve_quotes = True
    with open(path, "r", encoding="utf-8") as f:
        doc = yaml.load(f)
    if doc is None:
        from ruamel.yaml.comments import CommentedMap
        doc = CommentedMap()

    old_machine = doc.get("machine") or {}
    new_machine_for_write = _machine_for_write(new_machine, old_machine)

    diffs = list(_diff_machine(_ruamel_to_plain(old_machine), new_machine_for_write))
    if not diffs:
        print(f"config.{variant}.yaml already matches — nothing written.")
        return

    # Write through ruamel.yaml so comments / key ordering survive.
    from ruamel.yaml.comments import CommentedMap
    machine_cm = CommentedMap()
    for k, v in new_machine_for_write.items():
        if isinstance(v, dict):
            sub = CommentedMap()
            for sk, sv in v.items():
                sub[sk] = sv
            machine_cm[k] = sub
        else:
            machine_cm[k] = v
    doc["machine"] = machine_cm

    bak = path + ".bak"
    import shutil
    shutil.copy2(path, bak)
    with open(path, "w", encoding="utf-8") as f:
        yaml.dump(doc, f)
    print(f"Wrote {path} (backup: {bak})")


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


def _cmd_machine(args):
    """Handle 'machine' subcommand (also accessible as legacy alias 'cluster')."""
    from .config_manager import ConfigManager

    config_manager = ConfigManager()

    machine_values = {
        "variant": config_manager.get_variant(),
        "env_manager": config_manager.get_env_manager(),
        "scheduler": config_manager.get_scheduler(),
        "container_executor": config_manager.get_container_executor(),
        "slurm_modules": " ".join(config_manager.get_slurm_modules()),
        "shell_hook": config_manager.get_shell_hook_command(),
        "module_load": config_manager.get_module_load_line(),
        "activate": config_manager.get_activate_command("ENV_NAME"),
    }

    if not args:
        # List all machine settings
        for key, value in machine_values.items():
            print(f"{key}: {value}")
        return

    key = args[0]
    if key == "emails":
        emails = config_manager.get_emails()
        for user, email in emails.items():
            print(f"{user}: {email}")
    elif key in machine_values:
        print(machine_values[key])
    else:
        print(f"Unknown machine key: {key}", file=sys.stderr)
        print(f"Available keys: {', '.join(list(machine_values.keys()) + ['emails'])}", file=sys.stderr)
        sys.exit(1)
