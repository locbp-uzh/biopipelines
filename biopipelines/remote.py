"""Reach a cluster's outputs from the machine the agent runs on.

The cluster is the main platform, so the tools that read a run have to work against it, not
only against a local checkout. The server runs on the user's laptop and goes over the ssh alias
they already configured — no install rights on the cluster, no long-lived process on a shared
login node, credentials staying in their own agent.

**The remote resolves its own output root.** `biopipelines_output` is defined in the config on
each machine and expands `<username>` there, so asking the local config where the cluster keeps
its runs gives the wrong answer. `Ssh.output_root()` runs the resolution on the far side.

`LocalFS` and `Ssh` expose the same handful of operations, so `job_status` reads a run the same
way whether it is on this disk or on a login node.

Read-only queries are deliberately **not** written to a run's operations log. Recording every
status poll is what made `llm/log.sh` useless — a `squeue` and a six-hour submission became the
same kind of line. Only actions that change something are worth a record.
"""

import io
import json
import os
import pathlib
import re
import shlex
import subprocess

DEFAULT_TIMEOUT = 30
DEFAULT_REPO = "~/biopipelines"

# Resolve the remote's own `biopipelines_output`, creating nothing.
#
# Deliberately self-contained: it uses only `ConfigManager`, which every BioPipelines checkout
# has had for a long time. Calling `folders.output_root()` instead would make the server demand
# that the cluster run the same version as the laptop — the first live test against S3IT failed
# exactly that way, with `ImportError: cannot import name 'output_root'` on an older checkout.
#
# The sentinel prefix separates the answer from the config manager's variant banner.
_ROOT_SNIPPET = """
import getpass, os, re
import biopipelines
from biopipelines.config_manager import ConfigManager

resolved = {
    "username": getpass.getuser(),
    "cwd": os.path.dirname(os.path.dirname(os.path.abspath(biopipelines.__file__))),
}
# `base` is ordered so each entry only refers to earlier ones, which is why resolving in
# iteration order works without a second pass.
for key, template in ConfigManager().get_folder_config().get("base", {}).items():
    for _ in range(10):
        if "<" not in template:
            break
        template = re.sub(r"<([^>]+)>",
                          lambda m: resolved.get(m.group(1), m.group(0)), template)
    resolved[key] = template
print("BP_OUTPUT_ROOT=" + resolved.get("biopipelines_output", ""))
"""
_ROOT_SENTINEL = "BP_OUTPUT_ROOT="

# Stdlib only and no `biopipelines` import: ancestry must read a run on a cluster whose
# checkout is older than the caller's, and an import here is what makes the reader demand a
# version match. `csv` rather than cut/awk because column order varies per table and a quoted
# cell may contain commas.
_ID_COLUMNS_SNIPPET = r"""
import csv, fnmatch, json, os, sys

root, patterns = sys.argv[1], sys.argv[2].split("|")
extra = set(sys.argv[3].split("|")) if len(sys.argv) > 3 and sys.argv[3] else set()
out = {}
for dirpath, dirnames, filenames in os.walk(root):
    rel_dir = os.path.relpath(dirpath, root).replace("\\", "/")
    # Same reach as LocalFS: <step>/<file> and <step>/<stream>/<file>, never the job root.
    if rel_dir == ".":
        continue
    if rel_dir.count("/") >= 1:
        dirnames[:] = []
    for name in filenames:
        if not any(fnmatch.fnmatch(name, p) for p in patterns):
            continue
        rel = rel_dir + "/" + name
        try:
            with open(os.path.join(dirpath, name), newline="", encoding="utf-8",
                      errors="replace") as handle:
                reader = csv.DictReader(handle)
                cols = [c for c in (reader.fieldnames or [])
                        if c == "id" or c.endswith(".id") or c in extra]
                if "id" not in cols:
                    continue
                out[rel] = [{c: (row.get(c) or "") for c in cols} for row in reader]
        except (OSError, csv.Error):
            continue
print("BP_ID_COLUMNS=" + json.dumps(out))
"""
_ID_COLUMNS_SENTINEL = "BP_ID_COLUMNS="


def _decode_id_columns(text):
    """Pick the sentinel line out of remote stdout that may carry a config banner first."""
    for line in text.splitlines():
        if line.startswith(_ID_COLUMNS_SENTINEL):
            return json.loads(line[len(_ID_COLUMNS_SENTINEL):])
    return {}


class RemoteError(RuntimeError):
    pass


_LEADING_VARIABLE = re.compile(r"^\$(\{[A-Za-z_][A-Za-z0-9_]*\}|[A-Za-z_][A-Za-z0-9_]*)(?=/|$)")
_HOST = re.compile(r"^[A-Za-z0-9_][A-Za-z0-9._-]*(@[A-Za-z0-9_][A-Za-z0-9._-]*)?$")
# What an scp remote path may contain: a legacy-protocol scp hands the path to the remote shell.
_SCP_PATH = re.compile(r"^(\$\{?[A-Za-z_][A-Za-z0-9_]*\}?|~)?[A-Za-z0-9._/+@=,-]*$")

# BatchMode: a prompt the agent cannot answer must fail, not hang until the timeout.
SSH_OPTIONS = ["-o", "BatchMode=yes", "-o", "ConnectTimeout=15"]


def glob_escape(text):
    return re.sub(r"([*?\[])", r"[\1]", str(text))


def _decoded(data):
    if isinstance(data, bytes):
        return data.decode("utf-8", errors="replace")
    return data or ""


def quote(path):
    """Shell-quote a remote path, leaving a leading `~/` or `$VAR/` expandable.

    `shlex.quote("~/biopipelines")` returns it in single quotes, which stops the remote shell
    expanding the tilde — `cd '~/biopipelines'` fails with "No such file or directory". The same
    holds for Daint's `$SCRATCH/biopipelines`.
    """
    text = str(path)
    if text == "~":
        return "~"
    if text.startswith("~/"):
        return "~/" + shlex.quote(text[2:])
    variable = _LEADING_VARIABLE.match(text)
    if variable:
        rest = text[variable.end():]
        name = variable.group(1).strip("{}")
        # `:?` fails an unset variable instead of turning `$SCRATCH/x` into `/x`.
        return '"${' + name + ':?}"' + (shlex.quote(rest) if rest else "")
    return shlex.quote(text)


def interpreter(python):
    """The saved interpreter as a command word: a path is quoted, a hand-written `mamba run -n x python` is kept as typed."""
    text = str(python)
    return text if " " in text.strip() else quote(text)


def validate_host(host):
    """An ssh destination, never an option: `-oProxyCommand=...` as a host runs on this machine."""
    if not isinstance(host, str) or not _HOST.match(host):
        raise RemoteError(f"not a valid ssh host alias: {host!r}")
    return host


def validate_scp_path(path):
    """A remote path safe to hand to scp, whose legacy protocol expands it in the remote shell."""
    text = str(path)
    if not text or not _SCP_PATH.match(text) or ".." in text.split("/"):
        raise RemoteError(
            f"refusing remote path {text!r}: only letters, digits, '._/+@=,-', a leading ~ or "
            f"$VAR, and no '..' components")
    return text


class LocalFS:
    """The same operations as `Ssh`, against this machine's filesystem."""

    def join(self, *parts):
        return os.path.join(*parts)

    def is_dir(self, path):
        return pathlib.Path(path).is_dir()

    def is_file(self, path):
        return pathlib.Path(path).is_file()

    def listdir(self, path):
        target = pathlib.Path(path)
        if not target.is_dir():
            return []
        return sorted((p.name, p.is_dir()) for p in target.iterdir())

    def read_text(self, path):
        return pathlib.Path(path).read_text(encoding="utf-8", errors="replace")

    def append_text(self, path, text):
        """Append to a file, creating it and its parent. Returns True on success."""
        try:
            target = pathlib.Path(path)
            target.parent.mkdir(parents=True, exist_ok=True)
            with io.open(target, "a", encoding="utf-8") as handle:
                handle.write(text)
            return True
        except OSError:
            return False

    def tally(self, root, names):
        """{relative path: data row count} for every file under `root` matching `names`.

        One call, because on a cluster the alternative is a round trip per step and a large
        campaign has dozens. Row counts exclude the CSV header.
        """
        base = pathlib.Path(root)
        found = {}
        for name in names:
            # Map tables sit in the stream's own subfolder (`002_DSSP/dssp/dssp_map.csv`), not
            # beside it, so a one-level glob finds nothing and reports a run as traceless.
            for path in list(base.glob(f"*/{name}")) + list(base.glob(f"*/*/{name}")):
                try:
                    with io.open(path, encoding="utf-8", errors="replace") as handle:
                        found[str(path.relative_to(base)).replace("\\", "/")] = max(
                            sum(1 for _ in handle) - 1, 0)
                except OSError:
                    continue
        return found

    def id_columns(self, root, names, extra=()):
        """{relative path: [{column: value}]} for the `id` and `<axis>.id` columns, plus `extra`.

        Ancestry needs the provenance columns of every map table; the file-path columns beside
        them are the bulk of a map table and are never part of an edge.
        """
        import csv
        base = pathlib.Path(root)
        found = {}
        for name in names:
            for path in list(base.glob(f"*/{name}")) + list(base.glob(f"*/*/{name}")):
                try:
                    with io.open(path, newline="", encoding="utf-8", errors="replace") as handle:
                        reader = csv.DictReader(handle)
                        cols = [c for c in (reader.fieldnames or [])
                                if c == "id" or c.endswith(".id") or c in extra]
                        if "id" not in cols:
                            continue
                        rel = str(path.relative_to(base)).replace("\\", "/")
                        found[rel] = [{c: (row.get(c) or "") for c in cols} for row in reader]
                except (OSError, csv.Error):
                    continue
        return found

    def output_root(self, local_output=False):
        from biopipelines.folders import output_root as resolve
        return resolve(local_output=local_output)

    def find_job(self, root, job):
        """Every `<root>/<project>/<job>` directory."""
        base = pathlib.Path(root)
        if not base.is_dir():
            return []
        return sorted(str(p) for p in base.glob(f"*/{glob_escape(job)}") if p.is_dir())


class Ssh:
    """A cluster reached over an ssh alias from `~/.ssh/config`.

    Each call is a round trip, so the callers here batch what they can: reading a run's state
    costs three calls, not one per file.
    """

    def __init__(self, host, repo=DEFAULT_REPO, timeout=DEFAULT_TIMEOUT, python="python",
                 variant=None, prelude=None):
        self.host = validate_host(host)
        self.repo = repo
        self.timeout = timeout
        self.python = python
        # Shell to run before anything that needs the environment. Usually unnecessary: an
        # env's own binaries work by absolute path without activation, which is what `script()`
        # relies on. Kept for a site where they do not — an env behind a module, say.
        self.prelude = prelude
        # Which config variant this machine answers as. Daint is the case that needs it: the
        # resolver and `./submit` both call the framework, which otherwise takes whatever the
        # remote auto-detects, and generating a Daint run against the `cluster` variant produces
        # scripts for a scheduler layout that is not there.
        self.variant = variant

    @property
    def env_prefix(self):
        """What has to precede a command that calls the framework: prelude, then the variant.

        `ssh host "<command>"` is a non-interactive, non-login shell: it reads neither
        `~/.bashrc`'s interactive half nor `/etc/profile.d`, so it gets the bare system PATH.
        On S3IT that is `/usr/bin` and nothing else — no conda, no mamba, no console scripts.
        """
        parts = []
        if self.prelude:
            parts.append(self.prelude.rstrip("; &") + " && ")
        if self.variant:
            parts.append(f"BIOPIPELINES_CONFIG_VARIANT={quote(self.variant)} ")
        return "".join(parts)

    def script(self, name):
        """Where a console script lives for this host: beside its interpreter, if it is there.

        `bp-visualize` is installed into the environment, and a bare ssh cannot see it. An env's
        binaries do work by absolute path with no activation — verified on S3IT — so pointing a
        host's `python` at its env is enough to reach every script the env installed. Falls back
        to the bare name, which is right when `python` is the system one and the script was
        installed system-wide.
        """
        interpreter = str(self.python)
        if "/" not in interpreter:
            return name
        return interpreter.rsplit("/", 1)[0] + "/" + name

    def has_script(self, name):
        """Is that console script actually there? A missing one is a tool that cannot run."""
        resolved = self.script(name)
        if "/" not in resolved:
            code, _, _ = self.run(f"command -v {quote(resolved)} >/dev/null")
        else:
            code, _, _ = self.run(f"test -x {quote(resolved)}")
        return code == 0

    # --- transport -------------------------------------------------------

    def run(self, command, timeout=None, stdin=None):
        """Run a shell command on the remote. Returns (exit_code, stdout, stderr)."""
        try:
            # Bytes both ways: text-mode stdin on Windows writes CRLF into cluster files, and the
            # locale codec (cp1252) dies on the first non-ASCII byte of a cluster log.
            done = subprocess.run(["ssh", *SSH_OPTIONS, "--", self.host, command],
                                  input=stdin.encode("utf-8") if stdin is not None else None,
                                  capture_output=True, timeout=timeout or self.timeout)
        except FileNotFoundError:
            raise RemoteError("ssh is not on PATH")
        except subprocess.TimeoutExpired:
            raise RemoteError(f"ssh {self.host} timed out after {timeout or self.timeout}s")
        out, err = _decoded(done.stdout), _decoded(done.stderr)
        # 255 is ssh's own failure; letting it read as "absent" answers questions about a machine never reached.
        if done.returncode == 255:
            raise RemoteError(f"ssh {self.host} failed: {err.strip() or 'exit 255'}")
        return done.returncode, out, err

    def check(self):
        """Is the alias reachable? Returns (ok, detail) rather than raising, so a tool can say so."""
        try:
            code, out, err = self.run("echo ok", timeout=10)
        except RemoteError as exc:
            return False, str(exc)
        if code != 0:
            return False, (err or out).strip() or f"ssh exited {code}"
        return True, out.strip()

    # --- paths -----------------------------------------------------------

    def join(self, *parts):
        return "/".join(str(p).rstrip("/") for p in parts if str(p) != "")

    def output_root(self, local_output=False):
        """Ask the cluster where it writes runs, resolved by its own config."""
        command = (f"cd {quote(self.repo)} && {self.env_prefix}"
                   f"{interpreter(self.python)} -c {shlex.quote(_ROOT_SNIPPET)}")
        code, out, err = self.run(command)
        # The config manager prints a variant banner first, so pick the sentinel line rather
        # than trusting position.
        root = next((line[len(_ROOT_SENTINEL):].strip()
                     for line in out.splitlines() if line.startswith(_ROOT_SENTINEL)), "")
        if code != 0 or not root:
            raise RemoteError(
                f"could not resolve the output root on {self.host}: "
                f"{(err or out).strip() or f'exit {code}'}")
        return root

    # --- filesystem ------------------------------------------------------

    def is_dir(self, path):
        code, _, _ = self.run(f"test -d {quote(path)}")
        return code == 0

    def is_file(self, path):
        code, _, _ = self.run(f"test -f {quote(path)}")
        return code == 0

    def listdir(self, path):
        """(name, is_dir) for each entry. `ls -Ap` marks directories with a trailing slash."""
        code, out, _ = self.run(f"ls -Ap {quote(path)} 2>/dev/null")
        if code != 0:
            return []
        entries = []
        for line in out.splitlines():
            name = line.rstrip("\n")
            if not name:
                continue
            entries.append((name[:-1], True) if name.endswith("/") else (name, False))
        return sorted(entries)

    def read_text(self, path):
        code, out, err = self.run(f"cat {quote(path)}")
        if code != 0:
            raise RemoteError(f"cannot read {path} on {self.host}: {err.strip()}")
        return out

    def append_text(self, path, text):
        """Append to a file on the cluster, creating its parent.

        The operations log has to live with the run, and the run is on the cluster — writing it
        through local pathlib silently produced nothing, which is how the first real submission
        ended up with no record at all.
        """
        path = str(path)
        parent = path.rsplit("/", 1)[0] if "/" in path else "."
        code, _, _ = self.run(
            f"mkdir -p {quote(parent)} && cat >> {quote(path)}",
            stdin=text)
        return code == 0

    def tally(self, root, names):
        """{relative path: data row count} for matching files under `root`, in one round trip."""
        patterns = " -o ".join(f"-name {quote(n)}" for n in names)
        command = (f"cd {quote(root)} && find . -mindepth 2 -maxdepth 3 " + r"\( " + patterns + r" \) "
                   + "-exec wc -l {} + 2>/dev/null")
        code, out, _ = self.run(command, timeout=max(self.timeout, 120))
        if code != 0:
            return {}
        found = {}
        for line in out.splitlines():
            parts = line.strip().split(None, 1)
            if len(parts) != 2 or not parts[0].isdigit():
                continue
            count, path = int(parts[0]), parts[1].lstrip("./")
            if path == "total":
                continue
            found[path] = max(count - 1, 0)
        return found

    def id_columns(self, root, names, extra=()):
        """The provenance columns of every map table under `root`, plus `extra`, in one round trip.

        Sends a stdlib-only reader rather than shipping whole map tables back: the file-path
        column dominates their size and is never part of an edge.
        """
        command = (f"{self.env_prefix}{interpreter(self.python)} -c {shlex.quote(_ID_COLUMNS_SNIPPET)} "
                   f"{quote(root)} {shlex.quote('|'.join(names))} {shlex.quote('|'.join(extra))}")
        code, out, err = self.run(command, timeout=max(self.timeout, 120))
        if code != 0:
            raise RemoteError(
                f"could not read provenance columns under {root} on {self.host}: "
                f"{(err or out).strip() or f'exit {code}'}")
        return _decode_id_columns(out)

    def expand(self, path):
        """A leading `$VAR` resolved on the remote: scp's SFTP mode expands `~` but not variables."""
        text = validate_scp_path(path)
        if not _LEADING_VARIABLE.match(text):
            return text
        code, out, err = self.run(f"printf '%s' {quote(text)}")
        if code != 0 or not out.strip():
            raise RemoteError(f"could not resolve {text} on {self.host}: {err.strip() or 'empty'}")
        return validate_scp_path(out.strip())

    def find_job(self, root, job):
        """Every `<root>/<project>/<job>` directory, in one round trip rather than one per project."""
        code, out, _ = self.run(
            f"cd {quote(root)} && for d in */{shlex.quote(job)}; do [ -d \"$d\" ] && echo \"$d\"; done; true")
        if code != 0:
            return []
        return sorted(self.join(root, line.strip()) for line in out.splitlines() if line.strip())

    def upload(self, local, remote_path, timeout=None):
        """scp a local file to the cluster. Returns (exit_code, stdout, stderr) like `run`."""
        remote_path = self.expand(remote_path)
        try:
            done = subprocess.run(["scp", *SSH_OPTIONS, "--", str(local), f"{self.host}:{remote_path}"],
                                  capture_output=True, text=True,
                                  encoding="utf-8", errors="replace",
                                  timeout=timeout or max(self.timeout, 120))
        except FileNotFoundError:
            raise RemoteError("scp is not on PATH")
        except subprocess.TimeoutExpired:
            raise RemoteError(f"scp to {self.host} timed out")
        return done.returncode, done.stdout or "", done.stderr or ""

    def download(self, remote_path, local, timeout=None, recursive=None):
        """scp back from the cluster — one file, or a directory with `recursive`.

        `recursive=None` decides by asking the remote, so a caller does not have to know which
        it is. A run's structures live in a step folder, which is the common case for pulling
        results down to look at them.
        """
        remote_path = self.expand(remote_path)
        if recursive is None:
            recursive = self.is_dir(remote_path)
        argv = ["scp", *SSH_OPTIONS] + (["-r"] if recursive else []) + [
            "--", f"{self.host}:{remote_path}", str(local)]
        try:
            done = subprocess.run(argv, capture_output=True, text=True,
                                  encoding="utf-8", errors="replace",
                                  timeout=timeout or max(self.timeout, 600))
        except FileNotFoundError:
            raise RemoteError("scp is not on PATH")
        except subprocess.TimeoutExpired:
            raise RemoteError(f"scp from {self.host} timed out")
        return done.returncode, done.stdout or "", done.stderr or ""


def filesystem(host=None, repo=DEFAULT_REPO, **kwargs):
    """`LocalFS` when no host is given, an `Ssh` otherwise — what every caller should use."""
    return LocalFS() if not host else Ssh(host, repo=repo, **kwargs)


def _is_absolute(text):
    """A path the caller already resolved: POSIX root, `~`, or a Windows drive."""
    return text.startswith(("/", "~", "\\\\")) or re.match(r"^[A-Za-z]:[\\/]", text) is not None


def resolve_job_dir(fs, job, project=None, root=None, local_output=False):
    """Turn a job name into a full path, asking the machine that owns it where its outputs live.

    An absolute `job` is used as given. Anything else is taken as relative to the output root:
    `<project>/<job>` when it carries a separator, `<output_root>/<project>/<job>` when a
    project is named, otherwise the projects are searched one level down for that job name.
    """
    text = str(job)
    if _is_absolute(text):
        return text

    root = root or fs.output_root(local_output=local_output)
    if project:
        return fs.join(root, project, text)
    if "/" in text or "\\" in text:
        # `bp_runs` prints "2 runs in Template", so `Template/example_002` is what a caller
        # types next. Returning it unchanged made it a relative path that resolves to nothing,
        # and the run reported as missing while `bp_runs` had just listed it.
        return fs.join(root, *re.split(r"[\\/]+", text.strip("/\\")))

    matches = fs.find_job(root, text)
    if len(matches) > 1:
        # A cancel or resubmit on the wrong project's run is the failure this refuses.
        raise RemoteError(f"job {text!r} exists in {len(matches)} projects; name one with "
                          f"project= or pass <project>/{text}: " + ", ".join(matches))
    if matches:
        return matches[0]
    return fs.join(root, text)


# --- saved connections ----------------------------------------------------
#
# A first-time user with only the MCP tools has no way to discover their cluster's ssh alias or
# repo path, and repeating both on every call is noise. `bp_setup` verifies a connection and
# saves it here; every other tool falls back to it. Gitignored — it names one person's machines.
#
# Settings are keyed BY HOST, not global. A lab with access to two clusters is the ordinary
# case, and the flat version of this file carried one repo for whichever host you named — so
# `bp_status(job=..., host="daint")` sent S3IT's `~/biopipelines-locbp` to Daint, where the
# checkout lives on $SCRATCH. That resolves to nothing, or worse to something else's tree.
# Each host also needs its own interpreter (Daint's login nodes have no `python` at all) and
# its own config variant, neither of which a single slot can hold.

SETTINGS_FILE = pathlib.Path(__file__).resolve().parent.parent / ".bp-mcp.json"

PER_HOST_KEYS = ("repo", "python", "variant", "prelude")


def load_settings():
    """The saved connections as {"default": host, "hosts": {host: {...}}}.

    A file written by the flat version is read as one host so nobody has to run setup again.
    """
    import json
    try:
        raw = json.loads(SETTINGS_FILE.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {"default": None, "hosts": {}}
    if not isinstance(raw, dict):
        return {"default": None, "hosts": {}}
    if "hosts" in raw:
        raw.setdefault("default", next(iter(raw["hosts"]), None))
        return raw
    host = raw.get("host")
    if not host:
        return {"default": None, "hosts": {}}
    return {"default": host,
            "hosts": {host: {k: raw[k] for k in PER_HOST_KEYS if raw.get(k)}}}


def save_settings(host=None, make_default=True, **values):
    """Record what works for one host. Other hosts are left exactly as they were."""
    import json
    current = load_settings()
    if host:
        entry = dict(current["hosts"].get(host) or {})
        entry.update({k: v for k, v in values.items() if k in PER_HOST_KEYS and v is not None})
        current["hosts"][host] = entry
        if make_default or not current.get("default"):
            current["default"] = host
    try:
        SETTINGS_FILE.write_text(json.dumps(current, indent=2) + "\n", encoding="utf-8")
    except OSError:
        return None
    return current


def settings_for(host):
    """What is saved for one host: {} when nothing is."""
    return dict(load_settings()["hosts"].get(host) or {})


def forget_settings(host):
    """Drop one host. Returns True if it was there.

    Once settings hold several hosts they accumulate: a renamed alias or a typo would otherwise
    stay in the listing forever, and a stale entry is worse than none — it names a machine that
    may no longer answer.
    """
    import json
    current = load_settings()
    if host not in current["hosts"]:
        return False
    del current["hosts"][host]
    if current.get("default") == host:
        current["default"] = next(iter(current["hosts"]), None)
    try:
        SETTINGS_FILE.write_text(json.dumps(current, indent=2) + "\n", encoding="utf-8")
    except OSError:
        return False
    return True


def connection(host=None, repo=None, python=None, variant=None, prelude=None):
    """Everything needed to reach one machine: (host, repo, python, variant).

    Only `None` means "not given" for `host`. An explicit empty host is the escape hatch for
    reading this machine while a cluster is saved, so it must not fall through to the saved one.

    Nothing is inherited across hosts. Naming a host we know nothing about yields the default
    repo rather than another cluster's, because a wrong path that happens to exist is the one
    failure that produces a confident answer about the wrong machine.
    """
    saved = load_settings()
    host = saved.get("default") if host is None else host
    entry = dict(saved["hosts"].get(host) or {}) if host else {}
    return (host,
            repo or entry.get("repo") or DEFAULT_REPO,
            python or entry.get("python") or "python",
            variant or entry.get("variant") or None,
            prelude or entry.get("prelude") or None)


def defaults(host=None, repo=None):
    """(host, repo) for the callers that need no interpreter — see `connection`."""
    host, repo = connection(host, repo)[:2]
    return host, repo


def diagnose(host, repo=DEFAULT_REPO, python="python", variant=None, prelude=None):
    """Walk the chain a first run needs, stopping at the first thing that is not in place.

    Returns (ok, [(stage, ok, detail)]). The stages are ordered so the first failure is the
    thing to fix: a missing alias is reported as a missing alias, not as a config error three
    steps later.
    """
    ssh = Ssh(host, repo=repo, python=python, variant=variant, prelude=prelude)
    stages = []

    ok, detail = ssh.check()
    stages.append(("ssh alias reachable", ok, detail))
    if not ok:
        return False, stages

    code, out, _ = ssh.run(f"test -d {quote(repo)} && echo present")
    ok = code == 0 and "present" in out
    stages.append((f"repository at {repo}", ok,
                   "found" if ok else "not found — clone it there, or pass the right `repo`"))
    if not ok:
        return False, stages

    # Deliberately NOT from inside the repo. `cd <repo> && python -c "import biopipelines"`
    # puts the checkout on sys.path, so the check passed under /usr/bin/python importing the
    # source tree -- a green line that proved the clone exists and nothing about the
    # environment a run actually uses.
    code, out, err = ssh.run(
        f"cd ~ && {interpreter(python)} -c 'import biopipelines, sys; "
        f"print(biopipelines.__version__, sys.executable)'")
    ok = code == 0
    if ok:
        parts = out.strip().split()
        detail = f"{parts[0]} via {parts[-1]}" if len(parts) > 1 else out.strip()
    else:
        detail = ((err.strip().splitlines() or ["failed"])[-1]
                  + " — pass `python` for this host: the interpreter that has biopipelines "
                    "installed, e.g. <conda-root>/envs/biopipelines/bin/python")
    stages.append(("biopipelines importable", ok, detail))
    if not ok:
        return False, stages

    missing = [name for name in ("bp-visualize",) if not ssh.has_script(name)]
    stages.append(("console scripts on the interpreter's path", not missing,
                   "bp-visualize found" if not missing else
                   f"{', '.join(missing)} not beside {python} — bp_visualize cannot run. A bare "
                   "ssh gets only the system PATH, so the env's scripts are reachable only by "
                   "pointing `python` at that env."))

    try:
        root = ssh.output_root()
        stages.append(("config variant resolves", True, root))
    except RemoteError as exc:
        stages.append(("config variant resolves", False, str(exc).splitlines()[0]))
        return False, stages

    names = [n for n, is_dir in ssh.listdir(root) if is_dir]
    stages.append(("output root readable", True, f"{len(names)} projects"))
    return True, stages
