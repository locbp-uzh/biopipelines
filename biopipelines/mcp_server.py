"""`bp-mcp` — the MCP server that lets an agent operate BioPipelines instead of reading about it.

Runs on the user's machine over stdio and reaches a cluster over the configured ssh alias, so
it needs no install rights on the cluster and leaves no long-lived process on a shared login
node. Install the dependency with `pip install -e ".[mcp]"`.

`bp_tools` answers "what exists and what does it take", the question an agent otherwise answers
by reading 78k tokens of prose or by guessing a keyword that `**kwargs` will silently swallow.
`bp_status` and `bp_logs` answer "what happened to my run", which is where a half-failed
campaign otherwise turns into grepping. `bp_project` keeps the project folder's two documents.

The cluster is the main platform, so every run-reading tool takes a `host`: the ssh alias from
`~/.ssh/config`. With one, the tool asks that machine where its own `biopipelines_output` is —
the config resolves `<username>` on the machine that owns it, so the local config would answer
for the wrong user — and reads over ssh. Without one, it reads this machine's filesystem.
"""

import csv as _csv
import pathlib
import re
import sys
import tempfile

from biopipelines import (__version__, ancestry, job_status, lineage, manifest, project_docs,
                          remote, reproduce, submission, tool_docs, tool_schemas, tool_tags,
                          views)
from biopipelines.renderers import provenance_report
from biopipelines.tool_docs import collect, render, section, suggest

try:
    from mcp.server.mcpserver import MCPServer
    from mcp.server.mcpserver.utilities.types import Image
except ModuleNotFoundError:  # pragma: no cover - exercised only without the extra installed
    MCPServer = None
    Image = None

INSTRUCTIONS = """BioPipelines: 86 protein/ligand-modeling tools behind one declarative Pipeline API.

These tools ARE the interface to BioPipelines. Use them in preference to the shell and to the
repository's own files: do not `ssh` to submit a run, tail a log or list a directory, do not
`scp` a table back, do not open `docs/tool/*.md` to look up a signature, and do not run any
script under `skills/` or `versions/` — those are maintainer utilities, not the way in. Each
of those paths either loses the operations record these tools write, or re-derives something
they already report. Shell out only for what no tool here covers (`sinfo`, `squeue`), and only
after checking that none does.

Call `bp_tools` with no argument to get the catalog, then call it again with a tool name to get
that tool's full documentation before you write the stage that uses it. Do not guess a tool
name or a parameter: an unrecognized keyword is absorbed by **kwargs rather than rejected, so a
typo becomes a silently ignored argument instead of an error.

The cluster is the main platform. `bp_setup` connects to one and remembers it; after that
`bp_runs`, `bp_status` and `bp_logs` reach it without being told the host every time. Call
`bp_setup` first on a machine that has never been connected — it names whatever is missing
rather than failing three steps later.

`bp_submit` runs a pipeline there and `bp_fetch` pulls a file back. bp_submit spends compute:
call it only when the user has asked for the run. Nothing has executed when it returns — the
scheduler has merely accepted the job, so poll `bp_status` afterwards rather than reporting
success. After a run, `bp_status` reads its completion markers and `bp_logs` reads one step's log.
`bp_resubmit` resumes a run that failed partway and `bp_cancel` stops one; both spend or destroy
work, so both wait for the user to ask.

When a step finishes, call `bp_visualize` on it rather than describing its numbers in prose. It
returns the step's images inline and fetches the interactive page, and it needs no scheduler, so
a finished step can be shown while the rest of the campaign is still queued.

A project folder holds its job directories plus two documents: PROJECT.md (goal, inputs, and
the validated protocol — what we do and why we believe it) and HISTORY.md (append-only, dated:
runs, findings, restructures). `bp_project` surveys and scaffolds them, and appends history
entries. PROJECT.md carries judgment: draft and propose, but leave the content to the
scientist, and never rewrite what they wrote."""


def build_server():
    """The server object, so tests can exercise the tools without a transport."""
    if MCPServer is None:
        raise ModuleNotFoundError(
            "bp-mcp needs the MCP SDK: pip install -e \".[mcp]\"")

    server = MCPServer(name="biopipelines", version=__version__, instructions=INSTRUCTIONS)

    @server.tool(
        name="bp_tools",
        title="BioPipelines tool catalog",
        description=(
            "Find a BioPipelines tool, or read one's documentation. Call with no arguments for "
            "the catalog (one line per tool: version, CPU/GPU, verified platforms, purpose, and "
            "its upstream repo and paper). Call with `name` for that tool's complete entry, "
            "including its parameters — read this before writing a stage, because unrecognized "
            "keyword arguments are silently ignored rather than rejected. "
            "Call with `outputs` to find tools by what they RETURN — a stream, table or column "
            "name, case-insensitive substring, so `outputs=\"rmsd\"` finds `ligand_rmsd` and "
            "`RMSD_before`. Call with `inputs` to find tools by what they CONSUME, e.g. "
            "`inputs=\"structures\"`. Both are read from the tools' own source, so they are "
            "exact where a keyword search over prose is not. "
            "Call with `tags` to find tools by intent, from a closed 38-term vocabulary in four "
            "facets (action, subject, readout, capability): terms in the same facet are OR-ed, "
            "different facets AND-ed, so tags=[\"dock\", \"predict-structure\", \"covalent\"] "
            "means (docking or structure prediction) and covalent. `exclude` drops tools "
            "carrying any of its tags — exclude=[\"fetch\"] for what runs with no network, "
            "tags=[\"protein\", \"binding\"] exclude=[\"small-molecule\"] for protein-only "
            "binding tools. Pass an unknown tag to get the vocabulary printed."),
    )
    def bp_tools(name: str | None = None, outputs: str | None = None,
                 inputs: str | None = None, tags: list[str] | None = None,
                 exclude: list[str] | None = None) -> str:
        """The catalog, one tool's documentation, or a search by tag, output or input."""
        if outputs:
            return tool_schemas.summarize(outputs, where="outputs")
        if inputs:
            return tool_schemas.summarize(inputs, where="inputs")
        if tags or exclude:
            return _by_tags(tags or [], exclude or [])
        if name is None:
            return render(collect())
        found = section(name)
        if found is not None:
            return found
        close = suggest(name)
        hint = f" Closest documented names: {', '.join(close)}." if close else ""
        return (f"No BioPipelines tool named {name!r}.{hint}"
                " Call bp_tools with no arguments for the full catalog.")

    def _by_tags(tags, exclude):
        """Tools matching a tag query, or an error that teaches the vocabulary."""
        tags, renamed = tool_tags.resolve_query(tags)
        exclude, renamed_out = tool_tags.resolve_query(exclude)
        renamed.update(renamed_out)
        unknown = tool_tags.unknown_tags(list(tags) + list(exclude))
        if unknown:
            lines = [f"Not tags: {', '.join(unknown)}."]
            for tag in unknown:
                close = tool_tags.close_tags(tag)
                if close:
                    lines.append(f"  {tag} — did you mean {', '.join(close)}?")
            lines += ["", "The vocabulary is closed; these 38 terms are all of it:", "",
                      tool_tags.vocabulary_text(), "",
                      "Tags are the intent layer. For a named metric use bp_tools(outputs=\"rmsd\"), "
                      "which searches stream, table and column names from the tools' own source."]
            return "\n".join(lines)

        found = tool_tags.select(collect(), tags, exclude)
        query = " and ".join(filter(None, [
            f"tags {tags}" if tags else "", f"not {exclude}" if exclude else ""]))
        read_as = ("".join(f" (read {typed!r} as {meant!r})" for typed, meant in renamed.items())
                   if renamed else "")
        if not found:
            facets = len({tool_tags.FACET_OF[t] for t in tags})
            return (f"No tool matches {query}{read_as}. Terms in one facet are OR-ed and different "
                    f"facets AND-ed, so this query is the intersection of {facets} facet(s) — drop "
                    f"the least important term. The vocabulary:\n\n{tool_tags.vocabulary_text()}")
        header = (f"{len(found)} tool(s) matching {query}{read_as}. Tags are a discovery aid, not a "
                  f"guarantee: read the tool's own section with bp_tools(name=...) before using it.")
        return tool_docs.render_selection(found, header)

    def _fetch_to_scratch(fs, remote_path):
        """Copy one cluster file next to the others this server has pulled. None if it fails.

        A rendered page and its images have to exist locally before anything can open or embed
        them, and asking the user for a destination every time defeats the point of showing a
        result. They land under the OS temp dir, keyed by run, so nothing is written into the
        repository or into whatever directory the agent happened to start in.
        """
        if getattr(fs, "host", None) is None:
            return remote_path  # already this machine's filesystem
        # Keyed by host and the path's last three parts, so two runs' same-named views do not collide.
        parts = [re.sub(r"[^A-Za-z0-9._+-]", "_", p)
                 for p in str(remote_path).replace("\\", "/").split("/") if p not in ("", ".", "..")]
        target = pathlib.Path(tempfile.gettempdir()).joinpath("bp-mcp-views", fs.host, *parts[-3:])
        target.parent.mkdir(parents=True, exist_ok=True)
        try:
            code, _out, _err = fs.download(remote_path, str(target))
        except remote.RemoteError:
            return None
        return str(target) if code == 0 and target.exists() else None

    def _views_path(job_dir, suffix):
        """`<repo>/outputs/_views/<job><suffix>`: gitignored, findable, and independent of the server's cwd."""
        name = re.split(r"[\\/]+", str(job_dir).rstrip("/\\"))[-1]
        name = re.sub(r"[^A-Za-z0-9._+-]", "_", name) or "run"
        folder = pathlib.Path(__file__).resolve().parent.parent / "outputs" / "_views"
        folder.mkdir(parents=True, exist_ok=True)
        return folder / f"{name}{suffix}"

    def _connect(host, repo):
        """(host, filesystem) with that host's own repo, interpreter and config variant.

        Everything goes through here so a second cluster cannot be reached with the first's
        settings — the flat settings file carried one repo for whichever host was named.
        """
        host, repo, python, variant, prelude = remote.connection(host, repo)
        if not host:
            return host, remote.LocalFS()
        return host, remote.Ssh(host, repo=repo, python=python, variant=variant,
                                prelude=prelude)

    def _resolve(job, host, project, repo):
        """(filesystem, job path) for a run named or pathed, here or on a cluster."""
        _host, fs = _connect(host, repo)
        return fs, remote.resolve_job_dir(fs, job, project=project)

    @server.tool(
        name="bp_setup",
        title="Connect to a cluster",
        description=(
            "Check and remember how to reach a cluster. Pass `host` (the ssh alias from "
            "~/.ssh/config) and `repo` (the BioPipelines checkout there, default ~/biopipelines). "
            "Walks the chain a first run needs — alias reachable, repo present, biopipelines "
            "importable, config variant resolving, output root readable — and stops at the first "
            "thing that is not in place, so the report names the thing to fix. Settings are kept "
            "PER HOST — repo, interpreter and config variant — so a second cluster never "
            "replaces the first and nothing is inherited across them. Call with no arguments to "
            "see what is configured, `host` to add or re-check one, `make_default=False` to add "
            "one without changing which cluster unnamed calls use. Pass `variant` for a machine "
            "whose config variant is not what it auto-detects (Daint needs \"daint\"), and "
            "`python` where the interpreter is not simply `python` (Daint's login nodes have "
            "none) — point it at the environment's interpreter, not the system one, since that is "
            "also how the console scripts behind bp_visualize are found. `forget=True` with a "
            "host drops it, for a renamed alias or a typo. "
            "Call this first on a machine that has never been connected."),
    )
    def bp_setup(host: str | None = None, repo: str | None = None,
                 python: str | None = None, variant: str | None = None,
                 make_default: bool = True, save: bool = True, forget: bool = False) -> str:
        # `prelude` is raw shell run before every command, so it is set in .bp-mcp.json by hand, never through a tool call.
        prelude = None
        try:
            if host is not None:
                remote.validate_host(host)
            for label, value in (("python", python), ("repo", repo)):
                if value is not None:
                    remote.validate_scp_path(value)
            if variant is not None and not variant.replace("_", "").replace("-", "").isalnum():
                raise remote.RemoteError(f"not a config variant name: {variant!r}")
        except remote.RemoteError as exc:
            return f"Refused: {exc}"
        if forget and host:
            if remote.forget_settings(host):
                left = sorted(remote.load_settings()["hosts"]) or ["none"]
                return f"Forgot {host}. Still configured: {', '.join(left)}."
            return f"Nothing saved for {host}."
        saved = remote.load_settings()
        if host is None:
            if not saved["hosts"]:
                return ("No cluster is configured yet. Call bp_setup with `host` — the ssh alias "
                        "from ~/.ssh/config, e.g. host=\"cluster\" — and `repo` if the checkout "
                        "is not at ~/biopipelines.\n\nThe alias itself has to exist in "
                        "~/.ssh/config already; that file is outside what these tools can reach.")
            lines = ["Configured clusters:"]
            for name in sorted(saved["hosts"]):
                entry = saved["hosts"][name]
                mark = " (default)" if name == saved.get("default") else ""
                detail = f"repo={entry.get('repo', remote.DEFAULT_REPO)}"
                for key in ("python", "variant", "prelude"):
                    if entry.get(key):
                        detail += f" {key}={entry[key]}"
                lines.append(f"  {name}{mark}: {detail}")
            lines.append("\nEvery tool takes host= to reach one of these; without it they use the "
                         "default. Call bp_setup with host= to add or re-check a cluster.")
            return "\n".join(lines)

        host, repo, python, variant, prelude = remote.connection(
            host, repo, python, variant, prelude)
        try:
            ok, stages = remote.diagnose(host, repo=repo, python=python, variant=variant,
                                         prelude=prelude)
        except remote.RemoteError as exc:
            return f"Not ready: {host} ({repo})\n  [--] {exc}\n\nNothing was saved."
        lines = [f"{'Ready' if ok else 'Not ready'}: {host} ({repo})"]
        for name, good, detail in stages:
            lines.append(f"  [{'ok' if good else '--'}] {name}: {detail}")
        if ok and save:
            if remote.save_settings(host=host, repo=repo, python=python, variant=variant,
                                    prelude=prelude, make_default=make_default) is None:
                lines.append(f"\nReady, but NOT saved: {remote.SETTINGS_FILE} is not writable, "
                             f"so every call has to pass host= and repo= explicitly.")
                return "\n".join(lines)
            others = sorted(h for h in remote.load_settings()["hosts"] if h != host)
            note = f"\nSaved for {host}."
            if make_default:
                note += " Calls that name no host now use it."
            if others:
                plural = "is" if len(others) == 1 else "are"
                them = "it" if len(others) == 1 else "them"
                note += (f" {', '.join(others)} {plural} untouched — pass host= to reach {them}.")
            lines.append(note)
        elif not ok:
            lines.append("\nFix the first '--' line above, then call bp_setup again.")
        return "\n".join(lines)

    @server.tool(
        name="bp_runs",
        title="List runs",
        description=(
            "Projects and runs under the output root. Pass `host` (an ssh alias) for a cluster "
            "— that machine resolves its own output location. Omit `project` to list projects."),
    )
    def bp_runs(host: str | None = None, project: str | None = None,
                repo: str | None = None) -> str:
        try:
            host, fs = _connect(host, repo)
            root = fs.output_root()
        except Exception as exc:
            return f"Could not resolve the output root: {exc}"
        where = fs.join(root, project) if project else root
        try:
            names = [n for n, is_dir in fs.listdir(where) if is_dir]
        except remote.RemoteError as exc:
            return f"Could not list {where}: {exc}"
        if not names:
            return f"Nothing under {where}."
        label = f"runs in {project}" if project else "projects"
        return f"{len(names)} {label} under {where}:\n" + "\n".join("  " + n for n in names)

    @server.tool(
        name="bp_status",
        title="Run status",
        description=(
            "What happened to a run: every step with its completion marker, plus the run header "
            "if one was written. `job` is a job name (`<Job>_NNN`) or a full path; `host` is an "
            "ssh alias when the run is on a cluster. Read this before reporting a run as "
            "successful — a step can fail without halting the pipeline."),
    )
    def bp_status(job: str, host: str | None = None, project: str | None = None,
                  repo: str | None = None) -> str:
        try:
            fs, job_dir = _resolve(job, host, project, repo)
            return job_status.summarize(job_status.status(job_dir, fs=fs))
        except Exception as exc:
            return f"Could not read {job!r}{f' on {host}' if host else ''}: {exc}"

    @server.tool(
        name="bp_logs",
        title="Step log",
        description=(
            "The tail of one step's log. `step` is `001_ToolName`, or just the tool name when it "
            "appears once in the run. `host` is an ssh alias when the run is on a cluster. Use "
            "it on the failure bp_status names."),
    )
    def bp_logs(job: str, step: str, tail: int = job_status.LOG_TAIL_LINES,
                host: str | None = None, project: str | None = None,
                repo: str | None = None) -> str:
        try:
            fs, job_dir = _resolve(job, host, project, repo)
            found = job_status.log(job_dir, step, tail=tail, fs=fs)
            if found is None:
                available = [s["step"] for s in job_status.steps(job_dir, fs=fs) if s["has_log"]]
                listed = ", ".join(available) if available else "none"
                return f"No log for {step!r} in {job_dir}. Steps with logs: {listed}."
            head = f"{found['step']} — last {found['shown']} of {found['lines']} lines\n"
            return head + found["text"]
        except Exception as exc:
            return f"Could not read {job!r}{f' on {host}' if host else ''}: {exc}"

    @server.tool(
        name="bp_lineage",
        title="Where the designs went",
        description=(
            "Per step, how many IDs were produced per stream and how many were dropped, for a "
            "whole run. Answers 'I started with 500 and have 12 — which step removed the rest'. "
            "Read this before reporting a campaign's yield: the counts come from the map and "
            "missing tables the framework writes, not from arithmetic on result tables. "
            "Pass page=True to write a self-contained HTML page instead — provenance, attrition "
            "and every dropped ID with its reason — for someone who will not call these tools: "
            "a supervisor, a reviewer, a referee."),
    )
    def bp_lineage(job: str, page: bool = False, host: str | None = None,
                   project: str | None = None, repo: str | None = None) -> str:
        try:
            fs, job_dir = _resolve(job, host, project, repo)
            name = str(job_dir).rstrip("/").rsplit("/", 1)[-1]
            if page:
                out = str(_views_path(job_dir, "_provenance.html"))
                written = provenance_report.write(job_dir, out, fs=fs)
                return (f"Wrote {written}\n\nOpen it in a browser. It carries the run's recorded "
                        f"identity, the attrition step by step, the parameters as they resolved "
                        f"and every dropped ID with the reason its own step gave.")
            steps = lineage.collect(job_dir, fs=fs)
            return lineage.summarize(steps, job=name)
        except Exception as exc:
            return f"Could not read lineage for {job!r}: {exc}"

    @server.tool(
        name="bp_ancestry",
        title="Which design came from which input",
        description=(
            "The per-ID parent graph of a run. With `id`, walks one design back to the inputs "
            "it came from, step by step, showing every parent — a complex has both its designed "
            "chain and its tag, and only one of them is the binder. Without `id`, summarizes "
            "the whole graph: which step drew from which, and how each link was established. "
            "Use this for 'prove this design came from that structure', which bp_lineage cannot "
            "answer — it counts IDs, it does not connect them. An id that several steps carry "
            "unchanged is traced from the latest of them; pass `step` (a step folder) to start "
            "from another. Pass csv=True to write "
            "ancestry.csv, one row per link, for a supervisor or a referee."),
    )
    def bp_ancestry(job: str, id: str | None = None, step: str | None = None, csv: bool = False,
                    host: str | None = None, project: str | None = None,
                    repo: str | None = None) -> str:
        try:
            fs, job_dir = _resolve(job, host, project, repo)
            name = str(job_dir).rstrip("/").rsplit("/", 1)[-1]
            graph = ancestry.collect(job_dir, fs=fs)
            if id:
                found = ancestry.trace(graph, id, step=step)
                if not found:
                    return (f"{id!r} was not produced by any step of {name}. Call bp_ancestry "
                            f"without an id to see which IDs the run holds.")
                return "\n".join(ancestry.render_trace(found))
            if csv:
                out = _views_path(job_dir, "_ancestry.csv")
                rows = ancestry.to_rows(graph)
                with open(out, "w", newline="", encoding="utf-8") as handle:
                    writer = _csv.DictWriter(
                        handle, fieldnames=["child_step", "child", "parent_step", "parent", "tier"])
                    writer.writeheader()
                    writer.writerows(rows)
                return f"Wrote {out} — {len(rows)} parent link(s), one per row."
            return ancestry.summarize(graph, job=name)
        except Exception as exc:
            return f"Could not build ancestry for {job!r}: {exc}"

    @server.tool(
        name="bp_table",
        title="Read a result table",
        description=(
            "One step's CSV as rows — metrics, the dropped ids in missing.csv, a map table. "
            "Omit `table` to list what that step wrote. Truncated to `limit` rows: these run to "
            "thousands of lines, so ask for what you need rather than the whole file."),
    )
    def bp_table(job: str, step: str, table: str | None = None, limit: int = 20,
                 host: str | None = None, project: str | None = None,
                 repo: str | None = None) -> str:
        try:
            fs, job_dir = _resolve(job, host, project, repo)
            if table is None:
                names = lineage.tables(job_dir, step, fs=fs)
                return (f"{step} wrote: " + ", ".join(names)) if names else (
                    f"{step} wrote no CSV tables.")
            found = lineage.read_table(job_dir, step, table, limit=limit, fs=fs)
            if "error" in found:
                listed = ", ".join(found["available"]) or "none"
                return f"{found['error']}. Tables in that step: {listed}."
            head = (f"{found['step']}/{found['table']} — {found['shown']} of "
                    f"{found['rows']} rows\n")
            return head + "\n".join([found["header"]] + found["lines"])
        except Exception as exc:
            return f"Could not read {table or 'tables'} in {step!r}: {exc}"

    @server.tool(
        name="bp_submit",
        title="Submit a pipeline",
        description=(
            "Run a BioPipelines script on the cluster. `script` is its path relative to the "
            "repo there (e.g. my_pipelines/foo.py); pass `upload` (a local file path) to copy "
            "that file to `script` first, and `verbose=True` for ./submit -v. Returns the run directories created and the scheduler job ids. "
            "This is how a pipeline is submitted — do not run `./submit` over ssh instead, which "
            "skips the operations record written next to the run. "
            "THIS SPENDS COMPUTE — only call it when the user has asked for the run. Nothing "
            "has executed when it returns; poll bp_status afterwards."),
    )
    def bp_submit(script: str, upload: str | None = None, verbose: bool = False,
                  host: str | None = None, repo: str | None = None) -> str:
        try:
            host, ssh = _connect(host, repo)
            if not host:
                return ("bp_submit needs a cluster. Call bp_setup with your ssh alias first, "
                        "or pass host=.")
            result = submission.submit(ssh, script, upload=upload, verbose=verbose)
            return submission.summarize(result, host=host)
        except remote.RemoteError as exc:
            return f"Could not submit {script!r} on {host}: {exc}"

    @server.tool(
        name="bp_fetch",
        title="Pull results back",
        description=(
            "Copy a file or a whole step folder from the cluster to this machine — structures, "
            "tables, figures. `remote_path` is an absolute cluster path; `local` is where to "
            "put it. Directories are copied recursively without being asked. Use this when "
            "results exist only on the cluster and someone needs them locally."),
    )
    def bp_fetch(remote_path: str, local: str, overwrite: bool = False,
                 host: str | None = None, repo: str | None = None) -> str:
        try:
            host, ssh = _connect(host, repo)
            if not host:
                return "bp_fetch needs a cluster. Call bp_setup with your ssh alias first, or pass host=."
            if pathlib.Path(local).exists() and not overwrite and not pathlib.Path(local).is_dir():
                return (f"{local} already exists and would be overwritten. Pass overwrite=True "
                        f"if that is intended, or choose another `local`.")
            is_dir = ssh.is_dir(remote_path)
            code, out, err = ssh.download(remote_path, local, recursive=is_dir)
        except remote.RemoteError as exc:
            return f"Could not fetch {remote_path!r} from {host}: {exc}"
        if code != 0:
            return f"Could not fetch {remote_path!r} from {host}: {err.strip() or out.strip()}"
        kind = "folder" if is_dir else "file"
        return f"Fetched {kind} {remote_path} -> {local}"

    @server.tool(
        name="bp_resubmit",
        title="Resubmit a batch",
        description=(
            "Resume a run after a failure or a cancel by resubmitting one of its generated "
            "batch scripts. `script` is the file name inside <job>/RunTime — omit it when the "
            "run has only one, name it when it has several, because resubmitting the wrong "
            "batch spends compute on work that already succeeded; the error lists them. "
            "Dependency directives are stripped by default: they name the ORIGINAL run's job "
            "ids, which have finished or aged out, and keeping them leaves the job pending "
            "forever. Pass keep_dependencies only when the parent is still queued. "
            "THIS SPENDS COMPUTE — only call it when the user has asked to resume."),
    )
    def bp_resubmit(job: str, script: str | None = None, keep_dependencies: bool = False,
                    host: str | None = None, project: str | None = None,
                    repo: str | None = None) -> str:
        try:
            host, _fs = _connect(host, repo)
            if not host:
                return ("bp_resubmit needs a cluster: the batch scripts and the scheduler are "
                        "there. Call bp_setup first, or pass host=.")
            fs, job_dir = _resolve(job, host, project, repo)
            result = submission.resubmit(fs, job_dir, script=script,
                                         keep_dependencies=keep_dependencies)
            return submission.summarize_action(result, host=host)
        except Exception as exc:
            return f"Could not resubmit {job!r}: {exc}"

    @server.tool(
        name="bp_cancel",
        title="Cancel a run's jobs",
        description=(
            "Cancel the scheduler jobs of a run. With `job` alone it reads the ids out of that "
            "run's own operations log, so nobody has to find them first; `job_ids` overrides "
            "that. Only jobs the scheduler still holds are touched. Without confirm=True it "
            "cancels nothing and returns the ids it would cancel: show them to the user, then "
            "call again with confirm=True. THIS STOPS WORK IN PROGRESS and cannot be undone — "
            "only call it when the user has asked to cancel. A cancelled "
            "step leaves no FAILED marker, so bp_status will show it as pending, not failed."),
    )
    def bp_cancel(job: str | None = None, job_ids: list[str] | None = None,
                  confirm: bool = False, host: str | None = None, project: str | None = None,
                  repo: str | None = None) -> str:
        try:
            host, fs = _connect(host, repo)
            if not host:
                return "bp_cancel needs a cluster. Call bp_setup first, or pass host=."
            job_dir = remote.resolve_job_dir(fs, job, project=project) if job else None
            result = submission.cancel(fs, job_dir=job_dir, job_ids=job_ids, dry_run=not confirm)
            return submission.summarize_action(result, host=host)
        except Exception as exc:
            return f"Could not cancel: {exc}"

    @server.tool(
        name="bp_visualize",
        title="Show a step's results",
        description=(
            "Render one step's outputs as a page and bring it back, with any images the step "
            "produced returned inline. Use it the moment a step finishes rather than describing "
            "its numbers in prose — it runs on the login node, needs no scheduler, and works "
            "while the rest of the job is still queued. `descending`/`ascending` take "
            "TABLE.COLUMN to order by any column of any of the step's tables, and `max_items` "
            "caps how many are shown, which is how you answer 'show me the best five'. "
            "Structure views in the page are interactive and only render in a browser, so open "
            "the fetched file for those; plots and rendered images come back inline."),
    )
    def bp_visualize(job: str, step: str, descending: str | None = None,
                     ascending: str | None = None, max_items: int = 5,
                     ids: list[str] | None = None, fetch: bool = True,
                     host: str | None = None, project: str | None = None,
                     repo: str | None = None) -> list:
        try:
            host, _fs = _connect(host, repo)
            if not host:
                return ["bp_visualize needs a cluster. Call bp_setup first, or pass host=."]
            fs, job_dir = _resolve(job, host, project, repo)
            result = views.render(fs, job_dir, step, descending=descending, ascending=ascending,
                                  max_items=max_items, ids=ids)
            if not result["ok"]:
                return [views.summarize(result)]

            local_page = None
            if fetch:
                local_page = _fetch_to_scratch(fs, result["page"])

            blocks, shown, skipped = [], [], []
            for path, name in views.images(fs, job_dir, step):
                local = _fetch_to_scratch(fs, path)
                if local is None:
                    skipped.append(name)
                    continue
                data = pathlib.Path(local).read_bytes()
                if len(data) > views.MAX_IMAGE_BYTES:
                    skipped.append(name)
                    continue
                blocks.append(Image(data=data, format=name.rsplit(".", 1)[-1].lower()))
                shown.append(name)
            return [views.summarize(result, page_local=local_page, shown=shown,
                                    skipped=skipped)] + blocks
        except Exception as exc:
            return [f"Could not render {step!r} of {job!r}: {exc}"]

    @server.tool(
        name="bp_project",
        title="Project documents",
        description=(
            "Survey or scaffold a project folder's two documents, or append a dated HISTORY.md "
            "entry. action='survey' reports which documents and job folders exist; "
            "'scaffold' creates only what is missing and never overwrites; "
            "'history' appends one entry (needs `title`, optionally `body`). "
            "`project_dir` is a project name under the output root, or a full path; pass `host` "
            "for a cluster, as with the other tools. "
            "PROJECT.md holds the goal, inputs and validated protocol and carries judgment: "
            "propose changes to the scientist, never rewrite it yourself."),
    )
    def bp_project(project_dir: str, action: str = "survey",
                   title: str | None = None, body: str = "",
                   name: str | None = None, readme: bool = False,
                   host: str | None = None, repo: str | None = None) -> str:
        try:
            host, fs = _connect(host, repo)
            # A project name is what `bp_runs` prints, so accept it here too rather than making
            # the caller paste an absolute path it never showed them.
            project_dir = remote.resolve_job_dir(fs, project_dir)
        except remote.RemoteError as exc:
            return f"Could not reach {project_dir!r}: {exc}"

        if action == "survey":
            state = project_docs.survey(project_dir, fs=fs)
            if not state["exists"]:
                return (f"No project folder at {state['project_dir']}"
                        + (f" on {host}." if host else " on this machine.")
                        + " Call bp_runs to list the projects that do exist, or bp_project with "
                          "action='scaffold' to create this one.")
            docs = ", ".join(f"{k}={'yes' if v else 'NO'}" for k, v in state["documents"].items())
            jobs = ", ".join(state["jobs"]) or "none"
            recent = "; ".join(f"{e['date']} {e['title']}" for e in state["entries"][-5:]) or "none"
            return (f"{state['name']}\n  documents: {docs}\n  jobs: {jobs}\n"
                    f"  last history entries: {recent}")
        if action == "scaffold":
            result = project_docs.scaffold(project_dir, name=name, readme=readme, fs=fs)
            return "\n".join(f"{k}: {v}" for k, v in sorted(result.items()))
        if action == "history":
            if not title:
                return "action='history' needs a `title`."
            try:
                entry = project_docs.append_history(project_dir, title, body, fs=fs)
            except (OSError, remote.RemoteError) as exc:
                return f"Nothing was appended: {exc}"
            return "Appended to HISTORY.md:\n\n" + entry
        return f"Unknown action {action!r}. Use 'survey', 'scaffold' or 'history'."

    @server.tool(
        name="bp_provenance",
        title="What produced a run",
        description=(
            "The recorded identity of a run: BioPipelines version and commit, site variant, "
            "scheduler, and per step the tool version, environments, container image and the "
            "parameters as they resolved — including the defaults nobody passed, which are "
            "usually what decided the output. Pass `against` (a second job) to get the "
            "differences between two runs instead: tool versions, changed parameters and "
            "environment digests, named one by one. Use this to answer 'what produced this "
            "structure' and 'why did the same pipeline give different numbers' — it reads a "
            "file the run wrote, so it costs nothing and needs no scheduler."),
    )
    def bp_provenance(job: str, against: str | None = None, host: str | None = None,
                      project: str | None = None, repo: str | None = None) -> str:
        try:
            fs, job_dir = _resolve(job, host, project, repo)
            record = manifest.read(job_dir, fs=fs)
            if record is None:
                return (f"No {manifest.FILENAME} in {job_dir}. That run predates the manifest, so "
                        f"what produced it was never recorded.")
            if against is None:
                return _render_manifest(record, job_dir)
            _fs2, other_dir = _resolve(against, host, project, repo)
            other = manifest.read(other_dir, fs=fs)
            if other is None:
                return f"No {manifest.FILENAME} in {other_dir}, so there is nothing to compare."
            lines = manifest.compare(record, other)
            if not lines:
                if manifest.environments_were_compared(record, other):
                    return (f"{job} and {against} are the same run in every recorded respect: "
                            f"same versions, same parameters, same environments.")
                return (f"{job} and {against} agree on every recorded version and parameter. "
                        f"Their environments could NOT be compared — at least one recorded no "
                        f"digest, so this does not say the software was the same.")
            return (f"{job} -> {against} differs in {len(lines)} recorded respects:\n"
                    + "\n".join(f"  - {line}" for line in lines))
        except Exception as exc:
            return f"Could not read provenance for {job!r}: {exc}"

    def _render_manifest(record, job_dir):
        lines = [f"{record.get('project')}/{record.get('job')} — recorded {record.get('created')}",
                 f"  BioPipelines {record.get('biopipelines')} @ {record.get('commit')}"
                 f", variant {record.get('variant')}, {record.get('scheduler')}"
                 f", python {record.get('python')}",
                 f"  identity: {record.get('hash')}"]
        environments = record.get("environments_resolved") or {}
        if environments:
            lines.append("  environments: " + ", ".join(
                f"{name} @ {digest[:12]}" for name, digest in sorted(environments.items())))
        else:
            lines.append("  environments: not recorded (run predates the environment capture, or "
                         "the node had no sha256sum)")
        for tool in record.get("tools") or []:
            resolved = (tool.get("parameters") or {}).get("resolved") or {}
            passed = (tool.get("parameters") or {}).get("passed") or {}
            shown = ", ".join(f"{k}={v!r}" for k, v in list(resolved.items())[:8])
            lines.append(f"\n  {tool.get('order')}. {tool.get('tool')} v{tool.get('tool_version')}"
                         f"  env={','.join(tool.get('environments') or []) or '-'}"
                         + (f"  image={tool['container_image']}" if tool.get("container_image") else ""))
            lines.append(f"     passed: {', '.join(sorted(passed)) or 'nothing'}")
            lines.append(f"     resolved: {shown}"
                         + (f" … {len(resolved) - 8} more" if len(resolved) > 8 else ""))
        return "\n".join(lines)

    @server.tool(
        name="bp_reproduce",
        title="Plan a rerun of a recorded campaign",
        description=(
            "What it would take to run a past campaign again — here or on another cluster. "
            "Reads the run's manifest, finds the pipeline script the framework preserved in "
            "RunTime/, asks the target host what it would bring, and names every difference "
            "before anything is submitted: BioPipelines version and commit, site variant, "
            "environment manager, scheduler. Spends no compute and submits nothing; it ends with the exact "
            "bp_submit call to make, which is yours to decide on. Pass `target` to plan the "
            "rerun against a different host than the one that holds the run."),
    )
    def bp_reproduce(job: str, target: str | None = None, host: str | None = None,
                     project: str | None = None, repo: str | None = None) -> str:
        try:
            fs, job_dir = _resolve(job, host, project, repo)
            # A repo given for the run's host says nothing about where the target keeps its checkout.
            where, ssh = _connect(target or host, repo if not target or target == host else None)
            found = reproduce.plan(job_dir, fs, ssh=ssh if where else None,
                                   repo=getattr(ssh, "repo", "") or "")
            return reproduce.summarize(found, host=where or "")
        except Exception as exc:
            return f"Could not plan a rerun of {job!r}: {exc}"

    return server


def main():
    try:
        server = build_server()
    except ModuleNotFoundError as exc:
        print(exc, file=sys.stderr)
        return 1
    server.run(transport="stdio")
    return 0


if __name__ == "__main__":
    sys.exit(main())
