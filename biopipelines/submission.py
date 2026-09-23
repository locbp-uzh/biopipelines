"""Submit a pipeline to the cluster and record what was submitted, where it landed.

This is the first tool that *changes* something, so it is the first that writes to a run's
operations log. Read-only queries deliberately do not — recording every status poll is what
made `llm/log.sh` useless.

The ordering problem `RunLog`'s buffer exists for shows up here: the submission parameters are
known before `./submit` runs, but the job directory it will create is not known until the
output is parsed. Records are held and flushed once the directory is known, so the header and
the submission both land in the right run with the time they actually happened.

`./submit` prints what is needed to follow a job, in a stable enough shape to parse:
`Runtime directory: <path>/RunTime` per pipeline, and `Batch N submitted: Job ID <jid>` per
batch. Its failure modes are explicit too, so a failed submission is reported rather than
silently producing no job.
"""

import re

from biopipelines import run_log
from biopipelines.remote import RemoteError, quote, validate_scp_path

_RUNTIME = re.compile(r"Runtime directory:\s*(\S+)")
_JOB_ID = re.compile(r"Batch\s+(\S+)\s+submitted:\s*Job ID\s+(\S+)")
_SINGLE_JOB = re.compile(r"submitted:\s*Job ID\s+(\S+)")
_GEN_FAILED = re.compile(r"Pipeline generation failed with exit code\s+(\d+)")
_SUB_FAILED = re.compile(r"Submission failed \(([^)]*)\):\s*(.*)")


def parse(output):
    """What `./submit` said: job directories, job ids, and any failure it announced."""
    runtimes = _RUNTIME.findall(output)
    # `<Job>_NNN/RunTime` -> `<Job>_NNN`, the directory every other tool addresses.
    jobs = [r.rstrip("/").rsplit("/", 1)[0] for r in runtimes if r.rstrip("/").endswith("RunTime")]

    ids = [jid for _batch, jid in _JOB_ID.findall(output)]
    if not ids:
        ids = _SINGLE_JOB.findall(output)

    failures = []
    generation = _GEN_FAILED.search(output)
    if generation:
        failures.append(f"pipeline generation failed (exit {generation.group(1)})")
    for scheduler, detail in _SUB_FAILED.findall(output):
        failures.append(f"submission failed ({scheduler}): {detail.strip()}")
    if "No runtime directories found" in output:
        failures.append("no runtime directories in the pipeline output — nothing was generated")

    return {"jobs": jobs, "job_ids": ids, "failures": failures,
            "ok": bool(jobs) and not failures}


# `./submit` takes no other option; anything else here would be raw shell on the cluster.
SUBMIT_FLAGS = ("-v", "--verbose")


def _submit_flags(extra):
    tokens = str(extra or "").split()
    unknown = [t for t in tokens if t not in SUBMIT_FLAGS]
    if unknown:
        raise ValueError(f"./submit accepts only {', '.join(SUBMIT_FLAGS)}; refusing {unknown}")
    return " ".join(tokens)


def submit(ssh, script, upload=None, extra="", verbose=False):
    """Run `./submit <script>` on the cluster, recording it against the run it creates.

    `upload` is a local file to copy to `<repo>/<script>` first — the usual case, since the
    pipeline is authored on the machine the agent runs on.
    """
    try:
        flags = _submit_flags(("-v " if verbose else "") + str(extra or ""))
        validate_scp_path(script)
        if upload and str(script).startswith(("/", "~", "$")):
            raise ValueError("with upload=, script is the destination inside the repo and must be relative")
    except (ValueError, RemoteError) as exc:
        return {"ok": False, "jobs": [], "job_ids": [], "failures": [str(exc)], "output": ""}
    # The run lives on the cluster, so the record has to be written there too.
    log = run_log.RunLog(fs=ssh)
    log.header(submitted_from=True, pipeline=script, host=ssh.host, repo=ssh.repo)

    if upload:
        target = f"{ssh.repo}/{script}"
        code, out, err = ssh.upload(upload, target)
        log.command(["scp", str(upload), f"{ssh.host}:{target}"], exit_code=code,
                    output=(err or out)[:500] or None)
        if code != 0:
            return {"ok": False, "jobs": [], "job_ids": [],
                    "failures": [f"could not copy {upload} to {ssh.host}:{target}: "
                                 f"{err.strip() or out.strip() or f'exit {code}'}"],
                    "output": ""}

    prefix = getattr(ssh, "env_prefix", "")
    command = (f"cd {quote(ssh.repo)} && {prefix}./submit {flags} "
               f"{quote(script)}").replace("  ", " ")
    # Generation runs the pipeline script, which can take minutes before anything is submitted.
    try:
        code, out, err = ssh.run(command, timeout=max(ssh.timeout, 600))
    except RemoteError as exc:
        # The remote ./submit may still be queueing jobs; a retry would spend the compute twice.
        return {"ok": False, "jobs": [], "job_ids": [], "outcome_unknown": True, "output": "",
                "failures": [f"{exc}. The outcome is UNKNOWN: ./submit may still be queueing "
                             f"jobs. Check bp_runs and squeue before submitting again."]}
    combined = out + ("\n" + err if err.strip() else "")

    result = parse(combined)
    result["output"] = combined
    if code != 0 and not result["failures"]:
        result["failures"].append(f"./submit exited {code}")
        result["ok"] = False

    log.record("submitted", script=script, command=command, exit_code=code,
               job_ids=result["job_ids"] or None,
               failures=result["failures"] or None)

    # The job directory only becomes known here, which is what the buffer is for. One run is the
    # normal case; a script generating several gets the same record written into each.
    buffered = log.pending
    written = 0
    for job_dir in result["jobs"]:
        for action, fields in buffered:
            written += run_log.record(job_dir, action, fs=ssh, **fields) is not None
    result["recorded"] = written
    result["unrecorded"] = len(buffered) * len(result["jobs"]) - written
    return result


def summarize(result, host=None):
    """What an agent should read back after submitting."""
    where = f" on {host}" if host else ""
    if not result["ok"]:
        lines = [f"Submission failed{where}."]
        lines += [f"  - {f}" for f in result["failures"]] or ["  - no job directory was created"]
        tail = "\n".join(result.get("output", "").splitlines()[-15:])
        if tail:
            lines.append("\nLast lines of ./submit:\n" + tail)
        return "\n".join(lines)

    lines = [f"Submitted{where}: {len(result['jobs'])} run(s), "
             f"{len(result['job_ids'])} scheduler job(s)."]
    if result.get("unrecorded"):
        lines.append(f"  WARNING: {result['unrecorded']} operations-log record(s) could not be "
                     f"written into the run directory.")
    for job in result["jobs"]:
        lines.append(f"  {job}")
    if result["job_ids"]:
        lines.append("  job ids: " + ", ".join(result["job_ids"]))
    lines.append("\nNothing has run yet — poll with bp_status once the scheduler starts it.")
    return "\n".join(lines)


# --- resuming and cancelling -----------------------------------------------------------------
#
# A campaign that fails at step 4 of 9 is the ordinary case, and until now recovering from it
# meant the agent handing the user a shell command — which is the moment the premise of driving
# BioPipelines through tools visibly breaks. Both of these change scheduler state, so both write
# to the run's operations log the way `submit` does.

_JOB_SCRIPT = re.compile(r"^(slurm|lsf|pbs)(_batch\d+)?\.sh$")


def job_scripts(ssh, job_dir):
    """The generated batch scripts of a run, in batch order.

    Every `Resources()` call opens a batch and each batch is submitted as its own job, so a run
    has `slurm_batch<N>.sh` per batch — or a single `slurm.sh` when it ended up with one. Which
    exists is a property of the pipeline, not something to assume.
    """
    runtime = ssh.join(job_dir, "RunTime")
    if not ssh.is_dir(runtime):
        return []
    found = [name for name, is_dir in ssh.listdir(runtime)
             if not is_dir and _JOB_SCRIPT.match(name)]
    return sorted(found)


def resubmit(ssh, job_dir, script=None, keep_dependencies=False):
    """Resubmit one of a run's batch scripts. Returns what happened, for `summarize_action`.

    `script` is the file name inside `<job>/RunTime`. With one script it is optional; with
    several, naming one is required rather than guessed — resubmitting the wrong batch spends
    compute on work that already succeeded.
    """
    available = job_scripts(ssh, job_dir)
    if not available:
        return {"ok": False, "action": "resubmit", "available": [],
                "failures": [f"no batch script under {job_dir}/RunTime — nothing to resubmit"]}
    if script is None:
        if len(available) > 1:
            return {"ok": False, "action": "resubmit", "available": available,
                    "failures": [f"{len(available)} batch scripts in this run; name the one to "
                                 f"resubmit"]}
        script = available[0]
    if script not in available:
        return {"ok": False, "action": "resubmit", "available": available,
                "failures": [f"no script {script!r} in {job_dir}/RunTime"]}

    target = ssh.join(job_dir, "RunTime", script)
    # Dependency directives name the ORIGINAL run's job ids, which have since finished or aged
    # out; keeping them leaves the job pending on something that can never be satisfied.
    flag = "--keep-dependencies " if keep_dependencies else ""
    command = (f"cd {quote(ssh.repo)} && {getattr(ssh, 'env_prefix', '')}"
               f"./resubmit {flag}{quote(target)}")
    try:
        code, out, err = ssh.run(command, timeout=max(ssh.timeout, 300))
    except RemoteError as exc:
        # As with submit: ./resubmit may have queued the job before the connection dropped.
        return {"ok": False, "action": "resubmit", "script": script, "job_ids": [],
                "available": available, "output": "", "outcome_unknown": True,
                "failures": [f"{exc}. The outcome is UNKNOWN: check squeue before resubmitting."]}
    combined = out + ("\n" + err if err.strip() else "")

    ids = [jid for _batch, jid in _JOB_ID.findall(combined)] or _SINGLE_JOB.findall(combined)
    failures = [] if code == 0 else [f"./resubmit exited {code}"]
    run_log.record(job_dir, "resubmitted", fs=ssh, script=script, command=command,
                   exit_code=code, job_ids=ids or None,
                   keep_dependencies=keep_dependencies or None, failures=failures or None)
    return {"ok": code == 0, "action": "resubmit", "script": script, "job_ids": ids,
            "failures": failures, "output": combined, "available": available}


def recorded_job_ids(ssh, job_dir):
    """Scheduler ids this run's operations log has seen, newest last, de-duplicated."""
    ids = []
    for entry in run_log.read(job_dir, fs=ssh):
        for jid in entry.get("job_ids") or []:
            if jid not in ids:
                ids.append(jid)
    return ids


def live_job_ids(ssh, ids):
    """The subset of `ids` the scheduler still holds."""
    if not ids:
        return []
    code, out, err = ssh.run("squeue -h -o %i -j " + quote(",".join(str(i) for i in ids)))
    if code != 0 and not out.strip():
        # squeue exits non-zero once every id has left the queue; any other failure is not an answer.
        if "invalid job id" in err.lower():
            return []
        raise RemoteError(f"squeue failed: {err.strip() or f'exit {code}'}")
    held = set()
    for token in (line.strip() for line in out.splitlines() if line.strip()):
        # An array task (123_4, 123_[1-5]) or het component (123+0) is held under its own id and its job's.
        held.update({token, re.split(r"[_+]", token, maxsplit=1)[0]})
    return [i for i in ids if str(i) in held]


def cancel(ssh, job_dir=None, job_ids=None, dry_run=False):
    """Cancel scheduler jobs, by id or by reading the run's own log for them.

    Making the caller find the ids first is what pushes an agent back to a shell, so the ids of
    a named run are read from the record `submit` and `resubmit` already wrote. Only ids the
    scheduler still holds are cancelled, and `dry_run` names them without cancelling anything.
    """
    ids = list(job_ids or []) or (recorded_job_ids(ssh, job_dir) if job_dir else [])
    if not ids:
        where = f" in {job_dir}'s operations log" if job_dir else ""
        return {"ok": False, "action": "cancel", "job_ids": [],
                "failures": [f"no scheduler job ids{where} — pass job_ids explicitly"]}

    live = live_job_ids(ssh, ids)
    gone = [i for i in ids if i not in live]
    if dry_run:
        return {"ok": True, "action": "cancel", "dry_run": True, "job_ids": live,
                "finished": gone, "failures": []}
    if not live:
        return {"ok": False, "action": "cancel", "job_ids": [], "finished": gone,
                "failures": ["none of these jobs is still queued or running: "
                             + ", ".join(map(str, ids))]}
    ids = live

    command = "scancel " + " ".join(quote(str(i)) for i in ids)
    code, out, err = ssh.run(command)
    combined = (out + ("\n" + err if err.strip() else "")).strip()
    failures = [] if code == 0 else [f"scancel exited {code}: {combined or 'no output'}"]
    if job_dir:
        run_log.record(job_dir, "cancelled", fs=ssh, command=command, exit_code=code,
                       job_ids=ids, failures=failures or None)
    return {"ok": code == 0, "action": "cancel", "job_ids": ids, "finished": gone,
            "failures": failures, "output": combined}


def summarize_action(result, host=None):
    """One screen for a resubmit or a cancel, in the shape `summarize` uses for a submission."""
    where = f" on {host}" if host else ""
    verb = {"resubmit": "Resubmission", "cancel": "Cancellation"}[result["action"]]
    if not result["ok"]:
        lines = [f"{verb} failed{where}."]
        lines += [f"  - {f}" for f in result["failures"]]
        if result.get("available"):
            lines.append("  batch scripts in this run: " + ", ".join(result["available"]))
        tail = "\n".join(result.get("output", "").splitlines()[-10:])
        if tail:
            lines.append("\nLast lines:\n" + tail)
        return "\n".join(lines)

    if result["action"] == "cancel" and result.get("dry_run"):
        lines = [f"Would cancel{where}: " + (", ".join(map(str, result["job_ids"])) or "nothing — no job is still queued or running")]
        if result.get("finished"):
            lines.append("  already finished, left alone: " + ", ".join(map(str, result["finished"])))
        lines.append("\nNothing was cancelled. Call again with confirm=True to cancel these.")
        return "\n".join(lines)
    if result["action"] == "cancel":
        return (f"Cancelled{where}: " + ", ".join(result["job_ids"])
                + "\n\nThe scheduler may take a moment to release the allocation; bp_status "
                  "reads markers on disk, so a cancelled step stays 'pending' rather than "
                  "becoming 'failed'.")
    lines = [f"Resubmitted{where}: {result['script']}"]
    if result["job_ids"]:
        lines.append("  job ids: " + ", ".join(result["job_ids"]))
    lines.append("\nNothing has run yet — poll with bp_status once the scheduler starts it.")
    return "\n".join(lines)
