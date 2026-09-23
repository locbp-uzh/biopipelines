"""Run a recorded campaign again, and say up front what will not be the same.

Reproduction here is not reconstruction. The manifest records each tool's resolved parameters, but a parameter value that was a DataStream is recorded as a name, so the wiring — which step feeds which — is not recoverable from it. What is recoverable is the pipeline script itself: `save()` copies the calling `.py` into `RunTime/`, so the authored campaign survives next to its own record.

So `plan()` pairs the two. It reads the manifest to learn what the original run was, finds the preserved script, and asks the target host what it would bring to a rerun: which BioPipelines version and commit, which site variant, environment manager and scheduler. Environment digests are not probed: they exist only once a run has executed there. Then it names every difference before anything is submitted, because the useful moment to learn that a tool version moved is before the compute is spent, not after the numbers disagree.

A plan with differences is not a refusal. Reproducing a campaign on a second cluster differs by construction — different host, often a different variant. The point is that the differences are stated rather than discovered.
"""

import json
import posixpath
from typing import Any, Dict, List, Optional

from . import manifest

PROBE = (
    "import json,os,hashlib,sys\n"
    "sys.path.insert(0, os.path.expanduser(os.path.expandvars({repo!r})))\n"
    "out={{}}\n"
    "try:\n"
    "    from biopipelines import __version__ as v\n"
    "    out['biopipelines']=v\n"
    "except Exception as e:\n"
    "    out['biopipelines']=None\n"
    "try:\n"
    "    from biopipelines.config_manager import ConfigManager\n"
    "    cm=ConfigManager()\n"
    "    out['variant']=cm.get_variant()\n"
    "    out['env_manager']=cm.get_env_manager()\n"
    "    out['scheduler']=cm.get_scheduler()\n"
    "except Exception:\n"
    "    pass\n"
    "try:\n"
    "    from biopipelines import run_log\n"
    "    out['commit']=run_log._git_commit()\n"
    "except Exception:\n"
    "    pass\n"
    "print(json.dumps(out))\n"
)


def preserved_script(job_dir: str, fs) -> Optional[str]:
    """The authored pipeline `.py` that `save()` copied into `RunTime/`, if it is still there."""
    runtime = fs.join(job_dir, "RunTime")
    try:
        entries = fs.listdir(runtime)
    except Exception:
        return None
    scripts = [name for name, is_dir in entries if not is_dir and name.endswith(".py")]
    if not scripts:
        return None
    return fs.join(runtime, sorted(scripts)[0])


def probe(ssh, repo: str = "") -> Dict[str, Any]:
    """What the target host would bring to a rerun. Returns `{}` when it cannot be asked.

    Goes through `env_prefix` like everything else that calls the framework remotely: a bare non-interactive ssh has the system PATH only, and its `python` imports nothing.
    """
    from shlex import quote

    from .remote import interpreter

    code = PROBE.format(repo=repo or getattr(ssh, "repo", ""))
    try:
        exit_code, out, _err = ssh.run(
            f"{ssh.env_prefix}{interpreter(ssh.python)} -c {quote(code)}")
    except Exception:
        return {}
    if exit_code != 0:
        return {}
    for line in reversed(str(out).strip().splitlines()):
        try:
            parsed = json.loads(line)
        except Exception:
            continue
        if not isinstance(parsed, dict):
            continue
        # A target that could not import biopipelines reported nothing, and nothing must not read as a match.
        return parsed if parsed.get("biopipelines") else {}
    return {}


def differences(record: Dict[str, Any], target: Dict[str, Any]) -> List[str]:
    """What the target would change, field by field, worst first.

    Only fields the target actually reported are compared. A probe that could not answer says nothing rather than inventing a mismatch, which would be the same output as a genuine one.
    """
    lines = []
    for key, label in (("biopipelines", "BioPipelines version"),
                       ("commit", "commit"),
                       ("variant", "config variant"),
                       ("env_manager", "environment manager"),
                       ("scheduler", "scheduler")):
        if key not in target or target.get(key) is None:
            continue
        mine, theirs = record.get(key), target.get(key)
        # Records before 1.5.1 hold a short sha; a prefix of the full one is the same commit.
        if key == "commit" and mine and (str(theirs).startswith(str(mine))
                                         or str(mine).startswith(str(theirs))):
            continue
        if mine != theirs:
            lines.append(f"{label}: recorded {record.get(key)!r}, target has {target.get(key)!r}")
    return lines


def plan(job_dir: str, fs, ssh=None, repo: str = "") -> Dict[str, Any]:
    """Everything needed to decide whether to rerun, and nothing that spends compute."""
    record = manifest.read(job_dir, fs=fs)
    if record is None:
        return {"error": f"no {manifest.FILENAME} in {job_dir}; that run predates the manifest, "
                         f"so what produced it was never recorded"}
    script = preserved_script(job_dir, fs)
    target = probe(ssh, repo=repo) if ssh is not None else {}
    return {
        "job_dir": job_dir,
        "manifest": record,
        "script": script,
        "target": target,
        "differences": differences(record, target) if target else [],
        "environments": record.get("environments_resolved") or {},
    }


def summarize(found: Dict[str, Any], host: str = "") -> str:
    """The plan as an agent should read it: what ran, what would change, and the one command left."""
    if "error" in found:
        return found["error"]
    record = found["manifest"]
    tools = record.get("tools") or []
    where = f" on {host}" if host else ""
    lines = [
        f"{record.get('project')}/{record.get('job')} — {len(tools)} steps, "
        f"recorded {record.get('created')}",
        f"  BioPipelines {record.get('biopipelines')} @ {record.get('commit')}, "
        f"variant {record.get('variant')}, {record.get('scheduler')}",
    ]
    if tools:
        lines.append("  steps: " + " -> ".join(str(t.get("tool")) for t in tools))

    environments = found.get("environments") or {}
    if environments:
        lines.append("  environments: " + ", ".join(
            f"{name} @ {digest[:12]}" for name, digest in sorted(environments.items())))
    else:
        lines.append("  environments: not recorded — this run started before the environment "
                     "capture, or the node had no sha256sum")

    script = found.get("script")
    if script:
        lines.append(f"  script: {script}")
    else:
        lines.append("  script: NOT PRESERVED in RunTime/ — there is nothing to rerun from; "
                     "the campaign has to be reauthored")

    problems = found.get("differences") or []
    if not found.get("target"):
        lines.append(f"\nTarget{where} was not probed, so nothing is known about what a rerun "
                     f"would change.")
    elif problems:
        lines.append(f"\nA rerun{where} would differ:")
        lines += [f"  - {line}" for line in problems]
        lines.append("\nThat may be fine — a second cluster differs by construction. It is stated "
                     "so it is a decision rather than a surprise.")
    else:
        lines.append(f"\nThe target{where} matches the record on every field probed.")

    if script:
        # The preserved copy, never a same-named repo file that may have been edited since.
        name = script.replace("\\", "/").rsplit("/", 1)[-1]
        lines.append(f"\nTo rerun: bp_submit with script={script!r} on the host that holds it; on "
                     f"another host, bp_fetch it, then bp_submit with upload=<the local copy> and "
                     f"script=\"my_pipelines/{name}\". THIS SPENDS COMPUTE.")
    return "\n".join(lines)
