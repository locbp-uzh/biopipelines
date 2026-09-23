"""What produced a run, recorded on every run rather than on request.

`Pipeline(debug=True)` already snapshots the environment into `_debug_capture/`, and that snapshot is richer than this one: full env exports, `nvidia-smi`, the active config file. But it is opt-in, and a record you have to remember to switch on is not a record — by the time a result is surprising, the run that produced it is over. This writes the half that costs nothing: versions, resolved parameters, environment and container identity, hashed so two runs can be compared by a single value.

The hash deliberately covers only what would make a rerun differ. Timestamps, hostnames and the job's own name are recorded but excluded, so the same pipeline run twice hashes the same and a changed default cannot hide behind a changed date.
"""

import dataclasses
import hashlib
import inspect
import json
import os
from typing import Any, Dict, List, Optional

SCHEMA = 1
FILENAME = "manifest.json"

# Recorded for the reader, excluded from the hash: these differ between two runs that are otherwise identical.
VOLATILE = ("schema", "hash", "created", "host", "project", "job", "job_dir")


def path_for(job_dir: str, fs=None) -> str:
    join = fs.join if fs is not None else os.path.join
    return join(job_dir, FILENAME)


def _jsonable(value: Any) -> Any:
    """A JSON-safe, run-stable rendering of a parameter value.

    Stability is the point: a default `repr()` embeds a memory address, which would give every run a different hash and make comparison worthless. Anything that is not a primitive or a container collapses to a name or a type. A value that cannot be rendered collapses to its type rather than losing the whole manifest.
    """
    try:
        return _render(value)
    except Exception:
        return "<" + type(value).__name__ + ">"


def _path_tail(path: Any, parts: int = 3) -> str:
    """The step-relative end of an output path, so a moved checkout or output root hashes the same."""
    pieces = str(path).replace("\\", "/").rstrip("/").split("/")
    return "/".join(pieces[-parts:])


def _render(value: Any) -> Any:
    from .data_containers import TableInfo
    from .biopipelines_io import TableReference
    from .combinatorics import Grouped

    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, (set, frozenset)):
        # Iteration order is randomized per process; unordered in, ordered out.
        return sorted((_jsonable(v) for v in value), key=str)
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in sorted(value.items(), key=lambda kv: str(kv[0]))}
    # Before any duck-typed probe: a TableInfo answers every public attribute with a column reference.
    if isinstance(value, TableInfo):
        return "<table " + str(value.info.name) + " " + _path_tail(value.info.path) + ">"
    if isinstance(value, TableReference):
        return "<column " + str(value.column) + " of " + _path_tail(value.path) + ">"
    if isinstance(value, Grouped):
        return {"combinatorics": "Grouped", "source": _jsonable(value.source),
                "groups": _jsonable(value.groups)}
    tool_name = getattr(type(value), "TOOL_NAME", None)
    if tool_name:
        folder = getattr(value, "output_folder", "") or ""
        suffix = ":" + os.path.basename(folder) if folder else ""
        return "<tool " + str(tool_name) + suffix + ">"
    if hasattr(value, "ids") and hasattr(value, "format"):
        # The producer and the declared ids, not the stream's kind. Half the pipelines in the
        # repo have a `structures` stream, so the kind alone puts the wiring outside the hash.
        name = getattr(value, "name", "") or type(value).__name__
        return "<stream " + str(name) + _stream_origin(value) + ">"
    if hasattr(value, "streams") and hasattr(value, "tables"):
        # A StandardizedOutput is the commonest wiring value of all.
        inner = getattr(value, "streams", None)
        names = sorted(getattr(inner, "__dict__", {}) or {}) if inner is not None else []
        return "<output " + ",".join(str(n) for n in names) + _stream_origin(value) + ">"
    if hasattr(value, "sources"):
        return {"combinatorics": type(value).__name__,
                "sources": [_jsonable(s) for s in (getattr(value, "sources", None) or [])]}
    if hasattr(value, "__fspath__"):
        return str(value)
    # Last, so a DataStream (also a dataclass) keeps its path-free stream form above.
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        if _render.depth > 8:
            return "<" + type(value).__name__ + ">"
        _render.depth += 1
        try:
            # `__class__`, not `type`: Panda's Operation has a field named `type`.
            return {"__class__": type(value).__name__,
                    **{f.name: _jsonable(getattr(value, f.name)) for f in dataclasses.fields(value)}}
        finally:
            _render.depth -= 1
    return "<" + type(value).__name__ + ">"


_render.depth = 0


def _stream_origin(value: Any) -> str:
    """Which step produced this, and how wide it is — the part that identifies a connection."""
    producer = getattr(value, "_producer", None) or getattr(value, "producer", None)
    folder = getattr(producer, "output_folder", "") if producer is not None else ""
    where = os.path.basename(str(folder)) if folder else ""
    ids = getattr(value, "ids", None)
    try:
        width = len(ids) if ids is not None else None
    except TypeError:
        width = None
    parts = []
    if where:
        parts.append(" from " + where)
    if width is not None:
        parts.append(" n=%d" % width)
    return "".join(parts)


def _declared_parameters(cls: type) -> Dict[str, Any]:
    """Every named constructor parameter along the MRO, with its default, most-derived winning."""
    declared: Dict[str, Any] = {}
    for base in reversed(cls.__mro__):
        init = base.__dict__.get("__init__")
        if init is None:
            continue
        try:
            parameters = inspect.signature(init).parameters
        except (TypeError, ValueError):
            continue
        for name, parameter in parameters.items():
            if name == "self" or parameter.kind in (parameter.VAR_KEYWORD, parameter.VAR_POSITIONAL):
                continue
            declared[name] = None if parameter.default is inspect.Parameter.empty else parameter.default
    return declared


def parameters_of(tool) -> Dict[str, Any]:
    """What the user passed and what the tool ended up with.

    Both, because neither alone answers the question. `passed` is the three kwargs someone wrote; `resolved` is those plus the twenty defaults that actually determined the output. Resolution reads the attribute the wrapper bound, falling back to the signature default when a wrapper stores a parameter under another name.
    """
    # `params` holds only what BaseConfig received, which is the leftovers after the subclass bound
    # its named parameters -- empty for every typed argument. `_constructor_kwargs` is the call.
    written = getattr(tool, "_constructor_kwargs", None)
    if written is None:
        written = getattr(tool, "params", None) or {}
    passed = {k: _jsonable(v) for k, v in sorted(written.items())}
    resolved = {}
    for name, default in sorted(_declared_parameters(type(tool)).items()):
        # The bound attribute is what ran, including any normalization or override the
        # constructor applied. `passed` is the fallback only where the wrapper consumed the
        # parameter without binding it, which is the case that made this record empty before.
        missing = object()
        bound = getattr(tool, name, missing)
        if bound is not missing:
            resolved[name] = _jsonable(bound)
        elif name in passed:
            resolved[name] = passed[name]
        else:
            resolved[name] = _jsonable(default)
    record = {"passed": passed, "resolved": resolved}
    extra = getattr(tool, "extra_args", None)
    if extra:
        record["forwarded"] = {k: _jsonable(v) for k, v in sorted(extra.items())}
    inputs = _input_files(written)
    if inputs:
        record["input_files"] = inputs
    return record


def _input_files(written: Dict[str, Any]) -> Dict[str, str]:
    """`{parameter: sha256}` for every passed string that names an existing file: the path alone does not say what was in it."""
    found = {}
    # save() runs on a login node, so hashing stops at a budget; past it a file is identified by its size.
    budget = {"bytes": 1024 ** 3, "files": 200}
    for name, value in sorted(written.items()):
        candidates = value if isinstance(value, (list, tuple)) else [value]
        digests = []
        for item in candidates:
            if not isinstance(item, str) or len(item) > 4096:
                continue
            try:
                if not os.path.isfile(item):
                    continue
                size = os.path.getsize(item)
                if size > budget["bytes"] or budget["files"] <= 0:
                    digests.append(os.path.basename(item) + ":size=" + str(size))
                    continue
                budget["bytes"] -= size
                budget["files"] -= 1
                sha = hashlib.sha256()
                with open(item, "rb") as handle:
                    for block in iter(lambda: handle.read(1 << 20), b""):
                        sha.update(block)
                digests.append(os.path.basename(item) + ":" + sha.hexdigest())
            except OSError:
                continue
        if digests:
            found[name] = digests[0] if len(digests) == 1 and not isinstance(value, (list, tuple)) else digests
    return found


def _tool_entry(tool, order: int, folders: Dict[str, str]) -> Dict[str, Any]:
    entry = {
        "order": order,
        "tool": getattr(tool, "TOOL_NAME", type(tool).__name__),
        "class": type(tool).__name__,
        "tool_version": getattr(tool, "TOOL_VERSION", None),
        "environments": list(getattr(tool, "environments", []) or []),
        "container_image": "",
        "parameters": parameters_of(tool),
    }
    try:
        if tool.uses_container():
            entry["container_image"] = folders.get("container:" + str(tool.TOOL_NAME), "")
    except Exception:
        pass
    return entry


def _site(config_manager) -> Dict[str, Any]:
    out = {}
    for key, getter in (("variant", "get_variant"), ("env_manager", "get_env_manager"),
                        ("scheduler", "get_scheduler"), ("container_executor", "get_container_executor")):
        try:
            out[key] = getattr(config_manager, getter)()
        except Exception:
            out[key] = None
    return out


def build(pipeline) -> Dict[str, Any]:
    """The manifest for a pipeline, without touching the filesystem."""
    from . import run_log
    from .config_manager import ConfigManager

    header = run_log._header_fields()
    header.pop("config_variant", None)
    manifest: Dict[str, Any] = {"schema": SCHEMA, "created": run_log._now()}
    manifest.update(header)
    manifest.update(_site(ConfigManager()))
    manifest["project"] = getattr(pipeline, "project", None)
    manifest["job"] = getattr(pipeline, "job", None)
    folders = getattr(pipeline, "folders", {}) or {}
    manifest["tools"] = [_tool_entry(tool, order, folders)
                         for order, tool in enumerate(getattr(pipeline, "tools", []) or [], 1)]
    manifest["hash"] = digest(manifest)
    return manifest


def digest(manifest: Dict[str, Any]) -> str:
    """A hash of everything that would make a rerun differ."""
    material = {k: v for k, v in manifest.items() if k not in VOLATILE}
    blob = json.dumps(material, sort_keys=True, separators=(",", ":"), default=str)
    return "sha256:" + hashlib.sha256(blob.encode("utf-8")).hexdigest()


def write(pipeline, job_dir: Optional[str] = None, fs=None) -> Optional[str]:
    """Write the manifest beside the run. Never raises: a provenance record must not break a pipeline.

    The reason for a failure is kept on `write.last_error`, so the caller's warning can name it.
    """
    write.last_error = None
    try:
        manifest = build(pipeline)
        target = job_dir or (getattr(pipeline, "folders", {}) or {}).get("output")
        if not target:
            return None
        destination = path_for(target, fs=fs)
        text = json.dumps(manifest, indent=2, default=str)
        if fs is not None:
            fs.write_text(destination, text)
        else:
            with open(destination, "w", encoding="utf-8") as handle:
                handle.write(text)
        return destination
    except Exception as error:
        write.last_error = f"{type(error).__name__}: {error}"
        return None


write.last_error = None


def environments(job_dir: str, fs=None) -> Dict[str, str]:
    """`{environment name: sha256}` as the executing node recorded it.

    Written by the run itself, so it is absent until the job has started and absent for good if the node had no `sha256sum`. An empty map means unknown, never "identical".
    """
    join = fs.join if fs is not None else os.path.join
    folder = join(job_dir, "environments")
    found = {}
    try:
        lister = fs.listdir if fs is not None else (
            lambda p: sorted((n, False) for n in os.listdir(p)))
        entries = lister(folder)
    except Exception:
        return found
    for name, _is_dir in entries:
        if not name.endswith(".sha256"):
            continue
        try:
            path = join(folder, name)
            text = fs.read_text(path) if fs is not None else open(path, encoding="utf-8").read()
        except Exception:
            continue
        value = text.strip()
        if value and value != "unavailable":
            found[name[: -len(".sha256")]] = value
    return found


def read(job_dir: str, fs=None, environments_too: bool = True) -> Optional[Dict[str, Any]]:
    """The manifest, with the environment digests the run recorded merged in under `environments_resolved`."""
    destination = path_for(job_dir, fs=fs)
    # None means "no manifest", which callers report as a run that predates it; any other failure must say what it was.
    present = fs.is_file(destination) if fs is not None else os.path.isfile(destination)
    if not present:
        return None
    if fs is not None:
        text = fs.read_text(destination)
    else:
        with open(destination, encoding="utf-8") as handle:
            text = handle.read()
    try:
        record = json.loads(text)
    except ValueError as error:
        raise ValueError(f"{destination} is not valid JSON: {error}") from error
    # Off for callers that need only the identity: over ssh each digest is its own round trip.
    resolved = environments(job_dir, fs=fs) if environments_too else {}
    if resolved:
        record["environments_resolved"] = resolved
    return record


def environment_differences(left: Dict[str, Any], right: Dict[str, Any]) -> List[str]:
    """How the two runs' environments differ, saying so plainly when one of them did not record any.

    Silence here would claim the environments match, which is the one thing an absent record cannot support.
    """
    mine = left.get("environments_resolved") or {}
    theirs = right.get("environments_resolved") or {}
    if not mine or not theirs:
        if mine or theirs:
            return ["environment: one run recorded no environment digests, so the environments "
                    "cannot be compared"]
        # Neither recorded one. That is not a difference -- comparing a record with itself must
        # stay clean -- but it is not agreement either, so `compared()` reports it separately.
        return []
    lines = []
    for name in sorted(set(mine) | set(theirs)):
        if name not in mine or name not in theirs:
            lines.append("environment " + name + ": present in only one run")
        elif mine[name] != theirs[name]:
            lines.append("environment " + name + ": " + mine[name][:12] + " -> " + theirs[name][:12])
    return lines


def _declared_environments(record: Dict[str, Any]) -> set:
    return {env for tool in record.get("tools") or [] for env in tool.get("environments") or [] if env}


def environments_were_compared(left: Dict[str, Any], right: Dict[str, Any]) -> bool:
    """Did both runs record a digest for every environment their steps declare? Anything less leaves the comparison partial."""
    mine = left.get("environments_resolved") or {}
    theirs = right.get("environments_resolved") or {}
    if not mine or not theirs:
        return False
    declared = _declared_environments(left) | _declared_environments(right)
    return all(env in mine and env in theirs for env in declared)


def compare(left: Dict[str, Any], right: Dict[str, Any]) -> List[str]:
    """Why two runs would differ, in the order a reader would ask.

    An equal hash returns nothing. Otherwise every difference is named concretely, because "the manifests differ" is the answer that sends someone diffing two JSON files by hand.
    """
    differences = environment_differences(left, right)
    # The hash is sealed at save(); the environment digests arrive later, from the node that ran. So an equal hash means the same plan, never the same run.
    if not differences and left.get("hash") and left.get("hash") == right.get("hash"):
        return []
    for key in sorted(set(left) | set(right)):
        if key in VOLATILE or key in ("tools", "environments_resolved"):
            continue
        a, b = left.get(key), right.get(key)
        # A 1.5.0 manifest holds a short sha and no dirty flag; neither is a difference by itself.
        if key == "commit" and a and b and (str(a).startswith(str(b)) or str(b).startswith(str(a))):
            continue
        if key == "dirty" and (a is None or b is None):
            continue
        if a != b:
            differences.append(key + ": " + repr(a) + " -> " + repr(b))

    mine = {tool.get("order"): tool for tool in left.get("tools") or []}
    theirs = {tool.get("order"): tool for tool in right.get("tools") or []}
    for order in sorted(set(mine) | set(theirs), key=lambda o: (o is None, o)):
        a, b = mine.get(order), theirs.get(order)
        if a is None or b is None:
            present = b if a is None else a
            verb = "added" if a is None else "removed"
            differences.append("step " + str(order) + ": " + verb + " " + str(present.get("tool")))
            continue
        label = "step " + str(order) + " " + str(a.get("tool"))
        for key in ("tool", "class", "tool_version", "environments", "container_image"):
            if a.get(key) != b.get(key):
                differences.append(label + ": " + key + " " + repr(a.get(key)) + " -> " + repr(b.get(key)))
        left_record = a.get("parameters") or {}
        right_record = b.get("parameters") or {}
        left_params = left_record.get("resolved") or {}
        right_params = right_record.get("resolved") or {}
        for name in sorted(set(left_params) | set(right_params)):
            if left_params.get(name) != right_params.get(name):
                differences.append(label + ": " + name + " " + repr(left_params.get(name))
                                   + " -> " + repr(right_params.get(name)))
        for section in ("forwarded", "input_files", "passed"):
            mine_section = left_record.get(section) or {}
            theirs_section = right_record.get(section) or {}
            for name in sorted(set(mine_section) | set(theirs_section)):
                if section == "passed" and name in left_params and name in right_params:
                    continue
                if mine_section.get(name) != theirs_section.get(name):
                    differences.append(label + ": " + section + " " + name + " "
                                       + repr(mine_section.get(name)) + " -> "
                                       + repr(theirs_section.get(name)))
    if not differences and left.get("hash") != right.get("hash"):
        differences.append("hash: " + str(left.get("hash")) + " -> " + str(right.get("hash"))
                           + " (the difference is in a field this comparison does not list)")
    return differences
