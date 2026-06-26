"""Batch-scheduler backends.

Each backend is a pure string-generator: it turns the scheduler-agnostic
resource/dependency data carried by ``Pipeline`` into the native directive
text for one scheduler (SLURM ``#SBATCH``, LSF ``#BSUB``, PBS ``#PBS``).
``pipeline.py`` owns the surrounding batch-script template and delegates only
the directive fragments to the active backend, so the SLURM output stays
byte-for-byte identical to the pre-seam code.

The shared ``<JOBID_BATCH_NNN>`` placeholder scheme is unchanged across
schedulers; only the directive that wraps those placeholders differs. The
matching submit/job-id-capture logic lives in the bash ``submit`` wrapper.
"""

import re
from typing import Any, List, Optional, Tuple

_MEM_RE = re.compile(r"^\s*(\d+(?:\.\d+)?)\s*([KMGT]?B?)\s*$", re.IGNORECASE)
_MEM_UNIT_MB = {"": 1, "B": 1e-6, "KB": 1e-3, "K": 1e-3, "MB": 1, "M": 1,
                "GB": 1e3, "G": 1e3, "TB": 1e6, "T": 1e6}


def _memory_to_mb(memory: str) -> Optional[int]:
    """Parse a framework memory string ("16GB", "512MB", "8000") to integer MB.

    Returns None if the value cannot be parsed (caller passes it through).
    """
    if not isinstance(memory, str):
        return None
    m = _MEM_RE.match(memory)
    if not m:
        return None
    value, unit = float(m.group(1)), m.group(2).upper()
    factor = _MEM_UNIT_MB.get(unit)
    if factor is None:
        return None
    return int(round(value * factor))


def _hms_to_hhmm(t: str) -> str:
    """Convert SLURM-style walltime to LSF's [HH:]MM form (no seconds).

    "24:00:00" -> "24:00", "90:00" -> "90:00", "1-00:00:00" -> "24:00".
    Unparseable input is returned unchanged.
    """
    if not isinstance(t, str):
        return t
    days, rest = (t.split("-", 1) + [""])[:2] if "-" in t else ("0", t)
    parts = rest.split(":")
    try:
        nums = [int(p) for p in parts]
        d = int(days)
    except ValueError:
        return t
    if len(nums) == 3:
        h, m, _s = nums
    elif len(nums) == 2:
        h, m = nums
    elif len(nums) == 1:
        h, m = nums[0], 0
    else:
        return t
    h += d * 24
    return f"{h}:{m:02d}"


def _validate_directive_value(field: str, value: Any) -> None:
    """Reject values that would inject extra directive lines.

    A newline or carriage return in any interpolated resource value would
    smuggle an extra ``#SBATCH``/``#BSUB``/``#PBS`` directive into the
    generated script. Applies to every scheduler.
    """
    if value is None:
        return
    s = value if isinstance(value, str) else str(value)
    if "\n" in s or "\r" in s:
        raise ValueError(
            f"{field!r}={value!r} contains a newline or carriage return. "
        )


class SchedulerBackend:
    """Base class. Subclasses emit native directive fragments.

    ``gpu_setup`` is plain bash (``nvidia-smi``) and identical for every
    batch scheduler, so it lives here once.
    """

    name = "base"
    # The directive-option flag a sibling resource dict is keyed under
    # (``slurm_options`` / ``lsf_options`` / ``pbs_options``).
    options_key = "slurm_options"

    def header_directives(self, memory: str, time: str, output_name: str) -> str:
        raise NotImplementedError

    def gpu_directive(self, gpu_spec: Optional[str], gpus: int = 1) -> Tuple[str, List[str]]:
        """Return ``(directive_text, warnings)`` for a GPU request."""
        raise NotImplementedError

    def extra_options(self, options: dict) -> str:
        raise NotImplementedError

    def dependency_directive(self, afterok: List[str], after: List[str]) -> str:
        """Native dependency directive given placeholder/job-id terms.

        ``afterok``/``after`` are lists of ``<JOBID_BATCH_NNN>`` placeholders
        (or real ids for external deps). Returns "" when both are empty.
        """
        raise NotImplementedError

    def email_directive(self, email: str) -> str:
        raise NotImplementedError

    def gpu_setup(self, gpu_spec: Optional[str]) -> str:
        if gpu_spec is None or gpu_spec == "none" or gpu_spec == "":
            return ""
        return """
# Display GPU information
gpu_type=$(nvidia-smi --query-gpu=gpu_name --format=csv,noheader 2>/dev/null || echo "Unknown")
echo "GPU Type: $gpu_type"
"""


class SlurmBackend(SchedulerBackend):
    """Faithful extraction of the original ``#SBATCH`` generators.

    Output must remain byte-for-byte identical to the pre-seam code; the
    test-suite asserts exact strings.
    """

    name = "slurm"
    options_key = "slurm_options"

    def header_directives(self, memory: str, time: str, output_name: str) -> str:
        return (
            f"#SBATCH --mem={memory}\n"
            f"#SBATCH --time={time}\n"
            f"#SBATCH --output={output_name}\n"
            f"#SBATCH --begin=now+0hour"
        )

    def gpu_directive(self, gpu_spec: Optional[str], gpus: int = 1) -> Tuple[str, List[str]]:
        n = gpus if gpus else 1
        if gpu_spec is None or gpu_spec == "none" or gpu_spec == "":
            return "", []
        elif gpu_spec == "high-memory":
            return f"#SBATCH --gpus={n}\n#SBATCH --constraint=\"GPUMEM32GB|GPUMEM80GB|GPUMEM96GB\"", []
        elif gpu_spec == "gpu" or gpu_spec == "any":
            return f"#SBATCH --gpus={n}", []
        elif gpu_spec.startswith("!"):
            excluded_model = gpu_spec[1:]
            if excluded_model.upper() == "L4":
                return f"#SBATCH --gpus={n}\n#SBATCH --constraint=\"GPUMEM32GB|GPUMEM80GB|GPUMEM96GB\"", []
            else:
                return f"#SBATCH --gpus={n}\n#SBATCH --constraint=\"~GPU{excluded_model}\"", []
        elif gpu_spec in ["24GB", "32GB", "80GB", "96GB"] or "|" in gpu_spec:
            if "|" in gpu_spec:
                memory_options = gpu_spec.split("|")
                constraint_parts = [f"GPUMEM{mem}" for mem in memory_options]
                constraint = "|".join(constraint_parts)
            else:
                constraint = f"GPUMEM{gpu_spec}"
            return f"#SBATCH --gpus={n}\n#SBATCH --constraint=\"{constraint}\"", []
        else:
            return f"#SBATCH --gpus={gpu_spec}:{n}", []

    def extra_options(self, options: dict) -> str:
        if not options:
            return ""
        lines = []
        for param, value in options.items():
            _validate_directive_value(f"slurm_options[{param!r}]", value)
            if param == "cpus":
                slurm_param = "cpus-per-task"
            else:
                slurm_param = param.replace('_', '-')
            lines.append(f"#SBATCH --{slurm_param}={value}")
        return "\n" + "\n".join(lines)

    def dependency_directive(self, afterok: List[str], after: List[str]) -> str:
        dep_terms = []
        if afterok:
            dep_terms.append("afterok:" + ":".join(afterok))
        if after:
            dep_terms.append("after:" + ":".join(after))
        if not dep_terms:
            return ""
        return "\n#SBATCH --dependency=" + ",".join(dep_terms)

    def email_directive(self, email: str) -> str:
        if email == "":
            return ""
        return f"\n#SBATCH --mail-type=END,FAIL\n#SBATCH --mail-user={email}"


class LsfBackend(SchedulerBackend):
    """LSF (``bsub``/``#BSUB``).

    GPU/mem use best-effort translation: a GPU count maps to
    ``-gpu "num=<n>"``; a named model adds ``:gmodel=<model>``. SLURM
    constraint-style specs (``high-memory``, ``A|B``, ``!model``) have no clean
    LSF equivalent, so they fall back to the count-only directive plus a
    warning.
    """

    name = "lsf"
    options_key = "lsf_options"

    def header_directives(self, memory: str, time: str, output_name: str) -> str:
        # Translate the framework's canonical forms ("16GB", "HH:MM:SS") into
        # LSF-portable ones: memory as MB on -M plus a matching rusage request
        # (works regardless of LSF_UNIT_FOR_LIMITS), walltime as [HH:]MM (LSF -W
        # has no seconds field). Sites with different conventions override via
        # lsf_options. Unparseable values are passed through unchanged.
        mem_mb = _memory_to_mb(memory)
        mem_lines = (
            f'#BSUB -M {mem_mb}\n#BSUB -R "rusage[mem={mem_mb}]"'
            if mem_mb is not None else f"#BSUB -M {memory}"
        )
        walltime = _hms_to_hhmm(time)
        return (
            f"{mem_lines}\n"
            f"#BSUB -W {walltime}\n"
            f"#BSUB -o {output_name}"
        )

    def gpu_directive(self, gpu_spec: Optional[str], gpus: int = 1) -> Tuple[str, List[str]]:
        n = gpus if gpus else 1
        warnings: List[str] = []
        if gpu_spec is None or gpu_spec == "none" or gpu_spec == "":
            return "", warnings
        if gpu_spec in ("gpu", "any"):
            return f'#BSUB -gpu "num={n}"', warnings
        constraint_like = (
            gpu_spec == "high-memory"
            or gpu_spec.startswith("!")
            or "|" in gpu_spec
            or gpu_spec in ("24GB", "32GB", "80GB", "96GB")
        )
        if constraint_like:
            warnings.append(
                f"LSF: GPU spec {gpu_spec!r} has no LSF equivalent; "
                f"requesting {n} generic GPU(s). Use lsf_options for a native "
                f"-gpu/-R resource string."
            )
            return f'#BSUB -gpu "num={n}"', warnings
        # Named model (A100, V100, ...).
        return f'#BSUB -gpu "num={n}:gmodel={gpu_spec}"', warnings

    def extra_options(self, options: dict) -> str:
        if not options:
            return ""
        lines = []
        for param, value in options.items():
            _validate_directive_value(f"lsf_options[{param!r}]", value)
            if param == "cpus":
                lines.append(f"#BSUB -n {value}")
            else:
                lines.append(f"#BSUB -{param} {value}")
        return "\n" + "\n".join(lines)

    def dependency_directive(self, afterok: List[str], after: List[str]) -> str:
        # LSF has no job-id-list join: AND the per-job conditions instead.
        # afterok -> done(id) (exit 0); after (daemon started) -> started(id).
        conds = [f"done({jid})" for jid in afterok]
        conds += [f"started({jid})" for jid in after]
        if not conds:
            return ""
        return "\n#BSUB -w '" + " && ".join(conds) + "'"

    def email_directive(self, email: str) -> str:
        if email == "":
            return ""
        # -N notify at end, -B at begin; -u recipient. END+FAIL ~ -N.
        return f"\n#BSUB -N\n#BSUB -u {email}"


class PbsBackend(SchedulerBackend):
    """PBS/Torque (``qsub``/``#PBS``).

    PBS shares SLURM's ``:``-separated job-id join, so the dependency
    placeholders line up cleanly. mem/walltime/ncpus/ngpus are emitted as
    separate ``-l`` requests, which Torque and OpenPBS accept. Sites that
    require everything in one ``select`` chunk (some PBS Pro configs, e.g.
    ``select=1:ncpus=4:ngpus=1:mem=16gb``) should express that via
    ``pbs_options`` (a ``select=...`` value overrides the generated lines).
    GPU: ``-l select=1:ngpus=``, named model adds ``:gpu_model=``;
    constraint-style specs warn and fall back to a count-only ngpus request.
    """

    name = "pbs"
    options_key = "pbs_options"

    def header_directives(self, memory: str, time: str, output_name: str) -> str:
        return (
            f"#PBS -l mem={memory}\n"
            f"#PBS -l walltime={time}\n"
            f"#PBS -o {output_name}"
        )

    def gpu_directive(self, gpu_spec: Optional[str], gpus: int = 1) -> Tuple[str, List[str]]:
        n = gpus if gpus else 1
        warnings: List[str] = []
        if gpu_spec is None or gpu_spec == "none" or gpu_spec == "":
            return "", warnings
        if gpu_spec in ("gpu", "any"):
            return f"#PBS -l select=1:ngpus={n}", warnings
        constraint_like = (
            gpu_spec == "high-memory"
            or gpu_spec.startswith("!")
            or "|" in gpu_spec
            or gpu_spec in ("24GB", "32GB", "80GB", "96GB")
        )
        if constraint_like:
            warnings.append(
                f"PBS: GPU spec {gpu_spec!r} has no PBS equivalent; "
                f"requesting {n} generic GPU(s). Use pbs_options for a native "
                f"-l select resource string."
            )
            return f"#PBS -l select=1:ngpus={n}", warnings
        return f"#PBS -l select=1:ngpus={n}:gpu_model={gpu_spec}", warnings

    def extra_options(self, options: dict) -> str:
        if not options:
            return ""
        lines = []
        for param, value in options.items():
            _validate_directive_value(f"pbs_options[{param!r}]", value)
            if param == "cpus":
                lines.append(f"#PBS -l ncpus={value}")
            else:
                lines.append(f"#PBS -{param} {value}")
        return "\n" + "\n".join(lines)

    def dependency_directive(self, afterok: List[str], after: List[str]) -> str:
        # PBS separates dependency TYPES with ',' and the job ids WITHIN a
        # type with ':' (e.g. afterok:j1:j2,after:j3).
        dep_terms = []
        if afterok:
            dep_terms.append("afterok:" + ":".join(afterok))
        if after:
            dep_terms.append("after:" + ":".join(after))
        if not dep_terms:
            return ""
        return "\n#PBS -W depend=" + ",".join(dep_terms)

    def email_directive(self, email: str) -> str:
        if email == "":
            return ""
        return f"\n#PBS -m ae\n#PBS -M {email}"


_BACKENDS = {
    "slurm": SlurmBackend,
    "lsf": LsfBackend,
    "pbs": PbsBackend,
}

BATCH_SCHEDULERS = tuple(_BACKENDS.keys())


def get_backend(name: str) -> SchedulerBackend:
    """Return a backend instance for a batch-scheduler name.

    Raises for non-batch schedulers (``colab``/``none``); those are handled
    by the on-the-fly path in ``pipeline.py``, never via a backend.
    """
    try:
        return _BACKENDS[name]()
    except KeyError:
        raise ValueError(
            f"No batch-scheduler backend for {name!r}; "
            f"expected one of {BATCH_SCHEDULERS}."
        )
