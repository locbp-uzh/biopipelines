"""Single home for the framework's contract checks.

The contracts themselves are stated in `docs/developer_manual.md`. Historically they existed only as prose, which meant they were violated silently and — in at least one case — taught wrongly by the manual's own tool template. This module makes them executable without making the framework rigid: every check has a severity, and the default is `warn`, so a violation is reported without breaking a pipeline that would otherwise run.

Severity resolution, most specific first:

  1. `BIOPIPELINES_ENFORCE_<CHECK>` environment variable  (e.g. BIOPIPELINES_ENFORCE_COMPOUNDS_FORMAT=raise)
  2. `SEVERITY[<check>]` set in this module
  3. `DEFAULT_SEVERITY`

There is no global switch. `BIOPIPELINES_ENFORCE` promised two things and delivered neither: `=off` did not silence everything, because the rows reporting a working construction rather than a fault were skipped in both directions, and `=raise` — the setting the manual named for CI — refused every pipeline using `FORWARD_UNKNOWN_KWARGS`, a supported feature. The `INFORMATIONAL` set existed only to patch the second, and patching it in one direction is what broke the other. Name the checks you want fatal, one variable each.

Adding a contract: write a `check_*` function that returns a `Violation` or `None`, give it an entry in `SEVERITY`, and call it from the one place the contract applies. Keep the check pure — callers decide what to do with the result via `report()`.

Leftover constructor kwargs: a tool that types only a curated subset of its upstream options can set `FORWARD_UNKNOWN_KWARGS` (see `BaseConfig`) to render anything it does not recognise onto the upstream command line. `check_kwargs` is the one entry point `BaseConfig.__init__` calls; it splits the leftovers into probable misspellings (reported whether or not the tool forwards, because a near-miss on a typed name is a typo either way) and genuine passthrough (announced, so a forwarded flag is visible rather than silent). Whether the marker can work at all is not a matter of taste but of what the wrapper emits, so `assert_can_forward` derives it from the wrapper's own source instead of from a hand-kept list of tools.

Naming the offender: a `DataStream` does not know which tool built it and must stay constructible standalone (pipe scripts import it directly at runtime), so `BaseConfig` wraps each subclass's `get_output_files()` in `tool_context(TOOL_NAME)` and the stream checks read `where` from there. A tool whose stream names come from its user rather than from its wrapper (Mock(streams={...}), Scripting) sets `USER_STREAM_NAMES = True`, which drops the `stream_name` check for that tool only.
"""

import contextlib
import contextvars
import difflib
import inspect
import os
import sys
from dataclasses import dataclass
from typing import Any, Dict, Iterable, Iterator, List, Optional, Set, Union

OFF, WARN, RAISE = "off", "warn", "raise"
_VALID_SEVERITIES = (OFF, WARN, RAISE)

# These checks address tool authors, not the scientists running pipelines, so warning is the right
# default: a contract slip should not stop someone's run. `BIOPIPELINES_ENFORCE=1` promotes every
# check to raise for a debugging session, and the truthy spellings are accepted because a switch that
# silently ignores `True` is worse than one that rejects it.
_SEVERITY_ALIASES = {
    "1": RAISE, "true": RAISE, "yes": RAISE, "on": RAISE, "all": RAISE, "strict": RAISE,
    "0": OFF, "false": OFF, "no": OFF, "none": OFF,
}


_warned_unrecognized: Set[str] = set()


def _resolve_severity(value: Optional[str], variable: str = "") -> Optional[str]:
    """A severity name, or one of the truthy spellings a caller is likely to reach for.

    An unrecognized value says so once instead of being ignored, which is the same failure the aliases exist to prevent.
    """
    if value is None:
        return None
    text = value.strip().lower()
    if text in _VALID_SEVERITIES:
        return text
    resolved = _SEVERITY_ALIASES.get(text)
    if resolved is None and variable not in _warned_unrecognized:
        _warned_unrecognized.add(variable)
        accepted = ", ".join(_VALID_SEVERITIES + ("1", "0"))
        print(f"[contract] {variable}={value!r} is not a severity and was ignored; "
              f"accepted values are {accepted}.", file=sys.stderr)
    return resolved

DEFAULT_SEVERITY = WARN

# Per-check severity. Promote a row to RAISE once the warn phase shows the tree is clean.
SEVERITY: Dict[str, str] = {
    "stream_name": WARN,
    "compounds_format": WARN,
    "value_based_format": WARN,
    "unknown_kwargs": WARN,
    "forwarded_kwargs": WARN,
    "pattern_selection": WARN,
    "deprecated_alias": WARN,
    "code_only_ligand": WARN,
    "id_match_consistency": WARN,
}

# Every literal stream name declared across biopipelines/ and pipe_scripts/ as of
# the 1.4 audit; the developer manual documents five. Names computed at runtime
# (Consensus, Pool, Panda, ReMap) inherit a registered upstream name. An
# unregistered name is a warning, not an error: adding one here is how a new
# stream is declared, and nothing is forbidden.
KNOWN_STREAM_NAMES: Set[str] = {
    "structures", "sequences", "compounds", "msas",
    "accessibility", "aggregation", "annotated", "binding", "distances",
    "designs", "dssp", "fasta", "final_ranked_designs", "grids", "images",
    "intermediate_designs_inverse_folded", "movies", "plots",
    "precursors", "renders", "residues", "rmsf", "routes", "sessions",
    "ss", "trajectories",
}


# Constructor keys `BaseConfig.__init__` reads itself; every tool accepts them.
RESERVED_KWARGS: Set[str] = {"name", "pipeline", "resources", "dependencies", "_internal"}

# An entry here is a blind spot, not a fix: add one only for a tool that genuinely forwards keys nothing names.
UNKNOWN_KWARGS_EXEMPT_TOOLS: Set[str] = set()

# The `BaseConfig` accessors that put a forwarded token somewhere it can become argv. `extra_args_echo` is deliberately absent: it only names the tokens in the step log, so a wrapper calling it alone still drops them.
FORWARDING_ACCESSORS = ("extra_args_tokens", "extra_args_bash_tokens", "extra_args_bash")

# difflib ratio above which an unknown key is called a misspelling rather than an intentional upstream flag. Tuned so the genuine ProteinMPNN flag `omit_AAs` (0.76 against the typed `omit_AA_jsonl`) still forwards, while `num_desgins`/`num_designs` (0.91) and `symetry`/`symmetry` (0.93) are reported as typos. `code`/`codes` and `contig`/`contigs` scored above the cutoff too, but both are bound aliases now and resolve before this check sees them; `sampling_temp`/`temperature` scores 0.33 and is deliberately out of reach, since the two names mean the same thing but share almost no letters.
TYPO_SIMILARITY_CUTOFF = 0.80


class ContractViolation(Exception):
    """Raised when a contract check resolves to RAISE severity."""


@dataclass(frozen=True)
class Violation:
    check: str
    message: str
    hint: str = ""

    def render(self) -> str:
        text = f"[contract:{self.check}] {self.message}"
        return f"{text}\n    {self.hint}" if self.hint else text


def severity_for(check: str) -> str:
    variable = f"BIOPIPELINES_ENFORCE_{check.upper()}"
    env_specific = _resolve_severity(os.environ.get(variable), variable)
    if env_specific:
        return env_specific
    if check in SEVERITY:
        return SEVERITY[check]
    return DEFAULT_SEVERITY


# A violation inside a per-id loop would otherwise print hundreds of identical lines.
_reported: Set[str] = set()


def reset_reported() -> None:
    """Clear the dedup cache. For tests, and for long-lived sessions."""
    _reported.clear()


def report(violation: Optional[Violation]) -> None:
    """Act on a check result according to its severity. `None` is a pass."""
    if violation is None:
        return
    level = severity_for(violation.check)
    if level == OFF:
        return
    if level == RAISE:
        raise ContractViolation(violation.render())
    key = violation.render()
    if key not in _reported:
        _reported.add(key)
        print(violation.render(), file=sys.stderr)


def check_all(*violations: Optional[Violation]) -> None:
    """Report every violation in order. Raises on the first RAISE-severity one."""
    for violation in violations:
        report(violation)


# --- calling context -------------------------------------------------------

@dataclass(frozen=True)
class ToolContext:
    """Who is building the thing a check is about."""

    where: str = ""
    user_stream_names: bool = False


_NO_CONTEXT = ToolContext()

# A DataStream must stay usable standalone, so the tool name arrives out of band instead of through its signature.
_context: contextvars.ContextVar[ToolContext] = contextvars.ContextVar(
    "biopipelines_contract_context", default=_NO_CONTEXT
)


def current_context() -> ToolContext:
    """The innermost active `tool_context`, or an empty one."""
    return _context.get()


@contextlib.contextmanager
def tool_context(where: str, user_stream_names: bool = False) -> Iterator[None]:
    """Attribute every check inside the block to `where`.

    `user_stream_names=True` marks a tool whose stream names come from its user (Mock(streams={...}), Scripting), which no registry can know in advance.
    """
    token = _context.set(ToolContext(where=where, user_stream_names=user_stream_names))
    try:
        yield
    finally:
        _context.reset(token)


# --- contracts -------------------------------------------------------------

def check_stream_name(name: str, where: str = "") -> Optional[Violation]:
    """A stream name should be one the framework knows about."""
    if not name or name in KNOWN_STREAM_NAMES:
        return None
    context = f" (in {where})" if where else ""
    return Violation(
        check="stream_name",
        message=f"unregistered stream name {name!r}{context}",
        hint=("If this is a new stream, add it to KNOWN_STREAM_NAMES in "
              "contract_enforcement.py. If it is a typo, consumers reading it "
              "back by name will silently get nothing."),
    )


def check_compounds_format(
    name: str, files: Union[str, List[str]], format: str, where: str = ""
) -> Optional[Violation]:
    """The compounds stream carries chemistry as csv, never coordinate files.

    developer_manual.md:183. The manual's own tool template declared this as
    "sdf" for years, and five tools copied it.
    """
    if name != "compounds" or (format or "").lower() == "csv":
        return None
    context = f" (in {where})" if where else ""
    return Violation(
        check="compounds_format",
        message=f"compounds stream declared format={format!r}, expected 'csv'{context}",
        hint=("compounds carries ligand chemistry and identity in its map_table; "
              "coordinates belong on a structures stream. See the Ligand Contract "
              "in docs/developer_manual.md."),
    )


def check_value_based_format(
    name: str, files: Union[str, List[str]], format: str,
    ids: Optional[List[str]] = None, map_table: str = "", where: str = ""
) -> Optional[Violation]:
    """A stream with no files keeps its content in the map_table, so it is csv.

    developer_manual.md:167.
    """
    if isinstance(files, str) or files:
        return None
    if (format or "").lower() == "csv":
        return None
    # An empty placeholder means "this tool emits none of these", not "value-based".
    if not ids and not map_table:
        return None
    context = f" (in {where})" if where else ""
    return Violation(
        check="value_based_format",
        message=(f"value-based stream {name!r} (files=[]) declared "
                 f"format={format!r}, expected 'csv'{context}"),
        hint=("With no files the content lives in the map_table, which is csv. "
              "Consumers dispatch on format to decide how to read a stream."),
    )


def renders_forwarded_tokens(cls: type) -> bool:
    """Whether any wrapper class in `cls`'s ancestry reaches for a forwarded token.

    Read from source rather than from behavior because there is nothing to observe at class-creation time: the tokens only exist once a user passes an unrecognized kwarg. Every class up to the one defining the accessors (`BaseConfig`) is searched, so a subclass that inherits its script generation — SolubleMPNN from ProteinMPNN — answers on its parent's source. A class whose source cannot be read (a REPL, a zipped install) is given the benefit of the doubt, since an unreadable file is not evidence of a defect.
    """
    for klass in cls.__mro__:
        if any(name in vars(klass) for name in FORWARDING_ACCESSORS):
            return False
        try:
            source = inspect.getsource(klass)
        except (OSError, TypeError):
            return True
        if any(name in source for name in FORWARDING_ACCESSORS):
            return True
    return False


def assert_can_forward(cls: type) -> None:
    """Refuse a forwarding marker on a wrapper that never renders the forwarded tokens.

    Most wrappers emit `python pipe_<tool>.py --config <json>` and let the pipe script assemble argv, so a token the wrapper does not itself place has no command line to join and the marker is accepted and then silently dropped. What separates the two cases is whether the wrapper's own code calls one of `FORWARDING_ACCESSORS`, and that is what is checked here. This is the necessary condition, paid once per forwarding class at import; that the tokens actually reach the generated step script is the sufficient one, and `tests/test_forwarded_kwargs.py` checks that end to end for every tool the marker is found on.

    Not severity-gated on purpose: a marker that cannot work must fail at import rather than be reported and then ignored.
    """
    if renders_forwarded_tokens(cls):
        return
    tool = getattr(cls, "TOOL_NAME", cls.__name__)
    accessors = ", ".join(f"{name}()" for name in FORWARDING_ACCESSORS)
    raise ContractViolation(
        f"[contract:forwarded_kwargs] {tool} declares FORWARD_UNKNOWN_KWARGS, but "
        f"nothing in {cls.__name__} or its wrapper ancestry calls {accessors}, so the "
        f"forwarded tokens are rendered and then dropped. A tool whose script is "
        f"`python pipe_{tool.lower()}.py --config <json>` has no wrapper-written "
        f"command line for them to join: give the pipe script a typed parameter "
        f"instead, or emit the upstream command from the wrapper and interpolate "
        f"self.extra_args_bash() into it."
    )


def probable_typos(
    unknown: Iterable[str], known: Iterable[str], cutoff: float = TYPO_SIMILARITY_CUTOFF
) -> Dict[str, str]:
    """Map each unknown key that is lexically close to a real parameter name onto that name."""
    candidates = sorted(set(known) | RESERVED_KWARGS)
    matched: Dict[str, str] = {}
    for key in unknown:
        close = difflib.get_close_matches(key, candidates, n=1, cutoff=cutoff)
        if close:
            matched[key] = close[0]
    return matched


def _unknown_keys(kwargs: Iterable[str]) -> List[str]:
    return sorted(k for k in kwargs if k not in RESERVED_KWARGS)


def check_no_unknown_kwargs(
    tool: str, kwargs: Iterable[str], where: str = "", known: Iterable[str] = ()
) -> Optional[Violation]:
    """Every tool constructor ends in `**kwargs`, so a typo is silently accepted.

    Called from `BaseConfig.__init__`, i.e. after the subclass has bound its own named parameters, so anything still here is either a framework key or a mistake. For a tool that does not set `FORWARD_UNKNOWN_KWARGS` the value has no effect at all: `BaseConfig` stores it in `self.params` and never reads it again.

    `known` is the tool's real parameter names; when given, a key close to one of them is named in the message so a typo reads as a typo.
    """
    if tool in UNKNOWN_KWARGS_EXEMPT_TOOLS:
        return None
    unknown = _unknown_keys(kwargs)
    if not unknown:
        return None
    typos = probable_typos(unknown, known)
    context = where or tool
    listed = ", ".join(
        f"{k!r} (did you mean {typos[k]!r}?)" if k in typos else repr(k)
        for k in unknown
    )
    return Violation(
        check="unknown_kwargs",
        message=f"{context} got unknown constructor parameter(s): {listed}",
        hint=("A tool's **kwargs exists for the framework keys "
              f"({', '.join(sorted(RESERVED_KWARGS))}); anything else is parked in "
              ".params and never read, so the value you passed has no effect. Check "
              "the spelling against the constructor signature. If it is a real "
              f"upstream flag, {tool} would have to set FORWARD_UNKNOWN_KWARGS to "
              "pass it through."),
    )


def check_probable_typo(
    tool: str, kwargs: Iterable[str], where: str = "", known: Iterable[str] = ()
) -> Optional[Violation]:
    """A near-miss on a real parameter name, reported whether or not the tool forwards.

    A key one letter away from a typed parameter is almost certainly a misspelling rather than an upstream flag the wrapper does not know, so a forwarding tool must still say so — otherwise the typo silently becomes a bogus flag and only upstream complains, on a compute node.
    """
    if tool in UNKNOWN_KWARGS_EXEMPT_TOOLS:
        return None
    typos = probable_typos(_unknown_keys(kwargs), known)
    if not typos:
        return None
    context = where or tool
    listed = "; ".join(
        f"unknown parameter {k!r}; did you mean {v!r}?" for k, v in sorted(typos.items())
    )
    return Violation(
        check="unknown_kwargs",
        message=f"{context}: {listed}",
        hint=(f"{tool} forwards unrecognised keys to its upstream command line, so "
              "this one will be passed through as written; upstream will most likely "
              "reject it. Fix the spelling, or rename it if the near-match is a "
              "coincidence."),
    )


def check_forwarded_kwargs(
    tool: str, kwargs: Iterable[str], where: str = "", known: Iterable[str] = ()
) -> Optional[Violation]:
    """Announce the keys a forwarding tool is about to put on the upstream command line.

    Informational, not a fault: a forwarded argument is a wanted capability, and printing it is what keeps it from being an invisible one. Keys reported as probable typos are left to `check_probable_typo` so one mistake produces one line.
    """
    unknown = _unknown_keys(kwargs)
    typos = probable_typos(unknown, known)
    forwarded = [k for k in unknown if k not in typos]
    if not forwarded:
        return None
    context = where or tool
    listed = ", ".join(repr(k) for k in forwarded)
    return Violation(
        check="forwarded_kwargs",
        message=f"{context} is forwarding to {tool}: {listed}",
        hint=("These are not typed by the wrapper, so they are rendered onto the "
              "upstream command line as written and neither validated against "
              f"{tool}'s own options nor recorded in its output schema. Silence this "
              "line with BIOPIPELINES_ENFORCE_FORWARDED_KWARGS=off."),
    )


def check_kwargs(
    tool: str, kwargs: Iterable[str], known: Iterable[str] = (),
    forwards: bool = False, where: str = ""
) -> None:
    """Report the leftover-kwargs contracts for the one place a tool is constructed."""
    if forwards:
        check_all(
            check_probable_typo(tool, kwargs, where, known),
            check_forwarded_kwargs(tool, kwargs, where, known),
        )
    else:
        report(check_no_unknown_kwargs(tool, kwargs, where, known))


def check_code_only_ligand(codes, where: str = "") -> Optional[Violation]:
    """A `codes`-only Ligand names a residue but carries no chemistry.

    One `codes` parameter serves both modes, so the mode is decided by what else was passed. That is the simpler API, but it means a forgotten `smiles=` produces a chemistry-free stub instead of an error, and the failure would otherwise surface in whatever downstream tool needed the molecule.
    """
    shown = ", ".join(str(c) for c in codes) if isinstance(codes, (list, tuple)) else str(codes)
    return Violation(
        check="code_only_ligand",
        message=f"Ligand(codes={shown!r}) carries no chemistry, only a residue code{f' (in {where})' if where else ''}",
        hint="Pass lookup or smiles to retrieve chemical information.",
    )


def check_stream(
    name: str, files: Union[str, List[str]], format: str,
    ids: Optional[List[str]] = None, map_table: str = "", where: str = ""
) -> None:
    """Every stream-level contract, for the one place a DataStream is built."""
    context = current_context()
    where = where or context.where
    compounds = check_compounds_format(name, files, format, where)
    # A compounds stream is value-based by definition, so reporting both would
    # print two lines for one mistake.
    value_based = (None if compounds
                   else check_value_based_format(name, files, format, ids, map_table, where))
    stream_name = None if context.user_stream_names else check_stream_name(name, where)
    check_all(stream_name, compounds, value_based)


# --- pattern-selection contract --------------------------------------------

def check_pattern_selection(
    patterns: List[str], row_ids: List[str], selected: List[str], where: str = ""
) -> Optional[Violation]:
    """A pattern set that selects none of a non-empty row set has silently done nothing.

    `select_ids` returning `[]` is indistinguishable from "there was nothing to select", and exact matching makes a near-miss more likely, not less: a pattern for `design_1` against zero-padded rows `design_01…` correctly covers nothing at all.
    """
    if selected or not patterns or not row_ids:
        return None
    context = f" (in {where})" if where else ""
    shown_p = ", ".join(patterns[:4]) + (", ..." if len(patterns) > 4 else "")
    shown_r = ", ".join(row_ids[:4]) + (", ..." if len(row_ids) > 4 else "")
    return Violation(
        check="pattern_selection",
        message=(f"pattern set [{shown_p}] selected 0 of {len(row_ids)} "
                 f"map_table row(s) [{shown_r}]{context}"),
        hint=("Matching is exact outside the unresolved slots, so a shape "
              "mismatch (zero padding, a missing suffix delimiter, a renamed "
              "prefix) covers nothing. Compare the pattern against the "
              "map_table's id column."),
    )


# --- id-match consistency contract -----------------------------------------

# Ids named in one line before it degrades into a count; a pathological lookup must not print thousands of them.
MAX_NAMED_IDS = 8

# Mirrors id_map_utils.TIER_EXACT, copied so this module stays importable without the matcher; the shallowest tier reads as a plain fact, a deeper one names the route it took.
_EXACT_TIER = "exact"


def _render_tier_counts(tier_counts: Dict[str, int], unmatched: int) -> str:
    parts = [
        (f"{count} {tier}" if tier == _EXACT_TIER else f"{count} via {tier}")
        for tier, count in tier_counts.items()
    ]
    if unmatched:
        parts.append(f"{unmatched} unmatched")
    return ", ".join(parts)


def _name_ids(ids, limit: int = MAX_NAMED_IDS) -> str:
    named = list(ids)
    listed = ", ".join(str(i) for i in named[:limit])
    extra = len(named) - limit
    return f"{listed}, +{extra} more" if extra > 0 else listed


def check_id_match_consistency(
    tier_counts: Dict[str, int], minority_ids: Iterable[str], score: float,
    unmatched: int = 0, where: str = "", limit: int = MAX_NAMED_IDS
) -> Optional[Violation]:
    """Within one lookup, the ids should all resolve the same way.

    `get_mapped_ids` answers an id from the first of five tiers that fits, and the last three are string operations on the id, so they can return a different row's value: `design_9` can be handed `design_1`'s cell. The permissiveness is wanted — an upstream tool that renames ids relies on it — so this reports rather than refuses, and it reports on the shape of the whole lookup rather than on any single match, because that is the level where the two situations look different. Every row degrading identically is a systematic rename and passes; a few rows degrading differently from their neighbors is the case where a value probably came from the wrong row, and those ids are named so they can be checked.

    Log-only: matching behavior is unchanged, and nothing reads the score back. `score` and `minority_ids` come from `id_map_utils.score_id_match_tiers`, which defines the scoring; a uniform lookup scores 1.0 and passes here.
    """
    if not tier_counts or score >= 1.0:
        return None
    where = where or current_context().where
    context = f" (in {where})" if where else ""
    total = sum(tier_counts.values()) + unmatched
    named = _name_ids(minority_ids, limit)
    return Violation(
        check="id_match_consistency",
        message=(f"id matching{context} scored {score:.2f} — {total} rows: "
                 f"{_render_tier_counts(tier_counts, unmatched)}"
                 f"{f' ({named})' if named else ''}"),
        hint=("Uniform degradation is normally an upstream rename and is fine; these "
              "ids resolved through a different tier than the rest of the lookup, and "
              "the deeper tiers match on the id string alone, so such a row can carry "
              "another row's value. Compare the named ids against the table's id "
              "column. Nothing was gated on this: silence the line with "
              "BIOPIPELINES_ENFORCE_ID_MATCH_CONSISTENCY=off."),
    )
