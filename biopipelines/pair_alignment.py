# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Pair handling shared by the structure aligners that take a table of pairs (USalign, TMalign).

Reading the pair table, checking both ids of a pair against the structures stream, running the aligner once per pair and mode, parsing its compact output (`-outfmt 1 -m -`) and the rule that an unavailable value is NaN and a row is never dropped live here, so the two tools differ only in the binary, the modes and a few columns.

The module has two halves. `PairAligner` is the configuration-time base class both wrappers extend. Everything else is execution-time code that `pipe_scripts/pipe_pair_alignment.py` drives and the tests exercise directly.
"""

import math
import os
import re
import subprocess
from typing import Any, Callable, Dict, List, Optional, Tuple

import pandas as pd

try:
    from .base_config import BaseConfig, StandardizedOutput, resolve_table_reference
    from .biopipelines_io import TableReference, container_argv_prefix
    from .file_paths import Path
    from .datastream import DataStream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, resolve_table_reference
    from biopipelines_io import TableReference, container_argv_prefix
    from file_paths import Path
    from datastream import DataStream


MODE_SEQUENTIAL = "sequential"
MODE_CIRCULAR_PERMUTATION = "circular_permutation"
MODE_NON_SEQUENTIAL = "non_sequential"

# Compact fasta output, with the rotation matrix printed to stdout after it.
OUTPUT_ARGS = ["-outfmt", "1", "-m", "-"]

# The canonical column order of the `alignments` table. The wrappers repeat it as a literal in get_output_files() so versions/extract_output_schemas.py can read it, and tests/test_pair_alignment.py pins the two together.
ALIGNMENT_COLUMNS = [
    "id", "structure_1", "structure_2", "mode",
    "length_1", "length_2", "aligned_length", "rmsd", "tm_score_1", "tm_score_2",
    "rotation_11", "rotation_12", "rotation_13",
    "rotation_21", "rotation_22", "rotation_23",
    "rotation_31", "rotation_32", "rotation_33",
    "translation_1", "translation_2", "translation_3",
    "aligned_1", "aligned_2",
]
DELTA_COLUMNS = ["id", "structure_1", "structure_2", "delta_tm_1", "delta_tm_2"]
MISSING_COLUMNS = ["id", "removed_by", "kind", "cause"]

TEXT_COLUMNS = ("id", "structure_1", "structure_2", "mode", "aligned_1", "aligned_2")
NUMERIC_COLUMNS = [c for c in ALIGNMENT_COLUMNS if c not in TEXT_COLUMNS]
INTEGER_COLUMNS = ("length_1", "length_2", "aligned_length")

_HEADER = re.compile(
    r"^>(?P<name>.*)\tL=(?P<length>\d+)\td0=(?P<d0>\S+)\tseqID=(?P<seqid>\S+)\tTM-score=(?P<tm>\S+)\s*$")
_SUMMARY = re.compile(r"^#\s*Lali=(?P<aligned>\d+)\tRMSD=(?P<rmsd>\S+)\tseqID_ali=\S+\s*$")
_MATRIX_ROW = re.compile(
    r"^\s*(?P<row>[0-2])\s+(?P<t>-?\d+\.\d+)\s+(?P<u0>-?\d+\.\d+)\s+(?P<u1>-?\d+\.\d+)\s+(?P<u2>-?\d+\.\d+)\s*$")


class PairError(ValueError):
    """A pair that cannot be aligned for a reason that lies in its ids, not in the aligner."""


def pair_id(structure_1: str, structure_2: str) -> str:
    """Id of a pair: the two structure ids joined with '+', the framework's combinatorics separator."""
    return f"{structure_1}+{structure_2}"


def row_id(structure_1: str, structure_2: str, mode: str) -> str:
    """Id of one (pair, mode) row: the pair id with the mode appended after '_', as for any fanned-out id."""
    return f"{pair_id(structure_1, structure_2)}_{mode}"


# ── pair table ───────────────────────────────────────────────────────────────

def read_pairs(reference_1: str, reference_2: str) -> List[Tuple[str, str]]:
    """The distinct (structure_1, structure_2) pairs of the pair table, in first-seen order.

    Both references are `TABLE_REFERENCE:<path>:<column>` strings into one table. Ids are read as text, so an id such as "007" or "NA" survives, and surrounding blanks are stripped. A row repeated verbatim is one pair: its alignment would be identical. Two different pairs that would share an id, which '+' inside a structure id allows, are refused.
    """
    first = TableReference.from_string(reference_1)
    second = TableReference.from_string(reference_2)
    if first.path != second.path:
        raise ValueError(
            f"structure_1 and structure_2 must be columns of one table, got {first.path} and {second.path}")
    if not os.path.exists(first.path):
        raise FileNotFoundError(f"Pair table not found: {first.path}")
    table = pd.read_csv(first.path, dtype=str, keep_default_na=False)
    for column in (first.column, second.column):
        if column not in table.columns:
            raise ValueError(
                f"Column {column!r} is not in the pair table {first.path}; its columns are {list(table.columns)}")
    ids_1 = table[first.column].str.strip()
    ids_2 = table[second.column].str.strip()
    pairs = list(dict.fromkeys(zip(ids_1, ids_2)))
    owner: Dict[str, Tuple[str, str]] = {}
    for pair in pairs:
        earlier = owner.setdefault(pair_id(*pair), pair)
        if earlier != pair:
            raise ValueError(
                f"the pairs {earlier} and {pair} both have the id {pair_id(*pair)!r}, because a structure id "
                f"contains '+'; rename the structure so that the pair ids are distinct")
    return pairs


def resolve_structures(structure_1: str, structure_2: str, files: Dict[str, str],
                       declared: set) -> Tuple[str, str]:
    """The files of a pair, or a PairError naming the id that is not in the structures stream.

    `files` maps each id whose file exists to its path; `declared` is every id the stream declares. An id can be declared and still have no file when an upstream tool filtered it, which the cause says.
    """
    paths = []
    for label, sid in (("structure_1", structure_1), ("structure_2", structure_2)):
        if not sid:
            raise PairError(f"{label} is empty in the pair table")
        if sid in files:
            paths.append(files[sid])
        elif sid in declared:
            raise PairError(f"{label} {sid!r} has no file in the structures stream (removed upstream)")
        else:
            raise PairError(f"{label} {sid!r} is not in the structures stream")
    return paths[0], paths[1]


# ── running the aligner ──────────────────────────────────────────────────────

def build_command(binary: str, path_1: str, path_2: str, mode_args: List[str],
                  container_prefix: str = "") -> List[str]:
    """The argv of one alignment: the binary, the two files, the mode's own flags and the output flags."""
    return container_argv_prefix(container_prefix) + [binary, path_1, path_2] + list(mode_args) + OUTPUT_ARGS


def run_command(cmd: List[str]) -> Tuple[int, str, str]:
    """Run one command and return (exit code, stdout, stderr)."""
    result = subprocess.run(cmd, capture_output=True, text=True, errors="replace")
    return result.returncode, result.stdout, result.stderr


def _last_line(text: str) -> str:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    return lines[-1] if lines else ""


def _finite(value: str, what: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{what} is not finite: {value}")
    return number


def parse_alignment(text: str) -> Dict[str, Any]:
    """The values of one alignment from the compact output of US-align or TM-align (`-outfmt 1 -m -`).

    The output holds two `>` records (structure 1, then structure 2), each a header line and the aligned sequence, a `# Lali=... RMSD=...` line, `$$$$`, and the rotation matrix. Comment lines such as the circular-permutation note are ignored. The TM-score of a record is normalized by the length of its own structure. Raises ValueError when the text holds no alignment or lacks a piece of it, so a missing value is never read as zero.
    """
    records: List[Tuple[re.Match, str]] = []
    summary: Optional[re.Match] = None
    matrix: Dict[int, Tuple[float, float, float, float]] = {}
    lines = text.splitlines()
    index = 0
    while index < len(lines):
        header = _HEADER.match(lines[index])
        if header:
            if index + 1 >= len(lines):
                raise ValueError("an alignment record has no sequence line")
            records.append((header, lines[index + 1].strip()))
            index += 2
            continue
        found = _SUMMARY.match(lines[index])
        if found:
            summary = found
        row = _MATRIX_ROW.match(lines[index])
        if row:
            matrix[int(row.group("row"))] = (
                _finite(row.group("t"), "translation"), _finite(row.group("u0"), "rotation"),
                _finite(row.group("u1"), "rotation"), _finite(row.group("u2"), "rotation"))
        index += 1

    if not records:
        raise ValueError("the output holds no alignment")
    if len(records) != 2:
        raise ValueError(f"expected the alignment of one pair, found {len(records)} records")
    if summary is None:
        raise ValueError("the output has no 'Lali' summary line")
    if set(matrix) != {0, 1, 2}:
        raise ValueError("the output has no complete rotation matrix")
    (head_1, aligned_1), (head_2, aligned_2) = records
    if len(aligned_1) != len(aligned_2):
        raise ValueError("the two aligned sequences differ in length")

    values: Dict[str, Any] = {
        "length_1": int(head_1.group("length")),
        "length_2": int(head_2.group("length")),
        "aligned_length": int(summary.group("aligned")),
        "rmsd": _finite(summary.group("rmsd"), "RMSD"),
        "tm_score_1": _finite(head_1.group("tm"), "TM-score"),
        "tm_score_2": _finite(head_2.group("tm"), "TM-score"),
        "aligned_1": aligned_1,
        "aligned_2": aligned_2,
    }
    for k, (translation, *rotation) in sorted(matrix.items()):
        values[f"translation_{k + 1}"] = translation
        for j, value in enumerate(rotation):
            values[f"rotation_{k + 1}{j + 1}"] = value
    return values


def align_one(cmd: List[str], run: Callable[[List[str]], Tuple[int, str, str]] = run_command) -> Dict[str, Any]:
    """Run one alignment and return its parsed values; raise RuntimeError, with the aligner's own words, when it gave none.

    The binary exits 0 on an unreadable file and prints only a warning, so the output is parsed whatever the exit code says. A non-zero exit, a crash included, fails the alignment even if some output came out.
    """
    code, out, err = run(cmd)
    detail = _last_line(err) or _last_line(out)
    if code != 0:
        reason = f"killed by signal {-code}" if code < 0 else f"exited {code}"
        raise RuntimeError(f"aligner {reason}" + (f": {detail}" if detail else ""))
    try:
        return parse_alignment(out)
    except ValueError as problem:
        raise RuntimeError(f"{problem}" + (f": {detail}" if detail else "")) from problem


# ── rows ─────────────────────────────────────────────────────────────────────

def empty_row(structure_1: str, structure_2: str, mode: str) -> Dict[str, Any]:
    """A row of an alignment that gave nothing: every number NaN and both aligned sequences empty, the ids kept."""
    row: Dict[str, Any] = {column: math.nan for column in NUMERIC_COLUMNS}
    row.update({"id": row_id(structure_1, structure_2, mode), "structure_1": structure_1,
                "structure_2": structure_2, "mode": mode, "aligned_1": "", "aligned_2": ""})
    return row


def check_row_ids(pairs: List[Tuple[str, str]], modes: List[Dict[str, Any]]) -> None:
    """Refuse two (pair, mode) rows that would share an id.

    `read_pairs` already refuses two pairs with one id. A mode name holds the `_` that separates an id's suffix, so `non_sequential` of the pair (a, b) and `sequential` of the pair (a, b_non) are both `a+b_non_sequential`.
    """
    owner: Dict[str, Tuple[str, str, str]] = {}
    for structure_1, structure_2 in pairs:
        for mode in modes:
            row = (structure_1, structure_2, mode["label"])
            earlier = owner.setdefault(row_id(*row), row)
            if earlier != row:
                raise ValueError(
                    f"the alignments {earlier} and {row} both have the id {row_id(*row)!r}, because a mode name "
                    f"contains '_'; rename the structure so that the row ids are distinct")


def align_pairs(pairs: List[Tuple[str, str]], files: Dict[str, str], declared: set,
                modes: List[Dict[str, Any]], binary: str, step_id: str, container_prefix: str = "",
                run: Callable[[List[str]], Tuple[int, str, str]] = run_command,
                log: Callable[[str], None] = print, report_every: int = 100):
    """Align every pair in every mode; return (alignment rows, local missing rows).

    Each (pair, mode) gives exactly one row, and no two share an id. One that could not be aligned keeps its row with NaN values and adds a kind="failure" row, with the same id, to the missing rows.
    """
    check_row_ids(pairs, modes)
    rows: List[Dict[str, Any]] = []
    missing: List[Dict[str, Any]] = []
    declared = set(declared)
    for number, (structure_1, structure_2) in enumerate(pairs, 1):
        try:
            paths = resolve_structures(structure_1, structure_2, files, declared)
            problem = None
        except PairError as error:
            paths, problem = None, error
        for mode in modes:
            row = empty_row(structure_1, structure_2, mode["label"])
            try:
                if problem is not None:
                    raise problem
                cmd = build_command(binary, paths[0], paths[1], mode["args"], container_prefix)
                row.update(align_one(cmd, run))
            except Exception as error:  # noqa: BLE001 - one pair must not stop the others
                log(f"WARNING: {row['id']} failed: {error}")
                missing.append({"id": row["id"], "removed_by": step_id, "kind": "failure",
                                "cause": str(error)[:200]})
            rows.append(row)
        if number % report_every == 0 or number == len(pairs):
            log(f"  {number}/{len(pairs)} pair(s) aligned, {len(missing)} alignment(s) failed")
    return rows, missing


def alignments_frame(rows: List[Dict[str, Any]]) -> pd.DataFrame:
    """The `alignments` table: the canonical columns, with the length columns kept whole numbers where a failed row leaves a blank."""
    table = pd.DataFrame(rows, columns=ALIGNMENT_COLUMNS)
    for column in INTEGER_COLUMNS:
        table[column] = table[column].astype("Int64")
    return table


def delta_rows(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """One row per pair with the TM-score of the circular-permutation alignment minus the sequential one, for each normalization.

    NaN where either alignment of the pair gave no TM-score.
    """
    by_key = {(r["structure_1"], r["structure_2"], r["mode"]): r for r in rows}
    out = []
    for structure_1, structure_2 in dict.fromkeys((r["structure_1"], r["structure_2"]) for r in rows):
        sequential = by_key.get((structure_1, structure_2, MODE_SEQUENTIAL))
        permuted = by_key.get((structure_1, structure_2, MODE_CIRCULAR_PERMUTATION))
        row = {"id": pair_id(structure_1, structure_2), "structure_1": structure_1,
               "structure_2": structure_2, "delta_tm_1": math.nan, "delta_tm_2": math.nan}
        if sequential is not None and permuted is not None:
            for suffix in ("1", "2"):
                row[f"delta_tm_{suffix}"] = permuted[f"tm_score_{suffix}"] - sequential[f"tm_score_{suffix}"]
        out.append(row)
    return out


# ── wrapper base class ───────────────────────────────────────────────────────

class PairAligner(BaseConfig):
    """Configuration-time half shared by USalign and TMalign.

    A subclass declares its TOOL_NAME, BINARY and `modes()`, takes `structures`, `structure_1` and `structure_2` in its own `__init__` and hands them to `_set_inputs`, and writes its own `get_output_files()` with the table columns spelled out as literals. Everything else, the pair table checks, the config file, the emitted script and the pipe call, is here.
    """

    BINARY = ""
    # True when the tool also writes the per-pair `delta` table (TMalign).
    REPORTS_DELTA = False

    structures_json = Path(lambda self: self.configuration_path("structures.json"))
    config_yaml = Path(lambda self: self.configuration_path("config.yaml"))
    alignments_csv = Path(lambda self: self.table_path("alignments"))
    delta_csv = Path(lambda self: self.table_path("delta"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    local_missing_csv = Path(lambda self: self.execution_path("local_missing.csv"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_pair_alignment.py"))

    def modes(self) -> List[Dict[str, Any]]:
        """The alignments to run per pair, in order: [{"label": mode name, "args": [flags]}]."""
        raise NotImplementedError

    def _set_inputs(self, structures, structure_1, structure_2):
        self.structures = structures  # raw handle kept for missing-propagation
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")
        self.structure_1 = resolve_table_reference(structure_1, "structure_1")
        self.structure_2 = resolve_table_reference(structure_2, "structure_2")

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        for name in ("structure_1", "structure_2"):
            if not isinstance(getattr(self, name), TableReference):
                raise ValueError(
                    f"{name} must reference a column of the pair table, e.g. pairs.tables.data.<column> "
                    f"or (pairs.tables.data, '<column>'), got {getattr(self, name)!r}")
        if self.structure_1.path != self.structure_2.path:
            raise ValueError(
                "structure_1 and structure_2 must be columns of the same pair table, got "
                f"{self.structure_1.path} and {self.structure_2.path}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"PAIRS: {self.structure_1.column} -> {self.structure_2.column} (rotates structure_1 onto structure_2)")
        lines.append(f"ALIGNMENTS: {', '.join(mode['label'] for mode in self.modes())}")
        return lines

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        self._write_config_yaml()
        delta_flag = f'    --delta-csv "{self.delta_csv}" \\\n' if self.REPORTS_DELTA else ""
        script = "#!/bin/bash\n"
        script += f"# {self.TOOL_NAME} pairwise structure alignment script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += f"""echo "Aligning the pairs of the pair table with {self.TOOL_NAME}"
python "{self.helper_py}" \\
    --config-yaml "{self.config_yaml}" \\
    --structures-json "{self.structures_json}" \\
    --alignments-csv "{self.alignments_csv}" \\
{delta_flag}    --local-missing-csv "{self.local_missing_csv}" \\
    --container-prefix "{self.container_prefix()}"
"""
        script += self.generate_missing_propagation(
            self.structures, local_missing=self.local_missing_csv, missing_csv=self.missing_csv
        )
        script += self.generate_completion_check_footer()
        return script

    def _write_config_yaml(self):
        """Dump the parameters to configuration/config.yaml for the pipe to read."""
        import yaml
        with open(self.config_yaml, "w") as f:
            yaml.safe_dump({
                "binary": self.BINARY,
                "modes": self.modes(),
                "structure_1": str(self.structure_1),
                "structure_2": str(self.structure_2),
            }, f, sort_keys=False)
