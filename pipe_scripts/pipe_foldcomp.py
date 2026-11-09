#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Foldcomp runner. Decompresses every requested entry from a Foldcomp database in one `foldcomp decompress --id-list` pass, cuts each requested residue range out of the decompressed model, renumbers it from 1 and writes one PDB per domain with the mean pLDDT read from its B-factor column."""

import argparse
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, replace
from typing import Dict, List, Optional, Tuple

import pandas as pd
import yaml

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.biopipelines_io import (  # noqa: E402
    TableReference,
    load_table,
    step_id_from_table_path,
    container_argv_prefix,
)
from biopipelines.combinatorics import generate_multiplied_ids  # noqa: E402
from biopipelines.foldcomp import DOMAIN_COLUMNS  # noqa: E402
from biopipelines.pdb_parser import (  # noqa: E402
    field_b_factor,
    field_chain,
    field_res_name,
    field_res_seq,
    residue_shift,
)


ENTRY_NAME = re.compile(r"^AF-[A-Z0-9]+-F[0-9]+-model_v[0-9]+$")
ENTRY_NAME_ANYWHERE = re.compile(r"AF-[A-Z0-9]+-F[0-9]+-model_v[0-9]+")
UNIPROT_ACCESSION = re.compile(r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$")
INTEGER_COLUMNS = ["start", "end", "length", "offset", "entry_length"]


@dataclass
class Request:
    id: str
    entry: str
    name: str = ""
    start: Optional[int] = None
    end: Optional[int] = None
    error: str = ""


@dataclass
class Cut:
    text: str
    length: int
    offset: int
    entry_length: int
    plddt_mean: float
    clipped: bool


def parse_entry(entry: str, model_version: int) -> str:
    """Foldcomp database name of an entry given as a UniProt accession or as a full AFDB name."""
    if ENTRY_NAME.match(entry):
        return entry
    if UNIPROT_ACCESSION.match(entry):
        return f"AF-{entry}-F1-model_v{model_version}"
    raise ValueError(
        f"entry {entry!r} is neither a UniProt accession nor an AFDB entry name (AF-<accession>-F<n>-model_v<m>)")


def parse_position(value, label: str) -> int:
    """A residue position from a table cell, which may arrive as '12', 12 or 12.0."""
    if value is None or pd.isna(value) or str(value).strip() == "":
        raise ValueError(f"{label} is empty")
    try:
        number = float(str(value).strip())
    except ValueError:
        raise ValueError(f"{label} {value!r} is not an integer") from None
    if not number.is_integer():
        raise ValueError(f"{label} {value!r} is not an integer")
    return int(number)


def build_requests(entries, starts, ends, model_version: int) -> List[Request]:
    """One request per distinct domain id, in table order; a row that cannot be parsed keeps its place with its error."""
    requests: List[Request] = []
    seen = set()
    for index, (raw_entry, raw_start, raw_end) in enumerate(zip(entries, starts, ends)):
        entry = "" if pd.isna(raw_entry) else str(raw_entry).strip()
        request = Request(id="", entry=entry)
        try:
            if not entry:
                raise ValueError("entry is empty")
            request.name = parse_entry(entry, model_version)
            request.start = parse_position(raw_start, "start")
            request.end = parse_position(raw_end, "end")
            if request.start < 1:
                raise ValueError(f"start must be at least 1, got {request.start}")
            if request.start > request.end:
                raise ValueError(f"start {request.start} is after end {request.end}")
        except ValueError as e:
            request.error = str(e)
        start_text = request.start if request.start is not None else str(raw_start).strip()
        end_text = request.end if request.end is not None else str(raw_end).strip()
        # A row with no entry has no parent to multiply, and numbered_id would put an underscore in `row<n>`, so that id stays as it is.
        request.id = generate_multiplied_ids([entry], [f"{start_text}-{end_text}"])[0][0] if entry else f"row{index}"
        if request.id in seen:
            print(f"WARNING: {request.id} appears more than once in the table; its first row is used", file=sys.stderr)
            continue
        seen.add(request.id)
        requests.append(request)
    return requests


def check_database(database: str) -> None:
    """Foldcomp reads a prefix as a database only when <prefix>.dbtype exists, so every part is checked up front."""
    absent = [database + suffix for suffix in ("", ".dbtype", ".index", ".lookup") if not os.path.exists(database + suffix)]
    if absent:
        raise RuntimeError(f"not a Foldcomp database: {', '.join(absent)} not found")


def build_command(ids_file: str, database: str, out_dir: str, threads: int, container_prefix: str = "") -> List[str]:
    """The one `foldcomp decompress` argv that extracts every listed name from the database."""
    cmd = ["foldcomp", "decompress", "-t", str(threads), "--id-list", ids_file, database, out_dir]
    return container_argv_prefix(container_prefix) + cmd


def run_foldcomp(cmd: List[str]) -> str:
    """Run Foldcomp and return its stderr, which carries the not-found warnings; raise on a non-zero exit."""
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        tail = (result.stderr.strip() or result.stdout.strip())[-300:]
        raise RuntimeError(f"foldcomp exited {result.returncode}: {tail}")
    return result.stderr


def parse_not_found(stderr: str) -> set:
    """Names Foldcomp reported as absent from the database ('[Warning] <name> not found in database.')."""
    return set(re.findall(r"^\[Warning\] (\S+) not found in database\.$", stderr, re.MULTILINE))


def parse_mentions(stderr: str) -> Dict[str, str]:
    """The first stderr line that mentions each AFDB entry name, found in one pass."""
    mentions: Dict[str, str] = {}
    for line in stderr.splitlines():
        for name in ENTRY_NAME_ANYWHERE.findall(line):
            mentions.setdefault(name, line.strip())
    return mentions


def explain_absence(name: str, database: str, not_found: set, mentions: Dict[str, str]) -> str:
    """Why a requested name produced no decompressed file."""
    if name in not_found:
        return f"entry {name} is not in the Foldcomp database {os.path.basename(database)}"
    if name in mentions:
        return f"foldcomp wrote no structure for {name}: {mentions[name]}"
    return f"foldcomp wrote no structure for {name}"


def read_atoms(path: str) -> List[str]:
    with open(path) as f:
        return [line.rstrip("\n") for line in f if line.startswith("ATOM")]


def cut_domain(atoms: List[str], start: int, end: int, title: str) -> Cut:
    """Residues start..end of a one-chain model, clipped to what the model holds, renumbered from 1.

    Original number = new number + offset. Atom serials restart at 1. Raises when the range misses the model.
    """
    chains = {field_chain(line) for line in atoms}
    if len(chains) != 1:
        raise ValueError(f"the model holds {len(chains)} chains; a residue range needs exactly one")
    numbers = sorted({int(field_res_seq(line)) for line in atoms})
    first, last = numbers[0], numbers[-1]
    if start > last or end < first:
        raise ValueError(f"range {start}-{end} lies outside the entry (residues {first}-{last}, length {len(numbers)})")
    low, high = max(start, first), min(end, last)
    chosen = [line for line in atoms if low <= int(field_res_seq(line)) <= high]
    if not chosen:
        raise ValueError(f"no residue of the entry lies in range {start}-{end}")
    offset = low - 1
    lines = [f"TITLE     {title}"]
    plddt: Dict[int, float] = {}
    for serial, line in enumerate(chosen, 1):
        number = int(field_res_seq(line))
        plddt.setdefault(number, field_b_factor(line))
        # The serial and the residue number are rewritten in place, so these columns stay explicit.
        shift = residue_shift(line)
        lines.append(f"{line[:6]}{serial:5d}{line[11:22 + shift]}{number - offset:4d}{line[26 + shift:]}")
    final = chosen[-1]
    lines.append(f"TER   {len(chosen) + 1:5d}      {field_res_name(final):>3} "
                 f"{field_chain(final):1}{int(field_res_seq(final)) - offset:4d}")
    return Cut(
        text="\n".join(lines) + "\n",
        length=len(plddt),
        offset=offset,
        entry_length=len(numbers),
        plddt_mean=round(sum(plddt.values()) / len(plddt), 2),
        clipped=(low, high) != (start, end),
    )


def extract_entries(names: List[str], database: str, scratch_dir: str, threads: int,
                    container_prefix: str) -> Tuple[str, str]:
    """Decompress every name in one pass; return (folder of <name>.pdb files, Foldcomp's stderr)."""
    check_database(database)
    ids_file = os.path.join(scratch_dir, "foldcomp_ids.txt")
    out_dir = os.path.join(scratch_dir, "foldcomp_out")
    shutil.rmtree(out_dir, ignore_errors=True)
    with open(ids_file, "w") as f:
        f.write("".join(f"{name}\n" for name in names))
    log_path = os.path.join(scratch_dir, "foldcomp.log")
    try:
        stderr = run_foldcomp(build_command(ids_file, database, out_dir, threads, container_prefix))
    except RuntimeError as e:
        with open(log_path, "w") as f:
            f.write(str(e) + "\n")
        raise
    with open(log_path, "w") as f:
        f.write(stderr)
    return out_dir, stderr


def cut_entries(by_name: Dict[str, List[Request]], out_dir: str, structures_dir: str, database: str,
                stderr: str, database_error: str) -> Dict[str, Tuple[Optional[Cut], str]]:
    """{domain id: (Cut or None, failure cause)}; each domain is written the moment its entry is read, one entry in memory at a time."""
    not_found, mentions = parse_not_found(stderr), parse_mentions(stderr)
    outcomes: Dict[str, Tuple[Optional[Cut], str]] = {}
    for name, group in by_name.items():
        path = os.path.join(out_dir, f"{name}.pdb") if out_dir else ""
        atoms, reason = None, ""
        if database_error:
            reason = database_error
        elif os.path.exists(path):
            atoms = read_atoms(path)
        else:
            reason = explain_absence(name, database, not_found, mentions)
        for request in group:
            if atoms is None:
                outcomes[request.id] = (None, reason)
                continue
            try:
                cut = cut_domain(atoms, request.start, request.end, f"{name} residues {request.start}-{request.end}")
            except ValueError as e:
                outcomes[request.id] = (None, str(e))
                continue
            with open(os.path.join(structures_dir, f"{request.id}.pdb"), "w") as f:
                f.write(cut.text)
            outcomes[request.id] = (replace(cut, text=""), "")
        if atoms is not None:
            os.remove(path)
    return outcomes


def domain_row(request: Request, **values) -> dict:
    row = {"id": request.id, "entry": request.entry, "database_name": request.name or None,
           "start": request.start, "end": request.end, "length": None, "offset": None,
           "entry_length": None, "plddt_mean": None, "found": False, "clipped": None}
    row.update(values)
    return row


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config-yaml", required=True)
    p.add_argument("--structures-dir", required=True)
    p.add_argument("--structures-map-csv", required=True)
    p.add_argument("--domains-csv", required=True)
    p.add_argument("--scratch-dir", required=True)
    p.add_argument("--local-missing-csv", required=True)
    p.add_argument("--container-prefix", default="")
    args = p.parse_args()

    with open(args.config_yaml) as f:
        cfg = yaml.safe_load(f)
    step_id = step_id_from_table_path(args.local_missing_csv)
    database = cfg["database"]

    table, entry_column = load_table(cfg["entries"])
    start_column = TableReference.from_string(cfg["start"]).column
    end_column = TableReference.from_string(cfg["end"]).column
    for column in (entry_column, start_column, end_column):
        if column not in table.columns:
            print(f"ERROR: column {column!r} not in the entries table (columns: {list(table.columns)})", file=sys.stderr)
            sys.exit(1)
    requests = build_requests(table[entry_column].tolist(), table[start_column].tolist(),
                              table[end_column].tolist(), int(cfg["model_version"]))
    if not requests:
        print("ERROR: the entries table has no rows", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.structures_dir, exist_ok=True)
    os.makedirs(args.scratch_dir, exist_ok=True)
    names = list(dict.fromkeys(r.name for r in requests if not r.error))
    out_dir, stderr, database_error = "", "", ""
    if names:
        print(f"Extracting {len(names)} entr{'y' if len(names) == 1 else 'ies'} for {len(requests)} domain(s) from {database}")
        try:
            out_dir, stderr = extract_entries(names, database, args.scratch_dir, int(cfg["threads"]), args.container_prefix)
        except Exception as e:
            database_error = str(e)
            print(f"ERROR: {database_error}", file=sys.stderr)

    by_name: Dict[str, List[Request]] = {}
    for request in requests:
        if not request.error:
            by_name.setdefault(request.name, []).append(request)
    outcomes = cut_entries(by_name, out_dir, args.structures_dir, database, stderr, database_error)

    domain_rows, map_rows, local_missing = [], [], []
    for request in requests:
        cut, cause = outcomes.get(request.id, (None, request.error))
        if cut is None:
            print(f"WARNING: {request.id} failed: {cause}", file=sys.stderr)
            domain_rows.append(domain_row(request))
            local_missing.append({"id": request.id, "removed_by": step_id, "kind": "failure", "cause": cause[:200]})
            continue
        map_rows.append({"id": request.id, "file": os.path.join(args.structures_dir, f"{request.id}.pdb")})
        domain_rows.append(domain_row(request, length=cut.length, offset=cut.offset, entry_length=cut.entry_length,
                                      plddt_mean=cut.plddt_mean, found=True, clipped=cut.clipped))

    for path in (args.structures_map_csv, args.domains_csv, args.local_missing_csv):
        os.makedirs(os.path.dirname(path), exist_ok=True)
    domains = pd.DataFrame(domain_rows, columns=DOMAIN_COLUMNS)
    for column in INTEGER_COLUMNS:
        domains[column] = domains[column].astype("Int64")
    domains["clipped"] = domains["clipped"].astype("boolean")
    domains.to_csv(args.domains_csv, index=False)
    pd.DataFrame(map_rows, columns=["id", "file"]).to_csv(args.structures_map_csv, index=False)
    pd.DataFrame(local_missing, columns=["id", "removed_by", "kind", "cause"]).to_csv(args.local_missing_csv, index=False)
    if out_dir:
        shutil.rmtree(out_dir, ignore_errors=True)
    print(f"Domains: {args.domains_csv} ({len(domains)} rows, {len(map_rows)} extracted, {int(domains['clipped'].fillna(False).sum())} clipped)")
    print(f"Local failures: {len(local_missing)}")

    if len(local_missing) == len(requests):
        print(f"ERROR: all {len(requests)} attempted domain(s) failed", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
