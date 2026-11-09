# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Foldcomp tool: extract domains from a Foldcomp-compressed AlphaFold Database.

Takes a table of AFDB entries with residue ranges, decompresses every entry from a local Foldcomp database in one `foldcomp decompress --id-list` pass, cuts each range out of the decompressed model and emits one PDB per domain plus a `domains` table.

References: Kim et al., Bioinformatics, doi:10.1093/bioinformatics/btad153 (Foldcomp). https://github.com/steineggerlab/foldcomp
"""

import os
from typing import Dict, List, Any

try:
    from .base_config import BaseConfig, TableInfo, _validate_freeform_string, resolve_table_reference
    from .biopipelines_io import TableReference
    from .file_paths import Path
    from .datastream import DataStream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, TableInfo, _validate_freeform_string, resolve_table_reference
    from biopipelines_io import TableReference
    from file_paths import Path
    from datastream import DataStream


DOWNLOAD_URL = "https://opendata.mmseqs.org/foldcomp"
SWISSPROT_DATABASE = "afdb_swissprot_v4"
# Files of a database that `foldcomp decompress` reads; the `.source` file is not needed.
DATABASE_SUFFIXES = ("", ".index", ".dbtype", ".lookup")
DOMAIN_COLUMNS = ["id", "entry", "database_name", "start", "end", "length", "offset",
                  "entry_length", "plddt_mean", "found", "clipped"]
# Any id, including the underscore-separated range suffix: the ids come from the entries table at run time.
DOMAIN_ID_PATTERN = "[<?>][_<?>]"


class Foldcomp(BaseConfig):
    """
    Foldcomp: extract domains of the AlphaFold Database from a compressed Foldcomp database.

    Takes the AFDB entries of a table (a UniProt accession or a full AFDB entry name per row) with the start and end residue of each domain, and returns one structure per domain. The residue ranges are cut by this tool on the decompressed model: Foldcomp itself decompresses whole entries only.

    ENTRIES AND RANGES
    The AFDB clusters of Barrio-Hernandez et al. (Nature 2023, doi:10.1038/s41586-023-06510-w, cluster.foldseek.com) identify every member by its bare UniProt accession, for example `A0A009F5K6`, and cover whole chains: their membership files carry no residue ranges. The domain boundaries of that release are in `4-domain-clustering.zip` (`Structures_from_groups.csv`: `structure_ID` = accession + `_<n>`, `structure_start`, `structure_end`, `uniprot`). The release does not state the convention of the positions; they are read here as 1-based, inclusive positions on the UniProt sequence (the smallest start in the file is 1). Load such a file with `Table` and pass its columns: `uniprot` as `entries`, `structure_start` as `start`, `structure_end` as `end`. A bare accession is expanded to `AF-<accession>-F1-model_v<model_version>`, the name Foldcomp stores; a full AFDB name (`AF-P69905-F1-model_v4`) is used as given, which is how a fragment `F2` of a long protein is requested.

    DATABASE
    `database` is the path prefix of a Foldcomp database, the common stem of `<prefix>`, `<prefix>.index`, `<prefix>.dbtype` and `<prefix>.lookup`. Published databases (https://opendata.mmseqs.org/foldcomp/): `afdb_swissprot_v4` (3 GB), `afdb_rep_dark_v4` (2.3 GB), `afdb_rep_v4` (10 GB, the AFDB cluster representatives) and `afdb_uniprot_v4` (1.1 TB, the full AFDB, so the one to use for arbitrary cluster members). The tool never downloads them while running. `Foldcomp.install(swissprot_dir=...)` fetches `afdb_swissprot_v4` once; for the large ones download the files yourself and pass the prefix. Foldcomp reads the whole `.lookup` file into memory before it looks up a name (9 GB on disk for `afdb_uniprot_v4`), so request memory with Resources accordingly. An entry absent from the chosen database is a `failure` row in `missing`. The path is written into the generated script, so it may not contain `"`, a backtick, `$` or a backslash.

    NUMBERING
    Each domain is renumbered from 1 and keeps its chain letter; `offset` is recorded so that original residue number = residue number in the file + `offset`. Atom serials restart at 1. A range that runs past the model is clipped to the model and flagged in `clipped`.

    Inputs:
        entries: column of AFDB entries, a TableReference (`tool.tables.name.column`) or a `(TableInfo | path, column)` tuple.
        start: column of the first residue of each domain, from the same table as `entries`.
        end: column of the last residue of each domain (inclusive), from the same table.
        database: path prefix of a local Foldcomp database.
        model_version: the N of `model_vN` in the database's entry names, used to expand bare accessions (default 4, the version of every published AFDB database).
        threads: CPU threads Foldcomp uses (`-t`, default 1). Match the CPUs requested with Resources.

    Outputs:
        Streams:
            structures: one PDB per domain, ids `<entry>_<start>-<end>` as written in the table (for example `A0A009F5K6_12-150`), built by the shared multiplied-id rule with the range as the suffix, and `row<n>` for a row whose entry is empty; lazy, because the table is read at run time.
        Tables:
            domains: id | entry | database_name | start | end | length | offset | entry_length | plddt_mean | found | clipped. One row per distinct domain id of the table. `start`/`end` are the requested range, `length` the residues written, `offset` the shift from the file's numbering to the original one, `entry_length` the residues of the database model, `plddt_mean` the mean over the domain's residues of the B-factor column, `found` whether the domain was written and `clipped` whether the range was cut to the model. An id that failed carries NaN for the values it could not obtain.
            missing: id | removed_by | kind | cause. `kind="failure"` for an entry that is not in the database, a range outside the model (the cause gives the model's residues and length), an unparsable entry or position.

    Usage::

        members = Table("Structures_from_groups.csv", table_name="domains")
        domains = Foldcomp(entries=members.tables.domains.uniprot,
                           start=members.tables.domains.structure_start,
                           end=members.tables.domains.structure_end,
                           database="/data/foldcomp/afdb_rep_v4",
                           threads=4)
    """

    TOOL_NAME = "Foldcomp"
    TOOL_VERSION = "1.3"
    ENV_NAME = "foldcomp"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, swissprot_dir=None, **kwargs):
        """Install the foldcomp env; with `swissprot_dir`, also download `afdb_swissprot_v4` (about 3 GB) into it."""
        _validate_freeform_string("swissprot_dir", swissprot_dir)
        env = cls._install_env(env_manager)
        biopipelines = folders.get("biopipelines", "")
        env_check = cls._env_exists_check(env, env_manager)
        remove_block = cls._env_remove_block(env, env_manager) if force_reinstall else ""
        env_block = cls._env_install_block(env, env_manager, biopipelines)
        database_block = cls._swissprot_download_block(swissprot_dir) if swissprot_dir else ""
        return f"""echo "=== Installing Foldcomp ==="
{remove_block}
if {env_check}; then
    echo "foldcomp environment already present, skipping creation. Use force_reinstall=True to reinstall."
else
{env_block}
fi
{database_block}
if {cls._env_run(env, env_manager)}foldcomp --version >/dev/null 2>&1 && {cls._env_run(env, env_manager)}python -c "import pandas, yaml" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== Foldcomp installation complete ==="
else
    echo "ERROR: Foldcomp verification failed (foldcomp binary not on PATH, or pandas/yaml not importable)"
    exit 1
fi
"""

    @staticmethod
    def _swissprot_download_block(directory: str) -> str:
        """Bash that fetches the four files of `afdb_swissprot_v4`, resuming and skipping a file whose size already matches upstream."""
        suffixes = " ".join(f'"{s}"' for s in DATABASE_SUFFIXES)
        return f"""
echo "Downloading {SWISSPROT_DATABASE} into {directory}"
mkdir -p "{directory}"
for SUFFIX in {suffixes}; do
    URL="{DOWNLOAD_URL}/{SWISSPROT_DATABASE}$SUFFIX"
    FILE="{directory}/{SWISSPROT_DATABASE}$SUFFIX"
    REMOTE=$(curl -fsSIL "$URL" | tr -d '\\r' | awk 'tolower($1)=="content-length:" {{size=$2}} END {{print size}}')
    LOCAL=0
    [ -f "$FILE" ] && LOCAL=$(wc -c < "$FILE" | tr -d ' ')
    if [ -n "$REMOTE" ] && [ "$LOCAL" = "$REMOTE" ]; then
        echo "Already complete: $FILE"
    else
        curl -fL -C - -o "$FILE" "$URL" || exit 1
    fi
done
echo "Foldcomp database ready: pass database=\\"{directory}/{SWISSPROT_DATABASE}\\""
"""

    config_yaml = Path(lambda self: self.configuration_path("config.yaml"))
    structures_map = Path(lambda self: self.stream_map_path("structures"))
    domains_csv = Path(lambda self: self.table_path("domains"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    local_missing_csv = Path(lambda self: self.execution_path("local_missing.csv"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_foldcomp.py"))

    def __init__(self,
                 entries,
                 start,
                 end,
                 database: str,
                 model_version: int = 4,
                 threads: int = 1,
                 **kwargs):
        self.entries_ref = resolve_table_reference(entries, "entries")
        self.start_ref = resolve_table_reference(start, "start")
        self.end_ref = resolve_table_reference(end, "end")
        self.database = database
        self.model_version = int(model_version)
        self.threads = int(threads)
        super().__init__(**kwargs)

    def validate_params(self):
        references = (("entries", self.entries_ref), ("start", self.start_ref), ("end", self.end_ref))
        for label, reference in references:
            if not isinstance(reference, TableReference):
                raise ValueError(
                    f"{label} must be a column of a table, e.g. members.tables.domains.{label}; "
                    f"got {type(reference).__name__}")
        if len({str(reference.path) for _, reference in references}) != 1:
            raise ValueError("entries, start and end must be columns of the same table")
        if len({reference.column for _, reference in references}) != 3:
            raise ValueError("entries, start and end must be three different columns")
        if not isinstance(self.database, str) or not self.database.strip():
            raise ValueError("database must be the path prefix of a Foldcomp database")
        _validate_freeform_string("database", self.database)
        if self.model_version < 1:
            raise ValueError(f"model_version must be at least 1, got {self.model_version}")
        if self.threads < 1:
            raise ValueError("threads must be at least 1")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"DATABASE: {self.database}")
        lines.append(f"ENTRIES: {self.entries_ref.column} [{self.start_ref.column}-{self.end_ref.column}] of {self.entries_ref.path}")
        lines.append(f"MODEL VERSION: {self.model_version}")
        lines.append(f"THREADS: {self.threads}")
        return lines

    def generate_script(self, script_path: str) -> str:
        self._write_config_yaml()
        script = "#!/bin/bash\n"
        script += "# Foldcomp domain extraction script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += f"""echo "Extracting domains from {os.path.basename(self.database)}"
python "{self.helper_py}" \\
    --config-yaml "{self.config_yaml}" \\
    --structures-dir "{self.stream_folder('structures')}" \\
    --structures-map-csv "{self.structures_map}" \\
    --domains-csv "{self.domains_csv}" \\
    --scratch-dir "{self.execution_path()}" \\
    --local-missing-csv "{self.local_missing_csv}" \\
    --container-prefix "{self.container_prefix()}"
"""
        script += self.generate_missing_propagation(
            local_missing=self.local_missing_csv, missing_csv=self.missing_csv
        )
        script += self.generate_completion_check_footer()
        return script

    def _write_config_yaml(self):
        """Dump the parameters to configuration/config.yaml for the pipe to read."""
        import yaml
        with open(self.config_yaml, "w") as f:
            yaml.safe_dump({
                "entries": str(self.entries_ref),
                "start": str(self.start_ref),
                "end": str(self.end_ref),
                "database": self.database,
                "model_version": self.model_version,
                "threads": self.threads,
            }, f, sort_keys=False)

    def get_output_files(self) -> Dict[str, Any]:
        structures = DataStream(
            name="structures",
            ids=[DOMAIN_ID_PATTERN],
            files=[self.stream_path("structures", "<id>.pdb")],
            map_table=self.structures_map,
            format="pdb",
        )
        tables = {
            "domains": TableInfo(
                name="domains",
                path=self.domains_csv,
                columns=list(DOMAIN_COLUMNS),
                description="Requested range, written length, numbering offset and mean pLDDT of every domain",
            ),
            "missing": self.missing_table_info(self.missing_csv),
        }
        return {
            "structures": structures,
            "tables": tables,
            "output_folder": self.output_folder,
        }
