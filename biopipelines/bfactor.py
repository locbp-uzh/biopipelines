# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""BFactor tool: per-residue B-factor, summarized over named selections.

The B-factor column carries a structure predictor's per-residue confidence - Boltz2, AlphaFold and ESMFold all write pLDDT there on a 0-100 scale, while an experimental structure carries a real temperature factor in the same field. Every tool that reports confidence reports one number for the whole structure, so the per-residue vector already on disk has no reader. This is that reader.
"""

import os
import re
from typing import Dict, List, Any, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .data_containers import resolve_table_reference
    from .file_paths import Path
    from .datastream import DataStream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from data_containers import resolve_table_reference
    from file_paths import Path
    from datastream import DataStream

SUMMARY_BASE_COLUMNS = ("id", "n_residues", "all_mean", "all_sd")
SUMMARY_PER_SELECTION = ("mean", "sd", "min", "n", "delta")


class BFactor(BaseConfig):
    """
    BFactor: reads the B-factor column per residue and summarizes it over named selections.

    Values are reported as they appear in the file, never rescaled. Boltz2, AlphaFold and
    ESMFold write pLDDT in [0,100]; an experimental structure carries a temperature factor in
    Å², where higher means more disordered rather than more confident. Rescaling on a guessed
    range would corrupt one of them, so the caller keeps the units it supplied.

    One value per residue: the CA atom's B-factor when the residue has a CA, otherwise the
    mean over that residue's atoms.

    Inputs:
        structures: PDB/mmCIF structures.
        selections: optional {name: selection} mapping. Each selection is a framework
            selection string ("A75-77+A274"; unqualified "75-77" when only one chain has those
            residues; if several do, that id is recorded as a failure naming them) or a
            table column reference (``tool.tables.structures.designed``) resolved per
            structure id. One group of summary columns is emitted per name. A single-row
            table broadcasts to every structure.

    Outputs:
        Streams:
            bfactors: per-residue resi-csv (one <id>.csv per input) with columns
                      id | chain | resi | icode | bfactor. Consumable by Selection and Consensus,
                      e.g. ``Selection.add(bf.streams.bfactors, include="bfactor>=70")``.
        Tables:
            summary: id | n_residues | all_mean | all_sd, plus <name>_mean | <name>_sd |
                     <name>_min | <name>_n | <name>_delta for each named selection.
                     ``<name>_delta`` is the selection's mean minus the whole structure's -
                     negative means the region is less certain than the model around it.
            missing: id | removed_by | kind | cause
    """

    TOOL_NAME = "BFactor"
    TOOL_VERSION = "1.2"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== BFactor ==="
echo "Uses biopipelines environment (no additional installation needed)."
touch "$INSTALL_SUCCESS"
echo "=== BFactor ready ==="
"""

    structures_json = Path(lambda self: self.configuration_path("structures.json"))
    config_yaml = Path(lambda self: self.configuration_path("config.yaml"))
    bfactors_map = Path(lambda self: self.stream_map_path("bfactors"))
    summary_csv = Path(lambda self: self.table_path("summary"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    local_missing_csv = Path(lambda self: self.execution_path("local_missing.csv"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_bfactor.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 selections: Dict[str, Any] = None,
                 **kwargs):
        self.structures = structures
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")
        self.selections = {
            name: resolve_table_reference(value, f"selections[{name!r}]")
            for name, value in (selections or {}).items()
        }
        super().__init__(**kwargs)

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        for name in self.selections:
            # The name becomes a column prefix, so it has to survive a CSV header and a query.
            if not name or not name.replace("_", "").isalnum():
                raise ValueError(
                    f"selection name {name!r} must be alphanumeric with underscores; "
                    f"it is used as a prefix for the summary columns")
        clashes = set(self.selections) & {"all", "id", "n_residues"}
        if clashes:
            raise ValueError(
                f"selection name(s) {sorted(clashes)} collide with the summary table's own "
                f"columns; rename them")
        # A literal is checked here rather than on the GPU node, where a typo read as an empty selection.
        span = re.compile(r"^[A-Za-z]*\d+(-\d+)?$")
        for name, value in self.selections.items():
            if not isinstance(value, str) or value.startswith("TABLE_REFERENCE:"):
                continue
            bad = [p.strip() for p in value.replace(",", "+").split("+")
                   if p.strip() and not span.match(p.strip())]
            if bad:
                raise ValueError(f"selections[{name!r}]: not a residue span: {bad[0]!r} "
                                 f"(expected e.g. A75-77+A274)")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"SELECTIONS: {', '.join(self.selections) if self.selections else 'whole structure only'}")
        return lines

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        self._write_config_yaml()
        script = "#!/bin/bash\n"
        script += "# BFactor per-residue extraction script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += f"""echo "Running BFactor on {len(self.structures_stream)} structure(s)"
python "{self.helper_py}" \\
    --config-yaml "{self.config_yaml}" \\
    --structures-json "{self.structures_json}" \\
    --bfactors-dir "{self.stream_folder('bfactors')}" \\
    --bfactors-map-csv "{self.bfactors_map}" \\
    --summary-csv "{self.summary_csv}" \\
    --local-missing-csv "{self.local_missing_csv}"
"""
        script += self.generate_missing_propagation(
            self.structures, local_missing=self.local_missing_csv, missing_csv=self.missing_csv
        )
        script += self.generate_completion_check_footer()
        return script

    def _write_config_yaml(self):
        """A TableReference serializes to TABLE_REFERENCE:path:column, which the pipe resolves per id."""
        import yaml
        with open(self.config_yaml, "w") as f:
            yaml.safe_dump(
                {"selections": {name: str(value) for name, value in self.selections.items()}},
                f, sort_keys=False)

    def _summary_columns(self) -> List[str]:
        columns = list(SUMMARY_BASE_COLUMNS)
        for name in self.selections:
            columns += [f"{name}_{suffix}" for suffix in SUMMARY_PER_SELECTION]
        return columns

    def get_output_files(self) -> Dict[str, Any]:
        bfactors_stream = DataStream(
            name="bfactors",
            ids=self.structures_stream.ids,
            files=[self.stream_path("bfactors", "<id>.csv")],
            map_table=self.bfactors_map,
            format="resi-csv",
        )
        tables = {
            "summary": TableInfo(
                name="summary",
                path=self.summary_csv,
                columns=self._summary_columns(),
                description="Per-structure B-factor statistics, overall and per named selection",
            ),
            "missing": self.missing_table_info(self.missing_csv),
        }
        return {
            "bfactors": bfactors_stream,
            "tables": tables,
            "output_folder": self.output_folder,
        }
