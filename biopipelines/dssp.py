# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""DSSP tool: per-residue secondary-structure assignment.

For each input PDB, runs DSSP (`mkdssp`) and reports:
  * a per-residue table with the 8-state DSSP code,
  * a per-structure summary with helix/sheet/coil fractions.

Reference: https://github.com/PDB-REDO/dssp
"""

import os
from typing import Dict, List, Any, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream


class DSSP(BaseConfig):
    """
    DSSP: secondary-structure assignment from 3D coordinates.

    Inputs:
        structures: PDB structures.

    Outputs:
        Streams:
            dssp: one <id>.dssp file per input (raw mkdssp output, classic format).
            ss:   per-residue resi-csv (one <id>.csv per input) with columns
                  id | chain | resi | resn | ss_code | ss_simple | acc | rsa |
                  is_helix | is_sheet | is_coil. ``acc`` is absolute solvent
                  accessibility (Å²) and ``rsa`` is relative accessibility in
                  [0,1] (acc / Tien max-ACC per restype). The numeric acc/rsa
                  and 0/1 is_* columns let Selection threshold on accessibility
                  or SS class, e.g.
                  ``Selection.add(dssp.streams.ss, include="rsa>=0.25")`` or
                  ``Selection.add(dssp.streams.ss, include="is_helix==1")``.
        Tables:
            secondary_structure: id | chain | resnum | resi | resname | resn | ss_code | ss_simple
                                 (resi duplicates resnum, resn duplicates resname;
                                 PyMOL-style aliases for downstream compatibility)
            summary:             id | n_residues | helix_frac | sheet_frac | coil_frac
    """

    TOOL_NAME = "DSSP"
    TOOL_VERSION = "1.1"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        biopipelines = folders.get("biopipelines", "")
        env_check = cls._env_exists_check("dssp", env_manager)
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_check}; then
    echo "DSSP already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block("dssp", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("dssp", env_manager, biopipelines)
        return f"""echo "=== Installing DSSP ==="
{skip}{remove_block}
{env_block}

if {env_manager} run -n dssp mkdssp --help >/dev/null 2>&1 || {env_manager} run -n dssp dssp --help >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== DSSP installation complete ==="
else
    echo "ERROR: DSSP verification failed (mkdssp/dssp binary not on PATH)"
    exit 1
fi
"""

    structures_json = Path(lambda self: self.configuration_path("structures.json"))
    dssp_map = Path(lambda self: self.stream_map_path("dssp"))
    ss_map = Path(lambda self: self.stream_map_path("ss"))
    secondary_csv = Path(lambda self: self.table_path("secondary_structure"))
    summary_csv = Path(lambda self: self.table_path("summary"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_dssp.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 **kwargs):
        self.structures = structures
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")
        super().__init__(**kwargs)

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"STRUCTURES: {len(self.structures_stream)} files")
        return lines

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        upstream_missing_path = self._get_upstream_missing_table_path(self.structures)
        upstream_missing_flag = f' \\\n    --upstream-missing "{upstream_missing_path}"' if upstream_missing_path else ""
        script = "#!/bin/bash\n"
        script += "# DSSP secondary-structure script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += f"""echo "Running DSSP on {len(self.structures_stream)} structure(s)"
python "{self.helper_py}" \\
    --structures-json "{self.structures_json}" \\
    --dssp-dir "{self.stream_folder('dssp')}" \\
    --dssp-map-csv "{self.dssp_map}" \\
    --ss-dir "{self.stream_folder('ss')}" \\
    --ss-map-csv "{self.ss_map}" \\
    --secondary-csv "{self.secondary_csv}" \\
    --summary-csv "{self.summary_csv}" \\
    --container-prefix "{self.container_prefix()}" \\
    --missing-csv "{self.missing_csv}"{upstream_missing_flag}
"""
        script += self.generate_completion_check_footer()
        return script

    def get_output_files(self) -> Dict[str, Any]:
        dssp_stream = DataStream(
            name="dssp",
            ids=self.structures_stream.ids,
            files=[self.stream_path("dssp", "<id>.dssp")],
            map_table=self.dssp_map,
            format="dssp",
        )
        ss_stream = DataStream(
            name="ss",
            ids=self.structures_stream.ids,
            files=[self.stream_path("ss", "<id>.csv")],
            map_table=self.ss_map,
            format="resi-csv",
        )
        tables = {
            "secondary_structure": TableInfo(
                name="secondary_structure",
                path=self.secondary_csv,
                columns=["id", "chain", "resnum", "resi", "resname", "resn", "ss_code", "ss_simple"],
                description="DSSP 8-state secondary structure assignment per residue",
            ),
            "summary": TableInfo(
                name="summary",
                path=self.summary_csv,
                columns=["id", "n_residues", "helix_frac", "sheet_frac", "coil_frac"],
                description="DSSP per-structure SS fractions",
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "kind", "cause"],
                description="IDs removed (upstream or local failure) with removal reason",
            ),
        }
        return {
            "dssp": dssp_stream,
            "ss": ss_stream,
            "tables": tables,
            "output_folder": self.output_folder,
        }
