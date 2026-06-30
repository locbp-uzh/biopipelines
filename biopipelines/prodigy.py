# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Prodigy tool: binding-affinity prediction for protein-protein complexes.

Wraps `prodigy-prot` (https://github.com/haddocking/prodigy). For each complex
PDB, reports predicted dissociation constant Kd and binding free energy dG
between the specified chain groups.
"""

import os
from typing import Dict, List, Any, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream


class Prodigy(BaseConfig):
    """
    Prodigy: predict binding affinity for protein-protein complexes.

    Inputs:
        structures: PDB complexes.
        interface: space-separated chain groups, each group a comma-joined
                   list of chains treated as one molecule. Examples:
                   "A B"   -> chain A vs chain B
                   "A,B C" -> chains A+B together vs chain C
                   "A B C" -> all pairwise interfaces among A, B, C.
        temperature: temperature in Celsius (default: 25.0).

    Outputs:
        Streams: (none)
        Tables:
            affinity: id | interface | delta_g_kcal_mol | kd_M | n_intermol_contacts | percent_charged_nis | percent_apolar_nis
    """

    TOOL_NAME = "Prodigy"
    TOOL_VERSION = "1.1"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        biopipelines = folders.get("biopipelines", "")
        env_check = cls._env_exists_check("prodigy", env_manager)
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_check}; then
    echo "Prodigy already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block("prodigy", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("prodigy", env_manager, biopipelines)
        return f"""echo "=== Installing Prodigy ==="
{skip}{remove_block}
{env_block}

if {env_manager} run -n prodigy python -c "import prodigy_prot" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== Prodigy installation complete ==="
else
    echo "ERROR: Prodigy verification failed (cannot import prodigy_prot)"
    exit 1
fi
"""

    structures_json = Path(lambda self: self.configuration_path("structures.json"))
    affinity_csv = Path(lambda self: self.table_path("affinity"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_prodigy.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 interface: str = "A B",
                 temperature: float = 25.0,
                 **kwargs):
        self.structures = structures
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")
        self.interface = interface
        self.temperature = float(temperature)
        super().__init__(**kwargs)

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        if not self.interface:
            raise ValueError("interface is required (e.g. 'A B')")
        _validate_freeform_string("interface", self.interface)

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"STRUCTURES: {len(self.structures_stream)} complexes")
        lines.append(f"INTERFACE: {self.interface}")
        lines.append(f"TEMPERATURE: {self.temperature} C")
        return lines

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        script = "#!/bin/bash\n"
        script += "# Prodigy affinity script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += self.warn_container_unsupported()
        upstream_missing_flag = self.upstream_missing_flag(self.structures)
        script += f"""echo "Running Prodigy on {len(self.structures_stream)} complex(es)"
python "{self.helper_py}" \\
    --structures-json "{self.structures_json}" \\
    --interface "{self.interface}" \\
    --temperature {self.temperature} \\
    --affinity-csv "{self.affinity_csv}" \\
    --missing-csv "{self.missing_csv}"{upstream_missing_flag}
"""
        script += self.generate_completion_check_footer()
        return script

    def get_output_files(self) -> Dict[str, Any]:
        tables = {
            "affinity": TableInfo(
                name="affinity",
                path=self.affinity_csv,
                columns=["id", "interface", "delta_g_kcal_mol", "kd_M",
                         "n_intermol_contacts", "percent_charged_nis", "percent_apolar_nis"],
                description="Prodigy predicted binding affinity per complex",
            ),
            "missing": self.missing_table_info(self.missing_csv),
        }
        return {"tables": tables, "output_folder": self.output_folder}
