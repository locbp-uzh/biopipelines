# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Reduce tool: add explicit hydrogens to a structures stream via the
Richardson lab's `reduce` binary. Preserves protein residue topology and
ligand HETATM records, so the output is suitable as input to
interaction-fingerprint or MD-prep tools.

Reference: https://github.com/rlabduke/reduce
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


class Reduce(BaseConfig):
    """
    Reduce: add explicit hydrogens to a structures stream.

    Inputs:
        structures: PDB structures (proteins, complexes, ligand-bearing).

    Outputs:
        Streams:
            structures: one protonated <id>.pdb per input.
    """

    TOOL_NAME = "Reduce"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        biopipelines = folders.get("biopipelines", "")
        env_check = cls._env_exists_check("reduce", env_manager)
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_check}; then
    echo "Reduce already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block("reduce", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("reduce", env_manager, biopipelines)
        return f"""echo "=== Installing Reduce ==="
{skip}{remove_block}
{env_block}

if {env_manager} run -n reduce reduce -version 2>&1 | grep -q '^reduce'; then
    touch "$INSTALL_SUCCESS"
    echo "=== Reduce installation complete ==="
else
    echo "ERROR: Reduce verification failed (reduce binary not on PATH)"
    exit 1
fi
"""

    structures_json = Path(lambda self: self.configuration_path("structures.json"))
    structures_map = Path(lambda self: self.stream_map_path("structures"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_reduce.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 **kwargs):
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
        script = "#!/bin/bash\n"
        script += "# Reduce hydrogen-addition script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += f"""echo "Adding hydrogens to {len(self.structures_stream)} structure(s)"
python "{self.helper_py}" \\
    --structures-json "{self.structures_json}" \\
    --output-dir "{self.stream_folder('structures')}" \\
    --map-csv "{self.structures_map}"
"""
        script += self.generate_completion_check_footer()
        return script

    def get_output_files(self) -> Dict[str, Any]:
        structures = DataStream(
            name="structures",
            ids=self.structures_stream.ids,
            files=[self.stream_path("structures", "<id>.pdb")],
            map_table=self.structures_map,
            format="pdb",
        )
        return {
            "structures": structures,
            "sequences": DataStream.empty("sequences", "fasta"),
            "compounds": DataStream.empty("compounds", "sdf"),
            "tables": {},
            "output_folder": self.output_folder,
        }
