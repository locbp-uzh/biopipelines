# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""FPocket tool: binding-pocket detection.

For each input PDB, runs the fpocket binary (bioconda) and reports detected
pockets with their volume, druggability, and constituent residues.

Reference: https://github.com/Discngine/fpocket
"""

import os
from typing import Dict, List, Any, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
    from .biopipelines_io import Resolve
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream
    from biopipelines_io import Resolve


class FPocket(BaseConfig):
    """
    FPocket: detect ligand-binding pockets on protein structures.

    Inputs:
        structures: input PDB structures.
        min_alpha_spheres: minimum number of alpha-spheres per pocket
                           (fpocket `-i`, default: 35).
        min_radius: minimum alpha-sphere radius in Angstrom (fpocket `-m`, default: 3.0).
        max_radius: maximum alpha-sphere radius in Angstrom (fpocket `-M`, default: 6.0).
        clustering_distance: single-linkage clustering distance for grouping
                           alpha-spheres into pockets (fpocket `-D`, default: 2.4).

    Outputs:
        Streams:
            structures: one <id>.pdb per input — original structure annotated
                        with detected pocket alpha-spheres as HETATM STP records
                        (the `_out` suffix fpocket uses internally is dropped).
        Tables:
            pockets:  id | pocket_idx | druggability | volume | n_alpha_spheres | n_residues | residues | pocket_file
                      (pocket_file points to the per-pocket <id>_pocket<idx>_atm.pdb)
            summary:  id | n_pockets | top_druggability | top_volume | pymol_script
                      (pymol_script points to the .pml fpocket emits for visualisation)
    """

    TOOL_NAME = "FPocket"
    TOOL_VERSION = "1.1"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        biopipelines = folders.get("biopipelines", "")
        env_check = cls._env_exists_check("fpocket", env_manager)
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_check}; then
    echo "FPocket already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block("fpocket", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("fpocket", env_manager, biopipelines)
        return f"""echo "=== Installing FPocket ==="
{skip}{remove_block}
{env_block}

if {env_manager} run -n fpocket fpocket -h >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== FPocket installation complete ==="
else
    echo "ERROR: FPocket verification failed (fpocket binary not on PATH)"
    exit 1
fi
"""

    structures_json = Path(lambda self: self.configuration_path("structures.json"))
    structures_map = Path(lambda self: self.stream_map_path("structures"))
    pockets_csv = Path(lambda self: self.table_path("pockets"))
    summary_csv = Path(lambda self: self.table_path("summary"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_fpocket.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 min_alpha_spheres: int = 35,
                 min_radius: float = 3.0,
                 max_radius: float = 6.0,
                 clustering_distance: float = 2.4,
                 **kwargs):
        self.structures = structures
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")
        self.min_alpha_spheres = int(min_alpha_spheres)
        self.min_radius = float(min_radius)
        self.max_radius = float(max_radius)
        self.clustering_distance = float(clustering_distance)
        super().__init__(**kwargs)

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        if self.min_alpha_spheres < 1:
            raise ValueError("min_alpha_spheres must be >= 1")
        if self.min_radius <= 0 or self.max_radius <= 0:
            raise ValueError("min_radius and max_radius must be positive")
        if self.min_radius >= self.max_radius:
            raise ValueError("min_radius must be smaller than max_radius")
        if self.clustering_distance <= 0:
            raise ValueError("clustering_distance must be positive")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"STRUCTURES: {len(self.structures_stream)} files")
        lines.append(f"MIN ALPHA SPHERES: {self.min_alpha_spheres}")
        lines.append(f"ALPHA SPHERE RADIUS: {self.min_radius}-{self.max_radius} A")
        lines.append(f"CLUSTERING DISTANCE: {self.clustering_distance} A")
        return lines

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        upstream_missing_path = self._get_upstream_missing_table_path(self.structures)
        upstream_missing_flag = f' \\\n    --upstream-missing "{upstream_missing_path}"' if upstream_missing_path else ""
        script = "#!/bin/bash\n"
        script += "# FPocket pocket detection script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += f"""echo "Running FPocket on {len(self.structures_stream)} structure(s)"
python "{self.helper_py}" \\
    --structures-json "{self.structures_json}" \\
    --output-structures-dir "{self.stream_folder('structures')}" \\
    --structures-map-csv "{self.structures_map}" \\
    --scratch-dir "{self.extras_path()}" \\
    --min-alpha-spheres {self.min_alpha_spheres} \\
    --min-radius {self.min_radius} \\
    --max-radius {self.max_radius} \\
    --clustering-distance {self.clustering_distance} \\
    --pockets-csv "{self.pockets_csv}" \\
    --summary-csv "{self.summary_csv}" \\
    --container-prefix "{self.container_prefix()}" \\
    --missing-csv "{self.missing_csv}"{upstream_missing_flag}
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
        tables = {
            "pockets": TableInfo(
                name="pockets",
                path=self.pockets_csv,
                columns=["id", "pocket_idx", "druggability", "volume",
                         "n_alpha_spheres", "n_residues", "residues", "pocket_file"],
                description="FPocket detected pockets (one row per pocket)",
            ),
            "summary": TableInfo(
                name="summary",
                path=self.summary_csv,
                columns=["id", "n_pockets", "top_druggability", "top_volume", "pymol_script"],
                description="FPocket per-structure pocket summary",
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "kind", "cause"],
                description="IDs removed (upstream or local failure) with removal reason",
            ),
        }
        return {
            "structures": structures,
            "tables": tables,
            "output_folder": self.output_folder,
        }
