# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""APBS tool: electrostatic surface potential.

For each input PDB:
  1. Runs `pdb2pqr` (with PROPKA at physiological pH 7.4 by default) to assign protonation states.
  2. Runs `apbs` against the generated input file to compute the
     Poisson-Boltzmann electrostatic potential on a grid.
  3. Exposes the protonated/charged structure as a `structures` stream (PQR)
     and the electrostatic potential grid as a `grids` stream (DX).
  4. Reports per-structure metrics (net charge, n_basic, n_acidic, isoelectric_point)
     parsed from the pdb2pqr summary in `tables/electrostatics.csv`.
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


class APBS(BaseConfig):
    """
    APBS: Poisson-Boltzmann electrostatic surface potential.

    Inputs:
        structures: input PDB structures.
        ph: pH for pdb2pqr protonation assignment (default: 7.4, physiological).
        forcefield: pdb2pqr force field for charges/radii (`--ff`), one of
            "AMBER", "CHARMM", "PARSE", "TYL06", "PEOEPB", "SWANSON"
            (default "AMBER").
        ion_concentration: mobile-ion concentration in mol/L for the +1/-1
            species in the Poisson-Boltzmann calculation (default 0.150).
        grid_dim: cubic grid dimension (APBS `dime`, points per axis;
            default 65). Larger = finer grid, more memory.
        pdie: solute (protein) dielectric constant (default 2.0).
        sdie: solvent dielectric constant (default 78.5, water at 298 K).
        solver: Poisson-Boltzmann solver, "lpbe" (linearized, default) or
            "npbe" (nonlinear).

    Outputs:
        Streams:
            structures: one <id>.pqr per input (pdb2pqr protonated structure
                        with PROPKA-assigned charges and Amber radii).
            grids:      one <id>.dx per input (APBS Poisson-Boltzmann
                        electrostatic potential grid).
        Tables:
            electrostatics: id | net_charge | n_basic | n_acidic | pI | mean_potential
    """

    TOOL_NAME = "APBS"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        biopipelines = folders.get("biopipelines", "")
        env_check = cls._env_exists_check("apbs", env_manager)
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_check}; then
    echo "APBS already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block("apbs", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("apbs", env_manager, biopipelines)
        return f"""echo "=== Installing APBS ==="
{skip}{remove_block}
{env_block}

if {env_manager} run -n apbs apbs --version 2>&1 | grep -q "APBS" && {env_manager} run -n apbs pdb2pqr --help >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== APBS installation complete ==="
else
    echo "ERROR: APBS verification failed (apbs or pdb2pqr binary missing)"
    exit 1
fi
"""

    structures_json = Path(lambda self: self.configuration_path("structures.json"))
    pqr_map = Path(lambda self: self.stream_map_path("structures"))
    grid_map = Path(lambda self: self.stream_map_path("grids"))
    electrostatics_csv = Path(lambda self: self.table_path("electrostatics"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_apbs.py"))

    _FORCEFIELDS = ("AMBER", "CHARMM", "PARSE", "TYL06", "PEOEPB", "SWANSON")
    _SOLVERS = ("lpbe", "npbe")

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 ph: float = 7.4,
                 forcefield: str = "AMBER",
                 ion_concentration: float = 0.150,
                 grid_dim: int = 65,
                 pdie: float = 2.0,
                 sdie: float = 78.5,
                 solver: str = "lpbe",
                 **kwargs):
        self.structures = structures
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")
        self.ph = float(ph)
        self.forcefield = forcefield
        self.ion_concentration = float(ion_concentration)
        self.grid_dim = int(grid_dim)
        self.pdie = float(pdie)
        self.sdie = float(sdie)
        self.solver = solver
        super().__init__(**kwargs)

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        if not (0 <= self.ph <= 14):
            raise ValueError(f"ph must be in [0, 14], got {self.ph}")
        if self.forcefield not in self._FORCEFIELDS:
            raise ValueError(f"forcefield must be one of {self._FORCEFIELDS}, got '{self.forcefield}'")
        if self.solver not in self._SOLVERS:
            raise ValueError(f"solver must be one of {self._SOLVERS}, got '{self.solver}'")
        if self.ion_concentration < 0:
            raise ValueError("ion_concentration must be non-negative")
        if self.grid_dim < 1:
            raise ValueError("grid_dim must be a positive integer")
        if self.pdie <= 0 or self.sdie <= 0:
            raise ValueError("pdie and sdie must be positive")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"STRUCTURES: {len(self.structures_stream)} files")
        lines.append(f"pH: {self.ph}")
        lines.append(f"FORCEFIELD: {self.forcefield}")
        lines.append(f"ION CONCENTRATION: {self.ion_concentration} mol/L")
        lines.append(f"GRID DIM: {self.grid_dim}")
        lines.append(f"DIELECTRICS: pdie={self.pdie} sdie={self.sdie}")
        lines.append(f"SOLVER: {self.solver}")
        return lines

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        upstream_missing_path = self._get_upstream_missing_table_path(self.structures)
        upstream_missing_flag = f' \\\n    --upstream-missing "{upstream_missing_path}"' if upstream_missing_path else ""
        script = "#!/bin/bash\n"
        script += "# APBS electrostatic potential script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += f"""echo "Running APBS on {len(self.structures_stream)} structure(s)"
python "{self.helper_py}" \\
    --structures-json "{self.structures_json}" \\
    --pqr-dir "{self.stream_folder('structures')}" \\
    --pqr-map-csv "{self.pqr_map}" \\
    --grid-dir "{self.stream_folder('grids')}" \\
    --grid-map-csv "{self.grid_map}" \\
    --scratch-dir "{self.extras_path()}" \\
    --ph {self.ph} \\
    --forcefield "{self.forcefield}" \\
    --ion-concentration {self.ion_concentration} \\
    --grid-dim {self.grid_dim} \\
    --pdie {self.pdie} \\
    --sdie {self.sdie} \\
    --solver "{self.solver}" \\
    --electrostatics-csv "{self.electrostatics_csv}" \\
    --missing-csv "{self.missing_csv}"{upstream_missing_flag}
"""
        script += self.generate_completion_check_footer()
        return script

    def get_output_files(self) -> Dict[str, Any]:
        ids = self.structures_stream.ids
        structures = DataStream(
            name="structures",
            ids=ids,
            files=[self.stream_path("structures", "<id>.pqr")],
            map_table=self.pqr_map,
            format="pqr",
        )
        grids = DataStream(
            name="grids",
            ids=ids,
            files=[self.stream_path("grids", "<id>.dx")],
            map_table=self.grid_map,
            format="dx",
        )
        tables = {
            "electrostatics": TableInfo(
                name="electrostatics",
                path=self.electrostatics_csv,
                columns=["id", "net_charge", "n_basic", "n_acidic",
                         "pI", "mean_potential"],
                description="APBS electrostatic metrics per structure",
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
            "grids": grids,
            "tables": tables,
            "output_folder": self.output_folder,
        }
