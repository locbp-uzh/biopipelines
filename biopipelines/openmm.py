# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""OpenMM tool: energy minimisation of protein structures.

Scope is intentionally narrow — Amber14 + implicit solvent (GBn2) energy
minimisation only. Trajectory production is GPU-bound and out of scope here;
this tool exists to clean up bad clashes / bond lengths in pipeline outputs
before downstream metric calculation.

Reference: https://github.com/openmm/openmm
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


class OpenMM(BaseConfig):
    """
    OpenMM: Amber14-SB + GBn2 implicit solvent energy minimisation.

    Scope is intentionally a minimiser, not an MD engine — no trajectory
    production. The knobs below tune the force field, solvent model, compute
    platform, and an optional positional restraint.

    Inputs:
        structures: PDB structures.
        max_iterations: minimiser step cap (default: 1000; 0 = run to tolerance).
        tolerance_kj_per_mol_nm: convergence tolerance (default: 10).
        forcefield: protein force field, one of "amber14-all" (default),
            "amber99sb", "charmm36".
        solvent: implicit-solvent model, one of "implicit-gbn2" (default),
            "implicit-gbn", "implicit-obc2", "vacuum" (no solvent).
        platform: OpenMM compute platform, one of "auto" (default, fastest
            available), "CPU", "CUDA", "OpenCL".
        restraint_selection: optional chain-aware selection string (e.g.
            "A10-50") whose heavy atoms are harmonically restrained during
            minimisation. Empty (default) = no restraint.
        restraint_k: harmonic restraint force constant in kJ/mol/nm^2
            (default 1000.0), used only when restraint_selection is set.

    Outputs:
        Streams:
            structures: one minimised <id>.pdb per input.
        Tables:
            energies: id | energy_initial_kj_mol | energy_final_kj_mol | delta_kj_mol
    """

    TOOL_NAME = "OpenMM"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        biopipelines = folders.get("biopipelines", "")
        env_check = cls._env_exists_check("openmm", env_manager)
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_check}; then
    echo "OpenMM already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block("openmm", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("openmm", env_manager, biopipelines)
        return f"""echo "=== Installing OpenMM ==="
{skip}{remove_block}
{env_block}

if {env_manager} run -n openmm python -c "import openmm; import openmm.app" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== OpenMM installation complete ==="
else
    echo "ERROR: OpenMM verification failed (cannot import openmm)"
    exit 1
fi
"""

    structures_json = Path(lambda self: self.configuration_path("structures.json"))
    structures_map = Path(lambda self: self.stream_map_path("structures"))
    energies_csv = Path(lambda self: self.table_path("energies"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_openmm.py"))

    _FORCEFIELDS = ("amber14-all", "amber99sb", "charmm36")
    _SOLVENTS = ("implicit-gbn2", "implicit-gbn", "implicit-obc2", "vacuum")
    _PLATFORMS = ("auto", "CPU", "CUDA", "OpenCL")

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 max_iterations: int = 1000,
                 tolerance_kj_per_mol_nm: float = 10.0,
                 forcefield: str = "amber14-all",
                 solvent: str = "implicit-gbn2",
                 platform: str = "auto",
                 restraint_selection: str = "",
                 restraint_k: float = 1000.0,
                 **kwargs):
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")
        self.max_iterations = int(max_iterations)
        self.tolerance_kj_per_mol_nm = float(tolerance_kj_per_mol_nm)
        self.forcefield = forcefield
        self.solvent = solvent
        self.platform = platform
        self.restraint_selection = restraint_selection
        self.restraint_k = float(restraint_k)
        super().__init__(**kwargs)

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        if self.max_iterations < 0:
            raise ValueError("max_iterations must be >= 0")
        if self.tolerance_kj_per_mol_nm <= 0:
            raise ValueError("tolerance_kj_per_mol_nm must be positive")
        if self.forcefield not in self._FORCEFIELDS:
            raise ValueError(f"forcefield must be one of {self._FORCEFIELDS}, got '{self.forcefield}'")
        if self.solvent not in self._SOLVENTS:
            raise ValueError(f"solvent must be one of {self._SOLVENTS}, got '{self.solvent}'")
        if self.platform not in self._PLATFORMS:
            raise ValueError(f"platform must be one of {self._PLATFORMS}, got '{self.platform}'")
        if self.restraint_k <= 0:
            raise ValueError("restraint_k must be positive")
        if self.restraint_selection:
            _validate_freeform_string("restraint_selection", self.restraint_selection)

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"STRUCTURES: {len(self.structures_stream)} files")
        lines.append(f"MAX ITERATIONS: {self.max_iterations or 'until tolerance'}")
        lines.append(f"TOLERANCE: {self.tolerance_kj_per_mol_nm} kJ/mol/nm")
        lines.append(f"FORCEFIELD: {self.forcefield}")
        lines.append(f"SOLVENT: {self.solvent}")
        lines.append(f"PLATFORM: {self.platform}")
        if self.restraint_selection:
            lines.append(f"RESTRAINT: {self.restraint_selection} (k={self.restraint_k} kJ/mol/nm^2)")
        return lines

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        script = "#!/bin/bash\n"
        script += "# OpenMM energy minimisation script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += f"""echo "Minimising {len(self.structures_stream)} structure(s) with OpenMM"
python "{self.helper_py}" \\
    --structures-json "{self.structures_json}" \\
    --output-dir "{self.stream_folder('structures')}" \\
    --max-iterations {self.max_iterations} \\
    --tolerance {self.tolerance_kj_per_mol_nm} \\
    --forcefield "{self.forcefield}" \\
    --solvent "{self.solvent}" \\
    --platform "{self.platform}" \\
    --restraint-selection "{self.restraint_selection}" \\
    --restraint-k {self.restraint_k} \\
    --map-csv "{self.structures_map}" \\
    --energies-csv "{self.energies_csv}"
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
            "energies": TableInfo(
                name="energies",
                path=self.energies_csv,
                columns=["id", "energy_initial_kj_mol", "energy_final_kj_mol", "delta_kj_mol"],
                description="Energies before/after OpenMM minimisation",
            ),
        }
        return {
            "structures": structures,
            "sequences": DataStream.empty("sequences", "fasta"),
            "compounds": DataStream.empty("compounds", "sdf"),
            "tables": tables,
            "output_folder": self.output_folder,
        }
