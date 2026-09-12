# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""OpenMM tool: energy minimisation of protein structures.

Scope is intentionally narrow — Amber14 + implicit solvent (GBn2) energy
minimisation only. Trajectory production is GPU-bound and out of scope here;
this tool exists to clean up bad clashes / bond lengths in pipeline outputs
before downstream metric calculation.

An optional `ligand` parameterises a bound small molecule from its SMILES via
OpenFF, and `mobile_selection` / `frozen_selection` restrict which atoms are
allowed to move, so a pocket's side chains can be relaxed against a ligand while
the rest of the structure is held rigid.

Reference: https://github.com/openmm/openmm
"""

import os
from typing import Dict, List, Any, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string, resolve_table_reference
    from .biopipelines_io import TableReference
    from .file_paths import Path
    from .datastream import DataStream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string, resolve_table_reference
    from biopipelines_io import TableReference
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
        ligand: optional compounds stream naming a small molecule bound in the
            input structures. Its HETATM block is parameterised from the
            stream's SMILES via OpenFF, so no residue template is needed.
            Without it a structure carrying an unrecognised residue fails.
        ligand_charge_method: partial-charge model for the ligand, one of
            "am1bcc" (default, semi-empirical, minutes per molecule) or "nagl"
            (graph neural net, seconds). NAGL is much faster but covers only
            organic main-group elements — it rejects Si, P and metals, and falls
            back to AM1-BCC when it does.
        ligand_forcefield: small-molecule force field, "gaff-2.11" (default) or
            any GAFF/OpenFF name openmmforcefields accepts. GAFF is the default
            because SMIRNOFF/Sage has no silicon parameters; an "openff-*" value
            covers more common organic chemistry more accurately but fails
            outright on Si, B and metals. Note that solvent="implicit-gbn2"
            cannot solvate silicon either — pair such ligands with
            "implicit-obc2".
        covalent_anchor: protein atom name holding the ligand covalently, e.g.
            "SG" or "SG62" to pin the residue number too. Empty (default) treats
            the ligand as unbonded, which lets it drift out of the site during a
            long minimisation. Errors when no such atom is within
            covalent_max_distance — a missing link is a broken structure, not
            something to silently ignore.
            The link is held by a stiff harmonic bond: it cannot break, though the
            attachment stays free to swing. A real topology bond is not offered —
            it makes the anchor residue non-standard ("CYS ... has 1 S atom too
            many") and would need a custom residue template.
        covalent_k: restraint force constant in kJ/mol/nm^2 (default 300000).
        covalent_length: restraint equilibrium length in nm (default 0.18, a
            C-S single bond).
        covalent_max_distance: how far the anchor may be from the ligand and
            still count as bonded, in Angstroms (default 3.0).
        mobile_selection: chain-aware selection whose side chains are free to
            move; every other protein atom is frozen, and the ligand stays
            mobile. Backbone atoms of the selected residues are frozen too, so
            only rotamers relax. Empty (default) = whole structure mobile.
        frozen_selection: chain-aware selection frozen outright (all atoms).
            Mutually exclusive with mobile_selection.

    Frozen atoms have their mass set to zero, which makes them immovable rather
    than merely penalised — unlike restraint_selection, which lets the whole
    structure drift under a harmonic penalty.

    Outputs:
        Streams:
            structures: one minimised <id>.pdb per input.
            compounds: the input ligand stream, passed through unchanged (the
                minimiser reuses its residue codes), so downstream tools keep the
                ligand chemistry. Absent when no ligand was given.
        Tables:
            energies: id | energy_initial_kj_mol | energy_final_kj_mol | delta_kj_mol | n_mobile_atoms | n_frozen_atoms | charge_method
    """

    TOOL_NAME = "OpenMM"
    TOOL_VERSION = "3.2"
    ENV_NAME = "openmm"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, env_name=None, **kwargs):
        name = cls._install_env(env_manager, env_name)
        biopipelines = folders.get("biopipelines", "")

        remove_block = cls._env_remove_block(name, env_manager) if force_reinstall else ""
        env_block = cls._env_install_block(name, env_manager, biopipelines)
        install_block = f"{remove_block}\n{env_block}"
        if not force_reinstall:
            # An existing env still gets verified: the checks below are what a
            # half-built env fails, and skipping straight to the sentinel would
            # report that as a success.
            install_block = f"""if {cls._env_exists_check(name, env_manager)}; then
    echo "OpenMM env already present, verifying."
else
{install_block}
fi
"""
        verify = ("import openmm; import openmm.app; import openmmforcefields.generators; "
                  "import openff.toolkit")
        try:
            from .config_manager import ConfigManager as _CM
        except ImportError:
            from config_manager import ConfigManager as _CM
        # get_conda_env_root() is only defined for the venv manager; elsewhere the env run
        # prefix already puts the right bin on PATH.
        if _CM().get_env_manager() == "venv":
            conda_bin = f"{_CM().get_conda_env_root()}/{name}/bin"
            # Run antechamber rather than testing that the file exists: it is a wrapper
            # script sourcing ../amber.sh, so it can be present on PATH and still fail
            # on every invocation.
            amber_check = f"""
if ! PATH="{conda_bin}:$PATH" antechamber -h >/dev/null 2>&1; then
    echo "ERROR: antechamber is present but does not run — AM1-BCC charges would fail"
    PATH="{conda_bin}:$PATH" antechamber -h 2>&1 | tail -3
    exit 1
fi
"""
        else:
            amber_check = f"""
if ! {cls._env_run(name, env_manager)}antechamber -h >/dev/null 2>&1; then
    echo "ERROR: antechamber is present but does not run — AM1-BCC charges would fail"
    {cls._env_run(name, env_manager)}antechamber -h 2>&1 | tail -3
    exit 1
fi
"""
        return f"""echo "=== Installing OpenMM ({name}) ==="
{install_block}
{amber_check}
if {cls._env_run(name, env_manager)}python -c "{verify}" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== OpenMM installation complete ==="
else
    echo "ERROR: OpenMM verification failed"
    {cls._env_run(name, env_manager)}python -c "{verify}" 2>&1 | tail -5
    exit 1
fi
"""

    structures_json = Path(lambda self: self.configuration_path("structures.json"))
    ligand_json = Path(lambda self: self.configuration_path("ligand.json"))
    structures_map = Path(lambda self: self.stream_map_path("structures"))
    energies_csv = Path(lambda self: self.table_path("energies"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_openmm.py"))

    _FORCEFIELDS = ("amber14-all", "amber99sb", "charmm36")
    _SOLVENTS = ("implicit-gbn2", "implicit-gbn", "implicit-obc2", "vacuum")
    _PLATFORMS = ("auto", "CPU", "CUDA", "OpenCL")
    _CHARGE_METHODS = ("am1bcc", "nagl")

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 max_iterations: int = 1000,
                 tolerance_kj_per_mol_nm: float = 10.0,
                 forcefield: str = "amber14-all",
                 solvent: str = "implicit-gbn2",
                 platform: str = "auto",
                 restraint_selection: str = "",
                 restraint_k: float = 1000.0,
                 ligand: Union[DataStream, StandardizedOutput, None] = None,
                 ligand_charge_method: str = "am1bcc",
                 ligand_forcefield: str = "gaff-2.11",
                 covalent_anchor: str = "",
                 covalent_k: float = 300000.0,
                 covalent_length: float = 0.18,
                 covalent_max_distance: float = 3.0,
                 mobile_selection: str = "",
                 frozen_selection: str = "",
                 **kwargs):
        # Keep the original input for upstream missing-table detection.
        self.structures_input = structures
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
        self.ligand_input = ligand
        if ligand is None:
            self.ligand_stream = None
        elif isinstance(ligand, StandardizedOutput):
            self.ligand_stream = ligand.streams.compounds
        elif isinstance(ligand, DataStream):
            self.ligand_stream = ligand
        else:
            raise ValueError(f"ligand must be DataStream or StandardizedOutput, got {type(ligand)}")
        self.ligand_charge_method = ligand_charge_method
        self.ligand_forcefield = ligand_forcefield
        self.covalent_anchor = covalent_anchor
        self.covalent_k = float(covalent_k)
        self.covalent_length = float(covalent_length)
        self.covalent_max_distance = float(covalent_max_distance)
        self.mobile_selection = resolve_table_reference(mobile_selection, "mobile_selection")
        self.frozen_selection = resolve_table_reference(frozen_selection, "frozen_selection")
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
        if self.mobile_selection and self.frozen_selection:
            raise ValueError(
                "mobile_selection and frozen_selection are mutually exclusive: "
                "mobile_selection already freezes everything it does not name")
        # A TableReference resolves per-structure at runtime, so only a literal
        # string reaches the generated script and can be checked here.
        if self.mobile_selection and not isinstance(self.mobile_selection, TableReference):
            _validate_freeform_string("mobile_selection", self.mobile_selection)
        if self.frozen_selection and not isinstance(self.frozen_selection, TableReference):
            _validate_freeform_string("frozen_selection", self.frozen_selection)
        if self.ligand_charge_method not in self._CHARGE_METHODS:
            raise ValueError(
                f"ligand_charge_method must be one of {self._CHARGE_METHODS}, "
                f"got '{self.ligand_charge_method}'")
        if self.ligand_stream is not None and len(self.ligand_stream) == 0:
            raise ValueError("ligand stream is empty")
        if self.ligand_stream is not None and not self.ligand_forcefield:
            raise ValueError("ligand_forcefield must not be empty when ligand is set")
        if self.covalent_anchor:
            if self.ligand_stream is None:
                raise ValueError("covalent_anchor requires ligand to be set")
            _validate_freeform_string("covalent_anchor", self.covalent_anchor)
            if self.covalent_k <= 0:
                raise ValueError("covalent_k must be positive")
            if self.covalent_length <= 0:
                raise ValueError("covalent_length must be positive (nanometres)")
            if self.covalent_max_distance <= 0:
                raise ValueError("covalent_max_distance must be positive")

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
        if self.ligand_stream is not None:
            lines.append(f"LIGAND: {len(self.ligand_stream)} compound(s), charges={self.ligand_charge_method}")
        if self.covalent_anchor:
            lines.append(f"COVALENT: {self.covalent_anchor}-ligand restraint")
        if self.mobile_selection:
            lines.append(f"MOBILE: {self.mobile_selection} side chains (+ ligand); rest frozen")
        if self.frozen_selection:
            lines.append(f"FROZEN: {self.frozen_selection}")
        return lines

    def _ambertools_path_block(self) -> str:
        """Put the conda env's bin on PATH so AM1-BCC can reach AmberTools.

        antechamber is a wrapper script that sources ``../amber.sh`` relative to its
        own location, so a symlink into a venv shim breaks it — the binary must be
        called from the directory it ships in. No-op when the ligand path is unused
        or the env manager already puts it on PATH.
        """
        if self.ligand_stream is None:
            return ""
        try:
            from .config_manager import ConfigManager
        except ImportError:
            from config_manager import ConfigManager
        cm = ConfigManager()
        if cm.get_env_manager() != "venv":
            return ""
        env = cm.get_environment(self.TOOL_NAME) or "openmm"
        conda_bin = f"{cm.get_conda_env_root()}/{env}/bin"
        return (f'if [ -d "{conda_bin}" ]; then\n'
                f'    export PATH="{conda_bin}:$PATH"\n'
                f'    export AMBERHOME="{cm.get_conda_env_root()}/{env}"\n'
                f'fi\n')

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        if self.ligand_stream is not None:
            self.ligand_stream.save_json(self.ligand_json)
        script = "#!/bin/bash\n"
        script += "# OpenMM energy minimisation script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += self.warn_container_unsupported()
        script += self._ambertools_path_block()
        ligand_arg = f' \\\n    --ligand-json "{self.ligand_json}"' if self.ligand_stream is not None else ""
        mobile = str(self.mobile_selection) if self.mobile_selection else ""
        frozen = str(self.frozen_selection) if self.frozen_selection else ""
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
    --mobile-selection "{mobile}" \\
    --frozen-selection "{frozen}" \\
    --charge-method "{self.ligand_charge_method}" \\
    --ligand-forcefield "{self.ligand_forcefield}" \\
    --covalent-anchor "{self.covalent_anchor}" \\
    --covalent-k {self.covalent_k} \\
    --covalent-length {self.covalent_length} \\
    --covalent-max-distance {self.covalent_max_distance}{ligand_arg} \\
    --map-csv "{self.structures_map}" \\
    --energies-csv "{self.energies_csv}"
"""
        # Excuse ids an upstream filter removed, so they are not reported as our
        # failures. Panda over-declares in pool mode: it cannot know at config time
        # which rows survive, and writes the surviving set to its missing manifest.
        script += self.generate_missing_propagation(
            self.structures_input, self.ligand_input, missing_csv=self.missing_csv
        )
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
                columns=["id", "energy_initial_kj_mol", "energy_final_kj_mol", "delta_kj_mol",
                         "n_mobile_atoms", "n_frozen_atoms", "charge_method"],
                description="Energies before/after OpenMM minimisation",
            ),
        }
        if self._collect_upstream_missing_paths(self.structures_input, self.ligand_input):
            tables["missing"] = self.missing_table_info(self.missing_csv)
        return {
            "structures": structures,
            "sequences": DataStream.empty("sequences", "fasta"),
            # Pass the ligand chemistry through: the runtime reuses the input's
            # residue codes rather than assigning its own, so the stream is
            # unchanged and downstream tools would otherwise lose it.
            "compounds": (self.ligand_stream if self.ligand_stream is not None
                          else DataStream.empty("compounds", "csv")),
            "tables": tables,
            "output_folder": self.output_folder,
        }
