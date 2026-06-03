# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
PoseBusters validation for computationally generated molecule poses.

Validates docking/prediction outputs by checking bond lengths, angles,
steric clashes, protein interactions, stereochemistry, and more.
Supports ``dock`` mode (ligand + protein) and ``redock`` mode
(+ reference ligand for RMSD comparison).

Reference:
    Buttenschoen et al. (2024) PoseBusters: AI-based docking is not solved.
    Chem Sci 15, 3130. https://github.com/maabuu/posebusters
"""

import os
import json
from typing import Dict, List, Any, Optional, Union


# All boolean check columns PoseBusters emits, keyed by category.
# Names match what posebusters writes after rename_outputs (lowercased, spaces → "_").
_CHECK_GROUPS: Dict[str, List[str]] = {
    "loading": ["mol_pred_loaded", "mol_cond_loaded", "mol_true_loaded"],
    "sanity": ["sanitization", "inchi_convertible"],
    "radicals": ["no_radicals"],
    "connectivity": ["all_atoms_connected"],
    "geometry": [
        "bond_lengths", "bond_angles", "internal_steric_clash",
        "aromatic_ring_flatness", "non-aromatic_ring_non-flatness",
        "double_bond_flatness", "internal_energy",
    ],
    "intermolecular": [
        "protein-ligand_maximum_distance",
        "minimum_distance_to_protein",
        "minimum_distance_to_organic_cofactors",
        "minimum_distance_to_inorganic_cofactors",
        "minimum_distance_to_waters",
        "volume_overlap_with_protein",
        "volume_overlap_with_organic_cofactors",
        "volume_overlap_with_inorganic_cofactors",
        "volume_overlap_with_waters",
    ],
    "identity": [
        "molecular_formula", "molecular_bonds",
        "double_bond_stereochemistry", "tetrahedral_chirality",
        "stereochemistry_preserved",
    ],
    "rmsd": ["rmsd_≤_2å"],
}

# Presets selecting which check columns contribute to all_pass.
# "pose" (default) drops only no_radicals — predicted poses lack proper bond
# perception, so RDKit's check_radicals flags any unsaturated atom as a radical
# even when the pose itself is fine. Sanitization and inchi_convertible are
# kept since they catch genuine RDKit-parsing problems.
_CHECK_PRESETS: Dict[str, List[str]] = {
    "all": (
        _CHECK_GROUPS["loading"] + _CHECK_GROUPS["sanity"]
        + _CHECK_GROUPS["radicals"] + _CHECK_GROUPS["connectivity"]
        + _CHECK_GROUPS["geometry"] + _CHECK_GROUPS["intermolecular"]
        + _CHECK_GROUPS["identity"] + _CHECK_GROUPS["rmsd"]
    ),
    "pose": (
        _CHECK_GROUPS["loading"] + _CHECK_GROUPS["sanity"]
        + _CHECK_GROUPS["connectivity"] + _CHECK_GROUPS["geometry"]
        + _CHECK_GROUPS["intermolecular"]
    ),
    "loading": list(_CHECK_GROUPS["loading"]),
    "geometry": _CHECK_GROUPS["loading"] + _CHECK_GROUPS["geometry"],
    "chemistry": (
        _CHECK_GROUPS["loading"] + _CHECK_GROUPS["sanity"]
        + _CHECK_GROUPS["radicals"] + _CHECK_GROUPS["connectivity"]
    ),
}

_ALL_KNOWN_CHECKS = set(_CHECK_PRESETS["all"])

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .input_standardization import resolve_basic_input
    from .ligand import Ligand
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from input_standardization import resolve_basic_input
    from ligand import Ligand


class PoseBusters(BaseConfig):
    """
    Pipeline tool for validating computationally generated molecule poses.

    Takes protein-ligand complex structures and runs PoseBusters quality checks
    including bond lengths, bond angles, internal steric clashes, volume overlap
    with protein, and (in redock mode) RMSD to a reference ligand pose.

    Generates CSV with boolean pass/fail columns for each check and an
    ``all_pass`` summary column.

    Commonly used for:
    - Validating docking outputs from Gnina
    - Validating structure predictions from Boltz2
    - Quality-gating predicted poses before downstream analysis
    """

    TOOL_NAME = "PoseBusters"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        """Install PoseBusters in a dedicated conda environment.

        Creates a ``posebusters`` conda environment with Python 3.10 and
        installs the ``posebusters`` package via pip.
        """
        biopipelines = folders.get("biopipelines", "")
        env_check = cls._env_exists_check("posebusters", env_manager)
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_check}; then
    echo "PoseBusters already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block("posebusters", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("posebusters", env_manager, biopipelines)
        return f"""echo "=== Installing PoseBusters ==="
{skip}{remove_block}
{env_block}

# Verify installation. Check pandas too (not just posebusters): pipe_posebusters.py
# and biopipelines/base_config both import pandas at module load, so a pandas that
# can't import (e.g. missing python-dateutil) silently breaks every run while
# `import posebusters` still succeeds.
if {env_manager} run -n posebusters python -c "import posebusters, pandas, numpy, gemmi" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== PoseBusters installation complete ==="
else
    echo "ERROR: PoseBusters verification failed (import of posebusters/pandas/numpy/gemmi)"
    {env_manager} run -n posebusters python -c "import posebusters, pandas, numpy, gemmi" 2>&1 | tail -5
    exit 1
fi
"""

    # Lazy path descriptors
    analysis_csv = Path(lambda self: self.table_path("posebusters"))
    config_file = Path(lambda self: self.configuration_path("posebusters_config.json"))
    structures_ds_json = Path(lambda self: self.configuration_path("structures.json"))
    ligand_json = Path(lambda self: self.configuration_path("input_ligand.json"))
    posebusters_py = Path(lambda self: self.pipe_script_path("pipe_posebusters.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 ligand: Union[str, DataStream, StandardizedOutput],
                 reference_ligand: Union[DataStream, StandardizedOutput, None] = None,
                 mode: str = "dock",
                 check: Union[str, List[str]] = "pose",
                 exclude: Union[str, List[str], None] = None,
                 **kwargs):
        """
        Initialize PoseBusters validation tool.

        Args:
            structures: Input protein-ligand complexes as DataStream or StandardizedOutput
            ligand: Compounds stream (Ligand(code="LIG") or any compounds-producing
                    tool) naming the ligand. The residue `code` is read from the
                    stream's `code` column at runtime, and is also used as the
                    reference-structure residue code in redock mode.

                    IMPORTANT — provide a SMILES for non-trivial ligands. The ligand
                    is extracted from each complex and rebuilt as an SDF for validation.
                    When the stream carries a SMILES (e.g. ``Ligand(smiles=..., codes="LIG")``),
                    bond orders are templated from it; otherwise they are perceived from
                    coordinates alone (RDKit ``rdDetermineBonds``). Coordinate-only
                    perception FAILS on charged/conjugated ligands (dyes, metallo-organics,
                    zwitterions): the molecule won't build, so `mol_pred_loaded` fails and
                    EVERY check — including the coordinate-based distance/overlap checks —
                    reports False (the whole step yields all_pass=False). For such ligands
                    you MUST pass a SMILES-bearing Ligand whose ``codes`` matches the
                    residue code in the structures. Bare ``Ligand(code="LIG")`` (no SMILES)
                    is fine only for simple ligands RDKit can perceive from coordinates.
                    A bare code with no SMILES also triggers a (failing) RCSB SMILES
                    lookup unless the code is a real CCD code.
            reference_ligand: Reference ligand structure for redock mode (SDF or PDB)
            mode: 'dock' or 'redock' (auto-set to 'redock' if reference_ligand provided)
            check: Which checks contribute to ``all_pass``. Either a preset string or a
                list of column names.

                Presets:
                  - ``"pose"`` (default): everything except ``no_radicals``. That
                    one check fails spuriously on predicted poses because RDKit's
                    bond perception from coordinates flags unsaturated atoms as
                    radicals.
                  - ``"all"``: every available check (PoseBusters default).
                  - ``"loading"``: only that the molecule loaded.
                  - ``"geometry"``: loading + bond lengths/angles/clashes/flatness/energy.
                  - ``"chemistry"``: loading + sanity + radicals + connectivity.

                List form: ``check=["all_atoms_connected", "bond_lengths"]`` — only these
                columns determine ``all_pass``. All check columns are still written to
                the CSV; excluded columns are placed to the right of ``all_pass``.

            exclude: Check column(s) to drop from whatever ``check`` selected — a name
                or list of names. Easier than re-listing every check to keep when you
                only need to omit one or two. E.g. ``exclude="minimum_distance_to_protein"``
                keeps the ``pose`` preset but ignores the protein-distance check (which
                fails by construction on covalent-derived poses, where the attachment
                atom sits at bonding distance).

            **kwargs: Additional parameters

        Output:
            Streams: (none)
            Tables:
                posebusters: id | source_structure | <included checks...> | all_pass | <excluded checks...>
                In redock mode the column set additionally includes mol_true_loaded,
                identity columns, and rmsd_≤_2å.
        """
        # Resolve structures input to DataStream
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        # Ligand naming — a compounds stream; the residue code is resolved from
        # its `code` column at runtime. A bare string is shorthand for an
        # internal Ligand(code=...).
        self.ligand_stream: DataStream = resolve_basic_input(
            ligand, Ligand, "compounds", "code", allow_none=False)

        # Resolve reference_ligand input
        self.reference_ligand_stream: Optional[DataStream] = None
        if reference_ligand is not None:
            if isinstance(reference_ligand, StandardizedOutput):
                self.reference_ligand_stream = reference_ligand.streams.structures
            elif isinstance(reference_ligand, DataStream):
                self.reference_ligand_stream = reference_ligand
            else:
                raise ValueError(f"reference_ligand must be DataStream, StandardizedOutput, or None, got {type(reference_ligand)}")

        # Auto-set mode to redock if reference_ligand is provided
        if self.reference_ligand_stream is not None and mode == "dock":
            self.mode = "redock"
        else:
            self.mode = mode

        self.check = check
        self.exclude = exclude
        self._check_columns = self._resolve_check(check, exclude)

        super().__init__(**kwargs)

    @staticmethod
    def _resolve_check(check: Union[str, List[str]],
                       exclude: Union[str, List[str], None] = None) -> List[str]:
        """Resolve a check spec (preset name or list) to a list of column names,
        then drop any names in ``exclude``."""
        if isinstance(check, str):
            if check not in _CHECK_PRESETS:
                raise ValueError(
                    f"check preset '{check}' is unknown. "
                    f"Valid presets: {sorted(_CHECK_PRESETS)}, or pass a list of column names."
                )
            selected = list(_CHECK_PRESETS[check])
        elif isinstance(check, (list, tuple)):
            selected = []
            unknown = []
            for c in check:
                if not isinstance(c, str):
                    raise ValueError(f"check list entries must be strings, got {type(c)}")
                if c not in _ALL_KNOWN_CHECKS:
                    unknown.append(c)
                else:
                    selected.append(c)
            if unknown:
                raise ValueError(
                    f"check contains unknown column(s) {unknown}. "
                    f"Known: {sorted(_ALL_KNOWN_CHECKS)}"
                )
            if not selected:
                raise ValueError("check list is empty")
        else:
            raise ValueError(f"check must be a preset name or list of column names, got {type(check)}")

        # Drop any excluded checks from the selection.
        if exclude is not None:
            excl = [exclude] if isinstance(exclude, str) else list(exclude)
            unknown_excl = [c for c in excl if c not in _ALL_KNOWN_CHECKS]
            if unknown_excl:
                raise ValueError(
                    f"exclude contains unknown column(s) {unknown_excl}. "
                    f"Known: {sorted(_ALL_KNOWN_CHECKS)}"
                )
            excl_set = set(excl)
            selected = [c for c in selected if c not in excl_set]
            if not selected:
                raise ValueError("all checks were excluded; nothing left to determine all_pass")
        return selected

    def validate_params(self):
        """Validate PoseBusters parameters."""
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures cannot be empty")

        if not self.ligand_stream or len(self.ligand_stream) == 0:
            raise ValueError("ligand (a compounds stream, e.g. Ligand(code=...)) is required and must not be empty")

        if self.mode not in ("dock", "redock"):
            raise ValueError(f"mode must be 'dock' or 'redock', got '{self.mode}'")

        if self.mode == "redock" and self.reference_ligand_stream is None:
            raise ValueError("redock mode requires reference_ligand to be provided")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input structures."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        check_display = self.check if isinstance(self.check, str) else f"[{len(self._check_columns)} columns]"
        config_lines.extend([
            f"STRUCTURES: {len(self.structures_stream)} files",
            "LIGAND: (code resolved from compounds stream at runtime)",
            f"MODE: {self.mode}",
            f"CHECK: {check_display}",
        ])

        if self.reference_ligand_stream is not None:
            config_lines.append(f"REFERENCE LIGAND: {len(self.reference_ligand_stream)} files")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate PoseBusters execution script."""
        script_content = "#!/bin/bash\n"
        script_content += "# PoseBusters execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self.generate_script_run_posebusters()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_posebusters(self) -> str:
        """Generate the PoseBusters validation part of the script."""
        # Serialize structures + ligand DataStreams to JSON. The ligand residue
        # code is resolved from the compounds stream at runtime by the pipe
        # script (also used as the reference-structure code in redock mode).
        self.structures_stream.save_json(self.structures_ds_json)
        self.ligand_stream.save_json(self.ligand_json)

        # Build reference_pdb path for redock mode
        reference_pdb = None
        if self.reference_ligand_stream is not None:
            ref_ds_json = self.configuration_path("reference_ligand.json")
            self.reference_ligand_stream.save_json(ref_ds_json)
            reference_pdb = ref_ds_json

        config_data = {
            "structures_json": self.structures_ds_json,
            "ligand_json": self.ligand_json,
            "mode": self.mode,
            "reference_pdb": reference_pdb,
            "output_csv": self.analysis_csv,
            "execution_dir": self.execution_folder,
            "check_columns": self._check_columns,
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        return f"""echo "Running PoseBusters validation"
echo "Structures: {len(self.structures_stream)}"
echo "Mode: {self.mode}"
echo "Output: {self.analysis_csv}"

python "{self.posebusters_py}" --config "{self.config_file}"

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after PoseBusters validation."""
        # Full set of check columns PoseBusters can emit for this mode
        if self.mode == "redock":
            mode_checks = _CHECK_PRESETS["all"]
        else:
            mode_checks = [c for c in _CHECK_PRESETS["all"]
                           if c not in _CHECK_GROUPS["identity"]
                           and c not in _CHECK_GROUPS["rmsd"]
                           and c != "mol_true_loaded"]

        included = [c for c in self._check_columns if c in mode_checks]
        excluded = [c for c in mode_checks if c not in included]

        columns = ["id", "source_structure"] + included + ["all_pass"] + excluded

        tables = {
            "posebusters": TableInfo(
                name="posebusters",
                path=self.analysis_csv,
                columns=columns,
                description=f"PoseBusters validation ({self.mode} mode)"
            )
        }

        return {
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                "ligand_ids": list(self.ligand_stream.ids),
                "mode": self.mode,
                "check": self.check,
                "check_columns": self._check_columns,
            }
        })
        return base_dict
