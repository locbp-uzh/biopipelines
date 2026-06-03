# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""XTB tool: semi-empirical GFN2-xTB interaction-energy scoring.

For each complex PDB, splits the structure into a ligand fragment (selected
by 3-letter residue code) and the surrounding protein, runs single-point
GFN2-xTB calculations on protein, ligand, and complex, and reports the
interaction energy E_complex - E_protein - E_ligand. Useful as a
physics-grounded ranking signal complementary to ML scoring functions.

NOTE: a full GFN2-xTB single-point on a whole protein (~1500 atoms) takes >1 h; for routine use restrict to the binding pocket (a DistanceSelector `within` selection) or use gfnff.

Reference: https://github.com/grimme-lab/xtb
"""

import os
from typing import Dict, List, Any, Optional, Union

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


class XTB(BaseConfig):
    """
    XTB: GFN2-xTB interaction-energy scoring for protein-ligand complexes.

    Inputs:
        structures: PDB complexes (each must contain the ligand as HETATM
                    records with a 3-letter residue code matching the ligand).
        ligand: compounds stream (Ligand(code="LIG") or any compounds-producing
                tool) naming the ligand to isolate from each complex. The
                residue `code` is read from the stream's `code` column at
                runtime (a single distinct value is required).
        method: xtb level — "gfn2" (default), "gfn1", or "gfnff" (faster, less accurate).
        solvent: ALPB implicit solvent identifier (e.g. "water"). None = vacuum.
        charge: total complex charge (default 0). The ligand charge is
                inferred from the ligand fragment; protein charge is
                complex_charge - ligand_charge.
        opt: if True, run a constrained geometry optimisation before the
             single-point energies; default False (single-point only).

    Outputs:
        Streams: (none)
        Tables:
            interaction_energies: id | e_complex_kj | e_protein_kj | e_ligand_kj
                                | e_interaction_kj | e_interaction_kcal | charge_complex
            missing:              id | removed_by | kind | cause
    """

    TOOL_NAME = "XTB"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        biopipelines = folders.get("biopipelines", "")
        env_check = cls._env_exists_check("xtb", env_manager)
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_check}; then
    echo "XTB already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block("xtb", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("xtb", env_manager, biopipelines)
        return f"""echo "=== Installing XTB ==="
{skip}{remove_block}
{env_block}

if {env_manager} run -n xtb xtb --version >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== XTB installation complete ==="
else
    echo "ERROR: XTB verification failed (xtb binary not on PATH)"
    exit 1
fi
"""

    structures_json = Path(lambda self: self.configuration_path("structures.json"))
    ligand_json = Path(lambda self: self.configuration_path("ligand.json"))
    interaction_csv = Path(lambda self: self.table_path("interaction_energies"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_xtb.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 ligand: Union[DataStream, StandardizedOutput],
                 method: str = "gfn2",
                 solvent: Optional[str] = None,
                 charge: int = 0,
                 opt: bool = False,
                 **kwargs):
        self.structures = structures
        self.ligand = ligand
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        if isinstance(ligand, StandardizedOutput):
            self.ligand_stream: DataStream = ligand.streams.compounds
        elif isinstance(ligand, DataStream):
            self.ligand_stream = ligand
        else:
            raise ValueError(f"ligand must be DataStream or StandardizedOutput, got {type(ligand)}")

        self.method = method.lower()
        self.solvent = solvent
        self.charge = int(charge)
        self.opt = bool(opt)
        super().__init__(**kwargs)

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        if self.method not in ("gfn2", "gfn1", "gfnff"):
            raise ValueError(f"method must be one of gfn2, gfn1, gfnff; got {self.method!r}")
        if not self.ligand_stream or len(self.ligand_stream) == 0:
            raise ValueError("ligand (compounds stream with a `code` column) is required and must not be empty")
        if self.solvent:
            _validate_freeform_string("solvent", self.solvent)

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"STRUCTURES: {len(self.structures_stream)} complexes")
        lines.append("LIGAND CODE: (resolved from ligand stream at runtime)")
        lines.append(f"METHOD: {self.method}")
        lines.append(f"SOLVENT: {self.solvent or 'vacuum'}")
        lines.append(f"CHARGE: {self.charge}")
        lines.append(f"OPT: {self.opt}")
        return lines

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        self.ligand_stream.save_json(self.ligand_json)
        ligand_json_arg = f' --ligand-json "{self.ligand_json}"'
        solvent_arg = f' --solvent "{self.solvent}"' if self.solvent else ""
        opt_arg = " --opt" if self.opt else ""
        upstream_missing_flag = self.upstream_missing_flag(self.structures, self.ligand)
        script = "#!/bin/bash\n"
        script += "# XTB GFN2-xTB interaction-energy script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += f"""echo "Scoring {len(self.structures_stream)} complex(es) with xtb ({self.method})"
python "{self.helper_py}" \\
    --structures-json "{self.structures_json}"{ligand_json_arg} \\
    --method "{self.method}" \\
    --charge {self.charge}{solvent_arg}{opt_arg} \\
    --scratch-dir "{self.extras_path()}" \\
    --interaction-csv "{self.interaction_csv}" \\
    --missing-csv "{self.missing_csv}"{upstream_missing_flag}
"""
        script += self.generate_completion_check_footer()
        return script

    def get_output_files(self) -> Dict[str, Any]:
        tables = {
            "interaction_energies": TableInfo(
                name="interaction_energies",
                path=self.interaction_csv,
                columns=["id", "e_complex_kj", "e_protein_kj", "e_ligand_kj",
                         "e_interaction_kj", "e_interaction_kcal", "charge_complex"],
                description="GFN2-xTB interaction energy E_complex - E_protein - E_ligand per complex",
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "kind", "cause"],
                description="IDs removed (upstream or local failure) with removal reason",
            ),
        }
        return {"tables": tables, "output_folder": self.output_folder}
