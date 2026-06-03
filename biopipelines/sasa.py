# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
SASA (Solvent Accessible Surface Area) tool for calculating delta SASA.

Calculates the change in solvent accessible surface area when the binder protein
is present vs absent. This metric indicates how much of the ligand surface is
buried by protein binding.

delta_SASA = SASA_ligand_alone - SASA_ligand_in_complex
"""

import os
from typing import Dict, List, Any, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .biopipelines_io import Resolve
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from biopipelines_io import Resolve


_VALID_MODES = ("ligand", "residues")


class SASA(BaseConfig):
    """
    SASA tool for solvent accessible surface area, in two modes.

    Uses PyMOL's get_area command for SASA calculations.

    mode="ligand" (default): delta_SASA = SASA(ligand alone) - SASA(ligand in
    complex). A higher delta_SASA indicates more ligand surface is buried by the
    protein, suggesting tighter binding. Requires ``ligand``; emits the ``sasa``
    table.

    mode="residues": per-residue protein accessibility, ligand-independent.
    Emits a per-residue resi-csv ``accessibility`` stream (one <id>.csv per
    structure) with columns id | chain | resi | resn | sasa | rsa, where ``rsa``
    is relative accessibility in [0,1] (sasa / Tien max-ACC per restype). Lets
    Selection threshold on exposure, e.g.
    ``Selection.add(sasa.streams.accessibility, include="rsa>=0.25")``.
    """

    TOOL_NAME = "SASA"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        # SASA runs in ProteinEnv (PyMOL). Delegate the install so a notebook
        # calling SASA.install() doesn't have to know it depends on PyMOL.
        from .pymol import PyMOL
        return PyMOL._install_script(folders, env_manager=env_manager,
                                     force_reinstall=force_reinstall, **kwargs)

    # Lazy path descriptors
    results_csv = Path(lambda self: self.table_path("sasa_analysis"))
    structures_ds_json = Path(lambda self: self.configuration_path("structures.json"))
    ligand_json = Path(lambda self: self.configuration_path(".ligand_compounds.json"))
    accessibility_map = Path(lambda self: self.stream_map_path("accessibility"))
    helper_script = Path(lambda self: self.pipe_script_path("pipe_sasa.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 ligand: Union[DataStream, StandardizedOutput, None] = None,
                 mode: str = "ligand",
                 dot_density: int = 4,
                 **kwargs):
        """
        Initialize SASA configuration.

        Args:
            structures: Input structures as DataStream or StandardizedOutput
            ligand: The bound ligand whose buried surface is measured (required
                only when ``mode="ligand"``). Pass a ``Ligand`` / any
                compounds-producing tool output — the residue ``code`` is read
                from the compounds stream's map_table at runtime, per the Ligand
                Contract.
            mode: ``"ligand"`` (default) — ligand delta-SASA. ``"residues"`` —
                per-residue protein accessibility (ligand-independent).
            dot_density: Dot density for SASA calculation (1-4, higher = more accurate)
            **kwargs: Additional parameters

        Output:
            mode="ligand":
                Tables:
                    sasa: id | structure | sasa_ligand_alone | sasa_ligand_complex | delta_sasa
            mode="residues":
                Streams:
                    accessibility: per-residue resi-csv (id | chain | resi | resn | sasa | rsa)
        """
        # Resolve input to DataStream
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        # The ligand residue code is read from the compounds stream at runtime;
        # only needed in ligand mode.
        self.ligand_stream: Union[DataStream, None] = None
        if ligand is None:
            pass
        elif isinstance(ligand, StandardizedOutput):
            self.ligand_stream = ligand.streams.compounds
        elif isinstance(ligand, DataStream):
            self.ligand_stream = ligand
        else:
            raise ValueError("ligand must be a Ligand/compounds DataStream or "
                             f"StandardizedOutput, got {type(ligand)}")

        self.mode = mode
        self.dot_density = dot_density

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate SASA parameters."""
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        if self.mode not in _VALID_MODES:
            raise ValueError(f"mode must be one of {_VALID_MODES}, got {self.mode!r}")
        if self.mode == "ligand":
            if self.ligand_stream is None:
                raise ValueError("ligand parameter is required when mode='ligand'")
            if len(self.ligand_stream) == 0:
                raise ValueError("ligand compounds stream is empty")
        if self.dot_density < 1 or self.dot_density > 4:
            raise ValueError("dot_density must be between 1 and 4")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input structures."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()
        config_lines.append(f"INPUT STRUCTURES: {len(self.structures_stream)} files")
        config_lines.append(f"MODE: {self.mode}")
        if self.mode == "ligand":
            config_lines.append("LIGAND: (code from compounds stream at runtime)")
        config_lines.append(f"DOT DENSITY: {self.dot_density}")
        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate SASA execution script."""
        import json

        # Serialize structures DataStream to JSON for pipe_script to load
        self.structures_stream.save_json(self.structures_ds_json)

        script_content = "#!/bin/bash\n"
        script_content += "# SASA execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self._generate_script_run_sasa()
        script_content += self.generate_completion_check_footer()

        return script_content

    def _generate_script_run_sasa(self) -> str:
        """Generate the SASA execution part of the script."""
        if self.mode == "residues":
            return f"""echo "Running per-residue SASA analysis"
echo "Structures: {len(self.structures_stream)} files"

python "{self.helper_script}" \\
    --mode residues \\
    --structures "{self.structures_ds_json}" \\
    --accessibility-dir "{self.stream_folder('accessibility')}" \\
    --accessibility-map-csv "{self.accessibility_map}" \\
    --dot_density {self.dot_density}

if [ $? -eq 0 ]; then
    echo "Per-residue SASA analysis completed successfully"
else
    echo "Error: SASA analysis failed"
    exit 1
fi

"""

        # The ligand residue code is resolved at runtime from the compounds
        # stream's `code` column. Colons are stripped for the PyMOL selection.
        self.ligand_stream.save_json(self.ligand_json)
        lig_id = self.ligand_stream.ids[0]
        resolve_code_block = (
            f'LIG_CODE_RAW={Resolve.stream_item(self.ligand_json, lig_id, column="code")}\n'
            'LIGAND_RESN="${LIG_CODE_RAW//:/}"\n'
        )
        ligand_resn = "$LIGAND_RESN"
        ligand_echo = "$LIGAND_RESN"

        return f"""echo "Running SASA analysis"
echo "Structures: {len(self.structures_stream)} files"
{resolve_code_block}echo "Ligand: {ligand_echo}"
echo "Output: {self.results_csv}"

python "{self.helper_script}" \\
    --mode ligand \\
    --structures "{self.structures_ds_json}" \\
    --ligand "{ligand_resn}" \\
    --output_csv "{self.results_csv}" \\
    --dot_density {self.dot_density}

if [ $? -eq 0 ]; then
    echo "SASA analysis completed successfully"
    echo "Results written to: {self.results_csv}"
else
    echo "Error: SASA analysis failed"
    exit 1
fi

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after SASA execution."""
        if self.mode == "residues":
            accessibility_stream = DataStream(
                name="accessibility",
                ids=self.structures_stream.ids,
                files=[self.stream_path("accessibility", "<id>.csv")],
                map_table=self.accessibility_map,
                format="resi-csv",
            )
            return {
                "accessibility": accessibility_stream,
                "output_folder": self.output_folder,
            }

        tables = {
            "sasa": TableInfo(
                name="sasa",
                path=self.results_csv,
                columns=["id", "structure", "sasa_ligand_alone", "sasa_ligand_complex", "delta_sasa"],
                description="Solvent accessible surface area analysis"
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
            "sasa_params": {
                "mode": self.mode,
                "dot_density": self.dot_density
            }
        })
        return base_dict
