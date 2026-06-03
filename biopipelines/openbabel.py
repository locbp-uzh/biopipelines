# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""OpenBabel tool: convert a compounds or structures stream between chemical
file formats, optionally adding hydrogens (pH-aware) and generating 3D
coordinates."""

import os
from typing import Dict, List, Any, Optional, Union

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


# 3-D coordinate formats live on the `structures` stream; line/notation formats
# (no coordinates) become columns on the `compounds` stream.
THREE_D_FORMATS = {"sdf", "mol2", "mol", "pdb", "pdbqt", "xyz"}
LINE_FORMATS = {"smi", "inchi", "cml"}


class OpenBabel(BaseConfig):
    """
    OpenBabel: convert a compounds or structures stream between chemical
    file formats, optionally adding hydrogens.

    Inputs:
        compounds: source compounds stream (mutually exclusive with structures).
        structures: source structures stream (mutually exclusive with compounds).
        convert_3d: one 3-D coordinate format (sdf/mol2/mol/pdb/pdbqt/xyz) -> structures stream.
        convert_1d: line/notation formats (smi/inchi/cml) -> columns on the compounds stream.
        add_hydrogens: if True, add explicit hydrogens.
        pH: if set, add hydrogens at the given pH (implies add_hydrogens=True).
            Pass-through to `obabel -p <pH>`.
        gen3d: if True, generate 3D coordinates. Only valid with compounds input.
        gen3d_quality: 3-D embedding effort for gen3d/convert_3d, one of "fastest",
            "fast", "medium" (default), "better", "best".
        minimize: if True, run a force-field geometry minimization after embedding
            coordinates (uses `ff` and `minimize_steps`).
        ff: force field for minimization, one of "MMFF94" (default), "MMFF94s",
            "UFF", "GAFF", "Ghemical".
        minimize_steps: maximum minimization iterations (default 500).

    Outputs:
        Streams:
            - compounds input + a 3-D format (sdf/mol2/mol/pdb/pdbqt/xyz):
                emits a `structures` stream (the per-id coordinate files) AND
                passes the input `compounds` chemistry csv through unchanged.
                This is how a Ligand becomes a docking-ready 3-D ligand:
                `OpenBabel(compounds=Ligand("aspirin"), convert_3d="sdf")`.
            - compounds input + a line format (smi/inchi/cml): emits a
                `compounds` stream only (no coordinates).
            - structures input: emits a `structures` stream (one file per id).
    """

    TOOL_NAME = "OpenBabel"
    TOOL_VERSION = "1.1"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        # No-op: on the cluster/local the biopipelines conda env bundles
        # openbabel (pybel imports directly); on Colab the `colab` pip extra
        # (openbabel-wheel) puts pybel into base Python at `pip install -e
        # ".[colab]"` time. Either way nothing to install per-run here.
        return """echo "=== OpenBabel ==="
echo "Uses biopipelines environment (openbabel is bundled, no additional installation needed)."
touch "$INSTALL_SUCCESS"
echo "=== OpenBabel ready ==="
"""

    input_json = Path(lambda self: self.configuration_path("input.json"))
    coords_json = Path(lambda self: self.configuration_path("coords.json"))
    compounds_map = Path(lambda self: self.stream_map_path("compounds"))
    structures_map = Path(lambda self: self.stream_map_path("structures"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_openbabel.py"))

    def __init__(self,
                 compounds: Union[DataStream, StandardizedOutput, None] = None,
                 structures: Union[DataStream, StandardizedOutput, None] = None,
                 convert_3d: Optional[str] = None,
                 convert_1d: Optional[Union[str, List[str]]] = None,
                 add_hydrogens: bool = False,
                 pH: Optional[float] = None,
                 gen3d: bool = False,
                 gen3d_quality: str = "medium",
                 minimize: bool = False,
                 ff: str = "MMFF94",
                 minimize_steps: int = 500,
                 use_structure_template: bool = True,
                 **kwargs):
        if (compounds is None) == (structures is None):
            raise ValueError("OpenBabel requires exactly one of `compounds` or `structures`")

        if compounds is not None:
            self.input_kind = "compounds"
            src = compounds
        else:
            self.input_kind = "structures"
            src = structures

        if isinstance(src, StandardizedOutput):
            self.input_stream: DataStream = getattr(src.streams, self.input_kind)
        elif isinstance(src, DataStream):
            self.input_stream = src
        else:
            raise ValueError(f"{self.input_kind} must be DataStream or StandardizedOutput, got {type(src)}")

        # Posed-ligand path: a compounds StandardizedOutput may ALSO carry bound
        # coordinates on its structures stream (e.g. Ligand("STI", structures=
        # complex), which takes chemistry from the lookup and coords from the
        # bound HETATM). When use_structure_template is on and such coordinates
        # exist, a convert_3d="sdf" uses THOSE coordinates + the SMILES as a
        # bond-order template (no make3D) so the output SDF is posed, not a fresh
        # origin-centered embedding. self.coord_stream is None otherwise.
        self.use_structure_template = bool(use_structure_template)
        self.coord_stream = None
        if (self.input_kind == "compounds" and self.use_structure_template
                and isinstance(src, StandardizedOutput)):
            cand = getattr(src.streams, "structures", None)
            if cand is not None and len(cand) > 0:
                self.coord_stream = cand

        # convert_3d: one coordinate format -> structures stream.
        # convert_1d: line/notation formats -> columns on the compounds csv.
        self.convert_3d = convert_3d.lower() if convert_3d else None
        if convert_1d is None:
            self.convert_1d: List[str] = []
        elif isinstance(convert_1d, str):
            self.convert_1d = [convert_1d.lower()]
        else:
            self.convert_1d = [f.lower() for f in convert_1d]

        self.add_hydrogens = bool(add_hydrogens) or pH is not None
        self.pH = float(pH) if pH is not None else None
        self.gen3d = bool(gen3d)
        self.gen3d_quality = gen3d_quality.lower()
        self.minimize = bool(minimize)
        self.ff = ff
        self.minimize_steps = int(minimize_steps)
        super().__init__(**kwargs)

    _GEN3D_QUALITIES = ("fastest", "fast", "medium", "better", "best")
    _FORCEFIELDS = ("MMFF94", "MMFF94s", "UFF", "GAFF", "Ghemical")

    def validate_params(self):
        if not self.input_stream or len(self.input_stream) == 0:
            raise ValueError(f"{self.input_kind} parameter is required and must not be empty")
        if not self.convert_3d and not self.convert_1d:
            raise ValueError("OpenBabel requires at least one of `convert_3d` or `convert_1d`")
        if self.convert_3d is not None:
            if self.convert_3d not in THREE_D_FORMATS:
                raise ValueError(
                    f"convert_3d must be a 3-D format {sorted(THREE_D_FORMATS)}, got {self.convert_3d!r}"
                )
            _validate_freeform_string("convert_3d", self.convert_3d)
        for fmt in self.convert_1d:
            if fmt not in LINE_FORMATS:
                raise ValueError(
                    f"convert_1d entries must be line formats {sorted(LINE_FORMATS)}, got {fmt!r}"
                )
            _validate_freeform_string("convert_1d", fmt)
        if self.convert_1d and self.input_kind == "structures":
            raise ValueError("convert_1d is only valid with a compounds input (notations live on the compounds stream)")
        if self.gen3d and self.input_kind == "structures":
            raise ValueError("gen3d is only valid with a compounds input (structures already have 3D coordinates)")
        if self.gen3d_quality not in self._GEN3D_QUALITIES:
            raise ValueError(f"gen3d_quality must be one of {self._GEN3D_QUALITIES}, got '{self.gen3d_quality}'")
        if self.ff not in self._FORCEFIELDS:
            raise ValueError(f"ff must be one of {self._FORCEFIELDS}, got '{self.ff}'")
        if self.minimize_steps < 1:
            raise ValueError("minimize_steps must be at least 1")
        if self.minimize and not (self.gen3d or self.convert_3d):
            raise ValueError("minimize requires 3-D coordinates (set gen3d=True or convert_3d=...)")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"{self.input_kind.upper()}: {len(self.input_stream)} item(s)")
        if self.convert_3d:
            lines.append(f"CONVERT 3D: {self.convert_3d} -> structures")
        if self.convert_1d:
            lines.append(f"CONVERT 1D: {', '.join(self.convert_1d)} -> compounds columns")
        if self.pH is not None:
            lines.append(f"ADD HYDROGENS: yes (pH={self.pH})")
        elif self.add_hydrogens:
            lines.append("ADD HYDROGENS: yes")
        if self.gen3d:
            lines.append(f"GENERATE 3D: yes (quality={self.gen3d_quality})")
        if self.minimize:
            lines.append(f"MINIMIZE: {self.ff} ({self.minimize_steps} steps)")
        return lines

    def generate_script(self, script_path: str) -> str:
        self.input_stream.save_json(self.input_json)

        flags = [f"--input-kind {self.input_kind}"]
        if self.convert_3d:
            flags.append(f"--convert-3d {self.convert_3d}")
            flags.append(f"--structures-dir \"{self.stream_folder('structures')}\"")
            flags.append(f"--structures-map \"{self.structures_map}\"")
            # Posed-ligand path: hand the pipe the bound coordinate files so it
            # writes the SDF from THOSE coords (+ SMILES bond-order template),
            # not a fresh make3D embedding. Only when coords accompany compounds.
            if self.coord_stream is not None:
                self.coord_stream.save_json(self.coords_json)
                flags.append(f"--coords-json \"{self.coords_json}\"")
        if self.convert_1d:
            flags.append(f"--convert-1d {' '.join(self.convert_1d)}")
            flags.append(f"--compounds-map \"{self.compounds_map}\"")
        if self.add_hydrogens:
            flags.append("--add-hydrogens")
        if self.pH is not None:
            flags.append(f"--ph {self.pH}")
        if self.gen3d:
            flags.append("--gen3d")
        flags.append(f"--gen3d-quality {self.gen3d_quality}")
        if self.minimize:
            flags.append("--minimize")
            flags.append(f"--ff {self.ff}")
            flags.append(f"--minimize-steps {self.minimize_steps}")
        flag_str = " ".join(flags)

        targets = []
        if self.convert_3d:
            targets.append(f"3D:{self.convert_3d}")
        if self.convert_1d:
            targets.append(f"1D:{','.join(self.convert_1d)}")

        script = "#!/bin/bash\n"
        script += "# OpenBabel conversion script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += f"""echo "Converting {len(self.input_stream)} {self.input_kind} item(s) [{' '.join(targets)}]"
python "{self.helper_py}" \\
    --input-json "{self.input_json}" \\
    {flag_str}
"""
        script += self.generate_completion_check_footer()
        return script

    def get_output_files(self) -> Dict[str, Any]:
        result = {
            "compounds": DataStream.empty("compounds", "csv"),
            "structures": DataStream.empty("structures", "pdb"),
            "sequences": DataStream.empty("sequences", "fasta"),
            "tables": {},
            "output_folder": self.output_folder,
        }

        # 3-D conversion -> structures stream (per-id coordinate files).
        if self.convert_3d:
            result["structures"] = DataStream(
                name="structures",
                ids=self.input_stream.ids,
                files=[self.stream_path("structures", f"<id>.{self.convert_3d}")],
                map_table=self.structures_map,
                format=self.convert_3d,
            )

        # 1-D conversion -> compounds csv with one extra column per notation.
        # The chemistry passthrough is also emitted whenever the input is a
        # compounds stream so the SMILES/code travel with the 3-D structures.
        if self.input_kind == "compounds":
            result["compounds"] = DataStream(
                name="compounds",
                ids=self.input_stream.ids,
                files=[],
                map_table=self.compounds_map if self.convert_1d else self.input_stream.map_table,
                format="csv",
            )

        return result
