# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
DistanceSelector configuration for distance-based residue selection.

Analyzes protein structures to identify residues within or beyond a specified
distance from a reference (a ligand, a residue range, or any atom/atom-set
expressible in the framework's selection grammar), and emits a per-residue
distance resi-csv stream for downstream thresholding.
"""

import os
from typing import Dict, List, Any, Union, Optional

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .biopipelines_io import TableReference, Resolve
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from biopipelines_io import TableReference, Resolve


_VALID_MODES = ("min", "centroid")
_VALID_ATOM_CLASSES = ("all", "backbone", "CA", "sidechain")


class DistanceSelector(BaseConfig):
    """
    Distance-based residue selection.

    Computes the per-residue distance from each protein residue to a reference
    (ligand, residue range, or any selection grammar expression such as
    ``LIG.O132``, ``LIG.Cl+LIG.O3``, ``87.CA``, ``A141.CB``, ``first.CA``,
    ``D in IGDWG``) and partitions residues into within/beyond sets according
    to a distance cutoff, a top-K cap, or both.
    """

    TOOL_NAME = "DistanceSelector"
    TOOL_VERSION = "2.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== DistanceSelector ==="
echo "Uses biopipelines environment (no additional installation needed)."
touch "$INSTALL_SUCCESS"
echo "=== DistanceSelector ready ==="
"""

    selections_csv = Path(lambda self: self.table_path("selections"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    structures_json = Path(lambda self: self.configuration_path(".input_structures.json"))
    ligand_json = Path(lambda self: self.configuration_path(".ligand_compounds.json"))
    distances_map = Path(lambda self: self.stream_map_path("distances"))
    distance_selector_py = Path(lambda self: self.pipe_script_path("pipe_distance_selector.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 ligand: Union[DataStream, StandardizedOutput, None] = None,
                 distance: Optional[float] = 5.0,
                 reference: str = "ligand",
                 top_k: Optional[int] = None,
                 mode: str = "min",
                 atoms: str = "all",
                 restrict_to: Union[str, tuple, None] = None,
                 include_reference: bool = True,
                 **kwargs):
        """
        Args:
            structures: Input structures (DataStream / StandardizedOutput).
            ligand: Bound ligand (required only when ``reference="ligand"``).
                Pass a ``Ligand`` / any compounds-producing tool output; the
                residue ``code`` is read from the compounds stream's map_table
                at runtime (Ligand Contract).
            distance: Distance cutoff in Angstroms. Optional when ``top_k`` is
                set; required otherwise.
            reference: Reference specifier. Either the literal ``"ligand"``
                (resolves the ligand's residue code from the compounds stream
                at runtime) or any selection expression accepted by
                ``pdb_parser.resolve_selection``: ``"LIG.O132"``,
                ``"LIG.Cl+LIG.O3"``, ``"87.CA"``, ``"A141.CB"``,
                ``"first.CA"``, ``"D in IGDWG"``, ``"87-100"``,
                ``"10+15+20"``.
            top_k: Cap to the K nearest residues. Combinations:
                - distance only → all residues within cutoff (legacy);
                - top_k only → K nearest residues, no cutoff;
                - both → up to K nearest among residues within cutoff.
            mode: ``"min"`` (default) — distance = min over reference atoms
                × residue atoms. ``"centroid"`` — distance from each residue's
                atom set to the reference atom set's centroid.
            atoms: Which atoms enter the distance calc on BOTH protein and
                reference sides. One of ``"all"``, ``"backbone"`` (N, CA, C,
                O), ``"CA"``, ``"sidechain"``.
            restrict_to: Restrict the candidate residue set. Accepts:
                - residue selection string (``"10-20+30-40"``, ``"A10-20"``);
                - bare chain spec (``"B"`` or ``"chain B"``);
                - ``(TableInfo, column)`` per-structure column reference;
                - ``None`` — no restriction.
            include_reference: When the reference is itself protein residues
                (e.g. ``reference="87-100"``), whether to include them in
                ``within``. No-op for ligand/atom references.
        """
        # Keep original inputs for upstream missing-table detection.
        self.structures_input = structures
        self.ligand_input = ligand

        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

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

        self.distance = distance
        self.reference = reference
        self.top_k = top_k
        self.mode = mode
        self.atoms = atoms
        self.restrict_to_selection = restrict_to
        self.include_reference = include_reference

        super().__init__(**kwargs)

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")

        if not isinstance(self.reference, str) or not self.reference.strip():
            raise ValueError("reference must be a non-empty string")

        if self.reference == "ligand":
            if self.ligand_stream is None:
                raise ValueError("ligand parameter is required when reference is 'ligand'")
            if len(self.ligand_stream) == 0:
                raise ValueError("ligand compounds stream is empty")

        if self.distance is None and self.top_k is None:
            raise ValueError("at least one of `distance` or `top_k` must be set")
        if self.distance is not None and self.distance <= 0:
            raise ValueError("distance must be positive")
        if self.top_k is not None and (not isinstance(self.top_k, int) or self.top_k <= 0):
            raise ValueError("top_k must be a positive integer")

        if self.mode not in _VALID_MODES:
            raise ValueError(f"mode must be one of {_VALID_MODES}, got {self.mode!r}")
        if self.atoms not in _VALID_ATOM_CLASSES:
            raise ValueError(f"atoms must be one of {_VALID_ATOM_CLASSES}, got {self.atoms!r}")

        if self.restrict_to_selection is not None:
            if not isinstance(self.restrict_to_selection, (str, TableReference, tuple)):
                raise ValueError("restrict_to must be a string, (TableInfo, column) tuple, or None")
            if isinstance(self.restrict_to_selection, str):
                _validate_freeform_string("restrict_to", self.restrict_to_selection)

        _validate_freeform_string("reference", self.reference)

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"INPUT STRUCTURES: {len(self.structures_stream)} files")
        lines.append(f"REFERENCE: {self.reference}")
        if self.reference == "ligand":
            lines.append("LIGAND: (code from compounds stream at runtime)")
        if self.distance is not None:
            lines.append(f"DISTANCE: {self.distance}Å")
        if self.top_k is not None:
            lines.append(f"TOP_K: {self.top_k}")
        lines.append(f"MODE: {self.mode}")
        lines.append(f"ATOMS: {self.atoms}")
        lines.append(f"INCLUDE REFERENCE: {self.include_reference}")
        if self.restrict_to_selection is not None:
            if isinstance(self.restrict_to_selection, tuple):
                _, column = self.restrict_to_selection
                lines.append(f"RESTRICT TO: {column} from table")
            else:
                lines.append(f"RESTRICT TO: {self.restrict_to_selection}")
        return lines

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)

        # For the legacy `reference="ligand"` alias, resolve the ligand's residue
        # code from the compounds stream at runtime and inject it as the actual
        # selection string passed to the pipe script.
        resolve_code_block = ""
        if self.reference == "ligand":
            self.ligand_stream.save_json(self.ligand_json)
            lig_id = self.ligand_stream.ids[0]
            resolve_code_block = (
                f'LIG_CODE={Resolve.stream_item(self.ligand_json, lig_id, column="code")}\n'
            )
            reference_arg = "$LIG_CODE"
        else:
            reference_arg = self.reference

        restrict_spec = self.restrict_to_selection if self.restrict_to_selection is not None else ""
        restrict_echo = f'echo "Restricting to: {restrict_spec}"' if restrict_spec else ""

        include_reference_str = "true" if self.include_reference else "false"
        distance_arg = "" if self.distance is None else str(self.distance)
        top_k_arg = "" if self.top_k is None else str(self.top_k)

        script = "#!/bin/bash\n"
        script += "# DistanceSelector execution script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += f"""echo "Analyzing residue distances for {len(self.structures_stream)} structures"
{resolve_code_block}echo "Reference: {reference_arg}"
echo "Distance cutoff: {distance_arg}"
echo "Top-K: {top_k_arg}"
echo "Mode: {self.mode}"
echo "Atoms: {self.atoms}"
echo "Include reference: {self.include_reference}"
{restrict_echo}

python "{self.distance_selector_py}" \\
    --structures-json "{self.structures_json}" \\
    --reference "{reference_arg}" \\
    --distance "{distance_arg}" \\
    --top-k "{top_k_arg}" \\
    --mode "{self.mode}" \\
    --atoms "{self.atoms}" \\
    --restrict-to "{restrict_spec}" \\
    --include-reference "{include_reference_str}" \\
    --selections-csv "{self.selections_csv}" \\
    --distances-dir "{self.stream_folder('distances')}" \\
    --distances-map-csv "{self.distances_map}"

echo "Distance analysis completed"
echo "Selections saved to: {self.selections_csv}"
"""
        script += self.generate_missing_propagation(
            self.structures_input, self.ligand_input, missing_csv=self.missing_csv
        )
        script += self.generate_completion_check_footer()
        return script

    def get_output_files(self) -> Dict[str, Any]:
        distances_stream = DataStream(
            name="distances",
            ids=self.structures_stream.ids,
            files=[self.stream_path("distances", "<id>.csv")],
            map_table=self.distances_map,
            format="resi-csv",
        )

        tables = {
            "selections": TableInfo(
                name="selections",
                path=self.selections_csv,
                columns=["id", "pdb", "within", "beyond", "distance_cutoff", "top_k", "mode", "reference"],
                description="PyMOL-formatted residue selections based on distance to a reference",
            )
        }

        # Excuse upstream-filtered ids: declare a `missing` table when any input
        # axis carries one, so the completion check doesn't demand their files.
        if self._collect_upstream_missing_paths(self.structures_input, self.ligand_input):
            tables["missing"] = self.missing_table_info(self.missing_csv)

        return {
            "distances": distances_stream,
            "tables": tables,
            "output_folder": self.output_folder,
        }
