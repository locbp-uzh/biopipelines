# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
LigandAtomSelector configuration for distance-based residue selection
referenced to a SUBSET of atoms within a named ligand residue.

Use case: a chimeric ligand (e.g. dye+glutathione conjugate) sits in a single
ligand residue ("LIG"), but you want residues near only the glutathione atoms,
not the dye atoms. DistanceSelector takes the whole ligand residue as
reference; LigandAtomSelector takes a user-supplied atom-name subset.

Output schema is identical to DistanceSelector so the result drops into
LigandMPNN.redesigned= without changes.
"""

import os
from typing import Dict, List, Any, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string, resolve_table_reference
    from .file_paths import Path
    from .datastream import DataStream
    from .biopipelines_io import TableReference
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string, resolve_table_reference
    from file_paths import Path
    from datastream import DataStream
    from biopipelines_io import TableReference


class LigandAtomSelector(BaseConfig):
    """
    Distance-based residue selection referenced to a subset of atoms inside a
    named ligand residue.

    Parameters mirror DistanceSelector except that the reference is always
    ``ligand+atoms``: the user supplies both a ligand residue name and a
    PyMOL-style ``+``-joined list of atom names to use as the reference set.
    """

    TOOL_NAME = "LigandAtomSelector"
    TOOL_VERSION = "1.2"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== LigandAtomSelector ==="
echo "Uses biopipelines environment (no additional installation needed)."
touch "$INSTALL_SUCCESS"
echo "=== LigandAtomSelector ready ==="
"""

    # Lazy path descriptors
    selections_csv = Path(lambda self: self.table_path("selections"))
    structures_json = Path(lambda self: self.configuration_path(".input_structures.json"))
    ligand_atom_selector_py = Path(lambda self: self.pipe_script_path("pipe_ligand_atom_selector.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 ligand: str,
                 atoms: str,
                 distance: float = 5.0,
                 restrict_to: Union[str, tuple, None] = None,
                 include_reference: bool = True,
                 **kwargs):
        """
        Initialize LigandAtomSelector configuration.

        Args:
            structures: Input structures as DataStream or StandardizedOutput
            ligand: Ligand residue name (e.g. "LIG")
            atoms: ``+``-joined list of atom names within the ligand residue
                   to use as the distance reference set (e.g. "C61+C62+S57+O49").
                   All names must exist in the ligand residue of every input
                   structure or that structure is skipped.
            distance: Distance cutoff in Angstroms (default: 5.0)
            restrict_to: Optional selection to restrict the distance search to.
                         Accepts:
                         - Table reference tuple: (table, "column")
                         - Direct selection string: "10-20+30-40"
                         - None: Consider all protein residues (default)
            include_reference: Whether to include the ligand residue itself in
                               the "within" selection. Ligands are not protein
                               residues, so this has no effect on protein-residue
                               output (kept for API symmetry with DistanceSelector).
            **kwargs: Additional BaseConfig parameters.

        Output:
            Streams: (none)
            Tables:
                selections: id | pdb | within | beyond | distance_cutoff | reference_ligand
        """
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        self.ligand = ligand
        self.atoms = atoms
        self.distance = distance
        self.restrict_to_selection = resolve_table_reference(restrict_to, "restrict_to")
        self.include_reference = include_reference

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate LigandAtomSelector-specific parameters."""
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")

        if not self.ligand:
            raise ValueError("ligand parameter is required")

        if not self.atoms:
            raise ValueError("atoms parameter is required (e.g. 'C61+S57+O49')")

        if self.distance <= 0:
            raise ValueError("distance must be positive")

        if self.restrict_to_selection is not None:
            if not isinstance(self.restrict_to_selection, (str, TableReference)):
                raise ValueError("restrict_to_selection must be a string, TableReference, or None")

        _validate_freeform_string("ligand", self.ligand)
        _validate_freeform_string("atoms", self.atoms)

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input structures."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get LigandAtomSelector configuration display lines."""
        config_lines = super().get_config_display()
        atom_list = self.atoms.split("+")
        config_lines.extend([
            f"INPUT STRUCTURES: {len(self.structures_stream)} files",
            f"LIGAND: {self.ligand}",
            f"REFERENCE ATOMS ({len(atom_list)}): {self.atoms}",
            f"DISTANCE: {self.distance}Å",
            f"INCLUDE REFERENCE: {self.include_reference}",
        ])

        if self.restrict_to_selection is not None:
            if isinstance(self.restrict_to_selection, tuple):
                _, column = self.restrict_to_selection
                config_lines.append(f"RESTRICT TO: {column} from table")
            else:
                config_lines.append(f"RESTRICT TO: {self.restrict_to_selection}")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate execution script."""
        self.structures_stream.save_json(self.structures_json)

        if self.restrict_to_selection is not None:
            restrict_spec = self.restrict_to_selection
        else:
            restrict_spec = ""

        restrict_echo = f'echo "Restricting to selection: {restrict_spec}"' if restrict_spec else ""
        include_reference_str = "true" if self.include_reference else "false"

        script_content = "#!/bin/bash\n"
        script_content += "# LigandAtomSelector execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += f"""echo "Analyzing residue distances for {len(self.structures_stream)} structures"
echo "Ligand: {self.ligand}"
echo "Reference atoms: {self.atoms}"
echo "Distance cutoff: {self.distance}Å"
echo "Include reference: {self.include_reference}"
{restrict_echo}

# Run distance analysis with atom-subset reference
python {self.ligand_atom_selector_py} "{self.structures_json}" "{self.ligand}" "{self.atoms}" {self.distance} "{restrict_spec}" "{self.selections_csv}" {include_reference_str}

echo "Distance analysis completed"
echo "Selections saved to: {self.selections_csv}"

"""
        script_content += self.generate_completion_check_footer()
        return script_content

    def get_output_files(self) -> Dict[str, Any]:
        """Declare expected outputs."""
        tables = {
            "selections": TableInfo(
                name="selections",
                path=self.selections_csv,
                columns=["id", "pdb", "within", "beyond", "distance_cutoff", "reference_ligand"],
                description="PyMOL-formatted residue selections based on distance to a ligand-atom subset"
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
            "ligand_atom_selector_params": {
                "ligand": self.ligand,
                "atoms": self.atoms,
                "distance": self.distance,
                "restrict_to_selection": str(self.restrict_to_selection) if self.restrict_to_selection else None,
                "include_reference": self.include_reference
            }
        })
        return base_dict
