# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Contacts analysis for calculating contacts between selected protein regions and ligands.

Analyzes protein structures to calculate minimum distances between selected protein residues
and ligands, returning contact count and normalized distance sum.
Outputs CSV with contact metrics for all structures.
"""

import os
from typing import Dict, List, Any, Optional, Tuple, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .biopipelines_io import Resolve, TableReference
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from biopipelines_io import Resolve, TableReference


class Contacts(BaseConfig):
    """
    Pipeline tool for analyzing contacts between selected protein regions and ligands.

    Takes structures as input and calculates minimum distances between selected protein residues
    and ligands. For each residue, finds minimum distance to any ligand atom. Returns
    contact count (residues below threshold) and sum of distances normalized by √(number of residues).

    Generates CSV with contact metrics for all structures.

    Commonly used for:
    - Ligand-binding site contact analysis
    - Selected region-ligand coupling assessment
    - Binding affinity prediction features
    - Contact interface characterization
    """

    # Tool identification
    TOOL_NAME = "Contacts"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        # Contacts runs in ProteinEnv (PyMOL). Delegate the install so a notebook
        # calling Contacts.install() doesn't have to know it depends on PyMOL.
        from .pymol import PyMOL
        return PyMOL._install_script(folders, env_manager=env_manager,
                                     force_reinstall=force_reinstall, **kwargs)

    # Lazy path descriptors
    analysis_csv = Path(lambda self: self.table_path("contacts"))
    config_file = Path(lambda self: self.configuration_path("protein_ligand_config.json"))
    ligand_json = Path(lambda self: self.configuration_path(".ligand_compounds.json"))
    structures_ds_json = Path(lambda self: self.configuration_path("structures.json"))
    contacts_py = Path(lambda self: self.pipe_script_path("pipe_contacts.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 selections: Optional[Union[str, Tuple['TableInfo', str]]] = None,
                 ligand: Union[DataStream, StandardizedOutput, None] = None,
                 reference: Union[str, Tuple['TableInfo', str]] = "ligand",
                 contact_threshold: float = 5.0,
                 contact_metric_name: str = None,
                 **kwargs):
        """
        Initialize protein-ligand contact analysis tool.

        Args:
            structures: Input structures as DataStream or StandardizedOutput
            selections: Protein region selections. Accepts:
                       - None: analyze all protein residues (default)
                       - String: '10-20+30-40' (residue ranges, same for all structures)
                       - Table column reference (table, "column"): per-structure
                         selections, resolved by ID match at runtime. A single-row
                         table broadcasts to all structures; otherwise each structure
                         is matched by ID and any structure absent from the table is
                         skipped.
            ligand: Compounds stream naming the ligand whose contacts are counted.
                Required only when ``reference="ligand"``. The residue ``code`` is
                read from the stream's map_table at runtime (Ligand Contract).
            reference: What contacts are counted against. Accepts:
                - ``"ligand"`` (default): the ligand resolved from the ``ligand``
                  stream.
                - String residue selection accepted by
                  ``pdb_parser.resolve_selection`` (e.g. ``"84-182"``,
                  ``"A84-182"``, ``"10+15+20"``), same for all structures.
                - Table column reference (table, "column"): per-structure
                  residue references, resolved by ID match at runtime (single-row
                  table broadcasts to all structures; structures absent from the
                  table are skipped). Reference residues are excluded from the
                  protein side so a residue is not counted as contacting itself.
            contact_threshold: Distance threshold for counting contacts (default: 5.0 Å)
            contact_metric_name: Custom name for contact count column (default: "contacts")
            **kwargs: Additional parameters

        Output:
            Streams: (none)
            Tables:
                contacts: id | source_structure | selections | ligand | <contact_metric> | min_distance | max_distance | mean_distance | sum_distances_sqrt_normalized
        """
        # Resolve input to DataStream
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        self.reference = reference
        self.ligand_stream: Optional[DataStream] = None
        if reference == "ligand":
            if isinstance(ligand, StandardizedOutput):
                self.ligand_stream = ligand.streams.compounds
            elif isinstance(ligand, DataStream):
                self.ligand_stream = ligand
            else:
                raise ValueError("ligand must be a Ligand/compounds DataStream or "
                                 f"StandardizedOutput when reference='ligand', got {type(ligand)}")

        self.protein_selections = selections
        self.contact_threshold = contact_threshold
        self.custom_contact_metric_name = contact_metric_name

        super().__init__(**kwargs)

    def get_contact_metric_name(self) -> str:
        """Get the contact count metric name."""
        return self.custom_contact_metric_name if self.custom_contact_metric_name else "contacts"

    def validate_params(self):
        """Validate Contacts parameters."""
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures cannot be empty")

        if isinstance(self.protein_selections, str) and not self.protein_selections:
            raise ValueError("Protein selections specification cannot be empty string (use None for all protein)")

        if self.reference == "ligand" and (self.ligand_stream is None or len(self.ligand_stream) == 0):
            raise ValueError("ligand compounds stream is empty")

        if not isinstance(self.contact_threshold, (int, float)) or self.contact_threshold <= 0:
            raise ValueError("Contact threshold must be a positive number")

        if self.protein_selections is not None and not isinstance(self.protein_selections, (str, TableReference)):
            raise ValueError("selections must be a string, (TableInfo, column) tuple, or None")

        if not isinstance(self.reference, (str, TableReference)):
            raise ValueError("reference must be 'ligand', a residue selection string, or a (TableInfo, column) tuple")

        _validate_freeform_string("contact_metric_name", self.custom_contact_metric_name)
        if isinstance(self.protein_selections, str):
            _validate_freeform_string("selections", self.protein_selections)
        if isinstance(self.reference, str) and self.reference != "ligand":
            _validate_freeform_string("reference", self.reference)

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input structures."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        if self.protein_selections is None:
            selections_display = "All protein residues"
        elif isinstance(self.protein_selections, TableReference):
            selections_display = f"Column reference: {self.protein_selections.column}"
        else:
            selections_display = str(self.protein_selections)

        if self.reference == "ligand":
            ref_display = "ligand (code from compounds stream at runtime)"
        elif isinstance(self.reference, TableReference):
            ref_display = f"residues from column {self.reference.column}"
        else:
            ref_display = f"residues {self.reference}"
        config_lines.extend([
            f"STRUCTURES: {len(self.structures_stream)} files",
            f"PROTEIN SELECTIONS: {selections_display}",
            f"REFERENCE: {ref_display}",
            f"CONTACT THRESHOLD: {self.contact_threshold} Å",
            f"CONTACT METRIC: {self.get_contact_metric_name()}"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate Contacts execution script."""
        script_content = "#!/bin/bash\n"
        script_content += "# Contacts execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self.generate_script_run_contacts()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_contacts(self) -> str:
        """Generate the protein-ligand contact analysis part of the script."""
        import json

        # Serialize structures DataStream to JSON for pipe_script to load
        self.structures_stream.save_json(self.structures_ds_json)

        # Handle protein selections input
        if self.protein_selections is None:
            selections_config = {"type": "all_protein"}
        elif isinstance(self.protein_selections, TableReference):
            selections_config = {"type": "table_column",
                                 "table_path": self.protein_selections.path,
                                 "column_name": self.protein_selections.column}
        else:
            selections_config = {"type": "fixed", "value": self.protein_selections}

        # reference is "ligand", a static residue string, or a per-structure
        # table column. The string forms flow through "reference"; the column
        # form is carried separately as "reference_table" and resolved per
        # structure in the pipe script.
        if isinstance(self.reference, TableReference):
            reference_config_value = "table_column"
            reference_table = {"table_path": self.reference.path,
                               "column_name": self.reference.column}
        else:
            reference_config_value = self.reference
            reference_table = None

        config_data = {
            "structures_json": self.structures_ds_json,
            "protein_selections": selections_config,
            "ligand_name": "",   # filled at runtime via --ligand (ligand mode only)
            "reference": reference_config_value,
            "reference_table": reference_table,
            "contact_threshold": self.contact_threshold,
            "contact_metric_name": self.get_contact_metric_name(),
            "output_csv": self.analysis_csv
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        if self.protein_selections is None:
            selections_echo = 'All protein residues'
        elif isinstance(self.protein_selections, TableReference):
            selections_echo = f'per-structure from column {self.protein_selections.column}'
        else:
            selections_echo = self.protein_selections
        if self.reference == "ligand":
            self.ligand_stream.save_json(self.ligand_json)
            lig_id = self.ligand_stream.ids[0]
            resolve_code_block = (
                f'LIG_CODE_RAW={Resolve.stream_item(self.ligand_json, lig_id, column="code")}\n'
                'LIGAND_RESN="${LIG_CODE_RAW//:/}"\n'
            )
            return f"""echo "Running contact analysis (ligand reference)"
echo "Structures: {len(self.structures_stream)}"
echo "Protein selections: {selections_echo}"
{resolve_code_block}echo "Ligand: $LIGAND_RESN"
echo "Contact threshold: {self.contact_threshold} Å"
echo "Output: {self.analysis_csv}"

python "{self.contacts_py}" --config "{self.config_file}" --ligand "$LIGAND_RESN"

"""
        if isinstance(self.reference, TableReference):
            reference_echo = f'per-structure from column {self.reference.column}'
        else:
            reference_echo = self.reference
        return f"""echo "Running contact analysis (residue reference: {reference_echo})"
echo "Structures: {len(self.structures_stream)}"
echo "Protein selections: {selections_echo}"
echo "Reference residues: {reference_echo}"
echo "Contact threshold: {self.contact_threshold} Å"
echo "Output: {self.analysis_csv}"

python "{self.contacts_py}" --config "{self.config_file}"

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after protein-ligand contact analysis."""
        tables = {
            "contacts": TableInfo(
                name="contacts",
                path=self.analysis_csv,
                columns=["id", "source_structure", "selections", "ligand",
                        self.get_contact_metric_name(), "min_distance", "max_distance",
                        "mean_distance", "sum_distances_sqrt_normalized"],
                description="Protein-ligand contact analysis: ligand contacts with selected protein regions"
            )
        }

        return {
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        if self.protein_selections is None:
            selections_str = "all_protein"
        elif isinstance(self.protein_selections, TableReference):
            selections_str = f"table_column:{self.protein_selections.column}"
        else:
            selections_str = str(self.protein_selections)
        if isinstance(self.reference, TableReference):
            reference_str = f"table_column:{self.reference.column}"
        else:
            reference_str = str(self.reference)
        base_dict.update({
            "tool_params": {
                "protein_selections": selections_str,
                "reference": reference_str,
                "contact_threshold": self.contact_threshold,
                "contact_metric_name": self.get_contact_metric_name()
            }
        })
        return base_dict