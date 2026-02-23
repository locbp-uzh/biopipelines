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
from typing import Dict, List, Any, Optional, Union

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

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== Contacts ==="
echo "Requires ProteinEnv (installed with PyMOL.install())"
echo "No additional installation needed."
echo "=== Contacts ready ==="
"""

    # Lazy path descriptors
    analysis_csv = Path(lambda self: os.path.join(self.output_folder, "contacts.csv"))
    config_file = Path(lambda self: os.path.join(self.output_folder, "protein_ligand_config.json"))
    structures_ds_json = Path(lambda self: os.path.join(self.output_folder, "structures.json"))
    contacts_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_contacts.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 selections: Optional[str] = None,
                 ligand: str = None,
                 contact_threshold: float = 5.0,
                 contact_metric_name: str = None,
                 id_map: Dict[str, str] = {"*": "*_<S>"},
                 **kwargs):
        """
        Initialize protein-ligand contact analysis tool.

        Args:
            structures: Input structures as DataStream or StandardizedOutput
            selections: Protein region selections string or None for all protein
                       String format: '10-20+30-40' (residue ranges)
                       None: analyze all protein residues (default)
            ligand: Ligand residue name (3-letter code, e.g., 'LIG', 'ATP', 'GDP')
            contact_threshold: Distance threshold for counting contacts (default: 5.0 Å)
            contact_metric_name: Custom name for contact count column (default: "contacts")
            id_map: ID mapping pattern for matching structure IDs to table IDs
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

        self.protein_selections = selections
        self.ligand_name = ligand
        self.contact_threshold = contact_threshold
        self.custom_contact_metric_name = contact_metric_name
        self.id_map = id_map

        super().__init__(**kwargs)

    def get_contact_metric_name(self) -> str:
        """Get the contact count metric name."""
        return self.custom_contact_metric_name if self.custom_contact_metric_name else "contacts"

    def validate_params(self):
        """Validate Contacts parameters."""
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures cannot be empty")

        if self.protein_selections is not None and not self.protein_selections:
            raise ValueError("Protein selections specification cannot be empty string (use None for all protein)")

        if not self.ligand_name:
            raise ValueError("Ligand name cannot be empty")

        if not isinstance(self.contact_threshold, (int, float)) or self.contact_threshold <= 0:
            raise ValueError("Contact threshold must be a positive number")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input structures."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        selections_display = "All protein residues" if self.protein_selections is None else str(self.protein_selections)

        config_lines.extend([
            f"STRUCTURES: {len(self.structures_stream)} files",
            f"PROTEIN SELECTIONS: {selections_display}",
            f"LIGAND: {self.ligand_name}",
            f"CONTACT THRESHOLD: {self.contact_threshold} Å",
            f"CONTACT METRIC: {self.get_contact_metric_name()}"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate Contacts execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

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

        # Serialize structures DataStream to JSON for HelpScript to load
        with open(self.structures_ds_json, 'w') as f:
            json.dump(self.structures_stream.to_dict(), f, indent=2)

        # Handle protein selections input
        if self.protein_selections is None:
            selections_config = {"type": "all_protein"}
        else:
            selections_config = {"type": "fixed", "value": self.protein_selections}

        config_data = {
            "structures_json": self.structures_ds_json,
            "protein_selections": selections_config,
            "ligand_name": self.ligand_name,
            "contact_threshold": self.contact_threshold,
            "contact_metric_name": self.get_contact_metric_name(),
            "id_map": self.id_map,
            "output_csv": self.analysis_csv
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        return f"""echo "Running protein-ligand contact analysis"
echo "Structures: {len(self.structures_stream)}"
echo "Protein selections: {self.protein_selections if self.protein_selections is not None else 'All protein residues'}"
echo "Ligand: {self.ligand_name}"
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
                description=f"Protein-ligand contact analysis: {self.ligand_name} contacts with selected protein regions",
                count=len(self.structures_stream)
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
                "protein_selections": str(self.protein_selections) if self.protein_selections is not None else "all_protein",
                "ligand_name": self.ligand_name,
                "contact_threshold": self.contact_threshold,
                "contact_metric_name": self.get_contact_metric_name()
            }
        })
        return base_dict