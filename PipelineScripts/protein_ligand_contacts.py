"""
ProteinLigandContacts analysis for calculating contacts between selected protein regions and ligands.

Analyzes protein structures to calculate minimum distances between selected protein residues
and ligands, returning contact count and normalized distance sum.
Outputs CSV with contact metrics for all structures.
"""

import numpy as np
from typing import Dict, List, Any, Optional, Union, Tuple

import os

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo


class ProteinLigandContacts(BaseConfig):
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
    TOOL_NAME = "ProteinLigandContacts"
    DEFAULT_ENV = "ProteinEnv"
    DEFAULT_RESOURCES = {"gpu": "T4", "memory": "4GB", "time": "1:00:00"}

    def __init__(self,
                 structures: Union[ToolOutput, StandardizedOutput],
                 selections: Union[str, ToolOutput] = None,
                 ligand: str = None,
                 contact_threshold: float = 5.0,
                 contact_metric_name: str = None,
                 distance_metric_name: str = None,
                 **kwargs):
        """
        Initialize protein-ligand contact analysis tool.

        Args:
            structures: Input structures from previous tool (ToolOutput or StandardizedOutput)
            selections: Protein region selections - string, datasheet reference, or None for all protein
                       String format: '10-20+30-40' (residue ranges)
                       Datasheet format: <tool_output>.datasheets.<datasheet_name>.<column_name>
                       None: analyze all protein residues (default)
            ligand: Ligand residue name (3-letter code, e.g., 'LIG', 'ATP', 'GDP')
            contact_threshold: Distance threshold for counting contacts (default: 5.0 Å)
            contact_metric_name: Custom name for contact count column (default: "contacts")
            distance_metric_name: Custom name for distance column (default: "protein_ligand_distance")
            **kwargs: Additional parameters

        Selection Syntax:
            String format:
            - '10-20' → residues 10 to 20
            - '10-20+30-40' → residues 10-20 and 30-40
            - '145+147+150' → specific residues 145, 147, and 150

            Datasheet format:
            - rfdaa.datasheets.results.designed → use 'designed' column values
            - None → all protein residues (excluding ligands)

        Examples:
            # Analyze protein-ligand contacts with specific protein regions
            contact_analysis = pipeline.add(ProteinLigandContacts(
                structures=boltz_holo,
                selections='10-20+30-40',
                ligand='LIG',
                contact_threshold=4.0,
                contact_metric_name='close_contacts',
                distance_metric_name='binding_distance'
            ))

            # Use selections from datasheet (e.g., designed residues)
            contact_analysis = pipeline.add(ProteinLigandContacts(
                structures=boltz_holo,
                selections=rfdaa.datasheets.results.designed,
                ligand='ATP',
                contact_threshold=5.0
            ))

            # Analyze all protein residues (if selections not specified)
            contact_analysis = pipeline.add(ProteinLigandContacts(
                structures=boltz_holo,
                ligand='GDP'
            ))
        """
        self.structures_input = structures
        self.protein_selections = selections
        self.ligand_name = ligand
        self.contact_threshold = contact_threshold
        self.custom_contact_metric_name = contact_metric_name
        self.custom_distance_metric_name = distance_metric_name

        # Initialize base class
        super().__init__(**kwargs)

        # Set up dependencies
        if hasattr(structures, 'config'):
            self.dependencies.append(structures.config)
        if hasattr(selections, 'config'):
            self.dependencies.append(selections.config)

    def get_contact_metric_name(self) -> str:
        """Get the contact count metric name."""
        if self.custom_contact_metric_name:
            return self.custom_contact_metric_name
        return "contacts"

    def get_distance_metric_name(self) -> str:
        """Get the distance metric name."""
        if self.custom_distance_metric_name:
            return self.custom_distance_metric_name
        return "protein_ligand_distance"

    def get_analysis_csv_path(self) -> str:
        """Get the path for the analysis CSV file - defined once, used everywhere."""
        return os.path.join(self.output_folder, "protein_ligand_contacts.csv")

    def validate_params(self):
        """Validate ProteinLigandContacts parameters."""
        if not isinstance(self.structures_input, (ToolOutput, StandardizedOutput)):
            raise ValueError("Structures input must be a ToolOutput or StandardizedOutput object")

        # selections can be None (analyze all protein)
        if self.protein_selections is not None and not self.protein_selections:
            raise ValueError("Protein selections specification cannot be empty string (use None for all protein)")

        if not self.ligand_name:
            raise ValueError("Ligand name cannot be empty")

        if not isinstance(self.contact_threshold, (int, float)) or self.contact_threshold <= 0:
            raise ValueError("Contact threshold must be a positive number")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input structures from previous tool."""
        self.folders = pipeline_folders

        # Predict input structures paths
        self.input_structures = []
        if hasattr(self.structures_input, 'structures'):
            if isinstance(self.structures_input.structures, list):
                self.input_structures = self.structures_input.structures
            else:
                self.input_structures = [self.structures_input.structures]
        elif hasattr(self.structures_input, 'output_folder'):
            # Predict structure files in output folder
            output_folder = self.structures_input.output_folder
            predicted_structures = [
                os.path.join(output_folder, "predicted_structures.pdb"),
                os.path.join(output_folder, "structures.pdb"),
                os.path.join(output_folder, "output.pdb")
            ]
            self.input_structures = predicted_structures

        if not self.input_structures:
            raise ValueError(f"Could not predict input structure paths from: {self.structures_input}")

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        selections_display = "All protein residues" if self.protein_selections is None else str(self.protein_selections)
        if hasattr(self.protein_selections, 'datasheet_path'):
            selections_display = f"From datasheet: {self.protein_selections.datasheet_path}"

        config_lines.extend([
            f"STRUCTURES: {len(getattr(self, 'input_structures', []))} files",
            f"PROTEIN SELECTIONS: {selections_display}",
            f"LIGAND: {self.ligand_name}",
            f"CONTACT THRESHOLD: {self.contact_threshold} Å",
            f"CONTACT METRIC: {self.get_contact_metric_name()}",
            f"DISTANCE METRIC: {self.get_distance_metric_name()}"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """
        Generate protein-ligand contact analysis execution script.

        Args:
            script_path: Path where script should be written

        Returns:
            Script content as string
        """
        runtime_folder = os.path.dirname(script_path)
        output_folder = self.output_folder
        os.makedirs(output_folder, exist_ok=True)

        # Output CSV path - defined once in get_analysis_csv_path()
        analysis_csv = self.get_analysis_csv_path()

        # Create config file for protein-ligand contact calculation
        config_file = os.path.join(output_folder, "protein_ligand_config.json")

        # Handle protein selections input
        selections_config = None
        if self.protein_selections is None:
            selections_config = {"type": "all_protein"}
        elif isinstance(self.protein_selections, str):
            selections_config = {"type": "fixed", "value": self.protein_selections}
        elif isinstance(self.protein_selections, tuple) and len(self.protein_selections) == 2:
            # Tuple format: (DatasheetInfo, column_name)
            datasheet_info, column_name = self.protein_selections
            selections_config = {
                "type": "datasheet",
                "datasheet_path": datasheet_info.path if hasattr(datasheet_info, 'path') else '',
                "column_name": column_name
            }
        else:
            # Datasheet reference with attributes
            selections_config = {
                "type": "datasheet",
                "datasheet_path": getattr(self.protein_selections, 'datasheet_path', ''),
                "column_name": getattr(self.protein_selections, 'column_name', 'selections')
            }

        config_data = {
            "input_structures": self.input_structures,
            "protein_selections": selections_config,
            "ligand_name": self.ligand_name,
            "contact_threshold": self.contact_threshold,
            "contact_metric_name": self.get_contact_metric_name(),
            "distance_metric_name": self.get_distance_metric_name(),
            "output_csv": analysis_csv
        }

        import json
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        # Generate script content
        script_content = f"""#!/bin/bash
# ProteinLigandContacts execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Running protein-ligand contact analysis"
echo "Structures: {len(self.input_structures)}"
echo "Protein selections: {self.protein_selections if self.protein_selections is not None else 'All protein residues'}"
echo "Ligand: {self.ligand_name}"
echo "Contact threshold: {self.contact_threshold} Å"
echo "Output: {analysis_csv}"

# Run Python analysis script
python "{os.path.join(self.folders['HelpScripts'], 'pipe_protein_ligand_contacts.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Protein-ligand contact analysis completed successfully"
    echo "Results written to: {analysis_csv}"
else
    echo "Error: Protein-ligand contact analysis failed"
    exit 1
fi

{self.generate_completion_check_footer()}
"""

        return script_content

    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after protein-ligand contact analysis.

        Returns:
            Dictionary with output file paths
        """
        analysis_csv = self.get_analysis_csv_path()

        datasheets = {
            "contact_analysis": DatasheetInfo(
                name="contact_analysis",
                path=analysis_csv,
                columns=["id", "source_structure", "selections", "ligand",
                        self.get_contact_metric_name(), self.get_distance_metric_name()],
                description=f"Protein-ligand contact analysis: {self.ligand_name} contacts with selected protein regions",
                count=len(getattr(self, 'input_structures', []))
            )
        }

        return {
            "structures": [],
            "structure_ids": [],
            "compounds": [],
            "compound_ids": [],
            "sequences": [],
            "sequence_ids": [],
            "datasheets": datasheets,
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
                "contact_metric_name": self.get_contact_metric_name(),
                "distance_metric_name": self.get_distance_metric_name()
            }
        })
        return base_dict