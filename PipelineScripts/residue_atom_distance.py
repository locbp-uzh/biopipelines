"""
ResidueAtomDistance analysis for calculating distances between atoms and residues.

Analyzes protein structures to calculate distances between specific atoms and residues,
commonly used for ligand-protein binding site analysis and interaction validation.
Outputs CSV with distance metrics for all structures.
"""

import numpy as np
from typing import Dict, List, Any, Optional, Union, Tuple

import os

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo


class ResidueAtomDistance(BaseConfig):
    """
    Pipeline tool for analyzing structures to calculate distances between atoms and residues.
    
    Takes structures as input and outputs CSV with distance metrics for all structures.
    
    Generates CSV with distance metrics for all input structures.
    
    Commonly used for:
    - Ligand-protein binding site analysis
    - Metal coordination analysis  
    - Catalytic site geometry verification
    - Interaction distance measurements
    """
    
    # Tool identification
    TOOL_NAME = "ResidueAtomDistance"
    DEFAULT_ENV = None  # Loaded from config.yaml
    
    def __init__(self,
                 structures: Union[ToolOutput, StandardizedOutput],
                 atom: Union[str, List[str], None] = None,
                 residue: Union[str, List[str], None] = None,
                 method: str = "min",
                 metric_name: str = None,
                 **kwargs):
        """
        Initialize atom-residue distance analysis tool.

        Args:
            structures: Input structures from previous tool (ToolOutput or StandardizedOutput)
            atom: Atom selection string, list of two selections, or None
            residue: Residue selection string, list of two selections, or None
            method: How to calculate distance ("min", "max", "mean", "closest")
            metric_name: Custom name for the distance column (default: "distance")
            **kwargs: Additional parameters

        Selection Modes:
            1. Atom-Residue mode: atom=str, residue=str
            2. Atom-Atom mode: atom=[str, str], residue=None
            3. Residue-Residue mode: atom=None, residue=[str, str]

        Selection Syntax:
            Atom selections:
            - 'LIG.Cl' → ligand chlorine atoms
            - 'HAL.Br' → halogen bromine atoms
            - 'name CA' → all alpha carbon atoms

            Residue selections:
            - 'D in IGDWG' → aspartic acid in sequence context
            - '145' → residue number 145
            - '145-150' → residue range 145 to 150
            - '145+147+150' → specific residues 145, 147, and 150
            - '-1' → last residue (C-terminus)
            - '-2' → second-to-last residue
            - '1' → first residue (N-terminus)

        Examples:
            # Atom-Residue: ligand chlorine distance to specific aspartic acids
            distance_analysis = ResidueAtomDistance(
                structures=boltz_results,
                atom='LIG.Cl',
                residue='D in IGDWG',
                metric_name='chlorine_distance'
            )

            # Residue-Residue: distance between N and C termini
            termini_distance = ResidueAtomDistance(
                structures=af_results,
                atom=None,
                residue=['1', '-1'],  # First and last residue
                metric_name='termini_distance'
            )

            # Atom-Atom: distance between two ligand atoms
            ligand_internal = ResidueAtomDistance(
                structures=structure_tool,
                atom=['LIG.Cl', 'LIG.Br'],
                residue=None,
                metric_name='cl_br_distance'
            )
        """
        self.distance_input = structures
        self.atom_selection = atom
        self.residue_selection = residue
        self.distance_metric = method
        self.custom_metric_name = metric_name

        # Validate selection mode
        if atom is None and residue is None:
            raise ValueError("Both atom and residue cannot be None")

        # Validate that if one is a list, the other must be None
        if isinstance(atom, list) and residue is not None:
            raise ValueError("If atom is a list, residue must be None")
        if isinstance(residue, list) and atom is not None:
            raise ValueError("If residue is a list, atom must be None")

        # Validate list lengths
        if isinstance(atom, list) and len(atom) != 2:
            raise ValueError("atom list must contain exactly 2 selections")
        if isinstance(residue, list) and len(residue) != 2:
            raise ValueError("residue list must contain exactly 2 selections")

        # Validate distance metric
        if method not in ["min", "max", "mean", "closest"]:
            raise ValueError(f"Invalid method: {method}. Options: min, max, mean, closest")

        # Initialize base class
        super().__init__(**kwargs)

        # Set up dependency
        if hasattr(structures, 'config'):
            self.dependencies.append(structures.config)


    def get_metric_name(self) -> str:
        """Get the default metric name.""" 
        if self.custom_metric_name:
            return self.custom_metric_name
        return "distance"
    
    def get_analysis_csv_path(self) -> str:
        """Get the path for the analysis CSV file - defined once, used everywhere."""
        return os.path.join(self.output_folder, "analysis.csv")
    
    def validate_params(self):
        """Validate ResidueAtomDistance parameters."""
        if not isinstance(self.distance_input, (ToolOutput, StandardizedOutput)):
            raise ValueError("Input must be a ToolOutput or StandardizedOutput object")

        # At least one of atom or residue must be specified
        if self.atom_selection is None and self.residue_selection is None:
            raise ValueError("At least one of atom or residue must be specified")
    
    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input structures from previous tool."""
        self.folders = pipeline_folders
        
        # Predict input structures paths
        self.input_structures = []
        if hasattr(self.distance_input, 'structures'):
            if isinstance(self.distance_input.structures, list):
                self.input_structures = self.distance_input.structures
            else:
                self.input_structures = [self.distance_input.structures]
        elif hasattr(self.distance_input, 'output_folder'):
            # Predict structure files in output folder (don't check existence)
            output_folder = self.distance_input.output_folder
            # Predict common structure file patterns that tools would generate
            predicted_structures = [
                os.path.join(output_folder, "predicted_structures.pdb"),
                os.path.join(output_folder, "structures.pdb"),
                os.path.join(output_folder, "output.pdb")
            ]
            self.input_structures = predicted_structures
        
        if not self.input_structures:
            raise ValueError(f"Could not predict input structure paths from: {self.distance_input}")
    
    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        # Determine mode
        if isinstance(self.atom_selection, list):
            mode = "Atom-Atom"
            config_lines.append(f"MODE: {mode}")
            config_lines.append(f"ATOM 1: {self.atom_selection[0]}")
            config_lines.append(f"ATOM 2: {self.atom_selection[1]}")
        elif isinstance(self.residue_selection, list):
            mode = "Residue-Residue"
            config_lines.append(f"MODE: {mode}")
            config_lines.append(f"RESIDUE 1: {self.residue_selection[0]}")
            config_lines.append(f"RESIDUE 2: {self.residue_selection[1]}")
        else:
            mode = "Atom-Residue"
            config_lines.append(f"MODE: {mode}")
            config_lines.append(f"ATOM SELECTION: {self.atom_selection}")
            config_lines.append(f"RESIDUE SELECTION: {self.residue_selection}")

        config_lines.extend([
            f"DISTANCE METRIC: {self.distance_metric}",
            f"OUTPUT METRIC: {self.get_metric_name()}"
        ])

        return config_lines
    
    def generate_script(self, script_path: str) -> str:
        """
        Generate distance analysis execution script.
        
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
        
        # Create config file for distance calculation
        config_file = os.path.join(output_folder, "distance_config.json")
        config_data = {
            "input_structures": self.input_structures,
            "atom_selection": self.atom_selection,
            "residue_selection": self.residue_selection,
            "distance_metric": self.distance_metric,
            "metric_name": self.get_metric_name(),
            "output_csv": analysis_csv
        }
        
        import json
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)
        
        # Generate script content
        script_content = f"""#!/bin/bash
# ResidueAtomDistance execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Running distance analysis"
echo "Atom selection: {self.atom_selection}"
echo "Residue selection: {self.residue_selection}"
echo "Distance metric: {self.distance_metric}"
echo "Output: {analysis_csv}"

# Run Python analysis script
python "{os.path.join(self.folders['HelpScripts'], 'pipe_residue_atom_distance.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Distance analysis completed successfully"
    echo "Results written to: {analysis_csv}"
else
    echo "Error: Distance analysis failed"
    exit 1
fi

{self.generate_completion_check_footer()}
"""
        
        return script_content
    
    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after distance analysis.
        
        Returns:
            Dictionary with output file paths
        """
        analysis_csv = self.get_analysis_csv_path()
        
        tables = {
            "analysis": TableInfo(
                name="analysis", 
                path=analysis_csv,
                columns=["id", "source_structure", self.get_metric_name()],
                description=f"Distance analysis: {self.atom_selection} to {self.residue_selection}",
                count=len(self.input_structures) if hasattr(self, 'input_structures') else 0
            )
        }
        
        return {
            "structures": [],
            "structure_ids": [],
            "compounds": [],
            "compound_ids": [],
            "sequences": [],
            "sequence_ids": [],
            "tables": tables,
            "output_folder": self.output_folder
        }
    
    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()

        # Determine mode
        if isinstance(self.atom_selection, list):
            mode = "atom-atom"
        elif isinstance(self.residue_selection, list):
            mode = "residue-residue"
        else:
            mode = "atom-residue"

        base_dict.update({
            "tool_params": {
                "mode": mode,
                "atom_selection": self.atom_selection,
                "residue_selection": self.residue_selection,
                "distance_metric": self.distance_metric,
                "metric_name": self.get_metric_name()
            }
        })
        return base_dict