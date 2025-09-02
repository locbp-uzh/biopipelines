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
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput


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
    DEFAULT_ENV = "ProteinEnv" 
    COMPATIBLE_ENVS = ["ProteinEnv", "Boltz2Env"]
    DEFAULT_RESOURCES = {"gpu": "T4", "memory": "4GB", "time": "1:00:00"}
    
    def __init__(self,
                 input: Union[ToolOutput, StandardizedOutput],
                 atom: str,
                 residue: str,
                 method: str = "min",
                 metric_name: str = None,
                 **kwargs):
        """
        Initialize atom-residue distance analysis tool.
        
        Args:
            input: Input structures from previous tool (ToolOutput or StandardizedOutput)
            atom: Atom selection string (e.g., 'LIG.Cl', 'protein.CA', 'resname ZN')  
            residue: Residue selection string (e.g., 'protein.D in TRGDTGH', 'resid 145-150')
            method: How to calculate distance ("min", "max", "mean", "closest")
            metric_name: Custom name for the distance column (default: "distance")
            **kwargs: Additional parameters
            
        Examples:
            # Analyze ligand chlorine distance to specific aspartic acids
            distance_analysis = pipeline.add(ResidueAtomDistance(
                input=boltz_results.output,
                atom='LIG.Cl',
                residue='protein.D in TRGDTGH',
                metric_name='chlorine_distance'
            ))
            
            # Analyze protein CA atoms distance to ligand
            ca_analysis = pipeline.add(ResidueAtomDistance(
                input=structure_tool.output,
                atom='protein.CA',
                residue='ligand.all',
                distance_metric='mean'
            ))
        """
        self.distance_input = input
        self.atom_selection = atom
        self.residue_selection = residue  
        self.distance_metric = method
        self.custom_metric_name = metric_name
        
        # Validate distance metric
        if method not in ["min", "max", "mean", "closest"]:
            raise ValueError(f"Invalid method: {method}. Options: min, max, mean, closest")
        
        # Initialize base class
        super().__init__(**kwargs)
        
        # Set up dependency
        if hasattr(input, 'config'):
            self.dependencies.append(input.config)
    
    def get_metric_name(self) -> str:
        """Get the default metric name.""" 
        if self.custom_metric_name:
            return self.custom_metric_name
        return "distance"
    
    def validate_params(self):
        """Validate ResidueAtomDistance parameters."""
        if not isinstance(self.distance_input, (ToolOutput, StandardizedOutput)):
            raise ValueError("Input must be a ToolOutput or StandardizedOutput object")
        
        if not self.atom_selection:
            raise ValueError("atom selection cannot be empty")
        
        if not self.residue_selection:
            raise ValueError("residue selection cannot be empty")
    
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
        
        config_lines.extend([
            f"ATOM SELECTION: {self.atom_selection}",
            f"RESIDUE SELECTION: {self.residue_selection}",
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
        
        # Output CSV path
        analysis_csv = os.path.join(output_folder, "distance_analysis.csv")
        
        # Create config file for distance calculation
        config_file = os.path.join(runtime_folder, "distance_config.json")
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
        analysis_csv = os.path.join(self.output_folder, "distance_analysis.csv")
        
        datasheets = {
            "analysis": {
                "path": analysis_csv,
                "columns": ["id", "source_structure", self.get_metric_name()],
                "description": f"Distance analysis: {self.atom_selection} to {self.residue_selection}",
                "count": len(self.input_structures) if hasattr(self, 'input_structures') else "variable"
            }
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
                "atom_selection": self.atom_selection,
                "residue_selection": self.residue_selection,
                "distance_metric": self.distance_metric,
                "metric_name": self.get_metric_name()
            }
        })
        return base_dict