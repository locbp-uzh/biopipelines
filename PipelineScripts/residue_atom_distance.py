"""
ResidueAtomDistance criterion for filtering structures based on atom-residue distances.

Filters protein structures based on distances between specific atoms and residues,
commonly used for ligand-protein binding site analysis and interaction validation.
Uses 'distance' as the variable name in expressions.
"""

import numpy as np
from typing import Dict, List, Any, Optional, Union, Tuple

try:
    from .structure_criterion import StructureCriterion
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from structure_criterion import StructureCriterion


class ResidueAtomDistance(StructureCriterion):
    """
    Filter structures based on distances between specific atoms and residues.
    
    Uses 'distance' as the variable name in expressions.
    
    Commonly used for:
    - Ligand-protein binding site validation
    - Metal coordination analysis  
    - Catalytic site geometry verification
    - Interaction distance constraints
    """
    
    def __init__(self,
                 atom: str,
                 residue: str,
                 expression: str,
                 distance_metric: str = "min",
                 **kwargs):
        """
        Initialize atom-residue distance criterion.
        
        Args:
            atom: Atom selection string (e.g., 'ligand.Cl', 'protein.CA', 'resname ZN')
            residue: Residue selection string (e.g., 'protein.D in TRGDTGH', 'resid 145-150')
            expression: Distance constraint expression using 'distance' variable
                       (e.g., 'distance<=5.0', 'distance>2.5 and distance<8.0')
            distance_metric: How to calculate distance ("min", "max", "mean", "closest")
            **kwargs: Additional parameters
            
        Examples:
            # Ligand chlorine within 5Å of specific aspartic acids
            ResidueAtomDistance(
                atom='ligand.Cl',
                residue='protein.D in TRGDTGH', 
                expression='distance<=5.0'
            )
            
            # Protein CA atoms not too close to ligand
            ResidueAtomDistance(
                atom='protein.CA',
                residue='ligand.all',
                expression='distance>3.0'
            )
            
            # Metal coordination within specific range
            ResidueAtomDistance(
                atom='resname ZN',
                residue='protein.HIS.NE2',
                expression='distance>=1.9 and distance<=2.3'
            )
        """
        self.atom_selection = atom
        self.residue_selection = residue
        self.distance_metric = distance_metric
        
        # Validate distance metric
        if distance_metric not in ["min", "max", "mean", "closest"]:
            raise ValueError(f"Invalid distance_metric: {distance_metric}")
        
        # Initialize parent with expression and store parameters
        super().__init__(
            expression=expression,
            atom=atom,
            residue=residue,
            distance_metric=distance_metric,
            **kwargs
        )
    
    def get_variable_name(self) -> str:
        """Get the variable name used in expressions."""
        return "distance"
    
    def get_criterion_type(self) -> str:
        """Get the criterion type identifier."""
        return "residue_atom_distance"
    
    def get_runtime_script_path(self) -> str:
        """
        Get path to runtime script for distance calculation.
        
        Returns:
            Path to the pipe_criterion_residue_atom_distance.py script
        """
        return "pipe_criterion_residue_atom_distance.py"
    
    def get_parameter_summary(self) -> List[str]:
        """Get parameter summary for display."""
        summary = super().get_parameter_summary()
        
        # Add distance-specific parameters
        summary.extend([
            f"Atom selection: {self.atom_selection}",
            f"Residue selection: {self.residue_selection}",
            f"Distance metric: {self.distance_metric}"
        ])
        
        return summary
    
    def __str__(self) -> str:
        """String representation."""
        return f"ResidueAtomDistance({self.atom_selection} -> {self.residue_selection}: {self.expression})"