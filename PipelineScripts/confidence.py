"""
Confidence criterion for filtering structures based on confidence scores.

Filters protein structures based on predicted confidence scores like pLDDT,
commonly used for assessing structural quality and confidence in protein predictions.
Uses 'pLDDT' as the variable name in expressions.
"""

import os
import numpy as np
import pandas as pd
from typing import Dict, List, Any, Optional, Union, Tuple

try:
    from .structure_criterion import StructureCriterion
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from structure_criterion import StructureCriterion


class Confidence(StructureCriterion):
    """
    Filter structures based on confidence scores (pLDDT).
    
    Uses 'pLDDT' as the variable name in expressions.
    
    pLDDT (predicted Local Distance Difference Test) is a confidence measure
    for protein structure predictions, commonly used in AlphaFold and similar methods.
    Scores range from 0-100, where higher values indicate higher confidence.
    
    Common thresholds:
    - Very high confidence: >90
    - High confidence: 70-90  
    - Low confidence: 50-70
    - Very low confidence: <50
    """
    
    def __init__(self,
                 expression: str,
                 selection: Optional[Union[str, List[int], List[str]]] = None,
                 score_metric: str = "mean",
                 score_source: str = "bfactor",
                 **kwargs):
        """
        Initialize confidence criterion.
        
        Args:
            expression: pLDDT constraint expression using 'pLDDT' variable
                       (e.g., 'pLDDT>90', 'pLDDT>=70 and pLDDT<=90')
            selection: Region selection for scoring:
                - None: Use all residues
                - str: Datasheet reference like "input.datasheets.sequences.designed_residues"
                - List[int]: Residue numbers like [45, 46, 47, 48]
                - List[str]: Residue identifiers like ["A45", "A46", "B12"]
            score_metric: How to aggregate pLDDT scores ("mean", "min", "max", "median")
            score_source: Source of scores in PDB ("bfactor", "occupancy")
            **kwargs: Additional parameters
            
        Examples:
            # High confidence structures (mean pLDDT > 90)
            Confidence(
                expression='pLDDT>90',
                score_metric='mean'
            )
            
            # Filter based on designed regions only
            Confidence(
                expression='pLDDT>=70',
                selection='input.datasheets.sequences.designed_residues',
                score_metric='min'  # Ensure all designed residues are confident
            )
            
            # Specific residue range
            Confidence(
                expression='pLDDT>80',
                selection=[45, 46, 47, 48, 49, 50],
                score_metric='mean'
            )
        """
        self.selection = selection
        self.score_metric = score_metric
        self.score_source = score_source
        
        # Validate parameters
        if score_metric not in ["mean", "min", "max", "median"]:
            raise ValueError(f"Invalid score_metric: {score_metric}")
        
        if score_source not in ["bfactor", "occupancy"]:
            raise ValueError(f"Invalid score_source: {score_source}")
        
        # Initialize parent with expression and store parameters
        super().__init__(
            expression=expression,
            selection=selection,
            score_metric=score_metric,
            score_source=score_source,
            **kwargs
        )
    
    def get_variable_name(self) -> str:
        """Get the variable name used in expressions."""
        return "pLDDT"
    
    def get_criterion_type(self) -> str:
        """Get the criterion type identifier."""
        return "confidence"
    
    def get_runtime_script_path(self) -> str:
        """
        Get path to runtime script for confidence calculation.
        
        Returns:
            Path to the pipe_criterion_confidence.py script
        """
        return "pipe_criterion_confidence.py"
    
    # Confidence calculation methods moved to runtime helper script
    # These utility methods are included for reference in the runtime config
    
    def get_parameter_summary(self) -> List[str]:
        """Get parameter summary for display."""
        summary = super().get_parameter_summary()
        
        # Add confidence-specific parameters
        summary.extend([
            f"Selection: {self.selection if self.selection else 'all_residues'}",
            f"Score metric: {self.score_metric}",
            f"Score source: {self.score_source}"
        ])
        
        return summary
    
    def __str__(self) -> str:
        """String representation."""
        selection_str = str(self.selection) if self.selection else "all"
        return f"Confidence({selection_str}: {self.expression})"