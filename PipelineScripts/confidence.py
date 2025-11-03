"""
Confidence analysis for extracting confidence scores from protein structures.

Analyzes protein structures to extract predicted confidence scores like pLDDT,
commonly used for assessing structural quality and confidence in protein predictions.
Outputs CSV with confidence metrics for all structures.
"""

import os
import numpy as np
import pandas as pd
from typing import Dict, List, Any, Optional, Union, Tuple

try:
    from .analysis import Analysis
    from .base_config import BaseConfig
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from analysis import Analysis
    from base_config import BaseConfig


class Confidence(Analysis):
    """
    Analyze structures to extract confidence scores (pLDDT).
    
    Generates CSV with confidence metrics for all input structures.
    
    pLDDT (predicted Local Distance Difference Test) is a confidence measure
    for protein structure predictions, commonly used in AlphaFold and similar methods.
    Scores range from 0-100, where higher values indicate higher confidence.
    
    Common interpretation:
    - Very high confidence: >90
    - High confidence: 70-90  
    - Low confidence: 50-70
    - Very low confidence: <50
    """
    
    def __init__(self,
                 selection: Optional[Union[str, List[int], List[str]]] = None,
                 score_metric: str = "mean",
                 score_source: str = "bfactor",
                 metric_name: str = None,
                 **kwargs):
        """
        Initialize confidence analysis.
        
        Args:
            selection: Region selection for scoring:
                - None: Use all residues
                - str: Table reference like "sequences.designed_residues"
                - List[int]: Residue numbers like [45, 46, 47, 48]
                - List[str]: Residue identifiers like ["A45", "A46", "B12"]
            score_metric: How to aggregate pLDDT scores ("mean", "min", "max", "median")
            score_source: Source of scores in PDB ("bfactor", "occupancy")
            metric_name: Custom name for the pLDDT column (default: "pLDDT")
            **kwargs: Additional parameters
            
        Examples:
            # Analyze all residues with mean pLDDT
            Confidence(score_metric='mean')
            
            # Analyze designed regions only
            Confidence(
                selection='sequences.designed_residues',
                score_metric='min'
            )
            
            # Analyze specific residue range with custom metric name
            Confidence(
                selection=[45, 46, 47, 48, 49, 50],
                score_metric='mean',
                metric_name='binding_site_confidence'
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
        
        # Initialize parent with parameters
        super().__init__(
            metric_name=metric_name,
            selection=selection,
            score_metric=score_metric,
            score_source=score_source,
            **kwargs
        )
    
    def get_metric_name(self) -> str:
        """Get the default metric name."""
        return "pLDDT"
    
    def get_analysis_type(self) -> str:
        """Get the analysis type identifier."""
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
        metric = self.get_effective_metric_name()
        return f"Confidence({selection_str}: {metric}, {self.score_metric})"