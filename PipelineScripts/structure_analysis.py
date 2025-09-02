"""
StructureAnalysis base class for structure-based analysis tools.

Provides common functionality for analyzing protein structures, including
library imports, file parsing, and standardized interfaces for structure analysis.
"""

import os
import pandas as pd
from typing import Dict, List, Any, Optional, Union, Tuple
from abc import ABC, abstractmethod

try:
    from .analysis import Analysis
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from analysis import Analysis


class StructureAnalysis(Analysis):
    """
    Abstract base class for structure-based analysis tools.
    
    Provides common functionality for analyzing protein structures including:
    - Structure file parsing and validation
    - Selection string parsing
    - Coordinate extraction
    - Library availability checking
    """
    
    def __init__(self, metric_name: str = None, **kwargs):
        """
        Initialize structure analysis.
        
        Args:
            metric_name: Optional custom name for the metric column
            **kwargs: Structure-specific analysis parameters
        """
        super().__init__(metric_name, **kwargs)
        
        # Note: Library availability is checked at runtime by helper scripts
        # During pipeline compilation, we just store parameters
    
    def get_runtime_script_path(self) -> str:
        """
        Get path to the runtime helper script for structure analysis.
        
        Returns:
            Path relative to HelpScripts folder
        """
        analysis_name = self.get_analysis_type()
        return f"pipe_criterion_{analysis_name}.py"
    
    def parse_selection_string(self, selection: str) -> Dict[str, Any]:
        """
        Parse selection string into structured criteria.
        
        Handles various selection formats:
        - Entity.atom: 'ligand.Cl', 'protein.CA'
        - Sequence context: 'protein.D in TRGDTGH'
        - MDAnalysis style: 'resname LIG and name CL'
        - Residue IDs: 'resid 145-150'
        
        Args:
            selection: Selection string
            
        Returns:
            Dictionary with parsed selection criteria
        """
        if '.' in selection and ' in ' in selection:
            # Format: 'protein.D in TRGDTGH'
            parts = selection.split(' in ')
            if len(parts) == 2:
                entity_res, sequence = parts
                if '.' in entity_res:
                    entity, residue_type = entity_res.split('.', 1)
                    return {
                        "entity": entity.lower(),
                        "residue_type": residue_type,
                        "sequence_context": sequence,
                        "selection_type": "sequence_context"
                    }
        elif '.' in selection:
            # Format: 'ligand.Cl' or 'protein.CA'
            parts = selection.split('.', 1)
            if len(parts) == 2:
                entity, atom_type = parts
                return {
                    "entity": entity.lower(),
                    "atom_type": atom_type,
                    "selection_type": "entity_atom"
                }
        elif selection.lower().startswith('resname') or selection.lower().startswith('resid'):
            # MDAnalysis-style selection
            return {
                "raw_selection": selection,
                "selection_type": "mdanalysis"
            }
        
        # Default: treat as raw selection string
        return {
            "raw_selection": selection,
            "selection_type": "raw"
        }
    
    def resolve_datasheet_selection(self, selection_ref: str, context: Dict[str, Any]) -> Optional[List[int]]:
        """
        Resolve datasheet-based selection to list of residue numbers.
        
        Args:
            selection_ref: Reference like "input.datasheets.sequences.designed_residues"
            context: Context containing datasheets
            
        Returns:
            List of residue numbers or None if resolution fails
        """
        try:
            # Parse the reference
            if not selection_ref.startswith("input.datasheets."):
                return None
            
            parts = selection_ref.split(".")
            if len(parts) < 4:
                return None
            
            datasheet_name = parts[2]  # e.g., "sequences"
            column_name = parts[3]     # e.g., "designed_residues"
            
            # Look for datasheet in context
            datasheets = context.get('datasheets', {})
            
            if isinstance(datasheets, dict) and datasheet_name in datasheets:
                datasheet_info = datasheets[datasheet_name]
                
                # Get path from datasheet info
                if isinstance(datasheet_info, dict) and 'path' in datasheet_info:
                    datasheet_path = datasheet_info['path']
                else:
                    datasheet_path = str(datasheet_info)
                
                if os.path.exists(datasheet_path):
                    df = pd.read_csv(datasheet_path)
                    
                    if column_name in df.columns:
                        # Parse position specifications
                        positions = []
                        for value in df[column_name].dropna():
                            if isinstance(value, str):
                                positions.extend(self._parse_position_string(value))
                            elif isinstance(value, int):
                                positions.append(value)
                        return positions if positions else None
        
        except Exception as e:
            print(f"Warning: Could not resolve datasheet selection '{selection_ref}': {e}")
        
        return None
    
    def _parse_position_string(self, position_str: str) -> List[int]:
        """
        Parse position string like "45-50,52,55-57" into list of residue numbers.
        
        Args:
            position_str: Position specification string
            
        Returns:
            List of residue numbers
        """
        positions = []
        
        for part in str(position_str).split(','):
            part = part.strip()
            if '-' in part and not part.startswith('-'):
                # Range specification
                try:
                    start, end = map(int, part.split('-'))
                    positions.extend(range(start, end + 1))
                except ValueError:
                    continue
            else:
                # Single position
                try:
                    positions.append(int(part))
                except ValueError:
                    continue
        
        return positions
    
    # Structure analysis methods moved to runtime helper scripts
    # These utility methods are kept for reference and testing purposes
    
    def _get_structure_analysis_utilities(self) -> Dict[str, Any]:
        """
        Get structure analysis utility functions for runtime scripts.
        
        Returns:
            Dictionary containing utility functions that can be serialized
        """
        return {
            "parse_selection_string": self.parse_selection_string,
            "resolve_datasheet_selection": self.resolve_datasheet_selection,
            "_parse_position_string": self._parse_position_string
        }