"""
StructureFilter specialized tool class for structure-based filtering.

Provides a specialized interface for filtering protein structures with
structure-specific criteria, while maintaining all the flexibility of
the base Filter class.
"""

from typing import Dict, List, Any, Optional, Union

try:
    from .filter import Filter
    from .filter_criterion import FilterCriterion
    from .structure_criterion import StructureCriterion
    from .base_config import ToolOutput, StandardizedOutput
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from filter import Filter
    from filter_criterion import FilterCriterion
    from structure_criterion import StructureCriterion
    from base_config import ToolOutput, StandardizedOutput


class StructureFilter(Filter):
    """
    Specialized filter tool for protein structure filtering.
    
    Provides structure-specific functionality while maintaining all the
    capabilities of the base Filter class. Automatically sets filter_type
    to "structures" and provides enhanced validation for structure criteria.
    """
    
    # Tool identification
    TOOL_NAME = "StructureFilter"
    
    def __init__(self,
                 criteria: List[FilterCriterion],
                 input: Union[ToolOutput, StandardizedOutput, Dict[str, Any]],
                 combination: str = "AND",
                 score_weights: Optional[Dict[str, float]] = None,
                 max_items: Optional[int] = None,
                 **kwargs):
        """
        Initialize structure filter tool.
        
        Args:
            criteria: List of FilterCriterion objects (preferably StructureCriterion)
            input: Input from previous pipeline tool containing structures
            combination: How to combine criteria ("AND", "OR", "WEIGHTED")
            score_weights: Weights for WEIGHTED combination
            max_items: Maximum number of structures to keep after filtering
            **kwargs: Additional parameters for Filter base class
        """
        # Validate that we have structure-compatible criteria
        self._validate_structure_criteria(criteria)
        
        # Initialize base class with structure-specific settings
        super().__init__(
            criteria=criteria,
            input=input,
            combination=combination,
            score_weights=score_weights,
            max_items=max_items,
            filter_type="structures",  # Always filter structures
            **kwargs
        )
    
    def _validate_structure_criteria(self, criteria: List[FilterCriterion]):
        """
        Validate that criteria are compatible with structure filtering.
        
        Args:
            criteria: List of criteria to validate
        """
        if not criteria:
            raise ValueError("At least one FilterCriterion must be provided")
        
        # Check that criteria are preferably StructureCriterion
        non_structure_criteria = []
        for criterion in criteria:
            if not isinstance(criterion, StructureCriterion):
                non_structure_criteria.append(criterion.__class__.__name__)
        
        if non_structure_criteria:
            print(f"Warning: Non-structure criteria detected in StructureFilter: {non_structure_criteria}")
            print("Consider using the base Filter class or structure-specific criteria")
    
    def validate_params(self):
        """Validate structure filter-specific parameters."""
        super().validate_params()
        
        # Additional structure-specific validation
        if self.filter_type != "structures":
            raise ValueError(f"StructureFilter must have filter_type='structures', got '{self.filter_type}'")
        
        # Check that we have structure files in input
        structure_extensions = ['.pdb', '.cif', '.mmcif']
        valid_structures = []
        
        for item in self.input_items:
            if isinstance(item, str):
                if any(item.lower().endswith(ext) for ext in structure_extensions):
                    valid_structures.append(item)
        
        if not valid_structures and self.input_items:
            print(f"Warning: Input items may not be structure files. "
                  f"Expected extensions: {structure_extensions}")
    
    def get_structure_context(self) -> Dict[str, Any]:
        """
        Get enhanced context information for structure analysis.
        
        Returns:
            Enhanced context dictionary with structure-specific information
        """
        base_context = {
            'datasheets': self.input_datasheets,
            'filter_type': self.filter_type,
            'standardized_input': self.standardized_input
        }
        
        # Add structure-specific context
        base_context.update({
            'structure_files': self.input_items,
            'structure_count': len(self.input_items),
            # Could add more structure-specific context like:
            # - PDB metadata if available
            # - Chain information
            # - Ligand information
        })
        
        return base_context
    
    def apply_filtering(self, items: List[str]) -> 'CombinedFilterResult':
        """
        Apply structure filtering with enhanced context.
        
        Args:
            items: List of structure files to filter
            
        Returns:
            CombinedFilterResult with filtering outcome
        """
        # Validate structure files before processing
        valid_items = []
        invalid_items = []
        
        for item in items:
            if isinstance(item, str) and ('.pdb' in item.lower() or '.cif' in item.lower() or '.mmcif' in item.lower()):
                valid_items.append(item)
            else:
                invalid_items.append(item)
        
        if invalid_items:
            print(f"Warning: Skipping {len(invalid_items)} non-structure items:")
            for item in invalid_items[:3]:  # Show first 3
                print(f"  - {item}")
            if len(invalid_items) > 3:
                print(f"  ... and {len(invalid_items) - 3} more")
        
        if not valid_items:
            print("No valid structure files found for filtering")
            # Return empty result
            from filter import CombinedFilterResult
            return CombinedFilterResult([], items, {"error": "no_valid_structures"}, [])
        
        # Use the enhanced structure context
        print(f"Filtering {len(valid_items)} structure files")
        
        # Call parent method with valid structure files
        return super().apply_filtering(valid_items)
    
    def get_config_display(self) -> List[str]:
        """Get structure-specific configuration display."""
        config_lines = super().get_config_display()
        
        # Add structure-specific information
        structure_info = []
        
        # Count different structure types
        pdb_count = sum(1 for item in self.input_items if item.lower().endswith('.pdb'))
        cif_count = sum(1 for item in self.input_items if item.lower().endswith(('.cif', '.mmcif')))
        other_count = len(self.input_items) - pdb_count - cif_count
        
        if pdb_count > 0:
            structure_info.append(f"PDB files: {pdb_count}")
        if cif_count > 0:
            structure_info.append(f"CIF files: {cif_count}")
        if other_count > 0:
            structure_info.append(f"Other files: {other_count}")
        
        if structure_info:
            config_lines.extend(structure_info)
        
        # Count structure-specific criteria
        structure_criteria_count = sum(1 for criterion in self.criteria 
                                     if isinstance(criterion, StructureCriterion))
        if structure_criteria_count > 0:
            config_lines.append(f"Structure criteria: {structure_criteria_count}")
        
        return config_lines
    
    @staticmethod
    def create_distance_filter(atom: str, residue: str, expression: str, **kwargs):
        """
        Convenience method to create a structure filter with distance criterion.
        
        Args:
            atom: Atom selection string
            residue: Residue selection string  
            expression: Distance expression (e.g., 'distance<=5')
            **kwargs: Additional arguments for StructureFilter
            
        Returns:
            Configured StructureFilter instance
        """
        from residue_atom_distance import ResidueAtomDistance
        
        distance_criterion = ResidueAtomDistance(
            atom=atom,
            residue=residue, 
            expression=expression
        )
        
        return StructureFilter(
            criteria=[distance_criterion],
            **kwargs
        )
    
    @staticmethod
    def create_confidence_filter(expression: str, selection: Optional[str] = None, **kwargs):
        """
        Convenience method to create a structure filter with confidence criterion.
        
        Args:
            expression: Confidence expression (e.g., 'pLDDT>90')
            selection: Optional region selection
            **kwargs: Additional arguments for StructureFilter
            
        Returns:
            Configured StructureFilter instance
        """
        from confidence import Confidence
        
        confidence_criterion = Confidence(
            expression=expression,
            selection=selection
        )
        
        return StructureFilter(
            criteria=[confidence_criterion],
            **kwargs
        )