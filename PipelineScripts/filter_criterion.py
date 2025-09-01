"""
FilterCriterion base classes for protein design pipeline filtering.

Defines the abstract interface for filtering criteria that can be combined
within Filter tools to create sophisticated filtering workflows.
"""

import re
from abc import ABC, abstractmethod
from typing import Dict, List, Any, Optional, Union, Tuple


class FilterResult:
    """
    Result from applying a single filtering criterion.
    
    Contains information about which items passed/failed the criterion
    and their scores.
    """
    
    def __init__(self, 
                 kept_items: List[str], 
                 filtered_items: List[str],
                 criterion_info: Dict[str, Any],
                 item_scores: Dict[str, float] = None):
        """
        Initialize filter result.
        
        Args:
            kept_items: List of item IDs that passed the criterion
            filtered_items: List of item IDs that were filtered out
            criterion_info: Information about the criterion used
            item_scores: Optional scores for each item
        """
        self.kept_items = kept_items
        self.filtered_items = filtered_items
        self.criterion_info = criterion_info
        self.item_scores = item_scores or {}
        
        self.total_input = len(kept_items) + len(filtered_items)
        self.kept_count = len(kept_items)
        self.filtered_count = len(filtered_items)
        self.pass_rate = self.kept_count / self.total_input if self.total_input > 0 else 0.0
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "kept_items": self.kept_items,
            "filtered_items": self.filtered_items,
            "criterion_info": self.criterion_info,
            "item_scores": self.item_scores,
            "total_input": self.total_input,
            "kept_count": self.kept_count,
            "filtered_count": self.filtered_count,
            "pass_rate": self.pass_rate
        }
    
    def summary(self) -> str:
        """Get human-readable summary."""
        criterion_name = self.criterion_info.get('criterion_type', 'Filter')
        return (f"{criterion_name}: {self.kept_count}/{self.total_input} items kept "
                f"({self.pass_rate:.1%} pass rate)")


class FilterCriterion(ABC):
    """
    Abstract base class for all filtering criteria.
    
    A FilterCriterion defines a single filtering rule that can be applied
    to a list of items (structures, sequences, compounds, etc.) and determines
    which items pass or fail the criterion.
    """
    
    def __init__(self, expression: str, **kwargs):
        """
        Initialize base filtering criterion.
        
        Args:
            expression: Boolean expression defining the filtering rule
                       (e.g., 'distance<=5', 'pLDDT>90')
            **kwargs: Additional criterion-specific parameters
        """
        self.expression = expression
        self.parameters = kwargs
        
        # Validate expression format
        self._validate_expression()
    
    def _validate_expression(self):
        """Validate that the expression format is safe and well-formed."""
        # Check for basic safety (only allow mathematical expressions)
        allowed_pattern = r'^[a-zA-Z_][a-zA-Z0-9_]*\s*[<>=!]+\s*[0-9\.\s\+\-\*\/\(\)<=>\s\&\|and\sor\s]+$'
        
        if not re.match(allowed_pattern, self.expression.lower()):
            raise ValueError(f"Invalid expression format: {self.expression}")
        
        # Must contain the expected variable name
        expected_var = self.get_variable_name()
        if expected_var.lower() not in self.expression.lower():
            raise ValueError(f"Expression must contain '{expected_var}' variable, got: {self.expression}")
    
    @abstractmethod
    def get_variable_name(self) -> str:
        """
        Get the variable name used in expressions for this criterion.
        
        Returns:
            Variable name (e.g., 'distance', 'pLDDT', 'length')
        """
        pass
    
    @abstractmethod
    def get_criterion_type(self) -> str:
        """
        Get the type identifier for this criterion.
        
        Returns:
            Criterion type string (e.g., 'residue_atom_distance', 'confidence')
        """
        pass
    
    @abstractmethod
    def get_runtime_script_path(self) -> str:
        """
        Get the path to the runtime helper script for this criterion.
        
        Returns:
            Path to the pipe_criterion_*.py script that evaluates this criterion
        """
        pass
    
    def evaluate_expression(self, score: float) -> bool:
        """
        Evaluate the filtering expression with the calculated score.
        
        Args:
            score: Calculated score for the item
            
        Returns:
            True if the item passes the criterion
        """
        # Replace variable name with actual score
        variable_name = self.get_variable_name()
        expr = self.expression.lower().replace(variable_name.lower(), str(score))
        
        # Replace logical operators for Python evaluation
        expr = expr.replace(' and ', ' & ').replace(' or ', ' | ')
        
        try:
            # Use eval safely with restricted globals
            safe_globals = {
                "__builtins__": {},
                "__name__": "__main__",
                "__doc__": None,
            }
            return bool(eval(expr, safe_globals))
        except Exception as e:
            raise ValueError(f"Error evaluating expression '{self.expression}' with {variable_name}={score}: {e}")
    
    def generate_runtime_command(self, structures_dir: str, output_file: str, config_file: str) -> str:
        """
        Generate bash command to evaluate this criterion at runtime.
        
        Args:
            structures_dir: Directory containing structure files to analyze
            output_file: Path where criterion results should be written
            config_file: Path to JSON config file with criterion parameters
            
        Returns:
            Bash command string to execute this criterion
        """
        script_path = self.get_runtime_script_path()
        return f"python {script_path} --structures_dir {structures_dir} --output {output_file} --config {config_file}"
    
    def to_config_dict(self) -> Dict[str, Any]:
        """
        Convert criterion to configuration dictionary for runtime scripts.
        
        Returns:
            Dictionary with all criterion parameters for JSON serialization
        """
        config = {
            "criterion_type": self.get_criterion_type(),
            "criterion_class": self.__class__.__name__,
            "expression": self.expression,
            "variable_name": self.get_variable_name(),
            "parameters": self.parameters.copy()
        }
        return config
    
    def _get_item_display_name(self, item: str) -> str:
        """Get a short display name for an item."""
        import os
        return os.path.basename(item) if '/' in item or '\\' in item else item
    
    def get_parameter_summary(self) -> List[str]:
        """
        Get human-readable summary of criterion parameters.
        
        Returns:
            List of parameter description strings
        """
        summary = [f"Expression: {self.expression}"]
        
        # Add criterion-specific parameters (override in subclasses)
        for key, value in self.parameters.items():
            if key != 'expression':  # Don't duplicate expression
                summary.append(f"{key}: {value}")
        
        return summary
    
    def __str__(self) -> str:
        """String representation of the criterion."""
        return f"{self.__class__.__name__}({self.expression})"
    
    def __repr__(self) -> str:
        return self.__str__()