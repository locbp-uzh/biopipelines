"""
Analysis base classes for protein design pipeline analysis.

Defines the abstract interface for analysis tools that generate CSV tables
with metrics for structures, sequences, compounds, etc. Analysis tools preserve
all information and don't filter - filtering is done by separate Filter tools.
"""

import re
from abc import ABC, abstractmethod
from typing import Dict, List, Any, Optional, Union, Tuple


class AnalysisResult:
    """
    Result from applying a single analysis tool.
    
    Contains information about generated analysis data and CSV output files.
    """
    
    def __init__(self, 
                 output_csv: str,
                 item_count: int,
                 analysis_info: Dict[str, Any],
                 metrics_generated: List[str] = None):
        """
        Initialize analysis result.
        
        Args:
            output_csv: Path to generated CSV file with analysis results
            item_count: Number of items analyzed
            analysis_info: Information about the analysis performed
            metrics_generated: List of metric column names in the CSV
        """
        self.output_csv = output_csv
        self.item_count = item_count
        self.analysis_info = analysis_info
        self.metrics_generated = metrics_generated or []
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "output_csv": self.output_csv,
            "item_count": self.item_count,
            "analysis_info": self.analysis_info,
            "metrics_generated": self.metrics_generated
        }
    
    def summary(self) -> str:
        """Get human-readable summary."""
        analysis_name = self.analysis_info.get('analysis_type', 'Analysis')
        metrics_str = ", ".join(self.metrics_generated) if self.metrics_generated else "unknown"
        return (f"{analysis_name}: analyzed {self.item_count} items, "
                f"generated metrics: {metrics_str}")


class Analysis(ABC):
    """
    Abstract base class for all analysis tools.
    
    An Analysis tool generates CSV tables with metrics for items 
    (structures, sequences, compounds, etc.) without filtering them.
    All items are preserved and metrics are calculated for downstream filtering.
    """
    
    def __init__(self, metric_name: str = None, **kwargs):
        """
        Initialize base analysis tool.
        
        Args:
            metric_name: Optional name for the main metric column (for collision avoidance)
            **kwargs: Analysis-specific parameters
        """
        self.metric_name = metric_name
        self.parameters = kwargs
    
    @abstractmethod
    def get_metric_name(self) -> str:
        """
        Get the default metric column name for this analysis.
        
        Returns:
            Metric name (e.g., 'distance', 'pLDDT', 'length')
        """
        pass
    
    @abstractmethod
    def get_analysis_type(self) -> str:
        """
        Get the type identifier for this analysis.
        
        Returns:
            Analysis type string (e.g., 'residue_atom_distance', 'confidence')
        """
        pass
    
    @abstractmethod
    def get_runtime_script_path(self) -> str:
        """
        Get the path to the runtime helper script for this analysis.
        
        Returns:
            Path to the pipe_criterion_*.py script that performs this analysis
        """
        pass
    
    def get_effective_metric_name(self) -> str:
        """
        Get the effective metric name (custom name if provided, otherwise default).
        
        Returns:
            The metric column name to use in output CSV
        """
        return self.metric_name if self.metric_name else self.get_metric_name()
    
    def generate_runtime_command(self, input_dir: str, output_csv: str, config_file: str) -> str:
        """
        Generate bash command to run this analysis at runtime.
        
        Args:
            input_dir: Directory containing files to analyze
            output_csv: Path where analysis CSV should be written
            config_file: Path to JSON config file with analysis parameters
            
        Returns:
            Bash command string to execute this analysis
        """
        script_path = self.get_runtime_script_path()
        return f"python {script_path} --input_dir {input_dir} --output_csv {output_csv} --config {config_file}"
    
    def to_config_dict(self) -> Dict[str, Any]:
        """
        Convert analysis to configuration dictionary for runtime scripts.
        
        Returns:
            Dictionary with all analysis parameters for JSON serialization
        """
        config = {
            "analysis_type": self.get_analysis_type(),
            "analysis_class": self.__class__.__name__,
            "metric_name": self.get_effective_metric_name(),
            "default_metric_name": self.get_metric_name(),
            "parameters": self.parameters.copy()
        }
        return config
    
    def _get_item_display_name(self, item: str) -> str:
        """Get a short display name for an item."""
        import os
        return os.path.basename(item) if '/' in item or '\\' in item else item
    
    def get_parameter_summary(self) -> List[str]:
        """
        Get human-readable summary of analysis parameters.
        
        Returns:
            List of parameter description strings
        """
        summary = [f"Metric: {self.get_effective_metric_name()}"]
        
        # Add analysis-specific parameters (override in subclasses)
        for key, value in self.parameters.items():
            if key not in ['metric_name']:  # Don't duplicate metric name
                summary.append(f"{key}: {value}")
        
        return summary
    
    def __str__(self) -> str:
        """String representation of the analysis."""
        return f"{self.__class__.__name__}(metric={self.get_effective_metric_name()})"
    
    def __repr__(self) -> str:
        return self.__str__()