"""
Filter base tool class for protein design pipeline filtering.

Manages multiple FilterCriterion objects and combines their results with
sophisticated logic (AND, OR, WEIGHTED) while maintaining pipeline compatibility.
"""

import os
import json
from typing import Dict, List, Any, Optional, Union
from collections import defaultdict

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput
    from .filter_criterion import FilterCriterion, FilterResult
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput
    from filter_criterion import FilterCriterion, FilterResult


class CombinedFilterResult:
    """
    Result from combining multiple FilterCriterion results.
    
    Contains information about the final kept/filtered items after applying
    combination logic and any max_items limitation.
    """
    
    def __init__(self, 
                 kept_items: List[str],
                 filtered_items: List[str],
                 combination_info: Dict[str, Any],
                 individual_results: List[FilterResult],
                 final_scores: Dict[str, float] = None):
        """
        Initialize combined filter result.
        
        Args:
            kept_items: Final list of items that passed all criteria
            filtered_items: Final list of items that were filtered out
            combination_info: Information about how results were combined
            individual_results: Results from each individual criterion
            final_scores: Final combined scores for all items
        """
        self.kept_items = kept_items
        self.filtered_items = filtered_items
        self.combination_info = combination_info
        self.individual_results = individual_results
        self.final_scores = final_scores or {}
        
        self.total_input = len(kept_items) + len(filtered_items)
        self.kept_count = len(kept_items)
        self.filtered_count = len(filtered_items)
        self.pass_rate = self.kept_count / self.total_input if self.total_input > 0 else 0.0
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "kept_items": self.kept_items,
            "filtered_items": self.filtered_items,
            "combination_info": self.combination_info,
            "individual_results": [result.to_dict() for result in self.individual_results],
            "final_scores": self.final_scores,
            "total_input": self.total_input,
            "kept_count": self.kept_count,
            "filtered_count": self.filtered_count,
            "pass_rate": self.pass_rate
        }
    
    def summary(self) -> str:
        """Get human-readable summary."""
        method = self.combination_info.get('combination_method', 'unknown')
        return (f"Combined Filter ({method}): {self.kept_count}/{self.total_input} items kept "
                f"({self.pass_rate:.1%} pass rate)")


class Filter(BaseConfig):
    """
    Base filter tool for pipeline integration.
    
    Manages multiple FilterCriterion objects and combines their results using
    various combination methods (AND, OR, WEIGHTED). Provides pipeline integration
    through BaseConfig inheritance.
    """
    
    # Tool identification
    TOOL_NAME = "Filter"
    DEFAULT_ENV = "ProteinEnv"
    COMPATIBLE_ENVS = ["ProteinEnv", "Boltz2Env", "ligandmpnn_env"]
    DEFAULT_RESOURCES = {"gpu": "T4", "memory": "8GB", "time": "2:00:00"}
    
    def __init__(self,
                 criteria: List[FilterCriterion],
                 input: Union[ToolOutput, StandardizedOutput, Dict[str, Any]],
                 combination: str = "AND",
                 score_weights: Optional[Dict[str, float]] = None,
                 max_items: Optional[int] = None,
                 filter_type: str = "structures",
                 **kwargs):
        """
        Initialize base filter tool.
        
        Args:
            criteria: List of FilterCriterion objects to apply
            input: Input from previous pipeline tool
            combination: How to combine criteria ("AND", "OR", "WEIGHTED")
            score_weights: Weights for WEIGHTED combination (criterion class name -> weight)
            max_items: Maximum number of items to keep after filtering
            filter_type: Type of data to filter ("structures", "sequences", "compounds")
            **kwargs: Additional parameters for BaseConfig
        """
        self.criteria = criteria
        self.combination = combination
        self.score_weights = score_weights or {}
        self.max_items = max_items
        self.filter_type = filter_type
        
        # Validate parameters
        self._validate_parameters()
        
        # Initialize input handling
        if isinstance(input, StandardizedOutput):
            self.standardized_input = input
            self.input_datasheets = input.datasheets
            self.input_items = self._extract_input_items(input, filter_type)
        elif isinstance(input, ToolOutput):
            self.standardized_input = input.output
            self.input_datasheets = input.output.datasheets
            self.input_items = self._extract_input_items(input.output, filter_type)
        elif isinstance(input, dict):
            self.standardized_input = StandardizedOutput(input)
            self.input_datasheets = input.get('datasheets', {})
            self.input_items = input.get(filter_type, [])
        else:
            raise ValueError(f"Unsupported input type: {type(input)}")
        
        # Initialize file paths (set in _setup_file_paths)
        self.filter_manifest_file = None
        self.filtered_datasheet_file = None
        self.filter_report_file = None
        
        # Initialize base class
        super().__init__(**kwargs)
    
    def _validate_parameters(self):
        """Validate filter parameters."""
        if not self.criteria:
            raise ValueError("At least one FilterCriterion must be provided")
        
        if self.combination not in ["AND", "OR", "WEIGHTED"]:
            raise ValueError(f"Invalid combination method: {self.combination}")
        
        if self.combination == "WEIGHTED" and not self.score_weights:
            # Auto-assign equal weights if none provided
            self.score_weights = {criterion.__class__.__name__: 1.0 for criterion in self.criteria}
        
        if self.max_items is not None and self.max_items <= 0:
            raise ValueError("max_items must be positive if specified")
        
        if self.filter_type not in ["structures", "sequences", "compounds"]:
            raise ValueError(f"Unsupported filter_type: {self.filter_type}")
    
    def _extract_input_items(self, output: StandardizedOutput, filter_type: str) -> List[str]:
        """
        Extract the list of items to filter from standardized output.
        
        Args:
            output: Standardized output from previous tool
            filter_type: Type of data to extract
            
        Returns:
            List of item paths or IDs to filter
        """
        if filter_type == "structures":
            return output.structures or []
        elif filter_type == "sequences":
            return output.sequences or []
        elif filter_type == "compounds":
            return output.compounds or []
        else:
            raise ValueError(f"Unsupported filter type: {filter_type}")
    
    def _setup_file_paths(self):
        """Set up output file paths after output_folder is known."""
        self.filter_manifest_file = os.path.join(
            self.output_folder, f"{self.job_name}_filter_manifest.json"
        )
        self.filtered_datasheet_file = os.path.join(
            self.output_folder, f"{self.job_name}_filtered_{self.filter_type}.csv"
        )
        self.filter_report_file = os.path.join(
            self.output_folder, f"{self.job_name}_filter_report.txt"
        )
    
    def set_pipeline_context(self, pipeline_ref, execution_order: int, output_folder: str):
        """Set context when added to pipeline and initialize file paths."""
        super().set_pipeline_context(pipeline_ref, execution_order, output_folder)
        self._setup_file_paths()
    
    def validate_params(self):
        """Validate filter-specific parameters."""
        if not self.input_items:
            raise ValueError(f"No {self.filter_type} found in input for filtering")
    
    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input sources from pipeline context."""
        # Store pipeline folders for script generation
        self.folders = pipeline_folders
        
        # Input validation - items should already be set from initialization
        if not hasattr(self, 'input_items') or not self.input_items:
            raise ValueError(f"No {self.filter_type} available for filtering")
    
    # Direct filtering methods removed - filtering now happens at runtime via bash scripts
    # These methods are kept for reference and potential testing purposes
    
    def _legacy_apply_filtering(self, items: List[str]) -> CombinedFilterResult:
        """
        Legacy method for direct filtering (used only for testing).
        
        In the new architecture, filtering happens at runtime via bash scripts.
        This method is kept for reference and testing purposes only.
        
        Args:
            items: List of items to filter
            
        Returns:
            CombinedFilterResult with final filtering outcome
        """
        raise NotImplementedError(
            "Direct filtering has been replaced by runtime bash script execution. "
            "Use generate_script() to create the filter execution script."
        )
    
    def _combine_with_and(self, items: List[str], results: List[FilterResult]) -> CombinedFilterResult:
        """Combine results using AND logic (intersection)."""
        # Find items that passed ALL criteria
        kept_sets = [set(result.kept_items) for result in results]
        kept_intersection = set.intersection(*kept_sets) if kept_sets else set()
        
        # Items are filtered if they failed ANY criterion
        filtered_items = [item for item in items if item not in kept_intersection]
        
        # Combine scores (use minimum for AND logic)
        combined_scores = {}
        for item in items:
            scores = []
            for result in results:
                if item in result.item_scores:
                    scores.append(result.item_scores[item])
            
            if scores:
                combined_scores[item] = min(scores)  # Most restrictive
        
        combination_info = {
            "combination_method": "AND",
            "criteria_count": len(results),
            "criteria_types": [result.criterion_info.get('criterion_type', 'unknown') for result in results]
        }
        
        return CombinedFilterResult(
            kept_items=list(kept_intersection),
            filtered_items=filtered_items,
            combination_info=combination_info,
            individual_results=results,
            final_scores=combined_scores
        )
    
    def _combine_with_or(self, items: List[str], results: List[FilterResult]) -> CombinedFilterResult:
        """Combine results using OR logic (union)."""
        # Find items that passed ANY criterion
        kept_sets = [set(result.kept_items) for result in results]
        kept_union = set.union(*kept_sets) if kept_sets else set()
        
        # Items are filtered only if they failed ALL criteria
        filtered_items = [item for item in items if item not in kept_union]
        
        # Combine scores (use maximum for OR logic)
        combined_scores = {}
        for item in items:
            scores = []
            for result in results:
                if item in result.item_scores:
                    scores.append(result.item_scores[item])
            
            if scores:
                combined_scores[item] = max(scores)  # Most permissive
        
        combination_info = {
            "combination_method": "OR",
            "criteria_count": len(results),
            "criteria_types": [result.criterion_info.get('criterion_type', 'unknown') for result in results]
        }
        
        return CombinedFilterResult(
            kept_items=list(kept_union),
            filtered_items=filtered_items,
            combination_info=combination_info,
            individual_results=results,
            final_scores=combined_scores
        )
    
    def _combine_with_weights(self, items: List[str], results: List[FilterResult]) -> CombinedFilterResult:
        """Combine results using weighted scoring."""
        # Calculate weighted scores for each item
        combined_scores = {}
        
        for item in items:
            weighted_score = 0.0
            total_weight = 0.0
            
            for result in results:
                criterion_class = result.criterion_info.get('criterion_class', 'unknown')
                weight = self.score_weights.get(criterion_class, 1.0)
                
                if item in result.item_scores:
                    weighted_score += result.item_scores[item] * weight
                    total_weight += weight
            
            if total_weight > 0:
                combined_scores[item] = weighted_score / total_weight
            else:
                combined_scores[item] = 0.0
        
        # Keep items with positive weighted score
        kept_items = [item for item in items if combined_scores.get(item, 0.0) > 0]
        filtered_items = [item for item in items if item not in kept_items]
        
        combination_info = {
            "combination_method": "WEIGHTED",
            "score_weights": self.score_weights.copy(),
            "criteria_count": len(results),
            "criteria_types": [result.criterion_info.get('criterion_type', 'unknown') for result in results]
        }
        
        return CombinedFilterResult(
            kept_items=kept_items,
            filtered_items=filtered_items,
            combination_info=combination_info,
            individual_results=results,
            final_scores=combined_scores
        )
    
    def _apply_max_items_filter(self, result: CombinedFilterResult) -> CombinedFilterResult:
        """Apply max_items limitation to combined result."""
        if self.max_items is None or len(result.kept_items) <= self.max_items:
            return result
        
        # Sort kept items by score (descending - higher is better)
        if result.final_scores:
            scored_items = [(item, result.final_scores.get(item, 0.0)) for item in result.kept_items]
            scored_items.sort(key=lambda x: x[1], reverse=True)
            
            # Keep top max_items
            final_kept = [item for item, score in scored_items[:self.max_items]]
            additional_filtered = [item for item, score in scored_items[self.max_items:]]
        else:
            # No scores - just take first max_items
            final_kept = result.kept_items[:self.max_items]
            additional_filtered = result.kept_items[self.max_items:]
        
        # Update combination info
        updated_combination_info = result.combination_info.copy()
        updated_combination_info.update({
            "max_items_applied": self.max_items,
            "max_items_method": "score_based" if result.final_scores else "order_based"
        })
        
        return CombinedFilterResult(
            kept_items=final_kept,
            filtered_items=result.filtered_items + additional_filtered,
            combination_info=updated_combination_info,
            individual_results=result.individual_results,
            final_scores=result.final_scores
        )
    
    def generate_script(self, script_path: str) -> str:
        """
        Generate bash script for runtime filter execution.
        
        Args:
            script_path: Path where script should be written
            
        Returns:
            Script content as string
        """
        # Get helper script paths
        helpscripts_folder = self.folders.get("HelpScripts", "HelpScripts")
        filter_execution_script = os.path.join(helpscripts_folder, "pipe_filter_execution.py")
        
        # Create directory for storing structures
        structures_dir = os.path.join(self.output_folder, "structures_to_filter")
        
        # Create directory for criterion results
        criterion_results_dir = os.path.join(self.output_folder, "criterion_results")
        
        script_parts = []
        
        # Header
        script_parts.append(f"#!/bin/bash")
        script_parts.append(f"# {self.__class__.__name__} execution script")
        script_parts.append(f"# Generated by ProteinNotebooks pipeline system")
        script_parts.append("")
        script_parts.append(self.generate_completion_check_header())
        script_parts.append("")
        
        # Setup
        script_parts.append(f'echo "Starting {self.__class__.__name__} filtering..."')
        script_parts.append(f'echo "Criteria: {len(self.criteria)}"')
        script_parts.append(f'echo "Combination: {self.combination}"')
        script_parts.append(f'echo "Input {self.filter_type}: {len(self.input_items)} items"')
        script_parts.append("")
        
        # Create working directories
        script_parts.append(f'echo "Creating working directories..."')
        script_parts.append(f'mkdir -p "{structures_dir}"')
        script_parts.append(f'mkdir -p "{criterion_results_dir}"')
        script_parts.append("")
        
        # Copy/link structure files to working directory (if filter_type is structures)
        if self.filter_type == "structures":
            script_parts.append(f'echo "Preparing structure files..."')
            for item in self.input_items:
                item_name = os.path.basename(item)
                script_parts.append(f'cp "{item}" "{structures_dir}/{item_name}"')
            script_parts.append("")
        
        # Execute individual criteria
        criterion_result_files = []
        for i, criterion in enumerate(self.criteria):
            criterion_name = criterion.__class__.__name__
            script_parts.append(f'echo "Evaluating criterion {i+1}: {criterion_name}..."')
            
            # Create config file for this criterion
            config_file = os.path.join(self.output_folder, f"criterion_{i}_{criterion_name.lower()}.json")
            result_file = os.path.join(criterion_results_dir, f"criterion_{i}_{criterion_name.lower()}_results.json")
            criterion_result_files.append(result_file)
            
            # Generate config content
            config_content = json.dumps(criterion.to_config_dict(), indent=2)
            
            # Write config file and run criterion
            script_parts.append(f'cat > "{config_file}" << \'EOF\'')
            script_parts.append(config_content)
            script_parts.append("EOF")
            script_parts.append("")
            
            # Get runtime script path for this criterion
            criterion_script = os.path.join(helpscripts_folder, criterion.get_runtime_script_path())
            
            # Execute criterion
            if self.filter_type == "structures":
                cmd = f'python "{criterion_script}" --structures_dir "{structures_dir}" --output "{result_file}" --config "{config_file}"'
            else:
                # For other types, adapt as needed
                cmd = f'python "{criterion_script}" --input_dir "{structures_dir}" --output "{result_file}" --config "{config_file}"'
            
            script_parts.append(cmd)
            script_parts.append("")
        
        # Combine results
        script_parts.append(f'echo "Combining criterion results with {self.combination} method..."')
        
        # Build combine command
        combine_cmd = [
            "python", f'"{filter_execution_script}"',
            f'--input \'{json.dumps(self.standardized_input.to_dict())}\'',
            f'--filter-type "{self.filter_type}"',
            f'--output-folder "{self.output_folder}"',
            f'--job-name "{self.job_name}"',
            f'--combination "{self.combination}"',
            f'--criteria \'{json.dumps(criterion_result_files)}\''
        ]
        
        if self.score_weights:
            combine_cmd.append(f'--score-weights \'{json.dumps(self.score_weights)}\'')
        
        if self.max_items is not None:
            combine_cmd.append(f'--max-items {self.max_items}')
        
        script_parts.append(' '.join(combine_cmd) + f' | tee "{os.path.join(self.output_folder, f"{self.job_name}_filter.log")}"')
        script_parts.append("")
        
        # Cleanup
        script_parts.append(f'echo "Cleaning up temporary files..."')
        script_parts.append(f'rm -rf "{structures_dir}"')
        script_parts.append(f'rm -rf "{criterion_results_dir}"')
        script_parts.append(f'rm -f {os.path.join(self.output_folder, "criterion_*.json")}')
        script_parts.append("")
        
        script_parts.append(f'echo "Filtering completed"')
        script_parts.append("")
        script_parts.append(self.generate_completion_check_footer())
        
        return "\n".join(script_parts)
    
    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after filtering.
        
        Returns expected paths but cannot predict exact filtered content.
        The actual filtering happens at runtime.
        """
        # Ensure file paths are set up
        if not hasattr(self, 'filter_manifest_file') or self.filter_manifest_file is None:
            self._setup_file_paths()
        
        # Create output structure based on filter type
        if self.filter_type == "structures":
            return {
                "structures": self.input_items,  # Placeholder - actual filtering at runtime
                "structure_ids": [os.path.splitext(os.path.basename(f))[0] for f in self.input_items],
                "compounds": [],
                "compound_ids": [],
                "sequences": [],
                "sequence_ids": [],
                "datasheets": {
                    "filtered": {
                        "path": self.filtered_datasheet_file,
                        "columns": ["id", "file_path", "filter_passed", "final_score"],
                        "description": f"Filtering results for {self.filter_type}",
                        "count": len(self.input_items)
                    },
                    "manifest": {
                        "path": self.filter_manifest_file,
                        "columns": ["combination_method", "criteria_count", "kept_count", "filtered_count"],
                        "description": "Filter execution metadata and statistics",
                        "count": 1
                    }
                },
                "output_folder": self.output_folder,
                "filter_metadata": {
                    "filter_type": self.filter_type,
                    "input_count": len(self.input_items),
                    "max_items": self.max_items,
                    "is_filtered": True,
                    "combination_method": self.combination,
                    "criteria_count": len(self.criteria)
                }
            }
        else:
            # Generic structure for other filter types
            return {
                "structures": [],
                "structure_ids": [],
                "compounds": [],
                "compound_ids": [],
                "sequences": [],
                "sequence_ids": [],
                "datasheets": {
                    "filtered": {
                        "path": self.filtered_datasheet_file,
                        "columns": ["id", "data", "filter_passed", "final_score"],
                        "description": f"Filtering results for {self.filter_type}",
                        "count": len(self.input_items)
                    },
                    "manifest": {
                        "path": self.filter_manifest_file,
                        "columns": ["combination_method", "criteria_count", "kept_count"],
                        "description": "Filter execution metadata and statistics",
                        "count": 1
                    }
                },
                "output_folder": self.output_folder,
                "filter_metadata": {
                    "filter_type": self.filter_type,
                    "input_count": len(self.input_items),
                    "max_items": self.max_items,
                    "is_filtered": True,
                    "combination_method": self.combination,
                    "criteria_count": len(self.criteria)
                }
            }
    
    def get_config_display(self) -> List[str]:
        """Get configuration display lines for pipeline output."""
        config_lines = super().get_config_display()
        
        config_lines.extend([
            f"Filter type: {self.filter_type}",
            f"Input items: {len(self.input_items)}",
            f"Combination: {self.combination}",
            f"Criteria: {len(self.criteria)}"
        ])
        
        # List individual criteria
        for i, criterion in enumerate(self.criteria, 1):
            config_lines.append(f"  {i}. {criterion.__class__.__name__}: {criterion.expression}")
        
        if self.score_weights:
            config_lines.append("Score weights:")
            for name, weight in self.score_weights.items():
                config_lines.append(f"  {name}: {weight}")
        
        if self.max_items is not None:
            config_lines.append(f"Max items: {self.max_items}")
        
        return config_lines