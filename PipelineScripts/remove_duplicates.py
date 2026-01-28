"""
RemoveDuplicates tool for filtering duplicate sequences from tool outputs.

Removes duplicate protein sequences by comparing against a reference pool,
useful for iterative optimization cycles to avoid recomputing identical sequences.
"""

import os
import pandas as pd
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo


class RemoveDuplicates(BaseConfig):
    """
    Tool for removing duplicate sequences from new results against historical data.
    
    Compares protein sequences in new_pool against all sequences in reference_pool,
    returning only sequences from new_pool that don't exist in reference_pool.
    """
    
    TOOL_NAME = "RemoveDuplicates"
    
    
    def __init__(self,
                 pool: Union[ToolOutput, StandardizedOutput, TableInfo, str],
                 history: Optional[Union[ToolOutput, StandardizedOutput, TableInfo, str]] = None,
                 compare: str = "sequence",
                 **kwargs):
        """
        Initialize RemoveDuplicates tool.

        Args:
            pool: Tool output, direct table, or file path containing sequences to check for duplicates
            history: Tool output, direct table, or file path containing historical sequences to compare against (None for first cycle)
            compare: Column name to compare for duplicates (e.g. "sequence")
            **kwargs: Additional parameters

        Examples:
            # First cycle - remove self-duplicates only
            unique_sequences = pipeline.add(RemoveDuplicates(
                pool=lmpnn_current,
                history=None,
                compare="sequence"
            ))

            # Subsequent cycles - remove against history (tool output)
            unique_sequences = pipeline.add(RemoveDuplicates(
                pool=lmpnn_current,
                history=all_sequences_seen,
                compare="sequence"
            ))

            # Elegant direct table access (TableInfo objects)
            unique_sequences = pipeline.add(RemoveDuplicates(
                pool=composer.o.tables.sequences,
                history=all_sequences_seen.o.tables.concatenated,
                compare="sequence"
            ))

            # Direct file path access (most elegant)
            unique_sequences = pipeline.add(RemoveDuplicates(
                pool=composer.tables.sequences,
                history=all_sequences_seen.tables.concatenated,
                compare="sequence"
            ))
        """
        self.pool = pool
        self.history = history
        self.compare = compare
        
        # Validate inputs
        if not isinstance(pool, (ToolOutput, StandardizedOutput, TableInfo, str)):
            raise ValueError("pool must be a ToolOutput, StandardizedOutput, TableInfo object, or file path string")

        if history is not None and not isinstance(history, (ToolOutput, StandardizedOutput, TableInfo, str)):
            raise ValueError("history must be a ToolOutput, StandardizedOutput, TableInfo object, file path string, or None")

        # compare can be any column name - no validation needed for specific values

        # Initialize base class
        super().__init__(**kwargs)

        # Set up dependencies
        if hasattr(pool, 'config'):
            self.dependencies.append(pool.config)
        if history is not None and hasattr(history, 'config'):
            self.dependencies.append(history.config)
    
    def validate_params(self):
        """Validate RemoveDuplicates parameters."""
        if not isinstance(self.pool, (ToolOutput, StandardizedOutput, TableInfo, str)):
            raise ValueError("pool must be a ToolOutput, StandardizedOutput, TableInfo object, or file path string")

        if self.history is not None and not isinstance(self.history, (ToolOutput, StandardizedOutput, TableInfo, str)):
            raise ValueError("history must be a ToolOutput, StandardizedOutput, TableInfo object, file path string, or None")
        
        # compare can be any column name - no validation needed for specific values
    
    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure inputs from tool outputs."""
        self.folders = pipeline_folders
        
        # Extract pool configuration  
        self.pool_config = self._extract_pool_config(self.pool, "pool")
        
        # Extract history configuration (handle None case)
        if self.history is not None:
            self.history_config = self._extract_pool_config(self.history, "history")
        else:
            # Create empty history config for first cycle
            self.history_config = {
                "prefix": "history",
                "structures": [],
                "structure_ids": [],
                "sequences": [],
                "sequence_ids": [],
                "output_folder": None,
                "sequence_csv": None
            }
    
    def _extract_pool_config(self, pool: Union[ToolOutput, StandardizedOutput, TableInfo, str], prefix: str) -> Dict[str, Any]:
        """Extract configuration from a pool tool output, direct table, or file path."""
        config = {
            "prefix": prefix,
            "structures": [],
            "structure_ids": [],
            "sequences": [],
            "sequence_ids": [],
            "output_folder": None,
            "sequence_csv": None
        }

        # Handle direct file path strings
        if isinstance(pool, str):
            config["sequence_csv"] = pool
            config["output_folder"] = os.path.dirname(pool)
            return config

        # Handle direct TableInfo objects
        if isinstance(pool, TableInfo):
            config["sequence_csv"] = pool.path
            config["output_folder"] = os.path.dirname(pool.path)
            return config

        if hasattr(pool, 'output_folder'):
            config["output_folder"] = pool.output_folder
            
            # Predict structure locations
            structure_dirs = [
                pool.output_folder,
                os.path.join(pool.output_folder, 'structures'),
                os.path.join(pool.output_folder, 'pdbs')
            ]
            config["structure_dirs"] = structure_dirs
            
            # Predict sequence CSV location
            sequence_files = [
                os.path.join(pool.output_folder, 'sequences.csv'),
                os.path.join(pool.output_folder, 'table.csv'),
                os.path.join(pool.output_folder, 'results.csv')
            ]
            config["sequence_files"] = sequence_files
        
        # Extract tables if available
        if hasattr(pool, 'tables'):
            tables = pool.tables
            
            if hasattr(tables, '_tables'):
                # Standard BioPipelines format - look for sequences/concatenated table
                for name, info in tables._tables.items():
                    if 'sequence' in name.lower() or 'concatenated' in name.lower() or name == 'table':
                        if hasattr(info, 'path'):
                            config["sequence_csv"] = info.path
                        else:
                            config["sequence_csv"] = str(info)
                        break

                # Fallback: if no specific pattern found, use first available table
                if config["sequence_csv"] is None and tables._tables:
                    first_name, first_info = next(iter(tables._tables.items()))
                    if hasattr(first_info, 'path'):
                        config["sequence_csv"] = first_info.path
                    else:
                        config["sequence_csv"] = str(first_info)

            elif isinstance(tables, dict):
                # Dict format - look for sequences/concatenated
                for name, info in tables.items():
                    if 'sequence' in name.lower() or 'concatenated' in name.lower() or name == 'table':
                        if isinstance(info, dict) and 'path' in info:
                            config["sequence_csv"] = info['path']
                        else:
                            config["sequence_csv"] = str(info)
                        break

                # Fallback: if no specific pattern found, use first available table
                if config["sequence_csv"] is None and tables:
                    first_name, first_info = next(iter(tables.items()))
                    if isinstance(first_info, dict) and 'path' in first_info:
                        config["sequence_csv"] = first_info['path']
                    else:
                        config["sequence_csv"] = str(first_info)
        
        return config
    
    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()
        
        config_lines.extend([
            f"COMPARE COLUMN: {self.compare}"
        ])
        
        return config_lines
    
    def generate_script(self, script_path: str) -> str:
        """
        Generate script to remove duplicate sequences.
        
        Args:
            script_path: Path where script should be written
            
        Returns:
            Script content as string
        """
        runtime_folder = os.path.dirname(script_path)
        output_folder = self.output_folder
        os.makedirs(output_folder, exist_ok=True)
        
        # Preserve input filename - must have sequences
        if not (hasattr(self.pool, 'sequences') and self.pool.sequences):
            raise ValueError("Pool must have sequences attribute with valid sequence files")
        
        input_filename = os.path.basename(self.pool.sequences[0])
        output_csv = os.path.join(output_folder, input_filename)
        
        # Create config file for duplicate removal
        config_file = os.path.join(output_folder, "remove_duplicates_config.json")
        config_data = {
            "pool_config": self.pool_config,
            "history_config": self.history_config,
            "compare": self.compare,
            "output_csv": output_csv
        }
        
        import json
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)
        
        # Generate script content
        script_content = f"""#!/bin/bash
# RemoveDuplicates execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Removing duplicate sequences"
echo "Compare column: {self.compare}"
echo "Output: {output_folder}"

# Run Python duplicate removal script
python "{os.path.join(self.folders['HelpScripts'], 'pipe_remove_duplicates.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Successfully removed duplicates"
    echo "Unique sequences written to: {output_csv}"
else
    echo "Error: Failed to remove duplicates"
    exit 1
fi

{self.generate_completion_check_footer()}
"""
        
        return script_content
    
    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after duplicate removal.
        
        Returns:
            Dictionary with output file paths
        """
        # RemoveDuplicates operates on tables only - it filters sequences
        # We cannot predict exact outputs until runtime, so use the same structure as input
        # but indicate potentially filtered results
        
        # Output files - only sequences table since this is a table operation
        sequences_csv = None
        sequence_ids = []
        sequences_files = []
        
        # Preserve input filename - must have sequences or valid tables
        if hasattr(self.pool, 'sequences') and self.pool.sequences:
            input_filename = os.path.basename(self.pool.sequences[0])
            sequences_csv = os.path.join(self.output_folder, input_filename)
            sequences_files = [sequences_csv]
        elif hasattr(self.pool, 'tables') and hasattr(self.pool.tables, '_tables'):
            # Find the sequences table from input
            for name, info in self.pool.tables._tables.items():
                if 'sequence' in name.lower():
                    if hasattr(info, 'path'):
                        sequences_csv = os.path.join(self.output_folder, os.path.basename(info.path))
                        sequences_files = [sequences_csv]
                        break
            else:
                raise ValueError("No sequences table found in pool.tables")
        else:
            raise ValueError("Pool must have either sequences attribute or valid tables with sequences")
        
        # Extract sequence IDs from input (we don't know which will remain after filtering)
        if hasattr(self.pool, 'sequence_ids'):
            sequence_ids = self.pool.sequence_ids.copy() if self.pool.sequence_ids else []
        
        # Create tables matching input structure
        tables = {}
        if hasattr(self.pool, 'tables') and hasattr(self.pool.tables, '_tables'):
            for name, info in self.pool.tables._tables.items():
                if 'sequence' in name.lower():
                    # Copy the table structure but point to output location
                    output_path = os.path.join(self.output_folder, os.path.basename(info.path if hasattr(info, 'path') else str(info)))
                    tables[name] = TableInfo(
                        name=name,
                        path=output_path,
                        columns=info.columns if hasattr(info, 'columns') else ["id", "sequence"],
                        description=f"filtered {info.description if hasattr(info, 'description') else 'sequences'}",
                        count="variable"
                    )
                    break
        
        # If no sequences table found in input, create default
        if not tables:
            tables["sequences"] = TableInfo(
                name="sequences",
                path=sequences_csv,
                columns=["id", "sequence"],
                description="filtered unique sequences after duplicate removal",
                count="variable"
            )
        
        # Add missing table for filtered out sequences (consistent with Filter)
        missing_csv = os.path.join(self.output_folder, "missing.csv")
        tables["missing"] = TableInfo(
            name="missing",
            path=missing_csv,
            columns=["id", "structure", "msa"],
            description="sequences filtered out due to duplication",
            count="variable"
        )
        
        # RemoveDuplicates only works with tables, no structures or compounds
        return {
            "sequences": sequences_files,
            "sequence_ids": sequence_ids,
            "tables": tables,
            "output_folder": self.output_folder
        }
    
    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                "compare": self.compare
            }
        })
        return base_dict