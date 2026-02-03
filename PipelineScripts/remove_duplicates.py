"""
RemoveDuplicates tool for filtering duplicate sequences from tool outputs.

Removes duplicate protein sequences by comparing against a reference pool,
useful for iterative optimization cycles to avoid recomputing identical sequences.
"""

import os
import json
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream


class RemoveDuplicates(BaseConfig):
    """
    Tool for removing duplicate sequences from new results against historical data.

    Compares protein sequences in new_pool against all sequences in reference_pool,
    returning only sequences from new_pool that don't exist in reference_pool.
    """

    TOOL_NAME = "RemoveDuplicates"

    # Lazy path descriptors
    config_file = Path(lambda self: os.path.join(self.output_folder, "remove_duplicates_config.json"))
    missing_csv = Path(lambda self: os.path.join(self.output_folder, "missing.csv"))
    helper_script = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_remove_duplicates.py"))

    def __init__(self,
                 pool: Union[StandardizedOutput, TableInfo, str],
                 history: Optional[Union[StandardizedOutput, TableInfo, str]] = None,
                 compare: str = "sequence",
                 **kwargs):
        """
        Initialize RemoveDuplicates tool.

        Args:
            pool: StandardizedOutput, direct table, or file path containing sequences to check for duplicates
            history: StandardizedOutput, direct table, or file path containing historical sequences to compare against (None for first cycle)
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
                pool=composer.tables.sequences,
                history=all_sequences_seen.tables.concatenated,
                compare="sequence"
            ))
        """
        self.pool = pool
        self.history = history
        self.compare = compare

        # Initialize base class
        super().__init__(**kwargs)

    def validate_params(self):
        """Validate RemoveDuplicates parameters."""
        if not isinstance(self.pool, (StandardizedOutput, TableInfo, str)):
            raise ValueError("pool must be a StandardizedOutput, TableInfo object, or file path string")

        if self.history is not None and not isinstance(self.history, (StandardizedOutput, TableInfo, str)):
            raise ValueError("history must be a StandardizedOutput, TableInfo object, file path string, or None")

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
                "output_folder": None,
                "sequence_csv": None
            }

        # Get the output filename from input
        self._determine_output_filename()

    def _extract_pool_config(self, pool: Union[StandardizedOutput, TableInfo, str], prefix: str) -> Dict[str, Any]:
        """Extract configuration from a pool tool output, direct table, or file path."""
        config = {
            "prefix": prefix,
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

        # Extract tables if available
        if hasattr(pool, 'tables'):
            tables = pool.tables

            if hasattr(tables, '_tables'):
                # Standard BioPipelines format - look for sequences/concatenated table
                for name, info in tables._tables.items():
                    if 'sequence' in name.lower() or 'concatenated' in name.lower() or name == 'table':
                        if hasattr(info, 'path'):
                            config["sequence_csv"] = info.path
                        break

                # If no specific pattern found, use first available table
                if config["sequence_csv"] is None and tables._tables:
                    first_name, first_info = next(iter(tables._tables.items()))
                    if hasattr(first_info, 'path'):
                        config["sequence_csv"] = first_info.path

        # Also check sequences DataStream
        if config["sequence_csv"] is None and hasattr(pool, 'streams') and pool.streams.sequences:
            sequences = pool.streams.sequences
            if sequences and hasattr(sequences, 'map_table') and sequences.map_table:
                config["sequence_csv"] = sequences.map_table
            elif sequences and hasattr(sequences, 'files') and sequences.files:
                config["sequence_csv"] = sequences.files[0]

        return config

    def _determine_output_filename(self):
        """Determine the output filename from input."""
        # Try to get filename from pool config
        if self.pool_config.get("sequence_csv"):
            self.output_filename = os.path.basename(self.pool_config["sequence_csv"])
        elif hasattr(self.pool, 'sequences') and self.pool.sequences:
            if hasattr(self.pool.sequences, 'files') and self.pool.sequences.files:
                self.output_filename = os.path.basename(self.pool.sequences.files[0])
            elif hasattr(self.pool.sequences, 'map_table') and self.pool.sequences.map_table:
                self.output_filename = os.path.basename(self.pool.sequences.map_table)
            else:
                self.output_filename = "sequences.csv"
        else:
            self.output_filename = "sequences.csv"

    def _get_output_csv_path(self) -> str:
        """Get the path for the output CSV."""
        return os.path.join(self.output_folder, self.output_filename)

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"COMPARE COLUMN: {self.compare}"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate RemoveDuplicates execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        output_csv = self._get_output_csv_path()

        # Create config file for duplicate removal
        config_data = {
            "pool_config": self.pool_config,
            "history_config": self.history_config,
            "compare": self.compare,
            "output_csv": output_csv
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        script_content = "#!/bin/bash\n"
        script_content += "# RemoveDuplicates execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()

        script_content += f"""echo "Removing duplicate sequences"
echo "Compare column: {self.compare}"
echo "Output: {self.output_folder}"

python "{self.helper_script}" \\
  --config "{self.config_file}"

if [ $? -eq 0 ]; then
    echo "Successfully removed duplicates"
    echo "Unique sequences written to: {output_csv}"
else
    echo "Error: Failed to remove duplicates"
    exit 1
fi

"""
        script_content += self.generate_completion_check_footer()

        return script_content

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after duplicate removal."""
        output_csv = self._get_output_csv_path()

        # Get columns from input if available
        columns = ["id", "sequence"]
        if isinstance(self.pool, TableInfo) and hasattr(self.pool, 'columns'):
            columns = self.pool.columns
        elif hasattr(self.pool, 'tables') and hasattr(self.pool.tables, '_tables'):
            for name, info in self.pool.tables._tables.items():
                if 'sequence' in name.lower():
                    if hasattr(info, 'columns') and info.columns:
                        columns = info.columns
                    break

        # Create tables
        tables = {
            "sequences": TableInfo(
                name="sequences",
                path=output_csv,
                columns=columns,
                description="Filtered unique sequences after duplicate removal",
                count="variable"
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "structure", "msa"],
                description="Sequences filtered out due to duplication",
                count="variable"
            )
        }

        # Create sequences DataStream
        sequences = DataStream(
            name="sequences",
            ids=[],  # Will be populated at runtime
            files=[output_csv],
            map_table=output_csv,
            format="csv"
        )

        # RemoveDuplicates only works with tables, no structures or compounds
        return {
            "structures": DataStream.empty("structures", "pdb"),
            "sequences": sequences,
            "compounds": DataStream.empty("compounds", "sdf"),
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
