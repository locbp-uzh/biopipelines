"""
SliceTable tool for extracting a subset of rows from tool output tables.

Takes tool output and creates a new output containing only the first N rows
from each table, useful for sampling or testing subsets.
"""

import os
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


class SliceTable(BaseConfig):
    """
    Tool for extracting first N rows from tool output tables.

    Takes any tool output with tables and returns a new output containing
    only the first N rows from each table. Also copies associated structure
    files if they exist.
    """

    TOOL_NAME = "SliceTable"

    # Lazy path descriptors
    config_file = Path(lambda self: os.path.join(self.output_folder, "slice_config.json"))
    slice_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_slice_table.py"))

    def __init__(self,
                 input: StandardizedOutput,
                 n_rows: int = 10,
                 **kwargs):
        """
        Initialize SliceTable tool.

        Args:
            input: StandardizedOutput containing tables to slice
            n_rows: Number of rows to keep from the beginning of each table
            **kwargs: Additional parameters

        Examples:
            # Take first 10 sequences from LigandMPNN output
            first_ten = pipeline.add(SliceTable(
                input=lmpnn,
                n_rows=10
            ))

            # Take first 5 structures for quick testing
            sample = pipeline.add(SliceTable(
                input=boltz,
                n_rows=5
            ))
        """
        if not isinstance(input, StandardizedOutput):
            raise ValueError(f"input must be a StandardizedOutput object, got {type(input)}")

        if n_rows <= 0:
            raise ValueError("n_rows must be positive")

        self.input = input
        self.n_rows = n_rows

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate SliceTable parameters."""
        if not isinstance(self.input, StandardizedOutput):
            raise ValueError(f"input must be a StandardizedOutput object, got {type(self.input)}")

        if self.n_rows <= 0:
            raise ValueError("n_rows must be positive")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure inputs from tool output."""
        self.folders = pipeline_folders

        # Extract tables from input
        self.input_tables = {}
        if hasattr(self.input, 'tables') and hasattr(self.input.tables, '_tables'):
            self.input_tables = self.input.tables._tables

        if not self.input_tables:
            raise ValueError("Input tool must provide tables to slice")

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"N_ROWS: {self.n_rows}",
            f"INPUT DATASHEETS: {len(self.input_tables)}"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate SliceTable execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        script_content = "#!/bin/bash\n"
        script_content += "# SliceTable execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self.generate_script_run_slice_table()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_slice_table(self) -> str:
        """Generate the table slicing execution part of the script."""
        import json

        config_data = {
            "input_tables": {
                name: info.to_dict() if hasattr(info, 'to_dict') else {"path": str(info)}
                for name, info in self.input_tables.items()
            },
            "n_rows": self.n_rows,
            "output_folder": self.output_folder
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        return f"""echo "Slicing tables to first {self.n_rows} rows"
echo "Input tables: {len(self.input_tables)}"
echo "Output folder: {self.output_folder}"

python "{self.slice_py}" --config "{self.config_file}"

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after slicing."""
        tables = {}

        for name, info in self.input_tables.items():
            output_filename = os.path.basename(info.path if hasattr(info, 'path') else str(info))
            output_path = os.path.join(self.output_folder, output_filename)

            tables[name] = TableInfo(
                name=name,
                path=output_path,
                columns=info.columns if hasattr(info, 'columns') else [],
                description=f"First {self.n_rows} rows from {info.description if hasattr(info, 'description') else name}",
                count=self.n_rows
            )

        return {
            "structures": DataStream.empty("structures", "pdb"),
            "sequences": DataStream.empty("sequences", "fasta"),
            "compounds": DataStream.empty("compounds", "sdf"),
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                "n_rows": self.n_rows
            }
        })
        return base_dict