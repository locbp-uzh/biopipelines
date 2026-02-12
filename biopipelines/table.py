# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Table entity for loading existing CSV files into the pipeline.

Reads column headers from the CSV and provides a TableInfo for downstream tools.
"""

import os
import pandas as pd
from typing import Dict, List, Any, Optional

try:
    from .base_config import BaseConfig, TableInfo
    from .datastream import DataStream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, TableInfo
    from datastream import DataStream


class Table(BaseConfig):
    """
    Entity for loading existing CSV files into the pipeline.

    Reads the CSV at pipeline time to extract column headers, then provides
    a TableInfo that downstream tools can reference.

    Example:
        # Load existing metrics table
        metrics = Table("/path/to/metrics.csv")

        # Use column references in downstream tools
        Panda(
            tables=[metrics.tables.data],
            operations=[Panda.sort("score", ascending=False)]
        )

        # Access specific column for per-structure data
        ProteinMPNN(
            structures=proteins,
            redesigned=(metrics.tables.data, "designed_positions")
        )
    """

    TOOL_NAME = "Table"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== Table ==="
echo "Uses biopipelines environment (no additional installation needed)."
echo "=== Table ready ==="
"""

    def __init__(self,
                 path: str,
                 name: str = "data",
                 description: str = "",
                 **kwargs):
        """
        Initialize Table entity.

        Args:
            path: Path to existing CSV file
            name: Name for the table (default: "data")
            description: Description of the table contents
            **kwargs: Additional parameters
        """
        if not os.path.exists(path):
            raise FileNotFoundError(f"Table file not found: {path}")

        self.table_path = os.path.abspath(path)
        self.table_name = name
        self.table_description = description

        # Read CSV to get columns and count
        df = pd.read_csv(self.table_path)
        self.table_columns = list(df.columns)
        self.table_count = len(df)

        print(f"  Loaded table '{self.table_name}' with {self.table_count} rows and {len(self.table_columns)} columns")

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate Table parameters."""
        if not os.path.exists(self.table_path):
            raise FileNotFoundError(f"Table file not found: {self.table_path}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input parameters."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()
        config_lines.extend([
            f"PATH: {self.table_path}",
            f"NAME: {self.table_name}",
            f"ROWS: {self.table_count}",
            f"COLUMNS: {', '.join(self.table_columns)}"
        ])
        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate script - Table is a no-op since file already exists."""
        script_content = "#!/bin/bash\n"
        script_content += "# Table entity - no execution needed (file already exists)\n"
        script_content += self.generate_completion_check_header()
        script_content += f'echo "Table \'{self.table_name}\' already exists at: {self.table_path}"\n'
        script_content += f'echo "Columns: {", ".join(self.table_columns)}"\n'
        script_content += f'echo "Rows: {self.table_count}"\n'
        script_content += self.generate_completion_check_footer()
        return script_content

    def get_output_files(self) -> Dict[str, Any]:
        """Get output files - returns the loaded table as TableInfo."""
        tables = {
            self.table_name: TableInfo(
                name=self.table_name,
                path=self.table_path,
                columns=self.table_columns,
                description=self.table_description,
                count=self.table_count
            )
        }

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
            "table_params": {
                "path": self.table_path,
                "name": self.table_name,
                "columns": self.table_columns,
                "count": self.table_count,
                "description": self.table_description
            }
        })
        return base_dict
