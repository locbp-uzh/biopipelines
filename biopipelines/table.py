# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Table entity for loading existing CSV or Excel files into the pipeline.

Reads column headers from the file and provides a TableInfo for downstream tools.
Excel files (.xlsx/.xls) are converted to CSV internally before use.
"""

import os
import pandas as pd
from typing import Dict, List, Any, Optional

try:
    from .base_config import (
        BaseConfig, TableInfo,
        _validate_freeform_string, _escape_for_double_quotes,
    )
    from .datastream import DataStream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import (
        BaseConfig, TableInfo,
        _validate_freeform_string, _escape_for_double_quotes,
    )
    from datastream import DataStream


class Table(BaseConfig):
    """
    Entity for loading existing CSV files into the pipeline.

    Reads the CSV at configuration time to extract column headers, then provides
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
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== Table ==="
echo "Uses biopipelines environment (no additional installation needed)."
touch "$INSTALL_SUCCESS"
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
            path: Path to existing CSV or Excel (.xlsx/.xls) file.
                  Can be absolute, relative to the current directory, or a
                  filename inside the pipeline's tables/ folder.
                  Excel files are converted to CSV internally.
            name: Name for the table (default: "data")
            description: Description of the table contents
            **kwargs: Additional parameters
        """
        ext = os.path.splitext(path)[1].lower()
        if ext not in ('.csv', '.xlsx', '.xls'):
            raise ValueError(f"Table only accepts CSV or Excel files (.csv, .xlsx, .xls). Got: '{path}'")

        self.table_name = name
        self.table_description = description
        self._pending_filename = None
        # When the source is Excel, defer the CSV write until configure_inputs
        # (when tables_folder exists) so we never mutate the user's input dir.
        self._pending_excel_df = None
        self._excel_source_path = None

        if os.path.isfile(path):
            self._load_from_path(os.path.abspath(path))
        else:
            # Defer to configure_inputs to try joining with the tables/ folder
            self._pending_filename = path
            self.table_path = path
            self.table_columns = []
            self.table_count = 0

        super().__init__(**kwargs)

    def _load_from_path(self, abs_path: str):
        """Load the table from an absolute path; defer Excel-to-CSV until configure_inputs."""
        ext = os.path.splitext(abs_path)[1].lower()
        if ext in ('.xlsx', '.xls'):
            df = pd.read_excel(abs_path)
            self._pending_excel_df = df
            self._excel_source_path = abs_path
            # table_path will be assigned in configure_inputs once tables_folder
            # is known; provisional value lets get_config_display run pre-config.
            self.table_path = abs_path
        else:
            df = pd.read_csv(abs_path)
            self.table_path = abs_path

        self.table_columns = list(df.columns)
        self.table_count = len(df)

        print(f"  Loaded table '{self.table_name}' with {self.table_count} rows and {len(self.table_columns)} columns")

    def validate_params(self):
        """Validate Table parameters."""
        # Pending-filename tables are validated after configure_inputs resolves them
        if not self._pending_filename and not os.path.exists(self.table_path):
            raise FileNotFoundError(f"Table file not found: {self.table_path}")

        _validate_freeform_string("name", self.table_name)
        _validate_freeform_string("description", self.table_description)
        # table_path is a filesystem path, not a user-supplied free-form
        # string; Windows paths legitimately contain backslashes. We escape
        # it at bash-emission time instead of rejecting it outright.

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input parameters."""
        self.folders = pipeline_folders

        # Resolve a pending path against the tables/ folder
        if self._pending_filename:
            tables_folder = pipeline_folders.get("tables", "")
            candidate = os.path.join(tables_folder, self._pending_filename)
            if not os.path.isfile(candidate):
                raise FileNotFoundError(
                    f"Table file not found: '{self._pending_filename}'. "
                    f"Tried as-is and as '{candidate}'."
                )
            self._load_from_path(os.path.abspath(candidate))
            self._pending_filename = None

        # For Excel inputs: materialize the CSV inside the pipeline's tables/
        # folder so we don't write next to the user's input file.
        if self._pending_excel_df is not None:
            os.makedirs(self.tables_folder, exist_ok=True)
            base = os.path.splitext(os.path.basename(self._excel_source_path))[0]
            csv_path = os.path.join(self.tables_folder, f"{base}.csv")
            self._pending_excel_df.to_csv(csv_path, index=False)
            self.table_path = csv_path
            self._pending_excel_df = None
            print(f"  Converted Excel to CSV: {csv_path}")

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
        # Escape the filesystem path before interpolating into the echo
        # line (Windows paths contain backslashes that would otherwise be
        # interpreted as bash escape sequences).
        safe_path = _escape_for_double_quotes(self.table_path)
        script_content = "#!/bin/bash\n"
        script_content += "# Table entity - no execution needed (file already exists)\n"
        script_content += self.generate_completion_check_header()
        script_content += f'echo "Table \'{self.table_name}\' already exists at: {safe_path}"\n'
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
            )
        }

        return {
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
                "description": self.table_description
            }
        })
        return base_dict
