# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Plot tool with declarative operation-based API for data visualization.

Creates PNG plots from CSV tables using a sequence of operations like
Scatter, Histogram, Bar, Column. Supports multiple data sources for
comparison and table column references for data access.
"""

import os
import json
from typing import Dict, List, Any, Optional, Union, Tuple

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


class PlotOperation:
    """Base class for plot operations."""

    def __init__(self, op_type: str, **kwargs):
        self.op_type = op_type
        self.params = kwargs

    def to_dict(self) -> Dict[str, Any]:
        """Convert operation to dictionary for serialization."""
        return {"op": self.op_type, **self.params}


class Plot(BaseConfig):
    """
    Plot tool with declarative operation-based API for data visualization.

    Uses a sequence of operations (Scatter, Histogram, Bar, Column) to create
    plots from CSV data. Supports multiple data sources for
    comparison and table column references for data access.

    Example:
        # Single histogram
        Plot(
            Plot.Histogram(data=analysis.tables.results, x="pLDDT", bins=20, title="pLDDT Distribution")
        )

        # Multiple plots in one folder
        Plot(
            Plot.Histogram(data=analysis, x="pLDDT", title="pLDDT Distribution"),
            Plot.Histogram(data=analysis, x="affinity", title="Affinity Distribution"),
            Plot.Scatter(data=analysis, x="pLDDT", y="affinity", title="pLDDT vs Affinity"),
        )

        # Scatter with color grouping
        Plot(
            Plot.Scatter(data=analysis, x="pLDDT", y="affinity", color="source", title="By Source")
        )

        # Column plot (Prism-style) - comparing same metric across tools
        Plot(
            Plot.Column(
                data=[tool1.tables.results, tool2.tables.results, tool3.tables.results],
                y="pLDDT",
                labels=["RFdiffusion", "Boltz", "ESMFold"],
                title="pLDDT Comparison"
            )
        )

        # HeatMap - correlation matrix
        Plot(
            Plot.HeatMap(data=analysis, columns=["pLDDT", "affinity", "contacts"], title="Metric Correlations")
        )

        # HeatMap - pivot table
        Plot(
            Plot.HeatMap(data=analysis, x="model", y="condition", value="score", title="Model Performance")
        )
    """

    TOOL_NAME = "Plot"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== Plot ==="
echo "Uses biopipelines environment (no additional installation needed)."
echo "=== Plot ready ==="
"""

    # Lazy path descriptors
    config_file = Path(lambda self: os.path.join(self.output_folder, "plot_config.json"))
    metadata_csv = Path(lambda self: os.path.join(self.output_folder, "plot_metadata.csv"))
    plot_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_plot.py"))

    # --- Static methods for creating operations ---

    @staticmethod
    def Scatter(data: Union[StandardizedOutput, TableInfo],
                x: str,
                y: str,
                color: str = None,
                title: str = None,
                xlabel: str = None,
                ylabel: str = None,
                x_name: str = None,
                y_name: str = None,
                figsize: Tuple[int, int] = (800, 600),
                dpi: int = 100,
                x_tick_rotation: float = 0,
                y_tick_rotation: float = 0,
                grid: bool = True,
                color_legend_loc: str = "upper right",
                color_legend_outside: bool = False) -> PlotOperation:
        """
        Create X vs Y scatter plot, optionally colored by a third column.

        Args:
            data: Data source (tool output or table)
            x: Column name for X axis
            y: Column name for Y axis
            color: Column name for color grouping (optional)
            title: Plot title (optional)
            xlabel: X axis label (defaults to x_name or column name)
            ylabel: Y axis label (defaults to y_name or column name)
            x_name: Display name for X variable (alternative to xlabel)
            y_name: Display name for Y variable (alternative to ylabel)
            figsize: Figure size as (width, height) in pixels (default: (800, 600))
            dpi: Resolution in dots per inch (default: 100)
            x_tick_rotation: Rotation angle for x-axis tick labels in degrees
            y_tick_rotation: Rotation angle for y-axis tick labels in degrees
            grid: Show grid lines (default True)
            color_legend_loc: Color legend location ("upper right", "upper left", etc.)
            color_legend_outside: Place legend outside plot area to the right (default False)

        Returns:
            PlotOperation for scatter plot
        """
        return PlotOperation("scatter", data=data, x=x, y=y, color=color,
                            title=title, xlabel=xlabel, ylabel=ylabel,
                            x_name=x_name, y_name=y_name, figsize=figsize, dpi=dpi,
                            x_tick_rotation=x_tick_rotation, y_tick_rotation=y_tick_rotation,
                            grid=grid, color_legend_loc=color_legend_loc,
                            color_legend_outside=color_legend_outside)

    @staticmethod
    def Histogram(data: Union[StandardizedOutput, TableInfo],
                  x: str,
                  bins: int = 20,
                  title: str = None,
                  xlabel: str = None,
                  ylabel: str = None,
                  x_name: str = None,
                  y_name: str = None,
                  figsize: Tuple[int, int] = (800, 600),
                  dpi: int = 100) -> PlotOperation:
        """
        Create distribution histogram of a single column.

        Args:
            data: Data source (tool output or table)
            x: Column name to plot distribution of
            bins: Number of histogram bins
            title: Plot title (optional)
            xlabel: X axis label (defaults to x_name or column name)
            ylabel: Y axis label (defaults to y_name or "Count")
            x_name: Display name for X variable (alternative to xlabel)
            y_name: Display name for Y variable (alternative to ylabel)
            figsize: Figure size as (width, height) in pixels (default: (800, 600))
            dpi: Resolution in dots per inch (default: 100)

        Returns:
            PlotOperation for histogram
        """
        return PlotOperation("histogram", data=data, x=x, bins=bins,
                            title=title, xlabel=xlabel, ylabel=ylabel,
                            x_name=x_name, y_name=y_name, figsize=figsize, dpi=dpi)

    @staticmethod
    def Bar(data: Union[StandardizedOutput, TableInfo],
            x: str,
            y: str,
            y_right: str = None,
            title: str = None,
            xlabel: str = None,
            ylabel: str = None,
            ylabel_right: str = None,
            x_name: str = None,
            y_name: str = None,
            color_left: str = None,
            color_right: str = None,
            figsize: Tuple[int, int] = (800, 600),
            dpi: int = 100,
            x_tick_rotation: float = 0,
            y_tick_rotation: float = 0,
            grid: bool = True,
            legend_loc: str = "upper right",
            legend_outside: bool = False) -> PlotOperation:
        """
        Create bar chart for categorical X axis.

        Supports an optional second y column plotted on a secondary (right) y-axis,
        producing side-by-side bars at each x position with independent scales.

        Args:
            data: Data source (tool output or table)
            x: Column name for categories (X axis)
            y: Column name for values (left Y axis)
            y_right: Column name for values on secondary right Y axis (optional)
            title: Plot title (optional)
            xlabel: X axis label (defaults to x_name or column name)
            ylabel: Left Y axis label (defaults to y_name or column name)
            ylabel_right: Right Y axis label (defaults to y_right column name)
            x_name: Display name for X variable (alternative to xlabel)
            y_name: Display name for Y variable (alternative to ylabel)
            color_left: Color for left-axis bars (default: steelblue)
            color_right: Color for right-axis bars (default: coral)
            figsize: Figure size as (width, height) in pixels (default: (800, 600))
            dpi: Resolution in dots per inch (default: 100)
            x_tick_rotation: Rotation angle for x-axis tick labels in degrees
            y_tick_rotation: Rotation angle for y-axis tick labels in degrees
            grid: Show grid lines (default True)
            legend_loc: Legend location ("upper right", "upper left", etc.)
            legend_outside: Place legend outside plot area to the right (default False)

        Returns:
            PlotOperation for bar chart
        """
        return PlotOperation("bar", data=data, x=x, y=y, y_right=y_right,
                            title=title, xlabel=xlabel, ylabel=ylabel,
                            ylabel_right=ylabel_right,
                            x_name=x_name, y_name=y_name,
                            color_left=color_left, color_right=color_right,
                            figsize=figsize, dpi=dpi,
                            x_tick_rotation=x_tick_rotation, y_tick_rotation=y_tick_rotation,
                            grid=grid,
                            legend_loc=legend_loc, legend_outside=legend_outside)

    @staticmethod
    def Column(data: List[Union[StandardizedOutput, TableInfo]],
               y: str,
               labels: List[str] = None,
               title: str = None,
               xlabel: str = None,
               ylabel: str = None,
               x_name: str = None,
               y_name: str = None,
               style: str = "column",
               show_error: str = "sd",
               figsize: Tuple[int, int] = (800, 600),
               dpi: int = 100,
               color_groups: List[str] = None,
               colors: List[str] = None,
               color_legend_title: str = None,
               color_legend_loc: str = "upper right",
               color_legend_outside: bool = False,
               x_tick_rotation: float = 0,
               y_tick_rotation: float = 0,
               grid: bool = True) -> PlotOperation:
        """
        Create column plot comparing same metric across multiple data sources.

        Displays grouped data with configurable style, commonly used for comparing
        metrics between conditions or tools.

        Args:
            data: List of data sources [tool1.tables.x, tool2.tables.x, ...]
            y: Column name to compare across sources
            labels: Group labels for x-axis (defaults to source filenames)
            title: Plot title (optional)
            xlabel: X axis label (defaults to x_name)
            ylabel: Y axis label (defaults to y_name or column name)
            x_name: Display name for X axis (alternative to xlabel)
            y_name: Display name for Y variable (alternative to ylabel)
            style: Plot style - "column" (bars+points), "simple_bar" (bars only),
                   "scatter" (points only), "box" (box and whiskers), "floating_bar" (mean±error bars)
            show_error: Error bar type: "sd" (std dev), "sem" (std error), "ci" (95% CI), or None
            figsize: Figure size as (width, height) in pixels (default: (800, 600))
            dpi: Resolution in dots per inch (default: 100)
            color_groups: List of group names for coloring (e.g., protein names). Same color for same group.
            colors: Explicit list of colors (hex or named). If not provided, uses a color palette.
            color_legend_title: Title for the color legend (e.g., "Protein")
            color_legend_loc: Color legend location ("upper right", "upper left", "lower right", "lower left", etc.)
            color_legend_outside: Place legend outside plot area to the right (default False)
            x_tick_rotation: Rotation angle for x-axis tick labels in degrees
            y_tick_rotation: Rotation angle for y-axis tick labels in degrees
            grid: Show grid lines (default True)

        Returns:
            PlotOperation for column plot
        """
        return PlotOperation("column", data=data, y=y, labels=labels,
                            title=title, xlabel=xlabel, ylabel=ylabel,
                            x_name=x_name, y_name=y_name,
                            style=style, show_error=show_error, figsize=figsize, dpi=dpi,
                            color_groups=color_groups, colors=colors,
                            color_legend_title=color_legend_title, color_legend_loc=color_legend_loc,
                            color_legend_outside=color_legend_outside,
                            x_tick_rotation=x_tick_rotation, y_tick_rotation=y_tick_rotation,
                            grid=grid)

    @staticmethod
    def HeatMap(data: Union[StandardizedOutput, TableInfo],
                x: str = None,
                y: str = None,
                value: str = None,
                columns: List[str] = None,
                title: str = None,
                xlabel: str = None,
                ylabel: str = None,
                x_name: str = None,
                y_name: str = None,
                cmap: str = "viridis",
                annotate: bool = True,
                figsize: Tuple[int, int] = (1000, 800),
                dpi: int = 100) -> PlotOperation:
        """
        Create a heatmap visualization.

        Two modes:
        1. Pivot mode: Provide x, y, value columns to create a pivot table heatmap
        2. Correlation mode: Provide columns list to show correlation matrix

        Args:
            data: Data source (tool output or table)
            x: Column for X axis categories (pivot mode)
            y: Column for Y axis categories (pivot mode)
            value: Column for cell values (pivot mode)
            columns: List of column names for correlation matrix (correlation mode)
            title: Plot title (optional)
            xlabel: X axis label (defaults to x_name or column name)
            ylabel: Y axis label (defaults to y_name or column name)
            x_name: Display name for X axis (alternative to xlabel)
            y_name: Display name for Y axis (alternative to ylabel)
            cmap: Colormap name (default: "viridis")
            annotate: Show values in cells (default: True)
            figsize: Figure size as (width, height) in pixels (default: (1000, 800))
            dpi: Resolution in dots per inch (default: 100)

        Returns:
            PlotOperation for heatmap

        Examples:
            # Correlation matrix
            Plot.HeatMap(data=results, columns=["pLDDT", "affinity", "contacts"])

            # Pivot table heatmap
            Plot.HeatMap(data=results, x="model", y="condition", value="score")
        """
        return PlotOperation("heatmap", data=data, x=x, y=y, value=value,
                            columns=columns, title=title, xlabel=xlabel, ylabel=ylabel,
                            x_name=x_name, y_name=y_name, cmap=cmap,
                            annotate=annotate, figsize=figsize, dpi=dpi)

    # --- Instance methods ---

    def __init__(self, *args, **kwargs):
        """
        Initialize Plot tool with a sequence of plot operations.

        Args:
            *args: Sequence of plot operations (Scatter, Histogram, Bar, Column)
            **kwargs: Additional configuration parameters

        Example:
            Plot(
                Plot.Histogram(data=results, x="pLDDT"),
                Plot.Scatter(data=results, x="pLDDT", y="affinity")
            )
        """
        self.operations = list(args)

        # Extract data sources from operations
        self._data_sources = []
        self._extract_dependencies()

        super().__init__(**kwargs)

    def _extract_dependencies(self):
        """Extract data sources from operations."""
        for op in self.operations:
            data = op.params.get("data")

            if data is not None:
                if isinstance(data, list):
                    # Column plot with multiple data sources
                    for d in data:
                        if d not in self._data_sources:
                            self._data_sources.append(d)
                elif data not in self._data_sources:
                    self._data_sources.append(data)

    def validate_params(self):
        """Validate Plot parameters."""
        if not self.operations:
            raise ValueError("At least one plot operation must be provided")

        # Check that all operations are valid PlotOperation objects
        for i, op in enumerate(self.operations):
            if not isinstance(op, PlotOperation):
                raise ValueError(f"Operation {i} is not a PlotOperation: {type(op)}")

            # Validate operation-specific requirements
            if op.op_type == "scatter":
                if not op.params.get("x") or not op.params.get("y"):
                    raise ValueError(f"Scatter operation {i} requires both 'x' and 'y' parameters")
            elif op.op_type == "histogram":
                if not op.params.get("x"):
                    raise ValueError(f"Histogram operation {i} requires 'x' parameter")
            elif op.op_type == "bar":
                if not op.params.get("x") or not op.params.get("y"):
                    raise ValueError(f"Bar operation {i} requires both 'x' and 'y' parameters")
            elif op.op_type == "column":
                if not op.params.get("y"):
                    raise ValueError(f"Column operation {i} requires 'y' parameter")
                data = op.params.get("data")
                if not isinstance(data, list) or len(data) < 1:
                    raise ValueError(f"Column operation {i} requires 'data' to be a list with at least one source")
            elif op.op_type == "heatmap":
                # Either pivot mode (x, y, value) or correlation mode (columns)
                has_pivot = op.params.get("x") and op.params.get("y") and op.params.get("value")
                has_columns = op.params.get("columns")
                if not has_pivot and not has_columns:
                    raise ValueError(f"HeatMap operation {i} requires either (x, y, value) for pivot mode or 'columns' for correlation mode")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input sources from pipeline context."""
        self.folders = pipeline_folders
        self.input_sources = {}

        # Track data sources
        for i, source in enumerate(self._data_sources):
            if hasattr(source, 'info'):
                # TableInfo
                self.input_sources[f"data_{i}"] = {
                    'path': source.info.path,
                    'name': source.info.name
                }
            elif hasattr(source, 'output_folder'):
                # StandardizedOutput or ToolOutput
                self.input_sources[f"data_{i}"] = {
                    'output_folder': source.output_folder
                }

    def _get_table_path(self, data) -> str:
        """Extract table path from various data source types."""
        if isinstance(data, TableInfo):
            return data.info.path
        elif hasattr(data, 'tables'):
            tables = data.tables
            if hasattr(tables, '_tables'):
                # Get first table
                first_name, info = next(iter(tables._tables.items()))
                return info.info.path
            elif isinstance(tables, dict):
                first_name, info = next(iter(tables.items()))
                if isinstance(info, dict) and 'path' in info:
                    return info['path']
                elif hasattr(info, 'info'):
                    return info.info.path
                return str(info)
        raise ValueError(f"Cannot extract table path from: {type(data)}")

    def _get_table_name(self, data) -> str:
        """Extract table name from various data source types."""
        if isinstance(data, TableInfo):
            return data.info.name
        elif hasattr(data, 'tables'):
            tables = data.tables
            if hasattr(tables, '_tables'):
                first_name, _ = next(iter(tables._tables.items()))
                return first_name
        # Fallback to filename
        path = self._get_table_path(data)
        return os.path.splitext(os.path.basename(path))[0]

    def _serialize_data_source(self, data) -> Dict[str, Any]:
        """Serialize a data source to dictionary for JSON config."""
        if isinstance(data, TableInfo):
            return {
                "type": "table_info",
                "path": data.info.path,
                "name": data.info.name
            }
        elif hasattr(data, 'output_folder'):
            # StandardizedOutput or similar
            table_path = self._get_table_path(data)
            table_name = self._get_table_name(data)
            return {
                "type": "standardized_output",
                "path": table_path,
                "name": table_name,
                "output_folder": data.output_folder
            }
        else:
            raise ValueError(f"Unsupported data source type: {type(data)}")

    def _serialize_operation(self, op: PlotOperation) -> Dict[str, Any]:
        """Serialize an operation to a dictionary for JSON config."""
        result = {"op": op.op_type}

        for key, value in op.params.items():
            if key == "data":
                if isinstance(value, list):
                    # Multiple data sources (Column plot)
                    result[key] = [self._serialize_data_source(d) for d in value]
                elif value is not None:
                    result[key] = self._serialize_data_source(value)
            elif isinstance(value, tuple):
                # Convert figsize tuple to list for JSON
                result[key] = list(value)
            elif value is not None:
                result[key] = value

        return result

    def _generate_plot_filename(self, op: PlotOperation, index: int) -> str:
        """Generate output filename for a plot operation."""
        title = op.params.get("title")
        if title:
            # Sanitize title for filename
            safe_title = "".join(c if c.isalnum() or c in "_ -" else "_" for c in title)
            safe_title = safe_title.replace(" ", "_").lower()
            return f"{safe_title}.png"
        else:
            return f"{op.op_type}_{index + 1}.png"

    def generate_script(self, script_path: str) -> str:
        """Generate Plot execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        script_content = "#!/bin/bash\n"
        script_content += "# Plot execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self.generate_script_run_plot()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_plot(self) -> str:
        """Generate the plot creation part of the script."""
        config = {
            "operations": [self._serialize_operation(op) for op in self.operations],
            "output_folder": self.output_folder,
            "plot_filenames": [self._generate_plot_filename(op, i) for i, op in enumerate(self.operations)]
        }

        # Write config file at configuration time (not execution time)
        with open(self.config_file, 'w') as f:
            json.dump(config, f, indent=2)

        return f"""echo "Creating plots..."
echo "Output folder: {self.output_folder}"

python "{self.plot_py}" --config "{self.config_file}"

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files."""
        # Generate plot filenames
        plot_ids = []
        plot_files = []
        for i, op in enumerate(self.operations):
            plot_ids.append(f"plot_{i + 1}")
            filename = self._generate_plot_filename(op, i)
            plot_files.append(os.path.join(self.output_folder, filename))

        plots_stream = DataStream(
            name="plots",
            ids=plot_ids,
            files=plot_files,
            format="png"
        )

        tables = {
            "metadata": TableInfo(
                name="metadata",
                path=self.metadata_csv,
                columns=["filename", "type", "title", "x_column", "y_column", "data_sources"],
                description="Summary of all generated plots",
                count=len(self.operations)
            )
        }

        # Add a table entry for each plot's exported CSV data
        for i, op in enumerate(self.operations):
            png_filename = self._generate_plot_filename(op, i)
            csv_filename = os.path.splitext(png_filename)[0] + ".csv"
            csv_path = os.path.join(self.output_folder, csv_filename)
            table_name = os.path.splitext(png_filename)[0]
            tables[table_name] = TableInfo(
                name=table_name,
                path=csv_path,
                columns=[],
                description=f"Data for plot: {op.params.get('title', op.op_type)}",
                count="variable"
            )

        return {
            "plots": plots_stream,
            "tables": tables,
            "output_folder": self.output_folder,
        }

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"OPERATIONS: {len(self.operations)}",
            f"DATA SOURCES: {len(self._data_sources)}"
        ])

        # Show operation summary
        op_types = [op.op_type for op in self.operations]
        config_lines.append(f"PLOT TYPES: {', '.join(op_types)}")

        return config_lines

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration to dictionary."""
        base_dict = super().to_dict()
        base_dict.update({
            "plot_params": {
                "num_operations": len(self.operations),
                "operation_types": [op.op_type for op in self.operations]
            }
        })
        return base_dict

    def _repr_notebook_html(self, output: 'StandardizedOutput') -> str:
        """
        Render plots in notebook: each plot shown as a side-by-side layout
        with the companion CSV table on the left and the PNG image on the right.
        """
        import base64
        import html as html_module

        try:
            import pandas as pd
        except ImportError:
            return ""

        plots_stream = output.streams.get("plots")
        if not isinstance(plots_stream, DataStream) or len(plots_stream) == 0:
            return ""

        parts = []
        parts.append(
            '<div class="bp-section">'
            '<div class="bp-section-title">plots '
            f'<span style="font-weight: normal; color: #666;">'
            f'({plots_stream.format}, {len(plots_stream)} items)</span></div>'
        )

        for item_id, file_path in plots_stream:
            if not file_path or not os.path.isfile(file_path):
                continue

            # Read PNG as base64
            try:
                with open(file_path, "rb") as f:
                    img_data = base64.b64encode(f.read()).decode("utf-8")
            except Exception:
                continue

            # Find companion CSV (same name, .csv extension)
            csv_path = os.path.splitext(file_path)[0] + ".csv"
            table_html = ""
            if os.path.isfile(csv_path):
                try:
                    df = pd.read_csv(csv_path)
                    if len(df) > 0:
                        table_html = self._render_compact_table(df, html_module)
                except Exception:
                    pass

            # Title from filename
            title = os.path.splitext(os.path.basename(file_path))[0].replace("_", " ")

            # Side-by-side: table left, image right
            parts.append(
                f'<div style="margin: 8px 0;">'
                f'<div style="color: #666; font-size: 0.85em; margin-bottom: 4px;">{html_module.escape(title)}</div>'
                f'<div style="display: flex; align-items: flex-start; gap: 16px; flex-wrap: wrap;">'
            )

            if table_html:
                parts.append(
                    f'<div style="flex: 0 1 auto; max-width: 50%; overflow-x: auto;">'
                    f'{table_html}</div>'
                )

            parts.append(
                f'<div style="flex: 1 1 auto;">'
                f'<img src="data:image/png;base64,{img_data}" '
                f'style="max-width: 100%; height: auto;" /></div>'
            )

            parts.append('</div></div>')

        parts.append('</div>')
        return "\n".join(parts)

    @staticmethod
    def _render_compact_table(df, html_module) -> str:
        """Render a DataFrame as a compact HTML table with 2+...+2 row truncation."""
        import pandas as pd

        columns = list(df.columns)
        rows = []
        rows.append('<table class="bp-table"><tr>')
        for col in columns:
            rows.append(f'<th>{html_module.escape(str(col))}</th>')
        rows.append('</tr>')

        n_rows = len(df)
        if n_rows <= 4:
            display_rows = list(range(n_rows))
            ellipsis_after = None
        else:
            display_rows = list(range(2)) + list(range(n_rows - 2, n_rows))
            ellipsis_after = 2

        for row_i, idx in enumerate(display_rows):
            if ellipsis_after is not None and row_i == ellipsis_after:
                rows.append(
                    f'<tr class="bp-ellipsis"><td colspan="{len(columns)}">'
                    f'... {n_rows - 4} more ...</td></tr>'
                )
            row = df.iloc[idx]
            rows.append('<tr>')
            for col in columns:
                val = str(row[col]) if pd.notna(row[col]) else ''
                display_val = val if len(val) <= 60 else val[:57] + '...'
                rows.append(f'<td>{html_module.escape(display_val)}</td>')
            rows.append('</tr>')

        rows.append('</table>')
        return "\n".join(rows)

    def __str__(self) -> str:
        """String representation."""
        op_types = [op.op_type for op in self.operations]
        return f"Plot({', '.join(op_types)})"
