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
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo


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
    publication-ready plots from CSV data. Supports multiple data sources for
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
    DEFAULT_ENV = None

    # --- Static methods for creating operations ---

    @staticmethod
    def Scatter(data: Union[StandardizedOutput, ToolOutput, TableInfo],
                x: str,
                y: str,
                color: str = None,
                title: str = None,
                xlabel: str = None,
                ylabel: str = None,
                x_name: str = None,
                y_name: str = None,
                figsize: Tuple[float, float] = (8, 6),
                x_tick_rotation: float = 0,
                y_tick_rotation: float = 0,
                grid: bool = True,
                legend_loc: str = "upper right") -> PlotOperation:
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
            figsize: Figure size as (width, height) in inches
            x_tick_rotation: Rotation angle for x-axis tick labels in degrees
            y_tick_rotation: Rotation angle for y-axis tick labels in degrees
            grid: Show grid lines (default True)
            legend_loc: Legend location ("upper right", "upper left", etc.)

        Returns:
            PlotOperation for scatter plot
        """
        return PlotOperation("scatter", data=data, x=x, y=y, color=color,
                            title=title, xlabel=xlabel, ylabel=ylabel,
                            x_name=x_name, y_name=y_name, figsize=figsize,
                            x_tick_rotation=x_tick_rotation, y_tick_rotation=y_tick_rotation,
                            grid=grid, legend_loc=legend_loc)

    @staticmethod
    def Histogram(data: Union[StandardizedOutput, ToolOutput, TableInfo],
                  x: str,
                  bins: int = 20,
                  title: str = None,
                  xlabel: str = None,
                  ylabel: str = None,
                  x_name: str = None,
                  y_name: str = None,
                  figsize: Tuple[float, float] = (8, 6)) -> PlotOperation:
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
            figsize: Figure size as (width, height) in inches

        Returns:
            PlotOperation for histogram
        """
        return PlotOperation("histogram", data=data, x=x, bins=bins,
                            title=title, xlabel=xlabel, ylabel=ylabel,
                            x_name=x_name, y_name=y_name, figsize=figsize)

    @staticmethod
    def Bar(data: Union[StandardizedOutput, ToolOutput, TableInfo],
            x: str,
            y: str,
            title: str = None,
            xlabel: str = None,
            ylabel: str = None,
            x_name: str = None,
            y_name: str = None,
            figsize: Tuple[float, float] = (8, 6)) -> PlotOperation:
        """
        Create bar chart for categorical X axis.

        Args:
            data: Data source (tool output or table)
            x: Column name for categories (X axis)
            y: Column name for values (Y axis)
            title: Plot title (optional)
            xlabel: X axis label (defaults to x_name or column name)
            ylabel: Y axis label (defaults to y_name or column name)
            x_name: Display name for X variable (alternative to xlabel)
            y_name: Display name for Y variable (alternative to ylabel)
            figsize: Figure size as (width, height) in inches

        Returns:
            PlotOperation for bar chart
        """
        return PlotOperation("bar", data=data, x=x, y=y,
                            title=title, xlabel=xlabel, ylabel=ylabel,
                            x_name=x_name, y_name=y_name, figsize=figsize)

    @staticmethod
    def Column(data: List[Union[StandardizedOutput, ToolOutput, TableInfo]],
               y: str,
               labels: List[str] = None,
               title: str = None,
               xlabel: str = None,
               ylabel: str = None,
               x_name: str = None,
               y_name: str = None,
               show_points: bool = True,
               show_mean: bool = True,
               show_error: str = "sd",
               figsize: Tuple[float, float] = (8, 6),
               color_groups: List[str] = None,
               colors: List[str] = None,
               legend_title: str = None,
               legend_loc: str = "upper right",
               x_tick_rotation: float = 0,
               y_tick_rotation: float = 0,
               grid: bool = True) -> PlotOperation:
        """
        Create Prism-style column plot comparing same metric across multiple data sources.

        Displays grouped columns with error bars and optionally overlaid individual
        data points, commonly used for comparing metrics between conditions or tools.

        Args:
            data: List of data sources [tool1.tables.x, tool2.tables.x, ...]
            y: Column name to compare across sources
            labels: Group labels for x-axis (defaults to source filenames)
            title: Plot title (optional)
            xlabel: X axis label (defaults to x_name)
            ylabel: Y axis label (defaults to y_name or column name)
            x_name: Display name for X axis (alternative to xlabel)
            y_name: Display name for Y variable (alternative to ylabel)
            show_points: Overlay individual data points on bars
            show_mean: Show mean line/bar
            show_error: Error bar type: "sd" (std dev), "sem" (std error), "ci" (95% CI), or None
            figsize: Figure size as (width, height) in inches
            color_groups: List of group names for coloring (e.g., protein names). Same color for same group.
            colors: Explicit list of colors (hex or named). If not provided, uses a color palette.
            legend_title: Title for the color legend (e.g., "Protein")
            legend_loc: Legend location ("upper right", "upper left", "lower right", "lower left", etc.)
            x_tick_rotation: Rotation angle for x-axis tick labels in degrees
            y_tick_rotation: Rotation angle for y-axis tick labels in degrees
            grid: Show grid lines (default True)

        Returns:
            PlotOperation for column plot
        """
        return PlotOperation("column", data=data, y=y, labels=labels,
                            title=title, xlabel=xlabel, ylabel=ylabel,
                            x_name=x_name, y_name=y_name,
                            show_points=show_points, show_mean=show_mean,
                            show_error=show_error, figsize=figsize,
                            color_groups=color_groups, colors=colors,
                            legend_title=legend_title, legend_loc=legend_loc,
                            x_tick_rotation=x_tick_rotation, y_tick_rotation=y_tick_rotation,
                            grid=grid)

    @staticmethod
    def HeatMap(data: Union[StandardizedOutput, ToolOutput, TableInfo],
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
                figsize: Tuple[float, float] = (10, 8)) -> PlotOperation:
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
            figsize: Figure size as (width, height) in inches

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
                            annotate=annotate, figsize=figsize)

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
            if hasattr(source, 'path'):
                # TableInfo
                self.input_sources[f"data_{i}"] = {
                    'path': source.path,
                    'name': getattr(source, 'name', f'table_{i}')
                }
            elif hasattr(source, 'output_folder'):
                # StandardizedOutput or ToolOutput
                self.input_sources[f"data_{i}"] = {
                    'output_folder': source.output_folder
                }

    def _get_table_path(self, data) -> str:
        """Extract table path from various data source types."""
        if isinstance(data, TableInfo):
            return data.path
        elif hasattr(data, 'tables'):
            tables = data.tables
            if hasattr(tables, '_tables'):
                # Get first table
                first_name, info = next(iter(tables._tables.items()))
                return info.path
            elif isinstance(tables, dict):
                first_name, info = next(iter(tables.items()))
                if isinstance(info, dict) and 'path' in info:
                    return info['path']
                elif hasattr(info, 'path'):
                    return info.path
                return str(info)
        raise ValueError(f"Cannot extract table path from: {type(data)}")

    def _get_table_name(self, data) -> str:
        """Extract table name from various data source types."""
        if isinstance(data, TableInfo):
            return data.name
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
                "path": data.path,
                "name": data.name
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
        """Generate bash script for plot creation."""
        # Create config file with all operations
        config = {
            "operations": [self._serialize_operation(op) for op in self.operations],
            "output_folder": self.output_folder,
            "plot_filenames": [self._generate_plot_filename(op, i) for i, op in enumerate(self.operations)]
        }

        config_file = os.path.join(self.output_folder, "plot_config.json")

        # Get HelpScripts path
        help_scripts = self.folders.get("HelpScripts", "HelpScripts")
        plot_script = os.path.join(help_scripts, "pipe_plot.py")

        script_content = f"""#!/bin/bash
# Plot generation script
# Generated by BioPipelines

{self.generate_completion_check_header()}

echo "Creating plots..."
echo "Output folder: {self.output_folder}"

# Create output directory
mkdir -p "{self.output_folder}"

# Write configuration
cat > "{config_file}" << 'PLOT_CONFIG_EOF'
{json.dumps(config, indent=2)}
PLOT_CONFIG_EOF

# Run plot generation
python "{plot_script}" --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Plots created successfully"
else
    echo "ERROR: Failed to create plots"
    exit 1
fi

{self.generate_completion_check_footer()}
"""

        return script_content

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files."""
        # Generate plot filenames
        plot_files = []
        for i, op in enumerate(self.operations):
            filename = self._generate_plot_filename(op, i)
            plot_files.append(os.path.join(self.output_folder, filename))

        # Metadata table
        metadata_csv = os.path.join(self.output_folder, "plot_metadata.csv")

        # Build columns for metadata table
        metadata_columns = ["filename", "type", "title", "x_column", "y_column", "data_sources"]

        tables = {
            "metadata": TableInfo(
                name="metadata",
                path=metadata_csv,
                columns=metadata_columns,
                description="Summary of all generated plots",
                count=len(self.operations)
            )
        }

        return {
            "structures": [],
            "structure_ids": [],
            "compounds": [],
            "compound_ids": [],
            "sequences": [],
            "sequence_ids": [],
            "tables": tables,
            "output_folder": self.output_folder,
            "plots": plot_files
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

    def __str__(self) -> str:
        """String representation."""
        op_types = [op.op_type for op in self.operations]
        return f"Plot({', '.join(op_types)})"
