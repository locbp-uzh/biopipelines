#!/usr/bin/env python3
"""
Plot generation helper script for BioPipelines.

Executes a sequence of plot operations (Scatter, Histogram, Bar, Column)
based on a JSON configuration file. Produces publication-ready PNG figures
using matplotlib and seaborn.
"""

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Any, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Set matplotlib backend for headless rendering
plt.switch_backend('Agg')

# Apply publication-ready style
plt.style.use('seaborn-v0_8-whitegrid')


class PlotBuilder:
    """Builds plots from a sequence of operations."""

    def __init__(self, config: Dict[str, Any]):
        """
        Initialize plot builder.

        Args:
            config: Configuration dict with operations, output_folder, plot_filenames
        """
        self.config = config
        self.output_folder = config.get("output_folder", ".")
        self.plot_filenames = config.get("plot_filenames", [])

        # Track metadata for all plots
        self.metadata: List[Dict[str, Any]] = []

    def _load_table(self, path: str) -> pd.DataFrame:
        """Load a table CSV file."""
        if not os.path.exists(path):
            raise FileNotFoundError(f"Table not found: {path}")
        return pd.read_csv(path)

    def _resolve_data_source(self, data_ref: Dict[str, Any]) -> pd.DataFrame:
        """
        Resolve data reference to DataFrame.

        Args:
            data_ref: Serialized data source from config

        Returns:
            Loaded DataFrame
        """
        path = data_ref.get("path")
        if not path:
            raise ValueError(f"Data source has no path: {data_ref}")
        return self._load_table(path)

    def _get_data_source_name(self, data_ref: Dict[str, Any]) -> str:
        """Get display name for a data source."""
        name = data_ref.get("name")
        if name:
            return name
        path = data_ref.get("path", "")
        return os.path.splitext(os.path.basename(path))[0]

    def _apply_style(self, ax: plt.Axes, title: str = None,
                     xlabel: str = None, ylabel: str = None,
                     x_tick_rotation: float = 0, y_tick_rotation: float = 0,
                     grid: bool = True):
        """Apply consistent styling to axes."""
        if title:
            ax.set_title(title, fontsize=12, fontweight='bold', pad=10)
        if xlabel:
            ax.set_xlabel(xlabel, fontsize=10)
        if ylabel:
            ax.set_ylabel(ylabel, fontsize=10)

        ax.tick_params(axis='both', labelsize=9)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Apply tick rotation
        if x_tick_rotation != 0:
            ha = 'right' if x_tick_rotation > 0 else 'left'
            for label in ax.get_xticklabels():
                label.set_rotation(x_tick_rotation)
                label.set_horizontalalignment(ha)
        if y_tick_rotation != 0:
            for label in ax.get_yticklabels():
                label.set_rotation(y_tick_rotation)

        # Grid control
        if grid:
            ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
        else:
            ax.grid(False)

    def _resolve_label(self, explicit_label: str, name: str, column: str) -> str:
        """Resolve axis label from explicit label, name, or column."""
        if explicit_label:
            return explicit_label
        if name:
            return name
        return column

    def execute_scatter(self, op: Dict[str, Any], output_path: str):
        """
        Execute Scatter operation - create X vs Y scatter plot.

        Args:
            op: Operation dict with data, x, y, color, title, etc.
            output_path: Path to save the plot
        """
        data_ref = op.get("data")
        x_col = op.get("x")
        y_col = op.get("y")
        color_col = op.get("color")
        title = op.get("title")
        xlabel = self._resolve_label(op.get("xlabel"), op.get("x_name"), x_col)
        ylabel = self._resolve_label(op.get("ylabel"), op.get("y_name"), y_col)
        figsize = tuple(op.get("figsize", [8, 6]))
        x_tick_rotation = op.get("x_tick_rotation", 0)
        y_tick_rotation = op.get("y_tick_rotation", 0)
        grid = op.get("grid", True)
        color_legend_loc = op.get("color_legend_loc", "upper right")
        color_legend_outside = op.get("color_legend_outside", False)

        df = self._resolve_data_source(data_ref)

        # Validate columns exist
        if x_col not in df.columns:
            raise ValueError(f"Column '{x_col}' not found in data. Available: {list(df.columns)}")
        if y_col not in df.columns:
            raise ValueError(f"Column '{y_col}' not found in data. Available: {list(df.columns)}")

        fig, ax = plt.subplots(figsize=figsize)

        if color_col and color_col in df.columns:
            # Scatter with color grouping
            groups = df[color_col].unique()
            colors = plt.cm.tab10(np.linspace(0, 1, len(groups)))

            for group, c in zip(groups, colors):
                mask = df[color_col] == group
                ax.scatter(df.loc[mask, x_col], df.loc[mask, y_col],
                          c=[c], label=str(group), alpha=0.7, edgecolors='white', linewidth=0.5)

            if color_legend_outside:
                ax.legend(title=color_col, fontsize=9, title_fontsize=10,
                         loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=True)
            else:
                ax.legend(title=color_col, fontsize=9, title_fontsize=10, loc=color_legend_loc)
        else:
            # Simple scatter
            ax.scatter(df[x_col], df[y_col], c='#1f77b4', alpha=0.7,
                      edgecolors='white', linewidth=0.5)

        self._apply_style(ax, title, xlabel, ylabel,
                         x_tick_rotation=x_tick_rotation, y_tick_rotation=y_tick_rotation,
                         grid=grid)

        plt.tight_layout()
        plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close()

        print(f"  Created scatter plot: {output_path}")

        # Record metadata
        self.metadata.append({
            "filename": os.path.basename(output_path),
            "type": "scatter",
            "title": title or "",
            "x_column": x_col,
            "y_column": y_col,
            "data_sources": self._get_data_source_name(data_ref)
        })

    def execute_histogram(self, op: Dict[str, Any], output_path: str):
        """
        Execute Histogram operation - create distribution histogram.

        Args:
            op: Operation dict with data, x, bins, title, etc.
            output_path: Path to save the plot
        """
        data_ref = op.get("data")
        x_col = op.get("x")
        bins = op.get("bins", 20)
        title = op.get("title")
        xlabel = self._resolve_label(op.get("xlabel"), op.get("x_name"), x_col)
        ylabel = self._resolve_label(op.get("ylabel"), op.get("y_name"), "Count")
        figsize = tuple(op.get("figsize", [8, 6]))

        df = self._resolve_data_source(data_ref)

        # Validate column exists
        if x_col not in df.columns:
            raise ValueError(f"Column '{x_col}' not found in data. Available: {list(df.columns)}")

        fig, ax = plt.subplots(figsize=figsize)

        # Remove NaN values
        values = df[x_col].dropna()

        ax.hist(values, bins=bins, color='#1f77b4', edgecolor='white', alpha=0.8)

        # Add statistics annotation
        mean_val = values.mean()
        std_val = values.std()
        ax.axvline(mean_val, color='#d62728', linestyle='--', linewidth=1.5, label=f'Mean: {mean_val:.2f}')

        stats_text = f'n={len(values)}\nmean={mean_val:.2f}\nstd={std_val:.2f}'
        ax.text(0.95, 0.95, stats_text, transform=ax.transAxes, fontsize=9,
                verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        self._apply_style(ax, title, xlabel, ylabel, grid=True)

        plt.tight_layout()
        plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close()

        print(f"  Created histogram: {output_path}")

        # Record metadata
        self.metadata.append({
            "filename": os.path.basename(output_path),
            "type": "histogram",
            "title": title or "",
            "x_column": x_col,
            "y_column": "",
            "data_sources": self._get_data_source_name(data_ref)
        })

    def execute_bar(self, op: Dict[str, Any], output_path: str):
        """
        Execute Bar operation - create bar chart.

        Args:
            op: Operation dict with data, x, y, title, etc.
            output_path: Path to save the plot
        """
        data_ref = op.get("data")
        x_col = op.get("x")
        y_col = op.get("y")
        title = op.get("title")
        xlabel = self._resolve_label(op.get("xlabel"), op.get("x_name"), x_col)
        ylabel = self._resolve_label(op.get("ylabel"), op.get("y_name"), y_col)
        figsize = tuple(op.get("figsize", [8, 6]))

        df = self._resolve_data_source(data_ref)

        # Validate columns exist
        if x_col not in df.columns:
            raise ValueError(f"Column '{x_col}' not found in data. Available: {list(df.columns)}")
        if y_col not in df.columns:
            raise ValueError(f"Column '{y_col}' not found in data. Available: {list(df.columns)}")

        fig, ax = plt.subplots(figsize=figsize)

        # Aggregate by category if multiple values per category
        if df[x_col].duplicated().any():
            agg_df = df.groupby(x_col)[y_col].mean().reset_index()
        else:
            agg_df = df[[x_col, y_col]]

        x_positions = range(len(agg_df))
        ax.bar(x_positions, agg_df[y_col], color='#1f77b4', edgecolor='white', alpha=0.8)

        ax.set_xticks(x_positions)
        ax.set_xticklabels(agg_df[x_col], rotation=45, ha='right')

        self._apply_style(ax, title, xlabel, ylabel, grid=True)

        plt.tight_layout()
        plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close()

        print(f"  Created bar chart: {output_path}")

        # Record metadata
        self.metadata.append({
            "filename": os.path.basename(output_path),
            "type": "bar",
            "title": title or "",
            "x_column": x_col,
            "y_column": y_col,
            "data_sources": self._get_data_source_name(data_ref)
        })

    def execute_column(self, op: Dict[str, Any], output_path: str):
        """
        Execute Column operation - create column plot with configurable style.

        Supports multiple styles: column (bars+points), simple_bar, scatter, box, floating_bar.

        Args:
            op: Operation dict with data (list), y, labels, style, etc.
            output_path: Path to save the plot
        """
        data_refs = op.get("data", [])
        y_col = op.get("y")
        labels = op.get("labels")
        title = op.get("title")
        xlabel = self._resolve_label(op.get("xlabel"), op.get("x_name"), "")
        ylabel = self._resolve_label(op.get("ylabel"), op.get("y_name"), y_col)
        style = op.get("style", "column")
        show_error = op.get("show_error", "sd")
        figsize = tuple(op.get("figsize", [8, 6]))
        color_groups = op.get("color_groups")
        explicit_colors = op.get("colors")
        color_legend_title = op.get("color_legend_title")
        color_legend_loc = op.get("color_legend_loc", "upper right")
        color_legend_outside = op.get("color_legend_outside", False)
        x_tick_rotation = op.get("x_tick_rotation", 0)
        y_tick_rotation = op.get("y_tick_rotation", 0)
        grid = op.get("grid", True)

        # Load all data sources
        dataframes = []
        source_names = []

        for i, data_ref in enumerate(data_refs):
            df = self._resolve_data_source(data_ref)
            if y_col not in df.columns:
                raise ValueError(f"Column '{y_col}' not found in data source {i}. Available: {list(df.columns)}")
            dataframes.append(df)
            source_names.append(self._get_data_source_name(data_ref))

        # Use provided labels or source names
        if labels and len(labels) == len(dataframes):
            group_labels = labels
        else:
            group_labels = source_names

        fig, ax = plt.subplots(figsize=figsize)

        # Build color mapping based on color_groups
        if color_groups and len(color_groups) == len(dataframes):
            unique_groups = list(dict.fromkeys(color_groups))  # Preserve order, remove duplicates
            if explicit_colors and len(explicit_colors) >= len(unique_groups):
                # Use explicit colors for each unique group
                group_color_map = {g: explicit_colors[i] for i, g in enumerate(unique_groups)}
            else:
                # Generate colors from palette
                palette = plt.cm.tab10(np.linspace(0, 1, len(unique_groups)))
                group_color_map = {g: palette[i] for i, g in enumerate(unique_groups)}
            colors = [group_color_map[g] for g in color_groups]
        elif explicit_colors and len(explicit_colors) >= len(dataframes):
            colors = explicit_colors[:len(dataframes)]
            unique_groups = None
            group_color_map = None
        else:
            # Default color palette
            colors = plt.cm.tab10(np.linspace(0, 1, len(dataframes)))
            unique_groups = None
            group_color_map = None

        x_positions = np.arange(len(dataframes))
        bar_width = 0.6

        # Box plot style - handled separately
        if style == "box":
            box_data = [df[y_col].dropna().values for df in dataframes]
            bp = ax.boxplot(box_data, positions=x_positions, widths=bar_width * 0.8,
                           patch_artist=True, showfliers=True)
            for patch, color in zip(bp['boxes'], colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.6)
            for element in ['whiskers', 'caps', 'medians']:
                plt.setp(bp[element], color='black')
        else:
            # Other styles: column, simple_bar, scatter, floating_bar
            for i, (df, label, color) in enumerate(zip(dataframes, group_labels, colors)):
                values = df[y_col].dropna()
                mean_val = values.mean()
                n = len(values)

                # Calculate error
                if show_error == "sd":
                    error = values.std()
                elif show_error == "sem":
                    error = values.std() / np.sqrt(n) if n > 0 else 0
                elif show_error == "ci":
                    error = 1.96 * values.std() / np.sqrt(n) if n > 0 else 0
                else:
                    error = 0

                if style == "column":
                    # Bars with error bars and overlaid points
                    ax.bar(x_positions[i], mean_val, bar_width, color=color, alpha=0.6,
                          edgecolor='black', linewidth=1)
                    if error > 0:
                        ax.errorbar(x_positions[i], mean_val, yerr=error, fmt='none',
                                  color='black', capsize=5, capthick=1.5, linewidth=1.5)
                    jitter = (np.random.rand(len(values)) - 0.5) * 0.3
                    ax.scatter(x_positions[i] + jitter, values, color=color, alpha=0.8,
                              edgecolors='black', linewidth=0.5, s=30, zorder=5)

                elif style == "simple_bar":
                    # Just bars with error bars
                    ax.bar(x_positions[i], mean_val, bar_width, color=color, alpha=0.8,
                          edgecolor='black', linewidth=1)
                    if error > 0:
                        ax.errorbar(x_positions[i], mean_val, yerr=error, fmt='none',
                                  color='black', capsize=5, capthick=1.5, linewidth=1.5)

                elif style == "scatter":
                    # Just scattered points
                    jitter = (np.random.rand(len(values)) - 0.5) * 0.4
                    ax.scatter(x_positions[i] + jitter, values, color=color, alpha=0.7,
                              edgecolors='black', linewidth=0.5, s=40, zorder=5)
                    # Add mean line
                    ax.hlines(mean_val, x_positions[i] - 0.25, x_positions[i] + 0.25,
                             colors='black', linewidth=2, zorder=6)

                elif style == "floating_bar":
                    # Floating bars showing mean ± error (no baseline)
                    if error > 0:
                        ax.barh(x_positions[i], 2 * error, left=mean_val - error,
                               height=bar_width * 0.5, color=color, alpha=0.6,
                               edgecolor='black', linewidth=1)
                        ax.plot(mean_val, x_positions[i], 'k|', markersize=15, markeredgewidth=2)
                    else:
                        ax.plot(mean_val, x_positions[i], 'o', color=color, markersize=8)

        # For floating_bar, swap axes labels
        if style == "floating_bar":
            ax.set_yticks(x_positions)
            ax.set_yticklabels(group_labels)
            if xlabel:
                ax.set_ylabel(xlabel)
            if ylabel:
                ax.set_xlabel(ylabel)
        else:
            ax.set_xticks(x_positions)
            ax.set_xticklabels(group_labels)

            # Add sample size annotations (not for floating_bar)
            for i, df in enumerate(dataframes):
                n = len(df[y_col].dropna())
                ax.text(x_positions[i], ax.get_ylim()[0] - 0.05 * (ax.get_ylim()[1] - ax.get_ylim()[0]),
                       f'n={n}', ha='center', fontsize=8, color='gray')

        self._apply_style(ax, title,
                         xlabel if style != "floating_bar" else None,
                         ylabel if style != "floating_bar" else None,
                         x_tick_rotation=x_tick_rotation, y_tick_rotation=y_tick_rotation,
                         grid=grid)

        # Add legend for color groups
        if color_groups and group_color_map:
            legend_handles = [mpatches.Patch(color=group_color_map[g], label=g) for g in unique_groups]
            if color_legend_outside:
                ax.legend(handles=legend_handles, title=color_legend_title,
                         loc='center left', bbox_to_anchor=(1.02, 0.5),
                         fontsize=9, title_fontsize=10, frameon=True)
            else:
                ax.legend(handles=legend_handles, title=color_legend_title, loc=color_legend_loc,
                         fontsize=9, title_fontsize=10)

        # Add error bar type annotation if showing error (but no color group legend)
        if show_error and style in ["column", "simple_bar", "floating_bar"] and not color_groups:
            error_label = {"sd": "SD", "sem": "SEM", "ci": "95% CI"}.get(show_error, "")
            if error_label:
                ax.text(0.02, 0.98, f"Error bars: {error_label}",
                       transform=ax.transAxes, fontsize=8, verticalalignment='top')

        plt.tight_layout()
        plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close()

        print(f"  Created column plot ({style}): {output_path}")

        # Record metadata
        self.metadata.append({
            "filename": os.path.basename(output_path),
            "type": f"column ({style})",
            "title": title or "",
            "x_column": "",
            "y_column": y_col,
            "data_sources": "; ".join(source_names)
        })

    def execute_heatmap(self, op: Dict[str, Any], output_path: str):
        """
        Execute HeatMap operation - create heatmap visualization.

        Supports two modes:
        1. Pivot mode: Create heatmap from x, y, value columns
        2. Correlation mode: Create correlation matrix from selected columns

        Args:
            op: Operation dict with data, x, y, value, columns, etc.
            output_path: Path to save the plot
        """
        data_ref = op.get("data")
        x_col = op.get("x")
        y_col = op.get("y")
        value_col = op.get("value")
        columns = op.get("columns")
        title = op.get("title")
        xlabel = self._resolve_label(op.get("xlabel"), op.get("x_name"), x_col or "")
        ylabel = self._resolve_label(op.get("ylabel"), op.get("y_name"), y_col or "")
        cmap = op.get("cmap", "viridis")
        annotate = op.get("annotate", True)
        figsize = tuple(op.get("figsize", [10, 8]))

        df = self._resolve_data_source(data_ref)

        fig, ax = plt.subplots(figsize=figsize)

        if columns:
            # Correlation mode
            # Validate columns exist
            missing_cols = [c for c in columns if c not in df.columns]
            if missing_cols:
                raise ValueError(f"Columns not found: {missing_cols}. Available: {list(df.columns)}")

            # Compute correlation matrix
            corr_matrix = df[columns].corr()
            matrix_data = corr_matrix.values
            x_labels = columns
            y_labels = columns

            if not title:
                title = "Correlation Matrix"

            # Use diverging colormap for correlations
            cmap = "RdBu_r"
            vmin, vmax = -1, 1

        else:
            # Pivot mode
            if x_col not in df.columns:
                raise ValueError(f"Column '{x_col}' not found. Available: {list(df.columns)}")
            if y_col not in df.columns:
                raise ValueError(f"Column '{y_col}' not found. Available: {list(df.columns)}")
            if value_col not in df.columns:
                raise ValueError(f"Column '{value_col}' not found. Available: {list(df.columns)}")

            # Create pivot table
            pivot_df = df.pivot_table(index=y_col, columns=x_col, values=value_col, aggfunc='mean')
            matrix_data = pivot_df.values
            x_labels = pivot_df.columns.tolist()
            y_labels = pivot_df.index.tolist()

            vmin, vmax = None, None

        # Create heatmap
        im = ax.imshow(matrix_data, cmap=cmap, aspect='auto', vmin=vmin, vmax=vmax)

        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        if columns:
            cbar.set_label('Correlation', fontsize=10)
        elif value_col:
            cbar.set_label(value_col, fontsize=10)

        # Set ticks and labels
        ax.set_xticks(np.arange(len(x_labels)))
        ax.set_yticks(np.arange(len(y_labels)))
        ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=9)
        ax.set_yticklabels(y_labels, fontsize=9)

        # Add annotations
        if annotate:
            for i in range(len(y_labels)):
                for j in range(len(x_labels)):
                    value = matrix_data[i, j]
                    if not np.isnan(value):
                        # Choose text color based on background
                        text_color = 'white' if abs(value - (vmin or matrix_data.min())) > (((vmax or matrix_data.max()) - (vmin or matrix_data.min())) / 2) else 'black'
                        ax.text(j, i, f'{value:.2f}', ha='center', va='center',
                               color=text_color, fontsize=8)

        self._apply_style(ax, title, xlabel, ylabel, grid=False)

        plt.tight_layout()
        plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close()

        print(f"  Created heatmap: {output_path}")

        # Record metadata
        mode = "correlation" if columns else "pivot"
        self.metadata.append({
            "filename": os.path.basename(output_path),
            "type": f"heatmap ({mode})",
            "title": title or "",
            "x_column": x_col or ", ".join(columns or []),
            "y_column": y_col or ", ".join(columns or []),
            "data_sources": self._get_data_source_name(data_ref)
        })

    def execute_operation(self, op: Dict[str, Any], output_path: str):
        """
        Execute a single operation.

        Args:
            op: Operation dict with 'op' key indicating operation type
            output_path: Path to save the plot
        """
        op_type = op.get("op")

        if op_type == "scatter":
            self.execute_scatter(op, output_path)
        elif op_type == "histogram":
            self.execute_histogram(op, output_path)
        elif op_type == "bar":
            self.execute_bar(op, output_path)
        elif op_type == "column":
            self.execute_column(op, output_path)
        elif op_type == "heatmap":
            self.execute_heatmap(op, output_path)
        else:
            print(f"Unknown operation type: {op_type}")

    def build_plots(self):
        """Build all plots by executing all operations."""
        operations = self.config.get("operations", [])
        failed = []

        print(f"Building plots with {len(operations)} operations")
        print("=" * 60)

        for i, op in enumerate(operations):
            op_type = op.get("op", "unknown")
            output_filename = self.plot_filenames[i] if i < len(self.plot_filenames) else f"{op_type}_{i+1}.png"
            output_path = os.path.join(self.output_folder, output_filename)

            print(f"\n[{i+1}/{len(operations)}] Creating: {op_type}")
            print("-" * 40)

            try:
                self.execute_operation(op, output_path)
            except Exception as e:
                print(f"  ERROR creating {op_type} plot: {e}")
                failed.append((i + 1, op_type, str(e)))

        # Always write metadata CSV (even if some plots failed)
        self._write_metadata()

        print("\n" + "=" * 60)
        succeeded = len(operations) - len(failed)
        print(f"Plot generation complete: {succeeded}/{len(operations)} plots created")
        if failed:
            print("Failed plots:")
            for num, op_type, err in failed:
                print(f"  [{num}] {op_type}: {err}")
            sys.exit(1)

    def _write_metadata(self):
        """Write plot metadata to CSV file."""
        metadata_path = os.path.join(self.output_folder, "plot_metadata.csv")
        if self.metadata:
            df = pd.DataFrame(self.metadata)
        else:
            df = pd.DataFrame(columns=["filename", "type", "title", "x_column", "y_column", "data_sources"])
        df.to_csv(metadata_path, index=False)
        print(f"\nMetadata written to: {metadata_path}")


def main():
    parser = argparse.ArgumentParser(description='Create plots from pipeline operations')
    parser.add_argument('--config', required=True, help='JSON configuration file')

    args = parser.parse_args()

    # Load configuration
    with open(args.config, 'r') as f:
        config = json.load(f)

    # Create plot builder
    builder = PlotBuilder(config)

    # Build plots
    builder.build_plots()

    print("\nPlot generation completed successfully")


if __name__ == "__main__":
    main()
