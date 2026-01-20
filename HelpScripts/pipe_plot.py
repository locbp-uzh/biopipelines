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
                     xlabel: str = None, ylabel: str = None):
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
        xlabel = op.get("xlabel", x_col)
        ylabel = op.get("ylabel", y_col)
        figsize = tuple(op.get("figsize", [8, 6]))

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

            ax.legend(title=color_col, fontsize=9, title_fontsize=10)
        else:
            # Simple scatter
            ax.scatter(df[x_col], df[y_col], c='#1f77b4', alpha=0.7,
                      edgecolors='white', linewidth=0.5)

        self._apply_style(ax, title, xlabel, ylabel)

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
        xlabel = op.get("xlabel", x_col)
        ylabel = op.get("ylabel", "Count")
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

        self._apply_style(ax, title, xlabel, ylabel)

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
        xlabel = op.get("xlabel", x_col)
        ylabel = op.get("ylabel", y_col)
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

        self._apply_style(ax, title, xlabel, ylabel)

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
        Execute Column operation - create Prism-style column plot.

        Displays grouped columns with error bars and overlaid individual data points,
        commonly used for comparing metrics between conditions or tools.

        Args:
            op: Operation dict with data (list), y, labels, etc.
            output_path: Path to save the plot
        """
        data_refs = op.get("data", [])
        y_col = op.get("y")
        labels = op.get("labels")
        title = op.get("title")
        xlabel = op.get("xlabel", "")
        ylabel = op.get("ylabel", y_col)
        show_points = op.get("show_points", True)
        show_mean = op.get("show_mean", True)
        show_error = op.get("show_error", "sd")
        figsize = tuple(op.get("figsize", [8, 6]))

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

        # Color palette
        colors = plt.cm.tab10(np.linspace(0, 1, len(dataframes)))

        x_positions = np.arange(len(dataframes))
        bar_width = 0.6

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
                # 95% confidence interval
                error = 1.96 * values.std() / np.sqrt(n) if n > 0 else 0
            else:
                error = 0

            # Draw bar with error bars
            if show_mean:
                ax.bar(x_positions[i], mean_val, bar_width, color=color, alpha=0.6,
                      edgecolor='black', linewidth=1)

                if error > 0:
                    ax.errorbar(x_positions[i], mean_val, yerr=error, fmt='none',
                              color='black', capsize=5, capthick=1.5, linewidth=1.5)

            # Overlay individual points with jitter
            if show_points:
                jitter = (np.random.rand(len(values)) - 0.5) * 0.3
                ax.scatter(x_positions[i] + jitter, values, color=color, alpha=0.8,
                          edgecolors='black', linewidth=0.5, s=30, zorder=5)

        ax.set_xticks(x_positions)
        ax.set_xticklabels(group_labels, rotation=0)

        # Add sample size annotations
        for i, df in enumerate(dataframes):
            n = len(df[y_col].dropna())
            ax.text(x_positions[i], ax.get_ylim()[0] - 0.05 * (ax.get_ylim()[1] - ax.get_ylim()[0]),
                   f'n={n}', ha='center', fontsize=8, color='gray')

        self._apply_style(ax, title, xlabel, ylabel)

        # Add error bar type to legend if showing error
        if show_error and show_mean:
            error_label = {"sd": "SD", "sem": "SEM", "ci": "95% CI"}.get(show_error, "")
            if error_label:
                ax.text(0.02, 0.98, f"Error bars: {error_label}",
                       transform=ax.transAxes, fontsize=8, verticalalignment='top')

        plt.tight_layout()
        plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close()

        print(f"  Created column plot: {output_path}")

        # Record metadata
        self.metadata.append({
            "filename": os.path.basename(output_path),
            "type": "column",
            "title": title or "",
            "x_column": "",
            "y_column": y_col,
            "data_sources": "; ".join(source_names)
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
        else:
            print(f"Unknown operation type: {op_type}")

    def build_plots(self):
        """Build all plots by executing all operations."""
        operations = self.config.get("operations", [])

        print(f"Building plots with {len(operations)} operations")
        print("=" * 60)

        for i, op in enumerate(operations):
            op_type = op.get("op", "unknown")
            output_filename = self.plot_filenames[i] if i < len(self.plot_filenames) else f"{op_type}_{i+1}.png"
            output_path = os.path.join(self.output_folder, output_filename)

            print(f"\n[{i+1}/{len(operations)}] Creating: {op_type}")
            print("-" * 40)

            self.execute_operation(op, output_path)

        # Write metadata CSV
        self._write_metadata()

        print("\n" + "=" * 60)
        print(f"Plot generation complete: {len(operations)} plots created")

    def _write_metadata(self):
        """Write plot metadata to CSV file."""
        if not self.metadata:
            return

        metadata_path = os.path.join(self.output_folder, "plot_metadata.csv")
        df = pd.DataFrame(self.metadata)
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
