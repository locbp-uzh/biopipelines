#!/usr/bin/env python3
"""
PyMOL session creation helper script for BioPipelines.

Executes a sequence of PyMOL operations (Names, Load, Color, Align, etc.)
based on a JSON configuration file. Supports ID-based matching between
structures and table columns for per-structure selections and naming.
"""

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Any, Optional

import pandas as pd
import pymol
from pymol import cmd


class PyMOLSessionBuilder:
    """Builds PyMOL sessions from a sequence of operations."""

    def __init__(self, config: Dict[str, Any]):
        """
        Initialize session builder.

        Args:
            config: Configuration dict with operations, session_name, output_folder
        """
        self.config = config
        self.output_folder = config.get("output_folder", ".")
        self.session_name = config.get("session_name", "session")

        # Current naming map: id -> pymol_name
        self.naming_map: Dict[str, str] = {}

        # Track loaded objects: id -> pymol_name
        self.loaded_objects: Dict[str, str] = {}

        # Track first loaded object for alignment
        self.first_loaded_object: Optional[str] = None

    def setup_pymol(self):
        """Initialize PyMOL in headless mode with standard settings."""
        pymol.pymol_argv = ['pymol', '-c']  # -c for no GUI
        pymol.finish_launching()

        # Standard visualization settings
        cmd.set("seq_view", 1)
        cmd.set("cartoon_gap_cutoff", 20)
        cmd.set("cartoon_sampling", 20)
        cmd.set("sphere_scale", 0.2)
        cmd.set("ray_trace_mode", 1)
        cmd.set("ray_shadows", 0)
        cmd.set("spec_reflect", 0)
        cmd.set("ray_trace_frames", 1)
        cmd.set("ray_trace_color", "gray20")

        # Set up pLDDT colors
        self._setup_plddt_colors()

    def _setup_plddt_colors(self):
        """Set up standard pLDDT color scheme using hex colors."""
        # Colors match AlphaFold confidence visualization
        # Blue (high confidence), Cyan, Yellow, Orange (low confidence)
        pass  # Colors are applied directly via hex in coloraf

    def _load_table(self, table_path: str) -> pd.DataFrame:
        """Load a table CSV file."""
        if not os.path.exists(table_path):
            raise FileNotFoundError(f"Table not found: {table_path}")
        return pd.read_csv(table_path)

    def _get_pymol_name(self, structure_id: str) -> str:
        """Get PyMOL object name for a structure ID using current naming map."""
        if structure_id in self.naming_map:
            return self.naming_map[structure_id]
        return structure_id

    def _resolve_structures(self, structures_ref: Dict[str, Any]) -> List[Dict[str, str]]:
        """
        Resolve structure reference to list of {id, path} dicts.

        Args:
            structures_ref: Serialized structure reference from config

        Returns:
            List of dicts with 'id' and 'path' keys
        """
        ref_type = structures_ref.get("type", "")
        structures = structures_ref.get("structures", [])
        structure_ids = structures_ref.get("structure_ids", [])

        result = []
        for i, path in enumerate(structures):
            if i < len(structure_ids):
                struct_id = structure_ids[i]
            else:
                # Extract ID from filename
                struct_id = os.path.splitext(os.path.basename(path))[0]
            result.append({"id": struct_id, "path": path})

        return result

    def _resolve_table_column(self, column_ref: Dict[str, Any]) -> Dict[str, str]:
        """
        Resolve table column reference to id -> value mapping.

        Args:
            column_ref: Dict with table_path and column_name

        Returns:
            Dict mapping id to column value
        """
        table_path = column_ref.get("table_path")
        column_name = column_ref.get("column_name")

        if not table_path or not column_name:
            raise ValueError(f"Invalid table column reference: {column_ref}")

        df = self._load_table(table_path)

        if "id" not in df.columns:
            raise ValueError(f"Table {table_path} has no 'id' column")
        if column_name not in df.columns:
            raise ValueError(f"Column '{column_name}' not found in {table_path}")

        return dict(zip(df["id"].astype(str), df[column_name].astype(str)))

    def execute_names(self, op: Dict[str, Any]):
        """
        Execute Names operation - set up ID -> PyMOL name mapping.

        Args:
            op: Operation dict with prefix, basename, suffix
        """
        prefix = op.get("prefix", "")
        suffix = op.get("suffix", "")
        basename_ref = op.get("basename")

        # Clear current naming map
        self.naming_map = {}

        if basename_ref is None:
            print("Names: No basename provided, using IDs directly")
            return

        # Resolve basename to id -> value mapping
        if isinstance(basename_ref, dict) and basename_ref.get("type") == "table_column":
            id_to_basename = self._resolve_table_column(basename_ref)

            for struct_id, basename_value in id_to_basename.items():
                # Build PyMOL name: prefix_basename_suffix
                parts = []
                if prefix:
                    parts.append(prefix)
                parts.append(str(basename_value))
                if suffix:
                    parts.append(suffix)

                pymol_name = "_".join(parts)
                self.naming_map[struct_id] = pymol_name

            print(f"Names: Set up mapping for {len(self.naming_map)} structures")
            print(f"  Format: {prefix}_<{basename_ref.get('column_name')}>_{suffix}")
        else:
            print(f"Names: Unknown basename format: {basename_ref}")

    def execute_load(self, op: Dict[str, Any]):
        """
        Execute Load operation - load structures into PyMOL.

        Args:
            op: Operation dict with structures reference
        """
        structures_ref = op.get("structures")
        if not structures_ref:
            raise ValueError("Load operation requires structures")

        structures = self._resolve_structures(structures_ref)
        print(f"Load: Loading {len(structures)} structures")

        for struct in structures:
            struct_id = struct["id"]
            struct_path = struct["path"]

            if not os.path.exists(struct_path):
                print(f"  Warning: Structure file not found: {struct_path}")
                continue

            # Get PyMOL name from naming map
            pymol_name = self._get_pymol_name(struct_id)

            try:
                cmd.load(struct_path, pymol_name)
                self.loaded_objects[struct_id] = pymol_name

                if self.first_loaded_object is None:
                    self.first_loaded_object = pymol_name

                print(f"  Loaded: {pymol_name} (id: {struct_id})")
            except Exception as e:
                print(f"  Error loading {struct_path}: {e}")

    def execute_color(self, op: Dict[str, Any]):
        """
        Execute Color operation - color structures by selection.

        Args:
            op: Operation dict with structures, selection, color
        """
        structures_ref = op.get("structures")
        selection_ref = op.get("selection")
        color = op.get("color", "white")

        if not structures_ref:
            raise ValueError("Color operation requires structures")

        structures = self._resolve_structures(structures_ref)

        # Resolve selection - can be a fixed string or table column reference
        if isinstance(selection_ref, dict) and selection_ref.get("type") == "table_column":
            id_to_selection = self._resolve_table_column(selection_ref)
            print(f"Color: Applying color '{color}' with per-structure selections")
        elif isinstance(selection_ref, str):
            # Fixed selection for all structures
            id_to_selection = {s["id"]: selection_ref for s in structures}
            print(f"Color: Applying color '{color}' with fixed selection '{selection_ref}'")
        else:
            raise ValueError(f"Invalid selection format: {selection_ref}")

        for struct in structures:
            struct_id = struct["id"]

            if struct_id not in self.loaded_objects:
                print(f"  Skipping {struct_id} - not loaded")
                continue

            pymol_name = self.loaded_objects[struct_id]

            if struct_id not in id_to_selection:
                print(f"  Skipping {struct_id} - no selection value")
                continue

            selection_value = id_to_selection[struct_id]

            # Build PyMOL selection: object and resi selection
            # Selection value is in PyMOL format like "1-40" or "41" or "42-78"
            pymol_selection = f"{pymol_name} and resi {selection_value}"

            try:
                cmd.color(color, pymol_selection)
                print(f"  Colored {pymol_name} resi {selection_value} -> {color}")
            except Exception as e:
                print(f"  Error coloring {pymol_name}: {e}")

    def execute_coloraf(self, op: Dict[str, Any]):
        """
        Execute ColorAF operation - color structures by pLDDT (B-factor).

        Args:
            op: Operation dict with structures reference and optional upper value
        """
        structures_ref = op.get("structures")
        upper = op.get("upper", 100)

        if not structures_ref:
            raise ValueError("ColorAF operation requires structures")

        structures = self._resolve_structures(structures_ref)
        print(f"ColorAF: Applying pLDDT coloring to {len(structures)} structures (upper={upper})")

        # Calculate thresholds based on upper value
        if upper == 100:
            extremes = [int(float(upper) * x) for x in [0.5, 0.7, 0.9]]
        else:
            extremes = [float(upper) * x for x in [0.5, 0.7, 0.9]]

        for struct in structures:
            struct_id = struct["id"]

            if struct_id not in self.loaded_objects:
                print(f"  Skipping {struct_id} - not loaded")
                continue

            pymol_name = self.loaded_objects[struct_id]

            try:
                # Color by pLDDT ranges using B-factor with hex colors
                # Blue (high confidence >= 90%)
                cmd.color("0x126DFF", f"({pymol_name}) and (b > {extremes[2]} or b = {extremes[2]})")
                # Cyan (70-90%)
                cmd.color("0x0ECFF1", f"({pymol_name}) and ((b < {extremes[2]} and b > {extremes[1]}) or b = {extremes[1]})")
                # Yellow (50-70%)
                cmd.color("0xF6ED12", f"({pymol_name}) and ((b < {extremes[1]} and b > {extremes[0]}) or b = {extremes[0]})")
                # Orange (low confidence < 50%)
                cmd.color("0xEE831D", f"({pymol_name}) and b < {extremes[0]}")
                print(f"  Applied pLDDT coloring to {pymol_name}")
            except Exception as e:
                print(f"  Error coloring {pymol_name}: {e}")

    def execute_align(self, op: Dict[str, Any]):
        """
        Execute Align operation - align all loaded objects.

        Args:
            op: Operation dict with method and optional target
        """
        method = op.get("method", "align")
        target = op.get("target")

        if len(self.loaded_objects) <= 1:
            print("Align: Only one structure loaded, skipping alignment")
            return

        # Determine target for alignment
        if target:
            target_name = target
        elif self.first_loaded_object:
            target_name = self.first_loaded_object
        else:
            target_name = list(self.loaded_objects.values())[0]

        print(f"Align: Aligning to {target_name} using method '{method}'")

        for struct_id, pymol_name in self.loaded_objects.items():
            if pymol_name == target_name:
                continue

            try:
                if method == "align":
                    result = cmd.align(pymol_name, target_name)
                    rmsd = result[0] if isinstance(result, tuple) else 0
                    print(f"  Aligned {pymol_name}: RMSD = {rmsd:.2f}")
                elif method == "super":
                    result = cmd.super(pymol_name, target_name)
                    rmsd = result[0] if isinstance(result, tuple) else 0
                    print(f"  Super-aligned {pymol_name}: RMSD = {rmsd:.2f}")
                elif method == "cealign":
                    result = cmd.cealign(target_name, pymol_name)
                    rmsd = result.get("RMSD", 0) if isinstance(result, dict) else 0
                    print(f"  CE-aligned {pymol_name}: RMSD = {rmsd:.2f}")
                else:
                    print(f"  Unknown alignment method: {method}")
            except Exception as e:
                print(f"  Error aligning {pymol_name}: {e}")

    def execute_show(self, op: Dict[str, Any]):
        """
        Execute Show operation - show a representation.

        Args:
            op: Operation dict with representation and optional selection
        """
        representation = op.get("representation", "cartoon")
        selection = op.get("selection", "all")

        try:
            cmd.show(representation, selection)
            print(f"Show: {representation} for {selection}")
        except Exception as e:
            print(f"Error in show: {e}")

    def execute_hide(self, op: Dict[str, Any]):
        """
        Execute Hide operation - hide a representation.

        Args:
            op: Operation dict with representation and optional selection
        """
        representation = op.get("representation", "everything")
        selection = op.get("selection", "all")

        try:
            cmd.hide(representation, selection)
            print(f"Hide: {representation} for {selection}")
        except Exception as e:
            print(f"Error in hide: {e}")

    def execute_set(self, op: Dict[str, Any]):
        """
        Execute Set operation - set a PyMOL setting.

        Args:
            op: Operation dict with setting, value, and optional selection
        """
        setting = op.get("setting")
        value = op.get("value")
        selection = op.get("selection")

        try:
            if selection:
                cmd.set(setting, value, selection)
            else:
                cmd.set(setting, value)
            print(f"Set: {setting} = {value}")
        except Exception as e:
            print(f"Error in set: {e}")

    def execute_save(self, op: Dict[str, Any]):
        """
        Execute Save operation - save the session.

        Args:
            op: Operation dict with filename
        """
        filename = op.get("filename", f"{self.session_name}.pse")

        # Make path absolute if relative
        if not os.path.isabs(filename):
            filename = os.path.join(self.output_folder, filename)

        try:
            cmd.save(filename)
            print(f"Save: Saved session to {filename}")
        except Exception as e:
            print(f"Error saving session: {e}")

    def execute_operation(self, op: Dict[str, Any]):
        """
        Execute a single operation.

        Args:
            op: Operation dict with 'op' key indicating operation type
        """
        op_type = op.get("op")

        if op_type == "names":
            self.execute_names(op)
        elif op_type == "load":
            self.execute_load(op)
        elif op_type == "color":
            self.execute_color(op)
        elif op_type == "coloraf":
            self.execute_coloraf(op)
        elif op_type == "align":
            self.execute_align(op)
        elif op_type == "show":
            self.execute_show(op)
        elif op_type == "hide":
            self.execute_hide(op)
        elif op_type == "set":
            self.execute_set(op)
        elif op_type == "save":
            self.execute_save(op)
        else:
            print(f"Unknown operation type: {op_type}")

    def build_session(self):
        """Build the PyMOL session by executing all operations."""
        operations = self.config.get("operations", [])

        print(f"Building PyMOL session with {len(operations)} operations")
        print("=" * 60)

        for i, op in enumerate(operations):
            op_type = op.get("op", "unknown")
            print(f"\n[{i+1}/{len(operations)}] Executing: {op_type}")
            print("-" * 40)
            self.execute_operation(op)

        # Auto-save if no explicit Save operation
        has_save = any(op.get("op") == "save" for op in operations)
        if not has_save:
            print("\n" + "=" * 60)
            print("Auto-saving session...")
            session_path = os.path.join(self.output_folder, f"{self.session_name}.pse")
            cmd.save(session_path)
            print(f"Saved: {session_path}")

        print("\n" + "=" * 60)
        print(f"Session complete: {len(self.loaded_objects)} structures loaded")


def main():
    parser = argparse.ArgumentParser(description='Create PyMOL session from pipeline operations')
    parser.add_argument('--config', required=True, help='JSON configuration file')

    args = parser.parse_args()

    # Load configuration
    with open(args.config, 'r') as f:
        config = json.load(f)

    # Create session builder
    builder = PyMOLSessionBuilder(config)

    # Initialize PyMOL
    builder.setup_pymol()

    # Build session
    builder.build_session()

    # Quit PyMOL
    cmd.quit()
    print("\nPyMOL session creation completed successfully")


if __name__ == "__main__":
    main()
