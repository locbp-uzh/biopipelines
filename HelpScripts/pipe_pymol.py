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

# Set environment for headless PyMOL before importing
# This avoids libGL.so.1 errors on headless cluster nodes
os.environ.setdefault('PYMOL_SYMMETRY_VIEWER', '0')
os.environ.setdefault('DISPLAY', '')

import pandas as pd
import pymol
pymol.pymol_argv = ['pymol', '-cq']  # -c: no GUI, -q: quiet
pymol.finish_launching()
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
        """Initialize PyMOL with standard settings (already launched at import time)."""
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

        Note: When LoadOutput is used with validate_files=True (default), glob patterns
        are already resolved at pipeline runtime. This method handles:
        1. Standard: len(structures) == len(structure_ids), one file per ID
        2. Legacy glob: structures contains glob pattern(s) - fallback for validate_files=False
        3. Single file: one structure path with multiple IDs (expand by searching directory)

        Args:
            structures_ref: Serialized structure reference from config

        Returns:
            List of dicts with 'id' and 'path' keys
        """
        import glob as glob_module

        ref_type = structures_ref.get("type", "")
        structures = structures_ref.get("structures", [])
        structure_ids = structures_ref.get("structure_ids", [])

        result = []

        # Case 1: Standard format - one file per ID
        if len(structures) == len(structure_ids) and len(structures) > 0:
            for i, path in enumerate(structures):
                struct_id = structure_ids[i]
                result.append({"id": struct_id, "path": path})
            return result

        # Case 2: Legacy format - glob pattern(s) or fewer files than IDs
        # Need to match structure_ids to actual files in the directory
        if len(structures) > 0 and len(structure_ids) > 0:
            # Get directory from first structure path
            first_path = structures[0]

            # Check if it's a glob pattern (contains * or ?)
            if '*' in first_path or '?' in first_path:
                # Expand glob pattern
                expanded_files = sorted(glob_module.glob(first_path))
                print(f"  Expanded glob pattern: {len(expanded_files)} files found")

                # Try to match IDs to expanded files
                for struct_id in structure_ids:
                    matched_path = None
                    for file_path in expanded_files:
                        basename = os.path.splitext(os.path.basename(file_path))[0]
                        # Match if ID is in filename or filename contains ID
                        if struct_id in basename or basename in struct_id:
                            matched_path = file_path
                            break
                        # Also try matching rank numbers (e.g., rank0001 -> rankNNNN_something)
                        if struct_id.startswith("rank") and basename.startswith("rank"):
                            # Extract rank number from both
                            id_num = ''.join(filter(str.isdigit, struct_id))
                            file_num = ''.join(filter(str.isdigit, basename.split('_')[0]))
                            if id_num == file_num:
                                matched_path = file_path
                                break

                    if matched_path:
                        result.append({"id": struct_id, "path": matched_path})
                    else:
                        print(f"  Warning: No file found for ID '{struct_id}'")

            else:
                # Not a glob - try to find files in the same directory
                directory = os.path.dirname(first_path)
                if os.path.isdir(directory):
                    # List all structure files in directory
                    extensions = ['.pdb', '.cif', '.mmcif']
                    all_files = []
                    for ext in extensions:
                        all_files.extend(glob_module.glob(os.path.join(directory, f"*{ext}")))

                    # Match IDs to files
                    for struct_id in structure_ids:
                        matched_path = None
                        for file_path in all_files:
                            basename = os.path.splitext(os.path.basename(file_path))[0]
                            if struct_id in basename or basename == struct_id:
                                matched_path = file_path
                                break

                        if matched_path:
                            result.append({"id": struct_id, "path": matched_path})
                        else:
                            print(f"  Warning: No file found for ID '{struct_id}'")
                else:
                    # Just use what we have
                    for i, path in enumerate(structures):
                        if i < len(structure_ids):
                            result.append({"id": structure_ids[i], "path": path})

        # Case 3: Only structure_ids, no paths - try to infer from output_folder
        elif len(structure_ids) > 0 and len(structures) == 0:
            output_folder = structures_ref.get("output_folder", "")
            if output_folder and os.path.isdir(output_folder):
                extensions = ['.pdb', '.cif', '.mmcif']
                for struct_id in structure_ids:
                    for ext in extensions:
                        potential_path = os.path.join(output_folder, f"{struct_id}{ext}")
                        if os.path.exists(potential_path):
                            result.append({"id": struct_id, "path": potential_path})
                            break

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
            op: Operation dict with structures, representation, and optional selection
        """
        structures_ref = op.get("structures")
        representation = op.get("representation", "cartoon")
        selection_ref = op.get("selection")

        # If no structures specified, apply to all
        if not structures_ref:
            selection = selection_ref if selection_ref else "all"
            try:
                cmd.show(representation, selection)
                print(f"Show: {representation} for {selection}")
            except Exception as e:
                print(f"Error in show: {e}")
            return

        structures = self._resolve_structures(structures_ref)

        # Resolve selection if it's a table column reference
        if isinstance(selection_ref, dict) and selection_ref.get("type") == "table_column":
            id_to_selection = self._resolve_table_column(selection_ref)
        elif isinstance(selection_ref, str):
            id_to_selection = {s["id"]: selection_ref for s in structures}
        else:
            id_to_selection = None

        print(f"Show: {representation} for {len(structures)} structures")

        for struct in structures:
            struct_id = struct["id"]

            if struct_id not in self.loaded_objects:
                print(f"  Skipping {struct_id} - not loaded")
                continue

            pymol_name = self.loaded_objects[struct_id]

            if id_to_selection and struct_id in id_to_selection:
                pymol_selection = f"{pymol_name} and resi {id_to_selection[struct_id]}"
            else:
                pymol_selection = pymol_name

            try:
                cmd.show(representation, pymol_selection)
                print(f"  Show {representation} for {pymol_selection}")
            except Exception as e:
                print(f"  Error showing {pymol_name}: {e}")

    def execute_hide(self, op: Dict[str, Any]):
        """
        Execute Hide operation - hide a representation.

        Args:
            op: Operation dict with structures, representation, and optional selection
        """
        structures_ref = op.get("structures")
        representation = op.get("representation", "everything")
        selection_ref = op.get("selection")

        # If no structures specified, apply to all
        if not structures_ref:
            selection = selection_ref if selection_ref else "all"
            try:
                cmd.hide(representation, selection)
                print(f"Hide: {representation} for {selection}")
            except Exception as e:
                print(f"Error in hide: {e}")
            return

        structures = self._resolve_structures(structures_ref)

        # Resolve selection if it's a table column reference
        if isinstance(selection_ref, dict) and selection_ref.get("type") == "table_column":
            id_to_selection = self._resolve_table_column(selection_ref)
        elif isinstance(selection_ref, str):
            id_to_selection = {s["id"]: selection_ref for s in structures}
        else:
            id_to_selection = None

        print(f"Hide: {representation} for {len(structures)} structures")

        for struct in structures:
            struct_id = struct["id"]

            if struct_id not in self.loaded_objects:
                print(f"  Skipping {struct_id} - not loaded")
                continue

            pymol_name = self.loaded_objects[struct_id]

            if id_to_selection and struct_id in id_to_selection:
                pymol_selection = f"{pymol_name} and resi {id_to_selection[struct_id]}"
            else:
                pymol_selection = pymol_name

            try:
                cmd.hide(representation, pymol_selection)
                print(f"  Hide {representation} for {pymol_selection}")
            except Exception as e:
                print(f"  Error hiding {pymol_name}: {e}")

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

    def execute_orient(self, op: Dict[str, Any]):
        """
        Execute Orient operation - orient view to a selection.

        Args:
            op: Operation dict with selection and zoom flag
        """
        selection = op.get("selection", "all")
        zoom = op.get("zoom", True)

        try:
            cmd.orient(selection)
            if zoom:
                cmd.zoom(selection, buffer=5)
            print(f"Orient: Oriented view to '{selection}'")
        except Exception as e:
            print(f"Error orienting: {e}")

    def execute_ray(self, op: Dict[str, Any]):
        """
        Execute Ray operation - ray trace the current view.

        Args:
            op: Operation dict with width and height
        """
        width = op.get("width", 1920)
        height = op.get("height", 1080)

        try:
            cmd.ray(width, height)
            print(f"Ray: Ray traced at {width}x{height}")
        except Exception as e:
            print(f"Error ray tracing: {e}")

    def execute_png(self, op: Dict[str, Any]):
        """
        Execute PNG operation - save current view as PNG.

        Args:
            op: Operation dict with filename, dimensions, ray flag, and dpi
        """
        filename = op.get("filename", "render.png")
        width = op.get("width", 1920)
        height = op.get("height", 1080)
        ray = op.get("ray", True)
        dpi = op.get("dpi", 300)

        # Make path absolute if relative
        if not os.path.isabs(filename):
            filename = os.path.join(self.output_folder, filename)

        try:
            if ray:
                cmd.ray(width, height)
            cmd.png(filename, width=width, height=height, dpi=dpi)
            print(f"PNG: Saved image to {filename}")
        except Exception as e:
            print(f"Error saving PNG: {e}")

    def execute_render(self, op: Dict[str, Any]):
        """
        Execute Render operation - orient, ray trace, and save PNG.

        Args:
            op: Operation dict with structures, orient_selection, dimensions, filename
        """
        orient_selection = op.get("orient_selection", "all")
        width = op.get("width", 1920)
        height = op.get("height", 1080)
        filename = op.get("filename", "render.png")
        dpi = op.get("dpi", 300)

        # Make path absolute if relative
        if not os.path.isabs(filename):
            filename = os.path.join(self.output_folder, filename)

        try:
            cmd.orient(orient_selection)
            cmd.zoom(orient_selection, buffer=5)
            cmd.ray(width, height)
            cmd.png(filename, width=width, height=height, dpi=dpi)
            print(f"Render: Saved render to {filename}")
        except Exception as e:
            print(f"Error rendering: {e}")

    def execute_render_each(self, op: Dict[str, Any]):
        """
        Execute RenderEach operation - render each structure individually as PNG.

        For each structure:
        1. Loads the structure
        2. Colors protein and ligand
        3. Orients view towards ligand
        4. Adds title (top) and caption (bottom) with metrics from table
        5. Ray traces and saves PNG

        Args:
            op: Operation dict with structures, coloring options, title/caption format, table reference
        """
        structures_ref = op.get("structures")
        orient_selection = op.get("orient_selection", "resn LIG")
        color_protein = op.get("color_protein", "plddt")
        color_ligand = op.get("color_ligand", "byatom")
        ligand_selection = op.get("ligand_selection", "resn LIG")
        plddt_upper = op.get("plddt_upper", 1)  # Default to 1 for Boltz2 confidence
        title_format = op.get("title")
        caption_format = op.get("caption")
        data_table_ref = op.get("data_table")
        width = op.get("width", 1920)
        height = op.get("height", 1080)
        dpi = op.get("dpi", 300)
        background = op.get("background", "white")

        if not structures_ref:
            raise ValueError("RenderEach requires structures")

        structures = self._resolve_structures(structures_ref)
        print(f"RenderEach: Rendering {len(structures)} structures individually")

        # Load data table if provided
        data_table = None
        if data_table_ref and (title_format or caption_format):
            table_path = None
            if isinstance(data_table_ref, dict):
                table_path = data_table_ref.get("table_path")
            if table_path and os.path.exists(table_path):
                data_table = self._load_table(table_path)
                print(f"  Loaded data table: {table_path} ({len(data_table)} rows)")

        # Create renders subfolder
        renders_folder = os.path.join(self.output_folder, "renders")
        os.makedirs(renders_folder, exist_ok=True)

        # Set background color
        cmd.bg_color(background)

        for struct in structures:
            struct_id = struct["id"]
            struct_path = struct["path"]

            if not os.path.exists(struct_path):
                print(f"  Skipping {struct_id}: file not found ({struct_path})")
                continue

            print(f"\n  Rendering: {struct_id}")

            # Clear previous objects
            cmd.delete("all")

            # Load structure
            try:
                cmd.load(struct_path, struct_id)
            except Exception as e:
                print(f"    Error loading {struct_path}: {e}")
                continue

            # Show cartoon for protein
            cmd.show("cartoon", f"{struct_id} and polymer")
            cmd.hide("lines", f"{struct_id}")

            # Color protein
            if color_protein == "plddt":
                self._apply_plddt_coloring(struct_id, "polymer", upper=plddt_upper)
            else:
                cmd.color(color_protein, f"{struct_id} and polymer")

            # Show and color ligand
            cmd.show("sticks", f"{struct_id} and {ligand_selection}")
            if color_ligand == "byatom":
                cmd.util.cbag(f"{struct_id} and {ligand_selection}")  # Color by atom, gray carbons
            else:
                cmd.color(color_ligand, f"{struct_id} and {ligand_selection}")

            # Orient towards ligand/selection
            try:
                # First orient the whole structure
                cmd.orient(struct_id)
                # Then zoom to show the ligand prominently
                if cmd.count_atoms(f"{struct_id} and {orient_selection}") > 0:
                    cmd.zoom(f"{struct_id} and {orient_selection}", buffer=15)
            except Exception as e:
                print(f"    Warning: Could not orient to {orient_selection}: {e}")
                cmd.orient(struct_id)

            # Build text info from table
            text_info = self._get_text_info(struct_id, title_format, caption_format, data_table)

            # Ray trace and save
            output_file = os.path.join(renders_folder, f"{struct_id}.png")
            try:
                cmd.ray(width, height)
                cmd.png(output_file, width=width, height=height, dpi=dpi)
                print(f"    Saved: {output_file}")

                # Add title/caption overlay to the saved PNG
                if text_info:
                    self._add_text_overlay(output_file, text_info, width, height)

            except Exception as e:
                print(f"    Error saving PNG: {e}")

        print(f"\nRenderEach: Completed {len(structures)} renders")

    def _apply_plddt_coloring(self, obj_name: str, selection: str = "polymer", upper: float = 100):
        """
        Apply pLDDT/confidence coloring to a selection using B-factor values.

        Uses the AlphaFold color scheme:
        - Blue (high confidence >= 90%): #126DFF
        - Cyan (70-90%): #0ECFF1
        - Yellow (50-70%): #F6ED12
        - Orange (low confidence < 50%): #EE831D

        Args:
            obj_name: PyMOL object name
            selection: Selection within object (e.g., "polymer")
            upper: Upper bound for confidence values (100 for AlphaFold pLDDT, 1 for Boltz2)
        """
        full_sel = f"{obj_name} and {selection}"

        # Calculate thresholds based on upper value
        # For upper=100: extremes = [50, 70, 90]
        # For upper=1: extremes = [0.5, 0.7, 0.9]
        if upper == 100:
            extremes = [int(float(upper) * x) for x in [0.5, 0.7, 0.9]]
        else:
            extremes = [float(upper) * x for x in [0.5, 0.7, 0.9]]

        try:
            # Color by confidence ranges using B-factor
            # Blue (high confidence >= 90%)
            cmd.color("0x126DFF", f"({full_sel}) and (b > {extremes[2]} or b = {extremes[2]})")
            # Cyan (70-90%)
            cmd.color("0x0ECFF1", f"({full_sel}) and ((b < {extremes[2]} and b > {extremes[1]}) or b = {extremes[1]})")
            # Yellow (50-70%)
            cmd.color("0xF6ED12", f"({full_sel}) and ((b < {extremes[1]} and b > {extremes[0]}) or b = {extremes[0]})")
            # Orange (low confidence < 50%)
            cmd.color("0xEE831D", f"({full_sel}) and b < {extremes[0]}")
            print(f"    Applied pLDDT coloring (upper={upper}, thresholds={extremes})")
        except Exception as e:
            print(f"    Warning: pLDDT coloring failed: {e}")

    def _get_text_info(self, struct_id: str, title_format: Optional[str],
                       caption_format: Optional[str], data_table: Optional[pd.DataFrame]) -> Optional[Dict[str, Any]]:
        """
        Get title and caption information from table for a structure.

        Args:
            struct_id: Structure ID to look up in table
            title_format: Format string for title (top) with {column_name} placeholders
            caption_format: Format string for caption (bottom) with {column_name} placeholders
            data_table: DataFrame with 'id' column and value columns

        Returns:
            Dict with 'title', 'caption', 'values', 'struct_id', or None if no text needed
        """
        if not title_format and not caption_format:
            return None

        values = {'id': struct_id}

        if data_table is not None and 'id' in data_table.columns:
            row = data_table[data_table['id'] == struct_id]
            if row.empty:
                # Try matching - check if struct_id contains table id or vice versa
                for col_id in data_table['id']:
                    col_id_str = str(col_id)
                    # Match rank numbers (e.g., rank0001_xxx matches rank0001)
                    if struct_id.startswith("rank") and col_id_str.startswith("rank"):
                        struct_rank = struct_id.split('_')[0]
                        col_rank = col_id_str.split('_')[0] if '_' in col_id_str else col_id_str
                        if struct_rank == col_rank:
                            row = data_table[data_table['id'] == col_id]
                            break
                    elif struct_id in col_id_str or col_id_str in struct_id:
                        row = data_table[data_table['id'] == col_id]
                        break

            if not row.empty:
                values = row.iloc[0].to_dict()
                values['id'] = struct_id  # Use actual struct_id
            else:
                print(f"    Warning: ID '{struct_id}' not found in data table")

        return {
            'title_format': title_format,
            'caption_format': caption_format,
            'values': values,
            'struct_id': struct_id
        }

    def _add_text_overlay(self, image_path: str, text_info: Dict[str, Any], width: int, height: int):
        """
        Add title (top, centered) and caption (bottom, centered) overlay to a rendered PNG image.

        Args:
            image_path: Path to the PNG file to modify
            text_info: Dict with 'title_format', 'caption_format', 'values', 'struct_id'
            width: Image width
            height: Image height
        """
        try:
            from PIL import Image, ImageDraw, ImageFont
        except ImportError:
            print("    Warning: PIL not available, skipping text overlay")
            return

        try:
            # Open the image
            img = Image.open(image_path)
            draw = ImageDraw.Draw(img)

            # Font sizes based on image dimensions
            font_size_title = max(48, width // 30)
            font_size_caption = max(32, width // 45)

            # Try to load fonts
            font_title = None
            font_caption = None
            font_paths = [
                '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf',
                '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',
                'DejaVuSans-Bold.ttf',
                'DejaVuSans.ttf',
                'Arial.ttf',
                'arial.ttf',
                '/System/Library/Fonts/Helvetica.ttc',
                'C:\\Windows\\Fonts\\arial.ttf',
            ]
            for font_path in font_paths:
                try:
                    font_title = ImageFont.truetype(font_path, font_size_title)
                    font_caption = ImageFont.truetype(font_path, font_size_caption)
                    break
                except:
                    continue

            if font_title is None:
                font_title = ImageFont.load_default()
                font_caption = ImageFont.load_default()

            title_format = text_info.get('title_format')
            caption_format = text_info.get('caption_format')
            values = text_info.get('values', {})

            # Draw title at top center
            if title_format:
                try:
                    title_text = title_format.format(**values)
                except Exception as e:
                    title_text = title_format
                    print(f"    Warning: Could not format title: {e}")

                bbox = draw.textbbox((0, 0), title_text, font=font_title)
                text_width = bbox[2] - bbox[0]
                text_height = bbox[3] - bbox[1]

                x = (width - text_width) // 2
                y = max(30, height // 25)

                # Draw shadow for readability on any background
                shadow_offset = 3
                draw.text((x + shadow_offset, y + shadow_offset), title_text, font=font_title, fill='white')
                draw.text((x - 1, y - 1), title_text, font=font_title, fill='white')
                draw.text((x + 1, y - 1), title_text, font=font_title, fill='white')
                draw.text((x - 1, y + 1), title_text, font=font_title, fill='white')
                draw.text((x + 1, y + 1), title_text, font=font_title, fill='white')
                draw.text((x, y), title_text, font=font_title, fill='black')

            # Draw caption at bottom center
            if caption_format:
                try:
                    caption_text = caption_format.format(**values)
                except Exception as e:
                    caption_text = caption_format
                    print(f"    Warning: Could not format caption: {e}")

                bbox = draw.textbbox((0, 0), caption_text, font=font_caption)
                text_width = bbox[2] - bbox[0]
                text_height = bbox[3] - bbox[1]

                x = (width - text_width) // 2
                y = height - text_height - max(30, height // 25)

                # Draw shadow for readability
                shadow_offset = 2
                draw.text((x + shadow_offset, y + shadow_offset), caption_text, font=font_caption, fill='white')
                draw.text((x - 1, y - 1), caption_text, font=font_caption, fill='white')
                draw.text((x + 1, y - 1), caption_text, font=font_caption, fill='white')
                draw.text((x - 1, y + 1), caption_text, font=font_caption, fill='white')
                draw.text((x + 1, y + 1), caption_text, font=font_caption, fill='white')
                draw.text((x, y), caption_text, font=font_caption, fill='dimgray')

            # Save the modified image
            img.save(image_path)
            print(f"    Added text overlay")

        except Exception as e:
            print(f"    Warning: Could not add text overlay: {e}")

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
        elif op_type == "orient":
            self.execute_orient(op)
        elif op_type == "ray":
            self.execute_ray(op)
        elif op_type == "png":
            self.execute_png(op)
        elif op_type == "render":
            self.execute_render(op)
        elif op_type == "render_each":
            self.execute_render_each(op)
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
