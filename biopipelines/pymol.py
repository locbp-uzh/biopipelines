# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
PyMOL visualization tool with declarative operation-based API.

Creates PyMOL session files (.pse) from structure outputs using a sequence
of operations like Load, Color, Align, etc. Supports ID-based matching
between structures and table columns for per-structure selections.
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


class PyMOLOperation:
    """Base class for PyMOL operations."""

    def __init__(self, op_type: str, **kwargs):
        self.op_type = op_type
        self.params = kwargs

    def to_dict(self) -> Dict[str, Any]:
        """Convert operation to dictionary for serialization."""
        return {"op": self.op_type, **self.params}


class PyMOL(BaseConfig):
    """
    PyMOL visualization tool with declarative operation-based API.

    Uses a sequence of operations (Names, Load, Color, Align, etc.) to build
    PyMOL sessions. Supports ID-based matching between structures and table
    columns for per-structure selections and naming.

    Example:
        # Color by domain and align on calmodulin (D2)
        PyMOL(
            PyMOL.Names(prefix="sensor", basename=fuse.tables.sequences.lengths, suffix="parts"),
            PyMOL.Load(apo),
            PyMOL.Load(holo),
            PyMOL.Color("marine", selection=fuse.tables.sequences.D1),
            PyMOL.Color("gray50", selection=fuse.tables.sequences.L1),
            PyMOL.Color("white", selection=fuse.tables.sequences.D2),
            PyMOL.Color("orange"),  # color everything else
            PyMOL.Align(selection=fuse.tables.sequences.D2),
        )
    """

    TOOL_NAME = "PyMOL"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        biopipelines = folders.get("biopipelines", "")
        if env_manager == "pip":
            skip = "" if force_reinstall else """# Check if already installed
if python -c "import pymol" 2>/dev/null; then
    echo "PyMOL already installed, skipping. Use force_reinstall=True to reinstall."
    exit 0
fi
"""
            return f"""echo "=== Installing PyMOL (pip) ==="
{skip}pip install pymol-open-source

echo "=== PyMOL installation complete ==="
"""
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_manager} env list 2>/dev/null | grep -q "ProteinEnv"; then
    echo "PyMOL (ProteinEnv) already installed, skipping. Use force_reinstall=True to reinstall."
    exit 0
fi
"""
        return f"""echo "=== Installing PyMOL (ProteinEnv) ==="
{skip}{env_manager} env create -f {biopipelines}/Environments/ProteinEnv.yaml
if [ $? -ne 0 ]; then
    echo "ERROR: ProteinEnv creation failed."
    exit 1
fi
echo "=== PyMOL (ProteinEnv) installation complete ==="
"""

    def _repr_notebook_html(self, output) -> str:
        """Render py3Dmol viewer with structures from Load operations."""
        from .datastream import DataStream

        # Collect structures from all sources referenced in Load operations
        all_ids = []
        all_files = []
        for source in self._structure_sources:
            if hasattr(source, "streams"):
                structures_ds = source.streams.get("structures")
                if isinstance(structures_ds, DataStream) and len(structures_ds) > 0:
                    for struct_id, file_path in structures_ds:
                        all_ids.append(struct_id)
                        all_files.append(file_path)

        if not all_ids:
            return ""

        combined = DataStream(
            name="structures",
            ids=all_ids,
            files=all_files,
            format="pdb",
        )
        return self._build_py3dmol_html(combined)

    # Lazy path descriptors
    config_file = Path(lambda self: os.path.join(self.output_folder, "pymol_config.json"))
    session_file = Path(lambda self: os.path.join(self.output_folder, f"{self.session_name}.pse"))
    pymol_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_pymol.py"))

    # --- Static methods for creating operations ---

    @staticmethod
    def Names(prefix: str = "", basename=None, suffix: str = "") -> PyMOLOperation:
        """
        Set up ID → PyMOL name mapping for subsequent operations.

        For each ID in the basename table column, creates a mapping:
        ID → "{prefix}_{basename_value}_{suffix}"

        Args:
            prefix: Prefix for PyMOL object names
            basename: Table column reference for the variable part of the name
                      (e.g., fuse.tables.sequences.lengths)
            suffix: Suffix for PyMOL object names

        Returns:
            PyMOLOperation for names mapping
        """
        return PyMOLOperation("names", prefix=prefix, basename=basename, suffix=suffix)

    @staticmethod
    def Load(structures: Union[StandardizedOutput, DataStream]) -> PyMOLOperation:
        """
        Load structures into PyMOL with current naming.

        For each structure, uses its ID to look up the PyMOL name from the
        current naming map (set by Names). If no mapping exists, uses the ID
        as the PyMOL object name.

        Args:
            structures: Tool output containing structures to load

        Returns:
            PyMOLOperation for loading structures
        """
        return PyMOLOperation("load", structures=structures)

    @staticmethod
    def Color(color: str,
              structures: Union[StandardizedOutput, DataStream, None] = None,
              selection: Union[Tuple[TableInfo, str], str, None] = None) -> PyMOLOperation:
        """
        Color structures or selections.

        Can be used in several ways:
        - Color(color="marine") - color all loaded objects
        - Color(color="marine", selection="resi 1-40") - color a fixed selection on all objects
        - Color(color="marine", selection=fuse.tables.sequences.D1) - per-structure selection on all objects
        - Color(color="marine", structures=apo, selection=...) - restrict to specific structures

        Args:
            color: PyMOL color name (e.g., "white", "marine", "orange")
            structures: Tool output containing structures to color (optional, defaults to all loaded)
            selection: Table column reference for per-structure selection
                       (e.g., fuse.tables.sequences.D1) or fixed selection string (optional)

        Returns:
            PyMOLOperation for coloring
        """
        return PyMOLOperation("color", structures=structures, selection=selection, color=color)

    @staticmethod
    def ColorAF(structures: Union[StandardizedOutput, DataStream, None] = None,
                upper: float = 100) -> PyMOLOperation:
        """
        Color structures by AlphaFold pLDDT (B-factor spectrum).

        Applies spectrum coloring based on B-factor values, which contain
        pLDDT scores in AlphaFold-generated structures.

        Args:
            structures: Tool output containing structures (optional, defaults to all loaded)
            upper: Maximum B-factor value for scaling (default: 100 for pLDDT)

        Returns:
            PyMOLOperation for pLDDT coloring
        """
        return PyMOLOperation("coloraf", structures=structures, upper=upper)

    @staticmethod
    def ColorAlign(reference: Union[StandardizedOutput, DataStream],
                   targets: Union[StandardizedOutput, DataStream],
                   identical: str = "white",
                   similar: str = "wheat",
                   different: str = "wheat",
                   notcovered: str = "gray50",
                   show_mutations: bool = True) -> PyMOLOperation:
        """
        Color structures based on sequence alignment against a reference.

        Performs sequence alignment between reference and target structures,
        then colors residues based on alignment quality:
        - identical: residues that match exactly
        - similar: residues in the same amino acid group (e.g., aliphatic, aromatic)
        - different: mismatched residues not in similar groups
        - notcovered: gaps/unaligned regions

        Also performs structural alignment using PyMOL's align command.

        Args:
            reference: Reference structure(s) for alignment
            targets: Target structure(s) to color based on alignment
            identical: Color for identical residues (default: "white")
            similar: Color for similar group matches (default: "wheat")
            different: Color for non-similar mismatches (default: "wheat")
            notcovered: Color for gaps/unaligned regions (default: "gray50")
            show_mutations: Show sticks for mutated residues (default: True)

        Returns:
            PyMOLOperation for alignment-based coloring
        """
        return PyMOLOperation("coloralign",
                              reference=reference,
                              targets=targets,
                              identical=identical,
                              similar=similar,
                              different=different,
                              notcovered=notcovered,
                              show_mutations=show_mutations)

    @staticmethod
    def Align(method: str = "align", target: Optional[str] = None,
              selection: Optional[Union[Tuple[TableInfo, str], str]] = None) -> PyMOLOperation:
        """
        Align all loaded objects.

        If target is not specified, aligns all objects to the first loaded object.
        If selection is provided, alignment is restricted to matching residues
        (e.g., align only on a specific domain).

        Args:
            method: Alignment method - "align", "super", or "cealign"
            target: Target object name for alignment (default: first loaded)
            selection: Optional residue selection to restrict alignment to.
                       Table column reference (e.g., fuse.tables.sequences.D2)
                       or fixed selection string (e.g., "resi 100-250")

        Returns:
            PyMOLOperation for alignment
        """
        return PyMOLOperation("align", method=method, target=target, selection=selection)

    @staticmethod
    def Show(structures: Union[StandardizedOutput, DataStream, None] = None,
             representation: str = "cartoon",
             selection: Optional[Union[Tuple[TableInfo, str], str]] = None) -> PyMOLOperation:
        """
        Show a representation for structures.

        Args:
            structures: Tool output containing structures (optional, defaults to all)
            representation: PyMOL representation (cartoon, sticks, surface, etc.)
            selection: Optional selection - table column reference or fixed string

        Returns:
            PyMOLOperation for showing representation
        """
        return PyMOLOperation("show", structures=structures, representation=representation, selection=selection)

    @staticmethod
    def Hide(structures: Union[StandardizedOutput, DataStream, None] = None,
             representation: str = "everything",
             selection: Optional[Union[Tuple[TableInfo, str], str]] = None) -> PyMOLOperation:
        """
        Hide a representation for structures.

        Args:
            structures: Tool output containing structures (optional, defaults to all)
            representation: PyMOL representation to hide
            selection: Optional selection - table column reference or fixed string

        Returns:
            PyMOLOperation for hiding representation
        """
        return PyMOLOperation("hide", structures=structures, representation=representation, selection=selection)

    @staticmethod
    def Set(setting: str, value: Any, selection: Optional[str] = None) -> PyMOLOperation:
        """
        Set a PyMOL setting.

        Args:
            setting: PyMOL setting name
            value: Setting value
            selection: Optional selection to apply to

        Returns:
            PyMOLOperation for setting
        """
        return PyMOLOperation("set", setting=setting, value=value, selection=selection)

    @staticmethod
    def Save(filename: str = "session.pse") -> PyMOLOperation:
        """
        Save the PyMOL session.

        Args:
            filename: Output filename (relative to output folder)

        Returns:
            PyMOLOperation for saving
        """
        return PyMOLOperation("save", filename=filename)

    @staticmethod
    def Center(selection: str = "all") -> PyMOLOperation:
        """
        Center the view on a selection without changing orientation.

        Args:
            selection: PyMOL selection to center on (default: "all")

        Returns:
            PyMOLOperation for centering
        """
        return PyMOLOperation("center", selection=selection)

    @staticmethod
    def Zoom(selection: str = "all", buffer: float = 5.0) -> PyMOLOperation:
        """
        Zoom the view to fit a selection.

        Args:
            selection: PyMOL selection to zoom to (default: "all")
            buffer: Extra space around selection in Angstroms (default: 5.0)

        Returns:
            PyMOLOperation for zooming
        """
        return PyMOLOperation("zoom", selection=selection, buffer=buffer)

    @staticmethod
    def Orient(selection: str = "all") -> PyMOLOperation:
        """
        Orient the view to show a selection from the best angle.

        Rotates the camera to align the principal axes of the selection
        with the screen axes.

        Args:
            selection: PyMOL selection to orient towards (default: "all")

        Returns:
            PyMOLOperation for orienting
        """
        return PyMOLOperation("orient", selection=selection)

    @staticmethod
    def Ray(width: int = 1920, height: int = 1080) -> PyMOLOperation:
        """
        Ray trace the current view.

        Args:
            width: Image width in pixels
            height: Image height in pixels

        Returns:
            PyMOLOperation for ray tracing
        """
        return PyMOLOperation("ray", width=width, height=height)

    @staticmethod
    def PNG(filename: str = "render.png", width: int = 1920, height: int = 1080,
            ray: bool = True, dpi: int = 300) -> PyMOLOperation:
        """
        Save current view as PNG image.

        Args:
            filename: Output filename (relative to output folder, supports {id} placeholder)
            width: Image width in pixels
            height: Image height in pixels
            ray: Whether to ray trace before saving
            dpi: Image DPI

        Returns:
            PyMOLOperation for saving PNG
        """
        return PyMOLOperation("png", filename=filename, width=width, height=height, ray=ray, dpi=dpi)

    @staticmethod
    def Render(structures: Union[StandardizedOutput, DataStream],
               orient_selection: str = "all",
               width: int = 1920,
               height: int = 1080,
               filename: str = "render.png",
               dpi: int = 300) -> PyMOLOperation:
        """
        Render a single image of loaded structures.

        Orients, ray traces, and saves a PNG of the current view.

        Args:
            structures: Structures to render (must be loaded first)
            orient_selection: Selection to orient towards (e.g., "resn LIG")
            width: Image width in pixels
            height: Image height in pixels
            filename: Output filename
            dpi: Image DPI

        Returns:
            PyMOLOperation for rendering
        """
        return PyMOLOperation("render", structures=structures, orient_selection=orient_selection,
                              width=width, height=height, filename=filename, dpi=dpi)

    @staticmethod
    def RenderEach(structures: Union[StandardizedOutput, DataStream],
                   orient_selection: str = "resn LIG",
                   color_protein: str = "plddt",
                   color_ligand: str = "byatom",
                   ligand_selection: str = "resn LIG",
                   plddt_upper: float = 1,
                   width: int = 1920,
                   height: int = 1080,
                   dpi: int = 300,
                   background: str = "white") -> PyMOLOperation:
        """
        Render each structure individually as a separate PNG.

        For each structure:
        1. Loads the structure
        2. Colors protein (by pLDDT or specified color)
        3. Colors ligand (by atom or specified color)
        4. Orients view towards the ligand/selection
        5. Ray traces and saves PNG

        Args:
            structures: Structures to render
            orient_selection: Selection to orient towards (default: "resn LIG" for ligand)
            color_protein: Protein coloring - "plddt" for confidence coloring, or color name
            color_ligand: Ligand coloring - "byatom" for element colors, or color name
            ligand_selection: Selection for ligand (default: "resn LIG")
            plddt_upper: Upper bound for pLDDT values (default: 1 for Boltz2, use 100 for AlphaFold)
            width: Image width in pixels
            height: Image height in pixels
            dpi: Image DPI
            background: Background color (default: "white")

        Returns:
            PyMOLOperation for iterative rendering

        Example:
            PyMOL.RenderEach(
                structures=boltz_rifampicin,
                orient_selection="resn LIG",
                color_protein="plddt",
                color_ligand="byatom",
                ligand_selection="resn LIG",
                plddt_upper=1,  # Boltz2 uses 0-1 confidence scores
                width=1920,
                height=1080
            )
        """
        return PyMOLOperation("render_each",
                              structures=structures,
                              orient_selection=orient_selection,
                              color_protein=color_protein,
                              color_ligand=color_ligand,
                              ligand_selection=ligand_selection,
                              plddt_upper=plddt_upper,
                              width=width,
                              height=height,
                              dpi=dpi,
                              background=background)

    # --- Instance methods ---

    def __init__(self, *args, session: str = "session", **kwargs):
        """
        Initialize PyMOL tool with a sequence of operations.

        Args:
            session: Name for the output session file (without .pse extension)
            *args: Sequence of PyMOL operations (Names, Load, Color, etc.)
            **kwargs: Additional configuration parameters

        Example:
            PyMOL(session="my_session",
                  PyMOL.Load(structures),
                  PyMOL.ColorAF(structures))
        """
        self.operations = list(args)
        self.session_name = session

        # Extract structure sources and table references from operations
        self._structure_sources = []
        self._table_references = []
        self._extract_dependencies()

        super().__init__(**kwargs)

    def _extract_dependencies(self):
        """Extract structure sources and table references from operations."""
        for op in self.operations:
            if op.op_type == "load":
                structures = op.params.get("structures")
                if structures is not None:
                    self._structure_sources.append(structures)

            elif op.op_type in ("color", "coloraf", "show", "hide", "render"):
                structures = op.params.get("structures")
                if structures is not None and structures not in self._structure_sources:
                    self._structure_sources.append(structures)

                selection = op.params.get("selection")
                if isinstance(selection, tuple):
                    self._table_references.append(selection)

            elif op.op_type == "render_each":
                structures = op.params.get("structures")
                if structures is not None and structures not in self._structure_sources:
                    self._structure_sources.append(structures)

            elif op.op_type == "align":
                selection = op.params.get("selection")
                if isinstance(selection, tuple):
                    self._table_references.append(selection)

            elif op.op_type == "coloralign":
                reference = op.params.get("reference")
                if reference is not None and reference not in self._structure_sources:
                    self._structure_sources.append(reference)
                targets = op.params.get("targets")
                if targets is not None and targets not in self._structure_sources:
                    self._structure_sources.append(targets)

            elif op.op_type == "names":
                basename = op.params.get("basename")
                if isinstance(basename, tuple):
                    self._table_references.append(basename)

    def validate_params(self):
        """Validate PyMOL parameters."""
        if not self.operations:
            raise ValueError("At least one operation must be provided")

        # Check that all operations are valid PyMOLOperation objects
        for i, op in enumerate(self.operations):
            if not isinstance(op, PyMOLOperation):
                raise ValueError(f"Operation {i} is not a PyMOLOperation: {type(op)}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input sources from pipeline context."""
        self.folders = pipeline_folders
        self.input_sources = {}

        # Track structure sources
        for i, source in enumerate(self._structure_sources):
            if hasattr(source, 'output_folder'):
                self.input_sources[f"structures_{i}"] = {
                    'output_folder': source.output_folder,
                    'structures': getattr(source, 'structures', []),
                    'structure_ids': getattr(source, 'structure_ids', [])
                }

    def _serialize_operation(self, op: PyMOLOperation) -> Dict[str, Any]:
        """Serialize an operation to a dictionary for JSON config."""
        result = {"op": op.op_type}

        for key, value in op.params.items():
            if value is None:
                result[key] = None
            elif isinstance(value, DataStream):
                # Serialize DataStream - use files and ids directly
                result[key] = {
                    "type": "datastream",
                    "structures": value.files,
                    "structure_ids": value.ids,
                    "format": value.format
                }
            elif isinstance(value, StandardizedOutput):
                # Serialize StandardizedOutput reference
                # Check if structures is a DataStream
                structures_data = value.streams.structures
                if isinstance(structures_data, DataStream):
                    result[key] = {
                        "type": "standardized_output",
                        "output_folder": value.output_folder,
                        "structures": structures_data.files,
                        "structure_ids": structures_data.ids
                    }
                else:
                    result[key] = {
                        "type": "standardized_output",
                        "output_folder": value.output_folder,
                        "structures": structures_data if isinstance(structures_data, list) else [],
                        "structure_ids": value.structure_ids if hasattr(value, 'structure_ids') else []
                    }
            elif isinstance(value, TableInfo):
                # Direct TableInfo reference (for data_table in RenderEach)
                result[key] = {
                    "type": "table_info",
                    "table_path": value.info.path,
                    "columns": value.info.columns
                }
            elif isinstance(value, tuple) and len(value) == 2:
                # Table column reference: (TableInfo, column_name)
                table_info, column_name = value
                if hasattr(table_info, 'info'):
                    result[key] = {
                        "type": "table_column",
                        "table_path": table_info.info.path,
                        "column_name": column_name
                    }
                else:
                    result[key] = str(value)
            elif hasattr(value, 'output_folder'):
                # ToolOutput or similar - check for DataStream in structures
                structures_attr = getattr(value, 'structures', None)
                if isinstance(structures_attr, DataStream):
                    result[key] = {
                        "type": "tool_output",
                        "output_folder": value.output_folder,
                        "structures": structures_attr.files,
                        "structure_ids": structures_attr.ids
                    }
                else:
                    result[key] = {
                        "type": "tool_output",
                        "output_folder": value.output_folder,
                        "structures": structures_attr if isinstance(structures_attr, list) else [],
                        "structure_ids": getattr(value, 'structure_ids', [])
                    }
            else:
                result[key] = value

        return result

    def generate_script(self, script_path: str) -> str:
        """Generate PyMOL execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        script_content = "#!/bin/bash\n"
        script_content += "# PyMOL execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self.generate_script_run_pymol()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_pymol(self) -> str:
        """Generate the PyMOL session creation part of the script."""
        config = {
            "operations": [self._serialize_operation(op) for op in self.operations],
            "session_name": self.session_name,
            "output_folder": self.output_folder
        }

        # Write config file at pipeline time (not SLURM time)
        os.makedirs(self.output_folder, exist_ok=True)
        with open(self.config_file, 'w') as f:
            json.dump(config, f, indent=2)

        return f"""echo "Creating PyMOL session..."
echo "Output: {self.session_file}"

python "{self.pymol_py}" --config "{self.config_file}"

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files including PNG renders as DataStream."""
        # Predict PNG files from Render and RenderEach operations
        render_ids = []
        render_files = []

        for op in self.operations:
            if op.op_type == "render":
                # Single render - predict the output filename
                filename = op.params.get("filename", "render.png")
                if not os.path.isabs(filename):
                    filename = os.path.join(self.output_folder, filename)
                # Use filename without extension as ID
                render_id = os.path.splitext(os.path.basename(filename))[0]
                render_ids.append(render_id)
                render_files.append(filename)

            elif op.op_type == "render_each":
                # RenderEach - predict PNG for each structure
                structures = op.params.get("structures")
                if structures is not None:
                    # Get structure IDs to predict render filenames
                    structure_ids = []
                    if isinstance(structures, DataStream):
                        structure_ids = structures.ids
                    elif isinstance(structures, StandardizedOutput):
                        structures_ds = structures.streams.structures
                        if isinstance(structures_ds, DataStream):
                            structure_ids = structures_ds.ids

                    # RenderEach saves to renders/<id>.png
                    renders_folder = os.path.join(self.output_folder, "renders")
                    for struct_id in structure_ids:
                        render_ids.append(struct_id)
                        render_files.append(os.path.join(renders_folder, f"{struct_id}.png"))

            elif op.op_type == "png":
                # Direct PNG operation
                filename = op.params.get("filename", "render.png")
                if not os.path.isabs(filename):
                    filename = os.path.join(self.output_folder, filename)
                render_id = os.path.splitext(os.path.basename(filename))[0]
                render_ids.append(render_id)
                render_files.append(filename)

        # Create renders DataStream if there are any render operations
        if render_ids:
            renders = DataStream(
                name="renders",
                ids=render_ids,
                files=render_files,
                format="png"
            )
        else:
            renders = DataStream.empty("renders", "png")

        return {
            "structures": DataStream.empty("structures", "pdb"),
            "sequences": DataStream.empty("sequences", "fasta"),
            "compounds": DataStream.empty("compounds", "sdf"),
            "renders": renders,
            "tables": {},
            "output_folder": self.output_folder,
            "session_file": self.session_file
        }

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"OPERATIONS: {len(self.operations)}",
            f"SESSION: {self.session_name}.pse"
        ])

        # Show operation summary
        op_types = [op.op_type for op in self.operations]
        config_lines.append(f"SEQUENCE: {' → '.join(op_types)}")

        return config_lines

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration to dictionary."""
        base_dict = super().to_dict()
        base_dict.update({
            "pymol_params": {
                "session_name": self.session_name,
                "num_operations": len(self.operations),
                "operation_types": [op.op_type for op in self.operations]
            }
        })
        return base_dict

    def __str__(self) -> str:
        """String representation."""
        op_types = [op.op_type for op in self.operations]
        return f"PyMOL({' → '.join(op_types)})"
