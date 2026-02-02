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
        PyMOL(
            PyMOL.Names(prefix="sensor", basename=fuse.tables.sequences.lengths, suffix="parts"),
            PyMOL.Load(folded),
            PyMOL.Color(folded, selection=fuse.tables.sequences.D1, color="white"),
            PyMOL.Color(folded, selection=fuse.tables.sequences.L1, color="orange"),
            PyMOL.Align("align")
        )
    """

    TOOL_NAME = "PyMOL"

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
    def Color(structures: Union[StandardizedOutput, DataStream],
              selection: Union[Tuple[TableInfo, str], str],
              color: str) -> PyMOLOperation:
        """
        Color a selection on structures.

        For each structure, uses its ID to:
        1. Look up the PyMOL name from the current naming map
        2. Look up the selection value from the table column
        Then applies: color {color}, {pymol_name} and resi {selection}

        Args:
            structures: Tool output containing structures to color
            selection: Table column reference for per-structure selection
                       (e.g., fuse.tables.sequences.D1) or fixed selection string
            color: PyMOL color name (e.g., "white", "marine", "orange")

        Returns:
            PyMOLOperation for coloring
        """
        return PyMOLOperation("color", structures=structures, selection=selection, color=color)

    @staticmethod
    def ColorAF(structures: Union[StandardizedOutput, DataStream],
                upper: float = 100) -> PyMOLOperation:
        """
        Color structures by AlphaFold pLDDT (B-factor spectrum).

        Applies spectrum coloring based on B-factor values, which contain
        pLDDT scores in AlphaFold-generated structures.

        Args:
            structures: Tool output containing AlphaFold structures
            upper: Maximum B-factor value for scaling (default: 100 for pLDDT)

        Returns:
            PyMOLOperation for pLDDT coloring
        """
        return PyMOLOperation("coloraf", structures=structures, upper=upper)

    @staticmethod
    def Align(method: str = "align", target: Optional[str] = None) -> PyMOLOperation:
        """
        Align all loaded objects.

        If target is not specified, aligns all objects to the first loaded object.

        Args:
            method: Alignment method - "align", "super", or "cealign"
            target: Target object name for alignment (default: first loaded)

        Returns:
            PyMOLOperation for alignment
        """
        return PyMOLOperation("align", method=method, target=target)

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
    def Orient(selection: str = "all", zoom: bool = True) -> PyMOLOperation:
        """
        Orient the view to focus on a selection.

        Args:
            selection: PyMOL selection to orient towards (e.g., "resn LIG", "chain A")
            zoom: Whether to also zoom to the selection

        Returns:
            PyMOLOperation for orienting
        """
        return PyMOLOperation("orient", selection=selection, zoom=zoom)

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
                   orient_selection: str = "hetatm",
                   color_protein: str = "plddt",
                   color_ligand: str = "byatom",
                   ligand_selection: str = "hetatm",
                   title: Optional[str] = None,
                   title_table: Optional[TableInfo] = None,
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
        5. Adds title with optional metrics from table
        6. Ray traces and saves PNG

        Args:
            structures: Structures to render
            orient_selection: Selection to orient towards (default: "hetatm" for ligand)
            color_protein: Protein coloring - "plddt" for AlphaFold coloring, or color name
            color_ligand: Ligand coloring - "byatom" for element colors, or color name
            ligand_selection: Selection for ligand (default: "hetatm")
            title: Title format string with placeholders like {id}, {antibiotic},
                   {affinity_probability_binary:.2f}, etc. Values come from title_table.
            title_table: Table containing values for title placeholders (must have 'id' column)
            width: Image width in pixels
            height: Image height in pixels
            dpi: Image DPI
            background: Background color (default: "white")

        Returns:
            PyMOLOperation for iterative rendering

        Example:
            PyMOL.RenderEach(
                structures=boltzgen_out,
                orient_selection="hetatm",
                color_protein="plddt",
                color_ligand="byatom",
                title="Gentamicin - Affinity: {affinity_probability_binary:.2f}",
                title_table=boltzgen_out.tables.final_designs_metrics,
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
                              title=title,
                              title_table=title_table,
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

                title_table = op.params.get("title_table")
                if title_table is not None and hasattr(title_table, 'path'):
                    self._table_references.append((title_table, None))

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
            elif isinstance(value, StandardizedOutput):
                # Serialize StandardizedOutput reference
                result[key] = {
                    "type": "standardized_output",
                    "output_folder": value.output_folder,
                    "structures": value.structures,
                    "structure_ids": value.structure_ids
                }
            elif isinstance(value, TableInfo):
                # Direct TableInfo reference (for title_table in RenderEach)
                result[key] = {
                    "type": "table_info",
                    "table_path": value.path,
                    "columns": value.columns if hasattr(value, 'columns') else []
                }
            elif isinstance(value, tuple) and len(value) == 2:
                # Table column reference: (TableInfo, column_name)
                table_info, column_name = value
                if hasattr(table_info, 'path'):
                    result[key] = {
                        "type": "table_column",
                        "table_path": table_info.path,
                        "column_name": column_name
                    }
                else:
                    result[key] = str(value)
            elif hasattr(value, 'output_folder'):
                # ToolOutput or similar
                result[key] = {
                    "type": "tool_output",
                    "output_folder": value.output_folder,
                    "structures": getattr(value, 'structures', []),
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

        return f"""echo "Creating PyMOL session..."
echo "Output: {self.session_file}"

mkdir -p "{self.output_folder}"

cat > "{self.config_file}" << 'PYMOL_CONFIG_EOF'
{json.dumps(config, indent=2)}
PYMOL_CONFIG_EOF

python "{self.pymol_py}" --config "{self.config_file}"

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files."""
        return {
            "structures": DataStream.empty("structures", "pdb"),
            "sequences": DataStream.empty("sequences", "fasta"),
            "compounds": DataStream.empty("compounds", "sdf"),
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
