"""
Combinatorics module for Bundle/Each input handling in BioPipelines.

Provides wrapper classes and config generation for controlling how multiple
inputs are combined in tools.

Design:
- Each: Iterate over elements (cartesian product axis) - DEFAULT behavior
- Bundle: Group elements together as one entity
- Pipeline time: Generate a combinatorics config file describing groupings
- SLURM time: HelpScript reads config and generates appropriate configs

Example usage:
    outputA = Tool1()  # a1, a2
    outputB = Tool2()  # b1, b2

    # Default (Each): cartesian product a1+b1, a1+b2, a2+b1, a2+b2
    Tool3(input1=outputA, input2=outputB)

    # Bundled: all together a1+a2+b1+b2
    Tool3(input1=Bundle(outputA), input2=Bundle(outputB))

    # Mixed: a1+(b1+b2), a2+(b1+b2)
    Tool3(input1=outputA, input2=Bundle(outputB))
"""

import os
import json
from dataclasses import dataclass, field, asdict
from typing import List, Any, Union, Optional, Dict


class Each:
    """
    Marks input for iteration (cartesian product axis).
    This is the DEFAULT behavior for bare inputs.
    """
    def __init__(self, *sources):
        if not sources:
            raise ValueError("Each requires at least one source")
        self.sources = sources

    def __repr__(self) -> str:
        return f"Each({', '.join(repr(s) for s in self.sources)})"


class Bundle:
    """
    Groups multiple inputs together as one entity.
    All elements appear together in one config.
    """
    def __init__(self, *sources):
        if not sources:
            raise ValueError("Bundle requires at least one source")
        self.sources = sources

    def __repr__(self) -> str:
        return f"Bundle({', '.join(repr(s) for s in self.sources)})"


@dataclass
class AxisConfig:
    """
    Configuration for a single combinatorics axis.
    """
    name: str                                    # Axis name (e.g., "proteins", "ligands")
    mode: str = "each"                           # "each" or "bundle"
    sources: List[str] = field(default_factory=list)  # Paths to source CSVs

    def to_dict(self) -> Dict:
        return asdict(self)


@dataclass
class CombinatoricsConfig:
    """
    Configuration describing how multiple input axes should be combined.
    This gets written to a JSON file for HelpScript consumption.

    Generic design: supports any number of named axes.
    """
    axes: Dict[str, AxisConfig] = field(default_factory=dict)

    def to_dict(self) -> Dict:
        return {
            "axes": {name: axis.to_dict() for name, axis in self.axes.items()}
        }

    def save(self, path: str):
        """Save config to JSON file."""
        with open(path, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def load(cls, path: str) -> 'CombinatoricsConfig':
        """Load config from JSON file."""
        with open(path, 'r') as f:
            data = json.load(f)
        axes = {}
        for name, axis_data in data.get("axes", {}).items():
            axes[name] = AxisConfig(**axis_data)
        return cls(axes=axes)

    def get_axis(self, name: str) -> Optional[AxisConfig]:
        """Get axis config by name."""
        return self.axes.get(name)

    def get_mode(self, axis_name: str) -> str:
        """Get mode for an axis, defaults to 'each'."""
        axis = self.axes.get(axis_name)
        return axis.mode if axis else "each"

    def get_sources(self, axis_name: str) -> List[str]:
        """Get sources for an axis."""
        axis = self.axes.get(axis_name)
        return axis.sources if axis else []


def _extract_source_paths(source: Any, axis_name: str) -> List[str]:
    """
    Extract CSV file paths from a source (StandardizedOutput, ToolOutput, DataStream, etc.).

    Args:
        source: The source object to extract paths from
        axis_name: Which data type to extract ("proteins" -> sequences, "ligands" -> compounds)

    Returns list of paths to CSV files containing the data.
    """
    if source is None:
        return []

    # Direct path string
    if isinstance(source, str):
        if source.endswith('.csv'):
            return [source]
        return []

    # List of sources
    if isinstance(source, list):
        paths = []
        for item in source:
            paths.extend(_extract_source_paths(item, axis_name))
        return paths

    # DataStream - use map_table path
    if hasattr(source, 'map_table') and hasattr(source, 'ids'):
        if source.map_table:
            return [source.map_table]
        raise ValueError(f"DataStream for axis '{axis_name}' has no map_table")

    # StandardizedOutput - use axis_name to determine which DataStream to use
    if axis_name == "sequences":
        if hasattr(source, 'streams') and source.streams.sequences:
            if hasattr(source.streams.sequences, 'map_table') and source.streams.sequences.map_table:
                return [source.streams.sequences.map_table]
        raise ValueError(f"Source for 'sequences' axis must have sequences with map_table")
    elif axis_name == "compounds":
        if hasattr(source, 'streams') and source.streams.compounds:
            if hasattr(source.streams.compounds, 'map_table') and source.streams.compounds.map_table:
                return [source.streams.compounds.map_table]
        raise ValueError(f"Source for 'compounds' axis must have compounds with map_table")
    else:
        raise ValueError(f"Unknown axis name: {axis_name}. Must be 'sequences' or 'compounds'")


def _unwrap_sources(value: Any, axis_name: str) -> tuple:
    """
    Unwrap Bundle/Each wrappers and return (mode, sources).

    Args:
        value: The value to unwrap
        axis_name: Which data type to extract ("proteins" -> sequences, "ligands" -> compounds)

    Returns:
        (mode: str, sources: list) where mode is "bundle" or "each"
    """
    if isinstance(value, Bundle):
        # Collect all paths from bundle sources
        all_paths = []
        for src in value.sources:
            if isinstance(src, Each):
                # Each inside Bundle - flatten
                for sub_src in src.sources:
                    all_paths.extend(_extract_source_paths(sub_src, axis_name))
            else:
                all_paths.extend(_extract_source_paths(src, axis_name))
        return ("bundle", all_paths)

    elif isinstance(value, Each):
        # Each iterates
        all_paths = []
        for src in value.sources:
            all_paths.extend(_extract_source_paths(src, axis_name))
        return ("each", all_paths)

    elif isinstance(value, list):
        # Bare list defaults to Each
        all_paths = []
        for item in value:
            if isinstance(item, (Bundle, Each)):
                _, paths = _unwrap_sources(item, axis_name)
                all_paths.extend(paths)
            else:
                all_paths.extend(_extract_source_paths(item, axis_name))
        return ("each", all_paths)

    else:
        # Bare value defaults to Each
        paths = _extract_source_paths(value, axis_name)
        return ("each", paths)


def generate_combinatorics_config(
    output_path: str,
    **named_inputs: Any
) -> CombinatoricsConfig:
    """
    Generate combinatorics config file from named inputs.

    This is the main entry point called by tools.

    Args:
        output_path: Path to write the config JSON file
        **named_inputs: Named inputs using data type names (sequences, compounds)
                       e.g., sequences=tool1, compounds=Bundle(tool2)

    Returns:
        CombinatoricsConfig object (also saved to output_path)

    Example:
        generate_combinatorics_config(
            "config.json",
            sequences=lmpnn,
            compounds=Bundle(compounds)
        )
    """
    axes = {}
    for name, value in named_inputs.items():
        if value is None:
            continue
        mode, sources = _unwrap_sources(value, name)
        axes[name] = AxisConfig(name=name, mode=mode, sources=sources)

    config = CombinatoricsConfig(axes=axes)
    config.save(output_path)
    return config


def is_combinatorics_wrapper(value: Any) -> bool:
    """Check if a value is wrapped in Bundle or Each."""
    return isinstance(value, (Bundle, Each))


def contains_combinatorics_wrapper(value: Any) -> bool:
    """Check if a value contains any Bundle or Each wrappers."""
    if isinstance(value, (Bundle, Each)):
        return True
    if isinstance(value, list):
        return any(contains_combinatorics_wrapper(item) for item in value)
    return False


def get_mode(value: Any) -> str:
    """
    Get the combinatorics mode for a value.

    Returns "bundle" if wrapped in Bundle, "each" otherwise.
    """
    if isinstance(value, Bundle):
        return "bundle"
    return "each"


def load_ids_from_sources(sources: List[str]) -> List[str]:
    """
    Load IDs from CSV source files at pipeline time.

    Args:
        sources: List of CSV file paths

    Returns:
        List of IDs from the 'id' column of the CSVs
    """
    import pandas as pd

    all_ids = []
    for source_path in sources:
        if not source_path or not os.path.exists(source_path):
            continue
        try:
            df = pd.read_csv(source_path)
            if 'id' in df.columns:
                all_ids.extend(df['id'].tolist())
        except Exception as e:
            print(f"Warning: Could not load IDs from {source_path}: {e}")
    return all_ids


def _unwrap_to_first_source(value: Any) -> Any:
    """
    Recursively unwrap Bundle/Each wrappers to get the first actual source.

    Args:
        value: A value that may be wrapped in Bundle/Each

    Returns:
        The first non-Bundle/Each source found
    """
    while isinstance(value, (Bundle, Each)):
        if not value.sources:
            raise ValueError("Empty Bundle/Each wrapper")
        value = value.sources[0]
    return value


def predict_output_ids(
    bundled_name: str = "bundled_complex",
    **named_inputs: Any
) -> List[str]:
    """
    Predict the output IDs that will be generated by combinatorics.

    This function is called at pipeline time to predict what structure IDs
    will be generated at SLURM time based on the combinatorics modes.

    The ID generation follows these rules:
    - If all axes are bundled: single ID (bundled_name)
    - If one axis is "each" and others are "bundle": IDs from the "each" axis
    - If multiple axes are "each": cartesian product of IDs joined with "_"

    Args:
        bundled_name: Name to use when all axes are bundled (default: "bundled_complex")
        **named_inputs: Named inputs, each may be wrapped in Bundle/Each or bare

    Returns:
        List of expected output IDs
    """
    # Collect axis info: (name, mode, ids)
    axes_info = []
    for name, value in named_inputs.items():
        if value is None:
            continue
        mode, _ = _unwrap_sources(value, name)

        # Get IDs directly from source object - recursively unwrap to get actual source
        unwrapped = _unwrap_to_first_source(value)
        if name == "sequences":
            if hasattr(unwrapped, 'streams') and unwrapped.streams.sequences:
                ids = list(unwrapped.streams.sequences.ids)
            elif isinstance(unwrapped, str):
                # Direct sequence string - single ID
                ids = ["sequence"]
            else:
                raise ValueError(f"No sequences found in {name} input")
        elif name == "compounds":
            if hasattr(unwrapped, 'streams') and unwrapped.streams.compounds:
                ids = list(unwrapped.streams.compounds.ids)
            elif isinstance(unwrapped, str):
                # Direct SMILES string - single compound
                ids = ["compound"]
            else:
                raise ValueError(f"No compounds found in {name} input")
        else:
            raise ValueError(f"Unknown axis name: {name}. Must be 'sequences' or 'compounds'")

        axes_info.append((name, mode, ids))

    if not axes_info:
        return [bundled_name]

    # Separate into "each" and "bundle" axes
    each_axes = [(name, ids) for name, mode, ids in axes_info if mode == "each"]
    bundle_axes = [(name, ids) for name, mode, ids in axes_info if mode == "bundle"]

    # If all bundled, return single bundled name
    if not each_axes:
        return [bundled_name]

    # If single "each" axis, return its IDs
    if len(each_axes) == 1:
        return each_axes[0][1]

    # Multiple "each" axes: cartesian product
    # Special case: if one axis has only 1 element, use the other axis's IDs directly
    # This matches the behavior in pipe_boltz_config_unified.py generate_config_id()
    if len(each_axes) == 2:
        name0, ids0 = each_axes[0]
        name1, ids1 = each_axes[1]

        if len(ids0) == 1:
            # Single protein, multiple ligands: use ligand IDs
            return ids1
        elif len(ids1) == 1:
            # Multiple proteins, single ligand: use protein IDs
            return ids0
        else:
            # Multiple proteins, multiple ligands: cartesian product
            result_ids = []
            for id0 in ids0:
                for id1 in ids1:
                    result_ids.append(f"{id0}_{id1}")
            return result_ids

    # General case: cartesian product of all axes
    result_ids = each_axes[0][1]

    # Combine with remaining axes
    for _, ids in each_axes[1:]:
        new_result = []
        for existing_id in result_ids:
            for new_id in ids:
                new_result.append(f"{existing_id}_{new_id}")
        result_ids = new_result

    return result_ids
