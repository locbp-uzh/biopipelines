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

    Each source is a dict with:
    - path: Path to source CSV
    - iterate: True = part of iteration (Each), False = static/bundled with each iteration
    """
    name: str                                    # Axis name (e.g., "proteins", "ligands")
    mode: str = "each"                           # "each" or "bundle"
    sources: List[Dict] = field(default_factory=list)  # [{"path": str, "iterate": bool}]

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

    def get_sources(self, axis_name: str) -> List[Dict]:
        """Get sources for an axis (list of dicts with 'path' and 'iterate' keys)."""
        axis = self.axes.get(axis_name)
        return axis.sources if axis else []

    def get_source_paths(self, axis_name: str) -> List[str]:
        """Get just the source paths for an axis."""
        axis = self.axes.get(axis_name)
        if not axis:
            return []
        return [s["path"] if isinstance(s, dict) else s for s in axis.sources]


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
        sources is a list of dicts: [{"path": str, "iterate": bool}]

    For Bundle(Each(library), cofactor):
    - mode = "bundle" (preserves the bundle semantics)
    - Each sources marked iterate=True
    - Bare sources marked iterate=False
    """
    if isinstance(value, Bundle):
        sources_with_iterate = []
        for src in value.sources:
            if isinstance(src, Each):
                # Each inside Bundle - these are iterated
                for sub_src in src.sources:
                    for path in _extract_source_paths(sub_src, axis_name):
                        sources_with_iterate.append({"path": path, "iterate": True})
            else:
                # Bare source inside Bundle - static (bundled with each iteration)
                for path in _extract_source_paths(src, axis_name):
                    sources_with_iterate.append({"path": path, "iterate": False})

        return ("bundle", sources_with_iterate)

    elif isinstance(value, Each):
        # Each iterates - all sources are iterated
        sources_with_iterate = []
        for src in value.sources:
            for path in _extract_source_paths(src, axis_name):
                sources_with_iterate.append({"path": path, "iterate": True})
        return ("each", sources_with_iterate)

    elif isinstance(value, list):
        # Bare list defaults to Each - all iterated
        sources_with_iterate = []
        for item in value:
            if isinstance(item, (Bundle, Each)):
                _, sub_sources = _unwrap_sources(item, axis_name)
                sources_with_iterate.extend(sub_sources)
            else:
                for path in _extract_source_paths(item, axis_name):
                    sources_with_iterate.append({"path": path, "iterate": True})
        return ("each", sources_with_iterate)

    else:
        # Bare value defaults to Each - iterated
        sources_with_iterate = []
        for path in _extract_source_paths(value, axis_name):
            sources_with_iterate.append({"path": path, "iterate": True})
        return ("each", sources_with_iterate)


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


def load_ids_from_sources(sources: List[Union[str, Dict]], iterate_only: bool = False) -> List[str]:
    """
    Load IDs from CSV source files at pipeline time.

    Args:
        sources: List of CSV file paths or dicts with {"path": str, "iterate": bool}
        iterate_only: If True, only load IDs from sources with iterate=True

    Returns:
        List of IDs from the 'id' column of the CSVs
    """
    import pandas as pd

    all_ids = []
    for source in sources:
        # Handle both string and dict formats
        if isinstance(source, dict):
            source_path = source.get("path")
            is_iterate = source.get("iterate", True)
            if iterate_only and not is_iterate:
                continue
        else:
            source_path = source

        if not source_path or not os.path.exists(source_path):
            continue
        try:
            df = pd.read_csv(source_path)
            if 'id' in df.columns:
                all_ids.extend(df['id'].tolist())
        except Exception as e:
            print(f"Warning: Could not load IDs from {source_path}: {e}")
    return all_ids


def _collect_ids_from_value(value: Any, axis_name: str, iterate_only: bool = False) -> List[str]:
    """
    Collect all IDs from a value, handling Bundle/Each wrappers.

    For Each(a, b), collects IDs from both a and b.
    For Bundle(Each(library), cofactor), collects IDs from the Each (library) only.
    For Bundle(cofactor, Each(library)), finds the Each and collects its IDs only.
    For pure Bundle(a, b), collects IDs from first source only (single bundled config).

    Args:
        value: A value that may be wrapped in Bundle/Each
        axis_name: "sequences" or "compounds" to determine which stream to access
        iterate_only: If True, only collect IDs from iterated sources (used for Bundle mode)

    Returns:
        List of all IDs from the relevant sources
    """
    stream_name = "sequences" if axis_name == "sequences" else "compounds"

    def get_ids_from_source(src: Any) -> List[str]:
        """Get IDs from a single unwrapped source."""
        if hasattr(src, 'streams'):
            stream = getattr(src.streams, stream_name, None)
            if stream:
                return list(stream.ids)
        if isinstance(src, str):
            return ["sequence" if axis_name == "sequences" else "compound"]
        return []

    if isinstance(value, Bundle):
        # Check if there's an Each inside - if so, only collect IDs from the Each sources
        has_each_inside = any(isinstance(src, Each) for src in value.sources)

        if has_each_inside:
            # Collect IDs only from the Each sources (these are what we iterate over)
            all_ids = []
            for src in value.sources:
                if isinstance(src, Each):
                    for each_src in src.sources:
                        all_ids.extend(_collect_ids_from_value(each_src, axis_name))
            return all_ids
        else:
            # Pure Bundle - return IDs from first source (single bundled config)
            return _collect_ids_from_value(value.sources[0], axis_name)

    elif isinstance(value, Each):
        # Collect IDs from all sources in the Each
        all_ids = []
        for src in value.sources:
            all_ids.extend(_collect_ids_from_value(src, axis_name))
        return all_ids

    else:
        # Bare value - get IDs directly
        return get_ids_from_source(value)


def predict_output_ids(
    bundled_name: str = "bundled_complex",
    **named_inputs: Any
) -> List[str]:
    """
    Predict the output IDs that will be generated by combinatorics.

    This function is called at pipeline time to predict what structure IDs
    will be generated at SLURM time based on the combinatorics modes.

    The ID generation follows these rules:
    - If all axes are fully bundled (no iterate sources): single ID (bundled_name)
    - If one axis has iteration (each mode or bundle with iterate sources): IDs from that axis
    - If multiple axes have iteration: cartesian product of IDs joined with "_"

    Args:
        bundled_name: Name to use when all axes are bundled (default: "bundled_complex")
        **named_inputs: Named inputs, each may be wrapped in Bundle/Each or bare

    Returns:
        List of expected output IDs
    """
    # Collect axis info: (name, mode, sources, ids)
    axes_info = []
    for name, value in named_inputs.items():
        if value is None:
            continue
        mode, sources = _unwrap_sources(value, name)

        # Check if this axis has any iterated sources
        has_iteration = any(s.get("iterate", True) for s in sources)

        # Collect IDs from the value - for Bundle with nested Each, this gets only the Each IDs
        ids = _collect_ids_from_value(value, name)
        if not ids:
            raise ValueError(f"No {name} found in input")

        axes_info.append((name, mode, has_iteration, ids))

    if not axes_info:
        return [bundled_name]

    # Separate into axes with iteration and axes without iteration
    # An axis has iteration if: mode=="each" OR (mode=="bundle" AND has iterated sources)
    iterated_axes = [(name, ids) for name, mode, has_iteration, ids in axes_info
                     if mode == "each" or has_iteration]
    pure_bundle_axes = [(name, ids) for name, mode, has_iteration, ids in axes_info
                        if mode == "bundle" and not has_iteration]

    # If no iterated axes, return single bundled name
    if not iterated_axes:
        return [bundled_name]

    # If single iterated axis, return its IDs
    if len(iterated_axes) == 1:
        return iterated_axes[0][1]

    # Multiple iterated axes: cartesian product
    # Special case: if one axis has only 1 element, use the other axis's IDs directly
    # This matches the behavior in pipe_boltz_config_unified.py generate_config_id()
    if len(iterated_axes) == 2:
        name0, ids0 = iterated_axes[0]
        name1, ids1 = iterated_axes[1]

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

    # General case: cartesian product of all iterated axes
    result_ids = iterated_axes[0][1]

    # Combine with remaining axes
    for _, ids in iterated_axes[1:]:
        new_result = []
        for existing_id in result_ids:
            for new_id in ids:
                new_result.append(f"{existing_id}_{new_id}")
        result_ids = new_result

    return result_ids
