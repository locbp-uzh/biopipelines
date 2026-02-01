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


def _extract_source_paths(source: Any) -> List[str]:
    """
    Extract CSV file paths from a source (StandardizedOutput, ToolOutput, etc.).

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
            paths.extend(_extract_source_paths(item))
        return paths

    # StandardizedOutput - get sequences or compounds CSV
    if hasattr(source, 'sequences') and source.sequences:
        return source.sequences[:1]  # First sequence file
    if hasattr(source, 'compounds') and source.compounds:
        return source.compounds[:1]  # First compounds file

    # ToolOutput
    if hasattr(source, 'get_output_files'):
        outputs = source.get_output_files()
        if 'sequences' in outputs and outputs['sequences']:
            return outputs['sequences'][:1]
        if 'compounds' in outputs and outputs['compounds']:
            return outputs['compounds'][:1]

    return []


def _unwrap_sources(value: Any) -> tuple:
    """
    Unwrap Bundle/Each wrappers and return (mode, sources).

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
                    all_paths.extend(_extract_source_paths(sub_src))
            else:
                all_paths.extend(_extract_source_paths(src))
        return ("bundle", all_paths)

    elif isinstance(value, Each):
        # Each iterates
        all_paths = []
        for src in value.sources:
            all_paths.extend(_extract_source_paths(src))
        return ("each", all_paths)

    elif isinstance(value, list):
        # Bare list defaults to Each
        all_paths = []
        for item in value:
            if isinstance(item, (Bundle, Each)):
                _, paths = _unwrap_sources(item)
                all_paths.extend(paths)
            else:
                all_paths.extend(_extract_source_paths(item))
        return ("each", all_paths)

    else:
        # Bare value defaults to Each
        paths = _extract_source_paths(value)
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
        **named_inputs: Named inputs, each may be wrapped in Bundle/Each or bare
                       e.g., proteins=tool1, ligands=Bundle(tool2)

    Returns:
        CombinatoricsConfig object (also saved to output_path)

    Example:
        generate_combinatorics_config(
            "config.json",
            proteins=lmpnn,
            ligands=Bundle(compounds)
        )
    """
    axes = {}
    for name, value in named_inputs.items():
        if value is None:
            continue
        mode, sources = _unwrap_sources(value)
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
        mode, _ = _unwrap_sources(value)

        # Get IDs directly from source object
        unwrapped = value.sources[0] if isinstance(value, (Bundle, Each)) else value
        if name == "proteins":
            # DataStream refactor: sequences.ids instead of sequence_ids
            if hasattr(unwrapped, 'sequences') and unwrapped.sequences:
                ids = list(unwrapped.sequences.ids)
            elif isinstance(unwrapped, str):
                # Direct sequence string - single ID
                ids = ["protein"]
            else:
                raise ValueError(f"No sequences found in {name} input")
        elif name == "ligands":
            # DataStream refactor: compounds.ids instead of compound_ids
            if hasattr(unwrapped, 'compounds') and unwrapped.compounds:
                ids = list(unwrapped.compounds.ids)
            elif isinstance(unwrapped, str):
                # Direct SMILES string - single ligand, bundled by default
                ids = ["ligand"]
            else:
                raise ValueError(f"No compounds found in {name} input")
        else:
            raise ValueError(f"Unknown axis name: {name}")

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
