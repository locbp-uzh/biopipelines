# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

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
from typing import List, Any, Union, Optional, Dict, Tuple


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

    predicted_ids: List[str] = field(default_factory=list)
    provenance: Dict[str, List[str]] = field(default_factory=dict)

    def to_dict(self) -> Dict:
        result = {
            "axes": {name: axis.to_dict() for name, axis in self.axes.items()}
        }
        if self.predicted_ids:
            result["predicted_ids"] = self.predicted_ids
        if self.provenance:
            result["provenance"] = self.provenance
        return result

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
        return cls(
            axes=axes,
            predicted_ids=data.get("predicted_ids", []),
            provenance=data.get("provenance", {})
        )

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


def _unpack_input(name: str, value: Any) -> Tuple[Any, str]:
    """
    Unpack the (value, stream_name) tuple convention.

    Each named input can be either:
    - A bare value: the key serves as both alias and stream name
    - A (value, stream_name) tuple: key is alias, stream_name tells which stream to extract

    Args:
        name: The kwarg key (alias for provenance columns)
        value: The kwarg value — bare or (value, stream_name) tuple

    Returns:
        (actual_value, stream_name)
    """
    if isinstance(value, tuple) and len(value) == 2 and isinstance(value[1], str):
        # Check it's not a TableInfo tuple by verifying second element looks like a stream name
        actual_value, stream_name = value
        if not isinstance(actual_value, str):  # Not a (TableInfo, column_name) tuple
            return actual_value, stream_name
    return value, name


def _extract_source_paths(source: Any, stream_name: str) -> List[str]:
    """
    Extract CSV file paths from a source (StandardizedOutput, ToolOutput, DataStream, etc.).

    Args:
        source: The source object to extract paths from
        stream_name: Which stream to extract ("sequences", "compounds", "structures")

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
            paths.extend(_extract_source_paths(item, stream_name))
        return paths

    # DataStream - use map_table path
    if hasattr(source, 'map_table') and hasattr(source, 'ids'):
        if source.map_table:
            return [source.map_table]
        raise ValueError(f"DataStream for stream '{stream_name}' has no map_table")

    # StandardizedOutput - use stream_name to determine which DataStream to use
    if hasattr(source, 'streams'):
        stream = getattr(source.streams, stream_name, None)
        if stream and hasattr(stream, 'map_table') and stream.map_table:
            return [stream.map_table]
        raise ValueError(f"Source for '{stream_name}' stream must have {stream_name} with map_table")

    raise ValueError(f"Cannot extract source paths from {type(source)} for stream '{stream_name}'")


def _unwrap_sources(value: Any, stream_name: str) -> tuple:
    """
    Unwrap Bundle/Each wrappers and return (mode, sources).

    Args:
        value: The value to unwrap (already unpacked from tuple convention)
        stream_name: Which stream to extract ("sequences", "compounds", "structures")

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
        order = 0
        for src in value.sources:
            if isinstance(src, Each):
                # Each inside Bundle - these are iterated
                for sub_src in src.sources:
                    for path in _extract_source_paths(sub_src, stream_name):
                        sources_with_iterate.append({"path": path, "iterate": True, "order": order})
                        order += 1
            else:
                # Bare source inside Bundle - static (bundled with each iteration)
                for path in _extract_source_paths(src, stream_name):
                    sources_with_iterate.append({"path": path, "iterate": False, "order": order})
                    order += 1

        return ("bundle", sources_with_iterate)

    elif isinstance(value, Each):
        # Each iterates - all sources are iterated
        sources_with_iterate = []
        for src in value.sources:
            for path in _extract_source_paths(src, stream_name):
                sources_with_iterate.append({"path": path, "iterate": True})
        return ("each", sources_with_iterate)

    elif isinstance(value, list):
        # Bare list defaults to Each - all iterated
        sources_with_iterate = []
        for item in value:
            if isinstance(item, (Bundle, Each)):
                _, sub_sources = _unwrap_sources(item, stream_name)
                sources_with_iterate.extend(sub_sources)
            else:
                for path in _extract_source_paths(item, stream_name):
                    sources_with_iterate.append({"path": path, "iterate": True})
        return ("each", sources_with_iterate)

    else:
        # Bare value defaults to Each - iterated
        sources_with_iterate = []
        for path in _extract_source_paths(value, stream_name):
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
        **named_inputs: Named inputs — each value is either a bare value or
                       a (value, stream_name) tuple. The key is the alias
                       (used for provenance column naming).
                       e.g., proteins=(lmpnn, "sequences"), ligands=(compounds, "compounds")

    Returns:
        CombinatoricsConfig object (also saved to output_path)

    Example:
        generate_combinatorics_config(
            "config.json",
            proteins=(lmpnn, "sequences"),
            ligands=(Bundle(compounds), "compounds")
        )
    """
    axes = {}
    for name, raw_value in named_inputs.items():
        actual_value, stream_name = _unpack_input(name, raw_value)
        if actual_value is None:
            continue
        mode, sources = _unwrap_sources(actual_value, stream_name)
        axes[name] = AxisConfig(name=name, mode=mode, sources=sources)

    config = CombinatoricsConfig(axes=axes)

    # Compute and store predicted IDs and provenance
    try:
        predicted_ids, provenance = predict_output_ids_with_provenance(**named_inputs)
        config.predicted_ids = predicted_ids
        config.provenance = provenance
    except Exception:
        pass  # Non-critical: pipe script can still re-compute

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


def _collect_ids_from_value(value: Any, stream_name: str, iterate_only: bool = False) -> List[str]:
    """
    Collect all IDs from a value, handling Bundle/Each wrappers.

    For Each(a, b), collects IDs from both a and b.
    For Bundle(Each(library), cofactor), collects IDs from the Each (library) only.
    For Bundle(cofactor, Each(library)), finds the Each and collects its IDs only.
    For pure Bundle(a, b), collects IDs from first source only (single bundled config).

    Args:
        value: A value that may be wrapped in Bundle/Each (already unpacked from tuple)
        stream_name: Which stream to access ("sequences", "compounds", "structures")
        iterate_only: If True, only collect IDs from iterated sources (used for Bundle mode)

    Returns:
        List of all IDs from the relevant sources
    """
    def get_ids_from_source(src: Any) -> List[str]:
        """Get IDs from a single unwrapped source."""
        if hasattr(src, 'streams'):
            stream = getattr(src.streams, stream_name, None)
            if stream:
                return list(stream.ids)
        if isinstance(src, str):
            return [stream_name]
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
                        all_ids.extend(_collect_ids_from_value(each_src, stream_name))
            return all_ids
        else:
            # Pure Bundle - return IDs from first source (single bundled config)
            return _collect_ids_from_value(value.sources[0], stream_name)

    elif isinstance(value, Each):
        # Collect IDs from all sources in the Each
        all_ids = []
        for src in value.sources:
            all_ids.extend(_collect_ids_from_value(src, stream_name))
        return all_ids

    else:
        # Bare value - get IDs directly
        return get_ids_from_source(value)


def _collect_axes_info(
    bundled_name: str = "bundled_complex",
    **named_inputs: Any
) -> Tuple[List[tuple], List[tuple], List[tuple]]:
    """
    Collect axis info from named inputs and classify into iterated/pure_bundle.

    Each named_input value can be a bare value or a (value, stream_name) tuple.
    The key (name) is preserved as the provenance column alias.

    Returns:
        (axes_info, iterated_axes, pure_bundle_axes) where:
        - axes_info: [(name, mode, has_iteration, ids), ...]
        - iterated_axes: [(name, ids), ...] for axes with iteration
        - pure_bundle_axes: [(name, ids), ...] for pure bundle axes
    """
    axes_info = []
    for name, raw_value in named_inputs.items():
        actual_value, stream_name = _unpack_input(name, raw_value)
        if actual_value is None:
            continue
        mode, sources = _unwrap_sources(actual_value, stream_name)
        has_iteration = any(s.get("iterate", True) for s in sources)
        ids = _collect_ids_from_value(actual_value, stream_name)
        if not ids:
            raise ValueError(f"No {name} found in input")
        axes_info.append((name, mode, has_iteration, ids))

    iterated_axes = [(name, ids) for name, mode, has_iteration, ids in axes_info
                     if mode == "each" or has_iteration]
    pure_bundle_axes = [(name, ids) for name, mode, has_iteration, ids in axes_info
                        if mode == "bundle" and not has_iteration]

    return axes_info, iterated_axes, pure_bundle_axes


def predict_output_ids_with_provenance(
    bundled_name: str = "bundled_complex",
    **named_inputs: Any
) -> Tuple[List[str], Dict[str, List[str]]]:
    """
    Predict output IDs AND provenance mapping.

    ID generation always uses the full cartesian product of iterated axes
    joined with "_". Provenance columns track which input from each axis
    contributed to each output, enabling downstream joins on any axis.

    Each named_input value can be:
    - A bare value: the key serves as both alias and stream name
    - A (value, stream_name) tuple: key is alias (provenance column name),
      stream_name tells which stream to extract IDs from

    Returns:
        (output_ids, provenance) where provenance keys are the input aliases:
        {alias: [input_id_for_output_0, input_id_for_output_1, ...]}

    Example with 2 proteins x 3 ligands:
        predict_output_ids_with_provenance(
            proteins=(prot_source, "sequences"),
            ligands=(lig_source, "compounds")
        )
        output_ids = ["prot1_lig1", "prot1_lig2", "prot1_lig3",
                      "prot2_lig1", "prot2_lig2", "prot2_lig3"]
        provenance = {
            "proteins": ["prot1", "prot1", "prot1", "prot2", "prot2", "prot2"],
            "ligands": ["lig1", "lig2", "lig3", "lig1", "lig2", "lig3"]
        }
    """
    axes_info, iterated_axes, pure_bundle_axes = _collect_axes_info(
        bundled_name, **named_inputs
    )

    if not axes_info:
        return [bundled_name], {}

    # If no iterated axes, return single bundled name
    if not iterated_axes:
        provenance = {name: ids for name, ids in pure_bundle_axes}
        return [bundled_name], provenance

    # If single iterated axis, return its IDs directly (no join needed)
    if len(iterated_axes) == 1:
        iter_name, iter_ids = iterated_axes[0]
        provenance = {iter_name: list(iter_ids)}
        # Add pure bundle axes: repeat their IDs for every output
        for bundle_name, bundle_ids in pure_bundle_axes:
            if len(bundle_ids) == 1:
                provenance[bundle_name] = list(bundle_ids) * len(iter_ids)
            else:
                provenance[bundle_name] = list(bundle_ids)
        return list(iter_ids), provenance

    # Multiple iterated axes: always full cartesian product, no shortcuts
    # Start with first axis
    output_ids = list(iterated_axes[0][1])
    provenance = {iterated_axes[0][0]: list(iterated_axes[0][1])}

    for axis_name, axis_ids in iterated_axes[1:]:
        new_ids = []
        new_provenance = {k: [] for k in provenance}
        new_provenance[axis_name] = []
        for i, existing_id in enumerate(output_ids):
            for new_id in axis_ids:
                new_ids.append(f"{existing_id}_{new_id}")
                for k in provenance:
                    new_provenance[k].append(provenance[k][i])
                new_provenance[axis_name].append(new_id)
        output_ids = new_ids
        provenance = new_provenance

    # Add pure bundle axes: repeat their IDs for every output
    for bundle_name, bundle_ids in pure_bundle_axes:
        if len(bundle_ids) == 1:
            provenance[bundle_name] = list(bundle_ids) * len(output_ids)
        else:
            provenance[bundle_name] = list(bundle_ids)

    return output_ids, provenance


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
    output_ids, _ = predict_output_ids_with_provenance(bundled_name, **named_inputs)
    return output_ids


def predict_single_output_id(
    bundled_name: str = "bundled_complex",
    **axis_selections: Tuple[str, List[str], Optional[int]]
) -> str:
    """
    Predict the output ID for a single combination of axis selections.

    Each axis_selection is (mode, all_ids, selected_idx_or_None).
    Used by pipe scripts at SLURM time in their iteration loops.

    This mirrors predict_output_ids() logic but for a single output row,
    ensuring pipeline-time and SLURM-time ID generation always agree.

    Always uses full cartesian naming (all iterated axes joined with "_"),
    no shortcuts based on axis length.

    Args:
        bundled_name: Name for fully-bundled case
        **axis_selections: Named tuples of (mode, all_ids, selected_idx)
            where selected_idx is None for bundle mode (no iteration)

    Returns:
        Single output ID string

    Example:
        predict_single_output_id(
            bundled_name="bundled_complex",
            sequences=("each", ["prot1", "prot2"], 0),
            compounds=("each", ["lig1", "lig2", "lig3"], 1)
        )
        → "prot1_lig2"
    """
    # Collect selected IDs from iterated axes (in order)
    parts = []
    for name, (mode, all_ids, idx) in axis_selections.items():
        if mode == "each" or idx is not None:
            selected_idx = idx if idx is not None else 0
            parts.append(all_ids[selected_idx])

    if not parts:
        return bundled_name

    if len(parts) == 1:
        return parts[0]

    return "_".join(parts)


def generate_multiplied_ids(
    input_ids: List[str],
    suffixes: List[str],
    input_stream_name: str = ""
) -> Tuple[List[str], Dict[str, List[str]]]:
    """
    Generate output IDs by appending suffixes to each input ID.

    Args:
        input_ids: Parent IDs to multiply
        suffixes: Suffixes to append (e.g., ["1","2","3"] or ["50A","50V"])
        input_stream_name: Name of the input stream for provenance tracking

    Returns:
        (output_ids, provenance) where provenance maps input_stream_name
        to the parent ID for each output.

    Example:
        generate_multiplied_ids(
            ["prot_1", "prot_2"],
            ["1", "2", "3"],
            input_stream_name="structures"
        )
        -> (["prot_1_1", "prot_1_2", "prot_1_3", "prot_2_1", "prot_2_2", "prot_2_3"],
           {"structures": ["prot_1", "prot_1", "prot_1", "prot_2", "prot_2", "prot_2"]})
    """
    output_ids = []
    parent_ids = []
    for parent_id in input_ids:
        for suffix in suffixes:
            output_ids.append(f"{parent_id}_{suffix}")
            parent_ids.append(parent_id)

    provenance = {}
    if input_stream_name:
        provenance[input_stream_name] = parent_ids

    return output_ids, provenance
