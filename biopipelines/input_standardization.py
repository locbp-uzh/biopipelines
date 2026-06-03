# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Shared resolver that normalizes convenience inputs to canonical streams.

Many tools accept the same shape of input: a StandardizedOutput, a raw
DataStream, or an ergonomic shorthand (e.g. ``ligand="LIG"`` instead of
``ligand=Ligand(code="LIG")``). ``resolve_basic_input`` collapses the three
into the canonical stream so the rest of the framework only ever sees a
DataStream, preserving the entity-only invariant.
"""

from typing import Optional, Type

try:
    from .datastream import DataStream
    from .base_config import StandardizedOutput
except ImportError:
    from datastream import DataStream
    from base_config import StandardizedOutput


def resolve_basic_input(obj, cls: Type, stream: str, argument: str, *,
                        allow_none: bool = True) -> Optional[DataStream]:
    """Normalize a tool input to a named DataStream.

    - StandardizedOutput -> its ``streams.<stream>``.
    - DataStream         -> returned as-is.
    - shorthand value    -> ``cls(**{argument: obj}, _internal=True)`` then its
                            ``streams.<stream>``. Inside a Pipeline the entity
                            auto-registers as an internal tool, so it runs
                            before the consuming tool.

    Anything else (e.g. a raw tool instance) raises. ``None`` passes through
    when ``allow_none``.
    """
    if obj is None:
        if allow_none:
            return None
        raise ValueError(f"{argument} input is required")
    if isinstance(obj, StandardizedOutput):
        return getattr(obj.streams, stream)
    if isinstance(obj, DataStream):
        return obj
    if isinstance(obj, str):
        entity = cls(**{argument: obj}, _internal=True)
        # In a pipeline the entity auto-registers and returns a StandardizedOutput;
        # standalone it returns the raw tool instance.
        if isinstance(entity, StandardizedOutput):
            return getattr(entity.streams, stream)
        return getattr(StandardizedOutput(entity.get_output_files()).streams, stream)
    raise ValueError(
        f"input must be a string, DataStream, or StandardizedOutput, "
        f"got {type(obj).__name__}"
    )
