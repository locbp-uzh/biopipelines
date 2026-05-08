# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Simple lazy path descriptor for tools.
"""

import os


class Path:
    """
    Lazy path descriptor - computes path on first access.

    Example:
        class MyTool(BaseConfig):
            # Route artefacts through the canonical layout helpers on
            # BaseConfig — never use ``self.output_folder`` directly.
            main_table = Path(lambda self: self.table_path("results"))
            input_json = Path(lambda self: self.configuration_path("input.json"))
            raw_dump = Path(lambda self: self.execution_path("model_out"))
            structures_map = Path(lambda self: self.stream_map_path("structures"))
            helper_script = Path(lambda self: self.pipe_script_path("pipe_my_script.py"))
    """

    def __init__(self, compute_fn):
        """
        Args:
            compute_fn: Function taking instance and returning path string
        """
        self.compute_fn = compute_fn
        self._attr_name = None

    def __set_name__(self, owner, name):
        self._attr_name = f"_{name}"

    def __get__(self, instance, owner):
        if instance is None:
            return self
        if not hasattr(instance, self._attr_name):
            setattr(instance, self._attr_name, self.compute_fn(instance))
        return getattr(instance, self._attr_name)

    def __set__(self, instance, value):
        setattr(instance, self._attr_name, value)
