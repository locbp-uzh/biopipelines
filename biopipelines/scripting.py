# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Scripting - run a user-authored two-phase script as a typed pipeline step.

An escape hatch for one-off glue logic that doesn't justify a full tool
wrapper. The user writes a script implementing ``configuration(inputs)`` and
``execution(inputs, outputs)`` (see ``scripting_api``); ``Scripting`` wires its
declared streams/tables into the pipeline so downstream tools consume them like
any other tool's output.
"""

import importlib.util
import json
import os
from typing import Any, Dict, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .biopipelines_io import TableReference
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from biopipelines_io import TableReference


class _ConfigInput:
    """Config-time view of one input handed to the user's ``configuration``.

    Exposes the declared ids only — no disk access — so the script can derive
    output ids without reading files that don't exist yet.
    """

    def __init__(self, ids):
        self.ids = list(ids)


class Scripting(BaseConfig):
    """Run a user script as a typed, two-phase pipeline step.

    Args:
        script: Path to a .py file defining ``configuration(inputs)`` and
            ``execution(inputs, outputs)``.
        inputs: Dict mapping a name to a StandardizedOutput, DataStream, or a
            table-column reference (``tool.tables.x.col`` or ``(TableInfo, "col")``).
        env: Conda env the execution phase runs in. Defaults to the
            biopipelines env. The tool does not create it.

    Output:
        Streams/tables are whatever the script's ``configuration`` declares.
    """

    TOOL_NAME = "Scripting"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== Scripting ==="
echo "Uses an existing environment (no additional installation needed)."
touch "$INSTALL_SUCCESS"
echo "=== Scripting ready ==="
"""

    inputs_json = Path(lambda self: self.configuration_path("scripting_inputs.json"))
    outputs_json = Path(lambda self: self.configuration_path("scripting_outputs.json"))
    runner_py = Path(lambda self: self.pipe_script_path("pipe_scripting.py"))

    def __init__(self,
                 script: str,
                 inputs: Dict[str, Any],
                 env: Optional[str] = None,
                 **kwargs):
        self.script = script
        self.inputs_raw = inputs
        self.env = env
        self._user_module = None
        self._declared = None
        super().__init__(**kwargs)

    def validate_params(self):
        if not isinstance(self.script, str) or not self.script:
            raise ValueError("script must be a non-empty path string")
        # A path may carry backslashes on Windows; only reject chars that would
        # break the quoted bash interpolation of the path.
        if any(c in self.script for c in '"`$'):
            raise ValueError(f"script path contains an unsafe character (\" ` $): {self.script}")
        # A bare filename resolves against the configured scripts folder; an
        # existing path as-given always wins.
        self.script = self._resolve_script_path(self.script)
        if not isinstance(self.inputs_raw, dict) or not self.inputs_raw:
            raise ValueError("inputs must be a non-empty dict")
        for name in self.inputs_raw:
            if not isinstance(name, str) or not name:
                raise ValueError(f"input name must be a non-empty string, got {name!r}")
        if self.env is not None:
            _validate_freeform_string("env", self.env)
        self._load_user_module()
        for fn in ("configuration", "execution"):
            if not hasattr(self._user_module, fn):
                raise ValueError(f"script must define a '{fn}' function: {self.script}")

    def _resolve_script_path(self, script: str) -> str:
        """Resolve ``script`` to an existing file, searching the scripts folder.

        Order: (1) the path as given (absolute, or relative to cwd); (2) the
        configured ``folders.infrastructure.scripts`` folder. Raises with both
        candidates listed if neither exists, so a typo names the places looked.
        """
        if os.path.isfile(script):
            return os.path.abspath(script)

        tried = [script]
        try:
            from .config_manager import ConfigManager
        except ImportError:
            from config_manager import ConfigManager
        scripts_folder = ConfigManager().get_scripts_folder()
        if scripts_folder:
            candidate = os.path.join(scripts_folder, script)
            if os.path.isfile(candidate):
                return os.path.abspath(candidate)
            tried.append(candidate)

        raise ValueError(
            f"script not found: {script}. Looked in: "
            + ", ".join(repr(t) for t in tried)
        )

    def _load_user_module(self):
        if self._user_module is not None:
            return
        spec = importlib.util.spec_from_file_location("_bp_scripting_user", self.script)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self._user_module = module

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders
        # Resolve each input to either a serializable stream spec or a table ref.
        self.input_specs = {}
        self.input_streams = {}
        for name, obj in self.inputs_raw.items():
            kind, payload, stream = self._resolve_input(name, obj)
            self.input_specs[name] = {"kind": kind, "payload": payload}
            if stream is not None:
                self.input_streams[name] = stream

    def _resolve_input(self, name, obj):
        """Return (kind, payload, stream_or_None) for one input.

        kind="stream": payload is the input JSON path; stream is the DataStream.
        kind="table":  payload is a TABLE_REFERENCE string; stream is None.
        """
        if isinstance(obj, TableReference):
            return "table", str(obj), None
        if isinstance(obj, tuple) and len(obj) == 2 and isinstance(obj[1], str):
            table_info, column = obj
            path = table_info.info.path if hasattr(table_info, "info") else getattr(table_info, "path", None)
            if path is None:
                raise ValueError(f"input '{name}': cannot resolve table path from {table_info!r}")
            return "table", str(TableReference(path, column)), None
        if isinstance(obj, StandardizedOutput):
            stream = self._pick_stream(name, obj)
        elif isinstance(obj, DataStream):
            stream = obj
        else:
            raise ValueError(
                f"input '{name}' must be a StandardizedOutput, DataStream, or table reference, "
                f"got {type(obj).__name__}"
            )
        return "stream", self.configuration_path(f"input_{name}.json"), stream

    def _pick_stream(self, name, output: StandardizedOutput) -> DataStream:
        """Pick a non-empty stream of a StandardizedOutput input.

        When the input name matches a non-empty stream name, that stream wins — the key disambiguates. Otherwise fall back to the single non-empty stream, erroring only when the choice is genuinely ambiguous."""
        candidates = []
        for stream_name, stream in output.streams.items():
            if stream is not None and len(stream) > 0:
                if stream_name == name:
                    return stream
                candidates.append(stream)
        if not candidates:
            raise ValueError(f"input '{name}': StandardizedOutput has no non-empty stream")
        if len(candidates) > 1:
            raise ValueError(
                f"input '{name}': StandardizedOutput carries multiple streams "
                f"({[c.name for c in candidates]}); name the input after the "
                f"desired stream, or pass an explicit stream "
                f"(e.g. tool.streams.structures)"
            )
        return candidates[0]

    def _declare_outputs(self):
        """Call the user's configuration() with config-time input proxies."""
        if self._declared is not None:
            return self._declared
        self._load_user_module()
        proxies = {}
        for name, spec in self.input_specs.items():
            if spec["kind"] == "stream":
                proxies[name] = _ConfigInput(self.input_streams[name].ids)
            else:
                proxies[name] = _ConfigInput([])  # table inputs expose no ids
        try:
            declared = self._user_module.configuration(proxies)
        except Exception as e:
            raise ValueError(
                f"configuration() failed for script {self.script}: {e}. "
                f"The configuration phase must be pure prediction (no disk reads)."
            ) from e
        if not isinstance(declared, dict):
            raise ValueError(f"configuration() must return a dict, got {type(declared).__name__}")
        self._declared = declared
        return declared

    def _upstream_sources(self):
        """StandardizedOutput inputs — the only ones that can carry a `missing` table."""
        return [obj for obj in self.inputs_raw.values()
                if isinstance(obj, StandardizedOutput)]

    def generate_script(self, script_path: str) -> str:
        declared = self._declare_outputs()

        # Serialize stream inputs to JSON for the runner.
        manifest_inputs = {}
        for name, spec in self.input_specs.items():
            if spec["kind"] == "stream":
                self.input_streams[name].save_json(spec["payload"])
            manifest_inputs[name] = spec
        with open(self.inputs_json, "w") as f:
            json.dump(manifest_inputs, f, indent=2)

        # Serialize the declared output shape + resolved paths for the runner.
        out_manifest = self._build_output_manifest(declared)
        with open(self.outputs_json, "w") as f:
            json.dump(out_manifest, f, indent=2)

        script_content = "#!/bin/bash\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment(name=self.env)
        script_content += f"""
echo "Running Scripting: {os.path.basename(self.script)}"
python "{self.runner_py}" \\
    --script "{self.script}" \\
    --inputs "{self.inputs_json}" \\
    --outputs "{self.outputs_json}"
"""
        # Excuse upstream-filtered ids: merge any upstream `missing` manifests
        # into this tool's own tables/missing.csv so the completion check treats
        # those ids' absent outputs as expected, not failures.
        script_content += self.generate_missing_propagation(*self._upstream_sources())
        script_content += self.generate_completion_check_footer()
        return script_content

    def _build_output_manifest(self, declared):
        from .scripting_api import Stream as _Stream, Table as _Table
        streams = {}
        tables = {}
        for key, val in declared.items():
            if isinstance(val, _Stream):
                streams[key] = {
                    "format": val.format,
                    "map_table": self.stream_map_path(key),
                    "folder": self.stream_folder(key),
                }
            elif isinstance(val, _Table):
                tables[key] = {
                    "columns": val.columns,
                    "path": self.table_path(key),
                }
            else:
                raise ValueError(
                    f"configuration() returned {type(val).__name__} for '{key}'; "
                    f"expected Stream or Table"
                )
        return {"streams": streams, "tables": tables}

    def get_output_files(self) -> Dict[str, Any]:
        from .scripting_api import Stream as _Stream, Table as _Table
        declared = self._declare_outputs()

        result = {
            "structures": DataStream.empty("structures", "pdb"),
            "sequences": DataStream.empty("sequences", "fasta"),
            "compounds": DataStream.empty("compounds", "csv"),
            "tables": {},
            "output_folder": self.output_folder,
        }

        for key, val in declared.items():
            if isinstance(val, _Stream):
                is_value_based = val.format == "csv"
                files = [] if is_value_based else [self.stream_path(key, f"<id>.{val.format}")]
                result[key] = DataStream(
                    name=key,
                    ids=val.ids,
                    files=files,
                    map_table=self.stream_map_path(key),
                    format=val.format,
                )
            elif isinstance(val, _Table):
                if key == "missing":
                    raise ValueError(
                        "configuration() declared a table named 'missing', which is "
                        "reserved by the framework for upstream-filtered id propagation"
                    )
                result["tables"][key] = TableInfo(
                    name=key,
                    path=self.table_path(key),
                    columns=val.columns,
                    description=val.description,
                )
            else:
                raise ValueError(
                    f"configuration() returned {type(val).__name__} for '{key}'; "
                    f"expected Stream or Table"
                )

        # Excuse upstream-filtered ids: declare + own a `missing` table whenever
        # an upstream input carries one, so the completion check treats those
        # ids' absent outputs as expected. generate_missing_propagation fills it.
        if self._collect_upstream_missing_paths(*self._upstream_sources()):
            result["tables"]["missing"] = self.missing_table_info()
        return result
