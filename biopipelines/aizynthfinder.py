# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""AiZynthFinder: retrosynthesis route planning over a compounds stream.

Runs a Monte-Carlo tree search over a template-based expansion policy to break each target molecule down into purchasable building blocks, and emits the resulting routes and their precursors as first-class streams.

Reference:
    Genheden et al. (2020) AiZynthFinder: a fast, robust and flexible open-source software for retrosynthetic planning. J. Cheminform. 12, 70.
    https://github.com/MolecularAI/aizynthfinder
"""

import json
import os
from typing import Dict, List, Any, Union, Optional

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .combinatorics import generate_multiplied_ids_pattern
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from combinatorics import generate_multiplied_ids_pattern


TARGET_COLUMNS = [
    "id", "smiles", "is_solved", "top_score", "number_of_routes",
    "number_of_returned_routes", "number_of_solved_routes", "number_of_nodes",
    "search_time", "first_solution_time", "first_solution_iteration",
]

ROUTE_COLUMNS = [
    "id", "compounds.id", "route", "score", "is_solved", "number_of_steps",
    "number_of_precursors", "number_of_precursors_in_stock",
]

# Canonical compounds schema first (id | format | smiles | ccd) so this stream
# feeds any compounds consumer, then the route provenance this tool adds.
PRECURSOR_COLUMNS = [
    "id", "format", "smiles", "ccd", "compounds.id", "route", "precursor", "in_stock",
]


def _as_name_list(param: str, value) -> Optional[List[str]]:
    """Normalize a list-valued name parameter, rejecting a bare string."""
    if value is None:
        return None
    if isinstance(value, str):
        raise ValueError(
            f"{param} must be a list of names, not a string; "
            f"pass {param}=[{value!r}] for a single entry")
    if not isinstance(value, (list, tuple)):
        raise ValueError(
            f"{param} must be a list or tuple of names, got {type(value).__name__}")
    return list(value) or None


class AiZynthFinder(BaseConfig):
    """
    AiZynthFinder: retrosynthetic route search for a compounds stream.

    Each target SMILES is searched with a template-based expansion policy; the
    resulting routes are ranked by score (best first) and their leaf molecules
    become a compounds stream of purchasable precursors.

    Inputs:
        compounds:  target molecules. SMILES are read from the `smiles` column
                    of the compounds map_table.
        config:     path to an aizynthfinder config.yml, overriding the one
                    written by `download_public_data` at install time.
        stocks:     stock names to search against (default: all in the config).
        expansion_policy: expansion policy names (default: the first in the config).
        filter_policy:    filter policy names (default: all in the config).
        max_routes: most routes kept per target (RouteSelectionArguments.nmax).
        min_routes: fewest routes kept per target (RouteSelectionArguments.nmin).
        time_limit: per-target search cutoff in seconds.
        iteration_limit: per-target cap on tree-search iterations.
        nproc:      targets searched in parallel by aizynthcli.

    Outputs:
        Streams:
            routes:     one <id>.json per route, holding that route's reaction
                        tree. ids are <compound>_<route>.
            precursors: value-based compounds stream of route leaf molecules.
                        ids are <compound>_<route>_<precursor>.
        Tables:
            targets:    id | smiles | is_solved | top_score | number_of_routes |
                        number_of_returned_routes | number_of_solved_routes |
                        number_of_nodes | search_time | first_solution_time |
                        first_solution_iteration
                        (one row per input compound; unsolved targets included.
                        number_of_routes is what the search found;
                        number_of_returned_routes is what max_routes let through
                        and therefore how many route ids exist)
            routes:     id | compounds.id | route | score | is_solved |
                        number_of_steps | number_of_precursors |
                        number_of_precursors_in_stock
            precursors: id | format | smiles | ccd | compounds.id | route |
                        precursor | in_stock
                        (the precursors stream's map_table; filter it with Panda
                        and pass pool= to carry the stream along)
            missing:    id | removed_by | kind | cause

    Usage:
        with Pipeline(...):
            Resources(memory="16GB", time="4:00:00", cpus=8)
            library = CompoundLibrary("my_library.csv")
            azf = AiZynthFinder(compounds=library, nproc=8)

            buyable = Panda(tables=azf.tables.precursors,
                            operations=[Panda.filter("in_stock == True")],
                            pool=azf)
    """

    TOOL_NAME = "AiZynthFinder"
    TOOL_VERSION = "1.2"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        """Create the aizynthfinder env and download the public models.

        `download_public_data` fetches the expansion policy, templates and stock
        from figshare and writes the config.yml the CLI consumes, so the install
        is verified on that file's presence as well as on the import.
        """
        biopipelines = folders.get("biopipelines", "")
        data_dir = folders.get("AiZynthFinder", "")
        config_check = f'[ -f "{data_dir}/config.yml" ]'
        env_check = cls._env_exists_check("aizynthfinder", env_manager)
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_check} && {config_check} \\
   && {cls._env_run("aizynthfinder", env_manager)}python -c "import aizynthfinder" >/dev/null 2>&1; then
    echo "AiZynthFinder already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block("aizynthfinder", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("aizynthfinder", env_manager, biopipelines)
        return f"""echo "=== Installing AiZynthFinder ==="
{skip}{remove_block}
{env_block}
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create aizynthfinder environment."
    exit 1
fi

mkdir -p "{data_dir}"
if {config_check}; then
    echo "Public data already present in {data_dir}"
else
    echo "Downloading AiZynthFinder public models and templates"
    (cd "{data_dir}" && {cls._env_run("aizynthfinder", env_manager)}download_public_data .)
fi

# Verify: package importable AND the downloaded config the CLI needs is present.
if {config_check} \\
   && {cls._env_run("aizynthfinder", env_manager)}python -c "import aizynthfinder" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== AiZynthFinder installation complete ==="
else
    echo "ERROR: AiZynthFinder verification failed (cannot import aizynthfinder or config.yml missing)"
    exit 1
fi
"""

    config_json = Path(lambda self: self.configuration_path("aizynthfinder_config.json"))
    compounds_json = Path(lambda self: self.configuration_path("compounds_ds.json"))
    search_config_yaml = Path(lambda self: self.configuration_path("search_config.yml"))
    targets_smi = Path(lambda self: self.configuration_path("targets.smi"))
    cli_output = Path(lambda self: self.execution_path("aizynth_output.json.gz"))
    local_missing_csv = Path(lambda self: self.execution_path("local_missing.csv"))
    routes_map = Path(lambda self: self.stream_map_path("routes"))
    # Content-bearing: the precursors map_table IS the precursors table.
    precursors_csv = Path(lambda self: self.stream_path("precursors", "precursors.csv"))
    targets_csv = Path(lambda self: self.table_path("targets"))
    routes_csv = Path(lambda self: self.table_path("routes"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_aizynthfinder.py"))

    def __init__(self,
                 compounds: Union[DataStream, StandardizedOutput],
                 config: Optional[str] = None,
                 stocks: Optional[List[str]] = None,
                 expansion_policy: Optional[List[str]] = None,
                 filter_policy: Optional[List[str]] = None,
                 max_routes: int = 25,
                 min_routes: int = 5,
                 time_limit: Optional[int] = None,
                 iteration_limit: Optional[int] = None,
                 nproc: Optional[int] = None,
                 **kwargs):
        self.compounds = compounds  # raw handle kept for missing-propagation
        if isinstance(compounds, StandardizedOutput):
            self.compounds_stream: DataStream = compounds.streams.compounds
        elif isinstance(compounds, DataStream):
            self.compounds_stream = compounds
        else:
            raise ValueError(
                f"compounds must be DataStream or StandardizedOutput, got {type(compounds)}"
            )

        self.config = config
        # A bare string would list() into its characters and reach aizynthcli as
        # "--stocks s t o c k", failing deep inside the CLI instead of here.
        self.stocks = _as_name_list("stocks", stocks)
        self.expansion_policy = _as_name_list("expansion_policy", expansion_policy)
        self.filter_policy = _as_name_list("filter_policy", filter_policy)
        self.max_routes = max_routes
        self.min_routes = min_routes
        self.time_limit = time_limit
        self.iteration_limit = iteration_limit
        self.nproc = nproc
        super().__init__(**kwargs)

    def validate_params(self):
        if not self.compounds_stream or len(self.compounds_stream) == 0:
            raise ValueError("compounds parameter is required and must not be empty")

        if self.config is not None:
            _validate_freeform_string("config", self.config)

        for name, values in (("stocks", self.stocks),
                             ("expansion_policy", self.expansion_policy),
                             ("filter_policy", self.filter_policy)):
            if values is None:
                continue
            if not values:
                raise ValueError(f"{name} must be a non-empty list, or None to use the config default")
            for i, value in enumerate(values):
                if not isinstance(value, str):
                    raise ValueError(f"{name}[{i}] must be a string, got {type(value).__name__}")
                _validate_freeform_string(f"{name}[{i}]", value)

        if self.max_routes < 1:
            raise ValueError(f"max_routes must be >= 1, got {self.max_routes}")
        if self.min_routes < 1:
            raise ValueError(f"min_routes must be >= 1, got {self.min_routes}")
        if self.min_routes > self.max_routes:
            raise ValueError(
                f"min_routes ({self.min_routes}) cannot exceed max_routes ({self.max_routes})"
            )
        if self.time_limit is not None and self.time_limit <= 0:
            raise ValueError(f"time_limit must be positive, got {self.time_limit}")
        if self.iteration_limit is not None and self.iteration_limit <= 0:
            raise ValueError(f"iteration_limit must be positive, got {self.iteration_limit}")
        if self.nproc is not None and self.nproc < 1:
            raise ValueError(f"nproc must be >= 1, got {self.nproc}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"TARGETS: {len(self.compounds_stream)} molecules")
        lines.append(f"ROUTES PER TARGET: {self.min_routes}-{self.max_routes}")
        if self.stocks:
            lines.append(f"STOCKS: {', '.join(self.stocks)}")
        if self.expansion_policy:
            lines.append(f"EXPANSION POLICY: {', '.join(self.expansion_policy)}")
        if self.filter_policy:
            lines.append(f"FILTER POLICY: {', '.join(self.filter_policy)}")
        if self.time_limit is not None:
            lines.append(f"TIME LIMIT: {self.time_limit}s per target")
        if self.iteration_limit is not None:
            lines.append(f"ITERATION LIMIT: {self.iteration_limit}")
        return lines

    def generate_script(self, script_path: str) -> str:
        self.compounds_stream.save_json(self.compounds_json)

        # The installed config.yml is authoritative; the pipe copies it and
        # overlays only the search limits, so the models stay where they were
        # downloaded.
        base_config = self.config or os.path.join(self.folders.get("AiZynthFinder", ""), "config.yml")

        cfg = {
            "compounds_json": self.compounds_json,
            "base_config": base_config,
            "search_config_yaml": self.search_config_yaml,
            "targets_smi": self.targets_smi,
            "cli_output": self.cli_output,
            "routes_dir": self.stream_folder("routes"),
            "routes_map_csv": self.routes_map,
            "precursors_csv": self.precursors_csv,
            "targets_csv": self.targets_csv,
            "routes_csv": self.routes_csv,
            "local_missing_csv": self.local_missing_csv,
            "stocks": self.stocks,
            "expansion_policy": self.expansion_policy,
            "filter_policy": self.filter_policy,
            "max_routes": self.max_routes,
            "min_routes": self.min_routes,
            "time_limit": self.time_limit,
            "iteration_limit": self.iteration_limit,
            "nproc": self.nproc,
        }
        os.makedirs(os.path.dirname(self.config_json), exist_ok=True)
        with open(self.config_json, "w", encoding="utf-8") as f:
            json.dump(cfg, f, indent=2)

        script = "#!/bin/bash\n"
        script += "# AiZynthFinder retrosynthesis script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        # The helper imports biopipelines and stays on the host; only the
        # aizynthcli subprocess it launches enters the container, if configured.
        script += f"""echo "Running AiZynthFinder on {len(self.compounds_stream)} target(s)"
python "{self.helper_py}" \\
    --config "{self.config_json}" \\
    --container-prefix "{self.container_prefix()}"
"""
        script += self.generate_missing_propagation(
            self.compounds, local_missing=self.local_missing_csv, missing_csv=self.missing_csv
        )
        script += self.generate_completion_check_footer()
        return script

    def get_output_files(self) -> Dict[str, Any]:
        # Route and precursor counts are only known once the search runs, so both
        # suffixes stay in one lazy bracket.
        route_ids = generate_multiplied_ids_pattern(
            self.compounds_stream.ids, "[<N>]", input_stream_name="compounds"
        )
        precursor_ids = generate_multiplied_ids_pattern(
            self.compounds_stream.ids, "[<N>_<N>]", input_stream_name="compounds"
        )

        routes_stream = DataStream(
            name="routes",
            ids=route_ids,
            files=[self.stream_path("routes", "<id>.json")],
            map_table=self.routes_map,
            format="json",
        )
        precursors_stream = DataStream(
            name="precursors",
            ids=precursor_ids,
            files=[],
            map_table=self.precursors_csv,
            format="csv",
        )

        tables = {
            "targets": TableInfo(
                name="targets",
                path=self.targets_csv,
                columns=list(TARGET_COLUMNS),
                description="Per-target retrosynthesis search summary",
            ),
            "routes": TableInfo(
                name="routes",
                path=self.routes_csv,
                columns=list(ROUTE_COLUMNS),
                description="Per-route score and step count",
            ),
            "precursors": TableInfo(
                name="precursors",
                path=self.precursors_csv,
                columns=list(PRECURSOR_COLUMNS),
                description="Route leaf molecules (the precursors stream's content table)",
            ),
            "missing": self.missing_table_info(self.missing_csv),
        }

        return {
            "routes": routes_stream,
            "precursors": precursors_stream,
            "tables": tables,
            "output_folder": self.output_folder,
        }
