# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""AiZynthFinder: retrosynthesis route planning over a compounds stream.

Runs a tree search over an expansion policy (template-based by default) to break each target molecule down into purchasable building blocks, and emits the resulting routes and their precursors as first-class streams.

Reference:
    Genheden et al. (2020) AiZynthFinder: a fast, robust and flexible open-source software for retrosynthetic planning. J. Cheminform. 12, 70.
    Saigiridharan et al. (2024) AiZynthFinder 4.0: developments based on learnings from 3 years of industrial application. J. Cheminform. 16, 57.
    https://github.com/MolecularAI/aizynthfinder
"""

import json
import os
import re
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


SEARCH_ALGORITHMS = {
    "mcts": "mcts",
    "retrostar": "aizynthfinder.search.retrostar.search_tree.SearchTree",
    "breadth_first": "aizynthfinder.search.breadth_first.search_tree.SearchTree",
    "dfpn": "aizynthfinder.search.dfpn.search_tree.SearchTree",
}

# aizynthcli scores every route with these, so their values are always in the tree.
DEFAULT_SCORERS = (
    "state score", "number of reactions", "number of pre-cursors",
    "number of pre-cursors in stock",
)

# Scorers upstream builds from the config alone, no model file or other argument.
EXTRA_SCORERS = {
    "fraction in stock": "FractionInStockScorer",
    "average template occurrence": "AverageTemplateOccurrenceScorer",
    "sum of prices": "PriceSumScorer",
    "route cost": "RouteCostScorer",
}

STOP_CRITERIA_KEYS = ("weight", "price", "counts")


def _scorer_column(name: str) -> str:
    """Routes-table column for a scorer name: 'route cost' -> 'route_cost'."""
    return re.sub(r"[^0-9a-z]+", "_", name.lower()).strip("_")


def _is_number(value) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


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

    Each target SMILES is searched with an expansion policy (template-based in the public data); the resulting routes are ranked by score (best first) and their leaf molecules become a compounds stream of purchasable precursors.

    The parameters from search_algorithm onward were added for AiZynthFinder 4.x; each defaults to None, which leaves what the config file says untouched, and is overlaid onto a copy of the config at run time.

    Inputs:
        compounds:  target molecules. SMILES are read from the `smiles` column
                    of the compounds map_table.
        config:     path to an aizynthfinder config.yml, overriding the one
                    written by `download_public_data` at install time. This is
                    also how a user-supplied expansion model, filter model or
                    stock is used: define it in the file's `expansion:`,
                    `filter:` or `stock:` section and select it by name below.
        stocks:     stock names to search against (default: all in the config).
        expansion_policy: expansion policy names (default: the first in the config).
        filter_policy:    filter policy names (default: all in the config).
        max_routes: most routes kept per target (RouteSelectionArguments.nmax).
        min_routes: fewest routes kept per target (RouteSelectionArguments.nmin).
        time_limit: per-target search cutoff in seconds.
        iteration_limit: per-target cap on tree-search iterations.
        nproc:      targets searched in parallel by aizynthcli.
        search_algorithm: "mcts" (upstream default), "retrostar", "breadth_first" or "dfpn". Retro* runs with upstream's default zero molecule cost, so it is guided by the policy priors alone; a learned Retro* cost needs a config file. time_limit is checked between iterations and one breadth_first iteration expands a whole level, so breadth_first overruns it (298 s and 480 s on a 30 s limit).
        max_transforms: deepest route searched, in reaction steps (upstream default 6).
        route_scorer: scorer that ranks the routes and defines `top_score` (default "state score"). Any scorer in `scorers`, "number of reactions", "number of pre-cursors", "number of pre-cursors in stock", or a scorer named in your own config file when `config` is set. Scorers that ignore whether a route is solved can rank an unsolved route first, so filter the routes table on is_solved.
        scorers:    further scorers computed for every returned route, each adding a column to the routes table: "fraction in stock", "average template occurrence", "sum of prices", "route cost". Prices come from the stock and the public stock carries none, so the price scorers use default costs.
        filter_cutoff: feasibility threshold of the quick-filter, in (0, 1]; a reaction it scores below this is rejected (upstream default 0.05).
        frozen_substructures: SMARTS patterns that must survive every reaction; a reaction that splits one of them is rejected. Added as the `frozen_substructure` filter, which is always active. Not supported by "dfpn".
        stock_stop_criteria: extra conditions a stock molecule must meet to count as purchasable, {"weight": max exact mass in g/mol, "price": max price, "counts": {element: max atoms}}.

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
                        and therefore how many route ids exist.
                        top_score is the route_scorer's score of the best route.
                        first_solution_* are 0 when the search never flagged a solution, as upstream reports them (breadth_first can report 0 for a target it solves).
                        A target aizynthcli could not set up, such as an unsanitizable SMILES, keeps its row with everything but id and smiles empty and is listed in missing as a failure)
            routes:     id | compounds.id | route | score | is_solved |
                        number_of_steps | number_of_precursors |
                        number_of_precursors_in_stock | <one column per extra scorer>
                        (score is always the state score.
                        The extra columns are one per name in `scorers`, plus `route_scorer` when it is not one of the four always-computed scorers, spelled with underscores: "route cost" -> route_cost)
            precursors: id | format | smiles | ccd | compounds.id | route |
                        precursor | in_stock
                        (the precursors stream's map_table; filter it with Panda
                        and pass pool= to carry the stream along. A consumer reads the stream
                        called compounds, so hand it `azf.streams.precursors`)
            missing:    id | removed_by | kind | cause
                        (a target with no route is kind="filter"; a target with no SMILES, with a SMILES RDKit cannot sanitize, or that aizynthcli returned no row for is kind="failure")

    Usage:
        with Pipeline(...):
            Resources(memory="16GB", time="4:00:00", cpus=8)
            library = CompoundLibrary("my_library.csv")
            azf = AiZynthFinder(compounds=library, nproc=8)

            buyable = Panda(tables=azf.tables.precursors,
                            operations=[Panda.filter("in_stock == True")],
                            pool=azf)
            RDKit(compounds=buyable.streams.precursors)
    """

    TOOL_NAME = "AiZynthFinder"
    TOOL_VERSION = "2.6"
    ENV_NAME = "aizynthfinder"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        """Create the aizynthfinder env and download the public models.

        `download_public_data` fetches the expansion policies, templates and filter model from Zenodo and the stock from figshare, then writes the config.yml the CLI consumes, so the install is verified on that file's presence as well as on `aizynthcli` starting, which imports the whole package, and on the modules the pipe script imports. An existing env is kept, so a failed download is retried without rebuilding it; force_reinstall rebuilds both.
        """
        env = cls._install_env(env_manager)
        biopipelines = folders.get("biopipelines", "")
        data_dir = folders.get("AiZynthFinder", "")
        config_check = f'[ -f "{data_dir}/config.yml" ]'
        env_check = cls._env_exists_check(env, env_manager)
        env_run = cls._env_run(env, env_manager)
        cli_check = f'{env_run}aizynthcli --help >/dev/null 2>&1 && {env_run}python -c "import pandas, yaml, rdkit" >/dev/null 2>&1'
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_check} && {config_check} \\
   && {cli_check}; then
    echo "AiZynthFinder already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block(env, env_manager) if force_reinstall else ""
        env_block = cls._env_install_block(env, env_manager, biopipelines)
        return f"""echo "=== Installing AiZynthFinder ==="
{skip}{remove_block}
if {env_check}; then
    echo "{env} environment already present, not recreating it"
else
{env_block}
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create {env} environment."
    exit 1
fi
fi

mkdir -p "{data_dir}"
if {config_check}; then
    echo "Public data already present in {data_dir}"
else
    echo "Downloading AiZynthFinder public models and templates"
    (cd "{data_dir}" && {env_run}download_public_data .)
fi

# Verify: aizynthcli starts (imports the package) AND the downloaded config is present.
if {config_check} \\
   && {cli_check}; then
    touch "$INSTALL_SUCCESS"
    echo "=== AiZynthFinder installation complete ==="
else
    echo "ERROR: AiZynthFinder verification failed (aizynthcli does not start, pandas, yaml or rdkit does not import, or config.yml missing). Retry with force_reinstall=True to rebuild the env."
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
                 search_algorithm: Optional[str] = None,
                 max_transforms: Optional[int] = None,
                 route_scorer: Optional[str] = None,
                 scorers: Optional[List[str]] = None,
                 filter_cutoff: Optional[float] = None,
                 frozen_substructures: Optional[List[str]] = None,
                 stock_stop_criteria: Optional[Dict[str, Any]] = None,
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
        self.search_algorithm = search_algorithm
        self.max_transforms = max_transforms
        self.route_scorer = route_scorer
        self.scorers = _as_name_list("scorers", scorers)
        self.filter_cutoff = filter_cutoff
        self.frozen_substructures = _as_name_list("frozen_substructures", frozen_substructures)
        self.stock_stop_criteria = stock_stop_criteria
        # Scorers with a routes-table column: the requested ones plus a ranking scorer aizynthcli does not always compute.
        self.extra_scorer_names = list(self.scorers or [])
        if (isinstance(route_scorer, str) and route_scorer not in DEFAULT_SCORERS
                and route_scorer not in self.extra_scorer_names):
            self.extra_scorer_names.append(route_scorer)
        super().__init__(**kwargs)

    def _validate_scorers(self):
        known = DEFAULT_SCORERS + tuple(EXTRA_SCORERS)
        # A scorer defined in the user's own config file cannot be checked here.
        custom_ok = self.config is not None

        if self.route_scorer is not None:
            if not isinstance(self.route_scorer, str) or not self.route_scorer.strip():
                raise ValueError(f"route_scorer must be a scorer name, got {self.route_scorer!r}")
            _validate_freeform_string("route_scorer", self.route_scorer)
            if self.route_scorer not in known and not custom_ok:
                raise ValueError(
                    f"route_scorer {self.route_scorer!r} is not one of {list(known)}. A scorer "
                    f"defined in your own config file needs config= to be set as well")

        for i, name in enumerate(self.scorers or []):
            if not isinstance(name, str) or not name.strip():
                raise ValueError(f"scorers[{i}] must be a scorer name, got {name!r}")
            _validate_freeform_string(f"scorers[{i}]", name)
            if name in DEFAULT_SCORERS:
                raise ValueError(
                    f"scorers[{i}]={name!r} is always computed; its value is already in the "
                    f"routes table (score, number_of_steps, number_of_precursors, "
                    f"number_of_precursors_in_stock)")
            if name not in EXTRA_SCORERS and not custom_ok:
                raise ValueError(
                    f"scorers[{i}]={name!r} is not one of {list(EXTRA_SCORERS)}. A scorer "
                    f"defined in your own config file needs config= to be set as well")

        names = self.scorers or []
        if len(set(names)) != len(names):
            raise ValueError(f"scorers lists a scorer more than once: {names}")
        columns = [_scorer_column(n) for n in self.extra_scorer_names]
        for name, column in zip(self.extra_scorer_names, columns):
            if not column or column in ROUTE_COLUMNS:
                raise ValueError(
                    f"scorer {name!r} would write the routes column {column!r}, which "
                    f"{'is empty' if not column else 'already exists'}")
        if len(set(columns)) != len(columns):
            raise ValueError(f"scorer names {self.extra_scorer_names} map to colliding column names")

    def _validate_stop_criteria(self):
        criteria = self.stock_stop_criteria
        if criteria is None:
            return
        if not isinstance(criteria, dict) or not criteria:
            raise ValueError(
                f"stock_stop_criteria must be a non-empty dict with keys from "
                f"{list(STOP_CRITERIA_KEYS)}, got {criteria!r}")
        for key, value in criteria.items():
            if key not in STOP_CRITERIA_KEYS:
                raise ValueError(
                    f"stock_stop_criteria key {key!r} is not one of {list(STOP_CRITERIA_KEYS)}")
            if key == "counts":
                if not isinstance(value, dict) or not value:
                    raise ValueError("stock_stop_criteria['counts'] must be a non-empty dict {element: max atoms}")
                for element, count in value.items():
                    if not isinstance(element, str) or not element.isalpha():
                        raise ValueError(f"stock_stop_criteria['counts'] key {element!r} must be an element symbol")
                    if not isinstance(count, int) or isinstance(count, bool) or count < 0:
                        raise ValueError(
                            f"stock_stop_criteria['counts'][{element!r}] must be an integer >= 0, got {count!r}")
            elif not _is_number(value) or value <= 0:
                raise ValueError(f"stock_stop_criteria[{key!r}] must be a positive number, got {value!r}")

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

        if self.search_algorithm is not None and self.search_algorithm not in SEARCH_ALGORITHMS:
            raise ValueError(
                f"search_algorithm must be one of {list(SEARCH_ALGORITHMS)}, got {self.search_algorithm!r}")
        if self.max_transforms is not None and (
                not isinstance(self.max_transforms, int) or isinstance(self.max_transforms, bool)
                or self.max_transforms < 1):
            raise ValueError(f"max_transforms must be an integer >= 1, got {self.max_transforms!r}")

        self._validate_scorers()

        if self.filter_cutoff is not None and (
                not _is_number(self.filter_cutoff) or not 0 < self.filter_cutoff <= 1):
            raise ValueError(f"filter_cutoff must be in (0, 1], got {self.filter_cutoff!r}")
        if self.frozen_substructures is not None:
            for i, smarts in enumerate(self.frozen_substructures):
                if not isinstance(smarts, str) or not smarts.strip():
                    raise ValueError(f"frozen_substructures[{i}] must be a SMARTS string, got {smarts!r}")
        if self.search_algorithm == "dfpn" and (
                self.filter_policy or self.filter_cutoff is not None or self.frozen_substructures):
            raise ValueError(
                "search_algorithm='dfpn' does not support a filter policy upstream; drop "
                "filter_policy, filter_cutoff and frozen_substructures, or pick another algorithm")

        self._validate_stop_criteria()

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
        if self.search_algorithm is not None:
            lines.append(f"SEARCH ALGORITHM: {self.search_algorithm}")
        if self.max_transforms is not None:
            lines.append(f"MAX TRANSFORMS: {self.max_transforms}")
        if self.route_scorer is not None:
            lines.append(f"ROUTE SCORER: {self.route_scorer}")
        if self.extra_scorer_names:
            lines.append(f"EXTRA SCORERS: {', '.join(self.extra_scorer_names)}")
        if self.filter_cutoff is not None:
            lines.append(f"FILTER CUTOFF: {self.filter_cutoff}")
        if self.frozen_substructures:
            lines.append(f"FROZEN SUBSTRUCTURES: {len(self.frozen_substructures)} SMARTS")
        if self.stock_stop_criteria:
            lines.append(f"STOCK STOP CRITERIA: {self.stock_stop_criteria}")
        return lines

    def _pipe_config(self) -> Dict[str, Any]:
        """The JSON the pipe script reads; the new options reach it already resolved."""
        # The installed config.yml is authoritative; the pipe overlays only the options the user set.
        base_config = self.config or os.path.join(self.folders.get("AiZynthFinder", ""), "config.yml")

        return {
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
            "search_algorithm": SEARCH_ALGORITHMS.get(self.search_algorithm),
            "max_transforms": self.max_transforms,
            "route_scorer": self.route_scorer,
            "extra_scorers": [
                {"name": name, "column": _scorer_column(name), "class": EXTRA_SCORERS.get(name)}
                for name in self.extra_scorer_names
            ],
            "filter_cutoff": self.filter_cutoff,
            "frozen_substructures": self.frozen_substructures,
            "stock_stop_criteria": self.stock_stop_criteria,
        }

    def generate_script(self, script_path: str) -> str:
        self.compounds_stream.save_json(self.compounds_json)

        cfg = self._pipe_config()
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
            self.compounds_stream.ids, "[<#>]", input_stream_name="compounds"
        )
        precursor_ids = generate_multiplied_ids_pattern(
            self.compounds_stream.ids, "[<#>_<#>]", input_stream_name="compounds"
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
                columns=ROUTE_COLUMNS + [_scorer_column(n) for n in self.extra_scorer_names],
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
