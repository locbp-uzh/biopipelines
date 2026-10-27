# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
RFD4Proteina: all-atom flow-matching design with RFD4-Proteina (Baker Lab / IPD x NVIDIA GenAIR).

Sequence-structure co-design from a YAML condition spec: the wrapper turns each input structure plus a contig and a few conditions into the spec `rfd4-generate` reads, runs the model, and converts the designed annotated CIFs to PDB with their sequences and per-design tables.

Code: https://github.com/RosettaCommons/RFD4-Proteina (BSD-3-Clause). Docs: https://rfd4-proteina.org/latest/
"""

import copy
import os
import re
from typing import Any, Dict, List, NamedTuple, Optional, Union

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


TASKS = ("protein_binder", "small_molecule_binder", "motif_scaffold")
THINKING_LEVELS = ("low", "medium", "high")

# Upstream's task-aware recycling table (docs/how-to/run-the-model.md): PPI is front-loaded, small-molecule binders and M-CSA are midpeak, and the cookbook runs enzyme and motif scaffolding as medium/midpeak.
TASK_RECYCLE_SHAPE = {
    "protein_binder": "front_loaded",
    "small_molecule_binder": "midpeak",
    "motif_scaffold": "midpeak",
}

# Hydra keys the wrapper sets itself; forwarding one would conflict with it or break the id mapping to the produced files.
MANAGED_OVERRIDES = (
    "inputs", "out_dir", "num_replicates", "batch_size", "seed", "recycle_strategy",
    "hf_ckpt_path", "hf_ckpt_revision", "metrics_tags", "metrics", "compress_outputs",
    "global_prefix", "job_id", "gen_njobs", "subset_to_examples", "hydra.run.dir",
)

_PRESERVED_RE = re.compile(r"^([A-Za-z]+)(\d+)(?:-(\d+))?$")
_LENGTH_RE = re.compile(r"^(\d+)(?:-(\d+))?([PRD])?$")
_HOTSPOT_RESIDUE_RE = re.compile(r"^([A-Za-z]+)(\d+)$")
_ATOM_NAME_RE = re.compile(r"^[A-Za-z0-9]+$")
_TAG_RE = re.compile(r"^[A-Za-z0-9_:.-]+$")


class ContigToken(NamedTuple):
    kind: str  # "chainbreak", "preserved" or "expandable"
    chain: str = ""
    start: int = 0
    end: int = 0
    name: str = ""
    seg_min: int = 0
    seg_max: int = 0


def _parse_token(token: str) -> ContigToken:
    """One contig token, classified by shape the way upstream's contig_spec does."""
    token = token.strip()
    if not token:
        raise ValueError("empty contig token (check for stray or trailing commas)")
    if token == "/0":
        return ContigToken("chainbreak")
    if ":" in token:
        name, _, body = token.partition(":")
        if not name:
            raise ValueError(f"malformed named expandable segment (empty name): {token!r}")
        if name == "MASKED":
            raise ValueError("an expandable segment cannot be named 'MASKED'; upstream reserves it")
        match = _LENGTH_RE.match(body)
        if not match:
            raise ValueError(f"malformed expandable segment: {token!r}")
        return _expandable(match, name, token)
    if token[0].isalpha():
        match = _PRESERVED_RE.match(token)
        if not match:
            raise ValueError(f"malformed preserved segment: {token!r}")
        start = int(match.group(2))
        end = int(match.group(3)) if match.group(3) is not None else start
        if end < start:
            raise ValueError(f"preserved segment {token!r} has end < start")
        return ContigToken("preserved", chain=match.group(1), start=start, end=end)
    match = _LENGTH_RE.match(token)
    if not match:
        raise ValueError(f"malformed expandable segment: {token!r}")
    return _expandable(match, "", token)


def _expandable(match, name: str, token: str) -> ContigToken:
    seg_min = int(match.group(1))
    seg_max = int(match.group(2)) if match.group(2) is not None else seg_min
    if seg_max < seg_min:
        raise ValueError(f"expandable segment {token!r} has end < start")
    return ContigToken("expandable", name=name, seg_min=seg_min, seg_max=seg_max)


def contig_runs(tokens: List[ContigToken]) -> List[List[ContigToken]]:
    """The non-empty token runs between chain breaks; each run becomes one output chain."""
    runs: List[List[ContigToken]] = [[]]
    for token in tokens:
        if token.kind == "chainbreak":
            runs.append([])
        else:
            runs[-1].append(token)
    return [run for run in runs if run]


def parse_contig(contig: str) -> List[ContigToken]:
    """Parse and validate a contig string against the rules upstream enforces when it builds the structure.

    Token grammar and the checks below follow `rfproteina/inference_utils/contig_spec.py`: `/0` chain break, `A42-48` preserved segment, `E1:2-20` named and `20-30` unnamed expandable segment, a preserved segment required, no overlap on one chain, one input chain per run, unique names and unique output chains.
    """
    if not isinstance(contig, str) or not contig.strip():
        raise ValueError("contig is required, e.g. '50-120, /0, A17-131'")
    tokens = [_parse_token(piece) for piece in contig.split(",")]

    names = [t.name for t in tokens if t.kind == "expandable" and t.name]
    duplicates = sorted({n for n in names if names.count(n) > 1})
    if duplicates:
        raise ValueError(f"duplicate expandable-segment name(s) in contig: {duplicates}")

    runs = contig_runs(tokens)
    if not runs:
        raise ValueError("contig has no chain runs (only chain breaks)")

    by_chain: Dict[str, List[ContigToken]] = {}
    for token in tokens:
        if token.kind == "preserved":
            by_chain.setdefault(token.chain, []).append(token)
    for chain, segments in by_chain.items():
        segments = sorted(segments, key=lambda t: t.start)
        for previous, current in zip(segments, segments[1:]):
            if current.start <= previous.end:
                raise ValueError(
                    f"preserved segments on chain {chain!r} overlap: "
                    f"{chain}{previous.start}-{previous.end} and {chain}{current.start}-{current.end}")

    seen_chains = set()
    for run in runs:
        chains = {t.chain for t in run if t.kind == "preserved"}
        if len(chains) > 1:
            raise ValueError(
                f"one output chain mixes preserved segments of chains {sorted(chains)}; "
                "insert '/0' between them")
        for chain in chains:
            if chain in seen_chains:
                raise ValueError(f"contig produces two separate output chains both named {chain!r}")
            seen_chains.add(chain)

    if not by_chain:
        raise ValueError(
            "contig has no preserved segment: the spec's input structure exists only to be referenced "
            "by one (e.g. 'A17-131'), so upstream rejects a contig made only of expandable segments")
    return tokens


def preserved_chains(tokens: List[ContigToken]) -> List[str]:
    """Chains of the preserved segments, in order of first appearance."""
    chains: List[str] = []
    for token in tokens:
        if token.kind == "preserved" and token.chain not in chains:
            chains.append(token.chain)
    return chains


def _residue_select(chain: str, start: int, end: int) -> str:
    if start == end:
        return f"sel('{chain}/*/{start}')"
    return f"sel('{chain}/*/[{start}-{end}]')"


def hotspot_selects(hotspots: Union[str, Dict[str, str], None]) -> List[str]:
    """Selection strings for `C_HOT` from "A56:CG,OH;A115:CG,SD", "A67,A89" or {"A56": "CG,OH"}."""
    if hotspots is None or hotspots == "" or hotspots == {}:
        return []
    pairs: List[tuple] = []
    if isinstance(hotspots, dict):
        pairs = [(str(k).strip(), "" if v is None else str(v)) for k, v in hotspots.items()]
    elif isinstance(hotspots, str):
        for group in hotspots.split(";"):
            group = group.strip()
            if not group:
                continue
            if ":" in group:
                residue, atoms = group.split(":", 1)
                pairs.append((residue.strip(), atoms.strip()))
            else:
                pairs.extend((residue.strip(), "") for residue in group.split(",") if residue.strip())
    else:
        raise ValueError("hotspots must be a string like 'A56:CG,OH;A115' or a dict {'A56': 'CG,OH'}")

    selects = []
    for residue, atoms in pairs:
        match = _HOTSPOT_RESIDUE_RE.match(residue)
        if not match:
            raise ValueError(f"hotspot {residue!r} must be <chain><residue number>, e.g. 'A56'")
        chain, number = match.group(1), int(match.group(2))
        names = [a.strip() for a in atoms.split(",") if a.strip()]
        for name in names:
            if not _ATOM_NAME_RE.match(name):
                raise ValueError(f"hotspot atom name {name!r} must be alphanumeric")
        if not names:
            selects.append(f"sel('{chain}/*/{number}')")
        elif len(names) == 1:
            selects.append(f"sel('{chain}/*/{number}/{names[0]}')")
        else:
            selects.append(f"sel('{chain}/*/{number}/[{','.join(names)}]')")
    return selects


def normalize_conditions(conditions: Optional[Dict[Any, Any]]) -> Dict[Any, Any]:
    """Check the shape of an upstream `conditions:` mapping: annotation -> value -> [ {select, group?} ]."""
    if conditions is None:
        return {}
    if not isinstance(conditions, dict):
        raise ValueError("conditions must be a mapping {annotation: {value: [{'select': ...}]}}")
    normalized: Dict[Any, Any] = {}
    for annotation, values in conditions.items():
        if not isinstance(annotation, str) or not annotation:
            raise ValueError(f"condition annotation names must be non-empty strings, got {annotation!r}")
        if not isinstance(values, dict) or not values:
            raise ValueError(f"conditions[{annotation!r}] must map each value to a list of selections")
        normalized[annotation] = {}
        for value, entries in values.items():
            if value is not None and not isinstance(value, (bool, int, float, str)):
                raise ValueError(f"conditions[{annotation!r}] has a value key of type {type(value).__name__}")
            if not isinstance(entries, list) or not entries:
                raise ValueError(
                    f"conditions[{annotation!r}][{value!r}] must be a list of {{'select': ...}} entries, "
                    "even for one selection")
            for entry in entries:
                if not isinstance(entry, dict) or not isinstance(entry.get("select"), str) \
                        or not entry["select"].strip():
                    raise ValueError(
                        f"conditions[{annotation!r}][{value!r}] entries need a non-empty 'select' string")
                extra = set(entry) - {"select", "group"}
                if extra:
                    raise ValueError(f"conditions[{annotation!r}][{value!r}] has unknown key(s) {sorted(extra)}")
            normalized[annotation][value] = copy.deepcopy(entries)
    return normalized


def merge_conditions(base: Dict[Any, Any], extra: Dict[Any, Any]) -> Dict[Any, Any]:
    """Append `extra` after `base`; upstream resolves overlaps last write wins, in declaration order."""
    merged = copy.deepcopy(base)
    for annotation, values in extra.items():
        slot = merged.setdefault(annotation, {})
        for value, entries in values.items():
            slot.setdefault(value, []).extend(copy.deepcopy(entries))
    return merged


def build_conditions(task: str, tokens: List[ContigToken],
                     hotspots: Union[str, Dict[str, str], None] = None,
                     conditions: Optional[Dict[Any, Any]] = None) -> Dict[Any, Any]:
    """The `conditions:` block for one task, then the hotspots, then the user's raw conditions.

    protein_binder: the upstream PPI default (docs/how-to/input_prep/examples/condition_spec_ppi.yaml): target backbone as one internal-coordinate group, full target sequence.
    small_molecule_binder: nothing is added; upstream sequence-conditions atomized tokens (ligands) itself and leaves the pose free.
    motif_scaffold: every preserved segment is a motif, fixed in coordinates and identity (docs: a motif carries both C_CRD and C_SEQ).
    """
    built: Dict[Any, Any] = {}
    if task == "protein_binder":
        chains = preserved_chains(tokens)
        built["C_CRDGRP"] = {1: [
            {"select": f"sel('{chain}/*/*/[N,CA,C,O]') and not has_nan_coord()"} for chain in chains]}
        built["C_SEQ"] = {True: [{"select": f"sel('{chain}')"} for chain in chains]}
    elif task == "motif_scaffold":
        selects = [_residue_select(t.chain, t.start, t.end) for t in tokens if t.kind == "preserved"]
        built["C_CRD"] = {True: [{"select": s} for s in selects]}
        built["C_SEQ"] = {True: [{"select": s} for s in selects]}

    spots = hotspot_selects(hotspots)
    if spots:
        built = merge_conditions(built, {"C_HOT": {True: [{"select": s} for s in spots]}})
    return merge_conditions(built, normalize_conditions(conditions))


def recycle_strategy(task: str, thinking: str) -> str:
    """The upstream `recycle_strategy` group, `<low|medium|high>/<shape>`, for a task and thinking budget."""
    return f"{thinking}/{TASK_RECYCLE_SHAPE[task]}"


class RFD4Proteina(BaseConfig):
    """
    RFD4Proteina: sequence-structure co-design with RFD4-Proteina, driven by a YAML condition spec.

    Every input structure becomes one design problem: its `contig` keeps the residues named by preserved segments (`A17-131`) and inserts variable-length segments the model designs (`50-120`), and a few conditions say what is fixed, what is a hotspot and what is free. The model returns designed structures and their sequences in one pass.

    Tasks (`task=`):
        protein_binder: IMPLEMENTED. Binder to a protein target. Target backbone is conditioned as one internal-coordinate group (`C_CRDGRP`), the full target sequence is given (`C_SEQ`) and `hotspots` mark the interface (`C_HOT`): upstream's recommended PPI recipe. Every preserved chain is taken to be a protein target.
        small_molecule_binder: IMPLEMENTED. Protein pocket around a ligand present in the input structure. The ligand is kept by a preserved segment in the contig (for example `150, /0, B386`) and upstream conditions its identity on its own; its pose and conformer are left to the model unless `conditions` fix them.
        motif_scaffold: IMPLEMENTED for indexed motifs. Every preserved segment is fixed in coordinates and identity (`C_CRD` and `C_SEQ`) and the model designs the scaffold around them.
        NOT IMPLEMENTED (reachable only through `conditions` or forwarded overrides, untested): unindexed motifs and enzyme active-site scaffolding with tip atoms (`C_IDX`), DNA binders, antibody and nanobody design, partial diffusion and inverse folding, reward-guided search, refolding metrics.

    Inputs:
        structures: input structures (target, motif or ligand complex). CIF files are used as they are; a PDB is converted to CIF with AtomWorks first, because upstream refuses PDB (it does not encode the bonds and charges the model consumes). Residue numbers in `contig` and `hotspots` are those of the input structure. A `Ligand(structures=..., codes=...)` output is accepted: its `structures` stream is the input.
        contig: upstream contig string, one for all inputs, e.g. "50-120, /0, A17-131". `/0` breaks the chain, `A17-131` keeps a range of the input, `50-120` is a designed segment of variable length, `E1:2-20` a named one. At least one preserved segment is required. `contigs` is accepted as a synonym.
        task: "protein_binder", "small_molecule_binder" or "motif_scaffold" (default "protein_binder").
        hotspots: interface atoms on the target as `C_HOT`: "A56:CG,OH;A115:CG,SD" (atoms), "A67,A89" (whole residues) or {"A56": "CG,OH"}.
        conditions: raw upstream conditions appended after the task recipe, {annotation: {value: [{"select": "sel('A/*/[10-20]')"}]}}, for example {"C_RAS": {"BURIED": [{"select": "sel('B')"}]}}. Annotation names and selection syntax: https://rfd4-proteina.org/latest/reference/yaml-json.html
        num_designs: designs per input structure (default 1). Upstream calls it `num_replicates`, accepted as a synonym.
        thinking: "low", "medium" or "high" trunk-evaluation budget over 60 denoising steps (default "medium"). For protein_binder it is 120, 180 or 300 evaluations (front-loaded); for small_molecule_binder and motif_scaffold 120, 150 or 180 (midpeak).
        batch_size: designs generated in parallel on the GPU (default 1); raise it only after the longest input fits at 1.
        seed: random seed (default None, a fresh one each run).
        metrics_tags: upstream inline-metric tags that need no refolding oracle, e.g. ["design", "hotspot_conditioning"] (default None: no scoring, the `~metrics` smoke-test path). Tags that run an AlphaFold 2/3 or OpenFold3 oracle (`*_refolding`, `antibody`) are refused: the install does not set those environments up. Score designs with the Boltz2 or AlphaFold tools instead.
        checkpoint: directory of a model bundle (`config.json` + `model.safetensors`) to load instead of the one RFD4Proteina.install() downloaded (default None).
        Further keyword arguments are forwarded as Hydra overrides, for example `generation.args.nsteps=80` or `compile_lengths="null"` (skips the minutes-long torch.compile warm-up of a first run).

    Outputs:
        Streams:
            structures: designed structures, one PDB per design, ids <input id>_<1..num_designs>.
            sequences: value-based, one row per polymer chain of each design.
        Tables:
            structures: id | source_id | design | cif | pdb | n_chains | n_residues
            sequences:  id | source_id | source_pdb | chain | sequence | length
            metrics:    id | source_id, plus one column per upstream inline metric when `metrics_tags` is set (columns are prefixed by the metric key, e.g. `clash_metrics.frac_atoms_clashing_1p5A`); the run record (`seed`, `ckpt_path`, `ckpt_name`) is kept and the `generation.*` configuration dump is dropped.
            missing:    id | removed_by | kind | cause

    Requirements:
        Linux x86_64, Python 3.12, an NVIDIA GPU of Ampere generation or newer and a CUDA 13-compatible driver. Not available on aarch64 (Daint).
        RFD4Proteina.install() creates the `RFD4Proteina` env (`rfd4-proteina[all-cuda13]`) and downloads the weights bundle (about 1.7 GB, `config.json` and `model.safetensors`) of `nvidia/RFD4-Proteina` into `folders.repositories.RFD4Proteina`.

    Weights and license:
        The weights are pinned to the Hub commit that the installed package's own `configs/inference.yaml` names (`hf_ckpt_revision`), resolved to a commit hash at install time and written to `REVISION` beside the bundle. Runs load that local directory, so they are offline and a later move of the Hub's `main` cannot change a run; the step log prints the revision.
        License, verbatim from the Hugging Face model card (front matter `license: bsd-3-clause`, section "License"): "The RFD4-Proteína model and weights are released under the BSD 3-Clause License, Copyright (c) 2026, Institute for Protein Design, University of Washington. Optional third-party assets distributed alongside the model — the AlphaFold2 and OpenFold3/AlphaFold3 refolding images and the OpenFold3 payload — carry their own upstream terms; see the note under *Other repository files* and their provenance files. The training data has its own source terms, described on the dataset card." The GitHub page states only the code license (BSD-3-Clause), so the model card is the source for the weights.

    Examples:
        # Binder to PD-L1: keep chain A 17-131, design a 50-120 residue binder
        binder = RFD4Proteina(
            structures=target, contig="50-120, /0, A17-131",
            hotspots="A56:CG,OH;A115:CG,SD;A123:CD2,OH", num_designs=8, batch_size=8)

        # Small-molecule binder: ligand is residue 386 of chain B of the input complex
        pocket = RFD4Proteina(
            structures=complex_with_ligand, contig="150, /0, B386",
            task="small_molecule_binder", thinking="medium", num_designs=8)

        # Indexed motif scaffolding: keep A10-25 and A60-70, design 80-120 residues around them
        scaffold = RFD4Proteina(
            structures=motif, contig="A10-25, 20-40, A60-70, 20-40", task="motif_scaffold", num_designs=8)

    References:
        Code: https://github.com/RosettaCommons/RFD4-Proteina
        Docs: https://rfd4-proteina.org/latest/
        Weights: https://huggingface.co/nvidia/RFD4-Proteina
    """

    TOOL_NAME = "RFD4Proteina"
    TOOL_VERSION = "1.0"
    # Hydra application: untyped keyword arguments become extra `key=value` overrides.
    FORWARD_UNKNOWN_KWARGS = "hydra"
    # `num_replicates` is upstream's own name; forwarded it would collide with the managed override.
    PARAMETER_ALIASES = {"num_replicates": "num_designs", "contigs": "contig"}
    ENV_NAME = "RFD4Proteina"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        env = cls._install_env(env_manager)
        biopipelines = folders.get("biopipelines", "")
        weights_dir = folders.get("RFD4Proteina", "")
        if not weights_dir:
            raise ValueError(
                "RFD4Proteina: folders.repositories.RFD4Proteina is not set in the config, so there is "
                "nowhere to put the weights bundle. Add `RFD4Proteina: \"<data>/rfd4_proteina\"`.")
        env_check = cls._env_exists_check(env, env_manager)
        env_run = cls._env_run(env, env_manager)
        pipe_scripts = folders.get("pipe_scripts", f"{biopipelines}/pipe_scripts")
        weights_py = f"{pipe_scripts}/pipe_rfd4proteina_weights.py"
        # Entry point present in the env, without importing the model (a login node has no GPU).
        installed = (
            f'{env_run}python -c "import importlib.metadata as m; m.version(\'rfd4-proteina\'); '
            f'assert any(e.name == \'rfd4-generate\' for e in m.entry_points(group=\'console_scripts\'))" '
            '>/dev/null 2>&1'
        )
        weights_ok = f'[ -f "{weights_dir}/model.safetensors" ] && [ -f "{weights_dir}/config.json" ]'
        skip = "" if force_reinstall else f"""# Check if already fully installed (env + weights bundle)
if {installed} && {weights_ok}; then
    echo "RFD4Proteina already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block(env, env_manager) if force_reinstall else ""
        env_block = cls._env_install_block(env, env_manager, biopipelines)
        return f"""echo "=== Installing RFD4Proteina ({env}) ==="
{skip}{remove_block}
# Create {env} env (skip if it already exists)
if ! {env_check}; then
    {env_block}
else
    echo "{env} environment already exists, skipping creation."
fi

# Download the weights bundle pinned by the installed package (skipped when present)
mkdir -p "{weights_dir}"
if {weights_ok}; then
    echo "Weights already present in {weights_dir}, skipping download."
else
    {env_run}python "{weights_py}" --dest "{weights_dir}"
fi

# Verify installation
if {installed} && {weights_ok}; then
    touch "$INSTALL_SUCCESS"
    echo "=== RFD4Proteina installation complete ==="
else
    echo "ERROR: RFD4Proteina verification failed (rfd4-proteina missing from the env or weights bundle incomplete)"
    exit 1
fi
"""

    structures_json = Path(lambda self: self.configuration_path("structures.json"))
    config_yaml = Path(lambda self: self.configuration_path("config.yaml"))
    spec_yaml = Path(lambda self: self.configuration_path("condition_spec.yaml"))
    inputs_dir = Path(lambda self: self.execution_path("inputs"))
    manifest_json = Path(lambda self: self.execution_path("inputs_manifest.json"))
    run_dir = Path(lambda self: self.execution_path("rfd4"))
    hydra_dir = Path(lambda self: self.execution_path("hydra"))
    structures_map = Path(lambda self: self.stream_map_path("structures"))
    sequences_csv = Path(lambda self: self.stream_path("sequences", "sequences.csv"))
    structures_csv = Path(lambda self: self.table_path("structures"))
    metrics_csv = Path(lambda self: self.table_path("metrics"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    local_missing_csv = Path(lambda self: self.execution_path("local_missing.csv"))
    inputs_py = Path(lambda self: self.pipe_script_path("pipe_rfd4proteina_inputs.py"))
    postprocess_py = Path(lambda self: self.pipe_script_path("pipe_rfd4proteina_postprocess.py"))
    checkpoint_dir = Path(lambda self: self.checkpoint or self.folders.get("RFD4Proteina", ""))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 contig: str = None,
                 task: str = "protein_binder",
                 hotspots: Union[str, Dict[str, str]] = None,
                 conditions: Dict[Any, Any] = None,
                 num_designs: int = 1,
                 thinking: str = "medium",
                 batch_size: int = 1,
                 seed: int = None,
                 metrics_tags: List[str] = None,
                 checkpoint: str = None,
                 **kwargs):
        self.structures = structures  # raw handle kept for missing-propagation
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        self.contig = contig
        self.task = task
        self.hotspots = hotspots
        self.conditions = conditions
        self.num_designs = num_designs
        self.thinking = thinking
        self.batch_size = batch_size
        self.seed = seed
        self.metrics_tags = metrics_tags
        self.checkpoint = checkpoint
        super().__init__(**kwargs)

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        if self.task not in TASKS:
            raise ValueError(f"task must be one of {TASKS}, got {self.task!r}")
        if self.thinking not in THINKING_LEVELS:
            raise ValueError(f"thinking must be one of {THINKING_LEVELS}, got {self.thinking!r}")

        self._tokens = parse_contig(self.contig)
        self._conditions = build_conditions(self.task, self._tokens, self.hotspots, self.conditions)

        for name, value, minimum in (("num_designs", self.num_designs, 1), ("batch_size", self.batch_size, 1)):
            if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
                raise ValueError(f"{name} must be an integer >= {minimum}, got {value!r}")
        if self.seed is not None and (isinstance(self.seed, bool) or not isinstance(self.seed, int)
                                      or self.seed < 0):
            raise ValueError(f"seed must be a non-negative integer, got {self.seed!r}")

        if self.metrics_tags is not None:
            if isinstance(self.metrics_tags, str) or not isinstance(self.metrics_tags, (list, tuple)):
                raise ValueError("metrics_tags must be a list of tag strings, e.g. ['design']")
            for tag in self.metrics_tags:
                if not isinstance(tag, str) or not _TAG_RE.match(tag):
                    raise ValueError(f"metrics tag {tag!r} must be a string of letters, digits, '_', '.', ':' or '-'")
                if tag.endswith("_refolding") or tag == "antibody":
                    raise ValueError(
                        f"metrics tag {tag!r} runs an AlphaFold 2/3 or OpenFold3 refolding oracle, which "
                        "RFD4Proteina.install() does not set up. Refold the designs with Boltz2 or AlphaFold instead.")
        _validate_freeform_string("checkpoint", self.checkpoint)

        managed = {k.lstrip("+") for k in self.extra_args} & set(MANAGED_OVERRIDES)
        if managed:
            raise ValueError(
                f"{sorted(managed)} are set by RFD4Proteina itself and cannot be forwarded as overrides; "
                "use the matching named parameter (num_designs, batch_size, seed, thinking, metrics_tags, checkpoint).")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def _spec_template(self) -> Dict[str, Any]:
        """The per-design entry shared by every input: its contig and conditions, in upstream's key order."""
        template: Dict[str, Any] = {"contig": self.contig.strip()}
        if self._conditions:
            template["conditions"] = self._conditions
        return template

    def _hydra_overrides(self) -> List[str]:
        """Bash-ready `key=value` tokens for everything the wrapper manages."""
        tokens = [
            f'inputs="{self.spec_yaml}"',
            f'out_dir="{self.run_dir}"',
            f'hydra.run.dir="{self.hydra_dir}"',
            f'hf_ckpt_path="{self.checkpoint_dir}"',
            f"num_replicates={self.num_designs}",
            f"batch_size={self.batch_size}",
            f"recycle_strategy={recycle_strategy(self.task, self.thinking)}",
        ]
        if self.seed is not None:
            tokens.append(f"seed={self.seed}")
        if self.metrics_tags:
            tokens.append(f"metrics_tags='[{','.join(self.metrics_tags)}]'")
        else:
            tokens.append("'~metrics'")
        return tokens

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"TASK: {self.task}")
        lines.append(f"CONTIG: {self.contig.strip()}")
        if self.hotspots:
            lines.append(f"HOTSPOTS: {self.hotspots}")
        lines.append(f"THINKING: {self.thinking} ({recycle_strategy(self.task, self.thinking)})")
        lines.append(f"NUM DESIGNS PER INPUT: {self.num_designs}")
        lines.append(f"BATCH SIZE: {self.batch_size}")
        if self.seed is not None:
            lines.append(f"SEED: {self.seed}")
        lines.append(f"INLINE METRICS: {', '.join(self.metrics_tags) if self.metrics_tags else 'none'}")
        if self.checkpoint:
            lines.append(f"CHECKPOINT: {self.checkpoint}")
        return lines

    def _write_config_yaml(self):
        """Dump the shared design entry and the replicate count for the input builder and post-processing."""
        import yaml
        document = {"num_designs": self.num_designs, "spec": self._spec_template()}
        with open(self.config_yaml, "w") as f:
            yaml.safe_dump(document, f, sort_keys=False, width=1000)

    def _generate_script_inputs(self) -> str:
        return f"""echo "Preparing RFD4-Proteina inputs"
python "{self.inputs_py}" \\
    --config-yaml "{self.config_yaml}" \\
    --structures-json "{self.structures_json}" \\
    --inputs-dir "{self.inputs_dir}" \\
    --spec-yaml "{self.spec_yaml}" \\
    --manifest-json "{self.manifest_json}"
RFD4_INPUTS_RC=$?

"""

    def _generate_script_generate(self) -> str:
        if not self.checkpoint_dir:
            raise ValueError(
                "RFD4Proteina: no weights bundle configured. Add `RFD4Proteina: \"<data>/rfd4_proteina\"` under "
                "folders.repositories and run RFD4Proteina.install(), or pass checkpoint=<bundle directory>.")
        overrides = self._hydra_overrides() + self.extra_args_bash_tokens()
        command = "".join(f" \\\n    {token}" for token in overrides)
        return f"""RFD4_RC=1
if [ "$RFD4_INPUTS_RC" -ne 0 ]; then
    echo "ERROR: RFD4Proteina could not prepare its inputs; skipping generation"
elif [ ! -f "{self.checkpoint_dir}/model.safetensors" ] || [ ! -f "{self.checkpoint_dir}/config.json" ]; then
    echo "ERROR: RFD4-Proteina weights not found in {self.checkpoint_dir}"
    echo "Run RFD4Proteina.install() or pass checkpoint=<bundle directory>"
else
    echo "Weights bundle: {self.checkpoint_dir}"
    if [ -f "{self.checkpoint_dir}/REVISION" ]; then
        echo "Weights revision: $(cat "{self.checkpoint_dir}/REVISION")"
    fi
    echo "Starting RFD4-Proteina"
    {self.container_prefix()}rfd4-generate{command}
    RFD4_RC=$?
fi
echo "rfd4-generate exit status: $RFD4_RC"

"""

    def _generate_script_postprocess(self) -> str:
        return f"""echo "Post-processing RFD4-Proteina outputs"
{self.activate_environment(name="biopipelines")}
python "{self.postprocess_py}" \\
    --run-dir "{self.run_dir}" \\
    --manifest-json "{self.manifest_json}" \\
    --num-designs {self.num_designs} \\
    --structures-dir "{self.stream_folder('structures')}" \\
    --structures-map-csv "{self.structures_map}" \\
    --sequences-csv "{self.sequences_csv}" \\
    --structures-csv "{self.structures_csv}" \\
    --metrics-csv "{self.metrics_csv}" \\
    --local-missing-csv "{self.local_missing_csv}"
"""

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        self._write_config_yaml()

        script = "#!/bin/bash\n"
        script += "# RFD4Proteina execution script\n"
        script += "# Generated by BioPipelines\n\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += self.extra_args_echo()
        script += self._generate_script_inputs()
        script += self._generate_script_generate()
        script += self._generate_script_postprocess()
        script += self.generate_missing_propagation(
            self.structures, local_missing=self.local_missing_csv, missing_csv=self.missing_csv)
        script += self.generate_completion_check_footer()
        return script

    def _sequence_row_ids(self, structure_ids: List[str]) -> List[str]:
        """Declared ids of the `sequences` rows, which are per chain: `<structure>` or `<structure>_<chain>`.

        More than one run in the contig means a design can have several polymer chains, and the writer then appends `_<chain>`; the lazy bracket covers both shapes, as in RFdiffusion3.
        """
        if len(contig_runs(self._tokens)) <= 1:
            return list(structure_ids)
        return [f"{sid}[_<?>]" for sid in structure_ids]

    def get_output_files(self) -> Dict[str, Any]:
        structure_ids = generate_multiplied_ids_pattern(
            self.structures_stream.ids, f"<1..{self.num_designs}>", input_stream_name="structures")

        structures = DataStream(
            name="structures",
            ids=structure_ids,
            files=[self.stream_path("structures", "<id>.pdb")],
            map_table=self.structures_map,
            format="pdb",
        )
        sequences = DataStream(
            name="sequences",
            ids=self._sequence_row_ids(structure_ids),
            files=[],
            map_table=self.sequences_csv,
            format="csv",
        )
        tables = {
            "structures": TableInfo(
                name="structures",
                path=self.structures_csv,
                columns=["id", "source_id", "design", "cif", "pdb", "n_chains", "n_residues"],
                description="RFD4-Proteina designs: source structure, replicate, annotated CIF and PDB paths, size",
            ),
            "sequences": TableInfo(
                name="sequences",
                path=self.sequences_csv,
                columns=["id", "source_id", "source_pdb", "chain", "sequence", "length"],
                description="Designed sequence of each polymer chain",
            ),
            "metrics": TableInfo(
                name="metrics",
                path=self.metrics_csv,
                columns=["id", "source_id"],
                description="Upstream inline metrics per design (one column per metric when metrics_tags is set)",
            ),
            "missing": self.missing_table_info(self.missing_csv),
        }
        return {
            "structures": structures,
            "sequences": sequences,
            "tables": tables,
            "output_folder": self.output_folder,
        }

    def to_dict(self) -> Dict[str, Any]:
        base_dict = super().to_dict()
        base_dict.update({
            "rfd4proteina_params": {
                "task": self.task,
                "contig": self.contig,
                "hotspots": self.hotspots,
                "conditions": self._conditions,
                "num_designs": self.num_designs,
                "thinking": self.thinking,
                "recycle_strategy": recycle_strategy(self.task, self.thinking),
                "batch_size": self.batch_size,
                "seed": self.seed,
                "metrics_tags": self.metrics_tags,
                "checkpoint": self.checkpoint,
                "structures_ids": self.structures_stream.ids,
            }
        })
        return base_dict
