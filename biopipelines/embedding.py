# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Embedding tool: one fixed-length vector per input from a registered embedding model, and cosine similarities between chosen pairs.

The tool is model-agnostic. Every model it can run is one entry of `EMBEDDING_MODELS`, a frozen `EmbeddingModel` that declares everything the wrapper and the pipe script need to know about it, and nothing else in the tool names a model. Today the registry holds the three structure-graph models of the Progres family; a model of another input kind (a sequence embedder such as ESM2's) is one more entry plus an adapter class in `pipe_scripts/pipe_embedding.py`.

REGISTRY CONTRACT. Each `EmbeddingModel` declares:
    name             the value of the `model` parameter.
    input_kind       the stream the model reads, "structures" for both entries; the wrapper takes `inputs.streams.<input_kind>` from a standardized output.
    formats          file formats accepted on that stream; a stream advertising another format is refused at configuration time.
    adapter          key into `ADAPTERS` of `pipe_scripts/pipe_embedding.py`, the class that builds the model from the weights file and embeds one input.
    weights_env      the environment variable the generated script sets to the weights file; the adapter reads the path from it. A user-supplied file path is the only weights source: nothing is downloaded at run time.
    upstream_default where upstream ships a default weights file, as a note shown when `weights` is wrong; "" when upstream names none.
    dimension        length of one embedding vector; an adapter returning another length fails that id.
    node_module      name of the submodule whose output is the per-node embedding, "" when the model has no per-node output; `node_embeddings=True` needs a non-empty value.
    extra_args       keyword arguments the adapter passes to the upstream reader.

THE PROGRES FAMILY. Progres embeds a protein as a graph of its C-alpha atoms (an edge joins two atoms within 10 A) and runs six EGNN layers over it; a per-node decoder, a sum over nodes and a graph decoder give one 128-number vector, L2-normalized, so the dot product of two embeddings is their cosine. CIRPIN is the same network trained with synthetic circular permutations so that the embedding does not depend on where the termini sit. The two entries share one architecture and one checkpoint layout (a dict whose `model` key holds the state dict), and differ only in the weights file you pass. A checkpoint of another architecture fails the step at load with torch's own key or shape message, never silently.

References:
    Progres: Greener & Jamali (2025) Fast protein structure searching using structure graph embeddings, Bioinformatics Advances 5, vbaf042. https://github.com/greener-group/progres
    CIRPIN:  Kolodziej, Abulnaga & Ovchinnikov (2025) CIRPIN: Learning Circular Permutation-Invariant Representations to Uncover Putative Protein Homologs, bioRxiv 10.1101/2025.11.18.689110. https://github.com/aidenkoloj/CIRPIN
"""

import os
import shlex
from dataclasses import dataclass, field
from typing import Dict, List, Any, Mapping, Optional, Tuple, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
    from .data_containers import resolve_table_reference
    from .biopipelines_io import TableReference
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream
    from data_containers import resolve_table_reference
    from biopipelines_io import TableReference


@dataclass(frozen=True)
class EmbeddingModel:
    """One registered embedding model; see the registry contract in the module docstring."""
    name: str
    input_kind: str
    formats: Tuple[str, ...]
    adapter: str
    weights_env: str
    upstream_default: str
    dimension: int
    node_module: str = ""
    extra_args: Mapping[str, Any] = field(default_factory=dict)


_PROGRES_FAMILY = dict(
    input_kind="structures",
    formats=("pdb", "cif", "mmcif"),
    adapter="progres_graph",
    weights_env="TRAINED_MODEL_FP",
    dimension=128,
    node_module="node_dec",
    extra_args={"fileformat": "guess"},
)

EMBEDDING_MODELS: Dict[str, EmbeddingModel] = {
    entry.name: entry for entry in (
        EmbeddingModel(
            name="progres",
            upstream_default=(
                "the progres package downloads trained_model.pt from Zenodo (record 18245422) into "
                "$PROGRES_DATA_DIR/trained_models/v_0_2_0/ the first time its own command line runs, and the CIRPIN "
                "repository ships a copy at trained_models/Progres_model/Progres_trained_model.pt; "
                "Embedding.install() downloads neither"),
            **_PROGRES_FAMILY),
        EmbeddingModel(
            name="cirpin",
            upstream_default=(
                "the CIRPIN repository (github.com/aidenkoloj/CIRPIN) ships "
                "trained_models/CIRPIN_model/CIRPIN_model_5k_cp_epoch301.pt; Embedding.install() does not download it"),
            **_PROGRES_FAMILY),
    )
}

DEVICES = ("auto", "cpu", "cuda")

SIMILARITY_COLUMNS = ["id", "first", "second", "cosine_similarity"]


def registry_entry(model: str) -> EmbeddingModel:
    """The registry entry for `model`, or a ValueError naming the registered models."""
    if not isinstance(model, str) or model not in EMBEDDING_MODELS:
        raise ValueError(f"model must be one of {tuple(EMBEDDING_MODELS)}, got {model!r}")
    return EMBEDDING_MODELS[model]


def _pair_reference(value, parameter: str) -> Optional[TableReference]:
    """A column reference as a TableReference; None passes through and a literal string is refused."""
    resolved = resolve_table_reference(value, parameter)
    if resolved is not None and not isinstance(resolved, TableReference):
        raise ValueError(
            f"{parameter} must be a column of the pair table, e.g. pairs.tables.data.query, "
            f"or a (table, column) tuple; got {type(resolved).__name__}")
    return resolved


class Embedding(BaseConfig):
    """
    Embedding: one fixed-length vector per input from a registered embedding model, and cosine similarities between chosen pairs.

    Runs the model named by `model` (see the registry in this module) with the weights file you pass and writes one `.npy` vector per input. Output ids equal input ids. Today the models are `progres` and `cirpin`, both of which read structures.

    Inputs:
        inputs:          the structures to embed (a stream or a standardized output with a `structures` stream; PDB or mmCIF). Only the first chain of the first model is read, from its C-alpha atoms of ATOM records; HETATM records, other chains and other models are ignored without a warning, so pass one chain per id.
        model:           registered model, one of progres, cirpin. Required: the two share one network, so the name is what records which weights you meant.
        weights:         path to the checkpoint file of that model (a torch file whose `model` key is the state dict). Must exist when the pipeline is configured; the step exports it as TRAINED_MODEL_FP. No default is downloaded.
        first, second:   the two members of each pair to compare, as columns of one table, e.g. first=pairs.tables.data.query, second=pairs.tables.data.target. Both or neither; both must come from the same table (a pair per row). Ids are read as text, so `007` and `NA` stay as written, and compared with the ids of `inputs`; a pair listed twice is one pair, and two different pairs that would join to one id are refused. Omit both for no similarity table.
        node_embeddings: also write the `node_embeddings` stream, one `.npy` per input of shape (residues, dimension) (default: False). It is the output of the per-node decoder before the sum over residues, not normalized. The pip package of Progres has no such option; CIRPIN's own progres.py returns this tensor from `embed_nodes`.
        device:          "cpu", "cuda" or "auto" (cuda when available, else cpu) (default: "cpu", upstream's own default). "cuda" without a GPU fails the step rather than falling back.

    Embeddings: float32, shape (128,) for both registered models, L2-normalized by the network. Two structures with the same fold sit close together.

    Similarity: cosine similarity of the two embeddings, in [-1, 1]. It is NOT Progres's own score, which is (1 + cosine) / 2 and runs 0 to 1: convert with `(1 + cosine_similarity) / 2` downstream when comparing against a Progres threshold such as 0.8.

    Outputs:
        Streams:
            embeddings:      <id>.npy per input, the vector above (format npy).
            node_embeddings: <id>.npy per input, shape (residues, dimension) (format npy). Empty unless node_embeddings=True.
        Tables:
            similarities: id | first | second | cosine_similarity   (one row per distinct pair of the table, in first-seen order; id is `<first>+<second>`. NaN when either member has no embedding: empty, not in the input stream, dropped upstream, or its embedding failed. Absent without `first` and `second`.)
            missing:      id | removed_by | kind | cause   (an input that failed to embed, under its own id, kind="failure"; a pair that could not be compared, under its pair id: kind="unmatched" when a member is empty, not in the input stream or dropped upstream, kind="failure" when a member failed to embed)
    """

    TOOL_NAME = "Embedding"
    TOOL_VERSION = "1.3"
    ENV_NAME = "Embedding"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        env = cls._install_env(env_manager)
        biopipelines = folders.get("biopipelines", "")
        env_check = cls._env_exists_check(env, env_manager)
        env_run = cls._env_run(env, env_manager)
        imports_ok = (f'{env_run}python -c "import torch, torch_scatter, torch_geometric, einops, progres.progres, numpy, pandas, yaml" '
                      '>/dev/null 2>&1')
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_check} && {imports_ok}; then
    echo "Embedding already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block(env, env_manager) if force_reinstall else ""
        env_block = cls._env_install_block(env, env_manager, biopipelines)
        return f"""echo "=== Installing Embedding ({env}) ==="
{skip}{remove_block}
# Create {env} env (skip if it already exists)
if ! {env_check}; then
    {env_block}
else
    echo "{env} environment already exists, skipping creation."
fi

# The package is installed and imports; no weights are downloaded, the user supplies the checkpoint
if {imports_ok}; then
    touch "$INSTALL_SUCCESS"
    echo "=== Embedding installation complete ==="
else
    echo "ERROR: Embedding verification failed (torch, torch_scatter, torch_geometric, einops, progres, numpy, pandas or yaml missing from the env)"
    exit 1
fi
"""

    inputs_json = Path(lambda self: self.configuration_path("inputs.json"))
    config_yaml = Path(lambda self: self.configuration_path("config.yaml"))
    embeddings_map = Path(lambda self: self.stream_map_path("embeddings"))
    node_embeddings_map = Path(lambda self: self.stream_map_path("node_embeddings"))
    similarities_csv = Path(lambda self: self.table_path("similarities"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    local_missing_csv = Path(lambda self: self.execution_path("local_missing.csv"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_embedding.py"))

    def __init__(self,
                 inputs: Union[DataStream, StandardizedOutput],
                 model: str,
                 weights: str,
                 first: Optional[Union[TableReference, Tuple[Any, str]]] = None,
                 second: Optional[Union[TableReference, Tuple[Any, str]]] = None,
                 node_embeddings: bool = False,
                 device: str = "cpu",
                 **kwargs):
        entry = registry_entry(model)
        self.model = model
        self.inputs = inputs
        if isinstance(inputs, StandardizedOutput):
            stream = getattr(inputs.streams, entry.input_kind)
            if stream is None:
                raise ValueError(f"inputs has no '{entry.input_kind}' stream, which model {model!r} reads")
            self.inputs_stream: DataStream = stream
        elif isinstance(inputs, DataStream):
            self.inputs_stream = inputs
        else:
            raise ValueError(f"inputs must be DataStream or StandardizedOutput, got {type(inputs)}")
        self.weights = weights
        self.first = _pair_reference(first, "first")
        self.second = _pair_reference(second, "second")
        self.node_embeddings = node_embeddings
        self.device = device
        super().__init__(**kwargs)

    @property
    def entry(self) -> EmbeddingModel:
        return EMBEDDING_MODELS[self.model]

    def validate_params(self):
        entry = self.entry
        if not self.inputs_stream or len(self.inputs_stream) == 0:
            raise ValueError("inputs parameter is required and must not be empty")
        if not self.inputs_stream.has_only_formats(*entry.formats):
            raise ValueError(
                f"model {self.model!r} reads {entry.input_kind} in {entry.formats}, "
                f"got format {self.inputs_stream.format!r}")
        self._check_weights(entry)
        if (self.first is None) != (self.second is None):
            raise ValueError("first and second are given together or not at all")
        if self.first is not None and self.first.path != self.second.path:
            raise ValueError(
                f"first and second must be columns of the same table, got {self.first.path} and {self.second.path}")
        if not isinstance(self.node_embeddings, bool):
            raise ValueError(f"node_embeddings must be True or False, got {self.node_embeddings!r}")
        if self.node_embeddings and not entry.node_module:
            raise ValueError(f"model {self.model!r} has no per-node output, so node_embeddings=True is not available")
        if self.device not in DEVICES:
            raise ValueError(f"device must be one of {DEVICES}, got {self.device!r}")

    def _check_weights(self, entry: EmbeddingModel):
        """The weights are a local file that exists now; a URL, a directory or a missing file is refused with where upstream keeps a default."""
        hint = f" Upstream: {entry.upstream_default}." if entry.upstream_default else " Upstream names no default file."
        if not isinstance(self.weights, str) or not self.weights.strip():
            raise ValueError(f"weights must be the path to the {self.model} checkpoint file.{hint}")
        if "://" in self.weights:
            raise ValueError(f"weights must be a local file, not a URL ({self.weights}); download it first.{hint}")
        if not os.path.isfile(self._weights_path()):
            raise ValueError(f"weights file not found: {self._weights_path()}.{hint}")

    def _weights_path(self) -> str:
        return os.path.abspath(os.path.expanduser(self.weights))

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"MODEL: {self.model}")
        lines.append(f"WEIGHTS: {self._weights_path()}")
        lines.append(f"DEVICE: {self.device}")
        if self.first is not None:
            lines.append(f"PAIRS: {self.first.path} ({self.first.column}, {self.second.column})")
        if self.node_embeddings:
            lines.append("NODE EMBEDDINGS: yes")
        return lines

    def generate_script(self, script_path: str) -> str:
        entry = self.entry
        self.inputs_stream.save_json(self.inputs_json)
        self._write_config_yaml()
        script = "#!/bin/bash\n"
        script += "# Embedding script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += self.warn_container_unsupported()
        script += f"export {entry.weights_env}={shlex.quote(self._weights_path())}\n"
        script += f"""echo "Embedding {len(self.inputs_stream)} input(s) with {self.model}"
python "{self.helper_py}" \\
    --config-yaml "{self.config_yaml}" \\
    --inputs-json "{self.inputs_json}" \\
    --embeddings-dir "{self.stream_folder('embeddings')}" \\
    --embeddings-map-csv "{self.embeddings_map}" \\
    --node-embeddings-dir "{self.stream_folder('node_embeddings')}" \\
    --node-embeddings-map-csv "{self.node_embeddings_map}" \\
    --similarities-csv "{self.similarities_csv}" \\
    --local-missing-csv "{self.local_missing_csv}"
"""
        script += self.generate_missing_propagation(
            self.inputs, local_missing=self.local_missing_csv, missing_csv=self.missing_csv
        )
        script += self.generate_completion_check_footer()
        return script

    def _write_config_yaml(self):
        import yaml
        with open(self.config_yaml, "w") as f:
            yaml.safe_dump({
                "model": self.model,
                "node_embeddings": self.node_embeddings,
                "device": self.device,
                "first": str(self.first) if self.first is not None else None,
                "second": str(self.second) if self.second is not None else None,
            }, f, sort_keys=False)

    def get_output_files(self) -> Dict[str, Any]:
        embeddings = DataStream(
            name="embeddings",
            ids=self.inputs_stream.ids,
            files=[self.stream_path("embeddings", "<id>.npy")],
            map_table=self.embeddings_map,
            format="npy",
        )
        if self.node_embeddings:
            node_embeddings = DataStream(
                name="node_embeddings",
                ids=self.inputs_stream.ids,
                files=[self.stream_path("node_embeddings", "<id>.npy")],
                map_table=self.node_embeddings_map,
                format="npy",
            )
        else:
            node_embeddings = DataStream.empty("node_embeddings", "npy")
        tables = {"missing": self.missing_table_info(self.missing_csv)}
        if self.first is not None:
            tables["similarities"] = TableInfo(
                name="similarities",
                path=self.similarities_csv,
                columns=SIMILARITY_COLUMNS,
                description="Cosine similarity of the embeddings of each requested pair (NaN when a member has no embedding)",
            )
        return {
            "embeddings": embeddings,
            "node_embeddings": node_embeddings,
            "tables": tables,
            "output_folder": self.output_folder,
        }
