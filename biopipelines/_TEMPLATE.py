# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

# To add a tool: copy this file to biopipelines/<yourtool>.py, then walk the companion checklist in `_TEMPLATE_CHECKLIST.md`. Delete the annotation blocks as you go.

import os
from typing import Dict, List, Any, Union

# Dual import so the module works both as a package member (`from biopipelines import X`, in-env) and when a pipe_script adds biopipelines/ to sys.path (out-of-env execution). Both run modes are supported on purpose; do not collapse to a single import. Copy it as-is and only edit the symbol list.

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream


class TemplateTool(BaseConfig):
    # Document EVERY input param and EVERY output stream/table with its columns.
    """
    TemplateTool: counts residues per input structure.

    Inputs:
        structures: PDB structures.

    Outputs:
        Streams:
            annotated: one <id>.txt per input.
        Tables:
            counts:  id | n_residues
            missing: id | removed_by | kind | cause   (see MISSING block below)
    """

    """STREAM MODEL
    - Configuration: the wrapper declares outputs so downstream tools can wire up before execution. `files=` here is a PREDICTION of where each output will be ("<id>.pdb"), and ids= the predicted set. No map exists yet.
    - Execution: the map_table the pipe writes is the SOURCE OF TRUTH. The framework reads the real id set from the map (DataStream.ids_expanded reads `id` from map_table in runtime mode), so if the tool filtered/failed ids, the map has fewer rows than files= predicted - and the map wins. `files=` is never re-read as truth at runtime.
    Further: the ids represent a view over the entities, i.e., when passing a set of IDs smaller than what is contained in the map, the ID set wins. For this reason we can iterate a stream over its IDs, yielding DataStream objects with identical map table but reduced set of IDs.


    THREE STREAM SHAPES:
    1. per-id files: files=["<id>.ext"], one file each. map_table = id|file. Consume with iterate_files(ds) -> (id, path). (this template's `annotated`)
    2. shared file: files="one.fasta" (bare filename) - one artifact whose records correspond to all ids (e.g. a multi-record FASTA). map_table still lists the ids. iterate_files yields (id, same_path) for every id, so parse the artifact once and index by id rather than re-reading per id. The format must have a slicer registered in biopipelines/stream_slicers.py. 
    3. value-based: files=[] - the entity lives inline in map_table columns (a `sequences` stream's `sequence`, `compounds`' `smiles`). Consume with iterate_values(ds, columns=[...]) -> (id, {col: val}). Prefer this for short strings/numbers.
    In all three, the map is the runtime source of truth for which entities exist; the IDs for which entities to iterate over.
    """

    # TOOL_NAME is the public class handle and used as key in the config yaml files (containers, repos, environments, tool overrides)
    # Editing this wrapper (or its pipe scripts) requires bumping TOOL_VERSION and registering the tool in versions/tool_changelog.yaml and adding a CHANGELOG.md bullet.
    TOOL_NAME = "TemplateTool"
    TOOL_VERSION = "0.1"

    """INSTALLATION.
    Triggered by <Tool>.install() inside an install pipeline (see example_pipelines/install_tools.py), submitted like any pipeline. 
    Pick ONE variant:
    - dedicated env > keep this method. Create the env, git clone the repo, download weights. Ship a generic environments/<TOOL_NAME>.yaml and, if installation on Colab diverges, also <TOOL_NAME>.colab.yaml.
    - base-env tool > DELETE this method; runs in the biopipelines base env, no yaml.
    - shares another tool's env > DELEGATE: `return PyMOL._install_script(folders, env_manager=env_manager, force_reinstall=force_reinstall, **kwargs)`. A PyMOL-only tool (uses PyMOL/ProteinEnv, no deps of its own) should reuse that env instead of shipping its own yaml, so a user calling YourTool.install() doesn't need to know the dependency. Exemplars: contacts.py, sasa.py, conformational_change.py.
    - container/.sif > keep this method too: it still installs. Add the image SOURCE to environments/_containers.yaml (keyed by TOOL_NAME; a docker:// URI or a direct .sif URL) and call cls._container_pull_block(folders, force_reinstall) in the install; it pulls to the DESTINATION path the user configured under `containers:` in config.<variant>.yaml. Execution then flows through container_prefix() automatically (see generate_script). ⚠ Emit the pull BEFORE any early-exit "already installed" skip, or a user who enables a container after an env-mode install never pulls the image. (No tool wires this yet; rfdiffusion2.py is not a model - its image ships with the repo.)
    Colab fork: if either installation or execution diverges on Colab, guard on scheduler == "colab" rather than duplicating logic. 
    ⚠ Gate on _env_exists_check so re-installs are cheap; honor force_reinstall via _env_remove_block; verify the binary/imports and `touch "$INSTALL_SUCCESS"` or the installer reports failure.
    """
    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        biopipelines = folders.get("biopipelines", "")
        env_check = cls._env_exists_check("templatetool", env_manager)
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_check}; then
    echo "TemplateTool already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block("templatetool", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("templatetool", env_manager, biopipelines)
        return f"""echo "=== Installing TemplateTool ==="
{skip}{remove_block}
{env_block}

if {env_manager} run -n templatetool python -c "import sys" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== TemplateTool installation complete ==="
else
    echo "ERROR: TemplateTool verification failed"
    exit 1
fi
"""
    """PATH DESCRIPTORS
    `Path(lambda self: ...)` is a lazy descriptor: accessing self.counts_csv computes the path once the tool has pipeline context. We define them this way because paths depend on execution_order / output_folder, unknown at __init__. Declaring them here keeps every path in one place, reused by generate_script, get_output_files, or any other method.
    Helpers can be used that resolve under the tool's numbered output folder: 
    - configuration_path(f) > _configuration/f (JSON inputs, arg yaml); 
    - stream_folder(name) > <name>/ (a stream's file dir); 
    - stream_path(name,*p) > a file within it; 
    - stream_map_path(name) > the stream's id>file map CSV; 
    - table_path(name) > tables/<name>.csv; 
    - pipe_script_path(f) > the paired pipe_scripts/<f>; 
    - execution_path(*p) > scratch/raw dumps not part of the contract.
    ⚠ Every output stream returned in get_output_files must have a stream_map_path here, or downstream tools can't resolve its files.
    """
    structures_json = Path(lambda self: self.configuration_path("structures.json"))
    config_yaml = Path(lambda self: self.configuration_path("config.yaml"))
    annotated_map = Path(lambda self: self.stream_map_path("annotated"))
    counts_csv = Path(lambda self: self.table_path("counts"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    local_missing_csv = Path(lambda self: self.execution_path("local_missing.csv"))
    helper_py = Path(lambda self: self.pipe_script_path("_pipe_template.py"))

    """INPUT CONTRACT
    Normalization: an input may arrive as a StandardizedOutput (a whole upstream tool) or a bare DataStream; extract the stream you need from `.streams.<name>` in the former case. Keep the original handle too (self.structures) - the missing-propagation helpers need it. The isinstance chain below IS the convention for a stream-only input; it is byte-identical in dssp.py, apbs.py and ~40 other tools. Copy it.
    For basic input types - anything represented by streams - plan to receive a stream or standardized output rather than a bare value. For example, Tool(Sequence("MVLS...")) rather than Tool("MVLS..."). Sequence, PDB, Ligand, ... are themselves tools that turn a bare value into a stream; delegating to them means you debug generic-stream handling once instead of re-implementing bare-value conversion in every tool.
    Shorthand inputs (an input where a bare STRING is meaningful, e.g. a ligand code) use the shared helper instead of the chain: resolve_basic_input(obj, cls, stream, argument) collapses StandardizedOutput / DataStream / bare-string into one line. 
    Operation-DSL variant (a different input shape): some tools take op_type + a list of small serializable Operation objects instead of plain params, so a user composes a sequence of steps - Template(operations=[Template.operation(...), Template.operation(...)]). Each Operation carries its own to_dict for serialization into the config. Use it when the tool applies an ordered, user-composed set of transforms rather than fixed parameters. Exemplars: panda.py (Panda.concat/sort/head), pdb.py, selection.py.
    """
    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 mode: str = "default",
                 **kwargs):
        self.mode = mode
        self.structures = structures  # raw handle kept for missing-propagation
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")
        super().__init__(**kwargs)

    """PARAMETERS VALIDATION - fail fast and loud at pipeline-build time for bad inputs, before a cluster job is queued. Raise ValueError with an actionable message. Check emptiness, mutually-exclusive options, out-of-range numbers here.
    """
    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")

    # configure_inputs (ABSTRACT - required) - receive the pipeline's resolved folders (install paths, container images, caches) and stash what you need. Most tools just keep self.folders; container tools rely on it so container_prefix() can find the image.
    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    # get_config_display (optional override) - one line per notable setting, printed in the pipeline summary and the tool-metadata JSON. Always start from super().get_config_display(). Keep it to human-meaningful facts (modes, key params), not every kwarg.
    # ⚠ Runs at build time: do not do work here. Printing len(self.structures_stream) forces the stream to materialize just for a count - report static params instead.
    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"MODE: {self.mode}")
        return lines

    """SCRIPT GENERATION
    Emit the bash that runs the tool. The spine below is identical for ~89 tools:
    - header(): completion-check guard: if a COMPLETED marker is present, it skips the execution
    - activate_environment(): activates the configured environment.
    - <your command>
    - footer(): success marker the harness polls. Generates COMPLETED or FAILED markers.
    Style: keep it a short orchestrator. Consider spliting consecutive phases into a private _generate_script_<phase>(). Exemplar: protein_mpnn.py.
    Serialize inputs first: dump upstream streams to the JSON the pipe_script reads (save_json); the pipe_script re-loads them out-of-env.

    PASSING PARAMETERS - two styles, pick one:
    (1) config.yaml (preferred when there are more than a couple of params): dump the tool's parameters to a yaml under configuration_path and pass only that path; the pipe reads it. Keeps the command short and gives an editable, inspectable record of the run. Exemplars: boltz2.py, boltzgen.py, diffdock.py. This is the style shown below.
    (2) explicit --flags: pass each parameter as its own flag. Fine for few params. Whichever you pick, every value the wrapper sends MUST be read by the pipe (a flag needs a matching add_argument; a yaml key needs to be read).
    Paths (dirs, map csvs) are still passed as flags either way - they are runtime plumbing, not parameters.

    Container: container_prefix() returns "" for env tools and the `apptainer/singularity exec ... <image> ` prefix when a container is configured, so the same code runs both ways. Pick EXACTLY ONE of two architectures - never both, or you nest containers:
    (a) the heavy command is emitted here in bash and IS the containerized thing > prefix it directly, f"{self.container_prefix()}python run_inference.py ..." and do NOT pass --container-prefix anywhere (rfdiffusion.py, protein_mpnn.py).
    (b) a host helper (pipe script) dispatches the binary > run the helper UNPREFIXED (it imports biopipelines and must stay on the host), pass the prefix down as --container-prefix, and in the pipe build the call as container_argv_prefix(prefix) + cmd so only the binary enters the .sif (dssp.py / pipe_dssp.py; also xtb, fpocket, p2rank, reduce, apbs). This template uses (b).
    A tool that runs the model truly in-process can do neither: call self.warn_container_unsupported().

    SKIPPING IDS - three cases, all end up in missing.csv so pipe_check_completion can tell "correctly absent" from "silently lost". The `kind` column is a CONTROLLED value that decides excused-vs-flagged (kind="filter" is excused, absence expected; kind="failure" is NOT excused, absent output is flagged as a real problem). The human reason goes in the `cause` column, never in `kind`. For streams we have:
    - upstream-filtered: a declared id whose file is absent (a filtered Pool/Panda still declares every id but carries only survivors). The pipe's iterate_files() SKIPS these automatically, so the loop is already a view over surviving entities - nothing to write; they are carried by the upstream manifest and propagated.
    - local failure: an id whose processing raises. The pipe's try/except records kind="failure" (NOT excused).
    - deliberate drop: an id THIS tool chooses to exclude (below a threshold, wrong format, ...). Record kind="filter" (excused) and put the reason in cause.
    For tables: an id that failed or was dropped still gets its row in every result table this tool emits, with the unavailable values as NaN. The table records that the id was processed, the missing manifest records why it has no entity. Do not filter a table on a metric threshold: emit it complete and let the user apply Panda.filter, keeping selection an explicit pipeline step.
    (In bash, to sample one PRESENT id without tripping on a filtered one, use Resolve.stream_ids(..., valid_set=True) - see gnina.py autobox.)

    MISSING PROPAGATION - two patterns, pick one:
    (A) propagation step (shown below): the pipe writes only its OWN failures/drops to a separate local file (execution/local_missing.csv), then generate_missing_propagation(*inputs, local_missing=..., missing_csv=...) emits a second step that merges upstream manifests + that local file into tables/missing.csv. Modular; the pipe never touches upstream tables. Use local_missing=self.local_missing_csv so those rows survive the merge; drop it if the tool never skips per-id.
    (B) all-in-one: the pipe itself reads the upstream missing table (passed as --upstream-missing) and writes the merged tables/missing.csv in one pass, with no separate propagation step. Exemplars: dssp.py, gnina.py.
    Either way ids dropped upstream stay excused instead of counting as failures. Do NOT do both.
    """
    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        self._write_config_yaml()
        script = "#!/bin/bash\n"
        script += "# TemplateTool script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        # Architecture (b): the helper is a HOST process (it imports biopipelines), so it is NOT prefixed.
        # The prefix is passed down; the helper wraps only its binary subprocess with container_argv_prefix.
        script += f"""echo "Running TemplateTool"
python "{self.helper_py}" \\
    --config-yaml "{self.config_yaml}" \\
    --structures-json "{self.structures_json}" \\
    --annotated-dir "{self.stream_folder('annotated')}" \\
    --annotated-map-csv "{self.annotated_map}" \\
    --counts-csv "{self.counts_csv}" \\
    --local-missing-csv "{self.local_missing_csv}" \\
    --container-prefix "{self.container_prefix()}"
"""
        # Pattern (A): merge upstream manifests + this tool's local failures into tables/missing.csv.
        script += self.generate_missing_propagation(
            self.structures, local_missing=self.local_missing_csv, missing_csv=self.missing_csv
        )
        script += self.generate_completion_check_footer()
        return script

    def _write_config_yaml(self):
        """Dump the tool's parameters to configuration/config.yaml for the pipe to read."""
        import yaml
        with open(self.config_yaml, "w") as f:
            yaml.safe_dump({"mode": self.mode}, f, sort_keys=False)

    """OUTPUT ID SHAPE - two DISTINCT mechanisms, only if the tool is not 1:1 (this file is 1:1: output ids == input ids)

    (A) COMBINATORICS - combining MULTIPLE input axes with each other (e.g. N proteins x M ligands). Wrap inputs with Each (cartesian product, the default for bare inputs) or Bundle (group several inputs into one entity). generate_combinatorics_config writes the pairing to JSON for the pipe_script; predict_output_ids_with_provenance precomputes the combined ids + provenance:
      predict_output_ids_with_provenance(proteins=(prot_src,"sequences"), ligands=(lig_src,"compounds"))
        -> ids ["prot1+lig1","prot1+lig2","prot1+lig3","prot2+lig1",...]
           provenance {"proteins":["prot1","prot1","prot1","prot2",...], "ligands":["lig1","lig2","lig3","lig1",...]}
      Multi-axis ids join components with "+". Exemplars: boltz2.py, gnina.py, esmfold2.py.

    (B) ID MULTIPLICATION - ONE input axis fanning out to MANY outputs (e.g. num_sequences=3 per structure, a sweep). generate_multiplied_ids appends suffixes and records each child's parent:
      generate_multiplied_ids(["prot_1","prot_2"], ["1","2","3"], input_stream_name="structures")
        -> ids ["prot_1_1","prot_1_2","prot_1_3","prot_2_1",...]
           provenance {"structures":["prot_1","prot_1","prot_1","prot_2",...]}
      Suffixes join with "_" (reserved for parent->child, distinct from combinatorics "+"). For unexpanded patterns, generate_multiplied_ids_pattern(["5HG6_<0..4>"], "<1..3>") -> ["5HG6_<0..4>_<1..3>"]. Exemplars: protein_mpnn.py, ligand_mpnn.py.

    Feed the resulting ids into DataStream(ids=...) below; provenance lets downstream get_mapped_ids join back to any parent axis.
    """

    """OUTPUT PREDICTION
    Declare the tool's outputs so downstream tools can wire to them before anything runs. This is the CONFIG-TIME half of the contract (see STREAM MODEL block above); the return shape must be stable since ids are known from inputs. The pipe writes the actual map_table at runtime - that is the source of truth.
    For each stream choose one of the three shapes from STREAM MODEL (per-id files / shared file / value-based) and set files= accordingly.
    Stream fields: name, ids (from the input stream), files (per STREAM MODEL), map_table (the stream_map_path descriptor), and a `format` tag consumers key on ("pdb", "sdf", "fasta", "smiles", "csv", ...).
    Multiple formats: a stream can advertise several formats as a pipe-delimited string, e.g. format="pdb|cif". Consumers read them via stream.formats (tuple), stream.has_format("pdb"), and stream.has_only_formats("pdb", "cif") - the last returns True only if every advertised format is in the allowed set. Use has_only_formats in validate_params to reject inputs you can't handle.
    Table: TableInfo(name, path, columns, description): list every column expected.
    Missing table: include it whenever the tool emits at least one stream (safe to include always). Columns are fixed: id | removed_by | kind | cause. It lets the framework distinguish "excused upstream" from "this tool failed".
    ⚠ Return "output_folder": self.output_folder.
    """
    def get_output_files(self) -> Dict[str, Any]:
        annotated_stream = DataStream(
            name="annotated",
            ids=self.structures_stream.ids,
            files=[self.stream_path("annotated", "<id>.txt")],
            map_table=self.annotated_map,
            format="text",
        )
        tables = {
            "counts": TableInfo(
                name="counts",
                path=self.counts_csv,
                columns=["id", "n_residues"],
                description="Residue count per input structure",
            ),
            "missing": self.missing_table_info(self.missing_csv),  # standard id|removed_by|kind|cause TableInfo
        }
        return {
            "annotated": annotated_stream, # streams individually
            "tables": tables, # tables as one dictionarry
            "output_folder": self.output_folder, # always here
        }
