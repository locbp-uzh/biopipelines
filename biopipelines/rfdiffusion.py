# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
RFdiffusion configuration for protein backbone generation.

RFdiffusion is a diffusion-based generative model for designing protein backbones.
It supports unconditional generation, motif scaffolding, binder design, and partial
diffusion. See https://github.com/RosettaCommons/RFdiffusion for full documentation.
"""

import os
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .biopipelines_io import Resolve, TableReference
    from .combinatorics import generate_multiplied_ids_pattern
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from biopipelines_io import Resolve, TableReference
    from combinatorics import generate_multiplied_ids_pattern


def _normalize_selection_arg(name, value):
    """Normalize a per-PDB selection argument to a (kind, token) pair.

    RFdiffusion's ``contigs`` / ``inpaint`` / ``inpaint_str`` each accept either
    a plain string — broadcast to every input PDB — or a per-PDB column
    reference (a ``TableReference`` from ``tool.tables.X.col``, or a
    ``(TableInfo, "col")`` tuple). A reference is resolved per-id at runtime via
    ``Resolve.table_column``; a literal is emitted directly.

    Returns ``("literal", str)`` or ``("table", "TABLE_REFERENCE:path:col")``.
    A literal is the only kind validated by ``_validate_freeform_string`` (a
    table reference reaches bash only through a resolved subshell, not as raw
    interpolation).
    """
    if isinstance(value, TableReference):
        return ("table", str(value))
    # (TableInfo, "column") tuple — mirror the LigandMPNN/ThermoMPNN convention.
    if isinstance(value, tuple) and len(value) == 2 and isinstance(value[1], str):
        return ("table", str(TableReference(value[0].info.path, value[1])))
    if isinstance(value, str):
        return ("literal", value)
    raise ValueError(
        f"{name} must be a string or a (table, column) reference, got {type(value)}"
    )


def _selection_cli_arg(kind_token) -> str:
    """CLI arg for the contigs resolver: the token, or '-' when empty.

    ``kind_token`` is the ``(kind, token)`` pair from
    ``_normalize_selection_arg``. Shared by all three RFdiffusion tools.
    """
    kind, token = kind_token
    if kind == "literal" and not token:
        return "-"
    return token


class RFdiffusion(BaseConfig):
    """
    Configuration for RFdiffusion protein backbone generation.

    RFdiffusion generates novel protein backbones using a denoising diffusion
    probabilistic model conditioned on structural motifs, partial structures,
    or protein-protein interaction targets.

    Main workflows:
        - Unconditional generation: produce backbones of a given length range
          without any input structure (contigs only, no pdb).
        - Motif scaffolding: hold a fixed structural motif in place and design
          the surrounding scaffold (pdb + contigs specifying the motif).
        - Binder design: design a new protein that binds to a target (pdb of
          the target + contigs describing the binder length).
        - Partial diffusion: apply limited noise to an existing structure and
          re-diffuse, producing near-neighbour variants (partial_steps).

    Model checkpoints (see WEIGHTS / DEFAULT_WEIGHTS):
        By default only Base and Complex_base are downloaded (~2 GB). Add others
        via RFdiffusion.install(weights=[...]) when needed.

    Installation:
        with Pipeline(...):
            RFdiffusion.install()                         # default weights
            RFdiffusion.install(weights=["Base",
                                         "InpaintSeq"])   # custom selection
            rfd = RFdiffusion(contigs="50-100", ...)

    Reference:
        Watson et al. (2023) De novo design of protein structure and function
        with RFdiffusion. Nature 620, 1089-1100.
        https://github.com/RosettaCommons/RFdiffusion
    """

    TOOL_NAME = "RFdiffusion"
    TOOL_VERSION = "1.0"

    # Mapping of weight name -> (url_hash, filename)
    WEIGHTS = {
        "Base":              ("6f5902ac237024bdd0c176cb93063dc4", "Base_ckpt.pt"),
        "Complex_base":      ("e29311f6f1bf1af907f9ef9f44b8328b", "Complex_base_ckpt.pt"),
        "Complex_Fold_base": ("60f09a193fb5e5ccdc4980417708dbab", "Complex_Fold_base_ckpt.pt"),
        "InpaintSeq":        ("74f51cfb8b440f50d70878e05361d8f0", "InpaintSeq_ckpt.pt"),
        "InpaintSeq_Fold":   ("76d00716416567174cdb7ca96e208296", "InpaintSeq_Fold_ckpt.pt"),
        "ActiveSite":        ("5532d2e1f3a4738decd58b19d633b3c3", "ActiveSite_ckpt.pt"),
        "Base_epoch8":       ("12fc204edeae5b57713c5ad7dcb97d39", "Base_epoch8_ckpt.pt"),
        "RF_structure_prediction": ("1befcb9b28e2f778f53d47f18b7597fa", "RF_structure_prediction_weights.pt"),
    }
    # Weights downloaded by default — covers the two main workflows:
    #   Base          → unconditional generation, motif scaffolding
    #   Complex_base  → binder design / protein-protein interaction
    # Add more via RFdiffusion.install(weights=[...]):
    #   Complex_Fold_base  → binder design with fold/topology conditioning
    #   InpaintSeq         → auto-selected when contigmap.inpaint_seq is used
    #   InpaintSeq_Fold    → inpaint + fold conditioning together
    #   ActiveSite         → small active-site motifs (requires ckpt_override_path)
    #   Base_epoch8        → alternative base checkpoint
    #   RF_structure_prediction → legacy structure prediction (unrelated to diffusion)
    DEFAULT_WEIGHTS = ["Base", "Complex_base"]

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False,
                        weights=None, **kwargs):
        """
        Generate the bash installation script for RFdiffusion.

        Clones the RFdiffusion repository, downloads the requested model
        checkpoints, and creates (or reuses) the SE3nv conda/mamba environment.

        Args:
            folders: Resolved pipeline folder paths (must contain "RFdiffusion"
                     and "biopipelines" keys).
            env_manager: "mamba", "conda", or "micromamba".
            force_reinstall: If True, skip the already-installed early-exit check.
            weights: List of checkpoint names to download. Defaults to
                     DEFAULT_WEIGHTS (["Base", "Complex_base"]).
                     Valid names are the keys of WEIGHTS.
        """
        biopipelines = folders.get("biopipelines", "")
        repo_dir = folders.get("RFdiffusion", "")
        parent_dir = os.path.dirname(repo_dir)

        if weights is None:
            weights = cls.DEFAULT_WEIGHTS
        invalid = [w for w in weights if w not in cls.WEIGHTS]
        if invalid:
            raise ValueError(
                f"Unknown weight(s): {invalid}. "
                f"Valid options are: {list(cls.WEIGHTS.keys())}"
            )
        wget_lines = "\n".join(
            f"wget -nc http://files.ipd.uw.edu/pub/RFdiffusion/{cls.WEIGHTS[w][0]}/{cls.WEIGHTS[w][1]}"
            f" || echo \"WARNING: failed to download {cls.WEIGHTS[w][1]}, skipping\""
            for w in weights
        )

        # Skip only when the env exists AND `import rfdiffusion` actually
        # works inside it. Just the env name being present is not enough:
        # ProteinMPNN.install() may have created a leaner SE3nv (PyTorch
        # only) before us, in which case we must layer DGL, SE3Transformer,
        # and the rfdiffusion package on top instead of skipping.
        # The import probe is the load-bearing check (env + package); add the
        # weights-file test on top so we don't skip when checkpoints are missing.
        skip = "" if force_reinstall else f"""# Check if already installed
if [ -d "{repo_dir}/models" ] && [ -f "{repo_dir}/models/Base_ckpt.pt" ] \\
   && {env_manager} run -n SE3nv python -c "import rfdiffusion" >/dev/null 2>&1; then
    echo "RFdiffusion already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block("SE3nv", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("SE3nv", env_manager, biopipelines)
        return f"""echo "=== Installing RFdiffusion ==="
{skip}mkdir -p {parent_dir}
cd {parent_dir}
if [ ! -d "{repo_dir}" ]; then
    git clone https://github.com/RosettaCommons/RFdiffusion.git
fi
cd {repo_dir}

# Download model weights
mkdir -p models && cd models
{wget_lines}
cd ..

# Create SE3nv environment from BioPipelines specification
{remove_block}
{env_block}
if [ $? -ne 0 ]; then
    echo "WARNING: BioPipelines SE3nv env creation failed. Trying official RFdiffusion environment..."
    {env_manager} env create -f env/SE3nv.yml
    if [ $? -ne 0 ]; then
        echo "ERROR: SE3nv environment creation failed with both methods."
        echo "This is likely a CUDA version mismatch for your system."
        exit 1
    fi
fi

# Install DGL from pre-built wheel (special --find-links syntax)
{env_manager} run -n SE3nv pip install dgl -f https://data.dgl.ai/wheels/torch-2.4/cu124/repo.html --no-deps

# Install SE3Transformer and RFdiffusion (editable)
cd env/SE3Transformer
{env_manager} run -n SE3nv pip install --no-deps .
cd ../..
{env_manager} run -n SE3nv pip install -e .

# Verify installation
if {env_manager} run -n SE3nv python -c "import rfdiffusion" >/dev/null 2>&1 || [ -f "{repo_dir}/models/Base_ckpt.pt" ]; then
    touch "$INSTALL_SUCCESS"
    echo "=== RFdiffusion installation complete ==="
else
    echo "ERROR: RFdiffusion verification failed"
    exit 1
fi
"""

    # Lazy path descriptors
    #   main_table        — TableInfo CSV describing each design's provenance
    #                       (fixed/designed regions, pLDDT). Lives under tables/.
    #   table_py_file     — pipe script that builds main_table from .trb files.
    #   inference_py_file — RFdiffusion's own CLI entry point.
    #   pdb_ds_json       — serialized input DataStream (config-time artifact).
    #   update_map_py     — pipe script that rewrites structures_map.csv to
    #                       match actual files on disk after the run.
    main_table = Path(lambda self: self.table_path("structures"))
    table_py_file = Path(lambda self: self.pipe_script_path("pipe_rfdiffusion_table.py"))
    inference_py_file = Path(lambda self: os.path.join(self.folders["RFdiffusion"], "scripts", "run_inference.py"))
    pdb_ds_json = Path(lambda self: self.configuration_path("input_structures.json"))
    update_map_py = Path(lambda self: self.pipe_script_path("pipe_update_structures_map.py"))
    #   contigs_json       — per-PDB resolved {contigs, inpaint, inpaint_str},
    #                        written at runtime by pipe_rfdiffusion_contigs.py.
    #   contigs_py         — that resolver (one pass over all ids → JSON).
    #   contigs_reader_py  — per-id reader the bash loop calls for each PDB.
    contigs_json = Path(lambda self: self.configuration_path("contig_options.json"))
    contigs_py = Path(lambda self: self.pipe_script_path("pipe_rfdiffusion_contigs.py"))
    contigs_reader_py = Path(lambda self: self.pipe_script_path("resolve_rfdiffusion_contigs.py"))

    def __init__(self,
                 contigs: Union[str, TableReference, "tuple"],
                 pdb: Optional[Union[DataStream, StandardizedOutput]] = None,
                 inpaint: Union[str, TableReference, "tuple"] = "",
                 inpaint_str: Union[str, TableReference, "tuple"] = "",
                 num_designs: int = 1,
                 active_site: bool = False,
                 steps: int = 50,
                 partial_steps: int = 0,
                 reproducible: bool = False,
                 design_startnum: int = 1,
                 **kwargs):
        """
        Initialize RFdiffusion configuration.

        Args:
            contigs: Contig map describing which residues are fixed and which
                     are to be generated. Residues from the input PDB are
                     specified as chain+range (e.g. "A1-100"); new backbone
                     segments are specified as length ranges (e.g. "50-100").
                     Multiple segments are separated by "/" (e.g. "A1-50/30-50/A60-100").
                     For unconditional generation (no pdb) use a length range alone
                     (e.g. "100-200"). Either a plain string (broadcast to every
                     input PDB) or a per-PDB column reference
                     (``tool.tables.X.col`` / ``(TableInfo, "col")``) resolved by
                     id match at runtime — use the latter to give each input PDB
                     its own contigs. A column reference requires ``pdb``.
            pdb: Optional input structure(s). Required for motif scaffolding,
                 binder design, and partial diffusion. Accepts a DataStream or a
                 StandardizedOutput. When it carries multiple structures,
                 RFdiffusion runs once per input PDB (``num_designs`` designs
                 each); output ids are ``<pdb_id>_<n>``.
            inpaint: Residues whose sequence should be masked during diffusion
                     (same chain+range format as contigs). When set, the
                     InpaintSeq checkpoint is used automatically — ensure it
                     was downloaded at install time. Accepts the same
                     string-or-column-reference forms as ``contigs``.
            inpaint_str: Residues whose secondary structure should be masked
                     during diffusion (same chain+range format as contigs).
                     Wired to upstream ``contigmap.inpaint_str``. Empty string
                     disables (default). The overall length stays governed by
                     the ``contigs`` argument (use a range like ``"50-100"``
                     to control it). Accepts the same string-or-column-reference
                     forms as ``contigs``.
            num_designs: Number of independent backbone designs to generate.
            active_site: If True, use the ActiveSite checkpoint instead of
                         Base. Intended for scaffolding very small functional
                         motifs (< ~10 residues). Requires the ActiveSite
                         weight to have been downloaded at install time.
            steps: Number of denoising diffusion steps (default 50). Fewer
                   steps are faster but may reduce quality.
            partial_steps: Number of partial diffusion steps. When > 0, the
                           input structure is noised for this many steps and
                           then re-diffused, producing near-neighbour variants.
                           Must be less than steps.
            reproducible: If True, use deterministic (fixed-seed) sampling so
                          the same contigs always produce the same outputs.
            design_startnum: Integer appended to the pipeline name to number
                             output files (e.g. design_startnum=1 → name_1.pdb).
                             Useful when continuing a previous run.

        Output:
            Streams: structures (.pdb)
            Tables:
                structures: id | pdb | fixed | designed | source_fixed | plddt_mean | status
        """
        # Resolve optional pdb input — store stream for runtime resolution.
        # All input PDBs are iterated at execution time (see generate_script);
        # never index ids[0] at config time, which breaks under lazy IDs.
        self.pdb_stream: Optional[DataStream] = None
        if pdb is not None:
            if isinstance(pdb, StandardizedOutput):
                self.pdb_stream = pdb.streams.structures
            elif isinstance(pdb, DataStream):
                self.pdb_stream = pdb
            else:
                raise ValueError(f"pdb must be DataStream or StandardizedOutput, got {type(pdb)}")

        # Normalize the per-PDB selection args to (kind, token) pairs. A literal
        # is broadcast to every PDB; a table reference is resolved per-id at
        # runtime. Store both the original (for display) and the normalized form.
        self.contigs = contigs
        self.inpaint = inpaint
        self.inpaint_str = inpaint_str
        self._contigs_arg = _normalize_selection_arg("contigs", contigs)
        self._inpaint_arg = _normalize_selection_arg("inpaint", inpaint)
        self._inpaint_str_arg = _normalize_selection_arg("inpaint_str", inpaint_str)
        self.num_designs = num_designs
        self.active_site = active_site
        self.steps = steps
        self.partial_steps = partial_steps
        self.reproducible = reproducible
        self.design_startnum = design_startnum

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate RFdiffusion-specific parameters."""
        contigs_kind, contigs_token = self._contigs_arg
        if contigs_kind == "literal" and not contigs_token:
            raise ValueError("contigs parameter is required")

        if self.num_designs <= 0:
            raise ValueError("num_designs must be positive")

        if self.steps <= 0:
            raise ValueError("steps must be positive")

        if self.partial_steps < 0:
            raise ValueError("partial_steps cannot be negative")

        # A column reference is keyed by input-PDB id, so it needs PDBs to
        # match against.
        for name, (kind, _) in (("contigs", self._contigs_arg),
                                 ("inpaint", self._inpaint_arg),
                                 ("inpaint_str", self._inpaint_str_arg)):
            if kind == "table" and self.pdb_stream is None:
                raise ValueError(
                    f"{name} given as a table column reference requires a pdb input"
                )

        # Only literal selections are interpolated raw into bash; table
        # references reach the shell through a resolved subshell.
        for name, (kind, token) in (("contigs", self._contigs_arg),
                                     ("inpaint", self._inpaint_arg),
                                     ("inpaint_str", self._inpaint_str_arg)):
            if kind == "literal":
                _validate_freeform_string(name, token)

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input files."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get RFdiffusion configuration display lines."""
        config_lines = super().get_config_display()
        config_lines.extend([
            f"CONTIGS: {self.contigs}",
            f"NUM DESIGNS: {self.num_designs}",
            f"ACTIVE SITE: {self.active_site}",
            f"STEPS: {self.steps}"
        ])

        if self.pdb_stream:
            config_lines.append(f"PDB: {', '.join(self.pdb_stream.ids)}")
        if self._inpaint_arg[1]:
            config_lines.append(f"INPAINT: {self.inpaint}")
        if self._inpaint_str_arg[1]:
            config_lines.append(f"INPAINT_STR: {self.inpaint_str}")
        if self.partial_steps > 0:
            config_lines.append(f"PARTIAL STEPS: {self.partial_steps}")
        if self.reproducible:
            config_lines.append(f"REPRODUCIBLE: {self.reproducible}")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate RFdiffusion execution script."""
        # Serialize input DataStream to JSON for runtime file resolution.
        # The configuration/ folder was already created by the pipeline
        # after get_output_files() returned.
        if self.pdb_stream:
            self.pdb_stream.save_json(self.pdb_ds_json)

        script_content = "#!/bin/bash\n"
        script_content += "# RFdiffusion execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        # e3nn 0.3.3 uses torch.load() without weights_only=False,
        # which fails on PyTorch 2.6+ where the default flipped to True
        script_content += "export TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1\n" 
        script_content += self._generate_script_run_rfdiffusion()
        script_content += self._generate_script_create_table()
        script_content += self._generate_script_update_structures_map()
        script_content += self.generate_completion_check_footer()
        return script_content

    def _common_rfd_options(self) -> str:
        """Inference options shared by every input PDB (no contigs/pdb/prefix)."""
        opts = f"inference.num_designs={self.num_designs}"
        opts += f" inference.deterministic={self.reproducible}"
        opts += f" inference.design_startnum={self.design_startnum}"
        if self.steps != 50:
            opts += f" diffuser.T={self.steps}"
        if self.partial_steps > 0:
            opts += f" diffuser.partial_T={self.partial_steps}"
        if self.active_site:
            opts += " inference.ckpt_override_path=models/ActiveSite_ckpt.pt"
        return opts

    def _generate_script_run_rfdiffusion(self) -> str:
        """Generate the RFdiffusion execution part of the script.

        With a PDB input, RFdiffusion runs once per input structure inside a
        bash loop over the input ids (resolved at runtime — never ids[0] at
        config time, which breaks under lazy IDs). The per-PDB contigs /
        inpaint / inpaint_str selections (literal broadcast, or a table-column
        reference) are pre-resolved in one pass into contig_options.json by
        pipe_rfdiffusion_contigs.py, and the loop reads each id's values from
        it. Without a PDB it is a single unconditional run named after the
        pipeline.
        """
        common = self._common_rfd_options()
        structures_dir = self.stream_folder("structures")

        if not self.pdb_stream:
            # Unconditional generation — single run, contigs is a literal
            # (validate_params already rejected a table ref without pdb).
            contigs = self._contigs_arg[1]
            rfd_options = f"'contigmap.contigs=[{contigs}]'"
            if self._inpaint_arg[1]:
                rfd_options += f" 'contigmap.inpaint_seq=[{self._inpaint_arg[1]}]'"
            if self._inpaint_str_arg[1]:
                rfd_options += f" 'contigmap.inpaint_str=[{self._inpaint_str_arg[1]}]'"
            prefix = os.path.join(structures_dir, self.pipeline_name)
            rfd_options += f" inference.output_prefix={prefix} {common}"
            return f"""echo "Starting RFdiffusion (unconditional)"
echo "Output folder: {self.output_folder}"

cd {self.folders["RFdiffusion"]}
{self.container_prefix()}python {self.inference_py_file} {rfd_options}

"""

        # PDB-conditioned: pre-resolve per-PDB selections (one pass over all
        # ids) into a JSON, then loop over the input ids at runtime.
        contigs_arg = _selection_cli_arg(self._contigs_arg)
        inpaint_arg = _selection_cli_arg(self._inpaint_arg)
        inpaint_str_arg = _selection_cli_arg(self._inpaint_str_arg)

        return f"""echo "Resolving per-PDB contig options"
python {self.contigs_py} "{self.pdb_ds_json}" "{contigs_arg}" "{inpaint_arg}" "{inpaint_str_arg}" "{self.contigs_json}"

echo "Starting RFdiffusion"
echo "Output folder: {self.output_folder}"
cd {self.folders["RFdiffusion"]}

for STRUCT_ID in {Resolve.stream_ids(self.pdb_ds_json)}; do
    INPUT_PDB={Resolve.stream_item(self.pdb_ds_json, '$STRUCT_ID')}
    OPTS=$(python {self.contigs_reader_py} "{self.contigs_json}" "$STRUCT_ID")
    CONTIGS=$(echo "$OPTS" | sed -n '1p')
    INPAINT_SEL=$(echo "$OPTS" | sed -n '2p')
    INPAINT_STR_SEL=$(echo "$OPTS" | sed -n '3p')
    RFD_OPTIONS="'contigmap.contigs=[$CONTIGS]'"
    if [ -n "$INPAINT_SEL" ]; then
        RFD_OPTIONS="$RFD_OPTIONS 'contigmap.inpaint_seq=[$INPAINT_SEL]'"
    fi
    if [ -n "$INPAINT_STR_SEL" ]; then
        RFD_OPTIONS="$RFD_OPTIONS 'contigmap.inpaint_str=[$INPAINT_STR_SEL]'"
    fi
    RFD_OPTIONS="$RFD_OPTIONS inference.input_pdb=$INPUT_PDB"
    RFD_OPTIONS="$RFD_OPTIONS inference.output_prefix={structures_dir}/$STRUCT_ID {common}"
    echo "Running RFdiffusion for $STRUCT_ID"
    eval {self.container_prefix()}python {self.inference_py_file} $RFD_OPTIONS
done

"""

    def _generate_script_create_table(self) -> str:
        """Generate the table creation part of the script."""
        # pipe_rfdiffusion_table.py scans the structures/ stream folder for every
        # produced .pdb/.trb pair — one prefix per input PDB in the multi-PDB
        # case, one for the unconditional case.
        structures_dir = self.stream_folder("structures")
        return f"""echo "Creating results table"
python {self.table_py_file} "{structures_dir}" "{self.main_table}"

"""

    def _generate_script_update_structures_map(self) -> str:
        """Generate script to write structures_map.csv from the actual runtime PDBs."""
        structures_map = self.stream_map_path("structures")
        structures_dir = self.stream_folder("structures")
        # Each design id is "<pdb_id>_<n>"; its parent PDB is the id minus the
        # trailing "_<n>" suffix. Derive the `structures.id` provenance column
        # from that suffix at runtime (works for one or many parent PDBs). With
        # no PDB input there is no parent, so no provenance column.
        prov_arg = ' --provenance-from-suffix "structures.id"' if self.pdb_stream else ""
        return f"""echo "Writing structures map from actual output files"
python {self.update_map_py} --structures-map "{structures_map}" --output-folder "{structures_dir}"{prov_arg}

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after RFdiffusion execution."""
        # Pattern-based IDs — PDBs land under <output_folder>/structures/.
        # The per-design map_table is written at runtime by
        # _generate_script_update_structures_map(); here we only declare the
        # stream and its map_table path.
        start = self.design_startnum
        end = self.design_startnum + self.num_designs - 1
        suffix_pattern = f"<{start}..{end}>"
        if self.pdb_stream:
            # One <pdb_id>_<n> fan-out per input PDB. Keep parent ids compact /
            # lazy-safe via the shared multiplier helper.
            structure_ids = generate_multiplied_ids_pattern(
                self.pdb_stream.ids, suffix_pattern,
                input_stream_name="structures"
            )
        else:
            structure_ids = [f"{self.pipeline_name}_{suffix_pattern}"]
        structures_dir = self.stream_folder("structures")
        file_template = [os.path.join(structures_dir, "<id>.pdb")]
        structures_map = self.stream_map_path("structures")

        structures = DataStream(
            name="structures",
            ids=structure_ids,
            files=file_template,
            map_table=structures_map,
            format="pdb"
        )

        tables = {
            "structures": TableInfo(
                name="structures",
                path=self.main_table,
                columns=["id", "pdb", "fixed", "designed", "source_fixed", "plddt_mean", "status"],
                description="RFdiffusion structure generation results"
            )
        }

        return {
            "structures": structures,
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "rfd_params": {
                "pdb_input_ids": self.pdb_stream.ids if self.pdb_stream else None,
                "contigs": str(self.contigs),
                "inpaint": str(self.inpaint),
                "inpaint_str": str(self.inpaint_str),
                "num_designs": self.num_designs,
                "active_site": self.active_site,
                "steps": self.steps,
                "partial_steps": self.partial_steps,
                "reproducible": self.reproducible,
                "design_startnum": self.design_startnum
            }
        })
        return base_dict
