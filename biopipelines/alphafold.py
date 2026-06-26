# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
AlphaFold2/ColabFold configuration for protein structure prediction.

Handles sequence folding with AlphaFold2 using ColabFold implementation,
integrating with upstream sequence generation tools and providing comprehensive
ranking and analysis capabilities.
"""

import os
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
    from .datastream_resolver import resolve_input_to_datastream
    from .combinatorics import (
        Bundle, Each, get_mode, contains_combinatorics_wrapper,
        generate_combinatorics_config, predict_output_ids_with_provenance,
    )
    from ._weights_cache import link_weights_block
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream
    from datastream_resolver import resolve_input_to_datastream
    from combinatorics import (
        Bundle, Each, get_mode, contains_combinatorics_wrapper,
        generate_combinatorics_config, predict_output_ids_with_provenance,
    )
    from _weights_cache import link_weights_block


class AlphaFold(BaseConfig):
    """
    Configuration for AlphaFold2/ColabFold structure prediction.

    Predicts protein structures from amino acid sequences.

    By default each input sequence is folded as a separate monomer. Wrap the
    ``proteins`` input in ``Bundle(...)`` to fold the bundled sequences together
    as one multi-chain complex: ColabFold receives a single colon-joined query
    (``SEQ_A:SEQ_B``), auto-selects the multimer model, and runs its default
    paired+unpaired MSA pipeline. ``Bundle(static, Each(...))`` holds ``static``
    fixed and folds it against each iterated sequence (one complex per element).

    Example:
        # Using Sequence tool — two independent monomers
        proteins = Sequence(["MKTVRQ...", "AETGFT..."], ids=["p1", "p2"])
        af = AlphaFold(proteins=proteins)

        # Using output from another tool
        af = AlphaFold(proteins=mpnn_output)

        # Fold a single A:B complex (one prediction, id "p1+p2")
        af = AlphaFold(proteins=Bundle(seqs_a, seqs_b))

        # One complex per binder: fixed receptor + each binder
        af = AlphaFold(proteins=Bundle(receptor, Each(binders)))
    """

    TOOL_NAME = "AlphaFold"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        try:
            from .config_manager import ConfigManager
        except ImportError:
            from config_manager import ConfigManager
        scheduler = ConfigManager().get_scheduler()

        repo_dir = folders.get("AlphaFold", "")
        parent_dir = os.path.dirname(repo_dir)
        if scheduler == "colab":
            skip = "" if force_reinstall else """# Check if already installed
if /usr/bin/python3 -c "import colabfold" 2>/dev/null; then
    echo "AlphaFold (ColabFold) already installed, skipping pip install. Use force_reinstall=True to reinstall."
else
"""
            skip_end = "" if force_reinstall else "fi\n"
            return f"""echo "=== Installing AlphaFold (ColabFold on Colab) ==="
{skip}# Install into Colab's base Python to reuse JAX/GPU stack
/usr/bin/python3 -m pip install -q --no-warn-conflicts "colabfold[alphafold-minus-jax] @ git+https://github.com/sokrypton/ColabFold"
{skip_end}# Fix TF crashes (from the official ColabFold notebook) — always run,
# even if pip install was skipped, because a parallel session may have
# installed colabfold without applying these fixes.
#   * libtfkernel_sobol_op.so — older Sobol-op crash on TF import.
#   * tensorflow/lite/python/*/*.so — newer C-extension symbol mismatch
#     in tflite (e.g. _pywrap_tensorflow_lite_metrics_wrapper.so) hit
#     when ColabFold falls back to `import tensorflow` after failing to
#     find the optional `tpu_info` module on regular GPU runtimes.
rm -f /usr/local/lib/python3.*/dist-packages/tensorflow/core/kernels/libtfkernel_sobol_op.so \\
      /usr/local/lib/python3.*/dist-packages/tensorflow/lite/python/*/*.so

# Verify installation
if /usr/bin/python3 -c "import colabfold" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== AlphaFold installation complete ==="
else
    echo "ERROR: AlphaFold verification failed (cannot import colabfold)"
    exit 1
fi
"""
        # Skip only if the install dir AND a working colabfold-conda env are
        # present — a bare dir (or one whose env was clobbered, e.g. by a
        # downstream tool upgrading jax) must not read as "installed".
        skip = "" if force_reinstall else f"""# Check if already installed
if [ -d "{repo_dir}" ] && [ -x "{repo_dir}/colabfold-conda/bin/colabfold_batch" ] \\
   && "{repo_dir}/colabfold-conda/bin/python" -c "import colabfold" >/dev/null 2>&1; then
    echo "AlphaFold (LocalColabFold) already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        # On force_reinstall, remove the bundled Miniconda and the conda env —
        # install_colabbatch_linux.sh's Miniconda installer aborts if its target
        # dir already exists. Keep colabfold/params so the weights aren't re-downloaded.
        remove_block = f"""rm -rf "{repo_dir}/conda" "{repo_dir}/colabfold-conda"
""" if force_reinstall else ""
        return f"""echo "=== Installing AlphaFold (LocalColabFold) ==="
{skip}{remove_block}cd {parent_dir}
# install_colabbatch_linux.sh builds colabfold_batch into
# <dir>/localcolabfold/colabfold-conda via its own bundled Miniconda
# (python 3.10, jax[cuda11_pip]). Run it unmodified — it sources only the
# `conda` shell function, so rewriting conda->mamba would break its internal
# `conda activate` and leak installs into the active env. The subshell runs it
# with no host env active (stray writes can't touch a tool env) and accepts
# conda's ToS non-interactively (else `conda create` aborts).
wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/v1.5.5/install_colabbatch_linux.sh
# Ensure pip in the env it creates: `conda create ... python=3.10` (conda-forge
# python ships no pip) is followed by `colabfold-conda/bin/pip install`, which
# would otherwise be missing. Add pip to the create spec.
sed -i 's/git python=3.10/git python=3.10 pip/' install_colabbatch_linux.sh
(
export CONDA_PLUGINS_AUTO_ACCEPT_TOS=true
{env_manager} deactivate 2>/dev/null || true
bash install_colabbatch_linux.sh
)
rm install_colabbatch_linux.sh

# Verify installation
if [ -d "{repo_dir}" ] && [ -x "{repo_dir}/colabfold-conda/bin/colabfold_batch" ]; then
    touch "$INSTALL_SUCCESS"
    echo "=== AlphaFold installation complete ==="
else
    echo "ERROR: AlphaFold verification failed (colabfold_batch not found)"
    exit 1
fi
"""

    # Lazy path descriptors — routed through the canonical sub-layout.
    #   configuration/  — queries CSV/FASTA passed to colabfold_batch.
    #   execution/      — raw ColabFold Folding dumps.
    #   structures/     — best-rank PDBs + structures_map.csv.
    #   msas/           — per-ID .a3m files + msas_map.csv.
    #   tables/         — confidence, missing.
    queries_csv = Path(lambda self: self.configuration_path(f"{self.pipeline_name}_queries.csv"))
    queries_fasta = Path(lambda self: self.configuration_path(f"{self.pipeline_name}_queries.fasta"))
    combinatorics_config_file = Path(lambda self: self.configuration_path("combinatorics_config.json"))
    sequences_json = Path(lambda self: self.configuration_path(".input_sequences.json"))
    confidence_csv = Path(lambda self: self.table_path("confidence"))
    folding_folder = Path(lambda self: self.execution_folder)
    msas_folder = Path(lambda self: self.stream_folder("msas"))
    msa_csv = Path(lambda self: self.stream_map_path("msas"))
    structures_map = Path(lambda self: self.stream_map_path("structures"))
    missing_csv = Path(lambda self: self.table_path("missing"))

    # Helper script paths
    colabfold_batch = Path(lambda self: os.path.join(self.folders["AlphaFold"], "colabfold-conda/bin/colabfold_batch"))
    fa_to_csv_fasta_py = Path(lambda self: self.pipe_script_path("pipe_fa_to_csv_fasta.py"))
    alphafold_confidence_py = Path(lambda self: self.pipe_script_path("pipe_alphafold_confidence.py"))
    alphafold_msas_py = Path(lambda self: self.pipe_script_path("pipe_alphafold_msas.py"))
    alphafold_queries_py = Path(lambda self: self.pipe_script_path("pipe_alphafold_queries.py"))
    update_map_py = Path(lambda self: self.pipe_script_path("pipe_update_structures_map.py"))
    unsanitize_ids_py = Path(lambda self: self.pipe_script_path("pipe_unsanitize_ids.py"))
    msa_copy_py = Path(lambda self: self.pipe_script_path("pipe_boltz_msa_copy.py"))

    def __init__(self,
                 proteins: Union[DataStream, StandardizedOutput],
                 msas: Optional[StandardizedOutput] = None,
                 num_relax: int = 0,
                 num_recycle: int = 3,
                 rand_seed: int = 0,
                 **kwargs):
        """
        Initialize AlphaFold configuration.

        Args:
            proteins: Input protein sequences as DataStream or StandardizedOutput.
                      Use Sequence("MKTVRQ...") to create from raw sequence strings.
            msas: Pre-computed MSAs as StandardizedOutput with msas stream in A3M format.
                  When provided, these MSAs are copied into the Folding directory so
                  ColabFold skips MSA generation for matching queries.
                  Use MSA(source, convert="a3m") to convert from CSV if needed.
            num_relax: Number of best models to relax with AMBER
            num_recycle: Number of recycling iterations (default 3)
            rand_seed: Random seed for reproducible results (0 = random)

        Output:
            Streams: structures (.pdb), msas (.a3m)
            Tables:
                structures: id | file
                confidence: id | structure | plddt | max_pae | ptm
                msas: id | sequences.id | sequence | msa_file
                missing: id | removed_by | kind | cause
        """
        # Store original input for upstream missing table lookup and
        # combinatorics config generation (may be Bundle/Each-wrapped).
        self.proteins = proteins

        # When proteins is wrapped in Bundle/Each, sequences are grouped into
        # complexes via a combinatorics config + a queries-builder pipe script.
        # A bare input keeps the legacy one-monomer-per-sequence path.
        self._uses_combinatorics = contains_combinatorics_wrapper(proteins)
        self._proteins_mode = get_mode(proteins)

        # Resolve input to a DataStream for the msas/missing-table machinery.
        # resolve_input_to_datastream unwraps Bundle/Each to the sequences stream.
        if isinstance(proteins, StandardizedOutput):
            self.sequences_stream: DataStream = proteins.streams.sequences
        elif isinstance(proteins, DataStream):
            self.sequences_stream = proteins
        elif self._uses_combinatorics:
            self.sequences_stream = resolve_input_to_datastream(proteins, fallback_stream="sequences")
        else:
            raise ValueError(
                f"proteins must be DataStream, StandardizedOutput, or Bundle/Each, got {type(proteins)}"
            )

        # Pre-computed MSAs for recycling
        self.msas_input = msas
        if msas is not None:
            self.msas_stream_input = resolve_input_to_datastream(msas, fallback_stream="msas")
        else:
            self.msas_stream_input = None

        # Store AlphaFold-specific parameters
        self.num_relax = num_relax
        self.num_recycle = num_recycle
        self.rand_seed = rand_seed

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate AlphaFold-specific parameters."""
        if not self.sequences_stream or len(self.sequences_stream) == 0:
            raise ValueError("proteins parameter is required and must not be empty")

        if self.msas_stream_input is not None:
            fmt = self.msas_stream_input.format
            if fmt != "a3m":
                raise ValueError(
                    f"AlphaFold requires MSAs in A3M format, got '{fmt}'. "
                    f"Use MSA(source, convert=\"a3m\") to convert first."
                )
            if self._proteins_mode == "bundle":
                raise ValueError(
                    "Pre-computed msas= are not supported together with a Bundle "
                    "(multi-chain complex). Bundled complexes use ColabFold's own "
                    "paired+unpaired MSA pipeline; drop msas= for the bundle, or "
                    "fold the chains as separate monomers to reuse precomputed MSAs."
                )

        if self.num_relax < 0:
            raise ValueError("num_relax cannot be negative")

        if self.num_recycle < 1:
            raise ValueError("num_recycle must be at least 1")

        if self.rand_seed < 0:
            raise ValueError("rand_seed cannot be negative")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input sequences."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get AlphaFold configuration display lines."""
        config_lines = super().get_config_display()
        config_lines.extend([
            f"INPUT PROTEINS: {len(self.sequences_stream)} sequences",
            f"NUM RELAX: {self.num_relax}",
            f"NUM RECYCLE: {self.num_recycle}"
        ])

        if self.msas_stream_input is not None:
            config_lines.append(f"PRE-COMPUTED MSAs: {len(self.msas_stream_input)} files ({self.msas_stream_input.format})")

        if self.rand_seed > 0:
            config_lines.append(f"RAND SEED: {self.rand_seed}")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate AlphaFold execution script."""
        script_content = "#!/bin/bash\n"
        script_content += "# AlphaFold execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self._generate_script_prepare_sequences()
        script_content += self._generate_msa_copy_section()
        script_content += self._generate_script_run_alphafold()
        script_content += self._generate_script_extract_best_rank()
        script_content += self._generate_script_extract_confidence()
        script_content += self._generate_script_create_msas_table()
        script_content += self._generate_script_unsanitize_ids()
        script_content += self._generate_script_update_structures_map()
        script_content += self._generate_missing_table_propagation()
        script_content += self.generate_completion_check_footer()

        return script_content

    def _generate_script_prepare_sequences(self) -> str:
        """Generate the sequence preparation part of the script.

        Bundle/Each inputs go through the combinatorics queries-builder, which
        groups sequences into one colon-joined ColabFold query per complex. A
        bare input keeps the legacy copy/convert path (one monomer per row).
        """
        if self._uses_combinatorics:
            generate_combinatorics_config(
                self.combinatorics_config_file,
                proteins=(self.proteins, "sequences", "protein"),
            )
            return f"""echo "Building ColabFold queries from combinatorics config"
python {self.alphafold_queries_py} \\
    --combinatorics-config "{self.combinatorics_config_file}" \\
    --output-csv "{self.queries_csv}"

"""

        # Get source file from input DataStream
        source_file = self.sequences_stream.map_table

        # Check if we need to convert from FASTA folder to CSV
        # Only treat as ProteinMPNN seqs folder if it ends with /seqs or is a directory
        if source_file.endswith("seqs") and not source_file.endswith(".csv"):
            # ProteinMPNN seqs folder - convert .fa files to CSV. Pass the stream
            # JSON so an upstream id filter gates which .fa files are converted
            # (fa_to_csv builds the allowed-basename set from ids_expanded).
            self.sequences_stream.save_json(self.sequences_json)
            return f"""echo "Converting ProteinMPNN .fa files to queries CSV"
python {self.fa_to_csv_fasta_py} {source_file} {self.queries_csv} {self.queries_fasta} --ds-json "{self.sequences_json}"

"""
        else:
            # Direct CSV/FASTA stream - materialize a filtered queries CSV (honors
            # an upstream id filter), then re-activate the AlphaFold env.
            self.sequences_stream.save_json(self.sequences_json)
            block = self.generate_filtered_map_table_block(
                self.sequences_json, self.queries_csv, required_columns=["id", "sequence"]
            )
            block += self.activate_environment()
            return block

    def _generate_msa_copy_section(self) -> str:
        """Generate script section to copy pre-computed MSA files to Folding directory."""
        if self.msas_input is None:
            return ""

        msa_table_path = self.msas_input.tables.msas.info.path

        return f"""echo "Copying pre-computed MSA files to Folding (execution) directory"
python "{self.msa_copy_py}" \\
    --msa-table "{msa_table_path}" \\
    --output-folder "{self.folding_folder}"

"""

    def _generate_script_run_alphafold(self) -> str:
        """Generate the AlphaFold execution part of the script."""
        from .config_manager import ConfigManager
        cm = ConfigManager()
        scheduler = cm.get_scheduler()
        env_manager = cm.get_env_manager()

        # Build AlphaFold options
        af_options = ""
        if self.num_relax > 0:
            af_options += f" --amber --use-gpu-relax --num-relax {self.num_relax}"
        if self.num_recycle != 3:
            af_options += f" --num-recycle {self.num_recycle}"
        if self.rand_seed != 0:
            af_options += f" --random-seed {self.rand_seed}"

        # Determine colabfold_batch command
        if scheduler == "colab":
            # Colab: colabfold installed into system python; run via /usr/local/bin
            colabfold_cmd = "/usr/local/bin/colabfold_batch"
        elif env_manager == "pip":
            # Pure pip mode (no conda/mamba at all): rely on PATH
            colabfold_cmd = "colabfold_batch"
        else:
            # Cluster: absolute path to LocalColabFold binary
            colabfold_cmd = str(self.colabfold_batch)

        if scheduler == "colab":
            # Run in a clean subshell without conda env vars to avoid interference
            # with Colab's JAX/GPU/tensorflow stack
            run_colabfold = f"""(unset CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL; PATH="/usr/local/bin:/usr/bin:/bin:$PATH" {colabfold_cmd} {self.queries_csv} "{self.folding_folder}" {af_options})"""
        else:
            run_colabfold = f"""{self.container_prefix()}{colabfold_cmd} {self.queries_csv} "{self.folding_folder}" {af_options}"""

        # On Colab, colabfold_batch runs as the system Python and downloads its
        # AF2 params into $HOME/.cache/colabfold/params (i.e. /root/.cache/...),
        # which is ephemeral and not shared. Point that path at the shared
        # AlphaFoldParams cache (on Drive once the user repoints `cache`), so the
        # first run downloads into the shared cache and later sessions reuse it —
        # the same wiring AF2BIND/BioEmu use. On the cluster, LocalColabFold owns
        # its own params dir; leave it alone (don't regress that path).
        cache_wiring = ""
        if scheduler == "colab":
            af2_params_root = self.folders.get("AlphaFoldParams", "")
            if af2_params_root:
                params_dir = f"{af2_params_root}/params"
                cache_wiring = (
                    f'mkdir -p "{params_dir}"\n'
                    + link_weights_block(
                        target_dir=params_dir,
                        link_path="$HOME/.cache/colabfold/params",
                        label="AF2 params",
                    )
                    + "\n\n"
                )

        return f"""echo "Running AlphaFold2/ColabFold"
echo "Options: {af_options}"
echo "Output folder: {self.output_folder}"
echo "Folding (raw) folder: {self.folding_folder}"

# Share AF2 params via the AlphaFoldParams cache instead of ephemeral ~/.cache.
{cache_wiring}# execution/ auto-created by the pipeline.
# Run ColabFold batch
{run_colabfold}

"""

    def _generate_script_extract_best_rank(self) -> str:
        """Generate script to extract best rank structures.

        ColabFold dumps everything into execution/ (folding_folder). We
        relocate the best-rank PDB for each sequence into structures/ and
        copy any .a3m MSAs into msas/. Both destinations are pre-created
        by the pipeline layout step.
        """
        structures_dir = self.stream_folder("structures")
        msa_section = f"""
# Copy MSA files to msas/ stream folder
echo "Extracting MSA files"
for msa_file in *.a3m; do
    if [ -f "$msa_file" ]; then
        cp "$msa_file" "{self.msas_folder}/"
        echo "Copied MSA: $msa_file"
    fi
done
"""

        return f"""echo "Extracting best structures from Folding (execution) subfolder"
# AlphaFold creates files like:
# - sequenceid_unrelaxed_rank_001_alphafold2_ptm_model_N_seed_SSS.pdb
# - sequenceid_relaxed_rank_001_alphafold2_ptm_model_N_seed_SSS.pdb
# - sequenceid.a3m (MSA file, not generated in single_sequence mode)
# We route best-rank PDBs into structures/ and MSAs into msas/.

cd "{self.folding_folder}"

# Handle relaxed format (preferred if both exist)
for file in *_relaxed_rank_001_*.pdb; do
    if [ -f "$file" ]; then
        base=$(echo "$file" | sed 's/_relaxed_rank_001_.*/.pdb/')
        cp "$file" "{structures_dir}/$base"
        echo "Extracted: $file -> $base"
    fi
done

# Handle unrelaxed format (only if relaxed doesn't exist)
for file in *_unrelaxed_rank_001_*.pdb; do
    if [ -f "$file" ]; then
        base=$(echo "$file" | sed 's/_unrelaxed_rank_001_.*/.pdb/')
        if [ ! -f "{structures_dir}/$base" ]; then
            cp "$file" "{structures_dir}/$base"
            echo "Extracted: $file -> $base"
        else
            echo "Skipped $file (relaxed version already exists as $base)"
        fi
    fi
done
{msa_section}
cd - > /dev/null

"""

    def _generate_script_update_structures_map(self) -> str:
        """Generate script to update structures_map.csv with actual runtime output files."""
        structures_dir = self.stream_folder("structures")
        return f"""echo "Updating structures map with actual output files"
python {self.update_map_py} --structures-map "{self.structures_map}" --output-folder "{structures_dir}"

"""

    def _generate_script_extract_confidence(self) -> str:
        """Generate script section to extract confidence metrics from JSON files."""
        structures_dir = self.stream_folder("structures")
        return f"""echo "Extracting confidence metrics from JSON files"
python {self.alphafold_confidence_py} "{self.folding_folder}" "{structures_dir}" "{self.confidence_csv}"

"""

    def _generate_script_create_msas_table(self) -> str:
        """Generate script section to create MSAs CSV table."""
        return f"""echo "Creating MSAs table"
python {self.alphafold_msas_py} "{self.msas_folder}" "{self.queries_csv}" "{self.msa_csv}"

"""

    def _generate_script_unsanitize_ids(self) -> str:
        """Generate script section to restore original IDs after ColabFold sanitization.

        The expected ids come from the runtime queries CSV (written by
        _generate_script_prepare_sequences), not from a config-time map_table —
        the structures_map.csv is written later, from the actual PDBs.
        """
        structures_dir = self.stream_folder("structures")
        return f"""echo "Restoring original IDs (unsanitizing ColabFold output)"
python {self.unsanitize_ids_py} \\
    --predicted-ids "{self.queries_csv}" \\
    --rename-files "{structures_dir}" pdb \\
    --rename-files "{self.msas_folder}" a3m \\
    --fix-csv "{self.confidence_csv}" \\
    --fix-csv "{self.msa_csv}"

"""

    def _generate_missing_table_propagation(self) -> str:
        """Propagate the union of upstream `missing` manifests into tables/missing.csv."""
        return self.generate_missing_propagation(
            self.proteins, self.sequences_stream, missing_csv=self.missing_csv
        )

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after AlphaFold execution."""
        # Predict output IDs. With Bundle/Each, ids come from combinatorics
        # (e.g. a complex of p1+p2 yields id "p1+p2"); otherwise one id per
        # input sequence (monomer-per-row, the legacy path).
        if self._uses_combinatorics:
            sequence_ids, _provenance = predict_output_ids_with_provenance(
                proteins=(self.proteins, "sequences", "protein"),
            )
        else:
            sequence_ids = list(self.sequences_stream.ids)

        # Structure files land in structures/; map_table co-locates there.
        # The per-design map_table is written at runtime by
        # _generate_script_update_structures_map() from the actual PDBs; here
        # we only declare the stream and its map_table path.
        structure_files = [self.stream_path("structures", "<id>.pdb")]

        structures = DataStream(
            name="structures",
            ids=sequence_ids,
            files=structure_files,
            map_table=self.structures_map,
            format="pdb"
        )

        # MSA files (AlphaFold generates .a3m files)
        msa_files = [os.path.join(self.msas_folder, "<id>.a3m")]

        msas = DataStream(
            name="msas",
            ids=sequence_ids,
            files=msa_files,
            map_table=self.msa_csv,
            format="a3m"
        )

        # Tables
        tables = {
            "structures": TableInfo(
                name="structures",
                path=self.structures_map,
                columns=["id", "file"],
                description="AlphaFold predicted structures"
            ),
            "confidence": TableInfo(
                name="confidence",
                path=self.confidence_csv,
                columns=["id", "structure", "plddt", "max_pae", "ptm"],
                description="AlphaFold confidence metrics extracted from best rank models"
            ),
            "msas": TableInfo(
                name="msas",
                path=self.msa_csv,
                columns=["id", "sequences.id", "sequence", "msa_file"],
                description="MSA files for sequence recycling between predictions"
            )
        }

        # Declare a `missing` table when any input axis carries one.
        if self._collect_upstream_missing_paths(self.proteins, self.sequences_stream):
            tables["missing"] = self.missing_table_info(self.missing_csv)

        return {
            "structures": structures,
            "msas": msas,
            "tables": tables,
            "output_folder": self.output_folder,
            "rendering_parameters": {
                "structures": {
                    "color_by": "plddt",
                    "plddt_upper": 100,
                }
            }
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including AlphaFold-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "af_params": {
                "num_relax": self.num_relax,
                "num_recycle": self.num_recycle,
                "rand_seed": self.rand_seed,
                "msas": str(self.msas_input) if self.msas_input else None
            }
        })
        return base_dict
