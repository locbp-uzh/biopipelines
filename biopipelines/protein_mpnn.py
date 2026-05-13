# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
ProteinMPNN configuration for sequence design from protein structures.
"""

import os
from typing import Dict, List, Any, Union, Tuple

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .combinatorics import generate_multiplied_ids, generate_multiplied_ids_pattern
    from .biopipelines_io import Resolve
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from combinatorics import generate_multiplied_ids, generate_multiplied_ids_pattern
    from biopipelines_io import Resolve


class ProteinMPNN(BaseConfig):
    """
    Configuration for ProteinMPNN sequence design.
    """

    TOOL_NAME = "ProteinMPNN"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False,
                        cudatoolkit="11.3", **kwargs):
        try:
            from .config_manager import ConfigManager
        except ImportError:
            from config_manager import ConfigManager

        repo_dir = folders.get("ProteinMPNN", "")
        parent_dir = os.path.dirname(repo_dir)
        biopipelines = folders.get("biopipelines", "")

        pmpnn_env = ConfigManager().get_environment("ProteinMPNN") or "mlfold"
        rfd_env   = ConfigManager().get_environment("RFdiffusion")

        # Skip when the repo exists AND `import torch` works inside the
        # configured env. The env must be present too — when ProteinMPNN
        # is the first tool to bring up a shared SE3nv, neither side has
        # created it yet.
        # The torch import probe is the load-bearing check (env + package).
        skip = "" if force_reinstall else f"""# Check if already installed
if [ -d "{repo_dir}" ] \\
   && {env_manager} run -n {pmpnn_env} python -c "import torch" >/dev/null 2>&1; then
    echo "ProteinMPNN already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        clone_block = f"""mkdir -p {parent_dir}
cd {parent_dir}
if [ ! -d "{repo_dir}" ]; then
    git clone https://github.com/dauparas/ProteinMPNN.git
fi"""

        if rfd_env and pmpnn_env == rfd_env:
            # Same conda/mamba environment as RFdiffusion: ProteinMPNN
            # only needs PyTorch in there. If the env doesn't exist yet
            # (e.g. user runs ProteinMPNN.install() before RFdiffusion's),
            # create it from the BioPipelines spec; if it already exists
            # (RFdiffusion got there first), reuse it as-is. Either way the
            # repo gets cloned and the verification below confirms torch.
            env_block = cls._env_install_block(pmpnn_env, env_manager, biopipelines)
            shared_env_check = cls._env_exists_check(pmpnn_env, env_manager)
            return f"""echo "=== Installing ProteinMPNN ==="
{skip}{clone_block}

# Bring up the shared '{pmpnn_env}' env if RFdiffusion hasn't already done so.
if {shared_env_check}; then
    echo "Shared env '{pmpnn_env}' already exists — reusing it."
else
    echo "Shared env '{pmpnn_env}' missing — creating it from the BioPipelines spec."
    {env_block}
    if [ $? -ne 0 ]; then
        echo "ERROR: could not create shared env '{pmpnn_env}'."
        exit 1
    fi
fi

# Verify installation (repo cloned + torch importable in the shared env)
if [ -f "{repo_dir}/protein_mpnn_run.py" ] \\
   && {env_manager} run -n {pmpnn_env} python -c "import torch" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== ProteinMPNN installation complete ==="
else
    echo "ERROR: ProteinMPNN verification failed (script missing or torch not importable)"
    exit 1
fi
"""

        # Dedicated environment — create it and install PyTorch
        remove_block = cls._env_remove_block(pmpnn_env, env_manager) if force_reinstall else ""
        return f"""echo "=== Installing ProteinMPNN ==="
{skip}{clone_block}

{remove_block}
# Create dedicated ProteinMPNN environment
{env_manager} create --name {pmpnn_env} -y || echo "WARNING: environment '{pmpnn_env}' may already exist, continuing"
{env_manager} run -n {pmpnn_env} conda install pytorch torchvision torchaudio cudatoolkit={cudatoolkit} -c pytorch -y || echo "WARNING: PyTorch install failed, skipping"

# Verify installation
if [ -f "{repo_dir}/protein_mpnn_run.py" ] && {env_manager} run -n {pmpnn_env} python -c "import torch" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== ProteinMPNN installation complete ==="
else
    echo "ERROR: ProteinMPNN verification failed (script missing or torch not importable)"
    exit 1
fi
"""

    # Lazy path descriptors — routed through the canonical sub-layout.
    #   configuration/  — parsed_pdbs.jsonl, fixed_pos.jsonl, fixed_designed.csv,
    #                     .input_structures.json (all config-time inputs).
    #   execution/      — raw ProteinMPNN output; includes seqs/ (.fa files).
    #   sequences/      — content-bearing stream: sequences.csv IS the content
    #                     + map_table; sequences.fasta is its FASTA twin.
    #   tables/         — standalone TableInfo CSVs (missing, proteinmpnn_results).
    parsed_pdbs_jsonl = Path(lambda self: self.configuration_path("parsed_pdbs.jsonl"))
    fixed_jsonl = Path(lambda self: self.configuration_path("fixed_pos.jsonl"))
    sele_csv = Path(lambda self: self.configuration_path("fixed_designed.csv"))
    # ProteinMPNN's CLI creates seqs/<id>.fa inside --out_folder.
    # Point --out_folder at execution/ so the raw .fa dumps live there.
    pmpnn_out_folder = Path(lambda self: self.execution_folder)
    seqs_folder = Path(lambda self: self.execution_path("seqs"))
    main_table = Path(lambda self: self.table_path("proteinmpnn_results"))
    queries_csv = Path(lambda self: self.stream_path("sequences", "sequences.csv"))
    queries_fasta = Path(lambda self: self.stream_path("sequences", "sequences.fasta"))
    structures_json = Path(lambda self: self.configuration_path(".input_structures.json"))

    missing_csv = Path(lambda self: self.table_path("missing"))

    # Helper scripts
    fixed_py = Path(lambda self: self.pipe_script_path("pipe_pmpnn_fixed_positions.py"))
    table_py = Path(lambda self: self.pipe_script_path("pipe_pmpnn_table.py"))
    fa_to_csv_fasta_py = Path(lambda self: self.pipe_script_path("pipe_fa_to_csv_fasta.py"))

    # ProteinMPNN installation scripts
    parse_py = Path(lambda self: os.path.join(self.folders["ProteinMPNN"], "helper_scripts", "parse_multiple_chains.py"))
    pmpnn_py = Path(lambda self: os.path.join(self.folders["ProteinMPNN"], "protein_mpnn_run.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 num_sequences: int = 1,
                 fixed: Union[str, Tuple['TableInfo', str]] = "",
                 redesigned: Union[str, Tuple['TableInfo', str]] = "",
                 chain: str = "auto",
                 sampling_temp: float = 0.1,
                 model_name: str = "v_48_020",
                 soluble_model: bool = True,
                 remove_duplicates: bool = True,
                 fill_gaps: str = "G",
                 bias_AA_jsonl: str = "",
                 omit_AA_jsonl: str = "",
                 seed: int = 0,
                 ca_noise_std: float = 0.0,
                 **kwargs):
        """
        Initialize ProteinMPNN configuration.

        Args:
            structures: Input structures as DataStream or StandardizedOutput
            num_sequences: Number of sequences to generate per structure
            fixed: Fixed positions. Accepts:
                   - PyMOL-style selection string: "10-20+30-40"
                   - TableReference: table.column_name
            redesigned: Designed positions. Accepts:
                   - PyMOL-style selection string: "10-20+30-40"
                   - TableReference: table.column_name
            chain: Chain to apply fixed positions to ("auto" detects from input structure)
            sampling_temp: Sampling temperature for sequence generation
            model_name: ProteinMPNN model variant
            soluble_model: Use soluble protein model
            remove_duplicates: Remove duplicate sequences from output (default True)
            fill_gaps: Amino acid to replace X (unknown/gap residues) with (default "G" for glycine).
                       Empty string means no filling (X is kept as-is).
            bias_AA_jsonl: Path to a JSONL file biasing per-position amino-acid logits
                       (upstream --bias_AA_jsonl). Empty string disables.
            ca_noise_std: Backbone Cα Gaussian noise standard deviation injected during
                       inference (upstream --backbone_noise). 0.0 disables (default).
            omit_AA_jsonl: Path to a JSONL file omitting specific amino acids per position
                       (upstream --omit_AA_jsonl). Empty string disables.
            seed: Random seed for reproducible sampling (upstream --seed). 0 (default) lets
                       ProteinMPNN choose its own seed.

        Output:
            Streams: sequences (.csv), fasta (.fasta)
            Tables:
                sequences: id | structures.id | source_pdb | sequence | score | seq_recovery | gaps
                missing: id | removed_by | cause
        """
        # Resolve input to DataStream
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        # Store parameters
        self.num_sequences = num_sequences
        self.fixed = fixed
        self.redesigned = redesigned
        self.chain = chain
        self.sampling_temp = sampling_temp
        self.model_name = model_name
        self.soluble_model = soluble_model
        self.remove_duplicates = remove_duplicates
        self.fill_gaps = fill_gaps
        self.bias_AA_jsonl = bias_AA_jsonl
        self.omit_AA_jsonl = omit_AA_jsonl
        self.seed = seed
        self.ca_noise_std = ca_noise_std

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate ProteinMPNN-specific parameters."""
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")

        if self.num_sequences <= 0:
            raise ValueError("num_sequences must be positive")

        if self.sampling_temp <= 0:
            raise ValueError("sampling_temp must be positive")

        if self.ca_noise_std < 0:
            raise ValueError("ca_noise_std must be non-negative")

        if self.seed < 0:
            raise ValueError("seed must be non-negative")

        valid_models = ["v_48_002", "v_48_010", "v_48_020", "v_48_030"]
        if self.model_name not in valid_models:
            raise ValueError(f"model_name must be one of: {valid_models}")

        _validate_freeform_string("chain", self.chain)
        _validate_freeform_string("fill_gaps", self.fill_gaps)
        if isinstance(self.fixed, str):
            _validate_freeform_string("fixed", self.fixed)
        if isinstance(self.redesigned, str):
            _validate_freeform_string("redesigned", self.redesigned)

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input structures."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get ProteinMPNN configuration display lines."""
        config_lines = super().get_config_display()
        config_lines.extend([
            f"NUM SEQUENCES PER TARGET: {self.num_sequences}",
            f"FIXED: {self.fixed or 'None'}",
            f"REDESIGNED: {self.redesigned or 'None'}",
            f"CHAIN: {self.chain}",
            f"SAMPLING T: {self.sampling_temp}",
            f"MODEL: {self.model_name}",
            f"SOLUBLE: {self.soluble_model}"
        ])
        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate ProteinMPNN execution script."""
        script_content = "#!/bin/bash\n"
        script_content += "# ProteinMPNN execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self._generate_script_prepare_inputs()
        script_content += self._generate_script_run_proteinmpnn()
        script_content += self._generate_script_create_table()
        script_content += self.generate_completion_check_footer()
        return script_content

    def _generate_script_prepare_inputs(self) -> str:
        """Generate the input preparation part of the script."""

        # Serialize DataStream to JSON file (proper way to pass ids + files to pipe_script)
        self.structures_stream.save_json(self.structures_json)

        fixed_param = self.fixed if self.fixed else "-"
        designed_param = self.redesigned if self.redesigned else "-"

        # parse_multiple_chains.py wants a folder, not a file list — pick the
        # first expanded ID's file and take its dirname. The index=0 form is
        # the lazy-safe variant documented in developer_manual.md.
        return f"""FIRST_ID={Resolve.stream_ids(self.structures_json, index=0)}
FIRST_FILE={Resolve.stream_item(self.structures_json, "$FIRST_ID")}
INPUT_DIR=$(dirname "$FIRST_FILE")

echo "Determining fixed positions"
python {self.fixed_py} "{self.structures_json}" "{fixed_param}" "{designed_param}" "{self.chain}" "{self.fixed_jsonl}" "{self.sele_csv}"

echo "Parsing multiple PDBs"
python {self.parse_py} --input_path $INPUT_DIR --output_path {self.parsed_pdbs_jsonl}

"""

    def _generate_script_run_proteinmpnn(self) -> str:
        """Generate the ProteinMPNN execution part of the script.

        Iterates per-PDB so a structure that triggers a Python exception inside
        protein_mpnn_run.py (e.g. an unexpected chain in tied_featurize) does
        not abort the rest of the batch. Each iteration filters the merged
        parsed_pdbs.jsonl and fixed_pos.jsonl down to a single entry, runs the
        model, and on failure records the id + error in missing.csv.
        """
        pmpnn_options = f"--num_seq_per_target {self.num_sequences}"
        pmpnn_options += f" --sampling_temp {self.sampling_temp}"
        pmpnn_options += f" --model_name {self.model_name}"

        if self.soluble_model:
            pmpnn_options += " --use_soluble_model"

        if self.bias_AA_jsonl:
            pmpnn_options += f" --bias_AA_jsonl {self.bias_AA_jsonl}"

        if self.omit_AA_jsonl:
            pmpnn_options += f" --omit_AA_jsonl {self.omit_AA_jsonl}"

        if self.seed > 0:
            pmpnn_options += f" --seed {self.seed}"

        if self.ca_noise_std > 0:
            pmpnn_options += f" --backbone_noise {self.ca_noise_std}"

        per_pdb_dir = os.path.join(self.execution_folder, "per_pdb")

        return f"""echo "Running model (per-PDB loop)"
echo "Options: {pmpnn_options}"
echo "Output folder: {self.output_folder}"
echo "ProteinMPNN --out_folder (execution/): {self.pmpnn_out_folder}"

mkdir -p "{per_pdb_dir}"
# Seed missing.csv header now so partial failures still produce a readable file.
if [ ! -f "{self.missing_csv}" ]; then
    echo "id,removed_by,cause" > "{self.missing_csv}"
fi

for struct_id in {Resolve.stream_ids(self.structures_json)}; do
    SUBSET_JSONL="{per_pdb_dir}/${{struct_id}}.jsonl"
    SUBSET_FIXED="{per_pdb_dir}/${{struct_id}}_fixed.jsonl"
    SUBSET_LOG="{per_pdb_dir}/${{struct_id}}.log"

    python -c "
import json, sys
sid = sys.argv[1]
with open(sys.argv[2]) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        rec = json.loads(line)
        if rec.get('name') == sid:
            with open(sys.argv[3], 'w') as out:
                out.write(json.dumps(rec) + '\\n')
            break
    else:
        sys.exit(2)
with open(sys.argv[4]) as f:
    full_fixed = json.loads(f.read())
sub = {{sid: full_fixed[sid]}} if sid in full_fixed else {{}}
with open(sys.argv[5], 'w') as out:
    out.write(json.dumps(sub))
" "$struct_id" "{self.parsed_pdbs_jsonl}" "$SUBSET_JSONL" "{self.fixed_jsonl}" "$SUBSET_FIXED"

    if [ ! -s "$SUBSET_JSONL" ]; then
        echo "  [skip] $struct_id: not found in parsed_pdbs.jsonl"
        echo "${{struct_id}},ProteinMPNN,not_in_parsed_pdbs" >> "{self.missing_csv}"
        continue
    fi

    echo "  [run]  $struct_id"
    set +e
    {self.container_prefix()}python {self.pmpnn_py} --jsonl_path "$SUBSET_JSONL" --fixed_positions_jsonl "$SUBSET_FIXED" --out_folder {self.pmpnn_out_folder} {pmpnn_options} > "$SUBSET_LOG" 2>&1
    rc=$?
    set -e
    if [ $rc -ne 0 ]; then
        # Capture last non-empty line of the log as the cause; sanitize commas.
        CAUSE=$(grep -vE '^\\s*$' "$SUBSET_LOG" | tail -n1 | tr ',\\n' ' ')
        echo "  [fail] $struct_id (exit $rc): $CAUSE"
        echo "${{struct_id}},ProteinMPNN,\\"$CAUSE\\"" >> "{self.missing_csv}"
        cat "$SUBSET_LOG"
    fi
done

"""

    def _generate_script_create_table(self) -> str:
        """Generate the table creation part of the script."""
        duplicates_flag = " --duplicates" if not self.remove_duplicates else ""
        fill_gaps_flag = f' --fill-gaps "{self.fill_gaps}"' if self.fill_gaps else ""
        step_tool_name = os.path.basename(self.output_folder)

        # Check for upstream missing table
        upstream_missing_path = self._get_upstream_missing_table_path(
            self.structures_stream
        )
        upstream_missing_flag = f' --upstream-missing "{upstream_missing_path}"' if upstream_missing_path else ""

        return f"""echo "Creating results table and queries files"
python {self.table_py} {self.seqs_folder} {self.pipeline_name} "-" {self.main_table}

echo "Creating queries CSV and FASTA from results table"
python {self.fa_to_csv_fasta_py} {self.seqs_folder} {self.queries_csv} {self.queries_fasta} --ds-json "{self.structures_json}"{duplicates_flag}{fill_gaps_flag} --missing-csv "{self.missing_csv}" --step-tool-name "{step_tool_name}"{upstream_missing_flag}

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after ProteinMPNN execution."""
        # Raw ProteinMPNN .fa files land in execution/seqs/ (one per input
        # structure). These aren't a first-class pipeline stream — they're
        # an intermediate format the post-processing step reads.
        fasta_ids = self.structures_stream.ids
        fasta_files = [os.path.join(self.seqs_folder, "<id>.fa")]

        # Predict sequence IDs (stream_id + sequence number)
        suffix_pattern = f"<1..{self.num_sequences}>"
        sequence_ids = generate_multiplied_ids_pattern(
            self.structures_stream.ids, suffix_pattern,
            input_stream_name="structures"
        )

        # Content-bearing sequences stream: queries_csv lives under
        # sequences/ and doubles as the map_table.
        sequences = DataStream(
            name="sequences",
            ids=sequence_ids,
            files=[],
            map_table=self.queries_csv,
            format="csv"
        )

        # Fasta stream — points at the raw execution/seqs/<id>.fa dumps.
        fasta = DataStream(
            name="fasta",
            ids=fasta_ids,
            files=fasta_files,
            format="fasta",
            metadata={
                "sequences_per_file": self.num_sequences,
                "contains_original": False,
            }
        )

        tables = {
            "sequences": TableInfo(
                name="sequences",
                path=self.queries_csv,
                columns=["id", "structures.id", "source_pdb", "sequence", "score", "seq_recovery", "gaps"],
                description="ProteinMPNN sequence results"
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "cause"],
                description="IDs removed (duplicates or upstream) with removal reason"
            )
        }

        return {
            "sequences": sequences,
            "fasta": fasta,
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "mpnn_params": {
                "num_sequences": self.num_sequences,
                "fixed": self.fixed,
                "redesigned": self.redesigned,
                "chain": self.chain,
                "sampling_temp": self.sampling_temp,
                "model_name": self.model_name,
                "soluble_model": self.soluble_model,
                "remove_duplicates": self.remove_duplicates,
                "bias_AA_jsonl": self.bias_AA_jsonl,
                "omit_AA_jsonl": self.omit_AA_jsonl,
                "seed": self.seed,
                "ca_noise_std": self.ca_noise_std,
            }
        })
        return base_dict
