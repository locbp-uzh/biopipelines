# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
LigandMPNN configuration for ligand-aware sequence design.
"""

import os
import json
from typing import Dict, List, Any, Union, Tuple, Optional

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .combinatorics import generate_multiplied_ids, generate_multiplied_ids_pattern
    from .biopipelines_io import Resolve
    from .input_standardization import resolve_basic_input
    from .ligand import Ligand
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from combinatorics import generate_multiplied_ids, generate_multiplied_ids_pattern
    from biopipelines_io import Resolve
    from input_standardization import resolve_basic_input
    from ligand import Ligand


class LigandMPNN(BaseConfig):
    """
    LigandMPNN configuration for ligand-aware sequence design.
    """

    TOOL_NAME = "LigandMPNN"
    TOOL_VERSION = "1.1"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        repo_dir = folders.get("LigandMPNN", "")
        parent_dir = os.path.dirname(repo_dir)
        biopipelines = folders.get("biopipelines", "")
        env_check = cls._env_exists_check("ligandmpnn_env", env_manager)
        skip = "" if force_reinstall else f"""# Check if already installed
if [ -d "{repo_dir}" ] && [ -d "{repo_dir}/model_params" ] && {env_check}; then
    echo "LigandMPNN already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block("ligandmpnn_env", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("ligandmpnn_env", env_manager, biopipelines)
        return f"""echo "=== Installing LigandMPNN ==="
{skip}mkdir -p "{parent_dir}"
cd "{parent_dir}"
if [ ! -d "{repo_dir}" ]; then
    git clone https://github.com/dauparas/LigandMPNN.git
fi
cd "{repo_dir}"
bash get_model_params.sh "./model_params"

{remove_block}
{env_block}
# LigandMPNN's own pinned requirements (lives in the cloned repo)
{env_manager} run -n ligandmpnn_env pip install -r "{repo_dir}/requirements.txt"

# Verify installation
if [ -d "{repo_dir}/model_params" ] && {env_manager} run -n ligandmpnn_env python -c "import torch" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== LigandMPNN installation complete ==="
else
    echo "ERROR: LigandMPNN verification failed (model_params missing or torch not importable)"
    exit 1
fi
"""

    # Lazy path descriptors — same shape as ProteinMPNN on the new layout.
    #   configuration/  — structures JSON, positions JSON.
    #   execution/      — LigandMPNN CLI writes seqs/ here.
    #   sequences/      — content-bearing stream (sequences.csv + .fasta).
    #   tables/         — missing.
    lmpnn_out_folder = Path(lambda self: self.execution_folder)
    seqs_folder = Path(lambda self: self.execution_path("seqs"))
    queries_csv = Path(lambda self: self.stream_path("sequences", "sequences.csv"))
    queries_fasta = Path(lambda self: self.stream_path("sequences", "sequences.fasta"))
    structures_json = Path(lambda self: self.configuration_path(".input_structures.json"))
    positions_json = Path(lambda self: self.configuration_path("lmpnn_positions.json"))
    positions_args_json = Path(lambda self: self.configuration_path("lmpnn_positions_args.json"))

    missing_csv = Path(lambda self: self.table_path("missing"))

    # Helper script paths
    fa_to_csv_fasta_py = Path(lambda self: self.pipe_script_path("pipe_fa_to_csv_fasta.py"))
    lmpnn_folder = Path(lambda self: os.path.join(self.folders["data"], "LigandMPNN"))
    runtime_positions_py = Path(lambda self: self.pipe_script_path("pipe_lmpnn_runtime_positions.py"))
    resolve_positions_py = Path(lambda self: self.pipe_script_path("resolve_lmpnn_positions.py"))
    ligand_json = Path(lambda self: self.configuration_path("input_ligand.json"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 ligand: Optional[Union[str, DataStream, StandardizedOutput]] = None,
                 num_sequences: int = 1,
                 fixed: Union[str, Tuple['TableInfo', str]] = "",
                 redesigned: Union[str, Tuple['TableInfo', str]] = "",
                 design_within: float = 5.0,
                 chain: str = "A",
                 model: str = "v_32_010",
                 num_batches: int = 1,
                 remove_duplicates: bool = True,
                 fill_gaps: str = "G",
                 temperature: float = 0.0,
                 bias_AA_per_residue: str = "",
                 seed: int = 0,
                 **kwargs):
        """
        Initialize LigandMPNN configuration.

        Args:
            structures: Input structures as DataStream or StandardizedOutput
            ligand: Compounds stream (Ligand(code="LIG") or any
                    compounds-producing tool) naming the bound ligand to focus
                    the design around. A bare string ("LIG") is shorthand for
                    Ligand(code="LIG") and creates an internal code-only Ligand.
                    The residue `code` is read from the stream's `code` column
                    at runtime.
            num_sequences: Number of sequences per batch (maps to --batch_size in LigandMPNN)
            fixed: Fixed positions. Accepts:
                   - PyMOL selection string: "10-20+30-40"
                   - TableReference: table.column_name
            redesigned: Designed positions. Accepts:
                   - PyMOL selection string: "10-20+30-40"
                   - TableReference: table.column_name
            design_within: Distance in Angstrom from ligand to redesign (fallback if positions not specified)
            chain: Default chain ID for chainless position input (default "A")
            model: LigandMPNN model version to use
            num_batches: Number of batches to run
            remove_duplicates: Remove duplicate sequences from output (default True)
            fill_gaps: Amino acid to replace X (unknown/gap residues) with (default "G" for glycine).
                       Empty string means no filling (X is kept as-is).
            temperature: Sampling temperature (upstream --temperature). 0.0 (default) leaves
                       LigandMPNN at its built-in default.
            bias_AA_per_residue: Path to a JSON file biasing per-position amino-acid logits
                       (upstream --bias_AA_per_residue). Empty string disables. Format:
                       {'A12': {'G': -0.3, 'C': -2.0}, ...}.
            seed: Random seed for reproducible sampling (upstream --seed). 0 (default) lets
                       LigandMPNN choose its own seed.

        Output:
            Streams: sequences (.csv), fasta (.fasta)
            Tables:
                sequences: id | sequence | sample | T | seed | overall_confidence | ligand_confidence | seq_rec | gaps
                missing: id | removed_by | kind | cause
        """
        # Resolve input to DataStream
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        # Ligand binding-site focus — optional compounds stream; the residue
        # code is resolved from its `code` column at runtime. A bare string is
        # promoted to an internal code-only Ligand.
        self.ligand_stream: Optional[DataStream] = resolve_basic_input(
            ligand, Ligand, "compounds", "code")
        self.num_sequences = num_sequences
        self.fixed = fixed
        self.redesigned = redesigned
        self.design_within = design_within
        self.chain = chain
        self.model = model
        self.num_batches = num_batches
        self.remove_duplicates = remove_duplicates
        self.fill_gaps = fill_gaps
        self.temperature = temperature
        self.bias_AA_per_residue = bias_AA_per_residue
        self.seed = seed

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate LigandMPNN-specific parameters."""
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")

        if not self.ligand_stream or len(self.ligand_stream) == 0:
            raise ValueError("ligand (a compounds stream, e.g. Ligand(code=...)) is required and must not be empty")

        if self.num_sequences <= 0:
            raise ValueError("num_sequences must be positive")

        if self.num_batches <= 0:
            raise ValueError("num_batches must be positive")

        if self.design_within <= 0:
            raise ValueError("design_within must be positive")

        if self.temperature < 0:
            raise ValueError("temperature must be non-negative")

        if self.seed < 0:
            raise ValueError("seed must be non-negative")

        # Validate model name
        valid_models = ["v_32_005", "v_32_010", "v_32_020", "v_32_025"]
        if self.model not in valid_models:
            raise ValueError(f"model must be one of: {valid_models}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input structures."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get LigandMPNN configuration display lines."""
        config_lines = super().get_config_display()
        config_lines.extend([
            "LIGAND: (code resolved from compounds stream at runtime)",
            f"FIXED: {self.fixed or 'Auto (from table or ligand-based)'}",
            f"REDESIGNED: {self.redesigned or 'Auto (from table or ligand-based)'}",
            f"CHAIN: {self.chain}",
            f"DESIGN WITHIN: {self.design_within}A",
            f"NUM SEQUENCES: {self.num_sequences}",
            f"NUM BATCHES: {self.num_batches}",
            f"MODEL: {self.model}"
        ])
        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate LigandMPNN execution script."""

        # Serialize DataStream to JSON file (proper way to pass ids + files to pipe_script)
        self.structures_stream.save_json(self.structures_json)

        script_content = "#!/bin/bash\n"
        script_content += "# LigandMPNN execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self._generate_script_setup_positions()
        script_content += self._generate_script_run_ligandmpnn()
        script_content += self._generate_script_convert_outputs()
        script_content += self.generate_completion_check_footer()
        return script_content

    def _generate_script_setup_positions(self) -> str:
        """Generate the position setup part — produces positions JSON."""
        resolved_fixed = self.fixed if self.fixed else ""
        resolved_redesigned = self.redesigned if self.redesigned else ""

        # Determine input source for positions
        if resolved_fixed or resolved_redesigned:
            input_source = "selection"
            input_table = "-"
        else:
            input_source = "ligand"
            input_table = "-"

        fixed_param = resolved_fixed if resolved_fixed else "-"
        designed_param = resolved_redesigned if resolved_redesigned else "-"

        # The ligand `code` is read from the compounds stream inside the script.
        self.ligand_stream.save_json(self.ligand_json)

        with open(self.positions_args_json, "w") as f:
            json.dump({
                "structures_json": str(self.structures_json),
                "input_source": input_source,
                "input_table": input_table,
                "fixed_positions": fixed_param,
                "designed_positions": designed_param,
                "ligand": str(self.ligand_json),
                "design_within": self.design_within,
                "output_file": str(self.positions_json),
                "default_chain": self.chain,
            }, f, indent=2)

        return f"""echo "Setting up LigandMPNN position constraints"

# Create positions JSON for runtime lookup
echo "Computing position constraints..."
python {self.runtime_positions_py} "{self.positions_args_json}"

"""

    def _generate_script_run_ligandmpnn(self) -> str:
        """Generate the LigandMPNN execution part of the script."""
        # Build base LigandMPNN options
        base_options = f'--model_type "ligand_mpnn"'
        base_options += f' --checkpoint_ligand_mpnn "./model_params/ligandmpnn_{self.model}_25.pt"'
        base_options += f' --batch_size {self.num_sequences}'
        base_options += f' --number_of_batches {self.num_batches}'
        base_options += f' --ligand_mpnn_cutoff_for_score "{self.design_within}"'
        base_options += f' --out_folder "{self.lmpnn_out_folder}"'

        if self.temperature > 0:
            base_options += f' --temperature {self.temperature}'

        if self.bias_AA_per_residue:
            base_options += f' --bias_AA_per_residue "{self.bias_AA_per_residue}"'

        if self.seed > 0:
            base_options += f' --seed {self.seed}'

        return f"""echo "Executing LigandMPNN commands..."
cd {self.lmpnn_folder}
for struct_id in {Resolve.stream_ids(self.structures_json)}; do
    echo "Processing structure: $struct_id"
    PDB_FILE=$(resolve_stream_item "{self.structures_json}" "$struct_id")

    # Read position options from JSON
    POSITIONS=$(python "{self.resolve_positions_py}" "{self.positions_json}" "$struct_id")
    FIXED_OPTION=$(echo "$POSITIONS" | head -n1)
    REDESIGNED_OPTION=$(echo "$POSITIONS" | sed -n '2p')

    eval {self.container_prefix()}python run.py {base_options} --pdb_path '"$PDB_FILE"' $FIXED_OPTION $REDESIGNED_OPTION
done

"""

    def _generate_script_convert_outputs(self) -> str:
        """Generate the output conversion part of the script."""
        duplicates_flag = " --duplicates" if not self.remove_duplicates else ""
        fill_gaps_flag = f' --fill-gaps "{self.fill_gaps}"' if self.fill_gaps else ""
        step_tool_name = os.path.basename(self.output_folder)

        # Check for upstream missing table
        upstream_missing_path = self._get_upstream_missing_table_path(
            self.structures_stream
        )
        upstream_missing_flag = f' --upstream-missing "{upstream_missing_path}"' if upstream_missing_path else ""

        return f"""echo "Converting FASTA outputs to CSV format"
python {self.fa_to_csv_fasta_py} {self.seqs_folder} {self.queries_csv} {self.queries_fasta}{duplicates_flag}{fill_gaps_flag} --ds-json "{self.structures_json}" --missing-csv "{self.missing_csv}" --step-tool-name "{step_tool_name}"{upstream_missing_flag}

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after LigandMPNN execution."""
        # Expected FASTA files - one per input structure (keep compact/lazy patterns)
        fasta_ids = self.structures_stream.ids
        fasta_files = [os.path.join(self.seqs_folder, "<id>.fa")]

        # Predict sequence IDs (stream_id + sequence number)
        # Total sequences per structure = num_sequences (batch_size) * num_batches
        total_seqs = self.num_sequences * self.num_batches
        suffix_pattern = f"<1..{total_seqs}>"
        sequence_ids = generate_multiplied_ids_pattern(
            self.structures_stream.ids, suffix_pattern,
            input_stream_name="structures"
        )

        # Sequences stream - CSV-based with individual sequence IDs
        sequences = DataStream(
            name="sequences",
            ids=sequence_ids,
            files=[],
            map_table=self.queries_csv,
            format="csv"
        )

        # Fasta stream - file-based with structure IDs (one .fa file per structure)
        fasta = DataStream(
            name="fasta",
            ids=fasta_ids,
            files=fasta_files,
            format="fasta",
            metadata={
                "sequences_per_file": total_seqs,
                "contains_original": False,
            }
        )

        tables = {
            "sequences": TableInfo(
                name="sequences",
                path=self.queries_csv,
                columns=["id", "sequence", "sample", "T", "seed", "overall_confidence", "ligand_confidence", "seq_rec", "gaps"],
                description="LigandMPNN ligand-aware sequence generation results with binding scores"
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "kind", "cause"],
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
            "lmpnn_params": {
                "ligand_ids": list(self.ligand_stream.ids) if self.ligand_stream else [],
                "num_sequences": self.num_sequences,
                "num_batches": self.num_batches,
                "fixed": self.fixed,
                "redesigned": self.redesigned,
                "design_within": self.design_within,
                "chain": self.chain,
                "model": self.model,
                "remove_duplicates": self.remove_duplicates,
                "temperature": self.temperature,
                "bias_AA_per_residue": self.bias_AA_per_residue,
                "seed": self.seed,
            }
        })
        return base_dict
