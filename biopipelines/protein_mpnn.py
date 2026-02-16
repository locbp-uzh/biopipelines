# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
ProteinMPNN configuration for sequence design from protein structures.
"""

import os
from typing import Dict, List, Any, Union, Tuple

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
    from .combinatorics import generate_multiplied_ids
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream
    from combinatorics import generate_multiplied_ids


class ProteinMPNN(BaseConfig):
    """
    Configuration for ProteinMPNN sequence design.
    """

    TOOL_NAME = "ProteinMPNN"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        repo_dir = folders.get("ProteinMPNN", "")
        parent_dir = os.path.dirname(repo_dir)
        skip = "" if force_reinstall else f"""# Check if already installed
if [ -d "{repo_dir}" ]; then
    echo "ProteinMPNN already installed, skipping. Use force_reinstall=True to reinstall."
    exit 0
fi
"""
        # ProteinMPNN is git clone only — same for pip and conda
        return f"""echo "=== Installing ProteinMPNN ==="
{skip}mkdir -p {parent_dir}
cd {parent_dir}
git clone https://github.com/dauparas/ProteinMPNN.git

echo "=== ProteinMPNN installation complete ==="
"""

    # Lazy path descriptors
    parsed_pdbs_jsonl = Path(lambda self: os.path.join(self.output_folder, "parsed_pdbs.jsonl"))
    fixed_jsonl = Path(lambda self: os.path.join(self.output_folder, "fixed_pos.jsonl"))
    sele_csv = Path(lambda self: os.path.join(self.output_folder, "fixed_designed.csv"))
    seqs_folder = Path(lambda self: os.path.join(self.output_folder, "seqs"))
    main_table = Path(lambda self: os.path.join(self.output_folder, "proteinmpnn_results.csv"))
    queries_csv = Path(lambda self: os.path.join(self.output_folder, f"queries.csv"))
    queries_fasta = Path(lambda self: os.path.join(self.output_folder, f"queries.fasta"))
    structures_json = Path(lambda self: os.path.join(self.output_folder, ".input_structures.json"))
    id_map_json = Path(lambda self: os.path.join(self.output_folder, ".pdb_to_stream_id_map.json"))

    missing_csv = Path(lambda self: os.path.join(self.output_folder, "missing.csv"))

    # Helper scripts
    fixed_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_pmpnn_fixed_positions.py"))
    table_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_pmpnn_table.py"))
    fa_to_csv_fasta_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_fa_to_csv_fasta.py"))

    # ProteinMPNN installation scripts
    parse_py = Path(lambda self: os.path.join(self.folders["ProteinMPNN"], "helper_scripts", "parse_multiple_chains.py"))
    pmpnn_py = Path(lambda self: os.path.join(self.folders["ProteinMPNN"], "protein_mpnn_run.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 num_sequences: int = 1,
                 fixed: Union[str, Tuple['TableInfo', str]] = "",
                 redesigned: Union[str, Tuple['TableInfo', str]] = "",
                 fixed_chain: str = "A",
                 plddt_threshold: float = 100.0,
                 sampling_temp: float = 0.1,
                 model_name: str = "v_48_020",
                 soluble_model: bool = True,
                 remove_duplicates: bool = True,
                 fill_gaps: str = "G",
                 **kwargs):
        """
        Initialize ProteinMPNN configuration.

        Args:
            structures: Input structures as DataStream or StandardizedOutput
            num_sequences: Number of sequences to generate per structure
            fixed: Fixed positions. Accepts:
                   - PyMOL-style selection string: "10-20+30-40"
                   - Table column reference: (table, "column_name")
            redesigned: Designed positions. Accepts:
                   - PyMOL-style selection string: "10-20+30-40"
                   - Table column reference: (table, "column_name")
            fixed_chain: Chain to apply fixed positions to
            plddt_threshold: pLDDT threshold for automatic fixing (100 = no fixing)
            sampling_temp: Sampling temperature for sequence generation
            model_name: ProteinMPNN model variant
            soluble_model: Use soluble protein model
            remove_duplicates: Remove duplicate sequences from output (default True)
            fill_gaps: Amino acid to replace X (unknown/gap residues) with (default "G" for glycine).
                       Empty string means no filling (X is kept as-is).
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
        self.fixed_chain = fixed_chain
        self.plddt_threshold = plddt_threshold
        self.sampling_temp = sampling_temp
        self.model_name = model_name
        self.soluble_model = soluble_model
        self.remove_duplicates = remove_duplicates
        self.fill_gaps = fill_gaps

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate ProteinMPNN-specific parameters."""
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")

        if self.num_sequences <= 0:
            raise ValueError("num_sequences must be positive")

        if self.sampling_temp <= 0:
            raise ValueError("sampling_temp must be positive")

        valid_models = ["v_48_002", "v_48_010", "v_48_020", "v_48_030"]
        if self.model_name not in valid_models:
            raise ValueError(f"model_name must be one of: {valid_models}")

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
            f"FIXED CHAIN: {self.fixed_chain}",
            f"pLDDT THR: {self.plddt_threshold}",
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
        import json

        input_directory = os.path.dirname(self.structures_stream.files[0])

        # Serialize DataStream to JSON file (proper way to pass ids + files to HelpScript)
        datastream_dict = self.structures_stream.to_dict()
        with open(self.structures_json, 'w') as f:
            json.dump(datastream_dict, f, indent=2)

        # Write pdb_basename -> stream_id map for runtime ID remapping
        id_map = {}
        for struct_id, pdb_path in zip(self.structures_stream.ids, self.structures_stream.files):
            pdb_base = os.path.splitext(os.path.basename(pdb_path))[0]
            id_map[pdb_base] = struct_id
        with open(self.id_map_json, 'w') as f:
            json.dump(id_map, f, indent=2)

        # Resolve table references to DATASHEET_REFERENCE format
        resolved_fixed = self.resolve_table_reference(self.fixed) if self.fixed else ""
        resolved_redesigned = self.resolve_table_reference(self.redesigned) if self.redesigned else ""

        # Determine input source for fixed positions
        if resolved_fixed or resolved_redesigned:
            input_source = "selection"
        else:
            input_source = "plddt"

        fixed_param = resolved_fixed if resolved_fixed else "-"
        designed_param = resolved_redesigned if resolved_redesigned else "-"

        return f"""echo "Determining fixed positions"
python {self.fixed_py} "{self.structures_json}" "{input_source}" "-" {self.plddt_threshold} "{fixed_param}" "{designed_param}" "{self.fixed_chain}" "{self.fixed_jsonl}" "{self.sele_csv}"

echo "Parsing multiple PDBs"
python {self.parse_py} --input_path {input_directory} --output_path {self.parsed_pdbs_jsonl}

"""

    def _generate_script_run_proteinmpnn(self) -> str:
        """Generate the ProteinMPNN execution part of the script."""
        pmpnn_options = f"--num_seq_per_target {self.num_sequences}"
        pmpnn_options += f" --sampling_temp {self.sampling_temp}"
        pmpnn_options += f" --model_name {self.model_name}"

        if self.soluble_model:
            pmpnn_options += " --use_soluble_model"

        return f"""echo "Running model"
echo "Options: {pmpnn_options}"
echo "Output folder: {self.output_folder}"

python {self.pmpnn_py} --jsonl_path {self.parsed_pdbs_jsonl} --fixed_positions_jsonl {self.fixed_jsonl} --out_folder {self.output_folder} {pmpnn_options}

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
python {self.fa_to_csv_fasta_py} {self.seqs_folder} {self.queries_csv} {self.queries_fasta} --id-map {self.id_map_json}{duplicates_flag}{fill_gaps_flag} --missing-csv "{self.missing_csv}" --step-tool-name "{step_tool_name}"{upstream_missing_flag}

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after ProteinMPNN execution."""
        # Expected FASTA files - one per input structure
        fasta_files = []
        fasta_ids = []

        for struct_id, pdb_path in zip(self.structures_stream.ids, self.structures_stream.files):
            pdb_base = os.path.splitext(os.path.basename(pdb_path))[0]
            fasta_path = os.path.join(self.seqs_folder, f"{pdb_base}.fa")
            fasta_files.append(fasta_path)
            fasta_ids.append(struct_id)

        # Predict sequence IDs (stream_id + sequence number)
        suffixes = [str(i) for i in range(1, self.num_sequences + 1)]
        sequence_ids, provenance = generate_multiplied_ids(
            self.structures_stream.ids, suffixes,
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
            format="fasta"
        )

        tables = {
            "sequences": TableInfo(
                name="sequences",
                path=self.queries_csv,
                columns=["id", "structures.id", "source_pdb", "sequence", "score", "seq_recovery", "gaps"],
                description="ProteinMPNN sequence results",
                count=len(sequence_ids)
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "cause"],
                description="IDs removed (duplicates or upstream) with removal reason",
                count="variable"
            )
        }

        return {
            "structures": DataStream.empty("structures", "pdb"),
            "sequences": sequences,
            "fasta": fasta,
            "compounds": DataStream.empty("compounds", "sdf"),
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
                "fixed_chain": self.fixed_chain,
                "plddt_threshold": self.plddt_threshold,
                "sampling_temp": self.sampling_temp,
                "model_name": self.model_name,
                "soluble_model": self.soluble_model,
                "remove_duplicates": self.remove_duplicates
            }
        })
        return base_dict
