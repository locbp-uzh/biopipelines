# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
LigandMPNN configuration for ligand-aware sequence design.
"""

import os
from typing import Dict, List, Any, Union, Tuple

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


class LigandMPNN(BaseConfig):
    """
    LigandMPNN configuration for ligand-aware sequence design.
    """

    TOOL_NAME = "LigandMPNN"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        repo_dir = folders.get("LigandMPNN", "")
        parent_dir = os.path.dirname(repo_dir)
        skip = "" if force_reinstall else f"""# Check if already installed
if [ -d "{repo_dir}" ] && {env_manager} env list 2>/dev/null | grep -q "ligandmpnn_env"; then
    echo "LigandMPNN already installed, skipping. Use force_reinstall=True to reinstall."
    exit 0
fi
"""
        return f"""echo "=== Installing LigandMPNN ==="
{skip}cd {parent_dir}
git clone https://github.com/dauparas/LigandMPNN.git
cd {repo_dir}
bash get_model_params.sh "./model_params"

{env_manager} create -n ligandmpnn_env python=3.11 -y
{env_manager} activate ligandmpnn_env
pip install -r requirements.txt

echo "=== LigandMPNN installation complete ==="
"""

    # Lazy path descriptors
    seqs_folder = Path(lambda self: os.path.join(self.output_folder, "seqs"))
    queries_csv = Path(lambda self: os.path.join(self.output_folder, f"queries.csv"))
    queries_fasta = Path(lambda self: os.path.join(self.output_folder, f"queries.fasta"))
    structures_json = Path(lambda self: os.path.join(self.output_folder, ".input_structures.json"))
    commands_file = Path(lambda self: os.path.join(self.output_folder, "lmpnn_commands.sh"))
    replacement_script = Path(lambda self: os.path.join(self.output_folder, "lmpnn_positions_replacement.sh"))
    id_map_json = Path(lambda self: os.path.join(self.output_folder, ".pdb_to_stream_id_map.json"))

    missing_csv = Path(lambda self: os.path.join(self.output_folder, "missing.csv"))

    # Helper script paths
    fa_to_csv_fasta_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_fa_to_csv_fasta.py"))
    lmpnn_folder = Path(lambda self: os.path.join(self.folders["data"], "LigandMPNN"))
    runtime_positions_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_lmpnn_runtime_positions.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 ligand: str = "",
                 num_sequences: int = 1,
                 fixed: Union[str, Tuple['TableInfo', str]] = "",
                 redesigned: Union[str, Tuple['TableInfo', str]] = "",
                 design_within: float = 5.0,
                 model: str = "v_32_010",
                 num_batches: int = 1,
                 remove_duplicates: bool = True,
                 **kwargs):
        """
        Initialize LigandMPNN configuration.

        Args:
            structures: Input structures as DataStream or StandardizedOutput
            ligand: Ligand identifier for binding site focus
            num_sequences: Number of sequences per batch (maps to --batch_size in LigandMPNN)
            fixed: Fixed positions. Accepts:
                   - PyMOL selection string: "10-20+30-40"
                   - Table column reference: (table, "column_name")
            redesigned: Designed positions. Accepts:
                   - PyMOL selection string: "10-20+30-40"
                   - Table column reference: (table, "column_name")
            design_within: Distance in Angstrom from ligand to redesign (fallback if positions not specified)
            model: LigandMPNN model version to use
            num_batches: Number of batches to run
            remove_duplicates: Remove duplicate sequences from output (default True)
        """
        # Resolve input to DataStream
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        # Store LigandMPNN-specific parameters
        self.ligand = ligand
        self.num_sequences = num_sequences
        self.fixed = fixed
        self.redesigned = redesigned
        self.design_within = design_within
        self.model = model
        self.num_batches = num_batches
        self.remove_duplicates = remove_duplicates

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate LigandMPNN-specific parameters."""
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")

        if not self.ligand:
            raise ValueError("ligand parameter is required")

        if self.num_sequences <= 0:
            raise ValueError("num_sequences must be positive")

        if self.num_batches <= 0:
            raise ValueError("num_batches must be positive")

        if self.design_within <= 0:
            raise ValueError("design_within must be positive")

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
            f"LIGAND: {self.ligand}",
            f"FIXED: {self.fixed or 'Auto (from table or ligand-based)'}",
            f"REDESIGNED: {self.redesigned or 'Auto (from table or ligand-based)'}",
            f"DESIGN WITHIN: {self.design_within}A",
            f"NUM SEQUENCES: {self.num_sequences}",
            f"NUM BATCHES: {self.num_batches}",
            f"MODEL: {self.model}"
        ])
        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate LigandMPNN execution script."""
        import json

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
        """Generate the position setup and script modification part."""
        # Resolve table references to DATASHEET_REFERENCE format
        resolved_fixed = self.resolve_table_reference(self.fixed) if self.fixed else ""
        resolved_redesigned = self.resolve_table_reference(self.redesigned) if self.redesigned else ""

        # Determine input source for positions
        if resolved_fixed or resolved_redesigned:
            input_source = "selection"
            input_table = "-"
        else:
            input_source = "ligand"
            input_table = "-"

        fixed_param = resolved_fixed if resolved_fixed else "-"
        designed_param = resolved_redesigned if resolved_redesigned else "-"

        # Build base LigandMPNN options
        base_options = f'--model_type "ligand_mpnn"'
        base_options += f' --checkpoint_ligand_mpnn "./model_params/ligandmpnn_{self.model}_25.pt"'
        base_options += f' --batch_size {self.num_sequences}'
        base_options += f' --number_of_batches {self.num_batches}'
        base_options += f' --ligand_mpnn_cutoff_for_score "{self.design_within}"'
        base_options += f' --out_folder "{self.output_folder}"'

        # Generate commands for each structure using DataStream IDs
        commands = []
        for struct_id, pdb_path in zip(self.structures_stream.ids, self.structures_stream.files):
            pdb_name = os.path.basename(pdb_path)
            commands.append(f'echo "Processing {pdb_name} with positions for ID: {struct_id}"')
            commands.append(f'python run.py {base_options} --pdb_path "{pdb_path}" {struct_id}_FIXED_OPTION_PLACEHOLDER {struct_id}_REDESIGNED_OPTION_PLACEHOLDER')
            commands.append("")

        # Write commands file at pipeline time (not SLURM time)
        os.makedirs(os.path.dirname(self.commands_file), exist_ok=True)
        with open(self.commands_file, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write(f"cd {self.lmpnn_folder}\n\n")
            f.write("# LigandMPNN commands with placeholders - will be replaced by position script\n")
            f.write(chr(10).join(commands))
            f.write("\n")

        return f"""echo "Setting up LigandMPNN position constraints"
# Commands file: {self.commands_file}

# Make the commands file executable
chmod +x {self.commands_file}

# Use existing HelpScript to create position replacement script
echo "Creating position replacement script..."
python {self.runtime_positions_py} "{self.structures_json}" "{input_source}" "{input_table}" "{fixed_param}" "{designed_param}" "{self.ligand}" "{self.design_within}" "{self.replacement_script}"

"""

    def _generate_script_run_ligandmpnn(self) -> str:
        """Generate the LigandMPNN execution part of the script."""
        return f"""# Run the replacement script on the commands file (not this script)
echo "Running position replacement script on commands file: {self.commands_file}"
bash {self.replacement_script} {self.commands_file}

# Now execute the modified commands file
echo "Executing LigandMPNN commands..."
bash {self.commands_file}

"""

    def _generate_script_convert_outputs(self) -> str:
        """Generate the output conversion part of the script."""
        duplicates_flag = " --duplicates" if not self.remove_duplicates else ""
        step_tool_name = os.path.basename(self.output_folder)

        # Check for upstream missing table
        upstream_missing_path = self._get_upstream_missing_table_path(
            self.structures_stream
        )
        upstream_missing_flag = f' --upstream-missing "{upstream_missing_path}"' if upstream_missing_path else ""

        return f"""echo "Converting FASTA outputs to CSV format"
python {self.fa_to_csv_fasta_py} {self.seqs_folder} {self.queries_csv} {self.queries_fasta}{duplicates_flag} --id-map {self.id_map_json} --missing-csv "{self.missing_csv}" --step-tool-name "{step_tool_name}"{upstream_missing_flag}

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after LigandMPNN execution."""
        # Expected FASTA files - one per input structure
        fasta_files = []
        fasta_ids = []

        for struct_id, pdb_path in zip(self.structures_stream.ids, self.structures_stream.files):
            pdb_base = os.path.splitext(os.path.basename(pdb_path))[0]
            fasta_path = os.path.join(self.seqs_folder, f"{pdb_base}.fa")
            fasta_files.append(fasta_path)
            fasta_ids.append(struct_id)

        # Predict sequence IDs (stream_id + sequence number)
        # Total sequences per structure = num_sequences (batch_size) * num_batches
        total_seqs = self.num_sequences * self.num_batches
        sequence_ids = []
        for struct_id in self.structures_stream.ids:
            for seq_num in range(1, total_seqs + 1):
                sequence_ids.append(f"{struct_id}_{seq_num}")

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
                columns=["id", "sequence", "sample", "T", "seed", "overall_confidence", "ligand_confidence", "seq_rec"],
                description="LigandMPNN ligand-aware sequence generation results with binding scores",
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
            "lmpnn_params": {
                "ligand": self.ligand,
                "num_sequences": self.num_sequences,
                "num_batches": self.num_batches,
                "fixed": self.fixed,
                "redesigned": self.redesigned,
                "design_within": self.design_within,
                "model": self.model,
                "remove_duplicates": self.remove_duplicates
            }
        })
        return base_dict
