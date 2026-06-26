# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
RBSDesigner tool for designing synthetic ribosome binding sites.

Takes DNA sequences (e.g. from DNAEncoder) and designs RBS sequences using the
Salis thermodynamic model (Salis, Mirsky & Voigt, Nature Biotechnology 2009)
with simulated annealing. Outputs full genes ready for synthesis:
pre_sequence + rbs_sequence + cds_dna.
"""

import os
from typing import Dict, List, Any, Union

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


# TIR presets mapping
TIR_PRESETS = {
    "low": 100,
    "medium": 1000,
    "high": 10000,
    "maximum": 100000,
}


class RBSDesigner(BaseConfig):
    """
    Pipeline tool for designing synthetic ribosome binding sites (RBS).

    Uses the Salis thermodynamic model to predict translation initiation rates
    and designs RBS sequences via simulated annealing to match a target TIR.
    Requires ViennaRNA for RNA free energy calculations.

    Citation: Salis, Mirsky & Voigt, Nature Biotechnology 2009.
    """

    TOOL_NAME = "RBSDesigner"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        biopipelines = folders.get("biopipelines", "")
        env_check = cls._env_exists_check("rbs_designer", env_manager)
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_check}; then
    echo "RBSDesigner environment already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block("rbs_designer", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("rbs_designer", env_manager, biopipelines)
        return f"""echo "=== Installing RBSDesigner ==="
{skip}echo "Creating rbs_designer environment with ViennaRNA from bioconda..."
{remove_block}
{env_block}
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create rbs_designer environment."
    echo "ViennaRNA requires the bioconda channel with flexible channel priority."
    exit 1
fi

# Verify installation
if {env_manager} run -n rbs_designer python -c "import RNA" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== RBSDesigner installation complete ==="
else
    echo "ERROR: RBSDesigner verification failed (cannot import RNA / ViennaRNA bindings)"
    exit 1
fi
"""

    # Lazy path descriptors
    rbs_csv = Path(lambda self: self.table_path("rbs"))
    info_txt = Path(lambda self: os.path.join(self.extras_folder, "rbs_info.txt"))
    config_file = Path(lambda self: self.configuration_path("rbs_designer_config.json"))
    sequences_json = Path(lambda self: self.configuration_path(".input_sequences.json"))
    sequences_csv_path = Path(lambda self: self.configuration_path(".input_sequences.csv"))
    designer_py = Path(lambda self: self.pipe_script_path("pipe_rbs_designer.py"))

    def __init__(self,
                 sequences: Union[DataStream, StandardizedOutput],
                 tir: Union[str, int, float] = "medium",
                 pre_sequence: str = "",
                 add_start_codon: bool = False,
                 **kwargs):
        """
        Initialize RBS designer tool.

        Args:
            sequences: Input DNA sequences as DataStream or StandardizedOutput
                       (typically from DNAEncoder).
            tir: Target translation initiation rate. Either a numeric value or
                 a preset string: "low" (100), "medium" (1000), "high" (10000),
                 "maximum" (100000).
            pre_sequence: Optional fixed 5'UTR DNA sequence to prepend before
                          the designed RBS in the full gene output.
            add_start_codon: If True, prepend ATG to each input sequence before
                             RBS design. Use when sequences do not already begin
                             with a start codon.
            **kwargs: Additional parameters.

        Output:
            Streams: sequences (.dna)
            Tables:
                rbs: id | sequence | rbs_sequence | full_gene | converged | dg_total | tir_predicted | target_tir | target_dg | min_achievable_dg | spacing | dg_mrna_rrna | dg_start | dg_spacing | dg_mrna | dg_standby
        """
        # Resolve input to DataStream
        if isinstance(sequences, StandardizedOutput):
            self.sequences_stream: DataStream = sequences.streams.sequences
        elif isinstance(sequences, DataStream):
            self.sequences_stream = sequences
        else:
            raise ValueError(f"sequences must be DataStream or StandardizedOutput, got {type(sequences)}")

        # Resolve TIR preset or numeric value
        if isinstance(tir, str):
            tir_lower = tir.lower()
            if tir_lower in TIR_PRESETS:
                self.tir = TIR_PRESETS[tir_lower]
                self.tir_preset = tir_lower
            else:
                raise ValueError(
                    f"Invalid TIR preset '{tir}'. Must be one of: {list(TIR_PRESETS.keys())} "
                    f"or a numeric value."
                )
        else:
            self.tir = float(tir)
            self.tir_preset = None

        self.pre_sequence = pre_sequence.upper()
        self.add_start_codon = add_start_codon

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate RBSDesigner parameters."""
        if not self.sequences_stream or len(self.sequences_stream) == 0:
            raise ValueError("sequences parameter is required and must not be empty")

        if self.tir <= 0:
            raise ValueError(f"TIR must be positive, got {self.tir}")

        # Validate pre_sequence is valid DNA
        if self.pre_sequence:
            valid_bases = set("ATCG")
            invalid = set(self.pre_sequence) - valid_bases
            if invalid:
                raise ValueError(
                    f"pre_sequence contains invalid characters: {invalid}. "
                    f"Must be DNA (A, T, C, G only)."
                )

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input sequences."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        tir_display = f"{self.tir:.0f}"
        if self.tir_preset:
            tir_display += f" ({self.tir_preset})"

        config_lines.extend([
            f"SEQUENCES: {len(self.sequences_stream)} sequences",
            f"TARGET TIR: {tir_display}",
        ])

        if self.pre_sequence:
            config_lines.append(f"PRE_SEQUENCE: {self.pre_sequence}")

        if self.add_start_codon:
            config_lines.append("ADD_START_CODON: True (ATG prepended to each sequence)")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate script to perform RBS design."""
        self.sequences_stream.save_json(self.sequences_json)

        script_content = "#!/bin/bash\n"
        script_content += "# RBSDesigner execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.generate_filtered_map_table_block(
            self.sequences_json, self.sequences_csv_path, required_columns=["id", "sequence"]
        )
        script_content += self.activate_environment()
        script_content += self._generate_script_run_design()
        script_content += self.generate_completion_check_footer()

        return script_content

    def _generate_script_run_design(self) -> str:
        """Generate the RBS design part of the script."""
        import json

        config_data = {
            "sequences_csv": self.sequences_csv_path,
            "tir": self.tir,
            "pre_sequence": self.pre_sequence,
            "add_start_codon": self.add_start_codon,
            "rbs_output": self.rbs_csv,
            "info_output": self.info_txt,
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        return f"""echo "Designing RBS sequences"
echo "Target TIR: {self.tir}"
echo "Output folder: {self.output_folder}"

python "{self.designer_py}" --config "{self.config_file}"

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after RBS design."""
        sequence_ids = list(self.sequences_stream.ids)

        sequences = DataStream(
            name="sequences",
            ids=sequence_ids,
            files=[],
            map_table=self.rbs_csv,
            format="dna"
        )

        tables = {
            "rbs": TableInfo(
                name="rbs",
                path=self.rbs_csv,
                columns=[
                    "id", "sequence", "rbs_sequence", "full_gene",
                    "converged", "dg_total", "tir_predicted", "target_tir", "target_dg",
                    "min_achievable_dg", "spacing", "dg_mrna_rrna", "dg_start", "dg_spacing",
                    "dg_mrna", "dg_standby",
                ],
                description="RBS design results with thermodynamic parameters"
            )
        }

        return {
            "sequences": sequences,
            "tables": tables,
            "output_folder": self.output_folder,
            "info": self.info_txt,
        }
